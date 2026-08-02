#!/usr/bin/env python3
"""Classify generic six-stick polygons exactly using Calvo's invariants.

For an embedded hexagon H=(v1,...,v6), Calvo defines Delta_2, Delta_4,
and Delta_6 as algebraic intersections with the open triangular discs

    (v1,v2,v3), (v3,v4,v5), (v5,v6,v1).

The hexagon is a trefoil exactly when the three intersection numbers are all
+1 or all -1; it is an unknot exactly when at least one is zero.  Every
segment--open-triangle intersection is a Radon-partition predicate on five
points, so the test depends only on the 15 tetrahedral orientation signs of
the six vertices.  Conditional on Calvo's cited topological theorem, it gives
an exact finite classifier for every generic N=6 order-type bucket, with no
knot catalogue or pyknotid call.

The default input is the tracked 500,000-sample reference bucket table.  A stringent
reproduction command is:

    python3 scripts/repro/calvo_n6_exact.py \
      --expected-trefoils 1856 --validate-standard-direct \
      --audit-signs

The optional sign audit replays the exact seeded random-number stream, including
the reservoir-sampling draws that affect all later samples.  It checks every
bucket count, records the minimum scale-normalized Pluecker determinant, and
recomputes the nearest cases with 100-digit ``decimal`` arithmetic.  The full
500,000-sample audit takes roughly one to two minutes on a laptop; without
``--audit-signs`` the default CSV reclassification remains fast.

This certifies classification of the Monte Carlo sample; it does not turn the
Monte Carlo rate into an exact Haar chamber volume.
"""

from __future__ import annotations

import argparse
import csv
import heapq
import itertools
import math
import sys
import time
from collections import Counter
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
REFERENCE_DIR = ROOT / "data" / "reference" / "n6"
QUADS = tuple(itertools.combinations(range(6), 4))
QUAD_INDEX = {quad: index for index, quad in enumerate(QUADS)}


class DegenerateChirotope(ValueError):
    """Raised when an orientation sign is zero or absent."""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bucket-csv",
        type=Path,
        default=REFERENCE_DIR / "bucket_volume.csv",
        help="tracked N=6 order-type bucket table",
    )
    parser.add_argument(
        "--expected-trefoils",
        type=int,
        help="fail unless the exact Calvo total equals this count",
    )
    parser.add_argument(
        "--validate-standard-direct",
        action="store_true",
        help="reconstruct the saved 1000-sample Haar and Gaussian runs and compare every pyknotid label",
    )
    parser.add_argument(
        "--audit-signs",
        action="store_true",
        help="replay the seeded N=6 signature run, compare every bucket, and verify the nearest walls at high precision",
    )
    parser.add_argument(
        "--audit-samples",
        type=int,
        help="samples in the sign audit (default: infer the full total from --bucket-csv)",
    )
    parser.add_argument(
        "--audit-base-seed",
        type=int,
        default=20260604,
        help="base seed used by order_type_grouped_volume.py",
    )
    parser.add_argument(
        "--audit-reservoir-size",
        type=int,
        default=3,
        help="checks-per-group value whose RNG draws must be replayed exactly",
    )
    parser.add_argument(
        "--audit-nearest",
        type=int,
        default=100,
        help="number of nearest-wall samples to recompute with high precision",
    )
    parser.add_argument(
        "--audit-precision",
        type=int,
        default=100,
        help="decimal digits used to recompute nearest-wall determinants",
    )
    parser.add_argument(
        "--orientation-tolerance",
        type=float,
        default=1e-10,
        help="absolute 3x3 determinant tolerance used by the original signature generator",
    )
    return parser


def permutation_parity(values: tuple[int, int, int, int]) -> int:
    inversions = sum(
        values[left] > values[right]
        for left in range(4)
        for right in range(left + 1, 4)
    )
    return -1 if inversions % 2 else 1


def canonical_global(signature: str) -> str:
    flipped = "".join("-" if char == "+" else "+" if char == "-" else char for char in signature)
    return min(signature, flipped)


def orientation_sign(signature: str, a: int, b: int, c: int, d: int) -> int:
    """Return the alternating chirotope value chi(a,b,c,d)."""

    vertices = (a, b, c, d)
    char = signature[QUAD_INDEX[tuple(sorted(vertices))]]
    if char == "0":
        raise DegenerateChirotope(f"zero orientation on {tuple(sorted(vertices))}")
    if char not in "+-":
        raise DegenerateChirotope(f"invalid orientation character {char!r}")
    sorted_sign = 1 if char == "+" else -1
    return sorted_sign * permutation_parity(vertices)


def triangle_piercing_sign(
    signature: str,
    triangle: tuple[int, int, int],
    segment: tuple[int, int],
) -> int:
    """Return the oriented piercing sign, or zero when the open sets miss.

    For points (a,b,c,p,q), the signed maximal minors below are the
    coefficients of their unique affine dependence.  The open segment pq
    meets the open triangle abc exactly when the first three coefficients
    have one sign and the final two have the opposite sign.
    """

    a, b, c = triangle
    p, q = segment
    coefficients = (
        orientation_sign(signature, b, c, p, q),
        -orientation_sign(signature, a, c, p, q),
        orientation_sign(signature, a, b, p, q),
        -orientation_sign(signature, a, b, c, q),
        orientation_sign(signature, a, b, c, p),
    )
    triangle_same = coefficients[0] == coefficients[1] == coefficients[2]
    segment_same = coefficients[3] == coefficients[4]
    opposite_sides = coefficients[0] == -coefficients[3]
    if not (triangle_same and segment_same and opposite_sides):
        return 0

    # With the triangle oriented (a,b,c) and the segment directed p -> q,
    # this is sign(((b-a) x (c-a)) dot (q-p)).
    return orientation_sign(signature, a, b, c, q)


def calvo_deltas(signature: str) -> tuple[int, int, int]:
    if len(signature) != len(QUADS):
        raise DegenerateChirotope(f"expected {len(QUADS)} signs, got {len(signature)}")
    deltas: list[int] = []
    for start in (0, 2, 4):
        triangle = (start % 6, (start + 1) % 6, (start + 2) % 6)
        first_segment = ((start + 3) % 6, (start + 4) % 6)
        second_segment = ((start + 4) % 6, (start + 5) % 6)
        delta = triangle_piercing_sign(signature, triangle, first_segment)
        delta += triangle_piercing_sign(signature, triangle, second_segment)
        deltas.append(delta)
    return deltas[0], deltas[1], deltas[2]


def calvo_label(signature: str) -> tuple[str, tuple[int, int, int]]:
    deltas = calvo_deltas(signature)
    if deltas == (1, 1, 1) or deltas == (-1, -1, -1):
        return "3_1", deltas
    if 0 in deltas:
        return "0_1", deltas
    return "inconsistent_or_nongeneric", deltas


def signature_from_vertices(vertices: np.ndarray, tolerance: float = 1e-12) -> str:
    points = np.asarray(vertices, dtype=float)
    signs: list[str] = []
    for a, b, c, d in QUADS:
        determinant = float(
            np.linalg.det(
                np.column_stack(
                    (points[b] - points[a], points[c] - points[a], points[d] - points[a])
                )
            )
        )
        if abs(determinant) <= tolerance:
            signs.append("0")
        else:
            signs.append("+" if determinant > 0.0 else "-")
    return "".join(signs)


def classify_bucket_table(path: Path) -> tuple[int, int, Counter[tuple[int, int, int]], list[tuple[str, int]]]:
    total_samples = 0
    trefoil_samples = 0
    delta_counts: Counter[tuple[int, int, int]] = Counter()
    trefoil_buckets: list[tuple[str, int]] = []
    inconsistent: list[tuple[str, tuple[int, int, int]]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            signature = row["signature"]
            count = int(row["count"])
            label, deltas = calvo_label(signature)
            total_samples += count
            delta_counts[deltas] += count
            if label == "3_1":
                trefoil_samples += count
                trefoil_buckets.append((signature, count))
            elif label != "0_1":
                inconsistent.append((signature, deltas))
    if inconsistent:
        raise RuntimeError(f"Calvo's alternatives failed on {len(inconsistent)} buckets: {inconsistent[:3]}")
    return total_samples, trefoil_samples, delta_counts, trefoil_buckets


def read_bucket_counts(path: Path) -> dict[str, int]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["signature"]: int(row["count"]) for row in csv.DictReader(handle)}


def dual_quad_data() -> tuple[tuple[tuple[int, int], ...], np.ndarray]:
    """Return complementary pairs and Hodge-star parities for the 15 signs."""

    pairs: list[tuple[int, int]] = []
    parities: list[int] = []
    for quad in QUADS:
        complement = tuple(vertex for vertex in range(6) if vertex not in quad)
        sequence = quad + complement
        inversions = sum(
            sequence[left] > sequence[right]
            for left in range(6)
            for right in range(left + 1, 6)
        )
        pairs.append((complement[0], complement[1]))
        parities.append(-1 if inversions % 2 else 1)
    return tuple(pairs), np.asarray(parities, dtype=float)


def dual_signature_and_margin(
    gale_vectors: np.ndarray,
    first_indices: np.ndarray,
    second_indices: np.ndarray,
    parities: np.ndarray,
    orientation_tolerance: float,
) -> tuple[str, float, int, np.ndarray]:
    """Compute the affine chirotope through its rank-two Gale dual.

    The six-by-two Gale matrix need not have orthonormal columns.  All its
    two-by-two minors change by the same determinant under a basis change, so

        min_J |p_J| / sqrt(sum_J p_J^2)

    is a basis-free, scale-free distance-to-wall diagnostic.  For the centered
    row-orthonormal output frame used by the original pipeline, the affine
    three-by-three determinant is ``sqrt(6) * p_J / ||p||`` up to sign.
    """

    minors = (
        gale_vectors[first_indices, 0] * gale_vectors[second_indices, 1]
        - gale_vectors[first_indices, 1] * gale_vectors[second_indices, 0]
    )
    norm = math.sqrt(float(minors @ minors))
    if not math.isfinite(norm) or norm <= 0.0:
        raise RuntimeError("rank-deficient or non-finite Gale sample in sign audit")
    normalized = minors / norm
    closest_index = int(np.argmin(np.abs(normalized)))
    margin = float(abs(normalized[closest_index]))
    normalized_tolerance = orientation_tolerance / math.sqrt(6.0)
    if margin <= normalized_tolerance:
        chars = [
            "0"
            if abs(float(value)) <= normalized_tolerance
            else "+"
            if value * parity > 0.0
            else "-"
            for value, parity in zip(normalized, parities)
        ]
    else:
        chars = ["+" if value * parity > 0.0 else "-" for value, parity in zip(normalized, parities)]
    return canonical_global("".join(chars)), margin, closest_index, minors


def audit_sign_generation(
    path: Path,
    samples: int | None,
    base_seed: int,
    reservoir_size: int,
    nearest_count: int,
    precision: int,
    orientation_tolerance: float,
    expected_trefoils: int | None,
) -> None:
    """Replay and numerically certify the saved N=6 sign experiment."""

    if samples is not None and samples < 1:
        raise SystemExit("--audit-samples must be positive")
    if reservoir_size < 0:
        raise SystemExit("--audit-reservoir-size must be nonnegative")
    if nearest_count < 1:
        raise SystemExit("--audit-nearest must be positive")
    if precision < 30:
        raise SystemExit("--audit-precision must be at least 30 digits")
    if orientation_tolerance < 0.0:
        raise SystemExit("--orientation-tolerance must be nonnegative")

    sys.path.insert(0, str(ROOT))
    from gaussian_knots.general_chamber_geometry import (
        build_chamber_model,
        canonical_global as helper_canonical_global,
        d_orbit_representative,
        orientation_signature as helper_orientation_signature,
        vertices_from_kernel as helper_vertices_from_kernel,
    )

    # Guard against a future divergence between the local and shared helpers.
    if helper_canonical_global("-+") != canonical_global("-+"):
        raise RuntimeError("canonical-global conventions disagree")

    saved_counts = read_bucket_counts(path)
    saved_total = sum(saved_counts.values())
    audit_samples = saved_total if samples is None else samples
    model = build_chamber_model(6)
    pairs, parities = dual_quad_data()
    first_indices = np.asarray([pair[0] for pair in pairs], dtype=int)
    second_indices = np.asarray([pair[1] for pair in pairs], dtype=int)
    effective_seed = base_seed + 1009 * 6
    rng = np.random.default_rng(effective_seed)
    replay_counts: Counter[str] = Counter()
    canonical_cache: dict[str, str] = {}
    nearest: list[tuple[float, int, int, np.ndarray, str, np.ndarray]] = []
    started = time.perf_counter()

    for sample_index in range(audit_samples):
        # This is exactly the RNG-consuming matrix in random_kernel.  QR and
        # SVD consume no randomness, and duality lets us skip them in the main
        # loop without changing a sign.
        raw_kernel = rng.normal(size=(model.h_dimension, model.kernel_dimension))
        gale_vectors = model.simplex_vertices @ raw_kernel
        signature, margin, quad_index, minors = dual_signature_and_margin(
            gale_vectors,
            first_indices,
            second_indices,
            parities,
            orientation_tolerance,
        )
        canonical = canonical_cache.get(signature)
        if canonical is None:
            canonical = d_orbit_representative(
                signature,
                6,
                model.quads,
                model.quad_index,
            )
            canonical_cache[signature] = canonical
        replay_counts[canonical] += 1

        # order_type_grouped_volume.py uses the same RNG for reservoir draws.
        # Omitting these calls changes every later Gaussian sample.
        bucket_count = replay_counts[canonical]
        if reservoir_size > 0 and bucket_count > reservoir_size:
            rng.integers(bucket_count)

        candidate = (
            -margin,
            sample_index,
            quad_index,
            raw_kernel.copy(),
            canonical,
            minors.copy(),
        )
        if len(nearest) < nearest_count:
            heapq.heappush(nearest, candidate)
        elif margin < -nearest[0][0]:
            heapq.heapreplace(nearest, candidate)

        if (sample_index + 1) % 100000 == 0:
            current_minimum = min(-item[0] for item in nearest)
            elapsed = time.perf_counter() - started
            print(
                f"  audit progress={sample_index + 1} buckets={len(replay_counts)} "
                f"min_margin={current_minimum:.17g} elapsed={elapsed:.2f}s",
                flush=True,
            )

    all_signatures = set(saved_counts) | set(replay_counts)
    differences = [
        (signature, replay_counts.get(signature, 0), saved_counts.get(signature, 0))
        for signature in sorted(all_signatures)
        if replay_counts.get(signature, 0) != saved_counts.get(signature, 0)
    ]
    compare_full_table = audit_samples == saved_total
    if compare_full_table and differences:
        raise RuntimeError(f"sign-audit bucket mismatch: {differences[:5]}")

    zero_signatures = [signature for signature in replay_counts if "0" in signature]
    if zero_signatures:
        raise RuntimeError(
            f"the numerical tolerance produced {len(zero_signatures)} degenerate signatures"
        )
    replay_trefoils = sum(
        count
        for signature, count in replay_counts.items()
        if calvo_label(signature)[0] == "3_1"
    )
    if expected_trefoils is not None and compare_full_table and replay_trefoils != expected_trefoils:
        raise RuntimeError(
            f"sign audit expected {expected_trefoils} trefoils, found {replay_trefoils}"
        )

    nearest_rows = verify_nearest_signs(
        nearest,
        model,
        pairs,
        precision,
        orientation_tolerance,
        d_orbit_representative,
        helper_orientation_signature,
        helper_vertices_from_kernel,
    )
    closest = nearest_rows[0]
    elapsed = time.perf_counter() - started
    print("SIGN-GENERATION AUDIT")
    print(f"effective seed: {effective_seed} (= {base_seed} + 1009*6)")
    print(f"samples replayed: {audit_samples}")
    print(f"replayed buckets: {len(replay_counts)}")
    if compare_full_table:
        print(f"saved buckets: {len(saved_counts)}; bucket-count differences: 0")
    else:
        print(f"saved-table comparison skipped: audit is a {audit_samples}-sample prefix of {saved_total}")
    print(f"Calvo trefoils in replay: {replay_trefoils}")
    print(
        "minimum normalized Pluecker margin: "
        f"{closest['decimal_margin']:.17g} at zero-based sample "
        f"{closest['sample_index']}, quad {closest['quad']}"
    )
    print(
        "corresponding absolute affine determinant: "
        f"{closest['affine_determinant']:.17g} "
        f"({closest['tolerance_ratio']:.6g} times the configured tolerance)"
    )
    print(
        f"nearest {len(nearest_rows)}: high-precision sign mismatches="
        f"{sum(int(row['high_precision_sign_mismatches']) for row in nearest_rows)}, "
        f"original-pipeline signature mismatches="
        f"{sum(int(not row['pipeline_match']) for row in nearest_rows)}"
    )
    print(
        "closest determinant / [4 eps (|ad|+|bc|)] = "
        f"{closest['error_proxy_ratio']:.6g}"
    )
    print(f"audit elapsed: {elapsed:.2f}s")


def verify_nearest_signs(
    nearest: list[tuple[float, int, int, np.ndarray, str, np.ndarray]],
    model: object,
    pairs: tuple[tuple[int, int], ...],
    precision: int,
    orientation_tolerance: float,
    d_orbit_representative: object,
    helper_orientation_signature: object,
    helper_vertices_from_kernel: object,
) -> list[dict[str, object]]:
    """Recompute the nearest Pluecker coordinates from the original floats."""

    ordered = sorted(nearest, key=lambda item: -item[0])
    rows: list[dict[str, object]] = []
    with localcontext() as context:
        context.prec = precision
        decimal_basis = [
            [Decimal.from_float(float(value)) for value in row]
            for row in model.simplex_vertices
        ]
        for negative_margin, sample_index, quad_index, raw_kernel, canonical, float_minors in ordered:
            decimal_kernel = [
                [Decimal.from_float(float(value)) for value in row]
                for row in raw_kernel
            ]
            gale = [
                [
                    sum(
                        (decimal_basis[row][inner] * decimal_kernel[inner][column] for inner in range(5)),
                        Decimal(0),
                    )
                    for column in range(2)
                ]
                for row in range(6)
            ]
            decimal_minors = [
                gale[first][0] * gale[second][1]
                - gale[first][1] * gale[second][0]
                for first, second in pairs
            ]
            decimal_norm = sum((value * value for value in decimal_minors), Decimal(0)).sqrt()
            decimal_margin = min(abs(value) for value in decimal_minors) / decimal_norm
            sign_mismatches = sum(
                (float_value > 0.0) != (decimal_value > 0)
                for float_value, decimal_value in zip(float_minors, decimal_minors)
            )

            first, second = pairs[quad_index]
            left_product = gale[first][0] * gale[second][1]
            right_product = gale[first][1] * gale[second][0]
            cancellation_denominator = abs(left_product) + abs(right_product)
            cancellation_ratio = (
                abs(decimal_minors[quad_index]) / cancellation_denominator
                if cancellation_denominator
                else Decimal("Infinity")
            )

            q_matrix, r_matrix = np.linalg.qr(raw_kernel, mode="reduced")
            for column in range(model.kernel_dimension):
                if r_matrix[column, column] < 0.0:
                    q_matrix[:, column] *= -1.0
            vertices = helper_vertices_from_kernel(model, q_matrix)
            pipeline_signature = helper_orientation_signature(
                vertices,
                model.quads,
                tolerance=orientation_tolerance,
            )
            pipeline_canonical = d_orbit_representative(
                pipeline_signature,
                6,
                model.quads,
                model.quad_index,
            )

            double_gale = model.simplex_vertices @ raw_kernel
            ad = double_gale[first, 0] * double_gale[second, 1]
            bc = double_gale[first, 1] * double_gale[second, 0]
            error_proxy = 4.0 * np.finfo(float).eps * (abs(ad) + abs(bc))
            error_proxy_ratio = (
                abs(float_minors[quad_index]) / error_proxy
                if error_proxy > 0.0
                else math.inf
            )
            affine_determinant = math.sqrt(6.0) * float(decimal_margin)
            rows.append(
                {
                    "sample_index": sample_index,
                    "quad": QUADS[quad_index],
                    "float_margin": -negative_margin,
                    "decimal_margin": float(decimal_margin),
                    "affine_determinant": affine_determinant,
                    "tolerance_ratio": affine_determinant / orientation_tolerance
                    if orientation_tolerance > 0.0
                    else math.inf,
                    "cancellation_ratio": float(cancellation_ratio),
                    "error_proxy_ratio": error_proxy_ratio,
                    "high_precision_sign_mismatches": sign_mismatches,
                    "pipeline_match": pipeline_canonical == canonical,
                }
            )
    return rows


def validate_direct_csv(path: Path) -> tuple[int, int]:
    sys.path.insert(0, str(ROOT))
    from gaussian_knots.generation import projected_simplex_polygon

    checked = 0
    mismatches = 0
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if int(row["N"]) != 6:
                continue
            seed = int(row["sample_seed"])
            model = row["projection_model"]
            vertices = projected_simplex_polygon(6, np.random.default_rng(seed), projection_model=model)
            label, _ = calvo_label(signature_from_vertices(vertices))
            exact_nontrivial = label == "3_1"
            recorded = row["is_nontrivial"].strip().lower()
            recorded_nontrivial = recorded == "true"
            checked += 1
            mismatches += exact_nontrivial != recorded_nontrivial
    return checked, mismatches


def main() -> int:
    args = build_parser().parse_args()
    total, trefoils, delta_counts, trefoil_buckets = classify_bucket_table(args.bucket_csv)
    print("CALVO EXACT N=6 CLASSIFICATION")
    print(f"bucket table: {args.bucket_csv}")
    print(f"samples: {total}")
    print(f"trefoils: {trefoils}")
    print(f"trefoil rate: {trefoils / total:.9f}")
    print(f"trefoil buckets: {trefoil_buckets}")
    print("Calvo-delta sample masses:")
    for deltas, count in delta_counts.most_common():
        print(f"  {deltas}: {count}")

    if args.expected_trefoils is not None and trefoils != args.expected_trefoils:
        raise SystemExit(f"expected {args.expected_trefoils} trefoils, found {trefoils}")

    if args.validate_standard_direct:
        for model in ("haar", "gaussian"):
            path = REFERENCE_DIR / f"{model}_samples_N6.csv"
            checked, mismatches = validate_direct_csv(path)
            print(f"direct {model}: checked={checked}, Calvo/pyknotid mismatches={mismatches}")
            if mismatches:
                raise SystemExit(f"Calvo disagrees with {mismatches} rows in {path}")
    if args.audit_signs:
        audit_sign_generation(
            path=args.bucket_csv,
            samples=args.audit_samples,
            base_seed=args.audit_base_seed,
            reservoir_size=args.audit_reservoir_size,
            nearest_count=args.audit_nearest,
            precision=args.audit_precision,
            orientation_tolerance=args.orientation_tolerance,
            expected_trefoils=args.expected_trefoils,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
