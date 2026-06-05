#!/usr/bin/env python3
"""Targeted N=6 Grassmannian probe for the projected-simplex knot model.

For N=6, H = 1^perp has dimension 5.  Up to target rotations, a projection
H -> R^3 is determined by its kernel K in Gr_2(H).  The nine non-adjacent edge
pairs give determinant walls

    det[a, b, k_1, k_2, c] = 0,

where a and b are edge directions, c is the displacement between edge starts,
and k_1,k_2 span K.  This script samples K, classifies the resulting 6-stick
polygon, and records the determinant-wall sign signature.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from statistics import NormalDist

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
GAUSSIAN_KNOTS = ROOT
sys.path.insert(0, str(ROOT))

from gaussian_knots.generation import is_numerically_embedded
from gaussian_knots.pyknotid_adapter import identify_polygon, inspect_pyknotid_environment


VERTEX_COUNT = 6
H_DIMENSION = VERTEX_COUNT - 1
KERNEL_DIMENSION = 2
TARGET_DIMENSION = 3
NORMAL = NormalDist()
HALTON_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29)
EDGE_PAIRS = tuple(
    (first, second)
    for first in range(VERTEX_COUNT)
    for second in range(first + 1, VERTEX_COUNT)
    if (second - first) % VERTEX_COUNT not in (1, VERTEX_COUNT - 1)
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=2048, help="number of kernels to sample")
    parser.add_argument("--seed", type=int, default=20260604, help="seed for pseudo sampling or Halton shift")
    parser.add_argument(
        "--sampler",
        choices=("halton", "pseudo"),
        default="halton",
        help="kernel sampler in Gaussian coordinates before QR",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="output directory; default is GaussianKnots/results/n6_grassmannian_<sampler>_<samples>",
    )
    parser.add_argument(
        "--no-fast",
        action="store_true",
        help="do not request pyknotid fast/Cython helpers where available",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.samples < 1:
        raise SystemExit("--samples must be positive")

    output_dir = args.output_dir or GAUSSIAN_KNOTS / "results" / f"n6_grassmannian_{args.sampler}_{args.samples}"
    output_dir.mkdir(parents=True, exist_ok=True)

    environment = inspect_pyknotid_environment()
    if not environment.available:
        raise SystemExit(f"pyknotid is not available: {environment.error}")

    simplex_vertices = centered_simplex_vertices(VERTEX_COUNT) @ h_basis(VERTEX_COUNT)
    wall_data = determinant_wall_data(simplex_vertices)
    wall_normals = determinant_wall_normals(wall_data)
    arrangement_lifts = wall_automorphism_lifts(wall_normals)
    d6_lifts = wall_automorphism_lifts(wall_normals, allowed_permutations=dihedral_wall_permutations())
    sampler = kernel_sampler(args.sampler, args.samples, args.seed)

    records = []
    for sample_index, gaussian_vector in enumerate(sampler):
        kernel_basis = orthonormal_kernel(gaussian_vector)
        quotient_basis = orthogonal_complement(kernel_basis)
        vertices = simplex_vertices @ quotient_basis
        identification = identify_polygon(vertices, use_fast=not args.no_fast)
        wall_signature = determinant_wall_signature(wall_data, kernel_basis)
        d6_wall_signature = canonical_wall_signature(wall_signature, d6_lifts)
        arrangement_wall_signature = canonical_wall_signature(wall_signature, arrangement_lifts)
        records.append(
            {
                "N": VERTEX_COUNT,
                "sample_index": sample_index,
                "sampler": args.sampler,
                "seed": args.seed,
                "status": identification.status,
                "is_nontrivial": optional_bool(identification.is_nontrivial),
                "knot_label": knot_label(identification),
                "knot_types": ";".join(identification.knot_types),
                "determinant": identification.determinant or "",
                "vassiliev_2": identification.vassiliev_2 or "",
                "vassiliev_3": identification.vassiliev_3 or "",
                "crossing_count": optional_int(identification.crossing_count),
                "simplified_crossing_count": optional_int(identification.simplified_crossing_count),
                "wall_signature": wall_signature,
                "d6_wall_signature": d6_wall_signature,
                "arrangement_wall_signature": arrangement_wall_signature,
                "embedded_numerically": optional_bool(is_numerically_embedded(vertices)),
                "messages": " | ".join(identification.messages),
            }
        )

    write_csv(output_dir / "samples.csv", records)
    write_csv(output_dir / "type_counts.csv", type_counts(records, args.samples))
    write_csv(output_dir / "wall_signature_counts.csv", wall_signature_counts(records, args.samples))
    write_csv(
        output_dir / "d6_wall_signature_counts.csv",
        wall_signature_counts(records, args.samples, signature_key="d6_wall_signature"),
    )
    write_csv(
        output_dir / "arrangement_wall_signature_counts.csv",
        wall_signature_counts(records, args.samples, signature_key="arrangement_wall_signature"),
    )
    write_metadata(output_dir, args, environment, arrangement_lifts, d6_lifts)

    summary = summarize(records, args.samples)
    print(f"wrote outputs to {output_dir}")
    print(
        f"N=6 samples={args.samples} classified={summary['classified']} "
        f"nontrivial={summary['nontrivial']} trivial={summary['trivial']} "
        f"unknown={summary['unknown']} p_hat_3_1={summary['trefoil_rate']:.6f}"
    )
    return 0


def centered_simplex_vertices(vertex_count: int) -> np.ndarray:
    identity = np.eye(vertex_count)
    return identity - np.ones((vertex_count, vertex_count)) / vertex_count


def h_basis(vertex_count: int) -> np.ndarray:
    generators = np.eye(vertex_count)[:, : vertex_count - 1] - np.eye(vertex_count)[:, [vertex_count - 1]]
    q_matrix, _ = np.linalg.qr(generators)
    return q_matrix[:, : vertex_count - 1]


def nonadjacent_edge_pairs(vertex_count: int) -> list[tuple[int, int]]:
    return list(EDGE_PAIRS)


def determinant_wall_data(simplex_vertices: np.ndarray) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    data = []
    for first_edge, second_edge in nonadjacent_edge_pairs(len(simplex_vertices)):
        p0 = simplex_vertices[first_edge]
        p1 = simplex_vertices[(first_edge + 1) % len(simplex_vertices)]
        q0 = simplex_vertices[second_edge]
        q1 = simplex_vertices[(second_edge + 1) % len(simplex_vertices)]
        data.append((p1 - p0, q1 - q0, q0 - p0))
    return data


def kernel_sampler(sampler: str, samples: int, seed: int):
    if sampler == "pseudo":
        rng = np.random.default_rng(seed)
        for _ in range(samples):
            yield rng.normal(size=H_DIMENSION * KERNEL_DIMENSION)
        return

    rng = np.random.default_rng(seed)
    shift = rng.random(H_DIMENSION * KERNEL_DIMENSION)
    for index in range(samples):
        uniform = np.asarray(
            [
                (radical_inverse(index + 1, base) + shift[dimension]) % 1.0
                for dimension, base in enumerate(HALTON_PRIMES)
            ],
            dtype=float,
        )
        yield inverse_normal(uniform)


def radical_inverse(index: int, base: int) -> float:
    value = 0.0
    factor = 1.0 / base
    while index:
        index, digit = divmod(index, base)
        value += digit * factor
        factor /= base
    return value


def inverse_normal(values: np.ndarray) -> np.ndarray:
    clipped = np.clip(values, 1e-12, 1.0 - 1e-12)
    return np.asarray([NORMAL.inv_cdf(float(value)) for value in clipped], dtype=float)


def orthonormal_kernel(gaussian_vector: np.ndarray) -> np.ndarray:
    matrix = np.asarray(gaussian_vector, dtype=float).reshape(H_DIMENSION, KERNEL_DIMENSION)
    q_matrix, r_matrix = np.linalg.qr(matrix, mode="reduced")
    for column in range(KERNEL_DIMENSION):
        if r_matrix[column, column] < 0:
            q_matrix[:, column] *= -1.0
    return q_matrix


def orthogonal_complement(kernel_basis: np.ndarray) -> np.ndarray:
    _, _, vh_matrix = np.linalg.svd(kernel_basis.T, full_matrices=True)
    return vh_matrix[KERNEL_DIMENSION:].T


def determinant_wall_signature(
    wall_data: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
    kernel_basis: np.ndarray,
) -> str:
    k1 = kernel_basis[:, 0]
    k2 = kernel_basis[:, 1]
    signs = []
    for edge_a, edge_b, displacement in wall_data:
        value = float(np.linalg.det(np.column_stack([edge_a, edge_b, k1, k2, displacement])))
        if abs(value) < 1e-12:
            signs.append("0")
        elif value > 0:
            signs.append("+")
        else:
            signs.append("-")
    signature = "".join(signs)
    antipodal = "".join("+" if char == "-" else "-" if char == "+" else "0" for char in signature)
    return min(signature, antipodal)


def determinant_wall_normals(wall_data: list[tuple[np.ndarray, np.ndarray, np.ndarray]]) -> np.ndarray:
    normals = []
    for edge_a, edge_b, displacement in wall_data:
        matrix = np.column_stack([edge_a, edge_b, displacement])
        coordinates = []
        for first in range(H_DIMENSION):
            for second in range(first + 1, H_DIMENSION):
                complement = [index for index in range(H_DIMENSION) if index not in (first, second)]
                sequence = complement + [first, second]
                inversions = sum(
                    1
                    for left in range(H_DIMENSION)
                    for right in range(left + 1, H_DIMENSION)
                    if sequence[left] > sequence[right]
                )
                coordinates.append(((-1) ** inversions) * np.linalg.det(matrix[complement, :]))
        normal = np.asarray(coordinates, dtype=float)
        normal /= np.linalg.norm(normal)
        pivot = int(np.argmax(np.abs(normal)))
        if normal[pivot] < 0:
            normal *= -1.0
        normals.append(normal)
    return np.vstack(normals)


def wall_automorphism_lifts(
    wall_normals: np.ndarray,
    allowed_permutations: list[tuple[int, ...]] | None = None,
) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    gram = np.round(wall_normals @ wall_normals.T, 12)
    projective_gram = np.abs(gram)
    permutations = allowed_permutations or itertools.permutations(range(len(wall_normals)))
    lifts = []
    for permutation in permutations:
        if not all(
            projective_gram[row, column] == projective_gram[permutation[row], permutation[column]]
            for row in range(len(wall_normals))
            for column in range(len(wall_normals))
        ):
            continue
        signs = sign_lift(gram, permutation)
        if signs is not None:
            lifts.append((tuple(permutation), signs))
    return lifts


def sign_lift(gram: np.ndarray, permutation: tuple[int, ...]) -> tuple[int, ...] | None:
    signs: list[int | None] = [None] * len(permutation)
    signs[0] = 1
    stack = [0]
    while stack:
        row = stack.pop()
        for column in range(len(permutation)):
            if row == column or abs(gram[row, column]) < 1e-10:
                continue
            ratio = int(round(gram[row, column] / gram[permutation[row], permutation[column]]))
            value = signs[row] * ratio
            if signs[column] is None:
                signs[column] = value
                stack.append(column)
            elif signs[column] != value:
                return None
    if any(sign is None for sign in signs):
        return None
    resolved = tuple(int(sign) for sign in signs)
    if np.allclose(gram, np.outer(resolved, resolved) * gram[np.ix_(permutation, permutation)]):
        return resolved
    return None


def canonical_wall_signature(signature: str, lifts: list[tuple[tuple[int, ...], tuple[int, ...]]]) -> str:
    candidates = []
    for permutation, signs in lifts:
        transformed = [""] * len(signature)
        for source, target in enumerate(permutation):
            char = signature[source]
            if signs[source] < 0:
                char = flip_sign(char)
            transformed[target] = char
        candidate = "".join(transformed)
        candidates.append(min(candidate, flip_signature(candidate)))
    return min(candidates)


def flip_signature(signature: str) -> str:
    return "".join(flip_sign(char) for char in signature)


def flip_sign(char: str) -> str:
    if char == "+":
        return "-"
    if char == "-":
        return "+"
    return char


def dihedral_wall_permutations() -> list[tuple[int, ...]]:
    index = {pair: position for position, pair in enumerate(EDGE_PAIRS)}
    permutations = []
    for shift in range(VERTEX_COUNT):
        permutations.append(wall_permutation({i: (i + shift) % VERTEX_COUNT for i in range(VERTEX_COUNT)}, index))
        permutations.append(wall_permutation({i: (shift - i) % VERTEX_COUNT for i in range(VERTEX_COUNT)}, index))
    return permutations


def wall_permutation(vertex_permutation: dict[int, int], index: dict[tuple[int, int], int]) -> tuple[int, ...]:
    values = []
    for first, second in EDGE_PAIRS:
        image = tuple(sorted((vertex_permutation[first], vertex_permutation[second])))
        values.append(index[image])
    return tuple(values)


def knot_label(identification: object) -> str:
    knot_types = getattr(identification, "knot_types")
    is_nontrivial = getattr(identification, "is_nontrivial")
    if knot_types:
        return ";".join(knot_types)
    if is_nontrivial is False:
        return "0_1"
    if is_nontrivial is True:
        return "3_1_or_nontrivial"
    return "unknown"


def optional_bool(value: bool | None) -> str:
    if value is None:
        return ""
    return "true" if value else "false"


def optional_int(value: int | None) -> str:
    if value is None:
        return ""
    return str(value)


def type_counts(records: list[dict[str, object]], samples: int) -> list[dict[str, object]]:
    counts: Counter[str] = Counter(str(record["knot_label"]) for record in records)
    rows = []
    for label, count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
        rows.append({"knot_label": label, "count": count, "rate": count / samples})
    return rows


def wall_signature_counts(
    records: list[dict[str, object]],
    samples: int,
    signature_key: str = "wall_signature",
) -> list[dict[str, object]]:
    buckets: dict[str, Counter[str]] = defaultdict(Counter)
    for record in records:
        buckets[str(record[signature_key])][str(record["knot_label"])] += 1
    rows = []
    for signature, labels in sorted(buckets.items(), key=lambda item: (-sum(item[1].values()), item[0])):
        total = sum(labels.values())
        rows.append(
            {
                signature_key: signature,
                "count": total,
                "rate": total / samples,
                "labels": json.dumps(dict(sorted(labels.items())), sort_keys=True),
            }
        )
    return rows


def summarize(records: list[dict[str, object]], samples: int) -> dict[str, int | float]:
    nontrivial = sum(1 for record in records if record["is_nontrivial"] == "true")
    trivial = sum(1 for record in records if record["is_nontrivial"] == "false")
    unknown = samples - nontrivial - trivial
    return {
        "classified": nontrivial + trivial,
        "nontrivial": nontrivial,
        "trivial": trivial,
        "unknown": unknown,
        "trefoil_rate": nontrivial / samples,
    }


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("cannot write empty CSV")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_metadata(
    output_dir: Path,
    args: argparse.Namespace,
    environment: object,
    arrangement_lifts: list[tuple[tuple[int, ...], tuple[int, ...]]],
    d6_lifts: list[tuple[tuple[int, ...], tuple[int, ...]]],
) -> None:
    metadata = {
        "N": VERTEX_COUNT,
        "samples": args.samples,
        "seed": args.seed,
        "sampler": args.sampler,
        "quotient_model": "Gr_2(1^perp) kernel planes for N=6",
        "wall_count": len(nonadjacent_edge_pairs(VERTEX_COUNT)),
        "wall_pair_order": list(EDGE_PAIRS),
        "full_projective_wall_group_order": len(arrangement_lifts),
        "d6_wall_subgroup_order": len(d6_lifts),
        "wall_orbits_under_full_projective_group": [
            [list(pair) for pair in orbit]
            for orbit in wall_orbits_under_group([permutation for permutation, _ in arrangement_lifts])
        ],
        "wall_orbits_under_D6": [
            [list(pair) for pair in orbit] for orbit in wall_orbits_under_group(dihedral_wall_permutations())
        ],
        "pyknotid_environment": {
            "available": environment.available,
            "version": environment.version,
            "fast_backend_available": environment.fast_backend_available,
            "error": environment.error,
        },
    }
    with (output_dir / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")


def wall_orbits_under_group(permutations: object) -> list[list[tuple[int, int]]]:
    seen: set[tuple[int, int]] = set()
    orbits = []
    for pair in EDGE_PAIRS:
        if pair in seen:
            continue
        position = EDGE_PAIRS.index(pair)
        orbit_positions = sorted({permutation[position] for permutation in permutations})
        orbit = [EDGE_PAIRS[index] for index in orbit_positions]
        seen.update(orbit)
        orbits.append(orbit)
    return orbits


if __name__ == "__main__":
    raise SystemExit(main())
