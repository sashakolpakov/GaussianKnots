#!/usr/bin/env python3
"""Full N=6 order-type/symmetry analysis for the Stiefel knot model.

This script is intentionally focused on the first nontrivial case N=6.  It
combines:

* the 9 edge-pair determinant-wall symmetry group;
* its orthogonal lift to H = 1^perp;
* the safe D6 subgroup coming from cyclic Hamiltonian relabeling;
* the 15-sign order-type bucket proxy for six labeled points in R^3;
* classification of one representative per D6 order-type orbit;
* volumetric Monte Carlo estimates by order-type bucket; and
* comparison with the earlier direct pyknotid Monte Carlo table.

Run with the GaussianKnots virtual environment when classification is needed:

    python3 scripts/repro/n6_chamber_analysis.py
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from gaussian_knots.generation import projected_simplex_polygon  # noqa: E402
from gaussian_knots.pyknotid_adapter import identify_polygon, inspect_pyknotid_environment  # noqa: E402

from gaussian_knots.n6_geometry import (  # noqa: E402
    EDGE_PAIRS,
    H_DIMENSION,
    KERNEL_DIMENSION,
    QUADS,
    VERTEX_COUNT,
    centered_simplex_vertices,
    clique_action,
    common_clique_lines,
    d6_orbit,
    d6_orbit_representative,
    determinant_wall_data,
    determinant_wall_normals,
    dihedral_wall_permutations,
    h_basis,
    maximal_cliques,
    nonorthogonality_graph,
    orientation_signature,
    orthogonal_lift_errors,
    projective_automorphism_lifts,
    wilson_interval,
)


DEFAULT_TREFOIL_SEEDS = (1062611651, 696023607, 2485279349)


@dataclass
class OrbitRepresentative:
    signature: str
    kernel_basis: np.ndarray
    count: int = 0
    source: str = "sample"


@dataclass
class Classification:
    label: str
    status: str
    determinant: str
    vassiliev_2: str
    vassiliev_3: str


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def add(self, item: str) -> None:
        self.parent.setdefault(item, item)

    def find(self, item: str) -> str:
        self.add(item)
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root != right_root:
            self.parent[max(left_root, right_root)] = min(left_root, right_root)

    def classes(self) -> dict[str, list[str]]:
        result: dict[str, list[str]] = defaultdict(list)
        for item in sorted(self.parent):
            result[self.find(item)].append(item)
        return dict(result)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=200000, help="fast Haar samples for order-type volume estimates")
    parser.add_argument("--seed", type=int, default=20260604, help="seed for the order-type volume run")
    parser.add_argument(
        "--simple-mc",
        type=Path,
        default=ROOT / "results" / "haar_N5-12_1000" / "samples_N6.csv",
        help="existing direct classifier Monte Carlo CSV to compare against",
    )
    parser.add_argument(
        "--trefoil-seeds",
        type=int,
        nargs="*",
        default=list(DEFAULT_TREFOIL_SEEDS),
        help="known N=6 trefoil seeds to force into the representative set",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "n6_full_chamber_analysis.csv",
        help="CSV output for D6 order-type orbit representatives",
    )
    parser.add_argument("--no-classify", action="store_true", help="skip pyknotid classification of orbit representatives")
    parser.add_argument(
        "--max-closure-rounds",
        type=int,
        default=4,
        help="maximum rounds used to close sampled D6 order-type reps under the full 72-element wall group",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if args.samples < 1:
        raise SystemExit("--samples must be positive")

    simplex_vertices = centered_simplex_vertices(VERTEX_COUNT) @ h_basis(VERTEX_COUNT)
    wall_data = determinant_wall_data(simplex_vertices)
    wall_normals = determinant_wall_normals(wall_data)
    projective_gram = np.abs(wall_normals @ wall_normals.T)
    graph = nonorthogonality_graph(projective_gram)
    cliques = maximal_cliques(graph, size=3)
    common_lines = common_clique_lines(wall_data, cliques)
    full_lifts = projective_automorphism_lifts(wall_normals @ wall_normals.T)
    d6_wall_perms = set(dihedral_wall_permutations())
    d6_lifts = [lift for lift in full_lifts if lift[0] in d6_wall_perms]
    full_maps = orthogonal_maps(common_lines, cliques, full_lifts)
    full_errors = orthogonal_lift_errors(
        wall_normals,
        common_lines,
        cliques,
        full_lifts,
    )

    print_header("N=6 Wall Symmetry")
    print(f"wall order: {list(EDGE_PAIRS)}")
    print(f"projective Gram values: {sorted({round(float(value), 12) for value in projective_gram.ravel()})}")
    print(f"full wall group order: {len(full_lifts)}")
    print(f"D6 knot-model subgroup order: {len(d6_lifts)}")
    print(f"extra wall symmetries not in D6: {len(full_lifts) - len(d6_lifts)}")
    print("note: this is 60 extra symmetries, not 70")
    print(f"max orthogonal lift error ||A^T A-I||_F: {max(error[0] for error in full_errors):.3e}")
    print(f"max wall-normal image error: {max(error[1] for error in full_errors):.3e}")

    representatives = sample_order_type_representatives(args.samples, args.seed)
    add_seed_representatives(representatives, args.trefoil_seeds)
    close_under_full_wall_group(representatives, simplex_vertices, full_maps, args.max_closure_rounds)
    classifications = classify_representatives(representatives, simplex_vertices, skip=args.no_classify)
    full_classes = full_isometry_classes(representatives, simplex_vertices, full_maps)
    simple_mc = read_simple_mc(args.simple_mc)

    print_header("Order-Type Bucket Proxy")
    sampled_reps = [rep for rep in representatives.values() if rep.count > 0]
    raw_signature_count = sum(len(d6_orbit(rep.signature)) for rep in sampled_reps)
    print(f"15-sign order: {list(QUADS)}")
    print(f"sampled Haar points: {args.samples}")
    print(f"sampled D6 order-type orbit representatives: {len(sampled_reps)}")
    print(f"sampled raw 15-sign signatures represented by those D6 orbits: {raw_signature_count}")
    print(f"closure representatives under full wall group: {len(representatives)}")
    print(f"full orthogonal-isometry classes among representatives: {len(full_classes)}")

    label_summary = summarize_labels(representatives, classifications, args.samples)
    print_header("Knot Labels And Volumes")
    print_label_summary(label_summary, args.samples)

    print_header("Full Wall Group Versus Knot Type")
    print_full_class_summary(full_classes, representatives, classifications)
    print_trefoil_image_test(simplex_vertices, common_lines, cliques, full_lifts, d6_wall_perms)

    print_header("Comparison With Direct Monte Carlo")
    print_simple_mc_comparison(simple_mc, label_summary, args.samples)

    write_rows(args.output, representatives, classifications, full_classes, args.samples)
    print()
    print(f"wrote representative table: {args.output}")
    return 0


def sample_order_type_representatives(samples: int, seed: int) -> dict[str, OrbitRepresentative]:
    rng = np.random.default_rng(seed)
    representatives: dict[str, OrbitRepresentative] = {}
    for _ in range(samples):
        kernel = random_kernel(rng)
        vertices = vertices_from_kernel(centered_simplex_vertices(VERTEX_COUNT) @ h_basis(VERTEX_COUNT), kernel)
        signature = orientation_signature(vertices)
        representative = d6_orbit_representative(signature)
        if representative not in representatives:
            representatives[representative] = OrbitRepresentative(representative, kernel, 0, "sample")
        representatives[representative].count += 1
    return representatives


def add_seed_representatives(representatives: dict[str, OrbitRepresentative], seeds: list[int]) -> None:
    simplex_vertices = centered_simplex_vertices(VERTEX_COUNT) @ h_basis(VERTEX_COUNT)
    for seed in seeds:
        kernel = kernel_from_haar_seed(simplex_vertices, seed)
        signature = d6_orbit_representative(orientation_signature(vertices_from_kernel(simplex_vertices, kernel)))
        representatives.setdefault(signature, OrbitRepresentative(signature, kernel, 0, f"seed:{seed}"))


def close_under_full_wall_group(
    representatives: dict[str, OrbitRepresentative],
    simplex_vertices: np.ndarray,
    full_maps: list[np.ndarray],
    max_rounds: int,
) -> None:
    queue = deque(representatives.keys())
    rounds = 0
    while queue and rounds < max_rounds:
        rounds += 1
        current_keys = list(queue)
        queue.clear()
        for signature in current_keys:
            kernel = representatives[signature].kernel_basis
            for linear_map in full_maps:
                image_kernel = linear_map @ kernel
                image_signature = d6_orbit_representative(orientation_signature(vertices_from_kernel(simplex_vertices, image_kernel)))
                if image_signature not in representatives:
                    representatives[image_signature] = OrbitRepresentative(
                        image_signature,
                        image_kernel,
                        0,
                        f"full-image-of:{signature}",
                    )
                    queue.append(image_signature)


def classify_representatives(
    representatives: dict[str, OrbitRepresentative],
    simplex_vertices: np.ndarray,
    skip: bool,
) -> dict[str, Classification]:
    if skip:
        return {
            signature: Classification("unclassified", "skipped", "", "", "")
            for signature in representatives
        }

    environment = inspect_pyknotid_environment()
    if not environment.available:
        return {
            signature: Classification("unknown", f"pyknotid_unavailable:{environment.error}", "", "", "")
            for signature in representatives
        }

    result = {}
    for signature, representative in sorted(representatives.items(), key=lambda item: (-item[1].count, item[0])):
        vertices = vertices_from_kernel(simplex_vertices, representative.kernel_basis)
        identification = identify_polygon(vertices, use_fast=True)
        if identification.knot_types:
            label = ";".join(identification.knot_types)
        elif identification.is_nontrivial is True:
            label = "nontrivial"
        elif identification.is_nontrivial is False:
            label = "0_1"
        else:
            label = "unknown"
        result[signature] = Classification(
            label=label,
            status=identification.status,
            determinant=identification.determinant or "",
            vassiliev_2=identification.vassiliev_2 or "",
            vassiliev_3=identification.vassiliev_3 or "",
        )
    return result


def full_isometry_classes(
    representatives: dict[str, OrbitRepresentative],
    simplex_vertices: np.ndarray,
    full_maps: list[np.ndarray],
) -> dict[str, list[str]]:
    union_find = UnionFind()
    for signature in representatives:
        union_find.add(signature)
    for signature, representative in representatives.items():
        for linear_map in full_maps:
            image_signature = d6_orbit_representative(
                orientation_signature(vertices_from_kernel(simplex_vertices, linear_map @ representative.kernel_basis))
            )
            if image_signature in representatives:
                union_find.union(signature, image_signature)
    return union_find.classes()


def orthogonal_maps(
    common_lines: np.ndarray,
    cliques: list[tuple[int, ...]],
    full_lifts: list[tuple[tuple[int, ...], tuple[int, ...]]],
) -> list[np.ndarray]:
    source = common_lines.T
    source_pinv = np.linalg.pinv(source)
    maps = []
    for permutation, _ in full_lifts:
        action = clique_action(permutation, cliques)
        target = common_lines[list(action)].T
        maps.append(target @ source_pinv)
    return maps


def print_trefoil_image_test(
    simplex_vertices: np.ndarray,
    common_lines: np.ndarray,
    cliques: list[tuple[int, ...]],
    full_lifts: list[tuple[tuple[int, ...], tuple[int, ...]]],
    d6_wall_perms: set[tuple[int, ...]],
) -> None:
    kernel = kernel_from_haar_seed(simplex_vertices, DEFAULT_TREFOIL_SEEDS[0])
    source_label = classify_kernel(simplex_vertices, kernel).label
    counts_all: Counter[str] = Counter()
    counts_d6: Counter[str] = Counter()
    counts_extra: Counter[str] = Counter()
    for permutation, _ in full_lifts:
        linear_map = orthogonal_maps(common_lines, cliques, [(permutation, ())])[0]
        label = classify_kernel(simplex_vertices, linear_map @ kernel).label
        counts_all[label] += 1
        if permutation in d6_wall_perms:
            counts_d6[label] += 1
        else:
            counts_extra[label] += 1
    print(f"known trefoil seed {DEFAULT_TREFOIL_SEEDS[0]} source label: {source_label}")
    print(f"images under full 72-element wall group: {dict(counts_all)}")
    print(f"images under D6 subgroup: {dict(counts_d6)}")
    print(f"images under extra non-D6 wall symmetries: {dict(counts_extra)}")


def classify_kernel(simplex_vertices: np.ndarray, kernel_basis: np.ndarray) -> Classification:
    vertices = vertices_from_kernel(simplex_vertices, kernel_basis)
    identification = identify_polygon(vertices, use_fast=True)
    if identification.knot_types:
        label = ";".join(identification.knot_types)
    elif identification.is_nontrivial is True:
        label = "nontrivial"
    elif identification.is_nontrivial is False:
        label = "0_1"
    else:
        label = "unknown"
    return Classification(label, identification.status, identification.determinant or "", identification.vassiliev_2 or "", identification.vassiliev_3 or "")


def random_kernel(rng: np.random.Generator) -> np.ndarray:
    matrix = rng.normal(size=(H_DIMENSION, KERNEL_DIMENSION))
    q_matrix, r_matrix = np.linalg.qr(matrix, mode="reduced")
    for column in range(KERNEL_DIMENSION):
        if r_matrix[column, column] < 0:
            q_matrix[:, column] *= -1.0
    return q_matrix


def vertices_from_kernel(simplex_vertices: np.ndarray, kernel_basis: np.ndarray) -> np.ndarray:
    _, _, vh_matrix = np.linalg.svd(kernel_basis.T, full_matrices=True)
    quotient_basis = vh_matrix[KERNEL_DIMENSION:].T
    return simplex_vertices @ quotient_basis


def kernel_from_haar_seed(simplex_vertices: np.ndarray, seed: int) -> np.ndarray:
    vertices = projected_simplex_polygon(VERTEX_COUNT, np.random.default_rng(seed), projection_model="haar")
    quotient_basis = np.linalg.lstsq(simplex_vertices, vertices, rcond=None)[0]
    _, _, vh_matrix = np.linalg.svd(quotient_basis.T, full_matrices=True)
    return vh_matrix[3:].T


def summarize_labels(
    representatives: dict[str, OrbitRepresentative],
    classifications: dict[str, Classification],
    samples: int,
) -> dict[str, dict[str, float | int]]:
    summary: dict[str, dict[str, float | int]] = defaultdict(
        lambda: {"d6_orbits": 0, "raw_order_type_signatures": 0, "count": 0}
    )
    for signature, representative in representatives.items():
        label = classifications[signature].label
        summary[label]["d6_orbits"] += 1
        summary[label]["raw_order_type_signatures"] += len(d6_orbit(signature))
        summary[label]["count"] += representative.count
    for values in summary.values():
        count = int(values["count"])
        lower, upper = wilson_interval(count, samples)
        values["rate"] = count / samples
        values["wilson95_lower"] = lower
        values["wilson95_upper"] = upper
    return dict(summary)


def print_label_summary(summary: dict[str, dict[str, float | int]], samples: int) -> None:
    for label, values in sorted(summary.items(), key=lambda item: (-int(item[1]["count"]), item[0])):
        print(
            f"{label}: d6_orbits={values['d6_orbits']} "
            f"raw_15sign_signatures={values['raw_order_type_signatures']} "
            f"count={values['count']}/{samples} rate={values['rate']:.6f} "
            f"Wilson95=({values['wilson95_lower']:.6f}, {values['wilson95_upper']:.6f})"
        )


def print_full_class_summary(
    full_classes: dict[str, list[str]],
    representatives: dict[str, OrbitRepresentative],
    classifications: dict[str, Classification],
) -> None:
    mixed = []
    for root, signatures in full_classes.items():
        labels = Counter(classifications[signature].label for signature in signatures)
        sampled_count = sum(representatives[signature].count for signature in signatures)
        if len(labels) > 1:
            mixed.append((root, signatures, labels, sampled_count))
    print(f"full orthogonal-isometry classes: {len(full_classes)}")
    print(f"classes mixing knot labels: {len(mixed)}")
    for root, signatures, labels, sampled_count in sorted(mixed, key=lambda item: (-item[3], item[0]))[:12]:
        print(f"  class {root}: reps={len(signatures)} sampled_count={sampled_count} labels={dict(labels)}")


def read_simple_mc(path: Path) -> dict[str, int | float | tuple[float, float]]:
    if not path.exists():
        return {"available": 0}
    counts: Counter[str] = Counter()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            label = row["knot_label"] or ("unknown" if not row["is_nontrivial"] else "nontrivial")
            counts[label] += 1
    total = sum(counts.values())
    trefoil = counts.get("3_1", 0)
    lower, upper = wilson_interval(trefoil, total)
    return {
        "available": 1,
        "total": total,
        "counts": counts,
        "trefoil_rate": trefoil / total if total else math.nan,
        "trefoil_lower": lower,
        "trefoil_upper": upper,
    }


def print_simple_mc_comparison(
    simple_mc: dict[str, object],
    label_summary: dict[str, dict[str, float | int]],
    samples: int,
) -> None:
    if not simple_mc.get("available"):
        print("direct classifier Monte Carlo file not found")
        return
    trefoil_values = label_summary.get("3_1", {"count": 0, "rate": 0.0, "wilson95_lower": 0.0, "wilson95_upper": 0.0})
    print(
        f"direct classifier MC: {simple_mc['counts']} "
        f"trefoil_rate={simple_mc['trefoil_rate']:.6f} "
        f"Wilson95=({simple_mc['trefoil_lower']:.6f}, {simple_mc['trefoil_upper']:.6f})"
    )
    print(
        f"order-type volume MC: trefoil_count={trefoil_values['count']}/{samples} "
        f"rate={trefoil_values['rate']:.6f} "
        f"Wilson95=({trefoil_values['wilson95_lower']:.6f}, {trefoil_values['wilson95_upper']:.6f})"
    )
    print("interpretation: the direct 1000-sample run is noisy at N=6; the order-type run avoids classifier calls in the large loop")


def write_rows(
    output: Path,
    representatives: dict[str, OrbitRepresentative],
    classifications: dict[str, Classification],
    full_classes: dict[str, list[str]],
    samples: int,
) -> None:
    class_index = {}
    for index, signatures in enumerate(sorted(full_classes.values(), key=lambda values: (len(values), values[0]))):
        for signature in signatures:
            class_index[signature] = index
    rows = []
    for signature, representative in sorted(representatives.items(), key=lambda item: (-item[1].count, item[0])):
        classification = classifications[signature]
        lower, upper = wilson_interval(representative.count, samples)
        rows.append(
            {
                "d6_order_type_representative": signature,
                "full_isometry_class": class_index[signature],
                "d6_orbit_size": len(d6_orbit(signature)),
                "sample_count": representative.count,
                "sample_rate": representative.count / samples,
                "wilson95_lower": lower,
                "wilson95_upper": upper,
                "label": classification.label,
                "status": classification.status,
                "determinant": classification.determinant,
                "vassiliev_2": classification.vassiliev_2,
                "vassiliev_3": classification.vassiliev_3,
                "source": representative.source,
            }
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def print_header(title: str) -> None:
    print()
    print("=" * len(title))
    print(title)
    print("=" * len(title))


if __name__ == "__main__":
    raise SystemExit(main())
