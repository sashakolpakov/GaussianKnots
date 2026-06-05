#!/usr/bin/env python3
"""Two-phase order-type catalog and volume estimator for N=6,7,8.

Phase 1 samples Haar kernels, records D_N-canonical order-type buckets, and
classifies catalog representatives.

Phase 2 samples Haar kernels again and estimates the Haar mass of the classified
catalog buckets.  Samples whose bucket is not in the catalog are reported
explicitly as uncatalogued.  This is an order-type proxy, not a proof that the
catalog buckets are connected chambers.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from gaussian_knots.pyknotid_adapter import identify_polygon, inspect_pyknotid_environment  # noqa: E402
from gaussian_knots.general_chamber_geometry import (  # noqa: E402
    ChamberModel,
    build_chamber_model,
    d_orbit,
    d_orbit_representative,
    determinant_wall_signature,
    dihedral_wall_permutations,
    is_projective_automorphism,
    orientation_signature,
    random_kernel,
    vertices_from_kernel,
    wall_gram_automorphisms,
    wilson_interval,
)


@dataclass
class CatalogEntry:
    signature: str
    kernel_basis: np.ndarray
    wall_signature: str
    catalog_count: int = 0
    volume_count: int = 0
    label: str = "unclassified"
    status: str = ""
    determinant: str = ""
    vassiliev_2: str = ""
    vassiliev_3: str = ""


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vertices", default="6,7,8", help="comma-separated N values")
    parser.add_argument("--catalog-samples", type=int, default=1000, help="Haar samples used to discover the catalog")
    parser.add_argument("--volume-samples", type=int, default=500000, help="Haar samples used for the volume estimate")
    parser.add_argument("--seed", type=int, default=20260604, help="base random seed")
    parser.add_argument(
        "--max-catalog-classifications",
        type=int,
        default=None,
        help="classify only this many highest-count catalog buckets; default classifies every catalog bucket",
    )
    parser.add_argument("--no-classify", action="store_true", help="skip catalog representative classification")
    parser.add_argument("--no-fast", action="store_true", help="do not request pyknotid fast helpers")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results" / "order_type_volume_N6-8")
    parser.add_argument("--progress-every", type=int, default=50000, help="print volume progress every this many samples")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    vertex_counts = parse_vertices(args.vertices)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    for vertex_count in vertex_counts:
        run_vertex_count(vertex_count, args)
    return 0


def parse_vertices(text: str) -> list[int]:
    values = [int(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        raise SystemExit("--vertices must contain at least one N")
    for value in values:
        if value < 5:
            raise SystemExit("all N values must be at least 5")
    return values


def run_vertex_count(vertex_count: int, args: argparse.Namespace) -> None:
    out_dir = args.output_dir / f"N{vertex_count}"
    out_dir.mkdir(parents=True, exist_ok=True)
    model = build_chamber_model(vertex_count)

    print(f"N={vertex_count}: building catalog from {args.catalog_samples} Haar samples", flush=True)
    catalog = build_catalog(model, samples=args.catalog_samples, seed=args.seed + 1009 * vertex_count)
    classify_catalog(model, catalog, args)

    print(f"N={vertex_count}: estimating volumes with {args.volume_samples} Haar samples", flush=True)
    volume_stats = sample_volumes(
        model,
        catalog,
        samples=args.volume_samples,
        seed=args.seed + 7919 * vertex_count,
        progress_every=args.progress_every,
    )

    symmetry = symmetry_summary(model)
    write_catalog(out_dir / "catalog.csv", model, catalog, args.catalog_samples, args.volume_samples)
    write_volume_summary(out_dir / "volume_summary.csv", catalog, volume_stats, args.volume_samples)
    write_json(
        out_dir / "run_summary.json",
        {
            **symmetry,
            "catalog_samples": args.catalog_samples,
            "volume_samples": args.volume_samples,
            "seed": args.seed,
            "catalog_buckets": len(catalog),
            "classified_catalog_buckets": sum(1 for entry in catalog.values() if entry.label != "unclassified"),
            **volume_stats,
        },
    )
    print_summary(vertex_count, catalog, volume_stats, args.volume_samples, out_dir)


def build_catalog(model: ChamberModel, samples: int, seed: int) -> dict[str, CatalogEntry]:
    rng = np.random.default_rng(seed)
    catalog: dict[str, CatalogEntry] = {}
    for _ in range(samples):
        kernel = random_kernel(model, rng)
        vertices = vertices_from_kernel(model, kernel)
        signature = canonical_order_type(model, vertices)
        entry = catalog.get(signature)
        if entry is None:
            entry = CatalogEntry(
                signature=signature,
                kernel_basis=kernel,
                wall_signature=determinant_wall_signature(model.wall_data, kernel),
            )
            catalog[signature] = entry
        entry.catalog_count += 1
    return catalog


def classify_catalog(model: ChamberModel, catalog: dict[str, CatalogEntry], args: argparse.Namespace) -> None:
    if args.no_classify:
        return
    environment = inspect_pyknotid_environment()
    if not environment.available:
        for entry in catalog.values():
            entry.status = f"pyknotid_unavailable:{environment.error}"
        return

    entries = sorted(catalog.values(), key=lambda item: (-item.catalog_count, item.signature))
    if args.max_catalog_classifications is not None:
        entries = entries[: args.max_catalog_classifications]
    for index, entry in enumerate(entries, start=1):
        vertices = vertices_from_kernel(model, entry.kernel_basis)
        identification = identify_polygon(vertices, use_fast=not args.no_fast)
        entry.label = knot_label(identification)
        entry.status = identification.status
        entry.determinant = identification.determinant or ""
        entry.vassiliev_2 = identification.vassiliev_2 or ""
        entry.vassiliev_3 = identification.vassiliev_3 or ""
        if index % 100 == 0:
            print(f"  classified {index}/{len(entries)} catalog buckets", flush=True)


def sample_volumes(
    model: ChamberModel,
    catalog: dict[str, CatalogEntry],
    samples: int,
    seed: int,
    progress_every: int,
) -> dict[str, int]:
    rng = np.random.default_rng(seed)
    uncatalogued_count = 0
    uncatalogued_signatures: set[str] = set()
    for sample_index in range(samples):
        kernel = random_kernel(model, rng)
        vertices = vertices_from_kernel(model, kernel)
        signature = canonical_order_type(model, vertices)
        entry = catalog.get(signature)
        if entry is None:
            uncatalogued_count += 1
            uncatalogued_signatures.add(signature)
        else:
            entry.volume_count += 1
        if progress_every > 0 and (sample_index + 1) % progress_every == 0:
            print(f"  volume samples {sample_index + 1}/{samples}", flush=True)
    return {
        "catalogued_volume_count": samples - uncatalogued_count,
        "uncatalogued_volume_count": uncatalogued_count,
        "uncatalogued_unique_buckets": len(uncatalogued_signatures),
    }


def canonical_order_type(model: ChamberModel, vertices: np.ndarray) -> str:
    signature = orientation_signature(vertices, model.quads)
    return d_orbit_representative(signature, model.vertex_count, model.quads, model.quad_index)


def knot_label(identification: object) -> str:
    knot_types = getattr(identification, "knot_types", ())
    if knot_types:
        return ";".join(knot_types)
    is_nontrivial = getattr(identification, "is_nontrivial", None)
    if is_nontrivial is True:
        return "nontrivial"
    if is_nontrivial is False:
        return "0_1"
    return "unknown"


def symmetry_summary(model: ChamberModel) -> dict[str, object]:
    d_wall_perms = set(dihedral_wall_permutations(model))
    gram_group = wall_gram_automorphisms(model.wall_gram, max_count=200000)
    return {
        "vertex_count": model.vertex_count,
        "wall_count": len(model.edge_pairs),
        "order_type_sign_count": len(model.quads),
        "wall_gram_group_order": len(gram_group),
        "dihedral_subgroup_order": len(d_wall_perms),
        "extra_wall_gram_symmetries": len(set(gram_group) - d_wall_perms),
        "dihedral_preserves_wall_gram": all(is_projective_automorphism(model.wall_gram, perm) for perm in d_wall_perms),
    }


def write_catalog(
    path: Path,
    model: ChamberModel,
    catalog: dict[str, CatalogEntry],
    catalog_samples: int,
    volume_samples: int,
) -> None:
    fieldnames = [
        "rank",
        "signature",
        "wall_signature",
        "d_orbit_size",
        "catalog_count",
        "catalog_rate",
        "volume_count",
        "volume_rate",
        "label",
        "status",
        "determinant",
        "vassiliev_2",
        "vassiliev_3",
    ]
    rows = []
    for rank, entry in enumerate(sorted(catalog.values(), key=lambda item: (-item.volume_count, -item.catalog_count, item.signature)), start=1):
        rows.append(
            {
                "rank": rank,
                "signature": entry.signature,
                "wall_signature": entry.wall_signature,
                "d_orbit_size": len(d_orbit(entry.signature, model.vertex_count, model.quads, model.quad_index)),
                "catalog_count": entry.catalog_count,
                "catalog_rate": f"{entry.catalog_count / catalog_samples:.8f}",
                "volume_count": entry.volume_count,
                "volume_rate": f"{entry.volume_count / volume_samples:.8f}",
                "label": entry.label,
                "status": entry.status,
                "determinant": entry.determinant,
                "vassiliev_2": entry.vassiliev_2,
                "vassiliev_3": entry.vassiliev_3,
            }
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_volume_summary(
    path: Path,
    catalog: dict[str, CatalogEntry],
    volume_stats: dict[str, int],
    samples: int,
) -> None:
    label_counts = Counter()
    for entry in catalog.values():
        label_counts[entry.label] += entry.volume_count
    label_counts["uncatalogued"] += volume_stats["uncatalogued_volume_count"]

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["label", "count", "rate", "wilson95_lower", "wilson95_upper"])
        writer.writeheader()
        for label, count in sorted(label_counts.items(), key=lambda item: (-item[1], item[0])):
            low, high = wilson_interval(count, samples)
            writer.writerow(
                {
                    "label": label,
                    "count": count,
                    "rate": f"{count / samples:.8f}",
                    "wilson95_lower": f"{low:.8f}",
                    "wilson95_upper": f"{high:.8f}",
                }
            )


def write_json(path: Path, payload: dict[str, object]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def print_summary(
    vertex_count: int,
    catalog: dict[str, CatalogEntry],
    volume_stats: dict[str, int],
    samples: int,
    out_dir: Path,
) -> None:
    label_counts = Counter()
    for entry in catalog.values():
        label_counts[entry.label] += entry.volume_count
    label_counts["uncatalogued"] += volume_stats["uncatalogued_volume_count"]
    print(f"N={vertex_count}: catalog buckets={len(catalog)}", flush=True)
    for label, count in label_counts.most_common():
        low, high = wilson_interval(count, samples)
        print(f"  {label}: {count}/{samples} rate={count / samples:.6f} Wilson95=({low:.6f}, {high:.6f})", flush=True)
    print(f"  uncatalogued unique buckets: {volume_stats['uncatalogued_unique_buckets']}", flush=True)
    print(f"  wrote {out_dir}", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
