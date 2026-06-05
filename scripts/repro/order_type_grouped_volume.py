#!/usr/bin/env python3
"""Volume-first grouped Monte Carlo for Stiefel order-type buckets.

This is the order-type-assisted workflow:

1. sample Haar kernels;
2. bucket every sample by D_N-canonical order type;
3. keep a small reservoir of points from each bucket;
4. classify points inside the high-volume buckets;
5. report bucket volumes and within-bucket label checks.

It does not require a precomputed catalog to cover the volume sample, and it
does not enumerate connected chambers of the embedding discriminant complement.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from dataclasses import dataclass, field
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
    dihedral_wall_permutations,
    is_projective_automorphism,
    orientation_signature,
    random_kernel,
    vertices_from_kernel,
    wall_gram_automorphisms,
    wilson_interval,
)


@dataclass
class Bucket:
    signature: str
    count: int = 0
    reservoir: list[np.ndarray] = field(default_factory=list)
    label_counts: Counter[str] = field(default_factory=Counter)
    status_counts: Counter[str] = field(default_factory=Counter)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vertices", default="6,7,8", help="comma-separated N values")
    parser.add_argument("--samples", type=int, default=500000, help="Haar samples for bucket volumes")
    parser.add_argument("--seed", type=int, default=20260604, help="base random seed")
    parser.add_argument("--classify-top-groups", type=int, default=200, help="classify this many largest buckets")
    parser.add_argument("--checks-per-group", type=int, default=3, help="stored points to classify per selected bucket")
    parser.add_argument(
        "--direct-classify-samples",
        type=int,
        default=1000,
        help="also classify this many ordinary Monte Carlo samples directly",
    )
    parser.add_argument("--no-fast", action="store_true", help="do not request pyknotid fast helpers")
    parser.add_argument("--no-classify", action="store_true", help="only compute bucket volumes")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "results" / "order_type_grouped_volume_N6-8")
    parser.add_argument("--progress-every", type=int, default=100000)
    parser.add_argument(
        "--max-bucket-rows",
        type=int,
        default=50000,
        help="write only this many highest-volume bucket rows; use 0 to write all buckets",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for vertex_count in parse_vertices(args.vertices):
        run_vertex_count(vertex_count, args)
    return 0


def parse_vertices(text: str) -> list[int]:
    values = [int(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        raise SystemExit("--vertices must contain at least one N")
    return values


def run_vertex_count(vertex_count: int, args: argparse.Namespace) -> None:
    out_dir = args.output_dir / f"N{vertex_count}"
    out_dir.mkdir(parents=True, exist_ok=True)
    model = build_chamber_model(vertex_count)
    rng = np.random.default_rng(args.seed + 1009 * vertex_count)
    buckets: dict[str, Bucket] = {}
    direct_kernels: list[np.ndarray] = []

    print(f"N={vertex_count}: sampling {args.samples} Haar points into order-type buckets", flush=True)
    for sample_index in range(args.samples):
        kernel = random_kernel(model, rng)
        vertices = vertices_from_kernel(model, kernel)
        signature = canonical_order_type(model, vertices)
        bucket = buckets.get(signature)
        if bucket is None:
            bucket = Bucket(signature=signature)
            buckets[signature] = bucket
        add_to_bucket(bucket, kernel, args.checks_per_group, rng)
        if len(direct_kernels) < args.direct_classify_samples:
            direct_kernels.append(kernel.copy())
        if args.progress_every > 0 and (sample_index + 1) % args.progress_every == 0:
            print(f"  sampled {sample_index + 1}/{args.samples}", flush=True)

    direct_labels: Counter[str] = Counter()
    direct_statuses: Counter[str] = Counter()
    if not args.no_classify:
        classify_buckets(model, buckets, args)
        direct_labels, direct_statuses = classify_direct_samples(model, direct_kernels, args)

    write_direct_rows(out_dir / "direct_sample_summary.csv", direct_labels, direct_statuses, len(direct_kernels))
    write_json(
        out_dir / "run_summary.json",
        {
            **symmetry_summary(model),
            "samples": args.samples,
            "seed": args.seed,
            "bucket_count": len(buckets),
            "classify_top_groups": 0 if args.no_classify else args.classify_top_groups,
            "checks_per_group": args.checks_per_group,
            "direct_classify_samples": 0 if args.no_classify else len(direct_kernels),
            "direct_label_counts": dict(direct_labels),
            "direct_status_counts": dict(direct_statuses),
            "bucket_rows_written": min(len(buckets), args.max_bucket_rows) if args.max_bucket_rows > 0 else len(buckets),
        },
    )
    write_bucket_rows(out_dir / "bucket_volume.csv", model, buckets, args.samples, max_rows=args.max_bucket_rows)
    print_summary(vertex_count, buckets, direct_labels, args.samples, out_dir)


def canonical_order_type(model: ChamberModel, vertices: np.ndarray) -> str:
    signature = orientation_signature(vertices, model.quads)
    return d_orbit_representative(signature, model.vertex_count, model.quads, model.quad_index)


def add_to_bucket(bucket: Bucket, kernel: np.ndarray, reservoir_size: int, rng: np.random.Generator) -> None:
    bucket.count += 1
    if reservoir_size <= 0:
        return
    if len(bucket.reservoir) < reservoir_size:
        bucket.reservoir.append(kernel.copy())
        return
    draw = int(rng.integers(bucket.count))
    if draw < reservoir_size:
        bucket.reservoir[draw] = kernel.copy()


def classify_buckets(model: ChamberModel, buckets: dict[str, Bucket], args: argparse.Namespace) -> None:
    environment = inspect_pyknotid_environment()
    if not environment.available:
        print(f"pyknotid unavailable: {environment.error}", flush=True)
        return
    selected = sorted(buckets.values(), key=lambda item: (-item.count, item.signature))[: args.classify_top_groups]
    for index, bucket in enumerate(selected, start=1):
        for kernel in bucket.reservoir:
            label, status = classify_kernel(model, kernel, use_fast=not args.no_fast)
            bucket.label_counts[label] += 1
            bucket.status_counts[status] += 1
        if index % 100 == 0:
            print(f"  classified {index}/{len(selected)} high-volume buckets", flush=True)


def classify_direct_samples(
    model: ChamberModel,
    kernels: list[np.ndarray],
    args: argparse.Namespace,
) -> tuple[Counter[str], Counter[str]]:
    labels: Counter[str] = Counter()
    statuses: Counter[str] = Counter()
    for index, kernel in enumerate(kernels, start=1):
        label, status = classify_kernel(model, kernel, use_fast=not args.no_fast)
        labels[label] += 1
        statuses[status] += 1
        if index % 100 == 0:
            print(f"  directly classified {index}/{len(kernels)} samples", flush=True)
    return labels, statuses


def classify_kernel(model: ChamberModel, kernel: np.ndarray, use_fast: bool) -> tuple[str, str]:
    vertices = vertices_from_kernel(model, kernel)
    identification = identify_polygon(vertices, use_fast=use_fast)
    if identification.knot_types:
        label = ";".join(identification.knot_types)
    elif identification.is_nontrivial is True:
        label = "nontrivial"
    elif identification.is_nontrivial is False:
        label = "0_1"
    else:
        label = "unknown"
    return label, identification.status


def write_bucket_rows(
    path: Path,
    model: ChamberModel,
    buckets: dict[str, Bucket],
    samples: int,
    max_rows: int,
) -> None:
    fieldnames = [
        "rank",
        "signature",
        "d_orbit_size",
        "count",
        "rate",
        "wilson95_lower",
        "wilson95_upper",
        "check_count",
        "label_counts",
        "status_counts",
        "mixed_labels",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        ordered = sorted(buckets.values(), key=lambda item: (-item.count, item.signature))
        if max_rows > 0:
            ordered = ordered[:max_rows]
        orbit_size_cache: dict[str, int] = {}
        for rank, bucket in enumerate(ordered, start=1):
            low, high = wilson_interval(bucket.count, samples)
            orbit_size = orbit_size_cache.get(bucket.signature)
            if orbit_size is None:
                orbit_size = len(d_orbit(bucket.signature, model.vertex_count, model.quads, model.quad_index))
                orbit_size_cache[bucket.signature] = orbit_size
            writer.writerow(
                {
                    "rank": rank,
                    "signature": bucket.signature,
                    "d_orbit_size": orbit_size,
                    "count": bucket.count,
                    "rate": f"{bucket.count / samples:.8f}",
                    "wilson95_lower": f"{low:.8f}",
                    "wilson95_upper": f"{high:.8f}",
                    "check_count": sum(bucket.label_counts.values()),
                    "label_counts": ";".join(f"{label}:{count}" for label, count in sorted(bucket.label_counts.items())),
                    "status_counts": ";".join(f"{status}:{count}" for status, count in sorted(bucket.status_counts.items())),
                    "mixed_labels": len(bucket.label_counts) > 1,
                }
            )


def write_direct_rows(path: Path, labels: Counter[str], statuses: Counter[str], samples: int) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["kind", "label", "count", "rate", "wilson95_lower", "wilson95_upper"])
        writer.writeheader()
        for label, count in sorted(labels.items(), key=lambda item: (-item[1], item[0])):
            low, high = wilson_interval(count, samples)
            writer.writerow(
                {
                    "kind": "label",
                    "label": label,
                    "count": count,
                    "rate": f"{count / samples:.8f}" if samples else "",
                    "wilson95_lower": f"{low:.8f}",
                    "wilson95_upper": f"{high:.8f}",
                }
            )
        for status, count in sorted(statuses.items(), key=lambda item: (-item[1], item[0])):
            low, high = wilson_interval(count, samples)
            writer.writerow(
                {
                    "kind": "status",
                    "label": status,
                    "count": count,
                    "rate": f"{count / samples:.8f}" if samples else "",
                    "wilson95_lower": f"{low:.8f}",
                    "wilson95_upper": f"{high:.8f}",
                }
            )


def write_json(path: Path, payload: dict[str, object]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


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


def print_summary(
    vertex_count: int,
    buckets: dict[str, Bucket],
    direct_labels: Counter[str],
    samples: int,
    out_dir: Path,
) -> None:
    classified_bucket_mass = sum(bucket.count for bucket in buckets.values() if bucket.label_counts)
    mixed = sum(1 for bucket in buckets.values() if len(bucket.label_counts) > 1)
    print(f"N={vertex_count}: observed buckets={len(buckets)}", flush=True)
    print(f"  classified high-volume bucket mass={classified_bucket_mass}/{samples}", flush=True)
    print(f"  buckets with mixed labels among checks={mixed}", flush=True)
    if direct_labels:
        print(f"  direct sample labels={dict(sorted(direct_labels.items()))}", flush=True)
    print(f"  wrote {out_dir}", flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
