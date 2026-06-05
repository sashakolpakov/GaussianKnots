"""Shared runner for N=7 and N=8 wall/order-type probes."""

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
class SignatureRecord:
    signature: str
    kernel_basis: np.ndarray
    count: int = 0
    wall_signature: str = ""
    direct_labels: Counter[str] = field(default_factory=Counter)
    representative_label: str = "unclassified"
    representative_status: str = ""
    representative_determinant: str = ""
    representative_vassiliev_2: str = ""
    representative_vassiliev_3: str = ""


def main_for_vertex_count(vertex_count: int) -> int:
    args = build_parser(vertex_count).parse_args()
    if args.samples < 1:
        raise SystemExit("--samples must be positive")
    if args.no_classify:
        args.classify_reps = 0
        args.classify_samples = 0

    output_dir = args.output_dir or ROOT / "results" / f"n{vertex_count}_chamber_probe_{args.samples}"
    output_dir.mkdir(parents=True, exist_ok=True)

    model = build_chamber_model(vertex_count)
    symmetry_summary = analyze_symmetry(model, skip_full=args.skip_full_wall_symmetry, max_count=args.max_wall_automorphisms)
    representatives, sample_checks = sample_signatures(model, args)
    classify_representatives(model, representatives, args)
    fill_sample_matches(sample_checks, representatives)

    write_outputs(output_dir, model, symmetry_summary, representatives, sample_checks, args)
    print_report(model, symmetry_summary, representatives, sample_checks, args, output_dir)
    return 0


def build_parser(vertex_count: int) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            f"N={vertex_count} order-type/symmetry probe for the Stiefel projected-simplex knot model. "
            "This estimates signature sectors and checks direct-label consistency; it does not prove cell connectivity."
        )
    )
    parser.add_argument("--samples", type=int, default=5000, help="Haar kernel samples for signature-volume estimates")
    parser.add_argument("--seed", type=int, default=20260604, help="random seed")
    parser.add_argument(
        "--classify-reps",
        type=int,
        default=200,
        help="classify this many highest-mass D_N order-type representatives",
    )
    parser.add_argument(
        "--classify-samples",
        type=int,
        default=200,
        help="directly classify this many sampled polygons for consistency checks",
    )
    parser.add_argument("--no-classify", action="store_true", help="skip all pyknotid classification")
    parser.add_argument("--no-fast", action="store_true", help="do not request pyknotid fast helpers")
    parser.add_argument(
        "--skip-full-wall-symmetry",
        action="store_true",
        help="skip wall-normal Gram automorphism enumeration",
    )
    parser.add_argument(
        "--max-wall-automorphisms",
        type=int,
        default=100000,
        help="stop wall automorphism enumeration after this many automorphisms",
    )
    parser.add_argument("--output-dir", type=Path, default=None, help="output directory")
    return parser


def analyze_symmetry(model: ChamberModel, skip_full: bool, max_count: int) -> dict[str, object]:
    projective_gram = np.abs(model.wall_gram)
    d_wall_perms = set(dihedral_wall_permutations(model))
    d_valid = all(is_projective_automorphism(model.wall_gram, permutation) for permutation in d_wall_perms)
    automorphisms = None
    truncated = False
    if not skip_full:
        automorphisms = wall_gram_automorphisms(model.wall_gram, max_count=max_count)
        truncated = len(automorphisms) >= max_count
    full_set = set(automorphisms or [])
    return {
        "vertex_count": model.vertex_count,
        "h_dimension": model.h_dimension,
        "kernel_dimension": model.kernel_dimension,
        "wall_count": len(model.edge_pairs),
        "order_type_sign_count": len(model.quads),
        "wall_gram_values": sorted({round(float(value), 12) for value in projective_gram.ravel()}),
        "dihedral_subgroup_order": len(d_wall_perms),
        "dihedral_preserves_wall_gram": d_valid,
        "wall_gram_group_order": len(automorphisms) if automorphisms is not None else None,
        "wall_gram_group_truncated": truncated,
        "extra_wall_gram_symmetries": len(full_set - d_wall_perms) if automorphisms is not None else None,
    }


def sample_signatures(
    model: ChamberModel,
    args: argparse.Namespace,
) -> tuple[dict[str, SignatureRecord], list[dict[str, object]]]:
    rng = np.random.default_rng(args.seed)
    representatives: dict[str, SignatureRecord] = {}
    sample_checks: list[dict[str, object]] = []
    environment = inspect_pyknotid_environment() if args.classify_samples > 0 else None

    for sample_index in range(args.samples):
        kernel = random_kernel(model, rng)
        vertices = vertices_from_kernel(model, kernel)
        signature = orientation_signature(vertices, model.quads)
        representative_signature = d_orbit_representative(
            signature,
            model.vertex_count,
            model.quads,
            model.quad_index,
        )
        record = representatives.get(representative_signature)
        if record is None:
            record = SignatureRecord(
                signature=representative_signature,
                kernel_basis=kernel,
                wall_signature=determinant_wall_signature(model.wall_data, kernel),
            )
            representatives[representative_signature] = record
        record.count += 1

        if sample_index < args.classify_samples:
            label, status = classify_vertices(vertices, use_fast=not args.no_fast, environment_available=environment.available if environment else False)
            record.direct_labels[label] += 1
            sample_checks.append(
                {
                    "sample_index": sample_index,
                    "representative_signature": representative_signature,
                    "direct_label": label,
                    "direct_status": status,
                    "representative_label": "",
                    "match": "",
                }
            )

    return representatives, sample_checks


def classify_representatives(
    model: ChamberModel,
    representatives: dict[str, SignatureRecord],
    args: argparse.Namespace,
) -> None:
    if args.classify_reps <= 0:
        return
    environment = inspect_pyknotid_environment()
    if not environment.available:
        for record in representatives.values():
            record.representative_status = f"pyknotid_unavailable:{environment.error}"
        return

    for record in sorted(representatives.values(), key=lambda item: (-item.count, item.signature))[: args.classify_reps]:
        vertices = vertices_from_kernel(model, record.kernel_basis)
        identification = identify_polygon(vertices, use_fast=not args.no_fast)
        record.representative_label = knot_label(identification)
        record.representative_status = identification.status
        record.representative_determinant = identification.determinant or ""
        record.representative_vassiliev_2 = identification.vassiliev_2 or ""
        record.representative_vassiliev_3 = identification.vassiliev_3 or ""


def classify_vertices(vertices: np.ndarray, use_fast: bool, environment_available: bool) -> tuple[str, str]:
    if not environment_available:
        return "unknown", "pyknotid_unavailable"
    identification = identify_polygon(vertices, use_fast=use_fast)
    return knot_label(identification), identification.status


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


def fill_sample_matches(sample_checks: list[dict[str, object]], representatives: dict[str, SignatureRecord]) -> None:
    for row in sample_checks:
        record = representatives[str(row["representative_signature"])]
        representative_label = record.representative_label
        row["representative_label"] = representative_label
        if representative_label == "unclassified":
            row["match"] = ""
        else:
            row["match"] = str(row["direct_label"] == representative_label)


def write_outputs(
    output_dir: Path,
    model: ChamberModel,
    symmetry_summary: dict[str, object],
    representatives: dict[str, SignatureRecord],
    sample_checks: list[dict[str, object]],
    args: argparse.Namespace,
) -> None:
    write_json(
        output_dir / "symmetry_summary.json",
        {
            **symmetry_summary,
            "samples": args.samples,
            "seed": args.seed,
            "classify_reps": args.classify_reps,
            "classify_samples": args.classify_samples,
        },
    )
    write_signature_summary(output_dir / "signature_summary.csv", model, representatives, args.samples)
    write_sample_checks(output_dir / "sample_checks.csv", sample_checks)


def write_json(path: Path, payload: dict[str, object]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_signature_summary(
    path: Path,
    model: ChamberModel,
    representatives: dict[str, SignatureRecord],
    samples: int,
) -> None:
    rows = []
    for rank, record in enumerate(sorted(representatives.values(), key=lambda item: (-item.count, item.signature)), start=1):
        direct_mixed = len([label for label, count in record.direct_labels.items() if count > 0]) > 1
        comparable = record.representative_label != "unclassified"
        mismatches = sum(
            count for label, count in record.direct_labels.items() if comparable and label != record.representative_label
        )
        low, high = wilson_interval(record.count, samples)
        rows.append(
            {
                "rank": rank,
                "representative_signature": record.signature,
                "wall_signature": record.wall_signature,
                "d_orbit_size": len(d_orbit(record.signature, model.vertex_count, model.quads, model.quad_index)),
                "count": record.count,
                "rate": f"{record.count / samples:.8f}",
                "wilson95_lower": f"{low:.8f}",
                "wilson95_upper": f"{high:.8f}",
                "representative_label": record.representative_label,
                "representative_status": record.representative_status,
                "determinant": record.representative_determinant,
                "vassiliev_2": record.representative_vassiliev_2,
                "vassiliev_3": record.representative_vassiliev_3,
                "direct_sample_count": sum(record.direct_labels.values()),
                "direct_label_counts": ";".join(f"{label}:{count}" for label, count in sorted(record.direct_labels.items())),
                "direct_mixed": direct_mixed,
                "direct_rep_mismatches": mismatches,
            }
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["rank"])
        writer.writeheader()
        writer.writerows(rows)


def write_sample_checks(path: Path, sample_checks: list[dict[str, object]]) -> None:
    fieldnames = [
        "sample_index",
        "representative_signature",
        "direct_label",
        "direct_status",
        "representative_label",
        "match",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sample_checks)


def print_report(
    model: ChamberModel,
    symmetry_summary: dict[str, object],
    representatives: dict[str, SignatureRecord],
    sample_checks: list[dict[str, object]],
    args: argparse.Namespace,
    output_dir: Path,
) -> None:
    print(f"N={model.vertex_count}")
    print(f"H dimension: {model.h_dimension}")
    print(f"kernel dimension: {model.kernel_dimension}")
    print(f"edge-pair walls: {len(model.edge_pairs)}")
    print(f"order-type signs: {len(model.quads)}")
    print(f"wall-normal Gram values: {symmetry_summary['wall_gram_values']}")
    if symmetry_summary["wall_gram_group_order"] is None:
        print("wall-normal Gram group order: skipped")
    else:
        suffix = " (truncated)" if symmetry_summary["wall_gram_group_truncated"] else ""
        print(f"wall-normal Gram group order: {symmetry_summary['wall_gram_group_order']}{suffix}")
        print(f"extra Gram symmetries beyond D_N: {symmetry_summary['extra_wall_gram_symmetries']}")
    print(f"D_N wall subgroup order: {symmetry_summary['dihedral_subgroup_order']}")
    print(f"D_N preserves wall Gram: {symmetry_summary['dihedral_preserves_wall_gram']}")
    print(f"sampled Haar points: {args.samples}")
    print(f"observed D_N order-type buckets: {len(representatives)}")

    classified_records = [record for record in representatives.values() if record.representative_label != "unclassified"]
    classified_mass = sum(record.count for record in classified_records)
    print(f"classified representative buckets: {len(classified_records)}")
    print(f"classified representative sample mass: {classified_mass}/{args.samples}")
    if classified_mass < args.samples:
        print("classified label masses below omit unclassified representative buckets")
    for label, count in label_mass(classified_records).most_common():
        low, high = wilson_interval(count, args.samples)
        print(f"  {label}: {count}/{args.samples} rate={count / args.samples:.6f} Wilson95=({low:.6f}, {high:.6f})")

    mixed = [record for record in representatives.values() if len(record.direct_labels) > 1]
    comparable_checks = [row for row in sample_checks if row["match"] != ""]
    mismatches = [row for row in comparable_checks if row["match"] == "False"]
    direct_labels = Counter(str(row["direct_label"]) for row in sample_checks)
    print(f"directly classified sample checks: {len(sample_checks)}")
    if direct_labels:
        print(f"direct sample label counts: {dict(sorted(direct_labels.items()))}")
    print(f"buckets with mixed direct labels: {len(mixed)}")
    print(f"direct-vs-representative mismatches: {len(mismatches)}/{len(comparable_checks)}")
    print(f"wrote outputs to {output_dir}")


def label_mass(records: list[SignatureRecord]) -> Counter[str]:
    counter: Counter[str] = Counter()
    for record in records:
        counter[record.representative_label] += record.count
    return counter
