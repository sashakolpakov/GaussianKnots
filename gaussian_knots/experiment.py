"""Experiment runner for Gaussian projected stick knots."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import Iterable

import numpy as np

from .generation import (
    cycle_distortion,
    distance_deformation_stats,
    edge_lengths,
    is_numerically_embedded,
    minimum_nonadjacent_segment_distance,
    projected_simplex_polygon,
)
from .pyknotid_adapter import KnotIdentification, identify_polygon, inspect_pyknotid_environment


def parse_vertex_counts(raw: str) -> tuple[int, ...]:
    values = tuple(int(part.strip()) for part in raw.split(",") if part.strip())
    if not values:
        raise ValueError("at least one vertex count is required")
    if any(value < 5 for value in values):
        raise ValueError("all vertex counts must be at least 5 for the simplex experiment")
    return values


def run_experiment(
    vertex_counts: Iterable[int],
    samples: int,
    seed: int,
    output_dir: Path,
    use_fast: bool = True,
    allow_missing_pyknotid: bool = False,
    embedding_tolerance: float = 1e-9,
    projection_model: str = "haar",
) -> list[dict[str, object]]:
    """Run the experiment and write one CSV per N plus a summary CSV."""

    vertex_counts = tuple(vertex_counts)
    if samples < 1:
        raise ValueError("samples must be positive")

    environment = inspect_pyknotid_environment()
    if not environment.available and not allow_missing_pyknotid:
        raise RuntimeError(
            "pyknotid is not importable. Install requirements.txt, or rerun with "
            "--allow-missing-pyknotid to generate polygons without identification. "
            f"Import error: {environment.error}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    _write_metadata(output_dir, vertex_counts, samples, seed, use_fast, environment, projection_model)

    master_rng = np.random.default_rng(seed)
    summaries: list[dict[str, object]] = []

    for vertex_count in vertex_counts:
        records = []
        for sample_index in range(samples):
            sample_seed = int(master_rng.integers(0, np.iinfo(np.uint32).max))
            sample_rng = np.random.default_rng(sample_seed)
            vertices = projected_simplex_polygon(vertex_count, sample_rng, projection_model=projection_model)
            lengths = edge_lengths(vertices)
            deformation = distance_deformation_stats(vertices)
            min_nonadjacent_distance = (
                minimum_nonadjacent_segment_distance(vertices) if vertex_count >= 4 else float("inf")
            )
            identification = identify_polygon(vertices, use_fast=use_fast)
            record = _sample_record(
                vertex_count=vertex_count,
                sample_index=sample_index,
                sample_seed=sample_seed,
                projection_model=projection_model,
                vertices=vertices,
                lengths=lengths,
                deformation=deformation,
                min_nonadjacent_distance=min_nonadjacent_distance,
                embedded=is_numerically_embedded(vertices, tolerance=embedding_tolerance),
                identification=identification,
            )
            records.append(record)

        sample_path = output_dir / f"samples_N{vertex_count}.csv"
        _write_csv(sample_path, records)
        _write_csv(output_dir / f"type_counts_N{vertex_count}.csv", _type_counts(vertex_count, samples, records))
        summary = _summarize(vertex_count, samples, sample_path, records)
        summaries.append(summary)

    _write_csv(output_dir / "summary.csv", summaries)
    _write_csv(output_dir / "type_counts.csv", _all_type_counts(samples, summaries, output_dir))
    return summaries


def _sample_record(
    vertex_count: int,
    sample_index: int,
    sample_seed: int,
    projection_model: str,
    vertices: np.ndarray,
    lengths: np.ndarray,
    deformation: dict[str, float],
    min_nonadjacent_distance: float,
    embedded: bool,
    identification: KnotIdentification,
) -> dict[str, object]:
    return {
        "N": vertex_count,
        "sample_index": sample_index,
        "sample_seed": sample_seed,
        "projection_model": projection_model,
        "status": identification.status,
        "is_nontrivial": _optional_bool(identification.is_nontrivial),
        "knot_label": _knot_label(identification),
        "knot_types": ";".join(identification.knot_types),
        "determinant": identification.determinant or "",
        "alexander_roots": json.dumps(identification.alexander_roots, sort_keys=True),
        "vassiliev_2": identification.vassiliev_2 or "",
        "vassiliev_3": identification.vassiliev_3 or "",
        "gauss_code": identification.gauss_code or "",
        "simplified_gauss_code": identification.simplified_gauss_code or "",
        "crossing_count": _optional_int(identification.crossing_count),
        "simplified_crossing_count": _optional_int(identification.simplified_crossing_count),
        "edge_length_min": float(np.min(lengths)),
        "edge_length_mean": float(np.mean(lengths)),
        "edge_length_rms": float(np.sqrt(np.mean(lengths * lengths))),
        "edge_length_max": float(np.max(lengths)),
        "cycle_distortion": cycle_distortion(vertices),
        **deformation,
        "min_nonadjacent_segment_distance": min_nonadjacent_distance,
        "embedded_numerically": embedded,
        "fast_mode_requested": identification.fast_mode_requested,
        "fast_backend_available": _optional_bool(identification.fast_backend_available),
        "messages": " | ".join(identification.messages),
    }


def _summarize(
    vertex_count: int,
    samples: int,
    sample_path: Path,
    records: list[dict[str, object]],
) -> dict[str, object]:
    nontrivial = sum(1 for record in records if record["is_nontrivial"] == "true")
    trivial = sum(1 for record in records if record["is_nontrivial"] == "false")
    unknown = samples - nontrivial - trivial
    known = nontrivial + trivial
    numeric = _numeric_columns(records)
    return {
        "N": vertex_count,
        "samples": samples,
        "classified": known,
        "nontrivial": nontrivial,
        "trivial": trivial,
        "unknown": unknown,
        "nontrivial_rate_known": nontrivial / known if known else "",
        "nontrivial_lower_bound_rate": nontrivial / samples,
        "unknown_rate": unknown / samples,
        "cycle_distortion_mean": numeric["cycle_distortion_mean"],
        "cycle_distortion_median": numeric["cycle_distortion_median"],
        "pair_distance_distortion_mean": numeric["pair_distance_distortion_mean"],
        "pair_distance_distortion_median": numeric["pair_distance_distortion_median"],
        "pair_abs_ratio_rms_mean": numeric["pair_abs_ratio_rms_mean"],
        "pair_rms_scale_to_simplex_mean": numeric["pair_rms_scale_to_simplex_mean"],
        "pair_normalized_ratio_min_mean": numeric["pair_normalized_ratio_min_mean"],
        "pair_normalized_ratio_max_mean": numeric["pair_normalized_ratio_max_mean"],
        "sample_file": str(sample_path),
    }


def _type_counts(vertex_count: int, samples: int, records: list[dict[str, object]]) -> list[dict[str, object]]:
    counts: dict[str, dict[str, int]] = {}
    for record in records:
        label = str(record["knot_label"])
        bucket = counts.setdefault(label, {"count": 0, "nontrivial": 0, "trivial": 0, "unknown": 0})
        bucket["count"] += 1
        if record["is_nontrivial"] == "true":
            bucket["nontrivial"] += 1
        elif record["is_nontrivial"] == "false":
            bucket["trivial"] += 1
        else:
            bucket["unknown"] += 1

    rows = []
    for label, bucket in sorted(counts.items(), key=lambda item: (-item[1]["count"], item[0])):
        rows.append(
            {
                "N": vertex_count,
                "knot_label": label,
                "count": bucket["count"],
                "rate": bucket["count"] / samples,
                "nontrivial": bucket["nontrivial"],
                "trivial": bucket["trivial"],
                "unknown": bucket["unknown"],
            }
        )
    return rows


def _all_type_counts(
    samples: int,
    summaries: list[dict[str, object]],
    output_dir: Path,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for summary in summaries:
        path = output_dir / f"type_counts_N{summary['N']}.csv"
        with path.open(newline="", encoding="utf-8") as handle:
            rows.extend(csv.DictReader(handle))
    if not rows:
        return [{"N": "", "knot_label": "", "count": "", "rate": "", "nontrivial": "", "trivial": "", "unknown": ""}]
    return rows


def _write_csv(path: Path, records: list[dict[str, object]]) -> None:
    if not records:
        raise ValueError("cannot write empty CSV")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0].keys()))
        writer.writeheader()
        writer.writerows(records)


def _numeric_columns(records: list[dict[str, object]]) -> dict[str, float]:
    columns = {
        "cycle_distortion": [],
        "pair_distance_distortion": [],
        "pair_abs_ratio_rms": [],
        "pair_rms_scale_to_simplex": [],
        "pair_normalized_ratio_min": [],
        "pair_normalized_ratio_max": [],
    }
    for record in records:
        for key in columns:
            columns[key].append(float(record[key]))

    values: dict[str, float] = {}
    for key, column_values in columns.items():
        array = np.asarray(column_values, dtype=float)
        values[f"{key}_mean"] = float(np.mean(array))
        values[f"{key}_median"] = float(np.median(array))
    return values


def _write_metadata(
    output_dir: Path,
    vertex_counts: Iterable[int],
    samples: int,
    seed: int,
    use_fast: bool,
    environment: object,
    projection_model: str,
) -> None:
    metadata = {
        "vertex_counts": list(vertex_counts),
        "samples": samples,
        "seed": seed,
        "use_fast": use_fast,
        "projection_model": projection_model,
        "pyknotid_environment": asdict(environment),
    }
    with (output_dir / "run_metadata.json").open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")


def _optional_bool(value: bool | None) -> str:
    if value is None:
        return ""
    return "true" if value else "false"


def _optional_int(value: int | None) -> str:
    if value is None:
        return ""
    return str(value)


def _knot_label(identification: KnotIdentification) -> str:
    if identification.knot_types:
        return ";".join(identification.knot_types)
    if identification.is_nontrivial is False:
        return "unknot"
    if identification.is_nontrivial is True:
        return "nontrivial_detected"
    return "unknown"
