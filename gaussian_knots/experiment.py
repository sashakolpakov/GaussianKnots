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
    edge_lengths,
    gaussian_polygon,
    is_numerically_embedded,
    minimum_nonadjacent_segment_distance,
)
from .pyknotid_adapter import KnotIdentification, identify_polygon, inspect_pyknotid_environment


def parse_vertex_counts(raw: str) -> tuple[int, ...]:
    values = tuple(int(part.strip()) for part in raw.split(",") if part.strip())
    if not values:
        raise ValueError("at least one vertex count is required")
    if any(value < 3 for value in values):
        raise ValueError("all vertex counts must be at least 3")
    return values


def run_experiment(
    vertex_counts: Iterable[int],
    samples: int,
    seed: int,
    output_dir: Path,
    use_fast: bool = True,
    allow_missing_pyknotid: bool = False,
    embedding_tolerance: float = 1e-9,
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
    _write_metadata(output_dir, vertex_counts, samples, seed, use_fast, environment)

    master_rng = np.random.default_rng(seed)
    summaries: list[dict[str, object]] = []

    for vertex_count in vertex_counts:
        records = []
        for sample_index in range(samples):
            sample_seed = int(master_rng.integers(0, np.iinfo(np.uint32).max))
            sample_rng = np.random.default_rng(sample_seed)
            vertices = gaussian_polygon(vertex_count, sample_rng)
            lengths = edge_lengths(vertices)
            min_nonadjacent_distance = (
                minimum_nonadjacent_segment_distance(vertices) if vertex_count >= 4 else float("inf")
            )
            identification = identify_polygon(vertices, use_fast=use_fast)
            record = _sample_record(
                vertex_count=vertex_count,
                sample_index=sample_index,
                sample_seed=sample_seed,
                vertices=vertices,
                lengths=lengths,
                min_nonadjacent_distance=min_nonadjacent_distance,
                embedded=is_numerically_embedded(vertices, tolerance=embedding_tolerance),
                identification=identification,
            )
            records.append(record)

        sample_path = output_dir / f"samples_N{vertex_count}.csv"
        _write_csv(sample_path, records)
        summary = _summarize(vertex_count, samples, sample_path, records)
        summaries.append(summary)

    _write_csv(output_dir / "summary.csv", summaries)
    return summaries


def _sample_record(
    vertex_count: int,
    sample_index: int,
    sample_seed: int,
    vertices: np.ndarray,
    lengths: np.ndarray,
    min_nonadjacent_distance: float,
    embedded: bool,
    identification: KnotIdentification,
) -> dict[str, object]:
    return {
        "N": vertex_count,
        "sample_index": sample_index,
        "sample_seed": sample_seed,
        "status": identification.status,
        "is_nontrivial": _optional_bool(identification.is_nontrivial),
        "knot_types": ";".join(identification.knot_types),
        "determinant": identification.determinant or "",
        "alexander_roots": json.dumps(identification.alexander_roots, sort_keys=True),
        "vassiliev_2": identification.vassiliev_2 or "",
        "vassiliev_3": identification.vassiliev_3 or "",
        "crossing_count": _optional_int(identification.crossing_count),
        "simplified_crossing_count": _optional_int(identification.simplified_crossing_count),
        "edge_length_min": float(np.min(lengths)),
        "edge_length_max": float(np.max(lengths)),
        "cycle_distortion": cycle_distortion(vertices),
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
        "sample_file": str(sample_path),
    }


def _write_csv(path: Path, records: list[dict[str, object]]) -> None:
    if not records:
        raise ValueError("cannot write empty CSV")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0].keys()))
        writer.writeheader()
        writer.writerows(records)


def _write_metadata(
    output_dir: Path,
    vertex_counts: Iterable[int],
    samples: int,
    seed: int,
    use_fast: bool,
    environment: object,
) -> None:
    metadata = {
        "vertex_counts": list(vertex_counts),
        "samples": samples,
        "seed": seed,
        "use_fast": use_fast,
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
