#!/usr/bin/env python3
"""Summarize Monte Carlo knot-type probability estimates from GaussianKnots outputs."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


Z_95 = 1.959963984540054


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("result_dir", type=Path, help="GaussianKnots result directory")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="CSV output path; default is <result_dir>/probability_estimates.csv",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    result_dir = args.result_dir
    output = args.output or result_dir / "probability_estimates.csv"

    sample_counts = read_sample_counts(result_dir / "summary.csv")
    rows = summarize_type_counts(result_dir / "type_counts.csv", sample_counts)
    write_csv(output, rows)
    print(f"wrote {output}")
    for row in rows:
        print(
            f"N={row['N']} {row['knot_label']} "
            f"{row['count']}/{row['samples']} "
            f"p_hat={row['p_hat']} CI={row['wilson_95_low']}--{row['wilson_95_high']}"
        )
    return 0


def read_sample_counts(path: Path) -> dict[int, int]:
    counts: dict[int, int] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            counts[int(row["N"])] = int(row["samples"])
    return counts


def summarize_type_counts(path: Path, sample_counts: dict[int, int]) -> list[dict[str, object]]:
    buckets: dict[tuple[int, str, str], int] = defaultdict(int)
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            vertex_count = int(row["N"])
            count = int(row["count"])
            label, status = normalize_label(
                label=row["knot_label"],
                nontrivial=int(row["nontrivial"]),
                trivial=int(row["trivial"]),
                unknown=int(row["unknown"]),
            )
            buckets[(vertex_count, label, status)] += count

    rows: list[dict[str, object]] = []
    for (vertex_count, label, status), count in sorted(
        buckets.items(),
        key=lambda item: (item[0][0], status_order(item[0][2]), -item[1], item[0][1]),
    ):
        samples = sample_counts[vertex_count]
        low, high = wilson_interval(count, samples)
        rows.append(
            {
                "N": vertex_count,
                "knot_label": label,
                "status": status,
                "count": count,
                "samples": samples,
                "p_hat": f"{count / samples:.6f}",
                "wilson_95_low": f"{low:.6f}",
                "wilson_95_high": f"{high:.6f}",
            }
        )
    return rows


def normalize_label(label: str, nontrivial: int, trivial: int, unknown: int) -> tuple[str, str]:
    if unknown:
        return "unknown_or_ambiguous", "unknown"
    if trivial:
        return "0_1", "classified"
    if nontrivial and ";" in label:
        return "ambiguous_nontrivial_candidate_list", "ambiguous_nontrivial"
    if nontrivial:
        return label, "classified"
    return label or "unknown_or_ambiguous", "unknown"


def status_order(status: str) -> int:
    if status == "classified":
        return 0
    if status == "ambiguous_nontrivial":
        return 1
    return 2


def wilson_interval(successes: int, trials: int) -> tuple[float, float]:
    if trials <= 0:
        return (math.nan, math.nan)
    phat = successes / trials
    denominator = 1.0 + Z_95 * Z_95 / trials
    center = (phat + Z_95 * Z_95 / (2.0 * trials)) / denominator
    half_width = (
        Z_95
        * math.sqrt((phat * (1.0 - phat) + Z_95 * Z_95 / (4.0 * trials)) / trials)
        / denominator
    )
    return max(0.0, center - half_width), min(1.0, center + half_width)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "N",
                "knot_label",
                "status",
                "count",
                "samples",
                "p_hat",
                "wilson_95_low",
                "wilson_95_high",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    raise SystemExit(main())
