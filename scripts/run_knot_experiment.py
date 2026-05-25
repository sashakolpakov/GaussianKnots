#!/usr/bin/env python3
"""Run the Gaussian projected Hamiltonian-cycle knot experiment."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR))

from gaussian_knots.experiment import parse_vertex_counts, run_experiment
from gaussian_knots.pyknotid_adapter import inspect_pyknotid_environment


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--vertices",
        default="6,7,8,10",
        help="comma-separated vertex/stick counts to sample",
    )
    parser.add_argument("--samples", type=int, default=100, help="samples per vertex count")
    parser.add_argument("--seed", type=int, default=20260524, help="master random seed")
    parser.add_argument(
        "--projection-model",
        choices=("haar", "gaussian"),
        default="haar",
        help="projection model: Haar row-orthonormal simplex projection, or raw Gaussian coordinates",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=PROJECT_DIR / "results" / "knot_experiment",
        help="directory for CSV and metadata outputs",
    )
    parser.add_argument(
        "--allow-missing-pyknotid",
        action="store_true",
        help="write generated-polygon records even when pyknotid is unavailable",
    )
    parser.add_argument(
        "--no-fast",
        action="store_true",
        help="do not request pyknotid Cython/fast helpers where the installed API supports them",
    )
    parser.add_argument(
        "--embedding-tolerance",
        type=float,
        default=1e-9,
        help="tolerance for numerical segment-intersection sanity checks",
    )
    parser.add_argument(
        "--download-pyknotid-db",
        action="store_true",
        help="attempt to download the optional pyknotid catalogue database before running",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    try:
        vertex_counts = parse_vertex_counts(args.vertices)
    except ValueError as exc:
        parser.error(str(exc))

    if args.download_pyknotid_db:
        try:
            from pyknotid.catalogue.getdb import download_database

            download_database()
        except Exception as exc:
            parser.exit(2, f"Could not download pyknotid database: {exc}\n")

    environment = inspect_pyknotid_environment()
    if environment.available:
        print(
            "pyknotid available"
            f" version={environment.version or 'unknown'}"
            f" fast_backend={environment.fast_backend_available}"
        )
    else:
        print(f"pyknotid unavailable: {environment.error}", file=sys.stderr)

    try:
        summaries = run_experiment(
            vertex_counts=vertex_counts,
            samples=args.samples,
            seed=args.seed,
            output_dir=args.output_dir,
            use_fast=not args.no_fast,
            allow_missing_pyknotid=args.allow_missing_pyknotid,
            embedding_tolerance=args.embedding_tolerance,
            projection_model=args.projection_model,
        )
    except Exception as exc:
        parser.exit(2, f"{exc}\n")

    print(f"wrote outputs to {args.output_dir}")
    for summary in summaries:
        known_rate = summary["nontrivial_rate_known"]
        known_rate_text = "" if known_rate == "" else f"{known_rate:.6f}"
        print(
            f"N={summary['N']} samples={summary['samples']} "
            f"classified={summary['classified']} nontrivial={summary['nontrivial']} "
            f"trivial={summary['trivial']} unknown={summary['unknown']} "
            f"rate_known={known_rate_text} "
            f"lower_bound={summary['nontrivial_lower_bound_rate']:.6f}"
        )
    print(f"type counts written to {args.output_dir / 'type_counts.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
