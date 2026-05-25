#!/usr/bin/env python3
"""Small self-check for the GaussianKnots package."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR))

from gaussian_knots.generation import (
    cycle_distortion,
    distance_deformation_stats,
    edge_lengths,
    gaussian_polygon,
    haar_projected_simplex_polygon,
    is_numerically_embedded,
    minimum_nonadjacent_segment_distance,
    projected_simplex_polygon,
)
from gaussian_knots.pyknotid_adapter import identify_polygon, inspect_pyknotid_environment


def main() -> int:
    rng = np.random.default_rng(12345)
    vertices = projected_simplex_polygon(6, rng, projection_model="haar")
    lengths = edge_lengths(vertices)

    assert vertices.shape == (6, 3)
    assert np.allclose(vertices.sum(axis=0), 0.0)
    assert np.allclose(vertices.T @ vertices, np.eye(3))
    assert np.all(lengths > 0.0)
    assert cycle_distortion(vertices) >= 1.0
    deformation = distance_deformation_stats(vertices)
    assert deformation["pair_distance_distortion"] >= 1.0
    assert deformation["pair_normalized_ratio_min"] <= 1.0
    assert deformation["pair_normalized_ratio_max"] >= 1.0
    assert minimum_nonadjacent_segment_distance(vertices) >= 0.0
    assert isinstance(is_numerically_embedded(vertices), bool)

    environment = inspect_pyknotid_environment()
    identification = identify_polygon(vertices, use_fast=True)
    assert identification.status
    if not environment.available:
        assert identification.is_nontrivial is None
        print("generator checks passed; pyknotid is not installed, adapter fallback passed")
    else:
        print(
            "generator checks passed; pyknotid status="
            f"{identification.status}, fast_backend={identification.fast_backend_available}"
        )

    gaussian_vertices = gaussian_polygon(6, np.random.default_rng(23456))
    haar_vertices = haar_projected_simplex_polygon(6, np.random.default_rng(34567))
    assert gaussian_vertices.shape == haar_vertices.shape == (6, 3)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
