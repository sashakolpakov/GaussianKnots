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
    edge_lengths,
    gaussian_polygon,
    is_numerically_embedded,
    minimum_nonadjacent_segment_distance,
)
from gaussian_knots.pyknotid_adapter import identify_polygon, inspect_pyknotid_environment


def main() -> int:
    rng = np.random.default_rng(12345)
    vertices = gaussian_polygon(6, rng)
    lengths = edge_lengths(vertices)

    assert vertices.shape == (6, 3)
    assert np.all(lengths > 0.0)
    assert cycle_distortion(vertices) >= 1.0
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

