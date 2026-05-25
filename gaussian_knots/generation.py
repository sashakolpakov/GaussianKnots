"""Sampling and geometric checks for Gaussian Hamiltonian-cycle polygons."""

from __future__ import annotations

import math
from typing import Iterable, Iterator, Tuple

import numpy as np

DEFAULT_TARGET_DIM = 3


def gaussian_polygon(
    vertex_count: int,
    rng: np.random.Generator,
    target_dim: int = DEFAULT_TARGET_DIM,
    scale: float | None = None,
) -> np.ndarray:
    """Sample columns of a Gaussian projection matrix as a closed polygon.

    The paper's projection convention uses entries distributed as N(0, 1/n)
    for target dimension n.  For the knot experiment n=3, so the default
    standard deviation is 1/sqrt(3).  Multiplying all coordinates by a positive
    constant does not change knot type, but keeping the convention makes the
    generated data line up with the metric sections of the paper.
    """

    if vertex_count < 3:
        raise ValueError("a polygon needs at least 3 vertices")
    if target_dim != 3:
        raise ValueError("knot identification requires target_dim=3")
    if scale is None:
        scale = 1.0 / math.sqrt(target_dim)
    return rng.normal(loc=0.0, scale=scale, size=(vertex_count, target_dim))


def cycle_edges(vertex_count: int) -> Iterator[Tuple[int, int]]:
    """Yield edges in the Hamiltonian cycle 0-1-...-(N-1)-0."""

    for index in range(vertex_count):
        yield index, (index + 1) % vertex_count


def edge_lengths(vertices: np.ndarray) -> np.ndarray:
    """Return Euclidean lengths of consecutive cycle edges."""

    points = _as_vertices(vertices)
    diffs = np.roll(points, -1, axis=0) - points
    return np.linalg.norm(diffs, axis=1)


def cycle_distortion(vertices: np.ndarray) -> float:
    """Return max edge length divided by min edge length."""

    lengths = edge_lengths(vertices)
    minimum = float(np.min(lengths))
    if minimum <= 0.0:
        return math.inf
    return float(np.max(lengths) / minimum)


def minimum_nonadjacent_segment_distance(vertices: np.ndarray) -> float:
    """Return the minimum distance between nonconsecutive cycle edges."""

    points = _as_vertices(vertices)
    vertex_count = len(points)
    minimum = math.inf
    for first in range(vertex_count):
        p0 = points[first]
        p1 = points[(first + 1) % vertex_count]
        for second in range(first + 1, vertex_count):
            if _edges_are_adjacent(first, second, vertex_count):
                continue
            q0 = points[second]
            q1 = points[(second + 1) % vertex_count]
            distance = segment_segment_distance(p0, p1, q0, q1)
            if distance < minimum:
                minimum = distance
    return minimum


def is_numerically_embedded(vertices: np.ndarray, tolerance: float = 1e-9) -> bool:
    """Return whether the polygon passes simple numerical embedding checks."""

    lengths = edge_lengths(vertices)
    if np.min(lengths) <= tolerance:
        return False
    if len(vertices) < 4:
        return True
    return minimum_nonadjacent_segment_distance(vertices) > tolerance


def segment_segment_distance(
    p0: Iterable[float],
    p1: Iterable[float],
    q0: Iterable[float],
    q1: Iterable[float],
) -> float:
    """Distance between two line segments in R^3.

    This is the standard closest-points calculation, with explicit handling of
    near-degenerate segments.  It is only used as a numerical sanity check; the
    Gaussian model is embedded with probability one.
    """

    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)
    q0 = np.asarray(q0, dtype=float)
    q1 = np.asarray(q1, dtype=float)

    u = p1 - p0
    v = q1 - q0
    w = p0 - q0
    a = float(np.dot(u, u))
    b = float(np.dot(u, v))
    c = float(np.dot(v, v))
    d = float(np.dot(u, w))
    e = float(np.dot(v, w))
    eps = 1e-14

    if a <= eps and c <= eps:
        return float(np.linalg.norm(p0 - q0))
    if a <= eps:
        t = _clip01(e / c)
        return float(np.linalg.norm(p0 - (q0 + t * v)))
    if c <= eps:
        s = _clip01(-d / a)
        return float(np.linalg.norm((p0 + s * u) - q0))

    denominator = a * c - b * b
    if denominator <= eps:
        s = 0.0
    else:
        s = _clip01((b * e - c * d) / denominator)
    t = (b * s + e) / c

    if t < 0.0:
        t = 0.0
        s = _clip01(-d / a)
    elif t > 1.0:
        t = 1.0
        s = _clip01((b - d) / a)

    closest_p = p0 + s * u
    closest_q = q0 + t * v
    return float(np.linalg.norm(closest_p - closest_q))


def _as_vertices(vertices: np.ndarray) -> np.ndarray:
    points = np.asarray(vertices, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("vertices must have shape (N, 3)")
    if points.shape[0] < 3:
        raise ValueError("a polygon needs at least 3 vertices")
    return points


def _edges_are_adjacent(first: int, second: int, vertex_count: int) -> bool:
    return (
        first == second
        or (first + 1) % vertex_count == second
        or (second + 1) % vertex_count == first
    )


def _clip01(value: float) -> float:
    return max(0.0, min(1.0, value))

