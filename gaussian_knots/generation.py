"""Sampling and geometric checks for projected-simplex Hamiltonian cycles."""

from __future__ import annotations

import math
from typing import Iterable, Iterator, Tuple

import numpy as np

DEFAULT_TARGET_DIM = 3
SIMPLEX_DISTANCE = math.sqrt(2.0)
PROJECTION_MODELS = ("haar", "gaussian")


def projected_simplex_polygon(
    vertex_count: int,
    rng: np.random.Generator,
    projection_model: str = "haar",
) -> np.ndarray:
    """Sample the Hamiltonian cycle through projected simplex vertices.

    The manuscript model starts with simplex vertices e_1,...,e_N in R^N and
    projects them to R^3.  In the ``haar`` model this function samples a
    Haar-distributed 3-plane in the centered coordinate hyperplane 1^perp and
    returns the projected vertices in the cyclic order 1,2,...,N,1.

    The ``gaussian`` model samples a raw 3 x N Gaussian projection matrix.  Its
    knot-type law agrees with the Haar row-space model almost surely after
    centering and applying the polar row factor, but its metric edge lengths are
    sheared by the random singular values.
    """

    if projection_model == "haar":
        return haar_projected_simplex_polygon(vertex_count, rng)
    if projection_model == "gaussian":
        return gaussian_polygon(vertex_count, rng)
    raise ValueError(f"unknown projection model {projection_model!r}; expected one of {PROJECTION_MODELS}")


def haar_projected_simplex_polygon(vertex_count: int, rng: np.random.Generator) -> np.ndarray:
    """Project N simplex vertices by a Haar-distributed 3-plane.

    The returned array has shape (N, 3); row i is the projected point P e_i.
    The columns are centered and orthonormal:

        vertices.T @ vertices = I_3,   sum_i vertices[i] = 0.

    This is a concrete coordinate representative of a Haar point in
    Gr(3, 1^perp) or, before quotienting by target rotations, in the Stiefel
    manifold of row-orthonormal projections.
    """

    if vertex_count < 5:
        raise ValueError("the simplex experiment expects at least 5 vertices")

    gaussian_rows = rng.normal(size=(DEFAULT_TARGET_DIM, vertex_count))
    centered_rows = gaussian_rows - gaussian_rows.mean(axis=1, keepdims=True)
    q_matrix, _ = np.linalg.qr(centered_rows.T, mode="reduced")

    # QR is unique only up to signs.  The signs do not affect knot type, but
    # normalizing them makes seeded runs deterministic across LAPACK variants.
    for column_index in range(q_matrix.shape[1]):
        pivot = int(np.argmax(np.abs(q_matrix[:, column_index])))
        if q_matrix[pivot, column_index] < 0:
            q_matrix[:, column_index] *= -1.0
    return q_matrix


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

    if vertex_count < 5:
        raise ValueError("the simplex experiment expects at least 5 vertices")
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


def pairwise_distances(vertices: np.ndarray) -> np.ndarray:
    """Return all projected distances between distinct simplex vertices."""

    points = _as_vertices(vertices)
    values = []
    for first in range(len(points)):
        diffs = points[first + 1 :] - points[first]
        values.extend(np.linalg.norm(diffs, axis=1))
    return np.asarray(values, dtype=float)


def distance_deformation_stats(vertices: np.ndarray) -> dict[str, float]:
    """Return all-pair distance deformation statistics.

    The original regular simplex has pair distance sqrt(2) for every pair.
    Absolute ratios measure the literal metric deformation.  Normalized ratios
    divide projected distances by their sample RMS, removing one global scale;
    this isolates shape deformation and lets Haar and raw Gaussian projections
    be compared on the same footing.
    """

    distances = pairwise_distances(vertices)
    minimum = float(np.min(distances))
    maximum = float(np.max(distances))
    mean = float(np.mean(distances))
    rms = float(np.sqrt(np.mean(distances * distances)))
    std = float(np.std(distances))
    if minimum <= 0.0:
        distortion = math.inf
    else:
        distortion = maximum / minimum

    scale_to_simplex_rms = SIMPLEX_DISTANCE / rms if rms > 0.0 else math.inf
    normalized = distances / rms if rms > 0.0 else np.full_like(distances, math.inf)

    return {
        "pair_distance_min": minimum,
        "pair_distance_mean": mean,
        "pair_distance_rms": rms,
        "pair_distance_std": std,
        "pair_distance_max": maximum,
        "pair_distance_distortion": distortion,
        "pair_abs_ratio_min": minimum / SIMPLEX_DISTANCE,
        "pair_abs_ratio_mean": mean / SIMPLEX_DISTANCE,
        "pair_abs_ratio_rms": rms / SIMPLEX_DISTANCE,
        "pair_abs_ratio_max": maximum / SIMPLEX_DISTANCE,
        "pair_rms_scale_to_simplex": scale_to_simplex_rms,
        "pair_normalized_ratio_min": float(np.min(normalized)),
        "pair_normalized_ratio_mean": float(np.mean(normalized)),
        "pair_normalized_ratio_max": float(np.max(normalized)),
    }


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
