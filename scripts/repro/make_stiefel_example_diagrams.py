#!/usr/bin/env python3
"""Generate a reproducible unknot/trefoil diagram for the Stiefel manuscript."""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
GAUSSIAN_KNOTS = ROOT
sys.path.insert(0, str(ROOT))

from gaussian_knots.generation import projected_simplex_polygon


@dataclass(frozen=True)
class Example:
    title: str
    vertex_count: int
    sample_seed: int
    knot_label: str
    desired_crossings: int
    color: str


@dataclass(frozen=True)
class Crossing:
    edge_a: int
    edge_b: int
    t_a: float
    t_b: float
    point: np.ndarray
    edge_a_over: bool


EXAMPLES = (
    Example(
        title="Haar sample: unknot",
        vertex_count=5,
        sample_seed=896654480,
        knot_label="unknot",
        desired_crossings=0,
        color="blue!65!black",
    ),
    Example(
        title="Haar sample: trefoil",
        vertex_count=6,
        sample_seed=1062611651,
        knot_label=r"$3_1$",
        desired_crossings=3,
        color="red!70!black",
    ),
)


def main() -> int:
    out_dir = ROOT / "figures"
    out_dir.mkdir(exist_ok=True)
    panels = []
    for index, example in enumerate(EXAMPLES):
        vertices = sample_vertices(example)
        projection = choose_projection(vertices, example.desired_crossings, seed=example.sample_seed + 1729)
        panels.append(render_panel(example, vertices, projection, x_offset=4.8 * index))

    tex = standalone_document("\n".join(panels))
    output = out_dir / "stiefel_unknot_trefoil.tex"
    output.write_text(tex, encoding="utf-8")
    print(f"wrote {output}")
    return 0


def sample_vertices(example: Example) -> np.ndarray:
    rng = np.random.default_rng(example.sample_seed)
    return projected_simplex_polygon(example.vertex_count, rng, projection_model="haar")


def choose_projection(vertices: np.ndarray, desired_crossings: int, seed: int) -> np.ndarray:
    """Choose an orthonormal camera basis with a clear crossing diagram."""

    rng = np.random.default_rng(seed)
    best_basis = None
    best_score = -math.inf
    for _ in range(6000):
        basis = random_camera_basis(rng)
        projected, depth = project(vertices, basis)
        crossings = diagram_crossings(projected, depth)
        if len(crossings) != desired_crossings:
            continue
        score = projection_score(projected, crossings)
        if score > best_score:
            best_score = score
            best_basis = basis
    if best_basis is None:
        # Fall back to a deterministic oblique view if the requested crossing
        # count is not found.  The selected seeds currently find the target.
        best_basis = camera_basis_from_view(np.array([0.37, -0.58, 0.72]))
    return best_basis


def random_camera_basis(rng: np.random.Generator) -> np.ndarray:
    view = rng.normal(size=3)
    view /= np.linalg.norm(view)
    return camera_basis_from_view(view)


def camera_basis_from_view(view: np.ndarray) -> np.ndarray:
    w = np.asarray(view, dtype=float)
    w /= np.linalg.norm(w)
    helper = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(helper, w))) > 0.9:
        helper = np.array([0.0, 1.0, 0.0])
    u = np.cross(helper, w)
    u /= np.linalg.norm(u)
    v = np.cross(w, u)
    return np.vstack([u, v, w])


def project(vertices: np.ndarray, basis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    transformed = vertices @ basis.T
    return transformed[:, :2], transformed[:, 2]


def diagram_crossings(projected: np.ndarray, depth: np.ndarray) -> list[Crossing]:
    n = len(projected)
    crossings: list[Crossing] = []
    for a in range(n):
        a_next = (a + 1) % n
        for b in range(a + 1, n):
            if edges_adjacent(a, b, n):
                continue
            b_next = (b + 1) % n
            hit = segment_intersection_parameters(
                projected[a], projected[a_next], projected[b], projected[b_next]
            )
            if hit is None:
                continue
            t_a, t_b, point = hit
            z_a = (1.0 - t_a) * depth[a] + t_a * depth[a_next]
            z_b = (1.0 - t_b) * depth[b] + t_b * depth[b_next]
            if abs(z_a - z_b) < 1e-5:
                continue
            crossings.append(Crossing(a, b, t_a, t_b, point, edge_a_over=z_a > z_b))
    return crossings


def edges_adjacent(first: int, second: int, n: int) -> bool:
    return first == second or (first + 1) % n == second or (second + 1) % n == first


def segment_intersection_parameters(
    p: np.ndarray,
    p_next: np.ndarray,
    q: np.ndarray,
    q_next: np.ndarray,
) -> tuple[float, float, np.ndarray] | None:
    r = p_next - p
    s = q_next - q
    denom = cross2(r, s)
    if abs(denom) < 1e-10:
        return None
    qp = q - p
    t = cross2(qp, s) / denom
    u = cross2(qp, r) / denom
    if 1e-4 < t < 1.0 - 1e-4 and 1e-4 < u < 1.0 - 1e-4:
        return float(t), float(u), p + t * r
    return None


def cross2(a: np.ndarray, b: np.ndarray) -> float:
    return float(a[0] * b[1] - a[1] * b[0])


def projection_score(projected: np.ndarray, crossings: list[Crossing]) -> float:
    # Prefer views with well-spread vertices and separated crossings.
    normalized = normalize(projected)
    dists = []
    for i in range(len(normalized)):
        for j in range(i + 1, len(normalized)):
            dists.append(float(np.linalg.norm(normalized[i] - normalized[j])))
    min_vertex_distance = min(dists) if dists else 0.0
    if len(crossings) < 2:
        min_crossing_distance = 1.0
    else:
        min_crossing_distance = min(
            float(np.linalg.norm(normalized_crossing(a, projected, normalized) - normalized_crossing(b, projected, normalized)))
            for index, a in enumerate(crossings)
            for b in crossings[index + 1 :]
        )
    if crossings:
        min_endpoint_parameter = min(
            min(c.t_a, 1.0 - c.t_a, c.t_b, 1.0 - c.t_b) for c in crossings
        )
        min_angle = min(crossing_sine(projected, c) for c in crossings)
        min_crossing_vertex_distance = min(
            float(np.linalg.norm(normalized_crossing(c, projected, normalized) - point))
            for c in crossings
            for point in normalized
        )
    else:
        min_endpoint_parameter = 1.0
        min_angle = 1.0
        min_crossing_vertex_distance = 1.0
    return (
        min_vertex_distance
        + min_crossing_distance
        + min_crossing_vertex_distance
        + 1.25 * min_angle
        + 1.25 * min_endpoint_parameter
    )


def normalized_crossing(crossing: Crossing, original: np.ndarray, normalized: np.ndarray) -> np.ndarray:
    start = normalized[crossing.edge_a]
    end = normalized[(crossing.edge_a + 1) % len(normalized)]
    return (1.0 - crossing.t_a) * start + crossing.t_a * end


def crossing_sine(projected: np.ndarray, crossing: Crossing) -> float:
    n = len(projected)
    a = projected[(crossing.edge_a + 1) % n] - projected[crossing.edge_a]
    b = projected[(crossing.edge_b + 1) % n] - projected[crossing.edge_b]
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0.0:
        return 0.0
    return abs(cross2(a, b)) / denom


def render_panel(example: Example, vertices: np.ndarray, basis: np.ndarray, x_offset: float) -> str:
    projected, depth = project(vertices, basis)
    crossings = diagram_crossings(projected, depth)
    normalized = normalize(projected)
    crossing_data = rescale_crossings(crossings, projected, normalized)
    under_params = {edge: [] for edge in range(len(vertices))}
    for crossing in crossing_data:
        if crossing.edge_a_over:
            under_params[crossing.edge_b].append(crossing.t_b)
        else:
            under_params[crossing.edge_a].append(crossing.t_a)

    lines = [rf"\begin{{scope}}[xshift={x_offset:.3f}cm]"]
    lines.append(rf"\node[font=\small] at (0,1.85) {{{example.title}}};")
    lines.append(
        rf"\node[font=\scriptsize] at (0,1.58) "
        rf"{{$N={example.vertex_count}$, {example.knot_label}, seed {example.sample_seed}}};"
    )
    lines.append(rf"\draw[gray!35] (-1.75,-1.45) rectangle (1.75,1.45);")

    gap = 0.055
    for edge in range(len(vertices)):
        start = normalized[edge]
        end = normalized[(edge + 1) % len(vertices)]
        intervals = visible_intervals(sorted(under_params[edge]), gap)
        for left, right in intervals:
            a = (1.0 - left) * start + left * end
            b = (1.0 - right) * start + right * end
            lines.append(
                rf"\draw[very thick,{example.color}] "
                rf"({a[0]:.4f},{a[1]:.4f}) -- ({b[0]:.4f},{b[1]:.4f});"
            )
    for point in normalized:
        lines.append(rf"\fill[{example.color}] ({point[0]:.4f},{point[1]:.4f}) circle (0.035);")
    lines.append(r"\end{scope}")
    return "\n".join(lines)


def normalize(points: np.ndarray) -> np.ndarray:
    centered = points - np.mean(points, axis=0)
    scale = float(np.max(np.linalg.norm(centered, axis=1)))
    if scale <= 0:
        scale = 1.0
    return 1.25 * centered / scale


def rescale_crossings(
    crossings: list[Crossing],
    original: np.ndarray,
    normalized: np.ndarray,
) -> list[Crossing]:
    # Only the edge indices and parameters are used after normalization.
    return [
        Crossing(
            edge_a=c.edge_a,
            edge_b=c.edge_b,
            t_a=c.t_a,
            t_b=c.t_b,
            point=np.zeros(2),
            edge_a_over=c.edge_a_over,
        )
        for c in crossings
    ]


def visible_intervals(under_crossings: list[float], gap: float) -> list[tuple[float, float]]:
    intervals = []
    cursor = 0.0
    for t in under_crossings:
        left = max(0.0, t - gap)
        right = min(1.0, t + gap)
        if left > cursor:
            intervals.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1.0:
        intervals.append((cursor, 1.0))
    return intervals


def standalone_document(body: str) -> str:
    return rf"""\documentclass[tikz,border=4pt]{{standalone}}
\usepackage{{amsmath}}
\begin{{document}}
\begin{{tikzpicture}}[line cap=round,line join=round]
{body}
\end{{tikzpicture}}
\end{{document}}
"""


if __name__ == "__main__":
    raise SystemExit(main())
