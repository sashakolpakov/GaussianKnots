"""General geometry helpers for Stiefel wall and order-type probes."""

from __future__ import annotations

import itertools
import math
from collections import Counter
from dataclasses import dataclass
from statistics import NormalDist

import numpy as np


TARGET_DIMENSION = 3
NORMAL = NormalDist()


@dataclass(frozen=True)
class ChamberModel:
    vertex_count: int
    h_dimension: int
    kernel_dimension: int
    simplex_vertices: np.ndarray
    edge_pairs: tuple[tuple[int, int], ...]
    quads: tuple[tuple[int, int, int, int], ...]
    quad_index: dict[tuple[int, int, int, int], int]
    wall_data: list[tuple[np.ndarray, np.ndarray, np.ndarray]]
    wall_normals: np.ndarray
    wall_gram: np.ndarray


def build_chamber_model(vertex_count: int) -> ChamberModel:
    if vertex_count < 5:
        raise ValueError("the projected-simplex knot model expects N >= 5")
    h_dimension = vertex_count - 1
    kernel_dimension = h_dimension - TARGET_DIMENSION
    if kernel_dimension < 1:
        raise ValueError("target dimension 3 requires dim(1^perp) > 3")
    simplex_vertices = centered_simplex_vertices(vertex_count) @ h_basis(vertex_count)
    pairs = edge_pairs(vertex_count)
    quad_list = tuple(itertools.combinations(range(vertex_count), 4))
    wall_data = determinant_wall_data(simplex_vertices, pairs)
    wall_normals = determinant_wall_normals(wall_data, h_dimension, kernel_dimension)
    return ChamberModel(
        vertex_count=vertex_count,
        h_dimension=h_dimension,
        kernel_dimension=kernel_dimension,
        simplex_vertices=simplex_vertices,
        edge_pairs=pairs,
        quads=quad_list,
        quad_index={quad: index for index, quad in enumerate(quad_list)},
        wall_data=wall_data,
        wall_normals=wall_normals,
        wall_gram=wall_normals @ wall_normals.T,
    )


def centered_simplex_vertices(vertex_count: int) -> np.ndarray:
    identity = np.eye(vertex_count)
    return identity - np.ones((vertex_count, vertex_count)) / vertex_count


def h_basis(vertex_count: int) -> np.ndarray:
    generators = np.eye(vertex_count)[:, : vertex_count - 1] - np.eye(vertex_count)[:, [vertex_count - 1]]
    q_matrix, _ = np.linalg.qr(generators)
    return q_matrix[:, : vertex_count - 1]


def edge_pairs(vertex_count: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        (first, second)
        for first in range(vertex_count)
        for second in range(first + 1, vertex_count)
        if (second - first) % vertex_count not in (1, vertex_count - 1)
    )


def determinant_wall_data(
    simplex_vertices: np.ndarray,
    pairs: tuple[tuple[int, int], ...],
) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    data = []
    vertex_count = len(simplex_vertices)
    for first_edge, second_edge in pairs:
        p0 = simplex_vertices[first_edge]
        p1 = simplex_vertices[(first_edge + 1) % vertex_count]
        q0 = simplex_vertices[second_edge]
        q1 = simplex_vertices[(second_edge + 1) % vertex_count]
        data.append((p1 - p0, q1 - q0, q0 - p0))
    return data


def determinant_wall_normals(
    wall_data: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
    h_dimension: int,
    kernel_dimension: int,
) -> np.ndarray:
    coordinates = tuple(itertools.combinations(range(h_dimension), kernel_dimension))
    basis = np.eye(h_dimension)
    normals = []
    for edge_a, edge_b, displacement in wall_data:
        normal = []
        for coordinate in coordinates:
            matrix = np.column_stack([edge_a, edge_b, basis[:, coordinate], displacement])
            normal.append(float(np.linalg.det(matrix)))
        normal_array = np.asarray(normal, dtype=float)
        norm = float(np.linalg.norm(normal_array))
        if norm <= 0.0:
            raise RuntimeError("zero determinant-wall normal")
        normal_array /= norm
        pivot = int(np.argmax(np.abs(normal_array)))
        if normal_array[pivot] < 0.0:
            normal_array *= -1.0
        normals.append(normal_array)
    return np.vstack(normals)


def determinant_wall_signature(
    wall_data: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
    kernel_basis: np.ndarray,
    tolerance: float = 1e-10,
) -> str:
    signs = []
    for edge_a, edge_b, displacement in wall_data:
        value = float(np.linalg.det(np.column_stack([edge_a, edge_b, kernel_basis, displacement])))
        signs.append(sign_char(value, tolerance=tolerance))
    return "".join(signs)


def random_kernel(model: ChamberModel, rng: np.random.Generator) -> np.ndarray:
    matrix = rng.normal(size=(model.h_dimension, model.kernel_dimension))
    q_matrix, r_matrix = np.linalg.qr(matrix, mode="reduced")
    for column in range(model.kernel_dimension):
        if r_matrix[column, column] < 0.0:
            q_matrix[:, column] *= -1.0
    return q_matrix


def vertices_from_kernel(model: ChamberModel, kernel_basis: np.ndarray) -> np.ndarray:
    _, _, vh_matrix = np.linalg.svd(kernel_basis.T, full_matrices=True)
    quotient_basis = vh_matrix[model.kernel_dimension :].T
    return model.simplex_vertices @ quotient_basis


def orientation_signature(
    vertices: np.ndarray,
    quads: tuple[tuple[int, int, int, int], ...],
    tolerance: float = 1e-10,
) -> str:
    signs = []
    for first, second, third, fourth in quads:
        value = float(
            np.linalg.det(
                np.column_stack(
                    [
                        vertices[second] - vertices[first],
                        vertices[third] - vertices[first],
                        vertices[fourth] - vertices[first],
                    ]
                )
            )
        )
        signs.append(sign_char(value, tolerance=tolerance))
    return canonical_global("".join(signs))


def d_orbit_representative(
    signature: str,
    vertex_count: int,
    quads: tuple[tuple[int, int, int, int], ...],
    quad_index: dict[tuple[int, int, int, int], int],
) -> str:
    return min(d_orbit(signature, vertex_count, quads, quad_index))


def d_orbit(
    signature: str,
    vertex_count: int,
    quads: tuple[tuple[int, int, int, int], ...],
    quad_index: dict[tuple[int, int, int, int], int],
) -> set[str]:
    return {
        apply_vertex_permutation(signature, permutation, quads, quad_index)
        for permutation in dihedral_vertex_permutations(vertex_count)
    }


def dihedral_vertex_permutations(vertex_count: int) -> list[dict[int, int]]:
    permutations = []
    for shift in range(vertex_count):
        permutations.append({index: (index + shift) % vertex_count for index in range(vertex_count)})
        permutations.append({index: (shift - index) % vertex_count for index in range(vertex_count)})
    return permutations


def apply_vertex_permutation(
    signature: str,
    permutation: dict[int, int],
    quads: tuple[tuple[int, int, int, int], ...],
    quad_index: dict[tuple[int, int, int, int], int],
) -> str:
    transformed = [""] * len(quads)
    for index, quad in enumerate(quads):
        image = [permutation[vertex] for vertex in quad]
        sorted_image = tuple(sorted(image))
        char = signature[index]
        if permutation_parity(image) < 0:
            char = flip_sign(char)
        transformed[quad_index[sorted_image]] = char
    return canonical_global("".join(transformed))


def dihedral_wall_permutations(model: ChamberModel) -> list[tuple[int, ...]]:
    return [wall_permutation(model, permutation) for permutation in dihedral_vertex_permutations(model.vertex_count)]


def wall_permutation(model: ChamberModel, vertex_permutation: dict[int, int]) -> tuple[int, ...]:
    pair_index = {pair: position for position, pair in enumerate(model.edge_pairs)}
    image = []
    for first_edge, second_edge in model.edge_pairs:
        first_image = cycle_edge_image(first_edge, vertex_permutation, model.vertex_count)
        second_image = cycle_edge_image(second_edge, vertex_permutation, model.vertex_count)
        image.append(pair_index[tuple(sorted((first_image, second_image)))])
    return tuple(image)


def cycle_edge_image(edge_start: int, vertex_permutation: dict[int, int], vertex_count: int) -> int:
    left = vertex_permutation[edge_start]
    right = vertex_permutation[(edge_start + 1) % vertex_count]
    if (right - left) % vertex_count == 1:
        return left
    if (left - right) % vertex_count == 1:
        return right
    raise ValueError("vertex permutation does not preserve the Hamiltonian cycle")


def wall_gram_automorphisms(
    gram: np.ndarray,
    decimals: int = 12,
    max_count: int | None = None,
) -> list[tuple[int, ...]]:
    projective_gram = np.round(np.abs(gram), decimals)
    vertex_count = len(projective_gram)
    row_colors = [row_color(projective_gram[index]) for index in range(vertex_count)]
    color_buckets: dict[tuple[tuple[float, int], ...], list[int]] = {}
    for index, color in enumerate(row_colors):
        color_buckets.setdefault(color, []).append(index)

    mapping = [-1] * vertex_count
    used: set[int] = set()
    automorphisms: list[tuple[int, ...]] = []

    def candidates_for(source: int) -> list[int]:
        candidates = []
        for target in color_buckets[row_colors[source]]:
            if target in used:
                continue
            if all(
                projective_gram[source, mapped_source] == projective_gram[target, mapped_target]
                for mapped_source, mapped_target in enumerate(mapping)
                if mapped_target != -1
            ):
                candidates.append(target)
        return candidates

    def choose_source() -> tuple[int, list[int]]:
        best_source = -1
        best_candidates: list[int] | None = None
        for source in range(vertex_count):
            if mapping[source] != -1:
                continue
            candidates = candidates_for(source)
            if best_candidates is None or len(candidates) < len(best_candidates):
                best_source = source
                best_candidates = candidates
                if len(candidates) <= 1:
                    break
        return best_source, best_candidates or []

    def search() -> None:
        if max_count is not None and len(automorphisms) >= max_count:
            return
        if len(used) == vertex_count:
            automorphisms.append(tuple(mapping))
            return
        source, candidates = choose_source()
        for target in candidates:
            mapping[source] = target
            used.add(target)
            search()
            used.remove(target)
            mapping[source] = -1

    search()
    return automorphisms


def projective_automorphisms(
    gram: np.ndarray,
    decimals: int = 12,
    max_count: int | None = None,
) -> list[tuple[int, ...]]:
    """Backward-compatible alias for the wall-normal Gram automorphism group.

    This is not the full group of arbitrary projective-linear automorphisms of
    the normal configuration.  It is the group relevant to orthogonal wall
    symmetries and Haar-volume comparisons.
    """

    return wall_gram_automorphisms(gram, decimals=decimals, max_count=max_count)


def row_color(row: np.ndarray) -> tuple[tuple[float, int], ...]:
    return tuple(sorted(Counter(float(value) for value in row).items()))


def is_projective_automorphism(gram: np.ndarray, permutation: tuple[int, ...], decimals: int = 12) -> bool:
    projective_gram = np.round(np.abs(gram), decimals)
    return bool(np.all(projective_gram == projective_gram[np.ix_(permutation, permutation)]))


def permutation_parity(values: list[int]) -> int:
    inversions = sum(
        1 for first in range(len(values)) for second in range(first + 1, len(values)) if values[first] > values[second]
    )
    return -1 if inversions % 2 else 1


def canonical_global(signature: str) -> str:
    flipped = "".join(flip_sign(char) for char in signature)
    return min(signature, flipped)


def sign_char(value: float, tolerance: float = 1e-10) -> str:
    if abs(value) <= tolerance:
        return "0"
    return "+" if value > 0.0 else "-"


def flip_sign(char: str) -> str:
    if char == "+":
        return "-"
    if char == "-":
        return "+"
    return char


def wilson_interval(successes: int, trials: int, confidence: float = 0.95) -> tuple[float, float]:
    if trials <= 0:
        return math.nan, math.nan
    z = NORMAL.inv_cdf(0.5 + confidence / 2.0)
    phat = successes / trials
    denominator = 1.0 + z * z / trials
    center = (phat + z * z / (2.0 * trials)) / denominator
    radius = z * math.sqrt((phat * (1.0 - phat) + z * z / (4.0 * trials)) / trials) / denominator
    return max(0.0, center - radius), min(1.0, center + radius)
