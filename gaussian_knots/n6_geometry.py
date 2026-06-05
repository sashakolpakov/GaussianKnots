"""Shared N=6 geometry helpers for the Stiefel knot analysis."""

from __future__ import annotations

import itertools
import math
from statistics import NormalDist

import numpy as np


VERTEX_COUNT = 6
H_DIMENSION = VERTEX_COUNT - 1
KERNEL_DIMENSION = 2
EDGE_PAIRS = tuple(
    (first, second)
    for first in range(VERTEX_COUNT)
    for second in range(first + 1, VERTEX_COUNT)
    if (second - first) % VERTEX_COUNT not in (1, VERTEX_COUNT - 1)
)
QUADS = tuple(itertools.combinations(range(VERTEX_COUNT), 4))
QUAD_INDEX = {quad: index for index, quad in enumerate(QUADS)}
NORMAL = NormalDist()


def centered_simplex_vertices(vertex_count: int = VERTEX_COUNT) -> np.ndarray:
    identity = np.eye(vertex_count)
    return identity - np.ones((vertex_count, vertex_count)) / vertex_count


def h_basis(vertex_count: int = VERTEX_COUNT) -> np.ndarray:
    generators = np.eye(vertex_count)[:, : vertex_count - 1] - np.eye(vertex_count)[:, [vertex_count - 1]]
    q_matrix, _ = np.linalg.qr(generators)
    return q_matrix[:, : vertex_count - 1]


def simplex_vertices_in_h() -> np.ndarray:
    return centered_simplex_vertices(VERTEX_COUNT) @ h_basis(VERTEX_COUNT)


def determinant_wall_data(simplex_vertices: np.ndarray) -> list[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    data = []
    for first_edge, second_edge in EDGE_PAIRS:
        p0 = simplex_vertices[first_edge]
        p1 = simplex_vertices[(first_edge + 1) % len(simplex_vertices)]
        q0 = simplex_vertices[second_edge]
        q1 = simplex_vertices[(second_edge + 1) % len(simplex_vertices)]
        data.append((p1 - p0, q1 - q0, q0 - p0))
    return data


def determinant_wall_normals(wall_data: list[tuple[np.ndarray, np.ndarray, np.ndarray]]) -> np.ndarray:
    normals = []
    for edge_a, edge_b, displacement in wall_data:
        matrix = np.column_stack([edge_a, edge_b, displacement])
        coordinates = []
        for first in range(H_DIMENSION):
            for second in range(first + 1, H_DIMENSION):
                complement = [index for index in range(H_DIMENSION) if index not in (first, second)]
                sequence = complement + [first, second]
                inversions = sum(
                    1
                    for left in range(H_DIMENSION)
                    for right in range(left + 1, H_DIMENSION)
                    if sequence[left] > sequence[right]
                )
                coordinates.append(((-1) ** inversions) * np.linalg.det(matrix[complement, :]))
        normal = np.asarray(coordinates, dtype=float)
        normal /= np.linalg.norm(normal)
        pivot = int(np.argmax(np.abs(normal)))
        if normal[pivot] < 0:
            normal *= -1.0
        normals.append(normal)
    return np.vstack(normals)


def nonorthogonality_graph(projective_gram: np.ndarray) -> np.ndarray:
    graph = np.isclose(projective_gram, 0.25, atol=1e-10)
    np.fill_diagonal(graph, False)
    return graph


def projective_automorphism_lifts(gram: np.ndarray) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    projective_gram = np.round(np.abs(gram), 12)
    signed_gram = np.round(gram, 12)
    automorphisms = []
    for permutation in itertools.permutations(range(len(gram))):
        if not all(
            projective_gram[row, column] == projective_gram[permutation[row], permutation[column]]
            for row in range(len(gram))
            for column in range(len(gram))
        ):
            continue
        signs = sign_lift(signed_gram, permutation)
        if signs is not None:
            automorphisms.append((permutation, signs))
    return automorphisms


def sign_lift(gram: np.ndarray, permutation: tuple[int, ...]) -> tuple[int, ...] | None:
    signs: list[int | None] = [None] * len(permutation)
    signs[0] = 1
    stack = [0]
    while stack:
        row = stack.pop()
        for column in range(len(permutation)):
            if row == column or abs(gram[row, column]) < 1e-10:
                continue
            ratio = int(round(gram[row, column] / gram[permutation[row], permutation[column]]))
            value = signs[row] * ratio
            if signs[column] is None:
                signs[column] = value
                stack.append(column)
            elif signs[column] != value:
                return None
    if any(sign is None for sign in signs):
        return None
    resolved = tuple(int(sign) for sign in signs)
    if np.allclose(gram, np.outer(resolved, resolved) * gram[np.ix_(permutation, permutation)]):
        return resolved
    return None


def dihedral_wall_permutations() -> list[tuple[int, ...]]:
    index = {pair: position for position, pair in enumerate(EDGE_PAIRS)}
    permutations = []
    for shift in range(VERTEX_COUNT):
        permutations.append(wall_permutation({i: (i + shift) % VERTEX_COUNT for i in range(VERTEX_COUNT)}, index))
        permutations.append(wall_permutation({i: (shift - i) % VERTEX_COUNT for i in range(VERTEX_COUNT)}, index))
    return permutations


def wall_permutation(vertex_permutation: dict[int, int], index: dict[tuple[int, int], int]) -> tuple[int, ...]:
    values = []
    for first, second in EDGE_PAIRS:
        image = tuple(sorted((vertex_permutation[first], vertex_permutation[second])))
        values.append(index[image])
    return tuple(values)


def maximal_cliques(graph: np.ndarray, size: int) -> list[tuple[int, ...]]:
    cliques = []
    for subset in itertools.combinations(range(len(graph)), size):
        if all(graph[first, second] for first, second in itertools.combinations(subset, 2)):
            cliques.append(subset)
    return cliques


def common_clique_lines(
    wall_data: list[tuple[np.ndarray, np.ndarray, np.ndarray]],
    cliques: list[tuple[int, ...]],
) -> np.ndarray:
    projections = []
    for edge_a, edge_b, displacement in wall_data:
        normal = determinant_wall_normals([(edge_a, edge_b, displacement)])[0]
        skew = np.zeros((H_DIMENSION, H_DIMENSION))
        cursor = 0
        for first in range(H_DIMENSION):
            for second in range(first + 1, H_DIMENSION):
                skew[first, second] = normal[cursor]
                skew[second, first] = -normal[cursor]
                cursor += 1
        left, _, _ = np.linalg.svd(skew)
        basis = left[:, :2]
        projections.append(basis @ basis.T)

    lines = []
    for clique in cliques:
        projection_sum = sum(projections[index] for index in clique)
        _, eigenvectors = np.linalg.eigh(projection_sum)
        line = eigenvectors[:, -1]
        pivot = int(np.argmax(np.abs(line)))
        if line[pivot] < 0:
            line *= -1.0
        lines.append(line)
    return orient_regular_simplex_lines(np.vstack(lines))


def orient_regular_simplex_lines(lines: np.ndarray) -> np.ndarray:
    oriented = np.array(lines, copy=True)
    for index in range(1, len(oriented)):
        if float(oriented[0] @ oriented[index]) > 0:
            oriented[index] *= -1.0
    gram = oriented @ oriented.T
    expected = -0.2 * (np.ones_like(gram) - np.eye(len(oriented)))
    if not np.allclose(gram - np.eye(len(oriented)), expected):
        raise RuntimeError("common lines do not orient to a regular simplex")
    return oriented


def clique_action(permutation: tuple[int, ...], cliques: list[tuple[int, ...]]) -> tuple[int, ...]:
    clique_index = {tuple(sorted(clique)): index for index, clique in enumerate(cliques)}
    action = []
    for clique in cliques:
        image = tuple(sorted(permutation[index] for index in clique))
        action.append(clique_index[image])
    return tuple(action)


def orthogonal_maps(
    common_lines: np.ndarray,
    cliques: list[tuple[int, ...]],
    full_lifts: list[tuple[tuple[int, ...], tuple[int, ...]]],
) -> list[np.ndarray]:
    source = common_lines.T
    source_pinv = np.linalg.pinv(source)
    maps = []
    for permutation, _ in full_lifts:
        action = clique_action(permutation, cliques)
        target = common_lines[list(action)].T
        maps.append(target @ source_pinv)
    return maps


def orthogonal_lift_errors(
    normals: np.ndarray,
    common_lines: np.ndarray,
    cliques: list[tuple[int, ...]],
    automorphisms: list[tuple[tuple[int, ...], tuple[int, ...]]],
) -> list[tuple[float, float]]:
    maps = orthogonal_maps(common_lines, cliques, automorphisms)
    errors = []
    for linear_map, (wall_permutation, _) in zip(maps, automorphisms):
        orthogonal_error = float(np.linalg.norm(linear_map.T @ linear_map - np.eye(linear_map.shape[0])))
        wedge = wedge2_matrix(linear_map)
        wall_error = 0.0
        for source_index, target_index in enumerate(wall_permutation):
            image = wedge @ normals[source_index]
            image /= np.linalg.norm(image)
            target_normal = normals[target_index]
            wall_error = max(wall_error, min(np.linalg.norm(image - target_normal), np.linalg.norm(image + target_normal)))
        errors.append((orthogonal_error, wall_error))
    return errors


def wedge2_matrix(linear_map: np.ndarray) -> np.ndarray:
    basis_pairs = [(first, second) for first in range(H_DIMENSION) for second in range(first + 1, H_DIMENSION)]
    matrix = np.zeros((len(basis_pairs), len(basis_pairs)))
    for column, (first, second) in enumerate(basis_pairs):
        image = np.outer(linear_map[:, first], linear_map[:, second])
        skew = image - image.T
        for row, (target_first, target_second) in enumerate(basis_pairs):
            matrix[row, column] = skew[target_first, target_second]
    return matrix


def random_kernel(rng: np.random.Generator) -> np.ndarray:
    matrix = rng.normal(size=(H_DIMENSION, KERNEL_DIMENSION))
    q_matrix, r_matrix = np.linalg.qr(matrix, mode="reduced")
    for column in range(KERNEL_DIMENSION):
        if r_matrix[column, column] < 0:
            q_matrix[:, column] *= -1.0
    return q_matrix


def vertices_from_kernel(simplex_vertices: np.ndarray, kernel_basis: np.ndarray) -> np.ndarray:
    _, _, vh_matrix = np.linalg.svd(kernel_basis.T, full_matrices=True)
    quotient_basis = vh_matrix[KERNEL_DIMENSION:].T
    return simplex_vertices @ quotient_basis


def kernel_from_vertices(simplex_vertices: np.ndarray, vertices: np.ndarray) -> np.ndarray:
    quotient_basis = np.linalg.lstsq(simplex_vertices, vertices, rcond=None)[0]
    _, _, vh_matrix = np.linalg.svd(quotient_basis.T, full_matrices=True)
    return vh_matrix[3:].T


def orientation_signature(vertices: np.ndarray) -> str:
    signs = []
    for first, second, third, fourth in QUADS:
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
        signs.append("+" if value > 0.0 else "-")
    return canonical_global("".join(signs))


def d6_orbit_representative(signature: str) -> str:
    return min(d6_orbit(signature))


def d6_orbit(signature: str) -> set[str]:
    return {apply_vertex_permutation(signature, permutation) for permutation in d6_vertex_permutations()}


def d6_vertex_permutations() -> list[dict[int, int]]:
    permutations = []
    for shift in range(VERTEX_COUNT):
        permutations.append({index: (index + shift) % VERTEX_COUNT for index in range(VERTEX_COUNT)})
        permutations.append({index: (shift - index) % VERTEX_COUNT for index in range(VERTEX_COUNT)})
    return permutations


def apply_vertex_permutation(signature: str, permutation: dict[int, int]) -> str:
    transformed = [""] * len(QUADS)
    for index, quad in enumerate(QUADS):
        image = [permutation[vertex] for vertex in quad]
        sorted_image = tuple(sorted(image))
        char = signature[index]
        if permutation_parity(image) < 0:
            char = flip_sign(char)
        transformed[QUAD_INDEX[sorted_image]] = char
    return canonical_global("".join(transformed))


def permutation_parity(values: list[int]) -> int:
    inversions = sum(
        1 for first in range(len(values)) for second in range(first + 1, len(values)) if values[first] > values[second]
    )
    return -1 if inversions % 2 else 1


def canonical_global(signature: str) -> str:
    return min(signature, "".join(flip_sign(char) for char in signature))


def flip_sign(char: str) -> str:
    return "-" if char == "+" else "+"


def wilson_interval(successes: int, trials: int, confidence: float = 0.95) -> tuple[float, float]:
    if trials <= 0:
        return math.nan, math.nan
    z = NORMAL.inv_cdf(0.5 + confidence / 2.0)
    phat = successes / trials
    denominator = 1.0 + z * z / trials
    center = (phat + z * z / (2.0 * trials)) / denominator
    radius = z * math.sqrt((phat * (1.0 - phat) + z * z / (4.0 * trials)) / trials) / denominator
    return max(0.0, center - radius), min(1.0, center + radius)
