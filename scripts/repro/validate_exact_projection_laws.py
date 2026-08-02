#!/usr/bin/env python3
"""Reproduce two exact distributional laws for projected-simplex polygons.

The script checks:

1. For a fixed pair of nonincident edges in a planar Gaussian polygon,

       P(the edges cross) = (2 / pi) asin(1 / 3).

   Consequently the expected crossing count of the N-cycle is

       N(N - 3) / 2 * (2 / pi) asin(1 / 3).

   The Gaussian polygon and its centered Haar row factor are deliberately
   coupled below.  Centering and the invertible polar/QR factor do not change
   segment intersections, so their crossing counts must agree samplewise.

2. If Q is Haar on St_3(1^perp) and i != j, then

       ||Q(e_i - e_j)||^2 / 2 ~ Beta(3/2, (N - 4)/2),   N > 4.

   For N=4 the normalized squared distance is identically one.  The law is a
   one-pair marginal statement; the different pair distances are not
   independent (their squared sum is deterministically 3N).

SciPy is optional.  When it is available, a one-sample Kolmogorov--Smirnov
check is included for the beta law.
"""

from __future__ import annotations

import argparse
import math

import numpy as np


CROSSING_PROBABILITY = 2.0 / math.pi * math.asin(1.0 / 3.0)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crossing-n", default="6,10,20,40", help="comma-separated planar vertex counts")
    parser.add_argument("--crossing-samples", type=int, default=10000)
    parser.add_argument("--beta-n", default="5,6,8,12,20,50", help="comma-separated Haar vertex counts")
    parser.add_argument("--beta-samples", type=int, default=50000)
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20260802)
    return parser


def parse_counts(text: str, minimum: int) -> list[int]:
    counts = [int(part.strip()) for part in text.split(",") if part.strip()]
    if not counts or any(count < minimum for count in counts):
        raise SystemExit(f"all vertex counts must be at least {minimum}")
    return counts


def nonincident_edge_pairs(vertex_count: int) -> tuple[np.ndarray, np.ndarray]:
    first: list[int] = []
    second: list[int] = []
    for edge_a in range(vertex_count):
        for edge_b in range(edge_a + 1, vertex_count):
            if edge_b == edge_a + 1 or (edge_a == 0 and edge_b == vertex_count - 1):
                continue
            first.append(edge_a)
            second.append(edge_b)
    return np.asarray(first, dtype=int), np.asarray(second, dtype=int)


def cross2(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return left[..., 0] * right[..., 1] - left[..., 1] * right[..., 0]


def planar_crossing_counts(
    vertices: np.ndarray,
    first_edges: np.ndarray,
    second_edges: np.ndarray,
) -> np.ndarray:
    """Count proper crossings in a batch with shape (batch, N, 2)."""

    points = np.asarray(vertices, dtype=float)
    vertex_count = points.shape[1]
    a0 = points[:, first_edges]
    a1 = points[:, (first_edges + 1) % vertex_count]
    b0 = points[:, second_edges]
    b1 = points[:, (second_edges + 1) % vertex_count]

    orient_a_b0 = cross2(a1 - a0, b0 - a0)
    orient_a_b1 = cross2(a1 - a0, b1 - a0)
    orient_b_a0 = cross2(b1 - b0, a0 - b0)
    orient_b_a1 = cross2(b1 - b0, a1 - b0)
    separated_by_a = ((orient_a_b0 > 0.0) & (orient_a_b1 < 0.0)) | (
        (orient_a_b0 < 0.0) & (orient_a_b1 > 0.0)
    )
    separated_by_b = ((orient_b_a0 > 0.0) & (orient_b_a1 < 0.0)) | (
        (orient_b_a0 < 0.0) & (orient_b_a1 > 0.0)
    )
    return np.sum(separated_by_a & separated_by_b, axis=1)


def haar_frame_from_centered_gaussian(points: np.ndarray) -> np.ndarray:
    centered = points - np.mean(points, axis=1, keepdims=True)
    frame, _ = np.linalg.qr(centered, mode="reduced")
    return frame


def crossing_experiment(
    vertex_counts: list[int],
    samples: int,
    batch_size: int,
    rng: np.random.Generator,
) -> None:
    print("PLANAR CROSSINGS")
    print(f"fixed-pair exact probability = {CROSSING_PROBABILITY:.12f}")
    print(
        "McDiarmid: P(|C-E C| >= t) <= "
        "2 exp(-t^2 / (2 N (N-3)^2)); changing one vertex affects at most 2(N-3) indicators"
    )
    print("equivalently, for M=N(N-3)/2: P(|C-E C| >= epsilon M) <= 2 exp(-epsilon^2 N/8)")
    print("N,pairs,samples,exact_mean,gaussian_mean,haar_mean,empirical_sd,z_score,max_coupled_difference")
    for vertex_count in vertex_counts:
        edge_a, edge_b = nonincident_edge_pairs(vertex_count)
        values: list[np.ndarray] = []
        maximum_difference = 0
        for start in range(0, samples, batch_size):
            current = min(batch_size, samples - start)
            gaussian = rng.normal(size=(current, vertex_count, 2))
            haar = haar_frame_from_centered_gaussian(gaussian)
            gaussian_counts = planar_crossing_counts(gaussian, edge_a, edge_b)
            haar_counts = planar_crossing_counts(haar, edge_a, edge_b)
            maximum_difference = max(
                maximum_difference,
                int(np.max(np.abs(gaussian_counts - haar_counts))),
            )
            values.append(gaussian_counts)
        counts = np.concatenate(values).astype(float)
        exact_mean = len(edge_a) * CROSSING_PROBABILITY
        empirical_sd = float(np.std(counts))
        standard_error = empirical_sd / math.sqrt(samples)
        z_score = (float(np.mean(counts)) - exact_mean) / standard_error if standard_error else math.nan
        print(
            f"{vertex_count},{len(edge_a)},{samples},{exact_mean:.8f},"
            f"{np.mean(counts):.8f},{np.mean(counts):.8f},{empirical_sd:.8f},"
            f"{z_score:.4f},{maximum_difference}"
        )
        if maximum_difference != 0:
            raise RuntimeError("the coupled Gaussian and Haar crossing counts differ")


def beta_parameters(vertex_count: int) -> tuple[float, float]:
    return 1.5, (vertex_count - 4.0) / 2.0


def beta_exact_moments(vertex_count: int) -> tuple[float, float]:
    """Mean and variance of ||Q(e_i-e_j)||^2 / 2."""

    alpha, beta = beta_parameters(vertex_count)
    total = alpha + beta
    mean = alpha / total
    variance = alpha * beta / (total * total * (total + 1.0))
    return mean, variance


def beta_ks(values: np.ndarray, alpha: float, beta: float) -> tuple[float, float] | None:
    try:
        from scipy.stats import kstest
    except ImportError:
        return None
    statistic, p_value = kstest(values, "beta", args=(alpha, beta))
    return float(statistic), float(p_value)


def beta_experiment(
    vertex_counts: list[int],
    samples: int,
    batch_size: int,
    rng: np.random.Generator,
) -> None:
    print()
    print("HAAR SQUARED PAIR DISTANCES")
    print("law: ||Q(e_i-e_j)||^2 / 2 ~ Beta(3/2,(N-4)/2), for N>4")
    print("N,samples,empirical_mean,exact_mean,empirical_variance,exact_variance,KS_D,KS_p")
    for vertex_count in vertex_counts:
        if vertex_count == 4:
            print(f"4,{samples},1.00000000,1.00000000,0.00000000,0.00000000,degenerate,degenerate")
            continue
        values: list[np.ndarray] = []
        for start in range(0, samples, batch_size):
            current = min(batch_size, samples - start)
            gaussian = rng.normal(size=(current, vertex_count, 3))
            haar = haar_frame_from_centered_gaussian(gaussian)
            squared_distance = np.sum((haar[:, 0] - haar[:, 1]) ** 2, axis=1)
            values.append(squared_distance / 2.0)
        normalized = np.concatenate(values)
        exact_mean, exact_variance = beta_exact_moments(vertex_count)
        alpha, beta = beta_parameters(vertex_count)
        ks = beta_ks(normalized, alpha, beta)
        ks_statistic, ks_p = (math.nan, math.nan) if ks is None else ks
        print(
            f"{vertex_count},{samples},{np.mean(normalized):.8f},{exact_mean:.8f},"
            f"{np.var(normalized):.8f},{exact_variance:.8f},"
            f"{ks_statistic:.8f},{ks_p:.8f}"
        )


def main() -> int:
    args = build_parser().parse_args()
    if args.crossing_samples < 1 or args.beta_samples < 1 or args.batch_size < 1:
        raise SystemExit("sample counts and --batch-size must be positive")
    crossing_counts = parse_counts(args.crossing_n, minimum=4)
    beta_counts = parse_counts(args.beta_n, minimum=4)
    crossing_rng = np.random.default_rng(args.seed)
    beta_rng = np.random.default_rng(args.seed + 1)
    crossing_experiment(crossing_counts, args.crossing_samples, args.batch_size, crossing_rng)
    beta_experiment(beta_counts, args.beta_samples, args.batch_size, beta_rng)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
