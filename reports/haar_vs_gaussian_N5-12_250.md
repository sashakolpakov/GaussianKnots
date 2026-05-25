# Haar vs Gaussian Projection Experiment

Date: 2026-05-24

This is a first-pass Monte Carlo experiment for the Hamiltonian cycle
`1 -> 2 -> ... -> N -> 1` through projected simplex vertices.

The Haar run samples row-orthonormal projections, i.e. Haar points in the
Grassmann/Stiefel model governing knot type in the manuscript.  The Gaussian
run samples raw Gaussian vertices from a `3 x N` Gaussian projection matrix.
The manuscript predicts the same knot-type law for these two models; the
finite-sample tables below are consistent at this sample size but are not a
precision estimate.

Environment:

- `pyknotid` version: `0.5.4`
- fast backend: `True`
- pyknotid catalogue database: installed
- samples per `N`: `250`
- seed: `20260524`

Commands:

```sh
.venv/bin/python scripts/run_knot_experiment.py \
  --projection-model haar \
  --vertices 5,6,7,8,9,10,11,12 \
  --samples 250 \
  --seed 20260524 \
  --output-dir results/haar_N5-12_250

.venv/bin/python scripts/run_knot_experiment.py \
  --projection-model gaussian \
  --vertices 5,6,7,8,9,10,11,12 \
  --samples 250 \
  --seed 20260524 \
  --output-dir results/gaussian_N5-12_250
```

## Haar Summary

| N | classified | nontrivial | trivial | unknown | nontrivial rate |
|---:|---:|---:|---:|---:|---:|
| 5 | 250 | 0 | 250 | 0 | 0.000 |
| 6 | 250 | 1 | 249 | 0 | 0.004 |
| 7 | 250 | 1 | 249 | 0 | 0.004 |
| 8 | 250 | 8 | 242 | 0 | 0.032 |
| 9 | 250 | 9 | 241 | 0 | 0.036 |
| 10 | 250 | 27 | 223 | 0 | 0.108 |
| 11 | 249 | 29 | 220 | 1 | 0.116 |
| 12 | 248 | 36 | 212 | 2 | 0.145 |

Named nontrivial Haar detections:

| N | knot label | count | rate |
|---:|---|---:|---:|
| 6 | `3_1` | 1 | 0.004 |
| 7 | `3_1` | 1 | 0.004 |
| 8 | `3_1` | 7 | 0.028 |
| 8 | `4_1` | 1 | 0.004 |
| 9 | `3_1` | 7 | 0.028 |
| 9 | `4_1` | 2 | 0.008 |
| 10 | `3_1` | 23 | 0.092 |
| 10 | `4_1` | 2 | 0.008 |
| 10 | `5_2` | 2 | 0.008 |
| 11 | `3_1` | 21 | 0.084 |
| 11 | `5_2` | 4 | 0.016 |
| 11 | `4_1` | 2 | 0.008 |
| 11 | `5_1` | 1 | 0.004 |
| 11 | `8_19;K12n439` | 1 | 0.004 |
| 12 | `3_1` | 31 | 0.124 |
| 12 | `4_1` | 5 | 0.020 |

There were three ambiguous Haar samples where the available invariants matched
multiple catalogue entries, including `0_1`.  Those are counted as `unknown`.
The per-sample CSVs include the raw and simplified Gauss codes for later
manual inspection.

## Gaussian Summary

| N | classified | nontrivial | trivial | unknown | nontrivial rate |
|---:|---:|---:|---:|---:|---:|
| 5 | 250 | 0 | 250 | 0 | 0.000 |
| 6 | 250 | 1 | 249 | 0 | 0.004 |
| 7 | 250 | 6 | 244 | 0 | 0.024 |
| 8 | 250 | 7 | 243 | 0 | 0.028 |
| 9 | 250 | 20 | 230 | 0 | 0.080 |
| 10 | 250 | 24 | 226 | 0 | 0.096 |
| 11 | 250 | 24 | 226 | 0 | 0.096 |
| 12 | 250 | 30 | 220 | 0 | 0.120 |

Named nontrivial Gaussian detections:

| N | knot label | count | rate |
|---:|---|---:|---:|
| 6 | `3_1` | 1 | 0.004 |
| 7 | `3_1` | 6 | 0.024 |
| 8 | `3_1` | 7 | 0.028 |
| 9 | `3_1` | 19 | 0.076 |
| 9 | `4_1` | 1 | 0.004 |
| 10 | `3_1` | 19 | 0.076 |
| 10 | `4_1` | 2 | 0.008 |
| 10 | `5_2` | 2 | 0.008 |
| 10 | `6_1;9_46;K11n67;K11n97;K11n139` | 1 | 0.004 |
| 11 | `3_1` | 19 | 0.076 |
| 11 | `4_1` | 2 | 0.008 |
| 11 | `5_2` | 1 | 0.004 |
| 11 | `6_2` | 1 | 0.004 |
| 12 | `3_1` | 22 | 0.088 |
| 12 | `4_1` | 4 | 0.016 |
| 12 | `5_1` | 2 | 0.008 |
| 12 | `5_2` | 1 | 0.004 |
| 12 | `6_1;9_46` | 1 | 0.004 |

## Interpretation

The experiment supports the expected qualitative picture:

- `N = 5` is unknotted in this sample.
- Nontrivial knots are already detected by `N = 6`; the only possible
  six-stick nontrivial knot is the trefoil, and that is what appears.
- The nontrivial rate rises with `N` over this range.
- Trefoils dominate at these small stick counts, with figure-eight and
  five-crossing knots appearing less often.

The data are not yet a high-precision Haar estimate.  A next pass should use
larger sample sizes, binomial confidence intervals, and possibly a policy for
resolving ambiguous catalogue matches using stronger invariants or manual
Gauss-code reduction.
