# GaussianKnots

Reproducibility package for projected-simplex Hamiltonian-cycle stick-knot
experiments.  It projects the $N$ labelled simplex vertices to three-space,
orders them along the cycle $1\to2\to\cdots\to N\to1$, and estimates the
knot-type frequencies of the resulting closed polygon.

This package is intentionally script-first: no notebooks are required.

## Install

From this directory:

```sh
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
```

This repository uses the `sashakolpakov/pyknotid` fork and requests its
`performance` extra, which installs the Numba JIT backend advertised by that
fork.  If the fast backend is not available, the scripts still run through the
installed pyknotid path and report backend status in `run_metadata.json`.

`pyknotid` catalogue identification may require its separate knot database.  If
your run reports catalogue lookup failures, install the database with:

```sh
python3 scripts/run_knot_experiment.py --download-pyknotid-db --samples 1 --vertices 6
```

or directly:

```sh
python3 -c "from pyknotid.catalogue.getdb import download_database; download_database()"
```

## Run

Quick run:

```sh
python3 scripts/run_knot_experiment.py --vertices 6,7,8 --samples 25 --output-dir results/quick
```

Larger run:

```sh
python3 scripts/run_knot_experiment.py --vertices 6,7,8,10,12 --samples 1000 --seed 20260524 --output-dir results/main
```

Compare the Haar model with Gaussian vertices:

```sh
python3 scripts/run_knot_experiment.py --projection-model haar --vertices 6,7,8,10 --samples 250 --output-dir results/haar_250
python3 scripts/run_knot_experiment.py --projection-model gaussian --vertices 6,7,8,10 --samples 250 --output-dir results/gaussian_250
```

Current manuscript-scale runs use $N=5,\ldots,12$ and 1000 samples:

```sh
python3 scripts/run_knot_experiment.py \
  --projection-model haar \
  --vertices 5,6,7,8,9,10,11,12 \
  --samples 1000 \
  --seed 20260604 \
  --output-dir results/haar_N5-12_1000

python3 scripts/run_knot_experiment.py \
  --projection-model gaussian \
  --vertices 5,6,7,8,9,10,11,12 \
  --samples 1000 \
  --seed 20260604 \
  --output-dir results/gaussian_N5-12_1000
```

If you only want to verify polygon generation on a machine without pyknotid:

```sh
python3 scripts/run_knot_experiment.py --vertices 6 --samples 3 --allow-missing-pyknotid
```

## N=6 Wall/Order-Type Analysis

The first nontrivial stick count is $N=6$: only the unknot and trefoil can
occur.  The reproducibility scripts in `scripts/repro/` analyze this case more
finely than direct classifier Monte Carlo.

The current volume-first grouped entry point is:

```sh
python3 scripts/repro/order_type_grouped_volume.py \
  --vertices 6,7,8 \
  --samples 500000 \
  --seed 20260604 \
  --classify-top-groups 200 \
  --checks-per-group 3 \
  --direct-classify-samples 1000 \
  --output-dir results/order_type_grouped_volume_N6-8_500k
```

It records:

- the nine edge-pair determinant walls;
- the wall-normal Gram symmetry group, which has order 72 at `N=6`;
- the `D_6` subgroup of order 12 that comes from relabelling the Hamiltonian
  six-cycle and preserves knot type;
- the remaining 60 wall symmetries, which preserve the wall arrangement but can
  change knot type;
- the 15-sign order type of the six projected vertices, with one sign for each
  4-subset of vertices;
- classified checks in high-volume `D_N` order-type buckets; and
- volume Monte Carlo estimates over order-type buckets.

The 15 signs are

```math
\operatorname{sgn}\det(x_b-x_a,\;x_c-x_a,\;x_d-x_a),
\qquad \{a,b,c,d\}\subset\{0,\ldots,5\}.
```

There are $\binom64=15$ such signs.  They refine the nine non-adjacent
edge-pair signs: the nine signs are too coarse, while the 15-sign order type is
the current bucket proxy used for the $N=6$ volume estimate.

For `N=6`, this run checks all `155` observed order-type buckets.  The single
checked trefoil bucket gives

```text
p_6(3_1) ~= 1856/500000 = 0.003712
Wilson 95% interval: 0.003547--0.003884
```

This is consistent with the current direct Monte Carlo table.  The manuscript
Haar run has `2/1000 = 0.002` trefoils at `N=6`, with Wilson interval
`0.000549--0.007263`; the raw Gaussian comparison has `5/1000 = 0.005`, also
with a broad interval.  Both intervals contain the order-type estimate.  The
direct `N=6` run is therefore a useful sanity check, but it is too small to
measure `p_6(3_1)` sharply.

This should be read as a numerical order-type bucket estimate, not an exact
connected-chamber volume computation.  A proof-level exact computation would
still require proving that the relevant 15-sign cells are connected embedding
chambers, that segment-intersection inequalities introduce no knot-relevant
splitting or merging, and that the identified `D_6` trefoil buckets give the
complete trefoil sector.

The N=6-only structural diagnostic verifies the orthogonal lift of the 72 wall
symmetries:

```sh
python3 scripts/repro/n6_chamber_analysis.py \
  --samples 50000 \
  --seed 20260604 \
  --output results/n6_full_chamber_analysis_diagnostic.csv
```

The lower-level N=6 Grassmannian sampler is:

```sh
python3 scripts/repro/run_n6_grassmannian_probe.py \
  --samples 2048 \
  --sampler halton \
  --seed 20260604 \
  --output-dir results/n6_grassmannian_halton_2048
```

The exploratory N=7 and N=8 wall/order-type probes use the shared package helper
`gaussian_knots/general_chamber_geometry.py`:

```sh
python3 scripts/repro/n7_chamber_probe.py \
  --samples 5000 \
  --classify-reps 200 \
  --classify-samples 200 \
  --skip-full-wall-symmetry \
  --output-dir results/n7_chamber_probe_5000

python3 scripts/repro/n8_chamber_probe.py \
  --samples 5000 \
  --classify-reps 200 \
  --classify-samples 200 \
  --skip-full-wall-symmetry \
  --output-dir results/n8_chamber_probe_5000
```

These probes estimate `D_N` order-type bucket masses and compare limited
direct classifier samples with representative labels.  They are diagnostics,
not proof-level chamber decompositions.

In the same grouped run, the direct 1000-sample checks gave:
`N=6`: `0_1` 993, `3_1` 7;
`N=7`: `0_1` 982, `3_1` 16, `4_1` 2;
`N=8`: `0_1` 977, `3_1` 21, `4_1` 1, `5_1` 1.
The grouping showed 155, 54912, and 497097 observed order-type buckets,
respectively, in 500k samples.

Probability summaries with Wilson intervals can be regenerated from experiment
outputs:

```sh
python3 scripts/repro/summarize_knot_probabilities.py results/haar_N5-12_1000
python3 scripts/repro/summarize_knot_probabilities.py results/gaussian_N5-12_1000
```

## Outputs

Each run writes:

- `run_metadata.json`: parameters and pyknotid import/fast-backend status.
- `samples_N{N}.csv`: one record per sampled polygon.
- `summary.csv`: counts and rates by stick count.
- `type_counts.csv`: empirical knot-label frequencies by stick count.

Important summary columns:

- `classified`: samples with a detected nontrivial knot or a detected unknot.
- `nontrivial`: samples whose pyknotid type or invariants detect nontriviality.
- `trivial`: samples detected as unknotted.
- `unknown`: samples not classified by the available API/database/invariants.
- `nontrivial_rate_known`: `nontrivial / classified`.
- `nontrivial_lower_bound_rate`: `nontrivial / samples`.
- `cycle_distortion_mean`: average max/min length distortion along the
  Hamiltonian cycle.
- `pair_distance_distortion_mean`: average max/min distance distortion over all
  $N(N-1)/2$ projected simplex vertex pairs.
- `pair_abs_ratio_rms_mean`: average all-pair RMS projected distance divided by
  the original simplex distance $\sqrt{2}$.
- `pair_rms_scale_to_simplex_mean`: average global scale that would match each
  sample's all-pair RMS distance to $\sqrt{2}$.
- `pair_normalized_ratio_min_mean` and `pair_normalized_ratio_max_mean`: average
  scale-free all-pair distance range after RMS normalization.

Per-sample records also include `gauss_code`, `simplified_gauss_code`,
Alexander values at roots 2, 3, and 4, Vassiliev degree 2/3 values where
available, raw/simplified crossing counts, Hamiltonian-edge length statistics,
and all-pair distance deformation statistics.  These fields are useful when
catalogue identification is ambiguous or when comparing metric deformation.

When the catalogue database is unavailable, the code still tries pyknotid
invariants.  Nontriviality is counted only when an invariant detects it, so
`nontrivial_lower_bound_rate` remains conservative.  Samples with trivial
values for the computed invariants but no positive unknot identification are
reported as `unknown`.

## First Numeric Run

A first 250-sample run for $N=5,\ldots,12$, comparing Haar sampling with
Gaussian vertices, is summarized in:

- `reports/haar_vs_gaussian_N5-12_250.md`

## Sphinx Docs

The repository includes a small Sphinx site in `docs/` with MathJax enabled for
LaTeX rendering:

```sh
pip install -r docs/requirements.txt
sphinx-build -b html docs docs/_build/html
```

The GitHub Actions workflow in `.github/workflows/docs.yml` builds and deploys
the docs through GitHub Pages when Pages is enabled for the repository.

## Smoke Test

```sh
python3 scripts/smoke_test.py
```

The smoke test checks Gaussian polygon generation and the pyknotid adapter.  It
passes without pyknotid installed by verifying the adapter's missing-dependency
path.

## Model Notes

For a given $N$, the default `--projection-model haar` samples a
Haar-distributed row-orthonormal projection of the simplex:

```math
Q Q^T = I_3,\qquad \sum_i Q e_i = 0,\qquad x_i = Q e_i \in \mathbb{R}^3.
```

This directly samples the Grassmann/Stiefel model governing knot type in the
manuscript.  The optional `--projection-model gaussian` instead samples
Gaussian vertices

```math
x_i = P e_i \in \mathbb{R}^3,\qquad P_{ab}\sim\mathcal{N}(0,1/3).
```

The manuscript proves these two models have the same knot-type law: the
Gaussian map factors as $P=AQ$, and the invertible linear map $A$ preserves
ambient isotopy type.  They differ for metric statistics such as edge lengths.

The parameter space is cut by discriminant walls where non-adjacent edges meet.
Knot type is locally constant on each chamber, and a generic wall crossing
changes one diagram crossing.

Closed polygonal knots with fewer than six sticks are marked unknotted by the
stick-number obstruction.  For $N\ge 6$, classification is delegated to
pyknotid where available.

## pyknotid API Assumptions

The adapter imports `pyknotid.spacecurves.Knot`, constructs `Knot(points,
verbose=False)`, and treats it as a closed curve using pyknotid's default closed
`Knot` methods.  It dynamically inspects method signatures before passing
optional keywords such as `try_cython`, `roots`, `vassiliev_2`, and
`vassiliev_3`, because released pyknotid APIs differ.

The preferred fast path is the fork's Numba JIT support, installed through the
`performance` extra.  The adapter also requests older pyknotid fast-path
keywords such as `try_cython=True` when an installed API exposes them.  If those
backends are absent or an API does not accept an option, the run falls back to
the installed Python path and records messages in the sample CSV.
