# Reproducibility

Run commands from the repository root unless stated otherwise.  The manuscript
[`manuscript/knot_projections.tex`](manuscript/knot_projections.tex) is the
source of truth for the theorem statements and numerical tables.

## Python environment

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

The main requirements pin the exact `pyknotid` fork revision used by this
repository.  Install `requirements-analysis.txt` when SciPy's
Kolmogorov--Smirnov diagnostic is desired:

```sh
python -m pip install -r requirements-analysis.txt
```

## Exact-law checks

```sh
python scripts/repro/validate_exact_projection_laws.py \
  --crossing-samples 10000 \
  --beta-samples 50000 \
  --batch-size 500
```

The script checks the exact fixed-pair planar crossing probability
`(2/pi) asin(1/3)`, couples centered Gaussian samples to their Haar row
factor sample by sample, and checks the `Beta(3/2,(N-4)/2)` pair-distance
marginal. This is a
numerical audit of exact formulas proved in the manuscript.

## Exact N=6 finite classifier

The default command reads only tracked files under `data/reference/n6/`:

```sh
python scripts/repro/calvo_n6_exact.py \
  --expected-trefoils 1856 \
  --validate-standard-direct
```

It applies Calvo's finite sign rule to all 155 observed order-type buckets,
recovers 1,856 trefoils among 500,000 samples, and compares the rule against
the two tracked standard 1000-sample N=6 CSV files.

The stringent audit replays the seeded 500,000-sample stream and recomputes the
100 nearest-wall samples at 100-digit decimal precision:

```sh
python scripts/repro/calvo_n6_exact.py \
  --expected-trefoils 1856 \
  --validate-standard-direct \
  --audit-signs
```

The finite classification is exact for each generic sign vector, conditional
on Calvo's cited theorem. The proportion `1856/500000` remains a Monte Carlo
estimate of Haar mass.

## Manuscript Monte Carlo

```sh
python scripts/run_knot_experiment.py \
  --projection-model haar \
  --vertices 5,6,7,8,9,10,11,12 \
  --samples 1000 \
  --seed 20260604 \
  --output-dir results/haar_N5-12_1000

python scripts/run_knot_experiment.py \
  --projection-model gaussian \
  --vertices 5,6,7,8,9,10,11,12 \
  --samples 1000 \
  --seed 20260604 \
  --output-dir results/gaussian_N5-12_1000
```

The Haar and Gaussian runs have the same population knot-type law but are not
expected to agree sample by sample.  Outputs record seeds, classifier state,
unknown classifications, knot invariants, and metric summaries.

The volume-first order-type run is:

```sh
python scripts/repro/order_type_grouped_volume.py \
  --vertices 6,7,8 \
  --samples 500000 \
  --seed 20260604 \
  --classify-top-groups 200 \
  --checks-per-group 3 \
  --direct-classify-samples 1000 \
  --output-dir results/order_type_grouped_volume_N6-8_500k
```

For N=6, use `calvo_n6_exact.py` for final bucket labels. For N=7 and N=8,
the grouped labels remain exploratory classifier diagnostics and are not a
proof of connected chamber structure.

## Figure and manuscript

Regenerate the fixed diagram source:

```sh
python scripts/repro/make_stiefel_example_diagrams.py
pdflatex -interaction=nonstopmode -halt-on-error \
  -output-directory manuscript/figures \
  manuscript/figures/stiefel_unknot_trefoil.tex
```

Build the manuscript:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  -cd manuscript/knot_projections.tex
```

Both the TeX source and reviewed PDF are committed.

## Formal verification and documentation

```sh
cd formal
lake exe cache get
lake build GaussianKnots
```

See [`FORMAL_VERIFICATION.md`](FORMAL_VERIFICATION.md) for the checked and
external boundary.

The unfinished-proof guard is:

```sh
rg -n '(^|[[:space:]])(sorry|admit)([[:space:]]|$)|^[[:space:]]*axiom[[:space:]]|native_decide' \
  --glob '*.lean' formal
```

```sh
python -m pip install -r docs/requirements.txt
python -m sphinx -W --keep-going -b html docs docs/_build/html
```
