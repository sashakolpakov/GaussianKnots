# Random Knots via Stiefel Manifolds

GaussianKnots is the public manuscript and reproducibility repository for
*Random Knots via Stiefel Manifolds* by Alexander Kolpakov and Igor Rivin.

- [Compiled manuscript](manuscript/knot_projections.pdf)
- [LaTeX source](manuscript/knot_projections.tex)
- [Reproducibility guide](REPRODUCIBILITY.md)
- [Formal-verification scope](FORMAL_VERIFICATION.md)
- [Published documentation](https://sashakolpakov.github.io/GaussianKnots/)

The LaTeX manuscript is the source of truth for mathematical statements. The
README and Sphinx site provide navigation and reproducibility instructions;
if wording differs, defer to the manuscript.

## Main results

The model projects the $N$ labelled vertices of the standard simplex to
$\mathbb R^3$ and joins them in Hamiltonian order. The manuscript proves,
among other results:

- centered Gaussian projections and Haar row-orthonormal projections have the
  same knot-type law, although their metric shapes differ;
- the law has full support on precisely the knot types with stick number at
  most $N$;
- for $N>4$, $\lVert Q(e_i-e_j)\rVert^2/2 \sim \mathrm{Beta}(3/2,(N-4)/2)$;
- a generic planar shadow has expected crossing count
  $\mathbb E C_N=N(N-3)\arcsin(1/3)/\pi$;
- its crossing count satisfies
  $\Pr(|C_N-\mathbb E C_N|\geq t)\leq
  2\exp(-t^2/(2N(N-3)^2))$; and
- for a generic embedded hexagon, the 15 tetrahedral orientation signs determine
  Calvo's three oriented algebraic intersection numbers and hence give an exact
  unknot/right-trefoil/left-trefoil decision.

Applying the exact finite $N=6$ rule to the tracked 500,000-sample experiment
classifies 1,856 samples as trefoils:

$$
\widehat p_6(3_1)=1856/500000=0.003712.
$$

The classification of every generic sampled sign vector is exact, conditional
on Calvo's cited topological theorem. The frequency remains a Monte Carlo
estimate, not an exact Haar probability or chamber volume.

## Quick verification

Create an environment and install the computational dependencies:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

Run the fast exact $N=6$ reclassification against tracked reference data:

```sh
python scripts/repro/calvo_n6_exact.py \
  --expected-trefoils 1856 \
  --validate-standard-direct
```

Check the exact planar-crossing and Haar pair-distance laws numerically:

```sh
python -m pip install -r requirements-analysis.txt
python scripts/repro/validate_exact_projection_laws.py \
  --crossing-samples 10000 \
  --beta-samples 50000 \
  --batch-size 500
```

Build the Lean development:

```sh
cd formal
lake build GaussianKnots
```

The Lean development defines the actual iid Gaussian measure, crossing events,
crossing statistic, vertex replacement, and one 15-determinant hexagon
signature. It proves the metric identities, raw/centered crossing invariance,
the exact replacement sensitivity, and the five-point circuit criterion from
a concrete nonzero determinant. From a concrete four-point probability input
it derives the exact mean; from a typed Doob certificate it derives both tail
bounds using mathlib Azuma--Hoeffding. Gaussian polar/Wishart transfer, the
four-point formula, construction of the Doob certificate, PL isotopy, and
Calvo's topological theorem remain explicit external inputs. See
[FORMAL_VERIFICATION.md](FORMAL_VERIFICATION.md).

## Run a knot experiment

```sh
python scripts/run_knot_experiment.py \
  --projection-model haar \
  --vertices 5,6,7,8 \
  --samples 1000 \
  --seed 20260604 \
  --output-dir results/haar_example
```

Named knot identification uses the pinned
[`sashakolpakov/pyknotid`](https://github.com/sashakolpakov/pyknotid)
fork. Its optional catalogue database can be installed with:

```sh
python scripts/run_knot_experiment.py \
  --download-pyknotid-db \
  --samples 1 \
  --vertices 6
```

Generated result directories are intentionally ignored. Compact inputs needed
to audit the manuscript's $N=6$ calculation are versioned under
[`data/reference/n6/`](data/reference/n6/README.md).

## Repository map

- `manuscript/`: authoritative TeX source, compiled PDF, and figure sources;
- `formal/`: pinned Lean 4/mathlib project;
- `gaussian_knots/`: polygon generation, geometry, and classifier adapter;
- `scripts/repro/`: exact-law checks and manuscript experiments;
- `data/reference/`: compact tracked inputs for stated computational checks;
- `docs/`: Sphinx guide generated around the manuscript;
- [`CHAMBER_VOLUMES.md`](CHAMBER_VOLUMES.md): limits and prospects for chamber
  integration;
- [`STUDY_RANDOM.md`](STUDY_RANDOM.md): comparison with other random-knot
  models.

The earlier 250-sample pilot is retained as a clearly labelled historical
report in [`reports/`](reports/haar_vs_gaussian_N5-12_250.md); it is not the
source for the manuscript's 1000-sample tables.

## Authors and acknowledgment

The authors are Alexander Kolpakov and Igor Rivin. OpenAI's GPT-5.6-sol was
used to help improve the mathematical results and generate and refine code
used in the paper, including Lean files. See
[`ACKNOWLEDGMENTS.md`](ACKNOWLEDGMENTS.md) for the precise statement.

Citation metadata are available in [`CITATION.cff`](CITATION.cff).
