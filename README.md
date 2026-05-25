# GaussianKnots

Reproducibility package for the Gaussian projected Hamiltonian-cycle stick-knot
experiment.  It samples the columns of a `3 x N` Gaussian projection matrix,
orders them along the cycle `1-2-...-N-1`, and estimates how often the resulting
closed polygon is detected as a nontrivial knot.

This package is intentionally script-first: no notebooks are required.

## Install

From this directory:

```sh
python3 -m venv .venv
. .venv/bin/activate
pip install -r requirements.txt
```

`pyknotid` can use optional compiled helpers for some space-curve operations.
The requirements file lists `Cython` before `pyknotid` so a fresh install has a
chance to build those helpers.  If the compiled helpers are not available, the
scripts still request the fastest supported pyknotid API and report the backend
status in `run_metadata.json`.

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

If you only want to verify polygon generation on a machine without pyknotid:

```sh
python3 scripts/run_knot_experiment.py --vertices 6 --samples 3 --allow-missing-pyknotid
```

## Outputs

Each run writes:

- `run_metadata.json`: parameters and pyknotid import/fast-backend status.
- `samples_N{N}.csv`: one record per sampled polygon.
- `summary.csv`: counts and rates by stick count.

Important summary columns:

- `classified`: samples with a detected nontrivial knot or a detected unknot.
- `nontrivial`: samples whose pyknotid type or invariants detect nontriviality.
- `trivial`: samples detected as unknotted.
- `unknown`: samples not classified by the available API/database/invariants.
- `nontrivial_rate_known`: `nontrivial / classified`.
- `nontrivial_lower_bound_rate`: `nontrivial / samples`.

When the catalogue database is unavailable, the code still tries pyknotid
invariants.  Nontriviality is counted only when an invariant detects it, so
`nontrivial_lower_bound_rate` remains conservative.  Samples with trivial
values for the computed invariants but no positive unknot identification are
reported as `unknown`.

## Smoke Test

```sh
python3 scripts/smoke_test.py
```

The smoke test checks Gaussian polygon generation and the pyknotid adapter.  It
passes without pyknotid installed by verifying the adapter's missing-dependency
path.

## Model Notes

For a given `N`, the script samples vertices

```text
x_i = P e_i in R^3,   P_{ab} ~ N(0, 1/3),
```

then closes the Hamiltonian cycle through consecutive vertices.  The `1/sqrt(3)`
scale matches the projection convention in the paper; knot type is unchanged by
positive global rescaling.

Closed polygonal knots with fewer than six sticks are marked unknotted by the
stick-number obstruction.  For `N >= 6`, classification is delegated to
pyknotid where available.

## pyknotid API Assumptions

The adapter imports `pyknotid.spacecurves.Knot`, constructs `Knot(points,
verbose=False)`, and treats it as a closed curve using pyknotid's default closed
`Knot` methods.  It dynamically inspects method signatures before passing
optional keywords such as `try_cython`, `roots`, `vassiliev_2`, and
`vassiliev_3`, because released pyknotid APIs differ.

The preferred fast path is pyknotid's Cython-backed crossing/invariant support
where exposed through `try_cython=True`.  If that backend is absent or an API
does not accept the option, the run falls back to the installed Python path and
records messages in the sample CSV.

