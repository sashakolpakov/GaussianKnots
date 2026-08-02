# Lean verification for GaussianKnots

This pinned Lean 4.32/mathlib project gives conditional certification of the
principal algebraic and probability deductions in
`manuscript/knot_projections.tex`. Deep analytic and topological inputs are
ordinary typed theorem parameters; they are not project axioms.

## Build and source audit

```sh
lake exe cache get
lake build GaussianKnots

rg -n '(^|[[:space:]])(sorry|admit)([[:space:]]|$)|^[[:space:]]*axiom[[:space:]]|native_decide' \
  --glob '*.lean' .
```

The source scan is empty. `GaussianKnots.lean` imports the complete target and
runs `#print axioms` on the principal declarations. The expected dependencies
are only `propext`, `Classical.choice`, and `Quot.sound` from Lean/mathlib.

## Gaussian and crossing layer

`GaussianProbabilityCore.lean` defines the actual iid product law of a
`3 × N` matrix with coordinate distribution `N(0,1/3)`. Lean proves
row/column permutation invariance, coordinate means and variances,
independence, measurability of row centering, zero centered row sums,
pointwise preservation of inter-column distances, and

```text
E ||X_j-X_k||² = 2       (j ≠ k)
```

for both the raw and centered laws.

`CrossingProbabilityCertificate.lean` uses the edge set of the complement of
`cycleGraph N` as the actual type of nonincident edge pairs. It defines strict
proper crossing for the planar Gaussian vertices, proves measurability and the
exact `N(N-3)/2` slot count, and constructs the complete crossing statistic.
It proves pointwise invariance under row centering.

Lean also defines replacement of one matrix column and proves:

- every crossing indicator outside the two incident-edge families is
  unchanged; and
- the actual statistic satisfies
  `|C(x with vertex v replaced by p)-C(x)| ≤ 2(N-3)`.

From the concrete four-point crossing probability hypothesis, Lean derives
the exact raw and centered mean

```text
E C_N = N(N-3) / pi * arcsin(1/3).
```

A `CrossingDoobCertificate` contains an actual filtration and increments, the
identity expressing the centered crossing statistic as their sum, and the
conditional sub-Gaussian estimates. Mathlib's Azuma--Hoeffding theorem then
derives both the printed two-sided bound and
`2 exp(-epsilon² N/8)` after normalization. Constructing this certificate from
independent vertex exposure and the replacement theorem remains external; the
tail constants themselves are Lean deductions.

## Six-stick determinant classifier

`CircuitGenericity.lean` proves, from one nonzero affine determinant, that
every affine dependence of five points is a scalar multiple of the cofactor
circuit. Consequently, the strict three-versus-two sign split is equivalent
to actual open-segment/open-triangle intersection; no free
one-dimensionality conclusion is assumed.

`HexagonSignature.lean` starts from one genuine 15-determinant signature of six
vertices and derives all six Calvo records, including their oriented
contributions. Under the explicit nonvanishing condition, each executable
nonzero contribution is proved equivalent to its geometric piercing. Lean
also checks reflection, the canonical trefoil deltas, and the final finite
decision function.

## External theorem boundary

The following mathematics remains cited rather than re-proved:

- the Gaussian four-point/orthant formula;
- construction of the Doob certificate from iid vertex exposure;
- Gaussian polar/Wishart transfer from the centered raw model to Haar;
- the Haar beta marginal;
- PL stability and ambient-isotopy extension; and
- Calvo's topological theorem interpreting the three intersection numbers as
  knot type.

Thus the concrete determinant, crossing, sensitivity, expectation, tail, and
decision deductions are kernel checked. The analytic transfer theorems and
ambient-isotopy classification are not claimed to have been formalized.
