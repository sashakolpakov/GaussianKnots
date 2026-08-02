# Formal verification

The pinned Lean 4.32/mathlib project in [`formal/`](formal/) conditionally
certifies the manuscript's principal determinant, crossing, and probability
deductions on their actual mathematical objects.

## Checked content

| Layer | Lean status |
|---|---|
| Gaussian model | Defines the iid `N(0,1/3)` product measure on `3 × N` matrices; proves row/column invariance, means, variances, independence, centering identities, and expected squared distance `2` |
| Crossing slots | Uses actual unordered nonincident pairs from the complement of the cycle graph; proves exactly `N(N-3)/2` slots |
| Crossing events | Defines measurable strict segment crossing for the concrete planar Gaussian vertices and the full crossing statistic |
| Centering | Proves every crossing event and the full statistic are unchanged by row centering |
| Vertex replacement | Defines replacement of one column; proves all indicators outside the incident-edge families are unchanged and `|Delta C_N| ≤ 2(N-3)` |
| Exact mean | Derives `E C_N = N(N-3) arcsin(1/3)/pi` from the concrete four-point crossing-probability hypothesis |
| Concentration | From an actual typed Doob certificate, derives the exact printed two-sided tail and normalized `2 exp(-epsilon²N/8)` bound using mathlib Azuma--Hoeffding |
| Five-point circuit | Proves uniqueness of affine dependence from a nonzero determinant and equivalence between the circuit sign split and open segment--triangle intersection |
| Six-stick signature | Derives all six Calvo records and oriented contributions from one 15-determinant signature and checks the executable trefoil decision |

The exact mean and tail probabilities are downstream theorems, not free
numerical hypotheses. The mean theorem takes the four-point probability for
the actual measurable crossing events. The tail theorem takes a
`CrossingDoobCertificate` containing a filtration, increments, their sum
identity, and conditional sub-Gaussian estimates.

## External boundary

Lean does not re-prove:

- the Gaussian four-point/orthant formula;
- construction of the Doob certificate from independent vertex exposure and
  the now-certified replacement estimate;
- Gaussian polar/Wishart transfer to the Haar model;
- the Haar beta marginal;
- PL stability and ambient-isotopy extension; or
- Calvo's final topological classification theorem.

The six-stick determinant classifier is therefore kernel checked, conditional
on Calvo's interpretation as ambient-isotopy knot type. Likewise, the crossing
bookkeeping and exact constants are checked, conditional on the two stated
analytic bridges.

## Build and axiom audit

```sh
cd formal
lake exe cache get
lake build GaussianKnots

rg -n '(^|[[:space:]])(sorry|admit)([[:space:]]|$)|^[[:space:]]*axiom[[:space:]]|native_decide' \
  --glob '*.lean' .
```

The source scan is empty. The root module imports the complete project and
runs `#print axioms` on the principal declarations. Only the standard
Lean/mathlib foundations `propext`, `Classical.choice`, and `Quot.sound` are
reported; there are no project-specific axioms.

For the module-by-module inventory, see [`formal/README.md`](formal/README.md).
