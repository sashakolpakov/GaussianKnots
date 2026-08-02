# Chamber Volumes in the Stiefel Model

This note records what is feasible for probabilities and chamber volumes in the
Stiefel/simplex knot model.

## Basic Objects

For `H = 1^perp` with `dim H = N - 1`, the sampling space is

```text
St_3(H) = {Q: H -> R^3 with orthonormal rows}.
```

It has dimension

```text
dim St_3(R^{N-1}) = 3(N - 1) - 6 = 3N - 9.
```

There are three related decompositions, and they should not be conflated.

The embedding discriminant `Delta_N^emb` is the set where the polygon is not an
embedded polygon: an edge degenerates, or two non-adjacent segments meet.  The
complement

```text
St_3(H) \ Delta_N^emb
```

decomposes into connected embedding chambers.  The knot type is constant on
each such chamber.

The determinant-wall arrangement is related but not identical.  For each
non-adjacent edge pair, the condition that the two supporting lines are
coplanar is a determinant equation.  Actual segment intersection also requires
interval inequalities.  Therefore a determinant wall can contain false diagram
walls where the supporting lines are coplanar but the two segments do not meet,
and determinant-wall cells should not be identified with embedding chambers
without a separate connectivity and segment-intersection check.

The order-type buckets used in the scripts are finer still: for the projected
vertices, they record the signs of all four-point orientation determinants.
These buckets are useful for organizing Monte Carlo samples, but they have not
been proved to be connected embedding chambers.

For generic `N=6` configurations, connectedness is not needed for knot
classification: Calvo's theorem and the full 15-sign order type give an exact
finite unknot/trefoil rule.  Connectedness remains relevant if one insists on
enumerating individual chambers rather than integrating the whole knot sector.

For a knot type `K`,

```text
S_N(K) = tau_N^{-1}(K)
```

is the knot-type sector.  It is generally a disconnected union of embedding
chambers, and

```text
p_N(K) = Haar volume(S_N(K)).
```

Therefore `p_N(K)` is the sum of the Haar volumes of all embedding chambers
carrying `K`.  Our scripts estimate this sector mass directly, or through
order-type buckets, not by exact chamber-volume integration.

The total volume of the compact Stiefel or Grassmannian space is known for the
standard invariant Riemannian metrics; equivalently, one can normalize Haar
measure to total mass `1`.  This supplies the denominator of a probability.  The
hard part is the numerator: the volume of the semialgebraic region selected by
embedding conditions, determinant data, and the knot-type condition.

## What We Can Estimate Now

The direct Monte Carlo estimator is already valid:

```text
p_hat_N(K) = #{Q_r classified as K} / M.
```

This estimates the total sector volume without enumerating chambers.  It is the
right primary estimator for the manuscript, because the number of chambers is
large and a knot-type sector can be disconnected.

The current helper

```sh
python3 scripts/repro/summarize_knot_probabilities.py results/haar_N5-12_1000
```

writes

```text
results/haar_N5-12_1000/probability_estimates.csv
```

with per-type counts, `p_hat`, and Wilson 95 percent intervals.  Ambiguous
catalogue candidate lists are not assigned to a single knot type; they are kept
in explicit ambiguous or unknown buckets.

## Computing Individual Chamber Volumes

In principle, an embedding chamber is a semialgebraic subset of `St_3(H)` cut
out by determinant equations together with segment-parameter inequalities.  Its
Haar volume is therefore a well-defined semialgebraic integral.  Exact symbolic
computation is theoretically possible for very small `N`, but it is not a
realistic general method.

Westenberger's random-projection work is the right cautionary comparison:

- even the expectation of average inter-crossing number is obtained by reversing
  projection and configuration-space integration and then exploiting symmetries;
- the second moment of linking number requires a much more complicated
  configuration-space integral;
- in the Petaluma-inspired numerical chapter, the comparison between direct
  sampling and integral evaluation used `10^6` samples or integrand evaluations;
- small separation parameters lead to numerical instability because denominators
  approach zero;
- finite-type invariants are left as future directions precisely because the
  next configuration-space integrals are harder.

Reference:

- Christopher Westenberger, *Knots and Links from Random Projections*,
  arXiv:1602.01484, revised 2019:
  <https://arxiv.org/abs/1602.01484>

This means that "compute the volume of a chamber" should be treated as a
research problem, not as the next routine experiment.  The chamber-volume
integrand would be an indicator of a semialgebraic region in a high-dimensional
Stiefel manifold, not a smooth invariant like curvature or linking number.

For `N = 5`, every embedded polygon is unknotted, so the only meaningful sector
volume is `p_5(0_1) = 1`.  There is still a nontrivial wall arrangement: since
`dim(1^perp) = 4`, projecting to `R^3` is equivalent, up to target rotation, to
choosing a kernel line in `RP^3`.  The five non-adjacent edge pairs give five
determinant hyperplanes whose normals form a regular simplex in `R^4`.  Thus
the `N=5` determinant-wall cells are geometrically nontrivial, but every
embedded five-stick polygon lies in the unknot sector.

For `N = 6`, the first nontrivial sector appears.  This is the smallest plausible
case for attempting chamber-volume experiments.  Since target rotations do not
change the unoriented knot type, the relevant quotient of `St_3(R^5)` is the
Grassmannian

```text
Gr_3(R^5) ~= Gr_2(R^5),
```

which has dimension `6`.  There are

```text
N(N - 3)/2 = 9
```

non-adjacent edge pairs.  Each edge pair gives a determinant wall; in Pluecker
coordinates this is a hyperplane section of the Grassmannian.  Segment-parameter
inequalities then decide which parts of those walls correspond to actual
edge-through-edge events.

The cyclic Hamiltonian polygon has a dihedral `D_6` relabeling symmetry, but
that is not the full symmetry of the wall arrangement.  The N=6 structural
diagnostic

```sh
.venv/bin/python scripts/repro/n6_chamber_analysis.py \
  --samples 50000 \
  --seed 20260604 \
  --output results/n6_full_chamber_analysis_diagnostic.csv
```

computes the projective wall-normal Gram data.  With the wall order

```text
(0,2), (0,3), (0,4), (1,3), (1,4), (1,5), (2,4), (2,5), (3,5),
```

the nine projective walls have absolute pairwise inner products `0` or `1/4`,
and the nonorthogonality graph is the `3 x 3` rook graph.  Hence the full
projective wall permutation group has order `72`, naturally

```text
Aut(3 x 3 rook graph) ~= (S_3 x S_3) semidirect C_2.
```

Under this larger group all nine walls lie in one orbit.  The `D_6` relabeling
subgroup has order `12` and has two wall orbits: six edge pairs of cyclic
separation `2` or `4`, and three opposite-edge pairs of separation `3`.

This symmetry is useful, but it does not by itself compute `p_6(3_1)`.  It gives
equal-volume relations among determinant-wall or order-type cells in the same
orbit.  To turn that into an exact chamber probability, one would still need
either a complete connected-chamber enumeration with trefoil labels, or an
analytic integral over one representative of each chamber orbit.  The targeted
probe therefore records the raw wall signature, the `D_6`-canonical signature,
and the full arrangement-canonical signature as diagnostics.

The full `72`-element group is an actual orthogonal symmetry group of the wall
arrangement.  The consolidated script verifies this by constructing the six common lines
coming from the maximal cliques of the rook graph; these six lines form a
regular simplex in `H`, and every wall permutation lifts to an orthogonal map of
`H` with numerical error below `1e-14`.

However, the full wall-arrangement symmetry is larger than the symmetry of the
labeled Hamiltonian polygon.  The same script starts from the fixed `N=6`
trefoil seed `1062611651` and applies all `72`
orthogonal wall symmetries.  The result is:

```text
D_6 images:        12 trefoils
extra 60 images:   60 unknots
```

Thus the extra wall symmetries do give linear, volume-preserving maps between
wall-arrangement regions, but they do not preserve knot type.  They can relate
regions carrying different knots because they preserve the determinant-wall
geometry without preserving the labeled simplex vertices and cyclic Hamiltonian
order.  For knot-sector enumeration, only symmetries of the polygon model,
currently the `D_6` subgroup, are safe identifications.

There is another practical warning.  The nine determinant signs are not
connected embedding chambers.  Sampling points with the same raw sign pattern
`+-+---+--` produced both unknots and trefoils.  Hence "one point per sign
pattern" is not enough.  One point per connected embedding chamber would be
enough, but the current code does not construct those connected components.

A better finite proxy is the complete 15-sign order type of the six projected
vertices.  For every four-element subset `{a,b,c,d}` of `{0,...,5}`, record

```text
sign det(x_b - x_a, x_c - x_a, x_d - x_a).
```

There are `binom(6,4)=15` such signs.  A sign changes exactly when four
projected vertices become coplanar.  The previous nine signs are only those
four-vertex determinants attached to pairs of non-adjacent polygon edges; the
six remaining signs refine the vertex order type.  For a generic embedded
hexagon, the resulting sign vector determines Calvo's three intersection
numbers and hence the unknot or trefoil type.  This is stronger than the earlier
empirical observation that sampled 15-sign buckets did not mix labels.  The
nine determinant-wall signs alone are still too coarse.

The consolidated script uses fast Haar sampling with no classifier call in the
main loop.  The current grouped run samples `500000` Haar points and buckets
them by `D_6`-canonical 15-sign order type.  The separate exact Calvo script
classifies all `155` observed buckets and gives

```text
p_6(3_1) ~= 1856 / 500000 = 0.003712,
Wilson 95% interval: 0.003547--0.003884.
```

This is an exact classification of the finite generic sample and a numerical
estimate of trefoil-sector mass.  Calvo's criterion already identifies the
complete generic `N=6` trefoil sector from sign data; no connectivity assertion
about a sign realization set is required.  An exact probability would instead
require evaluating the invariant-measure integral of the resulting
semialgebraic sign predicate.  That integration problem remains open here.

If the chambers were ordinary Euclidean simplices in a flat affine chart and the
measure were the induced Lebesgue measure, then the volume problem would indeed
reduce to determinant comparisons.  The `N=6` situation is different: the walls
are linear in Pluecker coordinates, but the domain is the Grassmannian cut out
by Pluecker relations, and the invariant measure is the Grassmannian measure,
not flat ambient Lebesgue measure on the Pluecker coordinates.  Thus the cells
are better regarded as regions of a hyperplane arrangement restricted to a
curved algebraic manifold.  Determinants define the boundaries and symmetries;
they do not automatically make the chamber volumes determinant formulas.

Consequently,

```text
p_6(3_1) = integral_{Gr_2(R^5)} 1_{trefoil sector}(K) dmu(K)
```

is a well-defined six-dimensional indicator integral.  It is a legitimate small
case for numerical integration, stratified sampling, or diagram-signature
enumeration.  It is not obviously a tractable closed-form integral.

Stick number makes the `N=6` topology especially clean: the only possible knot
types are the unknot and the trefoil.  Thus

```text
p_6(0_1) + p_6(3_1) = 1.
```

The manuscript Haar Monte Carlo run gives `p_hat_6(3_1)=2/1000`, with a wide
Wilson interval.  The grouped order-type run above is the sharper numerical
estimate, but its frequency is not an exact chamber-volume calculation.

## Are These Integrals Algebraic Or Reconstructible?

The determinant description is valuable, but it does not usually make the volume
a rational or algebraic expression in the boundary data.

A chamber volume is an integral of the algebraic differential form giving Haar
measure over a semialgebraic region.  Such numbers are periods in the sense of
Kontsevich-Zagier when the defining data are algebraic.  Periods include simple
constants such as `pi`, but also logarithms, elliptic integrals, hypergeometric
integrals, and more complicated configuration-space integrals.

The determinant constraints give part of the boundary data, and the segment
inequalities decide which wall pieces correspond to actual edge intersections.
This boundary information does not, by itself, determine the volume by a
determinant formula.  Even in much
simpler settings:

- the area of a disk is a period involving `pi`;
- spherical polygon volumes involve angles such as `arccos` of algebraic data;
- spherical tetrahedron volumes can involve dilogarithmic-type expressions;
- elliptic periods arise from integrating algebraic forms over algebraic curves.

What can sometimes be reconstructed is a period function, not from zeros alone,
but from additional structure:

1. a known finite-dimensional ansatz, such as rational, algebraic, or
   hypergeometric;
2. high-precision values at enough parameter points;
3. known singular loci and monodromy or a Picard-Fuchs/GKZ differential system;
4. rigorous error bounds.

Without such prior structure, finite numerical values cannot determine an exact
period.  Zero sets alone are much weaker: many different analytic functions have
the same zeros, and a fixed chamber volume is just one real number.  For our
`N=6` problem, the determinant walls and Pluecker relations are enough to define
the period exactly, and they may lead to a useful Picard-Fuchs or hypergeometric
description after introducing parameters.  They do not imply a simple rational
function or determinant expression for the trefoil-sector volume.

## Practical Chamber-Volume Strategy

A feasible numerical strategy is:

1. Sample Haar points `Q` in `St_3(H)`.
2. Compute a complete diagram signature: planar crossing order, crossing signs,
   and over/under data from the omitted height coordinate.
3. Hash samples by this signature.
4. Classify one representative from each high-frequency signature.
5. Estimate the volume of each signature bucket by sample frequency.

This estimates diagram-signature cells, not necessarily connected chambers.  That
is acceptable for probabilities if the signature determines the knot diagram:
disconnected components with the same diagram have the same knot type and their
volumes should be added.

To estimate a genuinely connected chamber volume, one would need an additional
connectivity step:

1. Fix a sample point and its determinant inequalities.
2. Run a manifold MCMC sampler, such as geodesic hit-and-run or projected
   Langevin, constrained to those inequalities.
3. Test whether two sampled points in the same signature bucket can be connected
   without crossing a determinant wall.
4. Estimate relative chamber volumes by restricted sampling or by transition
   counts between adjacent walls.

This is substantially harder than estimating `p_N(K)`.

## Recommended Manuscript Position

The manuscript should claim:

- `p_N(K)` is a Haar volume of a knot-type sector.
- Monte Carlo directly estimates that sector volume.
- Chamber volumes are a possible refinement, but not needed for estimating
  `p_N(K)`.
- Exact chamber-volume computation is only plausible in very small cases.

The useful next experiment is not exact volume integration.  It is a
diagram-signature census for Haar samples, because it tells us whether the
observed knot-type mass is concentrated in a few large regions or spread across
many small regions.
