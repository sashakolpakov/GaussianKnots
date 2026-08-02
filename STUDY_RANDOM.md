# Study Note: Random-Knot Models and Comparisons

This note separates three issues that must stay distinct in the manuscript:

1. what we prove about the Stiefel/simplex model,
2. what we have actually computed,
3. which earlier random-knot models are close enough for meaningful comparison.

## Claims Ledger

### Mathematical Claims In The Manuscript

These are structural statements, not Monte Carlo conclusions.

- A generic projection of the labeled simplex vertices to `R^3`, joined in
  Hamiltonian order, is an embedded polygonal knot.
- The raw centered Gaussian projection and the Haar row-orthonormal Stiefel
  projection have the same exact knot-type law.  The polar factor removes
  singular-value shear by an invertible linear map in target space.
- The knot-type map on the embedded part of `St_3(1^perp)` is Borel measurable.
- The knot-type sectors
  `S_N(K) = tau_N^{-1}(K)` are unions of embedding chambers.
- The support is exactly controlled by stick number:
  `p_N(K) > 0` if and only if `stick(K) <= N`.

What this does not claim:

- It does not say the Stiefel law is uniform over knots with `stick(K) <= N`.
- It does not say finite Monte Carlo samples should see every knot with positive
  probability.
- It does not identify the Stiefel law with the standard random grid, random
  planar diagram, random polygon, or random jump laws.
- Outside the generic `N=6` Calvo sign classifier, it does not make a
  classifier-independent claim about numerical knot labels in the CSV files.

### Computations We Have Run

The code in this repository samples the
Hamiltonian cycle through projected simplex vertices and classifies the resulting
spatial polygon with `pyknotid`.

Historical 250-sample pilot:

- `results/haar_N5-12_250`
- `results/gaussian_N5-12_250`
- `N = 5,...,12`
- `M = 250` per `N`
- seed `20260524`

Current larger run:

- `results/haar_N5-12_1000`
- `results/gaussian_N5-12_1000`
- `N = 5,...,12`
- `M = 1000` per `N`
- seed `20260604`

These computations estimate the finite-`N` Stiefel probabilities

```text
p_N(K) = Haar{Q in St_3(1^perp): Gamma(Q) has knot type K}.
```

The estimator is

```text
p_hat_N(K) = (1 / M) * #{samples classified as K}.
```

The CSV data should be reported with sample size, random seed, projection model,
classifier/backend metadata, unknown/ambiguous classifications, and binomial
confidence intervals.  The raw Gaussian run is a numerical stability comparison,
not a different exact knot-type law.

### External Data We Use

Stick-number comparisons rely on external tables.  KnotInfo defines the stick
number as the minimum number of vertices in a polygonal description of a knot and
records sources for exact values and bounds:

- KnotInfo stick number table:
  <https://knotinfo.org/descriptions/polygon_index.html>
- Eddy-Shonkwiler used Monte Carlo sampling of confined polygons for stick-number
  bounds and generated 220 billion polygons:
  <https://arxiv.org/abs/1909.00917>
- Cantarella-Rechnitzer-Schumacher-Shonkwiler give 2025 upper-bound updates and
  a table of bounds through 13 crossings:
  <https://arxiv.org/abs/2508.18263>

These sources should be cited when the manuscript uses `stick(K)` data or proposes
a stick-number coverage benchmark.

## Relevant Random-Knot Models

### 1. Random Projections From High Dimension

Westenberger studies random knots by fixing a curve in `R^n`, `n > 3`, and
orthogonally projecting it to random 3-dimensional subspaces:

- Christopher Westenberger, *Knots and Links from Random Projections*:
  <https://arxiv.org/abs/1602.01484>

This is the closest conceptual predecessor.  Our model is a finite polygonal
specialization: the fixed curve is the Hamiltonian cycle through standard simplex
vertices, and the random projection is represented explicitly by a Haar point of
`St_3(1^perp)`.  The main difference is that we exploit the simplex and Stiefel
normalization to prove the Gaussian-Haar equivalence and exact support statement
in terms of stick number.

Westenberger should also calibrate our ambitions.  His computations are not
simple chamber-volume calculations.  The work derives projection-averaged
configuration-space integrals for curvature, average inter-crossing number, and
the second moment of linking number.  The second moment calculation requires
substantial symmetry reduction, and the numerical Petaluma-style comparisons use
`10^6` random samples or integrand evaluations.  The paper explicitly notes
numerical instability when strands are close and denominators become small, and
leaves more difficult finite-type-invariant integrals as future work.  This is
strong evidence that exact Stiefel chamber volumes are not a short-term
computational target.

Useful comparison:

- expected crossing count under random planar views,
- finite-type invariants as observables of the projected knot,
- how the choice of fixed high-dimensional curve changes the resulting knot law.
- direct sampling of knot-type sectors versus configuration-space integral
  computations of smoother invariants.

### 2. Polygonal Walks And Random Polygons In `R^3`

Classical models sample the spatial curve directly: self-avoiding lattice
polygons, equilateral polygons, Gaussian polygons, confined polygons, and related
Markov-chain or symplectic samplers.  Even-Zohar's survey is the best entry point:

- Chaim Even-Zohar, *Models of Random Knots*:
  <https://arxiv.org/abs/1711.10470>

For polygonal walks, Even-Zohar records Gaussian/equilateral models, the
Frisch-Wasserman-Delbruck unknotting question, and the role of local connected
summands.  The survey also notes empirical laws for fixed knot frequencies in
several polygonal/grid models.

We should compare to these models carefully.  Our points are not sampled as an
ordinary closed random walk.  They are columns of one random linear map applied to
simplex vertices, so the vertices are exchangeable but globally dependent.  The
Haar Stiefel model also enforces a compact normalization:

```text
sum_i x_i = 0,       X^T X = I_3.
```

Useful comparison:

- nontrivial rate by stick count `N`,
- knot-type frequency tables for small `N`,
- edge-length and all-pair distance deformation,
- finite-type invariants or determinant distributions,
- whether local summands are visible or rare.

### 3. Random Jump / Confined Polygon Models

The random jump model samples points independently in a bounded region, such as a
cube or ball, and joins them cyclically.  It is geometrically closer to our model
than small-step walks because typical edge lengths are comparable to the diameter
of the curve.  Even-Zohar places random jump, Petaluma, and grid diagrams in a
loose "3D" non-local class of models.

Our Stiefel model is also non-local: consecutive Hamiltonian vertices are
exchangeable projected simplex vertices, not short increments.  This suggests
comparisons with confined/random-jump experiments may be more meaningful than
comparisons with lattice self-avoiding polygons at the same number of segments.

Useful comparison:

- determinant and Alexander polynomial distributions,
- probability of prime vs composite types where classification is possible,
- crossing number in random planar views,
- hyperbolic volume for samples classified as hyperbolic.

### 4. Petaluma Model

The Petaluma model is diagrammatic and non-local.  A petal diagram has one
multi-crossing, and a uniformly random permutation gives the relative heights of
the strands.  Even-Zohar-Hass-Linial-Nowik prove and compute distributions of
finite type invariants in this model:

- *Invariants of Random Knots and Links*:
  <https://arxiv.org/abs/1411.3308>
- *The Distribution of Knots in the Petaluma Model*:
  <https://arxiv.org/abs/1706.06571>

This is not the same as our model, but the comparison may be productive because
both have long, global strands and height data determined by a compact set of
parameters.  Petaluma is permutation-based; Stiefel is continuous and geometric.

Useful comparison:

- finite-type invariants `c_2`, `v_3`,
- determinant/Alexander data,
- hyperbolic volume if we add a reliable volume pipeline,
- growth in crossing number of generic diagrams.

### 5. Random Grid Diagrams

A grid diagram has `n` horizontal and `n` vertical segments, with vertical segments
passing over horizontal segments.  It is encoded by two permutations, and the
random grid model samples those permutations uniformly.  This is not an
independent coin flip at each crossing.

Primary survey reference:

- Even-Zohar, *Models of Random Knots*, Section 3.6:
  <https://arxiv.org/abs/1711.10470>

Difference from our model:

- Grid diagrams prescribe over/under information by axis convention.
- Our planar diagram comes from projecting a 3D polygon to `R^2`; over/under is
  read from the third coordinate.
- If we project only to `R^2` and then assign crossing signs independently, that
  is a new diagrammatic shadow model, not the Stiefel knot law.

Useful comparison:

- crossing-count distribution after reducing/simplifying diagrams,
- finite-type invariant distributions,
- sign correlations at crossings,
- whether the marginal sign balance resembles the grid convention or random
  signs.

### 6. Random Planar Diagrams, Shadows, And Plane Curves

Random diagram models sample a 4-valent planar map or shadow with `n` crossings
and then decorate each crossing with over/under information.  Cantarella-Chapman-
Mastin tabulated diagrams through 10 crossings and computed exact knot
probabilities.  Chapman later proved asymptotic laws for random knot diagrams.

- Jason Cantarella, Harrison Chapman, Matt Mastin,
  *Knot Probabilities in Random Diagrams*:
  <https://arxiv.org/abs/1512.05749>
- Harrison Chapman, *Asymptotic Laws for Random Knot Diagrams*:
  <https://arxiv.org/abs/1608.02638>

This is the right family to compare with if we study the planar shadow of our
Stiefel polygon.  However, the shadow distribution is different:

- Random planar diagrams usually sample a planar map/shadow first.
- Our shadow is the projection of a specific cyclic polygon with vertices coming
  from a random Stiefel frame.
- Random diagram models often assign crossing signs independently.
- Our crossing signs are induced by one global height coordinate and are therefore
  correlated through the same sampled frame.

Useful comparison:

- number of crossings in a generic planar projection,
- shadow graph statistics,
- distribution of crossing signs and pairwise sign correlations,
- knot type conditioned on crossing count,
- frequency-rank behavior of knot types.

### 7. Large Random Knot Projections And Hamiltonian Plane Graphs

Diao-Ernst-Ziegler and related work study random knot projections and large
4-regular Hamiltonian plane graphs.  This is relevant because our cyclic order is
Hamiltonian, but the randomness enters differently.

- Yuanan Diao, Claus Ernst, Uta Ziegler,
  *Generating Large Random Knot Projections*, in *Physical and Numerical Models in
  Knot Theory*, 2005.
- O. Ascigil, Y. Diao, C. Ernst, D. High, U. Ziegler,
  *Generating 4-regular Hamiltonian Plane Graphs*:
  <https://combinatorialpress.com/article/jcmcc/Volume%20073/vol-073-paper%209.pdf>

The point to make in our manuscript is limited: these works randomize planar
projection combinatorics directly, whereas we randomize a Stiefel frame and obtain
both planar shadow and crossing data from geometry.

## The Plane Projection Issue

Given a Stiefel sample in `R^3`, a generic projection to `R^2` produces a knot
diagram:

```text
Q in St_3(1^perp)
    -> spatial polygon Gamma(Q) in R^3
    -> planar shadow pi(Gamma(Q)) in R^2
    -> over/under data from the omitted height coordinate.
```

This is the correct diagrammatic representation of our existing model.

A different model would be:

```text
Q_2 in St_2(1^perp)
    -> planar polygonal shadow in R^2
    -> independent fair coin at each crossing for over/under.
```

That model is closer to random planar diagrams or griddle models.  It is not
equivalent to the 3D Stiefel model unless one proves that the induced crossing-sign
law agrees with independent fair signs, and it almost certainly will not agree
jointly because all crossing signs are functions of the same hidden height vector.

This gives a concrete comparison project:

1. Sample `Q in St_3(1^perp)`.
2. Project to a fixed/generic plane and compute all crossings.
3. Record the sign at each crossing from the third coordinate.
4. Compare with the same shadow decorated by independent fair signs.
5. Compare both to standard random planar diagram statistics.

The statistic to check first is not knot type; it is sign correlation.  If signs are
substantially correlated, then the Stiefel shadow model is visibly different from a
coin-flip diagram model even before knot classification.

## Comparison Plan For The Manuscript

The manuscript should use three levels of comparison.

### Level A: Direct Finite-`N` Probability Tables

Primary object:

```text
p_N(K),  estimated by p_hat_N(K).
```

Use Haar Stiefel samples as the primary ensemble.  Report raw Gaussian only as a
classifier-conditioning check, because exact topology is invariant under the polar
factor but the classifier sees the sheared polygon.

### Level B: Stick-Number Coverage

For each `N`, use KnotInfo/database-knotinfo to list knots with known exact
`stick(K) <= N`.  Then report which types appeared in the Haar samples.  This is
not a proof of absence; it is a coverage diagnostic.

### Level C: Diagrammatic Comparison

For each Stiefel sample, produce a planar diagram and record:

- crossing count before simplification,
- simplified crossing count if the classifier reports one,
- crossing-sign balance,
- crossing-sign correlations,
- determinant/Alexander data,
- knot type when confidently classified.

Then compare against:

- independent-sign decoration of the same Stiefel shadows,
- random grid diagrams,
- random planar diagrams/shadows,
- random jump or confined polygon models at comparable vertex count.

This gives an honest answer to "compare how": compare observables shared by the
models, not the sampling laws themselves.

## Manuscript Edits Suggested By This Study

1. Add a references paragraph in the introduction separating:
   random spatial polygons, random projections, Petaluma/grid models, random
   planar diagrams, and stick-number computations.
2. State explicitly that random grid diagrams use a vertical-over-horizontal rule,
   while random planar diagram models use crossing decorations/signs.
3. Replace any wording suggesting equivalence to random diagram models with:
   "the Stiefel model induces a random diagram by planar projection, but its
   shadow and crossing-sign laws are geometrically constrained."
4. Add a "what was computed" paragraph before the tables: finite samples,
   classifier-dependent labels, Wilson intervals, unknown rates, and seed/backend
   metadata.
5. Keep the support theorem separate from the experiments.  The theorem says
   positive Haar measure exactly matches stick number; the experiments estimate
   visible mass in finite samples.
