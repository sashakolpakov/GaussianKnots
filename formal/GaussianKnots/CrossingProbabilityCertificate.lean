import Mathlib.Combinatorics.SimpleGraph.CycleGraph
import Mathlib.Combinatorics.SimpleGraph.DegreeSum
import Mathlib.Data.Nat.Choose.Cast
import GaussianKnots.GaussianProbabilityCore

/-!
# Conditional certificate for the planar crossing theorem

The edges of the cyclic polygon are the vertices of `cycleGraph N`.  Thus an unordered pair of
non-incident polygon edges is literally an edge of the complement of that graph.  This file uses
that concrete finite type as the index set of crossing indicators, proves its cardinality, and
performs the expectation and concentration bookkeeping.

The geometric four-point probability and the construction of the Doob increment certificates
remain explicit inputs.  The outputs are not probability hypotheses: Lean derives the exact mean
and the printed two-sided exponential tail from those structural inputs.
-/

open MeasureTheory ProbabilityTheory
open scoped BigOperators ENNReal NNReal

namespace GaussianKnots.CrossingProbabilityCertificate

/-- The graph whose vertices are polygon edges and whose edges are non-incident edge pairs. -/
def nonincidentEdgeGraph (N : ℕ) : SimpleGraph (Fin N) :=
  (SimpleGraph.cycleGraph N)ᶜ

/-- An unordered pair of non-incident edges of the cyclic `N`-gon. -/
abbrev CrossingSlot (N : ℕ) := (nonincidentEdgeGraph N).edgeSet

instance nonincidentEdgeGraphDecidableAdj (N : ℕ) :
    DecidableRel (nonincidentEdgeGraph N).Adj := by
  unfold nonincidentEdgeGraph
  infer_instance

instance crossingSlotFintype (N : ℕ) : Fintype (CrossingSlot N) :=
  by
    change Fintype (((⊤ : SimpleGraph (Fin N)) \ SimpleGraph.cycleGraph N).edgeSet)
    infer_instance

/-- A cycle with `n + 3` vertices has exactly `n + 3` edges. -/
theorem card_cycleGraph_edgeFinset (n : ℕ) :
    (SimpleGraph.cycleGraph (n + 3)).edgeFinset.card = n + 3 := by
  have h := SimpleGraph.sum_degrees_eq_twice_card_edges
    (SimpleGraph.cycleGraph (n + 3))
  simp_rw [SimpleGraph.cycleGraph_degree_three_le] at h
  simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at h
  apply Nat.eq_of_mul_eq_mul_left (by norm_num : 0 < 2)
  simpa [mul_comm] using h.symm

/-- The exact number of unordered non-incident edge pairs, first in subtraction form. -/
theorem card_crossingSlot_choose_sub (n : ℕ) :
    Fintype.card (CrossingSlot (n + 3)) = (n + 3).choose 2 - (n + 3) := by
  classical
  rw [← SimpleGraph.edgeFinset_card]
  have hedges : (nonincidentEdgeGraph (n + 3)).edgeFinset =
      (⊤ : SimpleGraph (Fin (n + 3))).edgeFinset \
        (SimpleGraph.cycleGraph (n + 3)).edgeFinset := by
    change ((⊤ : SimpleGraph (Fin (n + 3))) \
      SimpleGraph.cycleGraph (n + 3)).edgeFinset = _
    exact SimpleGraph.edgeFinset_sdiff
  rw [hedges,
    Finset.card_sdiff_of_subset (SimpleGraph.edgeFinset_mono le_top),
    SimpleGraph.card_edgeFinset_top_eq_card_choose_two,
    Fintype.card_fin]
  exact congrArg (fun z ↦ (n + 3).choose 2 - z)
    (card_cycleGraph_edgeFinset n)

/-- Twice the number of crossing slots is `N(N-3)`, with `N = n + 3`. -/
theorem two_mul_card_crossingSlot (n : ℕ) :
    2 * Fintype.card (CrossingSlot (n + 3)) = (n + 3) * n := by
  rw [card_crossingSlot_choose_sub, Nat.mul_sub_left_distrib]
  have hchoose : 2 * (n + 3).choose 2 = (n + 3) * (n + 2) := by
    rw [Nat.choose_two_right]
    have heven : Even ((n + 3) * ((n + 3) - 1)) := Nat.even_mul_pred_self (n + 3)
    rw [Nat.mul_div_cancel' heven.two_dvd]
    congr 1
  rw [hchoose]
  rw [mul_comm 2 (n + 3), ← Nat.mul_sub_left_distrib]
  congr 1

/-- The crossing-slot cardinality in the real-valued form used by the expectation formula. -/
theorem cast_card_crossingSlot (n : ℕ) :
    (Fintype.card (CrossingSlot (n + 3)) : ℝ) = ((n + 3 : ℕ) : ℝ) * n / 2 := by
  have h := congrArg (fun z : ℕ ↦ (z : ℝ)) (two_mul_card_crossingSlot n)
  norm_num at h ⊢
  linarith

/-- Every polygon edge has exactly `n` non-incident partners when `N = n + 3`. -/
theorem degree_nonincidentEdgeGraph (n : ℕ) (e : Fin (n + 3)) :
    (nonincidentEdgeGraph (n + 3)).degree e = n := by
  unfold nonincidentEdgeGraph
  rw [SimpleGraph.degree_compl, SimpleGraph.cycleGraph_degree_three_le]
  simp only [Fintype.card_fin]
  omega

/-- Slots that can involve one of the two polygon edges meeting a specified vertex. -/
def affectedSlotsAtVertex (n : ℕ) (v : Fin (n + 3)) : Finset (Sym2 (Fin (n + 3))) :=
  (nonincidentEdgeGraph (n + 3)).incidenceFinset (v - 1) ∪
    (nonincidentEdgeGraph (n + 3)).incidenceFinset v

/-- The actual affected-slot family has size at most `2(N-3)`. -/
theorem card_affectedSlotsAtVertex_le (n : ℕ) (v : Fin (n + 3)) :
    (affectedSlotsAtVertex n v).card ≤ 2 * n := by
  unfold affectedSlotsAtVertex
  calc
    ((nonincidentEdgeGraph (n + 3)).incidenceFinset (v - 1) ∪
      (nonincidentEdgeGraph (n + 3)).incidenceFinset v).card ≤
        ((nonincidentEdgeGraph (n + 3)).incidenceFinset (v - 1)).card +
          ((nonincidentEdgeGraph (n + 3)).incidenceFinset v).card :=
      Finset.card_union_le _ _
    _ = 2 * n := by
      rw [SimpleGraph.card_incidenceFinset_eq_degree,
        SimpleGraph.card_incidenceFinset_eq_degree,
        degree_nonincidentEdgeGraph, degree_nonincidentEdgeGraph]
      omega

/-- A measurable family of crossing events indexed by actual non-incident edge pairs. -/
structure CrossingProcess {Ω : Type*} [MeasurableSpace Ω] (N : ℕ) where
  crosses : CrossingSlot N → Ω → Prop
  measurable_crosses : ∀ s, MeasurableSet {ω | crosses s ω}

/-- The `0/1` indicator of one crossing slot. -/
noncomputable def CrossingProcess.indicator {Ω : Type*} [MeasurableSpace Ω] {N : ℕ}
    (X : CrossingProcess (Ω := Ω) N) (s : CrossingSlot N) : Ω → ℝ :=
  {ω | X.crosses s ω}.indicator (fun _ ↦ 1)

/-- The total number of crossing slots that occur. -/
noncomputable def CrossingProcess.statistic {Ω : Type*} [MeasurableSpace Ω] {N : ℕ}
    (X : CrossingProcess (Ω := Ω) N) : Ω → ℝ :=
  fun ω ↦ ∑ s : CrossingSlot N, X.indicator s ω

theorem CrossingProcess.measurable_indicator {Ω : Type*} [MeasurableSpace Ω] {N : ℕ}
    (X : CrossingProcess (Ω := Ω) N) (s : CrossingSlot N) : Measurable (X.indicator s) := by
  exact measurable_const.indicator (X.measurable_crosses s)

theorem CrossingProcess.measurable_statistic {Ω : Type*} [MeasurableSpace Ω] {N : ℕ}
    (X : CrossingProcess (Ω := Ω) N) : Measurable X.statistic := by
  unfold CrossingProcess.statistic
  exact Finset.measurable_fun_sum _ fun s _hs ↦ X.measurable_indicator s

section ActualGaussianCrossings

open GaussianKnots.GaussianProbabilityCore

/-- A point of the planar shadow, read from the first two rows of the raw matrix. -/
def rawPlanarPoint {N : ℕ} (x : RawMatrix N) (j : Fin N) : ℝ × ℝ :=
  (x (0, j), x (1, j))

/-- Translate a planar point by subtracting a common vector. -/
def subtractPoint (a t : ℝ × ℝ) : ℝ × ℝ :=
  (a.1 - t.1, a.2 - t.2)

/-- Row centering is one common translation of all planar vertices. -/
theorem rawPlanarPoint_centeredMatrix {N : ℕ} (x : RawMatrix N) (j : Fin N) :
    rawPlanarPoint (centeredMatrix N x) j =
      subtractPoint (rawPlanarPoint x j) (rowAverage N x 0, rowAverage N x 1) := by
  rfl

/-- Replace one three-dimensional matrix column, leaving every other vertex fixed. -/
def replaceColumn {N : ℕ} (x : RawMatrix N) (v : Fin N) (p : Fin 3 → ℝ) : RawMatrix N :=
  fun a ↦ if a.2 = v then p a.1 else x a

@[simp] theorem rawPlanarPoint_replaceColumn_of_ne {N : ℕ} (x : RawMatrix N)
    (v j : Fin N) (p : Fin 3 → ℝ) (hj : j ≠ v) :
    rawPlanarPoint (replaceColumn x v p) j = rawPlanarPoint x j := by
  simp [rawPlanarPoint, replaceColumn, hj]

/-- The signed planar area determinant. -/
def orient2 (a b c : ℝ × ℝ) : ℝ :=
  (b.1 - a.1) * (c.2 - a.2) - (b.2 - a.2) * (c.1 - a.1)

/-- The planar orientation determinant is invariant under common translation. -/
@[simp] theorem orient2_subtractPoint (a b c t : ℝ × ℝ) :
    orient2 (subtractPoint a t) (subtractPoint b t) (subtractPoint c t) = orient2 a b c := by
  unfold orient2 subtractPoint
  ring

/-- Strict segment crossing, excluding endpoints and collinear degeneracies. -/
def properSegmentsCross (a b c d : ℝ × ℝ) : Prop :=
  orient2 a b c * orient2 a b d < 0 ∧ orient2 c d a * orient2 c d b < 0

theorem properSegmentsCross_swap (a b c d : ℝ × ℝ) :
    properSegmentsCross a b c d = properSegmentsCross c d a b := by
  simp only [properSegmentsCross, and_comm]

/-- Proper segment crossing is invariant under common translation. -/
@[simp] theorem properSegmentsCross_subtractPoint (a b c d t : ℝ × ℝ) :
    properSegmentsCross (subtractPoint a t) (subtractPoint b t)
      (subtractPoint c t) (subtractPoint d t) = properSegmentsCross a b c d := by
  simp [properSegmentsCross]

/-- The crossing relation on two cyclic edge indices in an actual raw Gaussian matrix. -/
def rawSlotCrossesRel (n : ℕ) (x : RawMatrix (n + 3)) (i j : Fin (n + 3)) : Prop :=
  properSegmentsCross (rawPlanarPoint x i) (rawPlanarPoint x (i + 1))
    (rawPlanarPoint x j) (rawPlanarPoint x (j + 1))

theorem rawSlotCrossesRel_symm (n : ℕ) (x : RawMatrix (n + 3)) :
    Std.Symm (rawSlotCrossesRel n x) := by
  constructor
  intro i j
  unfold rawSlotCrossesRel
  rw [properSegmentsCross_swap]
  exact id

/-- Centering the Gaussian matrix does not change any concrete crossing relation. -/
@[simp] theorem rawSlotCrossesRel_centeredMatrix (n : ℕ) (x : RawMatrix (n + 3))
    (i j : Fin (n + 3)) :
    rawSlotCrossesRel n (centeredMatrix (n + 3) x) i j = rawSlotCrossesRel n x i j := by
  simp only [rawSlotCrossesRel, rawPlanarPoint_centeredMatrix,
    properSegmentsCross_subtractPoint]

private theorem rawSym2Crosses_replaceColumn_of_avoids (n : ℕ)
    (z : Sym2 (Fin (n + 3))) (x : RawMatrix (n + 3)) (v : Fin (n + 3))
    (p : Fin 3 → ℝ)
    (havoid : ∀ e ∈ z, e ≠ v ∧ e + 1 ≠ v) :
    (z ∈ Sym2.fromRel (rawSlotCrossesRel_symm n (replaceColumn x v p))) =
      (z ∈ Sym2.fromRel (rawSlotCrossesRel_symm n x)) := by
  induction z using Sym2.inductionOn with
  | _ i j =>
      simp only [Sym2.fromRel_prop]
      have hi := havoid i (Sym2.mem_mk_left i j)
      have hj := havoid j (Sym2.mem_mk_right i j)
      unfold rawSlotCrossesRel
      rw [rawPlanarPoint_replaceColumn_of_ne _ _ _ _ hi.1,
        rawPlanarPoint_replaceColumn_of_ne _ _ _ _ hi.2,
        rawPlanarPoint_replaceColumn_of_ne _ _ _ _ hj.1,
        rawPlanarPoint_replaceColumn_of_ne _ _ _ _ hj.2]

/-- The measurable crossing event belonging to an actual unordered crossing slot. -/
def rawSlotCrosses (n : ℕ) (s : CrossingSlot (n + 3))
    (x : RawMatrix (n + 3)) : Prop :=
  s.1 ∈ Sym2.fromRel (rawSlotCrossesRel_symm n x)

theorem measurableSet_rawSlotCrossesRel (n : ℕ) (i j : Fin (n + 3)) :
    MeasurableSet {x : RawMatrix (n + 3) | rawSlotCrossesRel n x i j} := by
  unfold rawSlotCrossesRel properSegmentsCross orient2 rawPlanarPoint
  measurability

private theorem measurableSet_rawSym2Crosses (n : ℕ) (z : Sym2 (Fin (n + 3))) :
    MeasurableSet {x : RawMatrix (n + 3) |
      z ∈ Sym2.fromRel (rawSlotCrossesRel_symm n x)} := by
  induction z using Sym2.inductionOn with
  | _ i j =>
      simpa only [Sym2.fromRel_prop] using measurableSet_rawSlotCrossesRel n i j

theorem measurableSet_rawSlotCrosses (n : ℕ) (s : CrossingSlot (n + 3)) :
    MeasurableSet {x : RawMatrix (n + 3) | rawSlotCrosses n s x} := by
  simpa only [rawSlotCrosses] using measurableSet_rawSym2Crosses n s.1

/-- Row centering leaves each actual unordered crossing event unchanged. -/
@[simp] theorem rawSlotCrosses_centeredMatrix (n : ℕ) (s : CrossingSlot (n + 3))
    (x : RawMatrix (n + 3)) :
    rawSlotCrosses n s (centeredMatrix (n + 3) x) = rawSlotCrosses n s x := by
  unfold rawSlotCrosses
  induction s.1 using Sym2.inductionOn with
  | _ i j =>
      simp only [Sym2.fromRel_prop, rawSlotCrossesRel_centeredMatrix]

/-- Replacing one vertex leaves every crossing slot outside the actual affected family unchanged. -/
theorem rawSlotCrosses_replaceColumn_of_not_mem_affected (n : ℕ)
    (s : CrossingSlot (n + 3)) (x : RawMatrix (n + 3)) (v : Fin (n + 3))
    (p : Fin 3 → ℝ) (hs : s.1 ∉ affectedSlotsAtVertex n v) :
    rawSlotCrosses n s (replaceColumn x v p) = rawSlotCrosses n s x := by
  let G := nonincidentEdgeGraph (n + 3)
  have hvnot : v ∉ s.1 := by
    intro hv
    apply hs
    rw [affectedSlotsAtVertex, Finset.mem_union]
    right
    apply (G.mem_incidenceFinset v s.1).2
    exact (G.edge_mem_incidenceSet_iff (e := s)).2 hv
  have hpnot : v - 1 ∉ s.1 := by
    intro hp
    apply hs
    rw [affectedSlotsAtVertex, Finset.mem_union]
    left
    apply (G.mem_incidenceFinset (v - 1) s.1).2
    exact (G.edge_mem_incidenceSet_iff (e := s)).2 hp
  have hav : ∀ e ∈ s.1, e ≠ v ∧ e + 1 ≠ v := by
    intro e he
    constructor
    · intro hev
      exact hvnot (hev ▸ he)
    · intro hev
      have hep : e = v - 1 := by
        calc
          e = (e + 1) - 1 := by simp
          _ = v - 1 := congrArg (fun z ↦ z - 1) hev
      exact hpnot (hep ▸ he)
  exact rawSym2Crosses_replaceColumn_of_avoids n s.1 x v p hav

/-- The actual proper-crossing process on the planar raw Gaussian polygon. -/
noncomputable def rawGaussianCrossingProcess (n : ℕ) :
    CrossingProcess (Ω := RawMatrix (n + 3)) (n + 3) where
  crosses := rawSlotCrosses n
  measurable_crosses := measurableSet_rawSlotCrosses n

/-- Its statistic is the manuscript's number of proper non-incident crossings. -/
noncomputable abbrev rawGaussianCrossingStatistic (n : ℕ) : RawMatrix (n + 3) → ℝ :=
  (rawGaussianCrossingProcess n).statistic

/-- The affected family as a finset of typed crossing slots. -/
def affectedCrossingSlots (n : ℕ) (v : Fin (n + 3)) :
    Finset (CrossingSlot (n + 3)) :=
  Finset.univ.filter fun s ↦ s.1 ∈ affectedSlotsAtVertex n v

theorem card_affectedCrossingSlots_le (n : ℕ) (v : Fin (n + 3)) :
    (affectedCrossingSlots n v).card ≤ (affectedSlotsAtVertex n v).card := by
  classical
  apply Finset.card_le_card_of_injOn (fun s : CrossingSlot (n + 3) ↦ s.1)
  · intro s hs
    simpa [affectedCrossingSlots] using hs
  · intro a _ha b _hb hab
    exact Subtype.ext hab

private theorem abs_rawCrossingIndicator_sub_le_one (n : ℕ)
    (s : CrossingSlot (n + 3)) (x y : RawMatrix (n + 3)) :
    |(rawGaussianCrossingProcess n).indicator s x -
      (rawGaussianCrossingProcess n).indicator s y| ≤ 1 := by
  classical
  unfold CrossingProcess.indicator rawGaussianCrossingProcess
  simp only [Set.indicator_apply, Set.mem_setOf_eq]
  split_ifs <;> norm_num

theorem rawCrossingIndicator_replaceColumn_of_not_mem_affected (n : ℕ)
    (s : CrossingSlot (n + 3)) (x : RawMatrix (n + 3)) (v : Fin (n + 3))
    (p : Fin 3 → ℝ) (hs : s ∉ affectedCrossingSlots n v) :
    (rawGaussianCrossingProcess n).indicator s (replaceColumn x v p) =
      (rawGaussianCrossingProcess n).indicator s x := by
  classical
  have hs' : s.1 ∉ affectedSlotsAtVertex n v := by
    simpa [affectedCrossingSlots] using hs
  unfold CrossingProcess.indicator rawGaussianCrossingProcess
  simp only [Set.indicator_apply, Set.mem_setOf_eq,
    rawSlotCrosses_replaceColumn_of_not_mem_affected n s x v p hs']

/-- The complete deterministic bounded-replacement estimate used by McDiarmid. -/
theorem abs_rawGaussianCrossingStatistic_replaceColumn_le (n : ℕ)
    (x : RawMatrix (n + 3)) (v : Fin (n + 3)) (p : Fin 3 → ℝ) :
    |rawGaussianCrossingStatistic n (replaceColumn x v p) -
      rawGaussianCrossingStatistic n x| ≤ (2 * n : ℕ) := by
  classical
  let A := affectedCrossingSlots n v
  let f : CrossingSlot (n + 3) → ℝ := fun s ↦
    (rawGaussianCrossingProcess n).indicator s (replaceColumn x v p) -
      (rawGaussianCrossingProcess n).indicator s x
  have hsum :
      rawGaussianCrossingStatistic n (replaceColumn x v p) -
          rawGaussianCrossingStatistic n x = ∑ s ∈ A, f s := by
    unfold rawGaussianCrossingStatistic CrossingProcess.statistic
    rw [← Finset.sum_sub_distrib]
    symm
    apply Finset.sum_subset (Finset.subset_univ A)
    intro s _hsu hsA
    exact sub_eq_zero.mpr
      (rawCrossingIndicator_replaceColumn_of_not_mem_affected n s x v p hsA)
  rw [hsum]
  calc
    |∑ s ∈ A, f s| ≤ ∑ s ∈ A, |f s| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _s ∈ A, (1 : ℝ) := by
      apply Finset.sum_le_sum
      intro s hs
      exact abs_rawCrossingIndicator_sub_le_one n s _ _
    _ = (A.card : ℝ) := by simp
    _ ≤ ((affectedSlotsAtVertex n v).card : ℝ) := by
      exact_mod_cast card_affectedCrossingSlots_le n v
    _ ≤ (2 * n : ℕ) := by
      exact_mod_cast card_affectedSlotsAtVertex_le n v

/-- Row centering leaves the complete proper-crossing statistic unchanged pointwise. -/
@[simp] theorem rawGaussianCrossingStatistic_centeredMatrix (n : ℕ)
    (x : RawMatrix (n + 3)) :
    rawGaussianCrossingStatistic n (centeredMatrix (n + 3) x) =
      rawGaussianCrossingStatistic n x := by
  classical
  unfold rawGaussianCrossingStatistic CrossingProcess.statistic CrossingProcess.indicator
  apply Finset.sum_congr rfl
  intro s _hs
  simp only [rawGaussianCrossingProcess, Set.indicator_apply, Set.mem_setOf_eq,
    rawSlotCrosses_centeredMatrix]

end ActualGaussianCrossings

section Mean

variable {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω) [IsFiniteMeasure μ]

omit [IsFiniteMeasure μ] in
theorem CrossingProcess.integral_indicator {N : ℕ} (X : CrossingProcess (Ω := Ω) N)
    (s : CrossingSlot N) :
    ∫ ω, X.indicator s ω ∂μ = μ.real {ω | X.crosses s ω} := by
  exact integral_indicator_one (X.measurable_crosses s)

/-- Linearity of expectation for the concrete finite crossing statistic. -/
theorem CrossingProcess.integral_statistic {N : ℕ} (X : CrossingProcess (Ω := Ω) N) :
    ∫ ω, X.statistic ω ∂μ = ∑ s : CrossingSlot N, μ.real {ω | X.crosses s ω} := by
  change (∫ ω, ∑ s : CrossingSlot N, X.indicator s ω ∂μ) = _
  rw [integral_finsetSum]
  · simp_rw [X.integral_indicator μ]
  · intro s _hs
    exact Integrable.indicator (integrable_const (1 : ℝ)) (X.measurable_crosses s)

/-- If every non-incident edge pair crosses with probability `p`, then the mean is `L_N p`. -/
theorem CrossingProcess.integral_statistic_of_constant_slot_probability {N : ℕ}
    (X : CrossingProcess (Ω := Ω) N) (p : ℝ)
    (hslot : ∀ s : CrossingSlot N, μ.real {ω | X.crosses s ω} = p) :
    ∫ ω, X.statistic ω ∂μ = Fintype.card (CrossingSlot N) * p := by
  rw [X.integral_statistic μ]
  simp [hslot]

/-- The exact manuscript mean for `N = n + 3`, conditional only on the four-point probability. -/
theorem CrossingProcess.exact_mean {n : ℕ}
    (X : CrossingProcess (Ω := Ω) (n + 3)) (p : ℝ)
    (hslot : ∀ s : CrossingSlot (n + 3), μ.real {ω | X.crosses s ω} = p) :
    ∫ ω, X.statistic ω ∂μ = (((n + 3 : ℕ) : ℝ) * n / 2) * p := by
  rw [X.integral_statistic_of_constant_slot_probability μ p hslot,
    cast_card_crossingSlot]

/-- The frequently used `p = 1/3` specialization gives `N(N-3)/6`. -/
theorem CrossingProcess.exact_mean_one_third {n : ℕ}
    (X : CrossingProcess (Ω := Ω) (n + 3))
    (hslot : ∀ s : CrossingSlot (n + 3), μ.real {ω | X.crosses s ω} = (1 / 3 : ℝ)) :
    ∫ ω, X.statistic ω ∂μ = ((n + 3 : ℕ) : ℝ) * n / 6 := by
  rw [X.exact_mean μ (1 / 3 : ℝ) hslot]
  ring

/-- The Gaussian four-point probability gives the exact mean printed in the manuscript. -/
theorem CrossingProcess.exact_gaussian_mean {n : ℕ}
    (X : CrossingProcess (Ω := Ω) (n + 3))
    (hslot : ∀ s : CrossingSlot (n + 3),
      μ.real {ω | X.crosses s ω} = 2 / Real.pi * Real.arcsin (1 / 3)) :
    ∫ ω, X.statistic ω ∂μ =
      (((n + 3 : ℕ) : ℝ) * n / Real.pi) * Real.arcsin (1 / 3) := by
  rw [X.exact_mean μ (2 / Real.pi * Real.arcsin (1 / 3)) hslot]
  field_simp [Real.pi_ne_zero]

/--
The exact mean for the actual raw planar Gaussian crossing statistic.  Its only remaining input is
the geometric four-point probability, now stated for each concrete non-incident edge-pair event.
-/
theorem rawGaussianCrossingStatistic_exact_mean (n : ℕ)
    (hfourPoint : ∀ s : CrossingSlot (n + 3),
      (GaussianProbabilityCore.rawGaussianLaw (n + 3)).real {x | rawSlotCrosses n s x} =
        2 / Real.pi * Real.arcsin (1 / 3)) :
    ∫ x, rawGaussianCrossingStatistic n x
      ∂GaussianProbabilityCore.rawGaussianLaw (n + 3) =
      (((n + 3 : ℕ) : ℝ) * n / Real.pi) * Real.arcsin (1 / 3) :=
  (rawGaussianCrossingProcess n).exact_gaussian_mean
    (GaussianProbabilityCore.rawGaussianLaw (n + 3)) hfourPoint

/-- The crossing statistic has exactly the same expectation before and after row centering. -/
theorem integral_rawGaussianCrossingStatistic_centered_eq_raw (n : ℕ) :
    ∫ x, rawGaussianCrossingStatistic n x
      ∂GaussianProbabilityCore.centeredGaussianLaw (n + 3) =
    ∫ x, rawGaussianCrossingStatistic n x
      ∂GaussianProbabilityCore.rawGaussianLaw (n + 3) := by
  rw [GaussianProbabilityCore.centeredGaussianLaw,
    integral_map (GaussianProbabilityCore.measurable_centeredMatrix (n + 3)).aemeasurable]
  · simp only [rawGaussianCrossingStatistic_centeredMatrix]
  · exact (rawGaussianCrossingProcess n).measurable_statistic.aestronglyMeasurable

/-- The same exact mean for the row-centered Gaussian polygon used downstream. -/
theorem centeredGaussianCrossingStatistic_exact_mean (n : ℕ)
    (hfourPoint : ∀ s : CrossingSlot (n + 3),
      (GaussianProbabilityCore.rawGaussianLaw (n + 3)).real {x | rawSlotCrosses n s x} =
        2 / Real.pi * Real.arcsin (1 / 3)) :
    ∫ x, rawGaussianCrossingStatistic n x
      ∂GaussianProbabilityCore.centeredGaussianLaw (n + 3) =
      (((n + 3 : ℕ) : ℝ) * n / Real.pi) * Real.arcsin (1 / 3) := by
  rw [integral_rawGaussianCrossingStatistic_centered_eq_raw]
  exact rawGaussianCrossingStatistic_exact_mean n hfourPoint

end Mean

section Concentration

variable {Ω : Type*} [MeasurableSpace Ω] [StandardBorelSpace Ω]
  {μ : Measure Ω} [IsProbabilityMeasure μ]

omit [StandardBorelSpace Ω] in
/-- A reusable two-sided Chernoff bound, obtained from mathlib's one-sided theorem and negation. -/
theorem measure_abs_ge_le_of_hasSubgaussianMGF {T : Ω → ℝ} {c : ℝ≥0}
    (hT : HasSubgaussianMGF T c μ) {t : ℝ} (ht : 0 ≤ t) :
    μ.real {ω | t ≤ |T ω|} ≤ 2 * Real.exp (-t ^ 2 / (2 * c)) := by
  let upper : Set Ω := {ω | t ≤ T ω}
  let lower : Set Ω := {ω | t ≤ -T ω}
  calc
    μ.real {ω | t ≤ |T ω|} ≤ μ.real (upper ∪ lower) := by
      refine measureReal_mono ?_ (measure_ne_top _ _)
      intro ω hω
      change t ≤ |T ω| at hω
      rw [Set.mem_union, Set.mem_setOf_eq, Set.mem_setOf_eq]
      exact (le_abs.mp hω).imp_right (by simp)
    _ ≤ μ.real upper + μ.real lower := measureReal_union_le _ _
    _ ≤ Real.exp (-t ^ 2 / (2 * c)) + Real.exp (-t ^ 2 / (2 * c)) := by
      apply add_le_add
      · exact hT.measure_ge_le ht
      · simpa [lower] using hT.neg.measure_ge_le ht
    _ = 2 * Real.exp (-t ^ 2 / (2 * c)) := by ring

/-- The sub-Gaussian proxy contributed by exposing one vertex. -/
def crossingIncrementProxy (N : ℕ) : ℝ≥0 :=
  (((N - 3 : ℕ) : ℝ≥0)) ^ 2

/--
The structural Doob-martingale input left by the bounded-replacement step.

This is substantially stronger than assuming a tail probability: it supplies an actual filtration,
an adapted increment process, an identity with the centered crossing statistic, and the conditional
MGF estimates with the manuscript's per-vertex proxy `(N-3)^2`.
-/
structure CrossingDoobCertificate (N : ℕ)
    (X : CrossingProcess (Ω := Ω) N) where
  filtration : Filtration ℕ ‹MeasurableSpace Ω›
  increments : ℕ → Ω → ℝ
  centered_eq_sum :
    (fun ω ↦ X.statistic ω - ∫ ω', X.statistic ω' ∂μ) =
      fun ω ↦ ∑ i ∈ Finset.range N, increments i ω
  stronglyAdapted : StronglyAdapted filtration increments
  initialSubgaussian : HasSubgaussianMGF (increments 0) (crossingIncrementProxy N) μ
  laterCondSubgaussian : ∀ i < N - 1,
    HasCondSubgaussianMGF (filtration i) (filtration.le i) (increments (i + 1))
      (crossingIncrementProxy N) μ

/-- Mathlib sums the certified Doob increments and obtains the exact total variance proxy. -/
theorem CrossingDoobCertificate.centered_hasSubgaussianMGF {N : ℕ}
    {X : CrossingProcess (Ω := Ω) N} (cert : CrossingDoobCertificate (μ := μ) N X) :
    HasSubgaussianMGF
      (fun ω ↦ X.statistic ω - ∫ ω', X.statistic ω' ∂μ)
      ((N : ℝ≥0) * crossingIncrementProxy N) μ := by
  rw [cert.centered_eq_sum]
  simpa [crossingIncrementProxy] using
    (HasSubgaussianMGF.sum_of_hasCondSubgaussianMGF
      (μ := μ) (Y := cert.increments)
      (cY := fun _ ↦ crossingIncrementProxy N) (ℱ := cert.filtration)
      cert.stronglyAdapted cert.initialSubgaussian N cert.laterCondSubgaussian)

/-- The manuscript's two-sided McDiarmid/Azuma bound, derived from the Doob certificate. -/
theorem CrossingDoobCertificate.two_sided_tail {N : ℕ}
    {X : CrossingProcess (Ω := Ω) N} (cert : CrossingDoobCertificate (μ := μ) N X)
    {t : ℝ} (ht : 0 ≤ t) :
    μ.real {ω |
      t ≤ |X.statistic ω - ∫ ω', X.statistic ω' ∂μ|} ≤
      2 * Real.exp
        (-t ^ 2 / (2 * (N : ℝ) *
          ((↑((N : ℝ≥0) - 3) : ℝ) ^ 2))) := by
  simpa [crossingIncrementProxy, mul_assoc] using
    measure_abs_ge_le_of_hasSubgaussianMGF cert.centered_hasSubgaussianMGF ht

/-- For `N ≥ 3`, the previous bound is exactly the printed real-arithmetic formula. -/
theorem CrossingDoobCertificate.two_sided_tail_printed {N : ℕ}
    {X : CrossingProcess (Ω := Ω) N} (cert : CrossingDoobCertificate (μ := μ) N X)
    (hN : 3 ≤ N) {t : ℝ} (ht : 0 ≤ t) :
    μ.real {ω |
      t ≤ |X.statistic ω - ∫ ω', X.statistic ω' ∂μ|} ≤
      2 * Real.exp
        (-t ^ 2 / (2 * (N : ℝ) * (((N : ℝ) - 3) ^ 2))) := by
  have hNnn : (3 : ℝ≥0) ≤ (N : ℝ≥0) := by exact_mod_cast hN
  simpa [NNReal.coe_sub hNnn] using cert.two_sided_tail ht

/-- The algebraically equivalent normalized tail bound printed after the main inequality. -/
theorem CrossingDoobCertificate.normalized_two_sided_tail {N : ℕ}
    {X : CrossingProcess (Ω := Ω) N} (cert : CrossingDoobCertificate (μ := μ) N X)
    (hN : 3 < N) (p : ℝ)
    (hmean : ∫ ω, X.statistic ω ∂μ = (N : ℝ) * ((N : ℝ) - 3) / 2 * p)
    {ε : ℝ} (hε : 0 ≤ ε) :
    μ.real {ω |
      ε ≤ |X.statistic ω / ((N : ℝ) * ((N : ℝ) - 3) / 2) - p|} ≤
      2 * Real.exp (-ε ^ 2 * (N : ℝ) / 8) := by
  let L : ℝ := (N : ℝ) * ((N : ℝ) - 3) / 2
  have hNR : (0 : ℝ) < N := by exact_mod_cast (by omega : 0 < N)
  have hNcast : (3 : ℝ) < (N : ℝ) := by exact_mod_cast hN
  have hNm3 : (0 : ℝ) < (N : ℝ) - 3 := by linarith
  have hL : 0 < L := by
    dsimp [L]
    positivity
  have hset :
      {ω | ε ≤ |X.statistic ω / L - p|} =
        {ω | ε * L ≤ |X.statistic ω - ∫ ω', X.statistic ω' ∂μ|} := by
    ext ω
    simp only [Set.mem_setOf_eq]
    rw [hmean]
    change ε ≤ |X.statistic ω / L - p| ↔
      ε * L ≤ |X.statistic ω - L * p|
    have hid : X.statistic ω / L - p = (X.statistic ω - L * p) / L := by
      field_simp [hL.ne']
    rw [hid, abs_div, abs_of_pos hL, le_div_iff₀ hL]
  change μ.real {ω | ε ≤ |X.statistic ω / L - p|} ≤ _
  rw [hset]
  calc
    μ.real {ω | ε * L ≤ |X.statistic ω - ∫ ω', X.statistic ω' ∂μ|} ≤
        2 * Real.exp
          (-(ε * L) ^ 2 / (2 * (N : ℝ) * (((N : ℝ) - 3) ^ 2))) :=
      cert.two_sided_tail_printed hN.le (mul_nonneg hε hL.le)
    _ = 2 * Real.exp (-ε ^ 2 * (N : ℝ) / 8) := by
      congr 2
      dsimp [L]
      field_simp [hNR.ne', hNm3.ne']
      ring

end Concentration

section RawGaussianConcentration

open GaussianKnots.GaussianProbabilityCore

/-- The printed unnormalized tail for the actual raw Gaussian crossing statistic. -/
theorem rawGaussianCrossingStatistic_two_sided_tail (n : ℕ)
    (cert : CrossingDoobCertificate
      (μ := rawGaussianLaw (n + 3)) (n + 3) (rawGaussianCrossingProcess n))
    {t : ℝ} (ht : 0 ≤ t) :
    (rawGaussianLaw (n + 3)).real {x |
      t ≤ |rawGaussianCrossingStatistic n x -
        ∫ x', rawGaussianCrossingStatistic n x' ∂rawGaussianLaw (n + 3)|} ≤
      2 * Real.exp (-t ^ 2 / (2 * ((n + 3 : ℕ) : ℝ) * (n : ℝ) ^ 2)) := by
  simpa only [Nat.cast_add, Nat.cast_ofNat, add_sub_cancel_right] using
    cert.two_sided_tail_printed (by omega) ht

/--
The normalized tail in the manuscript, now for the concrete raw Gaussian statistic.  The only
external inputs are the four-point crossing probability and the structural Doob certificate.
-/
theorem rawGaussianCrossingStatistic_normalized_tail (n : ℕ)
    (cert : CrossingDoobCertificate
      (μ := rawGaussianLaw (n + 3)) (n + 3) (rawGaussianCrossingProcess n))
    (hn : 0 < n)
    (hfourPoint : ∀ s : CrossingSlot (n + 3),
      (rawGaussianLaw (n + 3)).real {x | rawSlotCrosses n s x} =
        2 / Real.pi * Real.arcsin (1 / 3))
    {ε : ℝ} (hε : 0 ≤ ε) :
    (rawGaussianLaw (n + 3)).real {x |
      ε ≤ |rawGaussianCrossingStatistic n x /
          (((n + 3 : ℕ) : ℝ) * n / 2) -
        (2 / Real.pi * Real.arcsin (1 / 3))|} ≤
      2 * Real.exp (-ε ^ 2 * ((n + 3 : ℕ) : ℝ) / 8) := by
  have hmean :
      ∫ x, rawGaussianCrossingStatistic n x ∂rawGaussianLaw (n + 3) =
        (((n + 3 : ℕ) : ℝ) * (((n + 3 : ℕ) : ℝ) - 3) / 2) *
          (2 / Real.pi * Real.arcsin (1 / 3)) := by
    rw [rawGaussianCrossingStatistic_exact_mean n hfourPoint]
    have hn : (((n + 3 : ℕ) : ℝ) - 3) = (n : ℝ) := by
      push_cast
      ring
    rw [hn]
    field_simp [Real.pi_ne_zero]
  simpa only [Nat.cast_add, Nat.cast_ofNat, add_sub_cancel_right] using
    cert.normalized_two_sided_tail (show (3 : ℕ) < n + 3 by omega)
      (2 / Real.pi * Real.arcsin (1 / 3)) hmean hε

end RawGaussianConcentration

end GaussianKnots.CrossingProbabilityCertificate
