import Mathlib.MeasureTheory.Integral.Pi
import Mathlib.Probability.Distributions.Gaussian.Real
import Mathlib.Probability.Moments.SubGaussian
import Mathlib.Probability.Moments.Variance

/-!
# Kernel-checked probability core for the Gaussian polygon model

This file keeps the genuinely analytic or topological inputs of the manuscript out of the
probability bookkeeping.  It defines the finite iid Gaussian coordinate law used before
centering, proves its row/column permutation symmetry, computes the expected squared distance
between two distinct raw columns, and exposes the precise martingale hypotheses under which
mathlib's Azuma--Hoeffding theorem yields the concentration estimate.

The coordinate variance is `1 / 3`, so three independent coordinate differences have total
expected squared length `3 * (2 / 3) = 2`.
-/

open MeasureTheory ProbabilityTheory
open scoped BigOperators ENNReal NNReal

namespace GaussianKnots.GaussianProbabilityCore

/-- A coordinate of a `3 × N` matrix: a spatial row and a vertex column. -/
abbrev Coordinate (N : ℕ) := Fin 3 × Fin N

/-- A `3 × N` real matrix, represented as a function on its finite coordinate set. -/
abbrev RawMatrix (N : ℕ) := Coordinate N → ℝ

/-- The centered one-dimensional Gaussian law of variance `1 / 3`. -/
noncomputable def coordinateGaussian : Measure ℝ :=
  gaussianReal 0 (1 / 3)

noncomputable instance : IsProbabilityMeasure coordinateGaussian := by
  unfold coordinateGaussian
  infer_instance

/-- The law of a raw `3 × N` matrix with iid `N(0, 1/3)` entries. -/
noncomputable def rawGaussianLaw (N : ℕ) : Measure (RawMatrix N) :=
  Measure.pi fun _ : Coordinate N ↦ coordinateGaussian

instance (N : ℕ) : IsProbabilityMeasure (rawGaussianLaw N) := by
  unfold rawGaussianLaw
  infer_instance

/-- Reindex the rows and columns of a raw matrix. -/
noncomputable def reindex (N : ℕ) (σ : Equiv.Perm (Fin 3)) (τ : Equiv.Perm (Fin N)) :
    RawMatrix N ≃ᵐ RawMatrix N :=
  MeasurableEquiv.piCongrLeft (fun _ : Coordinate N ↦ ℝ) (Equiv.prodCongr σ τ)

/-- The iid raw Gaussian law is invariant under arbitrary row and column permutations. -/
theorem rawGaussianLaw_map_reindex (N : ℕ) (σ : Equiv.Perm (Fin 3))
    (τ : Equiv.Perm (Fin N)) :
    (rawGaussianLaw N).map (reindex N σ τ) = rawGaussianLaw N := by
  simpa [rawGaussianLaw, coordinateGaussian, reindex] using
    (MeasureTheory.measurePreserving_piCongrLeft
      (fun _ : Coordinate N ↦ coordinateGaussian) (Equiv.prodCongr σ τ)).map_eq

/-- In particular, the columns of the raw Gaussian matrix are exchangeable. -/
theorem rawGaussianLaw_map_columnPermutation (N : ℕ) (τ : Equiv.Perm (Fin N)) :
    (rawGaussianLaw N).map (reindex N (Equiv.refl (Fin 3)) τ) = rawGaussianLaw N :=
  rawGaussianLaw_map_reindex N (Equiv.refl (Fin 3)) τ

/-- The squared Euclidean distance between two columns of a raw matrix. -/
def columnDistanceSq {N : ℕ} (j k : Fin N) (x : RawMatrix N) : ℝ :=
  ∑ r : Fin 3, (x (r, j) - x (r, k)) ^ 2

theorem measurable_columnDistanceSq {N : ℕ} (j k : Fin N) :
    Measurable (columnDistanceSq j k) := by
  unfold columnDistanceSq
  fun_prop

/-- The average of one spatial row over all `N` columns. -/
noncomputable def rowAverage (N : ℕ) (x : RawMatrix N) (r : Fin 3) : ℝ :=
  (N : ℝ)⁻¹ * ∑ j : Fin N, x (r, j)

/-- Subtract the row-wise column average from every entry. -/
noncomputable def centeredMatrix (N : ℕ) (x : RawMatrix N) : RawMatrix N :=
  fun a ↦ x a - rowAverage N x a.1

theorem measurable_centeredMatrix (N : ℕ) : Measurable (centeredMatrix N) := by
  refine measurable_pi_lambda _ fun a ↦ ?_
  unfold centeredMatrix rowAverage
  fun_prop

/-- The singular Gaussian law obtained by projecting the raw matrix to the row-centered subspace. -/
noncomputable def centeredGaussianLaw (N : ℕ) : Measure (RawMatrix N) :=
  (rawGaussianLaw N).map (centeredMatrix N)

instance (N : ℕ) : IsProbabilityMeasure (centeredGaussianLaw N) := by
  unfold centeredGaussianLaw
  exact Measure.isProbabilityMeasure_map (measurable_centeredMatrix N).aemeasurable

/-- For a nonempty matrix, every spatial row of the centered matrix sums to zero. -/
theorem sum_centeredMatrix_row {N : ℕ} (hN : 0 < N) (x : RawMatrix N) (r : Fin 3) :
    ∑ j : Fin N, centeredMatrix N x (r, j) = 0 := by
  have hN0 : (N : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hN)
  simp only [centeredMatrix, rowAverage, Finset.sum_sub_distrib, Finset.sum_const,
    Finset.card_univ, Fintype.card_fin, nsmul_eq_mul]
  field_simp
  ring

/-- Centering cancels from every coordinate difference. -/
@[simp] theorem centeredMatrix_entry_sub {N : ℕ} (x : RawMatrix N)
    (r : Fin 3) (j k : Fin N) :
    centeredMatrix N x (r, j) - centeredMatrix N x (r, k) =
      x (r, j) - x (r, k) := by
  simp [centeredMatrix]

/-- Consequently, centering leaves every squared column distance unchanged. -/
@[simp] theorem columnDistanceSq_centeredMatrix {N : ℕ} (x : RawMatrix N)
    (j k : Fin N) :
    columnDistanceSq j k (centeredMatrix N x) = columnDistanceSq j k x := by
  simp [columnDistanceSq]

/-- Evaluation at one matrix coordinate preserves the specified Gaussian law. -/
theorem measurePreserving_entry {N : ℕ} (a : Coordinate N) :
    MeasurePreserving (fun x : RawMatrix N ↦ x a) (rawGaussianLaw N) coordinateGaussian := by
  simpa [rawGaussianLaw] using
    (MeasureTheory.measurePreserving_eval (fun _ : Coordinate N ↦ coordinateGaussian) a)

/-- Every raw matrix entry has mean zero. -/
@[simp] theorem integral_entry {N : ℕ} (a : Coordinate N) :
    ∫ x : RawMatrix N, x a ∂rawGaussianLaw N = 0 := by
  calc
    (∫ x : RawMatrix N, x a ∂rawGaussianLaw N) =
        ∫ y : ℝ, y ∂coordinateGaussian := by
      simpa [rawGaussianLaw] using
        (MeasureTheory.integral_comp_eval
          (μ := fun _ : Coordinate N ↦ coordinateGaussian) (i := a)
          (f := id) aestronglyMeasurable_id)
    _ = 0 := by simp [coordinateGaussian]

/-- Every raw matrix entry has variance `1 / 3`. -/
@[simp] theorem variance_entry {N : ℕ} (a : Coordinate N) :
    Var[fun x : RawMatrix N ↦ x a; rawGaussianLaw N] = (1 / 3 : ℝ) := by
  have h := (measurePreserving_entry a).variance_fun_comp
    (f := id) measurable_id.aemeasurable
  simpa [coordinateGaussian] using h

/-- Every raw entry has a finite second moment. -/
theorem memLp_two_entry {N : ℕ} (a : Coordinate N) :
    MemLp (fun x : RawMatrix N ↦ x a) 2 (rawGaussianLaw N) :=
  (memLp_id_gaussianReal (μ := 0) (v := (1 / 3)) 2).comp_measurePreserving
    (measurePreserving_entry a)

/-- Distinct entries of the raw matrix are independent. -/
theorem indepFun_entry {N : ℕ} {a b : Coordinate N} (hab : a ≠ b) :
    IndepFun (fun x : RawMatrix N ↦ x a) (fun x : RawMatrix N ↦ x b)
      (rawGaussianLaw N) := by
  have h : iIndepFun (fun a (x : RawMatrix N) ↦ x a) (rawGaussianLaw N) := by
    simpa [rawGaussianLaw] using
      (iIndepFun_pi (μ := fun _ : Coordinate N ↦ coordinateGaussian)
        (X := fun _ ↦ id) (fun _ ↦ measurable_id.aemeasurable))
  exact h.indepFun hab

/-- One spatial-coordinate difference between distinct columns has variance `2 / 3`. -/
theorem variance_entry_sub {N : ℕ} (r : Fin 3) {j k : Fin N} (hjk : j ≠ k) :
    Var[fun x : RawMatrix N ↦ x (r, j) - x (r, k); rawGaussianLaw N] = (2 / 3 : ℝ) := by
  have hcoord : (r, j) ≠ (r, k) := by
    intro h
    exact hjk (Prod.mk.inj h).2
  rw [variance_fun_sub (memLp_two_entry (r, j)) (memLp_two_entry (r, k)),
    variance_entry, variance_entry,
    (indepFun_entry hcoord).covariance_eq_zero (memLp_two_entry (r, j))
      (memLp_two_entry (r, k))]
  norm_num

/-- One spatial-coordinate squared difference between distinct columns has expectation `2 / 3`. -/
theorem integral_sq_entry_sub {N : ℕ} (r : Fin 3) {j k : Fin N} (hjk : j ≠ k) :
    ∫ x : RawMatrix N, (x (r, j) - x (r, k)) ^ 2 ∂rawGaussianLaw N = (2 / 3 : ℝ) := by
  have hvar := variance_entry_sub r hjk
  rw [variance_eq_integral (by fun_prop)] at hvar
  have hj_int := (memLp_two_entry (r, j)).integrable (by norm_num)
  have hk_int := (memLp_two_entry (r, k)).integrable (by norm_num)
  rw [integral_sub hj_int hk_int, integral_entry, integral_entry, sub_zero] at hvar
  simpa using hvar

/-- Distinct raw Gaussian columns have expected squared distance exactly `2`. -/
theorem integral_columnDistanceSq {N : ℕ} {j k : Fin N} (hjk : j ≠ k) :
    ∫ x : RawMatrix N, columnDistanceSq j k x ∂rawGaussianLaw N = 2 := by
  change (∫ x : RawMatrix N,
    ∑ r : Fin 3, (x (r, j) - x (r, k)) ^ 2 ∂rawGaussianLaw N) = 2
  rw [integral_finsetSum]
  · simp_rw [integral_sq_entry_sub _ hjk]
    norm_num
  · intro r _hr
    exact ((memLp_two_entry (r, j)).sub (memLp_two_entry (r, k))).integrable_sq

/-- The centered Gaussian model has the same expected squared inter-column distance `2`. -/
theorem integral_columnDistanceSq_centered {N : ℕ} {j k : Fin N} (hjk : j ≠ k) :
    ∫ x : RawMatrix N, columnDistanceSq j k x ∂centeredGaussianLaw N = 2 := by
  rw [centeredGaussianLaw, integral_map (measurable_centeredMatrix N).aemeasurable]
  · simpa using integral_columnDistanceSq hjk
  · exact (measurable_columnDistanceSq j k).aestronglyMeasurable

section AzumaHoeffding

variable {Ω : Type*} [mΩ : MeasurableSpace Ω] [StandardBorelSpace Ω]
  {μ : Measure Ω} [IsZeroOrProbabilityMeasure μ]

/--
A typed, kernel-checked concentration interface for the manuscript.

The inputs are not a probability bound: they are an adapted increment process, its initial
sub-Gaussian certificate, and conditional sub-Gaussian certificates for later increments.
Mathlib then proves the exponential upper-tail estimate.  A future bounded-differences layer may
construct these certificates from coordinate-replacement bounds without changing this theorem.
-/
theorem azumaHoeffding_of_increment_certificate
    {Y : ℕ → Ω → ℝ} {cY : ℕ → ℝ≥0}
    {ℱ : Filtration ℕ mΩ} (h_adapted : StronglyAdapted ℱ Y)
    (h0 : HasSubgaussianMGF (Y 0) (cY 0) μ) (n : ℕ)
    (h_subG : ∀ i < n - 1,
      HasCondSubgaussianMGF (ℱ i) (ℱ.le i) (Y (i + 1)) (cY (i + 1)) μ)
    {ε : ℝ} (hε : 0 ≤ ε) :
    μ.real {ω | ε ≤ ∑ i ∈ Finset.range n, Y i ω}
      ≤ Real.exp (-ε ^ 2 / (2 * ∑ i ∈ Finset.range n, cY i)) :=
  measure_sum_ge_le_of_hasCondSubgaussianMGF h_adapted h0 n h_subG hε

/-- The same Azuma--Hoeffding step with a named statistic identified with the increment sum. -/
theorem statistic_upperTail_of_increment_certificate
    {T : Ω → ℝ} {Y : ℕ → Ω → ℝ} {cY : ℕ → ℝ≥0}
    {ℱ : Filtration ℕ mΩ} (hT : T = fun ω ↦ ∑ i ∈ Finset.range n, Y i ω)
    (h_adapted : StronglyAdapted ℱ Y) (h0 : HasSubgaussianMGF (Y 0) (cY 0) μ)
    (h_subG : ∀ i < n - 1,
      HasCondSubgaussianMGF (ℱ i) (ℱ.le i) (Y (i + 1)) (cY (i + 1)) μ)
    {ε : ℝ} (hε : 0 ≤ ε) :
    μ.real {ω | ε ≤ T ω}
      ≤ Real.exp (-ε ^ 2 / (2 * ∑ i ∈ Finset.range n, cY i)) := by
  subst T
  exact azumaHoeffding_of_increment_certificate h_adapted h0 n h_subG hε

end AzumaHoeffding

end GaussianKnots.GaussianProbabilityCore
