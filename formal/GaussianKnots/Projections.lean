import GaussianKnots.Combinatorics
import GaussianKnots.Geometry

/-!
# Algebraic core for random polygon projections

This module checks identities used in `knot_projections.tex`: the all-pairs
RMS normalization, the signed `4 × 5` cofactor circuit, reflection symmetry,
and the executable Calvo intersection-number decision layer for generic
hexagons.
-/

namespace GaussianKnots

open scoped BigOperators

def normSq3 (u : Vec3) : ℝ := u 0 ^ 2 + u 1 ^ 2 + u 2 ^ 2

private theorem scalar_pairwise_sq_identity {N : ℕ} (x : Fin N → ℝ) :
    ∑ i, ∑ j, (x i - x j) ^ 2 =
      2 * N * ∑ i, (x i) ^ 2 - 2 * (∑ i, x i) ^ 2 := by
  have hmixed : (∑ i, ∑ j, x i * x j) = (∑ i, x i) ^ 2 := by
    calc
      (∑ i, ∑ j, x i * x j) = ∑ i, x i * (∑ j, x j) := by
        apply Finset.sum_congr rfl
        intro i _hi
        rw [Finset.mul_sum]
      _ = (∑ i, x i) * (∑ j, x j) := by rw [Finset.sum_mul]
      _ = (∑ i, x i) ^ 2 := by ring
  have hmixed2 : (∑ i, ∑ j, 2 * x i * x j) = 2 * (∑ i, x i) ^ 2 := by
    calc
      (∑ i, ∑ j, 2 * x i * x j) = 2 * (∑ i, ∑ j, x i * x j) := by
        rw [Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro i _hi
        rw [Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro j _hj
        ring
      _ = 2 * (∑ i, x i) ^ 2 := by rw [hmixed]
  have hsquares :
      (∑ i, (N : ℝ) * x i ^ 2) + (∑ i, (N : ℝ) * x i ^ 2) =
        ∑ i, 2 * (N : ℝ) * x i ^ 2 := by
    rw [← Finset.sum_add_distrib]
    apply Finset.sum_congr rfl
    intro i _hi
    ring
  simp_rw [sub_sq]
  simp only [Finset.sum_sub_distrib, Finset.sum_add_distrib]
  simp [Finset.mul_sum]
  rw [hmixed2]
  rw [← hsquares]
  ring

theorem scalar_pairwise_sq_of_centered {N : ℕ} (x : Fin N → ℝ)
    (hx : ∑ i, x i = 0) :
    ∑ i, ∑ j, (x i - x j) ^ 2 = 2 * N * ∑ i, (x i) ^ 2 := by
  rw [scalar_pairwise_sq_identity, hx]
  ring

/-- Centered all-pairs RMS identity in three dimensions.  In particular, a
Stiefel frame with total squared vertex norm `3` has total ordered-pair squared
distance `6N`. -/
theorem pairwise_normSq3_of_centered {N : ℕ} (p : Fin N → Vec3)
    (hp : ∀ k, ∑ i, p i k = 0) :
    ∑ i, ∑ j, normSq3 (p i - p j) = 2 * N * ∑ i, normSq3 (p i) := by
  have h0 := scalar_pairwise_sq_of_centered (fun i ↦ p i 0) (hp 0)
  have h1 := scalar_pairwise_sq_of_centered (fun i ↦ p i 1) (hp 1)
  have h2 := scalar_pairwise_sq_of_centered (fun i ↦ p i 2) (hp 2)
  simp only [normSq3, Pi.sub_apply, Finset.sum_add_distrib]
  linear_combination h0 + h1 + h2

theorem stiefel_pairwise_rms {N : ℕ} (p : Fin N → Vec3)
    (hp : ∀ k, ∑ i, p i k = 0)
    (hframe : ∑ i, normSq3 (p i) = 3) :
    ∑ i, ∑ j, normSq3 (p i - p j) = 6 * N := by
  rw [pairwise_normSq3_of_centered p hp, hframe]
  ring

/-- Oriented affine volume of four points. -/
def orient4 (a b c d : Vec3) : ℝ := det3 (b - a) (c - a) (d - a)

/-- Signed maximal minors of the augmented `4 × 5` affine-coordinate
matrix with columns `(a,b,c,p,q)`. -/
def circuit5 (a b c p q : Vec3) : Fin 5 → ℝ :=
  ![orient4 b c p q,
    -orient4 a c p q,
    orient4 a b p q,
    -orient4 a b c q,
    orient4 a b c p]

/-- First half of the rectangular cofactor identity: the signed maximal
minors sum to zero. -/
theorem circuit5_sum (a b c p q : Vec3) :
    circuit5 a b c p q 0 + circuit5 a b c p q 1 +
      circuit5 a b c p q 2 + circuit5 a b c p q 3 +
      circuit5 a b c p q 4 = 0 := by
  simp [circuit5, orient4, det3]
  ring

/-- Second half of the cofactor identity: the signed maximal minors give an
affine dependence of the five points. -/
theorem circuit5_affine_dependence (a b c p q : Vec3) :
    (circuit5 a b c p q 0) • a +
      (circuit5 a b c p q 1) • b +
      (circuit5 a b c p q 2) • c +
      (circuit5 a b c p q 3) • p +
      (circuit5 a b c p q 4) • q = 0 := by
  funext i
  fin_cases i <;> simp [circuit5, orient4, det3] <;> ring

/-- A strictly signed affine circuit normalizes to an equality between a
point in the open triangle and a point in the open segment.  This is the
sign-to-intersection half of the segment--triangle circuit criterion. -/
theorem positive_circuit_gives_open_segment_triangle
    (a b c p q : Vec3) (z : Fin 5 → ℝ)
    (hzsum : z 0 + z 1 + z 2 + z 3 + z 4 = 0)
    (hzvec : (z 0) • a + (z 1) • b + (z 2) • c +
        (z 3) • p + (z 4) • q = 0)
    (h0 : 0 < z 0) (h1 : 0 < z 1) (h2 : 0 < z 2)
    (h3 : z 3 < 0) (h4 : z 4 < 0) :
    ∃ α β γ u v : ℝ,
      0 < α ∧ 0 < β ∧ 0 < γ ∧ 0 < u ∧ 0 < v ∧
      α + β + γ = 1 ∧ u + v = 1 ∧
      α • a + β • b + γ • c = u • p + v • q := by
  let A : ℝ := z 0 + z 1 + z 2
  have hA : 0 < A := by dsimp [A]; linarith
  refine ⟨z 0 / A, z 1 / A, z 2 / A, -z 3 / A, -z 4 / A, ?_⟩
  constructor
  · exact div_pos h0 hA
  constructor
  · exact div_pos h1 hA
  constructor
  · exact div_pos h2 hA
  constructor
  · exact div_pos (neg_pos.mpr h3) hA
  constructor
  · exact div_pos (neg_pos.mpr h4) hA
  constructor
  · dsimp [A]
    field_simp
  constructor
  · dsimp [A]
    rw [← add_div]
    have heq : -z 3 + -z 4 = z 0 + z 1 + z 2 := by linarith
    rw [heq]
    exact div_self hA.ne'
  · funext i
    have hi := congrFun hzvec i
    change z 0 * a i + z 1 * b i + z 2 * c i + z 3 * p i + z 4 * q i = 0 at hi
    change z 0 / A * a i + z 1 / A * b i + z 2 / A * c i =
      -z 3 / A * p i + -z 4 / A * q i
    field_simp
    linarith

/-- The concrete cofactor-sign test used by the classifier: a `3 versus 2`
sign pattern in `circuit5` produces an open segment--open triangle
intersection. -/
theorem circuit5_signs_give_open_segment_triangle
    (a b c p q : Vec3)
    (h0 : 0 < circuit5 a b c p q 0)
    (h1 : 0 < circuit5 a b c p q 1)
    (h2 : 0 < circuit5 a b c p q 2)
    (h3 : circuit5 a b c p q 3 < 0)
    (h4 : circuit5 a b c p q 4 < 0) :
    ∃ α β γ u v : ℝ,
      0 < α ∧ 0 < β ∧ 0 < γ ∧ 0 < u ∧ 0 < v ∧
      α + β + γ = 1 ∧ u + v = 1 ∧
      α • a + β • b + γ • c = u • p + v • q := by
  exact positive_circuit_gives_open_segment_triangle
    a b c p q (circuit5 a b c p q)
    (circuit5_sum a b c p q) (circuit5_affine_dependence a b c p q)
    h0 h1 h2 h3 h4

/-! ## Executable Calvo decision layer -/

inductive OrientationSign where
  | neg
  | pos
  deriving DecidableEq, Repr

namespace OrientationSign

def flip : OrientationSign → OrientationSign
  | neg => pos
  | pos => neg

def toInt : OrientationSign → ℤ
  | neg => -1
  | pos => 1

@[simp] theorem flip_flip (s : OrientationSign) : flip (flip s) = s := by cases s <;> rfl

@[simp] theorem toInt_flip (s : OrientationSign) : toInt (flip s) = -toInt s := by
  cases s <;> rfl

end OrientationSign

structure PiercingSigns where
  c0 : OrientationSign
  c1 : OrientationSign
  c2 : OrientationSign
  c3 : OrientationSign
  c4 : OrientationSign
  orientedPiercing : OrientationSign
  deriving DecidableEq, Repr

def PiercingSigns.flip (x : PiercingSigns) : PiercingSigns where
  c0 := x.c0.flip
  c1 := x.c1.flip
  c2 := x.c2.flip
  c3 := x.c3.flip
  c4 := x.c4.flip
  orientedPiercing := x.orientedPiercing.flip

/-- The open segment pierces the open triangle exactly for the `3 versus 2`
Radon sign partition.  Its contribution is then the oriented determinant
sign; otherwise it contributes zero. -/
def piercingContribution (x : PiercingSigns) : ℤ :=
  if x.c0 = x.c1 ∧ x.c1 = x.c2 ∧ x.c3 = x.c4 ∧ x.c0 = x.c3.flip then
    x.orientedPiercing.toInt
  else 0

theorem piercingContribution_flip (x : PiercingSigns) :
    piercingContribution x.flip = -piercingContribution x := by
  rcases x with ⟨c0, c1, c2, c3, c4, o⟩
  cases c0 <;> cases c1 <;> cases c2 <;> cases c3 <;> cases c4 <;> cases o <;> decide

structure CalvoData where
  d2a : PiercingSigns
  d2b : PiercingSigns
  d4a : PiercingSigns
  d4b : PiercingSigns
  d6a : PiercingSigns
  d6b : PiercingSigns
  deriving DecidableEq, Repr

def CalvoData.flip (x : CalvoData) : CalvoData where
  d2a := x.d2a.flip
  d2b := x.d2b.flip
  d4a := x.d4a.flip
  d4b := x.d4b.flip
  d6a := x.d6a.flip
  d6b := x.d6b.flip

def calvoDeltas (x : CalvoData) : ℤ × ℤ × ℤ :=
  (piercingContribution x.d2a + piercingContribution x.d2b,
   piercingContribution x.d4a + piercingContribution x.d4b,
   piercingContribution x.d6a + piercingContribution x.d6b)

theorem calvoDeltas_flip (x : CalvoData) :
    calvoDeltas x.flip =
      (-(calvoDeltas x).1, -(calvoDeltas x).2.1, -(calvoDeltas x).2.2) := by
  simp only [calvoDeltas, CalvoData.flip, piercingContribution_flip]
  ring_nf

inductive CalvoDecision where
  | rightTrefoil
  | leftTrefoil
  | unknot
  | inadmissible
  deriving DecidableEq, Repr

/-- Calvo's topological classification theorem is external; this function is
the fully checked deterministic decision layer which consumes its three
intersection numbers. -/
def calvoDecision (x : CalvoData) : CalvoDecision :=
  let d := calvoDeltas x
  if d = (1, 1, 1) then .rightTrefoil
  else if d = (-1, -1, -1) then .leftTrefoil
  else if d.1 = 0 ∨ d.2.1 = 0 ∨ d.2.2 = 0 then .unknot
  else .inadmissible

private def miss : PiercingSigns :=
  ⟨.neg, .neg, .pos, .pos, .pos, .neg⟩

private def positivePiercing : PiercingSigns :=
  ⟨.pos, .pos, .pos, .neg, .neg, .pos⟩

/-- The six circuit tests obtained from the canonical 15-sign order-type
signature `+-++-++--+-++-+`. -/
def canonicalTrefoilCalvoData : CalvoData :=
  ⟨miss, positivePiercing, miss, positivePiercing, miss, positivePiercing⟩

example : calvoDeltas canonicalTrefoilCalvoData = (1, 1, 1) := by decide

example : calvoDecision canonicalTrefoilCalvoData = .rightTrefoil := by decide

end GaussianKnots
