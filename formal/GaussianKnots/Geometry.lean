import Mathlib

/-!
# Polynomial walls for polygonal knots

This file verifies the elementary geometric implication behind the wall
construction.  Segment intersection is encoded by an affine equation; no
genericity, topology, or appeal to a picture enters the proofs.
-/

namespace GaussianKnots

abbrev Vec3 := Fin 3 → ℝ

/-- The explicit scalar triple product.  This is the cubic used for a
nonincident pair of polygon edges. -/
def det3 (u v w : Vec3) : ℝ :=
  u 0 * (v 1 * w 2 - v 2 * w 1)
    - u 1 * (v 0 * w 2 - v 2 * w 0)
    + u 2 * (v 0 * w 1 - v 1 * w 0)

/-- The explicit cross product. -/
def cross3 (u v : Vec3) : Vec3 :=
  ![u 1 * v 2 - u 2 * v 1,
    u 2 * v 0 - u 0 * v 2,
    u 0 * v 1 - u 1 * v 0]

/-- A degree-four polynomial which vanishes precisely when two edge vectors
are collinear. -/
def crossSq (u v : Vec3) : ℝ :=
  (cross3 u v 0) ^ 2 + (cross3 u v 1) ^ 2 + (cross3 u v 2) ^ 2

theorem det3_smul_add_self (u v : Vec3) (a b : ℝ) :
    det3 u v (a • u + b • v) = 0 := by
  simp [det3]
  ring

/-- If the supporting lines of two directed segments meet, their edge
vectors and displacement are coplanar.  Segment-parameter inequalities are
unnecessary: the affine intersection equation alone implies the cubic wall. -/
theorem segment_intersection_implies_det_wall
    (a b c d : Vec3) (t s : ℝ)
    (h : a + t • (b - a) = c + s • (d - c)) :
    det3 (b - a) (d - c) (c - a) = 0 := by
  have hw : c - a = t • (b - a) + (-s) • (d - c) := by
    funext i
    have hi := congrFun h i
    change a i + t * (b i - a i) = c i + s * (d i - c i) at hi
    change c i - a i = t * (b i - a i) + (-s) * (d i - c i)
    linarith
  rw [hw]
  exact det3_smul_add_self (b - a) (d - c) t (-s)

private theorem cross3_eq_zero_of_smul_eq_smul
    (u v : Vec3) (a b : ℝ) (ha : a ≠ 0) (h : a • u = b • v) :
    cross3 u v = 0 := by
  have h0 := congrFun h 0
  have h1 := congrFun h 1
  have h2 := congrFun h 2
  change a * u 0 = b * v 0 at h0
  change a * u 1 = b * v 1 at h1
  change a * u 2 = b * v 2 at h2
  funext i
  fin_cases i
  · change u 1 * v 2 - u 2 * v 1 = 0
    have hz : a * (u 1 * v 2 - u 2 * v 1) = 0 := by
      calc
        a * (u 1 * v 2 - u 2 * v 1)
            = (a * u 1) * v 2 - (a * u 2) * v 1 := by ring
        _ = (b * v 1) * v 2 - (b * v 2) * v 1 := by rw [h1, h2]
        _ = 0 := by ring
    exact (mul_eq_zero.mp hz).resolve_left ha
  · change u 2 * v 0 - u 0 * v 2 = 0
    have hz : a * (u 2 * v 0 - u 0 * v 2) = 0 := by
      calc
        a * (u 2 * v 0 - u 0 * v 2)
            = (a * u 2) * v 0 - (a * u 0) * v 2 := by ring
        _ = (b * v 2) * v 0 - (b * v 0) * v 2 := by rw [h2, h0]
        _ = 0 := by ring
    exact (mul_eq_zero.mp hz).resolve_left ha
  · change u 0 * v 1 - u 1 * v 0 = 0
    have hz : a * (u 0 * v 1 - u 1 * v 0) = 0 := by
      calc
        a * (u 0 * v 1 - u 1 * v 0)
            = (a * u 0) * v 1 - (a * u 1) * v 0 := by ring
        _ = (b * v 0) * v 1 - (b * v 1) * v 0 := by rw [h0, h1]
        _ = 0 := by ring
    exact (mul_eq_zero.mp hz).resolve_left ha

/-- An intersection of consecutive edges away from their common endpoint
forces the adjacent-edge (quartic) wall. -/
theorem adjacent_improper_intersection_implies_cross_wall
    (a b c : Vec3) (t s : ℝ) (ht : t ≠ 1)
    (h : a + t • (b - a) = b + s • (c - b)) :
    cross3 (b - a) (c - b) = 0 := by
  have hsmul : (t - 1) • (b - a) = s • (c - b) := by
    funext i
    have hi := congrFun h i
    change a i + t * (b i - a i) = b i + s * (c i - b i) at hi
    change (t - 1) * (b i - a i) = s * (c i - b i)
    linarith
  apply cross3_eq_zero_of_smul_eq_smul (b - a) (c - b) (t - 1) s
  · exact sub_ne_zero.mpr ht
  · exact hsmul

/-- Over the reals the sum-of-squares quartic has no extra zeroes: it is
equivalent to the vector cross-product equation. -/
theorem crossSq_eq_zero_iff (u v : Vec3) :
    crossSq u v = 0 ↔ cross3 u v = 0 := by
  constructor
  · intro h
    have h0 : 0 ≤ (cross3 u v 0) ^ 2 := sq_nonneg _
    have h1 : 0 ≤ (cross3 u v 1) ^ 2 := sq_nonneg _
    have h2 : 0 ≤ (cross3 u v 2) ^ 2 := sq_nonneg _
    have hz0 : cross3 u v 0 = 0 := by
      unfold crossSq at h
      nlinarith
    have hz1 : cross3 u v 1 = 0 := by
      unfold crossSq at h
      nlinarith
    have hz2 : cross3 u v 2 = 0 := by
      unfold crossSq at h
      nlinarith
    funext i
    fin_cases i
    · exact hz0
    · exact hz1
    · exact hz2
  · intro h
    have h0 := congrFun h 0
    have h1 := congrFun h 1
    have h2 := congrFun h 2
    unfold crossSq
    rw [h0, h1, h2]
    norm_num

theorem adjacent_improper_intersection_implies_crossSq_wall
    (a b c : Vec3) (t s : ℝ) (ht : t ≠ 1)
    (h : a + t • (b - a) = b + s • (c - b)) :
    crossSq (b - a) (c - b) = 0 := by
  rw [crossSq_eq_zero_iff]
  exact adjacent_improper_intersection_implies_cross_wall a b c t s ht h

/-- Reflection in the first coordinate reverses every determinant sign. -/
def reflectX (u : Vec3) : Vec3 := ![-u 0, u 1, u 2]

theorem det3_reflectX (u v w : Vec3) :
    det3 (reflectX u) (reflectX v) (reflectX w) = -det3 u v w := by
  simp [det3, reflectX]
  ring

theorem det3_swap_first_two (u v w : Vec3) :
    det3 v u w = -det3 u v w := by
  simp [det3]
  ring

end GaussianKnots
