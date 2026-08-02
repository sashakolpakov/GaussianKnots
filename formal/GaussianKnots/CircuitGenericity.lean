import GaussianKnots.HexagonSignature

/-!
# Genericity closes the five-point circuit criterion

For five points in three-space, `CircuitConverse` previously accepted as an
input the assertion that every affine dependence is a multiple of the signed
cofactor circuit.  This file proves that assertion from the concrete
genericity condition actually used in the geometry: one of the five maximal
minors is nonzero.

No topological theorem is involved here.  The argument is elementary linear
algebra over `ℝ`, reduced to the explicit determinant already used by the
executable sign classifier.
-/

namespace GaussianKnots

/-- Three vectors with nonzero scalar triple product are linearly independent.
This coordinate lemma keeps the generic circuit proof tied to the exact
`det3` polynomial used by the manuscript and classifier. -/
theorem three_smul_eq_zero_of_det3_ne
    (u v w : Vec3) (x y z : ℝ) (hdet : det3 u v w ≠ 0)
    (hlin : x • u + y • v + z • w = 0) :
    x = 0 ∧ y = 0 ∧ z = 0 := by
  have h0 := congrFun hlin 0
  have h1 := congrFun hlin 1
  have h2 := congrFun hlin 2
  change x * u 0 + y * v 0 + z * w 0 = 0 at h0
  change x * u 1 + y * v 1 + z * w 1 = 0 at h1
  change x * u 2 + y * v 2 + z * w 2 = 0 at h2
  have hxmul : x * det3 u v w = 0 := by
    calc
      x * det3 u v w =
          (v 1 * w 2 - v 2 * w 1) *
              (x * u 0 + y * v 0 + z * w 0) -
            (v 0 * w 2 - v 2 * w 0) *
              (x * u 1 + y * v 1 + z * w 1) +
            (v 0 * w 1 - v 1 * w 0) *
              (x * u 2 + y * v 2 + z * w 2) := by
                simp only [det3]
                ring
      _ = 0 := by rw [h0, h1, h2]; ring
  have hymul : y * det3 u v w = 0 := by
    calc
      y * det3 u v w =
          -(u 1 * w 2 - u 2 * w 1) *
              (x * u 0 + y * v 0 + z * w 0) +
            (u 0 * w 2 - u 2 * w 0) *
              (x * u 1 + y * v 1 + z * w 1) -
            (u 0 * w 1 - u 1 * w 0) *
              (x * u 2 + y * v 2 + z * w 2) := by
                simp only [det3]
                ring
      _ = 0 := by rw [h0, h1, h2]; ring
  have hzmul : z * det3 u v w = 0 := by
    calc
      z * det3 u v w =
          (u 1 * v 2 - u 2 * v 1) *
              (x * u 0 + y * v 0 + z * w 0) -
            (u 0 * v 2 - u 2 * v 0) *
              (x * u 1 + y * v 1 + z * w 1) +
            (u 0 * v 1 - u 1 * v 0) *
              (x * u 2 + y * v 2 + z * w 2) := by
                simp only [det3]
                ring
      _ = 0 := by rw [h0, h1, h2]; ring
  exact ⟨(mul_eq_zero.mp hxmul).resolve_right hdet,
    (mul_eq_zero.mp hymul).resolve_right hdet,
    (mul_eq_zero.mp hzmul).resolve_right hdet⟩

/-- If `a,b,c,p` have nonzero oriented affine volume, every affine dependence
of `a,b,c,p,q` is a scalar multiple of the signed cofactor circuit. -/
theorem affineDependence_eq_smul_circuit5_of_orient4_ne
    (a b c p q : Vec3) (hgeneric : orient4 a b c p ≠ 0)
    (z : Fin 5 → ℝ)
    (hzsum : z 0 + z 1 + z 2 + z 3 + z 4 = 0)
    (hzvec : (z 0) • a + (z 1) • b + (z 2) • c +
        (z 3) • p + (z 4) • q = 0) :
    ∃ r : ℝ, z = r • circuit5 a b c p q := by
  let C := circuit5 a b c p q
  have hC4 : C 4 ≠ 0 := by
    simpa [C, circuit5] using hgeneric
  let r : ℝ := z 4 / C 4
  let e : Fin 5 → ℝ := fun i ↦ z i - r * C i
  have he4 : e 4 = 0 := by
    dsimp [e, r]
    rw [div_mul_cancel₀ _ hC4]
    ring
  have hesum : e 0 + e 1 + e 2 + e 3 + e 4 = 0 := by
    dsimp [e]
    have hCsum : C 0 + C 1 + C 2 + C 3 + C 4 = 0 := by
      simpa [C] using circuit5_sum a b c p q
    linear_combination hzsum - r * hCsum
  have hevec :
      (e 0) • a + (e 1) • b + (e 2) • c +
          (e 3) • p + (e 4) • q = 0 := by
    have hCvec := circuit5_affine_dependence a b c p q
    funext i
    have hzi := congrFun hzvec i
    have hCi := congrFun hCvec i
    change z 0 * a i + z 1 * b i + z 2 * c i +
        z 3 * p i + z 4 * q i = 0 at hzi
    change C 0 * a i + C 1 * b i + C 2 * c i +
        C 3 * p i + C 4 * q i = 0 at hCi
    change e 0 * a i + e 1 * b i + e 2 * c i +
        e 3 * p i + e 4 * q i = 0
    dsimp [e]
    linear_combination hzi - r * hCi
  have hrelative :
      (e 1) • (b - a) + (e 2) • (c - a) + (e 3) • (p - a) = 0 := by
    funext i
    have hvi := congrFun hevec i
    change e 0 * a i + e 1 * b i + e 2 * c i +
        e 3 * p i + e 4 * q i = 0 at hvi
    change e 1 * (b i - a i) + e 2 * (c i - a i) +
        e 3 * (p i - a i) = 0
    linear_combination hvi - a i * hesum - (q i - a i) * he4
  have hdet : det3 (b - a) (c - a) (p - a) ≠ 0 := by
    simpa [orient4] using hgeneric
  obtain ⟨he1, he2, he3⟩ := three_smul_eq_zero_of_det3_ne
    (b - a) (c - a) (p - a) (e 1) (e 2) (e 3) hdet hrelative
  have he0 : e 0 = 0 := by linarith
  refine ⟨r, ?_⟩
  funext i
  fin_cases i
  · change z 0 = r * C 0
    apply sub_eq_zero.mp
    simpa [e] using he0
  · change z 1 = r * C 1
    apply sub_eq_zero.mp
    simpa [e] using he1
  · change z 2 = r * C 2
    apply sub_eq_zero.mp
    simpa [e] using he2
  · change z 3 = r * C 3
    apply sub_eq_zero.mp
    simpa [e] using he3
  · change z 4 = r * C 4
    apply sub_eq_zero.mp
    simpa [e] using he4

/-- The generic five-point circuit test, with no conclusion-shaped
one-dimensionality hypothesis left external. -/
theorem circuitSignSplit_iff_openSegmentTriangle_of_orient4_ne
    (a b c p q : Vec3) (hgeneric : orient4 a b c p ≠ 0) :
    CircuitSignSplit a b c p q ↔ OpenSegmentTriangle a b c p q := by
  apply circuitSignSplit_iff_openSegmentTriangle_of_unique
  intro z hzsum hzvec
  exact affineDependence_eq_smul_circuit5_of_orient4_ne
    a b c p q hgeneric z hzsum hzvec

/-! ## The executable piercing predicate has the same geometry -/

/-- The Boolean-free proposition tested by `piercingContribution` before it
returns the oriented sign. -/
def PiercingDetected (x : PiercingSigns) : Prop :=
  x.c0 = x.c1 ∧ x.c1 = x.c2 ∧ x.c3 = x.c4 ∧ x.c0 = x.c3.flip

/-- `piercingContribution` is nonzero exactly when its executable circuit
test detects a piercing; the remaining field only chooses the orientation. -/
theorem piercingContribution_ne_zero_iff (x : PiercingSigns) :
    piercingContribution x ≠ 0 ↔ PiercingDetected x := by
  constructor
  · intro hn
    by_contra hdet
    have hif : ¬(x.c0 = x.c1 ∧ x.c1 = x.c2 ∧ x.c3 = x.c4 ∧
        x.c0 = x.c3.flip) := by simpa only [PiercingDetected] using hdet
    rw [piercingContribution, if_neg hif] at hn
    exact hn rfl
  · intro hdet
    have hif : x.c0 = x.c1 ∧ x.c1 = x.c2 ∧ x.c3 = x.c4 ∧
        x.c0 = x.c3.flip := by simpa only [PiercingDetected] using hdet
    rw [piercingContribution, if_pos hif]
    cases x.orientedPiercing <;> decide

/-- The geometric sign record before it is packaged as a labeled hexagon
query. -/
noncomputable def circuitPiercingSigns (a b c p q : Vec3) : PiercingSigns :=
  { c0 := OrientationSign.ofReal (circuit5 a b c p q 0)
    c1 := OrientationSign.ofReal (circuit5 a b c p q 1)
    c2 := OrientationSign.ofReal (circuit5 a b c p q 2)
    c3 := OrientationSign.ofReal (circuit5 a b c p q 3)
    c4 := OrientationSign.ofReal (circuit5 a b c p q 4)
    orientedPiercing := OrientationSign.ofReal (orient4 a b c q) }

@[simp] theorem circuitPiercingSigns_eq_geometricSigns
    (query : PiercingQuery) (v : Hexagon) :
    circuitPiercingSigns (v query.a) (v query.b) (v query.c)
      (v query.p) (v query.q) = query.geometricSigns v := rfl

private theorem ofReal_eq_pos_iff {x : ℝ} (hx : x ≠ 0) :
    OrientationSign.ofReal x = .pos ↔ 0 < x := by
  rcases lt_or_gt_of_ne hx with hx | hx
  · simp [OrientationSign.ofReal, hx, hx.not_gt]
  · simp [OrientationSign.ofReal, hx, hx.not_gt]

private theorem ofReal_eq_neg_iff {x : ℝ} (hx : x ≠ 0) :
    OrientationSign.ofReal x = .neg ↔ x < 0 := by
  rcases lt_or_gt_of_ne hx with hx | hx
  · simp [OrientationSign.ofReal, hx]
  · simp [OrientationSign.ofReal, hx.not_gt]

/-- Once the five cofactors are nonzero, equality of their two-valued signs
is exactly the strict `3 versus 2` circuit split. -/
theorem piercingDetected_circuitSigns_iff_circuitSignSplit
    (a b c p q : Vec3) (hcoeff : ∀ i, circuit5 a b c p q i ≠ 0) :
    PiercingDetected (circuitPiercingSigns a b c p q) ↔
      CircuitSignSplit a b c p q := by
  let C := circuit5 a b c p q
  have hn0 : C 0 ≠ 0 := hcoeff 0
  have hn1 : C 1 ≠ 0 := hcoeff 1
  have hn2 : C 2 ≠ 0 := hcoeff 2
  have hn3 : C 3 ≠ 0 := hcoeff 3
  have hn4 : C 4 ≠ 0 := hcoeff 4
  constructor
  · rintro ⟨h01, h12, h34, h03⟩
    change OrientationSign.ofReal (C 0) = OrientationSign.ofReal (C 1) at h01
    change OrientationSign.ofReal (C 1) = OrientationSign.ofReal (C 2) at h12
    change OrientationSign.ofReal (C 3) = OrientationSign.ofReal (C 4) at h34
    change OrientationSign.ofReal (C 0) =
      (OrientationSign.ofReal (C 3)).flip at h03
    unfold CircuitSignSplit
    change (0 < C 0 ∧ 0 < C 1 ∧ 0 < C 2 ∧ C 3 < 0 ∧ C 4 < 0) ∨
      (C 0 < 0 ∧ C 1 < 0 ∧ C 2 < 0 ∧ 0 < C 3 ∧ 0 < C 4)
    generalize hs0 : OrientationSign.ofReal (C 0) = s0 at h01 h03
    cases s0 with
    | neg =>
        right
        have hs1 : OrientationSign.ofReal (C 1) = .neg := by simpa [hs0] using h01.symm
        have hs2 : OrientationSign.ofReal (C 2) = .neg := h12.symm.trans hs1
        have hs3 : OrientationSign.ofReal (C 3) = .pos := by
          generalize hs3val : OrientationSign.ofReal (C 3) = s3 at h03
          cases s3 with
          | neg => simp [OrientationSign.flip] at h03
          | pos => rfl
        have hs4 : OrientationSign.ofReal (C 4) = .pos := h34.symm.trans hs3
        exact ⟨(ofReal_eq_neg_iff hn0).mp hs0,
          (ofReal_eq_neg_iff hn1).mp hs1,
          (ofReal_eq_neg_iff hn2).mp hs2,
          (ofReal_eq_pos_iff hn3).mp hs3,
          (ofReal_eq_pos_iff hn4).mp hs4⟩
    | pos =>
        left
        have hs1 : OrientationSign.ofReal (C 1) = .pos := by simpa [hs0] using h01.symm
        have hs2 : OrientationSign.ofReal (C 2) = .pos := h12.symm.trans hs1
        have hs3 : OrientationSign.ofReal (C 3) = .neg := by
          generalize hs3val : OrientationSign.ofReal (C 3) = s3 at h03
          cases s3 with
          | neg => rfl
          | pos => simp [OrientationSign.flip] at h03
        have hs4 : OrientationSign.ofReal (C 4) = .neg := h34.symm.trans hs3
        exact ⟨(ofReal_eq_pos_iff hn0).mp hs0,
          (ofReal_eq_pos_iff hn1).mp hs1,
          (ofReal_eq_pos_iff hn2).mp hs2,
          (ofReal_eq_neg_iff hn3).mp hs3,
          (ofReal_eq_neg_iff hn4).mp hs4⟩
  · intro hsplit
    rcases hsplit with hpos | hneg
    · simp [PiercingDetected, circuitPiercingSigns,
        OrientationSign.flip,
        OrientationSign.ofReal_of_pos hpos.1,
        OrientationSign.ofReal_of_pos hpos.2.1,
        OrientationSign.ofReal_of_pos hpos.2.2.1,
        OrientationSign.ofReal_of_neg hpos.2.2.2.1,
        OrientationSign.ofReal_of_neg hpos.2.2.2.2]
    · simp [PiercingDetected, circuitPiercingSigns,
        OrientationSign.flip,
        OrientationSign.ofReal_of_neg hneg.1,
        OrientationSign.ofReal_of_neg hneg.2.1,
        OrientationSign.ofReal_of_neg hneg.2.2.1,
        OrientationSign.ofReal_of_pos hneg.2.2.2.1,
        OrientationSign.ofReal_of_pos hneg.2.2.2.2]

theorem piercingDetected_geometricSigns_iff_circuitSignSplit
    (query : PiercingQuery) (v : Hexagon)
    (hcoeff : ∀ i, circuit5 (v query.a) (v query.b) (v query.c)
      (v query.p) (v query.q) i ≠ 0) :
    PiercingDetected (query.geometricSigns v) ↔
      CircuitSignSplit (v query.a) (v query.b) (v query.c)
        (v query.p) (v query.q) := by
  rw [← circuitPiercingSigns_eq_geometricSigns]
  exact piercingDetected_circuitSigns_iff_circuitSignSplit _ _ _ _ _ hcoeff

/-- Consequently, the sign record computed from a geometric query detects
the actual open segment--open triangle intersection. -/
theorem piercingDetected_geometricSigns_iff_openSegmentTriangle
    (query : PiercingQuery) (v : Hexagon)
    (hcoeff : ∀ i, circuit5 (v query.a) (v query.b) (v query.c)
      (v query.p) (v query.q) i ≠ 0) :
    PiercingDetected (query.geometricSigns v) ↔
      OpenSegmentTriangle (v query.a) (v query.b) (v query.c)
        (v query.p) (v query.q) := by
  rw [piercingDetected_geometricSigns_iff_circuitSignSplit query v hcoeff]
  apply circuitSignSplit_iff_openSegmentTriangle_of_orient4_ne
  simpa [circuit5] using hcoeff 4

/-! ## The six Calvo queries of a generic hexagon -/

def calvoD2aQuery : PiercingQuery := ⟨0, 1, 2, 3, 4⟩
def calvoD2bQuery : PiercingQuery := ⟨0, 1, 2, 4, 5⟩
def calvoD4aQuery : PiercingQuery := ⟨2, 3, 4, 5, 0⟩
def calvoD4bQuery : PiercingQuery := ⟨2, 3, 4, 0, 1⟩
def calvoD6aQuery : PiercingQuery := ⟨4, 5, 0, 1, 2⟩
def calvoD6bQuery : PiercingQuery := ⟨4, 5, 0, 2, 3⟩

theorem orient4_rotate_one_genericity (a b c d : Vec3) :
    orient4 b c d a = -orient4 a b c d := by
  simp [orient4, det3]
  ring

theorem orient4_rotate_two_genericity (a b c d : Vec3) :
    orient4 c d a b = orient4 a b c d := by
  simp [orient4, det3]
  ring

theorem orient4_rotate_three_genericity (a b c d : Vec3) :
    orient4 d a b c = -orient4 a b c d := by
  simp [orient4, det3]
  ring

private theorem genericHexagon_d2a_coefficients_ne (v : Hexagon)
    (hv : GenericHexagon v) :
    ∀ i, circuit5 (v 0) (v 1) (v 2) (v 3) (v 4) i ≠ 0 := by
  intro i
  fin_cases i <;> simp [circuit5] at *
  · exact hv 10
  · exact hv 6
  · exact hv 3
  · exact hv 1
  · exact hv 0

private theorem genericHexagon_d2b_coefficients_ne (v : Hexagon)
    (hv : GenericHexagon v) :
    ∀ i, circuit5 (v 0) (v 1) (v 2) (v 4) (v 5) i ≠ 0 := by
  intro i
  fin_cases i <;> simp [circuit5] at *
  · exact hv 12
  · exact hv 8
  · exact hv 5
  · exact hv 2
  · exact hv 1

private theorem genericHexagon_d4a_coefficients_ne (v : Hexagon)
    (hv : GenericHexagon v) :
    ∀ i, circuit5 (v 2) (v 3) (v 4) (v 5) (v 0) i ≠ 0 := by
  intro i
  fin_cases i
  · change orient4 (v 3) (v 4) (v 5) (v 0) ≠ 0
    rw [orient4_rotate_one_genericity]
    exact neg_ne_zero.mpr (by simpa [hexagonOrientationValues] using hv 9)
  · change -orient4 (v 2) (v 4) (v 5) (v 0) ≠ 0
    rw [orient4_rotate_one_genericity]
    simpa [hexagonOrientationValues] using hv 8
  · change orient4 (v 2) (v 3) (v 5) (v 0) ≠ 0
    rw [orient4_rotate_one_genericity]
    exact neg_ne_zero.mpr (by simpa [hexagonOrientationValues] using hv 7)
  · change -orient4 (v 2) (v 3) (v 4) (v 0) ≠ 0
    rw [orient4_rotate_one_genericity]
    simpa [hexagonOrientationValues] using hv 6
  · simpa [circuit5, hexagonOrientationValues] using hv 14

private theorem genericHexagon_d4b_coefficients_ne (v : Hexagon)
    (hv : GenericHexagon v) :
    ∀ i, circuit5 (v 2) (v 3) (v 4) (v 0) (v 1) i ≠ 0 := by
  intro i
  fin_cases i
  · change orient4 (v 3) (v 4) (v 0) (v 1) ≠ 0
    rw [orient4_rotate_two_genericity]
    simpa [hexagonOrientationValues] using hv 3
  · change -orient4 (v 2) (v 4) (v 0) (v 1) ≠ 0
    rw [orient4_rotate_two_genericity]
    exact neg_ne_zero.mpr (by simpa [hexagonOrientationValues] using hv 1)
  · change orient4 (v 2) (v 3) (v 0) (v 1) ≠ 0
    rw [orient4_rotate_two_genericity]
    simpa [hexagonOrientationValues] using hv 0
  · change -orient4 (v 2) (v 3) (v 4) (v 1) ≠ 0
    rw [orient4_rotate_one_genericity]
    simpa [hexagonOrientationValues] using hv 10
  · change orient4 (v 2) (v 3) (v 4) (v 0) ≠ 0
    rw [orient4_rotate_one_genericity]
    exact neg_ne_zero.mpr (by simpa [hexagonOrientationValues] using hv 6)

private theorem genericHexagon_d6a_coefficients_ne (v : Hexagon)
    (hv : GenericHexagon v) :
    ∀ i, circuit5 (v 4) (v 5) (v 0) (v 1) (v 2) i ≠ 0 := by
  intro i
  fin_cases i
  · change orient4 (v 5) (v 0) (v 1) (v 2) ≠ 0
    rw [orient4_rotate_three_genericity]
    exact neg_ne_zero.mpr (by simpa [hexagonOrientationValues] using hv 2)
  · change -orient4 (v 4) (v 0) (v 1) (v 2) ≠ 0
    rw [orient4_rotate_three_genericity]
    simpa [hexagonOrientationValues] using hv 1
  · change orient4 (v 4) (v 5) (v 1) (v 2) ≠ 0
    rw [orient4_rotate_two_genericity]
    simpa [hexagonOrientationValues] using hv 12
  · change -orient4 (v 4) (v 5) (v 0) (v 2) ≠ 0
    rw [orient4_rotate_two_genericity]
    exact neg_ne_zero.mpr (by simpa [hexagonOrientationValues] using hv 8)
  · change orient4 (v 4) (v 5) (v 0) (v 1) ≠ 0
    rw [orient4_rotate_two_genericity]
    simpa [hexagonOrientationValues] using hv 5

private theorem genericHexagon_d6b_coefficients_ne (v : Hexagon)
    (hv : GenericHexagon v) :
    ∀ i, circuit5 (v 4) (v 5) (v 0) (v 2) (v 3) i ≠ 0 := by
  intro i
  fin_cases i
  · change orient4 (v 5) (v 0) (v 2) (v 3) ≠ 0
    rw [orient4_rotate_three_genericity]
    exact neg_ne_zero.mpr (by simpa [hexagonOrientationValues] using hv 7)
  · change -orient4 (v 4) (v 0) (v 2) (v 3) ≠ 0
    rw [orient4_rotate_three_genericity]
    simpa [hexagonOrientationValues] using hv 6
  · change orient4 (v 4) (v 5) (v 2) (v 3) ≠ 0
    rw [orient4_rotate_two_genericity]
    simpa [hexagonOrientationValues] using hv 14
  · change -orient4 (v 4) (v 5) (v 0) (v 3) ≠ 0
    rw [orient4_rotate_two_genericity]
    exact neg_ne_zero.mpr (by simpa [hexagonOrientationValues] using hv 9)
  · change orient4 (v 4) (v 5) (v 0) (v 2) ≠ 0
    rw [orient4_rotate_two_genericity]
    simpa [hexagonOrientationValues] using hv 8

/-- For the first `d₂` query, the executable Calvo contribution is nonzero
exactly when the corresponding open segment pierces the open triangle. -/
theorem calvoD2a_contribution_ne_zero_iff_openSegmentTriangle
    (v : Hexagon) (hv : GenericHexagon v) :
    piercingContribution (calvoDataOfHexagon v).d2a ≠ 0 ↔
      OpenSegmentTriangle (v 0) (v 1) (v 2) (v 3) (v 4) := by
  have hd := congrArg (fun x : CalvoData ↦ x.d2a)
    (calvoDataOfHexagon_eq_geometric v hv)
  change (calvoDataOfHexagon v).d2a = calvoD2aQuery.geometricSigns v at hd
  rw [piercingContribution_ne_zero_iff, hd]
  exact piercingDetected_geometricSigns_iff_openSegmentTriangle
    calvoD2aQuery v (genericHexagon_d2a_coefficients_ne v hv)

/-- The second `d₂` query has the same certified interpretation. -/
theorem calvoD2b_contribution_ne_zero_iff_openSegmentTriangle
    (v : Hexagon) (hv : GenericHexagon v) :
    piercingContribution (calvoDataOfHexagon v).d2b ≠ 0 ↔
      OpenSegmentTriangle (v 0) (v 1) (v 2) (v 4) (v 5) := by
  have hd := congrArg (fun x : CalvoData ↦ x.d2b)
    (calvoDataOfHexagon_eq_geometric v hv)
  change (calvoDataOfHexagon v).d2b = calvoD2bQuery.geometricSigns v at hd
  rw [piercingContribution_ne_zero_iff, hd]
  exact piercingDetected_geometricSigns_iff_openSegmentTriangle
    calvoD2bQuery v (genericHexagon_d2b_coefficients_ne v hv)

/-- The first `d₄` query has the same certified interpretation. -/
theorem calvoD4a_contribution_ne_zero_iff_openSegmentTriangle
    (v : Hexagon) (hv : GenericHexagon v) :
    piercingContribution (calvoDataOfHexagon v).d4a ≠ 0 ↔
      OpenSegmentTriangle (v 2) (v 3) (v 4) (v 5) (v 0) := by
  have hd := congrArg (fun x : CalvoData ↦ x.d4a)
    (calvoDataOfHexagon_eq_geometric v hv)
  change (calvoDataOfHexagon v).d4a = calvoD4aQuery.geometricSigns v at hd
  rw [piercingContribution_ne_zero_iff, hd]
  exact piercingDetected_geometricSigns_iff_openSegmentTriangle
    calvoD4aQuery v (genericHexagon_d4a_coefficients_ne v hv)

/-- The second `d₄` query has the same certified interpretation. -/
theorem calvoD4b_contribution_ne_zero_iff_openSegmentTriangle
    (v : Hexagon) (hv : GenericHexagon v) :
    piercingContribution (calvoDataOfHexagon v).d4b ≠ 0 ↔
      OpenSegmentTriangle (v 2) (v 3) (v 4) (v 0) (v 1) := by
  have hd := congrArg (fun x : CalvoData ↦ x.d4b)
    (calvoDataOfHexagon_eq_geometric v hv)
  change (calvoDataOfHexagon v).d4b = calvoD4bQuery.geometricSigns v at hd
  rw [piercingContribution_ne_zero_iff, hd]
  exact piercingDetected_geometricSigns_iff_openSegmentTriangle
    calvoD4bQuery v (genericHexagon_d4b_coefficients_ne v hv)

/-- The first `d₆` query has the same certified interpretation. -/
theorem calvoD6a_contribution_ne_zero_iff_openSegmentTriangle
    (v : Hexagon) (hv : GenericHexagon v) :
    piercingContribution (calvoDataOfHexagon v).d6a ≠ 0 ↔
      OpenSegmentTriangle (v 4) (v 5) (v 0) (v 1) (v 2) := by
  have hd := congrArg (fun x : CalvoData ↦ x.d6a)
    (calvoDataOfHexagon_eq_geometric v hv)
  change (calvoDataOfHexagon v).d6a = calvoD6aQuery.geometricSigns v at hd
  rw [piercingContribution_ne_zero_iff, hd]
  exact piercingDetected_geometricSigns_iff_openSegmentTriangle
    calvoD6aQuery v (genericHexagon_d6a_coefficients_ne v hv)

/-- The second `d₆` query has the same certified interpretation. -/
theorem calvoD6b_contribution_ne_zero_iff_openSegmentTriangle
    (v : Hexagon) (hv : GenericHexagon v) :
    piercingContribution (calvoDataOfHexagon v).d6b ≠ 0 ↔
      OpenSegmentTriangle (v 4) (v 5) (v 0) (v 2) (v 3) := by
  have hd := congrArg (fun x : CalvoData ↦ x.d6b)
    (calvoDataOfHexagon_eq_geometric v hv)
  change (calvoDataOfHexagon v).d6b = calvoD6bQuery.geometricSigns v at hd
  rw [piercingContribution_ne_zero_iff, hd]
  exact piercingDetected_geometricSigns_iff_openSegmentTriangle
    calvoD6bQuery v (genericHexagon_d6b_coefficients_ne v hv)

end GaussianKnots
