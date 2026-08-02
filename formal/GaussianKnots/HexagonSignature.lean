import GaussianKnots.CircuitConverse

/-!
# Six-vertex orientation signatures and the Calvo input data

This file closes a data-flow gap in `GaussianKnots.Projections`.  A `CalvoData`
there contains six `PiercingSigns`, and each `PiercingSigns` has an
`orientedPiercing` field.  Taken in isolation, that field can be chosen
independently of the five circuit signs.  Here it is never an input: all six
records, including the oriented contribution, are derived from one typed
15-entry orientation signature on `Fin 6`.

The topological interpretation of `calvoDecision` remains external.  This
module only connects a six-vertex determinant signature to the already
checked circuit arithmetic and executable Calvo decision layer.
-/

namespace GaussianKnots

abbrev Hexagon := Fin 6 → Vec3

/-- The signs of the `15 = choose 6 4` increasing quadruples, in lexicographic
order

`0123, 0124, 0125, 0134, 0135, 0145, 0234, 0235, 0245, 0345,
 1234, 1235, 1245, 1345, 2345`.
-/
abbrev HexagonSignature := Fin 15 → OrientationSign

namespace OrientationSign

/-- The sign convention used by the finite classifier.  At zero this returns
`pos`; the geometric theorems use an explicit nonvanishing hypothesis. -/
noncomputable def ofReal (x : ℝ) : OrientationSign :=
  if x < 0 then .neg else .pos

@[simp] theorem ofReal_of_pos {x : ℝ} (hx : 0 < x) : ofReal x = .pos := by
  simp [ofReal, hx.not_gt]

@[simp] theorem ofReal_of_neg {x : ℝ} (hx : x < 0) : ofReal x = .neg := by
  simp [ofReal, hx]

theorem ofReal_neg {x : ℝ} (hx : x ≠ 0) : ofReal (-x) = (ofReal x).flip := by
  rcases lt_or_gt_of_ne hx with hxneg | hxpos
  · rw [ofReal_of_pos (neg_pos.mpr hxneg), ofReal_of_neg hxneg]
    rfl
  · rw [ofReal_of_neg (neg_lt_zero.mpr hxpos), ofReal_of_pos hxpos]
    rfl

end OrientationSign

/-- The fifteen increasing-quadruple determinants of a labeled hexagon. -/
def hexagonOrientationValues (v : Hexagon) : Fin 15 → ℝ :=
  ![orient4 (v 0) (v 1) (v 2) (v 3),
    orient4 (v 0) (v 1) (v 2) (v 4),
    orient4 (v 0) (v 1) (v 2) (v 5),
    orient4 (v 0) (v 1) (v 3) (v 4),
    orient4 (v 0) (v 1) (v 3) (v 5),
    orient4 (v 0) (v 1) (v 4) (v 5),
    orient4 (v 0) (v 2) (v 3) (v 4),
    orient4 (v 0) (v 2) (v 3) (v 5),
    orient4 (v 0) (v 2) (v 4) (v 5),
    orient4 (v 0) (v 3) (v 4) (v 5),
    orient4 (v 1) (v 2) (v 3) (v 4),
    orient4 (v 1) (v 2) (v 3) (v 5),
    orient4 (v 1) (v 2) (v 4) (v 5),
    orient4 (v 1) (v 3) (v 4) (v 5),
    orient4 (v 2) (v 3) (v 4) (v 5)]

/-- The genericity boundary needed to read all fifteen real determinants as
two-valued orientation signs. -/
def GenericHexagon (v : Hexagon) : Prop :=
  ∀ i, hexagonOrientationValues v i ≠ 0

/-- Extract the lexicographically ordered 15-sign signature of a labeled
hexagon.  The definition is total; `GenericHexagon v` records when no zero was
silently sent to `pos`. -/
noncomputable def hexagonSignature (v : Hexagon) : HexagonSignature :=
  fun i ↦ OrientationSign.ofReal (hexagonOrientationValues v i)

private theorem orient4_rotate_one (a b c d : Vec3) :
    orient4 b c d a = -orient4 a b c d := by
  simp [orient4, det3]
  ring

private theorem orient4_rotate_two (a b c d : Vec3) :
    orient4 c d a b = orient4 a b c d := by
  simp [orient4, det3]
  ring

private theorem orient4_rotate_three (a b c d : Vec3) :
    orient4 d a b c = -orient4 a b c d := by
  simp [orient4, det3]
  ring

/-- A bit-mask identifies the underlying four-element subset of `Fin 6`.
The Calvo queries below always contain four distinct indices. -/
def quadMask (a b c d : Fin 6) : ℕ :=
  2 ^ a.val + 2 ^ b.val + 2 ^ c.val + 2 ^ d.val

/-- Lexicographic index of a four-subset of `Fin 6`, encoded by `quadMask`.
The default branch is reached only for a tuple with repeated indices. -/
def quadIndex : ℕ → Fin 15
  | 15 => 0   -- 0123
  | 23 => 1   -- 0124
  | 39 => 2   -- 0125
  | 27 => 3   -- 0134
  | 43 => 4   -- 0135
  | 51 => 5   -- 0145
  | 29 => 6   -- 0234
  | 45 => 7   -- 0235
  | 53 => 8   -- 0245
  | 57 => 9   -- 0345
  | 30 => 10  -- 1234
  | 46 => 11  -- 1235
  | 54 => 12  -- 1245
  | 58 => 13  -- 1345
  | 60 => 14  -- 2345
  | _ => 0

private def inversion (a b : Fin 6) : ℕ :=
  if b < a then 1 else 0

/-- Inversion parity of an ordered quadruple. -/
def oddQuadruple (a b c d : Fin 6) : Prop :=
  (inversion a b + inversion a c + inversion a d +
    inversion b c + inversion b d + inversion c d) % 2 = 1

instance (a b c d : Fin 6) : Decidable (oddQuadruple a b c d) :=
  inferInstanceAs (Decidable (_ = _))

/-- Alternating extension of the 15 increasing-quadruple signs to an ordered
quadruple.  Repeated-index inputs use the documented default in `quadIndex`;
none occurs in the typed Calvo queries. -/
def orientationAt (s : HexagonSignature) (a b c d : Fin 6) : OrientationSign :=
  let base := s (quadIndex (quadMask a b c d))
  if oddQuadruple a b c d then base.flip else base

/-- A triangle/segment query.  The five labels are `Fin 6`; the six queries
used below contain five distinct labels by construction. -/
structure PiercingQuery where
  a : Fin 6
  b : Fin 6
  c : Fin 6
  p : Fin 6
  q : Fin 6
  deriving DecidableEq, Repr

/-- Derive every field of the old `PiercingSigns` record from one signature.
In particular, `orientedPiercing` is not supplied independently. -/
def PiercingQuery.signs (query : PiercingQuery) (s : HexagonSignature) : PiercingSigns where
  c0 := orientationAt s query.b query.c query.p query.q
  c1 := (orientationAt s query.a query.c query.p query.q).flip
  c2 := orientationAt s query.a query.b query.p query.q
  c3 := (orientationAt s query.a query.b query.c query.q).flip
  c4 := orientationAt s query.a query.b query.c query.p
  orientedPiercing := orientationAt s query.a query.b query.c query.q

/-- The same record read directly from the five real points.  This definition
is useful for stating the exact compatibility theorem for a vertex-derived
signature; unlike the old public structure constructor, it also derives the
oriented contribution rather than accepting it as a free argument. -/
noncomputable def PiercingQuery.geometricSigns (query : PiercingQuery)
    (v : Hexagon) : PiercingSigns :=
  let a := v query.a
  let b := v query.b
  let c := v query.c
  let p := v query.p
  let q := v query.q
  { c0 := OrientationSign.ofReal (circuit5 a b c p q 0)
    c1 := OrientationSign.ofReal (circuit5 a b c p q 1)
    c2 := OrientationSign.ofReal (circuit5 a b c p q 2)
    c3 := OrientationSign.ofReal (circuit5 a b c p q 3)
    c4 := OrientationSign.ofReal (circuit5 a b c p q 4)
    orientedPiercing := OrientationSign.ofReal (orient4 a b c q) }

@[simp] theorem PiercingQuery.signs_orientedPiercing (query : PiercingQuery)
    (s : HexagonSignature) :
    (query.signs s).orientedPiercing = orientationAt s query.a query.b query.c query.q := rfl

@[simp] theorem PiercingQuery.signs_c3 (query : PiercingQuery) (s : HexagonSignature) :
    (query.signs s).c3 = (orientationAt s query.a query.b query.c query.q).flip := rfl

/-- The compatibility relation missing from a freely constructed
`PiercingSigns`: the oriented contribution is forced by the same determinant
whose negated sign is the fourth circuit coefficient. -/
theorem PiercingQuery.signs_orientedPiercing_eq_flip_c3 (query : PiercingQuery)
    (s : HexagonSignature) :
    (query.signs s).orientedPiercing = (query.signs s).c3.flip := by
  simp

private def d2aQuery : PiercingQuery := ⟨0, 1, 2, 3, 4⟩
private def d2bQuery : PiercingQuery := ⟨0, 1, 2, 4, 5⟩
private def d4aQuery : PiercingQuery := ⟨2, 3, 4, 5, 0⟩
private def d4bQuery : PiercingQuery := ⟨2, 3, 4, 0, 1⟩
private def d6aQuery : PiercingQuery := ⟨4, 5, 0, 1, 2⟩
private def d6bQuery : PiercingQuery := ⟨4, 5, 0, 2, 3⟩

private theorem d2a_signs_hexagonSignature (v : Hexagon) (hv : GenericHexagon v) :
    d2aQuery.signs (hexagonSignature v) = d2aQuery.geometricSigns v := by
  have h1 :
      OrientationSign.ofReal (-orient4 (v 0) (v 1) (v 2) (v 4)) =
        (OrientationSign.ofReal (orient4 (v 0) (v 1) (v 2) (v 4))).flip :=
    OrientationSign.ofReal_neg (by simpa [hexagonOrientationValues] using hv 1)
  have h6 :
      OrientationSign.ofReal (-orient4 (v 0) (v 2) (v 3) (v 4)) =
        (OrientationSign.ofReal (orient4 (v 0) (v 2) (v 3) (v 4))).flip :=
    OrientationSign.ofReal_neg (by simpa [hexagonOrientationValues] using hv 6)
  simp [PiercingQuery.signs, PiercingQuery.geometricSigns, d2aQuery,
    orientationAt, oddQuadruple, inversion, quadMask, quadIndex,
    hexagonSignature, hexagonOrientationValues, circuit5, h1, h6]

private theorem d2b_signs_hexagonSignature (v : Hexagon) (hv : GenericHexagon v) :
    d2bQuery.signs (hexagonSignature v) = d2bQuery.geometricSigns v := by
  have h2 :
      OrientationSign.ofReal (-orient4 (v 0) (v 1) (v 2) (v 5)) =
        (OrientationSign.ofReal (orient4 (v 0) (v 1) (v 2) (v 5))).flip :=
    OrientationSign.ofReal_neg (by simpa [hexagonOrientationValues] using hv 2)
  have h8 :
      OrientationSign.ofReal (-orient4 (v 0) (v 2) (v 4) (v 5)) =
        (OrientationSign.ofReal (orient4 (v 0) (v 2) (v 4) (v 5))).flip :=
    OrientationSign.ofReal_neg (by simpa [hexagonOrientationValues] using hv 8)
  simp [PiercingQuery.signs, PiercingQuery.geometricSigns, d2bQuery,
    orientationAt, oddQuadruple, inversion, quadMask, quadIndex,
    hexagonSignature, hexagonOrientationValues, circuit5, h2, h8]

private theorem d4a_signs_hexagonSignature (v : Hexagon) (hv : GenericHexagon v) :
    d4aQuery.signs (hexagonSignature v) = d4aQuery.geometricSigns v := by
  have e0 : orient4 (v 3) (v 4) (v 5) (v 0) =
      -orient4 (v 0) (v 3) (v 4) (v 5) := orient4_rotate_one _ _ _ _
  have e1 : orient4 (v 2) (v 4) (v 5) (v 0) =
      -orient4 (v 0) (v 2) (v 4) (v 5) := orient4_rotate_one _ _ _ _
  have e2 : orient4 (v 2) (v 3) (v 5) (v 0) =
      -orient4 (v 0) (v 2) (v 3) (v 5) := orient4_rotate_one _ _ _ _
  have e3 : orient4 (v 2) (v 3) (v 4) (v 0) =
      -orient4 (v 0) (v 2) (v 3) (v 4) := orient4_rotate_one _ _ _ _
  have h6 :
      OrientationSign.ofReal (-orient4 (v 0) (v 2) (v 3) (v 4)) =
        (OrientationSign.ofReal (orient4 (v 0) (v 2) (v 3) (v 4))).flip :=
    OrientationSign.ofReal_neg (by simpa [hexagonOrientationValues] using hv 6)
  have h7 :
      OrientationSign.ofReal (-orient4 (v 0) (v 2) (v 3) (v 5)) =
        (OrientationSign.ofReal (orient4 (v 0) (v 2) (v 3) (v 5))).flip :=
    OrientationSign.ofReal_neg (by simpa [hexagonOrientationValues] using hv 7)
  have h9 :
      OrientationSign.ofReal (-orient4 (v 0) (v 3) (v 4) (v 5)) =
        (OrientationSign.ofReal (orient4 (v 0) (v 3) (v 4) (v 5))).flip :=
    OrientationSign.ofReal_neg (by simpa [hexagonOrientationValues] using hv 9)
  simp [PiercingQuery.signs, PiercingQuery.geometricSigns, d4aQuery,
    orientationAt, oddQuadruple, inversion, quadMask, quadIndex,
    hexagonSignature, hexagonOrientationValues, circuit5,
    e0, e1, e2, e3, h6, h7, h9]

private theorem d4b_signs_hexagonSignature (v : Hexagon) (hv : GenericHexagon v) :
    d4bQuery.signs (hexagonSignature v) = d4bQuery.geometricSigns v := by
  have e0 : orient4 (v 3) (v 4) (v 0) (v 1) =
      orient4 (v 0) (v 1) (v 3) (v 4) := orient4_rotate_two _ _ _ _
  have e1 : orient4 (v 2) (v 4) (v 0) (v 1) =
      orient4 (v 0) (v 1) (v 2) (v 4) := orient4_rotate_two _ _ _ _
  have e2 : orient4 (v 2) (v 3) (v 0) (v 1) =
      orient4 (v 0) (v 1) (v 2) (v 3) := orient4_rotate_two _ _ _ _
  have e3 : orient4 (v 2) (v 3) (v 4) (v 1) =
      -orient4 (v 1) (v 2) (v 3) (v 4) := orient4_rotate_one _ _ _ _
  have e4 : orient4 (v 2) (v 3) (v 4) (v 0) =
      -orient4 (v 0) (v 2) (v 3) (v 4) := orient4_rotate_one _ _ _ _
  have h1 :
      OrientationSign.ofReal (-orient4 (v 0) (v 1) (v 2) (v 4)) =
        (OrientationSign.ofReal (orient4 (v 0) (v 1) (v 2) (v 4))).flip :=
    OrientationSign.ofReal_neg (by simpa [hexagonOrientationValues] using hv 1)
  have h6 :
      OrientationSign.ofReal (-orient4 (v 0) (v 2) (v 3) (v 4)) =
        (OrientationSign.ofReal (orient4 (v 0) (v 2) (v 3) (v 4))).flip :=
    OrientationSign.ofReal_neg (by simpa [hexagonOrientationValues] using hv 6)
  have h10 :
      OrientationSign.ofReal (-orient4 (v 1) (v 2) (v 3) (v 4)) =
        (OrientationSign.ofReal (orient4 (v 1) (v 2) (v 3) (v 4))).flip :=
    OrientationSign.ofReal_neg (by simpa [hexagonOrientationValues] using hv 10)
  simp [PiercingQuery.signs, PiercingQuery.geometricSigns, d4bQuery,
    orientationAt, oddQuadruple, inversion, quadMask, quadIndex,
    hexagonSignature, hexagonOrientationValues, circuit5,
    e0, e1, e2, e3, e4, h1, h6, h10]

private theorem d6a_signs_hexagonSignature (v : Hexagon) (hv : GenericHexagon v) :
    d6aQuery.signs (hexagonSignature v) = d6aQuery.geometricSigns v := by
  have e0 : orient4 (v 5) (v 0) (v 1) (v 2) =
      -orient4 (v 0) (v 1) (v 2) (v 5) := orient4_rotate_three _ _ _ _
  have e1 : orient4 (v 4) (v 0) (v 1) (v 2) =
      -orient4 (v 0) (v 1) (v 2) (v 4) := orient4_rotate_three _ _ _ _
  have e2 : orient4 (v 4) (v 5) (v 1) (v 2) =
      orient4 (v 1) (v 2) (v 4) (v 5) := orient4_rotate_two _ _ _ _
  have e3 : orient4 (v 4) (v 5) (v 0) (v 2) =
      orient4 (v 0) (v 2) (v 4) (v 5) := orient4_rotate_two _ _ _ _
  have e4 : orient4 (v 4) (v 5) (v 0) (v 1) =
      orient4 (v 0) (v 1) (v 4) (v 5) := orient4_rotate_two _ _ _ _
  have h2 :
      OrientationSign.ofReal (-orient4 (v 0) (v 1) (v 2) (v 5)) =
        (OrientationSign.ofReal (orient4 (v 0) (v 1) (v 2) (v 5))).flip :=
    OrientationSign.ofReal_neg (by simpa [hexagonOrientationValues] using hv 2)
  have h8 :
      OrientationSign.ofReal (-orient4 (v 0) (v 2) (v 4) (v 5)) =
        (OrientationSign.ofReal (orient4 (v 0) (v 2) (v 4) (v 5))).flip :=
    OrientationSign.ofReal_neg (by simpa [hexagonOrientationValues] using hv 8)
  simp [PiercingQuery.signs, PiercingQuery.geometricSigns, d6aQuery,
    orientationAt, oddQuadruple, inversion, quadMask, quadIndex,
    hexagonSignature, hexagonOrientationValues, circuit5,
    e0, e1, e2, e3, e4, h2, h8]

private theorem d6b_signs_hexagonSignature (v : Hexagon) (hv : GenericHexagon v) :
    d6bQuery.signs (hexagonSignature v) = d6bQuery.geometricSigns v := by
  have e0 : orient4 (v 5) (v 0) (v 2) (v 3) =
      -orient4 (v 0) (v 2) (v 3) (v 5) := orient4_rotate_three _ _ _ _
  have e1 : orient4 (v 4) (v 0) (v 2) (v 3) =
      -orient4 (v 0) (v 2) (v 3) (v 4) := orient4_rotate_three _ _ _ _
  have e2 : orient4 (v 4) (v 5) (v 2) (v 3) =
      orient4 (v 2) (v 3) (v 4) (v 5) := orient4_rotate_two _ _ _ _
  have e3 : orient4 (v 4) (v 5) (v 0) (v 3) =
      orient4 (v 0) (v 3) (v 4) (v 5) := orient4_rotate_two _ _ _ _
  have e4 : orient4 (v 4) (v 5) (v 0) (v 2) =
      orient4 (v 0) (v 2) (v 4) (v 5) := orient4_rotate_two _ _ _ _
  have h7 :
      OrientationSign.ofReal (-orient4 (v 0) (v 2) (v 3) (v 5)) =
        (OrientationSign.ofReal (orient4 (v 0) (v 2) (v 3) (v 5))).flip :=
    OrientationSign.ofReal_neg (by simpa [hexagonOrientationValues] using hv 7)
  have h9 :
      OrientationSign.ofReal (-orient4 (v 0) (v 3) (v 4) (v 5)) =
        (OrientationSign.ofReal (orient4 (v 0) (v 3) (v 4) (v 5))).flip :=
    OrientationSign.ofReal_neg (by simpa [hexagonOrientationValues] using hv 9)
  simp [PiercingQuery.signs, PiercingQuery.geometricSigns, d6bQuery,
    orientationAt, oddQuadruple, inversion, quadMask, quadIndex,
    hexagonSignature, hexagonOrientationValues, circuit5,
    e0, e1, e2, e3, e4, h7, h9]

/-- The exact six Calvo queries, now generated from one 15-sign table. -/
def calvoDataOfSignature (s : HexagonSignature) : CalvoData where
  d2a := d2aQuery.signs s
  d2b := d2bQuery.signs s
  d4a := d4aQuery.signs s
  d4b := d4bQuery.signs s
  d6a := d6aQuery.signs s
  d6b := d6bQuery.signs s

/-- The six records read directly from the real circuit coefficients. -/
noncomputable def geometricCalvoData (v : Hexagon) : CalvoData where
  d2a := d2aQuery.geometricSigns v
  d2b := d2bQuery.geometricSigns v
  d4a := d4aQuery.geometricSigns v
  d4b := d4bQuery.geometricSigns v
  d6a := d6aQuery.geometricSigns v
  d6b := d6bQuery.geometricSigns v

/-- The vertex-to-signature-to-Calvo-data pipeline. -/
noncomputable def calvoDataOfHexagon (v : Hexagon) : CalvoData :=
  calvoDataOfSignature (hexagonSignature v)

/-- On the explicit nonvanishing boundary, the compact 15-sign pipeline is
exactly the data obtained by signing the six real `circuit5` coefficient
vectors and their oriented determinants. -/
theorem calvoDataOfHexagon_eq_geometric (v : Hexagon) (hv : GenericHexagon v) :
    calvoDataOfHexagon v = geometricCalvoData v := by
  simp [calvoDataOfHexagon, calvoDataOfSignature, geometricCalvoData,
    d2a_signs_hexagonSignature v hv, d2b_signs_hexagonSignature v hv,
    d4a_signs_hexagonSignature v hv, d4b_signs_hexagonSignature v hv,
    d6a_signs_hexagonSignature v hv, d6b_signs_hexagonSignature v hv]

/-- Every oriented-piercing field in data produced by the signature pipeline
is constrained by its corresponding circuit sign. -/
theorem calvoDataOfSignature_compatible (s : HexagonSignature) :
    let x := calvoDataOfSignature s
    x.d2a.orientedPiercing = x.d2a.c3.flip ∧
    x.d2b.orientedPiercing = x.d2b.c3.flip ∧
    x.d4a.orientedPiercing = x.d4a.c3.flip ∧
    x.d4b.orientedPiercing = x.d4b.c3.flip ∧
    x.d6a.orientedPiercing = x.d6a.c3.flip ∧
    x.d6b.orientedPiercing = x.d6b.c3.flip := by
  simp [calvoDataOfSignature]

/-- The canonical 15-sign literal recorded in the manuscript. -/
def canonicalTrefoilSignature : HexagonSignature :=
  ![.pos, .neg, .pos, .pos, .neg,
    .pos, .pos, .neg, .neg, .pos,
    .neg, .pos, .pos, .neg, .pos]

/-- The literal signature generates exactly the previously transcribed
`canonicalTrefoilCalvoData`; it is no longer a second independent datum. -/
theorem canonicalTrefoilSignature_calvoData :
    calvoDataOfSignature canonicalTrefoilSignature = canonicalTrefoilCalvoData := by
  rfl

theorem canonicalTrefoilSignature_deltas :
    calvoDeltas (calvoDataOfSignature canonicalTrefoilSignature) = (1, 1, 1) := by
  rfl

theorem canonicalTrefoilSignature_decision :
    calvoDecision (calvoDataOfSignature canonicalTrefoilSignature) = .rightTrefoil := by
  rfl

end GaussianKnots
