import GaussianKnots.Projections

/-!
# Generic converse for the segment--triangle circuit test

The converse needs only the usual rank/genericity fact that the space of
affine dependences of five points in three-space is one-dimensional.  That
fact is exposed as a hypothesis; all sign algebra is checked here.
-/

namespace GaussianKnots

def OpenSegmentTriangle (a b c p q : Vec3) : Prop :=
  ∃ α β γ u v : ℝ,
    0 < α ∧ 0 < β ∧ 0 < γ ∧ 0 < u ∧ 0 < v ∧
    α + β + γ = 1 ∧ u + v = 1 ∧
    α • a + β • b + γ • c = u • p + v • q

def CircuitSignSplit (a b c p q : Vec3) : Prop :=
  let C := circuit5 a b c p q
  (0 < C 0 ∧ 0 < C 1 ∧ 0 < C 2 ∧ C 3 < 0 ∧ C 4 < 0) ∨
  (C 0 < 0 ∧ C 1 < 0 ∧ C 2 < 0 ∧ 0 < C 3 ∧ 0 < C 4)

theorem negative_circuit5_signs_give_open_segment_triangle
    (a b c p q : Vec3)
    (h0 : circuit5 a b c p q 0 < 0)
    (h1 : circuit5 a b c p q 1 < 0)
    (h2 : circuit5 a b c p q 2 < 0)
    (h3 : 0 < circuit5 a b c p q 3)
    (h4 : 0 < circuit5 a b c p q 4) :
    OpenSegmentTriangle a b c p q := by
  let C := circuit5 a b c p q
  let z : Fin 5 → ℝ := fun i ↦ -C i
  have hzsum : z 0 + z 1 + z 2 + z 3 + z 4 = 0 := by
    dsimp [z, C]
    linarith [circuit5_sum a b c p q]
  have hzvec :
      (z 0) • a + (z 1) • b + (z 2) • c + (z 3) • p + (z 4) • q = 0 := by
    have hc := circuit5_affine_dependence a b c p q
    funext i
    have hci := congrFun hc i
    dsimp [z, C]
    change circuit5 a b c p q 0 * a i + circuit5 a b c p q 1 * b i +
      circuit5 a b c p q 2 * c i + circuit5 a b c p q 3 * p i +
      circuit5 a b c p q 4 * q i = 0 at hci
    change (-circuit5 a b c p q 0) * a i + (-circuit5 a b c p q 1) * b i +
      (-circuit5 a b c p q 2) * c i + (-circuit5 a b c p q 3) * p i +
      (-circuit5 a b c p q 4) * q i = 0
    linarith
  apply positive_circuit_gives_open_segment_triangle a b c p q z hzsum hzvec
  · dsimp [z, C]
    linarith
  · dsimp [z, C]
    linarith
  · dsimp [z, C]
    linarith
  · dsimp [z, C]
    linarith
  · dsimp [z, C]
    linarith

/-- Open intersection forces the `3 versus 2` cofactor sign split, assuming
the generic one-dimensionality of affine dependences. -/
theorem open_segment_triangle_implies_circuitSignSplit_of_unique
    (a b c p q : Vec3)
    (hunique : ∀ z : Fin 5 → ℝ,
      z 0 + z 1 + z 2 + z 3 + z 4 = 0 →
      (z 0) • a + (z 1) • b + (z 2) • c + (z 3) • p + (z 4) • q = 0 →
      ∃ r : ℝ, z = r • circuit5 a b c p q)
    (hopen : OpenSegmentTriangle a b c p q) :
    CircuitSignSplit a b c p q := by
  rcases hopen with ⟨α, β, γ, u, v, hα, hβ, hγ, hu, hv, habc, huv, hinter⟩
  let z : Fin 5 → ℝ := ![α, β, γ, -u, -v]
  have hzsum : z 0 + z 1 + z 2 + z 3 + z 4 = 0 := by
    change α + β + γ + -u + -v = 0
    linarith
  have hzvec :
      (z 0) • a + (z 1) • b + (z 2) • c + (z 3) • p + (z 4) • q = 0 := by
    funext i
    have hi := congrFun hinter i
    change α * a i + β * b i + γ * c i = u * p i + v * q i at hi
    change α * a i + β * b i + γ * c i + (-u) * p i + (-v) * q i = 0
    linarith
  obtain ⟨r, hr⟩ := hunique z hzsum hzvec
  have e0 : α = r * circuit5 a b c p q 0 := by
    have hi := congrFun hr 0
    simpa [z] using hi
  have e1 : β = r * circuit5 a b c p q 1 := by
    have hi := congrFun hr 1
    simpa [z] using hi
  have e2 : γ = r * circuit5 a b c p q 2 := by
    have hi := congrFun hr 2
    simpa [z] using hi
  have e3 : -u = r * circuit5 a b c p q 3 := by
    have hi := congrFun hr 3
    simpa [z] using hi
  have e4 : -v = r * circuit5 a b c p q 4 := by
    have hi := congrFun hr 4
    simpa [z] using hi
  have hr0 : r ≠ 0 := by
    intro hre
    rw [hre, zero_mul] at e0
    linarith
  have c0eq : circuit5 a b c p q 0 = α / r := by
    apply (eq_div_iff hr0).2
    nlinarith [e0]
  have c1eq : circuit5 a b c p q 1 = β / r := by
    apply (eq_div_iff hr0).2
    nlinarith [e1]
  have c2eq : circuit5 a b c p q 2 = γ / r := by
    apply (eq_div_iff hr0).2
    nlinarith [e2]
  have c3eq : circuit5 a b c p q 3 = (-u) / r := by
    apply (eq_div_iff hr0).2
    nlinarith [e3]
  have c4eq : circuit5 a b c p q 4 = (-v) / r := by
    apply (eq_div_iff hr0).2
    nlinarith [e4]
  unfold CircuitSignSplit
  rcases lt_trichotomy r 0 with hrneg | hre | hrpos
  · right
    exact ⟨c0eq ▸ div_neg_of_pos_of_neg hα hrneg,
      c1eq ▸ div_neg_of_pos_of_neg hβ hrneg,
      c2eq ▸ div_neg_of_pos_of_neg hγ hrneg,
      c3eq ▸ div_pos_of_neg_of_neg (neg_neg_of_pos hu) hrneg,
      c4eq ▸ div_pos_of_neg_of_neg (neg_neg_of_pos hv) hrneg⟩
  · exact (hr0 hre).elim
  · left
    exact ⟨c0eq ▸ div_pos hα hrpos,
      c1eq ▸ div_pos hβ hrpos,
      c2eq ▸ div_pos hγ hrpos,
      c3eq ▸ div_neg_of_neg_of_pos (neg_neg_of_pos hu) hrpos,
      c4eq ▸ div_neg_of_neg_of_pos (neg_neg_of_pos hv) hrpos⟩

/-- Under the explicit generic uniqueness hypothesis, the circuit sign test
is equivalent to open segment--open triangle intersection. -/
theorem circuitSignSplit_iff_openSegmentTriangle_of_unique
    (a b c p q : Vec3)
    (hunique : ∀ z : Fin 5 → ℝ,
      z 0 + z 1 + z 2 + z 3 + z 4 = 0 →
      (z 0) • a + (z 1) • b + (z 2) • c + (z 3) • p + (z 4) • q = 0 →
      ∃ r : ℝ, z = r • circuit5 a b c p q) :
    CircuitSignSplit a b c p q ↔ OpenSegmentTriangle a b c p q := by
  constructor
  · intro hsplit
    rcases hsplit with hpos | hneg
    · exact circuit5_signs_give_open_segment_triangle a b c p q
        hpos.1 hpos.2.1 hpos.2.2.1 hpos.2.2.2.1 hpos.2.2.2.2
    · exact negative_circuit5_signs_give_open_segment_triangle a b c p q
        hneg.1 hneg.2.1 hneg.2.2.1 hneg.2.2.2.1 hneg.2.2.2.2
  · exact open_segment_triangle_implies_circuitSignSplit_of_unique a b c p q hunique

end GaussianKnots
