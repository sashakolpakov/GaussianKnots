import Mathlib

/-!
# Exact crossing-slot arithmetic

This file isolates the finite arithmetic used by the crossing expectation and
bounded-difference estimates in the projection manuscript.
-/

namespace GaussianKnots

/-- Number of unordered pairs of nonadjacent edges in an `N`-cycle. -/
def nonincidentEdgePairCount (N : ℕ) : ℕ := N.choose 2 - N

private theorem choose_two_ge_self {N : ℕ} (hN : 3 ≤ N) : N ≤ N.choose 2 := by
  obtain ⟨k, rfl⟩ := Nat.exists_eq_add_of_le hN
  rw [show 3 + k = (k + 2) + 1 by omega, Nat.choose_succ_succ]
  simp only [Nat.choose_one_right]
  have hp : 0 < (k + 2).choose 2 := Nat.choose_pos (by omega)
  exact Nat.add_le_add_left hp (k + 2)

private theorem twice_choose_two (N : ℕ) :
    2 * N.choose 2 = N * (N - 1) := by
  rw [Nat.choose_two_right, Nat.mul_div_cancel' (Nat.even_mul_pred_self N).two_dvd]

/-- Exact cyclic count, in an integral form that avoids rational coercions. -/
theorem twice_nonincidentEdgePairCount {N : ℕ} (hN : 3 ≤ N) :
    2 * nonincidentEdgePairCount N = N * (N - 3) := by
  have hadd : nonincidentEdgePairCount N + N = N.choose 2 := by
    exact Nat.sub_add_cancel (choose_two_ge_self hN)
  have htwice := congrArg (fun x : ℕ ↦ 2 * x) hadd
  rw [mul_add, twice_choose_two] at htwice
  have hpred : N - 1 = (N - 3) + 2 := by omega
  rw [hpred, mul_add] at htwice
  have : 2 * nonincidentEdgePairCount N + 2 * N =
      N * (N - 3) + 2 * N := by
    simpa [mul_comm, mul_left_comm, mul_assoc] using htwice
  exact Nat.add_right_cancel this

theorem nonincidentEdgePairCount_eq {N : ℕ} (hN : 3 ≤ N) :
    nonincidentEdgePairCount N = N * (N - 3) / 2 := by
  exact Nat.eq_div_of_mul_eq_right two_ne_zero
    (twice_nonincidentEdgePairCount hN)

/-- Two incident edges each have at most `N-3` nonincident partners. -/
abbrev AffectedCrossingSlots (N : ℕ) := Fin 2 × Fin (N - 3)

theorem card_affectedCrossingSlots (N : ℕ) :
    Fintype.card (AffectedCrossingSlots N) = 2 * (N - 3) := by
  simp [AffectedCrossingSlots]

end GaussianKnots
