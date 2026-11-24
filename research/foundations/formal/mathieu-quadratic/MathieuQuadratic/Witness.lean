import MathieuQuadratic.Basic
import Mathlib.Algebra.Algebra.Rat
import Mathlib.RingTheory.SimpleRing.Basic
import Mathlib.Data.Nat.Choose.Sum
import Mathlib.Algebra.Order.Field.Rat

/-!
# The uniform normalizer witness from its two moments

The constants are positive rationals, and the marked binomial sum stays
nonzero after mapping to a characteristic-zero field.
-/

open scoped BigOperators

namespace MathieuQuadratic

variable {K A : Type*} [Field K] [CharZero K] [CommRing A] [Algebra K A]

/-- The phase from the uniform quadratic-complement construction. -/
def uniformPhase (M : ℕ) (P : A) : A := P ^ M * (1 + P)

theorem uniformPhase_pow (M m : ℕ) (P : A) :
    uniformPhase M P ^ m =
      ∑ r ∈ Finset.range (m + 1), (m.choose r : A) * P ^ (M * m + r) := by
  rw [uniformPhase, mul_pow, ← pow_mul, add_comm (1 : A), add_pow, Finset.mul_sum]
  refine Finset.sum_congr rfl fun r _ ↦ ?_
  rw [one_pow, mul_one, pow_add]
  ac_rfl

theorem uniformPhase_allOrderWitness_of_moments (M : ℕ)
    (trace : A →ₗ[K] K) (P H : A) (J : ℕ → ℚ)
    (hJ : ∀ k ≥ M, 0 < J k)
    (pureMoment : ∀ k ≥ M, trace (P ^ k) = 0)
    (markedMoment : ∀ k ≥ M, trace (H * P ^ k) = algebraMap ℚ K (J k)) :
    AllOrderWitness trace (uniformPhase M P) H ∧
      ∀ m : ℕ, 0 < m →
        trace (H * uniformPhase M P ^ m) =
          algebraMap ℚ K (∑ r ∈ Finset.range (m + 1), (m.choose r : ℚ) * J (M * m + r)) := by
  have bound (m r : ℕ) (hm : 0 < m) : M ≤ M * m + r := by
    calc
      M = M * 1 := (Nat.mul_one _).symm
      _ ≤ M * m := Nat.mul_le_mul_left _ hm
      _ ≤ M * m + r := Nat.le_add_right _ _
  have scale (n : ℕ) (a : A) : trace ((n : A) * a) = (n : K) * trace a := by
    rw [← map_natCast (algebraMap K A), ← Algebra.smul_def, map_smul, smul_eq_mul]
  have formula (m : ℕ) (hm : 0 < m) :
      trace (H * uniformPhase M P ^ m) =
        algebraMap ℚ K (∑ r ∈ Finset.range (m + 1), (m.choose r : ℚ) * J (M * m + r)) := by
    rw [uniformPhase_pow, Finset.mul_sum, map_sum, map_sum]
    refine Finset.sum_congr rfl fun r _ ↦ ?_
    rw [mul_left_comm H, scale, markedMoment _ (bound m r hm), map_mul, map_natCast]
  refine ⟨⟨?_, ?_⟩, formula⟩
  · intro m hm
    rw [uniformPhase_pow, map_sum]
    exact Finset.sum_eq_zero fun r _ ↦ by rw [scale, pureMoment _ (bound m r hm), mul_zero]
  · intro m hm
    rw [formula _ hm, ne_eq, map_eq_zero_iff _ (RingHom.injective (algebraMap ℚ K))]
    apply ne_of_gt
    apply Finset.sum_pos
    · intro r hr
      have hr' : r ≤ m := Nat.lt_succ_iff.mp (Finset.mem_range.mp hr)
      exact mul_pos (Nat.cast_pos.mpr (Nat.choose_pos hr')) (hJ _ (bound m r hm))
    · exact ⟨0, Finset.mem_range.mpr (Nat.zero_lt_succ _)⟩

end MathieuQuadratic
