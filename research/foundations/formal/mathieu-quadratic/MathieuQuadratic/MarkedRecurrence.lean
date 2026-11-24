import Mathlib.Algebra.Polynomial.Roots
import Mathlib.Tactic

/-!
# A finite eventual-zero test from a recurrence

The leading coefficient propagates a zero window forward. The
trailing coefficient propagates any eventual zero sequence backward
to the first index where the recurrence is regular.
-/

open scoped BigOperators

namespace MathieuQuadratic

variable {K : Type*} [Field K]

/-- Once both endpoint coefficients are nonzero, a recurrence of
order d gives a test using exactly the d terms in the initial window. -/
theorem eventual_zero_iff_recurrence_window (a : ℕ → K) (c : ℕ → ℕ → K)
    (d m₀ : ℕ)
    (recurrence : ∀ m ≥ m₀, ∑ j ∈ Finset.range (d + 1), c m j * a (m + j) = 0)
    (leading : ∀ m ≥ m₀, c m d ≠ 0) (trailing : ∀ m ≥ m₀, c m 0 ≠ 0) :
    (∃ start : ℕ, ∀ m ≥ start, a m = 0) ↔ ∀ j < d, a (m₀ + j) = 0 := by
  constructor
  · rintro ⟨start, hstart⟩
    have hback (n : ℕ) : ∀ m ≥ m₀, start ≤ m + n → a m = 0 := by
      induction n with
      | zero =>
        intro m _ hm
        exact hstart m (by omega)
      | succ n ih =>
        intro m hm hs
        have he := recurrence m hm
        rw [Finset.sum_range_succ', add_zero] at he
        have hz : ∑ j ∈ Finset.range d, c m (j + 1) * a (m + (j + 1)) = 0 := by
          apply Finset.sum_eq_zero
          intro j _
          rw [ih (m + (j + 1)) (by omega) (by omega), mul_zero]
        rw [hz, zero_add] at he
        exact (mul_eq_zero.mp he).resolve_left (trailing m hm)
    intro j _
    exact hback start (m₀ + j) (by omega) (by omega)
  · intro hwindow
    refine ⟨m₀, fun m hm ↦ ?_⟩
    have hforward (n : ℕ) : a (m₀ + n) = 0 := by
      induction n using Nat.strong_induction_on with
      | h n ih =>
        by_cases hn : n < d
        · exact hwindow n hn
        let k := n - d
        have hk : k + d = n := Nat.sub_add_cancel (by omega)
        have he := recurrence (m₀ + k) (by omega)
        rw [Finset.sum_range_succ] at he
        have hz : ∑ j ∈ Finset.range d, c (m₀ + k) j * a (m₀ + k + j) = 0 := by
          apply Finset.sum_eq_zero
          intro j hj
          have hjd := Finset.mem_range.mp hj
          rw [add_assoc, ih (k + j) (by omega), mul_zero]
        rw [hz, zero_add, add_assoc, hk] at he
        exact (mul_eq_zero.mp he).resolve_left (leading _ (by omega))
    simpa only [Nat.add_sub_of_le hm] using hforward (m - m₀)

variable [CharZero K]

private theorem polynomial_eventually_nonzero (p : Polynomial K) (hp : p ≠ 0) :
    ∃ r : ℕ, ∀ m ≥ r, p.eval (m : K) ≠ 0 := by
  have hset := (Polynomial.finite_setOf_isRoot hp).preimage (Nat.cast_injective (R := K)).injOn
  obtain ⟨r, hr⟩ := hset.bddAbove
  refine ⟨r + 1, fun m hm hz ↦ ?_⟩
  have hle : m ≤ r := hr (show (m : K) ∈ {x | p.IsRoot x} from hz)
  omega

/-- The finite window in the draft exists for every eventual
polynomial-coefficient recurrence with nonzero endpoint polynomials. -/
theorem eventual_zero_iff_polynomial_recurrence_window (a : ℕ → K)
    (p : ℕ → Polynomial K) (d m₁ : ℕ)
    (recurrence : ∀ m ≥ m₁, ∑ j ∈ Finset.range (d + 1), (p j).eval (m : K) * a (m + j) = 0)
    (leading : p d ≠ 0) (trailing : p 0 ≠ 0) :
    ∃ m₀ ≥ m₁, (∃ start : ℕ, ∀ m ≥ start, a m = 0) ↔ ∀ j < d, a (m₀ + j) = 0 := by
  obtain ⟨r, hr⟩ := polynomial_eventually_nonzero _ leading
  obtain ⟨s, hs⟩ := polynomial_eventually_nonzero _ trailing
  let m₀ := max m₁ (max r s)
  have hm₁ : m₁ ≤ m₀ := le_max_left _ _
  have hmr : r ≤ m₀ := (le_max_left _ _).trans (le_max_right _ _)
  have hms : s ≤ m₀ := (le_max_right _ _).trans (le_max_right _ _)
  exact ⟨m₀, hm₁, eventual_zero_iff_recurrence_window a
    (fun m j ↦ (p j).eval (m : K)) d m₀ (fun m hm ↦ recurrence m (hm₁.trans hm))
    (fun m hm ↦ hr m (hmr.trans hm)) (fun m hm ↦ hs m (hms.trans hm))⟩

end MathieuQuadratic
