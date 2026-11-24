import QuadraticNorms.MinusEightLattice
import QuadraticNorms.NormRepresentationBounds

open scoped BigOperators

namespace QuadraticNorms

def MinusEightSplitPrime (p : ℕ) : Prop := p.Prime ∧ (p % 8 = 1 ∨ p % 8 = 3)

instance : DecidablePred MinusEightSplitPrime := fun _ => inferInstanceAs (Decidable (_ ∧ _))

/-- The divisor product over the split primes of the norm. -/
def minusEightSplitDivisors (n : ℕ) : ℕ :=
  ∏ p ∈ n.primeFactors.filter MinusEightSplitPrime, (n.factorization p + 1)

/-- The largest split divisor product at a positive norm up to n. -/
def minusEightDivisorMaximum (n : ℕ) : ℕ :=
  (Finset.Icc 1 n).sup minusEightSplitDivisors

theorem minusEight_split_divisors_le_maximum {a n : ℕ} (ha : 0 < a) (han : a ≤ n) :
    minusEightSplitDivisors a ≤ minusEightDivisorMaximum n :=
  Finset.le_sup (f := minusEightSplitDivisors) (Finset.mem_Icc.mpr ⟨ha, han⟩)

private theorem minusEight_split_ne_two {p : ℕ} (hp : MinusEightSplitPrime p) : p ≠ 2 := by
  rintro rfl
  exact (by decide : ¬ (2 % 8 = 1 ∨ 2 % 8 = 3)) hp.2

private theorem minusEight_root_at_two_le_one : normRootCount 0 2 2 ≤ 1 := by
  rw [norm_root_count_apply]
  apply Finset.card_le_one.mpr
  have hex : ∀ t : ZMod 2, t ^ 2 - (0 : ZMod 2) * t + 2 = 0 → t = 0 := by decide
  exact fun x hx y hy => (hex x (mem_normRootResidues.mp hx)).trans
    (hex y (mem_normRootResidues.mp hy)).symm

private theorem minusEight_root_at_four_empty : ¬ HasNormRoot 0 2 4 := by
  unfold HasNormRoot
  decide

private theorem minusEight_nonsplit_root_eq_zero {p k : ℕ} (hp : p.Prime)
    (hp₂ : p ≠ 2) (hs : ¬ MinusEightSplitPrime p) (hk : 0 < k) :
    normRootCount 0 2 (p ^ k) = 0 := by
  apply norm_root_count_of_no_root (pow_pos hp.pos _)
  intro h
  letI : Fact p.Prime := ⟨hp⟩
  obtain ⟨t, ht⟩ := h.of_dvd (dvd_pow_self p hk.ne')
  have hsq : t ^ 2 = (-2 : ZMod p) := by
    simp only [Int.cast_zero, zero_mul, sub_zero, Int.cast_ofNat] at ht
    exact eq_neg_of_add_eq_zero_left ht
  exact hs ⟨hp, (ZMod.exists_sq_eq_neg_two_iff hp₂).mp ⟨t, by rw [← pow_two, hsq]⟩⟩

private theorem minusEight_split_root_le_two {p k : ℕ} (hp : MinusEightSplitPrime p) (hk : 0 < k) :
    normRootCount 0 2 (p ^ k) ≤ 2 := by
  apply norm_root_count_unramified_le_two 0 2 hp.1 _ hk
  intro hd
  have hh : p ∣ 2 ^ 3 := Int.natCast_dvd_natCast.mp (by norm_num at hd ⊢; exact hd)
  exact minusEight_split_ne_two hp ((Nat.dvd_prime Nat.prime_two).mp
    (hp.1.dvd_of_dvd_pow hh) |>.resolve_left hp.1.ne_one)

private theorem minusEight_norm_envelope_prime_power_le {p : ℕ} (hp : p.Prime) (k : ℕ) :
    normDivisorEnvelope (normRootCount 0 2) (p ^ k) ≤ if MinusEightSplitPrime p then k + 1 else 1 := by
  induction k using Nat.twoStepInduction with
  | zero =>
    rw [norm_divisor_envelope_prime_power hp]
    simp
  | one =>
    rw [norm_divisor_envelope_prime_power hp]
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.sub_zero, Nat.sub_self,
      show ¬ 1 % 2 = 0 by decide, Nat.zero_mod, if_true, if_false, zero_add, pow_one]
    by_cases hs : MinusEightSplitPrime p
    · simp only [hs, if_true]
      simpa only [pow_one] using minusEight_split_root_le_two hs (k := 1) (by decide)
    · rw [if_neg hs]
      by_cases h : p = 2
      · simpa only [h] using minusEight_root_at_two_le_one
      · have h₀ := minusEight_nonsplit_root_eq_zero hp h hs (k := 1) (by decide)
        simpa only [pow_one] using h₀.le.trans (by decide : 0 ≤ 1)
  | more k ih _ =>
    rw [norm_divisor_envelope_prime_power_step hp]
    by_cases hs : MinusEightSplitPrime p
    · rw [if_pos hs] at ih ⊢
      have h := minusEight_split_root_le_two hs (show 0 < k + 2 by omega)
      omega
    · rw [if_neg hs] at ih ⊢
      have h : normRootCount 0 2 (p ^ (k + 2)) = 0 := by
        by_cases h₂ : p = 2
        · subst p
          exact norm_root_count_of_no_root (by positivity) (fun h =>
            minusEight_root_at_four_empty (h.of_dvd (pow_dvd_pow 2 (by omega : 2 ≤ k + 2))))
        · exact minusEight_nonsplit_root_eq_zero hp h₂ hs (by omega)
      simpa only [h, add_zero] using ih

private theorem minusEight_norm_envelope_le_split_divisors {n : ℕ} (hn : 0 < n) :
    normDivisorEnvelope (normRootCount 0 2) n ≤ minusEightSplitDivisors n := by
  have hm : (normDivisorEnvelope (normRootCount 0 2)).IsMultiplicative :=
    squareIndicator_multiplicative.mul (norm_root_count_multiplicative 0 2)
  rw [hm.multiplicative_factorization _ hn.ne', Finsupp.prod, Nat.support_factorization]
  let f (p : ℕ) := if MinusEightSplitPrime p then n.factorization p + 1 else 1
  calc
    _ ≤ ∏ p ∈ n.primeFactors, f p := Finset.prod_le_prod (fun _ _ => zero_le)
      (fun _ hp => minusEight_norm_envelope_prime_power_le (Nat.prime_of_mem_primeFactors hp) _)
    _ = ∏ p ∈ n.primeFactors.filter MinusEightSplitPrime, f p := by
      symm
      exact Finset.prod_filter_of_ne (fun p _ hp => by by_contra h; exact hp (if_neg h))
    _ = minusEightSplitDivisors n := Finset.prod_congr rfl (fun _ hp => if_pos (Finset.mem_filter.mp hp).2)

/-- The signed representation bound r_Q(n) ≤ 2 times the split divisor product. -/
theorem minusEight_norm_fiber_card_le {n : ℕ} (hn : 0 < n) :
    (minusEightNormFiber n).card ≤ 2 * minusEightSplitDivisors n := by
  have h := binary_norm_fiber_card_le (b := 0) (c := 2) (by decide) minusEightNormFiber
    (fun n v => by
      rw [mem_minusEight_norm_fiber, ← Int.natCast_inj, minusEight_norm_cast]
      simp [form]) hn
  exact h.trans (Nat.mul_le_mul_left 2 (minusEight_norm_envelope_le_split_divisors hn))

end QuadraticNorms
