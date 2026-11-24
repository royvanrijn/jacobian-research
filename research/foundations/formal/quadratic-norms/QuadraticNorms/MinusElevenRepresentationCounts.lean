import QuadraticNorms.MinusElevenBadCounts
import QuadraticNorms.MinusElevenRootCounts

open scoped BigOperators

namespace QuadraticNorms

private noncomputable def minusElevenNormEnvelope : ArithmeticFunction ℕ :=
  normDivisorEnvelope minusElevenRootCount

private theorem norm_envelope_prime_power_le {p : ℕ} (hp : p.Prime) (k : ℕ) :
    minusElevenNormEnvelope (p ^ k) ≤ if MinusElevenSplitPrime p then k + 1 else 1 := by
  have h11 (hs : MinusElevenSplitPrime p) : p ≠ 11 := by
    rintro rfl
    exact (by decide : ¬ 11 % 11 ∈ [1, 3, 4, 5, 9]) hs.2
  induction k using Nat.twoStepInduction with
  | zero =>
    rw [minusElevenNormEnvelope, norm_divisor_envelope_prime_power hp]
    simp
  | one =>
    rw [minusElevenNormEnvelope, norm_divisor_envelope_prime_power hp]
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.sub_zero, Nat.sub_self,
      show ¬ 1 % 2 = 0 by decide, Nat.zero_mod, if_true, if_false, zero_add, pow_one]
    by_cases hs : MinusElevenSplitPrime p
    · simpa only [pow_one, hs, if_true] using minusEleven_unramified_root_count_le_two hp (h11 hs) (by decide)
        (k := 1)
    · rw [if_neg hs]
      by_cases h : p = 11
      · simpa only [h] using minusEleven_ramified_root_count_le_one
      · have h0 := minusEleven_nonsplit_root_count_eq_zero hp h hs (k := 1) (by decide)
        simpa only [pow_one, h0] using (show minusElevenRootCount (p ^ 1) ≤ 1 by rw [h0]; decide)
  | more k ih _ =>
    rw [minusElevenNormEnvelope, norm_divisor_envelope_prime_power_step hp]
    change minusElevenNormEnvelope (p ^ k) + _ ≤ _
    by_cases hs : MinusElevenSplitPrime p
    · rw [if_pos hs] at ih ⊢
      have h := minusEleven_unramified_root_count_le_two hp (h11 hs) (show 0 < k + 2 by omega)
      omega
    · rw [if_neg hs] at ih ⊢
      have h : minusElevenRootCount (p ^ (k + 2)) = 0 := by
        by_cases hp11 : p = 11
        · subst p
          exact minusEleven_ramified_root_count_eq_zero (by omega)
        · exact minusEleven_nonsplit_root_count_eq_zero hp hp11 hs (by omega)
      simpa only [h, add_zero] using ih

private theorem norm_envelope_le_split_divisors {n : ℕ} (hn : 0 < n) :
    minusElevenNormEnvelope n ≤ minusElevenSplitDivisors n := by
  have hm : minusElevenNormEnvelope.IsMultiplicative :=
    squareIndicator_multiplicative.mul minusEleven_root_count_multiplicative
  rw [hm.multiplicative_factorization _ hn.ne', Finsupp.prod, Nat.support_factorization]
  let f (p : ℕ) := if MinusElevenSplitPrime p then n.factorization p + 1 else 1
  calc
    _ ≤ ∏ p ∈ n.primeFactors, f p :=
      Finset.prod_le_prod (fun _ _ => zero_le)
        (fun _ hp => norm_envelope_prime_power_le (Nat.prime_of_mem_primeFactors hp) _)
    _ = ∏ p ∈ n.primeFactors.filter MinusElevenSplitPrime, f p := by
      symm
      exact Finset.prod_filter_of_ne (fun p _ hp => by by_contra h; exact hp (if_neg h))
    _ = minusElevenSplitDivisors n := Finset.prod_congr rfl (fun p hp => if_pos (Finset.mem_filter.mp hp).2)

private theorem minusEleven_norm_fiber_card_le_envelope {n : ℕ} (hn : 0 < n) :
    (minusElevenNormFiber n).card ≤ 2 * minusElevenNormEnvelope n :=
  binary_norm_fiber_card_le (b := 1) (c := 3) (by decide) minusElevenNormFiber
    (fun n v => by
      rw [mem_minusEleven_norm_fiber, ← Int.natCast_inj, minusEleven_norm_cast]) hn

/-- The signed representation count is bounded by twice the product of
one plus the exponent at every split prime. -/
theorem minusEleven_norm_fiber_card_le {n : ℕ} (hn : 0 < n) :
    (minusElevenNormFiber n).card ≤ 2 * minusElevenSplitDivisors n :=
  (minusEleven_norm_fiber_card_le_envelope hn).trans
    (Nat.mul_le_mul_left 2 (norm_envelope_le_split_divisors hn))

/-- The large inert-prime bound in the −11 sieve, with both counting
inputs proved for the signed-coordinate shell. -/
theorem minusEleven_inert_bad_shell_card_le {n p : ℕ} (hp : MinusElevenInertPrime p) :
    ((minusElevenBadShell n p).card : ℝ) ≤
      (48 / 5) * minusElevenDivisorMaximum n * n / (p : ℝ) ^ 2 :=
  minusEleven_inert_large_bad_count hp (fun _ ha han =>
    (minusEleven_norm_fiber_card_le ha).trans (Nat.mul_le_mul_left 2
      (minusEleven_split_divisors_le_maximum ha han)))

/-- The unrestricted −11 classification from the three spectral shell
bounds on the reduced targets above the finite scan. -/
theorem P11_sum_of_spectral_shell_bounds
    (hShell : ∀ t m : ℕ, 0 < t → Odd m → 2 ^ 25 < 2 ^ t * m → ¬ 11 ∣ 2 ^ t * m →
      IsReduced C11 (2 ^ t * m) →
      let v := 2 ^ t * m
      let E := minusElevenFilteredError t m
      let G := ((minusElevenOddShell v).card : ℝ)
      (18 / 5 : ℝ) * v * (1 - E) ≤ G ∧
        (1 - E) * (minusElevenBadShell v 11).card ≤ (1 / 3 + 2 * E) * G ∧
        ∀ p, MinusElevenInertPrime p →
          (1 - E) * (minusElevenBadShell v p).card ≤ (minusElevenLocalMass p + (5 / 2) * E) * G)
    (n : ℕ) :
    InSumset P11 P11 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 121 ∣ n := by
  apply P11_sum_of_shell_bounds ?_ n
  intro t m ht hm hcut h₁₁ hred
  obtain ⟨hLower, hRamified, hSmall⟩ := hShell t m ht hm hcut h₁₁ hred
  exact ⟨hLower, hRamified, hSmall, fun _ hp => minusEleven_inert_bad_shell_card_le hp⟩

end QuadraticNorms
