import QuadraticNorms.EulerDivisorBounds
import QuadraticNorms.MinusEightRepresentationCounts
import QuadraticNorms.MinusEightThetaNormalization
import Mathlib.Analysis.SpecialFunctions.Pow.Real

open scoped BigOperators ArithmeticFunction.sigma

namespace QuadraticNorms

private def smallSplitPrimes : Finset ℕ := (Finset.range 64).filter MinusEightSplitPrime
private def smallInertPrimes : Finset ℕ := (Finset.range 64).filter MinusEightInertPrime
private def smallOddPrimes : Finset ℕ := (Finset.range 64).filter (fun p => p.Prime ∧ p % 2 = 1)

private theorem divisor_cap_nonneg (k last p : ℕ) : 0 ≤ boundedEulerCap (divisorFactor k) last p :=
  (one_le_euler_cap (divisorFactor k) last p (by simp [divisorFactor])).trans' (by norm_num)

private theorem sigma_ratio_cap_nonneg (last p : ℕ) : 0 ≤ boundedEulerCap sigmaRatioFactor last p :=
  (one_le_euler_cap sigmaRatioFactor last p (by simp [sigmaRatioFactor])).trans' (by norm_num)

private def split_cap_sixthBound : ℕ → ℚ
  | 3 => 192902 / 1000
  | 11 => 6025 / 1000
  | 17 => 3765 / 1000
  | 19 => 3369 / 1000
  | 41 => 1561 / 1000
  | 43 => 1489 / 1000
  | 59 => 1085 / 1000
  | _ => 1

private theorem split_cap_sixthLocal :
    ∀ p ∈ ({3, 11, 17, 19, 41, 43, 59} : Finset ℕ),
      boundedEulerCap (divisorFactor 6) 10 p ≤ split_cap_sixthBound p := by
  decide +kernel

private theorem split_cap_sixth :
    (∏ p ∈ smallSplitPrimes, boundedEulerCap (divisorFactor 6) 10 p) < (6 : ℚ) ^ 6 := by
  rw [show smallSplitPrimes = {3, 11, 17, 19, 41, 43, 59} by decide +kernel]
  refine (Finset.prod_le_prod (f := fun p => boundedEulerCap (divisorFactor 6) 10 p)
    (g := split_cap_sixthBound) (fun p _ => ?_) split_cap_sixthLocal).trans_lt ?_
  · exact divisor_cap_nonneg 6 10 p
  · norm_num [split_cap_sixthBound, Finset.prod_insert, Finset.prod_singleton]

private def split_cap_fourthBound : ℕ → ℚ
  | 3 => 9482 / 1000
  | 11 => 1455 / 1000
  | 17 => 1000 / 1000
  | 19 => 1000 / 1000
  | 41 => 1000 / 1000
  | 43 => 1000 / 1000
  | 59 => 1000 / 1000
  | _ => 1

private theorem split_cap_fourthLocal :
    ∀ p ∈ ({3, 11, 17, 19, 41, 43, 59} : Finset ℕ),
      boundedEulerCap (divisorFactor 4) 10 p ≤ split_cap_fourthBound p := by
  decide +kernel

private theorem split_cap_fourth :
    (∏ p ∈ smallSplitPrimes, boundedEulerCap (divisorFactor 4) 10 p) < (2 : ℚ) ^ 4 := by
  rw [show smallSplitPrimes = {3, 11, 17, 19, 41, 43, 59} by decide +kernel]
  refine (Finset.prod_le_prod (f := fun p => boundedEulerCap (divisorFactor 4) 10 p)
    (g := split_cap_fourthBound) (fun p _ => ?_) split_cap_fourthLocal).trans_lt ?_
  · exact divisor_cap_nonneg 4 10 p
  · norm_num [split_cap_fourthBound, Finset.prod_insert, Finset.prod_singleton]

private def inert_cap_sixthBound : ℕ → ℚ
  | 5 => 32768 / 1000
  | 7 => 14878 / 1000
  | 13 => 4924 / 1000
  | 23 => 2783 / 1000
  | 29 => 2207 / 1000
  | 31 => 2065 / 1000
  | 37 => 1730 / 1000
  | 47 => 1362 / 1000
  | 53 => 1208 / 1000
  | 61 => 1050 / 1000
  | _ => 1

private theorem inert_cap_sixthLocal :
    ∀ p ∈ ({5, 7, 13, 23, 29, 31, 37, 47, 53, 61} : Finset ℕ),
      boundedEulerCap (divisorFactor 6) 10 p ≤ inert_cap_sixthBound p := by
  decide +kernel

private theorem inert_cap_sixth :
    (∏ p ∈ smallInertPrimes, boundedEulerCap (divisorFactor 6) 10 p) < (7 : ℚ) ^ 6 := by
  rw [show smallInertPrimes = {5, 7, 13, 23, 29, 31, 37, 47, 53, 61} by decide +kernel]
  refine (Finset.prod_le_prod (f := fun p => boundedEulerCap (divisorFactor 6) 10 p)
    (g := inert_cap_sixthBound) (fun p _ => ?_) inert_cap_sixthLocal).trans_lt ?_
  · exact divisor_cap_nonneg 6 10 p
  · norm_num [inert_cap_sixthBound, Finset.prod_insert, Finset.prod_singleton]

private theorem odd_cap_sixth :
    (∏ p ∈ smallOddPrimes, boundedEulerCap (divisorFactor 6) 10 p) < (42 : ℚ) ^ 6 := by
  have heq : smallOddPrimes = smallSplitPrimes ∪ smallInertPrimes := by decide +kernel
  have hdis : Disjoint smallSplitPrimes smallInertPrimes := by decide +kernel
  rw [heq, Finset.prod_union hdis, show (42 : ℚ) ^ 6 = (6 : ℚ) ^ 6 * (7 : ℚ) ^ 6 by ring]
  apply mul_lt_mul_of_nonneg_of_pos split_cap_sixth inert_cap_sixth.le
  · exact Finset.prod_nonneg fun p _ => divisor_cap_nonneg 6 10 p
  · norm_num

private theorem minusEight_split_divisors_power_le {n k : ℕ} (hn : 0 < n) (hk : k ≤ 6) :
    (minusEightSplitDivisors n : ℚ) ^ k ≤
      (∏ p ∈ smallSplitPrimes, boundedEulerCap (divisorFactor k) 10 p) * n := by
  unfold minusEightSplitDivisors
  simp only [Nat.cast_prod, Nat.cast_add, Nat.cast_one]
  exact divisor_product_power_le hn hk (Finset.filter_subset _ _)
    (fun p hp hlt => Finset.mem_filter.mpr ⟨Finset.mem_range.mpr hlt, (Finset.mem_filter.mp hp).2⟩)

theorem minusEight_split_divisors_sixth_power_lt {n : ℕ} (hn : 0 < n) :
    (minusEightSplitDivisors n : ℚ) ^ 6 < (6 : ℚ) ^ 6 * n :=
  (minusEight_split_divisors_power_le hn (by decide)).trans_lt
    (mul_lt_mul_of_pos_right split_cap_sixth (Nat.cast_pos.mpr hn))

theorem minusEight_split_divisors_fourth_power_lt {n : ℕ} (hn : 0 < n) :
    (minusEightSplitDivisors n : ℚ) ^ 4 < (2 : ℚ) ^ 4 * n :=
  (minusEight_split_divisors_power_le hn (by decide)).trans_lt
    (mul_lt_mul_of_pos_right split_cap_fourth (Nat.cast_pos.mpr hn))

private theorem minusEight_sigma_zero_power_le {m : ℕ} (hm : 0 < m) (small : Finset ℕ)
    (hsmall : ∀ p ∈ m.primeFactors, p < 64 → p ∈ small) :
    (σ 0 m : ℚ) ^ 6 ≤ (∏ p ∈ small, boundedEulerCap (divisorFactor 6) 10 p) * m := by
  have h := ArithmeticFunction.isMultiplicative_sigma.multiplicative_factorization (σ 0) hm.ne'
  rw [Finsupp.prod, Nat.support_factorization] at h
  have hτ : σ 0 m = ∏ p ∈ m.primeFactors, (m.factorization p + 1) :=
    h.trans (Finset.prod_congr rfl (fun _ hp =>
      ArithmeticFunction.sigma_zero_apply_prime_pow (Nat.prime_of_mem_primeFactors hp)))
  rw [hτ]
  simp only [Nat.cast_prod, Nat.cast_add, Nat.cast_one]
  exact divisor_product_power_le hm (by decide) (Finset.Subset.refl _) hsmall

theorem minusEight_inert_divisor_sixth_power_lt {m : ℕ} (hm : 0 < m)
    (hI : ∀ p ∈ m.primeFactors, MinusEightInertPrime p) :
    (σ 0 m : ℚ) ^ 6 < (7 : ℚ) ^ 6 * m := by
  refine (minusEight_sigma_zero_power_le hm smallInertPrimes
    (fun p hp hlt => Finset.mem_filter.mpr ⟨Finset.mem_range.mpr hlt, hI p hp⟩)).trans_lt ?_
  exact mul_lt_mul_of_pos_right inert_cap_sixth (Nat.cast_pos.mpr hm)

theorem minusEight_odd_divisor_sixth_power_lt {m : ℕ} (hm : Odd m) :
    (σ 0 m : ℚ) ^ 6 < (42 : ℚ) ^ 6 * m := by
  have hm₀ : 0 < m := by obtain ⟨k, hk⟩ := hm; omega
  refine (minusEight_sigma_zero_power_le hm₀ smallOddPrimes ?_).trans_lt
    (mul_lt_mul_of_pos_right odd_cap_sixth (Nat.cast_pos.mpr hm₀))
  intro p hp hlt
  have hodd := hm.of_dvd_nat (Nat.dvd_of_mem_primeFactors hp)
  exact Finset.mem_filter.mpr ⟨Finset.mem_range.mpr hlt, Nat.prime_of_mem_primeFactors hp, Nat.odd_iff.mp hodd⟩

private theorem minusEight_divisor_maximum_le_rpow {n k : ℕ} (hk : k ≠ 0) {C : ℚ}
    (hC : 0 ≤ C)
    (hpow : ∀ {m}, 0 < m → (minusEightSplitDivisors m : ℚ) ^ k ≤ C ^ k * m) :
    (minusEightDivisorMaximum n : ℝ) ≤ (C : ℝ) * (n : ℝ) ^ (1 / (k : ℝ)) := by
  have hC' : 0 ≤ (C : ℝ) := by exact_mod_cast hC
  have hr : 0 ≤ (C : ℝ) * (n : ℝ) ^ (1 / (k : ℝ)) := by positivity
  rw [← Nat.le_floor_iff hr]
  apply Finset.sup_le
  intro m hm
  apply Nat.le_floor
  obtain ⟨h₁, h₂⟩ := Finset.mem_Icc.mp hm
  apply (pow_le_pow_iff_left₀ (by positivity) hr hk).mp
  rw [mul_pow, ← Real.rpow_mul_natCast (Nat.cast_nonneg n),
    div_mul_cancel₀ _ (by exact_mod_cast hk), Real.rpow_one]
  have h : (minusEightSplitDivisors m : ℝ) ^ k ≤ (C : ℝ) ^ k * m := by
    exact_mod_cast hpow h₁
  exact h.trans (mul_le_mul_of_nonneg_left (by exact_mod_cast h₂) (by positivity))

/-- The two uniform envelopes for the split divisor maximum. -/
theorem minusEight_divisor_maximum_sixth_le (n : ℕ) :
    (minusEightDivisorMaximum n : ℝ) ≤ 6 * (n : ℝ) ^ (1 / 6 : ℝ) := by
  simpa using minusEight_divisor_maximum_le_rpow (n := n) (k := 6) (by decide) (C := 6) (by norm_num)
    (fun h => (minusEight_split_divisors_sixth_power_lt h).le)

theorem minusEight_divisor_maximum_fourth_le (n : ℕ) :
    (minusEightDivisorMaximum n : ℝ) ≤ 2 * (n : ℝ) ^ (1 / 4 : ℝ) := by
  simpa using minusEight_divisor_maximum_le_rpow (n := n) (k := 4) (by decide) (C := 2) (by norm_num)
    (fun h => (minusEight_split_divisors_fourth_power_lt h).le)

private def minusEight_inert_ratio_capBound : ℕ → ℚ
  | 5 => 8673 / 1000
  | 7 => 6005 / 1000
  | 13 => 3156 / 1000
  | 23 => 2156 / 1000
  | 29 => 1801 / 1000
  | 31 => 1707 / 1000
  | 37 => 1474 / 1000
  | 47 => 1201 / 1000
  | 53 => 1080 / 1000
  | 61 => 1000 / 1000
  | _ => 1

private def sigmaRatioGeometric (p a : ℕ) : ℚ :=
  ((a : ℚ) + 1) ^ 6 * (p : ℚ) ^ (5 * a) / ((∑ j ∈ Finset.range (a + 1), p ^ j : ℕ) : ℚ) ^ 6

private theorem minusEight_inert_ratio_capLocal :
    ∀ p ∈ ({5, 7, 13, 23, 29, 31, 37, 47, 53, 61} : Finset ℕ),
      boundedEulerCap sigmaRatioGeometric 3 p ≤ minusEight_inert_ratio_capBound p := by
  decide +kernel

private theorem minusEight_inert_ratio_cap :
    (∏ p ∈ smallInertPrimes, boundedEulerCap sigmaRatioFactor 3 p) < (18 / 5 : ℚ) ^ 6 := by
  have hps : smallInertPrimes = {5, 7, 13, 23, 29, 31, 37, 47, 53, 61} := by decide +kernel
  rw [hps]
  refine (Finset.prod_le_prod (f := fun p => boundedEulerCap sigmaRatioFactor 3 p)
    (g := minusEight_inert_ratio_capBound) (fun p _ => sigma_ratio_cap_nonneg 3 p) ?_).trans_lt ?_
  · intro p hp
    have hprime := ((Finset.mem_filter.mp (hps.symm ▸ hp)).2 : MinusEightInertPrime p).1
    have he : sigmaRatioFactor p = sigmaRatioGeometric p := by
      funext a
      rw [sigmaRatioFactor, sigmaRatioGeometric, ArithmeticFunction.sigma_one_apply_prime_pow hprime]
    have hcap : boundedEulerCap sigmaRatioFactor 3 p = boundedEulerCap sigmaRatioGeometric 3 p := by
      unfold boundedEulerCap
      rw [he]
    rw [hcap]
    exact minusEight_inert_ratio_capLocal p hp
  · norm_num [minusEight_inert_ratio_capBound, Finset.prod_insert, Finset.prod_singleton]

/-- The sharper Euler bound on the odd inert-supported shell error. -/
theorem minusEight_inert_shell_error_power_lt {m : ℕ} (hm : 0 < m)
    (hI : ∀ p ∈ m.primeFactors, MinusEightInertPrime p) :
    minusEightShellError m ^ 6 * (m : ℝ) ^ 2 < (18 / 5 : ℝ) ^ 6 := by
  have h := euler_product_le_small_primes sigmaRatioFactor 3 (small := smallInertPrimes) m.factorization
    (by intro p; simp [sigmaRatioFactor])
    (fun p hp hlt => Finset.mem_filter.mpr ⟨Finset.mem_range.mpr hlt, hI p hp⟩)
    (fun p hp => by
      unfold boundedEulerCap
      by_cases hlt : p < 64
      · rw [if_pos hlt]
        apply factor_le_finite_max
        intro k hk
        have hkQ : (3 : ℚ) ≤ k := by exact_mod_cast hk
        have hp₅ : (5 : ℚ) ≤ p := by
          have hmod := Nat.mod_le p 8
          exact_mod_cast (show 5 ≤ p by rcases (hI p hp).2 with h | h <;> omega)
        exact sigma_ratio_factor_step (hI p hp).1 (A := 5 / 4) (by linarith)
          ((by norm_num : (5 / 4 : ℚ) ^ 6 ≤ 5).trans hp₅)
      · rw [if_neg hlt]
        exact sigma_ratio_factor_large (hI p hp).1 (Nat.le_of_not_gt hlt))
    (fun _ _ => by unfold sigmaRatioFactor; positivity)
  have hcap := h.trans_lt minusEight_inert_ratio_cap
  rw [sigma_ratio_factor_product_eq hm] at hcap
  have hreal := (Rat.cast_lt (K := ℝ)).mpr hcap
  push_cast at hreal
  have hroot : Real.sqrt (m : ℝ) ^ 6 = (m : ℝ) ^ 3 := by
    rw [show (6 : ℕ) = 2 * 3 by decide, pow_mul, Real.sq_sqrt (by positivity)]
  unfold minusEightShellError
  rw [div_pow, mul_pow, hroot]
  rw [div_mul_eq_mul_div, mul_assoc, ← pow_add]
  exact hreal

end QuadraticNorms
