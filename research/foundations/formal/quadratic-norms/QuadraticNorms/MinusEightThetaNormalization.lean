import QuadraticNorms.MinusEightShell
import QuadraticNorms.SigmaBounds
import Mathlib.NumberTheory.ArithmeticFunction.Misc

open scoped BigOperators ArithmeticFunction.sigma

namespace QuadraticNorms

/-- The convention σ(n/d) = 0 at a nonintegral index. -/
noncomputable def minusEightSigmaQuotient (n d : ℕ) : ℝ :=
  if d ∣ n then (σ 1 (n / d) : ℝ) else 0

/-- The two modular coefficient identities used for the −8 shell,
together with the even-index vanishing of the level-32 eta product. -/
structure MinusEightThetaFormula (a : ℕ → ℝ) : Prop where
  full : ∀ n : ℕ, 0 < n → ((minusEightFullShell n).card : ℝ) =
    4 * minusEightSigmaQuotient n 1 - 4 * minusEightSigmaQuotient n 2 +
      8 * minusEightSigmaQuotient n 4 - 32 * minusEightSigmaQuotient n 8
  weighted : ∀ n : ℕ, 0 < n → ((minusEightWeightedShell 4 n).card : ℝ) =
    minusEightSigmaQuotient n 1 - minusEightSigmaQuotient n 2 +
      8 * minusEightSigmaQuotient n 16 - 32 * minusEightSigmaQuotient n 32 + a n
  even_zero : ∀ n : ℕ, 2 ∣ n → a n = 0

private theorem sigma_two_power_mul_odd (e m : ℕ) (hm : Odd m) :
    (σ 1 (2 ^ e * m) : ℝ) = ((2 : ℝ) ^ (e + 1) - 1) * (σ 1 m : ℝ) := by
  have hc := (Nat.Prime.coprime_iff_not_dvd Nat.prime_two).mpr hm.not_two_dvd_nat
  rw [ArithmeticFunction.isMultiplicative_sigma.map_mul_of_coprime (hc.pow_left e),
    ArithmeticFunction.sigma_one_apply_prime_pow Nat.prime_two, Nat.cast_mul, Nat.cast_sum]
  push_cast
  have h := geom_sum₂_mul_add (1 : ℝ) 1 (e + 1)
  simp only [one_add_one_eq_two, one_pow, mul_one] at h
  rw [← eq_sub_iff_add_eq] at h
  rw [h]

private theorem minusEight_pow_two_dvd_iff (e k m : ℕ) (hm : Odd m) :
    2 ^ k ∣ 2 ^ e * m ↔ k ≤ e := by
  have hc := (Nat.Prime.coprime_iff_not_dvd Nat.prime_two).mpr hm.not_two_dvd_nat
  exact ⟨fun h => (Nat.pow_dvd_pow_iff_le_right (by decide : 1 < 2)).mp
    ((hc.pow_left k).dvd_of_dvd_mul_right h),
    fun h => (Nat.pow_dvd_pow 2 h).trans (dvd_mul_right _ _)⟩

private theorem minusEight_sigma_quotient_two (e k m : ℕ) (hm : Odd m) :
    minusEightSigmaQuotient (2 ^ e * m) (2 ^ k) =
      if k ≤ e then ((2 : ℝ) ^ (e - k + 1) - 1) * (σ 1 m : ℝ) else 0 := by
  unfold minusEightSigmaQuotient
  simp only [minusEight_pow_two_dvd_iff _ _ _ hm]
  by_cases h : k ≤ e
  · rw [if_pos h, if_pos h]
    have hdiv : (2 ^ e * m) / 2 ^ k = 2 ^ (e - k) * m := by
      have hpow : (2 : ℕ) ^ e = 2 ^ k * 2 ^ (e - k) := by rw [← pow_add, Nat.add_sub_of_le h]
      rw [hpow, mul_assoc, Nat.mul_div_cancel_left _ (by positivity)]
    rw [hdiv, sigma_two_power_mul_odd _ _ hm]
  · rw [if_neg h, if_neg h]

/-- The coefficient c_e in the full shell count 4 c_e σ(m). -/
def minusEightShellFactor (e : ℕ) : ℕ :=
  if e = 0 then 1 else if e = 1 then 2 else 6

/-- The full shell counts in the four eligible dyadic strata. -/
theorem minusEight_full_shell_dyadic {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m : ℕ} (he : e ≤ 3) (hm : Odd m) :
    ((minusEightFullShell (2 ^ e * m)).card : ℝ) =
      4 * (minusEightShellFactor e : ℝ) * (σ 1 m : ℝ) := by
  have hm₀ : 0 < m := by obtain ⟨k, hk⟩ := hm; omega
  rw [hθ.full _ (by positivity)]
  have h₀ := minusEight_sigma_quotient_two e 0 m hm
  have h₁ := minusEight_sigma_quotient_two e 1 m hm
  have h₂ := minusEight_sigma_quotient_two e 2 m hm
  have h₃ := minusEight_sigma_quotient_two e 3 m hm
  norm_num at h₀ h₁ h₂ h₃
  rw [h₀, h₁, h₂, h₃]
  interval_cases e <;> norm_num [minusEightShellFactor] <;> ring

/-- The Q ⊕ 4Q shell counts needed for dyadic inclusion-exclusion. -/
theorem minusEight_weighted_shell_dyadic {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m : ℕ} (he : e ≤ 3) (hm : Odd m) :
    ((minusEightWeightedShell 4 (2 ^ e * m)).card : ℝ) =
      (2 : ℝ) ^ e * (σ 1 m : ℝ) + if e = 0 then a m else 0 := by
  have hm₀ : 0 < m := by obtain ⟨k, hk⟩ := hm; omega
  rw [hθ.weighted _ (by positivity)]
  have h₀ := minusEight_sigma_quotient_two e 0 m hm
  have h₁ := minusEight_sigma_quotient_two e 1 m hm
  have h₄ := minusEight_sigma_quotient_two e 4 m hm
  have h₅ := minusEight_sigma_quotient_two e 5 m hm
  norm_num at h₀ h₁ h₄ h₅
  rw [h₀, h₁, h₄, h₅, if_neg (by omega : ¬ 4 ≤ e), if_neg (by omega : ¬ 5 ≤ e)]
  by_cases he₀ : e = 0
  · subst e
    norm_num
  · rw [if_neg he₀, hθ.even_zero _ ((dvd_pow_self 2 he₀).trans (dvd_mul_right _ _))]
    interval_cases e <;> norm_num <;> ring

/-- The exact G₂ shell counts in the four eligible strata. The only
cusp correction occurs in the odd stratum. -/
theorem minusEight_dyadic_shell_theta {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m : ℕ} (he : e ≤ 3) (hm : Odd m) :
    ((minusEightDyadicShell (2 ^ e * m)).card : ℝ) =
      if e = 0 then 2 * (σ 1 m : ℝ) - 2 * a m
      else if e = 1 then 4 * (σ 1 m : ℝ)
      else if e = 2 then 20 * (σ 1 m : ℝ)
      else 16 * (σ 1 m : ℝ) := by
  have h : ((minusEightFullShell (2 ^ e * m)).card : ℝ) +
      ((if 4 ∣ 2 ^ e * m then (minusEightFullShell (2 ^ e * m / 4)).card else 0 : ℕ) : ℝ) =
      (minusEightDyadicShell (2 ^ e * m)).card + 2 * (minusEightWeightedShell 4 (2 ^ e * m)).card := by
    exact_mod_cast minusEight_dyadic_shell_card (2 ^ e * m)
  have hdvd : 4 ∣ 2 ^ e * m ↔ 2 ≤ e := minusEight_pow_two_dvd_iff e 2 m hm
  simp only [Nat.cast_ite, Nat.cast_zero, hdvd, minusEight_full_shell_dyadic hθ he hm,
    minusEight_weighted_shell_dyadic hθ he hm] at h
  by_cases h₂ : 2 ≤ e
  · have hq : 2 ^ e * m / 4 = 2 ^ (e - 2) * m := by
      conv_lhs => rw [show (2 : ℕ) ^ e = 2 ^ 2 * 2 ^ (e - 2) by rw [← pow_add, Nat.add_sub_of_le h₂]]
      rw [mul_assoc]
      exact Nat.mul_div_cancel_left _ (by decide : 0 < 4)
    have he' : e - 2 ≤ 3 := by omega
    rw [if_pos h₂, hq, minusEight_full_shell_dyadic hθ he' hm] at h
    interval_cases e <;> norm_num [minusEightShellFactor] at h ⊢ <;> linarith
  · rw [if_neg h₂] at h
    interval_cases e <;> norm_num [minusEightShellFactor] at h ⊢ <;> linarith

/-- The relative odd-part cusp error used in the −8 sieve. -/
noncomputable def minusEightShellError (m : ℕ) : ℝ :=
  (σ 0 m : ℝ) * Real.sqrt m / (σ 1 m : ℝ)

theorem minusEight_shell_error_nonneg (m : ℕ) : 0 ≤ minusEightShellError m := by
  unfold minusEightShellError
  positivity

theorem minusEight_shell_error_le_one (m : ℕ) : minusEightShellError m ≤ 1 := by
  by_cases hm : m = 0
  · simp [hm, minusEightShellError]
  · have hs : 0 < (σ 1 m : ℝ) := by exact_mod_cast ArithmeticFunction.sigma_pos 1 m hm
    apply (div_le_iff₀ hs).mpr
    simpa only [one_mul] using sigma_zero_mul_sqrt_le_sigma_one m

/-- The dyadic lower bound G₂/r ≥ (1−E)/2 used in the uniform sieve. -/
theorem minusEight_dyadic_shell_lower {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m : ℕ} (he : e ≤ 3) (hm : Odd m) (ha : |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m) :
    (1 / 2 - minusEightShellError m / 2) * (minusEightFullShell (2 ^ e * m)).card ≤
      (minusEightDyadicShell (2 ^ e * m)).card := by
  have hm₀ : 0 < m := by obtain ⟨k, hk⟩ := hm; omega
  have hs : 0 < (σ 1 m : ℝ) := by exact_mod_cast ArithmeticFunction.sigma_pos 1 m hm₀.ne'
  have hbound := mul_nonneg (minusEight_shell_error_nonneg m) hs.le
  have heq : minusEightShellError m * (σ 1 m : ℝ) = (σ 0 m : ℝ) * Real.sqrt m :=
    div_mul_cancel₀ _ hs.ne'
  rw [minusEight_full_shell_dyadic hθ he hm, minusEight_dyadic_shell_theta hθ he hm]
  interval_cases e <;> norm_num [minusEightShellFactor] <;>
    nlinarith only [hbound, hs, heq, le_abs_self (a m), ha]

/-- If an inert prime divides the target, the smaller-shell count gives
the required p⁻² bound without a Hecke operator at level p. -/
theorem minusEight_inert_target_bad_bound {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m p : ℕ} (he : e ≤ 3) (hm : Odd m) (hp : MinusEightInertPrime p) (hpn : p ∣ 2 ^ e * m) :
    ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
      (minusEightFullShell (2 ^ e * m)).card / (p : ℝ) ^ 2 := by
  have hp₂ : p ≠ 2 := by rcases hp.2 with h | h <;> omega
  have hco : (p ^ 2).Coprime (2 ^ e) := Nat.coprime_pow_primes _ _ hp.1 Nat.prime_two hp₂
  have hpR : 0 < (p : ℝ) := by exact_mod_cast hp.1.pos
  rw [minusEight_inert_target_bad_card hp hpn]
  split_ifs with hdiv
  · have hdivm := hco.dvd_of_dvd_mul_left hdiv
    have hodd : Odd (m / p ^ 2) := hm.of_dvd_nat (Nat.div_dvd_of_dvd hdivm)
    have hn : 0 < m := by obtain ⟨k, hk⟩ := hm; omega
    have hquot : 0 < m / p ^ 2 :=
      Nat.div_pos (Nat.le_of_dvd hn hdivm) (pow_pos hp.1.pos _)
    have heq : (2 ^ e * m) / p ^ 2 = 2 ^ e * (m / p ^ 2) := by
      conv_lhs => rw [← Nat.mul_div_cancel' hdivm, ← mul_assoc, mul_right_comm (2 ^ e)]
      exact Nat.mul_div_cancel _ (pow_pos hp.1.pos _)
    have hσ := pow_mul_sigma_le (k := 1) (pow_pos hp.1.pos 2) hquot
    rw [Nat.mul_div_cancel' hdivm, pow_one] at hσ
    have hσR : (p : ℝ) ^ 2 * (σ 1 (m / p ^ 2) : ℝ) ≤ σ 1 m := by exact_mod_cast hσ
    rw [heq, minusEight_full_shell_dyadic hθ he hodd, minusEight_full_shell_dyadic hθ he hm]
    apply (le_div_iff₀ (pow_pos hpR _)).mpr
    calc
      _ = 4 * (minusEightShellFactor e : ℝ) * ((p : ℝ) ^ 2 * (σ 1 (m / p ^ 2) : ℝ)) := by ring
      _ ≤ _ := mul_le_mul_of_nonneg_left hσR (by positivity)
  · simpa only [Nat.cast_zero] using div_nonneg (Nat.cast_nonneg _) (sq_nonneg (p : ℝ))

/-- The exact size of the selected odd-odd parts at 2m and 4m. -/
theorem minusEight_selected_even_shell_theta {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {m : ℕ} (hm : Odd m) :
    ((minusEightOddOddShell (2 * m)).card : ℝ) = 4 * (σ 1 m : ℝ) ∧
      ((minusEightOddOddShell (4 * m)).card : ℝ) = 16 * (σ 1 m : ℝ) := by
  have h₁ : ((minusEightFullShell (2 * m)).card : ℝ) =
      (minusEightOddOddShell (2 * m)).card + (minusEightFullShell m).card := by
    exact_mod_cast minusEight_odd_odd_shell_card m
  have h₂ : ((minusEightFullShell (2 * (2 * m))).card : ℝ) =
      (minusEightOddOddShell (2 * (2 * m))).card + (minusEightFullShell (2 * m)).card := by
    exact_mod_cast minusEight_odd_odd_shell_card (2 * m)
  have hz := minusEight_full_shell_dyadic hθ (e := 0) (by decide) hm
  have ho := minusEight_full_shell_dyadic hθ (e := 1) (by decide) hm
  have ht := minusEight_full_shell_dyadic hθ (e := 2) (by decide) hm
  norm_num [minusEightShellFactor, ← mul_assoc] at hz ho ht h₂
  exact ⟨by linarith only [h₁, hz, ho], by linarith only [h₂, ho, ht]⟩

end QuadraticNorms
