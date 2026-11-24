import QuadraticNorms.MinusElevenShell
import Mathlib.NumberTheory.ArithmeticFunction.Misc

open scoped BigOperators ArithmeticFunction.sigma

namespace QuadraticNorms

private theorem sigma_one_two_power (k : ℕ) : σ 1 (2 ^ k) = 2 ^ (k + 1) - 1 := by
  induction k with
  | zero => norm_num
  | succ k ih =>
    have hr := sigma_one_prime_pow_succ (a := k) Nat.prime_two
    rw [hr, ih, pow_succ 2 (k + 1)]
    have hpos : 0 < (2 : ℕ) ^ (k + 1) := by positivity
    omega

private theorem sigma_one_two_power_odd (k m : ℕ) (hm : Odd m) :
    σ 1 (2 ^ k * m) = (2 ^ (k + 1) - 1) * σ 1 m := by
  have hc := Nat.Prime.coprime_iff_not_dvd Nat.prime_two |>.mpr hm.not_two_dvd_nat
  rw [ArithmeticFunction.isMultiplicative_sigma.map_mul_of_coprime (hc.pow_left k), sigma_one_two_power]

noncomputable def minusElevenFilteredSigma (n : ℕ) : ℝ :=
  (σ 1 n : ℝ) - (if 4 ∣ n then (σ 1 (n / 4) : ℝ) else 0)

theorem minusEleven_filtered_sigma {t m : ℕ} (ht : 0 < t) (hm : Odd m) :
    minusElevenFilteredSigma (2 ^ t * m) = (3 / 2 : ℝ) * (2 : ℝ) ^ t * (σ 1 m : ℝ) := by
  unfold minusElevenFilteredSigma
  cases t with
  | zero => omega
  | succ t =>
    cases t with
    | zero =>
      have hn : ¬ 4 ∣ 2 ^ Nat.succ 0 * m := by
        obtain ⟨k, hk⟩ := hm
        norm_num
        omega
      rw [if_neg hn, sigma_one_two_power_odd _ _ hm]
      norm_num
    | succ k =>
      have hrewrite : 2 ^ (k + 1 + 1) * m = 4 * (2 ^ k * m) := by rw [pow_succ 2 (k+1), pow_succ 2 k]; ring
      change minusElevenFilteredSigma (2 ^ (k + 1 + 1) * m) = _
      unfold minusElevenFilteredSigma
      have hd : 4 ∣ 2 ^ (k + 1 + 1) * m := ⟨_, hrewrite⟩
      have hdiv : (2 ^ (k + 1 + 1) * m) / 4 = 2 ^ k * m := by
        rw [hrewrite, Nat.mul_div_cancel_left _ (by decide : 0 < 4)]
      rw [if_pos hd, hdiv, sigma_one_two_power_odd _ _ hm, sigma_one_two_power_odd _ _ hm]
      have h₁ : 1 ≤ (2 : ℕ) ^ (k + 1) := Nat.succ_le_of_lt (by positivity)
      have h₂ : 1 ≤ (2 : ℕ) ^ (k + 1 + 1 + 1) := Nat.succ_le_of_lt (by positivity)
      simp only [Nat.cast_mul, Nat.cast_sub h₁, Nat.cast_sub h₂, Nat.cast_pow, Nat.cast_ofNat, Nat.cast_one]
      rw [pow_succ (2 : ℝ) (k+1+1), pow_succ (2 : ℝ) (k+1), pow_succ (2 : ℝ) k]
      ring

noncomputable def minusElevenCuspErrorBound (t m : ℕ) : ℝ :=
  (σ 0 m : ℝ) * Real.sqrt m * Real.sqrt ((2 : ℝ) ^ t) * ((t : ℝ) + 3) / 2

theorem minusEleven_filtered_error_eq {t m : ℕ} (ht : 0 < t) (hm : Odd m) :
    minusElevenFilteredError t m =
      2 * minusElevenCuspErrorBound t m / (3 * minusElevenFilteredSigma (2 ^ t * m)) := by
  rw [minusEleven_filtered_sigma ht hm]
  unfold minusElevenFilteredError minusElevenCuspErrorBound
  have hm₀ : 0 < m := by
    obtain ⟨k, hk⟩ := hm
    omega
  have hs : (0 : ℝ) < σ 1 m := by exact_mod_cast ArithmeticFunction.sigma_pos 1 m hm₀.ne'
  have ht₀ : (0 : ℝ) < (2 : ℝ) ^ t := by positivity
  have hr : (0 : ℝ) < Real.sqrt ((2 : ℝ) ^ t) := Real.sqrt_pos.2 ht₀
  have hsq := Real.sq_sqrt ht₀.le
  field_simp [hs.ne', ht₀.ne', hr.ne']
  ring_nf
  rw [hsq]


theorem minusEleven_filtered_error_mul {t m : ℕ} (ht : 0 < t) (hm : Odd m) :
    3 * minusElevenFilteredSigma (2 ^ t * m) * minusElevenFilteredError t m =
      2 * minusElevenCuspErrorBound t m := by
  have hm₀ : 0 < m := by obtain ⟨k, hk⟩ := hm; omega
  have hσ : (0 : ℝ) < σ 1 m := by exact_mod_cast ArithmeticFunction.sigma_pos 1 m hm₀.ne'
  have hS : 0 < minusElevenFilteredSigma (2 ^ t * m) := by
    rw [minusEleven_filtered_sigma ht hm]
    positivity
  rw [mul_comm (3 * minusElevenFilteredSigma (2 ^ t * m))]
  exact (eq_div_iff (by positivity)).mp (minusEleven_filtered_error_eq ht hm)

/-- The cusp coefficient after subtracting the term at the integral index n/4. -/
noncomputable def minusElevenFilteredCusp (a : ℕ → ℝ) (n : ℕ) : ℝ :=
  a n - if 4 ∣ n then a (n / 4) else 0

/-- The theta coefficient identity on positive targets prime to eleven. -/
def MinusElevenThetaFormula (a : ℕ → ℝ) : Prop :=
  ∀ j : ℕ, 0 < j → ¬ 11 ∣ j →
    ((minusElevenFullShell j).card : ℝ) = (12 / 5) * (σ 1 j : ℝ) + (8 / 5) * a j

/-- The odd-shell theta identity follows from the theta coefficients
on positive eleven-unit targets. -/
theorem minusEleven_odd_shell_theta (a : ℕ → ℝ)
    (hTheta : MinusElevenThetaFormula a)
    {n : ℕ} (hn : 0 < n) (he : n % 2 = 0) (h₁₁ : ¬ 11 ∣ n) :
    ((minusElevenOddShell n).card : ℝ) =
      (12 / 5) * minusElevenFilteredSigma n + (8 / 5) * minusElevenFilteredCusp a n := by
  have h := minusEleven_full_shell_card_add he
  have ht := hTheta n hn h₁₁
  by_cases hd : 4 ∣ n
  · have hquarter : 0 < n / 4 := Nat.div_pos (Nat.le_of_dvd hn hd) (by decide)
    have hq₁₁ : ¬ 11 ∣ n / 4 := fun hdiv =>
      h₁₁ (hdiv.trans (Nat.div_dvd_of_dvd hd))
    rw [if_pos hd] at h
    have hcast : ((minusElevenFullShell n).card : ℝ) =
        ((minusElevenOddShell n).card : ℝ) + ((minusElevenFullShell (n / 4)).card : ℝ) := by
      exact_mod_cast h
    simp only [minusElevenFilteredSigma, minusElevenFilteredCusp, if_pos hd]
    linarith only [hcast, ht, hTheta (n / 4) hquarter hq₁₁]
  · rw [if_neg hd, Nat.add_zero] at h
    have hcast : ((minusElevenFullShell n).card : ℝ) = (minusElevenOddShell n).card := by
      exact_mod_cast h
    simp only [minusElevenFilteredSigma, minusElevenFilteredCusp, if_neg hd, sub_zero]
    exact hcast.symm.trans ht

/-- The odd-shell count lies in the interval supplied by the theta
identity and the bound on the filtered cusp term. -/
theorem minusEleven_odd_shell_bounds_of_theta (a : ℕ → ℝ)
    (hTheta : MinusElevenThetaFormula a) {t m : ℕ} (ht : 0 < t) (hm : Odd m)
    (h₁₁ : ¬ 11 ∣ 2 ^ t * m)
    (hcusp : |minusElevenFilteredCusp a (2 ^ t * m)| ≤ minusElevenCuspErrorBound t m) :
    let n := 2 ^ t * m
    let S := minusElevenFilteredSigma n
    let E := minusElevenFilteredError t m
    (12 / 5 : ℝ) * S * (1 - E) ≤ (minusElevenOddShell n).card ∧
      ((minusElevenOddShell n).card : ℝ) ≤ (12 / 5 : ℝ) * S + (8 / 5 : ℝ) * minusElevenCuspErrorBound t m := by
  let n := 2 ^ t * m
  let E := minusElevenFilteredError t m
  let S := minusElevenFilteredSigma n
  let D := minusElevenCuspErrorBound t m
  have hm₀ : 0 < m := by
    obtain ⟨k, hk⟩ := hm
    omega
  have hn : 0 < n := Nat.mul_pos (by positivity) hm₀
  have he : n % 2 = 0 := by
    apply Nat.dvd_iff_mod_eq_zero.mp
    exact (dvd_pow_self 2 ht.ne').trans (dvd_mul_right _ _)
  have hG := minusEleven_odd_shell_theta a hTheta hn he h₁₁
  have hpow : (3 * S) * E = 2 * D := minusEleven_filtered_error_mul ht hm
  obtain ⟨hClo, hChi⟩ := abs_le.mp hcusp
  change ((minusElevenOddShell n).card : ℝ) = (12 / 5) * S + (8 / 5) * minusElevenFilteredCusp a n at hG
  change (12 / 5 : ℝ) * S * (1 - E) ≤ (minusElevenOddShell n).card ∧
    ((minusElevenOddShell n).card : ℝ) ≤ (12 / 5 : ℝ) * S + (8 / 5 : ℝ) * D
  exact ⟨by nlinarith only [hG, hpow, hClo], by linarith only [hG, hChi]⟩

/-- The required shell lower bound follows from the exact theta formula
and the absolute bound on its filtered cusp coefficient. -/
theorem minusEleven_shell_lower_of_theta (a : ℕ → ℝ)
    (hTheta : MinusElevenThetaFormula a)
    {t m : ℕ} (ht : 0 < t) (hm : Odd m) (h₁₁ : ¬ 11 ∣ 2 ^ t * m)
    (hcusp : |minusElevenFilteredCusp a (2 ^ t * m)| ≤ minusElevenCuspErrorBound t m)
    (hE : minusElevenFilteredError t m ≤ 1) :
    (18 / 5 : ℝ) * (2 ^ t * m : ℕ) * (1 - minusElevenFilteredError t m) ≤
      (minusElevenOddShell (2 ^ t * m)).card := by
  let n := 2 ^ t * m
  let E := minusElevenFilteredError t m
  let S := minusElevenFilteredSigma n
  have hm₀ : 0 < m := by
    obtain ⟨k, hk⟩ := hm
    omega
  have hσ : m ≤ σ 1 m := by
    rw [ArithmeticFunction.sigma_one_apply]
    exact Finset.single_le_sum (f := id) (fun _ _ => Nat.zero_le _)
      (Nat.mem_divisors.mpr ⟨dvd_rfl, hm₀.ne'⟩)
  have hS : (3 / 2 : ℝ) * n ≤ S := by
    rw [show S = minusElevenFilteredSigma (2 ^ t * m) from rfl, minusEleven_filtered_sigma ht hm]
    push_cast [n]
    have hσR : (m : ℝ) ≤ σ 1 m := by exact_mod_cast hσ
    simpa only [mul_assoc] using mul_le_mul_of_nonneg_left hσR
      (by positivity : 0 ≤ (3 / 2 : ℝ) * (2 : ℝ) ^ t)
  calc
    (18 / 5 : ℝ) * n * (1 - E) = (12 / 5) * ((3 / 2) * n) * (1 - E) := by ring
    _ ≤ (12 / 5) * S * (1 - E) := mul_le_mul_of_nonneg_right
      (mul_le_mul_of_nonneg_left hS (by norm_num)) (sub_nonneg.mpr hE)
    _ ≤ _ := (minusEleven_odd_shell_bounds_of_theta a hTheta ht hm h₁₁ hcusp).1

end QuadraticNorms
