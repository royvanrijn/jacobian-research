import QuadraticNorms.MinusEightEulerBounds
import QuadraticNorms.MinusEightSieve

open scoped BigOperators ArithmeticFunction.sigma

namespace QuadraticNorms

private theorem minusEight_euler_error_le_scale {m : ℕ} {y C : ℝ} (hy : 0 < y) (hy12 : y ^ 12 = m)
    (hC : 0 ≤ C) (hτ : (σ 0 m : ℝ) ^ 6 ≤ C ^ 6 * m) :
    minusEightShellError m ≤ C / y ^ 4 := by
  have hm : (0 : ℝ) < m := hy12 ▸ pow_pos hy _
  have hn : 0 < m := by exact_mod_cast hm
  have hσN : m ≤ σ 1 m := by simpa using pow_mul_sigma_le (k := 1) hn Nat.one_pos
  have hσ : (m : ℝ) ≤ σ 1 m := by exact_mod_cast hσN
  have hτ' : (σ 0 m : ℝ) ≤ C * y ^ 2 := by
    apply (pow_le_pow_iff_left₀ (by positivity) (by positivity) (by decide : (6 : ℕ) ≠ 0)).mp
    rw [mul_pow, ← pow_mul, show 2 * 6 = 12 by decide, hy12]
    exact hτ
  have hroot : Real.sqrt (m : ℝ) = y ^ 6 := by
    rw [← hy12, show (12 : ℕ) = 6 * 2 by decide, pow_mul, Real.sqrt_sq (by positivity)]
  unfold minusEightShellError
  rw [hroot]
  calc
    _ ≤ (C * y ^ 2) * y ^ 6 / m := by
      exact div_le_div_of_nonneg_left (by positivity) hm hσ |>.trans
        (div_le_div_of_nonneg_right (mul_le_mul_of_nonneg_right hτ' (by positivity)) hm.le)
    _ = _ := by rw [← hy12]; field_simp

private theorem minusEight_uniform_tail_factor_le_scale {e m : ℕ} (he : e ≤ 3)
    {y : ℝ} (hy : 0 < y) (hy12 : y ^ 12 = m) :
    ((8 / 3 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum (2 ^ e * m)) ≤
      32 * y ^ 2 := by
  let n := 2 ^ e * m
  have hn : (n : ℝ) ≤ 8 * y ^ 12 := by
    have hpow : (2 : ℝ) ^ e ≤ 8 := (pow_le_pow_right₀ (by norm_num) he).trans_eq (by norm_num)
    dsimp [n]
    push_cast
    rw [← hy12]
    exact mul_le_mul_of_nonneg_right hpow (by positivity)
  have hM : (minusEightDivisorMaximum n : ℝ) ≤ 9 * y ^ 2 := by
    apply (pow_le_pow_iff_left₀ (by positivity) (by positivity) (by decide : (6 : ℕ) ≠ 0)).mp
    have hh := pow_le_pow_left₀ (show 0 ≤ (minusEightDivisorMaximum n : ℝ) by positivity)
      (minusEight_divisor_maximum_sixth_le n) 6
    rw [mul_pow, ← Real.rpow_mul_natCast (Nat.cast_nonneg n), show (1 / 6 : ℝ) * (6 : ℕ) = 1 by norm_num,
      Real.rpow_one] at hh
    refine hh.trans ((mul_le_mul_of_nonneg_left hn (by positivity)).trans ?_)
    rw [mul_pow, ← pow_mul, show 2 * 6 = 12 by decide]
    have hconst : (6 : ℝ) ^ 6 * 8 ≤ 9 ^ 6 := by norm_num
    simpa only [mul_assoc] using mul_le_mul_of_nonneg_right hconst (pow_nonneg hy.le 12)
  have hcoef : (8 / 3 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) ≤ 32 / 9 := by
    interval_cases e <;> norm_num [minusEightShellFactor]
  calc
    _ ≤ (32 / 9 : ℝ) * (9 * y ^ 2) := mul_le_mul hcoef hM (by positivity) (by norm_num)
    _ = _ := by ring

private theorem minusEight_margin_scaled {e m : ℕ} (he : e ≤ 3)
    {y C s : ℝ} (hy : 0 < y) (hy12 : y ^ 12 = m) (hC : 0 ≤ C) (hs : 0 < s)
    (hτ : (σ 0 m : ℝ) ^ 6 ≤ C ^ 6 * m) :
    3 / 16 + minusEightShellError m * ((s * y ^ 3) / 4 + 3 / 2) +
      ((8 / 3 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum (2 ^ e * m)) *
        (1 / (4 * (s * y ^ 3)) + 2 / (s * y ^ 3) ^ 2) ≤
      3 / 16 + (C * s / 4 + 8 / s) * (1 / y) + ((3 / 2) * C + 64 / s ^ 2) * (1 / y ^ 4) := by
  have hE := minusEight_euler_error_le_scale hy hy12 hC hτ
  have hL := minusEight_uniform_tail_factor_le_scale he hy hy12
  calc
    _ ≤ 3 / 16 + (C / y ^ 4) * ((s * y ^ 3) / 4 + 3 / 2) +
        (32 * y ^ 2) * (1 / (4 * (s * y ^ 3)) + 2 / (s * y ^ 3) ^ 2) := by
      exact add_le_add (add_le_add le_rfl (mul_le_mul_of_nonneg_right hE (by positivity)))
        (mul_le_mul_of_nonneg_right hL (by positivity))
    _ = _ := by field_simp; ring

/-- The uniform sieve has a positive arithmetic margin in every
eligible stratum once the odd part reaches 2⁷⁵. -/
theorem minusEight_general_large_margin_with_gap {e m : ℕ} (he : e ≤ 3) (hm : Odd m) (hlarge : 2 ^ 75 ≤ m) :
    ∃ X : ℝ, 0 < X ∧ 3 / 16 + minusEightShellError m * (X / 4 + 3 / 2) +
      ((8 / 3 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum (2 ^ e * m)) *
        (1 / (4 * X) + 2 / X ^ 2) < 1 / 2 - 1 / 100 := by
  let y := (m : ℝ) ^ (1 / 12 : ℝ)
  have hy12 : y ^ 12 = (m : ℝ) := by
    dsimp [y]
    rw [← Real.rpow_mul_natCast (Nat.cast_nonneg _) _ 12]
    norm_num
  have hy₀ : 0 ≤ y := by dsimp [y]; positivity
  have hy : 64 ≤ y := by
    apply (pow_le_pow_iff_left₀ (by norm_num) hy₀ (by decide : (12 : ℕ) ≠ 0)).mp
    rw [hy12]
    exact_mod_cast (by norm_num : (64 : ℕ) ^ 12 ≤ 2 ^ 75).trans hlarge
  have hy' : 0 < y := by linarith
  have hτ : (σ 0 m : ℝ) ^ 6 ≤ (42 : ℝ) ^ 6 * m := by
    exact_mod_cast (minusEight_odd_divisor_sixth_power_lt hm).le
  have h := minusEight_margin_scaled he hy' hy12 (C := 42) (s := 1) (by norm_num) (by norm_num) hτ
  refine ⟨y ^ 3, pow_pos hy' _, ?_⟩
  simp only [one_mul] at h
  refine h.trans_lt ?_
  have h₁ := one_div_le_one_div_of_le (by norm_num : (0 : ℝ) < 64) hy
  have h₂ := one_div_le_one_div_of_le (by norm_num : (0 : ℝ) < (64 : ℝ) ^ 4)
    (pow_le_pow_left₀ (by norm_num) hy 4)
  norm_num at h₁ h₂ ⊢
  linarith only [h₁, h₂]

theorem minusEight_general_large_margin {e m : ℕ} (he : e ≤ 3) (hm : Odd m) (hlarge : 2 ^ 75 ≤ m) :
    ∃ X : ℝ, 0 < X ∧ 3 / 16 + minusEightShellError m * (X / 4 + 3 / 2) +
      ((8 / 3 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum (2 ^ e * m)) *
        (1 / (4 * X) + 2 / X ^ 2) < 1 / 2 := by
  obtain ⟨X, hX, hg⟩ := minusEight_general_large_margin_with_gap he hm hlarge
  exact ⟨X, hX, hg.trans (by norm_num)⟩

/-- The uniform sieve has a positive arithmetic margin on inert cores
once the odd part reaches 2⁶⁰. -/
theorem minusEight_inert_large_margin {e m : ℕ} (he : e ≤ 3)
    (hI : ∀ p ∈ m.primeFactors, MinusEightInertPrime p) (hlarge : 2 ^ 60 ≤ m) :
    ∃ X : ℝ, 0 < X ∧ 3 / 16 + minusEightShellError m * (X / 4 + 3 / 2) +
      ((8 / 3 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum (2 ^ e * m)) *
        (1 / (4 * X) + 2 / X ^ 2) < 1 / 2 := by
  let y := (m : ℝ) ^ (1 / 12 : ℝ)
  have hy12 : y ^ 12 = (m : ℝ) := by
    dsimp [y]
    rw [← Real.rpow_mul_natCast (Nat.cast_nonneg _) _ 12]
    norm_num
  have hy₀ : 0 ≤ y := by dsimp [y]; positivity
  have hy : 32 ≤ y := by
    apply (pow_le_pow_iff_left₀ (by norm_num) hy₀ (by decide : (12 : ℕ) ≠ 0)).mp
    rw [hy12]
    exact_mod_cast (by norm_num : (32 : ℕ) ^ 12 ≤ 2 ^ 60).trans hlarge
  have hy' : 0 < y := by linarith
  have hτ : (σ 0 m : ℝ) ^ 6 ≤ (7 : ℝ) ^ 6 * m := by
    exact_mod_cast (minusEight_inert_divisor_sixth_power_lt (by omega : 0 < m) hI).le
  have h := minusEight_margin_scaled he hy' hy12 (C := 7) (s := 2) (by norm_num) (by norm_num) hτ
  refine ⟨2 * y ^ 3, by positivity, h.trans_lt ?_⟩
  have h₁ := one_div_le_one_div_of_le (by norm_num : (0 : ℝ) < 32) hy
  have h₂ := one_div_le_one_div_of_le (by norm_num : (0 : ℝ) < (32 : ℝ) ^ 4)
    (pow_le_pow_left₀ (by norm_num) hy 4)
  norm_num at h₁ h₂ ⊢
  linarith only [h₁, h₂]

private theorem minusEight_large_dyadic_decomposition {n : ℕ}
    (hn : 2 ^ 78 ≤ n) (h₁₆ : ¬ 16 ∣ n) :
    ∃ e m : ℕ, e ≤ 3 ∧ Odd m ∧ 2 ^ 75 ≤ m ∧ 2 ^ e * m = n := by
  obtain ⟨e, m, hm, rfl⟩ := Nat.exists_eq_two_pow_mul_odd (by omega : n ≠ 0)
  have he : e ≤ 3 := by
    by_contra! h
    exact h₁₆ ((Nat.pow_dvd_pow 2 h).trans (dvd_mul_right _ _))
  have hm' : 2 ^ 75 ≤ m := by
    have hpow : 2 ^ e ≤ 8 := (Nat.pow_le_pow_right (by decide) he).trans_eq (by decide)
    have h := hn.trans (Nat.mul_le_mul_right m hpow)
    omega
  exact ⟨e, m, he, hm, hm', rfl⟩

/-- The stated global onset follows from the two theta identities, the
odd cusp coefficient bound, and the unit-level small-prime estimate. -/
theorem minusEight_large_target_sum {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    (hcusp : ∀ m, Odd m → |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m)
    (hsmall : ∀ e m p : ℕ, e ≤ 3 → Odd m → p ≤ 2 ^ e * m → MinusEightInertPrime p →
      ¬ p ∣ 2 ^ e * m → ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card)
    {n : ℕ} (hn : 2 ^ 78 ≤ n) (h₁₆ : ¬ 16 ∣ n) : InSumset P8 P8 n := by
  obtain ⟨e, m, he, hm, hm', rfl⟩ := minusEight_large_dyadic_decomposition hn h₁₆
  obtain ⟨X, hX, hmargin⟩ := minusEight_general_large_margin he hm hm'
  exact minusEight_sum_of_arithmetic_margin hθ he hm (hcusp m hm) hX
    (fun p hp hi _ hd => hsmall e m p he hm hp hi hd) hmargin

/-- All three reduced core strata beyond the inert-part onset use the
same theta and unit-level small-prime inputs as the global bound. -/
theorem minusEight_large_inert_core_sums {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    (hcusp : ∀ m, Odd m → |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m)
    (hsmall : ∀ e m p : ℕ, e ≤ 3 → Odd m → p ≤ 2 ^ e * m → MinusEightInertPrime p →
      ¬ p ∣ 2 ^ e * m → ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card)
    {m : ℕ} (hm : Odd m) (hI : ∀ p ∈ m.primeFactors, MinusEightInertPrime p)
    (hlarge : 2 ^ 60 ≤ m) :
    InSumset C8 (Dilate 2 C8) m ∧ InSumset C8 C8 (2 * m) ∧ InSumset C8 C8 (4 * m) := by
  have hmixed : InSumset C8 (Dilate 2 C8) m := by
    obtain ⟨X, hX, hmargin⟩ := minusEight_inert_large_margin (e := 0) (by decide) hI hlarge
    exact minusEight_mixed_core_of_selected_margin hθ hm (hcusp m hm) hX
      (fun p hpn hi _ hd => by
        simpa only [pow_zero, one_mul] using hsmall 0 m p (by decide) hm (by simpa using hpn) hi
          (by simpa using hd))
      (by simpa [minusEightShellFactor] using hmargin)
  have heven (e : ℕ) (he : e = 1 ∨ e = 2) : InSumset C8 C8 (2 ^ e * m) := by
    have he₃ : e ≤ 3 := by omega
    obtain ⟨X, hX, hmargin⟩ := minusEight_inert_large_margin he₃ hI hlarge
    refine minusEight_core_sum_of_selected_margin hθ he hm hX
      (fun p hp hi _ hd => hsmall e m p he₃ hm hp hi hd) ?_
    have hE := minusEight_shell_error_nonneg m
    have hleft : 3 / 16 + minusEightShellError m * (X / 4 + 1) +
        ((2 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum (2 ^ e * m)) *
          (1 / (4 * X) + 2 / X ^ 2) ≤
        3 / 16 + minusEightShellError m * (X / 4 + 3 / 2) +
        ((8 / 3 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum (2 ^ e * m)) *
          (1 / (4 * X) + 2 / X ^ 2) := by
      gcongr <;> norm_num
    exact (hleft.trans_lt hmargin).trans_le (by split <;> norm_num)
  exact ⟨hmixed, by simpa using heven 1 (Or.inl rfl), heven 2 (Or.inr rfl)⟩

/-- At the prime-family boundary the mixed sieve has a positive
arithmetic margin, leaving only the analytic unit-level input. -/
theorem minusEight_large_prime_mixed_core {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {p : ℕ} (hp : p.Prime) (hcut : 2 ^ 22 ≤ p)
    (hcusp : |a p| ≤ (σ 0 p : ℝ) * Real.sqrt p)
    (hsmall : ∀ q ≤ p, MinusEightInertPrime q → ¬ q ∣ p →
      ((minusEightEitherDivisible p q).card : ℝ) ≤
        (minusEightLocalMass q + minusEightShellError p) * (minusEightFullShell p).card) :
    InSumset C8 (Dilate 2 C8) p := by
  have hodd := hp.odd_of_ne_two (by omega)
  let y := (p : ℝ) ^ (1 / 8 : ℝ)
  have hy₀ : 0 ≤ y := by dsimp [y]; positivity
  have hy8 : y ^ 8 = (p : ℝ) := by
    dsimp [y]
    rw [← Real.rpow_mul_natCast (Nat.cast_nonneg _) _ 8]
    norm_num
  have hy : 6 ≤ y := by
    apply (pow_le_pow_iff_left₀ (by norm_num) hy₀ (by decide : (8 : ℕ) ≠ 0)).mp
    rw [hy8]
    exact_mod_cast (by norm_num : (6 : ℕ) ^ 8 ≤ 2 ^ 22).trans hcut
  have hy' : 0 < y := by linarith
  have hp' : 0 < (p : ℝ) := by exact_mod_cast hp.pos
  have hpσ : (p : ℝ) ≤ (σ 1 p : ℝ) := by
    have h : p ≤ σ 1 p := by simpa using pow_mul_sigma_le (k := 1) hp.pos Nat.one_pos
    exact_mod_cast h
  have hτ : (σ 0 p : ℝ) = 2 := by
    have h : σ 0 p = 2 := by simpa using ArithmeticFunction.sigma_zero_apply_prime_pow (i := 1) hp
    exact_mod_cast h
  have hroot : Real.sqrt (p : ℝ) = y ^ 4 := by
    rw [← hy8, show (8 : ℕ) = 4 * 2 by decide, pow_mul, Real.sqrt_sq (by positivity)]
  have hE : minusEightShellError p ≤ 2 / y ^ 4 := by
    rw [minusEightShellError, hτ, hroot]
    refine (div_le_div_of_nonneg_left (by positivity) hp' hpσ).trans_eq ?_
    rw [← hy8]
    field_simp
  have hM : (minusEightDivisorMaximum p : ℝ) ≤ 2 * y ^ 2 := by
    refine (minusEight_divisor_maximum_fourth_le p).trans_eq ?_
    congr 1
    rw [← hy8, ← Real.rpow_natCast, ← Real.rpow_mul hy₀]
    norm_num
  let X := 2 * y ^ 3
  have hX : 0 < X := by dsimp [X]; positivity
  refine minusEight_mixed_core_of_selected_margin hθ hodd hcusp hX
    (fun q hq hi _ hd => hsmall q hq hi hd) ?_
  have hbound :
      3 / 16 + minusEightShellError p * (X / 4 + 3 / 2) +
        ((8 / 3 : ℝ) * minusEightDivisorMaximum p) * (1 / (4 * X) + 2 / X ^ 2) ≤
      3 / 16 + (5 / 3) * (1 / y) + (17 / 3) * (1 / y ^ 4) := by
    calc
      _ ≤ 3 / 16 + (2 / y ^ 4) * (X / 4 + 3 / 2) +
          ((8 / 3 : ℝ) * (2 * y ^ 2)) * (1 / (4 * X) + 2 / X ^ 2) := by
        gcongr
      _ = _ := by dsimp [X]; field_simp; ring
  refine hbound.trans_lt ?_
  have h₁ := one_div_le_one_div_of_le (by norm_num : (0 : ℝ) < 6) hy
  have h₂ := one_div_le_one_div_of_le (by norm_num : (0 : ℝ) < (6 : ℝ) ^ 4)
    (pow_le_pow_left₀ (by norm_num) hy 4)
  norm_num at h₁ h₂ ⊢
  linarith only [h₁, h₂]

/-- The ordered primitive value-pair count, indexed by the first value. -/
noncomputable def minusEightOrderedPairCount (n : ℕ) : ℕ := by
  classical
  exact pairCount P8 P8 n

private theorem minusEight_survivors_le_ordered_pair_count (n : ℕ) :
    (minusEightSieveSurvivors n (minusEightDyadicShell n)).card ≤
      (2 * minusEightDivisorMaximum n) ^ 2 * minusEightOrderedPairCount n := by
  classical
  let s := minusEightSieveSurvivors n (minusEightDyadicShell n)
  let f (v : MinusEightShellPoint) := minusEightNorm v.1
  have hvalues (v : MinusEightShellPoint) (hv : v ∈ s) :
      P8 (minusEightNorm v.1) ∧ P8 (minusEightNorm v.2) ∧
        minusEightNorm v.1 + minusEightNorm v.2 = n := by
    obtain ⟨hd, ha⟩ := mem_minusEight_sieve_survivors.mp hv
    obtain ⟨h₁, h₂⟩ := minusEight_values_of_unspoiled hd ha
    exact ⟨h₁, h₂, mem_minusEight_full_shell.mp (Finset.mem_filter.mp hd).1⟩
  change s.card ≤ _ * (pairFinset P8 P8 n).card
  apply Finset.card_le_mul_card_image_of_maps_to (f := f)
  · intro v hv
    have ⟨h₁, h₂, hsum⟩ := hvalues v hv
    exact mem_pairFinset.mpr ⟨by dsimp [f]; omega, h₁, by
      simpa only [f, ← hsum, Nat.add_sub_cancel_left] using h₂⟩
  · intro a ha
    obtain ⟨han, h₁, h₂⟩ := mem_pairFinset.mp ha
    have hsub : s.filter (fun v => f v = a) ⊆
        minusEightNormFiber a ×ˢ minusEightNormFiber (n - a) := by
      intro v hv
      obtain ⟨hvs, hvf⟩ := Finset.mem_filter.mp hv
      have hsum := (hvalues v hvs).2.2
      change minusEightNorm v.1 = a at hvf
      exact Finset.mem_product.mpr
        ⟨mem_minusEight_norm_fiber.mpr hvf, mem_minusEight_norm_fiber.mpr (by omega)⟩
    refine (Finset.card_le_card hsub).trans ?_
    rw [Finset.card_product, pow_two]
    exact Nat.mul_le_mul ((minusEight_norm_fiber_card_le h₁.1).trans
      (Nat.mul_le_mul_left _ (minusEight_split_divisors_le_maximum h₁.1 han)))
      ((minusEight_norm_fiber_card_le h₂.1).trans
        (Nat.mul_le_mul_left _ (minusEight_split_divisors_le_maximum h₂.1 (Nat.sub_le ..))))

private theorem minusEight_full_shell_ge_three_mul_target {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m : ℕ} (he : e ≤ 3) (hm : Odd m) :
    3 * (2 ^ e * m : ℕ) ≤ ((minusEightFullShell (2 ^ e * m)).card : ℝ) := by
  have hm₀ : 0 < m := by obtain ⟨k, hk⟩ := hm; omega
  have hσ : (m : ℝ) ≤ σ 1 m := by
    have h : m ≤ σ 1 m := by simpa using pow_mul_sigma_le (k := 1) hm₀ Nat.one_pos
    exact_mod_cast h
  rw [minusEight_full_shell_dyadic hθ he hm]
  have he' : (3 : ℝ) * 2 ^ e ≤ 4 * (minusEightShellFactor e : ℝ) := by
    interval_cases e <;> norm_num [minusEightShellFactor]
  push_cast
  rw [← mul_assoc]
  exact mul_le_mul he' hσ (by positivity) (by positivity)

/-- The quantitative −8 bound on ordered primitive value pairs. It uses
the same three analytic inputs as the classification. -/
theorem minusEight_large_ordered_pair_count {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    (hcusp : ∀ m, Odd m → |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m)
    (hsmall : ∀ e m p : ℕ, e ≤ 3 → Odd m → p ≤ 2 ^ e * m → MinusEightInertPrime p →
      ¬ p ∣ 2 ^ e * m → ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card)
    {n : ℕ} (hn : 2 ^ 78 ≤ n) (h₁₆ : ¬ 16 ∣ n) :
    (n : ℝ) ^ (2 / 3 : ℝ) / 4800 ≤ (minusEightOrderedPairCount n : ℝ) := by
  obtain ⟨e, m, he, hm, hmCut, rfl⟩ := minusEight_large_dyadic_decomposition hn h₁₆
  let n := 2 ^ e * m
  obtain ⟨X, hX, hmargin⟩ := minusEight_general_large_margin_with_gap he hm hmCut
  have hsurv := minusEight_uniform_survivor_card_le hθ he hm (hcusp m hm) hX
    (fun p hp hi _ hd => hsmall e m p he hm hp hi hd)
  have hcard : (3 / 100 : ℝ) * n ≤ (minusEightSieveSurvivors n (minusEightDyadicShell n)).card := by
    calc
      _ = (1 / 100 : ℝ) * (3 * n) := by ring
      _ ≤ (1 / 100 : ℝ) * (minusEightFullShell n).card :=
        mul_le_mul_of_nonneg_left (minusEight_full_shell_ge_three_mul_target hθ he hm) (by norm_num)
      _ ≤ _ := (mul_le_mul_of_nonneg_right (by linarith only [hmargin]) (by positivity)).trans hsurv
  have hcount : ((minusEightSieveSurvivors n (minusEightDyadicShell n)).card : ℝ) ≤
      (2 * (minusEightDivisorMaximum n : ℝ)) ^ 2 * minusEightOrderedPairCount n := by
    exact_mod_cast minusEight_survivors_le_ordered_pair_count n
  have hM : (2 * (minusEightDivisorMaximum n : ℝ)) ^ 2 ≤ 144 * (n : ℝ) ^ (1 / 3 : ℝ) := by
    have h := pow_le_pow_left₀ (show 0 ≤ 2 * (minusEightDivisorMaximum n : ℝ) by positivity)
      (mul_le_mul_of_nonneg_left (minusEight_divisor_maximum_sixth_le n) (by norm_num : (0 : ℝ) ≤ 2)) 2
    refine h.trans_eq ?_
    rw [show (2 : ℝ) * (6 * (n : ℝ) ^ (1 / 6 : ℝ)) = 12 * (n : ℝ) ^ (1 / 6 : ℝ) by ring,
      mul_pow, ← Real.rpow_mul_natCast (Nat.cast_nonneg n)]
    norm_num
  have hN : 0 < (n : ℝ) := by
    have hNN : 0 < n := (by decide : 0 < 2 ^ 78).trans_le hn
    exact Nat.cast_pos.mpr hNN
  have hpos := Real.rpow_pos_of_pos hN (1 / 3)
  change (n : ℝ) ^ (2 / 3 : ℝ) / 4800 ≤ (minusEightOrderedPairCount n : ℝ)
  apply (le_of_mul_le_mul_left · (show (0 : ℝ) < 144 * (n : ℝ) ^ (1 / 3 : ℝ) by positivity))
  calc
    _ = (3 / 100 : ℝ) * n := by
      have hpow : (n : ℝ) ^ (1 / 3 : ℝ) * (n : ℝ) ^ (2 / 3 : ℝ) = n := by
        rw [← Real.rpow_add hN, show (1 / 3 : ℝ) + 2 / 3 = 1 by norm_num, Real.rpow_one]
      calc
        _ = (3 / 100 : ℝ) * ((n : ℝ) ^ (1 / 3 : ℝ) * (n : ℝ) ^ (2 / 3 : ℝ)) := by ring
        _ = _ := by rw [hpow]
    _ ≤ _ := hcard.trans (hcount.trans (mul_le_mul_of_nonneg_right hM (by positivity)))

end QuadraticNorms
