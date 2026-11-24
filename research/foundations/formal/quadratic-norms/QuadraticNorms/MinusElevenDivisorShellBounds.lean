import QuadraticNorms.MinusElevenThetaNormalization
import QuadraticNorms.SigmaBounds

open scoped BigOperators ArithmeticFunction.sigma

namespace QuadraticNorms

private theorem two_power_sqrt_bound (j : ℕ) :
    (j : ℝ) + 4 ≤ 2 * Real.sqrt 2 * Real.sqrt ((2 : ℝ) ^ (j + 1)) := by
  have h₂ : Real.sqrt 2 ^ 2 = (2 : ℝ) := Real.sq_sqrt (by norm_num)
  induction j with
  | zero =>
    simp only [Nat.cast_zero, zero_add, pow_one]
    nlinarith only [h₂]
  | succ j ih =>
    have hsq : (5 / 4 : ℝ) ≤ Real.sqrt 2 :=
      (Real.le_sqrt (by norm_num) (by norm_num)).mpr (by norm_num)
    have hincr : (j : ℝ) + 5 ≤ ((j : ℝ) + 4) * Real.sqrt 2 := by
      have h := mul_le_mul_of_nonneg_left hsq (by positivity : 0 ≤ (j : ℝ) + 4)
      linarith only [h, show 0 ≤ (j : ℝ) by positivity]
    have hd := hincr.trans (mul_le_mul_of_nonneg_right ih (Real.sqrt_nonneg 2))
    rw [pow_succ (2 : ℝ) (j + 1), Real.sqrt_mul (by positivity)]
    push_cast
    linarith only [hd]

/-- The factor used in the shell estimates is uniformly less than one
on the even targets, independently of the prime support of m. -/
theorem minusEleven_filtered_error_lt_one {t m : ℕ} (ht : 0 < t) :
    minusElevenFilteredError t m < 1 := by
  by_cases hm : m = 0
  · simp [hm, minusElevenFilteredError]
  let W := Real.sqrt ((2 : ℝ) ^ t)
  let c := (4 / 9 : ℝ) * Real.sqrt 2
  have hs : (0 : ℝ) < σ 1 m := by exact_mod_cast ArithmeticFunction.sigma_pos 1 m hm
  have hw : 0 < W := Real.sqrt_pos.2 (by positivity)
  have hgrowth : (t : ℝ) + 3 ≤ 2 * Real.sqrt 2 * W := by
    obtain ⟨j, rfl⟩ := Nat.exists_eq_succ_of_ne_zero ht.ne'
    dsimp [W]
    push_cast
    linarith only [two_power_sqrt_bound j]
  have hc : c < 1 := by
    have hh : Real.sqrt 2 ≤ (3 / 2 : ℝ) :=
      Real.sqrt_le_iff.mpr (by norm_num)
    dsimp [c]
    linarith only [hh]
  refine lt_of_le_of_lt ?_ hc
  change 2 * ((t : ℝ) + 3) / 9 * (σ 0 m : ℝ) * Real.sqrt m / (W * σ 1 m) ≤ c
  apply (div_le_iff₀ (mul_pos hw hs)).mpr
  calc
    _ = (2 * ((t : ℝ) + 3) / 9) * ((σ 0 m : ℝ) * Real.sqrt m) := by ring
    _ ≤ (2 * ((t : ℝ) + 3) / 9) * (σ 1 m : ℝ) :=
      mul_le_mul_of_nonneg_left (sigma_zero_mul_sqrt_le_sigma_one m) (by positivity)
    _ ≤ (c * W) * (σ 1 m : ℝ) := by
      apply mul_le_mul_of_nonneg_right _ hs.le
      dsimp [c]
      linarith only [hgrowth]
    _ = _ := by ring

private theorem filtered_sigma_scale_le {p t m : ℕ} (hp : 0 < p) (ht : 0 < t) (hm : Odd m)
    (hdiv : p ^ 2 ∣ m) :
    (p : ℝ) ^ 2 * minusElevenFilteredSigma (2 ^ t * (m / p ^ 2)) ≤
      minusElevenFilteredSigma (2 ^ t * m) := by
  have hm₀ : 0 < m := by
    obtain ⟨k, hk⟩ := hm
    omega
  have hm' : 0 < m / p ^ 2 := Nat.div_pos (Nat.le_of_dvd hm₀ hdiv) (pow_pos hp _)
  rw [minusEleven_filtered_sigma ht (hm.of_dvd_nat (Nat.div_dvd_of_dvd hdiv)), minusEleven_filtered_sigma ht hm]
  have h := pow_mul_sigma_le (k := 1) (pow_pos hp 2) hm'
  rw [Nat.mul_div_cancel' hdiv, pow_one] at h
  have hR : (p : ℝ) ^ 2 * (σ 1 (m / p ^ 2) : ℝ) ≤ σ 1 m := by exact_mod_cast h
  calc
    _ = ((3 / 2 : ℝ) * (2 : ℝ) ^ t) * ((p : ℝ) ^ 2 * (σ 1 (m / p ^ 2) : ℝ)) := by ring
    _ ≤ _ := mul_le_mul_of_nonneg_left hR (by positivity)

private theorem cusp_error_scale_le {p t m : ℕ} (hp : 0 < p) (hm : 0 < m) (hdiv : p ^ 2 ∣ m) :
    (p : ℝ) * minusElevenCuspErrorBound t (m / p ^ 2) ≤ minusElevenCuspErrorBound t m := by
  let m' := m / p ^ 2
  have hm' : 0 < m' := Nat.div_pos (Nat.le_of_dvd hm hdiv) (pow_pos hp _)
  have hmEq : p ^ 2 * m' = m := Nat.mul_div_cancel' hdiv
  have h := pow_mul_sigma_le (k := 0) (pow_pos hp 2) hm'
  rw [hmEq, pow_zero, one_mul] at h
  have hτ : (σ 0 m' : ℝ) ≤ σ 0 m := by exact_mod_cast h
  have hs : Real.sqrt (m : ℝ) = (p : ℝ) * Real.sqrt m' := by
    have hmR : (m : ℝ) = (p : ℝ) ^ 2 * (m' : ℝ) := by exact_mod_cast hmEq.symm
    rw [hmR, Real.sqrt_mul (by positivity), Real.sqrt_sq (by positivity)]
  let c := Real.sqrt ((2 : ℝ) ^ t) * ((t : ℝ) + 3) / 2
  change (p : ℝ) * minusElevenCuspErrorBound t m' ≤ minusElevenCuspErrorBound t m
  calc
    _ = (Real.sqrt (m : ℝ) * c) * (σ 0 m' : ℝ) := by
      rw [hs]
      dsimp [c, minusElevenCuspErrorBound]
      ring
    _ ≤ (Real.sqrt (m : ℝ) * c) * (σ 0 m : ℝ) := mul_le_mul_of_nonneg_left hτ (by positivity)
    _ = _ := by dsimp [c, minusElevenCuspErrorBound]; ring

/-- At an inert prime dividing the target, the sharper small-prime
bound follows from the theta formula and the cusp coefficient bound.
The auxiliary Hecke level can therefore be restricted to primes
coprime to the target. -/
theorem minusEleven_inert_divisor_shell_bound (a : ℕ → ℝ)
    (hTheta : MinusElevenThetaFormula a)
    (hcusp : ∀ t m : ℕ, 0 < t → Odd m → ¬ 11 ∣ 2 ^ t * m →
      |minusElevenFilteredCusp a (2 ^ t * m)| ≤ minusElevenCuspErrorBound t m)
    {t m p : ℕ} (ht : 0 < t) (hm : Odd m) (h₁₁ : ¬ 11 ∣ 2 ^ t * m)
    (hp : MinusElevenInertPrime p) (hpN : p ∣ 2 ^ t * m) :
    let n := 2 ^ t * m
    let E := minusElevenFilteredError t m
    (1 - E) * (minusElevenBadShell n p).card ≤
      (1 / (p : ℝ) ^ 2 + E / p) * (minusElevenOddShell n).card := by
  let n := 2 ^ t * m
  let E := minusElevenFilteredError t m
  let S := minusElevenFilteredSigma n
  let D := minusElevenCuspErrorBound t m
  let G : ℝ := (minusElevenOddShell n).card
  let c := 1 / (p : ℝ) ^ 2 + E / p
  have hscale : 0 ≤ 1 - E := sub_nonneg.mpr (minusEleven_filtered_error_lt_one (m := m) ht).le
  have hE : 0 ≤ E := minusEleven_filtered_error_nonneg t m
  have hpR : 0 < (p : ℝ) := by exact_mod_cast hp.1.pos
  have hc : 0 ≤ c := by dsimp [c]; positivity
  change _ * ((minusElevenBadShell n p).card : ℝ) ≤ c * G
  rw [minusEleven_bad_shell_card_of_dvd hp hpN]
  by_cases hd : p ^ 2 ∣ n
  · rw [if_pos hd]
    have hcp : (p ^ 2).Coprime (2 ^ t) :=
      ((Nat.coprime_primes hp.1 Nat.prime_two).mpr (by have := hp.2.1; omega)).pow 2 t
    have hdm : p ^ 2 ∣ m := hcp.dvd_of_dvd_mul_left hd
    let m' := m / p ^ 2
    let j := 2 ^ t * m'
    have hm' : Odd m' := hm.of_dvd_nat (Nat.div_dvd_of_dvd hdm)
    have hnEq : p ^ 2 * j = n := by
      change p ^ 2 * (2 ^ t * m') = 2 ^ t * m
      rw [mul_left_comm, Nat.mul_div_cancel' hdm]
    have hj : n / p ^ 2 = j := by
      rw [← hnEq, Nat.mul_div_cancel_left _ (pow_pos hp.1.pos _)]
    have hj₁₁ : ¬ 11 ∣ j := fun h => h₁₁ (h.trans ⟨p ^ 2, hnEq.symm.trans (mul_comm ..)⟩)
    have hGold := (minusEleven_odd_shell_bounds_of_theta a hTheta ht hm h₁₁ (hcusp t m ht hm h₁₁)).1
    have hGnew := (minusEleven_odd_shell_bounds_of_theta a hTheta ht hm' hj₁₁
      (hcusp t m' ht hm' hj₁₁)).2
    have hS : minusElevenFilteredSigma j ≤ S / (p : ℝ) ^ 2 := by
      apply (le_div_iff₀ (by positivity)).mpr
      rw [mul_comm (minusElevenFilteredSigma j)]
      exact filtered_sigma_scale_le hp.1.pos ht hm hdm
    have hm₀ : 0 < m := by obtain ⟨k, hk⟩ := hm; omega
    have hD : minusElevenCuspErrorBound t m' ≤ D / (p : ℝ) := by
      apply (le_div_iff₀ hpR).mpr
      simpa only [mul_comm] using cusp_error_scale_le hp.1.pos hm₀ hdm (t := t)
    have hSD : (3 * S) * E = 2 * D := minusEleven_filtered_error_mul ht hm
    rw [hj]
    calc
      _ ≤ (1 - E) * ((12 / 5) * minusElevenFilteredSigma j + (8 / 5) * minusElevenCuspErrorBound t m') :=
        mul_le_mul_of_nonneg_left hGnew hscale
      _ ≤ (1 - E) * ((12 / 5) * (S / (p : ℝ) ^ 2) + (8 / 5) * (D / (p : ℝ))) :=
        mul_le_mul_of_nonneg_left (add_le_add (mul_le_mul_of_nonneg_left hS (by norm_num))
          (mul_le_mul_of_nonneg_left hD (by norm_num))) hscale
      _ = c * ((12 / 5) * S * (1 - E)) := by
        dsimp [c]
        have he : (8 / 5 : ℝ) * D = (12 / 5) * S * E := by nlinarith only [hSD]
        rw [← mul_div_assoc, ← mul_div_assoc, he]
        ring
      _ ≤ _ := mul_le_mul_of_nonneg_left hGold hc
  · rw [if_neg hd, Nat.cast_zero, mul_zero]
    exact mul_nonneg hc (by dsimp [G]; positivity)

theorem minusEleven_inert_divisor_small_bound (a : ℕ → ℝ)
    (hTheta : MinusElevenThetaFormula a)
    (hcusp : ∀ t m : ℕ, 0 < t → Odd m → ¬ 11 ∣ 2 ^ t * m →
      |minusElevenFilteredCusp a (2 ^ t * m)| ≤ minusElevenCuspErrorBound t m)
    {t m p : ℕ} (ht : 0 < t) (hm : Odd m) (h₁₁ : ¬ 11 ∣ 2 ^ t * m)
    (hp : MinusElevenInertPrime p) (hpN : p ∣ 2 ^ t * m) :
    let n := 2 ^ t * m
    let E := minusElevenFilteredError t m
    (1 - E) * (minusElevenBadShell n p).card ≤
      (minusElevenLocalMass p + (5 / 2) * E) * (minusElevenOddShell n).card := by
  let E := minusElevenFilteredError t m
  have hpR : 2 < (p : ℝ) := by exact_mod_cast hp.2.1
  have hE : 0 ≤ E := minusEleven_filtered_error_nonneg t m
  have hm₁ : 1 / (p : ℝ) ^ 2 ≤ minusElevenLocalMass p := by
    unfold minusElevenLocalMass
    have hpOne : 0 < (p : ℝ) - 1 := by linarith
    apply (div_le_div_iff₀ (by positivity) (by positivity)).mpr
    nlinarith only [hpR, sq_nonneg (p : ℝ)]
  have hm₂ : E / (p : ℝ) ≤ (5 / 2) * E := by
    calc
      _ ≤ E := div_le_self hE (by linarith)
      _ ≤ _ := by linarith
  exact (minusEleven_inert_divisor_shell_bound a hTheta hcusp ht hm h₁₁ hp hpN).trans
    (mul_le_mul_of_nonneg_right (add_le_add hm₁ hm₂) (by positivity))

end QuadraticNorms
