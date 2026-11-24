import QuadraticNorms.MinusElevenThetaNormalization
import QuadraticNorms.HeckeCoefficientBounds

open scoped BigOperators ArithmeticFunction.sigma

namespace QuadraticNorms

/-- Normalization, multiplicativity and the unramified Hecke coefficient
relations used for the weight-two level-eleven cusp term. -/
structure MinusElevenHeckeCoefficients (a : ℕ → ℝ) : Prop where
  normalized : a 1 = 1
  coprime_mul : ∀ {m n : ℕ}, m.Coprime n → a (m * n) = a m * a n
  prime_bound : ∀ p : ℕ, p.Prime → p ≠ 11 → |a p| ≤ 2 * Real.sqrt p
  prime_power : ∀ p : ℕ, p.Prime → p ≠ 11 → ∀ k : ℕ,
    a (p ^ (k + 2)) = a p * a (p ^ (k + 1)) - (p : ℝ) * a (p ^ k)

private theorem MinusElevenHeckeCoefficients.weightTwo {a : ℕ → ℝ}
    (ha : MinusElevenHeckeCoefficients a) : WeightTwoHeckeCoefficients a (· ≠ 11) :=
  ⟨ha.normalized, ha.coprime_mul, ha.prime_bound, ha.prime_power⟩

/-- The unfiltered cusp coefficient bound at an eleven-unit index. -/
theorem minusEleven_cusp_coefficient_bound {a : ℕ → ℝ} (ha : MinusElevenHeckeCoefficients a)
    {m : ℕ} (hm : 0 < m) (h₁₁ : ¬ 11 ∣ m) :
    |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m :=
  weight_two_hecke_coefficient_bound ha.weightTwo hm
    (fun _ _ hd he => h₁₁ (he ▸ hd))

private theorem minusEleven_two_power_filtered_bound {a : ℕ → ℝ}
    (ha : MinusElevenHeckeCoefficients a) (k : ℕ) :
    |a (2 ^ (k + 2)) - a (2 ^ k)| ≤ Real.sqrt ((2 : ℝ) ^ (k + 2)) * ((k : ℝ) + 5) / 2 := by
  let r := Real.sqrt (2 : ℝ)
  have hrSq : r ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  obtain ⟨u, v, hu, hv, huv, hpow⟩ := hecke_prime_power_roots ha.weightTwo Nat.prime_two (by decide)
  have he : ((a (2 ^ (k + 2)) - a (2 ^ k) : ℝ) : ℂ) =
      u ^ (k + 2) + v ^ (k + 2) + heckePowerSum u v k := by
    rw [Complex.ofReal_sub, hpow, hpow, hecke_power_sum_two_endpoints, huv]
    norm_num
    ring
  rw [← Real.norm_eq_abs, ← Complex.norm_real, he, sqrt_nat_power (by norm_num)]
  calc
    _ ≤ (‖u ^ (k + 2)‖ + ‖v ^ (k + 2)‖) + ‖heckePowerSum u v k‖ :=
      (norm_add_le ..).trans (add_le_add (norm_add_le ..) le_rfl)
    _ ≤ (r ^ (k + 2) + r ^ (k + 2)) + ((k : ℝ) + 1) * r ^ k := by
      rw [norm_pow, norm_pow, hu, hv]
      exact add_le_add le_rfl (norm_hecke_power_sum_le hu hv k)
    _ = r ^ (k + 2) * ((k : ℝ) + 5) / 2 := by rw [pow_add, hrSq]; ring

/-- The cancellation at two in the filtered Hecke recurrence yields
the paper's cusp error bound, with the factor (t+3)/2. -/
theorem minusEleven_filtered_cusp_bound {a : ℕ → ℝ} (ha : MinusElevenHeckeCoefficients a)
    {t m : ℕ} (ht : 0 < t) (hm : Odd m) (h₁₁ : ¬ 11 ∣ 2 ^ t * m) :
    |minusElevenFilteredCusp a (2 ^ t * m)| ≤ minusElevenCuspErrorBound t m := by
  have hm₀ : 0 < m := by obtain ⟨k, hk⟩ := hm; omega
  have h₁₁m : ¬ 11 ∣ m := fun h => h₁₁ (h.trans (dvd_mul_left ..))
  have hOdd := minusEleven_cusp_coefficient_bound ha hm₀ h₁₁m
  have hc (k : ℕ) : (2 ^ k).Coprime m :=
    ((Nat.Prime.coprime_iff_not_dvd Nat.prime_two).mpr hm.not_two_dvd_nat).pow_left k
  obtain ⟨j, rfl⟩ := Nat.exists_eq_succ_of_ne_zero ht.ne'
  cases j with
  | zero =>
    have hd : ¬ 4 ∣ 2 ^ (0 + 1) * m := by
      obtain ⟨k, h⟩ := hm
      norm_num
      omega
    change |minusElevenFilteredCusp a (2 ^ (0 + 1) * m)| ≤ minusElevenCuspErrorBound (0 + 1) m
    rw [minusElevenFilteredCusp, if_neg hd, sub_zero, ha.coprime_mul (hc _), abs_mul]
    have h₂ := ha.prime_bound 2 Nat.prime_two (by decide)
    calc
      _ ≤ (2 * Real.sqrt 2) * ((σ 0 m : ℝ) * Real.sqrt m) := by
        norm_num only [zero_add, pow_one]
        exact mul_le_mul h₂ hOdd (abs_nonneg _) (by positivity)
      _ = _ := by norm_num [minusElevenCuspErrorBound]; ring
  | succ k =>
    change |minusElevenFilteredCusp a (2 ^ (k + 2) * m)| ≤ minusElevenCuspErrorBound (k + 2) m
    have hn : 2 ^ (k + 2) * m = 4 * (2 ^ k * m) := by rw [pow_add]; ring
    have hdiv : 2 ^ (k + 2) * m / 4 = 2 ^ k * m := by
      rw [hn, Nat.mul_div_cancel_left _ (by decide : 0 < 4)]
    rw [minusElevenFilteredCusp, if_pos ⟨_, hn⟩, hdiv, ha.coprime_mul (hc _), ha.coprime_mul (hc _),
      ← sub_mul, abs_mul]
    calc
      _ ≤ (Real.sqrt ((2 : ℝ) ^ (k + 2)) * ((k : ℝ) + 5) / 2) * ((σ 0 m : ℝ) * Real.sqrt m) :=
        mul_le_mul (minusEleven_two_power_filtered_bound ha k) hOdd (abs_nonneg _) (by positivity)
      _ = _ := by dsimp [minusElevenCuspErrorBound]; push_cast; ring

end QuadraticNorms
