import QuadraticNorms.MinusElevenEulerBounds
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Data.Nat.Log

namespace QuadraticNorms

/-- The largest split-divisor product on positive integers up to X. -/
def minusElevenDivisorMaximum (X : ℕ) : ℕ :=
  (Finset.Icc 1 X).sup minusElevenSplitDivisors

theorem minusEleven_split_divisors_le_maximum {v X : ℕ} (hv : 0 < v) (hle : v ≤ X) :
    minusElevenSplitDivisors v ≤ minusElevenDivisorMaximum X :=
  Finset.le_sup (Finset.mem_Icc.mpr ⟨hv, hle⟩)

private theorem minusEleven_maximum_le_bin {row : MinusElevenBin} (hr : row ∈ minusElevenBins)
    {X : ℕ} (hX : X < 2 ^ (row.exponent + 1)) :
    minusElevenDivisorMaximum X ≤ row.divisorLimit := by
  apply Finset.sup_le
  intro v hv
  exact minusEleven_bin_divisor_limit hr (Finset.mem_Icc.mp hv).1 ((Finset.mem_Icc.mp hv).2.trans_lt hX)

private theorem minusEleven_maximum_le_scale {X : ℕ} {y : ℝ} (_hy : 0 ≤ y) (hX : (X : ℝ) ≤ y ^ 12) :
    (minusElevenDivisorMaximum X : ℝ) ≤ 7 * y ^ 2 := by
  rw [← Nat.le_floor_iff (by positivity : 0 ≤ 7 * y ^ 2)]
  apply Finset.sup_le
  intro v hv
  apply Nat.le_floor
  have ⟨hv₀, hv₁⟩ := Finset.mem_Icc.mp hv
  have hD : (minusElevenSplitDivisors v : ℝ) ^ 6 < (7 : ℝ) ^ 6 * v := by
    exact_mod_cast minusEleven_split_divisors_power_bound hv₀
  have hvR : (v : ℝ) ≤ X := by exact_mod_cast hv₁
  have hpow : (minusElevenSplitDivisors v : ℝ) ^ 6 ≤ (7 * y ^ 2) ^ 6 := by
    rw [mul_pow, ← pow_mul]
    exact hD.le.trans (mul_le_mul_of_nonneg_left (hvR.trans hX) (by positivity))
  exact (pow_le_pow_iff_left₀ (by positivity) (by positivity) (by decide : (6 : ℕ) ≠ 0)).mp hpow

private theorem error_le_scale {n E z : ℝ} (hE : 0 ≤ E) (hpow : E ^ 6 * n ^ 2 ≤ (11 / 4) ^ 6)
    (hz : 0 < z) (hzn : z ^ 3 ≤ n) : E ≤ 11 / (4 * z) := by
  have hleft : (E * z) ^ 6 ≤ (11 / 4 : ℝ) ^ 6 := by
    rw [mul_pow, show z ^ 6 = (z ^ 3) ^ 2 by rw [← pow_mul]]
    exact (mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (by positivity) hzn _) (by positivity)).trans hpow
  have := (pow_le_pow_iff_left₀ (by positivity) (by norm_num) (by decide : (6 : ℕ) ≠ 0)).mp hleft
  rw [← div_div]
  exact (le_div_iff₀ hz).mpr this

/-- The sufficient sieve expression from the −11 proof. A cutoff is
recorded by its floor, since all auxiliary primes are integers. -/
noncomputable def minusElevenSieveExpression (n P : ℕ) (E : ℝ) : ℝ :=
  minusElevenInertMass + (8 / 3) * minusElevenDivisorMaximum n * minusElevenInertTail P +
    E * (3 + (5 / 2) * minusElevenInertCount P)

private theorem sieve_expression_le {n P : ℕ} {E D T U A : ℝ}
    (hE : 0 ≤ E) (hD : (minusElevenDivisorMaximum n : ℝ) ≤ D)
    (hT : minusElevenInertTail P ≤ T) (hU : E ≤ U)
    (hA : (minusElevenInertCount P : ℝ) ≤ A) :
    minusElevenSieveExpression n P E ≤ 9 / 100 + (8 / 3) * D * T + U * (3 + (5 / 2) * A) := by
  unfold minusElevenSieveExpression
  have hD₀ : 0 ≤ D := (Nat.cast_nonneg _).trans hD
  have hpart₁ := mul_le_mul (mul_le_mul_of_nonneg_left hD (by norm_num : (0 : ℝ) ≤ 8 / 3))
    hT (minusEleven_inert_tail_nonneg _) (mul_nonneg (by norm_num) hD₀)
  have hpart₂ := mul_le_mul hU
    (add_le_add_right (mul_le_mul_of_nonneg_left hA (by norm_num : (0 : ℝ) ≤ 5 / 2)) 3)
    (by positivity) (hE.trans hU)
  exact add_le_add (add_le_add minusEleven_inert_mass_lt.le hpart₁) hpart₂

/-- Every inert-core target in a certified dyadic bin satisfies the
sieve inequality at that bin's cutoff. -/
theorem minusEleven_sieve_bin {row : MinusElevenBin} (hr : row ∈ minusElevenBins) (hv : row.valid)
    {n : ℕ} (hlo : 2 ^ row.exponent ≤ n) (hhi : n < 2 ^ (row.exponent + 1))
    {E : ℝ} (hE : 0 ≤ E) (hpow : E ^ 6 * (n : ℝ) ^ 2 ≤ (11 / 4 : ℝ) ^ 6) :
    minusElevenSieveExpression n row.primeCutoff E < 2 / 3 := by
  have hD : (minusElevenDivisorMaximum n : ℝ) ≤ row.divisorLimit := by
    exact_mod_cast minusEleven_maximum_le_bin hr hhi
  have hc : (0 : ℝ) < row.cubeRoot := by exact_mod_cast hv.1
  have hcube : (row.cubeRoot : ℝ) ^ 3 ≤ n := by exact_mod_cast hv.2.1.trans hlo
  have hE' : E ≤ (row.errorBound : ℝ) := by
    simpa only [MinusElevenBin.errorBound, Rat.cast_div, Rat.cast_ofNat, Rat.cast_mul, Rat.cast_natCast]
      using error_le_scale hE hpow hc hcube
  have hm := sieve_expression_le hE hD (MinusElevenBin.inertTail_le hv) hE'
    (A := row.inertCount) (by rw [MinusElevenBin.inertCount_eq hv])
  have hmargin := (Rat.cast_lt (K := ℝ)).mpr hv.2.2.2.2.2.2
  unfold MinusElevenBin.margin at hmargin
  push_cast at hmargin
  exact hm.trans_lt (hmargin.trans_le (by norm_num))

private theorem asymptotic_rational_margin {y : ℝ} (hy : 32 ≤ y) :
    9 / 100 + (8 / 3) * (7 * y ^ 2) * (5 / (22 * y ^ 3) + 5 / (y ^ 3) ^ 2) +
      (11 / (4 * y ^ 4)) * (3 + (5 / 2) * ((5 / 22) * y ^ 3 + 3 / 2)) < 2 / 3 := by
  have hy₀ : 0 < y := by linarith
  have hcalc :
      (8 / 3) * (7 * y ^ 2) * (5 / (22 * y ^ 3) + 5 / (y ^ 3) ^ 2) +
      (11 / (4 * y ^ 4)) * (3 + (5 / 2) * ((5 / 22) * y ^ 3 + 3 / 2)) =
      (3065 / 528) * (1 / y) + (5371 / 48) * (1 / y ^ 4) := by
    field_simp [hy₀.ne']
    ring
  rw [add_assoc, hcalc]
  have h₁ := one_div_le_one_div_of_le (by norm_num : (0 : ℝ) < 32) hy
  have h₂ := one_div_le_one_div_of_le (by norm_num : (0 : ℝ) < (32 : ℝ) ^ 4)
    (pow_le_pow_left₀ (by norm_num) hy 4)
  norm_num at h₁ h₂
  simp only [one_div]
  linarith only [h₁, h₂]

/-- The uniform large-target sieve margin at the floor of the fourth
root, using the sixth-power divisor and shell error bounds. -/
theorem minusEleven_sieve_large {n : ℕ} (hn : 2 ^ 60 ≤ n)
    {E : ℝ} (hE : 0 ≤ E) (hpow : E ^ 6 * (n : ℝ) ^ 2 ≤ (11 / 4 : ℝ) ^ 6) :
    ∃ P : ℕ, 0 < P ∧ minusElevenSieveExpression n P E < 2 / 3 := by
  let y := (n : ℝ) ^ (1 / 12 : ℝ)
  have hy12 : y ^ 12 = (n : ℝ) := by
    dsimp [y]
    rw [← Real.rpow_mul_natCast (Nat.cast_nonneg _) _ 12]
    norm_num
  have hy₀ : 0 ≤ y := Real.rpow_nonneg (Nat.cast_nonneg _) _
  have hy : 32 ≤ y := by
    apply (pow_le_pow_iff_left₀ (by norm_num) hy₀ (by decide : (12 : ℕ) ≠ 0)).mp
    rw [hy12]
    have hnR : (2 : ℝ) ^ 60 ≤ n := by exact_mod_cast hn
    norm_num at hnR ⊢
    exact hnR
  have hD := minusEleven_maximum_le_scale hy₀ hy12.ge
  have hE' : E ≤ 11 / (4 * y ^ 4) := by
    apply error_le_scale hE hpow
    · positivity
    · rw [← pow_mul]
      exact hy12.le
  let P := ⌊y ^ 3⌋₊
  have hP : 0 < P := by
    apply Nat.le_floor
    norm_num only [Nat.cast_one]
    exact one_le_pow₀ (by linarith)
  have hX : 0 < y ^ 3 := pow_pos (by linarith) _
  refine ⟨P, hP, ?_⟩
  have h := sieve_expression_le hE hD (minusEleven_inert_tail_le_real hX) hE'
    (minusEleven_inert_count_le_real hX.le)
  exact h.trans_lt (asymptotic_rational_margin hy)

/-- For every inert-core target above the finite scan, the exact
arithmetic estimates give a cutoff with a positive sieve margin. -/
theorem minusEleven_sieve_margin_above_scan {t m : ℕ} (hm : 0 < m)
    (hI : ∀ p ∈ m.primeFactors, MinusElevenInertPrime p)
    (hn : 2 ^ 25 ≤ 2 ^ t * m) :
    ∃ P : ℕ, 0 < P ∧ minusElevenSieveExpression (2 ^ t * m) P (minusElevenFilteredError t m) < 2 / 3 := by
  let n := 2 ^ t * m
  have hn₀ : 0 < n := Nat.mul_pos (by positivity) hm
  have hE := minusEleven_filtered_error_nonneg t m
  have hpow : minusElevenFilteredError t m ^ 6 * (n : ℝ) ^ 2 ≤ (11 / 4 : ℝ) ^ 6 := by
    simpa only [n, Nat.cast_mul, Nat.cast_pow, Nat.cast_ofNat] using
      (minusEleven_filtered_error_power_bound hm hI).le
  by_cases hlast : 2 ^ 60 ≤ n
  · exact minusEleven_sieve_large hlast hE hpow
  · let k := Nat.log 2 n
    have hk : 25 ≤ k := Nat.le_log_of_pow_le (by decide) hn
    have hlow : 2 ^ k ≤ n := Nat.pow_log_le_self _ hn₀.ne'
    have hk' : k < 60 := by
      by_contra! h
      exact hlast ((Nat.pow_le_pow_right (by decide : 0 < 2) h).trans hlow)
    have hhigh : n < 2 ^ (k + 1) := Nat.lt_pow_succ_log_self (by decide) n
    obtain ⟨row, hr, heq, hv⟩ := minusEleven_bin_certificate hk hk'
    rw [← heq] at hlow hhigh
    exact ⟨row.primeCutoff, hv.2.2.2.2.1, minusEleven_sieve_bin hr hv hlow hhigh hE hpow⟩

end QuadraticNorms
