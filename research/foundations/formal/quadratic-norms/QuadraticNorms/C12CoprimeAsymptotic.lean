import QuadraticNorms.C12CoprimeTransfer
import QuadraticNorms.PointwiseDivisorGrowth
import Mathlib.Analysis.SpecialFunctions.Pow.Asymptotics

open scoped Classical

namespace QuadraticNorms

private theorem c12_coprime_main_term_pos {A : ℝ} (hA : 0 < A) {n : ℕ} (hn : 1 < n) :
    0 < c12MainTerm A n := by
  have hJ : 0 < pointwiseTargetFactor c12ExceptionalPrime n :=
    zero_lt_one.trans_le (pointwise_target_factor_ge_one (fun p hp => by
      have := hp.1.two_le
      have := hp.2.1
      omega) n)
  have hnR : (1 : ℝ) < n := by exact_mod_cast hn
  unfold c12MainTerm
  exact div_pos (mul_pos (mul_pos hA hJ) (zero_lt_one.trans hnR))
    (Real.rpow_pos_of_pos (Real.log_pos hnR) _)

private theorem c12_main_term_normalization_le {A : ℝ} (hA : 0 < A) {n : ℕ} (hn : 1 < n) :
    (n : ℝ) / c12MainTerm A n ≤ (Real.log n) ^ (3 / 2 : ℝ) / A := by
  let J := pointwiseTargetFactor c12ExceptionalPrime n
  have hJ : 1 ≤ J := pointwise_target_factor_ge_one (fun p hp => by
    have := hp.1.two_le
    have := hp.2.1
    omega) n
  have hnR : (1 : ℝ) < n := by exact_mod_cast hn
  have hn₀ : (0 : ℝ) < n := zero_lt_one.trans hnR
  have hW : 0 < (Real.log n) ^ (3 / 2 : ℝ) := Real.rpow_pos_of_pos (Real.log_pos hnR) _
  have hP : 0 < J := zero_lt_one.trans_le hJ
  calc
    (n : ℝ) / c12MainTerm A n = (Real.log n) ^ (3 / 2 : ℝ) / (A * J) := by
      unfold c12MainTerm
      change (n : ℝ) / (A * J * n / (Real.log n) ^ (3 / 2 : ℝ)) = _
      field_simp [hA.ne', hn₀.ne', hP.ne', hW.ne']
    _ ≤ _ := div_le_div_of_nonneg_left hW.le hA (le_mul_of_one_le_right hA.le hJ)

private theorem c12_relative_error_budget {A C δ D : ℝ}
    (hA : 0 < A) (hC : 0 ≤ C) (hδ₁ : δ ≤ 1)
    (hD : 0 < D) {n : ℕ} (hn : 1 < n) (hL : 2 ≤ Real.log n)
    (hS₀ : c12AbsoluteDivisorMass n ≤ D * (Real.log n) ^ (δ / 4))
    (hS₁ : c12DivisorLogMoment n ≤ D * (Real.log n) ^ (δ / 4))
    (hPinv : (c12CoprimeFactor n)⁻¹ ≤ D * (Real.log n) ^ (δ / 4)) :
    (((8 * C / (Real.log n) ^ δ) * c12AbsoluteDivisorMass n +
        (6 / Real.log n) * c12DivisorLogMoment n +
        (1 + (n : ℝ) / c12MainTerm A n) * ((n.divisors.card : ℝ) / Real.sqrt n)) *
      (c12CoprimeFactor n)⁻¹) ≤
      (8 * C * D ^ 2 + 6 * D ^ 2 + (1 + A⁻¹) * ((7 : ℝ) ^ 16 * 16 ^ 3) * D) /
        (Real.log n) ^ (δ / 2) := by
  let L := Real.log n
  let g := δ / 4
  let r := δ / 2
  let P := c12CoprimeFactor n
  let S₀ := c12AbsoluteDivisorMass n
  let S₁ := c12DivisorLogMoment n
  let B := c12MainTerm A n
  let H := (7 : ℝ) ^ 16 * 16 ^ 3
  have hL₀ : 0 < L := by linarith
  have hL₁ : 1 ≤ L := by linarith
  have hP : 0 < P := c12_coprime_factor_pos n
  have hLg : 0 < L ^ g := Real.rpow_pos_of_pos hL₀ _
  have hLr : 0 < L ^ r := Real.rpow_pos_of_pos hL₀ _
  have hLδ : 0 < L ^ δ := Real.rpow_pos_of_pos hL₀ _
  have hg : L ^ g * L ^ g = L ^ r := by
    rw [← Real.rpow_add hL₀, show g + g = r by unfold g r; ring]
  have hr : L ^ r * L ^ r = L ^ δ := by
    rw [← Real.rpow_add hL₀, show r + r = δ by unfold r; ring]
  have hs₀ : 0 ≤ S₀ := by unfold S₀ c12AbsoluteDivisorMass; positivity
  have hs₁ : 0 ≤ S₁ := by
    unfold S₁ c12DivisorLogMoment
    apply Finset.sum_nonneg
    intro d _
    have hd : 0 ≤ Real.log (d : ℝ) := Real.log_natCast_nonneg d
    positivity
  have hfirst : ((8 * C / L ^ δ) * S₀) * P⁻¹ ≤ (8 * C * D ^ 2) / L ^ r := by
    calc
      _ ≤ ((8 * C / L ^ δ) * (D * L ^ g)) * (D * L ^ g) :=
        mul_le_mul (mul_le_mul_of_nonneg_left hS₀ (by positivity)) hPinv (by positivity) (by positivity)
      _ = _ := by
        rw [← hr]
        field_simp
        rw [pow_two (L ^ g), hg]
  have hsecond : ((6 / L) * S₁) * P⁻¹ ≤ (6 * D ^ 2) / L ^ r := by
    have hpowers : L ^ r / L ≤ 1 / L ^ r := by
      apply (div_le_div_iff₀ hL₀ hLr).mpr
      rw [hr, one_mul]
      simpa only [Real.rpow_one] using Real.rpow_le_rpow_of_exponent_le hL₁ hδ₁
    calc
      _ ≤ ((6 / L) * (D * L ^ g)) * (D * L ^ g) :=
        mul_le_mul (mul_le_mul_of_nonneg_left hS₁ (by positivity)) hPinv (by positivity) (by positivity)
      _ = (6 * D ^ 2) * (L ^ r / L) := by rw [← hg]; ring
      _ ≤ _ := by
        simpa only [mul_one_div] using mul_le_mul_of_nonneg_left hpowers (by positivity : 0 ≤ 6 * D ^ 2)
  have hthird : ((1 + (n : ℝ) / B) * ((n.divisors.card : ℝ) / Real.sqrt n)) * P⁻¹ ≤
      (1 + A⁻¹) * H * D / L ^ r := by
    have hnorm := c12_main_term_normalization_le hA hn
    change (n : ℝ) / B ≤ L ^ (3 / 2 : ℝ) / A at hnorm
    have hLpow : 1 ≤ L ^ (3 / 2 : ℝ) := by
      simpa only [Real.one_rpow] using Real.rpow_le_rpow (by norm_num : (0 : ℝ) ≤ 1) hL₁
        (by norm_num : (0 : ℝ) ≤ 3 / 2)
    have hnorm' : 1 + (n : ℝ) / B ≤ (1 + A⁻¹) * L ^ (3 / 2 : ℝ) := by
      calc
        _ ≤ L ^ (3 / 2 : ℝ) + L ^ (3 / 2 : ℝ) / A := add_le_add hLpow hnorm
        _ = _ := by ring
    have hτ := divisor_count_log_cube_over_sqrt hn
    change (n.divisors.card : ℝ) / Real.sqrt n ≤ H / L ^ 3 at hτ
    have hpowers : L ^ ((3 / 2 : ℝ) + g) / L ^ (3 : ℕ) ≤ 1 / L ^ r := by
      apply (div_le_div_iff₀ (pow_pos hL₀ _) hLr).mpr
      rw [one_mul, ← Real.rpow_add hL₀]
      have h : (3 / 2 : ℝ) + g + r ≤ (3 : ℕ) := by
        push_cast
        dsimp [g, r]
        linarith
      exact (Real.rpow_le_rpow_of_exponent_le hL₁ h).trans_eq (Real.rpow_natCast _ _)
    calc
      _ ≤ (((1 + A⁻¹) * L ^ (3 / 2 : ℝ)) * (H / L ^ 3)) * (D * L ^ g) := by
        apply mul_le_mul _ hPinv (by positivity) (by positivity)
        exact mul_le_mul hnorm' hτ (by positivity) (by positivity)
      _ = ((1 + A⁻¹) * H * D) * (L ^ ((3 / 2 : ℝ) + g) / L ^ (3 : ℕ)) := by
        rw [Real.rpow_add hL₀]
        ring
      _ ≤ _ := by
        simpa only [mul_one_div] using mul_le_mul_of_nonneg_left hpowers
          (by positivity : 0 ≤ (1 + A⁻¹) * H * D)
  change (((8 * C / L ^ δ) * S₀ + (6 / L) * S₁ +
    (1 + (n : ℝ) / B) * ((n.divisors.card : ℝ) / Real.sqrt n)) * P⁻¹) ≤ _
  rw [add_mul, add_mul]
  simpa only [add_div] using add_le_add (add_le_add hfirst hsecond) hthird

/-- The coprime corollary follows from the ordinary C₁₂ estimate
with the stated loss of half the logarithmic exponent. The divisor
growth and the large-divisor remainder require no extra hypotheses. -/
theorem c12_coprime_asymptotic_transfer {A C δ : ℝ} (hA : 0 < A) (hC : 0 ≤ C)
    (hδ : 0 < δ) (hδ₁ : δ ≤ 1) {N₀ : ℕ} (hN₀ : 2 ≤ N₀)
    (hCount : ∀ m ≥ N₀, m % 12 = 2 →
      |(c12OrderedPairs m : ℝ) - c12MainTerm A m| ≤ c12MainTerm A m * C / (Real.log m) ^ δ) :
    ∃ K : ℝ, 0 < K ∧ ∃ N : ℕ, ∀ n ≥ N, n % 12 = 2 →
      |(c12CoprimePairs n : ℝ) - c12MainTerm A n * c12CoprimeFactor n| ≤
        c12MainTerm A n * c12CoprimeFactor n * K / (Real.log n) ^ (δ / 2) := by
  obtain ⟨D, hD, hGrowth⟩ := c12_divisor_growth (ε := δ / 4) (by positivity)
  obtain ⟨N₁, hN₁⟩ := exists_nat_ge (Real.exp 2)
  let K := 8 * C * D ^ 2 + 6 * D ^ 2 + (1 + A⁻¹) * ((7 : ℝ) ^ 16 * 16 ^ 3) * D
  refine ⟨K, by dsimp [K]; positivity, max (N₀ ^ 2) N₁, fun n hn hclass => ?_⟩
  have hn₀ : N₀ ^ 2 ≤ n := (le_max_left ..).trans hn
  have hn₁ : 1 < n := by
    have h := (Nat.pow_le_pow_left hN₀ 2).trans hn₀
    omega
  have hLn : 2 ≤ Real.log n := by
    have hle : Real.exp 2 ≤ (n : ℝ) := hN₁.trans
      (by exact_mod_cast (le_max_right (N₀ ^ 2) N₁).trans hn)
    simpa only [Real.log_exp] using Real.log_le_log (Real.exp_pos _) hle
  obtain ⟨hS₀, hS₁, hPInv⟩ := hGrowth n hLn
  have hBudget := c12_relative_error_budget hA hC hδ₁ hD hn₁ hLn hS₀ hS₁ hPInv
  let B := c12MainTerm A n
  let P := c12CoprimeFactor n
  let L := Real.log n
  let U := (8 * C / L ^ δ) * c12AbsoluteDivisorMass n + (6 / L) * c12DivisorLogMoment n +
    (1 + (n : ℝ) / B) * ((n.divisors.card : ℝ) / Real.sqrt n)
  have hB : 0 < B := c12_coprime_main_term_pos hA hn₁
  have hP : 0 < P := c12_coprime_factor_pos n
  have hU : U * P⁻¹ ≤ K / L ^ (δ / 2) := hBudget
  calc
    _ ≤ B * U := by
      have h := c12_coprime_quantitative_transfer hA.le hC hδ₁ hN₀ hn₀ hclass hCount
      refine h.trans_eq ?_
      have hscale : B * (1 + (n : ℝ) / B) = B + n := by
        rw [mul_add, mul_one, mul_div_cancel₀ (n : ℝ) hB.ne']
      calc
        _ = B * ((8 * C / L ^ δ) * c12AbsoluteDivisorMass n +
            (6 / L) * c12DivisorLogMoment n) +
            (B + (n : ℝ)) * ((n.divisors.card : ℝ) / Real.sqrt n) := by
          rw [mul_div_assoc (c12MainTerm A n + (n : ℝ))]
        _ = B * U := by
          dsimp only [U]
          rw [← hscale, mul_assoc B (1 + (n : ℝ) / B), ← mul_add]
    _ = B * P * (U * P⁻¹) := by field_simp [hP.ne']
    _ ≤ _ := by
      rw [mul_div_assoc]
      exact mul_le_mul_of_nonneg_left hU (mul_nonneg hB.le hP.le)

theorem c12_coprime_pairs_pos_iff {n : ℕ} :
    0 < c12CoprimePairs n ↔ ∃ a b : ℕ, C12 a ∧ C12 b ∧ a + b = n ∧ a.Coprime b := by
  unfold c12CoprimePairs coprimePairFinset
  rw [Finset.card_pos]
  constructor
  · rintro ⟨a, ha⟩
    obtain ⟨hpair, hcoprime⟩ := Finset.mem_filter.mp ha
    obtain ⟨hle, ha₁, ha₂⟩ := mem_pairFinset.mp hpair
    exact ⟨a, n - a, ha₁, ha₂, by omega, Nat.coprime_iff_gcd_eq_one.mpr hcoprime⟩
  · rintro ⟨a, b, ha, hb, hab, hcoprime⟩
    have hsub : n - a = b := by omega
    refine ⟨a, Finset.mem_filter.mpr ⟨mem_pairFinset.mpr ⟨by omega, ha, hsub ▸ hb⟩, ?_⟩⟩
    rw [hsub]
    exact Nat.coprime_iff_gcd_eq_one.mp hcoprime

theorem c12_coprime_pairs_eq_zero_of_inadmissible {n : ℕ} (hn : n % 12 ≠ 2) :
    c12CoprimePairs n = 0 := by
  by_contra h
  obtain ⟨a, b, ha, hb, hab, _⟩ := c12_coprime_pairs_pos_iff.mp (Nat.pos_of_ne_zero h)
  exact hn (C12_sum_mod_twelve ⟨a, b, ha, hb, hab⟩)

/-- The ordinary estimate also gives coprime summands for every sufficiently large
admissible target. -/
theorem c12_eventually_coprime_summands {A C δ : ℝ} (hA : 0 < A) (hC : 0 ≤ C)
    (hδ : 0 < δ) (hδ₁ : δ ≤ 1) {N₀ : ℕ} (hN₀ : 2 ≤ N₀)
    (hCount : ∀ m ≥ N₀, m % 12 = 2 →
      |(c12OrderedPairs m : ℝ) - c12MainTerm A m| ≤ c12MainTerm A m * C / (Real.log m) ^ δ) :
    ∃ N : ℕ, ∀ n ≥ N, n % 12 = 2 →
      ∃ a b : ℕ, C12 a ∧ C12 b ∧ a + b = n ∧ a.Coprime b := by
  obtain ⟨K, hK, N₁, hbound⟩ := c12_coprime_asymptotic_transfer hA hC hδ hδ₁ hN₀ hCount
  have htrend : Filter.Tendsto (fun n : ℕ => (Real.log n) ^ (δ / 2)) Filter.atTop Filter.atTop :=
    (tendsto_rpow_atTop (by positivity : 0 < δ / 2)).comp
      (Real.tendsto_log_atTop.comp tendsto_natCast_atTop_atTop)
  obtain ⟨N₂, hN₂⟩ := Filter.eventually_atTop.mp (htrend.eventually (Filter.eventually_gt_atTop K))
  refine ⟨max (max N₁ N₂) 2, fun n hn hclass => c12_coprime_pairs_pos_iff.mp ?_⟩
  have hcommon : max N₁ N₂ ≤ n := (le_max_left ..).trans hn
  have hn₁ : 1 < n := (le_max_right ..).trans hn
  have hKpow : K < (Real.log n) ^ (δ / 2) := hN₂ n ((le_max_right ..).trans hcommon)
  have hpos := mul_pos (c12_coprime_main_term_pos hA hn₁) (c12_coprime_factor_pos n)
  have hr : c12MainTerm A n * c12CoprimeFactor n * K / (Real.log n) ^ (δ / 2) <
      c12MainTerm A n * c12CoprimeFactor n := by
    rw [div_lt_iff₀ (hK.trans hKpow)]
    exact mul_lt_mul_of_pos_left hKpow hpos
  have h := (hbound n ((le_max_left ..).trans hcommon) hclass).trans_lt hr
  have hpair : 0 < (c12CoprimePairs n : ℝ) := by
    have := (abs_lt.mp h).1
    linarith
  exact_mod_cast hpair

end QuadraticNorms
