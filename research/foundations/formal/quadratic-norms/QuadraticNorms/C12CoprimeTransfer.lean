import QuadraticNorms.C12CoprimeWeights
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.Real.Sqrt

open scoped BigOperators ArithmeticFunction.Moebius Classical

namespace QuadraticNorms

/-- Numerical pair count R₁₂ and its coprime restriction. -/
noncomputable def c12OrderedPairs (n : ℕ) : ℕ := pairCount C12 C12 n
noncomputable def c12CoprimePairs (n : ℕ) : ℕ :=
  (coprimePairFinset (fun p => p % 12 = 1) n).card

theorem c12_ordered_pairs_le (n : ℕ) : c12OrderedPairs n ≤ n := by
  unfold c12OrderedPairs pairCount
  have h : pairFinset C12 C12 n ⊆ Finset.Icc 1 n := by
    intro a ha
    obtain ⟨hle, hc, _⟩ := mem_pairFinset.mp ha
    exact Finset.mem_Icc.mpr ⟨hc.1, hle⟩
  simpa only [Nat.card_Icc, Nat.add_sub_cancel] using Finset.card_le_card h

/-- The pointwise main term before imposing coprimality. -/
noncomputable def c12MainTerm (A : ℝ) (n : ℕ) : ℝ :=
  A * pointwiseTargetFactor c12ExceptionalPrime n * n / (Real.log n) ^ (3 / 2 : ℝ)

private theorem c12_main_term_nonneg {A : ℝ} (hA : 0 ≤ A) (n : ℕ) : 0 ≤ c12MainTerm A n := by
  have hJ : 1 ≤ pointwiseTargetFactor c12ExceptionalPrime n :=
    pointwise_target_factor_ge_one (fun p hp => by
      have := hp.1.two_le
      have := hp.2.1
      omega) n
  unfold c12MainTerm
  positivity

private theorem log_ratio_estimate {L l : ℝ} (hl : 0 < l) (hle : l ≤ L) (hhalf : L ≤ 2 * l) :
    0 ≤ (L / l) ^ (3 / 2 : ℝ) - 1 ∧
      (L / l) ^ (3 / 2 : ℝ) - 1 ≤ 6 * (L - l) / L := by
  let r := L / l
  have hL : 0 < L := hl.trans_le hle
  have hr₁ : 1 ≤ r := (le_div_iff₀ hl).mpr (by simpa using hle)
  have hr₂ : r ≤ 2 := (div_le_iff₀ hl).mpr hhalf
  have hr₀ : 0 < r := zero_lt_one.trans_le hr₁
  have hdiff : (r - 1) * l = L - l := by
    unfold r
    rw [sub_mul, div_mul_cancel₀ _ hl.ne', one_mul]
  have h₁ : 1 ≤ r ^ (3 / 2 : ℝ) := by
    simpa only [Real.one_rpow] using
      Real.rpow_le_rpow (by norm_num : (0 : ℝ) ≤ 1) hr₁ (by norm_num : (0 : ℝ) ≤ 3 / 2)
  have h₂ : r ^ (3 / 2 : ℝ) ≤ r ^ 2 := by
    simpa only [Real.rpow_natCast] using
      Real.rpow_le_rpow_of_exponent_le hr₁ (by norm_num : (3 / 2 : ℝ) ≤ (2 : ℕ))
  refine ⟨sub_nonneg.mpr h₁, ?_⟩
  calc
    (L / l) ^ (3 / 2 : ℝ) - 1 ≤ (r - 1) * 3 := by
      change r ^ (3 / 2 : ℝ) - 1 ≤ _
      nlinarith only [h₂, mul_nonneg (sub_nonneg.mpr hr₁) (sub_nonneg.mpr hr₂)]
    _ ≤ _ := by
      apply (le_div_iff₀ hL).mpr
      linear_combination 3 * (mul_le_mul_of_nonneg_left hhalf (sub_nonneg.mpr hr₁)) + 6 * hdiff

private theorem c12_main_term_rescale {A : ℝ} {n d : ℕ} (hn : 0 < n) (hd : d ∣ n) (hc : C12 d)
    (hl : 0 < Real.log (n / d : ℕ)) (hL : 0 < Real.log n) :
    c12MainTerm A (n / d) =
      (c12MainTerm A n / d) * (Real.log n / Real.log (n / d : ℕ)) ^ (3 / 2 : ℝ) := by
  have hd0 : (0 : ℝ) < d := by exact_mod_cast hc.1
  have hq : ((n / d : ℕ) : ℝ) = (n : ℝ) / d := by
    apply (eq_div_iff hd0.ne').mpr
    rw [mul_comm]
    exact_mod_cast Nat.mul_div_cancel' hd
  unfold c12MainTerm
  rw [(c12_target_factor_div hn hd hc).2, Real.div_rpow hL.le hl.le]
  nth_rw 1 [hq]
  field_simp [ne_of_gt (Real.rpow_pos_of_pos hL (3 / 2 : ℝ)),
    ne_of_gt (Real.rpow_pos_of_pos hl (3 / 2 : ℝ)), hd0.ne']

private theorem c12_small_divisor_error {A C δ : ℝ} (hA : 0 ≤ A) (hC : 0 ≤ C) (hδ : δ ≤ 1)
    {n d : ℕ} (hn : 1 < n) (hd : d ∣ n) (hc : C12 d) (hsmall : d ^ 2 ≤ n)
    (hmain : |(c12OrderedPairs (n / d) : ℝ) - c12MainTerm A (n / d)| ≤
      c12MainTerm A (n / d) * C / (Real.log (n / d : ℕ)) ^ δ) :
    |(c12OrderedPairs (n / d) : ℝ) - c12MainTerm A n / d| ≤
      (c12MainTerm A n / d) * (8 * C / (Real.log n) ^ δ + 6 * Real.log d / Real.log n) := by
  let L := Real.log n
  let l := Real.log (n / d : ℕ)
  let r := L / l
  let B := c12MainTerm A n / (d : ℝ)
  have hd0 : (0 : ℝ) < d := by exact_mod_cast hc.1
  have hq0 : 0 < n / d := Nat.div_pos (Nat.le_of_dvd (by omega) hd) hc.1
  have hqR : (0 : ℝ) < (n / d : ℕ) := by exact_mod_cast hq0
  have hL : 0 < L := Real.log_pos (by exact_mod_cast hn)
  have hlog : L = Real.log d + l := by
    have h : (d : ℝ) * (n / d : ℕ) = n := by exact_mod_cast Nat.mul_div_cancel' hd
    dsimp [L, l]
    rw [← Real.log_mul hd0.ne' hqR.ne', h]
  have hdlog : 0 ≤ Real.log d := Real.log_nonneg (by exact_mod_cast hc.1)
  have hle : l ≤ L := by linarith only [hlog, hdlog]
  have hhalf : L ≤ 2 * l := by
    have hd2 : (d : ℝ) ^ 2 ≤ n := by exact_mod_cast hsmall
    have h := Real.log_le_log (pow_pos hd0 _) hd2
    rw [Real.log_pow] at h
    change (2 : ℝ) * Real.log d ≤ L at h
    linarith only [hlog, h]
  have hl : 0 < l := by linarith only [hhalf, hL]
  have hx₀ : 0 ≤ B := div_nonneg (c12_main_term_nonneg hA _) hd0.le
  have hscale := c12_main_term_rescale (A := A) (show 0 < n by omega) hd hc hl hL
  change c12MainTerm A (n / d) = B * r ^ (3 / 2 : ℝ) at hscale
  have hr₁ : 1 ≤ r := (le_div_iff₀ hl).mpr (by simpa using hle)
  have hr₂ : r ≤ 2 := (div_le_iff₀ hl).mpr hhalf
  have hr₀ : 0 < r := zero_lt_one.trans_le hr₁
  have h₃ : r ^ ((3 / 2 : ℝ) + δ) ≤ 8 := by
    apply (Real.rpow_le_rpow_of_exponent_le hr₁ (show (3 / 2 : ℝ) + δ ≤ (3 : ℕ) by push_cast; linarith)).trans
    rw [Real.rpow_natCast]
    calc
      r ^ 3 ≤ (2 : ℝ) ^ 3 := pow_le_pow_left₀ hr₀.le hr₂ _
      _ = 8 := by norm_num
  have hcount : |(c12OrderedPairs (n / d) : ℝ) - c12MainTerm A (n / d)| ≤ B * (8 * C / L ^ δ) := by
    apply hmain.trans
    rw [hscale]
    change B * r ^ (3 / 2 : ℝ) * C / l ^ δ ≤ _
    have heq : B * r ^ (3 / 2 : ℝ) * C / l ^ δ = (B * C / L ^ δ) * r ^ ((3 / 2 : ℝ) + δ) := by
      rw [Real.rpow_add hr₀]
      conv_rhs => rw [show r ^ δ = L ^ δ / l ^ δ from Real.div_rpow hL.le hl.le δ]
      field_simp [ne_of_gt (Real.rpow_pos_of_pos hL δ), ne_of_gt (Real.rpow_pos_of_pos hl δ)]
    rw [heq]
    calc
      _ ≤ (B * C / L ^ δ) * 8 := mul_le_mul_of_nonneg_left h₃ (by positivity)
      _ = _ := by ring
  have hlin := log_ratio_estimate hl hle hhalf
  have hlogerr : |c12MainTerm A (n / d) - B| ≤ B * (6 * Real.log d / L) := by
    rw [hscale, ← mul_sub_one, abs_mul, abs_of_nonneg hx₀, abs_of_nonneg hlin.1]
    apply mul_le_mul_of_nonneg_left _ hx₀
    simpa only [hlog, add_sub_cancel_right] using hlin.2
  change |(c12OrderedPairs (n / d) : ℝ) - B| ≤ B * (8 * C / L ^ δ + 6 * Real.log d / L)
  rw [mul_add]
  exact (abs_sub_le ..).trans (add_le_add hcount hlogerr)


/-- A uniform error bound for the ordinary pair count transfers to the
coprime pair count. The right side displays the two divisor moments and
the error from divisors larger than the square root of the target. -/
theorem c12_coprime_quantitative_transfer {A C δ : ℝ} (hA : 0 ≤ A) (hC : 0 ≤ C) (hδ : δ ≤ 1)
    {N₀ n : ℕ} (hN₀ : 2 ≤ N₀) (hn : N₀ ^ 2 ≤ n) (hclass : n % 12 = 2)
    (hCount : ∀ m ≥ N₀, m % 12 = 2 →
      |(c12OrderedPairs m : ℝ) - c12MainTerm A m| ≤ c12MainTerm A m * C / (Real.log m) ^ δ) :
    |(c12CoprimePairs n : ℝ) - c12MainTerm A n * c12CoprimeFactor n| ≤
      c12MainTerm A n *
        (8 * C / (Real.log n) ^ δ * c12AbsoluteDivisorMass n +
          (6 / Real.log n) * c12DivisorLogMoment n) +
        (c12MainTerm A n + n) * (n.divisors.card : ℝ) / Real.sqrt n := by
  let divs := c12AllowedDivisors n
  let B := c12MainTerm A n
  let L := Real.log n
  let T := (B + n) / Real.sqrt n
  let c₀ := 8 * C / L ^ δ
  let c₁ := 6 / L
  have hn₄ : 4 ≤ n := (Nat.pow_le_pow_left hN₀ 2).trans hn
  have hn₁ : 1 < n := by omega
  have hn₀ : 0 < n := by omega
  have hnR : (0 : ℝ) < n := by exact_mod_cast hn₀
  have hL : 0 < L := Real.log_pos (by exact_mod_cast hn₁)
  have hB : 0 ≤ B := c12_main_term_nonneg hA n
  have hT : 0 ≤ T := by positivity
  have h₁ (d : ℕ) (hd : d ∈ divs) : d ∣ n ∧ C12 d := by
    obtain ⟨hm, hc⟩ := Finset.mem_filter.mp hd
    exact ⟨Nat.dvd_of_mem_divisors hm, hc⟩
  have hterm (d : ℕ) (hd : d ∈ divs) :
      |(μ d : ℝ)| * |(c12OrderedPairs (n / d) : ℝ) - B / d| ≤
        B * (c₀ * (|(μ d : ℝ)| / d) + c₁ * (|(μ d : ℝ)| * Real.log d / d)) + T := by
    obtain ⟨hdn, hc⟩ := h₁ d hd
    have hd0 : (0 : ℝ) < d := by exact_mod_cast hc.1
    have hdlog : 0 ≤ Real.log d := Real.log_nonneg (by exact_mod_cast hc.1)
    have hC₀ : 0 ≤ c₀ := by positivity
    have hC₁ : 0 ≤ c₁ := by positivity
    have hw : |(μ d : ℝ)| ≤ 1 := by exact_mod_cast ArithmeticFunction.abs_moebius_le_one (n := d)
    by_cases hs : d ^ 2 ≤ n
    · have hq : N₀ ≤ n / d := by
        have hprod := Nat.mul_div_cancel' hdn
        have hdq : d ≤ n / d := by
          rw [pow_two, ← hprod] at hs
          exact Nat.le_of_mul_le_mul_left hs hc.1
        have hnq : n ≤ (n / d) ^ 2 := calc
          n = d * (n / d) := hprod.symm
          _ ≤ (n / d) * (n / d) := Nat.mul_le_mul_right _ hdq
          _ = _ := (pow_two _).symm
        exact (Nat.pow_le_pow_iff_left (by decide : (2 : ℕ) ≠ 0)).mp (hn.trans hnq)
      have hclass' : (n / d) % 12 = 2 := (c12_target_factor_div hn₀ hdn hc).1.trans hclass
      have hsmall := c12_small_divisor_error hA hC hδ hn₁ hdn hc hs (hCount _ hq hclass')
      change |(c12OrderedPairs (n / d) : ℝ) - B / d| ≤
        (B / d) * (c₀ + 6 * Real.log d / L) at hsmall
      calc
        _ ≤ |(μ d : ℝ)| * ((B / d) * (c₀ + 6 * Real.log d / L)) :=
          mul_le_mul_of_nonneg_left hsmall (abs_nonneg _)
        _ = B * (c₀ * (|(μ d : ℝ)| / d) + c₁ * (|(μ d : ℝ)| * Real.log d / d)) := by
          dsimp [c₁]
          ring
        _ ≤ _ := le_add_of_nonneg_right hT
    · have hlarge : Real.sqrt n < d :=
        (Real.sqrt_lt' hd0).mpr (by exact_mod_cast Nat.lt_of_not_ge hs)
      have hq : ((n / d : ℕ) : ℝ) = n / (d : ℝ) := by
        apply (eq_div_iff hd0.ne').mpr
        rw [mul_comm]
        exact_mod_cast Nat.mul_div_cancel' hdn
      have hnonneg : (0 : ℝ) ≤ (c12OrderedPairs (n / d) : ℝ) := Nat.cast_nonneg _
      have hcount : (c12OrderedPairs (n / d) : ℝ) ≤ n / (d : ℝ) := by
        rw [← hq]
        exact_mod_cast c12_ordered_pairs_le (n / d)
      calc
        _ ≤ |(c12OrderedPairs (n / d) : ℝ) - B / d| :=
          mul_le_of_le_one_left (abs_nonneg _) hw
        _ ≤ ((c12OrderedPairs (n / d) : ℝ) + B / d) := by
          simpa only [abs_of_nonneg hnonneg, abs_of_nonneg (div_nonneg hB hd0.le)] using
            (abs_sub (c12OrderedPairs (n / d) : ℝ) (B / d))
        _ ≤ n / (d : ℝ) + B / d := add_le_add_left hcount _
        _ = (B + n) / (d : ℝ) := by ring
        _ ≤ T := div_le_div_of_nonneg_left (by positivity) (Real.sqrt_pos.2 hnR) hlarge.le
        _ ≤ _ := le_add_of_nonneg_left (by positivity)
  have hcop : (c12CoprimePairs n : ℝ) =
      ∑ d ∈ divs, (μ d : ℝ) * (c12OrderedPairs (n / d) : ℝ) := by
    have h := coprime_pair_moebius (fun p => p % 12 = 1) hn₀
    exact_mod_cast h
  change |(c12CoprimePairs n : ℝ) - B * c12CoprimeFactor n| ≤
    B * (c₀ * c12AbsoluteDivisorMass n + c₁ * c12DivisorLogMoment n) +
      (B + n) * (n.divisors.card : ℝ) / Real.sqrt n
  rw [hcop, ← c12_signed_divisor_sum hn₀.ne', Finset.mul_sum, ← Finset.sum_sub_distrib]
  have hrewrite (d : ℕ) :
      (μ d : ℝ) * (c12OrderedPairs (n / d) : ℝ) - B * ((μ d : ℝ) / d) =
        (μ d : ℝ) * ((c12OrderedPairs (n / d) : ℝ) - B / d) := by ring
  simp_rw [hrewrite]
  refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
  calc
    _ ≤ ∑ d ∈ divs, (B * (c₀ * (|(μ d : ℝ)| / d) + c₁ * (|(μ d : ℝ)| * Real.log d / d)) + T) := by
      simp_rw [abs_mul]
      exact Finset.sum_le_sum hterm
    _ = B * (c₀ * c12AbsoluteDivisorMass n + c₁ * c12DivisorLogMoment n) +
        (divs.card : ℝ) * T := by
      unfold c12AbsoluteDivisorMass c12DivisorLogMoment
      change _ = B * (c₀ * (∑ d ∈ divs, |(μ d : ℝ)| / d) +
        c₁ * (∑ d ∈ divs, |(μ d : ℝ)| * Real.log d / d)) + (divs.card : ℝ) * T
      simp only [Finset.sum_add_distrib, Finset.mul_sum, Finset.sum_const, nsmul_eq_mul]
      rw [← Finset.mul_sum, Finset.sum_add_distrib, ← Finset.mul_sum, ← Finset.mul_sum]
    _ ≤ _ := by
      have hcard : (divs.card : ℝ) ≤ n.divisors.card := by
        exact_mod_cast Finset.card_le_card (Finset.filter_subset (s := n.divisors) C12)
      calc
        _ ≤ B * (c₀ * c12AbsoluteDivisorMass n + c₁ * c12DivisorLogMoment n) +
            (n.divisors.card : ℝ) * T := add_le_add_right (mul_le_mul_of_nonneg_right hcard hT) _
        _ = _ := by dsimp [T]; ring

end QuadraticNorms
