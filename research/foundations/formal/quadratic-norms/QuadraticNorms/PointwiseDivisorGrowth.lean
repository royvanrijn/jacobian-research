import QuadraticNorms.C12CoprimeWeights
import Mathlib.NumberTheory.Chebyshev
import Mathlib.NumberTheory.Harmonic.Bounds
import Mathlib.Data.Nat.Log

open scoped BigOperators Topology
open Filter

namespace QuadraticNorms

private theorem dyadic_prime_count_bound :
    ∃ K : ℝ, 0 < K ∧ ∀ j : ℕ,
      ((2 ^ (j + 1) : ℕ).primesLE.card : ℝ) ≤ K * (2 : ℝ) ^ j / ((j : ℝ) + 1) := by
  obtain ⟨Y, hY⟩ := eventually_atTop.mp (Chebyshev.eventually_primeCounting_le (ε := 1) (by norm_num))
  obtain ⟨J, hJ⟩ := exists_nat_ge Y
  let K := max (4 * ((J : ℝ) + 1)) (2 * (Real.log 4 + 1) / Real.log 2)
  have Kpos : 0 < K := lt_of_lt_of_le (by positivity) (le_max_left ..)
  refine ⟨K, Kpos, fun j => ?_⟩
  let P := (2 ^ (j + 1) : ℕ)
  change (P.primesLE.card : ℝ) ≤ _
  have hj : 0 < (j : ℝ) + 1 := by positivity
  have hlog2 : 0 < Real.log (2 : ℝ) := Real.log_pos (by norm_num)
  by_cases h : J ≤ j
  · have hJP : J ≤ P := h.trans ((Nat.lt_two_pow_self (n := j)).le.trans
      (Nat.pow_le_pow_right (by decide : 0 < 2) (Nat.le_succ _)))
    have hYP : Y ≤ (P : ℝ) := hJ.trans (by exact_mod_cast hJP)
    have hcount := hY P hYP
    rw [Nat.floor_natCast, ← Nat.primesLE_card_eq_primeCounting] at hcount
    have hpow : (P : ℝ) = 2 ^ (j + 1) := by simp [P]
    rw [hpow, Real.log_pow, Nat.cast_add, Nat.cast_one, pow_succ (2 : ℝ)] at hcount
    refine hcount.trans (le_trans (le_of_eq ?_) (div_le_div_of_nonneg_right
      (mul_le_mul_of_nonneg_right (le_max_right ..) (by positivity)) hj.le))
    field_simp
  · have hp : P.primesLE.card ≤ P + 1 := by
      rw [← Finset.card_range (P + 1)]
      exact Finset.card_le_card (fun p hp =>
        Finset.mem_range.mpr (Nat.lt_succ_of_le (Nat.mem_primesLE.mp hp).1))
    have hpR : (P.primesLE.card : ℝ) ≤ (P : ℝ) + 1 := by exact_mod_cast hp
    have hpow : (P : ℝ) = (2 : ℝ) ^ j * 2 := by simp [P, pow_succ]
    have hp2 : 1 ≤ (2 : ℝ) ^ j := one_le_pow₀ (by norm_num)
    have hjJ : (j : ℝ) + 1 ≤ (J : ℝ) + 1 :=
      add_le_add_left (by exact_mod_cast (Nat.le_of_lt (Nat.lt_of_not_ge h))) _
    apply (le_div_iff₀ hj).mpr
    calc
      (P.primesLE.card : ℝ) * ((j : ℝ) + 1) ≤ ((P : ℝ) + 1) * ((j : ℝ) + 1) :=
        mul_le_mul_of_nonneg_right hpR hj.le
      _ ≤ ((4 : ℝ) * 2 ^ j) * ((J : ℝ) + 1) := by
        apply mul_le_mul _ hjJ (by positivity) (by positivity)
        rw [hpow]
        linarith only [hp2]
      _ = (4 * ((J : ℝ) + 1)) * (2 : ℝ) ^ j := by ring
      _ ≤ _ := mul_le_mul_of_nonneg_right (le_max_left ..) (by positivity)

private theorem harmonic_as_range (n : ℕ) :
    (harmonic n : ℝ) = ∑ j ∈ Finset.range n, (1 : ℝ) / ((j : ℝ) + 1) := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [harmonic_succ, Finset.sum_range_succ]
    push_cast
    rw [ih, one_div]

/-- A reciprocal-prime bound obtained by summing Chebyshev's upper
bound over dyadic ranges. -/
theorem prime_reciprocal_loglog_bound :
    ∃ K : ℝ, 0 < K ∧ ∀ x : ℕ,
      (∑ p ∈ x.primesLE, (1 : ℝ) / p) ≤ K * (1 + Real.log ((Nat.log 2 x : ℝ) + 1)) := by
  obtain ⟨K, hK, hCount⟩ := dyadic_prime_count_bound
  refine ⟨K, hK, fun x => ?_⟩
  let b := Nat.log 2 x
  let bin (j : ℕ) := x.primesLE.filter (fun p => Nat.log 2 p = j)
  have hmap (p : ℕ) (hp : p ∈ x.primesLE) : Nat.log 2 p ∈ Finset.range (b + 1) :=
    Finset.mem_range.mpr (Nat.lt_succ_of_le (Nat.log_mono_right (Nat.mem_primesLE.mp hp).1))
  have hbin (j : ℕ) : (∑ p ∈ bin j, (1 : ℝ) / p) ≤ K / ((j : ℝ) + 1) := by
    have hp (p : ℕ) (hp : p ∈ bin j) : p.Prime ∧ 2 ^ j ≤ p ∧ p < 2 ^ (j + 1) := by
      obtain ⟨hpx, heq⟩ := Finset.mem_filter.mp hp
      have hp' := (Nat.mem_primesLE.mp hpx).2
      exact ⟨hp', heq ▸ Nat.pow_log_le_self 2 hp'.ne_zero,
        heq ▸ Nat.lt_pow_succ_log_self (by decide) p⟩
    have hcard : ((bin j).card : ℝ) ≤ ((2 ^ (j + 1) : ℕ).primesLE.card : ℝ) := by
      exact_mod_cast Finset.card_le_card (fun p h => Nat.mem_primesLE.mpr ⟨(hp p h).2.2.le, (hp p h).1⟩)
    calc
      _ ≤ ∑ p ∈ bin j, (1 : ℝ) / 2 ^ j := by
        refine Finset.sum_le_sum (fun p h => one_div_le_one_div_of_le (by positivity) ?_)
        exact_mod_cast (hp p h).2.1
      _ = ((bin j).card : ℝ) * (1 / 2 ^ j) := by rw [Finset.sum_const, nsmul_eq_mul]
      _ ≤ _ := by
        have h := mul_le_mul_of_nonneg_right (hcard.trans (hCount j)) (by positivity : (0 : ℝ) ≤ 1 / 2 ^ j)
        refine h.trans_eq ?_
        field_simp
  rw [← Finset.sum_fiberwise_of_maps_to hmap (fun p : ℕ => (1 : ℝ) / p)]
  change (∑ j ∈ Finset.range (b + 1), ∑ p ∈ bin j, (1 : ℝ) / p) ≤ _
  refine (Finset.sum_le_sum (fun j _ => hbin j)).trans ?_
  calc
    _ = K * (harmonic (b + 1) : ℝ) := by
      rw [harmonic_as_range, Finset.mul_sum]
      exact Finset.sum_congr rfl (fun j _ => by ring)
    _ ≤ _ := by
      apply mul_le_mul_of_nonneg_left _ hK.le
      simpa only [Nat.cast_add, Nat.cast_one] using harmonic_le_one_add_log (b + 1)

private theorem prime_divisor_log_sum_le {n : ℕ} (hn : 0 < n)
    {s : Finset ℕ} (hs : s ⊆ n.primeFactors) :
    (∑ p ∈ s, Real.log (p : ℝ)) ≤ Real.log n := by
  have hd : (∏ p ∈ s, p) ∣ n :=
    (Finset.prod_dvd_prod_of_subset s n.primeFactors id hs).trans (Nat.prod_primeFactors_dvd n)
  have hp (p : ℕ) (hp : p ∈ s) : (0 : ℝ) < p := by exact_mod_cast (Nat.prime_of_mem_primeFactors (hs hp)).pos
  have hprod : (∏ p ∈ s, (p : ℝ)) ≤ n := by
    rw [← Nat.cast_prod]
    exact_mod_cast Nat.le_of_dvd hn hd
  rw [← Real.log_prod (fun p hp' => (hp p hp').ne')]
  exact Real.log_le_log (Finset.prod_pos hp) hprod

/-- The reciprocal mass of the prime divisors is O(log log log n),
with an explicit logarithmic expression on the right. -/
theorem prime_divisor_reciprocal_growth :
    ∃ K : ℝ, 0 < K ∧ ∀ n : ℕ, 2 ≤ Real.log n → ∀ s : Finset ℕ, s ⊆ n.primeFactors →
      (∑ p ∈ s, (1 : ℝ) / p) ≤ K * (1 + Real.log (1 + Real.log (Real.log n))) := by
  obtain ⟨K, hK, hrecip⟩ := prime_reciprocal_loglog_bound
  let a := (Real.log 2)⁻¹
  let C := K * (1 + Real.log (1 + a)) + a
  have ha : 0 < a := inv_pos.mpr (Real.log_pos (by norm_num))
  have hC : 0 ≤ C := by
    have halog := Real.log_nonneg (show (1 : ℝ) ≤ 1 + a by linarith)
    dsimp [C]
    positivity
  refine ⟨C + K, by linarith, fun n hn s hs => ?_⟩
  let L := Real.log n
  let y := ⌊L⌋₊
  let ℓ := 1 + Real.log L
  let small := s.filter (fun p => p ≤ y)
  let large := s.filter (fun p => ¬ p ≤ y)
  have hL : 0 < L := by linarith
  have hn₀ : 0 < n := by
    by_contra! h
    obtain rfl : n = 0 := by omega
    norm_num at hn
  have hlogL : 0 ≤ Real.log L := Real.log_nonneg (by linarith)
  have hell : 1 ≤ ℓ := by dsimp [ℓ]; linarith
  have hy : 2 ≤ y := Nat.le_floor (by exact_mod_cast hn)
  have hyL : (y : ℝ) ≤ L := Nat.floor_le hL.le
  have hb : (Nat.log 2 y : ℝ) ≤ a * Real.log L := by
    have hp := Nat.pow_log_le_self 2 (show y ≠ 0 by omega)
    have hpR : (2 : ℝ) ^ Nat.log 2 y ≤ L := (by exact_mod_cast hp : (2 : ℝ) ^ Nat.log 2 y ≤ y).trans hyL
    have hl := Real.log_le_log (by positivity : (0 : ℝ) < 2 ^ Nat.log 2 y) hpR
    rw [Real.log_pow] at hl
    exact (le_div_iff₀ (Real.log_pos (by norm_num : (1 : ℝ) < 2))).mpr hl |>.trans_eq
      (by simp [a, div_eq_mul_inv, mul_comm])
  have hhead : (∑ p ∈ small, (1 : ℝ) / p) ≤ K * (1 + Real.log (1 + a) + Real.log ℓ) := by
    have hsub : small ⊆ y.primesLE := by
      intro p hp
      obtain ⟨hps, hpL⟩ := Finset.mem_filter.mp hp
      exact Nat.mem_primesLE.mpr ⟨hpL, Nat.prime_of_mem_primeFactors (hs hps)⟩
    refine (Finset.sum_le_sum_of_subset_of_nonneg hsub (fun _ _ _ => by positivity)).trans
      ((hrecip y).trans (mul_le_mul_of_nonneg_left ?_ hK.le))
    have hz : (Nat.log 2 y : ℝ) + 1 ≤ (1 + a) * ℓ := by
      dsimp [ℓ]
      nlinarith only [hb, ha.le, hlogL]
    have hl := Real.log_le_log (by positivity) hz
    rw [Real.log_mul (by positivity) (by linarith)] at hl
    linarith
  have htail : (∑ p ∈ large, (1 : ℝ) / p) ≤ a := by
    have hpoint (p : ℕ) (hp : p ∈ large) : (1 : ℝ) / p ≤ Real.log p / (L * Real.log 2) := by
      obtain ⟨hps, hp'⟩ := Finset.mem_filter.mp hp
      have hpP := Nat.prime_of_mem_primeFactors (hs hps)
      have hp2 : (2 : ℝ) ≤ p := by exact_mod_cast hpP.two_le
      have hpL : L ≤ p := by
        have h := Nat.lt_floor_add_one L
        have hle : (y : ℝ) + 1 ≤ p := by
          exact_mod_cast Nat.succ_le_of_lt (Nat.lt_of_not_ge hp')
        exact h.le.trans hle
      have hlogp := Real.log_le_log (show (0 : ℝ) < 2 by norm_num) hp2
      apply (div_le_div_iff₀ (by linarith) (mul_pos hL (Real.log_pos (by norm_num)))).mpr
      simpa only [one_mul, mul_comm] using
        mul_le_mul hpL hlogp (Real.log_pos (by norm_num : (1 : ℝ) < 2)).le (by linarith)
    refine (Finset.sum_le_sum hpoint).trans ?_
    rw [← Finset.sum_div]
    refine (div_le_div_of_nonneg_right (prime_divisor_log_sum_le hn₀
      ((Finset.filter_subset _ _).trans hs)) (by positivity)).trans_eq ?_
    change L / (L * Real.log 2) = a
    unfold a
    field_simp
  rw [← Finset.sum_filter_add_sum_filter_not s (fun p => p ≤ y)]
  change (∑ p ∈ small, (1 : ℝ) / p) + (∑ p ∈ large, (1 : ℝ) / p) ≤ (C + K) * (1 + Real.log ℓ)
  have hlogell := Real.log_nonneg hell
  have hsum : (∑ p ∈ small, (1 : ℝ) / p) + (∑ p ∈ large, (1 : ℝ) / p) ≤
      C + K * Real.log ℓ := by
    dsimp [C]
    nlinarith only [hhead, htail]
  exact hsum.trans (by nlinarith only [hK, mul_nonneg hC hlogell])

private theorem prime_divisor_logarithmic_weight_bound :
    ∃ K : ℝ, 0 < K ∧ ∀ n : ℕ, 2 ≤ Real.log n → ∀ s : Finset ℕ, s ⊆ n.primeFactors →
      (∑ p ∈ s, Real.log (p : ℝ) / ((p : ℝ) + 1)) ≤ K * (1 + Real.log (Real.log n)) ^ 2 := by
  obtain ⟨K, hK, hrecip⟩ := prime_divisor_reciprocal_growth
  refine ⟨K + 1, by linarith, fun n hn s hs => ?_⟩
  let L := Real.log n
  let y := ⌊L⌋₊
  let ℓ := 1 + Real.log L
  let small := s.filter (fun p => p ≤ y)
  let large := s.filter (fun p => ¬ p ≤ y)
  have hL : 0 < L := by linarith
  have hn₀ : 0 < n := Nat.pos_of_ne_zero (by rintro rfl; norm_num at hn)
  have hlog₀ : 0 ≤ Real.log L := Real.log_nonneg (by linarith)
  have hlogℓ : Real.log ℓ ≤ ℓ - 1 := Real.log_le_sub_one_of_pos (by dsimp [ℓ]; positivity)
  have hℓ : 1 ≤ ℓ := by dsimp [ℓ]; linarith
  have hsmall : (∑ p ∈ small, Real.log (p : ℝ) / ((p : ℝ) + 1)) ≤ K * ℓ ^ 2 := by
    have hp (p : ℕ) (hp : p ∈ small) : (0 : ℝ) < p ∧ (p : ℝ) ≤ L := by
      obtain ⟨hps, hpy⟩ := Finset.mem_filter.mp hp
      exact ⟨by exact_mod_cast (Nat.prime_of_mem_primeFactors (hs hps)).pos,
        (show (p : ℝ) ≤ y from by exact_mod_cast hpy).trans (Nat.floor_le hL.le)⟩
    calc
      _ ≤ ∑ p ∈ small, Real.log L * (1 / (p : ℝ)) := by
        apply Finset.sum_le_sum
        intro p hp'
        obtain ⟨hp₀, hpL⟩ := hp p hp'
        rw [mul_one_div]
        exact (div_le_div_of_nonneg_right (Real.log_le_log hp₀ hpL) (by positivity)).trans
          (div_le_div_of_nonneg_left hlog₀ hp₀ (by linarith))
      _ = Real.log L * (∑ p ∈ small, (1 : ℝ) / p) := (Finset.mul_sum ..).symm
      _ ≤ ℓ * (K * (1 + Real.log ℓ)) := by
        apply mul_le_mul (by dsimp [ℓ]; linarith)
          (hrecip n hn small ((Finset.filter_subset _ _).trans hs))
          (Finset.sum_nonneg (fun _ _ => by positivity)) (by linarith)
      _ ≤ _ := by
        calc
          ℓ * (K * (1 + Real.log ℓ)) ≤ ℓ * (K * ℓ) :=
            mul_le_mul_of_nonneg_left (mul_le_mul_of_nonneg_left (by linarith only [hlogℓ]) hK.le)
              (by linarith)
          _ = _ := by ring
  have hlarge : (∑ p ∈ large, Real.log (p : ℝ) / ((p : ℝ) + 1)) ≤ 1 := by
    have hpoint (p : ℕ) (hp : p ∈ large) : Real.log (p : ℝ) / ((p : ℝ) + 1) ≤ Real.log p / L := by
      obtain ⟨hps, hpy⟩ := Finset.mem_filter.mp hp
      have hpL : L ≤ (p : ℝ) := (Nat.lt_floor_add_one L).le.trans
        (by exact_mod_cast Nat.succ_le_of_lt (Nat.lt_of_not_ge hpy))
      have hp₂ : (2 : ℝ) ≤ p := by exact_mod_cast (Nat.prime_of_mem_primeFactors (hs hps)).two_le
      exact div_le_div_of_nonneg_left (Real.log_nonneg (by linarith)) hL (by linarith)
    refine (Finset.sum_le_sum hpoint).trans ?_
    rw [← Finset.sum_div]
    exact (div_le_one hL).mpr (prime_divisor_log_sum_le hn₀ ((Finset.filter_subset _ _).trans hs))
  rw [← Finset.sum_filter_add_sum_filter_not s (fun p => p ≤ y)]
  change (∑ p ∈ small, Real.log (p : ℝ) / ((p : ℝ) + 1)) +
    (∑ p ∈ large, Real.log (p : ℝ) / ((p : ℝ) + 1)) ≤ (K + 1) * ℓ ^ 2
  exact (add_le_add hsmall hlarge).trans (by nlinarith only [hℓ])

private theorem divisor_reciprocal_products_bound :
    ∃ K : ℝ, 0 < K ∧ ∀ n : ℕ, 2 ≤ Real.log n →
      c12AbsoluteDivisorMass n ≤ Real.exp (K * (1 + Real.log (1 + Real.log (Real.log n)))) ∧
      (c12CoprimeFactor n)⁻¹ ≤ Real.exp (K * (1 + Real.log (1 + Real.log (Real.log n)))) := by
  obtain ⟨K, hK, hrec⟩ := prime_divisor_reciprocal_growth
  refine ⟨2 * K, by linarith, fun n hn => ?_⟩
  let s := c12AllowedPrimeFactors n
  let L := Real.log n
  let ℓ := 1 + Real.log L
  have hn₀ : n ≠ 0 := by rintro rfl; norm_num at hn
  have hsub : s ⊆ n.primeFactors := Finset.filter_subset ..
  have hrec' : (∑ p ∈ s, (1 : ℝ) / p) ≤ K * (1 + Real.log ℓ) := hrec n hn s hsub
  have hℓ : 1 ≤ ℓ := by dsimp [ℓ]; linarith [Real.log_nonneg (show 1 ≤ L by linarith)]
  have hlogℓ := Real.log_nonneg hℓ
  have h2 (p : ℕ) (hp : p ∈ s) : (2 : ℝ) ≤ p := by
    exact_mod_cast (Nat.prime_of_mem_primeFactors (hsub hp)).two_le
  constructor
  · rw [c12_absolute_divisor_mass_product hn₀]
    change (∏ p ∈ s, (1 + (p : ℝ)⁻¹)) ≤ Real.exp (2 * K * (1 + Real.log ℓ))
    refine (Real.prod_one_add_le_exp_sum _ (fun _ => by positivity)).trans (Real.exp_le_exp.mpr ?_)
    simp_rw [← one_div]
    exact hrec'.trans (by nlinarith only [hK, hlogℓ])
  · let f (p : ℕ) : ℝ := 1 / max 1 ((p : ℝ) - 1)
    unfold c12CoprimeFactor
    rw [← Finset.prod_inv_distrib]
    have hterm (p : ℕ) (hp : p ∈ s) : (1 - (p : ℝ)⁻¹)⁻¹ = 1 + f p := by
      have hpR := h2 p hp
      have hp₀ : (p : ℝ) ≠ 0 := by linarith
      have hp₁ : (p : ℝ) - 1 ≠ 0 := by linarith
      have hd : 1 - (p : ℝ)⁻¹ ≠ 0 :=
        (sub_pos.mpr (inv_lt_one_of_one_lt₀ (by linarith))).ne'
      unfold f
      rw [max_eq_right (by linarith)]
      field_simp
      ring
    rw [Finset.prod_congr rfl hterm]
    refine (Real.prod_one_add_le_exp_sum _ (fun p => one_div_nonneg.mpr ((by norm_num : (0 : ℝ) ≤ 1).trans
      (le_max_left ..)))).trans (Real.exp_le_exp.mpr ?_)
    have hpoint (p : ℕ) (hp : p ∈ s) : f p ≤ 2 * (1 / (p : ℝ)) := by
      have hpR := h2 p hp
      unfold f
      rw [max_eq_right (by linarith), mul_one_div]
      exact (div_le_div_iff₀ (by linarith) (by linarith)).mpr (by linarith)
    exact (Finset.sum_le_sum hpoint).trans (by rw [← Finset.mul_sum]; linarith only [hrec'])

/-- For every ε>0, all three divisor factors grow at most as
(log n)^ε, uniformly in the target. This weaker form of the familiar
log-log estimates suffices for the exponent in the coprime corollary. -/
theorem c12_divisor_growth {ε : ℝ} (hε : 0 < ε) :
    ∃ K : ℝ, 0 < K ∧ ∀ n : ℕ, 2 ≤ Real.log n →
      c12AbsoluteDivisorMass n ≤ K * (Real.log n) ^ ε ∧
      c12DivisorLogMoment n ≤ K * (Real.log n) ^ ε ∧
      (c12CoprimeFactor n)⁻¹ ≤ K * (Real.log n) ^ ε := by
  obtain ⟨B, hB, hprod⟩ := divisor_reciprocal_products_bound
  obtain ⟨D, hD, hmoment⟩ := prime_divisor_logarithmic_weight_bound
  let P := B + 2
  let η := ε / P
  let K := (D + 1) * Real.exp B * (1 + 1 / η) ^ P
  have hP : 0 < P := by dsimp [P]; linarith
  have hη : 0 < η := div_pos hε hP
  refine ⟨K, by dsimp [K]; positivity, fun n hn => ?_⟩
  let L := Real.log n
  let ℓ := 1 + Real.log L
  let W := Real.exp B * ℓ ^ B
  have hL : 1 ≤ L := by linarith
  have hL₀ : 0 ≤ L := by linarith
  have hℓ : 1 ≤ ℓ := by
    dsimp [ℓ]
    linarith [Real.log_nonneg hL]
  have hℓ₀ : 0 < ℓ := zero_lt_one.trans_le hℓ
  have hlog : ℓ ≤ (1 + 1 / η) * L ^ η := by
    have h := Real.log_le_rpow_div hL₀ hη
    have h' : 1 ≤ L ^ η := by
      simpa only [Real.one_rpow] using
        Real.rpow_le_rpow (by norm_num : (0 : ℝ) ≤ 1) hL hη.le
    calc
      ℓ ≤ L ^ η + L ^ η / η := add_le_add h' h
      _ = _ := by ring
  have hpow : ℓ ^ P ≤ (1 + 1 / η) ^ P * L ^ ε := by
    refine (Real.rpow_le_rpow hℓ₀.le hlog hP.le).trans_eq ?_
    rw [Real.mul_rpow (by positivity) (by positivity), ← Real.rpow_mul hL₀]
    congr 1
    rw [show η * P = ε from div_mul_cancel₀ _ hP.ne']
  have hW : Real.exp (B * (1 + Real.log ℓ)) = W := by
    unfold W
    rw [mul_add, mul_one, Real.exp_add, Real.rpow_def_of_pos hℓ₀, mul_comm B]
  obtain ⟨hs₀, hpInv⟩ := hprod n hn
  change c12AbsoluteDivisorMass n ≤ Real.exp (B * (1 + Real.log ℓ)) at hs₀
  change (c12CoprimeFactor n)⁻¹ ≤ Real.exp (B * (1 + Real.log ℓ)) at hpInv
  rw [hW] at hs₀ hpInv
  have hcommon : W ≤ K * L ^ ε := by
    have hb : Real.exp B ≤ (D + 1) * Real.exp B :=
      le_mul_of_one_le_left (Real.exp_nonneg _) (by linarith)
    calc
      W ≤ ((D + 1) * Real.exp B) * ℓ ^ P :=
        mul_le_mul hb (Real.rpow_le_rpow_of_exponent_le hℓ (by dsimp [P]; linarith))
          (by positivity) (by positivity)
      _ ≤ _ := by
        dsimp [K]
        rw [mul_assoc ((D + 1) * Real.exp B)]
        exact mul_le_mul_of_nonneg_left hpow (by positivity)
  refine ⟨hs₀.trans hcommon, ?_, hpInv.trans hcommon⟩
  have hs : c12AllowedPrimeFactors n ⊆ n.primeFactors := Finset.filter_subset ..
  calc
    c12DivisorLogMoment n ≤ W * (D * ℓ ^ 2) := by
      rw [c12_divisor_log_moment_product (by rintro rfl; norm_num at hn)]
      apply mul_le_mul hs₀ (hmoment n hn _ hs) _ (by dsimp [W]; positivity)
      apply Finset.sum_nonneg
      intro p hp
      exact div_nonneg (Real.log_nonneg (by exact_mod_cast (Nat.prime_of_mem_primeFactors (hs hp)).one_lt.le))
        (by positivity)
    _ = (D * Real.exp B) * ℓ ^ P := by
      dsimp [W, P]
      rw [Real.rpow_add hℓ₀, show ℓ ^ (2 : ℝ) = ℓ ^ (2 : ℕ) from by
        exact Real.rpow_natCast ℓ 2]
      ring
    _ ≤ ((D + 1) * Real.exp B) * ℓ ^ P := by
      exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_right (by linarith) (Real.exp_nonneg _))
        (by positivity)
    _ ≤ _ := by
      dsimp [K]
      rw [mul_assoc ((D + 1) * Real.exp B)]
      exact mul_le_mul_of_nonneg_left hpow (by positivity)

private theorem fourth_power_factor_step {k : ℕ} (hk : 5 ≤ k) : (k + 2) ^ 4 ≤ 2 * (k + 1) ^ 4 := by
  have hratio : 6 * (k + 2) ≤ 7 * (k + 1) := by omega
  apply Nat.le_of_mul_le_mul_left (c := 6 ^ 4)
  · calc
      _ = (6 * (k + 2)) ^ 4 := (mul_pow ..).symm
      _ ≤ (7 * (k + 1)) ^ 4 := Nat.pow_le_pow_left hratio _
      _ = 7 ^ 4 * (k + 1) ^ 4 := mul_pow ..
      _ ≤ (6 ^ 4 * 2) * (k + 1) ^ 4 := Nat.mul_le_mul_right _ (by decide)
      _ = _ := by ring
  · decide

private theorem fourth_power_divisor_factor {p : ℕ} (hp : 2 ≤ p) (k : ℕ) :
    (k + 1) ^ 4 ≤ (if p < 16 then 7 ^ 4 else 1) * p ^ k := by
  by_cases hp16 : p < 16
  · rw [if_pos hp16]
    induction k with
    | zero => norm_num
    | succ k ih =>
      by_cases hk : 5 ≤ k
      · calc
          _ ≤ 2 * (k + 1) ^ 4 := fourth_power_factor_step hk
          _ ≤ p * (7 ^ 4 * p ^ k) := Nat.mul_le_mul hp ih
          _ = _ := by rw [pow_succ]; ring
      · calc
          _ ≤ (7 : ℕ) ^ 4 := Nat.pow_le_pow_left (by omega) _
          _ ≤ _ := Nat.le_mul_of_pos_right _ (pow_pos (by omega) _)
  · rw [if_neg hp16, one_mul]
    calc
      _ ≤ (2 ^ k : ℕ) ^ 4 := Nat.pow_le_pow_left (Nat.lt_two_pow_self (n := k)) _
      _ = ((2 : ℕ) ^ 4) ^ k := by rw [← pow_mul, mul_comm k, pow_mul]
      _ ≤ _ := Nat.pow_le_pow_left (Nat.le_of_not_gt hp16) _

/-- An elementary fourth-power divisor bound, with a fixed constant.
It uses only the small primes p<16 to absorb the possible factor losses. -/
theorem divisor_count_fourth_power_bound {n : ℕ} (hn : n ≠ 0) :
    n.divisors.card ^ 4 ≤ (7 ^ 4) ^ 16 * n := by
  let s := n.primeFactors
  let small := s.filter (fun p => p < 16)
  let K (p : ℕ) : ℕ := if p < 16 then 7 ^ 4 else 1
  have hprod := Nat.prod_factorization_pow_eq_self hn
  rw [Finsupp.prod, Nat.support_factorization] at hprod
  rw [Nat.card_divisors hn]
  change (∏ p ∈ s, (n.factorization p + 1)) ^ 4 ≤ _
  rw [← Finset.prod_pow]
  have hf (p : ℕ) (hp : p ∈ s) :
      (n.factorization p + 1) ^ 4 ≤ K p * p ^ n.factorization p :=
    fourth_power_divisor_factor (Nat.prime_of_mem_primeFactors hp).two_le _
  have hK : (∏ p ∈ s, K p) = (7 ^ 4) ^ small.card := by
    rw [← Finset.prod_filter_of_ne (p := fun p => p < 16)
      (fun p _ hp => by by_contra! h; exact hp (if_neg (Nat.not_lt_of_ge h)))]
    have he : (∏ p ∈ small, K p) = ∏ _p ∈ small, (7 ^ 4 : ℕ) :=
      Finset.prod_congr rfl (fun p hp => if_pos (Finset.mem_filter.mp hp).2)
    rw [he, Finset.prod_const]
  have hcard : small.card ≤ 16 := by
    rw [← Finset.card_range 16]
    exact Finset.card_le_card (fun p hp => Finset.mem_range.mpr (Finset.mem_filter.mp hp).2)
  calc
    _ ≤ ∏ p ∈ s, K p * p ^ n.factorization p := Finset.prod_le_prod (fun _ _ => zero_le) hf
    _ = (7 ^ 4) ^ small.card * n := by rw [Finset.prod_mul_distrib, hK, hprod]
    _ ≤ _ := Nat.mul_le_mul_right _ (Nat.pow_le_pow_right (by decide : 0 < 7 ^ 4) hcard)

/-- The bound used for the large-divisor remainder of the coprime
transfer. -/
theorem divisor_count_log_cube_over_sqrt {n : ℕ} (hn : 1 < n) :
    (n.divisors.card : ℝ) / Real.sqrt n ≤ (7 : ℝ) ^ 16 * 16 ^ 3 / (Real.log n) ^ 3 := by
  have hn0 : (0 : ℝ) < n := by exact_mod_cast (Nat.zero_lt_of_lt hn)
  have hn1 : (1 : ℝ) ≤ n := by exact_mod_cast hn.le
  have hlog : 0 < Real.log n := Real.log_pos (by exact_mod_cast hn)
  have hsqrt : 0 < Real.sqrt (n : ℝ) := Real.sqrt_pos.2 hn0
  have hτ : (n.divisors.card : ℝ) ^ 4 ≤ ((7 : ℝ) ^ 4) ^ 16 * n := by
    exact_mod_cast divisor_count_fourth_power_bound (Nat.zero_lt_of_lt hn).ne'
  have hlog_bound : Real.log n ≤ 16 * (n : ℝ) ^ (1 / 16 : ℝ) := by
    have hh := Real.log_natCast_le_rpow_div n (by norm_num : (0 : ℝ) < 1 / 16)
    convert hh using 1
    ring
  have hlog_pow : Real.log n ^ 12 ≤ (16 : ℝ) ^ 12 * n := by
    calc
      _ ≤ (16 * (n : ℝ) ^ (1 / 16 : ℝ)) ^ (12 : ℕ) :=
        pow_le_pow_left₀ hlog.le hlog_bound _
      _ = (16 : ℝ) ^ 12 * (n : ℝ) ^ (3 / 4 : ℝ) := by
        rw [mul_pow, ← Real.rpow_mul_natCast hn0.le]
        norm_num
      _ ≤ _ := by
        apply mul_le_mul_of_nonneg_left _ (by positivity)
        simpa only [Real.rpow_one] using Real.rpow_le_rpow_of_exponent_le hn1
          (by norm_num : (3 / 4 : ℝ) ≤ 1)
  have hp : (n.divisors.card : ℝ) * Real.log n ^ 3 ≤ (7 : ℝ) ^ 16 * 16 ^ 3 * Real.sqrt n := by
    apply (pow_le_pow_iff_left₀ (by positivity) (by positivity) (by decide : (4 : ℕ) ≠ 0)).mp
    rw [mul_pow, ← pow_mul]
    refine (mul_le_mul hτ hlog_pow (by positivity) (by positivity)).trans_eq ?_
    have hr : Real.sqrt (n : ℝ) ^ 4 = (n : ℝ) ^ 2 := by
      conv_lhs => rw [show (4 : ℕ) = 2 * 2 by decide, pow_mul, Real.sq_sqrt hn0.le]
    rw [mul_pow, mul_pow, hr]
    ring
  apply (div_le_div_iff₀ hsqrt (pow_pos hlog _)).mpr
  simpa only [mul_comm, mul_left_comm, mul_assoc] using hp

end QuadraticNorms
