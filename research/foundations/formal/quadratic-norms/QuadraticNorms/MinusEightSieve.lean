import QuadraticNorms.MinusEightBadCounts
import QuadraticNorms.MinusEightPrimeBounds
import QuadraticNorms.MinusEightThetaNormalization
import QuadraticNorms.SieveCriteria

open scoped BigOperators ArithmeticFunction.sigma

namespace QuadraticNorms

/-- The selected-shell inequality with the prime count and the inverse-square tail
bounded at the chosen cutoff. The three strata are m, 2m, and 4m. -/
def MinusEightSelectedMarginWithBounds (e m : ℕ) (count tail : ℝ) : Prop :=
  let lower : ℝ := if e = 2 then 2 / 3 else 1 / 2
  let δ : ℝ := if e = 0 then 1 / 2 else 0
  let L : ℝ := if e = 0 then 8 / 3 else if e = 1 then 2 else 4 / 3
  3 / 16 + minusEightShellError m * (count + δ) +
    L * minusEightDivisorMaximum (2 ^ e * m) * tail < lower

/-- Bounds for the two finite prime sums at a fixed sieve cutoff. -/
structure MinusEightCutoffBounds (X count tail : ℝ) : Prop where
  count_le : ∀ s : Finset ℕ, (∀ p ∈ s, MinusEightInertPrime p ∧ (p : ℝ) ≤ X) →
    (s.card : ℝ) ≤ count
  tail_le : ∀ s : Finset ℕ, (∀ p ∈ s, MinusEightInertPrime p ∧ X < (p : ℝ)) →
    (∑ p ∈ s, inverseSquare p) ≤ tail

theorem minusEight_cutoff_bounds_progressions {X : ℝ} (hX : 0 < X) :
    MinusEightCutoffBounds X (X / 4 + 1) (1 / (4 * X) + 2 / X ^ 2) := by
  constructor
  · intro s hs
    refine (minusEight_inert_card_le (P := ⌊X⌋₊) fun p h =>
      ⟨(hs p h).1, Nat.le_floor (hs p h).2⟩).trans ?_
    linarith [Nat.floor_le hX.le]
  · intro s hs
    apply minusEight_inert_tail_le (P := ⌊X⌋₊) hX (Nat.lt_floor_add_one _).le
    exact fun p h => ⟨(hs p h).1, (Nat.floor_lt hX.le).mpr (hs p h).2⟩

/-- The selected coordinate vectors surviving every inert prime dividing a positive block. -/
def minusEightSieveSurvivors (n : ℕ) (t : Finset MinusEightShellPoint) :
    Finset MinusEightShellPoint :=
  t \ (((Finset.range (n + 1)).filter MinusEightInertPrime).biUnion
    (fun p => t ∩ minusEightEitherDivisible n p))

theorem mem_minusEight_sieve_survivors {n : ℕ} {t : Finset MinusEightShellPoint}
    {v : MinusEightShellPoint} :
    v ∈ minusEightSieveSurvivors n t ↔ v ∈ t ∧
      ∀ p ≤ n, MinusEightInertPrime p → v ∉ minusEightEitherDivisible n p := by
  rw [minusEightSieveSurvivors, Finset.mem_sdiff]
  constructor
  · rintro ⟨ht, hbad⟩
    exact ⟨ht, fun p hpn hi hp =>
      hbad (Finset.mem_biUnion.mpr ⟨p,
        Finset.mem_filter.mpr ⟨Finset.mem_range.mpr (Nat.lt_succ_of_le hpn), hi⟩,
        Finset.mem_inter.mpr ⟨ht, hp⟩⟩)⟩
  · rintro ⟨ht, hgood⟩
    exact ⟨ht, fun hbad => by
      obtain ⟨p, hp, hbad⟩ := Finset.mem_biUnion.mp hbad
      obtain ⟨hpR, hi⟩ := Finset.mem_filter.mp hp
      exact hgood p (Nat.lt_succ_iff.mp (Finset.mem_range.mp hpR)) hi (Finset.mem_inter.mp hbad).2⟩

/-- The selected-shell sieve uses full-shell small-prime bounds and
restricted large-prime counts. Its survivors avoid every inert prime,
including the primes dividing the target. -/
theorem minusEight_sieve_survivor_card_le {n : ℕ} {X count tail lower δ E L : ℝ}
    (t : Finset MinusEightShellPoint)
    (hcut : MinusEightCutoffBounds X count tail) (hE : 0 ≤ E) (hL : 0 ≤ L)
    (hselected : (lower - δ * E) * (minusEightFullShell n).card ≤ t.card)
    (hsmall : ∀ p ≤ n, MinusEightInertPrime p → (p : ℝ) ≤ X →
      ((minusEightEitherDivisible n p).card : ℝ) ≤
        (minusEightLocalMass p + E) * (minusEightFullShell n).card)
    (hlarge : ∀ p ≤ n, MinusEightInertPrime p → X < (p : ℝ) →
      ((t ∩ minusEightEitherDivisible n p).card : ℝ) ≤
        L * inverseSquare p * (minusEightFullShell n).card) :
    (lower - (3 / 16 + E * (count + δ) + L * tail)) * (minusEightFullShell n).card ≤
      (minusEightSieveSurvivors n t).card := by
  let primes := (Finset.range (n + 1)).filter MinusEightInertPrime
  let small := primes.filter (fun (p : ℕ) => (p : ℝ) ≤ X)
  let large := primes.filter (fun (p : ℕ) => ¬ (p : ℝ) ≤ X)
  let bad (p : ℕ) := t ∩ minusEightEitherDivisible n p
  let R : ℝ := (minusEightFullShell n).card
  have hR : 0 ≤ R := by positivity
  have hp (p : ℕ) (hp : p ∈ primes) : p ≤ n ∧ MinusEightInertPrime p := by
    obtain ⟨h₁, h₂⟩ := Finset.mem_filter.mp hp
    exact ⟨Nat.lt_succ_iff.mp (Finset.mem_range.mp h₁), h₂⟩
  have hs (p : ℕ) (h : p ∈ small) : p ≤ n ∧ MinusEightInertPrime p ∧ (p : ℝ) ≤ X := by
    obtain ⟨h₁, h₂⟩ := Finset.mem_filter.mp h
    exact ⟨(hp p h₁).1, (hp p h₁).2, h₂⟩
  have hl (p : ℕ) (h : p ∈ large) : p ≤ n ∧ MinusEightInertPrime p ∧ X < (p : ℝ) := by
    obtain ⟨h₁, h₂⟩ := Finset.mem_filter.mp h
    exact ⟨(hp p h₁).1, (hp p h₁).2, lt_of_not_ge h₂⟩
  have hcount := hcut.count_le small (fun p hp => (hs p hp).2)
  have hmass := minusEight_inert_mass_lt (s := small) (fun p h => (hs p h).2.1)
  have htls := hcut.tail_le large (fun p hp => (hl p hp).2)
  have hsb : ∑ p ∈ small, ((bad p).card : ℝ) ≤ (3 / 16 + E * (count)) * R := by
    calc
      _ ≤ ∑ p ∈ small, (minusEightLocalMass p + E) * R := by
        apply Finset.sum_le_sum
        intro p h
        have hsub : (bad p).card ≤ (minusEightEitherDivisible n p).card :=
          Finset.card_le_card Finset.inter_subset_right
        exact (by exact_mod_cast hsub : ((bad p).card : ℝ) ≤ (minusEightEitherDivisible n p).card).trans
          (hsmall p (hs p h).1 (hs p h).2.1 (hs p h).2.2)
      _ = ((∑ p ∈ small, minusEightLocalMass p) + E * small.card) * R := by
        rw [← Finset.sum_mul, Finset.sum_add_distrib, Finset.sum_const, nsmul_eq_mul, mul_comm _ E]
      _ ≤ _ := mul_le_mul_of_nonneg_right
        (add_le_add hmass.le (mul_le_mul_of_nonneg_left hcount hE)) hR
  have hlb : ∑ p ∈ large, ((bad p).card : ℝ) ≤ (L * tail) * R := by
    calc
      _ ≤ ∑ p ∈ large, L * inverseSquare p * R :=
        Finset.sum_le_sum fun p h => hlarge p (hl p h).1 (hl p h).2.1 (hl p h).2.2
      _ = (L * ∑ p ∈ large, inverseSquare p) * R := by rw [← Finset.sum_mul, Finset.mul_sum]
      _ ≤ _ := mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left htls hL) hR
  let U := primes.biUnion bad
  have htotal : ∑ p ∈ primes, ((bad p).card : ℝ) ≤
      (3 / 16 + E * count + L * tail) * R := by
    rw [← Finset.sum_filter_add_sum_filter_not primes (fun (p : ℕ) => (p : ℝ) ≤ X)]
    exact (add_le_add hsb hlb).trans_eq (by ring)
  have hcard : (t.card : ℝ) ≤ (minusEightSieveSurvivors n t).card + ∑ p ∈ primes, ((bad p).card : ℝ) := by
    change (t.card : ℝ) ≤ (t \ U).card + _
    have hU : U.card ≤ ∑ p ∈ primes, (bad p).card := Finset.card_biUnion_le
    have hc := Finset.card_sdiff_add_card_inter t U
    have hi := Finset.card_le_card (s := t ∩ U) Finset.inter_subset_right
    have hnat : t.card ≤ (t \ U).card + ∑ p ∈ primes, (bad p).card := by omega
    exact_mod_cast hnat
  linarith only [htotal, hcard, hselected]

theorem minusEight_exists_unspoiled_of_cutoff_bounds {n : ℕ} {X count tail lower δ E L : ℝ}
    (t : Finset MinusEightShellPoint) (ht : t ⊆ minusEightDyadicShell n)
    (hn : 0 < (minusEightFullShell n).card) (hcut : MinusEightCutoffBounds X count tail) (hE : 0 ≤ E) (hL : 0 ≤ L)
    (hselected : (lower - δ * E) * (minusEightFullShell n).card ≤ t.card)
    (hsmall : ∀ p ≤ n, MinusEightInertPrime p → (p : ℝ) ≤ X →
      ((minusEightEitherDivisible n p).card : ℝ) ≤
        (minusEightLocalMass p + E) * (minusEightFullShell n).card)
    (hlarge : ∀ p ≤ n, MinusEightInertPrime p → X < (p : ℝ) →
      ((t ∩ minusEightEitherDivisible n p).card : ℝ) ≤
        L * inverseSquare p * (minusEightFullShell n).card)
    (hmargin : 3 / 16 + E * (count + δ) + L * tail < lower) :
    ∃ v ∈ t, v ∈ minusEightDyadicShell n ∧
      ∀ p ≤ n, MinusEightInertPrime p → v ∉ minusEightEitherDivisible n p := by
  have h := minusEight_sieve_survivor_card_le t hcut hE hL hselected hsmall hlarge
  have hpos : (0 : ℝ) < (minusEightSieveSurvivors n t).card :=
    (mul_pos (sub_pos.mpr hmargin) (Nat.cast_pos.mpr hn)).trans_le h
  obtain ⟨v, hv⟩ := Finset.card_pos.mp (Nat.cast_pos.mp hpos)
  obtain ⟨hvt, hvbad⟩ := mem_minusEight_sieve_survivors.mp hv
  exact ⟨v, hvt, ht hvt, hvbad⟩

theorem minusEight_exists_unspoiled_of_sieve {n : ℕ} {X lower δ E L : ℝ}
    (t : Finset MinusEightShellPoint) (ht : t ⊆ minusEightDyadicShell n)
    (hn : 0 < (minusEightFullShell n).card) (hX : 0 < X) (hE : 0 ≤ E) (hL : 0 ≤ L)
    (hselected : (lower - δ * E) * (minusEightFullShell n).card ≤ t.card)
    (hsmall : ∀ p ≤ n, MinusEightInertPrime p → (p : ℝ) ≤ X →
      ((minusEightEitherDivisible n p).card : ℝ) ≤
        (minusEightLocalMass p + E) * (minusEightFullShell n).card)
    (hlarge : ∀ p ≤ n, MinusEightInertPrime p → X < (p : ℝ) →
      ((t ∩ minusEightEitherDivisible n p).card : ℝ) ≤
        L * inverseSquare p * (minusEightFullShell n).card)
    (hmargin : 3 / 16 + E * (X / 4 + 1 + δ) + L * (1 / (4 * X) + 2 / X ^ 2) < lower) :
    ∃ v ∈ t, v ∈ minusEightDyadicShell n ∧
      ∀ p ≤ n, MinusEightInertPrime p → v ∉ minusEightEitherDivisible n p := by

  exact minusEight_exists_unspoiled_of_cutoff_bounds t ht hn
    (minusEight_cutoff_bounds_progressions hX) hE hL hselected hsmall hlarge hmargin

private theorem minusEight_small_prime_bad_bound {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m p : ℕ} (he : e ≤ 3) (hm : Odd m) (hp : MinusEightInertPrime p)
    (hsmall : ¬ p ∣ 2 ^ e * m →
      ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card) :
    ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
      (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card := by
  by_cases hd : p ∣ 2 ^ e * m
  · refine (minusEight_inert_target_bad_bound hθ he hm hp hd).trans ?_
    calc
      _ = inverseSquare p * (minusEightFullShell (2 ^ e * m)).card := by unfold inverseSquare; ring
      _ ≤ _ := mul_le_mul_of_nonneg_right ((minusEight_inverse_square_le_local_mass hp).trans
        (le_add_of_nonneg_right (minusEight_shell_error_nonneg m))) (by positivity)
  · exact hsmall hd

private theorem minusEight_full_shell_pos {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m : ℕ} (he : e ≤ 3) (hm : Odd m) : 0 < (minusEightFullShell (2 ^ e * m)).card := by
  have hm₀ : 0 < m := by obtain ⟨k, hk⟩ := hm; omega
  have hs : (0 : ℝ) < σ 1 m := by exact_mod_cast ArithmeticFunction.sigma_pos 1 m hm₀.ne'
  have hc : 0 < minusEightShellFactor e := by unfold minusEightShellFactor; split_ifs <;> decide
  apply (Nat.cast_pos (α := ℝ)).mp
  rw [minusEight_full_shell_dyadic hθ he hm]
  positivity

/-- The uniform sieve with the dyadic shell. Its large-prime input is
a bound on the actual selected-shell bad count. -/
theorem minusEight_sum_of_sieve_cutoff_bounds {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m : ℕ} (he : e ≤ 3) (hm : Odd m) (ha : |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m)
    {X count tail L : ℝ} (hcut : MinusEightCutoffBounds X count tail) (hL : 0 ≤ L)
    (hsmall : ∀ p ≤ 2 ^ e * m, MinusEightInertPrime p → (p : ℝ) ≤ X → ¬ p ∣ 2 ^ e * m →
      ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card)
    (hlarge : ∀ p ≤ 2 ^ e * m, MinusEightInertPrime p → X < (p : ℝ) →
      (((minusEightDyadicShell (2 ^ e * m)) ∩ minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        L * inverseSquare p * (minusEightFullShell (2 ^ e * m)).card)
    (hmargin : 3 / 16 + minusEightShellError m * (count + 1 / 2) +
      L * tail < 1 / 2) : InSumset P8 P8 (2 ^ e * m) := by
  let n := 2 ^ e * m
  obtain ⟨v, _, hv, huv⟩ := minusEight_exists_unspoiled_of_cutoff_bounds
    (lower := 1 / 2) (δ := 1 / 2) (minusEightDyadicShell n) (Finset.Subset.refl _)
    (minusEight_full_shell_pos hθ he hm) hcut (minusEight_shell_error_nonneg m) hL
    (by simpa only [div_eq_mul_inv, one_mul, mul_comm (2⁻¹ : ℝ) (minusEightShellError m)] using
      minusEight_dyadic_shell_lower hθ he hm ha)
    (fun p hpn hp hPX => minusEight_small_prime_bad_bound hθ he hm hp (hsmall p hpn hp hPX))
    hlarge (by linarith only [hmargin])
  exact minusEight_sum_of_unspoiled hv huv

theorem minusEight_sum_of_uniform_sieve {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m : ℕ} (he : e ≤ 3) (hm : Odd m) (ha : |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m)
    {X L : ℝ} (hX : 0 < X) (hL : 0 ≤ L)
    (hsmall : ∀ p ≤ 2 ^ e * m, MinusEightInertPrime p → (p : ℝ) ≤ X → ¬ p ∣ 2 ^ e * m →
      ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card)
    (hlarge : ∀ p ≤ 2 ^ e * m, MinusEightInertPrime p → X < (p : ℝ) →
      (((minusEightDyadicShell (2 ^ e * m)) ∩ minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        L * inverseSquare p * (minusEightFullShell (2 ^ e * m)).card)
    (hmargin : 3 / 16 + minusEightShellError m * (X / 4 + 3 / 2) +
      L * (1 / (4 * X) + 2 / X ^ 2) < 1 / 2) : InSumset P8 P8 (2 ^ e * m) := by

  exact minusEight_sum_of_sieve_cutoff_bounds hθ he hm ha (minusEight_cutoff_bounds_progressions hX)
    hL hsmall hlarge (by simpa only [add_assoc, show (1 : ℝ) + 1 / 2 = 3 / 2 by norm_num] using hmargin)

private theorem minusEight_bad_bound_normalize {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m p : ℕ} (he : e ≤ 3) (hm : Odd m) {K B : ℝ} (hK : 0 ≤ K)
    (hB : B ≤ (4 * K) * minusEightDivisorMaximum (2 ^ e * m) * (2 ^ e * m : ℕ) / (p : ℝ) ^ 2) :
    B ≤ (K * (2 : ℝ) ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum (2 ^ e * m)) *
      inverseSquare p * (minusEightFullShell (2 ^ e * m)).card := by
  have hm₀ : 0 < m := by obtain ⟨k, hk⟩ := hm; omega
  have hc : 0 < minusEightShellFactor e := by unfold minusEightShellFactor; split_ifs <;> decide
  have hcR : 0 < (minusEightShellFactor e : ℝ) := by exact_mod_cast hc
  have hσN : m ≤ σ 1 m := by
    simpa using pow_mul_sigma_le (k := 1) hm₀ Nat.one_pos
  have hσ : (m : ℝ) ≤ σ 1 m := by exact_mod_cast hσN
  rw [minusEight_full_shell_dyadic hθ he hm]
  refine hB.trans ?_
  push_cast
  unfold inverseSquare
  calc
    _ = (K * (2 : ℝ) ^ e / (minusEightShellFactor e : ℝ) *
        minusEightDivisorMaximum (2 ^ e * m)) * (1 / (p : ℝ) ^ 2) *
        (4 * (minusEightShellFactor e : ℝ) * m) := by field_simp [hcR.ne']
    _ ≤ _ := by gcongr

/-- The uniform sieve with the lattice and signed representation
inputs discharged. The remaining bad-count input is at unit levels. -/
theorem minusEight_sum_of_arithmetic_cutoff_bounds {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m : ℕ} (he : e ≤ 3) (hm : Odd m) (ha : |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m)
    {X count tail : ℝ} (hcut : MinusEightCutoffBounds X count tail)
    (hsmall : ∀ p ≤ 2 ^ e * m, MinusEightInertPrime p → (p : ℝ) ≤ X → ¬ p ∣ 2 ^ e * m →
      ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card)
    (hmargin : 3 / 16 + minusEightShellError m * (count + 1 / 2) +
      ((8 / 3 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum (2 ^ e * m)) *
        tail < 1 / 2) : InSumset P8 P8 (2 ^ e * m) := by
  apply minusEight_sum_of_sieve_cutoff_bounds hθ he hm ha hcut (by positivity) hsmall _ hmargin
  intro p _ hp _
  exact minusEight_bad_bound_normalize hθ he hm (by norm_num : (0 : ℝ) ≤ 8 / 3)
    (by simpa only [show (4 : ℝ) * (8 / 3) = 32 / 3 by norm_num] using minusEight_inert_dyadic_bad_card_le hp)

theorem minusEight_sum_of_arithmetic_margin {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m : ℕ} (he : e ≤ 3) (hm : Odd m) (ha : |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m)
    {X : ℝ} (hX : 0 < X)
    (hsmall : ∀ p ≤ 2 ^ e * m, MinusEightInertPrime p → (p : ℝ) ≤ X → ¬ p ∣ 2 ^ e * m →
      ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card)
    (hmargin : 3 / 16 + minusEightShellError m * (X / 4 + 3 / 2) +
      ((8 / 3 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum (2 ^ e * m)) *
        (1 / (4 * X) + 2 / X ^ 2) < 1 / 2) : InSumset P8 P8 (2 ^ e * m) := by

  exact minusEight_sum_of_arithmetic_cutoff_bounds hθ he hm ha
    (minusEight_cutoff_bounds_progressions hX) hsmall
    (by simpa only [add_assoc, show (1 : ℝ) + 1 / 2 = 3 / 2 by norm_num] using hmargin)

/-- The arithmetic sieve lower bound for the number of retained coordinate
vectors. Its only small-prime hypothesis concerns primes coprime to the target. -/
theorem minusEight_uniform_survivor_card_le {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m : ℕ} (he : e ≤ 3) (hm : Odd m) (ha : |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m)
    {X : ℝ} (hX : 0 < X)
    (hsmall : ∀ p ≤ 2 ^ e * m, MinusEightInertPrime p → (p : ℝ) ≤ X → ¬ p ∣ 2 ^ e * m →
      ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card) :
    let n := 2 ^ e * m
    let L := (8 / 3 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum n
    (1 / 2 - (3 / 16 + minusEightShellError m * (X / 4 + 3 / 2) +
      L * (1 / (4 * X) + 2 / X ^ 2))) * (minusEightFullShell n).card ≤
      (minusEightSieveSurvivors n (minusEightDyadicShell n)).card := by
  have h := minusEight_sieve_survivor_card_le (lower := 1 / 2) (δ := 1 / 2)
    (L := (8 / 3 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum (2 ^ e * m))
    (minusEightDyadicShell (2 ^ e * m)) (minusEight_cutoff_bounds_progressions hX)
    (minusEight_shell_error_nonneg m) (by positivity)
    (by simpa only [div_eq_mul_inv, one_mul, mul_comm (2⁻¹ : ℝ) (minusEightShellError m)] using
      minusEight_dyadic_shell_lower hθ he hm ha)
    (fun p hpn hp hPX => minusEight_small_prime_bad_bound hθ he hm hp (hsmall p hpn hp hPX))
    (fun p _ hp _ => minusEight_bad_bound_normalize hθ he hm (by norm_num : (0 : ℝ) ≤ 8 / 3)
      (by simpa only [show (4 : ℝ) * (8 / 3) = 32 / 3 by norm_num] using minusEight_inert_dyadic_bad_card_le hp))
  simpa only [add_assoc, show (1 : ℝ) + 1 / 2 = 3 / 2 by norm_num] using h

/-- The selected odd-odd sieve produces ordinary core pairs in the two
even strata. Its normalized large-prime coefficients are 2 and 4/3. -/
theorem minusEight_core_sum_of_selected_cutoff_bounds {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m : ℕ} (he : e = 1 ∨ e = 2) (hm : Odd m) {X count tail : ℝ} (hcut : MinusEightCutoffBounds X count tail)
    (hsmall : ∀ p ≤ 2 ^ e * m, MinusEightInertPrime p → (p : ℝ) ≤ X → ¬ p ∣ 2 ^ e * m →
      ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card)
    (hmargin : 3 / 16 + minusEightShellError m * (count) +
      ((2 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum (2 ^ e * m)) *
        tail < (if e = 1 then 1 / 2 else 2 / 3)) :
    InSumset C8 C8 (2 ^ e * m) := by
  let n := 2 ^ e * m
  have he₃ : e ≤ 3 := by omega
  have hselected : (if e = 1 then (1 / 2 : ℝ) else 2 / 3) * (minusEightFullShell n).card ≤
      (minusEightOddOddShell n).card := by
    have h := minusEight_full_shell_dyadic hθ he₃ hm
    obtain ⟨h₁, h₂⟩ := minusEight_selected_even_shell_theta hθ hm
    rcases he with rfl | rfl
    · dsimp only [n]
      rw [h]
      norm_num [minusEightShellFactor]
      linarith only [h₁]
    · dsimp only [n]
      rw [h]
      norm_num [minusEightShellFactor]
      linarith only [h₂]
  obtain ⟨v, hv, _, h⟩ := minusEight_exists_unspoiled_of_cutoff_bounds (δ := 0) (minusEightOddOddShell n)
    (minusEight_odd_odd_shell_subset n) (minusEight_full_shell_pos hθ he₃ hm)
    hcut (minusEight_shell_error_nonneg m)
    (L := (2 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum n)
    (by positivity) (by simpa only [zero_mul, sub_zero] using hselected)
    (fun p hpn hp hPX => minusEight_small_prime_bad_bound hθ he₃ hm hp (hsmall p hpn hp hPX))
    (fun p _ hp _ => minusEight_bad_bound_normalize hθ he₃ hm (by norm_num : (0 : ℝ) ≤ 2)
      (by simpa only [show (4 : ℝ) * 2 = 8 by norm_num] using minusEight_inert_odd_odd_bad_card_le hp))
    (by simpa only [add_zero] using hmargin)
  exact minusEight_core_sum_of_unspoiled_odd_odd hv h

theorem minusEight_core_sum_of_selected_margin {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {e m : ℕ} (he : e = 1 ∨ e = 2) (hm : Odd m) {X : ℝ} (hX : 0 < X)
    (hsmall : ∀ p ≤ 2 ^ e * m, MinusEightInertPrime p → (p : ℝ) ≤ X → ¬ p ∣ 2 ^ e * m →
      ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card)
    (hmargin : 3 / 16 + minusEightShellError m * (X / 4 + 1) +
      ((2 : ℝ) * 2 ^ e / (minusEightShellFactor e : ℝ) * minusEightDivisorMaximum (2 ^ e * m)) *
        (1 / (4 * X) + 2 / X ^ 2) < (if e = 1 then 1 / 2 else 2 / 3)) :
    InSumset C8 C8 (2 ^ e * m) := by

  exact minusEight_core_sum_of_selected_cutoff_bounds hθ he hm
    (minusEight_cutoff_bounds_progressions hX) hsmall hmargin

/-- The odd selected sieve produces a mixed core pair with the
large-prime coefficient 8/3. -/
theorem minusEight_mixed_core_of_selected_cutoff_bounds {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {m : ℕ} (hm : Odd m) (ha : |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m)
    {X count tail : ℝ} (hcut : MinusEightCutoffBounds X count tail)
    (hsmall : ∀ p ≤ m, MinusEightInertPrime p → (p : ℝ) ≤ X → ¬ p ∣ m →
      ((minusEightEitherDivisible m p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell m).card)
    (hmargin : 3 / 16 + minusEightShellError m * (count + 1 / 2) +
      ((8 / 3 : ℝ) * minusEightDivisorMaximum m) * tail < 1 / 2) :
    InSumset C8 (Dilate 2 C8) m := by
  have hsum : InSumset P8 P8 m := by
    simpa only [pow_zero, one_mul] using minusEight_sum_of_arithmetic_cutoff_bounds hθ
      (e := 0) (by decide) hm ha hcut (by simpa using hsmall) (by simpa [minusEightShellFactor] using hmargin)
  obtain ⟨x, y, hx, hy, he⟩ := hsum
  obtain ⟨j, hj⟩ := hm
  rcases (P8_iff_core x).mp hx with hx | ⟨u, hu, hxu⟩ <;>
    rcases (P8_iff_core y).mp hy with hy | ⟨v, hv, hyv⟩
  · have h₁ := C8_mod_eight hx
    have h₂ := C8_mod_eight hy
    omega
  · exact ⟨x, y, hx, ⟨v, hv, hyv⟩, he⟩
  · exact ⟨y, x, hy, ⟨u, hu, hxu⟩, by omega⟩
  · omega

theorem minusEight_mixed_core_of_selected_margin {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {m : ℕ} (hm : Odd m) (ha : |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m)
    {X : ℝ} (hX : 0 < X)
    (hsmall : ∀ p ≤ m, MinusEightInertPrime p → (p : ℝ) ≤ X → ¬ p ∣ m →
      ((minusEightEitherDivisible m p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell m).card)
    (hmargin : 3 / 16 + minusEightShellError m * (X / 4 + 3 / 2) +
      ((8 / 3 : ℝ) * minusEightDivisorMaximum m) * (1 / (4 * X) + 2 / X ^ 2) < 1 / 2) :
    InSumset C8 (Dilate 2 C8) m := by

  exact minusEight_mixed_core_of_selected_cutoff_bounds hθ hm ha
    (minusEight_cutoff_bounds_progressions hX) hsmall
    (by simpa only [add_assoc, show (1 : ℝ) + 1 / 2 = 3 / 2 by norm_num] using hmargin)

end QuadraticNorms
