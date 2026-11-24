import QuadraticNorms.EulerDivisorBounds
import QuadraticNorms.MinusElevenFrontier
import QuadraticNorms.MinusElevenPrimeBounds
import Mathlib.NumberTheory.ArithmeticFunction.Misc
import Mathlib.Analysis.Real.Sqrt
import Mathlib.Data.Finset.Lattice.Fold

open scoped BigOperators ArithmeticFunction.sigma

namespace QuadraticNorms

@[irreducible] private def smallSplitPrimes : Finset ℕ := (Finset.range 64).filter MinusElevenSplitPrime

private theorem mem_smallSplitPrimes {p : ℕ} : p ∈ smallSplitPrimes ↔ p < 64 ∧ MinusElevenSplitPrime p := by
  unfold smallSplitPrimes
  simp

private theorem divisor_cap_product :
    (∏ p ∈ smallSplitPrimes, boundedEulerCap (divisorFactor 6) 10 p) <
      (7 : ℚ) ^ 6 := by decide +kernel

private theorem split_factor_product_lt {n : ℕ} :
    (∏ p ∈ n.primeFactors.filter MinusElevenSplitPrime, divisorFactor 6 p (n.factorization p)) <
      (7 : ℚ) ^ 6 := by
  let s := n.primeFactors.filter MinusElevenSplitPrime
  have hSplit (p : ℕ) (hp : p ∈ s) : MinusElevenSplitPrime p := (Finset.mem_filter.mp hp).2
  have hprod := euler_product_le_small_primes (divisorFactor 6) 10 (support := s)
    (small := smallSplitPrimes) n.factorization
    (by intro p; simp [divisorFactor])
    (fun p hp hsmall => mem_smallSplitPrimes.mpr ⟨hsmall, hSplit p hp⟩)
    (fun p hp => divisor_factor_le_cap (by decide) (hSplit p hp).1.two_le _)
    (fun p _ => by unfold divisorFactor; positivity)
  exact lt_of_le_of_lt hprod divisor_cap_product

/-- The unrestricted split-divisor bound in sixth-power form. -/
theorem minusEleven_split_divisors_power_bound {n : ℕ} (hn : 0 < n) :
    (minusElevenSplitDivisors n : ℚ) ^ 6 < (7 : ℚ) ^ 6 * n := by
  let s := n.primeFactors.filter MinusElevenSplitPrime
  let e := n.factorization
  have h := split_factor_product_lt (n := n)
  unfold divisorFactor at h
  rw [Finset.prod_div_distrib, Finset.prod_pow] at h
  have hpos : 0 < ∏ p ∈ s, (p : ℚ) ^ e p := Finset.prod_pos (fun p hp => by
    have hp' := ((Finset.mem_filter.mp hp).2 : MinusElevenSplitPrime p).1.pos
    exact pow_pos (Nat.cast_pos.mpr hp') _)
  have hle : ∏ p ∈ s, (p : ℚ) ^ e p ≤ n := by
    have hall := Nat.prod_factorization_pow_eq_self hn.ne'
    rw [Finsupp.prod, Nat.support_factorization] at hall
    have hnat : ∏ p ∈ s, p ^ e p ≤ n := by
      rw [← hall]
      exact Finset.prod_le_prod_of_subset_of_one_le' (Finset.filter_subset _ _)
        (fun p hp _ => Nat.one_le_pow _ _ (Nat.prime_of_mem_primeFactors hp).pos)
    exact_mod_cast hnat
  have hv : (minusElevenSplitDivisors n : ℚ) = ∏ p ∈ s, (((e p : ℕ) : ℚ) + 1) := by
    unfold minusElevenSplitDivisors
    rw [Nat.cast_prod]
    exact Finset.prod_congr rfl (fun p _ => by simp only [e, Nat.cast_add, Nat.cast_one])
  rw [hv]
  change (∏ p ∈ s, (((e p : ℕ) : ℚ) + 1)) ^ 6 / (∏ p ∈ s, (p : ℚ) ^ e p) < (7 : ℚ) ^ 6 at h
  exact ((div_lt_iff₀ hpos).mp h).trans_le
    (mul_le_mul_of_nonneg_left hle (by positivity : 0 ≤ (7 : ℚ) ^ 6))

private theorem inert_factor_step {p a : ℕ} (hp : MinusElevenInertPrime p) (ha : 2 ≤ a) :
    sigmaRatioFactor p (a + 1) ≤ sigmaRatioFactor p a := by
  have hpN : 7 ≤ p := by
    by_contra! h
    have hex : ∀ q : Fin 7, ¬ MinusElevenInertPrime q.val := by decide
    exact hex ⟨p, h⟩ hp
  have hpQ : (7 : ℚ) ≤ p := by exact_mod_cast hpN
  have haQ : (2 : ℚ) ≤ a := by exact_mod_cast ha
  exact sigma_ratio_factor_step hp.1 (A := 4 / 3) (by linarith)
    ((by norm_num : (4 / 3 : ℚ) ^ 6 ≤ 7).trans hpQ)

@[irreducible] private def smallInertPrimes : Finset ℕ := (Finset.range 64).filter MinusElevenInertPrime

private theorem mem_smallInertPrimes {p : ℕ} : p ∈ smallInertPrimes ↔ p < 64 ∧ MinusElevenInertPrime p := by
  unfold smallInertPrimes
  simp

private theorem inert_cap_product :
    (∏ p ∈ smallInertPrimes, boundedEulerCap sigmaRatioFactor 2 p) < (11 / 4 : ℚ) ^ 6 := by
  decide +kernel

private theorem inert_factor_product_lt {m : ℕ}
    (hm : ∀ p ∈ m.primeFactors, MinusElevenInertPrime p) :
    (∏ p ∈ m.primeFactors, sigmaRatioFactor p (m.factorization p)) < (11 / 4 : ℚ) ^ 6 := by
  refine (euler_product_le_small_primes sigmaRatioFactor 2 (small := smallInertPrimes) m.factorization
    (by intro p; simp [sigmaRatioFactor])
    (fun p hp hsmall => mem_smallInertPrimes.mpr ⟨hsmall, hm p hp⟩) ?_ ?_).trans_lt inert_cap_product
  · intro p hp
    unfold boundedEulerCap
    by_cases h : p < 64
    · rw [if_pos h]
      exact factor_le_finite_max _ (fun _ ha => inert_factor_step (hm p hp) ha) _
    · rw [if_neg h]
      exact sigma_ratio_factor_large (hm p hp).1 (Nat.le_of_not_gt h)
  · intro p _
    unfold sigmaRatioFactor
    positivity


private def dyadicFactor (t : ℕ) : ℚ := (2 * ((t : ℚ) + 3) / 9) ^ 6 / (2 : ℚ) ^ t

private theorem dyadic_factor_step {t : ℕ} (ht : 6 ≤ t) :
    dyadicFactor (t + 1) ≤ dyadicFactor t := by
  have htQ : (6 : ℚ) ≤ t := by exact_mod_cast ht
  have hr : 2 * ((t : ℚ) + 4) / 9 ≤ (10 / 9) * (2 * ((t : ℚ) + 3) / 9) := by linarith
  have h := pow_le_pow_left₀ (show 0 ≤ 2 * ((t : ℚ) + 4) / 9 by positivity) hr 6
  rw [mul_pow] at h
  have h₂ := h.trans (mul_le_mul_of_nonneg_right (by norm_num : (10 / 9 : ℚ) ^ 6 ≤ 2)
    (by positivity : 0 ≤ (2 * ((t : ℚ) + 3) / 9) ^ 6))
  unfold dyadicFactor
  push_cast
  rw [show (t : ℚ) + 1 + 3 = (t : ℚ) + 4 by ring, pow_succ (2 : ℚ),
    div_le_div_iff₀ (by positivity) (by positivity)]
  simpa only [mul_left_comm, mul_comm, mul_assoc] using
    mul_le_mul_of_nonneg_right h₂ (by positivity : 0 ≤ (2 : ℚ) ^ t)

private theorem dyadic_factor_le_one (t : ℕ) : dyadicFactor t ≤ 1 := by
  induction t with
  | zero => norm_num [dyadicFactor]
  | succ t ih =>
    by_cases ht : 6 ≤ t
    · exact (dyadic_factor_step ht).trans ih
    · have h : ∀ a : Fin 7, dyadicFactor a.val ≤ 1 := by decide +kernel
      exact h ⟨t + 1, by omega⟩

/-- The error in the odd-shell estimate for the target 2^t * m.
The arithmetic functions σ 0 and σ 1 are the divisor count and divisor sum. -/
noncomputable def minusElevenFilteredError (t m : ℕ) : ℝ :=
  (2 * ((t : ℝ) + 3) / 9) * (σ 0 m : ℝ) * Real.sqrt m /
    (Real.sqrt ((2 : ℝ) ^ t) * (σ 1 m : ℝ))

theorem minusEleven_filtered_error_nonneg (t m : ℕ) : 0 ≤ minusElevenFilteredError t m := by
  unfold minusElevenFilteredError
  positivity

/-- The sharp Euler bound E^6 * n^2 < (11/4)^6 on inert cores. -/
theorem minusEleven_filtered_error_power_bound {t m : ℕ} (hm : 0 < m)
    (hI : ∀ p ∈ m.primeFactors, MinusElevenInertPrime p) :
    minusElevenFilteredError t m ^ 6 * ((2 : ℝ) ^ t * m) ^ 2 < (11 / 4 : ℝ) ^ 6 := by
  have hE := inert_factor_product_lt hI
  rw [sigma_ratio_factor_product_eq hm] at hE
  have hproduct : dyadicFactor t * ((σ 0 m : ℚ) ^ 6 * (m : ℚ) ^ 5 / (σ 1 m : ℚ) ^ 6) <
      (11 / 4 : ℚ) ^ 6 :=
    ((mul_le_mul_of_nonneg_right (dyadic_factor_le_one t) (by positivity)).trans_eq
      (one_mul _)).trans_lt hE
  have hcast := (Rat.cast_lt (K := ℝ)).mpr hproduct
  unfold dyadicFactor at hcast
  push_cast at hcast
  have hS : (0 : ℝ) < σ 1 m := by exact_mod_cast ArithmeticFunction.sigma_pos 1 m hm.ne'
  have hT : (0 : ℝ) < (2 : ℝ) ^ t := by positivity
  have hroot (x : ℝ) (hx : 0 ≤ x) : Real.sqrt x ^ 6 = x ^ 3 := by
    conv_lhs => rw [show (6 : ℕ) = 2 * 3 by decide, pow_mul, Real.sq_sqrt hx]
  have heq : minusElevenFilteredError t m ^ 6 * ((2 : ℝ) ^ t * m) ^ 2 =
      (2 * ((t : ℝ) + 3) / 9) ^ 6 / (2 : ℝ) ^ t *
        ((σ 0 m : ℝ) ^ 6 * (m : ℝ) ^ 5 / (σ 1 m : ℝ) ^ 6) := by
    unfold minusElevenFilteredError
    simp only [div_pow, mul_pow, hroot (m : ℝ) (Nat.cast_nonneg _), hroot _ hT.le]
    field_simp [hS.ne', hT.ne']
  rw [heq]
  exact hcast

end QuadraticNorms
