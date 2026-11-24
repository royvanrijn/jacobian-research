import MathieuQuadratic.Endpoint
import Mathlib.Algebra.Polynomial.Reverse

/-!
# The endpoint coefficients for the simple-pole phase

A polynomial differential identity computes the same first omitted
term as the reciprocal-square-root series in the draft.
-/

open Polynomial
open scoped BigOperators

namespace MathieuQuadratic

/-- The coefficient of degree k in the binomial series with exponent −1/2. -/
noncomputable def negativeHalfBinomial (k : ℕ) : ℚ :=
  (-1) ^ k * (ascPochhammer ℚ k).eval (1 / 2) / (k.factorial : ℚ)

private theorem negativeHalfBinomial_succ (k : ℕ) :
    ((k : ℚ) + 1) * negativeHalfBinomial (k + 1) =
      -((k : ℚ) + 1 / 2) * negativeHalfBinomial k := by
  rw [negativeHalfBinomial, negativeHalfBinomial, pow_succ, ascPochhammer_succ_eval,
    Nat.factorial_succ, Nat.cast_mul, Nat.cast_succ]
  have hk : (k.factorial : ℚ) ≠ 0 := by exact_mod_cast Nat.factorial_ne_zero k
  have hs : (k : ℚ) + 1 ≠ 0 := by positivity
  field_simp [hk, hs]
  ring

private noncomputable def halfBinomialTruncation (m : ℕ) : ℚ[X] :=
  ∑ k ∈ Finset.range (m + 1), monomial k (negativeHalfBinomial k)

private theorem halfBinomialTruncation_coeff (m k : ℕ) :
    (halfBinomialTruncation m).coeff k = if k ≤ m then negativeHalfBinomial k else 0 := by
  classical
  simp only [halfBinomialTruncation, finsetSum_coeff, coeff_monomial]
  by_cases hk : k ≤ m
  · rw [if_pos hk, Finset.sum_eq_single_of_mem k (Finset.mem_range.mpr (by omega)), if_pos rfl]
    intro j _ hj
    rw [if_neg hj]
  · rw [if_neg hk]
    apply Finset.sum_eq_zero
    intro j hj
    rw [if_neg (by
      have := Finset.mem_range.mp hj
      omega : j ≠ k)]

private theorem halfBinomialTruncation_derivative (m : ℕ) :
    (1 + X) * (halfBinomialTruncation m).derivative + C (1 / 2 : ℚ) * halfBinomialTruncation m =
      C (((m : ℚ) + 1 / 2) * negativeHalfBinomial m) * X ^ m := by
  let G := halfBinomialTruncation m
  have hx (k : ℕ) : (X * G.derivative).coeff k = (k : ℚ) * G.coeff k := by
    cases k with
    | zero => simp
    | succ k =>
      rw [coeff_X_mul, coeff_derivative, mul_comm, Nat.cast_succ]
  ext k
  change ((1 + X) * G.derivative + C (1 / 2 : ℚ) * G).coeff k = _
  rw [add_mul, one_mul, coeff_add, coeff_add, coeff_C_mul, coeff_C_mul, coeff_X_pow,
    hx, coeff_derivative]
  dsimp only [G]
  rw [halfBinomialTruncation_coeff, halfBinomialTruncation_coeff]
  by_cases hk : k < m
  · rw [if_pos (by omega : k + 1 ≤ m), if_pos hk.le, if_neg (by omega : k ≠ m), mul_zero]
    linear_combination negativeHalfBinomial_succ k
  · by_cases he : k = m
    · subst k
      rw [if_neg (by omega : ¬ m + 1 ≤ m), if_pos le_rfl, if_pos rfl]
      ring
    · rw [if_neg (by omega : ¬ k + 1 ≤ m), if_neg (by omega : ¬ k ≤ m), if_neg he]
      ring

/-- The reversed coefficient polynomial for the linear-marker phase. -/
noncomputable def linearPhaseEndpoint (m : ℕ) : ℚ[X] :=
  (1 + X) ^ m * (halfBinomialTruncation m).comp (X * (2 + X))

private noncomputable def linearPhaseAntiderivative (m : ℕ) : ℚ[X] :=
  (1 + X) * (halfBinomialTruncation m).comp (X * (2 + X))

private theorem linearPhaseAntiderivative_derivative (m : ℕ) :
    (linearPhaseAntiderivative m).derivative =
      X ^ m * (C ((2 * (m : ℚ) + 1) * negativeHalfBinomial m) * (2 + X) ^ m) := by
  let P : ℚ[X] := X * (2 + X)
  have hd : P.derivative = 2 * (1 + X) := by simp [P]; ring
  have hC : (1 + P) = (1 + X) ^ 2 := by dsimp [P]; ring
  have h := congrArg (fun p : ℚ[X] ↦ p.comp P) (halfBinomialTruncation_derivative m)
  simp only [add_comp, mul_comp, one_comp, X_comp, C_comp, pow_comp] at h
  unfold linearPhaseAntiderivative
  change derivative ((1 + X) * (halfBinomialTruncation m).comp P) = _
  rw [derivative_mul, derivative_comp, derivative_add, derivative_one, derivative_X, zero_add,
    one_mul, hd]
  calc
    _ = 2 * ((1 + P) * (halfBinomialTruncation m).derivative.comp P +
      C (1 / 2 : ℚ) * (halfBinomialTruncation m).comp P) := by
        have hc : (2 : ℚ[X]) * C (1 / 2 : ℚ) = 1 := by
          calc
            _ = C ((2 : ℚ) * (1 / 2)) := by rw [map_mul, map_ofNat]
            _ = _ := by norm_num
        rw [hC]
        symm
        rw [mul_add, ← mul_assoc 2 (C (1 / 2 : ℚ)), hc, one_mul]
        ring
    _ = _ := by
      rw [h, ← mul_assoc]
      have hc : (2 : ℚ[X]) * C (((m : ℚ) + 1 / 2) * negativeHalfBinomial m) =
          C ((2 * (m : ℚ) + 1) * negativeHalfBinomial m) := by
        rw [← map_ofNat C 2, ← map_mul]
        congr 1
        ring
      rw [hc]
      dsimp [P]
      rw [mul_pow]
      ring

private theorem linearPhaseAntiderivative_eval_zero (m : ℕ) :
    (linearPhaseAntiderivative m).eval 0 = 1 := by
  simp only [linearPhaseAntiderivative, eval_mul, eval_add, eval_one, eval_X, add_zero,
    eval_comp, zero_mul, one_mul]
  rw [← coeff_zero_eq_eval_zero, halfBinomialTruncation_coeff, if_pos (Nat.zero_le _)]
  simp [negativeHalfBinomial]

/-- The omitted part begins in degree m+1 and has the stated first coefficient. -/
private theorem linearPhaseAntiderivative_error (m : ℕ) :
    ∃ e : ℚ[X], linearPhaseAntiderivative m = 1 + X ^ (m + 1) * e ∧
      e.coeff 0 = (2 * (m : ℚ) + 1) * negativeHalfBinomial m * 2 ^ m / (m + 1 : ℕ) := by
  have hdiv : X ^ (m + 1) ∣ linearPhaseAntiderivative m - 1 := by
    apply X_pow_succ_dvd_of_derivative
    · rw [coeff_zero_eq_eval_zero, eval_sub, eval_one, linearPhaseAntiderivative_eval_zero, sub_self]
    · rw [derivative_sub, derivative_one, sub_zero, linearPhaseAntiderivative_derivative]
      exact dvd_mul_right _ _
  obtain ⟨e, he⟩ := hdiv
  have heq : linearPhaseAntiderivative m = 1 + X ^ (m + 1) * e := by rw [← he]; ring
  refine ⟨e, heq, ?_⟩
  have hc := congrArg (fun p : ℚ[X] ↦ p.coeff m) (linearPhaseAntiderivative_derivative m)
  rw [coeff_derivative, heq, coeff_add, coeff_one, if_neg (by omega), zero_add,
    coeff_X_pow_mul', if_pos le_rfl, Nat.sub_self, coeff_X_pow_mul', if_pos le_rfl, Nat.sub_self,
    coeff_C_mul] at hc
  have hcoef : ((2 + X : ℚ[X]) ^ m).coeff 0 = (2 : ℚ) ^ m := by
    rw [coeff_zero_eq_eval_zero, eval_pow, eval_add, eval_ofNat, eval_X, add_zero]
  rw [hcoef] at hc
  apply (eq_div_iff (by positivity)).2
  simpa only [Nat.cast_succ] using hc

/-- The two coefficients required for the simple-pole witness. -/
theorem linearPhaseEndpoint_coefficients (m : ℕ) (hm : 0 < m) :
    (linearPhaseEndpoint m).coeff m = 0 ∧
      (linearPhaseEndpoint m).coeff (m + 1) =
        (2 * (m : ℚ) + 1) * negativeHalfBinomial m * 2 ^ m / (m + 1 : ℕ) := by
  obtain ⟨e, he, hc⟩ := linearPhaseAntiderivative_error m
  have hp : linearPhaseEndpoint m = (1 + X) ^ (m - 1) * linearPhaseAntiderivative m := by
    rw [linearPhaseEndpoint, linearPhaseAntiderivative, ← mul_assoc, ← pow_succ,
      Nat.sub_add_cancel hm]
  have hdeg : ((1 + X : ℚ[X]) ^ (m - 1)).natDegree ≤ m - 1 := by
    refine natDegree_pow_le.trans ?_
    have hlin : (1 + X : ℚ[X]).natDegree ≤ 1 := by
      exact (natDegree_add_le _ _).trans (by simp)
    simpa using Nat.mul_le_mul_left (m - 1) hlin
  rw [hp, he, mul_add, mul_one, mul_left_comm _ (X ^ (m + 1))]
  constructor
  · rw [coeff_add, coeff_eq_zero_of_natDegree_lt (by omega), coeff_X_pow_mul',
      if_neg (by omega), add_zero]
  · rw [coeff_add, coeff_eq_zero_of_natDegree_lt (by omega), zero_add,
      coeff_X_pow_mul', if_pos le_rfl, Nat.sub_self, coeff_zero_eq_eval_zero, eval_mul,
      eval_pow, eval_add, eval_one, eval_X, add_zero, one_pow, one_mul,
      ← coeff_zero_eq_eval_zero, hc]

private theorem reflect_power_linear (p : ℚ[X]) (hp : p.natDegree ≤ 1) :
    ∀ j : ℕ, reflect j (p ^ j) = (reflect 1 p) ^ j := by
  intro j
  induction j with
  | zero => rw [pow_zero, pow_zero, ← C_1, reflect_C, pow_zero, mul_one]
  | succ j ih =>
    rw [pow_succ p, reflect_mul _ _ (natDegree_pow_le.trans (by simpa using Nat.mul_le_mul_left j hp)) hp,
      ih, ← pow_succ]

private theorem reflect_one_add_X : reflect 1 (1 + X : ℚ[X]) = 1 + X := by
  have hx : reflect 1 (X : ℚ[X]) = 1 := by
    simpa only [pow_one, revAt_le le_rfl, Nat.sub_self, pow_zero] using
      (reflect_monomial (R := ℚ) 1 1)
  rw [reflect_add, ← C_1, reflect_C, hx]
  simp [add_comm]

private theorem reflect_one_add_twoX : reflect 1 (1 + 2 * X : ℚ[X]) = 2 + X := by
  have hx : reflect 1 (X : ℚ[X]) = 1 := by
    simpa only [pow_one, revAt_le le_rfl, Nat.sub_self, pow_zero] using
      (reflect_monomial (R := ℚ) 1 1)
  rw [reflect_add, ← C_1, reflect_C, ← C_ofNat, reflect_C_mul, hx]
  simp [add_comm]

private theorem reversed_linearTerm (m k ℓ : ℕ) :
    reflect (m + k + ℓ) ((1 + X : ℚ[X]) ^ m * (1 + 2 * X) ^ k * X ^ ℓ) =
      (1 + X) ^ m * (2 + X) ^ k := by
  have h₁ : (1 + X : ℚ[X]).natDegree ≤ 1 := (natDegree_add_le _ _).trans (by simp)
  have h₂ : (1 + 2 * X : ℚ[X]).natDegree ≤ 1 :=
    (natDegree_add_le _ _).trans (by
      apply max_le (by simp)
      exact natDegree_mul_le.trans (by simp))
  have hdeg₁ := natDegree_pow_le.trans
    (show m * (1 + X : ℚ[X]).natDegree ≤ m by simpa using Nat.mul_le_mul_left m h₁)
  have hdeg₂ := natDegree_pow_le.trans
    (show k * (1 + 2 * X : ℚ[X]).natDegree ≤ k by simpa using Nat.mul_le_mul_left k h₂)
  have hdeg := (natDegree_mul_le (p := (1 + X : ℚ[X]) ^ m) (q := (1 + 2 * X) ^ k)).trans
    (Nat.add_le_add hdeg₁ hdeg₂)
  rw [reflect_mul _ _ hdeg (by simp : ((X : ℚ[X]) ^ ℓ).natDegree ≤ ℓ),
    reflect_mul _ _ hdeg₁ hdeg₂, reflect_power_linear _ h₁, reflect_power_linear _ h₂,
    reflect_one_add_X, reflect_one_add_twoX, reflect_monomial, revAt_le le_rfl, Nat.sub_self, pow_zero, mul_one]

/-- The coefficient reversal for the trace of x^ℓ p^m. -/
theorem linearPhaseEndpoint_coeff_eq_sum (m ℓ : ℕ) :
    (linearPhaseEndpoint m).coeff (m + ℓ) =
      ∑ k ∈ Finset.range (m + 1), negativeHalfBinomial k *
        (((1 + X : ℚ[X]) ^ m * (1 + 2 * X) ^ k * X ^ ℓ).coeff (2 * k)) := by
  rw [linearPhaseEndpoint, halfBinomialTruncation, sum_comp, Finset.mul_sum, finsetSum_coeff]
  refine Finset.sum_congr rfl fun k hk ↦ ?_
  have hkm : k ≤ m := by have := Finset.mem_range.mp hk; omega
  rw [← C_mul_X_pow_eq_monomial, mul_comp, C_comp, pow_comp, X_comp, mul_pow]
  have hmul : (1 + X : ℚ[X]) ^ m * (C (negativeHalfBinomial k) * (X ^ k * (2 + X) ^ k)) =
      C (negativeHalfBinomial k) * (X ^ k * ((1 + X) ^ m * (2 + X) ^ k)) := by ring
  rw [hmul, coeff_C_mul, coeff_X_pow_mul', if_pos (by omega), ← reversed_linearTerm m k ℓ,
    coeff_reflect, revAt_le (by omega)]
  congr 2
  omega

end MathieuQuadratic
