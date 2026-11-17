/-
Copyright (c) 2026 Jacobian Research contributors. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Jacobian Research contributors
-/

import Mathlib.Algebra.Group.ForwardDiff
import Mathlib.Algebra.Polynomial.Coeff
import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.Algebra.Polynomial.Div
import Mathlib.Algebra.Polynomial.Laurent
import Mathlib.Algebra.MvPolynomial.Eval
import Mathlib.Data.Nat.Factorial.DoubleFactorial
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.Ring

/-!
# The residual SIC2C4 finite sum

This module formalizes the discrete part of the direct coefficient-extraction
proof of SIC2C4: finite-difference vanishing, the general remainder principle
and its rank-one endpoint-residue specialization for alternating quotient
sums, the specialized product polynomials, the first-order beta recurrence,
its factorial/double-factorial evaluation, the triangular recurrence for
arbitrary repeated poles, the finite
remainder-to-jet identity, and the two normalized displayed binomial sums.
It also formalizes the coefficient functional and beta identity, the monomial
and balanced-array coefficient-extraction formulas, the scalar chart and
polynomial coefficient identities, and the formal integrals of the resulting
pure and mixed chart polynomials.  The literal Laurent chart witness, its
all-order binomial expansion, and the pure and mixed constant-term
identifications are formalized as well.  Finally, the module defines the
original four-variable `F,Q` as multivariate polynomials and proves that
their displayed chart substitution is exactly that Laurent witness.
-/

open scoped BigOperators
open scoped LaurentPolynomial
open Finset

namespace FiniteEtaleKeller.SIC2C4

/-- Algebraic chart identity underlying equation (4.3). -/
theorem chartIdentity {K : Type*} [Field K] [CharZero K]
    (x v : K) (hx : x ≠ 0) :
    (1 + x) * ((1 - v ^ 2) / (2 * x) - (2 + x) * v ^ 2 / 2) =
      (1 + x) / (2 * x) * (1 - v ^ 2 * (1 + x) ^ 2) := by
  field_simp [hx]
  ring

/-- The coefficient functional `I(t^j)=1/(j+1)` used as a purely formal
integral in the chart proof. -/
noncomputable def formalIntegral (P : Polynomial ℚ) : ℚ :=
  P.sum fun n a ↦ a / (n + 1 : ℚ)

lemma formalIntegral_add (P Q : Polynomial ℚ) :
    formalIntegral (P + Q) = formalIntegral P + formalIntegral Q := by
  unfold formalIntegral
  rw [Polynomial.sum_add_index]
  · intro n
    simp
  · intro n a b
    ring

lemma formalIntegral_monomial (n : ℕ) (a : ℚ) :
    formalIntegral (Polynomial.monomial n a) = a / (n + 1 : ℚ) := by
  simp [formalIntegral]

lemma formalIntegral_smul (c : ℚ) (P : Polynomial ℚ) :
    formalIntegral (c • P) = c * formalIntegral P := by
  induction P using Polynomial.induction_on' with
  | add P Q inductionP inductionQ =>
      rw [smul_add, formalIntegral_add, formalIntegral_add,
        inductionP, inductionQ]
      ring
  | monomial n a =>
      rw [Polynomial.smul_eq_C_mul, Polynomial.C_mul_monomial,
        formalIntegral_monomial, formalIntegral_monomial]
      ring

lemma formalIntegral_sub (P Q : Polynomial ℚ) :
    formalIntegral (P - Q) = formalIntegral P - formalIntegral Q := by
  rw [sub_eq_add_neg, show -Q = (-1 : ℚ) • Q by simp,
    formalIntegral_add, formalIntegral_smul]
  ring

/-- Formal fundamental theorem for the coefficient functional. -/
theorem formalIntegral_derivative (P : Polynomial ℚ) :
    formalIntegral P.derivative = P.eval 1 - P.eval 0 := by
  induction P using Polynomial.induction_on' with
  | add P Q inductionP inductionQ =>
      rw [Polynomial.derivative_add, formalIntegral_add, inductionP, inductionQ,
        Polynomial.eval_add, Polynomial.eval_add]
      ring
  | monomial n a =>
      cases n with
      | zero => simp [formalIntegral]
      | succ n =>
          rw [Polynomial.derivative_monomial_succ, formalIntegral_monomial]
          simp only [Polynomial.eval_monomial, one_pow, mul_one,
            zero_pow (Nat.succ_ne_zero _), mul_zero, sub_zero]
          have hnonzero : (n + 1 : ℚ) ≠ 0 := by positivity
          field_simp [hnonzero]

/-- The formal even moment used when converting the chart constant term to
the denominator `2k+1`. -/
theorem formalIntegral_two_mul_X_sub_one_even (k : ℕ) :
    formalIntegral
        ((Polynomial.C 2 * Polynomial.X - Polynomial.C 1) ^ (2 * k)) =
      1 / (2 * k + 1 : ℚ) := by
  let A : Polynomial ℚ :=
    Polynomial.C (1 / (2 * (2 * k + 1 : ℚ))) *
      (Polynomial.C 2 * Polynomial.X - Polynomial.C 1) ^ (2 * k + 1)
  have hderivative :
      A.derivative =
        (Polynomial.C 2 * Polynomial.X - Polynomial.C 1) ^ (2 * k) := by
    dsimp [A]
    rw [Polynomial.derivative_mul, Polynomial.derivative_C,
      zero_mul, zero_add, Polynomial.derivative_pow,
      Polynomial.derivative_sub, Polynomial.derivative_mul,
      Polynomial.derivative_C, Polynomial.derivative_X,
      Polynomial.derivative_C]
    simp only [zero_mul, zero_add, mul_one, Polynomial.C_1,
      sub_zero]
    rw [show 2 * k + 1 - 1 = 2 * k by omega]
    have hnonzero : (2 * k + 1 : ℚ) ≠ 0 := by positivity
    have hscalar :
        (1 / (2 * (2 * k + 1 : ℚ))) * (2 * k + 1 : ℚ) * 2 = 1 := by
      field_simp [hnonzero]
    push_cast
    have hC :
        Polynomial.C (1 / (2 * (2 * k + 1 : ℚ))) *
            Polynomial.C (2 * k + 1 : ℚ) * Polynomial.C 2 =
          1 := by
      rw [← Polynomial.C_mul, ← Polynomial.C_mul, hscalar,
        Polynomial.C_1]
    calc
      Polynomial.C (1 / (2 * (2 * k + 1 : ℚ))) *
          (Polynomial.C (2 * k + 1 : ℚ) *
            (Polynomial.C 2 * Polynomial.X - 1) ^ (2 * k) *
              Polynomial.C 2)
        = (Polynomial.C (1 / (2 * (2 * k + 1 : ℚ))) *
              Polynomial.C (2 * k + 1 : ℚ) * Polynomial.C 2) *
            (Polynomial.C 2 * Polynomial.X - 1) ^ (2 * k) := by ring
      _ = (Polynomial.C 2 * Polynomial.X - 1) ^ (2 * k) := by
        rw [hC, one_mul]
  rw [← hderivative, formalIntegral_derivative]
  dsimp [A]
  simp only [Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_pow,
    Polynomial.eval_sub, Polynomial.eval_X]
  have hnonzero : (2 * k + 1 : ℚ) ≠ 0 := by positivity
  field_simp [hnonzero]
  norm_num [pow_succ, pow_mul]

/-- The algebraic beta functional appearing in the monomial chart. -/
noncomputable def formalBeta (a b : ℕ) : ℚ :=
  formalIntegral
    (Polynomial.X ^ a * (1 - Polynomial.X) ^ b)

lemma formalBeta_zero (b : ℕ) :
    formalBeta 0 b = 1 / (b + 1 : ℚ) := by
  let A : Polynomial ℚ :=
    Polynomial.C (-1 / (b + 1 : ℚ)) *
      (1 - Polynomial.X) ^ (b + 1)
  have hderivative :
      A.derivative = (1 - Polynomial.X) ^ b := by
    dsimp [A]
    rw [Polynomial.derivative_mul, Polynomial.derivative_C,
      zero_mul, zero_add, Polynomial.derivative_pow,
      Polynomial.derivative_sub, Polynomial.derivative_one,
      Polynomial.derivative_X, zero_sub]
    rw [show b + 1 - 1 = b by omega]
    have hnonzero : (b + 1 : ℚ) ≠ 0 := by positivity
    have hscalar :
        (-1 / (b + 1 : ℚ)) * (b + 1 : ℚ) * (-1) = 1 := by
      field_simp [hnonzero]
    push_cast
    have hC :
        Polynomial.C (-1 / (b + 1 : ℚ)) *
            Polynomial.C (b + 1 : ℚ) * (-1) =
          1 := by
      rw [show (-1 : Polynomial ℚ) = Polynomial.C (-1) by norm_num,
        ← Polynomial.C_mul, ← Polynomial.C_mul, hscalar, Polynomial.C_1]
    calc
      Polynomial.C (-1 / (b + 1 : ℚ)) *
          (Polynomial.C (b + 1 : ℚ) * (1 - Polynomial.X) ^ b * (-1))
        = (Polynomial.C (-1 / (b + 1 : ℚ)) *
              Polynomial.C (b + 1 : ℚ) * (-1)) *
            (1 - Polynomial.X) ^ b := by ring
      _ = (1 - Polynomial.X) ^ b := by rw [hC, one_mul]
  rw [formalBeta, pow_zero, one_mul, ← hderivative,
    formalIntegral_derivative]
  dsimp [A]
  simp only [Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_pow,
    Polynomial.eval_sub, Polynomial.eval_one, Polynomial.eval_X]
  have hnonzero : (b + 1 : ℚ) ≠ 0 := by positivity
  field_simp [hnonzero]
  ring

lemma formalBeta_succ_recurrence (a b : ℕ) :
    (b + 1 : ℚ) * formalBeta (a + 1) b =
      (a + 1 : ℚ) * formalBeta a (b + 1) := by
  let H : Polynomial ℚ :=
    Polynomial.X ^ (a + 1) * (1 - Polynomial.X) ^ (b + 1)
  have hderivative :
      H.derivative =
        Polynomial.C (a + 1 : ℚ) *
            (Polynomial.X ^ a * (1 - Polynomial.X) ^ (b + 1)) -
          Polynomial.C (b + 1 : ℚ) *
            (Polynomial.X ^ (a + 1) * (1 - Polynomial.X) ^ b) := by
    dsimp [H]
    rw [Polynomial.derivative_mul, Polynomial.derivative_pow,
      Polynomial.derivative_X, mul_one, Polynomial.derivative_pow,
      Polynomial.derivative_sub, Polynomial.derivative_one,
      Polynomial.derivative_X, zero_sub]
    push_cast
    ring
  have hftc := formalIntegral_derivative H
  have hboundary : H.eval 1 - H.eval 0 = 0 := by
    dsimp [H]
    simp
  rw [hboundary, hderivative, formalIntegral_sub,
    ← Polynomial.smul_eq_C_mul, ← Polynomial.smul_eq_C_mul,
    formalIntegral_smul, formalIntegral_smul] at hftc
  simp only [formalBeta] at hftc ⊢
  change
    (a + 1 : ℚ) *
          formalIntegral (Polynomial.X ^ a * (1 - Polynomial.X) ^ (b + 1)) -
        (b + 1 : ℚ) *
          formalIntegral (Polynomial.X ^ (a + 1) * (1 - Polynomial.X) ^ b) =
      0 at hftc
  linarith

/-- Entirely algebraic beta identity (4.1). -/
theorem formalBeta_eq_factorial_ratio (a b : ℕ) :
    formalBeta a b =
      ((a.factorial : ℚ) * (b.factorial : ℚ)) /
        ((a + b + 1).factorial : ℚ) := by
  induction a generalizing b with
  | zero =>
      rw [formalBeta_zero]
      simp only [Nat.factorial_zero, Nat.cast_one, one_mul, zero_add]
      rw [Nat.factorial_succ]
      have hfactorial : (b.factorial : ℚ) ≠ 0 := by positivity
      have hnext : (b + 1 : ℚ) ≠ 0 := by positivity
      field_simp [hfactorial, hnext]
      norm_num
  | succ a inductionHypothesis =>
      have hrec := formalBeta_succ_recurrence a b
      rw [inductionHypothesis (b + 1)] at hrec
      have hb : (b + 1 : ℚ) ≠ 0 := by positivity
      apply (mul_left_cancel₀ hb)
      calc
        (b + 1 : ℚ) * formalBeta (a + 1) b
          = (a + 1 : ℚ) *
              (((a.factorial : ℚ) * ((b + 1).factorial : ℚ)) /
                ((a + (b + 1) + 1).factorial : ℚ)) := hrec
        _ = (b + 1 : ℚ) *
              (((a + 1).factorial : ℚ) * (b.factorial : ℚ) /
                ((a + 1 + b + 1).factorial : ℚ)) := by
              have hindex :
                  a + (b + 1) + 1 = a + 1 + b + 1 := by omega
              rw [Nat.factorial_succ a, Nat.factorial_succ b]
              rw [hindex]
              push_cast
              ring

/-- The two polynomial coefficients selected by the pure and mixed chart
constant terms. -/
theorem chartCoefficientIdentities (m k : ℕ) :
    (((1 + Polynomial.X) ^ (m + 2 * k) : Polynomial ℚ).coeff m =
        ((m + 2 * k).choose m : ℚ)) ∧
      (((1 + Polynomial.X) ^ (m + 2 * k) : Polynomial ℚ).coeff (m - 1) =
        ((m + 2 * k).choose (m - 1) : ℚ)) := by
  constructor <;> apply Polynomial.coeff_one_add_X_pow

/-- Multiplication by `x⁻ᵐ` turns the Laurent constant term into the
ordinary `x^m` coefficient. -/
theorem laurentConstantTerm_shift {R : Type*} [CommRing R]
    (P : Polynomial R) (m : ℕ) :
    ((LaurentPolynomial.T (-(m : ℤ)) : LaurentPolynomial R) *
        Polynomial.toLaurent P).coeff 0 =
      P.coeff m := by
  change
    (AddMonoidAlgebra.single (-(m : ℤ)) 1 *
      Polynomial.toLaurent P).coeff 0 = P.coeff m
  rw [AddMonoidAlgebra.coeff_single_mul_apply]
  simp only [neg_neg, add_zero, one_mul, LaurentPolynomial.coeff_toLaurent]
  change
    (Finsupp.mapDomain (⇑Nat.castEmbedding) P.toFinsupp.coeff)
        (Nat.castEmbedding m) = P.coeff m
  rw [Finsupp.mapDomain_apply Nat.castEmbedding.injective,
    Polynomial.toFinsupp_apply]

/-- Coefficient-valued version of `laurentConstantTerm_shift`. -/
theorem laurentConstantTerm_C_mul_shift {R : Type*} [CommRing R]
    (a : R) (P : Polynomial R) (m : ℕ) :
    (LaurentPolynomial.C a * LaurentPolynomial.T (-(m : ℤ)) *
        Polynomial.toLaurent P).coeff 0 =
      a * P.coeff m := by
  rw [mul_assoc]
  change
    (AddMonoidAlgebra.single 0 a *
      (LaurentPolynomial.T (-(m : ℤ)) *
        Polynomial.toLaurent P)).coeff 0 =
      a * P.coeff m
  rw [AddMonoidAlgebra.coeff_single_mul_apply]
  simp only [neg_zero, zero_add]
  rw [laurentConstantTerm_shift]

/-- The coefficient in `v` multiplying
`x⁻ᵐ(1+x)^(m+2k)` in the termwise Laurent expansion of `F^m`. -/
noncomputable def chartLaurentWeight (m k : ℕ) : Polynomial ℚ :=
  Polynomial.monomial (2 * k)
    ((1 / (2 : ℚ) ^ m) * (-1 : ℚ) ^ k * (m.choose k : ℚ))

lemma formalIntegral_sum {ι : Type*} (s : Finset ι)
    (f : ι → Polynomial ℚ) :
    formalIntegral (∑ i ∈ s, f i) = ∑ i ∈ s, formalIntegral (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [formalIntegral]
  | @insert a s ha inductionHypothesis =>
      rw [sum_insert ha, sum_insert ha, formalIntegral_add,
        inductionHypothesis]

/-- Pure chart constant term, viewed as a polynomial in `v`. -/
noncomputable def pureChartConstantTerm (m : ℕ) : Polynomial ℚ :=
  ∑ k ∈ range (m + 1),
    Polynomial.monomial (2 * k)
      ((1 / (2 : ℚ) ^ m) * (-1 : ℚ) ^ k * (m.choose k : ℚ) *
        ((m + 2 * k).choose m : ℚ))

/-- Mixed chart constant term, viewed as a polynomial in `v`. -/
noncomputable def mixedChartConstantTerm (m : ℕ) : Polynomial ℚ :=
  ∑ k ∈ range (m + 1),
    Polynomial.monomial (2 * k)
      ((1 / (2 : ℚ) ^ m) * (-1 : ℚ) ^ k * (m.choose k : ℚ) *
        ((m + 2 * k).choose (m - 1) : ℚ))

/-- The literal termwise Laurent expansion whose constant term is the pure
chart polynomial.  Its coefficient ring is `ℚ[v]` and its Laurent variable
is `x`. -/
noncomputable def pureChartLaurentExpansion (m : ℕ) :
    LaurentPolynomial (Polynomial ℚ) :=
  ∑ k ∈ range (m + 1),
    LaurentPolynomial.C (chartLaurentWeight m k) *
      LaurentPolynomial.T (-(m : ℤ)) *
        Polynomial.toLaurent
          ((1 + Polynomial.X) ^ (m + 2 * k) :
            Polynomial (Polynomial ℚ))

/-- The corresponding termwise Laurent expansion after multiplication by
the chart coordinate `x`. -/
noncomputable def mixedChartLaurentExpansion (m : ℕ) :
    LaurentPolynomial (Polynomial ℚ) :=
  ∑ k ∈ range (m + 1),
    LaurentPolynomial.C (chartLaurentWeight m k) *
      LaurentPolynomial.T (-((m - 1 : ℕ) : ℤ)) *
        Polynomial.toLaurent
          ((1 + Polynomial.X) ^ (m + 2 * k) :
            Polynomial (Polynomial ℚ))

/-- Binomial expansion in the form used by the Laurent chart witness. -/
lemma one_sub_pow_binomial {R : Type*} [CommRing R] (Y : R) (m : ℕ) :
    (1 - Y) ^ m =
      ∑ k ∈ range (m + 1),
        (-1 : R) ^ k * (m.choose k : R) * Y ^ k := by
  calc
    (1 - Y) ^ m = (-Y + 1) ^ m := by ring
    _ = ∑ k ∈ range (m + 1),
        (-Y) ^ k * 1 ^ (m - k) * (m.choose k : R) :=
          add_pow (-Y) 1 m
    _ = ∑ k ∈ range (m + 1),
        (-1 : R) ^ k * (m.choose k : R) * Y ^ k := by
      apply sum_congr rfl
      intro k _
      rw [neg_pow, one_pow]
      ring

/-- Equation (4.3) represented literally in `ℚ[v][x,x⁻¹]`. -/
noncomputable def chartLaurentWitness :
    LaurentPolynomial (Polynomial ℚ) :=
  LaurentPolynomial.C (Polynomial.C (1 / 2 : ℚ)) *
    LaurentPolynomial.T (-1) *
      (1 + LaurentPolynomial.T 1) *
        (1 - LaurentPolynomial.C (Polynomial.X ^ 2) *
          (1 + LaurentPolynomial.T 1) ^ 2)

/-- The coefficient type for the literal four-variable SIC2C4 witness.  The
coordinates are ordered as `ξ₁, ξ₂, z₁, z₂`. -/
abbrev WitnessPolynomial := MvPolynomial (Fin 4) ℚ

noncomputable def witnessXiOne : WitnessPolynomial := MvPolynomial.X 0
noncomputable def witnessXiTwo : WitnessPolynomial := MvPolynomial.X 1
noncomputable def witnessZOne : WitnessPolynomial := MvPolynomial.X 2
noncomputable def witnessZTwo : WitnessPolynomial := MvPolynomial.X 3

/-- The four displayed quadratic coordinates `R,Z,W,T`. -/
noncomputable def witnessR : WitnessPolynomial :=
  witnessXiOne * witnessZOne + witnessXiTwo * witnessZTwo

noncomputable def witnessZ : WitnessPolynomial :=
  witnessXiOne * witnessZTwo

noncomputable def witnessW : WitnessPolynomial :=
  2 * witnessXiTwo * witnessZOne

noncomputable def witnessT : WitnessPolynomial :=
  witnessXiOne * witnessZOne - witnessXiTwo * witnessZTwo

/-- The original four-variable polynomials `F,Q` from Theorem 1.1. -/
noncomputable def witnessF : WitnessPolynomial :=
  (witnessR + witnessZ) *
    (witnessR ^ 2 * witnessW -
      MvPolynomial.C (1 / 2 : ℚ) *
        (2 * witnessR + witnessZ) * witnessT ^ 2)

noncomputable def witnessQ : WitnessPolynomial := witnessZ

/-- The direct substitution
`(ξ₁,ξ₂,z₁,z₂)=(1,(1-v)/(2x),(1+v)/2,x)` into
`ℚ[v][x,x⁻¹]`. -/
noncomputable def chartSubstitutionPoint :
    Fin 4 → LaurentPolynomial (Polynomial ℚ)
  | 0 => 1
  | 1 =>
      LaurentPolynomial.C
          (Polynomial.C (1 / 2 : ℚ) * (1 - Polynomial.X)) *
        LaurentPolynomial.T (-1)
  | 2 =>
      LaurentPolynomial.C
        (Polynomial.C (1 / 2 : ℚ) * (1 + Polynomial.X))
  | 3 => LaurentPolynomial.T 1

noncomputable def chartSubstitution :
    WitnessPolynomial →+* LaurentPolynomial (Polynomial ℚ) :=
  MvPolynomial.eval₂Hom
    (algebraMap ℚ (LaurentPolynomial (Polynomial ℚ)))
    chartSubstitutionPoint

lemma chartT_inv_mul_chartT :
    (LaurentPolynomial.T (-1) :
        LaurentPolynomial (Polynomial ℚ)) * LaurentPolynomial.T 1 = 1 := by
  rw [← LaurentPolynomial.T_add]
  norm_num

lemma chartHalf_mul_two :
    LaurentPolynomial.C (Polynomial.C (1 / 2 : ℚ)) *
        (2 : LaurentPolynomial (Polynomial ℚ)) =
      1 := by
  change
    LaurentPolynomial.C (Polynomial.C (1 / 2 : ℚ)) *
        LaurentPolynomial.C (Polynomial.C 2) =
      1
  rw [← map_mul LaurentPolynomial.C, ← map_mul Polynomial.C]
  norm_num

lemma chartSubstitution_witnessR :
    chartSubstitution witnessR = 1 := by
  simp only [chartSubstitution, witnessR, witnessXiOne, witnessXiTwo,
    witnessZOne, witnessZTwo, map_add, map_mul,
    MvPolynomial.eval₂Hom_X']
  simp only [chartSubstitutionPoint, one_div, map_mul, map_add, map_one,
    one_mul, map_sub, Int.reduceNeg, LaurentPolynomial.mul_T_assoc,
    neg_add_cancel, LaurentPolynomial.T_zero, mul_one]
  calc
    _ = LaurentPolynomial.C (Polynomial.C (1 / 2 : ℚ)) * 2 := by ring_nf
    _ = 1 := chartHalf_mul_two

lemma chartSubstitution_witnessZ :
    chartSubstitution witnessZ = LaurentPolynomial.T 1 := by
  simp [chartSubstitution, witnessZ, witnessXiOne, witnessZTwo,
    chartSubstitutionPoint]

lemma chartSubstitution_witnessT :
    chartSubstitution witnessT =
      LaurentPolynomial.C Polynomial.X := by
  simp only [chartSubstitution, witnessT, witnessXiOne, witnessXiTwo,
    witnessZOne, witnessZTwo, map_sub, map_mul,
    MvPolynomial.eval₂Hom_X']
  simp only [chartSubstitutionPoint, one_div, map_mul, map_add, map_one,
    one_mul, map_sub, Int.reduceNeg, LaurentPolynomial.mul_T_assoc,
    neg_add_cancel, LaurentPolynomial.T_zero, mul_one]
  calc
    _ = (LaurentPolynomial.C (Polynomial.C (1 / 2 : ℚ)) * 2) *
        LaurentPolynomial.C Polynomial.X := by ring_nf
    _ = LaurentPolynomial.C Polynomial.X := by
      rw [chartHalf_mul_two, one_mul]

lemma chartSubstitution_witnessW :
    chartSubstitution witnessW =
      LaurentPolynomial.C
          (Polynomial.C (1 / 2 : ℚ) * (1 - Polynomial.X ^ 2)) *
        LaurentPolynomial.T (-1) := by
  simp only [chartSubstitution, witnessW, witnessXiTwo, witnessZOne,
    map_mul, map_ofNat, MvPolynomial.eval₂Hom_X']
  simp only [chartSubstitutionPoint, one_div, map_mul, map_sub, map_one,
    Int.reduceNeg, map_add, map_pow]
  have hsquare :
      2 * LaurentPolynomial.C (Polynomial.C (1 / 2 : ℚ)) ^ 2 =
        LaurentPolynomial.C (Polynomial.C (1 / 2 : ℚ)) := by
    calc
      _ = (LaurentPolynomial.C (Polynomial.C (1 / 2 : ℚ)) * 2) *
          LaurentPolynomial.C (Polynomial.C (1 / 2 : ℚ)) := by ring
      _ = _ := by rw [chartHalf_mul_two, one_mul]
  norm_num at hsquare ⊢
  linear_combination
    (LaurentPolynomial.T (-1) *
      (1 - LaurentPolynomial.C Polynomial.X ^ 2)) * hsquare

/-- Ring-theoretic form of (4.3), requiring only a specified inverse for
the Laurent coordinate. -/
lemma laurentChartIdentity_of_inverse
    (x xinv v : LaurentPolynomial (Polynomial ℚ))
    (h : xinv * x = 1) :
    (1 + x) *
        (algebraMap ℚ (LaurentPolynomial (Polynomial ℚ)) (1 / 2) *
            (1 - v ^ 2) * xinv -
          algebraMap ℚ (LaurentPolynomial (Polynomial ℚ)) (1 / 2) *
            (2 + x) * v ^ 2) =
      algebraMap ℚ (LaurentPolynomial (Polynomial ℚ)) (1 / 2) *
        (1 + x) * xinv * (1 - v ^ 2 * (1 + x) ^ 2) := by
  calc
    _ = algebraMap ℚ (LaurentPolynomial (Polynomial ℚ)) (1 / 2) *
          (1 + x) * xinv * (1 - v ^ 2 * (1 + x) ^ 2) +
        algebraMap ℚ (LaurentPolynomial (Polynomial ℚ)) (1 / 2) *
          v ^ 2 * (x + 1) * (x + 2) * (xinv * x - 1) := by ring
    _ = _ := by rw [h, sub_self, mul_zero, add_zero]

/-- The displayed four-variable `F` specializes exactly to the literal
Laurent witness used in the all-order proof. -/
theorem chartSubstitution_witnessF :
    chartSubstitution witnessF = chartLaurentWitness := by
  simp only [witnessF, map_mul, map_sub, map_pow, map_add, map_ofNat]
  rw [show chartSubstitution (MvPolynomial.C (1 / 2 : ℚ)) =
      algebraMap ℚ (LaurentPolynomial (Polynomial ℚ)) (1 / 2 : ℚ) by
    simp [chartSubstitution, MvPolynomial.eval₂Hom_C]]
  rw [chartSubstitution_witnessR, chartSubstitution_witnessZ,
    chartSubstitution_witnessT, chartSubstitution_witnessW]
  unfold chartLaurentWitness
  norm_num
  have h := laurentChartIdentity_of_inverse
    (LaurentPolynomial.T 1 :
      LaurentPolynomial (Polynomial ℚ))
    (LaurentPolynomial.T (-1) :
      LaurentPolynomial (Polynomial ℚ))
    (LaurentPolynomial.C Polynomial.X)
    chartT_inv_mul_chartT
  norm_num at h
  calc
    _ = LaurentPolynomial.C (Polynomial.C (1 / 2 : ℚ)) *
          (1 + LaurentPolynomial.T 1) * LaurentPolynomial.T (-1) *
            (1 - LaurentPolynomial.C Polynomial.X ^ 2 *
              (1 + LaurentPolynomial.T 1) ^ 2) := h
    _ = _ := by ring

theorem chartSubstitution_witnessQ :
    chartSubstitution witnessQ = LaurentPolynomial.T 1 := by
  rw [witnessQ, chartSubstitution_witnessZ]

/-- The substitution also intertwines all pure and mixed powers. -/
theorem chartSubstitution_witnessF_pow (m : ℕ) :
    chartSubstitution (witnessF ^ m) = chartLaurentWitness ^ m := by
  rw [map_pow, chartSubstitution_witnessF]

theorem chartSubstitution_witnessQ_mul_witnessF_pow (m : ℕ) :
    chartSubstitution (witnessQ * witnessF ^ m) =
      LaurentPolynomial.T 1 * chartLaurentWitness ^ m := by
  rw [map_mul, map_pow, chartSubstitution_witnessQ,
    chartSubstitution_witnessF]

/-- The all-order binomial expansion of the literal Laurent chart witness. -/
theorem chartLaurentWitness_pow_eq_pureExpansion (m : ℕ) :
    chartLaurentWitness ^ m = pureChartLaurentExpansion m := by
  unfold chartLaurentWitness pureChartLaurentExpansion chartLaurentWeight
  rw [mul_pow, mul_pow, mul_pow, one_sub_pow_binomial, mul_sum]
  apply sum_congr rfl
  intro k _
  have hhalfPoly :
      (Polynomial.C (1 / 2 : ℚ)) ^ m =
        Polynomial.C ((2 : ℚ)⁻¹ ^ m) := by
    rw [← map_pow]
    norm_num
  have hhalfLaurent :
      LaurentPolynomial.C (Polynomial.C (1 / 2 : ℚ)) ^ m =
        LaurentPolynomial.C (Polynomial.C ((2 : ℚ)⁻¹ ^ m)) := by
    rw [← map_pow, hhalfPoly]
  have hpowValue :
      LaurentPolynomial.C (Polynomial.C ((2 : ℚ)⁻¹ ^ m)) =
        LaurentPolynomial.C (Polynomial.C (((2 : ℚ) ^ m)⁻¹)) := by
    congr 2
    rw [inv_pow]
  simp only [map_pow, LaurentPolynomial.T_pow, mul_pow, map_add, map_one,
    Polynomial.toLaurent_X]
  rw [← pow_mul, ← pow_mul, ← Polynomial.C_mul_X_pow_eq_monomial]
  simp only [map_mul, map_pow]
  rw [show
    (1 + LaurentPolynomial.T (1 : ℤ)) ^ (m + 2 * k) =
      (1 + LaurentPolynomial.T (1 : ℤ)) ^ m *
        (1 + LaurentPolynomial.T (1 : ℤ)) ^ (2 * k) by
          rw [pow_add]]
  rw [show (m : ℤ) * -1 = -(m : ℤ) by omega]
  let Z : LaurentPolynomial (Polynomial ℚ) :=
    LaurentPolynomial.T (-(m : ℤ)) *
      (m.choose k : LaurentPolynomial (Polynomial ℚ)) *
        LaurentPolynomial.C Polynomial.X ^ (k * 2) * (-1) ^ k *
          (1 + LaurentPolynomial.T 1) ^ m *
            (1 + LaurentPolynomial.T 1) ^ (k * 2)
  calc
    _ = LaurentPolynomial.C (Polynomial.C (1 / 2 : ℚ)) ^ m * Z := by
      dsimp [Z]
      ring
    _ = LaurentPolynomial.C (Polynomial.C ((2 : ℚ)⁻¹ ^ m)) * Z := by
      rw [hhalfLaurent]
    _ = _ := by
      dsimp [Z]
      norm_num
      rw [hhalfLaurent, hpowValue]
      ring

/-- Multiplication by the chart coordinate `x` shifts the pure Laurent
expansion to the mixed one. -/
theorem chartX_mul_pureExpansion_eq_mixedExpansion
    (m : ℕ) (hm : 0 < m) :
    LaurentPolynomial.T 1 * pureChartLaurentExpansion m =
      mixedChartLaurentExpansion m := by
  unfold pureChartLaurentExpansion mixedChartLaurentExpansion
  rw [mul_sum]
  apply sum_congr rfl
  intro k _
  have hexponent :
      (1 : ℤ) + -(m : ℤ) = -((m - 1 : ℕ) : ℤ) := by omega
  calc
    LaurentPolynomial.T 1 *
        (LaurentPolynomial.C (chartLaurentWeight m k) *
          LaurentPolynomial.T (-(m : ℤ)) *
            Polynomial.toLaurent
              ((1 + Polynomial.X) ^ (m + 2 * k) :
                Polynomial (Polynomial ℚ)))
      = LaurentPolynomial.C (chartLaurentWeight m k) *
          (LaurentPolynomial.T 1 * LaurentPolynomial.T (-(m : ℤ))) *
            Polynomial.toLaurent
              ((1 + Polynomial.X) ^ (m + 2 * k) :
                Polynomial (Polynomial ℚ)) := by ring
    _ = LaurentPolynomial.C (chartLaurentWeight m k) *
          LaurentPolynomial.T (-((m - 1 : ℕ) : ℤ)) *
            Polynomial.toLaurent
              ((1 + Polynomial.X) ^ (m + 2 * k) :
                Polynomial (Polynomial ℚ)) := by
      rw [← LaurentPolynomial.T_add, hexponent]

/-- The literal mixed chart witness `x F^m` is the mixed termwise
expansion. -/
theorem chartX_mul_chartLaurentWitness_pow_eq_mixedExpansion
    (m : ℕ) (hm : 0 < m) :
    LaurentPolynomial.T 1 * chartLaurentWitness ^ m =
      mixedChartLaurentExpansion m := by
  rw [chartLaurentWitness_pow_eq_pureExpansion,
    chartX_mul_pureExpansion_eq_mixedExpansion m hm]

/-- Constant-term extraction from the termwise Laurent expansion produces
exactly `pureChartConstantTerm`. -/
theorem pureChartLaurentExpansion_constantTerm (m : ℕ) :
    (pureChartLaurentExpansion m).coeff 0 =
      pureChartConstantTerm m := by
  unfold pureChartLaurentExpansion pureChartConstantTerm
  rw [AddMonoidAlgebra.coeff_sum]
  simp only [Finset.sum_apply']
  apply sum_congr rfl
  intro k _
  rw [laurentConstantTerm_C_mul_shift,
    Polynomial.coeff_one_add_X_pow]
  unfold chartLaurentWeight
  change
    Polynomial.monomial (2 * k)
        ((1 / (2 : ℚ) ^ m) * (-1 : ℚ) ^ k * (m.choose k : ℚ)) *
        Polynomial.C ((m + 2 * k).choose m : ℚ) =
      _
  rw [Polynomial.monomial_mul_C]

/-- Constant-term extraction from the termwise mixed Laurent expansion
produces exactly `mixedChartConstantTerm`. -/
theorem mixedChartLaurentExpansion_constantTerm (m : ℕ) :
    (mixedChartLaurentExpansion m).coeff 0 =
      mixedChartConstantTerm m := by
  unfold mixedChartLaurentExpansion mixedChartConstantTerm
  rw [AddMonoidAlgebra.coeff_sum]
  simp only [Finset.sum_apply']
  apply sum_congr rfl
  intro k _
  rw [laurentConstantTerm_C_mul_shift,
    Polynomial.coeff_one_add_X_pow]
  unfold chartLaurentWeight
  change
    Polynomial.monomial (2 * k)
        ((1 / (2 : ℚ) ^ m) * (-1 : ℚ) ^ k * (m.choose k : ℚ)) *
        Polynomial.C ((m + 2 * k).choose (m - 1) : ℚ) =
      _
  rw [Polynomial.monomial_mul_C]

/-- The constant term of the literal `m`-th power of the chart witness is
the pure chart polynomial. -/
theorem chartLaurentWitness_pow_constantTerm (m : ℕ) :
    (chartLaurentWitness ^ m).coeff 0 =
      pureChartConstantTerm m := by
  rw [chartLaurentWitness_pow_eq_pureExpansion,
    pureChartLaurentExpansion_constantTerm]

/-- The constant term of the literal mixed chart witness `x F^m` is the
mixed chart polynomial. -/
theorem chartX_mul_chartLaurentWitness_pow_constantTerm
    (m : ℕ) (hm : 0 < m) :
    (LaurentPolynomial.T 1 * chartLaurentWitness ^ m).coeff 0 =
      mixedChartConstantTerm m := by
  rw [chartX_mul_chartLaurentWitness_pow_eq_mixedExpansion m hm,
    mixedChartLaurentExpansion_constantTerm]

/-- Applying the formal integral to the pure chart constant term produces
the normalized finite sum in (4.4). -/
theorem formalIntegral_pureChartConstantTerm (m : ℕ) :
    formalIntegral (pureChartConstantTerm m) =
      (1 / (2 : ℚ) ^ m) *
        (∑ k ∈ range (m + 1),
          (-1 : ℚ) ^ k * (m.choose k : ℚ) *
            (((m + 2 * k).choose m : ℚ) / (2 * k + 1))) := by
  rw [pureChartConstantTerm, formalIntegral_sum, mul_sum]
  apply sum_congr rfl
  intro k _
  rw [formalIntegral_monomial]
  push_cast
  ring

/-- Applying the formal integral to the mixed chart constant term produces
the normalized finite sum in (4.5). -/
theorem formalIntegral_mixedChartConstantTerm (m : ℕ) :
    formalIntegral (mixedChartConstantTerm m) =
      (1 / (2 : ℚ) ^ m) *
        (∑ k ∈ range (m + 1),
          (-1 : ℚ) ^ k * (m.choose k : ℚ) *
            (((m + 2 * k).choose (m - 1) : ℚ) / (2 * k + 1))) := by
  rw [mixedChartConstantTerm, formalIntegral_sum, mul_sum]
  apply sum_congr rfl
  intro k _
  rw [formalIntegral_monomial]
  push_cast
  ring

/-- Contraction of a balanced bidegree monomial after matching its two indices. -/
def monomialContraction (n a b : ℕ) : ℚ :=
  if a = b then (a.factorial : ℚ) * ((n - a).factorial : ℚ) else 0

/-- Constant-term and formal-beta value of the same monomial. -/
noncomputable def monomialChartValue (n a b : ℕ) : ℚ :=
  if a = b then
    ((n + 1).factorial : ℚ) *
      formalBeta a (n - a)
  else 0

/-- Monomial case of the contraction-to-chart formula (4.2). -/
theorem monomialCoefficientExtraction (n a b : ℕ) (ha : a ≤ n) :
    monomialChartValue n a b = monomialContraction n a b := by
  unfold monomialChartValue monomialContraction
  split_ifs
  · rw [formalBeta_eq_factorial_ratio]
    have hsum : a + (n - a) + 1 = n + 1 := by omega
    rw [hsum]
    have hfactorial : ((n + 1).factorial : ℚ) ≠ 0 := by positivity
    field_simp [hfactorial]
  · rfl

/-- Contraction of an arbitrary balanced bidegree-`(n,n)` coefficient
array. -/
noncomputable def balancedContractionValue
    (n : ℕ) (c : ℕ → ℕ → ℚ) : ℚ :=
  ∑ a ∈ range (n + 1), ∑ b ∈ range (n + 1),
    c a b * monomialContraction n a b

/-- Constant-term/formal-beta value of the same balanced coefficient
array. -/
noncomputable def balancedChartValue
    (n : ℕ) (c : ℕ → ℕ → ℚ) : ℚ :=
  ∑ a ∈ range (n + 1), ∑ b ∈ range (n + 1),
    c a b * monomialChartValue n a b

/-- Linear assembly of the monomial contraction-to-chart formula (4.2). -/
theorem balancedCoefficientExtraction
    (n : ℕ) (c : ℕ → ℕ → ℚ) :
    balancedChartValue n c = balancedContractionValue n c := by
  unfold balancedChartValue balancedContractionValue
  apply sum_congr rfl
  intro a ha
  apply sum_congr rfl
  intro b _
  rw [monomialCoefficientExtraction n a b (by
    simp only [mem_range] at ha
    omega)]

/-- The finite-difference cancellation used for the pure SIC2C4 sum. -/
theorem finiteDifference_vanishes (P : Polynomial ℚ) (m : ℕ)
    (hP : P.natDegree < m) :
    (∑ k ∈ range (m + 1),
      ((-1 : ℤ) ^ (m - k) * m.choose k) • P.eval (k : ℚ)) = 0 := by
  have hdiff := congr_fun (Polynomial.fwdDiff_iter_eq_zero_of_degree_lt hP) 0
  rw [Pi.zero_apply, fwdDiff_iter_eq_sum_shift] at hdiff
  simpa using hdiff

lemma negOnePow_sub_eq_mul (m k : ℕ) (hk : k ≤ m) :
    (-1 : ℚ) ^ (m - k) = (-1 : ℚ) ^ m * (-1 : ℚ) ^ k := by
  conv_rhs => lhs; rw [show m = (m - k) + k by omega, pow_add]
  rw [mul_assoc, ← mul_pow]
  norm_num

/-- The sign convention used in the displayed SIC2C4 sums. -/
theorem alternatingPolynomialSum_vanishes (P : Polynomial ℚ) (m : ℕ)
    (hP : P.natDegree < m) :
    (∑ k ∈ range (m + 1),
      (-1 : ℚ) ^ k * (m.choose k : ℚ) * P.eval (k : ℚ)) = 0 := by
  have hdiff := finiteDifference_vanishes P m hP
  have hq :
      (∑ k ∈ range (m + 1),
        (-1 : ℚ) ^ (m - k) * (m.choose k : ℚ) * P.eval (k : ℚ)) = 0 := by
    simpa [zsmul_eq_mul] using hdiff
  calc
    ∑ k ∈ range (m + 1),
        (-1 : ℚ) ^ k * (m.choose k : ℚ) * P.eval (k : ℚ)
      = (-1 : ℚ) ^ m *
          ∑ k ∈ range (m + 1),
            (-1 : ℚ) ^ (m - k) * (m.choose k : ℚ) * P.eval (k : ℚ) := by
              rw [mul_sum]
              apply sum_congr rfl
              intro k hk
              have hk_le : k ≤ m := by
                simp only [mem_range] at hk
                omega
              rw [negOnePow_sub_eq_mul m k hk_le]
              ring_nf
              rw [show m * 2 = 2 * m by omega, pow_mul]
              norm_num
    _ = 0 := by rw [hq, mul_zero]

/-- A summand in the residual beta sum. -/
def betaTerm (m k : ℕ) : ℚ :=
  (-1 : ℚ) ^ k * (m.choose k : ℚ) / (2 * k + 1)

/-- The residual sum left after polynomial finite-difference cancellation. -/
def betaSum (m : ℕ) : ℚ :=
  ∑ k ∈ range (m + 1), betaTerm m k

/-- Generalized beta summand for a pole of arbitrary positive order. -/
def repeatedBetaTerm (m s k : ℕ) : ℚ :=
  (-1 : ℚ) ^ k * (m.choose k : ℚ) / (2 * k + 1) ^ s

/-- Generalized beta sum governing a repeated pole at `X = -1/2`. -/
def repeatedBetaSum (m s : ℕ) : ℚ :=
  ∑ k ∈ range (m + 1), repeatedBetaTerm m s k

lemma repeatedBetaSum_zero_left (s : ℕ) : repeatedBetaSum 0 s = 1 := by
  simp [repeatedBetaSum, repeatedBetaTerm]

lemma repeatedBetaSum_zero_right (m : ℕ) (hm : 0 < m) :
    repeatedBetaSum m 0 = 0 := by
  have hdegree : (1 : Polynomial ℚ).natDegree < m := by
    simpa using hm
  simpa [repeatedBetaSum, repeatedBetaTerm] using
    alternatingPolynomialSum_vanishes (1 : Polynomial ℚ) m hdegree

lemma repeatedBetaSum_one (m : ℕ) :
    repeatedBetaSum m 1 = betaSum m := by
  simp [repeatedBetaSum, repeatedBetaTerm, betaSum, betaTerm]

/-- The mixed finite sum before removing its polynomial part. -/
def mixedTemplateSum (m : ℕ) (A : Polynomial ℚ) : ℚ :=
  ∑ k ∈ range (m + 1),
    (-1 : ℚ) ^ k * (m.choose k : ℚ) *
      (A.eval (k : ℚ) / (2 * k + 1))

/-- Alternating quotient transform with an arbitrary polynomial denominator. -/
def quotientTransform (m : ℕ) (L A : Polynomial ℚ) : ℚ :=
  ∑ k ∈ range (m + 1),
    (-1 : ℚ) ^ k * (m.choose k : ℚ) *
      (A.eval (k : ℚ) / L.eval (k : ℚ))

/-- The quotient transform ignores a denominator multiple whose quotient has
degree less than the finite-difference order. -/
theorem quotientTransform_add_mul_invariant
    (m : ℕ) (L A R D : Polynomial ℚ)
    (hfactor : A = R + L * D)
    (hD : D.natDegree < m)
    (hL : ∀ k ∈ range (m + 1), L.eval (k : ℚ) ≠ 0) :
    quotientTransform m L A = quotientTransform m L R := by
  have hvanish := alternatingPolynomialSum_vanishes D m hD
  rw [quotientTransform, quotientTransform]
  calc
    ∑ k ∈ range (m + 1),
        (-1 : ℚ) ^ k * (m.choose k : ℚ) *
          (A.eval (k : ℚ) / L.eval (k : ℚ))
      = ∑ k ∈ range (m + 1),
          ((-1 : ℚ) ^ k * (m.choose k : ℚ) *
              (R.eval (k : ℚ) / L.eval (k : ℚ)) +
            (-1 : ℚ) ^ k * (m.choose k : ℚ) * D.eval (k : ℚ)) := by
              apply sum_congr rfl
              intro k hk
              rw [hfactor]
              simp only [Polynomial.eval_add, Polynomial.eval_mul]
              have hden := hL k hk
              field_simp [hden]
    _ = (∑ k ∈ range (m + 1),
          (-1 : ℚ) ^ k * (m.choose k : ℚ) *
            (R.eval (k : ℚ) / L.eval (k : ℚ))) +
        ∑ k ∈ range (m + 1),
          (-1 : ℚ) ^ k * (m.choose k : ℚ) * D.eval (k : ℚ) := by
            rw [sum_add_distrib]
    _ = ∑ k ∈ range (m + 1),
          (-1 : ℚ) ^ k * (m.choose k : ℚ) *
            (R.eval (k : ℚ) / L.eval (k : ℚ)) := by
              rw [hvanish, add_zero]

/-- The quotient transform factors through Mathlib's monic remainder
operation whenever the corresponding quotient has degree below the
finite-difference order.  For a nonmonic denominator this statement is
tautological because `divByMonic` is zero. -/
theorem quotientTransform_eq_modByMonic_of_quotient_degree
    (m : ℕ) (L A : Polynomial ℚ)
    (hquotient : (A /ₘ L).natDegree < m)
    (hL : ∀ k ∈ range (m + 1), L.eval (k : ℚ) ≠ 0) :
    quotientTransform m L A = quotientTransform m L (A %ₘ L) := by
  apply quotientTransform_add_mul_invariant m L A (A %ₘ L) (A /ₘ L)
  · exact (Polynomial.modByMonic_add_div A L).symm
  · exact hquotient
  · exact hL

/-- A convenient degree criterion for the remainder principle. -/
theorem quotientTransform_eq_modByMonic_of_degree
    (m : ℕ) (L A : Polynomial ℚ)
    (hm : 0 < m)
    (hLmonic : L.Monic)
    (hAdegree : A.natDegree < m + L.natDegree)
    (hL : ∀ k ∈ range (m + 1), L.eval (k : ℚ) ≠ 0) :
    quotientTransform m L A = quotientTransform m L (A %ₘ L) := by
  apply quotientTransform_eq_modByMonic_of_quotient_degree m L A
  · rw [Polynomial.natDegree_divByMonic A hLmonic]
    omega
  · exact hL

theorem mixedTemplate_eq_endpoint_mul_betaSum
    (m : ℕ) (A D : Polynomial ℚ) (c : ℚ)
    (hfactor :
      A = Polynomial.C c +
        (Polynomial.C 2 * Polynomial.X + Polynomial.C 1) * D)
    (hD : D.natDegree < m) :
    mixedTemplateSum m A = c * betaSum m := by
  have hvanish := alternatingPolynomialSum_vanishes D m hD
  rw [mixedTemplateSum, betaSum]
  calc
    ∑ k ∈ range (m + 1),
        (-1 : ℚ) ^ k * (m.choose k : ℚ) *
          (A.eval (k : ℚ) / (2 * k + 1))
      = ∑ k ∈ range (m + 1),
          ((-1 : ℚ) ^ k * (m.choose k : ℚ) * D.eval (k : ℚ) +
            c * betaTerm m k) := by
              apply sum_congr rfl
              intro k _
              rw [hfactor]
              simp only [Polynomial.eval_add, Polynomial.eval_mul,
                Polynomial.eval_C, Polynomial.eval_X, betaTerm]
              have hden : (2 * (k : ℚ) + 1) ≠ 0 := by positivity
              field_simp [hden]
              ring
    _ = (∑ k ∈ range (m + 1),
          (-1 : ℚ) ^ k * (m.choose k : ℚ) * D.eval (k : ℚ)) +
        ∑ k ∈ range (m + 1), c * betaTerm m k := by
          rw [sum_add_distrib]
    _ = c * ∑ k ∈ range (m + 1), betaTerm m k := by
      rw [hvanish, zero_add, mul_sum]

theorem mixedTemplate_eq_betaSum (m : ℕ) (A D : Polynomial ℚ)
    (hfactor :
      A = Polynomial.C 1 +
        (Polynomial.C 2 * Polynomial.X + Polynomial.C 1) * D)
    (hD : D.natDegree < m) :
    mixedTemplateSum m A = betaSum m := by
  simpa using mixedTemplate_eq_endpoint_mul_betaSum m A D 1 hfactor hD

theorem mixedTemplate_eq_endpoint_mul_betaSum_of_degree
    (m : ℕ) (A : Polynomial ℚ) (hAdegree : A.natDegree < m) :
    mixedTemplateSum m A =
      A.eval (-1 / 2 : ℚ) * betaSum m := by
  obtain ⟨q, hq⟩ :=
    Polynomial.X_sub_C_dvd_sub_C_eval (p := A) (a := (-1 / 2 : ℚ))
  have hq_factor :
      A - Polynomial.C (A.eval (-1 / 2 : ℚ)) =
        (Polynomial.X - Polynomial.C (-1 / 2 : ℚ)) * q := by
    exact hq
  let D : Polynomial ℚ := Polynomial.C (1 / 2 : ℚ) * q
  have hfactor :
      A = Polynomial.C (A.eval (-1 / 2 : ℚ)) +
        (Polynomial.C 2 * Polynomial.X + Polynomial.C 1) * D := by
    dsimp [D]
    calc
      A = Polynomial.C (A.eval (-1 / 2 : ℚ)) +
          (A - Polynomial.C (A.eval (-1 / 2 : ℚ))) := by ring
      _ = Polynomial.C (A.eval (-1 / 2 : ℚ)) +
          (Polynomial.X - Polynomial.C (-1 / 2 : ℚ)) * q := by
            rw [hq_factor]
      _ = Polynomial.C (A.eval (-1 / 2 : ℚ)) +
          (Polynomial.C 2 * Polynomial.X + Polynomial.C 1) *
            (Polynomial.C (1 / 2 : ℚ) * q) := by
              have htwo :
                  Polynomial.C (2 : ℚ) * Polynomial.C (1 / 2 : ℚ) = 1 := by
                norm_num [← Polynomial.C_mul]
              have hminus :
                  Polynomial.C (-1 / 2 : ℚ) =
                    -Polynomial.C (1 / 2 : ℚ) := by norm_num
              simp only [Polynomial.C_1]
              rw [hminus]
              ring_nf
              have htwo' :
                  Polynomial.C (1 / 2 : ℚ) * Polynomial.C (2 : ℚ) = 1 := by
                rw [mul_comm, htwo]
              rw [show
                Polynomial.X * Polynomial.C (1 / 2 : ℚ) * q *
                    Polynomial.C (2 : ℚ) =
                  Polynomial.X * q *
                    (Polynomial.C (1 / 2 : ℚ) * Polynomial.C (2 : ℚ)) by
                      ring,
                htwo', mul_one]
  have hqdegree : q.natDegree ≤ A.natDegree := by
    by_cases hqzero : q = 0
    · simp [hqzero]
    · have hlinear :
          Polynomial.X - Polynomial.C (-1 / 2 : ℚ) ≠ 0 := by
        intro h
        have := congrArg (fun P : Polynomial ℚ ↦ P.coeff 1) h
        norm_num at this
      have hq_le_product :
          q.natDegree ≤
            ((Polynomial.X - Polynomial.C (-1 / 2 : ℚ)) * q).natDegree := by
        rw [Polynomial.natDegree_mul hlinear hqzero]
        omega
      have hproduct_le :
          ((Polynomial.X - Polynomial.C (-1 / 2 : ℚ)) * q).natDegree
            ≤ A.natDegree := by
        rw [← hq_factor]
        exact
          (Polynomial.natDegree_sub_le A
            (Polynomial.C (A.eval (-1 / 2 : ℚ)))).trans (by simp)
      exact hq_le_product.trans hproduct_le
  have hDdegree : D.natDegree < m := by
    apply lt_of_le_of_lt _ hAdegree
    exact (Polynomial.natDegree_C_mul_le (1 / 2 : ℚ) q).trans hqdegree
  exact mixedTemplate_eq_endpoint_mul_betaSum m A D
    (A.eval (-1 / 2 : ℚ)) hfactor hDdegree

theorem mixedTemplate_eq_betaSum_of_endpoint (m : ℕ) (A : Polynomial ℚ)
    (hAdegree : A.natDegree < m)
    (hAendpoint : A.eval (-1 / 2 : ℚ) = 1) :
    mixedTemplateSum m A = betaSum m := by
  rw [mixedTemplate_eq_endpoint_mul_betaSum_of_degree m A hAdegree,
    hAendpoint, one_mul]

/-- The normalized product `A_{n+1}` from the SIC2C4 finite sums. -/
noncomputable def coefficientPoly : ℕ → Polynomial ℚ
  | 0 => 1
  | n + 1 =>
      Polynomial.C (1 / (n + 1 : ℚ)) *
        (Polynomial.C 2 * Polynomial.X +
          Polynomial.C ((n + 2 : ℕ) : ℚ)) *
          coefficientPoly n

lemma coefficientPoly_degree (n : ℕ) :
    (coefficientPoly n).natDegree ≤ n := by
  induction n with
  | zero => simp [coefficientPoly]
  | succ n inductionHypothesis =>
      let L : Polynomial ℚ :=
        Polynomial.C 2 * Polynomial.X +
          Polynomial.C ((n + 2 : ℕ) : ℚ)
      rw [coefficientPoly, mul_assoc]
      change
        (Polynomial.C (1 / (n + 1 : ℚ)) *
          (L * coefficientPoly n)).natDegree ≤ n + 1
      have hL : L.natDegree ≤ 1 := by
        dsimp [L]
        exact (Polynomial.natDegree_add_le _ _).trans (by norm_num)
      calc
        (Polynomial.C (1 / (n + 1 : ℚ)) *
            (L * coefficientPoly n)).natDegree
          ≤ (L * coefficientPoly n).natDegree :=
            Polynomial.natDegree_C_mul_le _ _
        _ ≤ L.natDegree + (coefficientPoly n).natDegree :=
          Polynomial.natDegree_mul_le
        _ ≤ 1 + n := Nat.add_le_add hL inductionHypothesis
        _ = n + 1 := Nat.add_comm 1 n

lemma coefficientPoly_endpoint (n : ℕ) :
    (coefficientPoly n).eval (-1 / 2 : ℚ) = 1 := by
  induction n with
  | zero => simp [coefficientPoly]
  | succ n inductionHypothesis =>
      rw [coefficientPoly]
      simp only [Polynomial.eval_mul, Polynomial.eval_add, Polynomial.eval_C,
        Polynomial.eval_X, inductionHypothesis]
      have hnonzero : (n + 1 : ℚ) ≠ 0 := by positivity
      field_simp [hnonzero]
      push_cast
      ring

/-- At a nonnegative integer, the recursive product polynomial is exactly
the binomial coefficient occurring in the chart constant term. -/
theorem coefficientPoly_eval_nat (n k : ℕ) :
    (coefficientPoly n).eval (k : ℚ) =
      ((n + 1 + 2 * k).choose n : ℚ) := by
  induction n with
  | zero => simp [coefficientPoly]
  | succ n inductionHypothesis =>
      rw [coefficientPoly]
      simp only [Polynomial.eval_mul, Polynomial.eval_add, Polynomial.eval_C,
        Polynomial.eval_X, inductionHypothesis]
      have hden : (n + 1 : ℚ) ≠ 0 := by positivity
      have hchooseNat := Nat.add_one_mul_choose_eq (n + 1 + 2 * k) n
      have hchooseNat' :
          (n + 2 + 2 * k) * (n + 1 + 2 * k).choose n =
            (n + 2 + 2 * k).choose (n + 1) * (n + 1) := by
        simpa only [Nat.succ_eq_add_one,
          show n + 1 + 2 * k + 1 = n + 2 + 2 * k by omega] using hchooseNat
      have hchoose :
          (n + 2 + 2 * k : ℚ) *
              ((n + 1 + 2 * k).choose n : ℚ) =
            ((n + 2 + 2 * k).choose (n + 1) : ℚ) *
              (n + 1 : ℚ) := by
        exact_mod_cast hchooseNat'
      field_simp [hden]
      push_cast
      linear_combination hchoose

/-- The pure finite sum in product-polynomial form vanishes. -/
theorem pureProductSum_vanishes (m : ℕ) (hm : 0 < m) :
    (∑ k ∈ range (m + 1),
      (-1 : ℚ) ^ k * (m.choose k : ℚ) *
        (Polynomial.C (1 / (m : ℚ)) *
          coefficientPoly (m - 1)).eval (k : ℚ)) = 0 := by
  apply alternatingPolynomialSum_vanishes
  have hdegree :
      (Polynomial.C (1 / (m : ℚ)) *
        coefficientPoly (m - 1)).natDegree ≤ m - 1 := by
    exact (Polynomial.natDegree_C_mul_le
      (1 / (m : ℚ)) (coefficientPoly (m - 1))).trans
        (coefficientPoly_degree (m - 1))
  omega

/-- The pure displayed binomial quotient sum vanishes. -/
theorem pureBinomialSum_vanishes (m : ℕ) (hm : 0 < m) :
    (∑ k ∈ range (m + 1),
      (-1 : ℚ) ^ k * (m.choose k : ℚ) *
        (((m + 2 * k).choose m : ℚ) / (2 * k + 1))) = 0 := by
  rw [← pureProductSum_vanishes m hm]
  apply sum_congr rfl
  intro k _
  simp only [Polynomial.eval_mul, Polynomial.eval_C,
    coefficientPoly_eval_nat]
  have htop : m - 1 + 1 + 2 * k = m + 2 * k := by omega
  rw [htop]
  have hchooseNat := Nat.choose_succ_right_eq (m + 2 * k) (m - 1)
  have hpred : m - 1 + 1 = m := by omega
  have hsub : m + 2 * k - (m - 1) = 2 * k + 1 := by omega
  rw [hpred, hsub] at hchooseNat
  have hchoose :
      ((m + 2 * k).choose m : ℚ) * (m : ℚ) =
        ((m + 2 * k).choose (m - 1) : ℚ) * (2 * k + 1 : ℚ) := by
    exact_mod_cast hchooseNat
  have hmQ : (m : ℚ) ≠ 0 := by positivity
  have hodd : (2 * k + 1 : ℚ) ≠ 0 := by positivity
  have hratio :
      ((m + 2 * k).choose m : ℚ) / (2 * k + 1) =
        (1 / (m : ℚ)) * ((m + 2 * k).choose (m - 1) : ℚ) := by
    field_simp [hmQ, hodd]
    simpa [mul_comm] using hchoose
  rw [hratio]

/-- The mixed finite sum in product-polynomial form is the residual beta sum. -/
theorem mixedProductSum_eq_betaSum (m : ℕ) (hm : 0 < m) :
    mixedTemplateSum m (coefficientPoly (m - 1)) = betaSum m := by
  apply mixedTemplate_eq_betaSum_of_endpoint
  · exact lt_of_le_of_lt (coefficientPoly_degree (m - 1)) (by omega)
  · exact coefficientPoly_endpoint (m - 1)

/-- The mixed displayed binomial sum is the residual beta sum. -/
theorem mixedBinomialSum_eq_betaSum (m : ℕ) (hm : 0 < m) :
    (∑ k ∈ range (m + 1),
      (-1 : ℚ) ^ k * (m.choose k : ℚ) *
        (((m + 2 * k).choose (m - 1) : ℚ) / (2 * k + 1))) =
      betaSum m := by
  rw [← mixedProductSum_eq_betaSum m hm, mixedTemplateSum]
  apply sum_congr rfl
  intro k _
  rw [coefficientPoly_eval_nat]
  have htop : m - 1 + 1 + 2 * k = m + 2 * k := by omega
  rw [htop]

lemma betaSum_zero : betaSum 0 = 1 := by
  norm_num [betaSum, betaTerm]

lemma choose_predecessor_mul (n k : ℕ) :
    (n.choose k : ℚ) * (n + 1 : ℕ) =
      ((n + 1).choose k : ℚ) * ((n + 1 - k : ℕ) : ℚ) := by
  exact_mod_cast Nat.choose_mul_succ_eq n k

/-- Triangular creative-telescoping recurrence for every repeated-pole
beta sum; equation (4.8c) in the written proof. -/
theorem repeatedBetaSum_recurrence (n s : ℕ) :
    (2 * (n + 1) + 1 : ℚ) * repeatedBetaSum (n + 1) (s + 1) =
      repeatedBetaSum (n + 1) s +
        (2 * (n + 1) : ℚ) * repeatedBetaSum n (s + 1) := by
  have h_extend :
      repeatedBetaSum n (s + 1) =
        ∑ k ∈ range (n + 2), repeatedBetaTerm n (s + 1) k := by
    calc
      repeatedBetaSum n (s + 1) =
          ∑ k ∈ range (n + 1), repeatedBetaTerm n (s + 1) k := rfl
      _ = (∑ k ∈ range (n + 1), repeatedBetaTerm n (s + 1) k) +
          repeatedBetaTerm n (s + 1) (n + 1) := by
            simp [repeatedBetaTerm]
      _ = ∑ k ∈ range ((n + 1) + 1),
          repeatedBetaTerm n (s + 1) k :=
            (sum_range_succ
              (fun k ↦ repeatedBetaTerm n (s + 1) k) (n + 1)).symm
      _ = ∑ k ∈ range (n + 2), repeatedBetaTerm n (s + 1) k := by
        congr 2
  have h_term (k : ℕ) (hk : k ∈ range (n + 2)) :
      (2 * (n + 1) + 1 : ℚ) * repeatedBetaTerm (n + 1) (s + 1) k =
        repeatedBetaTerm (n + 1) s k +
          (2 * (n + 1) : ℚ) * repeatedBetaTerm n (s + 1) k := by
    have hk_le : k ≤ n + 1 := by
      simp only [mem_range] at hk
      omega
    have hchoose := choose_predecessor_mul n k
    rw [Nat.cast_sub hk_le] at hchoose
    have hinner :
        (2 * (n + 1) + 1 : ℚ) * ((n + 1).choose k : ℚ) -
            (2 * (n + 1) : ℚ) * (n.choose k : ℚ) =
          (2 * k + 1 : ℚ) * ((n + 1).choose k : ℚ) := by
      push_cast at hchoose ⊢
      linear_combination -2 * hchoose
    have hsplit :
        (2 * (n + 1) + 1 : ℚ) * ((n + 1).choose k : ℚ) =
          (2 * k + 1 : ℚ) * ((n + 1).choose k : ℚ) +
            (2 * (n + 1) : ℚ) * (n.choose k : ℚ) := by
      linarith [hinner]
    simp only [repeatedBetaTerm]
    calc
      (2 * (n + 1) + 1 : ℚ) *
          ((-1 : ℚ) ^ k * ((n + 1).choose k : ℚ) /
            (2 * k + 1) ^ (s + 1))
        = ((-1 : ℚ) ^ k / (2 * k + 1) ^ (s + 1)) *
            ((2 * (n + 1) + 1 : ℚ) *
              ((n + 1).choose k : ℚ)) := by ring
      _ = ((-1 : ℚ) ^ k / (2 * k + 1) ^ (s + 1)) *
            ((2 * k + 1 : ℚ) * ((n + 1).choose k : ℚ) +
              (2 * (n + 1) : ℚ) * (n.choose k : ℚ)) := by
                rw [hsplit]
      _ = (-1 : ℚ) ^ k * ((n + 1).choose k : ℚ) /
              (2 * k + 1) ^ s +
            (2 * (n + 1) : ℚ) *
              ((-1 : ℚ) ^ k * (n.choose k : ℚ) /
                (2 * k + 1) ^ (s + 1)) := by
                  have hden : (2 * k + 1 : ℚ) ≠ 0 := by positivity
                  rw [pow_succ]
                  field_simp [hden]
  rw [repeatedBetaSum, repeatedBetaSum, h_extend, mul_sum, mul_sum,
    ← sum_add_distrib]
  apply sum_congr rfl
  intro k hk
  exact h_term k hk

/-- A remainder expanded in powers of `2X+1` contributes exactly the
corresponding finite repeated-pole jet; the finite-sum part of (4.8d). -/
theorem repeatedPoleRemainderSum_eq_jet
    (m r : ℕ) (c : ℕ → ℚ) :
    (∑ k ∈ range (m + 1),
      (-1 : ℚ) ^ k * (m.choose k : ℚ) *
        ((∑ j ∈ range r, c j * (2 * k + 1 : ℚ) ^ j) /
          (2 * k + 1 : ℚ) ^ r)) =
      ∑ j ∈ range r, c j * repeatedBetaSum m (r - j) := by
  calc
    ∑ k ∈ range (m + 1),
        (-1 : ℚ) ^ k * (m.choose k : ℚ) *
          ((∑ j ∈ range r, c j * (2 * k + 1 : ℚ) ^ j) /
            (2 * k + 1 : ℚ) ^ r)
      = ∑ k ∈ range (m + 1),
          ∑ j ∈ range r, c j * repeatedBetaTerm m (r - j) k := by
            apply sum_congr rfl
            intro k _
            rw [div_eq_mul_inv, sum_mul, mul_sum]
            apply sum_congr rfl
            intro j hj
            have hj_lt : j < r := by simpa using hj
            have hden : (2 * k + 1 : ℚ) ≠ 0 := by positivity
            have hpow :
                (2 * k + 1 : ℚ) ^ r =
                  (2 * k + 1 : ℚ) ^ (r - j) *
                    (2 * k + 1 : ℚ) ^ j := by
              rw [← pow_add]
              congr 1
              omega
            rw [repeatedBetaTerm, hpow]
            field_simp [hden]
    _ = ∑ j ∈ range r,
          ∑ k ∈ range (m + 1), c j * repeatedBetaTerm m (r - j) k := by
            rw [sum_comm]
    _ = ∑ j ∈ range r, c j * repeatedBetaSum m (r - j) := by
      apply sum_congr rfl
      intro j _
      rw [repeatedBetaSum, mul_sum]

/-- Creative-telescoping recurrence for the residual beta sum. -/
theorem betaSum_recurrence (n : ℕ) :
    (2 * (n + 1) + 1 : ℚ) * betaSum (n + 1) =
      (2 * (n + 1) : ℚ) * betaSum n := by
  have h := repeatedBetaSum_recurrence n 0
  rw [repeatedBetaSum_one, repeatedBetaSum_one,
    repeatedBetaSum_zero_right (n + 1) (by omega), zero_add] at h
  simpa using h

/-- The product determined by the beta recurrence. -/
def betaProduct : ℕ → ℚ
  | 0 => 1
  | n + 1 =>
      (2 * (n + 1) : ℚ) / (2 * (n + 1) + 1) * betaProduct n

/-- Closed product evaluation of the residual finite sum. -/
theorem betaSum_eq_betaProduct (m : ℕ) : betaSum m = betaProduct m := by
  induction m with
  | zero => simp [betaSum_zero, betaProduct]
  | succ n inductionHypothesis =>
      rw [betaProduct, ← inductionHypothesis]
      have h := betaSum_recurrence n
      have hden : (2 * (n + 1) + 1 : ℚ) ≠ 0 := by positivity
      rw [div_mul_eq_mul_div]
      apply (eq_div_iff hden).2
      nlinarith [h]

/-- Factorial form of the closed beta product used in the mixed moment. -/
theorem betaProduct_eq_factorial_ratio (m : ℕ) :
    betaProduct m =
      ((2 : ℚ) ^ m * (m.factorial : ℚ)) /
        (Nat.doubleFactorial (2 * m + 1) : ℚ) := by
  induction m with
  | zero => norm_num [betaProduct, Nat.doubleFactorial]
  | succ n inductionHypothesis =>
      rw [betaProduct, inductionHypothesis, pow_succ, Nat.factorial_succ]
      rw [show 2 * (n + 1) + 1 = (2 * n + 1) + 2 by omega,
        Nat.doubleFactorial_add_two]
      have hodd :
          (Nat.doubleFactorial (2 * n + 1) : ℚ) ≠ 0 := by positivity
      have hnext :
          ((2 * n + 1 + 2 : ℕ) : ℚ) ≠ 0 := by positivity
      field_simp [hodd, hnext]
      push_cast
      ring

/-- Explicit all-order evaluation of the residual beta sum. -/
theorem betaSum_eq_factorial_ratio (m : ℕ) :
    betaSum m =
      ((2 : ℚ) ^ m * (m.factorial : ℚ)) /
        (Nat.doubleFactorial (2 * m + 1) : ℚ) := by
  rw [betaSum_eq_betaProduct, betaProduct_eq_factorial_ratio]

/-- Equation (4.4), including its normalization, has value zero. -/
theorem pureNormalizedBinomialSum_vanishes (m : ℕ) (hm : 0 < m) :
    (1 / (2 : ℚ) ^ m) *
      (∑ k ∈ range (m + 1),
        (-1 : ℚ) ^ k * (m.choose k : ℚ) *
          (((m + 2 * k).choose m : ℚ) / (2 * k + 1))) = 0 := by
  rw [pureBinomialSum_vanishes m hm, mul_zero]

/-- Equation (4.5), including its normalization and closed evaluation. -/
theorem mixedNormalizedBinomialSum_eq_factorial_ratio
    (m : ℕ) (hm : 0 < m) :
    (1 / (2 : ℚ) ^ m) *
      (∑ k ∈ range (m + 1),
        (-1 : ℚ) ^ k * (m.choose k : ℚ) *
          (((m + 2 * k).choose (m - 1) : ℚ) / (2 * k + 1))) =
      (m.factorial : ℚ) /
        (Nat.doubleFactorial (2 * m + 1) : ℚ) := by
  rw [mixedBinomialSum_eq_betaSum m hm, betaSum_eq_factorial_ratio]
  have hpow : (2 : ℚ) ^ m ≠ 0 := by positivity
  have hdouble :
      (Nat.doubleFactorial (2 * m + 1) : ℚ) ≠ 0 := by positivity
  field_simp [hpow, hdouble]

/-- The pure chart polynomial has zero formal integral in every positive
order. -/
theorem formalIntegral_pureChartConstantTerm_vanishes
    (m : ℕ) (hm : 0 < m) :
    formalIntegral (pureChartConstantTerm m) = 0 := by
  rw [formalIntegral_pureChartConstantTerm,
    pureNormalizedBinomialSum_vanishes m hm]

/-- The mixed chart polynomial has exactly the normalized detector value. -/
theorem formalIntegral_mixedChartConstantTerm_eq_factorial_ratio
    (m : ℕ) (hm : 0 < m) :
    formalIntegral (mixedChartConstantTerm m) =
      (m.factorial : ℚ) /
        (Nat.doubleFactorial (2 * m + 1) : ℚ) := by
  rw [formalIntegral_mixedChartConstantTerm,
    mixedNormalizedBinomialSum_eq_factorial_ratio m hm]

/-- End-to-end pure chart moment: the coefficientwise formal integral of
the literal Laurent witness power vanishes in every positive order. -/
theorem formalIntegral_chartLaurentWitness_pow_constantTerm_vanishes
    (m : ℕ) (hm : 0 < m) :
    formalIntegral ((chartLaurentWitness ^ m).coeff 0) = 0 := by
  rw [chartLaurentWitness_pow_constantTerm,
    formalIntegral_pureChartConstantTerm_vanishes m hm]

/-- End-to-end mixed chart moment: the coefficientwise formal integral of
`x F^m` has the required double-factorial value. -/
theorem formalIntegral_chartX_mul_chartLaurentWitness_pow_constantTerm
    (m : ℕ) (hm : 0 < m) :
    formalIntegral
        ((LaurentPolynomial.T 1 * chartLaurentWitness ^ m).coeff 0) =
      (m.factorial : ℚ) /
        (Nat.doubleFactorial (2 * m + 1) : ℚ) := by
  rw [chartX_mul_chartLaurentWitness_pow_constantTerm m hm,
    formalIntegral_mixedChartConstantTerm_eq_factorial_ratio m hm]

/-- The same two endpoints stated directly for the displayed four-variable
polynomials after the formal contraction-chart substitution. -/
theorem formalIntegral_chartSubstitution_witnessF_pow_constantTerm_vanishes
    (m : ℕ) (hm : 0 < m) :
    formalIntegral ((chartSubstitution (witnessF ^ m)).coeff 0) = 0 := by
  rw [chartSubstitution_witnessF_pow,
    formalIntegral_chartLaurentWitness_pow_constantTerm_vanishes m hm]

theorem
    formalIntegral_chartSubstitution_witnessQ_mul_witnessF_pow_constantTerm
    (m : ℕ) (hm : 0 < m) :
    formalIntegral
        ((chartSubstitution (witnessQ * witnessF ^ m)).coeff 0) =
      (m.factorial : ℚ) /
        (Nat.doubleFactorial (2 * m + 1) : ℚ) := by
  rw [chartSubstitution_witnessQ_mul_witnessF_pow,
    formalIntegral_chartX_mul_chartLaurentWitness_pow_constantTerm m hm]

end FiniteEtaleKeller.SIC2C4
