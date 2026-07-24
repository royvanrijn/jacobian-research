/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GMC2.LowerFace
import Mathlib.Algebra.MvPolynomial.Eval

/-!
# Ordinary bivariate polynomials

This file exposes the circular formalization through the usual polynomial
ring in two variables.  Variable `0` is `Z`, variable `1` is `W`, and the
substitution is

`Z ↦ T`, `W ↦ U T⁻¹`.

The Wick monomial formula follows directly from the definition.
-/

namespace GMC2

open Polynomial LaurentPolynomial MvPolynomial

abbrev BivariatePolynomial (K : Type*) [CommRing K] :=
  MvPolynomial (Fin 2) K

/-- Circular-coordinate substitution `Z ↦ T`, `W ↦ U T⁻¹`. -/
noncomputable def circularize
    (K : Type*) [CommRing K] :
    BivariatePolynomial K →+* CircularPolynomial K :=
  MvPolynomial.eval₂Hom
    ((LaurentPolynomial.C : K[X] →+* CircularPolynomial K).comp
      Polynomial.C)
    (Fin.cases (LaurentPolynomial.T 1)
      (fun _ ↦ LaurentPolynomial.C Polynomial.X * LaurentPolynomial.T (-1)))

@[simp] theorem circularize_Z
    {K : Type*} [CommRing K] :
    circularize K (MvPolynomial.X (0 : Fin 2)) =
      LaurentPolynomial.T 1 := by
  simp [circularize]

@[simp] theorem circularize_W
    {K : Type*} [CommRing K] :
    circularize K (MvPolynomial.X (1 : Fin 2)) =
      LaurentPolynomial.C Polynomial.X * LaurentPolynomial.T (-1) := by
  rw [circularize, MvPolynomial.eval₂Hom_X']
  rfl

/-- Gaussian expectation on ordinary bivariate polynomials. -/
noncomputable def bivariateGaussianExpectation
    {K : Type*} [CommRing K] (P : BivariatePolynomial K) : K :=
  circularExpectation (circularize K P)

def BivariatePureMomentsVanish
    {K : Type*} [CommRing K] (P : BivariatePolynomial K) : Prop :=
  ∀ m : ℕ, 0 < m → bivariateGaussianExpectation (P ^ m) = 0

def BivariateEventuallyMixedMomentsVanish
    {K : Type*} [CommRing K] (P : BivariatePolynomial K) : Prop :=
  ∀ Q : BivariatePolynomial K,
    ∃ M : ℕ, ∀ m ≥ M,
      bivariateGaussianExpectation (Q * P ^ m) = 0

/-- Wick's formula in the normalized circular variables. -/
theorem bivariateGaussianExpectation_monomial
    {K : Type*} [CommRing K] (a b : ℕ) :
    bivariateGaussianExpectation
        (MvPolynomial.X (0 : Fin 2) ^ a *
          MvPolynomial.X (1 : Fin 2) ^ b : BivariatePolynomial K) =
      if a = b then (Nat.factorial a : K) else 0 := by
  classical
  rw [bivariateGaussianExpectation, map_mul, map_pow, map_pow,
    circularize_Z, circularize_W, mul_pow]
  simp only [LaurentPolynomial.T_pow]
  rw [mul_one, mul_neg, mul_one]
  rw [← mul_assoc, LaurentPolynomial.T_mul]
  rw [← map_pow (LaurentPolynomial.C : K[X] →+* CircularPolynomial K)
    Polynomial.X b]
  rw [LaurentPolynomial.mul_T_assoc]
  by_cases hab : a = b
  · subst b
    simp only [add_neg_cancel, LaurentPolynomial.T_zero, mul_one, if_pos]
    rw [circularExpectation, constantTerm_eq_coeff,
      LaurentPolynomial.C_apply]
    rw [← Polynomial.monomial_one_right_eq_X_pow]
    rw [if_pos rfl]
    rw [factorialFunctional_monomial, one_mul]
  · have habZ : (a : ℤ) + -(b : ℤ) ≠ 0 := by
      exact sub_ne_zero.mpr (by exact_mod_cast hab)
    have hct :
        (LaurentPolynomial.C (Polynomial.X ^ b : K[X]) *
          (LaurentPolynomial.T ((a : ℤ) + -(b : ℤ)) :
            CircularPolynomial K)).coeff 0 = 0 := by
      change
        (AddMonoidAlgebra.single 0 (Polynomial.X ^ b : K[X]) *
          AddMonoidAlgebra.single ((a : ℤ) + -(b : ℤ)) 1).coeff 0 = 0
      rw [AddMonoidAlgebra.single_mul_single]
      simp only [zero_add, mul_one, AddMonoidAlgebra.coeff_single,
        Finsupp.single_apply]
      exact if_neg habZ
    rw [circularExpectation, constantTerm_eq_coeff, hct]
    simp [factorialFunctional, hab]

/-- The Gaussian Moments Conjecture for two ordinary polynomial
variables. -/
theorem gaussianMomentsConjecture_two_variables_bivariate
    {K : Type*} [Field K] [CharZero K]
    (P : BivariatePolynomial K)
    (hvan : BivariatePureMomentsVanish P) :
    BivariateEventuallyMixedMomentsVanish P := by
  have hcirc : PureMomentsVanish (circularize K P) := by
    intro m hm
    simpa [pureMoment, bivariateGaussianExpectation] using hvan m hm
  have hmain :=
    gaussianMomentsConjecture_two_variables (circularize K P) hcirc
  intro Q
  obtain ⟨M, hM⟩ := hmain (circularize K Q)
  refine ⟨M, ?_⟩
  intro m hm
  simpa [mixedMoment, bivariateGaussianExpectation] using hM m hm

end GMC2
