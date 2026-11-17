/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.Admissibility
import FiniteEtaleKeller.RealizationFiber

/-!
# The unchanged-parameter compiler for the universal fixed map

The universal relative quadratic-gauge map has parameters `u₄,...,u_N` and
target coordinates `π,b,c`.  Promoting the `u_j` to unchanged affine
coordinates gives the inverse polynomial

`S + b*S² + π*S³ + ∑ u_j*π^j*S^j - c/2`.

This module formalizes the coefficient compiler

`π=h₃, b=h₂, c=-2h₀, u_j=h_j/h₃^j`

and proves that it reconstructs every normalized polynomial
`H=h₀+S+h₂S²+...+h_NS^N` with `h₃ ≠ 0`.  It then applies the identity to the
normalized translate `P(a+S)/P'(a)`.

The polynomial identity and admissible-translation input are formalized here.
The promoted `N`-variable Jacobian, generic `S_N` monodromy, and
primitive-monodromy atomicity are separate theorem-level inputs.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

/-- Target and unchanged parameter values compiled from one normalized
inverse polynomial.  Values below index four are harmless and unused. -/
structure UniversalPromotedTarget (K : Type*) [Field K] where
  parameter : ℕ → K
  pi : K
  b : K
  c : K

/-- Compile the coefficients of a normalized polynomial into the target of
the promoted fixed map. -/
def compileUniversalPromotedTarget (H : K[X]) : UniversalPromotedTarget K where
  parameter j := H.coeff j / H.coeff 3 ^ j
  pi := H.coeff 3
  b := H.coeff 2
  c := -2 * H.coeff 0

/-- The inverse polynomial selected by a promoted target, truncated at `N`. -/
def universalPromotedTail
    (N : ℕ) (target : UniversalPromotedTarget K) : K[X] :=
  ∑ j ∈ Finset.Icc 4 N,
    C (target.parameter j * target.pi ^ j) * X ^ j

/-- The inverse polynomial selected by a promoted target, truncated at `N`. -/
def universalPromotedInversePolynomial
    (N : ℕ) (target : UniversalPromotedTarget K) : K[X] :=
  X + C target.b * X ^ 2 + C target.pi * X ^ 3 +
    universalPromotedTail N target -
    C (target.c / 2)

@[simp]
theorem compileUniversalPromotedTarget_pi (H : K[X]) :
    (compileUniversalPromotedTarget H).pi = H.coeff 3 := rfl

@[simp]
theorem compileUniversalPromotedTarget_b (H : K[X]) :
    (compileUniversalPromotedTarget H).b = H.coeff 2 := rfl

@[simp]
theorem compileUniversalPromotedTarget_c (H : K[X]) :
    (compileUniversalPromotedTarget H).c = -2 * H.coeff 0 := rfl

@[simp]
theorem compileUniversalPromotedTarget_parameter (H : K[X]) (j : ℕ) :
    (compileUniversalPromotedTarget H).parameter j =
      H.coeff j / H.coeff 3 ^ j := rfl

/-- On the nonzero-cubic chart, the compiled high-degree tail has exactly the
high-degree coefficients of `H`. -/
private theorem universalPromotedTail_compile_coeff
    (H : K[X]) (hthree : H.coeff 3 ≠ 0) (n : ℕ) :
    (universalPromotedTail H.natDegree
      (compileUniversalPromotedTarget H)).coeff n =
        if n ∈ Finset.Icc 4 H.natDegree then H.coeff n else 0 := by
  classical
  simp [universalPromotedTail, compileUniversalPromotedTarget,
    pow_ne_zero _ hthree]

/-- The promoted coefficient compiler reconstructs a normalized polynomial
coefficientwise. -/
theorem universalPromotedInversePolynomial_compile
    [CharZero K]
    (H : K[X]) (hone : H.coeff 1 = 1) (hthree : H.coeff 3 ≠ 0) :
    universalPromotedInversePolynomial H.natDegree
        (compileUniversalPromotedTarget H) = H := by
  ext n
  by_cases h0 : n = 0
  · subst n
    rw [universalPromotedInversePolynomial, Polynomial.coeff_sub,
      Polynomial.coeff_add, universalPromotedTail_compile_coeff H hthree]
    norm_num [compileUniversalPromotedTarget, Polynomial.coeff_X]
    field_simp
  by_cases h1 : n = 1
  · subst n
    rw [universalPromotedInversePolynomial, Polynomial.coeff_sub,
      Polynomial.coeff_add, universalPromotedTail_compile_coeff H hthree]
    simp [compileUniversalPromotedTarget, hone]
  by_cases h2 : n = 2
  · subst n
    rw [universalPromotedInversePolynomial, Polynomial.coeff_sub,
      Polynomial.coeff_add, universalPromotedTail_compile_coeff H hthree]
    norm_num [compileUniversalPromotedTarget, Polynomial.coeff_X]
  by_cases h3 : n = 3
  · subst n
    rw [universalPromotedInversePolynomial, Polynomial.coeff_sub,
      Polynomial.coeff_add, universalPromotedTail_compile_coeff H hthree]
    norm_num [compileUniversalPromotedTarget, Polynomial.coeff_X]
  have hn4 : 4 ≤ n := by omega
  by_cases hnN : n ≤ H.natDegree
  · have hnmem : n ∈ Finset.Icc 4 H.natDegree :=
      Finset.mem_Icc.mpr ⟨hn4, hnN⟩
    rw [universalPromotedInversePolynomial, Polynomial.coeff_sub,
      Polynomial.coeff_add, universalPromotedTail_compile_coeff H hthree]
    have hone' : (1 : ℕ) ≠ n := by omega
    simp [Polynomial.coeff_X, Polynomial.coeff_C, h0, hone', h2, h3, hnmem]
  · have hcoeff : H.coeff n = 0 := by
      by_contra hne
      have hle := Polynomial.le_natDegree_of_ne_zero hne
      omega
    have hnmem : n ∉ Finset.Icc 4 H.natDegree := by
      simp [Finset.mem_Icc, hn4, hnN]
    rw [universalPromotedInversePolynomial, Polynomial.coeff_sub,
      Polynomial.coeff_add, universalPromotedTail_compile_coeff H hthree]
    have hone' : (1 : ℕ) ≠ n := by omega
    simp [Polynomial.coeff_X, Polynomial.coeff_C, h0, hone', h2, h3,
      hnmem, hcoeff]

/-- The selected inverse polynomial has the full degree of the normalized
presentation. -/
theorem universalPromotedInversePolynomial_compile_natDegree
    [CharZero K]
    (H : K[X]) (hone : H.coeff 1 = 1) (hthree : H.coeff 3 ≠ 0) :
    (universalPromotedInversePolynomial H.natDegree
      (compileUniversalPromotedTarget H)).natDegree = H.natDegree := by
  rw [universalPromotedInversePolynomial_compile H hone hthree]

/-- A nonzero normalized coefficient compiles to a nonzero unchanged
parameter on the nonzero-cubic chart.  In particular the top coefficient of a
degree-`N` presentation gives `u_N ≠ 0`. -/
theorem compileUniversalPromotedTarget_parameter_ne_zero
    (H : K[X]) (j : ℕ) (hcoeff : H.coeff j ≠ 0)
    (hthree : H.coeff 3 ≠ 0) :
    (compileUniversalPromotedTarget H).parameter j ≠ 0 := by
  exact div_ne_zero hcoeff (pow_ne_zero j hthree)

/-- Normalize the translated polynomial by its nonzero linear coefficient. -/
def normalizedTranslatedPolynomial
    (P : K[X]) (a : K) : K[X] :=
  C (P.derivative.eval a)⁻¹ * translatePolynomial P a

@[simp]
theorem normalizedTranslatedPolynomial_coeff
    (P : K[X]) (a : K) (j : ℕ) :
    (normalizedTranslatedPolynomial P a).coeff j =
      (P.derivative.eval a)⁻¹ *
        (Polynomial.hasseDeriv j P).eval a := by
  change
    (C (P.derivative.eval a)⁻¹ * Polynomial.taylor a P).coeff j =
      (P.derivative.eval a)⁻¹ *
        (Polynomial.hasseDeriv j P).eval a
  simp [Polynomial.taylor_coeff]

@[simp]
theorem normalizedTranslatedPolynomial_coeff_one
    (P : K[X]) (a : K) (hone : P.derivative.eval a ≠ 0) :
    (normalizedTranslatedPolynomial P a).coeff 1 = 1 := by
  rw [normalizedTranslatedPolynomial_coeff]
  simpa using inv_mul_cancel₀ hone

theorem normalizedTranslatedPolynomial_coeff_three_ne_zero
    (P : K[X]) (a : K)
    (hone : P.derivative.eval a ≠ 0)
    (hthree : (Polynomial.hasseDeriv 3 P).eval a ≠ 0) :
    (normalizedTranslatedPolynomial P a).coeff 3 ≠ 0 := by
  rw [normalizedTranslatedPolynomial_coeff]
  exact mul_ne_zero (inv_ne_zero hone) hthree

/-- The fixed-map compiler reconstructs the normalized translated
presentation `P(a+S)/P'(a)` exactly. -/
theorem universalPromotedInversePolynomial_realization
    [CharZero K]
    (P : K[X]) (a : K)
    (hone : P.derivative.eval a ≠ 0)
    (hthree : (Polynomial.hasseDeriv 3 P).eval a ≠ 0) :
    universalPromotedInversePolynomial
        (normalizedTranslatedPolynomial P a).natDegree
        (compileUniversalPromotedTarget
          (normalizedTranslatedPolynomial P a)) =
      normalizedTranslatedPolynomial P a :=
  universalPromotedInversePolynomial_compile
    (normalizedTranslatedPolynomial P a)
    (normalizedTranslatedPolynomial_coeff_one P a hone)
    (normalizedTranslatedPolynomial_coeff_three_ne_zero
      P a hone hthree)

/-- The automatic admissible translation supplies the fixed-map compiler
without any remaining nonvanishing input. -/
theorem automaticUniversalPromotedInversePolynomial_realization
    [CharZero K]
    (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    universalPromotedInversePolynomial
        (normalizedTranslatedPolynomial P
          (chosenAdmissibleTranslation P hdeg)).natDegree
        (compileUniversalPromotedTarget
          (normalizedTranslatedPolynomial P
            (chosenAdmissibleTranslation P hdeg))) =
      normalizedTranslatedPolynomial P
        (chosenAdmissibleTranslation P hdeg) :=
  universalPromotedInversePolynomial_realization
    P (chosenAdmissibleTranslation P hdeg)
    (chosenAdmissibleTranslation_linear_ne_zero P hdeg)
    (chosenAdmissibleTranslation_cubic_ne_zero P hdeg)

#print axioms universalPromotedInversePolynomial_compile
#print axioms universalPromotedInversePolynomial_compile_natDegree
#print axioms universalPromotedInversePolynomial_realization
#print axioms automaticUniversalPromotedInversePolynomial_realization

end FiniteEtaleKeller
