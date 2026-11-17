/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.TranslationQuotient
import FiniteEtaleKeller.UniversalParameterCompiler

/-!
# Quotient algebra realized by the unchanged-parameter compiler

The universal compiler reconstructs the scalar-normalized translate

`P(a + S) / P'(a)`.

Multiplication by a nonzero scalar does not change the quotient algebra, and
translation identifies the resulting quotient with `K[T]/(P)`.  This module
formalizes both steps and composes them into the algebra equivalence needed by
the fixed-map realization argument.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

/-- Multiplication of a defining polynomial by a scalar induces a map from
the scaled quotient to the original quotient. -/
def scaledToOriginal (Q : K[X]) (r : K) :
    AdjoinRoot (C r * Q) →ₐ[K] AdjoinRoot Q :=
  AdjoinRoot.liftAlgHom (C r * Q)
    (Algebra.ofId K (AdjoinRoot Q))
    (AdjoinRoot.root Q) (by
      change Polynomial.aeval (AdjoinRoot.root Q) (C r * Q) = 0
      simp)

/-- When the scalar is nonzero, the original quotient maps back to the
scaled quotient. -/
def originalToScaled (Q : K[X]) (r : K) (hr : r ≠ 0) :
    AdjoinRoot Q →ₐ[K] AdjoinRoot (C r * Q) :=
  AdjoinRoot.liftAlgHom Q
    (Algebra.ofId K (AdjoinRoot (C r * Q)))
    (AdjoinRoot.root (C r * Q)) (by
      let A := AdjoinRoot (C r * Q)
      have hscaled :
          algebraMap K A r *
              Polynomial.aeval (AdjoinRoot.root (C r * Q)) Q = 0 := by
        rw [Polynomial.aeval_def]
        simpa [A] using
          (AdjoinRoot.eval₂_root (C r * Q))
      calc
        Polynomial.aeval (AdjoinRoot.root (C r * Q)) Q =
            1 * Polynomial.aeval (AdjoinRoot.root (C r * Q)) Q := by
              simp
        _ = (algebraMap K A r⁻¹ * algebraMap K A r) *
              Polynomial.aeval (AdjoinRoot.root (C r * Q)) Q := by
              rw [← map_mul]
              simp [hr]
        _ = algebraMap K A r⁻¹ *
              (algebraMap K A r *
                Polynomial.aeval (AdjoinRoot.root (C r * Q)) Q) := by
              rw [mul_assoc]
        _ = 0 := by rw [hscaled, mul_zero])

@[simp]
theorem scaledToOriginal_root (Q : K[X]) (r : K) :
    scaledToOriginal Q r (AdjoinRoot.root (C r * Q)) =
      AdjoinRoot.root Q := by
  simp [scaledToOriginal]

@[simp]
theorem originalToScaled_root (Q : K[X]) (r : K) (hr : r ≠ 0) :
    originalToScaled Q r hr (AdjoinRoot.root Q) =
      AdjoinRoot.root (C r * Q) := by
  simp [originalToScaled]

/-- Scaling a defining polynomial by a nonzero scalar does not change its
quotient algebra. -/
def scalingQuotientEquiv (Q : K[X]) (r : K) (hr : r ≠ 0) :
    AdjoinRoot (C r * Q) ≃ₐ[K] AdjoinRoot Q :=
  AlgEquiv.ofAlgHom (scaledToOriginal Q r) (originalToScaled Q r hr)
    (by
      apply AdjoinRoot.algHom_ext
      simp [scaledToOriginal, originalToScaled])
    (by
      apply AdjoinRoot.algHom_ext
      simp [scaledToOriginal, originalToScaled])

@[simp]
theorem scalingQuotientEquiv_root
    (Q : K[X]) (r : K) (hr : r ≠ 0) :
    scalingQuotientEquiv Q r hr (AdjoinRoot.root (C r * Q)) =
      AdjoinRoot.root Q := by
  simp [scalingQuotientEquiv]

/-- The compiler's normalized translated quotient is canonically equivalent
to the quotient by the original polynomial. -/
def normalizedTranslationQuotientEquiv
    (P : K[X]) (a : K) (hone : P.derivative.eval a ≠ 0) :
    AdjoinRoot (normalizedTranslatedPolynomial P a) ≃ₐ[K] AdjoinRoot P :=
  (scalingQuotientEquiv (translatePolynomial P a)
      (P.derivative.eval a)⁻¹ (inv_ne_zero hone)).trans
    (translationQuotientEquiv P a)

@[simp]
theorem normalizedTranslationQuotientEquiv_root
    (P : K[X]) (a : K) (hone : P.derivative.eval a ≠ 0) :
    normalizedTranslationQuotientEquiv P a hone
        (AdjoinRoot.root (normalizedTranslatedPolynomial P a)) =
      AdjoinRoot.root P - algebraMap K (AdjoinRoot P) a := by
  change
    translationQuotientEquiv P a
        (scalingQuotientEquiv (translatePolynomial P a)
          (P.derivative.eval a)⁻¹ (inv_ne_zero hone)
          (AdjoinRoot.root
            (C (P.derivative.eval a)⁻¹ * translatePolynomial P a))) =
      AdjoinRoot.root P - algebraMap K (AdjoinRoot P) a
  rw [scalingQuotientEquiv_root, translationQuotientEquiv_root]

/-- The automatically chosen admissible translation realizes the original
quotient algebra without any remaining nonvanishing hypothesis. -/
def automaticUniversalPromotedQuotientEquiv
    [CharZero K] (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    AdjoinRoot
        (normalizedTranslatedPolynomial P
          (chosenAdmissibleTranslation P hdeg)) ≃ₐ[K]
      AdjoinRoot P :=
  normalizedTranslationQuotientEquiv P
    (chosenAdmissibleTranslation P hdeg)
    (chosenAdmissibleTranslation_linear_ne_zero P hdeg)

#print axioms scalingQuotientEquiv
#print axioms normalizedTranslationQuotientEquiv
#print axioms automaticUniversalPromotedQuotientEquiv

end FiniteEtaleKeller
