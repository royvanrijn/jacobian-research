/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.Translation
import Mathlib.RingTheory.AdjoinRoot

/-!
# Translation equivalence of polynomial quotients

The inverse polynomial produced by the realization construction is
`P(a + S)`.  This module constructs the canonical algebra equivalence

`K[S]/(P(a+S)) ≃ K[T]/(P(T))`

by sending `S` to `T-a`, with inverse `T ↦ S+a`.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

/-- Map the translated quotient to the original quotient by `S ↦ T-a`. -/
def translatedToOriginal (P : K[X]) (a : K) :
    AdjoinRoot (translatePolynomial P a) →ₐ[K] AdjoinRoot P :=
  AdjoinRoot.liftAlgHom (translatePolynomial P a)
    (Algebra.ofId K (AdjoinRoot P))
    (AdjoinRoot.root P - algebraMap K (AdjoinRoot P) a) (by
      change Polynomial.aeval
        (AdjoinRoot.root P - algebraMap K (AdjoinRoot P) a)
        (translatePolynomial P a) = 0
      rw [translatePolynomial, Polynomial.aeval_comp]
      simp)

/-- Map the original quotient to the translated quotient by `T ↦ S+a`. -/
def originalToTranslated (P : K[X]) (a : K) :
    AdjoinRoot P →ₐ[K] AdjoinRoot (translatePolynomial P a) :=
  AdjoinRoot.liftAlgHom P
    (Algebra.ofId K (AdjoinRoot (translatePolynomial P a)))
    (AdjoinRoot.root (translatePolynomial P a)
      + algebraMap K (AdjoinRoot (translatePolynomial P a)) a) (by
      change Polynomial.aeval
        (AdjoinRoot.root (translatePolynomial P a)
          + algebraMap K (AdjoinRoot (translatePolynomial P a)) a) P = 0
      calc
        Polynomial.aeval
            (AdjoinRoot.root (translatePolynomial P a)
              + algebraMap K (AdjoinRoot (translatePolynomial P a)) a) P
          = Polynomial.aeval (AdjoinRoot.root (translatePolynomial P a))
              (P.comp (X + C a)) := by
                rw [Polynomial.aeval_comp]
                simp
        _ = Polynomial.aeval (AdjoinRoot.root (translatePolynomial P a))
              (translatePolynomial P a) := rfl
        _ = 0 := by simp)

@[simp]
theorem translatedToOriginal_root (P : K[X]) (a : K) :
    translatedToOriginal P a (AdjoinRoot.root (translatePolynomial P a)) =
      AdjoinRoot.root P - algebraMap K (AdjoinRoot P) a := by
  simp [translatedToOriginal]

@[simp]
theorem originalToTranslated_root (P : K[X]) (a : K) :
    originalToTranslated P a (AdjoinRoot.root P) =
      AdjoinRoot.root (translatePolynomial P a)
        + algebraMap K (AdjoinRoot (translatePolynomial P a)) a := by
  simp [originalToTranslated]

/-- Translation does not change the quotient algebra. -/
def translationQuotientEquiv (P : K[X]) (a : K) :
    AdjoinRoot (translatePolynomial P a) ≃ₐ[K] AdjoinRoot P :=
  AlgEquiv.ofAlgHom (translatedToOriginal P a) (originalToTranslated P a)
    (by
      apply AdjoinRoot.algHom_ext
      simp [translatedToOriginal, originalToTranslated])
    (by
      apply AdjoinRoot.algHom_ext
      simp [translatedToOriginal, originalToTranslated])

@[simp]
theorem translationQuotientEquiv_root (P : K[X]) (a : K) :
    translationQuotientEquiv P a
        (AdjoinRoot.root (translatePolynomial P a)) =
      AdjoinRoot.root P - algebraMap K (AdjoinRoot P) a := by
  simpa [translationQuotientEquiv] using translatedToOriginal_root P a

@[simp]
theorem translationQuotientEquiv_symm_root (P : K[X]) (a : K) :
    (translationQuotientEquiv P a).symm (AdjoinRoot.root P) =
      AdjoinRoot.root (translatePolynomial P a)
        + algebraMap K (AdjoinRoot (translatePolynomial P a)) a := by
  simpa [translationQuotientEquiv] using originalToTranslated_root P a

end FiniteEtaleKeller
