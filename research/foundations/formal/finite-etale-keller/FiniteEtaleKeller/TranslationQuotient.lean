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

variable {A : Type*} [CommRing A] [Algebra K A]

/-- Precomposition with translation identifies algebra maps out of the
original and translated polynomial quotients. -/
def translationQuotientHomEquiv (P : K[X]) (a : K) :
    (AdjoinRoot P →ₐ[K] A) ≃
      (AdjoinRoot (translatePolynomial P a) →ₐ[K] A) where
  toFun := fun φ => φ.comp (translationQuotientEquiv P a).toAlgHom
  invFun := fun ψ =>
    ψ.comp (translationQuotientEquiv P a).symm.toAlgHom
  left_inv := by
    intro φ
    apply DFunLike.ext _ _
    intro x
    exact congrArg φ ((translationQuotientEquiv P a).apply_symm_apply x)
  right_inv := by
    intro ψ
    apply DFunLike.ext _ _
    intro x
    exact congrArg ψ ((translationQuotientEquiv P a).symm_apply_apply x)

/-- Quotient translation commutes with postcomposition in the test
algebra. -/
theorem translationQuotientHomEquiv_natural
    {B : Type*} [CommRing B] [Algebra K B]
    (P : K[X]) (a : K) (f : A →ₐ[K] B)
    (φ : AdjoinRoot P →ₐ[K] A) :
    f.comp (translationQuotientHomEquiv (A := A) P a φ) =
      translationQuotientHomEquiv (A := B) P a (f.comp φ) := rfl

end FiniteEtaleKeller
