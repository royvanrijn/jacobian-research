/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.UniversalFiber
import FiniteEtaleKeller.TranslationQuotient
import FiniteEtaleKeller.Admissibility
import Mathlib.Algebra.Polynomial.Inductions

/-!
# The translated finite-étale realization theorem

This file instantiates the universal represented-fiber theorem with the inverse
polynomial

`E(S) = P(a + S)`

and transports the representing quotient across the canonical translation

`K[S]/(P(a+S)) ≃ K[T]/(P(T))`.

The final result is natural in every commutative test algebra: algebra maps out
of `K[T]/(P)` are equivalent to full quadratic-gauge source-fiber points, and
the equivalence commutes with postcomposition.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

/-- Differentiation commutes with translation of the polynomial variable. -/
theorem translatePolynomial_derivative (P : K[X]) (a : K) :
    (translatePolynomial P a).derivative =
      translatePolynomial P.derivative a := by
  simp [translatePolynomial, Polynomial.derivative_comp]

/-- The linear coefficient of `P(a+S)` is `P'(a)`. -/
@[simp]
theorem translatePolynomial_coeff_one (P : K[X]) (a : K) :
    (translatePolynomial P a).coeff 1 = P.derivative.eval a := by
  change (Polynomial.taylor a P).coeff 1 = P.derivative.eval a
  simp

/-- Translation preserves separability. -/
theorem translatePolynomial_separable (P : K[X]) (a : K)
    (hP : P.Separable) : (translatePolynomial P a).Separable := by
  rw [Polynomial.separable_def'] at hP ⊢
  rcases hP with ⟨u, v, h⟩
  refine ⟨translatePolynomial u a, translatePolynomial v a, ?_⟩
  rw [translatePolynomial_derivative]
  have ht := congrArg (fun q : K[X] => translatePolynomial q a) h
  simpa [translatePolynomial] using ht

/-- The part of the normalized translated derivative left after removing the
constant and prescribed quadratic chart terms. -/
def realizationDerivativeRemainder (P : K[X]) (a : K) : K[X] :=
  C ((P.derivative.eval a)⁻¹) * (translatePolynomial P a).derivative - 1 - X ^ 2

/-- The marked polynomial obtained by dividing the zero-constant remainder by
`S`.  This is the polynomial form of
`(E'(S)/P'(a) - 1 - S²)/S`. -/
def realizationBeta (P : K[X]) (a : K) : K[X] :=
  (realizationDerivativeRemainder P a).divX

@[simp]
theorem realizationDerivativeRemainder_coeff_zero
    (P : K[X]) (a : K) (h₁ : P.derivative.eval a ≠ 0) :
    (realizationDerivativeRemainder P a).coeff 0 = 0 := by
  simp [realizationDerivativeRemainder, Polynomial.coeff_derivative, h₁]

/-- The defining identity for the marked polynomial. -/
theorem X_mul_realizationBeta
    (P : K[X]) (a : K) (h₁ : P.derivative.eval a ≠ 0) :
    X * realizationBeta P a = realizationDerivativeRemainder P a := by
  simpa [realizationBeta, realizationDerivativeRemainder_coeff_zero P a h₁]
    using Polynomial.X_mul_divX_add (realizationDerivativeRemainder P a)

/-- The translated derivative has exactly the factorization required by the
universal quadratic-gauge datum at first target coordinate `1` and marked
target coordinate `0`. -/
theorem translatePolynomial_derivative_eq_markedChartPolynomial
    (P : K[X]) (a : K) (h₁ : P.derivative.eval a ≠ 0) :
    (translatePolynomial P a).derivative =
      C (P.derivative.eval a) *
        markedChartPolynomial 1 0 (realizationBeta P a) := by
  rw [markedChartPolynomial]
  simp only [map_zero, zero_mul, sub_zero, map_one, one_mul]
  have hβ := X_mul_realizationBeta P a h₁
  rw [realizationDerivativeRemainder] at hβ
  have hg :
      C (P.derivative.eval a) * C ((P.derivative.eval a)⁻¹) =
        (1 : K[X]) := by
    rw [← C_mul]
    simp [h₁]
  calc
    (translatePolynomial P a).derivative =
        (C (P.derivative.eval a) * C ((P.derivative.eval a)⁻¹)) *
          (translatePolynomial P a).derivative := by rw [hg, one_mul]
    _ = C (P.derivative.eval a) *
          (C ((P.derivative.eval a)⁻¹) *
            (translatePolynomial P a).derivative) := by ring
    _ = C (P.derivative.eval a) *
          (1 + X * realizationBeta P a + X ^ 2) := by rw [hβ]; ring

/-- The nonzero cubic Taylor coefficient, carried as a unit. -/
def realizationCubicUnit (P : K[X]) (a : K)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0) : Kˣ :=
  Units.mk0 ((Polynomial.hasseDeriv 3 P).eval a) h₃

/-- The source coefficient `g₁/g₃` used by the quadratic-gauge reconstruction. -/
def realizationGaugeCoefficient (P : K[X]) (a : K)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0) : K :=
  P.derivative.eval a / (realizationCubicUnit P a h₃ : K)

section CharacteristicZero

variable [CharZero K]

/-- The universal quadratic-gauge datum attached to a squarefree polynomial and
an admissible translation parameter.  Its inverse polynomial is `P(a+S)`. -/
def realizationDatum
    (P : K[X]) (a : K) (hP : Squarefree P)
    (h₁ : P.derivative.eval a ≠ 0)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0) :
    QuadraticGaugeFiberDatum K where
  E := translatePolynomial P a
  β := realizationBeta P a
  pi := 1
  b := 0
  a := realizationGaugeCoefficient P a h₃
  g₁ := Units.mk0 (P.derivative.eval a) h₁
  separable := translatePolynomial_separable P a
    ((PerfectField.separable_iff_squarefree).2 hP)
  derivative_eq := by
    change (translatePolynomial P a).derivative =
      C (P.derivative.eval a) *
        markedChartPolynomial 1 0 (realizationBeta P a)
    exact translatePolynomial_derivative_eq_markedChartPolynomial P a h₁

@[simp]
theorem realizationDatum_E
    (P : K[X]) (a : K) (hP : Squarefree P)
    (h₁ : P.derivative.eval a ≠ 0)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0) :
    (realizationDatum P a hP h₁ h₃).E = translatePolynomial P a := rfl

/-- At the prescribed third target coordinate, the inverse polynomial used by
the datum is exactly `P(a+S)`. -/
theorem realizationDatum_inversePolynomial
    (P : K[X]) (a : K) (hP : Squarefree P)
    (h₁ : P.derivative.eval a ≠ 0)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0) :
    rootedTranslate P a -
        C (P.derivative.eval a / 2 *
          realizationTargetC P a (P.derivative.eval a)) =
      (realizationDatum P a hP h₁ h₃).E := by
  change rootedTranslate P a -
      C (P.derivative.eval a / 2 *
        realizationTargetC P a (P.derivative.eval a)) =
    translatePolynomial P a
  exact rootedTranslate_inverse_at_target P a (P.derivative.eval a) h₁

end CharacteristicZero

variable {A B : Type*}
variable [CommRing A] [Algebra K A]
variable [CommRing B] [Algebra K B]

/-- Precomposition with the translation quotient equivalence identifies maps
out of the original and translated quotients. -/
def translatedQuotientHomEquiv (P : K[X]) (a : K) :
    (AdjoinRoot P →ₐ[K] A) ≃
      (AdjoinRoot (translatePolynomial P a) →ₐ[K] A) where
  toFun := fun φ => φ.comp (translationQuotientEquiv P a).toAlgHom
  invFun := fun ψ => ψ.comp (translationQuotientEquiv P a).symm.toAlgHom
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

section RepresentedRealization

variable [CharZero K]

/-- The quotient-hom equivalence with its translated source type exposed through
the realization datum.  Keeping this dependent type explicit makes subsequent
naturality statements robust under elaboration. -/
def realizationQuotientHomEquiv
    (P : K[X]) (a : K) (hP : Squarefree P)
    (h₁ : P.derivative.eval a ≠ 0)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0) :
    (AdjoinRoot P →ₐ[K] A) ≃
      (AdjoinRoot (realizationDatum P a hP h₁ h₃).E →ₐ[K] A) := by
  change (AdjoinRoot P →ₐ[K] A) ≃
    (AdjoinRoot (translatePolynomial P a) →ₐ[K] A)
  exact translatedQuotientHomEquiv P a

/-- Translation of the representing quotient commutes with postcomposition. -/
theorem realizationQuotientHomEquiv_natural
    (P : K[X]) (a : K) (hP : Squarefree P)
    (h₁ : P.derivative.eval a ≠ 0)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0)
    (f : A →ₐ[K] B) (φ : AdjoinRoot P →ₐ[K] A) :
    f.comp (realizationQuotientHomEquiv (A := A) P a hP h₁ h₃ φ) =
      realizationQuotientHomEquiv (A := B) P a hP h₁ h₃ (f.comp φ) := by
  apply DFunLike.ext _ _
  intro x
  rfl

/-- The translated finite-étale realization theorem in represented
functor-of-points form.  The source fiber is represented directly by
`K[T]/(P)`, not merely by the translated quotient. -/
def translatedFiberRepresentingEquiv
    (P : K[X]) (a : K) (hP : Squarefree P)
    (h₁ : P.derivative.eval a ≠ 0)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0) :
    (AdjoinRoot P →ₐ[K] A) ≃
      (realizationDatum P a hP h₁ h₃).Point A :=
  (realizationQuotientHomEquiv (A := A) P a hP h₁ h₃).trans
    ((realizationDatum P a hP h₁ h₃).representingEquiv A)

/-- The represented realization equivalence is natural under every morphism of
commutative test `K`-algebras. -/
theorem translatedFiberRepresentingEquiv_natural
    (P : K[X]) (a : K) (hP : Squarefree P)
    (h₁ : P.derivative.eval a ≠ 0)
    (h₃ : (Polynomial.hasseDeriv 3 P).eval a ≠ 0)
    (f : A →ₐ[K] B) (φ : AdjoinRoot P →ₐ[K] A) :
    GaugeFiberPoint.map f
        (translatedFiberRepresentingEquiv (A := A) P a hP h₁ h₃ φ) =
      translatedFiberRepresentingEquiv (A := B) P a hP h₁ h₃ (f.comp φ) := by
  change
    GaugeFiberPoint.map f
        ((realizationDatum P a hP h₁ h₃).representingEquiv A
          (realizationQuotientHomEquiv (A := A) P a hP h₁ h₃ φ)) =
      (realizationDatum P a hP h₁ h₃).representingEquiv B
        (realizationQuotientHomEquiv (A := B) P a hP h₁ h₃ (f.comp φ))
  rw [(realizationDatum P a hP h₁ h₃).representingEquiv_natural]
  rw [realizationQuotientHomEquiv_natural]

#print axioms translatedFiberRepresentingEquiv_natural

end RepresentedRealization

end FiniteEtaleKeller
