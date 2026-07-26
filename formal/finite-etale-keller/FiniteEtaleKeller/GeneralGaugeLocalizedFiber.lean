/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.LocalizedGaugeFiberPoints
import FiniteEtaleKeller.GeneralGaugeNormalization

/-!
# The localized fiber of the actual all-degree Keller map

For a target with invertible first coordinate, the literal polynomial-map
fiber is represented by the inverse-polynomial quotient localized away from
the derivative.  This theorem has no separability hypothesis.  The earlier
finite étale theorem is the special case in which the derivative class is
already a unit modulo the inverse polynomial.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K A B : Type*} [Field K] [CharZero K]
variable [CommRing A] [Algebra K A]
variable [CommRing B] [Algebra K B]

/-- Roots of the inverse polynomial with invertible derivative are equivalent
to literal points of the unnormalized all-degree map fiber. -/
def generalGaugeLocalizedRawFiberEquiv
    (G : K[X]) (pi : Kˣ) (b c : K)
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (A : Type*) [CommRing A] [Algebra K A] :
    LocalizedPolynomialRoot
        (generalGaugeInversePolynomial G (pi : K) b c) A ≃
      GeneralGaugeRawFiberPoint G pi b c A :=
  (localizedRootEquivGaugeFiberPoint (A := A)
      (G.coeff 1 / G.coeff 3)
      (Units.mk0 (G.coeff 1) h₁)
      (generalGaugeInversePolynomial_derivative G (pi : K) b c h₁)).trans
    ((generalGaugeDisplayedFiberEquiv (A := A) h₁ h₃).symm.trans
      (generalGaugeRawFiberEquivDisplayed (A := A)).symm)

/-- The localized-root/literal-fiber equivalence is natural in the test
algebra. -/
theorem generalGaugeLocalizedRawFiberEquiv_natural
    (G : K[X]) (pi : Kˣ) (b c : K)
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (f : A →ₐ[K] B)
    (s : LocalizedPolynomialRoot
      (generalGaugeInversePolynomial G (pi : K) b c) A) :
    GeneralGaugeRawFiberPoint.map f
        (generalGaugeLocalizedRawFiberEquiv G pi b c h₁ h₃ A s) =
      generalGaugeLocalizedRawFiberEquiv G pi b c h₁ h₃ B (s.map f) := by
  apply GeneralGaugeRawFiberPoint.ext
  funext i
  have h := localizedRootEquivGaugeFiberPoint_natural
    (A := A) (B := B)
    (G.coeff 1 / G.coeff 3)
    (Units.mk0 (G.coeff 1) h₁)
    (generalGaugeInversePolynomial_derivative G (pi : K) b c h₁)
    f s
  fin_cases i
  · exact congrArg (fun p => p.source.x) h
  · exact congrArg (fun p => p.source.y) h
  · exact congrArg (fun p => p.source.z) h

/-- The localized inverse quotient represents the literal map fiber at every
target whose first coordinate is a unit. -/
def generalGaugeLocalizedRawRepresentingEquiv
    (G : K[X]) (pi : Kˣ) (b c : K)
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (A : Type*) [CommRing A] [Algebra K A] :
    (LocalizedAdjoinRoot
        (generalGaugeInversePolynomial G (pi : K) b c) →ₐ[K] A) ≃
      GeneralGaugeRawFiberPoint G pi b c A :=
  (LocalizedPolynomialRoot.localizedAlgHomEquiv
      (generalGaugeInversePolynomial G (pi : K) b c) A).trans
    (generalGaugeLocalizedRawFiberEquiv G pi b c h₁ h₃ A)

/-- Naturality of the localized literal-fiber representation. -/
theorem generalGaugeLocalizedRawRepresentingEquiv_natural
    (G : K[X]) (pi : Kˣ) (b c : K)
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (f : A →ₐ[K] B)
    (φ : LocalizedAdjoinRoot
      (generalGaugeInversePolynomial G (pi : K) b c) →ₐ[K] A) :
    GeneralGaugeRawFiberPoint.map f
        (generalGaugeLocalizedRawRepresentingEquiv G pi b c h₁ h₃ A φ) =
      generalGaugeLocalizedRawRepresentingEquiv G pi b c h₁ h₃ B (f.comp φ) := by
  change GeneralGaugeRawFiberPoint.map f
      (generalGaugeLocalizedRawFiberEquiv G pi b c h₁ h₃ A
        (LocalizedPolynomialRoot.localizedAlgHomEquiv
          (generalGaugeInversePolynomial G (pi : K) b c) A φ)) =
    generalGaugeLocalizedRawFiberEquiv G pi b c h₁ h₃ B
      (LocalizedPolynomialRoot.localizedAlgHomEquiv
        (generalGaugeInversePolynomial G (pi : K) b c) B (f.comp φ))
  rw [generalGaugeLocalizedRawFiberEquiv_natural]
  congr 1
  exact (LocalizedPolynomialRoot.localizedAlgHomEquiv_natural f φ).symm

/-- At a zero second target coordinate, the same localized quotient represents
the literal fiber of the determinant-one normalized map. -/
def generalGaugeLocalizedJacobianOneRepresentingEquiv
    (G : K[X]) (pi : Kˣ) (c : K)
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (A : Type*) [CommRing A] [Algebra K A] :
    (LocalizedAdjoinRoot
        (generalGaugeInversePolynomial G (pi : K) 0 c) →ₐ[K] A) ≃
      GeneralGaugeJacobianOneFiberPoint G pi c A :=
  (generalGaugeLocalizedRawRepresentingEquiv G pi 0 c h₁ h₃ A).trans
    (generalGaugeRawFiberEquivJacobianOne (A := A))

/-- Naturality of the determinant-one localized-fiber theorem. -/
theorem generalGaugeLocalizedJacobianOneRepresentingEquiv_natural
    (G : K[X]) (pi : Kˣ) (c : K)
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (f : A →ₐ[K] B)
    (φ : LocalizedAdjoinRoot
      (generalGaugeInversePolynomial G (pi : K) 0 c) →ₐ[K] A) :
    GeneralGaugeJacobianOneFiberPoint.map f
        (generalGaugeLocalizedJacobianOneRepresentingEquiv
          G pi c h₁ h₃ A φ) =
      generalGaugeLocalizedJacobianOneRepresentingEquiv
        G pi c h₁ h₃ B (f.comp φ) := by
  apply GeneralGaugeJacobianOneFiberPoint.ext
  funext i
  have h := generalGaugeLocalizedRawRepresentingEquiv_natural
    G pi 0 c h₁ h₃ f φ
  fin_cases i
  · exact congrArg (fun p => p.point 0) h
  · exact congrArg (fun p => p.point 1) h
  · exact congrArg (fun p => p.point 2) h

#print axioms generalGaugeLocalizedRawRepresentingEquiv
#print axioms generalGaugeLocalizedRawRepresentingEquiv_natural
#print axioms generalGaugeLocalizedJacobianOneRepresentingEquiv
#print axioms generalGaugeLocalizedJacobianOneRepresentingEquiv_natural

end FiniteEtaleKeller
