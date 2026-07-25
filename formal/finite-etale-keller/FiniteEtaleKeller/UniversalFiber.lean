/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.FiberNaturality

/-!
# The universal finite-étale fiber theorem

A quadratic-gauge fiber datum consists only of a separable inverse polynomial,
the marked polynomial, the target parameters, and the derivative
factorization.  The resulting full source fiber functor is naturally
represented by the polynomial quotient `K[S]/(E)`.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K A B : Type*} [Field K]
variable [CommRing A] [Algebra K A]
variable [CommRing B] [Algebra K B]

/-- The minimal data needed for the universal quadratic-gauge fiber theorem. -/
structure QuadraticGaugeFiberDatum (K : Type*) [Field K] where
  E : K[X]
  β : K[X]
  pi : K
  b : K
  a : K
  g₁ : Kˣ
  separable : E.Separable
  derivative_eq :
    E.derivative = C (g₁ : K) * markedChartPolynomial pi b β

namespace QuadraticGaugeFiberDatum

/-- Full source fiber points over a commutative test algebra. -/
def Point (D : QuadraticGaugeFiberDatum K)
    (A : Type*) [CommRing A] [Algebra K A] :=
  GaugeFiberPoint D.E D.β D.pi D.b D.a A

/-- Root form of the fiber theorem over a test algebra. -/
def rootEquivPoint (D : QuadraticGaugeFiberDatum K)
    (A : Type*) [CommRing A] [Algebra K A] :
    PolynomialRoot D.E A ≃ D.Point A :=
  rootEquivGaugeFiberPoint D.a D.separable D.g₁ D.derivative_eq

/-- The full functor-of-points theorem: the source fiber is represented by
`K[S]/(E)`. -/
def representingEquiv (D : QuadraticGaugeFiberDatum K)
    (A : Type*) [CommRing A] [Algebra K A] :
    (AdjoinRoot D.E →ₐ[K] A) ≃ D.Point A :=
  (PolynomialRoot.algHomEquiv D.E A).trans (D.rootEquivPoint A)

/-- The quotient universal property is natural under postcomposition. -/
theorem rootOfAlgHom_natural (D : QuadraticGaugeFiberDatum K)
    (f : A →ₐ[K] B) (φ : AdjoinRoot D.E →ₐ[K] A) :
    (PolynomialRoot.ofAlgHom φ).map f =
      PolynomialRoot.ofAlgHom (f.comp φ) := by
  apply PolynomialRoot.ext
  rfl

/-- Naturality of the represented fiber equivalence. -/
theorem representingEquiv_natural (D : QuadraticGaugeFiberDatum K)
    (f : A →ₐ[K] B) (φ : AdjoinRoot D.E →ₐ[K] A) :
    GaugeFiberPoint.map f (D.representingEquiv A φ) =
      D.representingEquiv B (f.comp φ) := by
  change
    GaugeFiberPoint.map f
        (D.rootEquivPoint A (PolynomialRoot.ofAlgHom φ)) =
      D.rootEquivPoint B (PolynomialRoot.ofAlgHom (f.comp φ))
  rw [rootEquivGaugeFiberPoint_natural]
  rw [D.rootOfAlgHom_natural f φ]

/-- Construct the broad datum from the characteristic-zero squarefree
hypothesis used in the paper. -/
def ofSquarefree [CharZero K]
    (E β : K[X]) (pi b a : K) (g₁ : Kˣ)
    (hE : Squarefree E)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β) :
    QuadraticGaugeFiberDatum K where
  E := E
  β := β
  pi := pi
  b := b
  a := a
  g₁ := g₁
  separable := (PerfectField.separable_iff_squarefree).2 hE
  derivative_eq := hderiv

end QuadraticGaugeFiberDatum

end FiniteEtaleKeller
