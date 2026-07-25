/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeRawFiber

/-!
# The determinant-one target normalization on the literal fiber

The main theorem uses the fixed output scaling `diag(1,-1/2,1)`.  At the
chosen target the second coordinate is zero, so the normalized and unnormalized
literal fibers are naturally equivalent over every commutative test algebra.
-/

noncomputable section

open Polynomial
open MvPolynomial

namespace FiniteEtaleKeller

variable {K A B : Type*} [Field K] [CharZero K]
variable [CommRing A] [Algebra K A]
variable [CommRing B] [Algebra K B]

/-- A literal point of the determinant-one normalized map over target
`(pi,0,c)`. -/
@[ext]
structure GeneralGaugeJacobianOneFiberPoint
    (G : K[X]) (pi : Kˣ) (c : K)
    (A : Type*) [CommRing A] [Algebra K A] where
  point : Fin 3 → A
  pi_eq : MvPolynomial.eval₂ (algebraMap K A) point
      (generalGaugeJacobianOneMap G 0) = algebraMap K A (pi : K)
  b_eq : MvPolynomial.eval₂ (algebraMap K A) point
      (generalGaugeJacobianOneMap G 1) = 0
  c_eq : MvPolynomial.eval₂ (algebraMap K A) point
      (generalGaugeJacobianOneMap G 2) = algebraMap K A c

namespace GeneralGaugeRawFiberPoint

variable {G : K[X]} {pi : Kˣ} {c : K}

/-- Apply the target-preserving determinant-one normalization to a raw fiber
with second target coordinate zero. -/
def toJacobianOne
    (p : GeneralGaugeRawFiberPoint G pi 0 c A) :
    GeneralGaugeJacobianOneFiberPoint G pi c A where
  point := p.point
  pi_eq := by
    simpa [generalGaugeJacobianOneMap, scaleOutput, generalGaugeMap] using p.pi_eq
  b_eq := by
    simp [generalGaugeJacobianOneMap, scaleOutput, generalGaugeMap, p.b_eq]
  c_eq := by
    simpa [generalGaugeJacobianOneMap, scaleOutput, generalGaugeMap] using p.c_eq

end GeneralGaugeRawFiberPoint

namespace GeneralGaugeJacobianOneFiberPoint

variable {G : K[X]} {pi : Kˣ} {c : K}

/-- Undo the output normalization on the zero second-coordinate fiber.  The
proof uses the inverse scalar explicitly and therefore remains valid over a
commutative test ring with zero divisors. -/
def toRaw
    (p : GeneralGaugeJacobianOneFiberPoint G pi c A) :
    GeneralGaugeRawFiberPoint G pi 0 c A where
  point := p.point
  pi_eq := by
    simpa [generalGaugeJacobianOneMap, scaleOutput, generalGaugeMap] using p.pi_eq
  b_eq := by
    have hb := p.b_eq
    change algebraMap K A (-1 / 2 : K) *
      MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeB G) = 0 at hb
    have hscale :
        algebraMap K A (-2 : K) * algebraMap K A (-1 / 2 : K) = 1 := by
      rw [← map_mul]
      norm_num
    calc
      MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeB G) =
          1 * MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeB G) := by
            rw [one_mul]
      _ = (algebraMap K A (-2 : K) * algebraMap K A (-1 / 2 : K)) *
          MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeB G) := by
            rw [hscale]
      _ = algebraMap K A (-2 : K) *
          (algebraMap K A (-1 / 2 : K) *
            MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeB G)) := by
              ring
      _ = 0 := by rw [hb]; ring
  c_eq := by
    simpa [generalGaugeJacobianOneMap, scaleOutput, generalGaugeMap] using p.c_eq

/-- Normalized literal fibers are functorial in the test algebra. -/
def map (f : A →ₐ[K] B)
    (p : GeneralGaugeJacobianOneFiberPoint G pi c A) :
    GeneralGaugeJacobianOneFiberPoint G pi c B where
  point := fun i => f (p.point i)
  pi_eq := by
    have h := congrArg f p.pi_eq
    rw [MvPolynomial.hom_eval₂] at h
    simpa [RingHom.comp_apply] using h
  b_eq := by
    have h := congrArg f p.b_eq
    rw [MvPolynomial.hom_eval₂] at h
    simpa [RingHom.comp_apply] using h
  c_eq := by
    have h := congrArg f p.c_eq
    rw [MvPolynomial.hom_eval₂] at h
    simpa [RingHom.comp_apply] using h

end GeneralGaugeJacobianOneFiberPoint

/-- The zero second-coordinate fiber is unchanged by the determinant-one
output normalization. -/
def generalGaugeRawFiberEquivJacobianOne
    {G : K[X]} {pi : Kˣ} {c : K} :
    GeneralGaugeRawFiberPoint G pi 0 c A ≃
      GeneralGaugeJacobianOneFiberPoint G pi c A where
  toFun := GeneralGaugeRawFiberPoint.toJacobianOne
  invFun := GeneralGaugeJacobianOneFiberPoint.toRaw
  left_inv := by
    intro p
    apply GeneralGaugeRawFiberPoint.ext
    rfl
  right_inv := by
    intro p
    apply GeneralGaugeJacobianOneFiberPoint.ext
    rfl

/-- The generic inverse quotient represents the literal determinant-one fiber
at a zero second target coordinate. -/
def generalGaugeJacobianOneRepresentingEquiv
    (G : K[X]) (pi : Kˣ) (c : K)
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (hE : (generalGaugeInversePolynomial G (pi : K) 0 c).Separable)
    (A : Type*) [CommRing A] [Algebra K A] :
    (AdjoinRoot (generalGaugeInversePolynomial G (pi : K) 0 c) →ₐ[K] A) ≃
      GeneralGaugeJacobianOneFiberPoint G pi c A :=
  (generalGaugeRawRepresentingEquiv G pi 0 c h₁ h₃ hE A).trans
    (generalGaugeRawFiberEquivJacobianOne (A := A))

/-- Naturality of the determinant-one literal represented fiber. -/
theorem generalGaugeJacobianOneRepresentingEquiv_natural
    (G : K[X]) (pi : Kˣ) (c : K)
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (hE : (generalGaugeInversePolynomial G (pi : K) 0 c).Separable)
    (f : A →ₐ[K] B)
    (φ : AdjoinRoot (generalGaugeInversePolynomial G (pi : K) 0 c) →ₐ[K] A) :
    GeneralGaugeJacobianOneFiberPoint.map f
        (generalGaugeJacobianOneRepresentingEquiv G pi c h₁ h₃ hE A φ) =
      generalGaugeJacobianOneRepresentingEquiv G pi c h₁ h₃ hE B (f.comp φ) := by
  apply GeneralGaugeJacobianOneFiberPoint.ext
  funext i
  have h := generalGaugeRawRepresentingEquiv_natural
    G pi 0 c h₁ h₃ hE f φ
  fin_cases i
  · exact congrArg (fun p => p.point 0) h
  · exact congrArg (fun p => p.point 1) h
  · exact congrArg (fun p => p.point 2) h

#print axioms generalGaugeRawFiberEquivJacobianOne
#print axioms generalGaugeJacobianOneRepresentingEquiv
#print axioms generalGaugeJacobianOneRepresentingEquiv_natural

end FiniteEtaleKeller
