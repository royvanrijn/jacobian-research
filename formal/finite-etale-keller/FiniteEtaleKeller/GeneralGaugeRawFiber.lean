/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeDisplayedFiber

/-!
# The literal fiber of the all-degree polynomial map

The previous module identifies the displayed coordinate equations after a
source-chart unit has been supplied.  Here a unit first target coordinate is
used to construct that chart directly from a raw triple satisfying the three
`MvPolynomial` equations.  This produces the literal functor-of-points fiber
of the polynomial map and removes the final auxiliary source-chart wrapper.
-/

noncomputable section

open Polynomial
open MvPolynomial

namespace FiniteEtaleKeller

variable {K A B : Type*} [Field K] [CharZero K]
variable [CommRing A] [Algebra K A]
variable [CommRing B] [Algebra K B]

/-- A point of the literal fiber of the all-degree quadratic-gauge map.  The
first target coordinate is stored as a unit because this is exactly what makes
the reciprocal source chart global on the whole fiber. -/
@[ext]
structure GeneralGaugeRawFiberPoint
    (G : K[X]) (pi : Kˣ) (b c : K)
    (A : Type*) [CommRing A] [Algebra K A] where
  point : Fin 3 → A
  pi_eq : MvPolynomial.eval₂ (algebraMap K A) point (generalGaugePi G) =
    algebraMap K A (pi : K)
  b_eq : MvPolynomial.eval₂ (algebraMap K A) point (generalGaugeB G) =
    algebraMap K A b
  c_eq : MvPolynomial.eval₂ (algebraMap K A) point (generalGaugeC G) =
    algebraMap K A c

/-- If a product is a supplied scalar unit, then its first factor is a unit.
This is valid over an arbitrary commutative ring and contains no cancellation
or domain hypothesis. -/
def firstFactorUnitOfMulEqScalarUnit
    (t q : A) (pi : Kˣ)
    (h : t * q = algebraMap K A (pi : K)) : Aˣ where
  val := t
  inv := q * algebraMap K A (↑pi⁻¹ : K)
  val_inv := by
    calc
      t * (q * algebraMap K A (↑pi⁻¹ : K)) =
          (t * q) * algebraMap K A (↑pi⁻¹ : K) := by ring
      _ = algebraMap K A (pi : K) * algebraMap K A (↑pi⁻¹ : K) := by rw [h]
      _ = algebraMap K A ((pi : K) * (↑pi⁻¹ : K)) := by rw [map_mul]
      _ = 1 := by simp
  inv_val := by
    calc
      (q * algebraMap K A (↑pi⁻¹ : K)) * t =
          (t * q) * algebraMap K A (↑pi⁻¹ : K) := by ring
      _ = algebraMap K A (pi : K) * algebraMap K A (↑pi⁻¹ : K) := by rw [h]
      _ = algebraMap K A ((pi : K) * (↑pi⁻¹ : K)) := by rw [map_mul]
      _ = 1 := by simp

namespace GeneralGaugeRawFiberPoint

variable {G : K[X]} {pi : Kˣ} {b c : K}

/-- The recurrent values attached to a raw source triple. -/
def t (p : GeneralGaugeRawFiberPoint G pi b c A) : A :=
  1 + p.point 0 * p.point 1

def q (p : GeneralGaugeRawFiberPoint G pi b c A) : A :=
  MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeQ G)

/-- The first coordinate equation is precisely `t*q=pi`. -/
theorem t_mul_q (p : GeneralGaugeRawFiberPoint G pi b c A) :
    p.t * p.q = algebraMap K A (pi : K) := by
  simpa [t, q, generalGaugePi] using p.pi_eq

/-- The raw first coordinate equation supplies the global source-chart unit. -/
def tUnit (p : GeneralGaugeRawFiberPoint G pi b c A) : Aˣ :=
  firstFactorUnitOfMulEqScalarUnit p.t p.q pi p.t_mul_q

@[simp]
theorem tUnit_val (p : GeneralGaugeRawFiberPoint G pi b c A) :
    (p.tUnit : A) = p.t := rfl

/-- Convert a literal polynomial-map fiber point into the displayed source
fiber with its reciprocal chart. -/
def toDisplayed (p : GeneralGaugeRawFiberPoint G pi b c A) :
    GeneralGaugeDisplayedFiberPoint G (pi : K) b c A where
  source :=
    { t := p.tUnit
      x := p.point 0
      y := p.point 1
      z := p.point 2
      t_eq := rfl
      pi_eq := by
        change p.t *
          (p.t ^ 2 * p.point 2 +
            algebraMap K A (G.coeff 1 / G.coeff 3) * p.point 1 ^ 2 *
              (1 + 3 * p.t)) = algebraMap K A (pi : K)
        have hq :
            p.q = p.t ^ 2 * p.point 2 +
              algebraMap K A (G.coeff 1 / G.coeff 3) * p.point 1 ^ 2 *
                (1 + 3 * p.t) := by
          change MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeQ G) = _
          simpa only [t] using
            (FiniteEtaleKeller.eval₂_generalGaugeQ (A := A) G p.point)
        calc
          p.t *
              (p.t ^ 2 * p.point 2 +
                algebraMap K A (G.coeff 1 / G.coeff 3) * p.point 1 ^ 2 *
                  (1 + 3 * p.t)) =
            p.t * p.q := congrArg (fun u : A => p.t * u) hq.symm
          _ = algebraMap K A (pi : K) := p.t_mul_q }
  b_eq := p.b_eq
  c_eq := p.c_eq

end GeneralGaugeRawFiberPoint

namespace GeneralGaugeDisplayedFiberPoint

variable {G : K[X]} {pi : Kˣ} {b c : K}

/-- Forget the explicit source-chart witness and retain the literal polynomial
map fiber point. -/
def toRaw (p : GeneralGaugeDisplayedFiberPoint G (pi : K) b c A) :
    GeneralGaugeRawFiberPoint G pi b c A where
  point := p.source.point
  pi_eq := p.source.eval₂_generalGaugePi
  b_eq := p.b_eq
  c_eq := p.c_eq

end GeneralGaugeDisplayedFiberPoint

/-- A unit first target coordinate makes the literal map fiber and the
source-chart displayed fiber equivalent over every commutative test algebra. -/
def generalGaugeRawFiberEquivDisplayed
    {G : K[X]} {pi : Kˣ} {b c : K} :
    GeneralGaugeRawFiberPoint G pi b c A ≃
      GeneralGaugeDisplayedFiberPoint G (pi : K) b c A where
  toFun := GeneralGaugeRawFiberPoint.toDisplayed
  invFun := GeneralGaugeDisplayedFiberPoint.toRaw
  left_inv := by
    intro p
    apply GeneralGaugeRawFiberPoint.ext
    funext i
    fin_cases i <;> rfl
  right_inv := by
    intro p
    apply GeneralGaugeDisplayedFiberPoint.ext
    apply GaugeSource.ext
    · apply Units.ext
      change 1 + p.source.x * p.source.y = (p.source.t : A)
      exact p.source.t_eq.symm
    · rfl
    · rfl
    · rfl

private theorem eval₂_map_algHom
    (f : A →ₐ[K] B) (point : Fin 3 → A) (P : GaugePolynomial K) :
    MvPolynomial.eval₂ (algebraMap K B) (fun i => f (point i)) P =
      f (MvPolynomial.eval₂ (algebraMap K A) point P) := by
  have hcomp : f.toRingHom.comp (algebraMap K A) = algebraMap K B := by
    ext r
    exact f.commutes r
  calc
    MvPolynomial.eval₂ (algebraMap K B) (fun i => f (point i)) P =
        MvPolynomial.eval₂ (f.toRingHom.comp (algebraMap K A))
          (fun i => f (point i)) P := by rw [hcomp]
    _ = f (MvPolynomial.eval₂ (algebraMap K A) point P) :=
      (MvPolynomial.hom_eval₂ P (algebraMap K A) f.toRingHom point).symm

namespace GeneralGaugeRawFiberPoint

variable {G : K[X]} {pi : Kˣ} {b c : K}

/-- Literal polynomial-map fibers are functorial in the commutative test
algebra. -/
def map (f : A →ₐ[K] B)
    (p : GeneralGaugeRawFiberPoint G pi b c A) :
    GeneralGaugeRawFiberPoint G pi b c B where
  point := fun i => f (p.point i)
  pi_eq := by
    calc
      MvPolynomial.eval₂ (algebraMap K B) (fun i => f (p.point i))
          (generalGaugePi G) =
        f (MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugePi G)) :=
          eval₂_map_algHom f p.point (generalGaugePi G)
      _ = f (algebraMap K A (pi : K)) := congrArg f p.pi_eq
      _ = algebraMap K B (pi : K) := f.commutes _
  b_eq := by
    calc
      MvPolynomial.eval₂ (algebraMap K B) (fun i => f (p.point i))
          (generalGaugeB G) =
        f (MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeB G)) :=
          eval₂_map_algHom f p.point (generalGaugeB G)
      _ = f (algebraMap K A b) := congrArg f p.b_eq
      _ = algebraMap K B b := f.commutes _
  c_eq := by
    calc
      MvPolynomial.eval₂ (algebraMap K B) (fun i => f (p.point i))
          (generalGaugeC G) =
        f (MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeC G)) :=
          eval₂_map_algHom f p.point (generalGaugeC G)
      _ = f (algebraMap K A c) := congrArg f p.c_eq
      _ = algebraMap K B c := f.commutes _

@[simp]
theorem map_point (f : A →ₐ[K] B)
    (p : GeneralGaugeRawFiberPoint G pi b c A) (i : Fin 3) :
    (p.map f).point i = f (p.point i) := rfl

end GeneralGaugeRawFiberPoint

/-- The polynomial quotient of the generic inverse equation naturally
represents the literal fiber of the actual all-degree polynomial map. -/
def generalGaugeRawRepresentingEquiv
    (G : K[X]) (pi : Kˣ) (b c : K)
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (hE : (generalGaugeInversePolynomial G (pi : K) b c).Separable)
    (A : Type*) [CommRing A] [Algebra K A] :
    (AdjoinRoot (generalGaugeInversePolynomial G (pi : K) b c) →ₐ[K] A) ≃
      GeneralGaugeRawFiberPoint G pi b c A :=
  (generalGaugeRepresentingEquiv G (pi : K) b c h₁ hE A).trans
    ((generalGaugeDisplayedFiberEquiv (A := A) h₁ h₃).symm.trans
      (generalGaugeRawFiberEquivDisplayed (A := A)).symm)

/-- Naturality of the literal represented fiber equivalence under
postcomposition of test-algebra maps. -/
theorem generalGaugeRawRepresentingEquiv_natural
    (G : K[X]) (pi : Kˣ) (b c : K)
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (hE : (generalGaugeInversePolynomial G (pi : K) b c).Separable)
    (f : A →ₐ[K] B)
    (φ : AdjoinRoot (generalGaugeInversePolynomial G (pi : K) b c) →ₐ[K] A) :
    GeneralGaugeRawFiberPoint.map f
        (generalGaugeRawRepresentingEquiv G pi b c h₁ h₃ hE A φ) =
      generalGaugeRawRepresentingEquiv G pi b c h₁ h₃ hE B (f.comp φ) := by
  apply GeneralGaugeRawFiberPoint.ext
  funext i
  have h :=
    (generalGaugeDatum G (pi : K) b c h₁ hE).representingEquiv_natural f φ
  fin_cases i
  · exact congrArg (fun p => p.source.x) h
  · exact congrArg (fun p => p.source.y) h
  · exact congrArg (fun p => p.source.z) h

#print axioms generalGaugeRawFiberEquivDisplayed
#print axioms generalGaugeRawRepresentingEquiv
#print axioms generalGaugeRawRepresentingEquiv_natural

end FiniteEtaleKeller
