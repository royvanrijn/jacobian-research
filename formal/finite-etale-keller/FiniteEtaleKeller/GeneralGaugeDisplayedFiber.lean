/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeInverse
import FiniteEtaleKeller.GaugeInverseAssembly

/-!
# The displayed all-degree map and the represented source fiber

This module connects the actual `MvPolynomial` coordinates to the abstract
source-fiber datum.  The first coordinate is already one of the defining
relations of `GaugeSource`; the two theorems below identify the evaluated
second and third coordinates with the marked equations used by the universal
represented-fiber theorem.
-/

noncomputable section

open Polynomial
open MvPolynomial

namespace FiniteEtaleKeller

variable {K A : Type*} [Field K] [CharZero K]
variable [CommRing A] [Algebra K A]

namespace GaugeSource

/-- The underlying source triple of a globally charted source point. -/
def point {pi a : A} (p : GaugeSource A pi a) : Fin 3 → A :=
  ![p.x, p.y, p.z]

@[simp]
theorem point_zero {pi a : A} (p : GaugeSource A pi a) : p.point 0 = p.x := rfl

@[simp]
theorem point_one {pi a : A} (p : GaugeSource A pi a) : p.point 1 = p.y := rfl

@[simp]
theorem point_two {pi a : A} (p : GaugeSource A pi a) : p.point 2 = p.z := rfl

end GaugeSource

@[simp]
theorem aeval_generalGaugeBeta (G : K[X]) (pi : K) (s : A) :
    Polynomial.aeval s (generalGaugeBeta G pi) =
      algebraMap K A (2 * (G.coeff 2 / G.coeff 1) * pi) +
        algebraMap K A ((3 * (G.coeff 3 / G.coeff 1) - 1) * pi) * s +
        ∑ k ∈ Finset.Icc 4 G.natDegree,
          algebraMap K A ((k : K) * (G.coeff k / G.coeff 1) * pi ^ k) *
            s ^ (k - 2) := by
  simp [generalGaugeBeta]

@[simp]
theorem aeval_generalGaugeSeedPolynomial (G : K[X]) (pi : K) (s : A) :
    Polynomial.aeval s (generalGaugeSeedPolynomial G pi) =
      algebraMap K A (G.coeff 1) * s +
        algebraMap K A pi *
          (algebraMap K A (G.coeff 2) * s ^ 2 +
            algebraMap K A (G.coeff 3) * s ^ 3) +
        ∑ k ∈ Finset.Icc 4 G.natDegree,
          algebraMap K A (G.coeff k * pi ^ k) * s ^ k := by
  simp [generalGaugeSeedPolynomial]

namespace GaugeSource

variable {G : K[X]} {pi : K}

/-- The recurrent polynomial `q` evaluates to the source's derived `q`. -/
@[simp]
theorem eval₂_generalGaugeQ
    (p : GaugeSource A (algebraMap K A pi)
      (algebraMap K A (G.coeff 1 / G.coeff 3))) :
    MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeQ G) = p.q := by
  simp [point, generalGaugeQ, generalGaugeT, GaugeSource.q, p.t_eq]

/-- The first displayed coordinate is exactly the first source-fiber equation. -/
@[simp]
theorem eval₂_generalGaugePi
    (p : GaugeSource A (algebraMap K A pi)
      (algebraMap K A (G.coeff 1 / G.coeff 3))) :
    MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugePi G) =
      algebraMap K A pi := by
  rw [FiniteEtaleKeller.eval₂_generalGaugePi]
  rw [eval₂_generalGaugeQ]
  rw [← p.t_eq]
  exact p.t_mul_q

/-- On the global source chart, the evaluated second displayed coordinate is
exactly `Q + β(π,S)`. -/
theorem eval₂_generalGaugeB_eq_marked
    (p : GaugeSource A (algebraMap K A pi)
      (algebraMap K A (G.coeff 1 / G.coeff 3))) :
    MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeB G) =
      p.Q + Polynomial.aeval p.S (generalGaugeBeta G pi) := by
  let ι : K →+* A := algebraMap K A
  have h := quadraticGauge_fullBIdentity
    p.t p.x p.y p.q
    (ι (G.coeff 2 / G.coeff 1))
    (ι (G.coeff 3 / G.coeff 1))
    (fun k => ι ((k : K) * (G.coeff k / G.coeff 1)))
    G.natDegree
  calc
    MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeB G) =
        p.y + 3 * ι (G.coeff 3 / G.coeff 1) * p.x * p.q +
          2 * ι (G.coeff 2 / G.coeff 1) * (p.t : A) * p.q +
          ∑ k ∈ Finset.Icc 4 G.natDegree,
            ι ((k : K) * (G.coeff k / G.coeff 1)) *
              (p.t : A) ^ 2 * p.x ^ (k - 2) * p.q ^ k := by
      rw [FiniteEtaleKeller.eval₂_generalGaugeB]
      rw [eval₂_generalGaugeQ]
      simp only [point_zero, point_one]
      rw [← p.t_eq]
      simp [ι, map_mul]
    _ = (p.y + p.x * p.q) +
          2 * ι (G.coeff 2 / G.coeff 1) * ((p.t : A) * p.q) +
          (3 * ι (G.coeff 3 / G.coeff 1) - 1) *
            ((p.t : A) * p.q) * (p.x * (↑p.t⁻¹ : A)) +
          ∑ k ∈ Finset.Icc 4 G.natDegree,
            ι ((k : K) * (G.coeff k / G.coeff 1)) *
              ((((p.t : A) * p.q) ^ k) *
                (p.x * (↑p.t⁻¹ : A)) ^ (k - 2)) := h.symm
    _ = p.Q + Polynomial.aeval p.S (generalGaugeBeta G pi) := by
      simp only [GaugeSource.Q, GaugeSource.S, aeval_generalGaugeBeta]
      rw [p.t_mul_q]
      simp [ι, map_mul, map_sub, map_pow]

/-- On the same chart, the third displayed coordinate is the inverse-equation
expression `2*G_π/g₁ - B*S²`. -/
theorem eval₂_generalGaugeC_eq_inverse
    (p : GaugeSource A (algebraMap K A pi)
      (algebraMap K A (G.coeff 1 / G.coeff 3)))
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeC G) =
      algebraMap K A (2 / G.coeff 1) *
          Polynomial.aeval p.S (generalGaugeSeedPolynomial G pi) -
        MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeB G) * p.S ^ 2 := by
  let ι : K →+* A := algebraMap K A
  have hraK : (G.coeff 3 / G.coeff 1) * (G.coeff 1 / G.coeff 3) = 1 := by
    field_simp [h₁, h₃]
  have hra :
      ι (G.coeff 3 / G.coeff 1) * ι (G.coeff 1 / G.coeff 3) = 1 := by
    rw [← map_mul, hraK, map_one]
  have hdisplay := quadraticGauge_fullCIdentity
    p.t p.x p.y p.z p.q
    (ι (G.coeff 1 / G.coeff 3))
    (ι (G.coeff 3 / G.coeff 1))
    (fun k => ι (((k - 2 : ℕ) : K) * (G.coeff k / G.coeff 1)))
    G.natDegree p.t_eq rfl hra
  have hexpand := quadraticGauge_inverseCExpansion
    p.S p.Q (ι pi)
    (ι (G.coeff 2 / G.coeff 1))
    (ι (G.coeff 3 / G.coeff 1))
    (fun k => ι (G.coeff k / G.coeff 1))
    G.natDegree
  calc
    MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeC G) =
        p.x * (5 - 3 * (p.t : A)) -
          ι (G.coeff 3 / G.coeff 1) * p.x ^ 3 * p.z -
          ∑ k ∈ Finset.Icc 4 G.natDegree,
            ι (((k - 2 : ℕ) : K) * (G.coeff k / G.coeff 1)) *
              (p.x * p.q) ^ k := by
      rw [FiniteEtaleKeller.eval₂_generalGaugeC]
      rw [eval₂_generalGaugeQ]
      simp only [point_zero, point_two]
      rw [← p.t_eq]
      simp [ι, map_mul]
    _ = 2 * p.S - p.Q * p.S ^ 2 +
          (1 - ι (G.coeff 3 / G.coeff 1)) * ι pi * p.S ^ 3 -
          ∑ k ∈ Finset.Icc 4 G.natDegree,
            (((k - 2 : ℕ) : A) * ι (G.coeff k / G.coeff 1)) *
              ι pi ^ k * p.S ^ k := by
      rw [← hdisplay]
      rw [p.t_mul_q]
      simp [ι, map_mul, map_pow]
    _ = 2 *
          (p.S + ι (G.coeff 2 / G.coeff 1) * ι pi * p.S ^ 2 +
            ι (G.coeff 3 / G.coeff 1) * ι pi * p.S ^ 3 +
            ∑ k ∈ Finset.Icc 4 G.natDegree,
              ι (G.coeff k / G.coeff 1) * ι pi ^ k * p.S ^ k) -
          (p.Q + 2 * ι (G.coeff 2 / G.coeff 1) * ι pi +
            (3 * ι (G.coeff 3 / G.coeff 1) - 1) * ι pi * p.S +
            ∑ k ∈ Finset.Icc 4 G.natDegree,
              (k : A) * ι (G.coeff k / G.coeff 1) *
                ι pi ^ k * p.S ^ (k - 2)) * p.S ^ 2 := hexpand.symm
    _ = algebraMap K A (2 / G.coeff 1) *
          Polynomial.aeval p.S (generalGaugeSeedPolynomial G pi) -
        (p.Q + Polynomial.aeval p.S (generalGaugeBeta G pi)) * p.S ^ 2 := by
      simp only [aeval_generalGaugeSeedPolynomial, aeval_generalGaugeBeta]
      simp [ι, map_mul, map_sub, map_pow]
      field_simp [h₁]
      ring
    _ = algebraMap K A (2 / G.coeff 1) *
          Polynomial.aeval p.S (generalGaugeSeedPolynomial G pi) -
        MvPolynomial.eval₂ (algebraMap K A) p.point (generalGaugeB G) * p.S ^ 2 := by
      rw [eval₂_generalGaugeB_eq_marked]

end GaugeSource

/-- Points of the actual displayed map fiber, with the global first-coordinate
chart carried explicitly. -/
@[ext]
structure GeneralGaugeDisplayedFiberPoint
    (G : K[X]) (pi b c : K) (A : Type*) [CommRing A] [Algebra K A] where
  source : GaugeSource A (algebraMap K A pi)
    (algebraMap K A (G.coeff 1 / G.coeff 3))
  b_eq : MvPolynomial.eval₂ (algebraMap K A) source.point (generalGaugeB G) =
    algebraMap K A b
  c_eq : MvPolynomial.eval₂ (algebraMap K A) source.point (generalGaugeC G) =
    algebraMap K A c

namespace GeneralGaugeDisplayedFiberPoint

variable {G : K[X]} {pi b c : K}

/-- A displayed-map fiber point determines the abstract represented source
fiber point. -/
def toGaugeFiberPoint
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (p : GeneralGaugeDisplayedFiberPoint G pi b c A) :
    GaugeFiberPoint
      (generalGaugeInversePolynomial G pi b c)
      (generalGaugeBeta G pi) pi b (G.coeff 1 / G.coeff 3) A where
  source := p.source
  marked_eq := by
    rw [← p.source.eval₂_generalGaugeB_eq_marked]
    exact p.b_eq
  root_eq := by
    have hc := p.source.eval₂_generalGaugeC_eq_inverse h₁ h₃
    rw [p.b_eq, p.c_eq] at hc
    simp only [generalGaugeInversePolynomial, map_sub, map_mul,
      Polynomial.aeval_sub, Polynomial.aeval_mul, Polynomial.aeval_C,
      Polynomial.aeval_add, Polynomial.aeval_X, Polynomial.aeval_pow]
    rw [hc]
    have hscale :
        algebraMap K A (G.coeff 1 / 2) *
          algebraMap K A (2 / G.coeff 1) = 1 := by
      rw [← map_mul]
      field_simp [h₁]
    rw [← map_mul]
    ring_nf at hscale ⊢
    rw [hscale]
    ring

end GeneralGaugeDisplayedFiberPoint

namespace GaugeFiberPoint

variable {G : K[X]} {pi b c : K}

/-- An abstract represented source point satisfies the actual displayed map's
second and third coordinate equations. -/
def toGeneralGaugeDisplayedFiberPoint
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (p : GaugeFiberPoint
      (generalGaugeInversePolynomial G pi b c)
      (generalGaugeBeta G pi) pi b (G.coeff 1 / G.coeff 3) A) :
    GeneralGaugeDisplayedFiberPoint G pi b c A where
  source := p.source
  b_eq := by
    rw [p.source.eval₂_generalGaugeB_eq_marked]
    exact p.marked_eq
  c_eq := by
    have hB :
        MvPolynomial.eval₂ (algebraMap K A) p.source.point (generalGaugeB G) =
          algebraMap K A b := by
      rw [p.source.eval₂_generalGaugeB_eq_marked]
      exact p.marked_eq
    have hC := p.source.eval₂_generalGaugeC_eq_inverse h₁ h₃
    rw [hB]
    have hroot := p.root_eq
    simp only [generalGaugeInversePolynomial, map_sub, map_mul,
      Polynomial.aeval_sub, Polynomial.aeval_mul, Polynomial.aeval_C,
      Polynomial.aeval_add, Polynomial.aeval_X, Polynomial.aeval_pow] at hroot
    have hscale :
        algebraMap K A (2 / G.coeff 1) *
          algebraMap K A (G.coeff 1 / 2) = 1 := by
      rw [← map_mul]
      field_simp [h₁]
    have hseed :
        algebraMap K A (2 / G.coeff 1) *
            Polynomial.aeval p.source.S (generalGaugeSeedPolynomial G pi) =
          algebraMap K A b * p.source.S ^ 2 + algebraMap K A c := by
      calc
        algebraMap K A (2 / G.coeff 1) *
            Polynomial.aeval p.source.S (generalGaugeSeedPolynomial G pi) =
          algebraMap K A (2 / G.coeff 1) *
            (algebraMap K A (G.coeff 1 / 2) *
              (algebraMap K A b * p.source.S ^ 2 + algebraMap K A c)) := by
                congr 1
                exact sub_eq_zero.mp hroot
        _ = algebraMap K A b * p.source.S ^ 2 + algebraMap K A c := by
          rw [← mul_assoc, hscale, one_mul]
    rw [hC, hseed]
    ring

end GaugeFiberPoint

/-- The actual displayed source-fiber equations and the abstract represented
source-fiber datum are equivalent over every commutative test algebra. -/
def generalGaugeDisplayedFiberEquiv
    {G : K[X]} {pi b c : K}
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    GeneralGaugeDisplayedFiberPoint G pi b c A ≃
      GaugeFiberPoint
        (generalGaugeInversePolynomial G pi b c)
        (generalGaugeBeta G pi) pi b (G.coeff 1 / G.coeff 3) A where
  toFun := GeneralGaugeDisplayedFiberPoint.toGaugeFiberPoint h₁ h₃
  invFun := GaugeFiberPoint.toGeneralGaugeDisplayedFiberPoint h₁ h₃
  left_inv := by
    intro p
    apply GeneralGaugeDisplayedFiberPoint.ext
    rfl
  right_inv := by
    intro p
    apply GaugeFiberPoint.ext
    rfl

#print axioms GaugeSource.eval₂_generalGaugeB_eq_marked
#print axioms GaugeSource.eval₂_generalGaugeC_eq_inverse
#print axioms generalGaugeDisplayedFiberEquiv

end FiniteEtaleKeller
