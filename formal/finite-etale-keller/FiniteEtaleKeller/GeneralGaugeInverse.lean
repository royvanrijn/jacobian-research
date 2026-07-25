/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeMap
import FiniteEtaleKeller.UniversalFiber

/-!
# The generic inverse polynomial of the all-degree quadratic gauge

This module formalizes the paper's polynomials `G_π`, `β(π,S)`, and
`E_{π,b,c}` for an arbitrary finite seed.  It proves the exact derivative
factorization required by the represented-fiber theorem, including the entire
finite high-degree sum.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

/-- The inverse polynomial before imposing the `B,C` target coordinates. -/
def generalGaugeSeedPolynomial (G : K[X]) (pi : K) : K[X] :=
  C (G.coeff 1) * X +
    C pi * (C (G.coeff 2) * X ^ 2 + C (G.coeff 3) * X ^ 3) +
    ∑ k ∈ Finset.Icc 4 G.natDegree,
      C (G.coeff k * pi ^ k) * X ^ k

/-- The marked polynomial
`β(π,S) = (G_π'(S)/g₁ - 1 - πS²)/S`, written without polynomial division. -/
def generalGaugeBeta (G : K[X]) (pi : K) : K[X] :=
  C (2 * (G.coeff 2 / G.coeff 1) * pi) +
    C ((3 * (G.coeff 3 / G.coeff 1) - 1) * pi) * X +
    ∑ k ∈ Finset.Icc 4 G.natDegree,
      C ((k : K) * (G.coeff k / G.coeff 1) * pi ^ k) * X ^ (k - 2)

/-- The generic inverse equation from the paper. -/
def generalGaugeInversePolynomial
    (G : K[X]) (pi b c : K) : K[X] :=
  generalGaugeSeedPolynomial G pi -
    C (G.coeff 1 / 2) * (C b * X ^ 2 + C c)

private theorem generalGauge_tail_derivative
    (G : K[X]) (pi : K) (h₁ : G.coeff 1 ≠ 0) :
    (∑ k ∈ Finset.Icc 4 G.natDegree,
      C (G.coeff k * pi ^ k) * X ^ k).derivative =
      C (G.coeff 1) * X *
        (∑ k ∈ Finset.Icc 4 G.natDegree,
          C ((k : K) * (G.coeff k / G.coeff 1) * pi ^ k) * X ^ (k - 2)) := by
  rw [Polynomial.derivative_sum]
  rw [Finset.mul_sum]
  apply Finset.sum_congr rfl
  intro k hk
  have hk4 : 4 ≤ k := (Finset.mem_Icc.mp hk).1
  have hpow :
      (X : K[X]) * X ^ (k - 2) = X ^ (k - 1) := by
    calc
      (X : K[X]) * X ^ (k - 2) = X ^ 1 * X ^ (k - 2) := by simp
      _ = X ^ (1 + (k - 2)) := (pow_add X 1 (k - 2)).symm
      _ = X ^ (k - 1) := by congr; omega
  have hcoeff :
      G.coeff k * pi ^ k * (k : K) =
        G.coeff 1 * ((k : K) * (G.coeff k / G.coeff 1) * pi ^ k) := by
    field_simp [h₁]
    ring
  rw [Polynomial.derivative_C_mul_X_pow, ← hpow, hcoeff]
  simp only [C_mul]
  ring

private theorem generalGauge_low_derivative
    (G : K[X]) (pi : K) (h₁ : G.coeff 1 ≠ 0) :
    (C (G.coeff 1) * X +
      C pi * (C (G.coeff 2) * X ^ 2 + C (G.coeff 3) * X ^ 3)).derivative =
      C (G.coeff 1) *
        (1 + X *
          (C (2 * (G.coeff 2 / G.coeff 1) * pi) +
            C ((3 * (G.coeff 3 / G.coeff 1) - 1) * pi) * X) +
          C pi * X ^ 2) := by
  simp only [Polynomial.derivative_add, Polynomial.derivative_C_mul,
    Polynomial.derivative_X, Polynomial.derivative_X_pow,
    Polynomial.derivative_C, zero_mul, add_zero, mul_one]
  field_simp [h₁]
  ring_nf

/-- The explicit `β` has exactly the normalized-derivative relation stated in
the paper. -/
theorem generalGaugeSeedPolynomial_derivative
    (G : K[X]) (pi : K) (h₁ : G.coeff 1 ≠ 0) :
    (generalGaugeSeedPolynomial G pi).derivative =
      C (G.coeff 1) *
        (1 + X * generalGaugeBeta G pi + C pi * X ^ 2) := by
  rw [generalGaugeSeedPolynomial, generalGaugeBeta]
  rw [Polynomial.derivative_add]
  rw [generalGauge_low_derivative G pi h₁]
  rw [generalGauge_tail_derivative G pi h₁]
  ring

/-- The generic inverse polynomial has the derivative factorization required
by the universal represented-fiber theorem. -/
theorem generalGaugeInversePolynomial_derivative [CharZero K]
    (G : K[X]) (pi b c : K) (h₁ : G.coeff 1 ≠ 0) :
    (generalGaugeInversePolynomial G pi b c).derivative =
      C (G.coeff 1) * markedChartPolynomial pi b (generalGaugeBeta G pi) := by
  rw [generalGaugeInversePolynomial, Polynomial.derivative_sub,
    generalGaugeSeedPolynomial_derivative G pi h₁]
  simp only [Polynomial.derivative_C_mul, Polynomial.derivative_add,
    Polynomial.derivative_C, Polynomial.derivative_X_pow, zero_mul, add_zero]
  rw [markedChartPolynomial]
  field_simp [h₁]
  ring_nf

section RepresentedFiber

variable [CharZero K]

/-- The abstract source-fiber datum attached directly to the displayed generic
inverse polynomial. -/
def generalGaugeDatum
    (G : K[X]) (pi b c : K)
    (h₁ : G.coeff 1 ≠ 0)
    (hE : (generalGaugeInversePolynomial G pi b c).Separable) :
    QuadraticGaugeFiberDatum K where
  E := generalGaugeInversePolynomial G pi b c
  β := generalGaugeBeta G pi
  pi := pi
  b := b
  a := G.coeff 1 / G.coeff 3
  g₁ := Units.mk0 (G.coeff 1) h₁
  separable := hE
  derivative_eq := by
    change (generalGaugeInversePolynomial G pi b c).derivative =
      C (G.coeff 1) * markedChartPolynomial pi b (generalGaugeBeta G pi)
    exact generalGaugeInversePolynomial_derivative G pi b c h₁

/-- The generic inverse quotient naturally represents the abstract complete
source fiber over every commutative test algebra. -/
def generalGaugeRepresentingEquiv
    (G : K[X]) (pi b c : K)
    (h₁ : G.coeff 1 ≠ 0)
    (hE : (generalGaugeInversePolynomial G pi b c).Separable)
    (A : Type*) [CommRing A] [Algebra K A] :
    (AdjoinRoot (generalGaugeInversePolynomial G pi b c) →ₐ[K] A) ≃
      (generalGaugeDatum G pi b c h₁ hE).Point A :=
  (generalGaugeDatum G pi b c h₁ hE).representingEquiv A

#print axioms generalGaugeInversePolynomial_derivative
#print axioms generalGaugeRepresentingEquiv

end RepresentedFiber

end FiniteEtaleKeller
