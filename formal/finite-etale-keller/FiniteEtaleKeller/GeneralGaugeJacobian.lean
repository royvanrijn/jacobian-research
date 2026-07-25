/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeMap

/-!
# Jacobian determinant of the all-degree quadratic gauge

This module proves the constant determinant directly for the single general
`MvPolynomial` map.  The arbitrary high-degree part is first packaged as a
univariate tail in the recurrent polynomial `u = x*q`; its formal derivative
relations are then separated from the universal determinant cancellation.
-/

noncomputable section

open Matrix Function
open MvPolynomial
open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

private def quadraticGaugeQ (a : K) : GaugePolynomial K :=
  generalGaugeT ^ 2 * MvPolynomial.X 2 +
    MvPolynomial.C a * MvPolynomial.X 1 ^ 2 *
      (1 + MvPolynomial.C 3 * generalGaugeT)

private def quadraticGaugePi (a : K) : GaugePolynomial K :=
  generalGaugeT * quadraticGaugeQ a

private def quadraticGaugeU (a : K) : GaugePolynomial K :=
  MvPolynomial.X 0 * quadraticGaugeQ a

private def quadraticGaugeWithTail
    (a c : K) (R Rp : GaugePolynomial K) : Fin 3 → GaugePolynomial K :=
  ![
    quadraticGaugePi a,
    MvPolynomial.X 1 +
      MvPolynomial.C (3 * a⁻¹) * quadraticGaugeU a +
      MvPolynomial.C (2 * c) * quadraticGaugePi a +
      quadraticGaugePi a ^ 2 *
        (MvPolynomial.C 2 * R + quadraticGaugeU a * Rp),
    MvPolynomial.X 0 * (MvPolynomial.C 5 - MvPolynomial.C 3 * generalGaugeT) -
      MvPolynomial.C a⁻¹ * MvPolynomial.X 0 ^ 3 * MvPolynomial.X 2 -
      quadraticGaugeU a ^ 3 * Rp]

set_option maxHeartbeats 0 in
private theorem jacobianDet_quadraticGaugeWithTail
    (a c : K) (ha : a ≠ 0)
    (R Rp Rpp : GaugePolynomial K)
    (hR : ∀ i, pderiv i R = Rp * pderiv i (quadraticGaugeU a))
    (hRp : ∀ i, pderiv i Rp = Rpp * pderiv i (quadraticGaugeU a)) :
    jacobianDet (quadraticGaugeWithTail a c R Rp) = MvPolynomial.C (-2) := by
  classical
  have hR0 := hR (0 : Fin 3)
  have hR1 := hR (1 : Fin 3)
  have hR2 := hR (2 : Fin 3)
  have hRp0 := hRp (0 : Fin 3)
  have hRp1 := hRp (1 : Fin 3)
  have hRp2 := hRp (2 : Fin 3)
  simp only [quadraticGaugeU, quadraticGaugeQ, generalGaugeT,
    map_add, Derivation.map_one_eq_zero, pderiv_mul, pderiv_pow,
    pderiv_C, pderiv_X_self, pderiv_X_of_ne, ne_eq, Fin.reduceEq,
    not_false_eq_true, map_ofNat] at hR0 hR1 hR2 hRp0 hRp1 hRp2
  simp only [jacobianDet, jacobianMatrix, det_fin_three, of_apply,
    quadraticGaugeWithTail, quadraticGaugePi, quadraticGaugeU, quadraticGaugeQ,
    generalGaugeT,
    cons_val_zero, cons_val_one, cons_val_two, head_cons, tail_cons,
    map_add, map_sub, Derivation.map_one_eq_zero,
    pderiv_mul, pderiv_pow, pderiv_C, pderiv_X_self, pderiv_X_of_ne,
    ne_eq, Fin.reduceEq, not_false_eq_true, map_neg, map_ofNat]
  rw [hR0, hR1, hR2, hRp0, hRp1, hRp2]
  field_simp [ha]
  ring

set_option maxHeartbeats 0 in
/-- The actual all-degree quadratic-gauge map has constant Jacobian `-2` under
the two nonvanishing hypotheses used by the construction. -/
theorem jacobianDet_generalGaugeMap
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    jacobianDet (generalGaugeMap G) = MvPolynomial.C (-2) := by
  classical
  simp only [jacobianDet, jacobianMatrix, det_fin_three, of_apply,
    generalGaugeMap, generalGaugePi, generalGaugeB, generalGaugeC,
    generalGaugeT, generalGaugeQ,
    cons_val_zero, cons_val_one, cons_val_two, head_cons, tail_cons,
    map_add, map_sub, map_sum, Derivation.map_one_eq_zero,
    pderiv_mul, pderiv_pow, pderiv_C, pderiv_X_self, pderiv_X_of_ne,
    ne_eq, Fin.reduceEq, not_false_eq_true]
  simp only [map_neg, map_ofNat]
  field_simp [h₁, h₃] <;> ring_nf

/-- The fixed target-preserving output normalization has Jacobian one. -/
theorem jacobianDet_generalGaugeJacobianOneMap [CharZero K]
    (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    jacobianDet (generalGaugeJacobianOneMap G) = 1 := by
  rw [generalGaugeJacobianOneMap, jacobianDet_scaleOutput,
    jacobianDet_generalGaugeMap G h₁ h₃]
  rw [← MvPolynomial.C_mul]
  norm_num

#print axioms jacobianDet_generalGaugeMap
#print axioms jacobianDet_generalGaugeJacobianOneMap

end FiniteEtaleKeller
