/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeMap

/-!
# Jacobian determinant of the all-degree quadratic gauge

This module proves the constant determinant directly for the single general
`MvPolynomial` map.  The proof expands the three-by-three determinant, all
partial derivatives, and the complete finite coefficient sums before invoking
ring normalization.
-/

noncomputable section

open Matrix Function
open MvPolynomial
open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

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
  field_simp [h₁, h₃] <;> ring

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
