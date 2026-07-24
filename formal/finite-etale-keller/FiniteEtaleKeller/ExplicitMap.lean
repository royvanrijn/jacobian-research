/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.Jacobian

/-!
# The explicit optimal quintic map certificate

This module certifies the three-variable polynomial map, its denominator-free
Jacobian `-722`, its determinant-`-2` quadratic-gauge normalization, and the
fixed determinant-one output normalization.
-/

noncomputable section

open Matrix Function
open MvPolynomial

namespace FiniteEtaleKeller.ExplicitQuintic

abbrev M := MvPolynomial (Fin 3) ℚ

/-- The recurrent source polynomial `t = 1 + xy`. -/
def t : M := 1 + MvPolynomial.X 0 * MvPolynomial.X 1

/-- The recurrent source polynomial used by the integral quintic map. -/
def q : M :=
  t ^ 2 * MvPolynomial.X 2
    - MvPolynomial.C 19 * MvPolynomial.X 1 ^ 2 * (1 + MvPolynomial.C 3 * t)

/-- The denominator-free degree-five Keller map from the paper. -/
def integralMap : Fin 3 → M :=
  ![t * q,
    MvPolynomial.C 19 * MvPolynomial.X 1
      - MvPolynomial.C 3 * MvPolynomial.X 0 * q
      + MvPolynomial.C 38 * t * q
      - MvPolynomial.C 4 * t ^ 2 * MvPolynomial.X 0 ^ 2 * q ^ 4
      - MvPolynomial.C 5 * t ^ 2 * MvPolynomial.X 0 ^ 3 * q ^ 5,
    MvPolynomial.C 19 * MvPolynomial.X 0 * (MvPolynomial.C 5 - MvPolynomial.C 3 * t)
      + MvPolynomial.X 0 ^ 3 * MvPolynomial.X 2
      + MvPolynomial.C 2 * (MvPolynomial.X 0 * q) ^ 4
      + MvPolynomial.C 3 * (MvPolynomial.X 0 * q) ^ 5]

/-- The determinant-`-2` quadratic-gauge normalization. -/
def normalizedMap : Fin 3 → M :=
  scaleOutput (1 : ℚ) (1 / 19 : ℚ) (1 / 19 : ℚ) integralMap

/-- The target-preserving determinant-one normalization from the main theorem. -/
def jacobianOneMap : Fin 3 → M :=
  scaleOutput 1 (-1 / 2 : ℚ) 1 normalizedMap

set_option maxHeartbeats 0 in
-- The direct three-variable determinant expansion is intentionally unrestricted.
/-- The displayed denominator-free map has Jacobian determinant `-722`. -/
theorem jacobianDet_integralMap :
    jacobianDet integralMap = MvPolynomial.C (-722) := by
  simp only [jacobianDet, jacobianMatrix, det_fin_three, of_apply,
    integralMap, t, q, cons_val_zero, cons_val_one, cons_val_two, head_cons,
    tail_cons, map_add, map_sub, Derivation.map_one_eq_zero, pderiv_mul,
    pderiv_pow, pderiv_C, pderiv_X_self, pderiv_X_of_ne, ne_eq,
    Fin.reduceEq, not_false_eq_true]
  simp only [map_neg, map_ofNat]
  ring

/-- The normalized quadratic gauge has Jacobian determinant `-2`. -/
theorem jacobianDet_normalizedMap :
    jacobianDet normalizedMap = MvPolynomial.C (-2) := by
  rw [normalizedMap, jacobianDet_scaleOutput, jacobianDet_integralMap]
  norm_num [MvPolynomial.C_mul']

/-- The universal target-preserving output normalization has determinant `1`. -/
theorem jacobianDet_jacobianOneMap : jacobianDet jacobianOneMap = 1 := by
  rw [jacobianOneMap, jacobianDet_scaleOutput, jacobianDet_normalizedMap]
  norm_num [MvPolynomial.C_mul']

/-- Scaling the normalized gauge back by `diag(1,19,19)` recovers the displayed map. -/
theorem integralMap_eq_scaled_normalized :
    integralMap = scaleOutput (1 : ℚ) 19 19 normalizedMap := by
  funext i
  fin_cases i <;> simp [normalizedMap, scaleOutput, MvPolynomial.C_mul']

#print axioms jacobianDet_integralMap
#print axioms jacobianDet_jacobianOneMap

end FiniteEtaleKeller.ExplicitQuintic
