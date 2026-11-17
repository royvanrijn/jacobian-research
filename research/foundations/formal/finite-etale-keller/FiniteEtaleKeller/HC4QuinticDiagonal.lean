/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib

/-!
# The scalar core of the diagonal HC(4) quintic obstruction

This file formalizes the ring identities and the characteristic-zero
fifth-power conclusion used by `HC4QD1`.  The extraction of these scalar
coefficients from the generic Hessian determinant is certified separately by
`scripts/verify_hc4_quintic_diagonal_schur.py`.
-/

namespace FiniteEtaleKeller.HC4QuinticDiagonal

variable {K : Type*} [Field K] [CharZero K]

omit [CharZero K] in
theorem diagonalSchurNorm
    (a b c x y z : K) :
    16 *
        (a ^ 2 * x ^ 6 * y ^ 4 * z ^ 4
          + b ^ 2 * x ^ 4 * y ^ 6 * z ^ 4
          + c ^ 2 * x ^ 4 * y ^ 4 * z ^ 6) =
      (x ^ 4 * y ^ 4 * z ^ 4) *
        (16 * (a ^ 2 * x ^ 2 + b ^ 2 * y ^ 2 + c ^ 2 * z ^ 2)) := by
  ring

theorem diagonalFacesForceZero
    (a b c δ : K)
    (face13 :
      -2 * (32 * a ^ 3 + 32 * b ^ 3 + 32 * c ^ 3 - 3 * δ) = 0)
    (face11x :
      -32 * a ^ 2 * (32 * b ^ 3 + 32 * c ^ 3 - 3 * δ) = 0)
    (face11y :
      -32 * b ^ 2 * (32 * a ^ 3 + 32 * c ^ 3 - 3 * δ) = 0)
    (face11z :
      -32 * c ^ 2 * (32 * a ^ 3 + 32 * b ^ 3 - 3 * δ) = 0) :
    a = 0 ∧ b = 0 ∧ c = 0 := by
  have ha5 : (1024 : K) * a ^ 5 = 0 := by
    linear_combination face11x + (-16 * a ^ 2) * face13
  have hb5 : (1024 : K) * b ^ 5 = 0 := by
    linear_combination face11y + (-16 * b ^ 2) * face13
  have hc5 : (1024 : K) * c ^ 5 = 0 := by
    linear_combination face11z + (-16 * c ^ 2) * face13
  have haPow : a ^ 5 = 0 :=
    (mul_eq_zero.mp ha5).resolve_left (by norm_num)
  have hbPow : b ^ 5 = 0 :=
    (mul_eq_zero.mp hb5).resolve_left (by norm_num)
  have hcPow : c ^ 5 = 0 :=
    (mul_eq_zero.mp hc5).resolve_left (by norm_num)
  exact
    ⟨(pow_eq_zero_iff (by norm_num : (5 : ℕ) ≠ 0)).mp haPow,
      (pow_eq_zero_iff (by norm_num : (5 : ℕ) ≠ 0)).mp hbPow,
      (pow_eq_zero_iff (by norm_num : (5 : ℕ) ≠ 0)).mp hcPow⟩

omit [CharZero K] in
theorem coefficientFamilyZeroOfFourthPowers
    {n : ℕ} (coefficients : Fin n → K)
    (fourthPowers : ∀ i, coefficients i ^ 4 = 0) :
    coefficients = 0 := by
  funext i
  exact
    (pow_eq_zero_iff (by norm_num : (4 : ℕ) ≠ 0)).mp
      (fourthPowers i)

omit [CharZero K] in
theorem coefficientFamilyZeroOfCubes
    {n : ℕ} (coefficients : Fin n → K)
    (cubes : ∀ i, coefficients i ^ 3 = 0) :
    coefficients = 0 := by
  funext i
  exact
    (pow_eq_zero_iff (by norm_num : (3 : ℕ) ≠ 0)).mp
      (cubes i)

end FiniteEtaleKeller.HC4QuinticDiagonal
