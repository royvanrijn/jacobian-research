/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib
import FiniteEtaleKeller.Bezout

/-!
# The explicit optimal quintic polynomial certificate

This lightweight module certifies the univariate part of the smallest Hasse
fiber: expansion, inverse-polynomial identity, derivative, a constructive
Bézout identity, and the derivative inverse in the quotient algebra.
-/

noncomputable section

namespace FiniteEtaleKeller.ExplicitQuintic

/-- The classical minimal intersective quintic. -/
def p5 : Polynomial ℚ :=
  (Polynomial.X ^ 3 - 19) * (Polynomial.X ^ 2 + Polynomial.X + 1)

/-- Expanded form of the classical quintic. -/
theorem p5_expanded :
    p5 = Polynomial.X ^ 5 + Polynomial.X ^ 4 + Polynomial.X ^ 3
      - 19 * Polynomial.X ^ 2 - 19 * Polynomial.X - 19 := by
  unfold p5
  ring

/-- The rooted quadratic-gauge seed. -/
def g5 : Polynomial ℚ :=
  Polynomial.X ^ 5 + Polynomial.X ^ 4 + Polynomial.X ^ 3
    - 19 * Polynomial.X ^ 2 - 19 * Polynomial.X

/-- At normalized target `C = -2`, the inverse polynomial is exactly `p5`. -/
theorem inversePolynomial_eq_p5 : g5 - 19 = p5 := by
  rw [p5_expanded]
  unfold g5
  rfl

/-- Explicit derivative of the quintic. -/
theorem p5_derivative :
    p5.derivative =
      5 * Polynomial.X ^ 4 + 4 * Polynomial.X ^ 3
        + 3 * Polynomial.X ^ 2 - 38 * Polynomial.X - 19 := by
  rw [p5_expanded]
  simp
  have h2 : Polynomial.C (2 : ℚ) = (2 : Polynomial ℚ) := by norm_num
  have h3 : Polynomial.C (3 : ℚ) = (3 : Polynomial ℚ) := by norm_num
  have h4 : Polynomial.C (4 : ℚ) = (4 : Polynomial ℚ) := by norm_num
  rw [h2, h3, h4]
  ring

/-- Integral numerator of the first Bézout coefficient. -/
def bezoutUInt : Polynomial ℚ :=
  100 * Polynomial.X ^ 3 + 30 * Polynomial.X ^ 2
    - 75 * Polynomial.X - 676

/-- Integral numerator of the inverse derivative class. -/
def bezoutVInt : Polynomial ℚ :=
  -20 * Polynomial.X ^ 4 - 10 * Polynomial.X ^ 3
    + 9 * Polynomial.X ^ 2 + 371 * Polynomial.X + 190

/-- Denominator-cleared Bézout identity. -/
theorem p5_bezout_integral :
    bezoutUInt * p5 + bezoutVInt * p5.derivative = 9234 := by
  rw [p5_derivative, p5_expanded]
  unfold bezoutUInt bezoutVInt
  ring

/-- First coefficient in the normalized Bézout identity. -/
def bezoutU : Polynomial ℚ :=
  Polynomial.C (1 / 9234) * bezoutUInt

/-- Explicit inverse of the derivative class modulo `p5`. -/
def bezoutV : Polynomial ℚ :=
  Polynomial.C (1 / 9234) * bezoutVInt

/-- Constructive squarefreeness certificate for the quintic. -/
theorem p5_bezout : bezoutU * p5 + bezoutV * p5.derivative = 1 := by
  calc
    bezoutU * p5 + bezoutV * p5.derivative
        = Polynomial.C (1 / 9234) *
            (bezoutUInt * p5 + bezoutVInt * p5.derivative) := by
              unfold bezoutU bezoutV
              ring
    _ = Polynomial.C (1 / 9234) * (9234 : Polynomial ℚ) := by
          rw [p5_bezout_integral]
    _ = Polynomial.C ((1 / 9234 : ℚ) * 9234) := by
          change Polynomial.C (1 / 9234) * Polynomial.C (9234 : ℚ) =
            Polynomial.C ((1 / 9234 : ℚ) * 9234)
          exact Polynomial.C_mul.symm
    _ = 1 := by norm_num

/-- The class of `p5'` has the displayed inverse in `ℚ[X]/(p5)`. -/
theorem p5_derivative_inverse :
    AdjoinRoot.mk p5 bezoutV * AdjoinRoot.mk p5 p5.derivative = 1 :=
  adjoinRoot_derivative_inverse p5 bezoutU bezoutV p5_bezout

/-- The normalized target `(1,0,-2)` scales to the integral target `(1,0,-38)`. -/
theorem target_scaling : (19 : ℚ) * (-2) = -38 := by
  norm_num

#print axioms inversePolynomial_eq_p5
#print axioms p5_derivative_inverse

end FiniteEtaleKeller.ExplicitQuintic
