/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.Bezout
import Mathlib.Tactic.NativeDecide

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
  (Polynomial.X ^ 3 - Polynomial.C 19)
    * (Polynomial.X ^ 2 + Polynomial.X + 1)

/-- Expanded form of the classical quintic. -/
theorem p5_expanded :
    p5 = Polynomial.X ^ 5 + Polynomial.X ^ 4 + Polynomial.X ^ 3
      - Polynomial.C 19 * Polynomial.X ^ 2
      - Polynomial.C 19 * Polynomial.X - Polynomial.C 19 := by
  native_decide

/-- The rooted quadratic-gauge seed. -/
def g5 : Polynomial ℚ :=
  Polynomial.X ^ 5 + Polynomial.X ^ 4 + Polynomial.X ^ 3
    - Polynomial.C 19 * Polynomial.X ^ 2
    - Polynomial.C 19 * Polynomial.X

/-- At normalized target `C = -2`, the inverse polynomial is exactly `p5`. -/
theorem inversePolynomial_eq_p5 : g5 - Polynomial.C 19 = p5 := by
  native_decide

/-- Explicit derivative of the quintic. -/
theorem p5_derivative :
    p5.derivative =
      Polynomial.C 5 * Polynomial.X ^ 4
        + Polynomial.C 4 * Polynomial.X ^ 3
        + Polynomial.C 3 * Polynomial.X ^ 2
        - Polynomial.C 38 * Polynomial.X
        - Polynomial.C 19 := by
  native_decide

/-- First coefficient in an explicit Bézout identity for `p5` and `p5'`. -/
def bezoutU : Polynomial ℚ :=
  Polynomial.C (50 / 4617) * Polynomial.X ^ 3
    + Polynomial.C (5 / 1539) * Polynomial.X ^ 2
    - Polynomial.C (25 / 3078) * Polynomial.X
    - Polynomial.C (338 / 4617)

/-- Second coefficient in an explicit Bézout identity for `p5` and `p5'`. -/
def bezoutV : Polynomial ℚ :=
  -Polynomial.C (10 / 4617) * Polynomial.X ^ 4
    - Polynomial.C (5 / 4617) * Polynomial.X ^ 3
    + Polynomial.C (1 / 1026) * Polynomial.X ^ 2
    + Polynomial.C (371 / 9234) * Polynomial.X
    + Polynomial.C (5 / 243)

/-- Constructive squarefreeness certificate for the quintic. -/
theorem p5_bezout : bezoutU * p5 + bezoutV * p5.derivative = 1 := by
  native_decide

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
