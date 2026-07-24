import FiniteEtaleKeller.Bezout
import FiniteEtaleKeller.Jacobian

/-!
# The explicit optimal quintic certificate

This first formalization stage certifies:

* the denominator-free degree-five map has Jacobian determinant `-722`;
* it is `diag(1,19,19)` of the determinant-`-2` quadratic gauge;
* a further fixed output scaling gives Jacobian determinant `1` while fixing
  the distinguished target because its second coordinate is zero;
* the inverse polynomial at the distinguished target is the Berend--Bilu
  quintic;
* an explicit Bezout identity supplies the inverse of its derivative in the
  quotient algebra.
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

/-- The fixed determinant-one normalization used by the main realization theorem. -/
def jacobianOneMap : Fin 3 → M :=
  scaleOutput 1 (-1 / 2 : ℚ) 1 normalizedMap

set_option maxHeartbeats 0 in
/-- The displayed denominator-free map has Jacobian determinant `-722`. -/
theorem jacobianDet_integralMap :
    jacobianDet integralMap = MvPolynomial.C (-722) := by
  simp only [jacobianDet, jacobianMatrix, det_fin_three, Matrix.of_apply,
    integralMap, t, q, Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.cons_val_two, Matrix.head_cons, Matrix.tail_cons, map_add, map_sub,
    Derivation.map_one_eq_zero, pderiv_mul, pderiv_pow, pderiv_C,
    pderiv_X_self, pderiv_X_of_ne, ne_eq, Fin.reduceEq, not_false_eq_true]
  simp only [map_neg, map_ofNat]
  ring

/-- The normalized quadratic gauge has Jacobian determinant `-2`. -/
theorem jacobianDet_normalizedMap :
    jacobianDet normalizedMap = MvPolynomial.C (-2) := by
  rw [normalizedMap, jacobianDet_scaleOutput, jacobianDet_integralMap]
  norm_num

/-- The universal output normalization gives Jacobian determinant `1`. -/
theorem jacobianDet_jacobianOneMap : jacobianDet jacobianOneMap = 1 := by
  rw [jacobianOneMap, jacobianDet_scaleOutput, jacobianDet_normalizedMap]
  norm_num

/-- Scaling the normalized gauge back by `diag(1,19,19)` recovers the displayed map. -/
theorem integralMap_eq_scaled_normalized :
    integralMap = scaleOutput (1 : ℚ) 19 19 normalizedMap := by
  funext i
  fin_cases i <;> simp [normalizedMap, scaleOutput] <;> ring

/-- The classical minimal intersective quintic. -/
def p5 : ℚ[Polynomial.X] :=
  (Polynomial.X ^ 3 - Polynomial.C 19)
    * (Polynomial.X ^ 2 + Polynomial.X + 1)

/-- Expanded form of the classical quintic. -/
theorem p5_expanded :
    p5 = Polynomial.X ^ 5 + Polynomial.X ^ 4 + Polynomial.X ^ 3
      - Polynomial.C 19 * Polynomial.X ^ 2
      - Polynomial.C 19 * Polynomial.X - Polynomial.C 19 := by
  simp [p5]
  ring

/-- The rooted quadratic-gauge seed. -/
def g5 : ℚ[Polynomial.X] :=
  Polynomial.X ^ 5 + Polynomial.X ^ 4 + Polynomial.X ^ 3
    - Polynomial.C 19 * Polynomial.X ^ 2
    - Polynomial.C 19 * Polynomial.X

/-- At normalized target `C = -2`, the inverse polynomial is exactly `p5`. -/
theorem inversePolynomial_eq_p5 : g5 - Polynomial.C 19 = p5 := by
  rw [p5_expanded]
  simp [g5]

/-- Explicit derivative of the quintic. -/
theorem p5_derivative :
    p5.derivative =
      Polynomial.C 5 * Polynomial.X ^ 4
        + Polynomial.C 4 * Polynomial.X ^ 3
        + Polynomial.C 3 * Polynomial.X ^ 2
        - Polynomial.C 38 * Polynomial.X
        - Polynomial.C 19 := by
  rw [p5_expanded]
  simp
  ring

/-- First coefficient in an explicit Bezout identity for `p5` and `p5'`. -/
def bezoutU : ℚ[Polynomial.X] :=
  Polynomial.C (50 / 4617) * Polynomial.X ^ 3
    + Polynomial.C (5 / 1539) * Polynomial.X ^ 2
    - Polynomial.C (25 / 3078) * Polynomial.X
    - Polynomial.C (338 / 4617)

/-- Second coefficient in an explicit Bezout identity for `p5` and `p5'`. -/
def bezoutV : ℚ[Polynomial.X] :=
  -Polynomial.C (10 / 4617) * Polynomial.X ^ 4
    - Polynomial.C (5 / 4617) * Polynomial.X ^ 3
    + Polynomial.C (1 / 1026) * Polynomial.X ^ 2
    + Polynomial.C (371 / 9234) * Polynomial.X
    + Polynomial.C (5 / 243)

/-- Constructive squarefreeness certificate for the quintic. -/
theorem p5_bezout : bezoutU * p5 + bezoutV * p5.derivative = 1 := by
  rw [p5_expanded, p5_derivative]
  simp [bezoutU, bezoutV]
  ring

/-- The class of `p5'` has the displayed inverse in `ℚ[X]/(p5)`. -/
theorem p5_derivative_inverse :
    AdjoinRoot.mk p5 bezoutV * AdjoinRoot.mk p5 p5.derivative = 1 :=
  adjoinRoot_derivative_inverse p5 bezoutU bezoutV p5_bezout

/-- The normalized target `(1,0,-2)` scales to the integral target `(1,0,-38)`. -/
theorem target_scaling : (19 : ℚ) * (-2) = -38 := by
  norm_num

#print axioms jacobianDet_integralMap
#print axioms jacobianDet_jacobianOneMap
#print axioms inversePolynomial_eq_p5
#print axioms p5_derivative_inverse

end FiniteEtaleKeller.ExplicitQuintic
