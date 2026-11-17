/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeJacobian
import FiniteEtaleKeller.ExplicitPolynomial

/-!
# The explicit optimal quintic map certificate

This module certifies the three-variable polynomial map, its denominator-free
Jacobian `-722`, its determinant-`-2` quadratic-gauge normalization, and the
fixed determinant-one output normalization.  The determinant calculation is a
specialization of the all-degree theorem.  The displayed integral map is then
identified with the output scaling `diag(1,19,19)` of that specialization,
without expanding the high powers of the recurrent polynomial `q`.
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

/-- The determinant-`-2` quadratic-gauge normalization.  Defining it from the
universal map keeps the expensive Jacobian calculation entirely general. -/
def normalizedMap : Fin 3 → M := generalGaugeMap g5

/-- The target-preserving determinant-one normalization from the main theorem. -/
def jacobianOneMap : Fin 3 → M :=
  scaleOutput 1 (-1 / 2 : ℚ) 1 normalizedMap

@[simp]
private theorem t_eq_generalGaugeT : t = (generalGaugeT : M) := rfl

@[simp]
private theorem q_eq_generalGaugeQ : q = generalGaugeQ g5 := by
  simp [q, generalGaugeQ, t, generalGaugeT, g5, Polynomial.coeff_X]
  ring

private theorem g5_natDegree : g5.natDegree = 5 := by
  unfold g5
  compute_degree!

private theorem scale_coefficient_two (a b : M) (r : ℚ) :
    MvPolynomial.C 19 * a * b * MvPolynomial.C r =
      MvPolynomial.C (19 * r) * a * b := by
  calc
    _ = MvPolynomial.C 19 * MvPolynomial.C r * a * b := by ring
    _ = MvPolynomial.C (19 * r) * a * b := by
      rw [← MvPolynomial.C_mul]

private theorem scale_coefficient_three (a b c : M) (r : ℚ) :
    MvPolynomial.C 19 * a * b * c * MvPolynomial.C r =
      MvPolynomial.C (19 * r) * a * b * c := by
  calc
    _ = MvPolynomial.C 19 * MvPolynomial.C r * a * b * c := by ring
    _ = MvPolynomial.C (19 * r) * a * b * c := by
      rw [← MvPolynomial.C_mul]

/-- The explicitly displayed normalized quintic map is literally the
all-degree quadratic gauge attached to the seed `g5`. -/
theorem normalizedMap_eq_generalGaugeMap :
    normalizedMap = generalGaugeMap g5 := rfl

/-- The normalized quadratic gauge has Jacobian determinant `-2`. -/
theorem jacobianDet_normalizedMap :
    jacobianDet normalizedMap = MvPolynomial.C (-2) := by
  rw [normalizedMap_eq_generalGaugeMap]
  apply jacobianDet_generalGaugeMap g5
  · norm_num [g5, Polynomial.coeff_X]
  · norm_num [g5, Polynomial.coeff_X]

/-- Scaling the normalized gauge by `diag(1,19,19)` recovers the displayed
integer-coefficient map. -/
theorem integralMap_eq_scaled_normalized :
    integralMap = scaleOutput (1 : ℚ) 19 19 normalizedMap := by
  rw [normalizedMap_eq_generalGaugeMap]
  funext i
  fin_cases i
  · simp [integralMap, scaleOutput, generalGaugeMap, generalGaugePi]
  · change
      MvPolynomial.C 19 * MvPolynomial.X 1
          - MvPolynomial.C 3 * MvPolynomial.X 0 * q
          + MvPolynomial.C 38 * t * q
          - MvPolynomial.C 4 * t ^ 2 * MvPolynomial.X 0 ^ 2 * q ^ 4
          - MvPolynomial.C 5 * t ^ 2 * MvPolynomial.X 0 ^ 3 * q ^ 5 =
        MvPolynomial.C 19 * generalGaugeB g5
    rw [q_eq_generalGaugeQ, t_eq_generalGaugeT, generalGaugeB, g5_natDegree]
    rw [show Finset.Icc 4 5 = ({4, 5} : Finset ℕ) by decide]
    norm_num [g5, Polynomial.coeff_X]
    ring_nf
    simp_rw [scale_coefficient_two, scale_coefficient_three]
    norm_num
    simp only [map_ofNat]
    ring
  · change
      MvPolynomial.C 19 * MvPolynomial.X 0 *
            (MvPolynomial.C 5 - MvPolynomial.C 3 * t)
          + MvPolynomial.X 0 ^ 3 * MvPolynomial.X 2
          + MvPolynomial.C 2 * (MvPolynomial.X 0 * q) ^ 4
          + MvPolynomial.C 3 * (MvPolynomial.X 0 * q) ^ 5 =
        MvPolynomial.C 19 * generalGaugeC g5
    rw [q_eq_generalGaugeQ, t_eq_generalGaugeT, generalGaugeC, g5_natDegree]
    rw [show Finset.Icc 4 5 = ({4, 5} : Finset ℕ) by decide]
    norm_num [g5, Polynomial.coeff_X]
    ring_nf
    simp_rw [scale_coefficient_two]
    norm_num
    simp only [map_ofNat]
    ring

/-- The displayed denominator-free map has Jacobian determinant `-722`. -/
theorem jacobianDet_integralMap :
    jacobianDet integralMap = MvPolynomial.C (-722) := by
  rw [integralMap_eq_scaled_normalized, jacobianDet_scaleOutput,
    jacobianDet_normalizedMap]
  rw [← MvPolynomial.C_mul]
  norm_num

/-- The universal target-preserving output normalization has determinant `1`. -/
theorem jacobianDet_jacobianOneMap : jacobianDet jacobianOneMap = 1 := by
  rw [jacobianOneMap, jacobianDet_scaleOutput, jacobianDet_normalizedMap]
  rw [← MvPolynomial.C_mul]
  norm_num

#print axioms normalizedMap_eq_generalGaugeMap
#print axioms integralMap_eq_scaled_normalized
#print axioms jacobianDet_integralMap
#print axioms jacobianDet_jacobianOneMap

end FiniteEtaleKeller.ExplicitQuintic
