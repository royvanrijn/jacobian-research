/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.Jacobian
import FiniteEtaleKeller.UniversalPromotedBlock
import Mathlib.Algebra.MvPolynomial.Equiv
import Mathlib.Algebra.MvPolynomial.PDeriv

/-!
# Literal unchanged-parameter promotion of a polynomial map

Let `F` be a `q`-variable polynomial map whose coefficients are polynomials
in a finite parameter type `p`.  Mathlib's `sumAlgEquiv` flattens those nested
polynomials to polynomials in the literal variable type `q ⊕ p`.

This module defines the promoted self-map

`(x, u) ↦ (F_u(x), u)`

and identifies its actual `MvPolynomial` Jacobian with the block matrix

`[[D_x F, D_u F], [0, I]]`.

It then proves that flattening preserves the vertical determinant, so a
determinant-one nested family promotes to a determinant-one literal
polynomial self-map.
-/

noncomputable section

open Matrix
open MvPolynomial

namespace FiniteEtaleKeller

variable {R : Type*} [CommRing R]
variable {p q : Type*}

/-- Flatten `q` polynomial variables with polynomial coefficients in the
parameter variables `p` into the literal variable type `q ⊕ p`. -/
def flattenParameterPolynomial :
    MvPolynomial q (MvPolynomial p R) →+*
      MvPolynomial (q ⊕ p) R :=
  (MvPolynomial.sumAlgEquiv R q p).symm.toRingEquiv.toRingHom

@[simp]
theorem flattenParameterPolynomial_C_C (r : R) :
    flattenParameterPolynomial (p := p) (q := q)
        (C (C r)) = C r := by
  simp [flattenParameterPolynomial]

@[simp]
theorem flattenParameterPolynomial_X (i : q) :
    flattenParameterPolynomial (p := p) (q := q)
        (X i : MvPolynomial q (MvPolynomial p R)) =
      (X (Sum.inl i) : MvPolynomial (q ⊕ p) R) := by
  simp [flattenParameterPolynomial]

@[simp]
theorem flattenParameterPolynomial_C_X (j : p) :
    flattenParameterPolynomial (p := p) (q := q)
        (C (X j : MvPolynomial p R)) =
      (X (Sum.inr j) : MvPolynomial (q ⊕ p) R) := by
  simp [flattenParameterPolynomial]

/-- Flattening commutes with differentiation in a vertical variable. -/
theorem pderiv_flattenParameterPolynomial
    (j : q) (f : MvPolynomial q (MvPolynomial p R)) :
    pderiv (Sum.inl j)
        (flattenParameterPolynomial (p := p) (q := q) f) =
      flattenParameterPolynomial (p := p) (q := q) (pderiv j f) := by
  apply (MvPolynomial.sumAlgEquiv R q p).injective
  rw [← MvPolynomial.pderiv_sumAlgEquiv]
  simp [flattenParameterPolynomial]

/-- A coefficient-ring homomorphism commutes with the Jacobian
determinant. -/
theorem jacobianDet_map
    [Fintype q] [DecidableEq q]
    {S : Type*} [CommRing S] (φ : R →+* S)
    (F : q → MvPolynomial q R) :
    jacobianDet (fun i => MvPolynomial.map φ (F i)) =
      MvPolynomial.map φ (jacobianDet F) := by
  rw [jacobianDet, jacobianDet]
  have hmatrix :
      jacobianMatrix (fun i => MvPolynomial.map φ (F i)) =
        (jacobianMatrix F).map (MvPolynomial.map φ) := by
    apply Matrix.ext
    intro i j
    simpa [jacobianMatrix] using
      (MvPolynomial.pderiv_map
        (R := R) (S := S) (σ := q) (φ := φ) (f := F i) (i := j))
  rw [hmatrix]
  exact ((MvPolynomial.map φ).map_det (jacobianMatrix F)).symm

/-- Literal promotion of a polynomial family with polynomial parameters.
Vertical coordinates come first and parameter coordinates are unchanged. -/
def unchangedParameterPromotion
    (F : q → MvPolynomial q (MvPolynomial p R)) :
    q ⊕ p → MvPolynomial (q ⊕ p) R :=
  Sum.elim
    (fun i => flattenParameterPolynomial (p := p) (F i))
    (fun j => X (Sum.inr j))

/-- The vertical derivative block of the literal promoted map. -/
def promotedVerticalJacobian
    (F : q → MvPolynomial q (MvPolynomial p R)) :
    Matrix q q (MvPolynomial (q ⊕ p) R) :=
  fun i j =>
    flattenParameterPolynomial (p := p) (pderiv j (F i))

/-- The parameter derivative block of the literal promoted map. -/
def promotedParameterJacobian
    (F : q → MvPolynomial q (MvPolynomial p R)) :
    Matrix q p (MvPolynomial (q ⊕ p) R) :=
  fun i j =>
    pderiv (Sum.inr j)
      (flattenParameterPolynomial (p := p) (F i))

/-- The actual Jacobian of the promoted `MvPolynomial (q ⊕ p)` map is the
expected unchanged-coordinate block matrix. -/
theorem jacobianMatrix_unchangedParameterPromotion
    [DecidableEq p] [DecidableEq q]
    (F : q → MvPolynomial q (MvPolynomial p R)) :
    jacobianMatrix (unchangedParameterPromotion F) =
      Matrix.fromBlocks
        (promotedVerticalJacobian F)
        (promotedParameterJacobian F)
        0 1 := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · simp [jacobianMatrix, unchangedParameterPromotion,
      promotedVerticalJacobian, pderiv_flattenParameterPolynomial]
  · rfl
  · simp [jacobianMatrix, unchangedParameterPromotion]
  · classical
    simp [jacobianMatrix, unchangedParameterPromotion, Matrix.one_apply,
      MvPolynomial.pderiv_X, Pi.single_apply]

/-- The promoted determinant is the determinant of its vertical block. -/
theorem jacobianDet_unchangedParameterPromotion_eq_vertical
    [Fintype p] [DecidableEq p] [Fintype q] [DecidableEq q]
    (F : q → MvPolynomial q (MvPolynomial p R)) :
    jacobianDet (unchangedParameterPromotion F) =
      (promotedVerticalJacobian F).det := by
  rw [jacobianDet, jacobianMatrix_unchangedParameterPromotion,
    Matrix.det_fromBlocks_zero₂₁, Matrix.det_one, mul_one]

/-- The determinant of the vertical block is the flattened nested Jacobian
determinant. -/
theorem det_promotedVerticalJacobian
    [Fintype q] [DecidableEq q]
    (F : q → MvPolynomial q (MvPolynomial p R)) :
    (promotedVerticalJacobian F).det =
      flattenParameterPolynomial (p := p) (jacobianDet F) := by
  rw [jacobianDet]
  have hmatrix :
      promotedVerticalJacobian F =
        (jacobianMatrix F).map
          (flattenParameterPolynomial (p := p) (q := q)) := by
    ext i j
    rfl
  rw [hmatrix]
  exact
    ((flattenParameterPolynomial (p := p) (q := q)).map_det
      (jacobianMatrix F)).symm

/-- Unchanged-coordinate promotion preserves determinant one. -/
theorem jacobianDet_unchangedParameterPromotion
    [Fintype p] [DecidableEq p] [Fintype q] [DecidableEq q]
    (F : q → MvPolynomial q (MvPolynomial p R))
    (hF : jacobianDet F = 1) :
    jacobianDet (unchangedParameterPromotion F) = 1 := by
  rw [jacobianDet_unchangedParameterPromotion_eq_vertical,
    det_promotedVerticalJacobian, hF, map_one]

#print axioms jacobianMatrix_unchangedParameterPromotion
#print axioms jacobianDet_unchangedParameterPromotion_eq_vertical
#print axioms det_promotedVerticalJacobian
#print axioms jacobianDet_unchangedParameterPromotion
#print axioms jacobianDet_map

end FiniteEtaleKeller
