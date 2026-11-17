/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

/-!
# The determinant block used by unchanged-coordinate promotion

Promoting parameters to unchanged coordinates gives a Jacobian matrix

`[[I,0],[C,D]]`.

This module formalizes the determinant step for arbitrary finite parameter and
vertical index types.  Identifying the literal promoted map's Jacobian with
this block matrix remains a separate `MvPolynomial (Fin N)` formalization
task.
-/

namespace FiniteEtaleKeller

open Matrix

variable {R : Type*} [CommRing R]
variable {p q : Type*}
variable [Fintype p] [DecidableEq p] [Fintype q] [DecidableEq q]

/-- An unchanged-coordinate block promotion has the determinant of its
vertical block. -/
theorem det_unchangedCoordinateBlock
    (C : Matrix q p R) (D : Matrix q q R) :
    (Matrix.fromBlocks (1 : Matrix p p R) 0 C D).det = D.det := by
  rw [Matrix.det_fromBlocks_zero₁₂, Matrix.det_one, one_mul]

/-- If the vertical determinant is one, so is the promoted determinant. -/
theorem det_unchangedCoordinateBlock_eq_one
    (C : Matrix q p R) (D : Matrix q q R) (hD : D.det = 1) :
    (Matrix.fromBlocks (1 : Matrix p p R) 0 C D).det = 1 := by
  rw [det_unchangedCoordinateBlock C D, hD]

#print axioms det_unchangedCoordinateBlock
#print axioms det_unchangedCoordinateBlock_eq_one

end FiniteEtaleKeller
