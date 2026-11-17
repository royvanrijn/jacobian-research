/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib

/-!
# Jacobian helpers

Minimal multivariate-polynomial definitions used by the finite-etale Keller
fiber certificates.
-/

noncomputable section

open Matrix Function
open MvPolynomial

namespace FiniteEtaleKeller

variable {R : Type*} {σ : Type*}

/-- Jacobian matrix of a family of multivariate polynomials. -/
def jacobianMatrix [CommSemiring R] [DecidableEq σ]
    (F : σ → MvPolynomial σ R) : Matrix σ σ (MvPolynomial σ R) :=
  Matrix.of fun i j => pderiv j (F i)

/-- Jacobian determinant of a family of multivariate polynomials. -/
def jacobianDet [CommRing R] [Fintype σ] [DecidableEq σ]
    (F : σ → MvPolynomial σ R) : MvPolynomial σ R :=
  (jacobianMatrix F).det

/-- Polynomial self-map induced by a family of multivariate polynomials. -/
def evalMap [CommSemiring R] (F : σ → MvPolynomial σ R) (p : σ → R) : σ → R :=
  fun i => eval p (F i)

/-- Scale the three output coordinates of a polynomial map. -/
def scaleOutput [CommSemiring R] (a b c : R)
    (F : Fin 3 → MvPolynomial (Fin 3) R) :
    Fin 3 → MvPolynomial (Fin 3) R :=
  ![C a * F 0, C b * F 1, C c * F 2]

/-- Scaling the three output rows scales the Jacobian determinant by `a*b*c`. -/
theorem jacobianDet_scaleOutput [CommRing R] (a b c : R)
    (F : Fin 3 → MvPolynomial (Fin 3) R) :
    jacobianDet (scaleOutput a b c F) = C (a * b * c) * jacobianDet F := by
  simp only [jacobianDet, jacobianMatrix, scaleOutput, det_fin_three, of_apply,
    cons_val_zero, cons_val_one, cons_val_two, head_cons, tail_cons, pderiv_mul,
    pderiv_C, zero_mul, zero_add, map_mul]
  ring

end FiniteEtaleKeller
