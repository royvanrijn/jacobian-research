/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.RingTheory.AdjoinRoot

/-!
# Bezout inversion in a polynomial quotient

The scheme reconstruction uses no hidden localization: a Bezout identity
`U * E + V * E' = 1` gives the inverse of the derivative class explicitly in
`K[S]/(E)`.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [CommRing K]

/-- A Bezout identity gives an explicit inverse after quotienting by `E`. -/
theorem adjoinRoot_bezout_inverse (E U V W : K[X])
    (h : U * E + V * W = 1) :
    AdjoinRoot.mk E V * AdjoinRoot.mk E W = 1 := by
  have hq := congrArg (AdjoinRoot.mk E) h
  simpa only [map_add, map_mul, map_one, AdjoinRoot.mk_self, mul_zero,
    zero_add] using hq

/-- Package the two-sided quotient inverse supplied by a Bezout identity as an
actual unit. -/
def adjoinRootUnitOfBezout (E U V W : K[X])
    (h : U * E + V * W = 1) : (AdjoinRoot E)ˣ where
  val := AdjoinRoot.mk E W
  inv := AdjoinRoot.mk E V
  val_inv := by
    simpa [mul_comm] using adjoinRoot_bezout_inverse E U V W h
  inv_val := adjoinRoot_bezout_inverse E U V W h

/-- Specialization to the derivative class used in the fiber reconstruction. -/
theorem adjoinRoot_derivative_inverse (E U V : K[X])
    (h : U * E + V * E.derivative = 1) :
    AdjoinRoot.mk E V * AdjoinRoot.mk E E.derivative = 1 :=
  adjoinRoot_bezout_inverse E U V E.derivative h

/-- The derivative class as an explicit unit of `K[S]/(E)`. -/
def adjoinRootDerivativeUnit (E U V : K[X])
    (h : U * E + V * E.derivative = 1) : (AdjoinRoot E)ˣ :=
  adjoinRootUnitOfBezout E U V E.derivative h

@[simp]
theorem adjoinRootDerivativeUnit_val (E U V : K[X])
    (h : U * E + V * E.derivative = 1) :
    (adjoinRootDerivativeUnit E U V h : AdjoinRoot E) =
      AdjoinRoot.mk E E.derivative := rfl

@[simp]
theorem adjoinRootDerivativeUnit_inv_val (E U V : K[X])
    (h : U * E + V * E.derivative = 1) :
    (↑((adjoinRootDerivativeUnit E U V h)⁻¹) : AdjoinRoot E) =
      AdjoinRoot.mk E V := rfl

/-- In particular, the derivative class is a unit in the quotient ring. -/
theorem adjoinRoot_derivative_isUnit (E U V : K[X])
    (h : U * E + V * E.derivative = 1) :
    IsUnit (AdjoinRoot.mk E E.derivative) := by
  simpa using (adjoinRootDerivativeUnit E U V h).isUnit

end FiniteEtaleKeller
