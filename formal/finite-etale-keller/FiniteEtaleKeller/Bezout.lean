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

/-- Specialization to the derivative class used in the fiber reconstruction. -/
theorem adjoinRoot_derivative_inverse (E U V : K[X])
    (h : U * E + V * E.derivative = 1) :
    AdjoinRoot.mk E V * AdjoinRoot.mk E E.derivative = 1 :=
  adjoinRoot_bezout_inverse E U V E.derivative h

end FiniteEtaleKeller
