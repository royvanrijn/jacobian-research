/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib

/-!
# The tangent-pencil parametrization

This file formalizes equations (1.2), (1.3), and the differential identities
(3.1), (3.2) from *Generic Discriminants of Polynomial Tangent Pencils*.
-/

open Polynomial

namespace DiscriminantPencils

variable {𝕜 : Type*} [Field 𝕜]

/-- The polynomial pencil `H(W) - sW + t`. -/
noncomputable def pencil (H : 𝕜[X]) (s t : 𝕜) : 𝕜[X] :=
  H - C s * X + C t

/-- The affine tangent-discriminant parametrization
`r ↦ (H'(r), r H'(r) - H(r))`. -/
noncomputable def tangentMap (H : 𝕜[X]) (r : 𝕜) : 𝕜 × 𝕜 :=
  (H.derivative.eval r, r * H.derivative.eval r - H.eval r)

/-- A root detected together with the vanishing of the first derivative.

Over a field this is equivalent to root multiplicity at least two, but this
equational form is exactly what the tangent-pencil calculation needs.
-/
def IsRepeatedRoot (p : 𝕜[X]) (r : 𝕜) : Prop :=
  p.eval r = 0 ∧ p.derivative.eval r = 0

@[simp] theorem eval_pencil (H : 𝕜[X]) (s t r : 𝕜) :
    (pencil H s t).eval r = H.eval r - s * r + t := by
  simp [pencil]

@[simp] theorem derivative_pencil (H : 𝕜[X]) (s t : 𝕜) :
    (pencil H s t).derivative = H.derivative - C s := by
  simp [pencil]

@[simp] theorem eval_derivative_pencil (H : 𝕜[X]) (s t r : 𝕜) :
    (pencil H s t).derivative.eval r = H.derivative.eval r - s := by
  simp

/-- Equation (1.2): a parameter is a repeated root exactly at its tangent
parameter pair. -/
theorem isRepeatedRoot_pencil_iff (H : 𝕜[X]) (s t r : 𝕜) :
    IsRepeatedRoot (pencil H s t) r ↔ tangentMap H r = (s, t) := by
  simp only [IsRepeatedRoot, eval_pencil, eval_derivative_pencil, tangentMap,
    Prod.mk.injEq]
  constructor
  · rintro ⟨hp, hd⟩
    constructor
    · exact sub_eq_zero.mp hd
    · apply sub_eq_zero.mp
      calc
        r * H.derivative.eval r - H.eval r - t =
            -(H.eval r - H.derivative.eval r * r + t) := by ring
        _ = 0 := by rw [sub_eq_zero.mp hd, hp]; ring
  · rintro ⟨hs, ht⟩
    subst s
    subst t
    constructor
    · ring
    · ring

theorem repeatedRoot_at_tangentMap (H : 𝕜[X]) (r : 𝕜) :
    IsRepeatedRoot (pencil H (tangentMap H r).1 (tangentMap H r).2) r :=
  (isRepeatedRoot_pencil_iff H _ _ r).2 rfl

/-- Coordinatewise polynomial model of the tangent map. -/
noncomputable def tangentMapPoly (H : 𝕜[X]) : 𝕜[X] × 𝕜[X] :=
  (H.derivative, X * H.derivative - H)

@[simp] theorem tangentMapPoly_eval (H : 𝕜[X]) (r : 𝕜) :
    ((tangentMapPoly H).1.eval r, (tangentMapPoly H).2.eval r) = tangentMap H r := by
  simp [tangentMapPoly, tangentMap]

/-- The derivative of the second coordinate is `X H''`. -/
theorem derivative_tangentMapPoly_snd (H : 𝕜[X]) :
    (tangentMapPoly H).2.derivative = X * H.derivative.derivative := by
  simp [tangentMapPoly]

/-- Polynomial form of equation (3.1). -/
theorem tangent_velocity (H : 𝕜[X]) (r : 𝕜) :
    ((tangentMapPoly H).1.derivative.eval r,
      (tangentMapPoly H).2.derivative.eval r) =
      (H.derivative.derivative.eval r,
        r * H.derivative.derivative.eval r) := by
  simp [tangentMapPoly]

/-- The determinant of the tangent directions `(1,r)` and `(1,u)`. -/
theorem tangent_direction_det (r u : 𝕜) :
    (1 : 𝕜) * u - r * 1 = u - r := by
  ring

theorem tangent_direction_det_ne_zero {r u : 𝕜} (hru : r ≠ u) :
    (1 : 𝕜) * u - r * 1 ≠ 0 := by
  simpa [tangent_direction_det] using sub_ne_zero.mpr hru.symm

/-- Algebraic content of equation (3.2).  If `a = H'''(r)`, the determinant
of the displayed second and third jets is `2a²`. -/
theorem cusp_jet_det (r a b : 𝕜) :
    a * (2 * a + r * b) - (r * a) * b = 2 * a ^ 2 := by
  ring

theorem cusp_jet_det_ne_zero [CharZero 𝕜] {r a b : 𝕜} (ha : a ≠ 0) :
    a * (2 * a + r * b) - (r * a) * b ≠ 0 := by
  rw [cusp_jet_det]
  exact mul_ne_zero (by norm_num) (pow_ne_zero 2 ha)

end DiscriminantPencils
