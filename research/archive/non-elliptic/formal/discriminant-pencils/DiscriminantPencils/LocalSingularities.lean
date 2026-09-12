/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import DiscriminantPencils.TangentPencil

/-!
# Local singularity calculations

This file proves the jet and transversality calculations used in Section 3
of *Generic Discriminants of Polynomial Tangent Pencils*.
-/

open Polynomial

namespace DiscriminantPencils

variable {𝕜 : Type*} [Field 𝕜]

/-- The tangent-map parametrization is unramified at `r` exactly when the
second derivative of `H` does not vanish there. -/
def IsUnramified (H : 𝕜[X]) (r : 𝕜) : Prop :=
  H.derivative.derivative.eval r ≠ 0

/-- A parameter at which the first nonzero velocity jet is the quadratic
one: `H''(r) = 0` but `H'''(r) ≠ 0`. -/
def IsOrdinaryCuspParameter (H : 𝕜[X]) (r : 𝕜) : Prop :=
  H.derivative.derivative.eval r = 0 ∧
    H.derivative.derivative.derivative.eval r ≠ 0

/-- Contact of the graph with its tangent line is always at least two. -/
theorem tangent_has_contact_two (H : 𝕜[X]) (r : 𝕜) :
    IsRepeatedRoot
      (pencil H (H.derivative.eval r)
        (r * H.derivative.eval r - H.eval r)) r := by
  simpa [tangentMap] using repeatedRoot_at_tangentMap H r

/-- The first two derivatives of the tangent-map coordinates. -/
theorem tangent_second_jet (H : 𝕜[X]) (r : 𝕜) :
    ((tangentMapPoly H).1.derivative.derivative.eval r,
      (tangentMapPoly H).2.derivative.derivative.eval r) =
      (H.derivative.derivative.derivative.eval r,
        H.derivative.derivative.eval r +
          r * H.derivative.derivative.derivative.eval r) := by
  simp [tangentMapPoly]

/-- The first three derivatives of the tangent-map coordinates. -/
theorem tangent_third_jet (H : 𝕜[X]) (r : 𝕜) :
    ((tangentMapPoly H).1.derivative.derivative.derivative.eval r,
      (tangentMapPoly H).2.derivative.derivative.derivative.eval r) =
      (H.derivative.derivative.derivative.derivative.eval r,
        2 * H.derivative.derivative.derivative.eval r +
          r * H.derivative.derivative.derivative.derivative.eval r) := by
  simp [tangentMapPoly]
  ring

/-- At a simple zero of `H''`, the determinant of the quadratic and cubic
jets is `2 H'''(r)^2`, equation (3.2). -/
theorem ordinary_cusp_jet_det [CharZero 𝕜] (H : 𝕜[X]) (r : 𝕜)
    (hr : IsOrdinaryCuspParameter H r) :
    let s₂ := (tangentMapPoly H).1.derivative.derivative.eval r
    let t₂ := (tangentMapPoly H).2.derivative.derivative.eval r
    let s₃ := (tangentMapPoly H).1.derivative.derivative.derivative.eval r
    let t₃ := (tangentMapPoly H).2.derivative.derivative.derivative.eval r
    s₂ * t₃ - t₂ * s₃ =
      2 * (H.derivative.derivative.derivative.eval r) ^ 2 := by
  dsimp
  change H.derivative.derivative.eval r = 0 ∧
    H.derivative.derivative.derivative.eval r ≠ 0 at hr
  simp [tangentMapPoly]
  rw [hr.1]
  ring

theorem ordinary_cusp_jet_det_ne_zero [CharZero 𝕜] (H : 𝕜[X]) (r : 𝕜)
    (hr : IsOrdinaryCuspParameter H r) :
    let s₂ := (tangentMapPoly H).1.derivative.derivative.eval r
    let t₂ := (tangentMapPoly H).2.derivative.derivative.eval r
    let s₃ := (tangentMapPoly H).1.derivative.derivative.derivative.eval r
    let t₃ := (tangentMapPoly H).2.derivative.derivative.derivative.eval r
    s₂ * t₃ - t₂ * s₃ ≠ 0 := by
  dsimp
  rw [ordinary_cusp_jet_det H r hr]
  change H.derivative.derivative.eval r = 0 ∧
    H.derivative.derivative.derivative.eval r ≠ 0 at hr
  exact mul_ne_zero (by norm_num) (pow_ne_zero 2 hr.2)

/-- The determinant of the two actual velocity vectors at parameters `r`
and `u`. -/
theorem tangent_velocity_det (H : 𝕜[X]) (r u : 𝕜) :
    let ar := H.derivative.derivative.eval r
    let au := H.derivative.derivative.eval u
    ar * (u * au) - (r * ar) * au = ar * au * (u - r) := by
  dsimp
  ring

/-- Two distinct unramified normalization parameters have transverse image
branches. This is the node transversality assertion in Section 3. -/
theorem distinct_unramified_branches_transverse (H : 𝕜[X]) {r u : 𝕜}
    (hru : r ≠ u) (hr : IsUnramified H r) (hu : IsUnramified H u) :
    let ar := H.derivative.derivative.eval r
    let au := H.derivative.derivative.eval u
    ar * (u * au) - (r * ar) * au ≠ 0 := by
  dsimp
  rw [tangent_velocity_det]
  change H.derivative.derivative.eval r ≠ 0 at hr
  change H.derivative.derivative.eval u ≠ 0 at hu
  exact mul_ne_zero (mul_ne_zero hr hu) (sub_ne_zero.mpr hru.symm)

end DiscriminantPencils
