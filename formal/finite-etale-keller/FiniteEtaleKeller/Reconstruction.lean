/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib

/-!
# Algebraic reconstruction identities

These identities are the field-valued model for the scheme-theoretic inverse.
The later quotient-ring stage replaces `d⁻¹` by the explicit unit supplied by
a Bézout identity for the inverse polynomial and its derivative.
-/

noncomputable section

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

/-- Once `d = 1 - S*Q + Π*S^2` is invertible, the reconstruction formulas
recover the source chart and its defining equation exactly. -/
theorem reconstruction_identities
    (S Q Π d a : K) (hd : d ≠ 0)
    (hD : d = 1 - S * Q + Π * S ^ 2) :
    let t := d⁻¹
    let x := S * d⁻¹
    let y := Q - Π * S
    let q := Π * d
    let z := d ^ 2 * (q - a * y ^ 2 * (1 + 3 * t))
    (1 + x * y = t)
      ∧ (t ^ 2 * z + a * y ^ 2 * (1 + 3 * t) = q)
      ∧ (x / t = S)
      ∧ (y + x * q = Q)
      ∧ (t * q = Π) := by
  dsimp
  constructor
  · field_simp [hd]
    linear_combination hD
  constructor
  · field_simp [hd]
    ring
  constructor
  · field_simp [hd]
  constructor
  · field_simp [hd]
    ring
  · field_simp [hd]

/-- The same reconstruction also recovers the marked coordinate `D = t⁻¹`. -/
theorem reconstruction_recovers_D (d : K) (hd : d ≠ 0) :
    (d⁻¹)⁻¹ = d := by
  exact inv_inv d

end FiniteEtaleKeller
