/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib

/-!
# Algebraic reconstruction identities

These identities are the algebraic core of the scheme-theoretic inverse.  The
field-valued theorem mirrors the rational reconstruction formulas; the unit
version works in an arbitrary commutative ring and is the form needed inside a
polynomial quotient algebra.
-/

noncomputable section

namespace FiniteEtaleKeller

section FieldReconstruction

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

end FieldReconstruction

section UnitReconstruction

variable {R : Type*} [CommRing R]

/-- Ring-level reconstruction using an explicit unit.  This is the form used
in `K[S]/(E)`, where Bézout makes the derivative class a unit without any
localization. -/
theorem unitReconstruction_identities
    (S Q Π a : R) (d : Rˣ)
    (hD : (d : R) = 1 - S * Q + Π * S ^ 2) :
    let t : R := ↑d⁻¹
    let x := S * t
    let y := Q - Π * S
    let q := Π * (d : R)
    let z := (d : R) ^ 2 * (q - a * y ^ 2 * (1 + 3 * t))
    (1 + x * y = t)
      ∧ (t ^ 2 * z + a * y ^ 2 * (1 + 3 * t) = q)
      ∧ (x * (d : R) = S)
      ∧ (y + x * q = Q)
      ∧ (t * q = Π)
      ∧ (t * (d : R) = 1) := by
  dsimp
  have hunit : (↑d⁻¹ : R) * (d : R) = 1 := by
    simp
  have hsum : (d : R) + S * Q - Π * S ^ 2 = 1 := by
    rw [hD]
    ring
  have hsq : (↑d⁻¹ : R) ^ 2 * (d : R) ^ 2 = 1 := by
    rw [← mul_pow]
    simp
  constructor
  · calc
      1 + S * (↑d⁻¹ : R) * (Q - Π * S)
          = (↑d⁻¹ : R) * ((d : R) + S * Q - Π * S ^ 2) := by
              rw [hunit]
              ring
      _ = (↑d⁻¹ : R) := by rw [hsum, mul_one]
  constructor
  · calc
      (↑d⁻¹ : R) ^ 2
            * ((d : R) ^ 2
              * (Π * (d : R)
                - a * (Q - Π * S) ^ 2 * (1 + 3 * (↑d⁻¹ : R))))
          + a * (Q - Π * S) ^ 2 * (1 + 3 * (↑d⁻¹ : R))
          = ((↑d⁻¹ : R) ^ 2 * (d : R) ^ 2)
              * (Π * (d : R)
                - a * (Q - Π * S) ^ 2 * (1 + 3 * (↑d⁻¹ : R)))
            + a * (Q - Π * S) ^ 2 * (1 + 3 * (↑d⁻¹ : R)) := by ring
      _ = Π * (d : R) := by rw [hsq]; ring
  constructor
  · calc
      S * (↑d⁻¹ : R) * (d : R) = S * ((↑d⁻¹ : R) * (d : R)) := by ring
      _ = S := by rw [hunit, mul_one]
  constructor
  · calc
      Q - Π * S + S * (↑d⁻¹ : R) * (Π * (d : R))
          = Q - Π * S + S * Π * ((↑d⁻¹ : R) * (d : R)) := by ring
      _ = Q := by rw [hunit]; ring
  constructor
  · calc
      (↑d⁻¹ : R) * (Π * (d : R))
          = Π * ((↑d⁻¹ : R) * (d : R)) := by ring
      _ = Π := by rw [hunit, mul_one]
  · exact hunit

end UnitReconstruction

end FiniteEtaleKeller
