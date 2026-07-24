/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib

/-!
# Universal marked-line identities

This file formalizes the algebraic core of the quadratic-gauge proof before
encoding the full polynomial family.  It contains the reciprocal source-chart
identity, the inverse-equation cancellation, and the two-by-two Jacobian
cancellation responsible for the constant determinant.
-/

noncomputable section

namespace FiniteEtaleKeller

section SourceChart

variable {K : Type*} [Field K]

/-- On the chart `t = 1 + x*y ≠ 0`, the source factor
`D = 1 - S*Q + pi*S^2` is exactly `t⁻¹`. -/
theorem sourceChart_reciprocal (x y q : K) (ht : 1 + x * y ≠ 0) :
    let t := 1 + x * y
    let S := x / t
    let pi := t * q
    let Q := y + x * q
    1 - S * Q + pi * S ^ 2 = 1 / t := by
  dsimp
  field_simp [ht]
  ring

/-- If `t*q = pi` and `pi` is nonzero, then both chart factors are nonzero.
This is the pointwise algebra behind the scheme-theoretic observation that
`t` and `q` are units on a fiber with nonzero first target coordinate. -/
theorem sourceChart_factors_ne_zero {t q pi : K} (hpi : pi ≠ 0) (htq : t * q = pi) :
    t ≠ 0 ∧ q ≠ 0 := by
  constructor
  · intro ht
    apply hpi
    simpa [ht] using htq.symm
  · intro hq
    apply hpi
    simpa [hq] using htq.symm

end SourceChart

section MarkedLine

variable {R : Type*} [CommRing R]

/-- The defining relation for `β` converts the normalized derivative into the
same factor `D` that appears in the source chart. -/
theorem normalizedDerivative_eq_chartFactor
    (S Q pi β h : R) (hβ : S * β = h - 1 - pi * S ^ 2) :
    h - (Q + β) * S = 1 - S * Q + pi * S ^ 2 := by
  have hh : h = 1 + pi * S ^ 2 + S * β := by
    linear_combination hβ
  rw [hh]
  ring

/-- The determinant of the marked-line differential is `-2` times the
normalized derivative factor.  The terms involving `β'` cancel identically. -/
theorem markedLine_planeJacobian
    (S B h β' : R) :
    β' * (-S ^ 2) - (2 * h - β' * S ^ 2 - 2 * B * S)
      = -2 * (h - B * S) := by
  ring

/-- Combining the `β` relation with the differential cancellation identifies
the marked-line Jacobian with `-2D`, ready to cancel the reciprocal source
Jacobian. -/
theorem markedLine_planeJacobian_eq_chartFactor
    (S Q pi β h β' : R) (hβ : S * β = h - 1 - pi * S ^ 2) :
    β' * (-S ^ 2)
        - (2 * h - β' * S ^ 2 - 2 * (Q + β) * S)
      = -2 * (1 - S * Q + pi * S ^ 2) := by
  calc
    β' * (-S ^ 2)
          - (2 * h - β' * S ^ 2 - 2 * (Q + β) * S)
        = -2 * (h - (Q + β) * S) := by ring
    _ = -2 * (1 - S * Q + pi * S ^ 2) := by
      rw [normalizedDerivative_eq_chartFactor S Q pi β h hβ]

end MarkedLine

section InverseEquation

variable {K : Type*} [Field K] [CharZero K]

/-- The displayed definition of `C` makes the inverse equation vanish
identically. -/
theorem markedLine_inverseEquation_zero
    (g₁ g S B : K) (hg₁ : g₁ ≠ 0) :
    g - (g₁ / 2) * (B * S ^ 2 + (2 * g / g₁ - B * S ^ 2)) = 0 := by
  field_simp [hg₁]
  ring

end InverseEquation

end FiniteEtaleKeller
