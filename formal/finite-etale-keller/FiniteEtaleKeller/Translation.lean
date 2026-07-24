/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib

/-!
# Translation and the realization target

This module formalizes the passage

`P(T) ↦ G(S) = P(a + S) - P(a)`

and the target value which makes the quadratic-gauge inverse equation equal to
the translated polynomial `P(a + S)`.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

/-- Translation of the polynomial variable by `a`. -/
def translatePolynomial (P : K[X]) (a : K) : K[X] :=
  P.comp (X + C a)

/-- The translated polynomial with its constant term removed. -/
def rootedTranslate (P : K[X]) (a : K) : K[X] :=
  translatePolynomial P a - C (P.eval a)

/-- Evaluating the translated polynomial at `s` is the same as evaluating the
original polynomial at `s + a`. -/
theorem translatePolynomial_eval (P : K[X]) (a s : K) :
    (translatePolynomial P a).eval s = P.eval (s + a) := by
  simp [translatePolynomial]

/-- The rooted translation has zero constant term. -/
@[simp]
theorem rootedTranslate_eval_zero (P : K[X]) (a : K) :
    (rootedTranslate P a).eval 0 = 0 := by
  simp [rootedTranslate, translatePolynomial]

/-- Restoring the removed constant recovers the translated polynomial. -/
theorem rootedTranslate_add_constant (P : K[X]) (a : K) :
    rootedTranslate P a + C (P.eval a) = translatePolynomial P a := by
  unfold rootedTranslate
  ring

/-- The third target coordinate used by the realization theorem. -/
def realizationTargetC (P : K[X]) (a g₁ : K) : K :=
  -2 * P.eval a / g₁

/-- At the realization target, the inverse polynomial is exactly the translated
input polynomial. -/
theorem rootedTranslate_inverse_at_target [CharZero K]
    (P : K[X]) (a g₁ : K) (hg₁ : g₁ ≠ 0) :
    rootedTranslate P a
        - C (g₁ / 2 * realizationTargetC P a g₁)
      = translatePolynomial P a := by
  have hscalar : g₁ / 2 * realizationTargetC P a g₁ = -P.eval a := by
    unfold realizationTargetC
    field_simp [hg₁]
    ring
  rw [hscalar]
  simp [rootedTranslate]

/-- Combining evaluation with the target identity gives the exact root
correspondence `S ↔ T = S + a`. -/
theorem realizationTarget_eval [CharZero K]
    (P : K[X]) (a g₁ s : K) (hg₁ : g₁ ≠ 0) :
    (rootedTranslate P a
        - C (g₁ / 2 * realizationTargetC P a g₁)).eval s
      = P.eval (s + a) := by
  rw [rootedTranslate_inverse_at_target P a g₁ hg₁]
  exact translatePolynomial_eval P a s

end FiniteEtaleKeller
