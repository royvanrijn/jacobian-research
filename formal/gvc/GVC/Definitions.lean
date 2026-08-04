/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.Algebra.MvPolynomial.PDeriv
import Mathlib.Data.Nat.Factorial.Basic

/-!
# Constant-coefficient differential operators and GVC

The action is defined coefficientwise.  A symbol monomial `X^α` sends a
polynomial monomial `x^β` to

`(∏ i, βᵢ.descFactorial αᵢ) x^(β-α)`.

This formula also handles `α ≰ β`, since the descending factorial is then
zero.
-/

namespace GVC

open MvPolynomial

/-- The coefficient multiplier produced by differentiating `x^β` by the
multi-index `α`. -/
noncomputable def multiDescFactorial
    {σ : Type*} [Fintype σ] (β α : σ →₀ ℕ) : ℕ :=
  ∏ i, (β i).descFactorial (α i)

/-- Coefficientwise action of a constant-coefficient differential symbol. -/
noncomputable def differentialAction
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (symbol p : MvPolynomial σ K) : MvPolynomial σ K :=
  Finsupp.sum (AddMonoidAlgebra.coeff symbol) fun α a ↦
    Finsupp.sum (AddMonoidAlgebra.coeff p) fun β b ↦
      monomial (β - α) (a * b * (multiDescFactorial β α : K))

@[simp] theorem differentialAction_zero_right
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (symbol : MvPolynomial σ K) :
    differentialAction symbol 0 = 0 := by
  simp [differentialAction]

def PurePowersVanish
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (symbol p : MvPolynomial σ K) : Prop :=
  ∀ m : ℕ, 0 < m → differentialAction (symbol ^ m) (p ^ m) = 0

def EventuallyMixedPowersVanish
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (symbol p : MvPolynomial σ K) : Prop :=
  ∀ q : MvPolynomial σ K, ∃ M : ℕ, ∀ m ≥ M,
    differentialAction (symbol ^ m) (q * p ^ m) = 0

/-- GVC for one symbol/polynomial pair. -/
def GeneralizedVanishingFor
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (symbol p : MvPolynomial σ K) : Prop :=
  PurePowersVanish symbol p → EventuallyMixedPowersVanish symbol p

/-- GVC for all constant-coefficient symbols and polynomials with variable
index type `σ`. -/
def GeneralizedVanishingConjecture
    (σ K : Type*) [Fintype σ] [CommSemiring K] : Prop :=
  ∀ symbol p : MvPolynomial σ K,
    GeneralizedVanishingFor symbol p

end GVC
