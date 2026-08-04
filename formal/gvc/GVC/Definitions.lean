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

@[simp] theorem differentialAction_add_right
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (symbol p q : MvPolynomial σ K) :
    differentialAction symbol (p + q) =
      differentialAction symbol p + differentialAction symbol q := by
  classical
  simp [differentialAction, Finsupp.sum_add_index, add_mul, mul_add]

/-- The coefficientwise semantics agrees with Mathlib's formal partial
derivative for a degree-one symbol. -/
theorem differentialAction_X
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (i : σ) (p : MvPolynomial σ K) :
    differentialAction (X i) p = pderiv i p := by
  classical
  induction p using MvPolynomial.induction_on' with
  | monomial β b =>
      have hprod :
          (∏ j, ((β j).descFactorial ((Finsupp.single i 1) j) : K)) =
            (β i : K) := by
        rw [← Nat.cast_prod]
        congr
        rw [Finset.prod_eq_single i]
        · simp
        · intro j _hj hji
          simp [Finsupp.single_eq_of_ne hji]
        · simp
      simpa [differentialAction, multiDescFactorial, X,
        pderiv_monomial] using
        congrArg
          (fun z ↦ monomial (β - Finsupp.single i 1) (b * z)) hprod
  | add p q hp hq =>
      rw [differentialAction_add_right, map_add, hp, hq]

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
