/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GVC.Definitions
import Mathlib.RingTheory.MvPolynomial.Homogeneous

/-!
# Top homogeneous contractions

Equal-degree homogeneous symbols and polynomials contract to a constant.
This file proves the exact coefficientwise formula used by the algebraic
Reynolds functional in the ternary counterexample.
-/

namespace GVC

open MvPolynomial

/-- The scalar obtained by pairing equal multi-indices with their full
descending-factorial weight. -/
noncomputable def topContraction
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (symbol p : MvPolynomial σ K) : K :=
  Finsupp.sum (AddMonoidAlgebra.coeff symbol) fun α a ↦
    a * p.coeff α * (multiDescFactorial α α : K)

@[simp] theorem differentialAction_zero_left
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (p : MvPolynomial σ K) :
    differentialAction 0 p = 0 := by
  simp [differentialAction]

@[simp] theorem differentialAction_add_left
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (symbol₁ symbol₂ p : MvPolynomial σ K) :
    differentialAction (symbol₁ + symbol₂) p =
      differentialAction symbol₁ p + differentialAction symbol₂ p := by
  classical
  simp [differentialAction, Finsupp.sum_add_index, add_mul]

@[simp] theorem topContraction_zero_left
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (p : MvPolynomial σ K) :
    topContraction 0 p = 0 := by
  simp [topContraction]

@[simp] theorem topContraction_zero_right
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (symbol : MvPolynomial σ K) :
    topContraction symbol 0 = 0 := by
  simp [topContraction]

@[simp] theorem topContraction_add_left
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (symbol₁ symbol₂ p : MvPolynomial σ K) :
    topContraction (symbol₁ + symbol₂) p =
      topContraction symbol₁ p + topContraction symbol₂ p := by
  classical
  simp [topContraction, Finsupp.sum_add_index, add_mul]

@[simp] theorem topContraction_add_right
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (symbol p₁ p₂ : MvPolynomial σ K) :
    topContraction symbol (p₁ + p₂) =
      topContraction symbol p₁ + topContraction symbol p₂ := by
  classical
  simp [topContraction, add_mul, mul_add]

@[simp] theorem topContraction_monomial_left
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (α : σ →₀ ℕ) (a : K) (p : MvPolynomial σ K) :
    topContraction (monomial α a) p =
      a * p.coeff α * (multiDescFactorial α α : K) := by
  classical
  by_cases ha : a = 0
  · subst a
    simp
  · simp [topContraction]

theorem multiDescFactorial_eq_zero_of_not_le
    {σ : Type*} [Fintype σ] (β α : σ →₀ ℕ) (hnot : ¬α ≤ β) :
    multiDescFactorial β α = 0 := by
  simp only [Finsupp.le_iff] at hnot
  push Not at hnot
  obtain ⟨i, hi⟩ := hnot
  rw [multiDescFactorial, Finset.prod_eq_zero (Finset.mem_univ i)]
  exact Nat.descFactorial_eq_zero_iff_lt.mpr hi.2

@[simp] theorem multiDescFactorial_self
    {σ : Type*} [Fintype σ] (α : σ →₀ ℕ) :
    multiDescFactorial α α = ∏ i, Nat.factorial (α i) := by
  simp only [multiDescFactorial, Nat.descFactorial_self]

private theorem finsupp_eq_of_le_of_degree_eq
    {σ : Type*} [Finite σ] {α β : σ →₀ ℕ}
    (hle : α ≤ β) (hdegree : α.degree = β.degree) : α = β := by
  classical
  letI := Fintype.ofFinite σ
  apply Finsupp.ext
  intro i
  have hpoint : α i ≤ β i := hle i
  by_contra hne
  have hstrict : α i < β i := lt_of_le_of_ne hpoint hne
  have hsum : (∑ j, α j) < ∑ j, β j :=
    Finset.sum_lt_sum (fun j _ ↦ hle j) ⟨i, Finset.mem_univ i, hstrict⟩
  have halpha : α.degree = ∑ j, α j := by
    change α.sum (fun _ c ↦ c) = _
    exact Finsupp.sum_fintype _ _ (fun _ ↦ rfl)
  have hbeta : β.degree = ∑ j, β j := by
    change β.sum (fun _ c ↦ c) = _
    exact Finsupp.sum_fintype _ _ (fun _ ↦ rfl)
  omega

theorem differentialAction_monomial_monomial
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (α β : σ →₀ ℕ) (a b : K) :
    differentialAction (monomial α a) (monomial β b) =
      monomial (β - α) (a * b * (multiDescFactorial β α : K)) := by
  classical
  simp [differentialAction]

theorem multiDescFactorial_add
    {σ : Type*} [Fintype σ] (β γ α : σ →₀ ℕ) :
    multiDescFactorial β (γ + α) =
      multiDescFactorial β γ * multiDescFactorial (β - γ) α := by
  classical
  rw [multiDescFactorial, multiDescFactorial, multiDescFactorial,
    ← Finset.prod_mul_distrib]
  apply Finset.prod_congr rfl
  intro i _hi
  simpa [mul_comm] using
    (Nat.descFactorial_mul_descFactorial
      (n := β i) (k := γ i) (m := γ i + α i) (by omega)).symm

private theorem finsupp_sub_sub
    {σ : Type*} (β γ α : σ →₀ ℕ) :
    (β - γ) - α = β - (γ + α) := by
  ext i
  change (β i - γ i) - α i = β i - (γ i + α i)
  omega

theorem differentialAction_monomial_mul_comp
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (α γ β : σ →₀ ℕ) (a c b : K) :
    differentialAction (monomial α a * monomial γ c) (monomial β b) =
      differentialAction (monomial α a)
        (differentialAction (monomial γ c) (monomial β b)) := by
  classical
  rw [monomial_mul, differentialAction_monomial_monomial,
    differentialAction_monomial_monomial,
    differentialAction_monomial_monomial, finsupp_sub_sub]
  have hexponent : β - (α + γ) = β - (γ + α) := by
    rw [add_comm α γ]
  have hfactor :
      multiDescFactorial β (α + γ) =
        multiDescFactorial β γ * multiDescFactorial (β - γ) α := by
    rw [add_comm α γ, multiDescFactorial_add]
  rw [hexponent, hfactor]
  push_cast
  congr 1
  ring

/-- Multiplication of constant-coefficient symbols agrees with composition
of their differential actions. -/
theorem differentialAction_mul_left
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (symbol₁ symbol₂ p : MvPolynomial σ K) :
    differentialAction (symbol₁ * symbol₂) p =
      differentialAction symbol₁ (differentialAction symbol₂ p) := by
  classical
  induction symbol₁ using MvPolynomial.induction_on' with
  | add symbol₁ symbol₂ ih₁ ih₂ =>
      rw [add_mul, differentialAction_add_left,
        differentialAction_add_left, ih₁, ih₂]
  | monomial α a =>
      induction symbol₂ using MvPolynomial.induction_on' with
      | add symbol₂ symbol₃ ih₂ ih₃ =>
          rw [mul_add, differentialAction_add_left,
            differentialAction_add_left, differentialAction_add_right,
            ih₂, ih₃]
      | monomial γ c =>
          induction p using MvPolynomial.induction_on' with
          | add p q ihp ihq =>
              rw [differentialAction_add_right, differentialAction_add_right,
                differentialAction_add_right, ihp, ihq]
          | monomial β b =>
              exact differentialAction_monomial_mul_comp α γ β a c b

theorem differentialAction_monomial_monomial_of_degree_eq
    {σ K : Type*} [Fintype σ] [Field K] [CharZero K]
    (α β : σ →₀ ℕ) (a b : K) (hdegree : α.degree = β.degree) :
    differentialAction (monomial α a) (monomial β b) =
      C (topContraction (monomial α a) (monomial β b)) := by
  classical
  by_cases hab : α = β
  · subst β
    rw [differentialAction_monomial_monomial,
      topContraction_monomial_left, coeff_monomial]
    simp
  · have hnle : ¬α ≤ β := fun hle ↦
      hab (finsupp_eq_of_le_of_degree_eq hle hdegree)
    have hzero := multiDescFactorial_eq_zero_of_not_le β α hnle
    rw [differentialAction_monomial_monomial,
      topContraction_monomial_left, coeff_monomial]
    have hba : β ≠ α := Ne.symm hab
    simp [hba, hzero]

/-- An equal-degree homogeneous constant-coefficient contraction is the
constant obtained by pairing equal monomial exponents. -/
theorem differentialAction_eq_C_topContraction_of_homogeneous
    {σ K : Type*} [Fintype σ] [Field K] [CharZero K]
    {symbol p : MvPolynomial σ K} {d : ℕ}
    (hsymbol : symbol.IsHomogeneous d) (hp : p.IsHomogeneous d) :
    differentialAction symbol p = C (topContraction symbol p) := by
  induction hsymbol using IsWeightedHomogeneous.induction_on with
  | zero => simp
  | add symbol₁ symbol₂ hs₁ hs₂ ih₁ ih₂ =>
      rw [differentialAction_add_left, ih₁, ih₂,
        topContraction_add_left, C_add]
  | monomial α a hα =>
      induction hp using IsWeightedHomogeneous.induction_on with
      | zero => simp
      | add p₁ p₂ hp₁ hp₂ ih₁ ih₂ =>
          rw [differentialAction_add_right, ih₁, ih₂,
            topContraction_add_right, C_add]
      | monomial β b hβ =>
          apply differentialAction_monomial_monomial_of_degree_eq
          calc
            α.degree = (Finsupp.weight (1 : σ → ℕ)) α :=
              congrArg (fun f ↦ f α) Finsupp.degree_eq_weight_one
            _ = d := hα
            _ = (Finsupp.weight (1 : σ → ℕ)) β := hβ.symm
            _ = β.degree :=
              (congrArg (fun f ↦ f β) Finsupp.degree_eq_weight_one).symm

end GVC
