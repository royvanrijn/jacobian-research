/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GVC.ConcreteWitness

/-!
# Unused-variable padding

Injectively renaming the variables of both a constant-coefficient symbol and
its input polynomial commutes with the differential action.  Consequently a
counterexample witnessed by nonvanishing mixed powers remains a
counterexample after adjoining unused variables.
-/

namespace GVC

open MvPolynomial

@[simp] theorem differentialAction_C
    {σ K : Type*} [Fintype σ] [CommSemiring K]
    (a : K) (p : MvPolynomial σ K) :
    differentialAction (C a) p = C a * p := by
  classical
  induction p using MvPolynomial.induction_on' with
  | add p q hp hq =>
      rw [differentialAction_add_right, hp, hq, mul_add]
  | monomial β b =>
      rw [show C a = monomial 0 a by rfl,
        differentialAction_monomial_monomial]
      simpa [multiDescFactorial] using
        (MvPolynomial.C_mul_monomial (a := a) (a' := b) (s := β)).symm

/-- Injectively adjoining or renaming variables commutes with the full
coefficientwise differential action. -/
theorem differentialAction_rename
    {σ τ K : Type*} [Fintype σ] [Fintype τ] [CommSemiring K]
    (f : σ → τ) (hf : Function.Injective f)
    (symbol p : MvPolynomial σ K) :
    differentialAction (rename f symbol) (rename f p) =
      rename f (differentialAction symbol p) := by
  induction symbol using MvPolynomial.induction_on generalizing p with
  | C a => simp
  | add symbol₁ symbol₂ ih₁ ih₂ =>
      rw [map_add, differentialAction_add_left, differentialAction_add_left,
        ih₁, ih₂, map_add]
  | mul_X symbol i ih =>
      rw [map_mul, rename_X, differentialAction_mul_left,
        differentialAction_X, pderiv_rename hf, ih,
        differentialAction_mul_left, differentialAction_X]

/-- Pure-power vanishing is preserved when unused variables are adjoined. -/
theorem purePowersVanish_rename
    {σ τ K : Type*} [Fintype σ] [Fintype τ] [CommSemiring K]
    (f : σ → τ) (hf : Function.Injective f)
    {symbol p : MvPolynomial σ K}
    (hvanish : PurePowersVanish symbol p) :
    PurePowersVanish (rename f symbol) (rename f p) := by
  intro m hm
  calc
    differentialAction ((rename f symbol) ^ m) ((rename f p) ^ m) =
        differentialAction (rename f (symbol ^ m)) (rename f (p ^ m)) := by
          rw [map_pow, map_pow]
    _ = rename f (differentialAction (symbol ^ m) (p ^ m)) :=
      differentialAction_rename f hf _ _
    _ = 0 := by rw [hvanish m hm, map_zero]

/-- GVC descends from a variable type to any type embedded in it. -/
theorem generalizedVanishingConjecture_of_injective
    {σ τ K : Type*} [Fintype σ] [Fintype τ] [CommSemiring K]
    (f : σ → τ) (hf : Function.Injective f)
    (hGVC : GeneralizedVanishingConjecture τ K) :
    GeneralizedVanishingConjecture σ K := by
  intro symbol p hvanish q
  obtain ⟨M, hM⟩ :=
    hGVC (rename f symbol) (rename f p)
      (purePowersVanish_rename f hf hvanish) (rename f q)
  refine ⟨M, fun m hm ↦ ?_⟩
  apply rename_injective f hf
  rw [← differentialAction_rename f hf]
  simpa only [map_pow, map_mul, map_zero] using hM m hm

/-- In particular, the binary theorem contains the one-variable theorem. -/
theorem unary_gvc_of_binary
    {K : Type*} [CommSemiring K]
    (hGVC : GeneralizedVanishingConjecture (Fin 2) K) :
    GeneralizedVanishingConjecture (Fin 1) K :=
  generalizedVanishingConjecture_of_injective
    (Fin.castLE (by omega : 1 ≤ 2))
    (Fin.castLE_injective (by omega : 1 ≤ 2)) hGVC

/-- A mixed-power sequence that is nonzero at every positive index remains
nonzero after unused-variable padding. -/
theorem mixedPowers_ne_zero_rename
    {σ τ K : Type*} [Fintype σ] [Fintype τ] [CommSemiring K]
    (f : σ → τ) (hf : Function.Injective f)
    {symbol p q : MvPolynomial σ K}
    (hne : ∀ m : ℕ, 0 < m →
      differentialAction (symbol ^ m) (q * p ^ m) ≠ 0) :
    ∀ m : ℕ, 0 < m →
      differentialAction ((rename f symbol) ^ m)
        (rename f q * (rename f p) ^ m) ≠ 0 := by
  intro m hm hzero
  have hpadded :
      rename f (differentialAction (symbol ^ m) (q * p ^ m)) = 0 := by
    rw [← differentialAction_rename f hf]
    simpa only [map_pow, map_mul] using hzero
  exact hne m hm ((rename_injective f hf) hpadded)

/-- An explicit pure-vanishing/mixed-nonvanishing counterexample remains a
counterexample after injectively adjoining variables. -/
theorem not_generalizedVanishingFor_rename
    {σ τ K : Type*} [Fintype σ] [Fintype τ] [CommSemiring K]
    (f : σ → τ) (hf : Function.Injective f)
    {symbol p q : MvPolynomial σ K}
    (hvanish : PurePowersVanish symbol p)
    (hne : ∀ m : ℕ, 0 < m →
      differentialAction (symbol ^ m) (q * p ^ m) ≠ 0) :
    ¬ GeneralizedVanishingFor (rename f symbol) (rename f p) := by
  intro hGVC
  obtain ⟨M, hM⟩ := hGVC (purePowersVanish_rename f hf hvanish) (rename f q)
  let m := max M 1
  have hmM : M ≤ m := le_max_left _ _
  have hm : 0 < m := lt_of_lt_of_le Nat.zero_lt_one (le_max_right _ _)
  exact mixedPowers_ne_zero_rename f hf hne m hm (hM m hmM)

/-- The manuscript's ternary witness, conditional only on phase extraction,
continues to refute GVC after any injective variable embedding. -/
theorem gvc3_padded_not_generalizedVanishingFor
    (B : ConcreteCounterexampleBridge)
    {τ : Type*} [Fintype τ] (f : Fin 3 → τ)
    (hf : Function.Injective f) :
    ¬ GeneralizedVanishingFor (rename f gvcLambda) (rename f gvcP) :=
  not_generalizedVanishingFor_rename f hf (gvc3_purePowersVanish B)
    (gvc3_mixed_ne_zero B)

/-- In particular, the same conditional witness exists in every finite
dimension `n ≥ 3`. -/
theorem gvc3_fin_padded_not_generalizedVanishingFor
    (B : ConcreteCounterexampleBridge) {n : ℕ} (hn : 3 ≤ n) :
    ¬ GeneralizedVanishingFor
      (rename (Fin.castLE hn) gvcLambda)
      (rename (Fin.castLE hn) gvcP) :=
  gvc3_padded_not_generalizedVanishingFor B (Fin.castLE hn)
    (Fin.castLE_injective hn)

/-- Thus, conditional only on phase extraction, GVC itself fails over `ℚ`
in every dimension `n ≥ 3`. -/
theorem gvc3_fin_not_generalizedVanishingConjecture
    (B : ConcreteCounterexampleBridge) {n : ℕ} (hn : 3 ≤ n) :
    ¬ GeneralizedVanishingConjecture (Fin n) ℚ := by
  intro hGVC
  exact gvc3_fin_padded_not_generalizedVanishingFor B hn
    (hGVC _ _)

end GVC
