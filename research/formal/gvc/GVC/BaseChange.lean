/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GVC.Padding

/-!
# Base change of differential counterexamples

Mapping coefficients through a ring homomorphism commutes with the
coefficientwise differential action.  An injective coefficient map therefore
preserves both pure vanishing and mixed nonvanishing.  This formalizes the
manuscript's descent of the rational ternary identities to every
characteristic-zero field.
-/

namespace GVC

open MvPolynomial

/-- Coefficient base change commutes with constant-coefficient
differentiation. -/
theorem differentialAction_map
    {σ K L : Type*} [Fintype σ] [CommSemiring K] [CommSemiring L]
    (φ : K →+* L) (symbol p : MvPolynomial σ K) :
    differentialAction (map φ symbol) (map φ p) =
      map φ (differentialAction symbol p) := by
  induction symbol using MvPolynomial.induction_on generalizing p with
  | C a => simp
  | add symbol₁ symbol₂ ih₁ ih₂ =>
      rw [map_add, differentialAction_add_left, differentialAction_add_left,
        ih₁, ih₂, map_add]
  | mul_X symbol i ih =>
      rw [map_mul, map_X, differentialAction_mul_left,
        differentialAction_X, pderiv_map, ih,
        differentialAction_mul_left, differentialAction_X]

/-- Pure-power vanishing is preserved by any coefficient base change. -/
theorem purePowersVanish_map
    {σ K L : Type*} [Fintype σ] [CommSemiring K] [CommSemiring L]
    (φ : K →+* L) {symbol p : MvPolynomial σ K}
    (hvanish : PurePowersVanish symbol p) :
    PurePowersVanish (map φ symbol) (map φ p) := by
  intro m hm
  calc
    differentialAction ((map φ symbol) ^ m) ((map φ p) ^ m) =
        differentialAction (map φ (symbol ^ m)) (map φ (p ^ m)) := by
          rw [map_pow, map_pow]
    _ = map φ (differentialAction (symbol ^ m) (p ^ m)) :=
      differentialAction_map φ _ _
    _ = 0 := by rw [hvanish m hm, map_zero]

/-- Everywhere mixed nonvanishing is preserved by an injective coefficient
base change. -/
theorem mixedPowers_ne_zero_map
    {σ K L : Type*} [Fintype σ] [CommSemiring K] [CommSemiring L]
    (φ : K →+* L) (hφ : Function.Injective φ)
    {symbol p q : MvPolynomial σ K}
    (hne : ∀ m : ℕ, 0 < m →
      differentialAction (symbol ^ m) (q * p ^ m) ≠ 0) :
    ∀ m : ℕ, 0 < m →
      differentialAction ((map φ symbol) ^ m)
        (map φ q * (map φ p) ^ m) ≠ 0 := by
  intro m hm hzero
  have hmapped :
      map φ (differentialAction (symbol ^ m) (q * p ^ m)) = 0 := by
    rw [← differentialAction_map φ]
    simpa only [map_pow, map_mul] using hzero
  apply hne m hm
  apply MvPolynomial.map_injective φ hφ
  simpa using hmapped

/-- An explicit counterexample remains one after an injective coefficient
base change. -/
theorem not_generalizedVanishingFor_map
    {σ K L : Type*} [Fintype σ] [CommSemiring K] [CommSemiring L]
    (φ : K →+* L) (hφ : Function.Injective φ)
    {symbol p q : MvPolynomial σ K}
    (hvanish : PurePowersVanish symbol p)
    (hne : ∀ m : ℕ, 0 < m →
      differentialAction (symbol ^ m) (q * p ^ m) ≠ 0) :
    ¬ GeneralizedVanishingFor (map φ symbol) (map φ p) := by
  intro hGVC
  obtain ⟨M, hM⟩ := hGVC (purePowersVanish_map φ hvanish) (map φ q)
  let m := max M 1
  have hmM : M ≤ m := le_max_left _ _
  have hm : 0 < m := lt_of_lt_of_le Nat.zero_lt_one (le_max_right _ _)
  exact mixedPowers_ne_zero_map φ hφ hne m hm (hM m hmM)

section CharacteristicZero

variable {K : Type*} [Field K] [CharZero K]

private noncomputable abbrev rationalBaseChange (K : Type*)
    [Field K] [CharZero K] : ℚ →+* K :=
  Rat.castHom K

private theorem rationalBaseChange_injective (K : Type*)
    [Field K] [CharZero K] :
    Function.Injective (rationalBaseChange K) :=
  RingHom.injective _

/-- The exact nonzero scalar identity of Theorem 8.1 after base change from
`ℚ` to an arbitrary characteristic-zero field. -/
theorem gvc3_exact_next_mixed_charZero
    (B : ConcreteCounterexampleBridge) (m : ℕ) (hm : 0 < m) :
    differentialAction (map (rationalBaseChange K) gvcDelta)
      (differentialAction ((map (rationalBaseChange K) gvcLambda) ^ m)
        (map (rationalBaseChange K) gvcQ *
          (map (rationalBaseChange K) gvcP) ^ m)) =
      C (rationalBaseChange K (mixedDerivativeValue m)) := by
  let φ := rationalBaseChange K
  have hinner :
      differentialAction ((map φ gvcLambda) ^ m)
          (map φ gvcQ * (map φ gvcP) ^ m) =
        map φ (differentialAction (gvcLambda ^ m) (gvcQ * gvcP ^ m)) := by
    simpa only [map_pow, map_mul] using
      differentialAction_map φ (gvcLambda ^ m) (gvcQ * gvcP ^ m)
  rw [hinner, differentialAction_map,
    gvc3_exact_next_mixed B m hm, map_C]

/-- The pure identity of Theorem 8.1 after base change from `ℚ` to an
arbitrary characteristic-zero field. -/
theorem gvc3_pure_identity_charZero
    (B : ConcreteCounterexampleBridge) (m : ℕ) (hm : 0 < m) :
    differentialAction ((map (rationalBaseChange K) gvcLambda) ^ m)
      ((map (rationalBaseChange K) gvcP) ^ m) = 0 := by
  let φ := rationalBaseChange K
  calc
    differentialAction ((map φ gvcLambda) ^ m) ((map φ gvcP) ^ m) =
        differentialAction (map φ (gvcLambda ^ m))
          (map φ (gvcP ^ m)) := by rw [map_pow, map_pow]
    _ = map φ (differentialAction (gvcLambda ^ m) (gvcP ^ m)) :=
      differentialAction_map φ _ _
    _ = 0 := by rw [gvc3_pure_identity B m hm, map_zero]

/-- The manuscript's mixed output is nonzero over every
characteristic-zero field. -/
theorem gvc3_mixed_ne_zero_charZero
    (B : ConcreteCounterexampleBridge) (m : ℕ) (hm : 0 < m) :
    differentialAction ((map (rationalBaseChange K) gvcLambda) ^ m)
      (map (rationalBaseChange K) gvcQ *
        (map (rationalBaseChange K) gvcP) ^ m) ≠ 0 :=
  mixedPowers_ne_zero_map (rationalBaseChange K)
    (rationalBaseChange_injective K) (gvc3_mixed_ne_zero B) m hm

/-- The ternary counterexample descends from `ℚ` to every
characteristic-zero field. -/
theorem gvc3_charZero_not_generalizedVanishingFor
    (B : ConcreteCounterexampleBridge) :
    ¬ GeneralizedVanishingFor
      (map (rationalBaseChange K) gvcLambda)
      (map (rationalBaseChange K) gvcP) :=
  not_generalizedVanishingFor_map (rationalBaseChange K)
    (rationalBaseChange_injective K) (gvc3_purePowersVanish B)
    (gvc3_mixed_ne_zero B)

/-- Combining coefficient base change and unused-variable padding transports
any concrete bridge to every characteristic-zero field and every finite
dimension `n ≥ 3`. -/
theorem gvc3_charZero_fin_padded_not_generalizedVanishingFor
    (B : ConcreteCounterexampleBridge) {n : ℕ} (hn : 3 ≤ n) :
    ¬ GeneralizedVanishingFor
      (rename (Fin.castLE hn) (map (rationalBaseChange K) gvcLambda))
      (rename (Fin.castLE hn) (map (rationalBaseChange K) gvcP)) :=
  not_generalizedVanishingFor_rename (Fin.castLE hn)
    (Fin.castLE_injective hn)
    (purePowersVanish_map (rationalBaseChange K) (gvc3_purePowersVanish B))
    (mixedPowers_ne_zero_map (rationalBaseChange K)
      (rationalBaseChange_injective K) (gvc3_mixed_ne_zero B))

/-- Generic bridge-parametrized form of the negative half of the paper's
dimension classification.  `GVC.VerifiedCounterexample` instantiates the
bridge unconditionally. -/
theorem gvc3_charZero_fin_not_generalizedVanishingConjecture
    (B : ConcreteCounterexampleBridge) {n : ℕ} (hn : 3 ≤ n) :
    ¬ GeneralizedVanishingConjecture (Fin n) K := by
  intro hGVC
  exact gvc3_charZero_fin_padded_not_generalizedVanishingFor B hn
    (hGVC _ _)

end CharacteristicZero

end GVC
