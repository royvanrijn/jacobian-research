/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn, OpenAI Codex
-/

import Mathlib.Algebra.Module.Torsion.PrimaryComponent
import Mathlib.RingTheory.Ideal.AssociatedPrime.Finiteness
import Mathlib.RingTheory.Ideal.AssociatedPrime.Localization
import Mathlib.RingTheory.Regular.IsSMulRegular

/-!
# Support saturation

This file formalizes the module-theoretic kernel of the repository's
support-saturation theorem.  Mathlib's `Ideal.primaryComponent M I` is

`⋃ n, 0 :_M I ^ n`,

the zeroth local-cohomology module `H⁰_I(M)`.

The geometric `S₁` hypothesis is represented by its associated-prime
consequence: every associated prime is minimal over `ann(M)`.  This avoids
introducing a local-depth API solely for this theorem.
-/

open Set Submodule Module

namespace SupportSaturation

variable {R M : Type*} [CommRing R] [AddCommGroup M] [Module R M]

/-- A module has no embedded associated primes when every associated prime
is minimal over its annihilator.  For finite modules over Noetherian rings,
this is the associated-prime consequence of Serre's condition `S₁`. -/
def NoEmbeddedAssociatedPrimes : Prop :=
  associatedPrimes R M ⊆ (Module.annihilator R M).minimalPrimes

/-- The support ideal avoids every irreducible component of `Supp(M)`. -/
def AvoidsMinimalSupport (I : Ideal R) : Prop :=
  ∀ p ∈ (Module.annihilator R M).minimalPrimes, ¬I ≤ p

/-- A regular element in `I` kills all `I`-power torsion. -/
theorem primaryComponent_eq_bot_of_regular (I : Ideal R) {r : R}
    (hrI : r ∈ I) (hr : IsSMulRegular M r) :
    I.primaryComponent M = ⊥ := by
  rw [Submodule.eq_bot_iff]
  intro x hx
  rw [Ideal.primaryComponent_mem] at hx
  obtain ⟨n, hn⟩ := hx
  have hrpow : r ^ n ∈ I ^ n := Ideal.pow_mem_pow hrI n
  rw [Submodule.mem_torsionBySet_iff] at hn
  have hzero : r ^ n • x = 0 := hn ⟨r ^ n, hrpow⟩
  exact (hr.pow n).right_eq_zero_of_smul hzero

/-- Avoidance of all associated primes produces an `M`-regular element of
the ideal.  This is the prime-avoidance step in support saturation. -/
theorem exists_regular_of_avoids_associatedPrimes [IsNoetherianRing R]
    [Module.Finite R M] (I : Ideal R)
    (hAvoid : ∀ p ∈ associatedPrimes R M, ¬I ≤ p) :
    ∃ r ∈ I, IsSMulRegular M r := by
  cases subsingleton_or_nontrivial M with
  | inl hM =>
      letI : Subsingleton M := hM
      exact ⟨0, I.zero_mem, IsSMulRegular.zero⟩
  | inr hM =>
      letI : Nontrivial M := hM
      by_contra! hregular
      have hsubset : (I : Set R) ⊆ ⋃ p ∈ associatedPrimes R M, p := by
        rw [biUnion_associatedPrimes_eq_compl_regular R M]
        exact fun r hr ↦ hregular r hr
      obtain ⟨p, hp, hIp⟩ :=
        (I.subset_union_prime_finite (associatedPrimes.finite R M)
          (f := id) 0 0 (fun _ hp _ _ ↦ hp.isPrime)).mp hsubset
      exact hAvoid p hp hIp

/-- Associated-prime avoidance is exactly the hypothesis needed for the
zeroth local-cohomology module to vanish. -/
theorem primaryComponent_eq_bot_of_avoids_associatedPrimes
    [IsNoetherianRing R] [Module.Finite R M] (I : Ideal R)
    (hAvoid : ∀ p ∈ associatedPrimes R M, ¬I ≤ p) :
    I.primaryComponent M = ⊥ := by
  obtain ⟨r, hrI, hr⟩ := exists_regular_of_avoids_associatedPrimes I hAvoid
  exact primaryComponent_eq_bot_of_regular I hrI hr

/-- `S₁` in associated-prime form plus positive relative support height
forces support saturation. -/
theorem primaryComponent_eq_bot_of_noEmbeddedAssociatedPrimes
    [IsNoetherianRing R] [Module.Finite R M] (I : Ideal R)
    (hnoEmbedded : NoEmbeddedAssociatedPrimes (R := R) (M := M))
    (hheight : AvoidsMinimalSupport (R := R) (M := M) I) :
    I.primaryComponent M = ⊥ := by
  apply primaryComponent_eq_bot_of_avoids_associatedPrimes I
  intro p hp
  exact hheight p (hnoEmbedded hp)

/-- A nonzero vector annihilated by every element of `I` is a concrete
obstruction to support saturation.  This is the formal kernel of the
counterexample `R \oplus R/I`: the class of `1` in the quotient summand is
nonzero and is annihilated by `I`. -/
theorem primaryComponent_ne_bot_of_annihilated
    (I : Ideal R) {x : M} (hx : x ≠ 0)
    (hIx : ∀ a ∈ I, a • x = 0) :
    I.primaryComponent M ≠ ⊥ := by
  intro hbot
  have hxPrimary : x ∈ I.primaryComponent M := by
    rw [Ideal.primaryComponent_mem]
    refine ⟨1, ?_⟩
    rw [Submodule.mem_torsionBySet_iff]
    intro a
    exact hIx a.1 (by simpa using a.2)
  rw [hbot] at hxPrimary
  exact hx hxPrimary

/-- Vanishing of the primary component is equivalent to avoidance of every
associated prime.  Together with finite prime avoidance, this is the exact
associated-prime form of the support-saturation principle. -/
theorem primaryComponent_eq_bot_iff_avoids_associatedPrimes
    [IsNoetherianRing R] [Module.Finite R M] (I : Ideal R) :
    I.primaryComponent M = ⊥ ↔
      ∀ p ∈ associatedPrimes R M, ¬I ≤ p := by
  constructor
  · intro hprimary p hp hIp
    rw [AssociatedPrimes.mem_iff, isAssociatedPrime_iff] at hp
    obtain ⟨hpPrime, x, hpx⟩ := hp
    have hx : x ≠ 0 := by
      intro hx
      apply hpPrime.ne_top
      rw [hpx, hx]
      simp
    have hIx : ∀ a ∈ I, a • x = 0 := fun a ha ↦ by
      have haColon : a ∈ (⊥ : Submodule R M).colon {x} := by
        rw [← hpx]
        exact hIp ha
      simpa [Submodule.mem_colon_singleton] using haColon
    exact (primaryComponent_ne_bot_of_annihilated I hx hIx) hprimary
  · exact primaryComponent_eq_bot_of_avoids_associatedPrimes I

/-- Equivalently, the primary component vanishes exactly when the ideal
contains a module-regular element. -/
theorem primaryComponent_eq_bot_iff_exists_regular
    [IsNoetherianRing R] [Module.Finite R M] (I : Ideal R) :
    I.primaryComponent M = ⊥ ↔
      ∃ r ∈ I, IsSMulRegular M r := by
  constructor
  · intro hprimary
    exact exists_regular_of_avoids_associatedPrimes I
      ((primaryComponent_eq_bot_iff_avoids_associatedPrimes I).mp hprimary)
  · rintro ⟨r, hrI, hr⟩
    exact primaryComponent_eq_bot_of_regular I hrI hr

end SupportSaturation
