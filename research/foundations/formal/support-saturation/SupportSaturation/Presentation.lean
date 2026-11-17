/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn, OpenAI Codex
-/

import SupportSaturation.Core

/-!
# Presentation saturation

This file identifies vanishing of the primary component of a quotient
`F / N` with the elementwise saturation statement used by finite free
presentations.
-/

open Submodule

namespace SupportSaturation

variable {R F : Type*} [CommRing R] [AddCommGroup F] [Module R F]

/-- Elementwise saturation of a submodule by all powers of an ideal. -/
def IsSaturated (I : Ideal R) (N : Submodule R F) : Prop :=
  ∀ x : F, (∃ n : ℕ, ∀ a ∈ I ^ n, a • x ∈ N) → x ∈ N

/-- If the `I`-primary component of `F / N` vanishes, every representative
whose `I`-power multiples lie in `N` already lies in `N`. -/
theorem mem_of_pow_smul_mem_of_primaryComponent_eq_bot
    (I : Ideal R) (N : Submodule R F) {x : F}
    (hprimary : I.primaryComponent (F ⧸ N) = ⊥)
    (hcolon : ∃ n : ℕ, ∀ a ∈ I ^ n, a • x ∈ N) :
    x ∈ N := by
  obtain ⟨n, hn⟩ := hcolon
  have hx : N.mkQ x ∈ I.primaryComponent (F ⧸ N) := by
    rw [Ideal.primaryComponent_mem]
    refine ⟨n, ?_⟩
    rw [Submodule.mem_torsionBySet_iff]
    intro a
    change N.mkQ ((a : R) • x) = 0
    rw [Submodule.mkQ_apply, Submodule.Quotient.mk_eq_zero]
    exact hn a a.2
  rw [hprimary] at hx
  exact (Submodule.Quotient.mk_eq_zero N).mp hx

/-- Quotient primary-component vanishing is exactly submodule saturation. -/
theorem primaryComponent_eq_bot_iff_isSaturated
    (I : Ideal R) (N : Submodule R F) :
    I.primaryComponent (F ⧸ N) = ⊥ ↔ IsSaturated I N := by
  constructor
  · intro hprimary x hcolon
    exact mem_of_pow_smul_mem_of_primaryComponent_eq_bot I N hprimary hcolon
  · intro hSaturated
    rw [Submodule.eq_bot_iff]
    intro y hy
    obtain ⟨x, rfl⟩ := N.mkQ_surjective y
    apply (Submodule.Quotient.mk_eq_zero N).mpr
    apply hSaturated x
    rw [Ideal.primaryComponent_mem] at hy
    obtain ⟨n, hn⟩ := hy
    refine ⟨n, fun a ha ↦ ?_⟩
    rw [Submodule.mem_torsionBySet_iff] at hn
    have hzero : (a : R) • N.mkQ x = 0 := hn ⟨a, ha⟩
    change N.mkQ ((a : R) • x) = 0 at hzero
    exact (Submodule.Quotient.mk_eq_zero N).mp hzero

/-- The associated-prime form of `S₁`, together with support avoidance,
implies presentation saturation. -/
theorem mem_of_pow_smul_mem_of_noEmbeddedAssociatedPrimes
    [IsNoetherianRing R] (N : Submodule R F) [Module.Finite R (F ⧸ N)]
    (I : Ideal R) {x : F}
    (hnoEmbedded : NoEmbeddedAssociatedPrimes (R := R) (M := F ⧸ N))
    (hheight : AvoidsMinimalSupport (R := R) (M := F ⧸ N) I)
    (hcolon : ∃ n : ℕ, ∀ a ∈ I ^ n, a • x ∈ N) :
    x ∈ N := by
  apply mem_of_pow_smul_mem_of_primaryComponent_eq_bot I N
    (primaryComponent_eq_bot_of_noEmbeddedAssociatedPrimes I hnoEmbedded hheight)
    hcolon

end SupportSaturation
