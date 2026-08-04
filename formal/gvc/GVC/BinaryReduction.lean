/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GVC.Definitions

/-!
# Binary GVC reduction surface

This file records the exact logical interface between Sections 4--6 of the
manuscript.  It does not assume the binary theorem as one opaque proposition:
the two remaining obligations are separately named as envelope closure and
terminality of a common threshold.
-/

namespace GVC

open MvPolynomial

abbrev BinaryPolynomial (K : Type*) [CommSemiring K] :=
  MvPolynomial (Fin 2) K

def binaryWeight (w : Fin 2 → ℕ) (α : Fin 2 →₀ ℕ) : ℕ :=
  ∑ i, w i * α i

def PositiveUnequalWeight (w : Fin 2 → ℕ) : Prop :=
  (∀ i, 0 < w i) ∧ w 0 ≠ w 1

/-- The support inequalities at the first common envelope threshold. -/
def HasCommonThreshold
    {K : Type*} [CommSemiring K]
    (symbol p : BinaryPolynomial K) (w : Fin 2 → ℕ) (W : ℕ) : Prop :=
  PositiveUnequalWeight w ∧
    (∀ α ∈ symbol.support, W ≤ binaryWeight w α) ∧
    (∀ β ∈ p.support, binaryWeight w β ≤ W)

/-- The two proof obligations left by the current Lean development of the
binary theorem. -/
structure BinaryEnvelopeBridge
    (K : Type*) [Field K] [CharZero K] where
  /-- Sections 3, 4, and 6: Hall localization plus shifted-ray separation
  force the moving envelopes to reach a common threshold. -/
  envelope_closure : ∀ symbol p : BinaryPolynomial K,
    PurePowersVanish symbol p →
    ∃ w W, HasCommonThreshold symbol p w W
  /-- Section 5: bounded weight defect makes a common unequal threshold
  terminal for every fixed multiplier. -/
  common_threshold_terminal : ∀ symbol p : BinaryPolynomial K,
    PurePowersVanish symbol p →
    (∃ w W, HasCommonThreshold symbol p w W) →
    EventuallyMixedPowersVanish symbol p

/-- The binary theorem follows from precisely the envelope and terminality
obligations displayed above. -/
theorem binary_gvc_of_envelope_bridge
    {K : Type*} [Field K] [CharZero K]
    (B : BinaryEnvelopeBridge K) :
    GeneralizedVanishingConjecture (Fin 2) K := by
  intro symbol p hpure
  exact B.common_threshold_terminal symbol p hpure
    (B.envelope_closure symbol p hpure)

end GVC
