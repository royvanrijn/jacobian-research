/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import DiscriminantPencils.LocalSingularities
import Mathlib.FieldTheory.IsAlgClosed.Basic
import Mathlib.FieldTheory.Perfect

/-!
# Counting cusp parameters

Over an algebraically closed characteristic-zero field, a separable second
derivative of a degree-`n` polynomial has exactly `n - 2` distinct roots.
Each root is an ordinary cusp parameter.
-/

open Polynomial

namespace DiscriminantPencils

variable {𝕜 : Type*} [Field 𝕜] [CharZero 𝕜] [IsAlgClosed 𝕜]

/-- The multiset of potential cusp parameters, namely the roots of `H''`.
Under the separability hypothesis below, it has no repetitions. -/
noncomputable def cuspParameters (H : 𝕜[X]) : Multiset 𝕜 :=
  H.derivative.derivative.roots

omit [CharZero 𝕜] [IsAlgClosed 𝕜] in
theorem ordinaryCuspParameter_iff (H : 𝕜[X])
    (hsep : H.derivative.derivative.Separable) (r : 𝕜) :
    IsOrdinaryCuspParameter H r ↔ H.derivative.derivative.eval r = 0 := by
  constructor
  · exact fun hr ↦ hr.1
  · intro hr
    refine ⟨hr, ?_⟩
    have hr' : H.derivative.derivative.eval₂ (RingHom.id 𝕜) r = 0 := by
      simpa using hr
    simpa using hsep.eval₂_derivative_ne_zero (RingHom.id 𝕜) hr'

omit [CharZero 𝕜] [IsAlgClosed 𝕜] in
theorem cuspParameters_nodup (H : 𝕜[X])
    (hsep : H.derivative.derivative.Separable) :
    (cuspParameters H).Nodup :=
  nodup_roots hsep

/-- A degree-`n` polynomial has `n - 2` second-derivative roots, counted
with multiplicity. -/
theorem card_cuspParameters (H : 𝕜[X]) (n : ℕ) (hdeg : H.natDegree = n)
    : (cuspParameters H).card = n - 2 := by
  calc
    (cuspParameters H).card = H.derivative.derivative.natDegree :=
      (IsAlgClosed.splits H.derivative.derivative).natDegree_eq_card_roots.symm
    _ = n - 2 := by
      rw [natDegree_derivative, natDegree_derivative, hdeg]
      omega

/-- Under the paper's squarefreeness hypothesis, the `n - 2` parameters are
distinct and each is an ordinary cusp parameter. -/
theorem cuspParameters_count_and_nodup_of_squarefree (H : 𝕜[X]) (n : ℕ)
    (hdeg : H.natDegree = n) (hsq : Squarefree H.derivative.derivative) :
    (cuspParameters H).card = n - 2 ∧ (cuspParameters H).Nodup := by
  refine ⟨card_cuspParameters H n hdeg, ?_⟩
  exact cuspParameters_nodup H
    (PerfectField.separable_iff_squarefree.mpr hsq)

end DiscriminantPencils
