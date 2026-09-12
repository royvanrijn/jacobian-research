/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GVC.BaseChange
import GVC.PhaseKernel

/-!
# Unconditional characteristic-zero counterexamples

The concrete phase bridge is constructed in `GVC.PhaseKernel`.  This file
feeds it through the already proved coefficient base-change and padding
theorems, removing the last hypothesis from the negative half of the
dimension classification.
-/

namespace GVC

/-- The ternary GVC counterexample exists over every characteristic-zero
field. -/
theorem verified_gvc3_charZero_not_generalizedVanishingConjecture
    {K : Type*} [Field K] [CharZero K] :
    ¬ GeneralizedVanishingConjecture (Fin 3) K :=
  gvc3_charZero_fin_not_generalizedVanishingConjecture
    verifiedConcreteCounterexampleBridge (K := K) (by omega)

/-- Unconditional negative half of the paper's dimension classification:
GVC fails over every characteristic-zero field in every finite dimension
at least three. -/
theorem verified_gvc3_charZero_fin_not_generalizedVanishingConjecture
    {K : Type*} [Field K] [CharZero K] {n : ℕ} (hn : 3 ≤ n) :
    ¬ GeneralizedVanishingConjecture (Fin n) K :=
  gvc3_charZero_fin_not_generalizedVanishingConjecture
    verifiedConcreteCounterexampleBridge (K := K) hn

end GVC
