/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib

/-!
# The numerical genus calculation

This file isolates the final arithmetic step of the paper: once the genus
formula supplies the total delta invariant and the cusp count is `n - 2`,
the remaining number of delta-one singularities is
`(n - 2)(n - 3) / 2`.
-/

namespace DiscriminantPencils

/-- Twice the node-count identity, avoiding division while exposing the
exact algebraic content of the genus calculation. -/
theorem twice_node_count_from_genus (n N : ℕ) (hn : 3 ≤ n)
    (hgenus : (n - 1) * (n - 2) = 2 * ((n - 2) + N)) :
    2 * N = (n - 2) * (n - 3) := by
  have h₁ : n - 1 = (n - 2) + 1 := by omega
  have h₂ : n - 2 = (n - 3) + 1 := by omega
  rw [h₁] at hgenus
  nlinarith

/-- The node count in the form stated in the paper. -/
theorem node_count_from_genus (n N : ℕ) (hn : 3 ≤ n)
    (hgenus : (n - 1) * (n - 2) = 2 * ((n - 2) + N)) :
    N = (n - 2) * (n - 3) / 2 := by
  have htwo := twice_node_count_from_genus n N hn hgenus
  omega

end DiscriminantPencils
