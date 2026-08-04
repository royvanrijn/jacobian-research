/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Topology.Algebra.Ring.Real

/-!
# The continuous envelope closure

This isolates the final topological step in Section 6 of the manuscript.
The finite-support Newton argument supplies continuity, a positive starting
gap, and a later nonpositive gap.  The intermediate value theorem then
produces a common threshold.
-/

namespace GVC

/-- A continuous envelope gap which starts positive and later becomes
nonpositive reaches a common threshold. -/
theorem exists_common_threshold
    (gap : ℝ → ℝ) (s₀ s₁ : ℝ)
    (hs : s₀ ≤ s₁)
    (hcont : ContinuousOn gap (Set.Icc s₀ s₁))
    (hstart : 0 < gap s₀) (hend : gap s₁ ≤ 0) :
    ∃ s ∈ Set.Icc s₀ s₁, gap s = 0 := by
  have hz : (0 : ℝ) ∈ Set.Icc (gap s₁) (gap s₀) := ⟨hend, hstart.le⟩
  obtain ⟨s, hsIcc, hs0⟩ := intermediate_value_Icc' hs hcont hz
  exact ⟨s, hsIcc, hs0⟩

end GVC
