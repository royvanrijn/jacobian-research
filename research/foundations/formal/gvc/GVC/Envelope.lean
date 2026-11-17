/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.Topology.Order.IntermediateValue
import Mathlib.Topology.Algebra.Ring.Real
import Mathlib.Tactic.Linarith

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

/-- Once the finite envelopes are on their final affine pieces, strict
right-to-left separation makes the gap's final slope negative.  Positivity
cannot then persist forever. -/
theorem exists_common_threshold_of_negative_affine_tail
    (gap : ℝ → ℝ) (s₀ constant lowerSlope upperSlope : ℝ)
    (hcont : Continuous gap)
    (hstart : 0 < gap s₀)
    (hslope : upperSlope < lowerSlope)
    (htail : ∀ s, s₀ ≤ s →
      gap s = constant + s * (upperSlope - lowerSlope)) :
    ∃ s, s₀ ≤ s ∧ gap s = 0 := by
  let slopeGap := lowerSlope - upperSlope
  have hslopeGap : 0 < slopeGap := by
    dsimp [slopeGap]
    exact sub_pos.mpr hslope
  let s₁ := max s₀ ((constant + 1) / slopeGap)
  have hs₀s₁ : s₀ ≤ s₁ := le_max_left _ _
  have hratio : (constant + 1) / slopeGap ≤ s₁ := le_max_right _ _
  have hscaled : constant + 1 ≤ s₁ * slopeGap :=
    (div_le_iff₀ hslopeGap).mp hratio
  have hend : gap s₁ ≤ 0 := by
    rw [htail s₁ hs₀s₁]
    dsimp [slopeGap] at hscaled
    nlinarith
  obtain ⟨s, hs, hzero⟩ :=
    exists_common_threshold gap s₀ s₁ hs₀s₁
      hcont.continuousOn hstart hend
  exact ⟨s, hs.1, hzero⟩

end GVC
