import QuadraticNorms.MinusEightScan

namespace QuadraticNorms

open MinusEightScan

private theorem decodeBit {values bound a : ℕ} (hzero : values.testBit 0 = false)
    (hsound : ∀ j, 0 < j → j ≤ bound → values.testBit j = true → C8 j)
    (hle : a ≤ bound) (h : values.testBit a = true) : C8 a := by
  by_cases hpos : 0 < a
  · exact hsound a hpos hle h
  · obtain rfl : a = 0 := by omega
    simp [hzero] at h

private theorem ordinary_sound {values bound n : ℕ}
    (hsound : ∀ j ≤ bound, values.testBit j = true → C8 j)
    {anchors : List ℕ} (hanchor : ∀ j ∈ anchors, values.testBit j = true)
    {targets : ℕ} (hcheck : targets &&& packedShiftUnion anchors values = targets)
    (hle : n ≤ bound) (htarget : targets.testBit n = true) : InSumset C8 C8 n := by
  obtain ⟨a, b, ha, hb, hab⟩ := packedShiftUnion_covers
    (P := fun a => a ≤ bound → C8 a) (Q := fun b => b ≤ bound → C8 b)
    (fun a ha hle' => hsound a hle' (hanchor a ha)) (fun b hb hle' => hsound b hle' hb) hcheck htarget
  exact ⟨a, b, ha (by omega), hb (by omega), hab⟩

private theorem mixed_sound {values bound n : ℕ}
    (hsound : ∀ j ≤ bound, values.testBit j = true → C8 j)
    {anchors : List ℕ} (hanchor : ∀ j ∈ anchors, values.testBit j = true)
    {targets : ℕ} (hcheck : targets &&& packedShiftUnion (anchors.map (2 * ·)) values = targets)
    (hle : n ≤ bound) (htarget : targets.testBit n = true) : InSumset C8 (Dilate 2 C8) n := by
  have hsum := packedShiftUnion_covers
    (P := fun a => a ≤ bound → Dilate 2 C8 a) (Q := fun b => b ≤ bound → C8 b)
    (fun a ha hale => by
      obtain ⟨x, hx, rfl⟩ := List.mem_map.mp ha
      exact ⟨x, hsound x (by omega) (hanchor x hx), rfl⟩)
    (fun b hb hle' => hsound b hle' hb) hcheck htarget
  obtain ⟨a, b, ha, hb, hab⟩ := hsum
  exact ⟨b, a, hb (by omega), ha (by omega), by omega⟩

private theorem eightValues_sound (a : ℕ) :
    a ≤ eightBound → eightMask.testBit a = true → C8 a :=
  decodeBit eightMask_zero (fun _ ha hl => (eightMask_bit ha hl).mp)

theorem C8_finite_ordinary_coverage {n : ℕ} (hle : n ≤ 2 ^ 27)
    (he : n % 2 = 0) (h₈ : ¬ 8 ∣ n) : InSumset C8 C8 n :=
  ordinary_sound (bound := eightBound) eightValues_sound (scanAnchors_sound eightMask)
    (dualChecks_eq_true eightScanPasses_checked).1 hle
    (ordinaryTarget_has_bit (bound := eightBound) hle he h₈)

theorem C8_finite_mixed_coverage {n : ℕ} (hle : n ≤ 2 ^ 27)
    (h₃ : 3 ≤ n) (ho : n % 2 = 1) : InSumset C8 (Dilate 2 C8) n :=
  mixed_sound (bound := eightBound) eightValues_sound (scanAnchors_sound eightMask)
    (dualChecks_eq_true eightScanPasses_checked).2 hle
    (mixedTarget_has_bit (bound := eightBound) hle h₃ ho)

end QuadraticNorms
