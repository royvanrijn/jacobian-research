import QuadraticNorms.PeriodicCoreCertificate
import Mathlib.Data.ZMod.Basic

namespace QuadraticNorms

private def C11Residues : List ℕ := [1, 3, 4, 5, 9]
private def C11ScanBound : ℕ := 2 ^ 25
private def C11ScanMask : ℕ := periodicCoreMask 11 C11ScanBound 5793 C11Residues

private theorem C11ScanMask_bit {n : ℕ} (hn : 0 < n) (hle : n ≤ C11ScanBound) :
    C11ScanMask.testBit n = true ↔ C11 n := by
  rw [C11ScanMask, periodicCoreMask_bit hn hle (hle.trans_lt (by decide))]
  · simp [C11, C11Residues, List.mem_cons]
  · decide
  · decide
  · intro a b ha hb
    have rule : ∀ x y : ZMod 11, x.val ∈ C11Residues → y.val ∈ C11Residues →
        (x * y).val ∈ C11Residues := by decide
    simpa only [← Nat.cast_mul, ZMod.val_natCast] using rule a b
      (by simpa only [ZMod.val_natCast] using ha) (by simpa only [ZMod.val_natCast] using hb)
  · intro a b hab hb
    have rule : ∀ x y : ZMod 11, (x * y).val ∈ C11Residues → y.val ∈ C11Residues →
        x.val ∈ C11Residues := by decide
    simpa only [ZMod.val_natCast] using rule a b
      (by simpa only [← Nat.cast_mul, ZMod.val_natCast] using hab)
      (by simpa only [ZMod.val_natCast] using hb)

private theorem C11ScanMask_zero : C11ScanMask.testBit 0 = false :=
  periodicCoreMask_zero (by simp [C11Residues]) (by decide)

private def C11RequestedMask (bound : ℕ) : ℕ :=
  let clipped := 2 ^ (bound + 1) - 1
  let evenPos := progressionMask 2 (bound / 2) <<< 2
  let forbidden := progressionMask 22 (bound / 22 + 1)
  evenPos &&& (clipped ^^^ forbidden)

private theorem C11RequestedMask_bit {bound n : ℕ} (hn : n ≤ bound) :
    (C11RequestedMask bound).testBit n = true ↔ 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 11 ∣ n := by
  unfold C11RequestedMask
  rw [Nat.testBit_and, Nat.testBit_shiftLeft, Nat.testBit_xor, Nat.testBit_two_pow_sub_one,
    decide_eq_true (by omega : n < bound + 1), Bool.true_xor, Bool.and_eq_true,
    Bool.and_eq_true, decide_eq_true_eq, progressionMask_bit]
  have hneg (b : Bool) : (!b) = true ↔ ¬ b = true := by cases b <;> decide
  rw [hneg, progressionMask_bit]
  constructor
  · rintro ⟨⟨h₂, j, _, hj⟩, hnot⟩
    refine ⟨h₂, by omega, fun h => ?_⟩
    have hmod : n % 22 = 0 := by omega
    apply hnot
    exact ⟨n / 22, by have := Nat.div_le_div_right (c := 22) hn; omega, by omega⟩
  · rintro ⟨h₂, he, h₁₁⟩
    refine ⟨⟨h₂, (n - 2) / 2, by omega, by omega⟩, ?_⟩
    rintro ⟨k, _, rfl⟩
    exact h₁₁ ⟨2 * k, by omega⟩

private def C11ShortAnchors : List ℕ := (List.range 5001).filter C11ScanMask.testBit

private theorem C11FullScanPasses :
    C11RequestedMask C11ScanBound &&& packedShiftUnion C11ShortAnchors C11ScanMask = C11RequestedMask C11ScanBound := by
  native_decide

/-- Full exact C₁₁ binary coverage through the published finite cutoff. -/
theorem C11_finite_core_coverage {n : ℕ} (hle : n ≤ 2 ^ 25) (h₂ : 2 ≤ n)
    (he : n % 2 = 0) (h₁₁ : ¬ 11 ∣ n) : InSumset C11 C11 n := by
  have hl : n ≤ C11ScanBound := hle
  have hn := C11RequestedMask_bit hl |>.mpr ⟨h₂, he, h₁₁⟩
  have hbit := congrArg (fun bits : ℕ => bits.testBit n) C11FullScanPasses
  rw [Nat.testBit_and, hn, Bool.true_and] at hbit
  obtain ⟨a, ha, hal, hb⟩ := packedShiftUnion_bit.mp hbit
  obtain ⟨_, hab⟩ := List.mem_filter.mp ha
  have hap : 0 < a := by
    by_contra! h
    obtain rfl : a = 0 := by omega
    rw [C11ScanMask_zero] at hab
    exact (by decide : (false : Bool) ≠ true) hab
  have hbp : 0 < n - a := by
    by_contra! h
    rw [show n - a = 0 by omega, C11ScanMask_zero] at hb
    exact (by decide : (false : Bool) ≠ true) hb
  exact ⟨a, n - a, (C11ScanMask_bit hap (hal.trans hl)).mp hab,
    (C11ScanMask_bit hbp ((Nat.sub_le ..).trans hl)).mp hb, Nat.add_sub_of_le hal⟩

end QuadraticNorms
