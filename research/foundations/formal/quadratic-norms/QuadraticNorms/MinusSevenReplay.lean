import QuadraticNorms.PeriodicCoreCertificate
import Mathlib.Data.ZMod.Basic

namespace QuadraticNorms

private def sevenResidues : List ℕ := [1, 2, 4]
private def sevenBound : ℕ := 2 ^ 24
private def sevenMask : ℕ := periodicCoreMask 7 sevenBound 4097 sevenResidues

private theorem sevenMask_bit {n : ℕ} (hn : 0 < n) (hle : n ≤ sevenBound) :
    sevenMask.testBit n = true ↔ C7 n := by
  rw [sevenMask, periodicCoreMask_bit hn hle (hle.trans_lt (by decide))]
  · simp [C7, PrimeCore, sevenResidues]
  · decide
  · decide
  · intro a b ha hb
    have rule : ∀ x y : ZMod 7, x.val ∈ sevenResidues → y.val ∈ sevenResidues →
        (x * y).val ∈ sevenResidues := by decide
    simpa only [← Nat.cast_mul, ZMod.val_natCast] using rule a b
      (by simpa only [ZMod.val_natCast] using ha) (by simpa only [ZMod.val_natCast] using hb)
  · intro a b hab hb
    have rule : ∀ x y : ZMod 7, (x * y).val ∈ sevenResidues → y.val ∈ sevenResidues →
        x.val ∈ sevenResidues := by decide
    simpa only [ZMod.val_natCast] using rule a b
      (by simpa only [← Nat.cast_mul, ZMod.val_natCast] using hab)
      (by simpa only [ZMod.val_natCast] using hb)

private theorem sevenMask_zero : sevenMask.testBit 0 = false :=
  periodicCoreMask_zero (by simp [sevenResidues]) (by decide)

private def sevenRequestedMask (bound : ℕ) : ℕ :=
  let all := 2 ^ (bound + 1) - 1
  let forbidden := progressionMask 7 (bound / 7 + 1) ||| packedValues [1, 11, 29]
  all &&& (all ^^^ forbidden)

private theorem sevenRequestedMask_has_bit {bound n : ℕ}
    (hle : n ≤ bound) (h₂ : 2 ≤ n) (h₇ : ¬ 7 ∣ n) (h₁₁ : n ≠ 11) (h₂₉ : n ≠ 29) :
    (sevenRequestedMask bound).testBit n = true := by
  have hleft : (progressionMask 7 (bound / 7 + 1)).testBit n = false := by
    apply Bool.eq_false_iff.mpr
    intro hbit
    obtain ⟨j, _, rfl⟩ := (progressionMask_bit ..).mp hbit
    exact h₇ (dvd_mul_right ..)
  have hright : (packedValues [1, 11, 29]).testBit n = false := by
    apply Bool.eq_false_iff.mpr
    intro hbit
    rcases (show n = 1 ∨ n = 11 ∨ n = 29 by simpa using packedValues_bit.mp hbit) with h | h | h <;> omega
  simp only [sevenRequestedMask, Nat.testBit_and, Nat.testBit_two_pow_sub_one,
    decide_eq_true (by omega : n < bound + 1), Nat.testBit_xor, Nat.testBit_or, hleft, hright]
  decide

private def sevenAnchors : List ℕ := (List.range 5001).filter sevenMask.testBit

private theorem sevenScanPasses_checked :
    sevenRequestedMask sevenBound &&& packedShiftUnion sevenAnchors sevenMask =
      sevenRequestedMask sevenBound := by native_decide

private theorem seven_bit_below {a : ℕ} (hle : a ≤ sevenBound) (h : sevenMask.testBit a = true) :
    C7 a := by
  have hpos : 0 < a := by
    by_contra! hno
    obtain rfl : a = 0 := by omega
    rw [sevenMask_zero] at h
    exact Bool.noConfusion h
  exact (sevenMask_bit hpos hle).mp h

theorem C7_finite_core_coverage {n : ℕ} (hle : n ≤ 2 ^ 24)
    (h₂ : 2 ≤ n) (h₇ : ¬ 7 ∣ n) (h₁₁ : n ≠ 11) (h₂₉ : n ≠ 29) : InSumset C7 C7 n := by
  have hl : n ≤ sevenBound := hle
  have htarget := sevenRequestedMask_has_bit hl h₂ h₇ h₁₁ h₂₉
  have hbit := congrArg (fun bits : ℕ => bits.testBit n) sevenScanPasses_checked
  rw [Nat.testBit_and, htarget, Bool.true_and] at hbit
  obtain ⟨a, ha, hal, hb⟩ := packedShiftUnion_bit.mp hbit
  have hab := (List.mem_filter.mp ha).2
  exact ⟨a, n - a, seven_bit_below (hal.trans hl) hab,
    seven_bit_below ((Nat.sub_le ..).trans hl) hb, Nat.add_sub_of_le hal⟩

end QuadraticNorms
