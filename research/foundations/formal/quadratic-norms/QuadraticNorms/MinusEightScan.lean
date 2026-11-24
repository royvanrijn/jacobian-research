import QuadraticNorms.PeriodicCoreCertificate
import Mathlib.Data.ZMod.Basic

namespace QuadraticNorms.MinusEightScan

private def eightResidues : List ℕ := [1, 3]
def eightBound : ℕ := 2 ^ 27
@[irreducible] def eightMask : ℕ := periodicCoreMask 8 eightBound 11586 eightResidues

private theorem coreMask_bit_aux {n bound : ℕ} (hn : 0 < n) (hle : n ≤ bound)
    (hbound : bound < 11586 ^ 2) :
    (periodicCoreMask 8 bound 11586 eightResidues).testBit n = true ↔ C8 n := by
  refine (periodicCoreMask_bit hn hle (hle.trans_lt hbound) ?_ ?_ ?_ ?_).trans ?_
  · decide
  · decide
  · intro a b ha hb
    have rule : ∀ x y : ZMod 8, x.val ∈ eightResidues → y.val ∈ eightResidues →
        (x * y).val ∈ eightResidues := by decide
    simpa only [← Nat.cast_mul, ZMod.val_natCast] using rule a b
      (by simpa only [ZMod.val_natCast] using ha) (by simpa only [ZMod.val_natCast] using hb)
  · intro a b hab hb
    have rule : ∀ x y : ZMod 8, (x * y).val ∈ eightResidues → y.val ∈ eightResidues →
        x.val ∈ eightResidues := by decide
    simpa only [ZMod.val_natCast] using rule a b
      (by simpa only [← Nat.cast_mul, ZMod.val_natCast] using hab)
      (by simpa only [ZMod.val_natCast] using hb)

  · simp [C8, PrimeCore, eightResidues]

private theorem eightMask_eq : eightMask = periodicCoreMask 8 eightBound 11586 eightResidues := by
  unfold eightMask
  rfl

theorem eightMask_bit {n : ℕ} (hn : 0 < n) (hle : n ≤ eightBound) :
    eightMask.testBit n = true ↔ C8 n :=
  Eq.mpr (congrArg (fun values => values.testBit n = true ↔ C8 n) eightMask_eq)
    (coreMask_bit_aux hn hle (by decide))

theorem eightMask_zero : eightMask.testBit 0 = false :=
  Eq.mpr (congrArg (fun values => values.testBit 0 = false) eightMask_eq)
    (periodicCoreMask_zero (by simp [eightResidues]) (by decide))

def eightOrdinaryMask (bound : ℕ) : ℕ :=
  packedShiftUnion [2, 4, 6] (progressionMask 8 (bound / 8 + 1))

theorem eightOrdinaryMask_has_bit {bound n : ℕ}
    (hle : n ≤ bound) (h₂ : n % 2 = 0) (h₈ : ¬ 8 ∣ n) :
    (eightOrdinaryMask bound).testBit n = true := by
  have hres : n % 8 = 2 ∨ n % 8 = 4 ∨ n % 8 = 6 := by omega
  have hmem : n % 8 ∈ [2, 4, 6] := by simpa using hres
  have hshift : n - n % 8 = 8 * (n / 8) := by omega
  apply packedShiftUnion_bit.mpr
  exact ⟨n % 8, hmem, Nat.mod_le _ _, (progressionMask_bit ..).mpr
    ⟨n / 8, Nat.lt_succ_of_le (Nat.div_le_div_right hle), hshift⟩⟩

def eightMixedMask (bound : ℕ) : ℕ := progressionMask 2 (bound / 2) <<< 3

theorem eightMixedMask_has_bit {bound n : ℕ}
    (hle : n ≤ bound) (h₃ : 3 ≤ n) (hodd : n % 2 = 1) :
    (eightMixedMask bound).testBit n = true := by
  rw [eightMixedMask, Nat.testBit_shiftLeft, Bool.and_eq_true, decide_eq_true_eq, progressionMask_bit]
  exact ⟨h₃, (n - 3) / 2, by omega, by omega⟩

@[irreducible] def scanAnchors (values : ℕ) : List ℕ := (List.range 5001).filter values.testBit

theorem scanAnchors_sound (values a : ℕ) (ha : a ∈ scanAnchors values) :
    values.testBit a = true := by
  unfold scanAnchors at ha
  exact (List.mem_filter.mp ha).2

@[irreducible] def ordinaryTarget (bound : ℕ) := eightOrdinaryMask bound &&& (2 ^ (bound + 1) - 1)
@[irreducible] def mixedTarget (bound : ℕ) := eightMixedMask bound &&& (2 ^ (bound + 1) - 1)

private theorem target_clip_has_bit {bound target n : ℕ} (hn : n ≤ bound)
    (ht : target.testBit n = true) : (target &&& (2 ^ (bound + 1) - 1)).testBit n = true := by
  rw [Nat.testBit_and, ht, Nat.testBit_two_pow_sub_one, decide_eq_true (Nat.lt_succ_of_le hn)]
  decide

theorem ordinaryTarget_has_bit {bound n : ℕ}
    (hle : n ≤ bound) (h₂ : n % 2 = 0) (h₈ : ¬ 8 ∣ n) :
    (ordinaryTarget bound).testBit n = true := by
  unfold ordinaryTarget
  exact target_clip_has_bit hle (eightOrdinaryMask_has_bit hle h₂ h₈)

theorem mixedTarget_has_bit {bound n : ℕ}
    (hle : n ≤ bound) (h₃ : 3 ≤ n) (hodd : n % 2 = 1) :
    (mixedTarget bound).testBit n = true := by
  unfold mixedTarget
  exact target_clip_has_bit hle (eightMixedMask_has_bit hle h₃ hodd)

@[irreducible] def dualChecks (ordinary mixed values : ℕ) (anchors : List ℕ) : Bool :=
  ((ordinary &&& packedShiftUnion anchors values) == ordinary) &&
    ((mixed &&& packedShiftUnion (anchors.map (2 * ·)) values) == mixed)

theorem dualChecks_eq_true {ordinary mixed values : ℕ} {anchors : List ℕ}
    (h : dualChecks ordinary mixed values anchors = true) :
    ordinary &&& packedShiftUnion anchors values = ordinary ∧
      mixed &&& packedShiftUnion (anchors.map (2 * ·)) values = mixed := by
  simpa only [dualChecks, Bool.and_eq_true, beq_iff_eq] using h

theorem eightScanPasses_checked :
    dualChecks (ordinaryTarget eightBound) (mixedTarget eightBound)
      eightMask (scanAnchors eightMask) = true := by native_decide

end QuadraticNorms.MinusEightScan
