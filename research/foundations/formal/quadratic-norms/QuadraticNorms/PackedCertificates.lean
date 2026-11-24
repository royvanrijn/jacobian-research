import QuadraticNorms.FiniteValues
import Mathlib.Data.Nat.Factors

namespace QuadraticNorms

/-- Executable membership in a prime core. The empty factorization
represents one, and zero is excluded explicitly. -/
def primeCoreTest (allowed : ℕ → Bool) (n : ℕ) : Bool :=
  (n != 0) && n.primeFactorsList.all allowed

theorem primeCoreTest_spec {allowed : ℕ → Bool} {n : ℕ} :
    primeCoreTest allowed n = true ↔ PrimeCore (fun p => allowed p = true) n := by
  simp only [primeCoreTest, Bool.and_eq_true, bne_iff_ne, List.all_eq_true,
    Nat.mem_primeFactorsList', PrimeCore, Nat.pos_iff_ne_zero]
  exact ⟨fun ⟨hn, h⟩ => ⟨hn, fun p hp hd => h p ⟨hp, hd, hn⟩⟩,
    fun ⟨hn, h⟩ => ⟨hn, fun p hp => h p hp.1 hp.2.1⟩⟩

/-- The bit at position n records membership of n among the values. -/
def packedValues (values : List ℕ) : ℕ :=
  values.foldr (fun a acc => 2 ^ a ||| acc) 0

theorem packedValues_bit {values : List ℕ} {n : ℕ} :
    (packedValues values).testBit n = true ↔ n ∈ values := by
  induction values with
  | nil => simp [packedValues]
  | cons a values ih =>
    change (2 ^ a ||| packedValues values).testBit n = true ↔ _
    rw [Nat.testBit_or, Nat.testBit_two_pow, Bool.or_eq_true, decide_eq_true_eq, ih]
    simp [eq_comm]

/-- The packed progression with bits at 0, step, …, (count−1)*step.
The doubling step constructs it in logarithmically many shifts. -/
def progressionMask (step count : ℕ) : ℕ :=
  if count = 0 then 0 else
    let half := count / 2
    let initial := progressionMask step half
    let doubled := initial ||| (initial <<< (step * half))
    if count % 2 = 1 then doubled ||| (1 <<< (step * (2 * half))) else doubled
termination_by count
decreasing_by exact Nat.div_lt_self (by omega) (by decide)

theorem progressionMask_bit (step count n : ℕ) :
    (progressionMask step count).testBit n = true ↔ ∃ i < count, n = step * i := by
  induction count using Nat.strong_induction_on generalizing n with
  | h count ih =>
    by_cases hc : count = 0
    · subst count
      rw [progressionMask]
      simp
    have hhalf : count / 2 < count := Nat.div_lt_self (Nat.pos_of_ne_zero hc) (by decide)
    let half := count / 2
    have htwice : 2 * half + count % 2 = count := by dsimp [half]; omega
    have hdouble (n : ℕ) :
        ((progressionMask step half) ||| ((progressionMask step half) <<< (step * half))).testBit n
          = true ↔ ∃ i < 2 * half, n = step * i := by
      rw [Nat.testBit_or, Nat.testBit_shiftLeft, Bool.or_eq_true, Bool.and_eq_true, decide_eq_true_eq,
        ih half hhalf, ih half hhalf]
      constructor
      · rintro (⟨i, hi, rfl⟩ | ⟨hlo, i, hi, heq⟩)
        · exact ⟨i, by omega, rfl⟩
        · refine ⟨half + i, by omega, ?_⟩
          rw [Nat.mul_add]
          omega
      · rintro ⟨i, hi, rfl⟩
        by_cases hlow : i < half
        · exact Or.inl ⟨i, hlow, rfl⟩
        · have hle : half ≤ i := by omega
          exact Or.inr ⟨Nat.mul_le_mul_left _ hle, i - half, by omega,
            (Nat.mul_sub_left_distrib ..).symm⟩
    rw [progressionMask, if_neg hc]
    dsimp only
    change (if count % 2 = 1 then
      (progressionMask step half ||| (progressionMask step half <<< (step * half))) |||
        (1 <<< (step * (2 * half)))
      else progressionMask step half ||| (progressionMask step half <<< (step * half))).testBit n = true ↔ _
    by_cases hone : count % 2 = 1
    · rw [if_pos hone, Nat.testBit_or, Bool.or_eq_true, hdouble, Nat.testBit_shiftLeft,
        Bool.and_eq_true, decide_eq_true_eq]
      have hsingle : step * (2 * half) ≤ n ∧ (1 : ℕ).testBit (n - step * (2 * half)) = true ↔
          n = step * (2 * half) := by
        rw [show (1 : ℕ) = 2 ^ 0 by decide, Nat.testBit_two_pow, decide_eq_true_eq]
        omega
      rw [hsingle]
      constructor
      · rintro (⟨i, hi, heq⟩ | heq)
        · exact ⟨i, by omega, heq⟩
        · exact ⟨2 * half, by omega, heq⟩
      · rintro ⟨i, hi, heq⟩
        by_cases hlt : i < 2 * half
        · exact Or.inl ⟨i, hlt, heq⟩
        · obtain rfl : i = 2 * half := by omega
          exact Or.inr heq
    · rw [if_neg hone, hdouble]
      obtain heq : 2 * half = count := by omega
      rw [heq]

/-- The union of shifts of a membership mask. A shift by an anchor
corresponds to adding that anchor to every represented value. -/
def packedShiftUnion (anchors : List ℕ) (values : ℕ) : ℕ :=
  anchors.foldr (fun a acc => (values <<< a) ||| acc) 0

theorem packedShiftUnion_bit {anchors : List ℕ} {values n : ℕ} :
    (packedShiftUnion anchors values).testBit n = true ↔
      ∃ a ∈ anchors, a ≤ n ∧ values.testBit (n - a) = true := by
  induction anchors with
  | nil => simp [packedShiftUnion]
  | cons a anchors ih =>
    change ((values <<< a) ||| packedShiftUnion anchors values).testBit n = true ↔ _
    rw [Nat.testBit_or, Nat.testBit_shiftLeft, Bool.or_eq_true, Bool.and_eq_true, decide_eq_true_eq, ih]
    simp only [List.mem_cons, exists_eq_or_imp]

/-- A packed inclusion certificate gives a genuine decomposition at
every marked target when the anchor and second-summand masks are sound. -/
theorem packedShiftUnion_covers {P Q : ℕ → Prop} {anchors : List ℕ}
    {values targets n : ℕ}
    (hanchors : ∀ a ∈ anchors, P a)
    (hvalues : ∀ b, values.testBit b = true → Q b)
    (hcheck : targets &&& packedShiftUnion anchors values = targets)
    (hn : targets.testBit n = true) : InSumset P Q n := by
  have hbit := congrArg (fun bits : ℕ => bits.testBit n) hcheck
  rw [Nat.testBit_and, hn, Bool.true_and] at hbit
  obtain ⟨a, ha, hal, hb⟩ := packedShiftUnion_bit.mp hbit
  exact ⟨a, n - a, hanchors a ha, hvalues _ hb, Nat.add_sub_of_le hal⟩

/-- In the coordinate-based scanner, each marked position is an
authenticated primitive value, and every such value in the box occurs. -/
theorem packedPrimitiveValues_bit {b c : ℤ} {n bound : ℕ}
    (hb : b = -1 ∨ b = 0 ∨ b = 1) (hc : 1 ≤ c) :
    (packedValues (finitePrimitiveValues b c bound)).testBit n = true ↔
      PrimitiveValues b c n ∧ n ≤ bound :=
  packedValues_bit.trans (mem_finitePrimitiveValues hb hc)

private def boundedPrimitiveSumMask (b c : ℤ) (bound : ℕ) : ℕ :=
  let vals := finitePrimitiveValues b c bound
  packedShiftUnion vals (packedValues vals)

private def targets7 (bound : ℕ) : List ℕ :=
  (List.range (bound + 1)).filter fun n =>
    2 ≤ n ∧ ¬ 49 ∣ n ∧ n ≠ 7 ∧ n ≠ 77 ∧ n ≠ 203

set_option maxRecDepth 4096 in
/-- The full −7 classification through the largest exceptional
target, using exact coordinate masks and a kernel-checked inclusion. -/
theorem P7_sum_through_203 {n : ℕ} (hn : n ≤ 203) :
    InSumset P7 P7 n ↔ 2 ≤ n ∧ ¬ 49 ∣ n ∧ n ≠ 7 ∧ n ≠ 77 ∧ n ≠ 203 := by
  refine ⟨P7_sum_necessary, fun h => ?_⟩
  have hinclusion : packedValues (targets7 203) &&& boundedPrimitiveSumMask 1 2 203 =
      packedValues (targets7 203) := by decide +kernel
  refine packedShiftUnion_covers (P := P7) (Q := P7)
    (fun a ha => ((mem_finitePrimitiveValues (by decide) (by decide)).mp ha).1)
    (fun _ hb => ((packedPrimitiveValues_bit (by decide) (by decide)).mp hb).1)
    hinclusion (packedValues_bit.mpr ?_)
  exact List.mem_filter.mpr ⟨List.mem_range.mpr (by omega), by simpa using h⟩

end QuadraticNorms
