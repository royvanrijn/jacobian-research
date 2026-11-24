import QuadraticNorms.MinusEightEulerBounds
import QuadraticNorms.MinusEightSieve
import QuadraticNorms.RoundedPrimeTails
import QuadraticNorms.SplitDivisorFrontier
import Mathlib.Data.Nat.Sqrt
import Mathlib.Data.Nat.Log

open scoped BigOperators

namespace QuadraticNorms

/-- The input for one selected-shell inequality on a dyadic bin in the odd part. -/
structure MinusEightIntervalBin where
  dyadicExponent : ℕ
  exponent : ℕ
  divisorLimit : ℕ
  cubeRoot : ℕ
  primeCutoff : ℕ
  deriving DecidableEq

/-- The divisor, cube root, and cutoff inputs from the manuscript's
three groups of interval rows. The prime counts and tails are recomputed. -/
def minusEightIntervalBins : List MinusEightIntervalBin := [
  ⟨0, 27, 128, 512, 220⟩,
  ⟨0, 28, 160, 645, 276⟩,
  ⟨0, 29, 160, 812, 310⟩,
  ⟨0, 30, 192, 1024, 381⟩,
  ⟨0, 31, 192, 1290, 428⟩,
  ⟨0, 32, 240, 1625, 537⟩,
  ⟨0, 33, 288, 2048, 660⟩,
  ⟨0, 34, 320, 2580, 782⟩,
  ⟨0, 35, 336, 3250, 899⟩,
  ⟨0, 36, 384, 4096, 1079⟩,
  ⟨0, 37, 448, 5160, 1308⟩,
  ⟨0, 38, 480, 6501, 1520⟩,
  ⟨0, 39, 576, 8192, 1869⟩,
  ⟨0, 40, 640, 10321, 2211⟩,
  ⟨0, 41, 672, 13003, 2544⟩,
  ⟨0, 42, 768, 16384, 3052⟩,
  ⟨0, 43, 864, 20642, 3634⟩,
  ⟨0, 44, 960, 26007, 4300⟩,
  ⟨0, 45, 1152, 32768, 5287⟩,
  ⟨0, 46, 1152, 41285, 5935⟩,
  ⟨0, 47, 1344, 52015, 7196⟩,
  ⟨0, 48, 1536, 65536, 8635⟩,
  ⟨0, 49, 1728, 82570, 10280⟩,
  ⟨0, 50, 1920, 104031, 12163⟩,
  ⟨0, 51, 2048, 131072, 14101⟩,
  ⟨0, 52, 2304, 165140, 16788⟩,
  ⟨0, 53, 2688, 208063, 20353⟩,
  ⟨0, 54, 2880, 262144, 23648⟩,
  ⟨0, 55, 3072, 330280, 27414⟩,
  ⟨0, 56, 3840, 416127, 34404⟩,
  ⟨0, 57, 4032, 524288, 39571⟩,
  ⟨0, 58, 4608, 660561, 47483⟩,
  ⟨0, 59, 5120, 832255, 56181⟩,
  ⟨1, 26, 128, 406, 169⟩,
  ⟨1, 27, 160, 512, 213⟩,
  ⟨1, 28, 160, 645, 239⟩,
  ⟨1, 29, 192, 812, 294⟩,
  ⟨1, 30, 192, 1024, 330⟩,
  ⟨1, 31, 240, 1290, 414⟩,
  ⟨1, 32, 288, 1625, 509⟩,
  ⟨1, 33, 320, 2048, 603⟩,
  ⟨1, 34, 336, 2580, 693⟩,
  ⟨1, 35, 384, 3250, 832⟩,
  ⟨1, 36, 448, 4096, 1009⟩,
  ⟨1, 37, 480, 5160, 1173⟩,
  ⟨1, 38, 576, 6501, 1442⟩,
  ⟨1, 39, 640, 8192, 1706⟩,
  ⟨1, 40, 672, 10321, 1962⟩,
  ⟨1, 41, 768, 13003, 2355⟩,
  ⟨1, 42, 864, 16384, 2804⟩,
  ⟨1, 43, 960, 20642, 3317⟩,
  ⟨1, 44, 1152, 26007, 4079⟩,
  ⟨1, 45, 1152, 32768, 4579⟩,
  ⟨1, 46, 1344, 41285, 5552⟩,
  ⟨1, 47, 1536, 52015, 6662⟩,
  ⟨1, 48, 1728, 65536, 7931⟩,
  ⟨1, 49, 1920, 82570, 9384⟩,
  ⟨1, 50, 2048, 104031, 10879⟩,
  ⟨1, 51, 2304, 131072, 12952⟩,
  ⟨1, 52, 2688, 165140, 15703⟩,
  ⟨1, 53, 2880, 208063, 18245⟩,
  ⟨1, 54, 3072, 262144, 21151⟩,
  ⟨1, 55, 3840, 330280, 26544⟩,
  ⟨1, 56, 4032, 416127, 30530⟩,
  ⟨1, 57, 4608, 524288, 36635⟩,
  ⟨1, 58, 5120, 660561, 43346⟩,
  ⟨1, 59, 5760, 832255, 51606⟩,
  ⟨2, 20, 64, 101, 48⟩,
  ⟨2, 21, 80, 128, 61⟩,
  ⟨2, 22, 84, 161, 70⟩,
  ⟨2, 23, 96, 203, 84⟩,
  ⟨2, 24, 120, 256, 106⟩,
  ⟨2, 25, 128, 322, 123⟩,
  ⟨2, 26, 160, 406, 155⟩,
  ⟨2, 27, 160, 512, 174⟩,
  ⟨2, 28, 192, 645, 214⟩,
  ⟨2, 29, 192, 812, 240⟩,
  ⟨2, 30, 240, 1024, 301⟩,
  ⟨2, 31, 288, 1290, 370⟩,
  ⟨2, 32, 320, 1625, 438⟩,
  ⟨2, 33, 336, 2048, 504⟩,
  ⟨2, 34, 384, 2580, 605⟩,
  ⟨2, 35, 448, 3250, 734⟩,
  ⟨2, 36, 480, 4096, 853⟩,
  ⟨2, 37, 576, 5160, 1049⟩,
  ⟨2, 38, 640, 6501, 1241⟩,
  ⟨2, 39, 672, 8192, 1427⟩,
  ⟨2, 40, 768, 10321, 1713⟩,
  ⟨2, 41, 864, 13003, 2039⟩,
  ⟨2, 42, 960, 16384, 2413⟩,
  ⟨2, 43, 1152, 20642, 2967⟩,
  ⟨2, 44, 1152, 26007, 3331⟩,
  ⟨2, 45, 1344, 32768, 4038⟩,
  ⟨2, 46, 1536, 41285, 4846⟩,
  ⟨2, 47, 1728, 52015, 5769⟩,
  ⟨2, 48, 1920, 65536, 6826⟩,
  ⟨2, 49, 2048, 82570, 7913⟩,
  ⟨2, 50, 2304, 104031, 9421⟩,
  ⟨2, 51, 2688, 131072, 11423⟩,
  ⟨2, 52, 2880, 165140, 13272⟩,
  ⟨2, 53, 3072, 208063, 15386⟩,
  ⟨2, 54, 3840, 262144, 19308⟩,
  ⟨2, 55, 4032, 330280, 22208⟩,
  ⟨2, 56, 4608, 416127, 26649⟩,
  ⟨2, 57, 5120, 524288, 31531⟩,
  ⟨2, 58, 5760, 660561, 37539⟩,
  ⟨2, 59, 6144, 832255, 43518⟩]

def minusEightSieveBound : ℕ := 65536

/-- Every inert prime up to the sieve limit, in increasing order. -/
def minusEightIntervalPrimes : List ℕ :=
  (List.range (minusEightSieveBound + 1)).filter (fun p => decide (MinusEightInertPrime p))

theorem mem_minusEight_intervalPrimes {p : ℕ} :
    p ∈ minusEightIntervalPrimes ↔ p ≤ minusEightSieveBound ∧ MinusEightInertPrime p := by
  simp [minusEightIntervalPrimes]

theorem minusEight_intervalPrimes_nodup : minusEightIntervalPrimes.Nodup := by
  unfold minusEightIntervalPrimes
  exact List.nodup_range.filter _

attribute [irreducible] minusEightIntervalPrimes

private def minusEightSplitPrefix : List ℕ := [3, 11, 17, 19, 41, 43, 59, 67, 73, 83, 89, 97, 107]

private theorem minusEight_splitPrefix_covers {p : ℕ} (h : p ≤ 107) (hp : MinusEightSplitPrime p) :
    p ∈ minusEightSplitPrefix := by
  have hs : ∀ q : Fin 108, MinusEightSplitPrime q.val → q.val ∈ minusEightSplitPrefix := by
    decide +kernel
  exact hs ⟨p, by omega⟩ hp

private theorem minusEight_split_divisors_le_envelope {n X : ℕ} (hn : 0 < n)
    (hnX : n ≤ X) (hX : X ≤ 2 ^ 62) :
    minusEightSplitDivisors n ≤ splitDivisorEnvelope 39 minusEightSplitPrefix X := by
  apply split_divisors_le_prefix_envelope MinusEightSplitPrime minusEightSplitPrefix (base := 3)
    (last := 107) hn hnX (by decide +kernel) (by decide +kernel)
    (fun _ hp hsplit => minusEight_splitPrefix_covers hp hsplit) (by decide)
  · intro p hp hsplit
    have hp₀ := hp.two_le
    have h := hsplit.2
    omega
  · exact hX.trans_lt (by decide +kernel)
  · exact hX.trans_lt (by decide +kernel)

def MinusEightIntervalBin.targetUpper (row : MinusEightIntervalBin) : ℕ :=
  2 ^ row.dyadicExponent * (2 ^ (row.exponent + 1) - 1)

def MinusEightIntervalBin.coefficientNumerator (row : MinusEightIntervalBin) : ℕ :=
  if row.dyadicExponent = 0 then 8 else if row.dyadicExponent = 1 then 2 else 4

def MinusEightIntervalBin.coefficientDenominator (row : MinusEightIntervalBin) : ℕ :=
  if row.dyadicExponent = 1 then 1 else 3

def MinusEightIntervalBin.coefficient (row : MinusEightIntervalBin) : ℚ :=
  (row.coefficientNumerator : ℚ) / row.coefficientDenominator

def MinusEightIntervalBin.errorBound (row : MinusEightIntervalBin) : ℚ :=
  18 / (5 * (row.cubeRoot : ℚ))

def MinusEightIntervalBin.allowance (row : MinusEightIntervalBin) : ℚ :=
  if row.dyadicExponent = 0 then 1 / 2 else 0

def MinusEightIntervalBin.lower (row : MinusEightIntervalBin) : ℚ :=
  if row.dyadicExponent = 2 then 2 / 3 else 1 / 2

def MinusEightIntervalBin.inertCount (row : MinusEightIntervalBin) (primes : List ℕ) : ℕ :=
  (primes.filter fun p => decide (p ≤ row.primeCutoff)).length

def MinusEightIntervalBin.tailPrimes (row : MinusEightIntervalBin) (primes : List ℕ) : List ℕ :=
  primes.filter fun p => decide (row.primeCutoff < p)

def MinusEightIntervalBin.inertTailBound (row : MinusEightIntervalBin) (primes : List ℕ) : ℚ :=
  let rounded := ((row.tailPrimes primes).map roundedInverseSquareWeight).sum
  (rounded : ℚ) / inverseSquareFixedPointScale + 1 / (4 * (minusEightSieveBound : ℚ)) +
    2 / (minusEightSieveBound : ℚ) ^ 2

def MinusEightIntervalBin.sieveExpression (row : MinusEightIntervalBin) (primes : List ℕ) : ℚ :=
  3 / 16 + row.errorBound * (row.inertCount primes + row.allowance) +
    row.coefficient * row.divisorLimit * row.inertTailBound primes

/-- The complete arithmetic check on a row, including the searched
divisor envelope and the exact margin computed from the bounded prime list. -/
structure MinusEightIntervalBin.Valid (row : MinusEightIntervalBin) (primes : List ℕ) : Prop where
  dyadic : row.dyadicExponent ≤ 2
  upperExponent : row.exponent < 60
  cubePos : 0 < row.cubeRoot
  cubeLower : row.cubeRoot ^ 3 ≤ 2 ^ row.exponent
  cubeUpper : 2 ^ row.exponent < (row.cubeRoot + 1) ^ 3
  cutoffExact : row.primeCutoff = min (minusEightSieveBound - 1)
      (Nat.sqrt (row.coefficientNumerator * row.divisorLimit * (5 * row.cubeRoot) /
        (row.coefficientDenominator * 18)))
  cutoffPos : 0 < row.primeCutoff
  cutoffBound : row.primeCutoff < minusEightSieveBound
  divisors : splitDivisorEnvelope 39 minusEightSplitPrefix row.targetUpper ≤ row.divisorLimit
  margin : row.sieveExpression primes < row.lower - 1 / 256

instance instDecidableMinusEightIntervalBinValid (row : MinusEightIntervalBin) (primes : List ℕ) : Decidable (row.Valid primes) :=
  let check := row.dyadicExponent ≤ 2 ∧ row.exponent < 60 ∧
    0 < row.cubeRoot ∧ row.cubeRoot ^ 3 ≤ 2 ^ row.exponent ∧
    2 ^ row.exponent < (row.cubeRoot + 1) ^ 3 ∧
    row.primeCutoff = min (minusEightSieveBound - 1)
      (Nat.sqrt (row.coefficientNumerator * row.divisorLimit * (5 * row.cubeRoot) /
        (row.coefficientDenominator * 18))) ∧
    0 < row.primeCutoff ∧ row.primeCutoff < minusEightSieveBound ∧
    splitDivisorEnvelope 39 minusEightSplitPrefix row.targetUpper ≤ row.divisorLimit ∧
    row.sieveExpression primes < row.lower - 1 / 256
  have h : check ↔ row.Valid primes := by
    constructor
    · rintro ⟨h₁, h₂, h₃, h₄, h₅, h₆, h₇, h₈, h₉, h₁₀⟩
      exact ⟨h₁, h₂, h₃, h₄, h₅, h₆, h₇, h₈, h₉, h₁₀⟩
    · rintro ⟨h₁, h₂, h₃, h₄, h₅, h₆, h₇, h₈, h₉, h₁₀⟩
      exact ⟨h₁, h₂, h₃, h₄, h₅, h₆, h₇, h₈, h₉, h₁₀⟩
  decidable_of_iff check h

private theorem minusEight_interval_bins_pass :
    let primes := minusEightIntervalPrimes
    minusEightIntervalBins.all (fun row => decide (row.Valid primes)) = true := by
  native_decide

/-- All 107 interval inequalities pass, with the divisor bounds and
integer square and cube root checks included. -/
theorem minusEight_interval_row_valid {row : MinusEightIntervalBin} (hrow : row ∈ minusEightIntervalBins) :
    row.Valid minusEightIntervalPrimes :=
  of_decide_eq_true (List.all_eq_true.mp minusEight_interval_bins_pass row hrow)

/-- The rows cover all dyadic bins above the finite scan and below the
inert-core onset. -/
theorem minusEight_interval_bin_certificate {e k : ℕ} (he : e ≤ 2)
    (hk : 27 ≤ e + k) (hk' : k < 60) :
    ∃ row ∈ minusEightIntervalBins, row.dyadicExponent = e ∧ row.exponent = k ∧
      row.Valid minusEightIntervalPrimes := by
  have hspan : (minusEightIntervalBins.map fun row => (row.dyadicExponent, row.exponent)) =
      ((List.range' 27 33).map fun k => (0, k)) ++
        ((List.range' 26 34).map fun k => (1, k)) ++
        ((List.range' 20 40).map fun k => (2, k)) := by decide +kernel
  have hm : (e, k) ∈ (minusEightIntervalBins.map fun row => (row.dyadicExponent, row.exponent)) := by
    rw [hspan]
    interval_cases e
    · exact List.mem_append_left _ (List.mem_append_left _ (List.mem_map.mpr
        ⟨k, List.mem_range'_1.mpr ⟨by omega, by omega⟩, rfl⟩))
    · exact List.mem_append_left _ (List.mem_append_right _ (List.mem_map.mpr
        ⟨k, List.mem_range'_1.mpr ⟨by omega, by omega⟩, rfl⟩))
    · exact List.mem_append_right _ (List.mem_map.mpr
        ⟨k, List.mem_range'_1.mpr ⟨by omega, by omega⟩, rfl⟩)
  obtain ⟨row, hr, heq⟩ := List.mem_map.mp hm
  exact ⟨row, hr, congrArg Prod.fst heq, congrArg Prod.snd heq, minusEight_interval_row_valid hr⟩

/-- A row's divisor envelope controls every target at most its upper endpoint. -/
theorem MinusEightIntervalBin.divisor_maximum_le {row : MinusEightIntervalBin}
    (hr : row.Valid minusEightIntervalPrimes) {n : ℕ} (hn : n ≤ row.targetUpper) :
    minusEightDivisorMaximum n ≤ row.divisorLimit := by
  have hX : row.targetUpper ≤ 2 ^ 62 := by
    unfold MinusEightIntervalBin.targetUpper
    have he : 2 ^ row.dyadicExponent ≤ 2 ^ 2 := Nat.pow_le_pow_right (by decide) hr.dyadic
    have hk : 2 ^ (row.exponent + 1) - 1 ≤ 2 ^ 60 :=
      (Nat.sub_le ..).trans (Nat.pow_le_pow_right (by decide) (by have := hr.upperExponent; omega))
    exact (Nat.mul_le_mul he hk).trans_eq (by decide +kernel)
  apply Finset.sup_le
  intro m hm
  obtain ⟨hm₀, hm₁⟩ := Finset.mem_Icc.mp hm
  exact (minusEight_split_divisors_le_envelope hm₀ (hm₁.trans hn) hX).trans hr.divisors

private theorem MinusEightIntervalBin.tail_primes_pos (row : MinusEightIntervalBin) :
    ∀ p ∈ row.tailPrimes minusEightIntervalPrimes, 0 < p := by
  intro p hp
  exact ((mem_minusEight_intervalPrimes.mp (List.mem_filter.mp hp).1).2.1).pos

private theorem rounded_partial_inverse_square_le_real (primes : List ℕ)
    (hp : ∀ p ∈ primes, 0 < p) :
    (primes.map inverseSquare).sum ≤
      (((primes.map roundedInverseSquareWeight).sum : ℕ) : ℝ) / inverseSquareFixedPointScale := by
  have h := (Rat.cast_le (K := ℝ)).mpr (inverse_square_sum_le_rounded primes hp)
  simp only [Rat.cast_list_sum, List.map_map, Function.comp_def, Rat.cast_div, Rat.cast_natCast,
    Rat.cast_pow, Rat.cast_one] at h
  exact h

/-- The computed prime count and fixed-point tail bound the prime
sums in the actual shell at any finite target. -/
theorem MinusEightIntervalBin.cutoff_bounds {row : MinusEightIntervalBin} (hr : row.Valid minusEightIntervalPrimes) :
    MinusEightCutoffBounds (row.primeCutoff : ℝ) (row.inertCount minusEightIntervalPrimes : ℝ)
      (row.inertTailBound minusEightIntervalPrimes : ℝ) := by
  constructor
  · intro s hs
    let primes := minusEightIntervalPrimes.filter fun p => decide (p ≤ row.primeCutoff)
    have hsub : s ⊆ primes.toFinset := by
      intro p hp
      obtain ⟨hi, hP⟩ := hs p hp
      have hle : p ≤ row.primeCutoff := by exact_mod_cast hP
      exact List.mem_toFinset.mpr (List.mem_filter.mpr
        ⟨mem_minusEight_intervalPrimes.mpr ⟨hle.trans hr.cutoffBound.le, hi⟩, decide_eq_true hle⟩)
    have hcard := Finset.card_le_card hsub
    rw [List.toFinset_card_of_nodup (minusEight_intervalPrimes_nodup.filter _)] at hcard
    exact_mod_cast hcard
  · intro s hs
    let p₀ := s.filter (fun p => p ≤ minusEightSieveBound)
    let p₁ := s.filter (fun p => ¬ p ≤ minusEightSieveBound)
    let primes := row.tailPrimes minusEightIntervalPrimes
    have hsub : p₀ ⊆ primes.toFinset := by
      intro p hp
      obtain ⟨h₁, h₂⟩ := Finset.mem_filter.mp hp
      have hlt : row.primeCutoff < p := by exact_mod_cast (hs p h₁).2
      exact List.mem_toFinset.mpr (List.mem_filter.mpr
        ⟨mem_minusEight_intervalPrimes.mpr ⟨h₂, (hs p h₁).1⟩, decide_eq_true hlt⟩)
    have hhead : ∑ p ∈ p₀, inverseSquare p ≤
        (((primes.map roundedInverseSquareWeight).sum : ℕ) : ℝ) / inverseSquareFixedPointScale := by
      refine (Finset.sum_le_sum_of_subset_of_nonneg hsub
        (fun p _ _ => inverseSquare_nonneg p)).trans ?_
      have hnd : primes.Nodup := minusEight_intervalPrimes_nodup.filter _
      rw [List.sum_toFinset inverseSquare hnd]
      exact rounded_partial_inverse_square_le_real primes row.tail_primes_pos
    have htail : ∑ p ∈ p₁, inverseSquare p ≤
        1 / (4 * (minusEightSieveBound : ℝ)) + 2 / (minusEightSieveBound : ℝ) ^ 2 := by
      refine minusEight_inert_tail_le (P := minusEightSieveBound) (X := minusEightSieveBound)
        (by norm_num [minusEightSieveBound]) (by linarith) ?_
      intro p hp
      obtain ⟨h₁, h₂⟩ := Finset.mem_filter.mp hp
      exact ⟨(hs p h₁).1, Nat.lt_of_not_ge h₂⟩
    rw [← Finset.sum_filter_add_sum_filter_not s (fun p => p ≤ minusEightSieveBound)]
    simp only [MinusEightIntervalBin.inertTailBound, Rat.cast_add, Rat.cast_div, Rat.cast_mul,
      Rat.cast_ofNat, Rat.cast_natCast, Rat.cast_pow, Rat.cast_one]
    exact (add_le_add hhead htail).trans_eq (add_assoc ..).symm

private theorem MinusEightIntervalBin.shell_error_le {row : MinusEightIntervalBin}
    (hr : row.Valid minusEightIntervalPrimes) {m : ℕ} (hm : 0 < m)
    (hI : ∀ p ∈ m.primeFactors, MinusEightInertPrime p) (hk : 2 ^ row.exponent ≤ m) :
    minusEightShellError m ≤ (row.errorBound : ℝ) := by
  let c : ℝ := row.cubeRoot
  have hc : 0 < c := Nat.cast_pos.mpr hr.cubePos
  have hcm : c ^ 3 ≤ (m : ℝ) := by
    change (row.cubeRoot : ℝ) ^ 3 ≤ (m : ℝ)
    exact_mod_cast hr.cubeLower.trans hk
  have hE := minusEight_inert_shell_error_power_lt hm hI
  have hE₀ := minusEight_shell_error_nonneg m
  have hpow : (minusEightShellError m * c) ^ 6 ≤ (18 / 5 : ℝ) ^ 6 := by
    rw [mul_pow, show (6 : ℕ) = 3 * 2 by decide, pow_mul c 3 2]
    exact (mul_le_mul_of_nonneg_left (pow_le_pow_left₀ (pow_nonneg hc.le 3) hcm 2)
      (pow_nonneg hE₀ 6)).trans hE.le
  have h := (pow_le_pow_iff_left₀ (mul_nonneg hE₀ hc.le) (by norm_num)
    (by decide : (6 : ℕ) ≠ 0)).mp hpow
  have hle := (le_div_iff₀ hc).mpr h
  simpa only [MinusEightIntervalBin.errorBound, Rat.cast_div, Rat.cast_mul, Rat.cast_ofNat,
    Rat.cast_natCast, ← div_div] using hle

private theorem MinusEightIntervalBin.cast_margin {row : MinusEightIntervalBin}
    (hr : row.Valid minusEightIntervalPrimes) :
    3 / 16 + (row.errorBound : ℝ) * ((row.inertCount minusEightIntervalPrimes : ℝ) +
      (row.allowance : ℝ)) + (row.coefficient : ℝ) * (row.divisorLimit : ℝ) *
        (row.inertTailBound minusEightIntervalPrimes : ℝ) < (row.lower : ℝ) := by
  have h := (Rat.cast_lt (K := ℝ)).mpr hr.margin
  simp only [MinusEightIntervalBin.sieveExpression, Rat.cast_sub, Rat.cast_add, Rat.cast_mul,
    Rat.cast_div, Rat.cast_ofNat, Rat.cast_natCast] at h
  linarith only [h]

/-- A checked row supplies the exact selected-shell margin at every
inert-supported odd part in its dyadic bin. -/
theorem MinusEightIntervalBin.selected_margin {row : MinusEightIntervalBin}
    (hr : row.Valid minusEightIntervalPrimes) {m : ℕ} (hm : 0 < m)
    (hI : ∀ p ∈ m.primeFactors, MinusEightInertPrime p)
    (hlo : 2 ^ row.exponent ≤ m) (hhi : m < 2 ^ (row.exponent + 1)) :
    MinusEightSelectedMarginWithBounds row.dyadicExponent m (row.inertCount minusEightIntervalPrimes)
      (row.inertTailBound minusEightIntervalPrimes : ℝ) := by
  have hD : (minusEightDivisorMaximum (2 ^ row.dyadicExponent * m) : ℝ) ≤ row.divisorLimit := by
    exact_mod_cast MinusEightIntervalBin.divisor_maximum_le hr (Nat.mul_le_mul_left _ (Nat.le_sub_one_of_lt hhi))
  have hE := MinusEightIntervalBin.shell_error_le hr hm hI hlo
  have hE₀ := minusEight_shell_error_nonneg m
  have hC : (0 : ℝ) ≤ row.coefficient := by unfold MinusEightIntervalBin.coefficient; positivity
  have hδ : (0 : ℝ) ≤ row.allowance := by unfold MinusEightIntervalBin.allowance; split <;> norm_num
  have hT : (0 : ℝ) ≤ row.inertTailBound minusEightIntervalPrimes := by
    unfold MinusEightIntervalBin.inertTailBound
    positivity
  have hmargin :
      3 / 16 + minusEightShellError m * ((row.inertCount minusEightIntervalPrimes : ℝ) +
        (row.allowance : ℝ)) + (row.coefficient : ℝ) *
        minusEightDivisorMaximum (2 ^ row.dyadicExponent * m) *
          (row.inertTailBound minusEightIntervalPrimes : ℝ) < (row.lower : ℝ) := by
    refine lt_of_le_of_lt ?_ (MinusEightIntervalBin.cast_margin hr)
    gcongr
  unfold MinusEightSelectedMarginWithBounds
  unfold MinusEightIntervalBin.allowance MinusEightIntervalBin.coefficient
    MinusEightIntervalBin.coefficientNumerator MinusEightIntervalBin.coefficientDenominator
    MinusEightIntervalBin.lower at hmargin
  by_cases h₀ : row.dyadicExponent = 0
  · norm_num [h₀] at hmargin ⊢
    exact hmargin
  · by_cases h₁ : row.dyadicExponent = 1
    · norm_num [h₀, h₁] at hmargin ⊢
      exact hmargin
    · have h₂ : row.dyadicExponent = 2 := by have := hr.dyadic; omega
      norm_num [h₂] at hmargin ⊢
      exact hmargin

/-- The bounded arithmetic gap between the finite scan and the inert-core
onset is covered by the verified interval rows. -/
theorem minusEight_certified_interval_cutoff {e m : ℕ} (he : e ≤ 2)
    (hm : 0 < m) (hI : ∀ p ∈ m.primeFactors, MinusEightInertPrime p)
    (hscan : 2 ^ 27 < 2 ^ e * m) (hlarge : m < 2 ^ 60) :
    ∃ X count tail : ℝ, MinusEightCutoffBounds X count tail ∧
      MinusEightSelectedMarginWithBounds e m count tail := by
  let k := Nat.log 2 m
  have hk₀ : 2 ^ k ≤ m := Nat.pow_log_le_self _ hm.ne'
  have hk₁ : m < 2 ^ (k + 1) := Nat.lt_pow_succ_log_self (by decide) m
  have hk : k < 60 := by
    have h := hk₀.trans_lt hlarge
    exact (Nat.pow_lt_pow_iff_right (by decide : 1 < 2)).mp h
  have hstart : 27 ≤ e + k := by
    by_contra! h
    have huv : 2 ^ e * m < 2 ^ (e + (k + 1)) := by
      rw [pow_add]
      exact Nat.mul_lt_mul_of_pos_left hk₁ (by positivity)
    exact (not_lt_of_ge (Nat.pow_le_pow_right (by decide : 0 < 2) (by omega : e + (k + 1) ≤ 27)))
      (hscan.trans huv)
  obtain ⟨row, _, heq, hkq, hvr⟩ := minusEight_interval_bin_certificate he hstart hk
  refine ⟨_, _, _, MinusEightIntervalBin.cutoff_bounds hvr, ?_⟩
  rw [← heq]
  exact MinusEightIntervalBin.selected_margin hvr hm hI (hkq ▸ hk₀) (hkq ▸ hk₁)

end QuadraticNorms
