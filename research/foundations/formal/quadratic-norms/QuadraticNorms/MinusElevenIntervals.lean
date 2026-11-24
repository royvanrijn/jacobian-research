import QuadraticNorms.SieveCriteria
import QuadraticNorms.RoundedPrimeTails
import Mathlib.Data.Nat.Prime.Basic
import Mathlib.Data.Nat.Sqrt

namespace QuadraticNorms

/-- One dyadic bin in the −11 certificate. -/
structure MinusElevenBin where
  exponent : ℕ
  divisorLimit : ℕ
  cubeRoot : ℕ
  primeCutoff : ℕ
  deriving DecidableEq

/-- The certificate supplies the divisor limits and the proposed integer
roots. The prime counts and fixed-point tails are computed from primes. -/
def minusElevenBins : List MinusElevenBin := [
  ⟨25, 120, 322, 122⟩,
  ⟨26, 128, 406, 141⟩,
  ⟨27, 160, 512, 178⟩,
  ⟨28, 168, 645, 205⟩,
  ⟨29, 192, 812, 245⟩,
  ⟨30, 200, 1024, 281⟩,
  ⟨31, 256, 1290, 357⟩,
  ⟨32, 288, 1625, 426⟩,
  ⟨33, 320, 2048, 504⟩,
  ⟨34, 336, 2580, 579⟩,
  ⟨35, 400, 3250, 710⟩,
  ⟨36, 480, 4096, 873⟩,
  ⟨37, 512, 5160, 1012⟩,
  ⟨38, 576, 6501, 1205⟩,
  ⟨39, 640, 8192, 1426⟩,
  ⟨40, 768, 10321, 1753⟩,
  ⟨41, 800, 13003, 2008⟩,
  ⟨42, 960, 16384, 2469⟩,
  ⟨43, 1024, 20642, 2863⟩,
  ⟨44, 1152, 26007, 3408⟩,
  ⟨45, 1280, 32768, 4033⟩,
  ⟨46, 1536, 41285, 4959⟩,
  ⟨47, 1600, 52015, 5681⟩,
  ⟨48, 1920, 65536, 6986⟩,
  ⟨49, 2048, 82570, 8098⟩,
  ⟨50, 2304, 104031, 9642⟩,
  ⟨51, 2560, 131072, 11408⟩,
  ⟨52, 3072, 165140, 14027⟩,
  ⟨53, 3200, 208063, 16070⟩,
  ⟨54, 3584, 262144, 19089⟩,
  ⟨55, 4096, 330280, 22907⟩,
  ⟨56, 4608, 416127, 27272⟩,
  ⟨57, 5120, 524288, 32267⟩,
  ⟨58, 6144, 660561, 39676⟩,
  ⟨59, 6400, 832255, 45453⟩]

def minusElevenInertPrimesThrough (bound : ℕ) : List ℕ :=
  (List.range (bound + 1)).filter (fun p =>
    decide (p.Prime ∧ p > 2 ∧ p % 11 ∈ [2, 6, 7, 8, 10]))

theorem mem_minusEleven_inertPrimes {p bound : ℕ} :
    p ∈ minusElevenInertPrimesThrough bound ↔
      p ≤ bound ∧ p.Prime ∧ p > 2 ∧ p % 11 ∈ [2, 6, 7, 8, 10] := by
  simp [minusElevenInertPrimesThrough, Nat.lt_succ_iff]

def minusElevenSieveBound : ℕ := 65536
private def countThrough (primes : List ℕ) (cutoff : ℕ) : ℕ :=
  (primes.filter (fun p => decide (p ≤ cutoff))).length

private def tailBound (primes : List ℕ) (cutoff : ℕ) : ℚ :=
  let tail := primes.filter (fun p => decide (cutoff < p))
  let sum : ℕ := (tail.map roundedInverseSquareWeight).sum
  (sum : ℚ) / (inverseSquareFixedPointScale : ℚ) + 5 / (22 * (minusElevenSieveBound : ℚ)) + 5 / (minusElevenSieveBound : ℚ) ^ 2

def MinusElevenBin.inertCount (row : MinusElevenBin) : ℕ :=
  countThrough (minusElevenInertPrimesThrough minusElevenSieveBound) row.primeCutoff

def MinusElevenBin.inertTailBound (row : MinusElevenBin) : ℚ :=
  tailBound (minusElevenInertPrimesThrough minusElevenSieveBound) row.primeCutoff

/-- Rounding up each term bounds the contribution through the sieve
limit. The added allowance is the manuscript's bound for the later primes. -/
theorem MinusElevenBin.rounded_tail_bound (row : MinusElevenBin) :
    let primes := (minusElevenInertPrimesThrough minusElevenSieveBound).filter
      (fun p => decide (row.primeCutoff < p))
    (primes.map (fun (p : ℕ) => (1 : ℚ) / (p : ℚ) ^ 2)).sum +
      5 / (22 * (minusElevenSieveBound : ℚ)) + 5 / (minusElevenSieveBound : ℚ) ^ 2 ≤
      row.inertTailBound := by
  intro primes
  have hpos (p : ℕ) (hp : p ∈ primes) : 0 < p :=
    ((mem_minusEleven_inertPrimes.mp (List.mem_filter.mp hp).1).2.1).pos
  exact add_le_add_left (add_le_add_left (inverse_square_sum_le_rounded primes hpos) _) _

def MinusElevenBin.errorBound (row : MinusElevenBin) : ℚ := 11 / (4 * (row.cubeRoot : ℚ))

def MinusElevenBin.margin (row : MinusElevenBin) : ℚ :=
  9 / 100 + (8 / 3) * row.divisorLimit * row.inertTailBound +
    row.errorBound * (3 + (5 / 2) * row.inertCount)

def MinusElevenBin.valid (row : MinusElevenBin) : Prop :=
  0 < row.cubeRoot ∧
    row.cubeRoot ^ 3 ≤ 2 ^ row.exponent ∧ 2 ^ row.exponent < (row.cubeRoot + 1) ^ 3 ∧
    row.primeCutoff = Nat.sqrt (64 * row.divisorLimit * row.cubeRoot / 165) ∧
    0 < row.primeCutoff ∧ row.primeCutoff < minusElevenSieveBound ∧
    row.margin < 2 / 3 - 1 / 1024

instance (row : MinusElevenBin) : Decidable row.valid := by unfold MinusElevenBin.valid; infer_instance

private theorem bins_pass :
    minusElevenBins.all (fun row => decide row.valid) = true := by
  native_decide

/-- Every dyadic bin from 25 through 59 is present and passes the
rational inequality using the stated divisor limit. -/
theorem minusEleven_bin_certificate {k : ℕ} (hk : 25 ≤ k) (hk' : k < 60) :
    ∃ row ∈ minusElevenBins, row.exponent = k ∧ row.valid := by
  have hspan : minusElevenBins.map MinusElevenBin.exponent = List.range' 25 35 := by decide
  have hmem : k ∈ List.range' 25 35 := by
    apply List.mem_range'_1.mpr
    omega
  rw [← hspan] at hmem
  obtain ⟨row, hr, heq⟩ := List.mem_map.mp hmem
  exact ⟨row, hr, heq, of_decide_eq_true (List.all_eq_true.mp bins_pass row hr)⟩

/-- The dyadic interval margin implies the sieve inequality once the
arithmetic divisor, error and inert-tail estimates are supplied. -/
theorem MinusElevenBin.sieve_margin {row : MinusElevenBin} (hr : row.valid)
    {base divisors tail error count : ℚ}
    (hbase : base ≤ 9 / 100)
    (hdivisors' : divisors ≤ row.divisorLimit)
    (htail : 0 ≤ tail) (htail' : tail ≤ row.inertTailBound)
    (herror : 0 ≤ error) (herror' : error ≤ row.errorBound)
    (hcount : 0 ≤ count) (hcount' : count ≤ row.inertCount) :
    base + (8 / 3) * divisors * tail + error * (3 + (5 / 2) * count) < 2 / 3 := by
  have hdnonneg : (0 : ℚ) ≤ row.divisorLimit := Nat.cast_nonneg _
  have hcountnonneg : (0 : ℚ) ≤ 3 + (5 / 2) * count := by positivity
  have hfirst := mul_le_mul (mul_le_mul_of_nonneg_left hdivisors' (by norm_num : (0 : ℚ) ≤ 8 / 3))
    htail' htail (mul_nonneg (by norm_num) hdnonneg)
  have hsecond := mul_le_mul herror'
    (add_le_add_right (mul_le_mul_of_nonneg_left hcount' (by norm_num : (0 : ℚ) ≤ 5 / 2)) 3)
    hcountnonneg (le_trans herror herror')
  have h := hr.2.2.2.2.2.2
  unfold MinusElevenBin.margin at h
  linarith

/-- The rational comparison used on the range starting at 2^60,
once the two negative powers have been bounded. -/
theorem minusEleven_large_tail_margin {u v : ℚ} (hu : u ≤ 1 / 32) (hv : v ≤ 1 / (2 ^ 20)) :
    9 / 100 + (3065 / 528) * u + (5371 / 48) * v < 2 / 3 := by
  norm_num at hv ⊢
  linarith

end QuadraticNorms
