import QuadraticNorms.MinusElevenIntervals
import QuadraticNorms.SplitDivisorFrontier

namespace QuadraticNorms

def MinusElevenSplitPrime (p : ℕ) : Prop :=
  p.Prime ∧ p % 11 ∈ [1, 3, 4, 5, 9]

instance : DecidablePred MinusElevenSplitPrime := fun p => by unfold MinusElevenSplitPrime; infer_instance

private def minusElevenSplitPrefix : List ℕ :=
  [3, 5, 23, 31, 37, 47, 53, 59, 67, 71, 89, 97]

private theorem minusElevenSplitPrefix_covers {p : ℕ} (hp : MinusElevenSplitPrime p) (hle : p ≤ 97) :
    p ∈ minusElevenSplitPrefix := by
  have h : ∀ q : Fin 98, MinusElevenSplitPrime q.val → q.val ∈ minusElevenSplitPrefix := by decide +kernel
  exact h ⟨p, by omega⟩ hp

/-- The product of one plus the valuation at each split prime dividing n. -/
def minusElevenSplitDivisors (n : ℕ) : ℕ :=
  ∏ p ∈ n.primeFactors.filter MinusElevenSplitPrime, (n.factorization p + 1)

private theorem split_divisors_le_envelope {n X : ℕ} (hn : 0 < n)
    (hnX : n ≤ X) (hX : X ≤ 2 ^ 61) :
    minusElevenSplitDivisors n ≤ splitDivisorEnvelope 61 minusElevenSplitPrefix X := by
  apply split_divisors_le_prefix_envelope MinusElevenSplitPrime minusElevenSplitPrefix (base := 2)
    (last := 97) hn hnX (by decide +kernel) (by decide +kernel)
    (fun _ hp hsplit => minusElevenSplitPrefix_covers hsplit hp) (by decide)
    (fun _ hp _ => hp.two_le)
  · exact hX.trans_lt (by decide +kernel)
  · exact hX.trans_lt (by decide +kernel)

private def verifyDivisorLimit (row : MinusElevenBin) : Bool :=
  decide (splitDivisorEnvelope 61 minusElevenSplitPrefix (2 ^ (row.exponent + 1) - 1) ≤ row.divisorLimit)

private theorem divisorLimits_pass : minusElevenBins.all verifyDivisorLimit = true := by
  native_decide

/-- The divisor frontier verifies the stated upper bound on the split
divisor product at every positive integer in the bin's range. -/
theorem minusEleven_bin_divisor_limit {row : MinusElevenBin} (hr : row ∈ minusElevenBins)
    {n : ℕ} (hn : 0 < n) (hn' : n < 2 ^ (row.exponent + 1)) :
    minusElevenSplitDivisors n ≤ row.divisorLimit := by
  have hcheck := of_decide_eq_true (List.all_eq_true.mp divisorLimits_pass row hr)
  have hrows : minusElevenBins.map MinusElevenBin.exponent = List.range' 25 35 := by decide
  have hk : row.exponent < 60 := by
    have hmem : row.exponent ∈ minusElevenBins.map MinusElevenBin.exponent :=
      List.mem_map.mpr ⟨row, hr, rfl⟩
    rw [hrows, List.mem_range'_1] at hmem
    exact hmem.2
  have hbound : 2 ^ (row.exponent + 1) - 1 ≤ 2 ^ 61 :=
    (Nat.sub_le ..).trans (Nat.pow_le_pow_right (by decide) (by omega))
  exact (split_divisors_le_envelope hn (by omega) hbound).trans hcheck

/-- The interval sieve margin with the divisor-product bound discharged. -/
theorem minusEleven_bin_sieve_margin {row : MinusElevenBin}
    (hr : row ∈ minusElevenBins) (hv : row.valid)
    {n : ℕ} (hn : 0 < n) (hn' : n < 2 ^ (row.exponent + 1))
    {base tail error count : ℚ}
    (hbase : base ≤ 9 / 100) (htail : 0 ≤ tail) (htail' : tail ≤ row.inertTailBound)
    (herr : 0 ≤ error) (herr' : error ≤ row.errorBound)
    (hcount : 0 ≤ count) (hcount' : count ≤ row.inertCount) :
    base + (8 / 3) * minusElevenSplitDivisors n * tail + error * (3 + (5 / 2) * count) < 2 / 3 := by
  apply MinusElevenBin.sieve_margin hv hbase _ htail htail' herr herr' hcount hcount'
  exact_mod_cast minusEleven_bin_divisor_limit hr hn hn'

end QuadraticNorms
