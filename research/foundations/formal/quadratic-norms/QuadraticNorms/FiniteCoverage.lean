import QuadraticNorms.EisensteinReplay
import QuadraticNorms.MinusSevenReplay
import QuadraticNorms.MinusEightReplay
import QuadraticNorms.MinusElevenReplay

namespace QuadraticNorms

/-- The Eisenstein primitive sumset on the complete finite range of the paper. -/
theorem P3_sum_below_cutoff {n : ℕ} (hn : n < 100000000) :
    InSumset P3 P3 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ n % 18 ≠ 0 ∧ n % 18 ≠ 12 :=
  P3_sum_classification_up_to n restoration3
    (fun _ hm h => eisensteinOrdinary_finite (hm.trans_lt hn) h)
    (fun _ hm h => eisensteinMixed_finite (hm.trans_lt hn) h) n le_rfl

/-- The discriminant −7 primitive sumset through its scan cutoff. -/
theorem P7_sum_through_cutoff {n : ℕ} (hn : n ≤ 2 ^ 24) :
    InSumset P7 P7 n ↔ 2 ≤ n ∧ ¬ 49 ∣ n ∧ n ≠ 7 ∧ n ≠ 77 ∧ n ≠ 203 :=
  P7_sum_classification_up_to (2 ^ 24) restoration7
    (fun _ => C7_finite_core_coverage) n hn

/-- The discriminant −8 primitive sumset through its scan cutoff. -/
theorem P8_sum_through_cutoff {n : ℕ} (hn : n ≤ 2 ^ 27) :
    InSumset P8 P8 n ↔ 2 ≤ n ∧ ¬ 16 ∣ n :=
  P8_sum_classification_up_to (2 ^ 27) restoration8
    (fun _ hm _ => C8_finite_ordinary_coverage hm)
    (fun _ => C8_finite_mixed_coverage) n hn

/-- The discriminant −11 primitive sumset through its scan cutoff. -/
theorem P11_sum_through_cutoff {n : ℕ} (hn : n ≤ 2 ^ 25) :
    InSumset P11 P11 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 121 ∣ n :=
  P11_sum_classification_up_to (2 ^ 25) restoration11
    (fun _ => C11_finite_core_coverage) n hn

/-- Only these reduced targets and split-prime boundaries remain after the
complete finite scans. Each cutoff matches the corresponding scan. -/
structure BeyondFiniteCoreCoverage : Prop where
  eisensteinOrdinary : ∀ n ≥ 100000000, n % 6 = 2 → IsReduced C3 n → InSumset C3 C3 n
  eisensteinMixed : ∀ n ≥ 100000000, n % 6 = 4 → IsReduced C3 n → InSumset C3 (Dilate 3 C3) n
  seven : ∀ n > 2 ^ 24, ¬ 7 ∣ n → IsReduced C7 n → InSumset C7 C7 n
  sevenPrimes : ∀ p > 2 ^ 24, p.Prime → C7 p → InSumset C7 C7 p
  eightOrdinary : ∀ n > 2 ^ 27, n % 2 = 0 → ¬ 8 ∣ n → IsReduced C8 n → InSumset C8 C8 n
  eightMixed : ∀ n > 2 ^ 27, n % 2 = 1 → IsReduced C8 n → InSumset C8 (Dilate 2 C8) n
  eightPrimes : ∀ p > 2 ^ 27, p.Prime → C8 p → InSumset C8 (Dilate 2 C8) p
  eleven : ∀ n > 2 ^ 25, n % 2 = 0 → ¬ 11 ∣ n → IsReduced C11 n → InSumset C11 C11 n

theorem reducedCoreCoverage_of_beyondFinite (h : BeyondFiniteCoreCoverage) : ReducedCoreCoverage := by
  constructor
  · intro n _ hclass hred
    exact if hn : n < 100000000 then eisensteinOrdinary_finite hn hclass
      else h.eisensteinOrdinary n (Nat.le_of_not_gt hn) hclass hred
  · intro n _ hclass hred
    exact if hn : n < 100000000 then eisensteinMixed_finite hn hclass
      else h.eisensteinMixed n (Nat.le_of_not_gt hn) hclass hred
  · intro n h₂ h₇ h₁₁ h₂₉ hred
    exact if hn : n ≤ 2 ^ 24 then C7_finite_core_coverage hn h₂ h₇ h₁₁ h₂₉
      else h.seven n (Nat.lt_of_not_ge hn) h₇ hred
  · intro p hp hc h₁₁ h₂₉
    by_cases hn : p ≤ 2 ^ 24
    · have h₇ : ¬ 7 ∣ p := fun hd => by
        have he := (Nat.dvd_prime hp).mp hd
        rcases he with h | rfl
        · omega
        · exact (by decide : ¬ (7 % 7 = 1 ∨ 7 % 7 = 2 ∨ 7 % 7 = 4)) (hc.2 7 (by decide) dvd_rfl)
      exact C7_finite_core_coverage hn hp.two_le h₇ h₁₁ h₂₉
    · exact h.sevenPrimes p (Nat.lt_of_not_ge hn) hp hc
  · intro n _ he h₈ hred
    exact if hn : n ≤ 2 ^ 27 then C8_finite_ordinary_coverage hn he h₈
      else h.eightOrdinary n (Nat.lt_of_not_ge hn) he h₈ hred
  · intro n h₃ ho hred
    exact if hn : n ≤ 2 ^ 27 then C8_finite_mixed_coverage hn h₃ ho
      else h.eightMixed n (Nat.lt_of_not_ge hn) ho hred
  · intro p hp hc
    by_cases hn : p ≤ 2 ^ 27
    · have hres := C8_mod_eight hc
      have h₂ := hp.two_le
      exact C8_finite_mixed_coverage hn (by omega) (by omega)
    · exact h.eightPrimes p (Nat.lt_of_not_ge hn) hp hc
  · intro n h₂ he h₁₁ hred
    exact if hn : n ≤ 2 ^ 25 then C11_finite_core_coverage hn h₂ he h₁₁
      else h.eleven n (Nat.lt_of_not_ge hn) he h₁₁ hred

/-- The four exact primitive classifications from coverage above the scan cutoffs. -/
theorem exact_overview_of_beyond_finite (h : BeyondFiniteCoreCoverage) (n : ℕ) :
    (InSumset P3 P3 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ n % 18 ≠ 0 ∧ n % 18 ≠ 12) ∧
    (InSumset P7 P7 n ↔ 2 ≤ n ∧ ¬ 49 ∣ n ∧ n ≠ 7 ∧ n ≠ 77 ∧ n ≠ 203) ∧
    (InSumset P8 P8 n ↔ 2 ≤ n ∧ ¬ 16 ∣ n) ∧
    (InSumset P11 P11 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 121 ∣ n) :=
  exact_overview_of_reduced_coverage (reducedCoreCoverage_of_beyondFinite h) n

end QuadraticNorms
