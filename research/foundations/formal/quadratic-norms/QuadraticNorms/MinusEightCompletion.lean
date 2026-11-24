import QuadraticNorms.MinusEightAsymptotic
import QuadraticNorms.MinusEightIntervals
import QuadraticNorms.MinusEightReplay
import QuadraticNorms.CoreDescent

open scoped ArithmeticFunction.sigma

namespace QuadraticNorms

/-- The three finite interval inequalities for the selected shells of
norm m, 2m, and 4m. -/
def MinusEightSelectedMargin (e m : ℕ) (X : ℝ) : Prop :=
  MinusEightSelectedMarginWithBounds e m (X / 4 + 1) (1 / (4 * X) + 2 / X ^ 2)

private theorem minusEight_reduced_odd_part_inert {e m : ℕ} (hm : Odd m)
    (hred : IsReduced C8 (2 ^ e * m)) :
    ∀ p ∈ m.primeFactors, MinusEightInertPrime p := by
  intro p hp
  have hprime := Nat.prime_of_mem_primeFactors hp
  have hd := Nat.dvd_of_mem_primeFactors hp
  have hodd := hm.of_dvd_nat hd
  have h₂ := Nat.odd_iff.mp hodd
  have hnot : ¬ (p % 8 = 1 ∨ p % 8 = 3) :=
    fun h => hred p hprime (hd.trans (dvd_mul_left _ _)) ((primeCore_prime hprime).mpr h)
  have hmod := Nat.mod_mod_of_dvd p (by decide : 2 ∣ 8)
  exact ⟨hprime, by omega⟩

/-- The unrestricted −8 sumset follows from the two theta identities,
the cusp coefficient bound, the unit-level Hecke estimate, and the
selected-shell margins in the interval between the finite scan and
the inert-core onset. The split-prime boundary is included. -/
private theorem minusEight_sumsets_of_theta_unit_and_cutoff_bounds {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    (hcusp : ∀ m, Odd m → |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m)
    (hsmall : ∀ e m p : ℕ, e ≤ 3 → Odd m → p ≤ 2 ^ e * m → MinusEightInertPrime p →
      ¬ p ∣ 2 ^ e * m → ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card)
    (hIntervals : ∀ e m : ℕ, e ≤ 2 → Odd m → (∀ p ∈ m.primeFactors, MinusEightInertPrime p) →
      2 ^ 27 < 2 ^ e * m → m < 2 ^ 60 → ∃ X count tail : ℝ, MinusEightCutoffBounds X count tail ∧
        MinusEightSelectedMarginWithBounds e m count tail)
    (n : ℕ) :
    (InSumset C8 C8 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 8 ∣ n) ∧
      (InSumset C8 (Dilate 2 C8) n ↔ 3 ≤ n ∧ n % 2 = 1) ∧
      (InSumset P8 P8 n ↔ 2 ≤ n ∧ ¬ 16 ∣ n) := by
  have hordinary (u : ℕ) (hu : 2 ≤ u) (heven : u % 2 = 0)
      (h₈ : ¬ 8 ∣ u) (hred : IsReduced C8 u) : InSumset C8 C8 u := by
    by_cases hscan : u ≤ 2 ^ 27
    · exact C8_finite_ordinary_coverage hscan heven h₈
    obtain ⟨e, m, hm, rfl⟩ := Nat.exists_eq_two_pow_mul_odd (by omega : u ≠ 0)
    have he : e = 1 ∨ e = 2 := by
      have he₀ : e ≠ 0 := by
        rintro rfl
        rw [pow_zero, one_mul, Nat.odd_iff.mp hm] at heven
        omega
      have he₃ : e < 3 := by
        by_contra! h
        exact h₈ ((Nat.pow_dvd_pow 2 h).trans (dvd_mul_right _ _))
      omega
    have hI := minusEight_reduced_odd_part_inert hm hred
    by_cases hlarge : 2 ^ 60 ≤ m
    · obtain ⟨_, h₁, h₂⟩ := minusEight_large_inert_core_sums hθ hcusp hsmall hm hI hlarge
      rcases he with rfl | rfl
      · simpa using h₁
      · exact h₂
    obtain ⟨X, count, tail, hcut, hmargin⟩ := hIntervals e m (by omega) hm hI (by omega) (by omega)
    apply minusEight_core_sum_of_selected_cutoff_bounds hθ he hm hcut
      (fun p hp hi _ hd => hsmall e m p (by omega) hm hp hi hd)
    rcases he with rfl | rfl <;> norm_num [MinusEightSelectedMarginWithBounds, minusEightShellFactor] at hmargin ⊢ <;>
      exact hmargin
  have hmixed (u : ℕ) (hu : 3 ≤ u) (hodd : u % 2 = 1)
      (hred : IsReduced C8 u) : InSumset C8 (Dilate 2 C8) u := by
    by_cases hscan : u ≤ 2 ^ 27
    · exact C8_finite_mixed_coverage hscan hu hodd
    have hm := Nat.odd_iff.mpr hodd
    have hI := minusEight_reduced_odd_part_inert (e := 0) hm (by simpa using hred)
    by_cases hlarge : 2 ^ 60 ≤ u
    · exact (minusEight_large_inert_core_sums hθ hcusp hsmall hm hI hlarge).1
    obtain ⟨X, count, tail, hcut, hmargin⟩ := hIntervals 0 u (by decide) hm hI
      (by simpa using Nat.lt_of_not_ge hscan) (by omega)
    apply minusEight_mixed_core_of_selected_cutoff_bounds hθ hm (hcusp u hm) hcut
    · intro p hp hi _ hd
      simpa only [pow_zero, one_mul] using hsmall 0 u p (by decide) hm
        (by simpa using hp) hi (by simpa using hd)
    · simpa [MinusEightSelectedMarginWithBounds] using hmargin
  have hprime (p : ℕ) (hp : p.Prime) (hcore : C8 p) : InSumset C8 (Dilate 2 C8) p := by
    by_cases hscan : p ≤ 2 ^ 27
    · have hmod := C8_mod_eight hcore
      have hp₂ := hp.two_le
      exact C8_finite_mixed_coverage hscan (by omega) (by omega)
    · have hpOdd := hp.odd_of_ne_two (by omega)
      apply minusEight_large_prime_mixed_core hθ hp (by omega) (hcusp p hpOdd)
      intro q hq hi hd
      simpa only [pow_zero, one_mul] using hsmall 0 p q (by decide) hpOdd
        (by simpa using hq) hi (by simpa using hd)
  have ho (v : ℕ) (hn : 2 ≤ v) (he : v % 2 = 0) (h₈ : ¬ 8 ∣ v) :=
    C8_ordinary_core_coverage hordinary hn he h₈
  have hm (v : ℕ) (hn : 3 ≤ v) (ho : v % 2 = 1) :=
    C8_mixed_core_coverage hmixed hprime hn ho
  exact ⟨⟨C8_ordinary_necessary, fun ⟨hn, he, h₈⟩ => ho n hn he h₈⟩,
    ⟨C8_mixed_necessary, fun ⟨hn, hodd⟩ => hm n hn hodd⟩,
    P8_sum_classification restoration8 ho hm n⟩

/-- The unrestricted −8 classification using the progression-only form
of the remaining selected-shell inequalities. -/
theorem P8_sum_of_theta_unit_and_intervals {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    (hcusp : ∀ m, Odd m → |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m)
    (hsmall : ∀ e m p : ℕ, e ≤ 3 → Odd m → p ≤ 2 ^ e * m → MinusEightInertPrime p →
      ¬ p ∣ 2 ^ e * m → ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card)
    (hIntervals : ∀ e m : ℕ, e ≤ 2 → Odd m → (∀ p ∈ m.primeFactors, MinusEightInertPrime p) →
      2 ^ 27 < 2 ^ e * m → m < 2 ^ 60 → ∃ X : ℝ, 0 < X ∧ MinusEightSelectedMargin e m X)
    (n : ℕ) : InSumset P8 P8 n ↔ 2 ≤ n ∧ ¬ 16 ∣ n := by
  refine (minusEight_sumsets_of_theta_unit_and_cutoff_bounds hθ hcusp hsmall ?_ n).2.2
  intro e m he hm hI hscan hlarge
  obtain ⟨X, hX, hmargin⟩ := hIntervals e m he hm hI hscan hlarge
  exact ⟨_, _, _, minusEight_cutoff_bounds_progressions hX, hmargin⟩

/-- The unrestricted −8 classification with all arithmetic and
certificate checks discharged. The three inputs state the modular
coefficient identities, the cusp bound, and the unit-level Hecke bound. -/
theorem minusEight_exact_sumsets_of_theta_and_unit {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    (hcusp : ∀ m, Odd m → |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m)
    (hsmall : ∀ e m p : ℕ, e ≤ 3 → Odd m → p ≤ 2 ^ e * m → MinusEightInertPrime p →
      ¬ p ∣ 2 ^ e * m → ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card)
    (n : ℕ) :
    (InSumset C8 C8 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 8 ∣ n) ∧
      (InSumset C8 (Dilate 2 C8) n ↔ 3 ≤ n ∧ n % 2 = 1) ∧
      (InSumset P8 P8 n ↔ 2 ≤ n ∧ ¬ 16 ∣ n) := by
  refine minusEight_sumsets_of_theta_unit_and_cutoff_bounds hθ hcusp hsmall ?_ n
  intro e m he hm hI hscan hlarge
  have hm₀ : 0 < m := by obtain ⟨k, hk⟩ := hm; omega
  exact minusEight_certified_interval_cutoff he hm₀ hI hscan hlarge

/-- The unrestricted primitive row with the arithmetic inputs discharged. -/
theorem P8_sum_of_theta_and_unit {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    (hcusp : ∀ m, Odd m → |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m)
    (hsmall : ∀ e m p : ℕ, e ≤ 3 → Odd m → p ≤ 2 ^ e * m → MinusEightInertPrime p →
      ¬ p ∣ 2 ^ e * m → ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card)
    (n : ℕ) : InSumset P8 P8 n ↔ 2 ≤ n ∧ ¬ 16 ∣ n := by
  exact (minusEight_exact_sumsets_of_theta_and_unit hθ hcusp hsmall n).2.2

end QuadraticNorms
