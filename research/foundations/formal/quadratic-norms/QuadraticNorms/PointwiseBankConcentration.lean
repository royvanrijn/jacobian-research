import QuadraticNorms.PointwiseBankSuccess

open scoped BigOperators Classical

namespace QuadraticNorms

/-- The mass of one set of successes in independent Bernoulli trials. -/
noncomputable def bankPatternMass {ι : Type*} [DecidableEq ι] (s : Finset ι) (q : ι → ℝ)
    (A : Finset ι) : ℝ :=
  (∏ i ∈ A, q i) * ∏ i ∈ s \ A, (1 - q i)

theorem bank_pattern_mass_sum {ι : Type*} [DecidableEq ι] (s : Finset ι) (q : ι → ℝ) :
    ∑ A ∈ s.powerset, bankPatternMass s q A = 1 := by
  unfold bankPatternMass
  rw [← Finset.prod_add]
  simp only [add_sub_cancel, Finset.prod_const_one]

private theorem bank_pattern_mass_nonneg {ι : Type*} [DecidableEq ι] {s A : Finset ι} {q : ι → ℝ}
    (hA : A ⊆ s) (hq : ∀ i ∈ s, 0 ≤ q i ∧ q i ≤ 1) : 0 ≤ bankPatternMass s q A := by
  unfold bankPatternMass
  exact mul_nonneg (Finset.prod_nonneg (fun i hi => (hq i (hA hi)).1))
    (Finset.prod_nonneg (fun i hi => sub_nonneg.mpr (hq i (Finset.mem_sdiff.mp hi).1).2))

/-- Mass of the records with fewer than one sixteenth of the banks
successful. -/
noncomputable def bankDiscardMass {ι : Type*} [DecidableEq ι] (s : Finset ι) (q : ι → ℝ) : ℝ :=
  ∑ A ∈ s.powerset, if (A.card : ℝ) < (s.card : ℝ) / 16 then bankPatternMass s q A else 0

private theorem bank_pattern_exp_sum {ι : Type*} [DecidableEq ι] (s : Finset ι) (q : ι → ℝ) :
    (∑ A ∈ s.powerset, bankPatternMass s q A * Real.exp (-(A.card : ℝ))) =
      ∏ i ∈ s, (q i * Real.exp (-1) + (1 - q i)) := by
  rw [Finset.prod_add]
  refine Finset.sum_congr rfl (fun A _ => ?_)
  have he : Real.exp (-(A.card : ℝ)) = ∏ _i ∈ A, Real.exp (-1) := by
    rw [Finset.prod_const, ← Real.exp_nat_mul]
    congr 1
    ring
  rw [he, Finset.prod_mul_distrib, bankPatternMass]
  ring

/-- The discarded-bank bound with the exponent in the paper. -/
theorem bank_discard_mass_le {ι : Type*} [DecidableEq ι] (s : Finset ι) {q : ι → ℝ}
    (hq : ∀ i ∈ s, 1 / 8 ≤ q i ∧ q i ≤ 1) :
    bankDiscardMass s q ≤ Real.exp (-(s.card : ℝ) / 64) := by
  have hq₀ (i : ι) (hi : i ∈ s) : 0 ≤ q i ∧ q i ≤ 1 :=
    ⟨by linarith [(hq i hi).1], (hq i hi).2⟩
  have hExp : Real.exp (-1) ≤ 3 / 8 := by
    rw [Real.exp_neg, ← one_div]
    exact (div_le_iff₀ (Real.exp_pos 1)).mpr (by linarith [Real.exp_one_gt_d9])
  have hpoint (i : ι) (hi : i ∈ s) : q i * Real.exp (-1) + (1 - q i) ≤ Real.exp (-5 / 64) := by
    have h₁ := mul_le_mul_of_nonneg_left hExp (hq₀ i hi).1
    have h₂ := Real.add_one_le_exp (-5 / 64 : ℝ)
    linarith only [h₁, h₂, (hq i hi).1]
  have hm :
      (∑ A ∈ s.powerset, bankPatternMass s q A * Real.exp (-(A.card : ℝ))) ≤
        Real.exp (-5 / 64 * (s.card : ℝ)) := by
    rw [bank_pattern_exp_sum]
    calc
      _ ≤ ∏ _i ∈ s, Real.exp (-5 / 64) := Finset.prod_le_prod (fun i hi =>
        add_nonneg (mul_nonneg (hq₀ i hi).1 (Real.exp_nonneg _)) (sub_nonneg.mpr (hq₀ i hi).2)) hpoint
      _ = _ := by rw [Finset.prod_const, ← Real.exp_nat_mul, mul_comm]
  unfold bankDiscardMass
  calc
    _ ≤ ∑ A ∈ s.powerset, Real.exp ((s.card : ℝ) / 16) *
        (bankPatternMass s q A * Real.exp (-(A.card : ℝ))) := by
      apply Finset.sum_le_sum
      intro A hA
      have hW := bank_pattern_mass_nonneg (Finset.mem_powerset.mp hA) hq₀
      rw [← mul_assoc, mul_comm (Real.exp ((s.card : ℝ) / 16)), mul_assoc, ← Real.exp_add]
      split
      · rename_i hb
        exact le_mul_of_one_le_right hW (Real.one_le_exp (by linarith only [hb]))
      · exact mul_nonneg hW (Real.exp_nonneg _)
    _ ≤ Real.exp ((s.card : ℝ) / 16) * Real.exp (-5 / 64 * (s.card : ℝ)) := by
      rw [← Finset.mul_sum]
      exact mul_le_mul_of_nonneg_left hm (Real.exp_nonneg _)
    _ = _ := by rw [← Real.exp_add]; congr 1; ring

/-- The discarded mass when either of the two independent bank
ranges has too few successful banks. -/
noncomputable def twoRangeBankDiscardMass {ι κ : Type*} [DecidableEq ι] [DecidableEq κ]
    (s : Finset ι) (t : Finset κ) (q : ι → ℝ) (r : κ → ℝ) : ℝ :=
  ∑ A ∈ s.powerset, ∑ B ∈ t.powerset,
    if (A.card : ℝ) < (s.card : ℝ) / 16 ∨ (B.card : ℝ) < (t.card : ℝ) / 16 then
      bankPatternMass s q A * bankPatternMass t r B
    else 0

theorem two_range_bank_discard_mass_le {ι κ : Type*} [DecidableEq ι] [DecidableEq κ]
    (s : Finset ι) (t : Finset κ) {q : ι → ℝ} {r : κ → ℝ}
    (hq : ∀ i ∈ s, 1 / 8 ≤ q i ∧ q i ≤ 1) (hr : ∀ j ∈ t, 1 / 8 ≤ r j ∧ r j ≤ 1) :
    twoRangeBankDiscardMass s t q r ≤
      Real.exp (-(s.card : ℝ) / 64) + Real.exp (-(t.card : ℝ) / 64) := by
  have hq₀ i hi : 0 ≤ q i ∧ q i ≤ 1 := ⟨by linarith [(hq i hi).1], (hq i hi).2⟩
  have hr₀ j hj : 0 ≤ r j ∧ r j ≤ 1 := ⟨by linarith [(hr j hj).1], (hr j hj).2⟩
  calc
    _ ≤ ∑ A ∈ s.powerset, ∑ B ∈ t.powerset,
        ((if (A.card : ℝ) < (s.card : ℝ) / 16 then bankPatternMass s q A else 0) * bankPatternMass t r B +
          bankPatternMass s q A * (if (B.card : ℝ) < (t.card : ℝ) / 16 then bankPatternMass t r B else 0)) := by
      refine Finset.sum_le_sum (fun A hA => Finset.sum_le_sum (fun B hB => ?_))
      have hW := mul_nonneg (bank_pattern_mass_nonneg (Finset.mem_powerset.mp hA) hq₀)
        (bank_pattern_mass_nonneg (Finset.mem_powerset.mp hB) hr₀)
      by_cases hs : (A.card : ℝ) < (s.card : ℝ) / 16 <;>
        by_cases ht : (B.card : ℝ) < (t.card : ℝ) / 16 <;> simp [hs, ht, hW]
    _ = bankDiscardMass s q + bankDiscardMass t r := by
      simp_rw [Finset.sum_add_distrib, ← Finset.mul_sum, bank_pattern_mass_sum, mul_one,
        ← Finset.sum_mul, bank_pattern_mass_sum, one_mul]
      rfl
    _ ≤ _ := add_le_add (bank_discard_mass_le s hq) (bank_discard_mass_le t hr)

/-- Substitution of the number of banks gives the paper's logarithmic
saving for the discarded part of the comparison. -/
theorem two_range_bank_discard_log_bound {ι κ : Type*} [DecidableEq ι] [DecidableEq κ]
    (s : Finset ι) (t : Finset κ) {q : ι → ℝ} {r : κ → ℝ}
    (hq : ∀ i ∈ s, 1 / 8 ≤ q i ∧ q i ≤ 1) (hr : ∀ j ∈ t, 1 / 8 ≤ r j ∧ r j ≤ 1)
    {L φ : ℝ} (hL : 0 < L) (hφ : 0 < φ)
    (hs : Real.log L / (64 * φ) ≤ s.card) (ht : Real.log L / (64 * φ) ≤ t.card) :
    twoRangeBankDiscardMass s t q r ≤ 2 * L ^ (-(1 / (4096 * φ))) := by
  have h (b : ℕ) (hb : Real.log L / (64 * φ) ≤ b) :
      Real.exp (-(b : ℝ) / 64) ≤ L ^ (-(1 / (4096 * φ))) := by
    rw [Real.rpow_def_of_pos hL]
    apply Real.exp_le_exp.mpr
    have he : Real.log L * (-(1 / (4096 * φ))) = -(Real.log L / (64 * φ)) / 64 := by
      field_simp
      ring
    rw [he]
    linarith only [hb]
  exact (two_range_bank_discard_mass_le s t hq hr).trans
    (by rw [two_mul]; exact add_le_add (h s.card hs) (h t.card ht))

section IndependentTrials

variable {ι : Type*} [Fintype ι] [DecidableEq ι] {α : ι → Type*} [∀ i, Fintype (α i)]

/-- The set of successful indices in one joint trial outcome. -/
def successfulTrials (success : ∀ i, α i → Bool) (z : ∀ i, α i) : Finset ι :=
  Finset.univ.filter (fun i => success i (z i) = true)

/-- The sum of the outcome weights for a successful trial. -/
noncomputable def trialSuccessMass (w : ∀ i, α i → ℝ) (success : ∀ i, α i → Bool) (i : ι) : ℝ :=
  ∑ x ∈ Finset.univ.filter (fun x => success i x = true), w i x

/-- Grouping independent finite trial outcomes by their success set
gives exactly the Bernoulli pattern mass. -/
theorem independent_trial_success_pattern (w : ∀ i, α i → ℝ) (success : ∀ i, α i → Bool)
    (hnorm : ∀ i, ∑ z : α i, w i z = 1) (A : Finset ι) :
    (∑ z : (∀ i, α i), if successfulTrials success z = A then ∏ i, w i (z i) else 0) =
      bankPatternMass Finset.univ (trialSuccessMass w success) A := by
  let cells (i : ι) := Finset.univ.filter (fun x : α i => success i x = true ↔ i ∈ A)
  have heq : Finset.univ.filter (fun z : ∀ i, α i => successfulTrials success z = A) =
      Fintype.piFinset cells := by
    ext z
    rw [Fintype.mem_piFinset]
    simp only [Finset.mem_filter, Finset.mem_univ, true_and, cells]
    constructor
    · intro hz i
      rw [← hz]
      simp [successfulTrials]
    · intro hz
      ext i
      simpa only [successfulTrials, Finset.mem_filter, Finset.mem_univ, true_and] using hz i
  have hcell (i : ι) : ∑ x ∈ cells i, w i x =
      if i ∈ A then trialSuccessMass w success i else 1 - trialSuccessMass w success i := by
    by_cases hi : i ∈ A
    · simp [cells, hi, trialSuccessMass]
    · have hs := Finset.sum_filter_add_sum_filter_not (Finset.univ : Finset (α i))
        (fun x => success i x = true) (w i)
      rw [hnorm i] at hs
      have hset : cells i = Finset.univ.filter (fun x : α i => ¬ success i x = true) := by
        ext x
        simp [cells, hi]
      rw [hset, if_neg hi]
      exact eq_sub_of_add_eq' hs
  rw [← Finset.sum_filter, heq, ← Finset.prod_univ_sum]
  simp_rw [hcell]
  rw [← Finset.prod_filter_mul_prod_filter_not Finset.univ (fun i => i ∈ A),
    Finset.filter_mem_eq_inter, Finset.univ_inter, Finset.filter_notMem_eq_sdiff]
  unfold bankPatternMass
  congr 1
  · exact Finset.prod_congr rfl (fun i hi => if_pos hi)
  · exact Finset.prod_congr rfl (fun i hi => if_neg (Finset.mem_sdiff.mp hi).2)

/-- Any event in the bank success set has the same mass under the
independent finite trial law and the Bernoulli law. -/
theorem independent_trial_success_event (w : ∀ i, α i → ℝ) (success : ∀ i, α i → Bool)
    (hnorm : ∀ i, ∑ z : α i, w i z = 1) (event : Finset ι → Prop) [DecidablePred event] :
    (∑ z : (∀ i, α i), if event (successfulTrials success z) then ∏ i, w i (z i) else 0) =
      ∑ A ∈ (Finset.univ : Finset ι).powerset,
        if event A then bankPatternMass Finset.univ (trialSuccessMass w success) A else 0 := by
  rw [← Finset.sum_fiberwise_of_maps_to (g := successfulTrials success)
    (fun _ _ => Finset.mem_powerset.mpr (Finset.subset_univ _))]
  refine Finset.sum_congr rfl (fun A _ => ?_)
  by_cases he : event A
  · rw [if_pos he, ← independent_trial_success_pattern w success hnorm A,
      ← Finset.sum_filter (fun z => successfulTrials success z = A) (fun z => ∏ i, w i (z i))]
    refine Finset.sum_congr rfl (fun z hz => ?_)
    exact if_pos ((Finset.mem_filter.mp hz).2.symm ▸ he)
  · rw [if_neg he]
    refine Finset.sum_eq_zero (fun z hz => ?_)
    exact if_neg ((Finset.mem_filter.mp hz).2.symm ▸ he)

/-- The rare-bank estimate for independent finite samples, before
forgetting their actual records. -/
theorem independent_trial_discard_mass_le (w : ∀ i, α i → ℝ) (success : ∀ i, α i → Bool)
    (hnorm : ∀ i, ∑ z : α i, w i z = 1)
    (hq : ∀ i, 1 / 8 ≤ trialSuccessMass w success i ∧ trialSuccessMass w success i ≤ 1) :
    (∑ z : (∀ i, α i), if ((successfulTrials success z).card : ℝ) < (Fintype.card ι : ℝ) / 16 then
      ∏ i, w i (z i) else 0) ≤ Real.exp (-(Fintype.card ι : ℝ) / 64) := by
  rw [independent_trial_success_event w success hnorm (fun A => (A.card : ℝ) < (Fintype.card ι : ℝ) / 16)]
  have hb := bank_discard_mass_le (Finset.univ : Finset ι) (fun i _ => hq i)
  convert hb using 1 <;> simp only [bankDiscardMass, Finset.card_univ]
  refine Finset.sum_congr rfl (fun A _ => ?_)
  by_cases h : (A.card : ℝ) < (Fintype.card ι : ℝ) / 16 <;> simp [h]

end IndependentTrials

end QuadraticNorms
