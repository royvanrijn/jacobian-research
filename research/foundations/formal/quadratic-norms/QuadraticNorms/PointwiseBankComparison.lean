import QuadraticNorms.PointwiseBankReconstruction
import QuadraticNorms.PointwiseBankConcentration

open scoped BigOperators Classical

namespace QuadraticNorms

/-- Failure of the success threshold in either size range. -/
def twoRangeBankFailure {ι κ X : Type*} (small : Finset ι) (large : Finset κ)
    (smallSuccess : X → ι → Bool) (largeSuccess : X → κ → Bool) (e : X) : Prop :=
  ((small.filter (fun i => smallSuccess e i = true)).card : ℝ) < (small.card : ℝ) / 16 ∨
    ((large.filter (fun j => largeSuccess e j = true)).card : ℝ) < (large.card : ℝ) / 16

/-- The bank reconstruction and discarded-bank estimates combine to
bound the full comparison. The active-set correlation and positive
record envelope are the hypotheses supplied by the analytic argument. -/
theorem two_range_model_comparison_of_bank_estimates {ι κ X : Type*} [DecidableEq ι] [DecidableEq κ]
    (edges : Finset X) (small : Finset ι) (large : Finset κ)
    (smallSuccess : X → ι → Bool) (largeSuccess : X → κ → Bool)
    (f₁ f₂ d₁ d₂ : X → ℝ) (hf₁ : ∀ e ∈ edges, 0 ≤ f₁ e) (hf₂ : ∀ e ∈ edges, 0 ≤ f₂ e)
    (hd₁ : ∀ e ∈ edges, 0 ≤ d₁ e) (hd₂ : ∀ e ∈ edges, 0 ≤ d₂ e)
    {q : ι → ℝ} {r : κ → ℝ} (hq : ∀ i ∈ small, 1 / 8 ≤ q i ∧ q i ≤ 1)
    (hr : ∀ j ∈ large, 1 / 8 ≤ r j ∧ r j ≤ 1)
    {L φ S R ε : ℝ} (hL : 0 < L) (hφ : 0 < φ) (hS : 0 ≤ S) (hε : 0 ≤ ε)
    (hs : Real.log L / (64 * φ) ≤ small.card) (ht : Real.log L / (64 * φ) ≤ large.card)
    (hActive : ∀ s ∈ small.powerset, ∀ t ∈ large.powerset,
      ⌈(small.card : ℝ) / 16⌉₊ ≤ s.card → ⌈(large.card : ℝ) / 16⌉₊ ≤ t.card →
      |∑ e ∈ edges, (f₁ e * f₂ e - d₁ e * d₂ e) *
        bankSelectionWeight smallSuccess s e * bankSelectionWeight largeSuccess t e| ≤ ε)
    (hEnvelope : ∀ a b : X → ℝ, (a = f₁ ∨ a = d₁) → (b = f₂ ∨ b = d₂) →
      (∑ e ∈ edges, if twoRangeBankFailure small large smallSuccess largeSuccess e then a e * b e else 0) ≤
        S * twoRangeBankDiscardMass small large q r + R) :
    |∑ e ∈ edges, (f₁ e * f₂ e - d₁ e * d₂ e)| ≤
      (3 : ℝ) ^ (small.card + large.card) * ε + 4 * S * L ^ (-(1 / (4096 * φ))) + 2 * R := by
  let w (e : X) := f₁ e * f₂ e - d₁ e * d₂ e
  let fail (e : X) := twoRangeBankFailure small large smallSuccess largeSuccess e
  let kept (e : X) := if ¬ fail e then w e else 0
  let lost (e : X) := if fail e then w e else 0
  have hKept : |∑ e ∈ edges, kept e| ≤ (3 : ℝ) ^ (small.card + large.card) * ε := by
    have h := two_range_banked_correlation_bound edges small large
      ⌈(small.card : ℝ) / 16⌉₊ ⌈(large.card : ℝ) / 16⌉₊ smallSuccess largeSuccess w hε hActive
    simpa only [kept, fail, twoRangeBankFailure, not_or, not_lt, Nat.ceil_le] using h
  have hLost : |∑ e ∈ edges, lost e| ≤ 2 * (S * twoRangeBankDiscardMass small large q r + R) := by
    calc
      _ ≤ ∑ e ∈ edges, |lost e| := Finset.abs_sum_le_sum_abs ..
      _ ≤ ∑ e ∈ edges, ((if fail e then f₁ e * f₂ e else 0) + (if fail e then d₁ e * d₂ e else 0)) := by
        refine Finset.sum_le_sum (fun e he => ?_)
        by_cases h : fail e
        · simp only [lost, if_pos h, w]
          apply abs_le.mpr
          have hf : 0 ≤ f₁ e * f₂ e := mul_nonneg (hf₁ e he) (hf₂ e he)
          have hd : 0 ≤ d₁ e * d₂ e := mul_nonneg (hd₁ e he) (hd₂ e he)
          exact ⟨by linarith only [hf], by linarith only [hd]⟩
        · simp only [lost, if_neg h, abs_zero, add_zero, le_refl]
      _ ≤ _ := by
        rw [Finset.sum_add_distrib, two_mul]
        exact add_le_add (hEnvelope f₁ f₂ (Or.inl rfl) (Or.inl rfl)) (hEnvelope d₁ d₂ (Or.inr rfl) (Or.inr rfl))
  have hSplit : (∑ e ∈ edges, w e) = (∑ e ∈ edges, kept e) + ∑ e ∈ edges, lost e := by
    rw [← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun e _ => ?_)
    by_cases h : fail e <;> simp [kept, lost, h]
  have hRare := mul_le_mul_of_nonneg_left (two_range_bank_discard_log_bound small large hq hr hL hφ hs ht) hS
  change |∑ e ∈ edges, w e| ≤ _
  rw [hSplit]
  exact (abs_add_le ..).trans ((add_le_add hKept hLost).trans (by linarith only [hRare]))

end QuadraticNorms
