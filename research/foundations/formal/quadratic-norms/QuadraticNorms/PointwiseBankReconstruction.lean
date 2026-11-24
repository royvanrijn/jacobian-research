import Mathlib.Data.Nat.Choose.Sum
import Mathlib.Algebra.Order.BigOperators.GroupWithZero.Finset
import Mathlib.Algebra.BigOperators.Group.Finset.Powerset
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic.Ring

open scoped BigOperators

namespace QuadraticNorms

private def subsetCoefficient {ι : Type*} (g : Finset ι → ℤ) (s : Finset ι) : ℤ :=
  ∑ t ∈ s.powerset, (-1) ^ (s.card - t.card) * g t

private theorem supersets_sign_sum {ι : Type*} [DecidableEq ι] {s u : Finset ι} (hu : u ⊆ s) :
    (∑ t ∈ s.powerset.filter (fun t => u ⊆ t), ((-1) ^ (t.card - u.card) : ℤ)) =
      if u = s then 1 else 0 := by
  have heq : (∑ t ∈ (s \ u).powerset, ((-1) ^ t.card : ℤ)) =
      ∑ t ∈ s.powerset.filter (fun t => u ⊆ t), ((-1) ^ (t.card - u.card) : ℤ) := by
    refine Finset.sum_bij' (fun t _ => t ∪ u) (fun t _ => t \ u) ?_ ?_ ?_ ?_ ?_
    · intro t ht
      have hsub := Finset.mem_powerset.mp ht
      exact Finset.mem_filter.mpr ⟨Finset.mem_powerset.mpr
        (Finset.union_subset (hsub.trans Finset.sdiff_subset) hu), Finset.subset_union_right⟩
    · intro t ht
      have ⟨hsub, _⟩ := Finset.mem_filter.mp ht
      apply Finset.mem_powerset.mpr
      intro p hp
      obtain ⟨hpt, hpu⟩ := Finset.mem_sdiff.mp hp
      exact Finset.mem_sdiff.mpr ⟨Finset.mem_powerset.mp hsub hpt, hpu⟩
    · intro t ht
      ext a
      have hsub := Finset.mem_powerset.mp ht
      have he (ha : a ∈ t) : a ∉ u := (Finset.mem_sdiff.mp (hsub ha)).2
      simp only [Finset.mem_sdiff, Finset.mem_union]
      exact ⟨fun h => h.1.resolve_right h.2, fun h => ⟨Or.inl h, he h⟩⟩
    · intro t ht
      exact Finset.sdiff_union_of_subset (Finset.mem_filter.mp ht).2
    · intro t ht
      have hd : Disjoint t u := Finset.disjoint_left.mpr
        (fun _ ha => (Finset.mem_sdiff.mp (Finset.mem_powerset.mp ht ha)).2)
      rw [Finset.card_union_of_disjoint hd, Nat.add_sub_cancel]
  rw [← heq, Finset.sum_powerset_neg_one_pow_card]
  have h : s \ u = ∅ ↔ u = s := by
    rw [Finset.sdiff_eq_empty_iff_subset]
    exact ⟨fun h => Finset.Subset.antisymm hu h, fun h => h ▸ Finset.Subset.rfl⟩
  simp only [h]

private theorem sum_subset_coefficients {ι : Type*} [DecidableEq ι]
    (g : Finset ι → ℤ) (s : Finset ι) :
    ∑ t ∈ s.powerset, subsetCoefficient g t = g s := by
  classical
  calc
    _ = ∑ t ∈ s.powerset, ∑ u ∈ s.powerset,
        if u ⊆ t then (-1) ^ (t.card - u.card) * g u else 0 := by
      apply Finset.sum_congr rfl
      intro t ht
      unfold subsetCoefficient
      rw [← Finset.sum_filter]
      congr 1
      ext u
      simp only [Finset.mem_filter, Finset.mem_powerset]
      exact ⟨fun h => ⟨h.trans (Finset.mem_powerset.mp ht), h⟩, And.right⟩
    _ = ∑ u ∈ s.powerset, g u * (∑ t ∈ s.powerset.filter (fun t => u ⊆ t),
          ((-1) ^ (t.card - u.card) : ℤ)) := by
      rw [Finset.sum_comm]
      apply Finset.sum_congr rfl
      intro u _
      rw [Finset.mul_sum, Finset.sum_filter]
      exact Finset.sum_congr rfl (fun t _ => by split <;> simp [mul_comm])
    _ = g s := by
      rw [Finset.sum_eq_single_of_mem s (Finset.mem_powerset.mpr Finset.Subset.rfl)]
      · rw [supersets_sign_sum Finset.Subset.rfl]
        simp
      · intro u hu hne
        rw [supersets_sign_sum (Finset.mem_powerset.mp hu)]
        simp [hne]

/-- The coefficient of a bank set in the expansion of the condition
that at least k banks succeed. -/
def bankThresholdCoefficient {ι : Type*} (k : ℕ) (s : Finset ι) : ℤ :=
  subsetCoefficient (fun t => if k ≤ t.card then 1 else 0) s

private theorem coefficient_complement_sum {ι : Type*} [DecidableEq ι]
    (k : ℕ) (s : Finset ι) :
    bankThresholdCoefficient k s =
      ∑ u ∈ s.powerset, if u.card + k ≤ s.card then ((-1) ^ u.card : ℤ) else 0 := by
  unfold bankThresholdCoefficient subsetCoefficient
  refine Finset.sum_bij' (fun t _ => s \ t) (fun u _ => s \ u) ?_ ?_ ?_ ?_ ?_
  · intro t ht
    exact Finset.mem_powerset.mpr Finset.sdiff_subset
  · intro u hu
    exact Finset.mem_powerset.mpr Finset.sdiff_subset
  · intro t ht
    ext p
    have h (hp : p ∈ t) : p ∈ s := Finset.mem_powerset.mp ht hp
    simp only [Finset.mem_sdiff, not_and, Classical.not_not]
    exact ⟨fun ⟨hps, hc⟩ => hc hps, fun ht' => ⟨h ht', fun _ => ht'⟩⟩
  · intro t ht
    ext p
    have h (hp : p ∈ t) : p ∈ s := Finset.mem_powerset.mp ht hp
    simp only [Finset.mem_sdiff, not_and, Classical.not_not]
    exact ⟨fun ⟨hps, hc⟩ => hc hps, fun ht' => ⟨h ht', fun _ => ht'⟩⟩
  · intro t ht
    rw [Finset.card_sdiff_of_subset (Finset.mem_powerset.mp ht)]
    have hle := Finset.card_le_card (Finset.mem_powerset.mp ht)
    have hiff : s.card - t.card + k ≤ s.card ↔ k ≤ t.card := by omega
    simp only [hiff]
    split <;> simp

/-- The signed binomial coefficient in the reconstruction. -/
theorem bank_threshold_coefficient_formula {ι : Type*} [DecidableEq ι]
    (s : Finset ι) {k : ℕ} (hk : 0 < k) :
    bankThresholdCoefficient k s =
      if k ≤ s.card then (-1) ^ (s.card - k) * ((s.card - 1).choose (k - 1) : ℤ) else 0 := by
  rw [coefficient_complement_sum, Finset.sum_powerset_apply_card
    (fun j => if j + k ≤ s.card then ((-1) ^ j : ℤ) else 0)]
  by_cases hks : k ≤ s.card
  · rw [if_pos hks]
    let f (j : ℕ) : ℤ := s.card.choose j • (if j + k ≤ s.card then (-1) ^ j else 0)
    change (∑ j ∈ Finset.range (s.card + 1), f j) = _
    rw [← Finset.sum_subset (s₁ := Finset.range (s.card - k + 1))
      (Finset.range_mono (by omega : s.card - k + 1 ≤ s.card + 1))
      (f := f) (fun j hj hcut => by
        have hcut' : ¬ j < s.card - k + 1 := fun h => hcut (Finset.mem_range.mpr h)
        simp [f, show ¬ j + k ≤ s.card by omega])]
    have hterm (j : ℕ) (hj : j ∈ Finset.range (s.card - k + 1)) :
        f j = (-1) ^ j * (s.card.choose j : ℤ) := by
      have hj' := Finset.mem_range.mp hj
      simp [f, show j + k ≤ s.card by omega, mul_comm]
    rw [Finset.sum_congr rfl hterm]
    have hcard : s.card - 1 + 1 = s.card := by omega
    have hAlt := Int.alternating_sum_range_choose_eq_choose (n := s.card - 1) (m := s.card - k)
    rw [hcard, Nat.choose_symm_of_eq_add (by omega : s.card - 1 = (s.card - k) + (k - 1))] at hAlt
    exact hAlt
  · simp [hks, show ∀ j : ℕ, ¬ j + k ≤ s.card by omega]

/-- The total absolute coefficient mass for one bank range. -/
theorem bank_threshold_coefficient_mass {ι : Type*} (s : Finset ι) (k : ℕ) :
    (∑ t ∈ s.powerset, |bankThresholdCoefficient k t|) ≤ (3 : ℤ) ^ s.card := by
  have hpoint (t : Finset ι) : |bankThresholdCoefficient k t| ≤ (2 : ℤ) ^ t.card := by
    unfold bankThresholdCoefficient subsetCoefficient
    refine (Finset.abs_sum_le_sum_abs _ _).trans ?_
    calc
      _ ≤ ∑ _ ∈ t.powerset, (1 : ℤ) := by
        apply Finset.sum_le_sum
        intro u _
        dsimp only
        split <;> simp [abs_pow]
      _ = _ := by rw [Finset.sum_const, nsmul_eq_mul, mul_one, Finset.card_powerset, Nat.cast_pow, Nat.cast_ofNat]
  refine (Finset.sum_le_sum (fun t _ => hpoint t)).trans_eq ?_
  simpa only [Finset.prod_const, show (1 : ℤ) + 2 = 3 by decide] using
    (Finset.prod_one_add (f := fun _ : ι => (2 : ℤ)) s).symm

/-- The threshold indicator is the exact polynomial in the bank
success indicators. -/
theorem bank_threshold_polynomial {ι : Type*} [DecidableEq ι] (s : Finset ι) (k : ℕ)
    (success : ι → Prop) [DecidablePred success] :
    (∑ t ∈ s.powerset, bankThresholdCoefficient k t * ∏ i ∈ t, (if success i then 1 else 0 : ℤ)) =
      if k ≤ (s.filter success).card then 1 else 0 := by
  have hp (t : Finset ι) (ht : t ∈ s.powerset) :
      (∏ i ∈ t, (if success i then 1 else 0 : ℤ)) = if t ⊆ s.filter success then 1 else 0 := by
    by_cases h : t ⊆ s.filter success
    · rw [if_pos h]
      exact Finset.prod_eq_one (fun i hi => if_pos (Finset.mem_filter.mp (h hi)).2)
    · rw [if_neg h]
      have hfail : ∃ i ∈ t, ¬ success i := by
        by_contra! hgood
        exact h (fun i hi => Finset.mem_filter.mpr ⟨Finset.mem_powerset.mp ht hi, hgood i hi⟩)
      obtain ⟨i, hi, hno⟩ := hfail
      exact Finset.prod_eq_zero hi (if_neg hno)
  calc
    _ = ∑ t ∈ s.powerset, if t ⊆ s.filter success then bankThresholdCoefficient k t else 0 := by
      apply Finset.sum_congr rfl
      intro t ht
      rw [hp _ ht]
      split <;> simp
    _ = ∑ t ∈ (s.filter success).powerset, bankThresholdCoefficient k t := by
      rw [← Finset.sum_filter]
      congr 1
      ext t
      simp only [Finset.mem_filter, Finset.mem_powerset]
      exact ⟨And.right, fun h => ⟨h.trans (Finset.filter_subset ..), h⟩⟩
    _ = _ := sum_subset_coefficients _ _


private theorem coefficient_zero_below {ι : Type*} {s : Finset ι} {k : ℕ} (hs : s.card < k) :
    bankThresholdCoefficient k s = 0 := by
  unfold bankThresholdCoefficient subsetCoefficient
  refine Finset.sum_eq_zero (fun t ht => ?_)
  have hlt := (Finset.card_le_card (Finset.mem_powerset.mp ht)).trans_lt hs
  simp [Nat.not_le_of_lt hlt]

/-- A product of success indicators for one active set of banks. -/
noncomputable def bankSelectionWeight {ι X : Type*}
    (success : X → ι → Bool) (s : Finset ι) (e : X) : ℝ :=
  ∏ i ∈ s, if success e i then 1 else 0

private theorem threshold_weight_polynomial {ι X : Type*} [DecidableEq ι]
    (edges : Finset X) (banks : Finset ι) (k : ℕ) (success : X → ι → Bool) (w : X → ℝ) :
    (∑ e ∈ edges, if k ≤ (banks.filter (fun i => success e i = true)).card then w e else 0) =
      ∑ t ∈ banks.powerset, (bankThresholdCoefficient k t : ℝ) *
        ∑ e ∈ edges, w e * bankSelectionWeight success t e := by
  have h (e : X) :
      (∑ t ∈ banks.powerset, (bankThresholdCoefficient k t : ℝ) * bankSelectionWeight success t e) =
        if k ≤ (banks.filter (fun i => success e i = true)).card then 1 else 0 := by
    have heq := congrArg (fun a : ℤ => (a : ℝ))
      (bank_threshold_polynomial banks k (fun i => success e i = true))
    push_cast at heq
    simpa [bankSelectionWeight] using heq
  calc
    _ = ∑ e ∈ edges, w e *
          (∑ t ∈ banks.powerset, (bankThresholdCoefficient k t : ℝ) * bankSelectionWeight success t e) := by
      apply Finset.sum_congr rfl
      intro e _
      rw [h e]
      split <;> simp
    _ = _ := by
      simp_rw [Finset.mul_sum]
      rw [Finset.sum_comm]
      apply Finset.sum_congr rfl
      intro t _
      apply Finset.sum_congr rfl
      intro e _
      ring

/-- Reconstruction of a weighted correlation from bounds on its
active-bank correlations in one size range. -/
theorem banked_correlation_bound {ι X : Type*} [DecidableEq ι] (edges : Finset X)
    (banks : Finset ι) (k : ℕ) (success : X → ι → Bool) (w : X → ℝ)
    {ε : ℝ} (hε : 0 ≤ ε)
    (hActive : ∀ t ∈ banks.powerset, k ≤ t.card →
      |∑ e ∈ edges, w e * bankSelectionWeight success t e| ≤ ε) :
    |∑ e ∈ edges, if k ≤ (banks.filter (fun i => success e i = true)).card then w e else 0| ≤
      (3 : ℝ) ^ banks.card * ε := by
  rw [threshold_weight_polynomial]
  apply (Finset.abs_sum_le_sum_abs _ _).trans
  calc
    _ ≤ ∑ t ∈ banks.powerset, |(bankThresholdCoefficient k t : ℝ)| * ε := by
      apply Finset.sum_le_sum
      intro t ht
      by_cases h : k ≤ t.card
      · rw [abs_mul]
        exact mul_le_mul_of_nonneg_left (hActive t ht h) (abs_nonneg _)
      · simp [coefficient_zero_below (Nat.lt_of_not_ge h)]
    _ = (∑ t ∈ banks.powerset, |(bankThresholdCoefficient k t : ℝ)|) * ε := (Finset.sum_mul ..).symm
    _ ≤ _ := by
      apply mul_le_mul_of_nonneg_right _ hε
      have hcast := (Int.cast_le (R := ℝ)).mpr (bank_threshold_coefficient_mass banks k)
      push_cast at hcast
      exact hcast

/-- Reconstruct the part of a comparison with enough successful small
and large banks. The two coefficient-mass costs multiply. -/
theorem two_range_banked_correlation_bound {ι κ X : Type*} [DecidableEq ι] [DecidableEq κ]
    (edges : Finset X) (small : Finset ι) (large : Finset κ) (k l : ℕ)
    (smallSuccess : X → ι → Bool) (largeSuccess : X → κ → Bool) (w : X → ℝ)
    {ε : ℝ} (hε : 0 ≤ ε)
    (hActive : ∀ s ∈ small.powerset, ∀ t ∈ large.powerset, k ≤ s.card → l ≤ t.card →
      |∑ e ∈ edges, w e * bankSelectionWeight smallSuccess s e * bankSelectionWeight largeSuccess t e| ≤ ε) :
    |∑ e ∈ edges,
      if k ≤ (small.filter (fun i => smallSuccess e i = true)).card ∧
        l ≤ (large.filter (fun i => largeSuccess e i = true)).card then w e else 0| ≤
      (3 : ℝ) ^ (small.card + large.card) * ε := by
  let v (e : X) := if k ≤ (small.filter (fun i => smallSuccess e i = true)).card then w e else 0
  have hActiveLarge (t : Finset κ) (ht : t ∈ large.powerset) (htl : l ≤ t.card) :
      |∑ e ∈ edges, v e * bankSelectionWeight largeSuccess t e| ≤ (3 : ℝ) ^ small.card * ε := by
    have heq (e : X) : v e * bankSelectionWeight largeSuccess t e =
        if k ≤ (small.filter (fun i => smallSuccess e i = true)).card then
          w e * bankSelectionWeight largeSuccess t e else 0 := by
      dsimp [v]
      split <;> simp
    simp_rw [heq]
    refine banked_correlation_bound edges small k smallSuccess
      (fun e => w e * bankSelectionWeight largeSuccess t e) hε ?_
    intro s hs hsk
    simpa only [mul_right_comm] using hActive s hs t ht hsk htl
  have h := banked_correlation_bound edges large l largeSuccess v (by positivity) hActiveLarge
  have heq (e : X) :
      (if l ≤ (large.filter (fun i => largeSuccess e i = true)).card then v e else 0) =
      if k ≤ (small.filter (fun i => smallSuccess e i = true)).card ∧
        l ≤ (large.filter (fun i => largeSuccess e i = true)).card then w e else 0 := by
    by_cases hs : k ≤ (small.filter (fun i => smallSuccess e i = true)).card <;>
      by_cases ht : l ≤ (large.filter (fun i => largeSuccess e i = true)).card <;> simp [v, hs, ht]
  simp_rw [heq] at h
  simpa only [pow_add, mul_assoc, mul_left_comm, mul_comm] using h

end QuadraticNorms
