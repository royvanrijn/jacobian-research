import QuadraticNorms.PointwiseBankValuations
import Mathlib.Analysis.Complex.ExponentialBounds

open scoped BigOperators Topology Classical
open Filter

namespace QuadraticNorms

/-- A successful bank has exactly one left hit and one right hit. -/
def PrimeBankSuccess {ι : Type*} (z : ι → BankPrimeOutcome) : Prop :=
  (∃! i, z i = .left) ∧ (∃! i, z i = .right)

section Success

variable {ι : Type*} [Fintype ι] [DecidableEq ι]

/-- The mass of successful records in one bank with reciprocal weights r. -/
noncomputable def primeBankSuccessMass (r : ι → ℝ) : ℝ :=
  ∑ z : ι → BankPrimeOutcome, if PrimeBankSuccess z then ∏ i, bankPrimeOutcomeWeight (r i) (z i) else 0

private theorem bank_outcome_weights_sum (r : ℝ) :
    ∑ z : BankPrimeOutcome, bankPrimeOutcomeWeight r z = 1 := by
  have heq : (Finset.univ : Finset BankPrimeOutcome) = {.neither, .left, .right} := by
    ext z
    cases z <;> simp
  rw [heq]
  simp [Finset.sum_insert, bankPrimeOutcomeWeight]
  ring

theorem prime_bank_success_mass_mem_Icc (r : ι → ℝ) (hr : ∀ i, 0 ≤ r i ∧ r i ≤ 1 / 2) :
    primeBankSuccessMass r ∈ Set.Icc 0 1 := by
  classical
  have hw (i : ι) (z : BankPrimeOutcome) : 0 ≤ bankPrimeOutcomeWeight (r i) z := by
    cases z <;> simp [bankPrimeOutcomeWeight] <;> linarith [(hr i).1, (hr i).2]
  have ht : (∑ z : ι → BankPrimeOutcome, ∏ i, bankPrimeOutcomeWeight (r i) (z i)) = 1 := by
    simpa only [Fintype.piFinset_univ, bank_outcome_weights_sum, Finset.prod_const_one] using
      Finset.sum_prod_piFinset Finset.univ (fun i z => bankPrimeOutcomeWeight (r i) z)
  refine ⟨Finset.sum_nonneg (fun z _ => ?_), ?_⟩
  · split
    · exact Finset.prod_nonneg (fun i _ => hw i _)
    · exact le_rfl
  · rw [← ht]
    exact Finset.sum_le_sum (fun z _ => by split; exact le_rfl; exact Finset.prod_nonneg (fun i _ => hw i _))

private def bankSuccessRecord (i j : ι) (k : ι) : BankPrimeOutcome :=
  if k = i then .left else if k = j then .right else .neither

omit [Fintype ι] in
private theorem bankSuccessRecord_left (i j k : ι) :
    bankSuccessRecord i j k = .left ↔ k = i := by
  unfold bankSuccessRecord
  by_cases h : k = i <;> by_cases h' : k = j <;> simp [h, h']

omit [Fintype ι] in
private theorem bankSuccessRecord_right {i j : ι} (hij : i ≠ j) (k : ι) :
    bankSuccessRecord i j k = .right ↔ k = j := by
  unfold bankSuccessRecord
  by_cases h : k = i
  · subst k
    simp [hij]
  · simp [h]

omit [Fintype ι] in
private theorem bank_success_iff_exists_pair (z : ι → BankPrimeOutcome) :
    PrimeBankSuccess z ↔ ∃ i j : ι, i ≠ j ∧ z = bankSuccessRecord i j := by
  constructor
  · rintro ⟨⟨i, hi, hiu⟩, ⟨j, hj, hju⟩⟩
    have hij : i ≠ j := by rintro rfl; rw [hi] at hj; contradiction
    refine ⟨i, j, hij, funext (fun k => ?_)⟩
    by_cases h : k = i
    · subst k
      simpa [bankSuccessRecord] using hi
    · by_cases h' : k = j
      · subst k
        simpa [bankSuccessRecord, h] using hj
      · have hz : z k = .neither := by
          cases hz' : z k
          · rfl
          · exact (h (hiu k hz')).elim
          · exact (h' (hju k hz')).elim
        simpa [bankSuccessRecord, h, h'] using hz
  · rintro ⟨i, j, hij, rfl⟩
    exact ⟨⟨i, (bankSuccessRecord_left i j i).mpr rfl, fun k h => (bankSuccessRecord_left i j k).mp h⟩,
      ⟨j, (bankSuccessRecord_right hij j).mpr rfl, fun k h => (bankSuccessRecord_right hij k).mp h⟩⟩

private theorem bank_success_record_weight (r : ι → ℝ)
    (hr : ∀ i, 1 - 2 * r i ≠ 0) {i j : ι} (hij : i ≠ j) :
    (∏ k, bankPrimeOutcomeWeight (r k) (bankSuccessRecord i j k)) =
      (∏ k, (1 - 2 * r k)) * ((r i / (1 - 2 * r i)) * (r j / (1 - 2 * r j))) := by
  let w := fun k => r k / (1 - 2 * r k)
  have h (k : ι) : bankPrimeOutcomeWeight (r k) (bankSuccessRecord i j k) =
      (1 - 2 * r k) * ((if k = i then w k else 1) * (if k = j then w k else 1)) := by
    by_cases hi : k = i
    · subst k
      simp only [bankSuccessRecord, if_true, hij, if_false, mul_one, bankPrimeOutcomeWeight, w]
      rw [← mul_div_assoc, mul_div_cancel_left₀ _ (hr i)]
    · by_cases hj : k = j
      · subst k
        simp only [bankSuccessRecord, hi, if_false, if_true, one_mul, bankPrimeOutcomeWeight, w]
        rw [← mul_div_assoc, mul_div_cancel_left₀ _ (hr j)]
      · simp [bankSuccessRecord, hi, hj, bankPrimeOutcomeWeight]
  simp_rw [h, Finset.prod_mul_distrib, Finset.prod_ite_eq']
  simp [w]

private theorem bank_distinct_pair_sum (w : ι → ℝ) :
    (∑ ij ∈ Finset.univ.filter (fun ij : ι × ι => ij.1 ≠ ij.2), w ij.1 * w ij.2) =
      (∑ i, w i) ^ 2 - ∑ i, w i ^ 2 := by
  have hsplit := Finset.sum_filter_add_sum_filter_not (Finset.univ : Finset (ι × ι))
    (fun ij => ij.1 ≠ ij.2) (fun ij => w ij.1 * w ij.2)
  have hdiag : (∑ ij ∈ Finset.univ.filter (fun ij : ι × ι => ¬ ij.1 ≠ ij.2), w ij.1 * w ij.2) =
      ∑ i, w i ^ 2 := by
    symm
    refine Finset.sum_bij (fun i _ => (i, i)) ?_ ?_ ?_ ?_
    · intro i _
      simp
    · intro a _ b _ h
      exact congrArg Prod.fst h
    · rintro ⟨i, j⟩ h
      have he := not_not.mp (Finset.mem_filter.mp h).2
      exact ⟨i, Finset.mem_univ _, Prod.ext rfl he⟩
    · intro i _
      exact pow_two _
  have hfull : (∑ ij : ι × ι, w ij.1 * w ij.2) = (∑ i, w i) ^ 2 := by
    rw [Fintype.sum_prod_type (fun ij : ι × ι => w ij.1 * w ij.2), pow_two, Finset.sum_mul]
    simp_rw [Finset.mul_sum]
  rw [hdiag, hfull] at hsplit
  linarith only [hsplit]

/-- The exact probability of one hit on each side in a prime bank.
The square minus the diagonal is the ordered sum over distinct primes. -/
theorem prime_bank_success_mass_exact (r : ι → ℝ)
    (hr : ∀ i, 1 - 2 * r i ≠ 0) :
    primeBankSuccessMass r =
      (∏ i, (1 - 2 * r i)) *
        ((∑ i, r i / (1 - 2 * r i)) ^ 2 - ∑ i, (r i / (1 - 2 * r i)) ^ 2) := by
  let pairs := Finset.univ.filter (fun ij : ι × ι => ij.1 ≠ ij.2)
  let P := ∏ i, (1 - 2 * r i)
  let w := fun i => r i / (1 - 2 * r i)
  unfold primeBankSuccessMass
  rw [← Finset.sum_filter]
  calc
    _ = ∑ ij ∈ pairs, P * (w ij.1 * w ij.2) := by
      symm
      refine Finset.sum_bij (fun ij _ => bankSuccessRecord ij.1 ij.2) ?_ ?_ ?_ ?_
      · intro ij hij
        exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, (bank_success_iff_exists_pair _).mpr
          ⟨ij.1, ij.2, (Finset.mem_filter.mp hij).2, rfl⟩⟩
      · intro a ha b hb h
        have hl : bankSuccessRecord b.1 b.2 a.1 = .left := by
          rw [← h]
          exact (bankSuccessRecord_left a.1 a.2 a.1).mpr rfl
        have hr' : bankSuccessRecord b.1 b.2 a.2 = .right := by
          rw [← h]
          exact (bankSuccessRecord_right (Finset.mem_filter.mp ha).2 a.2).mpr rfl
        exact Prod.ext ((bankSuccessRecord_left b.1 b.2 a.1).mp hl)
          ((bankSuccessRecord_right (Finset.mem_filter.mp hb).2 a.2).mp hr')
      · intro z hz
        obtain ⟨i, j, hij, rfl⟩ := (bank_success_iff_exists_pair z).mp (Finset.mem_filter.mp hz).2
        exact ⟨(i, j), Finset.mem_filter.mpr ⟨Finset.mem_univ _, hij⟩, rfl⟩
      · intro ij hij
        exact (bank_success_record_weight r hr (Finset.mem_filter.mp hij).2).symm
    _ = _ := by rw [← Finset.mul_sum, bank_distinct_pair_sum]

private theorem prod_sub_prod_le_sum_sub {α : Type*} (s : Finset α) (a b : α → ℝ)
    (ha : ∀ i ∈ s, 0 ≤ a i) (hab : ∀ i ∈ s, a i ≤ b i) (hb : ∀ i ∈ s, b i ≤ 1) :
    (∏ i ∈ s, b i) - (∏ i ∈ s, a i) ≤ ∑ i ∈ s, (b i - a i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert i s hi ih =>
    have ha' j (hj : j ∈ s) := ha j (Finset.mem_insert_of_mem hj)
    have hab' j (hj : j ∈ s) := hab j (Finset.mem_insert_of_mem hj)
    have hb' j (hj : j ∈ s) := hb j (Finset.mem_insert_of_mem hj)
    have hB : ∏ j ∈ s, b j ≤ 1 := by
      simpa only [Finset.prod_const_one] using Finset.prod_le_prod
        (fun j hj => (ha' j hj).trans (hab' j hj)) hb'
    have hprod := sub_nonneg.mpr (Finset.prod_le_prod ha' hab')
    have hai := ha i (Finset.mem_insert_self ..)
    have hbi := (hab i (Finset.mem_insert_self ..)).trans (hb i (Finset.mem_insert_self ..))
    have hdi := sub_nonneg.mpr (hab i (Finset.mem_insert_self ..))
    rw [Finset.prod_insert hi, Finset.prod_insert hi, Finset.sum_insert hi]
    calc
      _ = (b i - a i) * (∏ j ∈ s, b j) + a i * ((∏ j ∈ s, b j) - ∏ j ∈ s, a j) := by ring
      _ ≤ (b i - a i) * 1 + 1 * ((∏ j ∈ s, b j) - ∏ j ∈ s, a j) :=
        add_le_add (mul_le_mul_of_nonneg_left hB hdi) (mul_le_mul_of_nonneg_right hbi hprod)
      _ ≤ _ := by simpa only [mul_one, one_mul] using add_le_add (le_refl (b i - a i)) (ih ha' hab' hb')

omit [DecidableEq ι] in
private theorem bank_no_hit_product_approx (r : ι → ℝ) (hr : ∀ i, 0 ≤ r i ∧ r i ≤ 1 / 4) :
    0 ≤ Real.exp (-2 * ∑ i, r i) - ∏ i, (1 - 2 * r i) ∧
      Real.exp (-2 * ∑ i, r i) - ∏ i, (1 - 2 * r i) ≤ 3 * ∑ i, r i ^ 2 := by
  have ha (i : ι) : 0 ≤ 1 - 2 * r i := by linarith [(hr i).2]
  have hab (i : ι) : 1 - 2 * r i ≤ Real.exp (-2 * r i) := by
    linarith only [Real.add_one_le_exp (-2 * r i)]
  have hb (i : ι) : Real.exp (-2 * r i) ≤ 1 := Real.exp_le_one_iff.mpr (by linarith [(hr i).1])
  have herr (i : ι) : Real.exp (-2 * r i) - (1 - 2 * r i) ≤ 3 * r i ^ 2 := by
    have harg : |-2 * r i| ≤ 1 := by
      rw [abs_of_nonpos (by linarith [(hr i).1])]
      linarith [(hr i).2]
    have h := (abs_le.mp (Real.exp_bound harg (n := 2) (by decide))).2
    simp only [Finset.sum_range_succ, Finset.sum_range_zero, zero_add, pow_zero, pow_one,
      Nat.factorial_zero, Nat.factorial_one, Nat.cast_one, div_one] at h
    rw [sq_abs] at h
    norm_num at h
    rw [← neg_mul] at h
    nlinarith only [h]
  rw [Finset.mul_sum, Real.exp_sum]
  refine ⟨sub_nonneg.mpr (Finset.prod_le_prod (fun i _ => ha i) (fun i _ => hab i)), ?_⟩
  refine (prod_sub_prod_le_sum_sub Finset.univ _ _ (fun i _ => ha i)
    (fun i _ => hab i) (fun i _ => hb i)).trans ?_
  rw [Finset.mul_sum]
  exact Finset.sum_le_sum (fun i _ => herr i)

/-- Uniform error in the approximation to the success probability.
The constant forty is sufficient when the total reciprocal mass is at
most two and each reciprocal weight is at most one quarter. -/
theorem prime_bank_success_mass_approx (r : ι → ℝ)
    (hr : ∀ i, 0 ≤ r i ∧ r i ≤ 1 / 4) (hV : ∑ i, r i ≤ 2) :
    |primeBankSuccessMass r - (∑ i, r i) ^ 2 * Real.exp (-2 * ∑ i, r i)| ≤ 40 * ∑ i, r i ^ 2 := by
  let V := ∑ i, r i
  let T := ∑ i, r i ^ 2
  let P := ∏ i, (1 - 2 * r i)
  let w := fun i => r i / (1 - 2 * r i)
  let A := ∑ i, w i
  let B := ∑ i, w i ^ 2
  have hden (i : ι) : 0 < 1 - 2 * r i := by linarith [(hr i).2]
  have hw (i : ι) : r i ≤ w i ∧ w i ≤ 2 * r i ∧ w i - r i ≤ 4 * r i ^ 2 := by
    have h := hr i
    dsimp only [w]
    refine ⟨(le_div_iff₀ (hden i)).mpr (by nlinarith only [sq_nonneg (r i)]), ?_,
      (sub_le_iff_le_add).mpr ((div_le_iff₀ (hden i)).mpr ?_)⟩
    · exact (div_le_iff₀ (hden i)).mpr (by
        have hh := mul_nonneg h.1 (show 0 ≤ 1 - 4 * r i by linarith [h.2])
        nlinarith only [hh])
    · have hh := mul_nonneg (sq_nonneg (r i)) (show 0 ≤ 2 - 8 * r i by linarith [h.2])
      nlinarith only [hh]
  have hV0 : 0 ≤ V := Finset.sum_nonneg (fun i _ => (hr i).1)
  have hAV : V ≤ A := Finset.sum_le_sum (fun i _ => (hw i).1)
  have hA : A ≤ 2 * V := by
    unfold A V
    rw [Finset.mul_sum]
    exact Finset.sum_le_sum (fun i _ => (hw i).2.1)
  have hAerr : A - V ≤ 4 * T := by
    unfold A V T
    rw [← Finset.sum_sub_distrib, Finset.mul_sum]
    exact Finset.sum_le_sum (fun i _ => (hw i).2.2)
  have hB0 : 0 ≤ B := Finset.sum_nonneg (fun _ _ => sq_nonneg _)
  have hB : B ≤ 4 * T := by
    unfold B T
    rw [Finset.mul_sum]
    refine Finset.sum_le_sum (fun i _ => ?_)
    have h := pow_le_pow_left₀ ((hr i).1.trans (hw i).1) (hw i).2.1 2
    nlinarith only [h]
  have hA2 : 0 ≤ A ^ 2 - V ^ 2 := by
    linarith only [pow_le_pow_left₀ hV0 hAV 2]
  have hA2' : A ^ 2 - V ^ 2 ≤ 24 * T := by
    have hsum : A + V ≤ 6 := by linarith only [hA, hV]
    have h := mul_le_mul hAerr hsum (by linarith only [hAV, hV0]) (by
      unfold T
      positivity : 0 ≤ 4 * T)
    nlinarith only [h]
  have hT : 0 ≤ T := Finset.sum_nonneg (fun _ _ => sq_nonneg _)
  have hdiff : |(A ^ 2 - B) - V ^ 2| ≤ 28 * T := by
    apply abs_le.mpr
    constructor <;> linarith only [hA2, hA2', hB, hB0, hT]
  have hP0 : 0 ≤ P := Finset.prod_nonneg (fun i _ => (hden i).le)
  have hP1 : P ≤ 1 := by
    simpa only [Finset.prod_const_one] using Finset.prod_le_prod (s := Finset.univ) (g := fun _ : ι => (1 : ℝ))
      (fun i _ => (hden i).le) (fun i _ => by linarith [(hr i).1])
  have ⟨hPE, hPE'⟩ : 0 ≤ Real.exp (-2 * V) - P ∧ Real.exp (-2 * V) - P ≤ 3 * T :=
    bank_no_hit_product_approx r hr
  have hPEabs : |P - Real.exp (-2 * V)| ≤ 3 * T := by
    rw [abs_sub_comm, abs_of_nonneg (show 0 ≤ Real.exp (-2 * V) - P from hPE)]
    exact hPE'
  have hV2 : V ^ 2 ≤ 4 := by
    have h := pow_le_pow_left₀ hV0 hV 2
    norm_num at h
    exact h
  rw [prime_bank_success_mass_exact r (fun i => (hden i).ne')]
  change |P * (A ^ 2 - B) - V ^ 2 * Real.exp (-2 * V)| ≤ 40 * T
  calc
    _ = |P * ((A ^ 2 - B) - V ^ 2) + V ^ 2 * (P - Real.exp (-2 * V))| := by
      congr 1
      ring
    _ ≤ |P| * |(A ^ 2 - B) - V ^ 2| + |V ^ 2| * |P - Real.exp (-2 * V)| := by
      rw [← abs_mul, ← abs_mul]
      exact abs_add_le ..
    _ ≤ 1 * (28 * T) + 4 * (3 * T) := by
      rw [abs_of_nonneg hP0, abs_of_nonneg (sq_nonneg V)]
      exact add_le_add (mul_le_mul hP1 hdiff (abs_nonneg _) (by norm_num))
        (mul_le_mul hV2 hPEabs (abs_nonneg _) (by norm_num))
    _ = _ := by ring

/-- The success mass is the proportion of actual residues whose CRT
record has exactly one hit on each side in the bank. -/
theorem prime_bank_crt_success_proportion (p k : ι → ℕ)
    (hp : ∀ i, (p i).Prime) (hDistinct : Function.Injective p) (hk : ∀ i, 0 < k i)
    (u : ∀ i, (ZMod (p i ^ k i))ˣ) (d : ∀ i, ZMod (p i ^ k i))
    (hd : ∀ i, ¬ p i ∣ (d i).val) :
    letI : NeZero (∏ i, p i ^ k i) := ⟨Finset.prod_ne_zero_iff.mpr
      (fun i _ => pow_ne_zero _ (hp i).ne_zero)⟩
    ((Finset.univ.filter (fun x : ZMod (∏ i, p i ^ k i) =>
      PrimeBankSuccess (fun i => bankPrimeOutcome (p i) (u i) (d i)
        (bankCRTCoordinate (fun j => p j ^ k j) x i)))).card : ℝ) / (∏ i, (p i : ℝ) ^ k i) =
      primeBankSuccessMass (fun i => 1 / (p i : ℝ)) :=
  bank_joint_outcome_event_proportion p k hp hDistinct hk u d hd PrimeBankSuccess

end Success

/-- A family of banks with reciprocal mass tending to one and
square-reciprocal mass tending to zero eventually has success
probability above one eighth. -/
theorem prime_bank_success_eventually_gt_eighth {α : Type*} {l : Filter α}
    {q V T : α → ℝ} (hV : Tendsto V l (𝓝 1)) (hT : Tendsto T l (𝓝 0))
    (hq : ∀ᶠ x in l, |q x - V x ^ 2 * Real.exp (-2 * V x)| ≤ 40 * T x) :
    ∀ᶠ x in l, 1 / 8 < q x := by
  have hExp : Real.exp (-2) > (1 / 8 : ℝ) := by
    have he : Real.exp 2 = Real.exp 1 ^ 2 := by
      simpa only [Nat.cast_ofNat, mul_one] using Real.exp_nat_mul (1 : ℝ) 2
    rw [Real.exp_neg, he, ← one_div]
    have h := pow_lt_pow_left₀
      (Real.exp_one_lt_d9.trans (by norm_num : (2.7182818286 : ℝ) < 14 / 5)) (by positivity : 0 ≤ Real.exp 1) (by decide : 2 ≠ 0)
    have h' : (Real.exp 1) ^ 2 < 8 := h.trans (by norm_num)
    exact one_div_lt_one_div_of_lt (by positivity) h'
  have hF : Tendsto (fun x => V x ^ 2 * Real.exp (-2 * V x) - 40 * T x) l (𝓝 (Real.exp (-2))) := by
    have hNeg : Tendsto (fun x => -2 * V x) l (𝓝 (-2)) := by
      simpa only [mul_one] using
        (show Tendsto (fun _ : α => (-2 : ℝ)) l (𝓝 (-2)) from tendsto_const_nhds).mul hV
    have hE : Tendsto (fun x => Real.exp (-2 * V x)) l (𝓝 (Real.exp (-2))) :=
      (Real.continuous_exp.tendsto _).comp hNeg
    have hErr : Tendsto (fun x => 40 * T x) l (𝓝 0) := by
      simpa only [mul_zero] using
        (show Tendsto (fun _ : α => (40 : ℝ)) l (𝓝 40) from tendsto_const_nhds).mul hT
    simpa using ((hV.pow 2).mul hE).sub hErr
  filter_upwards [hF.eventually (eventually_gt_nhds hExp), hq] with x hx hqx
  have hh := (abs_le.mp hqx).1
  linarith only [hx, hh]

end QuadraticNorms
