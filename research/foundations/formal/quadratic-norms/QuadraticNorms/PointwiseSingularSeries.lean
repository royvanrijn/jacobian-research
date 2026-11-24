import QuadraticNorms.PointwiseCounts
import Mathlib.Analysis.PSeries
import Mathlib.Analysis.SpecialFunctions.Log.Summable
import Mathlib.Data.Nat.Factorization.Basic

open scoped BigOperators Topology

namespace QuadraticNorms

/-- The factor at a prime forbidden in both coordinates before
accounting for the merger of the excluded roots. -/
noncomputable def commonForbiddenFactor (p : ℕ) : ℝ :=
  1 - 1 / ((p : ℝ) - 1) ^ 2

/-- The correction when the two excluded roots coincide. -/
noncomputable def rootMergerFactor (p : ℕ) : ℝ :=
  ((p : ℝ) - 1) / ((p : ℝ) - 2)

private theorem commonForbiddenFactor_pos {p : ℕ} (hp : 3 ≤ p) :
    0 < commonForbiddenFactor p := by
  have h : (3 : ℝ) ≤ p := by exact_mod_cast hp
  unfold commonForbiddenFactor
  have hsq : (1 : ℝ) < ((p : ℝ) - 1) ^ 2 := by nlinarith
  exact sub_pos.mpr ((div_lt_one (by positivity)).mpr hsq)

private theorem commonForbiddenFactor_le_one (p : ℕ) : commonForbiddenFactor p ≤ 1 :=
  sub_le_self _ (by positivity)

private theorem rootMergerFactor_eq {p : ℕ} (hp : 3 ≤ p) :
    rootMergerFactor p = 1 + 1 / ((p : ℝ) - 2) := by
  have h : (3 : ℝ) ≤ p := by exact_mod_cast hp
  unfold rootMergerFactor
  field_simp [show (p : ℝ) - 2 ≠ 0 by linarith]
  ring

private theorem rootMergerFactor_ge_one {p : ℕ} (hp : 3 ≤ p) : 1 ≤ rootMergerFactor p := by
  have h : (3 : ℝ) ≤ p := by exact_mod_cast hp
  rw [rootMergerFactor_eq hp]
  exact le_add_of_nonneg_right (div_nonneg (by norm_num) (by linarith))

/-- The fraction of the line residues surviving the restrictions,
divided by the two marginal fractions. The parameter is the normalized
first coordinate, so the other coordinate is zero exactly at N. -/
noncomputable def normalizedGoodPrimeFactor (p : ℕ) [NeZero p] (N : ZMod p)
    (restrictLeft restrictRight : Bool) : ℝ :=
  let live := Finset.univ.filter (fun x : ZMod p =>
    (restrictLeft = true → x ≠ 0) ∧ (restrictRight = true → N - x ≠ 0))
  ((live.card : ℝ) / p) /
    ((if restrictLeft then 1 - 1 / (p : ℝ) else 1) *
      (if restrictRight then 1 - 1 / (p : ℝ) else 1))

/-- At a good prime, either single restriction cancels with its
marginal. The two-sided restriction supplies exactly the factors
appearing in the singular series. -/
theorem normalized_good_prime_factor {p : ℕ} [NeZero p] (hp : 3 ≤ p) (N : ZMod p)
    (restrictLeft restrictRight : Bool) :
    normalizedGoodPrimeFactor p N restrictLeft restrictRight =
      if restrictLeft && restrictRight then
        commonForbiddenFactor p * (if N = 0 then rootMergerFactor p else 1)
      else 1 := by
  have hpR : (3 : ℝ) ≤ p := by exact_mod_cast hp
  have hp0 : (p : ℝ) ≠ 0 := by linarith
  have hpm1 : (p : ℝ) - 1 ≠ 0 := by linarith
  have hpm2 : (p : ℝ) - 2 ≠ 0 := by linarith
  have one_excluded (a : ZMod p) : (Finset.univ.filter (fun x : ZMod p => x ≠ a)).card = p - 1 := by
    have heq : Finset.univ.filter (fun x : ZMod p => x ≠ a) = Finset.univ.erase a := by ext x; simp
    rw [heq, Finset.card_erase_of_mem (Finset.mem_univ _), Finset.card_univ, ZMod.card]
  cases restrictLeft <;> cases restrictRight
  · simp [normalizedGoodPrimeFactor, ZMod.card, hp0]
  · simp only [normalizedGoodPrimeFactor, Bool.false_eq_true, false_implies, true_and,
      forall_const, Bool.false_and, Bool.false_eq_true, if_false, if_true, one_mul]
    simp_rw [sub_ne_zero, ne_comm (a := N)]
    rw [one_excluded]
    push_cast [Nat.cast_sub (by omega : 1 ≤ p)]
    field_simp
  · simp only [normalizedGoodPrimeFactor, Bool.false_eq_true, false_implies, and_true,
      forall_const, Bool.true_and, Bool.false_eq_true, if_false, if_true, mul_one]
    rw [one_excluded]
    push_cast [Nat.cast_sub (by omega : 1 ≤ p)]
    field_simp
  · simp only [normalizedGoodPrimeFactor, forall_const, Bool.true_and, if_true]
    rw [goodPrime_unit_pair_count]
    by_cases hN : N = 0 <;> simp only [hN, if_true, if_false, mul_one]
    · push_cast [Nat.cast_sub (by omega : 1 ≤ p)]
      unfold commonForbiddenFactor rootMergerFactor
      field_simp
      ring
    · push_cast [Nat.cast_sub (by omega : 2 ≤ p)]
      unfold commonForbiddenFactor
      field_simp
      ring

/-- The finite target factor J_E(N). -/
noncomputable def pointwiseTargetFactor (bad : ℕ → Prop) [DecidablePred bad] (n : ℕ) : ℝ :=
  ∏ p ∈ n.primeFactors.filter bad, rootMergerFactor p

/-- The finite target factor is at least one. -/
theorem pointwise_target_factor_ge_one {bad : ℕ → Prop} [DecidablePred bad]
    (hbad : ∀ p, bad p → 3 ≤ p) (n : ℕ) :
    1 ≤ pointwiseTargetFactor bad n :=
  Finset.one_le_prod (fun p hp => rootMergerFactor_ge_one (hbad p (Finset.mem_filter.mp hp).2))

/-- The finite local product factors into the constant part and the
target correction. This also authenticates the finite truncation. -/
theorem finite_good_prime_product {bad : ℕ → Prop} [DecidablePred bad]
    (hbad : ∀ p, bad p → p.Prime ∧ 3 ≤ p)
    (s : Finset ℕ) {n : ℕ} (hn : n ≠ 0) :
    (∏ p ∈ s.filter bad, commonForbiddenFactor p * (if p ∣ n then rootMergerFactor p else 1)) =
      (∏ p ∈ s.filter bad, commonForbiddenFactor p) *
        (∏ p ∈ s ∩ n.primeFactors.filter bad, rootMergerFactor p) := by
  rw [Finset.prod_mul_distrib, ← Finset.prod_filter]
  have heq : ((s.filter bad).filter (fun p => p ∣ n)) = s ∩ n.primeFactors.filter bad := by
    ext p
    simp only [Finset.mem_filter, Finset.mem_inter, Nat.mem_primeFactors]
    exact ⟨fun ⟨⟨hs, hb⟩, hd⟩ => ⟨hs, ⟨(hbad p hb).1, hd, hn⟩, hb⟩,
      fun ⟨hs, ⟨_, hd, _⟩, hb⟩ => ⟨⟨hs, hb⟩, hd⟩⟩
  rw [heq]

/-- The local deficit from one in the ordinary singular product. -/
noncomputable def pointwiseExcludedMass (bad : ℕ → Prop) [DecidablePred bad] (p : ℕ) : ℝ :=
  if bad p then 1 / ((p : ℝ) - 1) ^ 2 else 0

private theorem pointwiseExcludedMass_nonneg (bad : ℕ → Prop) [DecidablePred bad] (p : ℕ) :
    0 ≤ pointwiseExcludedMass bad p := by unfold pointwiseExcludedMass; split <;> positivity

private theorem pointwiseExcludedMass_summable {bad : ℕ → Prop} [DecidablePred bad]
    (hbad : ∀ p, bad p → 3 ≤ p) : Summable (pointwiseExcludedMass bad) := by
  apply Summable.of_nonneg_of_le (pointwiseExcludedMass_nonneg bad)
    (f := fun p : ℕ => 4 * (1 / (p : ℝ) ^ 2)) _
    ((Real.summable_one_div_nat_pow.mpr (by decide : 1 < (2 : ℕ))).mul_left 4)
  intro p
  unfold pointwiseExcludedMass
  split
  · rename_i hp
    have hpR : (3 : ℝ) ≤ p := by exact_mod_cast hbad p hp
    have hsq : (p : ℝ) ^ 2 ≤ 4 * ((p : ℝ) - 1) ^ 2 := by nlinarith
    rw [mul_one_div]
    have hp1 : (0 : ℝ) < (p : ℝ) - 1 := by linarith
    have hp0 : (0 : ℝ) < p := by linarith
    exact (div_le_div_iff₀ (pow_pos hp1 _) (pow_pos hp0 _)).mpr (by simpa using hsq)
  · positivity

private theorem excludedFactor_pos {bad : ℕ → Prop} [DecidablePred bad]
    (hbad : ∀ p, bad p → 3 ≤ p) (p : ℕ) : 0 < 1 - pointwiseExcludedMass bad p := by
  unfold pointwiseExcludedMass
  split
  · rename_i hp
    exact commonForbiddenFactor_pos (hbad p hp)
  · norm_num

/-- The ordinary, target-independent part C_E of the singular series. -/
noncomputable def pointwiseSingularConstant (bad : ℕ → Prop) [DecidablePred bad] : ℝ :=
  ∏' p : ℕ, (1 - pointwiseExcludedMass bad p)

/-- The ordinary product converges absolutely to a positive value. -/
theorem pointwise_singular_constant_converges {bad : ℕ → Prop} [DecidablePred bad]
    (hbad : ∀ p, bad p → 3 ≤ p) :
    Summable (fun p : ℕ => ‖(1 - pointwiseExcludedMass bad p) - 1‖) ∧
      HasProd (fun p : ℕ => 1 - pointwiseExcludedMass bad p) (pointwiseSingularConstant bad) ∧
      0 < pointwiseSingularConstant bad := by
  have hs := pointwiseExcludedMass_summable hbad
  have hlog : Summable (fun p : ℕ => Real.log (1 - pointwiseExcludedMass bad p)) := by
    simpa only [sub_eq_add_neg] using Real.summable_log_one_add_of_summable hs.neg
  have hprod := Real.hasProd_of_hasSum_log (excludedFactor_pos hbad) hlog.hasSum
  refine ⟨?_, ?_, ?_⟩
  · simpa only [sub_sub_cancel_left, norm_neg, Real.norm_of_nonneg (pointwiseExcludedMass_nonneg bad _)]
      using hs
  · unfold pointwiseSingularConstant
    rw [hprod.tprod_eq]
    exact hprod
  · unfold pointwiseSingularConstant
    rw [hprod.tprod_eq]
    exact Real.exp_pos _

private theorem one_sub_finite_prod_le_sum {ι : Type*} (s : Finset ι) (u : ι → ℝ)
    (h0 : ∀ i, 0 ≤ u i) (h1 : ∀ i, u i ≤ 1) :
    1 - (∏ i ∈ s, (1 - u i)) ≤ ∑ i ∈ s, u i := by
  classical
  induction s using Finset.induction with
  | empty => simp
  | insert a s ha ih =>
    rw [Finset.prod_insert ha, Finset.sum_insert ha]
    have hl : ∏ i ∈ s, (1 - u i) ≤ 1 :=
      Finset.prod_le_one (fun i _ => sub_nonneg.mpr (h1 i)) (fun i _ => sub_le_self _ (h0 i))
    have hm := mul_le_mul_of_nonneg_left hl (h0 a)
    nlinarith only [hm, ih]

private theorem reciprocal_square_tail_le {y : ℕ} (hy : 2 ≤ y) :
    (∑' j : ℕ, 1 / (((y : ℝ) + j) ^ 2)) ≤ 1 / ((y : ℝ) - 1) := by
  apply Real.tsum_le_of_sum_range_le (fun _ => by positivity)
  intro n
  let F (j : ℕ) : ℝ := -(1 / ((y : ℝ) + j - 1))
  have hyR : (2 : ℝ) ≤ y := by exact_mod_cast hy
  have hp (j : ℕ) : 1 / (((y : ℝ) + j) ^ 2) ≤ F (j + 1) - F j := by
    have hj : (0 : ℝ) ≤ j := Nat.cast_nonneg _
    unfold F
    push_cast
    have h₁ : 0 < (y : ℝ) + j - 1 := by linarith
    have h₂ : 0 < (y : ℝ) + j := by linarith
    rw [show (y : ℝ) + (j + 1) - 1 = (y : ℝ) + j by ring]
    field_simp
    nlinarith
  calc
    _ ≤ ∑ j ∈ Finset.range n, (F (j + 1) - F j) := Finset.sum_le_sum (fun j _ => hp j)
    _ = F n - F 0 := Finset.sum_range_sub F n
    _ ≤ _ := by
      unfold F
      have hpos : 0 ≤ (y : ℝ) + n - 1 := by linarith [Nat.cast_nonneg (α := ℝ) n]
      simp only [Nat.cast_zero, add_zero]
      linarith only [div_nonneg (by norm_num : (0 : ℝ) ≤ 1) hpos]

/-- The ordinary singular product truncated at y. -/
noncomputable def pointwiseSingularHead (bad : ℕ → Prop) [DecidablePred bad] (y : ℕ) : ℝ :=
  ∏ p ∈ Finset.range (y + 1), (1 - pointwiseExcludedMass bad p)

/-- The ordinary tail contributes at most 1/(y-1) in absolute error. -/
theorem pointwise_singular_head_error {bad : ℕ → Prop} [DecidablePred bad]
    (hbad : ∀ p, bad p → 3 ≤ p) {y : ℕ} (hy : 2 ≤ y) :
    0 ≤ pointwiseSingularHead bad y - pointwiseSingularConstant bad ∧
      pointwiseSingularHead bad y - pointwiseSingularConstant bad ≤ 1 / ((y : ℝ) - 1) := by
  let s := Finset.range (y + 1)
  let f (p : ℕ) := 1 - pointwiseExcludedMass bad p
  let tail (p : ℕ) := if y < p then pointwiseExcludedMass bad p else 0
  let Q := ∏' p : ℕ, (1 - tail p)
  have hTail : Summable tail :=
    Summable.of_nonneg_of_le (fun p => by dsimp [tail]; split <;> simp [pointwiseExcludedMass_nonneg])
      (fun p => by dsimp [tail]; split <;> simp [pointwiseExcludedMass_nonneg]) (pointwiseExcludedMass_summable hbad)
  have hfactor (p) : 0 < 1 - tail p := by
    dsimp [tail]
    split
    · exact excludedFactor_pos hbad p
    · norm_num
  have hmul : Multipliable (fun p => 1 - tail p) :=
    Real.multipliable_of_summable_log hfactor
      (by simpa only [sub_eq_add_neg] using Real.summable_log_one_add_of_summable hTail.neg)
  have htail_le : Q ≤ 1 := by
    apply le_of_tendsto' (f := fun t : Finset ℕ => ∏ p ∈ t, (1 - tail p)) hmul.hasProd
    intro t
    apply Finset.prod_le_one (fun p _ => (hfactor p).le)
    intro p _
    exact sub_le_self _ (by dsimp [tail]; split <;> simp [pointwiseExcludedMass_nonneg])
  have heq : pointwiseSingularConstant bad = pointwiseSingularHead bad y * Q := by
    let head (p : ℕ) := if p ≤ y then f p else 1
    have hf (p : ℕ) : f p = head p * (1 - tail p) := by
      by_cases hp : p ≤ y
      · simp [head, tail, hp, Nat.not_lt_of_ge hp]
      · simp [f, head, tail, hp, Nat.lt_of_not_ge hp]
    have hout (p : ℕ) (hp : p ∉ s) : head p = 1 :=
      if_neg (fun hp' => hp (Finset.mem_range.mpr (Nat.lt_succ_of_le hp')))
    unfold pointwiseSingularConstant
    change (∏' p : ℕ, f p) = _
    simp_rw [hf]
    rw [(multipliable_of_ne_finset_one hout).tprod_mul hmul, tprod_eq_prod hout]
    congr 1
    exact Finset.prod_congr rfl (fun p hp => if_pos (Nat.lt_succ_iff.mp (Finset.mem_range.mp hp)))
  have htail_sum : (∑' p : ℕ, tail p) ≤ 1 / ((y : ℝ) - 1) := by
    let g (j : ℕ) := y + 1 + j
    have hinj : Function.Injective g := by intro a b hab; dsimp [g] at hab; omega
    have hsupp : Function.support tail ⊆ Set.range g := by
      intro p hp
      have hpy : y < p := by by_contra h; exact hp (if_neg h)
      exact ⟨p - (y + 1), by dsimp [g]; omega⟩
    have hsq : Summable (fun j : ℕ => 1 / (((y : ℝ) + j) ^ 2)) := by
      simpa only [Function.comp_def, Nat.cast_add] using
        (Real.summable_one_div_nat_pow.mpr (by decide : 1 < (2 : ℕ))).comp_injective
          (fun a b h => Nat.add_left_cancel h : Function.Injective (fun j : ℕ => y + j))
    rw [← hinj.tsum_eq hsupp]
    have hpoint (j : ℕ) : tail (g j) ≤ 1 / (((y : ℝ) + j) ^ 2) := by
      have hcut : y < g j := by dsimp [g]; omega
      rw [show tail (g j) = pointwiseExcludedMass bad (g j) from if_pos hcut]
      unfold pointwiseExcludedMass
      split
      · dsimp [g]; push_cast
        rw [show (y : ℝ) + 1 + j - 1 = (y : ℝ) + j by ring]
      · positivity
    exact ((hTail.comp_injective hinj).tsum_le_tsum hpoint hsq).trans (reciprocal_square_tail_le hy)
  have hQ : 1 - 1 / ((y : ℝ) - 1) ≤ Q := by
    apply ge_of_tendsto' (f := fun t : Finset ℕ => ∏ p ∈ t, (1 - tail p)) hmul.hasProd
    intro t
    have ht₀ (p : ℕ) : 0 ≤ tail p := by dsimp [tail]; split <;> simp [pointwiseExcludedMass_nonneg]
    have ht₁ (p : ℕ) : tail p ≤ 1 := by linarith only [hfactor p]
    have h := (one_sub_finite_prod_le_sum t tail ht₀ ht₁).trans
      ((hTail.sum_le_tsum _ (fun p _ => ht₀ p)).trans htail_sum)
    linarith only [h]
  have hhead₀ : 0 ≤ pointwiseSingularHead bad y := Finset.prod_nonneg (fun p _ => (excludedFactor_pos hbad p).le)
  have hhead₁ : pointwiseSingularHead bad y ≤ 1 :=
    Finset.prod_le_one (fun p _ => (excludedFactor_pos hbad p).le)
      (fun p _ => sub_le_self _ (pointwiseExcludedMass_nonneg bad p))
  rw [heq, ← mul_one_sub]
  refine ⟨mul_nonneg hhead₀ (sub_nonneg.mpr htail_le), ?_⟩
  exact (mul_le_of_le_one_left (sub_nonneg.mpr htail_le) hhead₁).trans (by linarith only [hQ])


/-- The finite target factor truncated to primes at most y. -/
noncomputable def pointwiseTargetHead (bad : ℕ → Prop) [DecidablePred bad] (n y : ℕ) : ℝ :=
  ∏ p ∈ (n.primeFactors.filter bad).filter (fun p => p ≤ y), rootMergerFactor p

/-- The missing target-prime factor is bounded uniformly in the
valuations by exp(2 log n / (y log y)). -/
theorem pointwise_target_head_bounds {bad : ℕ → Prop} [DecidablePred bad]
    (hbad : ∀ p, bad p → 3 ≤ p) {n y : ℕ} (hn : 0 < n) (hy : 4 ≤ y) :
    pointwiseTargetHead bad n y ≤ pointwiseTargetFactor bad n ∧
      pointwiseTargetFactor bad n ≤ pointwiseTargetHead bad n y *
        Real.exp (2 * Real.log n / ((y : ℝ) * Real.log y)) := by
  let s := n.primeFactors.filter bad
  let tail := s.filter (fun p => ¬ p ≤ y)
  let f (p : ℕ) : ℝ := if y < p then 1 / ((p : ℝ) - 2) else 0
  have hyR : (4 : ℝ) ≤ y := by exact_mod_cast hy
  have hy0 : (0 : ℝ) < y := by linarith
  have hylog : 0 < Real.log (y : ℝ) := Real.log_pos (by linarith)
  have hf (p : ℕ) : 0 ≤ f p := by
    unfold f
    split
    · rename_i h
      have hp : (y : ℝ) < p := by exact_mod_cast h
      exact div_nonneg (by norm_num) (by linarith)
    · norm_num
  have hcut (p : ℕ) (hp : p ∈ tail) : y < p :=
    Nat.lt_of_not_ge (Finset.mem_filter.mp hp).2
  have htailprod : (∏ p ∈ tail, rootMergerFactor p) ≤ Real.exp (∑ p ∈ tail, f p) := by
    rw [show (∏ p ∈ tail, rootMergerFactor p) = ∏ p ∈ tail, (1 + f p) from
      Finset.prod_congr rfl (fun p hp => by
        have hp3 : 3 ≤ p := by have := hcut p hp; omega
        rw [rootMergerFactor_eq hp3, show f p = 1 / ((p : ℝ) - 2) from if_pos (hcut p hp)])]
    exact Real.prod_one_add_le_exp_sum tail hf
  have hcard : (tail.card : ℝ) * Real.log y ≤ Real.log n := by
    have hpmin (p : ℕ) (hp : p ∈ tail) : (y : ℝ) < p := by exact_mod_cast hcut p hp
    have hpos (p : ℕ) (hp : p ∈ tail) : (0 : ℝ) < p := hy0.trans (hpmin p hp)
    have hProdDvd : (∏ p ∈ tail, p) ∣ n :=
      (Finset.prod_dvd_prod_of_subset _ _ id
        ((Finset.filter_subset _ _).trans (Finset.filter_subset _ _))).trans (Nat.prod_primeFactors_dvd n)
    have hProd : ∏ p ∈ tail, p ≤ n := Nat.le_of_dvd hn hProdDvd
    have hNN : (∏ p ∈ tail, (p : ℝ)) ≤ n := by
      rw [← Nat.cast_prod]
      exact_mod_cast hProd
    calc
      (tail.card : ℝ) * Real.log y = ∑ p ∈ tail, Real.log (y : ℝ) := by
        rw [Finset.sum_const, nsmul_eq_mul]
      _ ≤ ∑ p ∈ tail, Real.log (p : ℝ) := Finset.sum_le_sum
        (fun p hp => Real.log_le_log hy0 (hpmin p hp).le)
      _ = Real.log (∏ p ∈ tail, (p : ℝ)) := (Real.log_prod (fun p hp => (hpos p hp).ne')).symm
      _ ≤ Real.log n := Real.log_le_log (Finset.prod_pos hpos) hNN
  have hsum : (∑ p ∈ tail, f p) ≤ 2 * Real.log n / ((y : ℝ) * Real.log y) := by
    calc
      _ ≤ ∑ p ∈ tail, (2 : ℝ) / y := Finset.sum_le_sum (fun p hp => by
        have hpy : (y : ℝ) < p := by exact_mod_cast hcut p hp
        rw [show f p = 1 / ((p : ℝ) - 2) from if_pos (hcut p hp)]
        apply (div_le_div_iff₀ (by linarith) hy0).mpr
        linarith)
      _ = (tail.card : ℝ) * (2 / y) := by rw [Finset.sum_const, nsmul_eq_mul]
      _ ≤ _ := by
        apply (le_div_iff₀ (mul_pos hy0 hylog)).mpr
        calc
          _ = 2 * ((tail.card : ℝ) * Real.log y) := by field_simp
          _ ≤ _ := mul_le_mul_of_nonneg_left hcard (by norm_num)
  have heq : pointwiseTargetFactor bad n =
      pointwiseTargetHead bad n y * ∏ p ∈ tail, rootMergerFactor p :=
    (Finset.prod_filter_mul_prod_filter_not s _ _).symm
  have hhead : 0 ≤ pointwiseTargetHead bad n y := Finset.prod_nonneg (fun p hp => by
    have hb : bad p := (Finset.mem_filter.mp (Finset.mem_filter.mp hp).1).2
    exact (by norm_num : (0 : ℝ) ≤ 1).trans (rootMergerFactor_ge_one (hbad p hb)))
  have hOne : 1 ≤ ∏ p ∈ tail, rootMergerFactor p := Finset.one_le_prod (fun p hp =>
    rootMergerFactor_ge_one (hbad p (Finset.mem_filter.mp (Finset.mem_filter.mp hp).1).2))
  rw [heq]
  refine ⟨?_, mul_le_mul_of_nonneg_left (htailprod.trans (Real.exp_le_exp.mpr hsum)) hhead⟩
  simpa only [mul_one] using mul_le_mul_of_nonneg_left hOne hhead

end QuadraticNorms
