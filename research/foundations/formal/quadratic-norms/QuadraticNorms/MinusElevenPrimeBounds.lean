import QuadraticNorms.MinusElevenIntervals
import QuadraticNorms.ProgressionTail
import Mathlib.Algebra.Order.Floor.Semiring

open scoped BigOperators

namespace QuadraticNorms

/-- Odd inert primes for the discriminant −11 core. -/
def MinusElevenInertPrime (p : ℕ) : Prop :=
  p.Prime ∧ p > 2 ∧ p % 11 ∈ [2, 6, 7, 8, 10]

instance : DecidablePred MinusElevenInertPrime := fun p => by
  unfold MinusElevenInertPrime
  infer_instance

private def inertResidues : Finset ℕ := {7, 13, 17, 19, 21}

private theorem inert_residue {p : ℕ} (hp : MinusElevenInertPrime p) : p % 22 ∈ inertResidues := by
  have hodd := (hp.1.mod_two_eq_one_iff_ne_two).mpr (by have := hp.2.1; omega)
  have hmod : (p % 22) % 11 = p % 11 := Nat.mod_mod_of_dvd p (by decide)
  have hpar : (p % 22) % 2 = 1 := (Nat.mod_mod_of_dvd p (by decide)).trans hodd
  have hfinite : ∀ n : Fin 22, n.val % 2 = 1 → n.val % 11 ∈ [2, 6, 7, 8, 10] →
      n.val ∈ inertResidues := by decide
  exact hfinite ⟨p % 22, Nat.mod_lt p (by decide)⟩ hpar (hmod ▸ hp.2.2)

private theorem inert_at_least_seven {p : ℕ} (hp : MinusElevenInertPrime p) : 7 ≤ p := by
  have hshort : ∀ r ∈ inertResidues, 7 ≤ r := by decide
  exact (hshort _ (inert_residue hp)).trans (Nat.mod_le _ _)

private def nextInProgression (P r : ℕ) : ℕ := P + 1 + ((r + 22 - (P + 1) % 22) % 22)

private theorem nextInProgression_spec (P : ℕ) {r : ℕ} (hr : r < 22) :
    P < nextInProgression P r ∧ nextInProgression P r ≤ P + 22 ∧ nextInProgression P r % 22 = r := by
  unfold nextInProgression
  omega

private theorem mem_progression {P r p : ℕ} (hr : r < 22) (hP : P < p) (hp : p % 22 = r) :
    ∃ j : ℕ, p = nextInProgression P r + 22 * j := by
  obtain ⟨_, hle, hres⟩ := nextInProgression_spec P hr
  exact ⟨p / 22 - nextInProgression P r / 22, by omega⟩

private noncomputable def residueTail (P r : ℕ) : ℝ :=
  ∑' p : ℕ, if P < p ∧ p % 22 = r then inverseSquare p else 0

private theorem sum_residueTail_le (P : ℕ) {r : ℕ} (hr : r < 22) :
    residueTail P r ≤ inverseSquare (nextInProgression P r) + 1 / (22 * (nextInProgression P r : ℝ)) := by
  let a := nextInProgression P r
  let f (p : ℕ) : ℝ := if P < p ∧ p % 22 = r then inverseSquare p else 0
  let g (j : ℕ) := a + 22 * j
  obtain ⟨hpa, _, hamod⟩ := nextInProgression_spec P hr
  have hinj : Function.Injective g := by intro x y h; dsimp [g] at h; omega
  have hsupport : Function.support f ⊆ Set.range g := by
    intro p hp
    have hcase : P < p ∧ p % 22 = r := by
      by_contra hneg
      exact hp (if_neg hneg)
    obtain ⟨j, heq⟩ := mem_progression hr hcase.1 hcase.2
    exact ⟨j, heq.symm⟩
  have heq (j : ℕ) : f (g j) = inverseSquare (g j) := by
    apply if_pos
    dsimp [g, a]
    omega
  unfold residueTail
  change tsum f ≤ _
  rw [← hinj.tsum_eq hsupport]
  simp_rw [heq]
  exact progression_tail_bound (by omega : 0 < a) (by decide : 0 < 22)

noncomputable def minusElevenInertTail (P : ℕ) : ℝ :=
  ∑' p : ℕ, if MinusElevenInertPrime p ∧ P < p then inverseSquare p else 0

private theorem inert_tail_le_progressions (P : ℕ) :
    minusElevenInertTail P ≤ ∑ r ∈ inertResidues,
      (inverseSquare (nextInProgression P r) + 1 / (22 * (nextInProgression P r : ℝ))) := by
  let f (r p : ℕ) := if P < p ∧ p % 22 = r then inverseSquare p else 0
  have hsum (r : ℕ) : Summable (f r) := conditional_square_summable _
  have htotal : Summable (fun p => ∑ r ∈ inertResidues, f r p) :=
    (hasSum_sum (fun r _ => (hsum r).hasSum)).summable
  calc
    minusElevenInertTail P ≤ ∑' p, ∑ r ∈ inertResidues, f r p := by
      apply (conditional_square_summable (fun p => MinusElevenInertPrime p ∧ P < p)).tsum_le_tsum _ htotal
      intro p
      split
      · rename_i h
        have hp := inert_residue h.1
        exact (show inverseSquare p = f (p % 22) p by simp [f, h.2]) ▸
          Finset.single_le_sum (f := fun r => f r p)
            (fun r _ => by dsimp [f]; split <;> simp [inverseSquare_nonneg]) hp
      · exact Finset.sum_nonneg (fun r _ => by dsimp [f]; split <;> simp [inverseSquare_nonneg])
    _ = ∑ r ∈ inertResidues, residueTail P r := Summable.tsum_finsetSum (fun r _ => hsum r)
    _ ≤ _ := Finset.sum_le_sum (fun r hr =>
      sum_residueTail_le P ((show ∀ r ∈ inertResidues, r < 22 by decide) r hr))

private theorem tail_le_at_cutoff {P : ℕ} {X : ℝ} (hX : 0 < X) (hPX : X ≤ (P : ℝ) + 1) :
    minusElevenInertTail P ≤ 5 / (22 * X) + 5 / X ^ 2 := by
  have hpoint (r : ℕ) (hr : r ∈ inertResidues) :
      inverseSquare (nextInProgression P r) + 1 / (22 * (nextInProgression P r : ℝ)) ≤
        1 / X ^ 2 + 1 / (22 * X) := by
    have hres : r < 22 := (show ∀ r ∈ inertResidues, r < 22 by decide) r hr
    have hle : (P : ℝ) + 1 ≤ nextInProgression P r := by
      exact_mod_cast Nat.succ_le_of_lt (nextInProgression_spec P hres).1
    have hle' := hPX.trans hle
    apply add_le_add
    · exact one_div_le_one_div_of_le (pow_pos hX _) (pow_le_pow_left₀ hX.le hle' _)
    · exact one_div_le_one_div_of_le (by positivity) (mul_le_mul_of_nonneg_left hle' (by norm_num))
  refine (inert_tail_le_progressions P).trans ((Finset.sum_le_sum hpoint).trans ?_)
  have hcard : inertResidues.card = 5 := by decide
  rw [Finset.sum_const, hcard, nsmul_eq_mul]
  ring_nf
  rfl


theorem minusEleven_inert_tail_nonneg (P : ℕ) : 0 ≤ minusElevenInertTail P :=
  tsum_nonneg (fun p => by split <;> simp [inverseSquare_nonneg])

/-- The inverse-square tail at an integer cutoff. -/
theorem minusEleven_inert_tail_le {P : ℕ} (hP : 0 < P) :
    minusElevenInertTail P ≤ 5 / (22 * (P : ℝ)) + 5 / (P : ℝ) ^ 2 :=
  tail_le_at_cutoff (by exact_mod_cast hP) (by linarith)

/-- Primes above a real cutoff are exactly the primes above its natural floor. -/
theorem minusEleven_inert_tail_le_real {P : ℝ} (hP : 0 < P) :
    minusElevenInertTail ⌊P⌋₊ ≤ 5 / (22 * P) + 5 / P ^ 2 :=
  tail_le_at_cutoff hP (Nat.lt_floor_add_one _).le

def minusElevenInertCount (P : ℕ) : ℕ :=
  ((Finset.range (P + 1)).filter MinusElevenInertPrime).card

private theorem residue_count_bound (P r : ℕ) (hr : r < 22) (s : Finset ℕ)
    (hs : ∀ p ∈ s, p ≤ P ∧ p % 22 = r) :
    (s.card : ℝ) ≤ (P : ℝ) / 22 + 1 - (r : ℝ) / 22 := by
  have hcard : s.card ≤ (P + 22 - r) / 22 := by
    rw [← Finset.card_range ((P + 22 - r) / 22)]
    apply Finset.card_le_card_of_injOn (fun p => p / 22)
    · intro p hp
      rw [Finset.mem_coe, Finset.mem_range]
      change p / 22 < (P + 22 - r) / 22
      have ⟨hle, heq⟩ := hs p hp
      omega
    · intro p hp q hq hdiv
      change p / 22 = q / 22 at hdiv
      have := (hs p hp).2
      have := (hs q hq).2
      omega
  calc
    (s.card : ℝ) ≤ (((P + 22 - r) / 22 : ℕ) : ℝ) := by exact_mod_cast hcard
    _ ≤ (((P + 22 - r : ℕ) : ℝ) / 22) := Nat.cast_div_le
    _ = _ := by
      rw [Nat.cast_sub (by omega), Nat.cast_add]
      ring

/-- The count bound for odd inert primes, by the five residue classes. -/
theorem minusEleven_inert_count_le (P : ℕ) :
    (minusElevenInertCount P : ℝ) ≤ (5 / 22) * P + 3 / 2 := by
  let s := (Finset.range (P + 1)).filter MinusElevenInertPrime
  have hfiber := Finset.card_eq_sum_card_fiberwise (s := s) (t := inertResidues) (f := (· % 22))
    (fun p hp => inert_residue (Finset.mem_filter.mp hp).2)
  unfold minusElevenInertCount
  change (s.card : ℝ) ≤ _
  rw [hfiber, Nat.cast_sum]
  refine (Finset.sum_le_sum (fun r hr => residue_count_bound P r ?_ _ ?_)).trans ?_
  · exact (show ∀ r ∈ inertResidues, r < 22 by decide) r hr
  · intro p hp
    obtain ⟨hps, heq⟩ := Finset.mem_filter.mp hp
    exact ⟨Nat.lt_succ_iff.mp (Finset.mem_range.mp (Finset.mem_filter.mp hps).1), heq⟩
  · simp [inertResidues]
    linarith

theorem minusEleven_inert_count_le_real {P : ℝ} (hP : 0 ≤ P) :
    (minusElevenInertCount ⌊P⌋₊ : ℝ) ≤ (5 / 22) * P + 3 / 2 :=
  (minusEleven_inert_count_le _).trans
    (add_le_add_left (mul_le_mul_of_nonneg_left (Nat.floor_le hP) (by norm_num)) _)

private theorem tsum_split_at (f : ℕ → ℝ) (hf : Summable f) (P : ℕ) :
    (∑' p, f p) =
      (∑ p ∈ Finset.range (P + 1), f p) + ∑' p : ℕ, if P < p then f p else 0 := by
  let g (p : ℕ) := p + (P + 1)
  let fTail (p : ℕ) := if P < p then f p else 0
  have hg : Function.Injective g := by intro a b h; dsimp [g] at h; omega
  have hsupport : Function.support fTail ⊆ Set.range g := by
    intro p hp
    have hle : P + 1 ≤ p := by
      by_contra! h
      exact hp (if_neg (by omega))
    exact ⟨p - (P + 1), Nat.sub_add_cancel hle⟩
  have htail (p : ℕ) : fTail (g p) = f (p + (P + 1)) := if_pos (by dsimp [g]; omega)
  change _ = _ + tsum fTail
  rw [← hg.tsum_eq hsupport]
  simp_rw [htail]
  exact (hf.sum_add_tsum_nat_add (P + 1)).symm

private theorem sum_inert_primes_through (f : ℕ → ℝ) (P : ℕ) :
    (∑ p ∈ Finset.range (P + 1), if MinusElevenInertPrime p then f p else 0) =
      ((minusElevenInertPrimesThrough P).map f).sum := by
  have hnodup : (minusElevenInertPrimesThrough P).Nodup :=
    List.Nodup.filter _ List.nodup_range
  rw [← List.sum_toFinset f hnodup, ← Finset.sum_filter]
  congr 1
  ext p
  simp [mem_minusEleven_inertPrimes, MinusElevenInertPrime]

noncomputable def minusElevenLocalMass (p : ℕ) : ℝ := 2 / ((p : ℝ) * ((p : ℝ) - 1))

noncomputable def minusElevenInertMass : ℝ :=
  ∑' p : ℕ, if MinusElevenInertPrime p then minusElevenLocalMass p else 0

private theorem inert_mass_term_bounds {p : ℕ} (hp : MinusElevenInertPrime p) :
    0 ≤ minusElevenLocalMass p ∧ minusElevenLocalMass p ≤ (7 / 3) * inverseSquare p := by
  have hp' : (7 : ℝ) ≤ p := by exact_mod_cast inert_at_least_seven hp
  have hzero : (0 : ℝ) < p := by linarith
  have hsub : (0 : ℝ) < (p : ℝ) - 1 := by linarith
  unfold minusElevenLocalMass inverseSquare
  constructor
  · exact div_nonneg (by norm_num) (mul_nonneg hzero.le hsub.le)
  · rw [mul_one_div]
    apply (div_le_div_iff₀ (mul_pos hzero hsub) (pow_pos hzero _)).mpr
    nlinarith only [mul_nonneg hzero.le (show 0 ≤ (p : ℝ) - 7 by linarith)]

private theorem inert_mass_summable :
    Summable (fun p : ℕ => if MinusElevenInertPrime p then minusElevenLocalMass p else 0) := by
  let f (p : ℕ) : ℝ := if MinusElevenInertPrime p then minusElevenLocalMass p else 0
  have hnonneg (p : ℕ) : 0 ≤ f p := by
    by_cases hp : MinusElevenInertPrime p
    · simpa [f, hp] using (inert_mass_term_bounds hp).1
    · simp [f, hp]
  have hmajor (p : ℕ) : f p ≤ (7 / 3) * inverseSquare p := by
    by_cases hp : MinusElevenInertPrime p
    · simpa [f, hp] using (inert_mass_term_bounds hp).2
    · simpa [f, hp] using mul_nonneg (by norm_num : (0 : ℝ) ≤ 7 / 3) (inverseSquare_nonneg p)
  exact Summable.of_nonneg_of_le hnonneg hmajor (inverseSquare_summable.mul_left (7 / 3))

theorem minusEleven_inert_mass_nonneg : 0 ≤ minusElevenInertMass := by
  apply tsum_nonneg
  intro p
  by_cases hp : MinusElevenInertPrime p
  · simpa [hp] using (inert_mass_term_bounds hp).1
  · simp [hp]

theorem minusEleven_inert_mass_sum_le {s : Finset ℕ} (hs : ∀ p ∈ s, MinusElevenInertPrime p) :
    (∑ p ∈ s, minusElevenLocalMass p) ≤ minusElevenInertMass := by
  have h := inert_mass_summable.sum_le_tsum s
    (fun p _ => by by_cases hp : MinusElevenInertPrime p <;>
      simp [hp, inert_mass_term_bounds])
  rwa [Finset.sum_congr rfl (fun p hp => if_pos (hs p hp))] at h

theorem minusEleven_inert_tail_sum_le {s : Finset ℕ} {P : ℕ}
    (hs : ∀ p ∈ s, MinusElevenInertPrime p ∧ P < p) :
    (∑ p ∈ s, (1 : ℝ) / (p : ℝ) ^ 2) ≤ minusElevenInertTail P := by
  have h := (conditional_square_summable (fun p => MinusElevenInertPrime p ∧ P < p)).sum_le_tsum s
    (fun p _ => by split <;> simp [inverseSquare_nonneg])
  rwa [Finset.sum_congr rfl (fun p hp => if_pos (hs p hp))] at h

theorem minusEleven_inert_card_le {s : Finset ℕ} {P : ℕ}
    (hs : ∀ p ∈ s, MinusElevenInertPrime p ∧ p ≤ P) :
    s.card ≤ minusElevenInertCount P :=
  Finset.card_le_card (fun p hp =>
    Finset.mem_filter.mpr ⟨Finset.mem_range.mpr (Nat.lt_succ_of_le (hs p hp).2), (hs p hp).1⟩)

/-- The total local inert mass is below 9/100. The finite head is exact
and each progression tail is bounded from its first possible term. -/
theorem minusEleven_inert_mass_lt : minusElevenInertMass < 9 / 100 := by
  let f (p : ℕ) : ℝ := if MinusElevenInertPrime p then minusElevenLocalMass p else 0
  have hnonneg (p : ℕ) : 0 ≤ f p := by
    by_cases hp : MinusElevenInertPrime p
    · simpa [f, hp] using (inert_mass_term_bounds hp).1
    · simp [f, hp]
  have hf : Summable f := inert_mass_summable
  let tail (p : ℕ) : ℝ := if 64 < p then f p else 0
  have htailS := conditional_square_summable (fun p => MinusElevenInertPrime p ∧ 64 < p)
  have hcmp (p : ℕ) : tail p ≤ (7 / 3) * (if MinusElevenInertPrime p ∧ 64 < p then inverseSquare p else 0) := by
    by_cases hp : MinusElevenInertPrime p <;> by_cases hc : 64 < p <;>
      simp [tail, f, hp, hc, inert_mass_term_bounds]
  have htail := (Summable.of_nonneg_of_le (fun p => by dsimp [tail]; split <;> simp [hnonneg])
    hcmp (htailS.mul_left (7 / 3))).tsum_le_tsum hcmp (htailS.mul_left (7 / 3))
  rw [htailS.tsum_mul_left] at htail
  change tsum tail ≤ (7 / 3) * minusElevenInertTail 64 at htail
  have hsum : minusElevenInertPrimesThrough 64 = [7, 13, 17, 19, 29, 41, 43, 61] := by decide +kernel
  change tsum f < _
  rw [tsum_split_at f hf 64, sum_inert_primes_through, hsum]
  apply lt_of_le_of_lt (add_le_add_right (htail.trans (mul_le_mul_of_nonneg_left
    (inert_tail_le_progressions 64) (by norm_num))) _)
  norm_num [minusElevenLocalMass, inverseSquare, inertResidues, nextInProgression]


private theorem sum_after_eq_filtered_sum (f : ℕ → ℝ) (P : ℕ) (primes : List ℕ) :
    (primes.map (fun p => if P < p then f p else 0)).sum =
      ((primes.filter (fun p => decide (P < p))).map f).sum := by
  induction primes with
  | nil => rfl
  | cons p rest ih =>
    by_cases h : P < p <;> simp [h, ih]

/-- The sieve count in a verified row is the exact inert-prime count. -/
theorem MinusElevenBin.inertCount_eq {row : MinusElevenBin} (hr : row.valid) :
    row.inertCount = minusElevenInertCount row.primeCutoff := by
  let primes := (minusElevenInertPrimesThrough minusElevenSieveBound).filter
    (fun p => decide (p ≤ row.primeCutoff))
  change primes.length = _
  have hnd : primes.Nodup := (List.nodup_range.filter _).filter _
  have hset : primes.toFinset = (Finset.range (row.primeCutoff + 1)).filter MinusElevenInertPrime := by
    ext p
    have hle : row.primeCutoff ≤ minusElevenSieveBound := hr.2.2.2.2.2.1.le
    simp only [primes, List.mem_toFinset, List.mem_filter, mem_minusEleven_inertPrimes,
      decide_eq_true_eq, Finset.mem_filter, Finset.mem_range, Nat.lt_succ_iff, MinusElevenInertPrime]
    constructor
    · rintro ⟨⟨_, hp, htwo, hmod⟩, hP⟩
      exact ⟨hP, hp, htwo, hmod⟩
    · rintro ⟨hP, hp, htwo, hmod⟩
      exact ⟨⟨hP.trans hle, hp, htwo, hmod⟩, hP⟩
  rw [← List.toFinset_card_of_nodup hnd, hset]
  rfl

private theorem cast_tail_inequality (primes : List ℕ) (R : ℕ) (bound : ℚ)
    (h : (primes.map (fun (p : ℕ) => (1 : ℚ) / (p : ℚ) ^ 2)).sum +
      5 / (22 * (R : ℚ)) + 5 / (R : ℚ) ^ 2 ≤ bound) :
    (primes.map inverseSquare).sum +
      (5 / (22 * (R : ℝ)) + 5 / (R : ℝ) ^ 2) ≤ (bound : ℝ) := by
  have hreal := (Rat.cast_le (K := ℝ)).mpr h
  push_cast [Rat.cast_list_sum] at hreal
  rw [show inverseSquare = (fun p : ℕ => ((p : ℝ) ^ 2)⁻¹) by funext p; simp [inverseSquare]]
  simpa [List.map_map, Function.comp_def, ← add_assoc] using hreal

/-- The fixed-point tail in a verified row bounds the entire inverse-square
series beyond its cutoff, including primes above the sieve limit. -/
theorem MinusElevenBin.inertTail_le {row : MinusElevenBin} (hr : row.valid) :
    minusElevenInertTail row.primeCutoff ≤ (row.inertTailBound : ℝ) := by
  let P := row.primeCutoff
  let R := minusElevenSieveBound
  let f (p : ℕ) : ℝ := if MinusElevenInertPrime p ∧ P < p then inverseSquare p else 0
  let primes := (minusElevenInertPrimesThrough R).filter (fun p => decide (P < p))
  have hPR : P < R := hr.2.2.2.2.2.1
  have hf : Summable f := conditional_square_summable _
  have htail (p : ℕ) : (if R < p then f p else 0) =
      if MinusElevenInertPrime p ∧ R < p then inverseSquare p else 0 := by
    by_cases h : R < p
    · simp [f, h, hPR.trans h]
    · simp [h]
  have hhead (p : ℕ) : f p =
      if MinusElevenInertPrime p then (if P < p then inverseSquare p else 0) else 0 := by
    by_cases hi : MinusElevenInertPrime p <;> simp [f, hi]
  have hheadSum : (∑ p ∈ Finset.range (R + 1), f p) = (primes.map inverseSquare).sum := by
    simp_rw [hhead]
    rw [sum_inert_primes_through, sum_after_eq_filtered_sum]
  change tsum f ≤ _
  rw [tsum_split_at f hf R, hheadSum]
  simp_rw [htail]
  refine (add_le_add_right (minusEleven_inert_tail_le (by decide : 0 < R)) _).trans ?_
  exact cast_tail_inequality primes R row.inertTailBound row.rounded_tail_bound

end QuadraticNorms
