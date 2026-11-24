import QuadraticNorms.MinusEightShell
import QuadraticNorms.ProgressionTail

open scoped BigOperators

namespace QuadraticNorms

/-- The nonsplit torus normalizer's proportion of the projective quotient. -/
noncomputable def minusEightLocalMass (p : ℕ) : ℝ := 2 / ((p : ℝ) * ((p : ℝ) - 1))

private theorem minusEight_inert_at_least_five {p : ℕ} (hp : MinusEightInertPrime p) : 5 ≤ p := by
  have h := Nat.mod_le p 8
  rcases hp.2 with h' | h' <;> omega

theorem minusEight_local_mass_nonneg (p : ℕ) : 0 ≤ minusEightLocalMass p := by
  rcases p with _ | p
  · norm_num [minusEightLocalMass]
  · unfold minusEightLocalMass
    apply div_nonneg (by norm_num)
    push_cast
    simpa only [add_sub_cancel_right] using mul_nonneg (show 0 ≤ (p : ℝ) + 1 by positivity)
      (show 0 ≤ (p : ℝ) by positivity)

theorem minusEight_inverse_square_le_local_mass {p : ℕ} (hp : MinusEightInertPrime p) :
    inverseSquare p ≤ minusEightLocalMass p := by
  have hp' : (5 : ℝ) ≤ p := by exact_mod_cast minusEight_inert_at_least_five hp
  have hpos : (0 : ℝ) < p := by linarith
  unfold inverseSquare minusEightLocalMass
  apply (div_le_div_iff₀ (sq_pos_of_pos hpos) (mul_pos hpos (by linarith))).mpr
  nlinarith only [sq_nonneg (p : ℝ), hpos.le]

private theorem minusEight_local_mass_le {p : ℕ} (hp : MinusEightInertPrime p) :
    minusEightLocalMass p ≤ (5 / 2) * inverseSquare p := by
  have hp' : (5 : ℝ) ≤ p := by exact_mod_cast minusEight_inert_at_least_five hp
  have hpos : (0 : ℝ) < p := by linarith
  unfold inverseSquare minusEightLocalMass
  rw [mul_one_div]
  apply (div_le_div_iff₀ (mul_pos hpos (by linarith)) (sq_pos_of_pos hpos)).mpr
  nlinarith only [mul_nonneg hpos.le (show 0 ≤ (p : ℝ) - 5 by linarith)]

private theorem minusEight_residue_card_le (s : Finset ℕ) (P r : ℕ) (hr : r < 8)
    (hs : ∀ p ∈ s, p ≤ P ∧ p % 8 = r) :
    (s.card : ℝ) ≤ (P : ℝ) / 8 + 1 - (r : ℝ) / 8 := by
  have hcard : s.card ≤ (P + 8 - r) / 8 := by
    rw [← Finset.card_range ((P + 8 - r) / 8)]
    apply Finset.card_le_card_of_injOn (fun p => p / 8)
    · intro p hp
      have ⟨hP, hm⟩ := hs p hp
      rw [Finset.mem_coe, Finset.mem_range]
      change p / 8 < (P + 8 - r) / 8
      omega
    · intro p hp q hq hd
      have hp' := (hs p hp).2
      have hq' := (hs q hq).2
      change p / 8 = q / 8 at hd
      omega
  calc
    (s.card : ℝ) ≤ (((P + 8 - r) / 8 : ℕ) : ℝ) := by exact_mod_cast hcard
    _ ≤ (((P + 8 - r : ℕ) : ℝ) / 8) := Nat.cast_div_le
    _ = _ := by rw [Nat.cast_sub (by omega), Nat.cast_add]; norm_num; ring

/-- The inert primes up to P occupy only the classes five and seven modulo eight. -/
theorem minusEight_inert_card_le {s : Finset ℕ} {P : ℕ}
    (hs : ∀ p ∈ s, MinusEightInertPrime p ∧ p ≤ P) :
    (s.card : ℝ) ≤ (P : ℝ) / 4 + 1 := by
  let c₅ := s.filter (fun p => p % 8 = 5)
  let c₇ := s.filter (fun p => p % 8 ≠ 5)
  have h₅ := minusEight_residue_card_le c₅ P 5 (by decide)
    (fun p hp => ⟨(hs p (Finset.mem_filter.mp hp).1).2, (Finset.mem_filter.mp hp).2⟩)
  have h₇ := minusEight_residue_card_le c₇ P 7 (by decide) (fun p hp => by
    obtain ⟨h₁, h₂⟩ := Finset.mem_filter.mp hp
    exact ⟨(hs p h₁).2, (hs p h₁).1.2.resolve_left h₂⟩)
  have hc : (c₅.card : ℝ) + c₇.card = s.card := by
    exact_mod_cast Finset.card_filter_add_card_filter_not (s := s) (fun p => p % 8 = 5)
  linarith

private def minusEightNextResidue (P r : ℕ) : ℕ :=
  P + 1 + ((r + 8 - (P + 1) % 8) % 8)

private theorem minusEightResidue_tail_le (s : Finset ℕ) {P r : ℕ} (hr : r < 8)
    (hs : ∀ p ∈ s, P < p ∧ p % 8 = r) :
    (∑ p ∈ s, inverseSquare p) ≤
      inverseSquare (minusEightNextResidue P r) + 1 / (8 * (minusEightNextResidue P r : ℝ)) := by
  apply inverseSquare_progression_sum_le s (by unfold minusEightNextResidue; omega) (by decide : 0 < 8)
  intro p hp
  obtain ⟨hP, hm⟩ := hs p hp
  unfold minusEightNextResidue
  omega

private theorem minusEight_inert_tail_head {s : Finset ℕ} {P : ℕ}
    (hs : ∀ p ∈ s, MinusEightInertPrime p ∧ P < p) :
    (∑ p ∈ s, inverseSquare p) ≤
      (inverseSquare (minusEightNextResidue P 5) + 1 / (8 * (minusEightNextResidue P 5 : ℝ))) +
        (inverseSquare (minusEightNextResidue P 7) + 1 / (8 * (minusEightNextResidue P 7 : ℝ))) := by
  rw [← Finset.sum_filter_add_sum_filter_not s (fun p => p % 8 = 5) inverseSquare]
  apply add_le_add
  · exact minusEightResidue_tail_le _ (by decide) fun p hp =>
      ⟨(hs p (Finset.mem_filter.mp hp).1).2, (Finset.mem_filter.mp hp).2⟩
  · exact minusEightResidue_tail_le _ (by decide) fun p hp =>
      ⟨(hs p (Finset.mem_filter.mp hp).1).2, (hs p (Finset.mem_filter.mp hp).1).1.2.resolve_left
        (Finset.mem_filter.mp hp).2⟩

/-- The inverse-square tail in the sieve, also valid for a real cutoff X
with the finite prime filter taken above its natural floor. -/
theorem minusEight_inert_tail_le {s : Finset ℕ} {P : ℕ} {X : ℝ} (hX : 0 < X)
    (hPX : X ≤ (P : ℝ) + 1) (hs : ∀ p ∈ s, MinusEightInertPrime p ∧ P < p) :
    (∑ p ∈ s, inverseSquare p) ≤ 1 / (4 * X) + 2 / X ^ 2 := by
  have hterm (r : ℕ) :
      inverseSquare (minusEightNextResidue P r) + 1 / (8 * (minusEightNextResidue P r : ℝ)) ≤
        1 / X ^ 2 + 1 / (8 * X) := by
    have hle : (P : ℝ) + 1 ≤ (minusEightNextResidue P r : ℝ) := by
      exact_mod_cast (show P + 1 ≤ minusEightNextResidue P r by unfold minusEightNextResidue; omega)
    have h := hPX.trans hle
    exact add_le_add (one_div_le_one_div_of_le (sq_pos_of_pos hX)
      (pow_le_pow_left₀ hX.le h _)) (one_div_le_one_div_of_le (by positivity)
        (mul_le_mul_of_nonneg_left h (by norm_num)))
  exact (minusEight_inert_tail_head hs).trans
    ((add_le_add (hterm 5) (hterm 7)).trans_eq (by ring))

private def minusEightInertThrough (P : ℕ) : Finset ℕ :=
  (Finset.range (P + 1)).filter MinusEightInertPrime

/-- The local inert mass is bounded by the exact finite head and the
two progression tails starting at 69 and 71. -/
theorem minusEight_inert_mass_lt {s : Finset ℕ} (hs : ∀ p ∈ s, MinusEightInertPrime p) :
    (∑ p ∈ s, minusEightLocalMass p) < 3 / 16 := by
  have hhead : ∑ p ∈ s.filter (fun p => p ≤ 64), minusEightLocalMass p ≤
      ∑ p ∈ minusEightInertThrough 64, minusEightLocalMass p := by
    apply Finset.sum_le_sum_of_subset_of_nonneg _ (fun p _ _ => minusEight_local_mass_nonneg p)
    intro p hp
    obtain ⟨h, hP⟩ := Finset.mem_filter.mp hp
    exact Finset.mem_filter.mpr ⟨Finset.mem_range.mpr (Nat.lt_succ_of_le hP), hs p h⟩
  let tail := s.filter (fun p => ¬ p ≤ 64)
  have htail : ∑ p ∈ tail, minusEightLocalMass p ≤
      (5 / 2) * ((inverseSquare 69 + 1 / (8 * 69)) + (inverseSquare 71 + 1 / (8 * 71))) := by
    calc
      _ ≤ ∑ p ∈ tail, (5 / 2) * inverseSquare p :=
        Finset.sum_le_sum (fun p hp => minusEight_local_mass_le (hs p (Finset.mem_filter.mp hp).1))
      _ = (5 / 2) * ∑ p ∈ tail, inverseSquare p := (Finset.mul_sum ..).symm
      _ ≤ _ := by
        have h := minusEight_inert_tail_head (s := tail) (P := 64) fun p hp =>
          ⟨hs p (Finset.mem_filter.mp hp).1, Nat.lt_of_not_ge (Finset.mem_filter.mp hp).2⟩
        have h₅ : minusEightNextResidue 64 5 = 69 := by decide
        have h₇ : minusEightNextResidue 64 7 = 71 := by decide
        rw [h₅, h₇] at h
        exact mul_le_mul_of_nonneg_left h (by norm_num)
  rw [← Finset.sum_filter_add_sum_filter_not s (fun p => p ≤ 64) minusEightLocalMass]
  refine (add_le_add hhead htail).trans_lt ?_
  rw [show minusEightInertThrough 64 = {5, 7, 13, 23, 29, 31, 37, 47, 53, 61} by decide +kernel]
  norm_num [minusEightLocalMass, inverseSquare, Finset.sum_insert, Finset.sum_singleton]

end QuadraticNorms
