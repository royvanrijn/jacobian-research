import Mathlib.Analysis.PSeries

open scoped BigOperators

namespace QuadraticNorms

noncomputable def inverseSquare (p : ℕ) : ℝ := 1 / (p : ℝ) ^ 2

theorem inverseSquare_nonneg (p : ℕ) : 0 ≤ inverseSquare p := by
  unfold inverseSquare
  positivity

theorem inverseSquare_summable : Summable inverseSquare :=
  Real.summable_one_div_nat_pow.mpr (by decide : 1 < (2 : ℕ))

theorem progression_tail_bound {a d : ℕ} (ha : 0 < a) (hd : 0 < d) :
    (∑' j : ℕ, inverseSquare (a + d * j)) ≤
      inverseSquare a + 1 / ((d : ℝ) * (a : ℝ)) := by
  have hinj : Function.Injective (fun j : ℕ => a + d * j) := by
    intro x y h
    exact (Nat.mul_left_cancel_iff hd).mp (Nat.add_left_cancel h)
  have hsum : Summable (fun j : ℕ => inverseSquare (a + d * j)) :=
    inverseSquare_summable.comp_injective hinj
  rw [← hsum.sum_add_tsum_nat_add 1, Finset.sum_range_one, Nat.mul_zero, Nat.add_zero]
  apply add_le_add_right
  apply Real.tsum_le_of_sum_range_le (fun j => inverseSquare_nonneg _)
  intro n
  let F (j : ℕ) : ℝ := -(1 / ((d : ℝ) * (a + d * j)))
  have haR : (0 : ℝ) < a := by exact_mod_cast ha
  have hdR : (0 : ℝ) < d := by exact_mod_cast hd
  have hpoint (j : ℕ) : inverseSquare (a + d * (j + 1)) ≤ F (j + 1) - F j := by
    unfold inverseSquare F
    push_cast
    have hprior : 0 < (a : ℝ) + d * j := by positivity
    have hnext : 0 < (a : ℝ) + d * ((j : ℝ) + 1) := by positivity
    field_simp [hprior.ne', hnext.ne', hdR.ne']
    nlinarith
  calc
    _ ≤ ∑ j ∈ Finset.range n, (F (j + 1) - F j) := Finset.sum_le_sum (fun j _ => hpoint j)
    _ = F n - F 0 := Finset.sum_range_sub F n
    _ ≤ _ := by dsimp [F]; push_cast; simp; positivity

theorem conditional_square_summable (pred : ℕ → Prop) [DecidablePred pred] :
    Summable (fun p : ℕ => if pred p then inverseSquare p else 0) :=
  Summable.of_nonneg_of_le (fun p => by split <;> simp [inverseSquare_nonneg])
    (fun p => by split <;> simp [inverseSquare_nonneg]) inverseSquare_summable

/-- A finite subset of a progression above its first term is bounded
by the inverse-square tail of the whole progression. -/
theorem inverseSquare_progression_sum_le (s : Finset ℕ) {a d : ℕ}
    (ha : 0 < a) (hd : 0 < d) (hs : ∀ p ∈ s, a ≤ p ∧ p % d = a % d) :
    (∑ p ∈ s, inverseSquare p) ≤ inverseSquare a + 1 / ((d : ℝ) * a) := by
  let j (p : ℕ) := p / d - a / d
  have hj (p : ℕ) (hp : p ∈ s) : a + d * j p = p := by
    have hle : a / d ≤ p / d := Nat.div_le_div_right (hs p hp).1
    unfold j
    conv_lhs => enter [1]; rw [← Nat.mod_add_div a d]
    rw [add_assoc, ← Nat.mul_add, Nat.add_sub_of_le hle,
      ← (hs p hp).2, Nat.mod_add_div]
  have hinj : Set.InjOn j (s : Set ℕ) := by
    intro p hp q hq h
    rw [← hj p hp, h, hj q hq]
  have hseries : Summable (fun x : ℕ => inverseSquare (a + d * x)) :=
    inverseSquare_summable.comp_injective (fun _ _ h =>
      (Nat.mul_left_cancel_iff hd).mp (Nat.add_left_cancel h))
  calc
    _ = ∑ x ∈ s.image j, inverseSquare (a + d * x) := by
      rw [Finset.sum_image hinj]
      exact Finset.sum_congr rfl (fun p hp => congrArg inverseSquare (hj p hp).symm)
    _ ≤ ∑' x : ℕ, inverseSquare (a + d * x) :=
      hseries.sum_le_tsum _ (fun x _ => inverseSquare_nonneg _)
    _ ≤ _ := progression_tail_bound ha hd

end QuadraticNorms
