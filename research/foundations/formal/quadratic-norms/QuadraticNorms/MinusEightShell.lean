import QuadraticNorms.PrimitiveClassification
import Mathlib.Data.Nat.Sqrt

namespace QuadraticNorms

/-- A signed, ordered point in a shell for the form x² + 2y². -/
abbrev MinusEightShellPoint := (ℤ × ℤ) × (ℤ × ℤ)

def minusEightNorm (v : ℤ × ℤ) : ℕ := (v.1 ^ 2 + 2 * v.2 ^ 2).toNat

theorem minusEight_norm_cast (v : ℤ × ℤ) :
    (minusEightNorm v : ℤ) = v.1 ^ 2 + 2 * v.2 ^ 2 :=
  Int.toNat_of_nonneg (by positivity)

def minusEightVectorBox (n : ℕ) : Finset (ℤ × ℤ) :=
  let coordinates := (Finset.range (2 * n.sqrt + 1)).map
    ⟨fun k : ℕ => (k : ℤ) - (n.sqrt : ℤ), by
      intro a b h
      change (a : ℤ) - (n.sqrt : ℤ) = (b : ℤ) - (n.sqrt : ℤ) at h
      omega⟩
  coordinates ×ˢ coordinates

theorem mem_minusEight_vector_box {n : ℕ} {v : ℤ × ℤ} (h : minusEightNorm v ≤ n) :
    v ∈ minusEightVectorBox n := by
  have he := minusEight_norm_cast v
  have hn : (minusEightNorm v : ℤ) ≤ n := by exact_mod_cast h
  have hmem (z : ℤ) (hz : z ^ 2 ≤ n) :
      ∃ j ∈ Finset.range (2 * n.sqrt + 1), (j : ℤ) - (n.sqrt : ℤ) = z := by
    rw [← Int.natAbs_sq] at hz
    have hnz : (z.natAbs : ℤ) ≤ (n.sqrt : ℤ) := by
      exact_mod_cast Nat.le_sqrt'.mpr (by exact_mod_cast hz)
    rw [Int.natCast_natAbs, abs_le] at hnz
    refine ⟨(z + (n.sqrt : ℤ)).toNat, Finset.mem_range.mpr ?_, ?_⟩ <;> omega
  exact Finset.mem_product.mpr
    ⟨Finset.mem_map.mpr (hmem v.1 (by nlinarith [sq_nonneg v.2])),
      Finset.mem_map.mpr (hmem v.2 (by nlinarith [sq_nonneg v.1]))⟩

/-- The shell for Q ⊕ wQ. The two theta series use w = 1 and w = 4. -/
def minusEightWeightedShell (w n : ℕ) : Finset MinusEightShellPoint :=
  (minusEightVectorBox n ×ˢ minusEightVectorBox n).filter
    (fun v => minusEightNorm v.1 + w * minusEightNorm v.2 = n)

abbrev minusEightFullShell := minusEightWeightedShell 1

theorem mem_minusEight_weighted_shell {w n : ℕ} (hw : 0 < w) {v : MinusEightShellPoint} :
    v ∈ minusEightWeightedShell w n ↔ minusEightNorm v.1 + w * minusEightNorm v.2 = n := by
  constructor
  · exact fun hv => (Finset.mem_filter.mp hv).2
  · intro hv
    have h₁ : minusEightNorm v.1 ≤ n := by omega
    have h₂ : minusEightNorm v.2 ≤ w * minusEightNorm v.2 :=
      Nat.le_mul_of_pos_left _ hw
    exact Finset.mem_filter.mpr ⟨Finset.mem_product.mpr
      ⟨mem_minusEight_vector_box h₁, mem_minusEight_vector_box (by omega)⟩, hv⟩

theorem mem_minusEight_full_shell {n : ℕ} {v : MinusEightShellPoint} :
    v ∈ minusEightFullShell n ↔ minusEightNorm v.1 + minusEightNorm v.2 = n := by
  simpa only [one_mul] using mem_minusEight_weighted_shell (n := n) (v := v) Nat.one_pos

/-- Coordinate divisibility is the local condition being sieved. -/
abbrev MinusEightDivisible (d : ℕ) (v : ℤ × ℤ) : Prop :=
  (d : ℤ) ∣ v.1 ∧ (d : ℤ) ∣ v.2

def minusEightScaleVector (d : ℕ) (v : ℤ × ℤ) : ℤ × ℤ :=
  ((d : ℤ) * v.1, (d : ℤ) * v.2)

theorem minusEight_norm_scale (d : ℕ) (v : ℤ × ℤ) :
    minusEightNorm (minusEightScaleVector d v) = d ^ 2 * minusEightNorm v := by
  have he : (minusEightNorm (minusEightScaleVector d v) : ℤ) =
      ((d ^ 2 * minusEightNorm v : ℕ) : ℤ) := by
    push_cast
    rw [minusEight_norm_cast, minusEight_norm_cast]
    dsimp [minusEightScaleVector]
    ring
  exact_mod_cast he

private def minusEightVectorEmbedding (d : ℕ) (hd : 0 < d) : (ℤ × ℤ) ↪ (ℤ × ℤ) where
  toFun := minusEightScaleVector d
  inj' := by
    intro v w h
    have hd' : (d : ℤ) ≠ 0 := by exact_mod_cast hd.ne'
    exact Prod.ext (mul_left_cancel₀ hd' (congrArg Prod.fst h))
      (mul_left_cancel₀ hd' (congrArg Prod.snd h))

theorem minusEight_divisible_iff_scaled (d : ℕ) (v : ℤ × ℤ) :
    MinusEightDivisible d v ↔ ∃ w, minusEightScaleVector d w = v := by
  constructor
  · rcases v with ⟨x, y⟩
    rintro ⟨⟨x', rfl⟩, ⟨y', rfl⟩⟩
    exact ⟨(x', y'), rfl⟩
  · rintro ⟨w, rfl⟩
    exact ⟨dvd_mul_right _ _, dvd_mul_right _ _⟩

private def minusEightFirstDivisible (n d : ℕ) : Finset MinusEightShellPoint :=
  (minusEightFullShell n).filter (fun v => MinusEightDivisible d v.1)

private def minusEightSecondDivisible (n d : ℕ) : Finset MinusEightShellPoint :=
  (minusEightFullShell n).filter (fun v => MinusEightDivisible d v.2)

/-- The shell with both coordinate blocks divisible by d. -/
def minusEightBothDivisible (n d : ℕ) : Finset MinusEightShellPoint :=
  (minusEightFullShell n).filter (fun v => MinusEightDivisible d v.1 ∧ MinusEightDivisible d v.2)

/-- The shell with either block divisible by d. -/
def minusEightEitherDivisible (n d : ℕ) : Finset MinusEightShellPoint :=
  (minusEightFullShell n).filter (fun v => MinusEightDivisible d v.1 ∨ MinusEightDivisible d v.2)

/-- The G₂ shell retains exactly the points with primitive dyadic blocks. -/
def minusEightDyadicShell (n : ℕ) : Finset MinusEightShellPoint :=
  (minusEightFullShell n).filter (fun v => ¬ (MinusEightDivisible 2 v.1 ∨ MinusEightDivisible 2 v.2))

private theorem minusEight_second_divisible_card (n d : ℕ) (hd : 0 < d) :
    (minusEightSecondDivisible n d).card = (minusEightWeightedShell (d ^ 2) n).card := by
  let e : MinusEightShellPoint ↪ MinusEightShellPoint :=
    Function.Embedding.prodMap (Function.Embedding.refl _) (minusEightVectorEmbedding d hd)
  suffices h : minusEightSecondDivisible n d = (minusEightWeightedShell (d ^ 2) n).map e by
    rw [h, Finset.card_map]
  ext v
  simp only [minusEightSecondDivisible, Finset.mem_filter, mem_minusEight_full_shell,
    Finset.mem_map, mem_minusEight_weighted_shell (pow_pos hd 2)]
  constructor
  · rintro ⟨hv, hdv⟩
    obtain ⟨b, hb⟩ := (minusEight_divisible_iff_scaled d v.2).mp hdv
    have hb' : minusEightNorm v.2 = d ^ 2 * minusEightNorm b := hb ▸ minusEight_norm_scale d b
    exact ⟨(v.1, b), by simpa only [hb'] using hv, Prod.ext rfl hb⟩
  · rintro ⟨v, hv, rfl⟩
    change minusEightNorm v.1 + minusEightNorm (minusEightScaleVector d v.2) = n ∧
      MinusEightDivisible d (minusEightScaleVector d v.2)
    exact ⟨by simpa only [minusEight_norm_scale] using hv,
      (minusEight_divisible_iff_scaled _ _).mpr ⟨_, rfl⟩⟩

private theorem minusEight_first_divisible_card (n d : ℕ) :
    (minusEightFirstDivisible n d).card = (minusEightSecondDivisible n d).card := by
  let e := (Equiv.prodComm (ℤ × ℤ) (ℤ × ℤ)).toEmbedding
  suffices h : minusEightFirstDivisible n d = (minusEightSecondDivisible n d).map e by
    rw [h, Finset.card_map]
  ext v
  simp only [minusEightFirstDivisible, minusEightSecondDivisible, Finset.mem_filter,
    mem_minusEight_full_shell, Finset.mem_map]
  constructor
  · intro ⟨h, hp⟩
    exact ⟨v.swap, ⟨by change minusEightNorm v.2 + minusEightNorm v.1 = n
                       simpa only [add_comm] using h, hp⟩, rfl⟩
  · rintro ⟨w, ⟨h, hp⟩, rfl⟩
    change minusEightNorm w.2 + minusEightNorm w.1 = n ∧ MinusEightDivisible d w.2
    exact ⟨by simpa only [add_comm] using h, hp⟩

/-- Cancelling the same integral dilation in both blocks gives the
smaller full shell, with an empty shell when d² does not divide n. -/
theorem minusEight_both_divisible_card (n d : ℕ) (hd : 0 < d) :
    (minusEightBothDivisible n d).card =
      if d ^ 2 ∣ n then (minusEightFullShell (n / d ^ 2)).card else 0 := by
  let e := Function.Embedding.prodMap (minusEightVectorEmbedding d hd) (minusEightVectorEmbedding d hd)
  suffices h : minusEightBothDivisible n d =
      if d ^ 2 ∣ n then (minusEightFullShell (n / d ^ 2)).map e else ∅ by
    rw [h]
    split <;> simp
  ext v
  simp only [minusEightBothDivisible, Finset.mem_filter, mem_minusEight_full_shell]
  constructor
  · rintro ⟨hv, h₁, h₂⟩
    obtain ⟨a, ha⟩ := (minusEight_divisible_iff_scaled d v.1).mp h₁
    obtain ⟨b, hb⟩ := (minusEight_divisible_iff_scaled d v.2).mp h₂
    have hn : d ^ 2 * (minusEightNorm a + minusEightNorm b) = n := by
      rw [Nat.mul_add, ← minusEight_norm_scale, ← minusEight_norm_scale, ha, hb, hv]
    rw [if_pos ⟨_, hn.symm⟩]
    apply Finset.mem_map.mpr
    exact ⟨(a, b), mem_minusEight_full_shell.mpr
      (by rw [← hn, Nat.mul_div_cancel_left _ (pow_pos hd 2)]), Prod.ext ha hb⟩
  · by_cases hdiv : d ^ 2 ∣ n
    · rw [if_pos hdiv]
      intro hv
      obtain ⟨w, hw, rfl⟩ := Finset.mem_map.mp hv
      refine ⟨?_, (minusEight_divisible_iff_scaled _ _).mpr ⟨_, rfl⟩,
        (minusEight_divisible_iff_scaled _ _).mpr ⟨_, rfl⟩⟩
      change minusEightNorm (minusEightScaleVector d w.1) +
        minusEightNorm (minusEightScaleVector d w.2) = n
      rw [minusEight_norm_scale, minusEight_norm_scale, ← Nat.mul_add,
        mem_minusEight_full_shell.mp hw, Nat.mul_div_cancel' hdiv]
    · simp [hdiv]

/-- Inclusion-exclusion on the two even coordinate blocks. -/
theorem minusEight_dyadic_shell_card (n : ℕ) :
    (minusEightFullShell n).card + (if 4 ∣ n then (minusEightFullShell (n / 4)).card else 0) =
      (minusEightDyadicShell n).card + 2 * (minusEightWeightedShell 4 n).card := by
  have h := Finset.card_filter_add_card_filter_not
    (s := minusEightFullShell n) (fun v => MinusEightDivisible 2 v.1 ∨ MinusEightDivisible 2 v.2)
  change (minusEightEitherDivisible n 2).card + (minusEightDyadicShell n).card =
    (minusEightFullShell n).card at h
  have h' := Finset.card_union_add_card_inter (minusEightFirstDivisible n 2) (minusEightSecondDivisible n 2)
  unfold minusEightFirstDivisible minusEightSecondDivisible at h'
  rw [← Finset.filter_or, ← Finset.filter_and] at h'
  change (minusEightEitherDivisible n 2).card + (minusEightBothDivisible n 2).card =
    (minusEightFirstDivisible n 2).card + (minusEightSecondDivisible n 2).card at h'
  rw [minusEight_first_divisible_card, minusEight_second_divisible_card n 2 (by decide),
    minusEight_both_divisible_card n 2 (by decide)] at h'
  change _ + (if 4 ∣ n then (minusEightFullShell (n / 4)).card else 0) =
    (minusEightWeightedShell 4 n).card + (minusEightWeightedShell 4 n).card at h'
  omega

/-- The inert primes for x² + 2y², all of which are odd. -/
def MinusEightInertPrime (p : ℕ) : Prop :=
  p.Prime ∧ (p % 8 = 5 ∨ p % 8 = 7)

instance : DecidablePred MinusEightInertPrime :=
  fun _ => inferInstanceAs (Decidable (_ ∧ _))

/-- At an inert prime a zero of the norm is a zero of both coordinates. -/
theorem minusEight_inert_dvd_norm_iff {p : ℕ} (hp : MinusEightInertPrime p) (v : ℤ × ℤ) :
    p ∣ minusEightNorm v ↔ MinusEightDivisible p v := by
  letI : Fact p.Prime := ⟨hp.1⟩
  let x : ZMod p := v.1
  let y : ZMod p := v.2
  have hv : (minusEightNorm v : ZMod p) = x ^ 2 + 2 * y ^ 2 := by
    rw [← Int.cast_natCast, minusEight_norm_cast]
    simp [x, y]
  have hp₂ : p ≠ 2 := by rcases hp.2 with h | h <;> omega
  have hns : ¬ IsSquare (-2 : ZMod p) := by
    rw [ZMod.exists_sq_eq_neg_two_iff hp₂]
    rcases hp.2 with h | h <;> omega
  change p ∣ minusEightNorm v ↔ (p : ℤ) ∣ v.1 ∧ (p : ℤ) ∣ v.2
  rw [← ZMod.natCast_eq_zero_iff, hv, ← ZMod.intCast_zmod_eq_zero_iff_dvd,
    ← ZMod.intCast_zmod_eq_zero_iff_dvd]
  change x ^ 2 + 2 * y ^ 2 = 0 ↔ x = 0 ∧ y = 0
  constructor
  · intro h
    by_cases hy : y = 0
    · exact ⟨by simpa [hy] using h, hy⟩
    · have hr : (x / y) ^ 2 = -2 := by
        field_simp [hy]
        linear_combination h
      exact (hns ⟨x / y, by rw [← pow_two, hr]⟩).elim
  · rintro ⟨hx, hy⟩
    simp [hx, hy]

/-- At an inert prime dividing the target, a bad block forces both blocks
to be divisible by the prime. The smaller-shell formula retains the full valuation. -/
theorem minusEight_inert_target_bad_card {n p : ℕ} (hp : MinusEightInertPrime p) (hpn : p ∣ n) :
    (minusEightEitherDivisible n p).card =
      if p ^ 2 ∣ n then (minusEightFullShell (n / p ^ 2)).card else 0 := by
  rw [← minusEight_both_divisible_card n p hp.1.pos]
  congr 1
  ext v
  simp only [minusEightEitherDivisible, minusEightBothDivisible, Finset.mem_filter,
    ← minusEight_inert_dvd_norm_iff hp]
  constructor
  · rintro ⟨hv, h | h⟩
    · exact ⟨hv, h, (Nat.dvd_add_iff_right h).mpr (mem_minusEight_full_shell.mp hv ▸ hpn)⟩
    · exact ⟨hv, (Nat.dvd_add_iff_left h).mpr (mem_minusEight_full_shell.mp hv ▸ hpn), h⟩
  · rintro ⟨hv, h, _⟩
    exact ⟨hv, Or.inl h⟩

private theorem minusEight_val_even_iff {q : ℕ} [NeZero q] (hq : 2 ∣ q) (z : ℤ) :
    (z : ZMod q).val % 2 = 0 ↔ (2 : ℤ) ∣ z := by
  have hq' : (2 : ℤ) ∣ (q : ℤ) := by exact_mod_cast hq
  rw [Int.dvd_iff_emod_eq_zero, ← Int.emod_emod_of_dvd z hq', ← ZMod.val_intCast (n := q)]
  norm_cast

/-- A norm is a multiple of four exactly when both coordinates are even. -/
theorem minusEight_four_dvd_norm_iff (v : ℤ × ℤ) :
    4 ∣ minusEightNorm v ↔ MinusEightDivisible 2 v := by
  have hcases : ∀ x y : ZMod 4, x ^ 2 + 2 * y ^ 2 = 0 ↔ x.val % 2 = 0 ∧ y.val % 2 = 0 := by decide
  have hv : (minusEightNorm v : ZMod 4) = (v.1 : ZMod 4) ^ 2 + 2 * (v.2 : ZMod 4) ^ 2 := by
    rw [← Int.cast_natCast, minusEight_norm_cast]
    push_cast
    rfl
  change 4 ∣ minusEightNorm v ↔ (2 : ℤ) ∣ v.1 ∧ (2 : ℤ) ∣ v.2
  simp only [← ZMod.natCast_eq_zero_iff, hv, hcases, minusEight_val_even_iff (by decide : 2 ∣ 4)]

theorem minusEight_norm_pos_of_not_divisible_two (v : ℤ × ℤ) (hv : ¬ MinusEightDivisible 2 v) :
    0 < minusEightNorm v := by
  by_contra! hn
  have hz : minusEightNorm v = 0 := by omega
  exact hv ((minusEight_four_dvd_norm_iff v).mp (hz ▸ dvd_zero 4))

private theorem minusEight_value_of_avoids (v : ℤ × ℤ) (h₂ : ¬ MinusEightDivisible 2 v)
    (hinert : ∀ p, MinusEightInertPrime p → ¬ MinusEightDivisible p v) :
    P8 (minusEightNorm v) := by
  let n := minusEightNorm v
  have hn : 0 < n := minusEight_norm_pos_of_not_divisible_two v h₂
  have h₄ : ¬ 4 ∣ n := mt (minusEight_four_dvd_norm_iff v).mp h₂
  have hno (p : ℕ) (hp : MinusEightInertPrime p) : ¬ p ∣ n :=
    mt (minusEight_inert_dvd_norm_iff hp v).mp (hinert p hp)
  have hcore {a : ℕ} (ha : 0 < a) (ha₂ : ¬ 2 ∣ a) (hd : a ∣ n) : C8 a := by
    refine ⟨ha, ?_⟩
    intro p hp hpa
    rcases hp.eq_two_or_odd with rfl | hp₂
    · exact (ha₂ hpa).elim
    · have hi : ¬ (p % 8 = 5 ∨ p % 8 = 7) :=
        fun h => hno p ⟨hp, h⟩ (hpa.trans hd)
      have hmod := Nat.mod_mod_of_dvd p (by decide : 2 ∣ 8)
      omega
  apply (P8_iff_core _).mpr
  by_cases hdn : 2 ∣ n
  · refine Or.inr ⟨n / 2, hcore (Nat.div_pos (Nat.le_of_dvd hn hdn) (by decide)) ?_
      (Nat.div_dvd_of_dvd hdn), Nat.mul_div_cancel' hdn⟩
    intro hd
    exact h₄ ((Nat.dvd_div_iff_mul_dvd hdn).mp hd)
  · exact Or.inl (hcore hn hdn dvd_rfl)

/-- A surviving G₂-shell point gives two primitive numerical values. -/
theorem minusEight_values_of_unspoiled {n : ℕ} {v : MinusEightShellPoint}
    (hv : v ∈ minusEightDyadicShell n) (h : ∀ p ≤ n, MinusEightInertPrime p →
      v ∉ minusEightEitherDivisible n p) : P8 (minusEightNorm v.1) ∧ P8 (minusEightNorm v.2) := by
  obtain ⟨hfull, h₂⟩ := Finset.mem_filter.mp hv
  have hs := mem_minusEight_full_shell.mp hfull
  obtain ⟨h₁, h₂⟩ := not_or.mp h₂
  have hainert (p : ℕ) (hp : MinusEightInertPrime p) : ¬ MinusEightDivisible p v.1 := by
    intro hdiv
    have hle := Nat.le_of_dvd (minusEight_norm_pos_of_not_divisible_two v.1 h₁) ((minusEight_inert_dvd_norm_iff hp _).mpr hdiv)
    exact h p (hle.trans (by omega)) hp (Finset.mem_filter.mpr ⟨hfull, Or.inl hdiv⟩)
  have hbinert (p : ℕ) (hp : MinusEightInertPrime p) : ¬ MinusEightDivisible p v.2 := by
    intro hdiv
    have hle := Nat.le_of_dvd (minusEight_norm_pos_of_not_divisible_two v.2 h₂) ((minusEight_inert_dvd_norm_iff hp _).mpr hdiv)
    exact h p (hle.trans (by omega)) hp (Finset.mem_filter.mpr ⟨hfull, Or.inr hdiv⟩)
  exact ⟨minusEight_value_of_avoids v.1 h₁ hainert, minusEight_value_of_avoids v.2 h₂ hbinert⟩

theorem minusEight_sum_of_unspoiled {n : ℕ} {v : MinusEightShellPoint}
    (hv : v ∈ minusEightDyadicShell n) (h : ∀ p ≤ n, MinusEightInertPrime p →
      v ∉ minusEightEitherDivisible n p) : InSumset P8 P8 n :=
  ⟨_, _, (minusEight_values_of_unspoiled hv h).1, (minusEight_values_of_unspoiled hv h).2,
    mem_minusEight_full_shell.mp (Finset.mem_filter.mp hv).1⟩

private theorem minusEight_norm_mod_sixteen (v : ℤ × ℤ) (hv : ¬ MinusEightDivisible 2 v) :
    minusEightNorm v % 16 ∈ [1, 2, 3, 6, 9, 11] := by
  have hcases : ∀ x y : ZMod 16,
      ¬ (x.val % 2 = 0 ∧ y.val % 2 = 0) →
        (x ^ 2 + 2 * y ^ 2).val ∈ [1, 2, 3, 6, 9, 11] := by decide +kernel
  have he : (minusEightNorm v : ZMod 16) = (v.1 : ZMod 16) ^ 2 + 2 * (v.2 : ZMod 16) ^ 2 := by
    rw [← Int.cast_natCast, minusEight_norm_cast]
    push_cast
    rfl
  rw [← ZMod.val_natCast, he]
  exact hcases _ _ (by
    simpa only [minusEight_val_even_iff (by decide : 2 ∣ 16), MinusEightDivisible,
      Nat.cast_ofNat] using hv)

/-- The dyadic shell is empty at all targets divisible by sixteen. -/
theorem minusEight_dyadic_shell_eq_empty {n : ℕ} (hn : 16 ∣ n) :
    minusEightDyadicShell n = ∅ := by
  apply Finset.filter_eq_empty_iff.mpr
  intro v hv
  rw [not_not]
  by_contra h
  obtain ⟨h₁, h₂⟩ := not_or.mp h
  have ha := minusEight_norm_mod_sixteen v.1 h₁
  have hb := minusEight_norm_mod_sixteen v.2 h₂
  have hsum := mem_minusEight_full_shell.mp hv
  have hmod := Nat.add_mod (minusEightNorm v.1) (minusEightNorm v.2) 16
  rw [hsum, Nat.dvd_iff_mod_eq_zero.mp hn] at hmod
  simp only [List.mem_cons, List.not_mem_nil, or_false] at ha hb
  omega

/-- The smaller selected shell used for even core targets. -/
def minusEightOddOddShell (n : ℕ) : Finset MinusEightShellPoint :=
  (minusEightFullShell n).filter (fun v => ¬ 2 ∣ minusEightNorm v.1 ∧ ¬ 2 ∣ minusEightNorm v.2)

theorem minusEight_even_norm_iff (v : ℤ × ℤ) :
    2 ∣ minusEightNorm v ↔ (2 : ℤ) ∣ v.1 := by
  have hcases : ∀ x : ZMod 2, x ^ 2 = 0 ↔ x = 0 := by decide
  have hv : (minusEightNorm v : ZMod 2) = (v.1 : ZMod 2) ^ 2 := by
    rw [← Int.cast_natCast, minusEight_norm_cast]
    simp [show (2 : ZMod 2) = 0 by decide]
  rw [← ZMod.natCast_eq_zero_iff, hv, hcases]
  exact ZMod.intCast_zmod_eq_zero_iff_dvd v.1 2

private def minusEightRamifiedMap (v : ℤ × ℤ) : ℤ × ℤ := (-2 * v.2, v.1)

private def minusEightRamifiedEmbedding : (ℤ × ℤ) ↪ (ℤ × ℤ) where
  toFun := minusEightRamifiedMap
  inj' := by
    intro a b h
    have h₁ : (-2 : ℤ) * a.2 = (-2 : ℤ) * b.2 := congrArg Prod.fst h
    exact Prod.ext (congrArg Prod.snd h) (mul_left_cancel₀ (by decide : (-2 : ℤ) ≠ 0) h₁)

private theorem minusEight_ramified_norm (v : ℤ × ℤ) :
    minusEightNorm (minusEightRamifiedMap v) = 2 * minusEightNorm v := by
  have h : (minusEightNorm (minusEightRamifiedMap v) : ℤ) = 2 * (minusEightNorm v : ℤ) := by
    rw [minusEight_norm_cast, minusEight_norm_cast]
    dsimp [minusEightRamifiedMap]
    ring
  exact_mod_cast h

theorem minusEight_even_norm_iff_ramified (v : ℤ × ℤ) :
    2 ∣ minusEightNorm v ↔ ∃ a, minusEightRamifiedMap a = v := by
  rw [minusEight_even_norm_iff]
  constructor
  · intro hv
    obtain ⟨x, hx⟩ := hv
    exact ⟨(v.2, -x), Prod.ext (by dsimp [minusEightRamifiedMap]; rw [hx]; ring) rfl⟩
  · rintro ⟨a, rfl⟩
    exact ⟨-a.2, by dsimp [minusEightRamifiedMap]; ring⟩

/-- Multiplying both blocks by √−2 accounts for the entire even-even
part of an even target shell. Its complement consists of two odd norms. -/
theorem minusEight_odd_odd_shell_card (n : ℕ) :
    (minusEightFullShell (2 * n)).card =
      (minusEightOddOddShell (2 * n)).card + (minusEightFullShell n).card := by
  let e := Function.Embedding.prodMap minusEightRamifiedEmbedding minusEightRamifiedEmbedding
  have heq : ((minusEightFullShell (2 * n)).filter (fun v =>
      ¬ (¬ 2 ∣ minusEightNorm v.1 ∧ ¬ 2 ∣ minusEightNorm v.2))) =
        (minusEightFullShell n).map e := by
    ext v
    simp only [Finset.mem_filter, Finset.mem_map, mem_minusEight_full_shell, not_and_or, not_not]
    constructor
    · rintro ⟨h, h₁ | h₂⟩
      all_goals
        have hn : 2 ∣ minusEightNorm v.1 + minusEightNorm v.2 := h ▸ dvd_mul_right 2 n
        have hboth : (2 ∣ minusEightNorm v.1) ∧ (2 ∣ minusEightNorm v.2) := by
          first
          | exact ⟨h₁, (Nat.dvd_add_iff_right h₁).mpr hn⟩
          | exact ⟨(Nat.dvd_add_iff_left h₂).mpr hn, h₂⟩
        obtain ⟨x, hx⟩ := (minusEight_even_norm_iff_ramified _).mp hboth.1
        obtain ⟨y, hy⟩ := (minusEight_even_norm_iff_ramified _).mp hboth.2
        rw [← hx, ← hy, minusEight_ramified_norm, minusEight_ramified_norm, ← Nat.mul_add] at h
        exact ⟨(x, y), Nat.eq_of_mul_eq_mul_left (by decide) h, Prod.ext hx hy⟩
    · rintro ⟨w, hw, rfl⟩
      change minusEightNorm (minusEightRamifiedMap w.1) +
        minusEightNorm (minusEightRamifiedMap w.2) = 2 * n ∧ _
      rw [minusEight_ramified_norm, minusEight_ramified_norm, ← Nat.mul_add, hw]
      exact ⟨rfl, Or.inl ((minusEight_even_norm_iff_ramified _).mpr ⟨_, rfl⟩)⟩
  have h := Finset.card_filter_add_card_filter_not (s := minusEightFullShell (2 * n))
    (fun v => ¬ 2 ∣ minusEightNorm v.1 ∧ ¬ 2 ∣ minusEightNorm v.2)
  rw [heq, Finset.card_map] at h
  exact h.symm

theorem minusEight_odd_odd_shell_subset (n : ℕ) :
    minusEightOddOddShell n ⊆ minusEightDyadicShell n := by
  intro v hv
  obtain ⟨hs, h₁, h₂⟩ := Finset.mem_filter.mp hv
  exact Finset.mem_filter.mpr ⟨hs, not_or.mpr
    ⟨fun h => h₁ ((show 2 ∣ 4 by decide).trans ((minusEight_four_dvd_norm_iff _).mpr h)),
     fun h => h₂ ((show 2 ∣ 4 by decide).trans ((minusEight_four_dvd_norm_iff _).mpr h))⟩⟩

/-- Odd norm blocks that survive the sieve give the ordinary core pair. -/
theorem minusEight_core_sum_of_unspoiled_odd_odd {n : ℕ} {v : MinusEightShellPoint}
    (hv : v ∈ minusEightOddOddShell n) (h : ∀ p ≤ n, MinusEightInertPrime p →
      v ∉ minusEightEitherDivisible n p) : InSumset C8 C8 n := by
  obtain ⟨hfull, ha, hb⟩ := Finset.mem_filter.mp hv
  obtain ⟨h₁, h₂⟩ := minusEight_values_of_unspoiled (minusEight_odd_odd_shell_subset n hv) h
  have hcore {a : ℕ} (ha : ¬ 2 ∣ a) (hP : P8 a) : C8 a :=
    ((P8_iff_core a).mp hP).resolve_right (fun ⟨_, _, he⟩ => ha ⟨_, he.symm⟩)
  exact ⟨_, _, hcore ha h₁, hcore hb h₂, mem_minusEight_full_shell.mp hfull⟩

end QuadraticNorms
