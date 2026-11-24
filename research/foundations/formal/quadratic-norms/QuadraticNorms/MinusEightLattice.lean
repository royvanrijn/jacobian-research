import QuadraticNorms.MinusEightShell
import Mathlib.Analysis.Real.Sqrt
import Mathlib.Data.Int.Interval

open scoped BigOperators

namespace QuadraticNorms

/-- Nonzero signed vectors of norm at most n. -/
def minusEightLatticeBall (n : ℕ) : Finset (ℤ × ℤ) :=
  (minusEightVectorBox n).filter fun v => 0 < minusEightNorm v ∧ minusEightNorm v ≤ n

theorem mem_minusEight_lattice_ball {n : ℕ} {v : ℤ × ℤ} :
    v ∈ minusEightLatticeBall n ↔ 0 < minusEightNorm v ∧ minusEightNorm v ≤ n :=
  ⟨fun h => (Finset.mem_filter.mp h).2,
    fun h => Finset.mem_filter.mpr ⟨mem_minusEight_vector_box h.2, h⟩⟩

/-- Signed vectors of odd norm at most n. -/
def minusEightOddLatticeBall (n : ℕ) : Finset (ℤ × ℤ) :=
  (minusEightLatticeBall n).filter fun v => ¬ 2 ∣ minusEightNorm v

/-- All signed representations of the numerical norm n. -/
def minusEightNormFiber (n : ℕ) : Finset (ℤ × ℤ) :=
  (minusEightVectorBox n).filter fun v => minusEightNorm v = n

theorem mem_minusEight_norm_fiber {n : ℕ} {v : ℤ × ℤ} :
    v ∈ minusEightNormFiber n ↔ minusEightNorm v = n :=
  ⟨fun h => (Finset.mem_filter.mp h).2,
    fun h => Finset.mem_filter.mpr ⟨mem_minusEight_vector_box h.le, h⟩⟩

private theorem minusEight_norm_real (v : ℤ × ℤ) :
    (minusEightNorm v : ℝ) = (v.1 : ℝ) ^ 2 + 2 * (v.2 : ℝ) ^ 2 := by
  rw [← Int.cast_natCast, minusEight_norm_cast]
  push_cast
  rfl

private theorem int_card_le_of_bounds (s : Finset ℤ) {lo hi : ℝ} (hle : lo ≤ hi)
    (hs : ∀ x ∈ s, lo ≤ (x : ℝ) ∧ (x : ℝ) ≤ hi) : (s.card : ℝ) ≤ hi - lo + 1 := by
  by_cases hne : s.Nonempty
  · let a := s.min' hne
    let b := s.max' hne
    have ha : a ∈ s := Finset.min'_mem ..
    have hb : b ∈ s := Finset.max'_mem ..
    have hab : a ≤ b := Finset.min'_le _ _ hb
    have hsub : s ⊆ Finset.Icc a b := fun _ h =>
      Finset.mem_Icc.mpr ⟨Finset.min'_le _ _ h, Finset.le_max' _ _ h⟩
    have hc : (s.card : ℤ) ≤ ((Finset.Icc a b).card : ℤ) := by exact_mod_cast Finset.card_le_card hsub
    rw [Int.card_Icc_of_le _ _ (by omega)] at hc
    have hcr : (s.card : ℝ) ≤ (b : ℝ) + 1 - a := by exact_mod_cast hc
    linarith only [hcr, (hs a ha).1, (hs b hb).2]
  · rw [Finset.not_nonempty_iff_eq_empty.mp hne, Finset.card_empty, Nat.cast_zero]
    linarith

private theorem signed_int_quad_sum (u v : ℝ) (m : ℕ) :
    ∑ y ∈ Finset.Icc (-(m : ℤ)) m, (u - v * (y : ℝ) ^ 2) =
      (2 * (m : ℝ) + 1) * u - (v / 3) * (m : ℝ) * (m + 1) * (2 * m + 1) := by
  induction m with
  | zero => simp
  | succ m ih =>
    have heq : Finset.Icc (-((m + 1 : ℕ) : ℤ)) ((m + 1 : ℕ) : ℤ) =
        insert (-((m + 1 : ℕ) : ℤ)) (insert ((m + 1 : ℕ) : ℤ) (Finset.Icc (-(m : ℤ)) m)) := by
      ext y
      simp only [Finset.mem_Icc, Finset.mem_insert, Nat.cast_add, Nat.cast_one]
      omega
    have hlo : -((m + 1 : ℕ) : ℤ) ∉
        insert ((m + 1 : ℕ) : ℤ) (Finset.Icc (-(m : ℤ)) m) := by
      simp only [Finset.mem_Icc, Finset.mem_insert, Nat.cast_add, Nat.cast_one]
      omega
    have hhi : ((m + 1 : ℕ) : ℤ) ∉ Finset.Icc (-(m : ℤ)) m := by simp [Finset.mem_Icc]
    rw [heq, Finset.sum_insert hlo, Finset.sum_insert hhi, ih]
    push_cast
    ring

private theorem minusEight_lattice_row_bound {s : Finset (ℤ × ℤ)} {n : ℕ} (hn : 0 < n)
    (hs : ∀ v ∈ s, minusEightNorm v ≤ n)
    {d : ℕ} (hd : 0 < d) (r y : ℤ)
    (hres : ∀ v ∈ s, (d : ℤ) ∣ v.1 - r) (hy : 2 * (y : ℝ) ^ 2 ≤ n) :
    (((s.filter fun v => v.2 = y).card : ℝ)) ≤
      (2 * Real.sqrt n / d + 1) - (2 / ((d : ℝ) * Real.sqrt n)) * (y : ℝ) ^ 2 := by
  let R := Real.sqrt (n : ℝ)
  let W := Real.sqrt ((n : ℝ) - 2 * (y : ℝ) ^ 2)
  let row := s.filter (fun v => v.2 = y)
  let idx (v : ℤ × ℤ) := (v.1 - r) / (d : ℤ)
  have hidx (v : ℤ × ℤ) (hv : v ∈ row) : (d : ℤ) * idx v + r = v.1 := by
    dsimp [idx]
    rw [Int.mul_ediv_cancel' (hres v (Finset.mem_filter.mp hv).1), sub_add_cancel]
  have hi : Set.InjOn idx (row : Set (ℤ × ℤ)) := by
    intro a ha b hb h
    apply Prod.ext
    · rw [← hidx a ha, ← hidx b hb, h]
    · exact (Finset.mem_filter.mp ha).2.trans (Finset.mem_filter.mp hb).2.symm
  have hdR : 0 < (d : ℝ) := by exact_mod_cast hd
  have hR : 0 < R := Real.sqrt_pos.2 (by exact_mod_cast hn)
  have hRsq : R ^ 2 = (n : ℝ) := Real.sq_sqrt (by positivity)
  have hW : 0 ≤ W := Real.sqrt_nonneg _
  have hWsq : W ^ 2 = (n : ℝ) - 2 * (y : ℝ) ^ 2 := Real.sq_sqrt (sub_nonneg.mpr hy)
  have habs (v : ℤ × ℤ) (hv : v ∈ row) : |(v.1 : ℝ)| ≤ W := by
    obtain ⟨hv', h₂⟩ := Finset.mem_filter.mp hv
    have hN : (minusEightNorm v : ℝ) ≤ n := by exact_mod_cast hs v hv'
    rw [minusEight_norm_real, h₂] at hN
    apply (sq_le_sq₀ (abs_nonneg _) hW).mp
    rw [sq_abs, hWsq]
    linarith only [hN]
  have hwidth : W ≤ R - (y : ℝ) ^ 2 / R := by
    have hb : (y : ℝ) ^ 2 / R ≤ R - W :=
      (div_le_iff₀ hR).mpr (by nlinarith only [hRsq, hWsq, sq_nonneg (R - W)])
    linarith
  have hcard : ((row.image idx).card : ℝ) ≤ ((W - r) / d - ((-W - r) / d) + 1) := by
    apply int_card_le_of_bounds _ ((div_le_div_iff_of_pos_right hdR).mpr (by linarith))
    intro x hx
    obtain ⟨v, hv, rfl⟩ := Finset.mem_image.mp hx
    have he : (d : ℝ) * (idx v : ℝ) + (r : ℝ) = (v.1 : ℝ) := by exact_mod_cast hidx v hv
    obtain ⟨ha, hb⟩ := abs_le.mp (habs v hv)
    exact ⟨(div_le_iff₀ hdR).mpr (by nlinarith), (le_div_iff₀ hdR).mpr (by nlinarith)⟩
  rw [Finset.card_image_of_injOn hi] at hcard
  change (row.card : ℝ) ≤ 2 * R / d + 1 - _
  calc
    _ ≤ 2 * W / d + 1 := hcard.trans_eq (by ring)
    _ ≤ _ := by
      rw [show 2 * R / (d : ℝ) + 1 - (2 / ((d : ℝ) * R)) * (y : ℝ) ^ 2 =
        (2 * (R - (y : ℝ) ^ 2 / R)) / d + 1 by field_simp; ring]
      exact add_le_add (div_le_div_of_nonneg_right
        (mul_le_mul_of_nonneg_left hwidth (by norm_num : 0 ≤ (2 : ℝ))) hdR.le) le_rfl

private theorem minusEight_lattice_bound {s : Finset (ℤ × ℤ)} {n : ℕ} (hn : 0 < n)
    (hs : ∀ v ∈ s, minusEightNorm v ≤ n)
    {d : ℕ} (hd : 0 < d) (r : ℤ) (hres : ∀ v ∈ s, (d : ℤ) ∣ v.1 - r) :
    (s.card : ℝ) ≤ (50 / (21 * (d : ℝ))) * n + (2 / (d : ℝ) + 10 / 7) * Real.sqrt n + 1 := by
  let R := Real.sqrt (n : ℝ)
  let B := R / Real.sqrt 2
  let m := ⌊B⌋₊
  let ys := Finset.Icc (-(m : ℤ)) m
  let u := 2 * R / (d : ℝ) + 1
  let v := 2 / ((d : ℝ) * R)
  let F (t : ℝ) := (2 * t + 1) * u - (2 * v / 3) * t ^ 3
  have hdR : 0 < (d : ℝ) := by exact_mod_cast hd
  have hR : 0 < R := Real.sqrt_pos.2 (by exact_mod_cast hn)
  have hRsq : R ^ 2 = (n : ℝ) := Real.sq_sqrt (by positivity)
  have hB₀ : 0 ≤ B := by dsimp [B]; positivity
  have hm : (m : ℝ) ≤ B := Nat.floor_le hB₀
  have hBsq : 2 * B ^ 2 = n := by
    dsimp [B]
    rw [div_pow, Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
    linarith only [hRsq]
  have hmap : Set.MapsTo (Prod.snd : ℤ × ℤ → ℤ) ↑s ↑ys := by
    intro w hw
    have hN : (minusEightNorm w : ℝ) ≤ n := by exact_mod_cast hs w hw
    rw [minusEight_norm_real] at hN
    have hbound : |(w.2 : ℝ)| ≤ B := by
      apply (sq_le_sq₀ (abs_nonneg _) hB₀).mp
      rw [sq_abs]
      nlinarith only [hN, hBsq, sq_nonneg (w.1 : ℝ)]
    have habs : |w.2| ≤ (m : ℤ) := by
      have hnat : (w.2.natAbs : ℝ) ≤ B := by
        rw [← Int.cast_natCast, Int.natCast_natAbs, Int.cast_abs]
        exact hbound
      rw [← Int.natCast_natAbs]
      exact_mod_cast (Nat.le_floor hnat)
    exact Finset.mem_Icc.mpr (abs_le.mp habs)
  have hrow (y : ℤ) (hy : y ∈ ys) :
      (((s.filter fun w => w.2 = y).card : ℝ)) ≤ u - v * (y : ℝ) ^ 2 := by
    obtain ⟨hy₁, hy₂⟩ := Finset.mem_Icc.mp hy
    have hyabs : |(y : ℝ)| ≤ (m : ℝ) :=
      abs_le.mpr ⟨by exact_mod_cast hy₁, by exact_mod_cast hy₂⟩
    have hySq := pow_le_pow_left₀ (abs_nonneg (y : ℝ)) (hyabs.trans hm) 2
    rw [sq_abs] at hySq
    exact minusEight_lattice_row_bound hn hs hd r y hres (by linarith only [hySq, hBsq])
  have hF : (∑ y ∈ ys, (u - v * (y : ℝ) ^ 2)) ≤ F m := by
    dsimp only [ys]
    rw [signed_int_quad_sum]
    have hm₀ : (0 : ℝ) ≤ m := by positivity
    have hpol : 2 * (m : ℝ) ^ 3 ≤ (m : ℝ) * (m + 1) * (2 * m + 1) := by
      nlinarith only [sq_nonneg (m : ℝ), hm₀]
    have hh := mul_le_mul_of_nonneg_left hpol (show 0 ≤ v / 3 by dsimp [v]; positivity)
    dsimp [F]
    nlinarith only [hh]
  have hscale : (2 * v / 3) * (3 * B ^ 2) = 2 * R / (d : ℝ) := by
    dsimp [v]
    field_simp
    linarith only [hBsq, hRsq]
  have hmono : F m ≤ F B := by
    have hq : B ^ 2 + B * (m : ℝ) + (m : ℝ) ^ 2 ≤ 3 * B ^ 2 := by
      nlinarith only [mul_le_mul_of_nonneg_left hm hB₀,
        pow_le_pow_left₀ (show (0 : ℝ) ≤ m by positivity) hm 2]
    have hv : 0 ≤ 2 * v / 3 := by dsimp [v]; positivity
    have hh := mul_le_mul_of_nonneg_left hq hv
    rw [hscale] at hh
    have hpos : 0 ≤ 2 * u - (2 * v / 3) * (B ^ 2 + B * m + (m : ℝ) ^ 2) := by
      apply sub_nonneg.mpr (hh.trans ?_)
      have h' : 0 ≤ 2 * R / (d : ℝ) := div_nonneg (by positivity) hdR.le
      dsimp only [u]
      linarith only [h']
    have hh' := mul_nonneg (sub_nonneg.mpr hm) hpos
    dsimp [F]
    nlinarith only [hh']
  have hB_le : B ≤ (5 / 7) * R := by
    have hsqrt : (7 / 5 : ℝ) ≤ Real.sqrt 2 :=
      (Real.le_sqrt (by norm_num) (by norm_num)).mpr (by norm_num)
    dsimp [B]
    apply (div_le_iff₀ (Real.sqrt_pos.2 (by norm_num))).mpr
    nlinarith only [mul_le_mul_of_nonneg_right hsqrt hR.le]
  have hF_eq : F B = 2 * R / d + 1 + ((10 / 3) * R / (d : ℝ) + 2) * B := by
    have h : (2 * v / 3) * B ^ 3 = (2 / 3) * R / (d : ℝ) * B := by
      calc
        _ = ((2 * v / 3) * (3 * B ^ 2)) * B / 3 := by ring
        _ = _ := by rw [hscale]; ring
    dsimp only [F, u]
    rw [h]
    ring
  calc
    _ = ∑ y ∈ ys, (((s.filter fun w => w.2 = y).card : ℝ)) := by
      exact_mod_cast Finset.card_eq_sum_card_fiberwise hmap
    _ ≤ ∑ y ∈ ys, (u - v * (y : ℝ) ^ 2) := Finset.sum_le_sum hrow
    _ ≤ F m := hF
    _ ≤ F B := hmono
    _ = _ := hF_eq
    _ ≤ _ := by
      have h := mul_le_mul_of_nonneg_left hB_le
        (show 0 ≤ (10 / 3 : ℝ) * R / (d : ℝ) + 2 by positivity)
      have he : ((10 / 3 : ℝ) * R / (d : ℝ) + 2) * ((5 / 7) * R) =
          (50 / (21 * (d : ℝ))) * n + (10 / 7) * R := by rw [← hRsq]; ring
      rw [he] at h
      change 2 * R / (d : ℝ) + 1 + _ ≤ _ + (2 / (d : ℝ) + 10 / 7) * R + 1
      calc
        _ ≤ 2 * R / (d : ℝ) + 1 + ((50 / (21 * (d : ℝ))) * n + (10 / 7) * R) :=
          add_le_add (le_refl _) h
        _ = _ := by ring

set_option maxRecDepth 16000 in
set_option maxHeartbeats 3000000 in
private theorem minusEight_lattice_ball_small :
    ∀ j : Fin 169, 3 * (minusEightLatticeBall j.val).card ≤ 8 * j.val := by
  decide +kernel

/-- The sharp positive lattice count used for the odd selected shell. -/
theorem minusEight_lattice_ball_card_le (n : ℕ) :
    ((minusEightLatticeBall n).card : ℝ) ≤ (8 / 3) * n := by
  by_cases hn : 169 ≤ n
  · have h := minusEight_lattice_bound (s := minusEightLatticeBall n) (by omega : 0 < n)
      (fun _ hv => (mem_minusEight_lattice_ball.mp hv).2)
      (d := 1) (by decide) 0 (by intro v _; simp)
    have hR : (13 : ℝ) ≤ Real.sqrt n :=
      (Real.le_sqrt (by norm_num) (by positivity)).mpr (by exact_mod_cast hn)
    have hs := Real.sq_sqrt (show 0 ≤ (n : ℝ) by positivity)
    have hlin := mul_le_mul_of_nonneg_right hR (Real.sqrt_nonneg (n : ℝ))
    norm_num at h
    linarith only [h, hs, hlin, show (169 : ℝ) ≤ n by exact_mod_cast hn]
  · have h := minusEight_lattice_ball_small ⟨n, by omega⟩
    have hR : 3 * ((minusEightLatticeBall n).card : ℝ) ≤ 8 * n := by exact_mod_cast h
    linarith only [hR]

set_option maxRecDepth 8000 in
private theorem minusEight_odd_ball_small :
    ∀ j : Fin 25, (minusEightOddLatticeBall j.val).card ≤ 2 * j.val := by
  decide +kernel

/-- The sharp odd lattice count used for the two even selected shells. -/
theorem minusEight_odd_lattice_ball_card_le (n : ℕ) :
    ((minusEightOddLatticeBall n).card : ℝ) ≤ 2 * n := by
  by_cases hn : 25 ≤ n
  · have hres (w : ℤ × ℤ) (hw : w ∈ minusEightOddLatticeBall n) : (2 : ℤ) ∣ w.1 - 1 := by
      have ho := mt (minusEight_even_norm_iff w).mpr (Finset.mem_filter.mp hw).2
      rw [Int.dvd_iff_emod_eq_zero] at ho ⊢
      omega
    have h := minusEight_lattice_bound (s := minusEightOddLatticeBall n) (by omega : 0 < n)
      (fun _ hv => (mem_minusEight_lattice_ball.mp (Finset.mem_filter.mp hv).1).2)
      (d := 2) (by decide) 1 (by exact hres)
    have hR : (5 : ℝ) ≤ Real.sqrt n :=
      (Real.le_sqrt (by norm_num) (by positivity)).mpr (by exact_mod_cast hn)
    have hs := Real.sq_sqrt (show 0 ≤ (n : ℝ) by positivity)
    have hlin := mul_le_mul_of_nonneg_right hR (Real.sqrt_nonneg (n : ℝ))
    norm_num at h
    linarith only [h, hs, hlin, show (25 : ℝ) ≤ n by exact_mod_cast hn]
  · exact_mod_cast minusEight_odd_ball_small ⟨n, by omega⟩

end QuadraticNorms
