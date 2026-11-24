import QuadraticNorms.MinusElevenShell
import Mathlib.Data.Int.Interval

open scoped BigOperators

namespace QuadraticNorms

/-- The nonzero signed vectors of norm at most n. -/
def minusElevenLatticeBall (n : ℕ) : Finset MinusElevenVector :=
  (minusElevenCoordinateBox n).filter fun v => 0 < minusElevenNorm v ∧ minusElevenNorm v ≤ n

theorem mem_minusEleven_lattice_ball {n : ℕ} {v : MinusElevenVector} :
    v ∈ minusElevenLatticeBall n ↔ 0 < minusElevenNorm v ∧ minusElevenNorm v ≤ n := by
  exact ⟨fun hv => (Finset.mem_filter.mp hv).2, fun hv =>
    Finset.mem_filter.mpr ⟨mem_minusEleven_coordinate_box hv.2, hv⟩⟩

/-- All signed representations of a given norm. -/
def minusElevenNormFiber (a : ℕ) : Finset MinusElevenVector :=
  (minusElevenCoordinateBox a).filter fun v => minusElevenNorm v = a

theorem mem_minusEleven_norm_fiber {a : ℕ} {v : MinusElevenVector} :
    v ∈ minusElevenNormFiber a ↔ minusElevenNorm v = a := by
  exact ⟨fun hv => (Finset.mem_filter.mp hv).2, fun hv =>
    Finset.mem_filter.mpr ⟨mem_minusEleven_coordinate_box hv.le, hv⟩⟩

private theorem minusEleven_norm_real (v : MinusElevenVector) :
    (minusElevenNorm v : ℝ) =
      (v.1 : ℝ) ^ 2 + (v.1 : ℝ) * (v.2 : ℝ) + 3 * (v.2 : ℝ) ^ 2 := by
  rw [← Int.cast_natCast, minusEleven_norm_cast]
  simp [form]

private noncomputable def latticeRowBound (R : ℝ) (y : ℤ) : ℝ := 2 * R + 1 - 11 * (y : ℝ) ^ 2 / (4 * R)

private theorem minusEleven_lattice_ball_row_le {n : ℕ} (hn : 0 < n) (y : ℤ)
    (hy : 11 * (y : ℝ) ^ 2 ≤ 4 * n) :
    (((minusElevenLatticeBall n).filter fun v => v.2 = y).card : ℝ) ≤
      latticeRowBound (Real.sqrt n) y := by
  let row := (minusElevenLatticeBall n).filter (fun v => v.2 = y)
  let xs := row.image Prod.fst
  let R := Real.sqrt (n : ℝ)
  have hR : 0 < R := Real.sqrt_pos.2 (by exact_mod_cast hn)
  have hRsq : R ^ 2 = (n : ℝ) := Real.sq_sqrt (by positivity)
  have hrow (x : ℤ) (hx : x ∈ xs) : (2 * (x : ℝ) + y) ^ 2 ≤ 4 * n - 11 * (y : ℝ) ^ 2 := by
    obtain ⟨v, hv, rfl⟩ := Finset.mem_image.mp hx
    obtain ⟨hball, hvy⟩ := Finset.mem_filter.mp hv
    have hvalue : ((minusElevenNorm v : ℕ) : ℝ) ≤ n := by
      exact_mod_cast (mem_minusEleven_lattice_ball.mp hball).2
    have hnorm := minusEleven_norm_real v
    rw [hvy] at hnorm
    nlinarith only [hnorm, hvalue]
  have hc : row.card = xs.card := Finset.card_image_of_injOn (fun a ha b hb hab =>
    Prod.ext hab ((Finset.mem_filter.mp ha).2.trans (Finset.mem_filter.mp hb).2.symm)) |>.symm
  let W := Real.sqrt (4 * n - 11 * (y : ℝ)^2)
  have hWsq : W ^ 2 = 4 * n - 11 * (y : ℝ)^2 := Real.sq_sqrt (sub_nonneg.mpr hy)
  have hW : W ≤ 2 * R - 11 * (y : ℝ)^2 / (4 * R) := by
    have h := sq_nonneg (W - 2 * R)
    have hd : 0 < 4 * R := by positivity
    have hb : 11 * (y : ℝ)^2 / (4 * R) ≤ 2 * R - W := by
      apply (div_le_iff₀ hd).2
      nlinarith only [h, hWsq, hRsq]
    linarith only [hb]
  change (row.card : ℝ) ≤ _
  rw [hc]
  by_cases hx : xs.Nonempty
  · let a := xs.min' hx
    let b := xs.max' hx
    have ha : a ∈ xs := Finset.min'_mem _ _
    have hb : b ∈ xs := Finset.max'_mem _ _
    have hab : a ≤ b := Finset.min'_le _ _ hb
    have hcard := Finset.card_le_card (s := xs) (t := Finset.Icc a b)
      (fun x h => Finset.mem_Icc.mpr ⟨Finset.min'_le _ _ h, Finset.le_max' _ _ h⟩)
    have hcardI : ((Finset.Icc a b).card : ℤ) = b + 1 - a :=
      Int.card_Icc_of_le _ _ (by omega)
    have hcardR : (xs.card : ℝ) ≤ (b : ℝ) + 1 - a := by
      have hcardZ : (xs.card : ℤ) ≤ ((Finset.Icc a b).card : ℤ) := by exact_mod_cast hcard
      rw [hcardI] at hcardZ
      exact_mod_cast hcardZ
    have hnorm (x : ℤ) (hx : x ∈ xs) : |2 * (x : ℝ) + y| ≤ W := by
      apply (Real.le_sqrt (abs_nonneg _) (sub_nonneg.mpr hy)).2
      simpa only [sq_abs] using hrow x hx
    have ha' := (abs_le.mp (hnorm a ha)).1
    have hb' := (abs_le.mp (hnorm b hb)).2
    change (xs.card : ℝ) ≤ 2 * R + 1 - 11 * (y : ℝ)^2 / (4 * R)
    linarith only [hcardR, ha', hb', hW]
  · rw [Finset.not_nonempty_iff_eq_empty.mp hx]
    simp only [Finset.card_empty, Nat.cast_zero]
    change 0 ≤ 2 * R + 1 - 11 * (y : ℝ)^2 / (4 * R)
    linarith only [hW, Real.sqrt_nonneg (4 * n - 11 * (y : ℝ)^2)]

private theorem lattice_row_bound_sum (R : ℝ) (m : ℕ) :
    ∑ y ∈ Finset.Icc (-(m : ℤ)) m, latticeRowBound R y =
      (2 * (m : ℝ) + 1) * (2 * R + 1) - (11 / (12 * R)) * (m : ℝ) * (m + 1) * (2 * m + 1) := by
  induction m with
  | zero => simp [latticeRowBound]
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
    have hhi : ((m + 1 : ℕ) : ℤ) ∉ Finset.Icc (-(m : ℤ)) m := by
      simp [Finset.mem_Icc]
    rw [heq, Finset.sum_insert hlo, Finset.sum_insert hhi, ih]
    simp only [latticeRowBound, Nat.cast_succ, Int.cast_neg, Int.cast_add, Int.cast_one, Int.cast_natCast]
    by_cases hR : R = 0
    · simp [hR]
      ring
    · field_simp
      ring

private theorem minusEleven_lattice_ball_large {n : ℕ} (hn : 81 ≤ n) :
    ((minusElevenLatticeBall n).card : ℝ) ≤ (12 / 5) * n := by
  let R := Real.sqrt (n : ℝ)
  let B := 2 * R / Real.sqrt 11
  let m := ⌊B⌋₊
  let ys := Finset.Icc (-(m : ℤ)) m
  let F (t : ℝ) := (2 * t + 1) * (2 * R + 1) - (11 / (6 * R)) * t ^ 3
  have hnR : (81 : ℝ) ≤ n := by exact_mod_cast hn
  have hR : 0 < R := Real.sqrt_pos.2 (by linarith)
  have hR9 : 9 ≤ R := (Real.le_sqrt (by norm_num) (by positivity)).2 (by linarith)
  have hRsq : R ^ 2 = (n : ℝ) := Real.sq_sqrt (by positivity)
  have hB₀ : 0 ≤ B := by dsimp [B]; positivity
  have hm : (m : ℝ) ≤ B := Nat.floor_le hB₀
  have hBsq : 11 * B ^ 2 = 4 * n := by
    dsimp [B]
    rw [div_pow, Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 11)]
    linarith only [hRsq]
  have hmap : Set.MapsTo (Prod.snd : MinusElevenVector → ℤ) ↑(minusElevenLatticeBall n) ↑ys := by
    intro v hv
    have hvN : minusElevenNorm v ≤ n := (mem_minusEleven_lattice_ball.mp hv).2
    have hcast : ((minusElevenNorm v : ℕ) : ℝ) ≤ n := by exact_mod_cast hvN
    have hnorm := minusEleven_norm_real v
    have habs : |(v.2 : ℝ)| ≤ B := by
      apply (sq_le_sq₀ (abs_nonneg _) hB₀).mp
      rw [sq_abs]
      nlinarith only [hcast, hnorm, hBsq, sq_nonneg (2 * (v.1 : ℝ) + v.2)]
    have hnat : (v.2.natAbs : ℝ) ≤ B := by
      rw [← Int.cast_natCast, Int.natCast_natAbs, Int.cast_abs]
      exact habs
    have habsZ : |v.2| ≤ (m : ℤ) := by
      rw [← Int.natCast_natAbs]
      exact_mod_cast Nat.le_floor hnat
    exact Finset.mem_Icc.mpr (abs_le.mp habsZ)
  have hrowBound (y : ℤ) (hy : y ∈ ys) :
      (((minusElevenLatticeBall n).filter fun v => v.2 = y).card : ℝ) ≤ latticeRowBound R y := by
    obtain ⟨hy₁, hy₂⟩ := Finset.mem_Icc.mp hy
    have hyR : |(y : ℝ)| ≤ (m : ℝ) := by
      rw [abs_le]
      exact ⟨by exact_mod_cast hy₁, by exact_mod_cast hy₂⟩
    have hsq := (sq_le_sq₀ (abs_nonneg (y : ℝ)) hB₀).mpr (hyR.trans hm)
    rw [sq_abs] at hsq
    exact minusEleven_lattice_ball_row_le (by omega) y (by linarith only [hBsq, hsq])
  have hcalc : ∑ y ∈ ys, latticeRowBound R y ≤ F m := by
    rw [lattice_row_bound_sum R m]
    have hpoly : 2 * (m : ℝ) ^ 3 ≤ (m : ℝ) * (m + 1) * (2 * m + 1) := by
      calc
        _ = (m : ℝ) * m * (2 * m) := by ring
        _ ≤ _ := mul_le_mul (mul_le_mul_of_nonneg_left (by linarith) (by positivity))
          (by linarith) (by positivity) (by positivity)
    have h := mul_le_mul_of_nonneg_left hpoly (show 0 ≤ 11 / (12 * R) by positivity)
    have he : (11 / (12 * R)) * (2 * (m : ℝ) ^ 3) = (11 / (6 * R)) * (m : ℝ) ^ 3 := by
      field_simp
      ring
    rw [he] at h
    simpa only [F, mul_assoc] using sub_le_sub_left h ((2 * (m : ℝ) + 1) * (2 * R + 1))
  have h_comparison : F m ≤ F B := by
    have hz : (0 : ℝ) ≤ m := by positivity
    have hcross := mul_le_mul_of_nonneg_left hm hB₀
    have hsquare := pow_le_pow_left₀ hz hm 2
    have hquad : B ^ 2 + B * (m : ℝ) + (m : ℝ) ^ 2 ≤ 3 * B ^ 2 := by
      nlinarith only [hcross, hsquare]
    have hcoeff : 11 / (6 * R) * (3 * B ^ 2) = 2 * R := by
      field_simp
      nlinarith only [hRsq, hBsq]
    have hnonneg : 0 ≤ 4 * R + 2 - 11 / (6 * R) * (B ^ 2 + B * (m : ℝ) + (m : ℝ) ^ 2) := by
      have h := mul_le_mul_of_nonneg_left hquad (by positivity : 0 ≤ 11 / (6 * R))
      rw [hcoeff] at h
      linarith only [h, hR]
    have h := mul_nonneg (sub_nonneg.mpr hm) hnonneg
    dsimp [F]
    nlinarith only [h]
  have hB_le : B ≤ (20 / 33) * R := by
    have h11 : (33 / 10 : ℝ) ≤ Real.sqrt 11 :=
      (Real.le_sqrt (by norm_num) (by norm_num)).mpr (by norm_num)
    dsimp [B]
    apply (div_le_iff₀ (Real.sqrt_pos.2 (by norm_num))).mpr
    nlinarith only [mul_le_mul_of_nonneg_left h11 hR.le]
  have hSimplify : F B = 2 * R + 1 + ((10 / 3) * R + 2) * B := by
    have hcoeff : 11 / (6 * R) * B ^ 2 = (2 / 3) * R := by
      field_simp
      nlinarith only [hRsq, hBsq]
    dsimp [F]
    rw [pow_succ, ← mul_assoc, hcoeff]
    ring
  calc
    ((minusElevenLatticeBall n).card : ℝ) =
        ∑ y ∈ ys, (((minusElevenLatticeBall n).filter fun v => v.2 = y).card : ℝ) := by
      exact_mod_cast Finset.card_eq_sum_card_fiberwise hmap
    _ ≤ ∑ y ∈ ys, latticeRowBound R y := Finset.sum_le_sum hrowBound
    _ ≤ F m := hcalc
    _ ≤ F B := h_comparison
    _ ≤ (200 / 99) * n + (106 / 33) * R + 1 := by
      rw [hSimplify]
      have h := mul_le_mul_of_nonneg_left hB_le (show 0 ≤ (10 / 3) * R + 2 by positivity)
      nlinarith only [h, hRsq]
    _ ≤ _ := by
      have h := mul_le_mul_of_nonneg_right hR9 hR.le
      rw [← pow_two] at h
      linarith only [h, hRsq, hnR]

set_option maxRecDepth 9000 in
set_option maxHeartbeats 1000000 in
private theorem minusEleven_lattice_ball_small :
    ∀ j : Fin 81, 5 * (minusElevenLatticeBall j.val).card ≤ 12 * j.val := by
  decide +kernel

/-- The lattice point estimate for the −11 norm in the large-prime sieve.
The row count uses a quadratic upper bound for the square root, and
kernel computation handles the initial range. -/
theorem minusEleven_lattice_ball_card_le (n : ℕ) :
    ((minusElevenLatticeBall n).card : ℝ) ≤ (12 / 5) * n := by
  by_cases hn : 81 ≤ n
  · exact minusEleven_lattice_ball_large hn
  · have h := minusEleven_lattice_ball_small ⟨n, Nat.lt_of_not_ge hn⟩
    have hR : 5 * ((minusElevenLatticeBall n).card : ℝ) ≤ 12 * n := by exact_mod_cast h
    linarith only [hR]

end QuadraticNorms
