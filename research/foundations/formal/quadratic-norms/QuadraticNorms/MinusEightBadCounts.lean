import QuadraticNorms.MinusEightRepresentationCounts

open scoped BigOperators

namespace QuadraticNorms

private theorem minusEight_scaled_first_bad_card_le {n p : ℕ}
    (s : Finset MinusEightShellPoint) (ball : Finset (ℤ × ℤ))
    (hfull : s ⊆ minusEightFullShell n) (hpos : ∀ v ∈ s, 0 < minusEightNorm v.2)
    (hball : ∀ u w, (minusEightScaleVector p u, w) ∈ s → u ∈ ball) :
    (s.filter (fun v => MinusEightDivisible p v.1)).card ≤
      ball.card * (2 * minusEightDivisorMaximum n) := by
  let base := ball.filter (fun u => p ^ 2 * minusEightNorm u < n)
  let fiber (u : ℤ × ℤ) := minusEightNormFiber (n - p ^ 2 * minusEightNorm u)
  let reconstruct (v : Σ _ : ℤ × ℤ, ℤ × ℤ) : MinusEightShellPoint :=
    (minusEightScaleVector p v.1, v.2)
  have hsubset : (s.filter fun v => MinusEightDivisible p v.1) ⊆
      (base.sigma fiber).image reconstruct := by
    rintro ⟨x, w⟩ hv
    obtain ⟨hs, hd⟩ := Finset.mem_filter.mp hv
    obtain ⟨u, rfl⟩ := (minusEight_divisible_iff_scaled p x).mp hd
    have h := mem_minusEight_full_shell.mp (hfull hs)
    change minusEightNorm (minusEightScaleVector p u) + minusEightNorm w = n at h
    rw [minusEight_norm_scale] at h
    have hw : 0 < minusEightNorm w := hpos _ hs
    have hlt : p ^ 2 * minusEightNorm u < n := by omega
    refine Finset.mem_image.mpr ⟨⟨u, w⟩, Finset.mem_sigma.mpr
      ⟨Finset.mem_filter.mpr ⟨hball u w hs, hlt⟩,
        mem_minusEight_norm_fiber.mpr (by dsimp [fiber]; omega)⟩, rfl⟩
  calc
    _ ≤ ((base.sigma fiber).image reconstruct).card := Finset.card_le_card hsubset
    _ ≤ (base.sigma fiber).card := Finset.card_image_le
    _ = ∑ a ∈ base, (fiber a).card := Finset.card_sigma base fiber
    _ ≤ ∑ _a ∈ base, 2 * minusEightDivisorMaximum n := by
      refine Finset.sum_le_sum fun a ha => ?_
      have hpos := Nat.sub_pos_of_lt (Finset.mem_filter.mp ha).2
      exact (minusEight_norm_fiber_card_le hpos).trans
        (Nat.mul_le_mul_left _ (minusEight_split_divisors_le_maximum hpos (Nat.sub_le ..)))
    _ = base.card * (2 * minusEightDivisorMaximum n) := by simp
    _ ≤ _ := Nat.mul_le_mul_right _ (Finset.card_filter_le ..)

private theorem minusEight_selected_bad_card_le_twice_first (n p : ℕ)
    (s : Finset MinusEightShellPoint) (hfull : s ⊆ minusEightFullShell n)
    (hswap : ∀ v : MinusEightShellPoint, v.swap ∈ s ↔ v ∈ s) :
    (s ∩ minusEightEitherDivisible n p).card ≤ 2 * (s.filter (fun v => MinusEightDivisible p v.1)).card := by
  let e := (Equiv.prodComm (ℤ × ℤ) (ℤ × ℤ)).toEmbedding
  let f := s.filter (fun v => MinusEightDivisible p v.1)
  have he : s.filter (fun v => MinusEightDivisible p v.2) = f.map e := by
    ext v
    simp only [f, Finset.mem_filter, Finset.mem_map]
    constructor
    · intro ⟨h, hd⟩
      exact ⟨v.swap, ⟨(hswap v).mpr h, hd⟩, rfl⟩
    · rintro ⟨w, ⟨h, hd⟩, rfl⟩
      exact ⟨(hswap w).mpr h, hd⟩
  have ht : s ∩ minusEightEitherDivisible n p =
      s.filter (fun v => MinusEightDivisible p v.1 ∨ MinusEightDivisible p v.2) := by
    ext v
    simp only [minusEightEitherDivisible, Finset.mem_inter, Finset.mem_filter]
    exact ⟨fun h => ⟨h.1, h.2.2⟩, fun h => ⟨h.1, hfull h.1, h.2⟩⟩
  rw [ht, Finset.filter_or, he]
  exact (Finset.card_union_le ..).trans_eq (by rw [Finset.card_map, two_mul])

private theorem minusEight_inert_large_bad_count {n p : ℕ} (hp : MinusEightInertPrime p)
    {s : Finset MinusEightShellPoint} (hfull : s ⊆ minusEightFullShell n)
    (hpos : ∀ v ∈ s, 0 < minusEightNorm v.2)
    (hswap : ∀ v : MinusEightShellPoint, v.swap ∈ s ↔ v ∈ s)
    (ball : Finset (ℤ × ℤ)) {K : ℝ} (hK : 0 ≤ K)
    (hcard : (ball.card : ℝ) ≤ K * ((n / p ^ 2 : ℕ) : ℝ))
    (hball : ∀ u w, (minusEightScaleVector p u, w) ∈ s → u ∈ ball) :
    ((s ∩ minusEightEitherDivisible n p).card : ℝ) ≤
      (4 * K) * minusEightDivisorMaximum n * n / (p : ℝ) ^ 2 := by
  have h := (minusEight_selected_bad_card_le_twice_first n p s hfull hswap).trans
    (Nat.mul_le_mul_left 2 (minusEight_scaled_first_bad_card_le (p := p) s ball hfull hpos hball))
  have hR : ((s ∩ minusEightEitherDivisible n p).card : ℝ) ≤
      2 * (ball.card : ℝ) * (2 * minusEightDivisorMaximum n) := by
    rw [mul_assoc]
    exact_mod_cast h
  have hpR : 0 < (p : ℝ) ^ 2 := by exact_mod_cast pow_pos hp.1.pos 2
  have hdiv : ((n / p ^ 2 : ℕ) : ℝ) ≤ (n : ℝ) / (p : ℝ) ^ 2 := by
    apply (le_div_iff₀ hpR).mpr
    exact_mod_cast Nat.div_mul_le_self n (p ^ 2)
  calc
    _ ≤ _ := hR
    _ ≤ 2 * (K * ((n : ℝ) / (p : ℝ) ^ 2)) * (2 * minusEightDivisorMaximum n) := by
      exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left
        (hcard.trans (mul_le_mul_of_nonneg_left hdiv hK)) (by norm_num)) (by positivity)
    _ = _ := by ring

private theorem minusEight_scaled_mem_lattice_ball {n p : ℕ} (hp : 0 < p) {u w : ℤ × ℤ}
    (hs : (minusEightScaleVector p u, w) ∈ minusEightFullShell n)
    (hu : 0 < minusEightNorm (minusEightScaleVector p u)) :
    u ∈ minusEightLatticeBall (n / p ^ 2) := by
  have h := mem_minusEight_full_shell.mp hs
  change minusEightNorm (minusEightScaleVector p u) + minusEightNorm w = n at h
  rw [minusEight_norm_scale] at hu h
  refine mem_minusEight_lattice_ball.mpr ⟨Nat.pos_of_mul_pos_left hu, ?_⟩
  apply (Nat.le_div_iff_mul_le (pow_pos hp 2)).mpr
  rw [mul_comm]
  omega

/-- The large-prime bad count for the dyadic shell, with the sharp
lattice coefficient needed for the odd selected shell. -/
theorem minusEight_inert_dyadic_bad_card_le {n p : ℕ} (hp : MinusEightInertPrime p) :
    ((minusEightDyadicShell n ∩ minusEightEitherDivisible n p).card : ℝ) ≤
      (32 / 3) * minusEightDivisorMaximum n * n / (p : ℝ) ^ 2 := by
  have hfull : minusEightDyadicShell n ⊆ minusEightFullShell n := Finset.filter_subset _ _
  have h := minusEight_inert_large_bad_count hp hfull
      (s := minusEightDyadicShell n)
      (fun v hv => minusEight_norm_pos_of_not_divisible_two v.2 (not_or.mp (Finset.mem_filter.mp hv).2).2)
      (fun v => by
        simp only [minusEightDyadicShell, Finset.mem_filter, mem_minusEight_full_shell, Prod.swap]
        rw [add_comm (minusEightNorm v.2)]
        tauto)
      (minusEightLatticeBall (n / p ^ 2)) (by norm_num : (0 : ℝ) ≤ 8 / 3)
      (minusEight_lattice_ball_card_le _) (fun u w h =>
        minusEight_scaled_mem_lattice_ball hp.1.pos (hfull h)
          (minusEight_norm_pos_of_not_divisible_two _
            (not_or.mp (Finset.mem_filter.mp h).2).1))
  norm_num at h ⊢
  exact h

/-- On the odd-odd shell the dilated block and its cofactor have odd
norm, so the sharper odd lattice count applies. -/
theorem minusEight_inert_odd_odd_bad_card_le {n p : ℕ} (hp : MinusEightInertPrime p) :
    ((minusEightOddOddShell n ∩ minusEightEitherDivisible n p).card : ℝ) ≤
      8 * minusEightDivisorMaximum n * n / (p : ℝ) ^ 2 := by
  have hfull : minusEightOddOddShell n ⊆ minusEightFullShell n := Finset.filter_subset _ _
  have h := minusEight_inert_large_bad_count hp hfull
    (fun v hv => minusEight_norm_pos_of_not_divisible_two v.2
      (fun h => (Finset.mem_filter.mp hv).2.2
        ((show 2 ∣ 4 by decide).trans ((minusEight_four_dvd_norm_iff _).mpr h))))
    (fun v => by
      simp only [minusEightOddOddShell, Finset.mem_filter, mem_minusEight_full_shell, Prod.swap]
      rw [add_comm (minusEightNorm v.2)]
      tauto)
    (minusEightOddLatticeBall (n / p ^ 2)) (by norm_num : (0 : ℝ) ≤ 2)
    (minusEight_odd_lattice_ball_card_le _) (fun u w h => by
      obtain ⟨hs, hu, _⟩ := Finset.mem_filter.mp h
      have hu' : 0 < minusEightNorm (minusEightScaleVector p u) := by
        by_contra! hz
        exact hu (by rw [Nat.eq_zero_of_le_zero hz]; exact dvd_zero 2)
      exact Finset.mem_filter.mpr ⟨minusEight_scaled_mem_lattice_ball hp.1.pos hs hu',
        fun hd => hu (by rw [minusEight_norm_scale]; exact hd.trans (dvd_mul_left _ _))⟩)
  norm_num at h
  exact h

end QuadraticNorms
