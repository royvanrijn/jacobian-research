import QuadraticNorms.MinusElevenLattice

open scoped BigOperators

namespace QuadraticNorms

private noncomputable def minusElevenBadFirst (n p : ℕ) : Finset MinusElevenShellPoint :=
  (minusElevenOddShell n).filter fun v => p ∣ minusElevenNorm v.1

private theorem minusEleven_bad_first_card_le {n p : ℕ} (hp : MinusElevenInertPrime p)
    (hReps : ∀ a, 0 < a → a ≤ n → (minusElevenNormFiber a).card ≤ 2 * minusElevenDivisorMaximum n) :
    (minusElevenBadFirst n p).card ≤
      (minusElevenLatticeBall (n / p ^ 2)).card * (2 * minusElevenDivisorMaximum n) := by
  let base := (minusElevenLatticeBall (n / p ^ 2)).filter (fun u => p ^ 2 * minusElevenNorm u < n)
  let fiber (u : MinusElevenVector) := minusElevenNormFiber (n - p ^ 2 * minusElevenNorm u)
  let reconstruct (v : Σ _ : MinusElevenVector, MinusElevenVector) : MinusElevenShellPoint :=
    (minusElevenScaleVector p v.1, v.2)
  have hsubset : minusElevenBadFirst n p ⊆ (base.sigma fiber).image reconstruct := by
    rintro ⟨⟨x, y⟩, v⟩ hbad
    obtain ⟨hshell, hd⟩ := Finset.mem_filter.mp hbad
    obtain ⟨⟨a, rfl⟩, ⟨b, rfl⟩⟩ :=
      (minusEleven_inert_dvd_norm_iff hp (x, y)).mp hd
    let u : MinusElevenVector := (a, b)
    have h := mem_minusEleven_odd_shell.mp hshell
    change minusElevenNorm (minusElevenScaleVector p u) % 2 = 1 ∧ minusElevenNorm v % 2 = 1 ∧
      minusElevenNorm (minusElevenScaleVector p u) + minusElevenNorm v = n at h
    rw [minusEleven_norm_scale] at h
    have hu_pos : 0 < minusElevenNorm u := by
      by_contra! hu
      have hu_zero : minusElevenNorm u = 0 := by omega
      simp [hu_zero] at h
    have hlt : p ^ 2 * minusElevenNorm u < n := by omega
    have hu_le : minusElevenNorm u ≤ n / p ^ 2 := by
      apply (Nat.le_div_iff_mul_le (pow_pos hp.1.pos _)).mpr
      simpa only [mul_comm] using hlt.le
    exact Finset.mem_image.mpr
      ⟨⟨u, v⟩, Finset.mem_sigma.mpr
        ⟨Finset.mem_filter.mpr ⟨mem_minusEleven_lattice_ball.mpr ⟨hu_pos, hu_le⟩, hlt⟩,
          mem_minusEleven_norm_fiber.mpr (by dsimp [fiber]; omega)⟩, rfl⟩
  calc
    _ ≤ ((base.sigma fiber).image reconstruct).card := Finset.card_le_card hsubset
    _ ≤ (base.sigma fiber).card := Finset.card_image_le
    _ = ∑ a ∈ base, (fiber a).card := Finset.card_sigma base fiber
    _ ≤ ∑ _a ∈ base, 2 * minusElevenDivisorMaximum n := by
      refine Finset.sum_le_sum (fun a ha => hReps _ ?_ (Nat.sub_le ..))
      exact Nat.sub_pos_of_lt (Finset.mem_filter.mp ha).2
    _ = base.card * (2 * minusElevenDivisorMaximum n) := by
      simp only [Finset.sum_const, nsmul_eq_mul, Nat.cast_id]
    _ ≤ _ := Nat.mul_le_mul_right _ (Finset.card_filter_le ..)

private theorem minusEleven_bad_card_le_twice_first (n p : ℕ) :
    (minusElevenBadShell n p).card ≤ 2 * (minusElevenBadFirst n p).card := by
  let swap := (Equiv.prodComm MinusElevenVector MinusElevenVector).toEmbedding
  have heq : (minusElevenOddShell n).filter (fun v => p ∣ minusElevenNorm v.2) =
      (minusElevenBadFirst n p).map swap := by
    have hswap (v : MinusElevenShellPoint) : swap v ∈ minusElevenOddShell n ↔
        v ∈ minusElevenOddShell n := by
      simp only [mem_minusEleven_odd_shell, swap, Equiv.toEmbedding_apply, Equiv.prodComm_apply]
      rw [add_comm]
      tauto
    ext v
    constructor
    · intro hv
      obtain ⟨h, hd⟩ := Finset.mem_filter.mp hv
      exact Finset.mem_map.mpr ⟨swap v, Finset.mem_filter.mpr ⟨(hswap v).mpr h, hd⟩, rfl⟩
    · intro hv
      obtain ⟨w, hw, rfl⟩ := Finset.mem_map.mp hv
      obtain ⟨h, hd⟩ := Finset.mem_filter.mp hw
      exact Finset.mem_filter.mpr ⟨(hswap w).mpr h, hd⟩
  unfold minusElevenBadShell
  rw [Finset.filter_or]
  change (minusElevenBadFirst n p ∪
    (minusElevenOddShell n).filter (fun v => p ∣ minusElevenNorm v.2)).card ≤ _
  rw [heq]
  exact (Finset.card_union_le ..).trans_eq (by rw [Finset.card_map, two_mul])

/-- The large inert-prime bad-count estimate follows from the number of
representations of each positive norm. The lattice-ball estimate and
the two orientations are included in the proof. -/
theorem minusEleven_inert_large_bad_count {n p : ℕ} (hp : MinusElevenInertPrime p)
    (hReps : ∀ a, 0 < a → a ≤ n → (minusElevenNormFiber a).card ≤ 2 * minusElevenDivisorMaximum n) :
    ((minusElevenBadShell n p).card : ℝ) ≤
      (48 / 5) * minusElevenDivisorMaximum n * n / (p : ℝ) ^ 2 := by
  have h := (minusEleven_bad_card_le_twice_first n p).trans
    (Nat.mul_le_mul_left 2 (minusEleven_bad_first_card_le hp hReps))
  have hR : ((minusElevenBadShell n p).card : ℝ) ≤
      2 * ((minusElevenLatticeBall (n / p ^ 2)).card : ℝ) * (2 * minusElevenDivisorMaximum n) := by
    rw [mul_assoc]
    exact_mod_cast h
  have hpR : 0 < (p : ℝ) ^ 2 := by exact_mod_cast pow_pos hp.1.pos 2
  have hdiv : ((n / p ^ 2 : ℕ) : ℝ) ≤ (n : ℝ) / (p : ℝ) ^ 2 := by
    apply (le_div_iff₀ hpR).mpr
    exact_mod_cast Nat.div_mul_le_self n (p ^ 2)
  calc
    _ ≤ _ := hR
    _ ≤ 2 * ((12 / 5) * ((n : ℝ) / (p : ℝ) ^ 2)) * (2 * minusElevenDivisorMaximum n) := by
      have hBall := (minusEleven_lattice_ball_card_le (n / p ^ 2)).trans
        (mul_le_mul_of_nonneg_left hdiv (by norm_num))
      exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hBall (by norm_num)) (by positivity)
    _ = _ := by ring

end QuadraticNorms
