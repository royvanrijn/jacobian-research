import MathieuQuadratic.Draft
import Mathlib.Algebra.Order.Floor.Div

/-!
# Sharp multiplier bounds in the product-phase theorem

The trace formula gives both the last detecting order below the degree
threshold and the detecting norm power at the threshold.
-/

namespace MathieuQuadratic

variable {K A : Type*} [Field K] [CharZero K] [CommRing A] [Algebra K A]

private theorem choose_ne_zero_of_neg {u : ℤ} (hu : u < 0) (n : ℕ) :
    Ring.choose u n ≠ 0 := by
  let a := (-u + n - 1 : ℤ).toNat
  have ha : (a : ℤ) = -u + n - 1 := Int.toNat_of_nonneg (by omega)
  have hn : n ≤ a := by omega
  rw [← neg_neg u, Ring.choose_neg, ← ha, Ring.choose_natCast, Units.smul_def, zsmul_eq_mul]
  exact mul_ne_zero (Int.negOnePow (n : ℤ)).ne_zero (by exact_mod_cast Nat.choose_ne_zero_iff.mpr hn)

private theorem choose_eq_zero_of_nonneg_lt {u : ℤ} (hu : 0 ≤ u) {n : ℕ} (hn : u < n) :
    Ring.choose u n = 0 := by
  rw [← Int.toNat_of_nonneg hu, Ring.choose_natCast, Nat.choose_eq_zero_of_lt (by omega), Nat.cast_zero]

variable (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
variable {M γ : ℕ} (hM : 0 < M) (hγ : M ≤ γ)
variable (weight : ∀ (a : ℤ) (i : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ i) = 0)
variable (height : ∀ i : ℕ, trace (S ^ i) = (M : K) / (↑(M + 2 * i) : K))

include hM hγ weight height in
/-- Every monomial in the block average below the degree threshold
vanishes from the claimed cutoff onward. -/
theorem product_cutoff_vanishing {J ℓ j m : ℕ} (hJ : J < M)
    (hj : j ≤ ℓ) (hℓ : ℓ ≤ J) (hm : (2 * J ⌈/⌉ M) + 1 ≤ m) :
    trace ((U : A) ^ ℓ * S ^ j * productPhase U S M M γ ^ m) = 0 := by
  have hm₀ : 0 < m := by omega
  have hprod := Nat.mul_le_mul_right m hγ
  have hsmall := Nat.le_mul_of_pos_right M hm₀
  have hcut : 2 * J ≤ M * (m - 1) := (ceilDiv_le_iff_le_mul hM).1 (by omega)
  rw [Nat.mul_sub, mul_one] at hcut
  have hup : 0 ≤ ((M * m : ℕ) - (M + 2 * j : ℕ) : ℤ) := by omega
  have hz : (U : A) ^ ℓ = (↑(U ^ (ℓ : ℤ)) : A) := by simp
  rw [hz, product_mixedMoment trace U S M M M γ m j hM ℓ (by omega) weight height,
    if_pos (by omega : (ℓ : ℤ) ≤ (M * m : ℕ)),
    choose_eq_zero_of_nonneg_lt hup (by omega)]
  simp

include hM hγ weight height in
/-- The block power Bᴶ detects at the last possible order below cutoff. -/
theorem product_cutoff_sharp {J : ℕ} (hJ₀ : 0 < J) (hJ : J < M) :
    let m := 2 * J ⌈/⌉ M
    trace ((U : A) ^ J * S ^ J * productPhase U S M M γ ^ m) ≠ 0 := by
  let m := 2 * J ⌈/⌉ M
  change trace ((U : A) ^ J * S ^ J * productPhase U S M M γ ^ m) ≠ 0
  have hc : 2 * J ≤ M * m := (ceilDiv_le_iff_le_mul hM).1 le_rfl
  have hm₀ : 0 < m := by
    by_contra! h
    simp [show m = 0 by omega] at hc
    omega
  have hsmall := Nat.le_mul_of_pos_right M hm₀
  have hprod := Nat.mul_le_mul_right m hγ
  have hcut : ¬ 2 * J ≤ M * (m - 1) := by
    rw [← ceilDiv_le_iff_le_mul hM]
    change ¬ m ≤ m - 1
    omega
  rw [Nat.mul_sub, mul_one] at hcut
  have hup : ((M * m : ℕ) - (M + 2 * J : ℕ) : ℤ) < 0 := by omega
  have hz : (U : A) ^ J = (↑(U ^ (J : ℤ)) : A) := by simp
  rw [hz, product_mixedMoment trace U S M M M γ m J hM J (by omega) weight height,
    if_pos (by omega : (J : ℤ) ≤ (M * m : ℕ))]
  apply (map_ne_zero_iff (algebraMap ℚ K) (RingHom.injective (algebraMap ℚ K))).mpr
  exact mul_ne_zero (mixedBetaMoment_pos _ _ _ hM).ne'
    (by exact_mod_cast choose_ne_zero_of_neg hup _)

include hM weight height in
/-- In the product family with γ ≥ α, the purity inequalities are necessary
and sufficient. -/
theorem product_pure_iff {α β : ℕ} (hαγ : α ≤ γ) :
    Pure trace (productPhase U S α β γ) ↔ M ≤ β ∧ β ≤ α := by
  have hval (m : ℕ) :
      trace (productPhase U S α β γ ^ m) = algebraMap ℚ K
        (mixedBetaMoment M (γ * m) 0 *
          ((Ring.choose ((β * m : ℕ) - M : ℤ) (α * m) : ℤ) : ℚ)) := by
    have hc := Nat.mul_le_mul_right m hαγ
    have hh := product_mixedMoment trace U S M α β γ m 0 hM 0 (by omega) weight height
    rw [if_pos (by positivity : (0 : ℤ) ≤ (α * m : ℕ)), sub_zero, Int.toNat_natCast] at hh
    simpa using hh
  have hzero (m : ℕ) :
      trace (productPhase U S α β γ ^ m) = 0 ↔
        Ring.choose ((β * m : ℕ) - M : ℤ) (α * m) = 0 := by
    rw [hval, map_eq_zero_iff _ (RingHom.injective (algebraMap ℚ K)), mul_eq_zero]
    simp [(mixedBetaMoment_pos _ _ _ hM).ne']
  constructor
  · intro hp
    have hβ : M ≤ β := by
      by_contra! h
      have hn := choose_ne_zero_of_neg (show ((β * 1 : ℕ) - M : ℤ) < 0 by omega) (α * 1)
      exact hn ((hzero 1).1 (hp 1 Nat.one_pos))
    refine ⟨hβ, ?_⟩
    by_contra! h
    have hmul : (α + 1) * M ≤ β * M := Nat.mul_le_mul_right M h
    rw [Nat.add_mul, one_mul] at hmul
    have hu : 0 ≤ ((β * M : ℕ) - M : ℤ) := by omega
    have hn : α * M ≤ (((β * M : ℕ) - M : ℤ).toNat) := by omega
    have hz := (hzero M).1 (hp M hM)
    rw [← Int.toNat_of_nonneg hu, Ring.choose_natCast, Nat.cast_eq_zero] at hz
    exact (Nat.choose_ne_zero_iff.mpr hn) hz
  · rintro ⟨hβ, hαβ⟩ m hm
    apply (hzero m).2
    have hmβ := Nat.le_mul_of_pos_right β hm
    have hmul := Nat.mul_le_mul_right m hαβ
    exact choose_eq_zero_of_nonneg_lt (by omega) (by omega)

include hM hγ weight height in
/-- The least positive cutoff on the span of the monomials of degree at most 2J. -/
theorem product_cutoff_span {J : ℕ} (hJ : J < M) {m₀ : ℕ} (hm₀ : 0 < m₀) :
    (∀ h ∈ Submodule.span K {v : A | ∃ ℓ j : ℕ, j ≤ ℓ ∧ ℓ ≤ J ∧ v = (U : A) ^ ℓ * S ^ j},
      ∀ m ≥ m₀, trace (h * productPhase U S M M γ ^ m) = 0) ↔
        (2 * J ⌈/⌉ M) + 1 ≤ m₀ := by
  constructor
  · intro hz
    by_contra! h
    let m := 2 * J ⌈/⌉ M
    have hm : m₀ ≤ m := by omega
    have hJ₀ : 0 < J := by
      by_contra! hj
      have hzero : m = 0 := by simp [m, show J = 0 by omega, zero_ceilDiv]
      omega
    exact product_cutoff_sharp trace U S hM hγ weight height hJ₀ hJ
      (hz _ (Submodule.subset_span ⟨J, J, le_rfl, le_rfl, rfl⟩) _ hm)
  · intro h hmul hmem m hm
    induction hmem using Submodule.span_induction with
    | mem x hx =>
        obtain ⟨ℓ, j, hj, hℓ, rfl⟩ := hx
        exact product_cutoff_vanishing trace U S hM hγ weight height hJ hj hℓ (h.trans hm)
    | zero => simp
    | add x y _ _ hx hy => rw [add_mul, map_add, hx, hy, add_zero]
    | smul c x _ hx => rw [smul_mul_assoc, map_smul, hx, smul_zero]

include hM weight height in
/-- When γ = α > β, the inverse norm selects precisely the first discrepancy. -/
theorem product_inverse_norm_moment {α β : ℕ} (hβ : M ≤ β) (hαβ : β < α)
    (m : ℕ) (hm : 0 < m) :
    trace ((↑(U⁻¹) : A) * productPhase U S α β α ^ m) =
      algebraMap ℚ K ((M : ℚ) * (-2) ^ (α * m) / (α * m + 1 : ℕ)) := by
  have he := product_mixedMoment_first_error trace U S M α β α m 0 hM weight height
  simp only [sub_self, zero_sub, zpow_neg_one, pow_zero, mul_one, mul_zero, add_zero] at he
  rw [he]
  have hmβ := Nat.le_mul_of_pos_right β hm
  have hc := Nat.mul_le_mul_right m hαβ.le
  rw [choose_eq_zero_of_nonneg_lt (by omega) (by omega)]
  simp

include hM weight height in
/-- Eventual vanishing of a fixed Laurent monomial in the product family
when α > β, away from the first omitted coefficient. -/
theorem product_laurent_monomial_eventually_zero {α β : ℕ} (hβ : M ≤ β) (hαβ : β < α)
    (hαγ : α ≤ γ) (ℓ : ℤ) (j : ℕ) (hℓ : 0 ≤ ℓ ∨ α < γ) :
    ∃ m₀ : ℕ, ∀ m ≥ m₀,
      trace ((↑(U ^ ℓ) : A) * S ^ j * productPhase U S α β γ ^ m) = 0 := by
  refine ⟨ℓ.natAbs + M + 2 * j + 1, fun m hm ↦ ?_⟩
  have hm₀ : 0 < m := by omega
  have hmul := Nat.mul_le_mul_right m (show β + 1 ≤ α by omega)
  have hβm := Nat.mul_le_mul_right m (show 1 ≤ β by omega)
  rw [Nat.add_mul, one_mul] at hmul
  rw [one_mul] at hβm
  have hrange : (α : ℤ) * m - (γ : ℤ) * m ≤ ℓ := by
    have hc := Nat.mul_le_mul_right m hαγ
    rcases hℓ with hℓ | h
    · omega
    · have hγm := Nat.mul_le_mul_right m (show α + 1 ≤ γ by omega)
      rw [Nat.add_mul, one_mul] at hγm
      omega
  rw [product_mixedMoment trace U S M α β γ m j hM ℓ hrange weight height]
  split_ifs with h
  · rw [choose_eq_zero_of_nonneg_lt (by omega) (by omega)]
    simp
  · rfl

include hM weight height in
/-- The marker qᴹ works on the entire boundary α = β of the pure family. -/
theorem product_equal_marked_norm {α : ℕ} (hα : M ≤ α) (hαγ : α ≤ γ)
    (m : ℕ) (hm : 0 < m) :
    trace ((U : A) ^ M * productPhase U S α α γ ^ m) =
      algebraMap ℚ K (mixedBetaMoment M (γ * m) 0) := by
  have hmα := Nat.le_mul_of_pos_right α hm
  have hprod := Nat.mul_le_mul_right m hαγ
  have hsub : M ≤ α * m := by omega
  have he : (U : A) ^ M * productPhase U S α α γ ^ m =
      (↑(U ^ (M : ℤ)) : A) * S ^ 0 * productPhase U S α α γ ^ m := by simp
  rw [he, product_mixedMoment trace U S M α α γ m 0 hM M (by omega) weight height,
    if_pos (by exact_mod_cast hsub),
    show ((α * m : ℕ) - (M + 2 * 0 : ℕ) : ℤ) = (↑(α * m - M) : ℤ) by
      simpa using (Int.ofNat_sub hsub).symm,
    ← Int.ofNat_sub hsub, Int.toNat_natCast, Ring.choose_natCast, Nat.choose_self]
  simp

include hM weight height in
/-- The detecting-multiplier alternatives for the product family, given
the finite block expansion of each marked trace. -/
theorem product_family_detects_iff {α β : ℕ} (hβ : M ≤ β) (hβα : β ≤ α) (hαγ : α ≤ γ)
    (average : ∀ h : A, ∃ (support : Finset (ℤ × ℕ)) (c : ℤ × ℕ → K),
      ∀ m : ℕ, 0 < m → trace (h * productPhase U S α β γ ^ m) =
        trace ((∑ p ∈ support, c p • ((↑(U ^ p.1) : A) * S ^ p.2)) *
          productPhase U S α β γ ^ m)) :
    (∃ h : A, ¬ ∃ m₀ : ℕ, ∀ m ≥ m₀,
        trace (h * productPhase U S α β γ ^ m) = 0) ↔
      α = β ∨ α = γ := by
  constructor
  · rintro ⟨h, hh⟩
    by_contra! hn
    obtain ⟨support, c, havg⟩ := average h
    have hp (p : ℤ × ℕ) := product_laurent_monomial_eventually_zero trace U S hM weight height
      hβ (by omega : β < α) hαγ p.1 p.2 (Or.inr (by omega))
    choose n hn' using hp
    refine hh ⟨support.sup n + 1, fun m hm ↦ ?_⟩
    rw [havg m (by omega), Finset.sum_mul, map_sum]
    apply Finset.sum_eq_zero
    intro p hp
    rw [smul_mul_assoc, map_smul, hn' p m ((Finset.le_sup hp).trans (by omega)), smul_zero]
  · intro h
    rcases h with rfl | rfl
    · refine ⟨(U : A) ^ M, ?_⟩
      rintro ⟨m₀, h⟩
      have hz := h (m₀ + 1) (by omega)
      rw [product_equal_marked_norm trace U S hM weight height hβ hαγ _ (by omega)] at hz
      exact (map_ne_zero_iff (algebraMap ℚ K) (RingHom.injective (algebraMap ℚ K))).mpr
        (mixedBetaMoment_pos M _ 0 hM).ne' hz
    · by_cases h : α = β
      · subst β
        refine ⟨(U : A) ^ M, ?_⟩
        rintro ⟨m₀, h⟩
        have hz := h (m₀ + 1) (by omega)
        rw [product_equal_marked_norm trace U S hM weight height hβ le_rfl _ (by omega)] at hz
        exact (map_ne_zero_iff (algebraMap ℚ K) (RingHom.injective (algebraMap ℚ K))).mpr
          (mixedBetaMoment_pos M _ 0 hM).ne' hz
      · refine ⟨(↑(U⁻¹) : A), ?_⟩
        rintro ⟨m₀, hz⟩
        have hh := hz (m₀ + 1) (by omega)
        rw [product_inverse_norm_moment trace U S hM weight height hβ (by omega) _ (by omega)] at hh
        apply (map_ne_zero_iff (algebraMap ℚ K) (RingHom.injective (algebraMap ℚ K))).mpr _ hh
        exact div_ne_zero (mul_ne_zero (by exact_mod_cast hM.ne') (pow_ne_zero _ (by norm_num)))
          (by positivity)

include hM weight height in
/-- The complementary index and rising-factorial factor in the draft's
fixed-multiplier formula, at orders beyond the finite-support cutoff. -/
theorem product_tail_monomial (ℓ : ℤ) (j m : ℕ) (hγ : M < γ)
    (hsize : ℓ.natAbs + M + 2 * j ≤ m) :
    trace ((↑(U ^ ℓ) : A) * S ^ j * productPhase U S M M γ ^ m) =
      algebraMap ℚ K (mixedBetaMoment M (γ * m) 0 *
        (if (M + 2 * j : ℕ) ≤ ℓ then
          ((ascPochhammer ℚ j).eval ((M : ℚ) / 2) /
            (ascPochhammer ℚ j).eval ((M : ℚ) / 2 + (γ * m : ℕ) + 1)) *
          ((Ring.choose ((M * m : ℕ) - (M + 2 * j : ℕ) : ℤ)
            (ℓ - (M + 2 * j : ℕ)).toNat : ℤ) : ℚ) else 0)) := by
  have hm := Nat.mul_le_mul_right m (show 1 ≤ M by omega)
  rw [one_mul] at hm
  have hprod := Nat.mul_le_mul_right m (show M + 1 ≤ γ by omega)
  rw [Nat.add_mul, one_mul] at hprod
  have hu : 0 ≤ ((M * m : ℕ) - (M + 2 * j : ℕ) : ℤ) := by omega
  have hl : ℓ ≤ (M * m : ℕ) := by omega
  rw [product_mixedMoment trace U S M M M γ m j hM ℓ (by omega) weight height, if_pos hl]
  split_ifs with hcut
  · rw [mixedBetaMoment_height_factor M (γ * m) j hM, mul_assoc]
    let n := ((M * m : ℕ) - (M + 2 * j : ℕ) : ℤ).toNat
    have hn : (n : ℤ) = ((M * m : ℕ) - (M + 2 * j : ℕ) : ℤ) := Int.toNat_of_nonneg hu
    rw [← hn, Ring.choose_natCast, Ring.choose_natCast,
      Nat.choose_symm_of_eq_add (show n = ((M * m : ℕ) - ℓ : ℤ).toNat +
        (ℓ - (M + 2 * j : ℕ)).toNat by omega)]
  · rw [choose_eq_zero_of_nonneg_lt hu (by omega)]
    simp

end MathieuQuadratic
