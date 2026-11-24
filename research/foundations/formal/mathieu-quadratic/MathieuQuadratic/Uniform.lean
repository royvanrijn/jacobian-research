import MathieuQuadratic.Moments

/-!
# Uniform consequences of the mixed moments

This proves the pure-power threshold, both fixed radial marker thresholds
and the all-order witness, from the two trace identities used in the proof.
-/

open scoped BigOperators

namespace MathieuQuadratic

variable {K A : Type*} [Field K] [CharZero K] [CommRing A] [Algebra K A]

/-- The two trace identities to be obtained from the scalar action and
rank-two rotations. -/
structure NormalizerTrace (M : ℕ) (trace : A →ₗ[K] K) (U : Aˣ) (S : A) : Prop where
  weight : ∀ (a : ℤ) (j : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ j) = 0
  height : ∀ j : ℕ, trace (S ^ j) = (M : K) / (↑(M + 2 * j) : K)

variable {M : ℕ} {trace : A →ₗ[K] K} {U : Aˣ} {S : A}

theorem normalizer_marked_nonzero_iff (h : NormalizerTrace M trace U S)
    (hM : 0 < M) {k : ℕ} (hk : M ≤ k) {ℓ : ℤ} (hlower : -(k : ℤ) ≤ ℓ) :
    trace ((↑(U ^ ℓ) : A) * normalizerPhase U S ^ k) ≠ 0 ↔
      (M : ℤ) ≤ ℓ ∧ ℓ ≤ (k : ℤ) := by
  rw [normalizer_mixedMoment trace U S M k hM ℓ hlower h.weight h.height]
  by_cases hℓ : ℓ ≤ (k : ℤ)
  · have hd : (k : ℤ) - M = ((k - M : ℕ) : ℤ) := by omega
    simp only [hℓ, if_true, hd, Ring.choose_natCast, Int.cast_natCast, ne_eq,
      map_eq_zero_iff _ (RingHom.injective (algebraMap ℚ K)), mul_eq_zero,
      (betaMoment_pos M k hM).ne', false_or, Nat.cast_eq_zero, Nat.choose_ne_zero_iff, and_true]
    omega
  · simp [hℓ]

/-- The signed binomial evaluation of the first moment below the
pure-power threshold, as in equation (15). -/
theorem normalizer_unmarked_eq_of_lt (h : NormalizerTrace M trace U S)
    {k : ℕ} (hk : k < M) :
    trace (normalizerPhase U S ^ k) =
      algebraMap ℚ K (betaMoment M k * ((-1 : ℚ) ^ k * ((M - 1).choose k : ℚ))) := by
  have formula := normalizer_mixedMoment trace U S M k (by omega) 0 (by omega) h.weight h.height
  rw [zpow_zero, Units.val_one, one_mul, if_pos (Int.natCast_nonneg k), sub_zero,
    Int.toNat_natCast] at formula
  have hn : (k : ℤ) - M = -((M : ℤ) - k) := by omega
  have hi : (M : ℤ) - k + k - 1 = ((M - 1 : ℕ) : ℤ) := by omega
  rw [formula, hn, Ring.choose_neg, hi, Ring.choose_natCast, Units.smul_def, smul_eq_mul,
    Int.cast_mul, Int.cast_negOnePow_natCast, Int.cast_natCast]

theorem normalizer_unmarked_ne_zero_of_lt (h : NormalizerTrace M trace U S)
    {k : ℕ} (hk : k < M) :
    trace (normalizerPhase U S ^ k) ≠ 0 := by
  rw [normalizer_unmarked_eq_of_lt h hk, ne_eq,
    map_eq_zero_iff _ (RingHom.injective (algebraMap ℚ K))]
  exact mul_ne_zero (betaMoment_pos M k (by omega)).ne' (mul_ne_zero (pow_ne_zero _ (by norm_num))
    (by exact_mod_cast Nat.choose_ne_zero (by omega : k ≤ M - 1)))

theorem normalizer_pure_pow_iff (h : NormalizerTrace M trace U S)
    (hM : 0 < M) {a : ℕ} (_ha : 0 < a) :
    Pure trace (normalizerPhase U S ^ a) ↔ M ≤ a := by
  constructor
  · intro hp
    by_contra hn
    have hf := hp 1 Nat.one_pos
    rw [pow_one] at hf
    exact normalizer_unmarked_ne_zero_of_lt h (by omega) hf
  · intro hMa m hm
    rw [← pow_mul]
    have hk : M ≤ a * m := hMa.trans (Nat.le_mul_of_pos_right a hm)
    by_contra hn
    have marked := (normalizer_marked_nonzero_iff h hM hk (ℓ := 0) (by omega)).1
      (by simpa only [zpow_zero, Units.val_one, one_mul] using hn)
    omega

private theorem marker_large_order (h : NormalizerTrace M trace U S)
    (hM : 0 < M) {a : ℕ} (hMa : M ≤ a) (ℓ : ℤ) {m : ℕ}
    (hm : ℓ.natAbs + 1 ≤ m) :
    trace ((↑(U ^ ℓ) : A) * (normalizerPhase U S ^ a) ^ m) ≠ 0 ↔ (M : ℤ) ≤ ℓ := by
  have ha : 0 < a := hM.trans_le hMa
  have hm' : m ≤ a * m := by
    simpa using Nat.mul_le_mul_right m ha
  rw [← pow_mul, normalizer_marked_nonzero_iff h hM (hMa.trans
    (Nat.le_mul_of_pos_right a (by omega))) (by omega), and_iff_left (by omega)]

theorem normalizer_marker_eventually_nonzero_iff (h : NormalizerTrace M trace U S)
    (hM : 0 < M) {a : ℕ} (hMa : M ≤ a) (ℓ : ℤ) :
    (∃ cutoff : ℕ, ∀ m ≥ cutoff,
      trace ((↑(U ^ ℓ) : A) * (normalizerPhase U S ^ a) ^ m) ≠ 0) ↔ (M : ℤ) ≤ ℓ := by
  constructor
  · rintro ⟨cutoff, hc⟩
    exact (marker_large_order h hM hMa ℓ (le_max_right cutoff (ℓ.natAbs + 1))).1
      (hc _ (le_max_left _ _))
  · intro hℓ
    exact ⟨ℓ.natAbs + 1, fun m hm ↦ (marker_large_order h hM hMa ℓ hm).2 hℓ⟩

theorem normalizer_marker_eventually_zero_iff (h : NormalizerTrace M trace U S)
    (hM : 0 < M) {a : ℕ} (hMa : M ≤ a) (ℓ : ℤ) :
    (∃ cutoff : ℕ, ∀ m ≥ cutoff,
      trace ((↑(U ^ ℓ) : A) * (normalizerPhase U S ^ a) ^ m) = 0) ↔ ℓ < (M : ℤ) := by
  constructor
  · rintro ⟨cutoff, hc⟩
    apply lt_of_not_ge
    intro hℓ
    exact (marker_large_order h hM hMa ℓ (le_max_right cutoff (ℓ.natAbs + 1))).2 hℓ
      (hc _ (le_max_left _ _))
  · intro hℓ
    refine ⟨ℓ.natAbs + 1, fun m hm ↦ ?_⟩
    by_contra hne
    exact (not_le_of_gt hℓ) ((marker_large_order h hM hMa ℓ hm).1 hne)

theorem normalizer_marker_all_orders_iff (h : NormalizerTrace M trace U S)
    (hM : 0 < M) {a : ℕ} (hMa : M ≤ a) (ℓ : ℤ) :
    (∀ m : ℕ, 0 < m → trace ((↑(U ^ ℓ) : A) * (normalizerPhase U S ^ a) ^ m) ≠ 0) ↔
      (M : ℤ) ≤ ℓ ∧ ℓ ≤ (a : ℤ) := by
  constructor
  · intro hm
    have hl := (marker_large_order h hM hMa ℓ (le_refl (ℓ.natAbs + 1))).1
      (hm _ (Nat.zero_lt_succ _))
    have hfirst := hm 1 Nat.one_pos
    rw [pow_one, normalizer_mixedMoment trace U S M a hM ℓ (by omega) h.weight h.height] at hfirst
    exact ⟨hl, by
      by_contra hn
      exact hfirst (if_neg hn)⟩
  · rintro ⟨hl, hu⟩ m hm
    have ham := Nat.le_mul_of_pos_right a hm
    rw [← pow_mul]
    exact (normalizer_marked_nonzero_iff h hM (hMa.trans ham) (by omega)).2 ⟨hl, by omega⟩

theorem normalizer_marked_normPower (h : NormalizerTrace M trace U S) (hM : 0 < M)
    {k : ℕ} (hk : M ≤ k) :
    trace ((U : A) ^ M * normalizerPhase U S ^ k) = algebraMap ℚ K (betaMoment M k) := by
  have hm := normalizer_mixedMoment trace U S M k hM (M : ℤ) (by omega) h.weight h.height
  have hd : (k : ℤ) - (M : ℤ) = ((k - M : ℕ) : ℤ) := by omega
  simpa [show (M : ℤ) ≤ (k : ℤ) by omega, hd, Ring.choose_natCast] using hm

/-- The two exact moments in (6), including the nonzero marked sum. -/
theorem uniformPhase_allOrderWitness (h : NormalizerTrace M trace U S) (hM : 0 < M) :
    AllOrderWitness trace (uniformPhase M (normalizerPhase U S)) ((U : A) ^ M) ∧
      ∀ m : ℕ, 0 < m →
        trace ((U : A) ^ M * uniformPhase M (normalizerPhase U S) ^ m) =
          algebraMap ℚ K (∑ r ∈ Finset.range (m + 1),
            (m.choose r : ℚ) * betaMoment M (M * m + r)) := by
  apply uniformPhase_allOrderWitness_of_moments
  · intro k _
    have j := betaMoment_pos M k hM
    exact j
  · intro k hk
    have hp := (normalizer_pure_pow_iff h hM (a := k) (hM.trans_le hk)).2 hk
    simpa using hp 1 Nat.one_pos
  · exact fun _ hk ↦ normalizer_marked_normPower h hM hk

theorem normalizer_direct_ternary (h : NormalizerTrace 1 trace U S) :
    AllOrderWitness trace (normalizerPhase U S) (U : A) ∧
      ∀ m : ℕ, 0 < m →
        trace ((U : A) * normalizerPhase U S ^ m) = algebraMap ℚ K (betaMoment 1 m) := by
  have hformula (m : ℕ) (hm : 0 < m) :
      trace ((U : A) * normalizerPhase U S ^ m) = algebraMap ℚ K (betaMoment 1 m) := by
    rw [← pow_one (U : A)]
    exact normalizer_marked_normPower h Nat.one_pos hm
  refine ⟨⟨?_, fun m hm ↦ ?_⟩, hformula⟩
  · simpa using (normalizer_pure_pow_iff h Nat.one_pos (a := 1) Nat.one_pos).2 le_rfl
  · rw [hformula m hm, ne_eq, map_eq_zero_iff _ (RingHom.injective (algebraMap ℚ K))]
    exact (betaMoment_pos 1 m Nat.one_pos).ne'

end MathieuQuadratic
