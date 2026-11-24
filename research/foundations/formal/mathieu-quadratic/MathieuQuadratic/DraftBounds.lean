import MathieuQuadratic.BlockPolynomials
import MathieuQuadratic.PoleOrder
import Mathlib.LinearAlgebra.Finsupp.LinearCombination

/-!
# The polynomial degree statements on the quadratic complement
-/

open MvPolynomial
open scoped BigOperators

namespace MathieuQuadratic

variable {K : Type*} [Field K] [CharZero K] {M : ℕ}
variable (a : Fin 2 → K) (b : Fin M → K)

local notation "R" => QuadPolynomial K M
local notation "A" => NormComplement a b

omit [CharZero K] in
private theorem heightNorm_totalDegree_le : (heightNorm b).totalDegree ≤ 2 := by
  rw [heightNorm]
  apply totalDegree_finsetSum_le
  intro i _
  exact (totalDegree_mul _ _).trans (by simp [totalDegree_X_pow])

omit [CharZero K] in
private theorem quadNorm_totalDegree_le : (quadNorm a b).totalDegree ≤ 2 := by
  rw [quadNorm]
  apply (totalDegree_add _ _).trans
  refine max_le ?_ (heightNorm_totalDegree_le b)
  rw [planeNorm]
  apply totalDegree_finsetSum_le
  intro i _
  exact (totalDegree_mul _ _).trans (by simp [totalDegree_X_pow])

omit [CharZero K] in
private theorem heightNorm_image :
    algebraMap R A (heightNorm b) = (normUnit a b : A) * height a b := by
  rw [height, mul_left_comm, Units.mul_inv, mul_one]

/-- The polynomial part of the block average theorem, as an input
common to the two detecting-degree results. -/
def HasPolynomialBlockAverages (trace : A →ₗ[K] K) (phase : A) : Prop :=
  ∀ D : ℕ, ∀ p : R, p.totalDegree ≤ D →
    ∃ w ∈ blockMonomialSpan (K := K) (normUnit a b) (height a b) (D / 2),
      ∀ m : ℕ, 0 < m → trace (algebraMap R A p * phase ^ m) = trace (w * phase ^ m)

omit [CharZero K] in
private theorem shiftedBlockSpan_has_expansion {J n : ℕ} {w : A}
    (hw : w ∈ blockMonomialSpan (K := K) (normUnit a b) (height a b) J) :
    ∃ (support : Finset (ℤ × ℕ)) (c : ℤ × ℕ → K),
      (↑((normUnit a b)⁻¹) : A) ^ n * w =
        ∑ p ∈ support, c p • ((↑(normUnit a b ^ p.1) : A) * height a b ^ p.2) := by
  let U := normUnit a b
  let S := height a b
  let basis (p : ℤ × ℕ) : A := (↑(U ^ p.1) : A) * S ^ p.2
  let W := Submodule.span K (Set.range basis)
  let shift := LinearMap.mulLeft K ((↑(U⁻¹) : A) ^ n)
  have hsub : blockMonomialSpan (K := K) U S J ≤ W.comap shift := by
    apply Submodule.span_le.2
    rintro _ ⟨ℓ, j, _, _, rfl⟩
    change (↑(U⁻¹) : A) ^ n * ((U : A) ^ ℓ * S ^ j) ∈ W
    have he : (↑(U⁻¹) : A) ^ n * ((U : A) ^ ℓ * S ^ j) = basis ((ℓ : ℤ) - n, j) := by
      change _ = (↑(U ^ ((ℓ : ℤ) - n)) : A) * S ^ j
      rw [zpow_sub, zpow_natCast, zpow_natCast, ← inv_pow,
        Units.val_mul, Units.val_pow_eq_pow_val, Units.val_pow_eq_pow_val]
      ring
    rw [he]
    exact Submodule.subset_span ⟨_, rfl⟩
  obtain ⟨c, hc⟩ := Finsupp.mem_span_range_iff_exists_finsupp.1 (hsub hw)
  exact ⟨c.support, c, hc.symm⟩

omit [CharZero K] in
/-- Clearing the norm denominator reduces every regular multiplier to
the polynomial block averaging statement. -/
theorem draft_block_average_regular (trace : A →ₗ[K] K) (phase h : A)
    (average : ∀ r : ℕ, HasPolynomialBlockAverages a b
      (trace.comp (LinearMap.mulLeft K ((↑((normUnit a b)⁻¹) : A) ^ r))) phase) :
    ∃ (support : Finset (ℤ × ℕ)) (c : ℤ × ℕ → K),
      ∀ m : ℕ, 0 < m → trace (h * phase ^ m) =
        trace ((∑ p ∈ support, c p • ((↑(normUnit a b ^ p.1) : A) * height a b ^ p.2)) *
          phase ^ m) := by
  let p := (IsLocalization.Away.sec (quadNorm a b) h).1
  let n := (IsLocalization.Away.sec (quadNorm a b) h).2
  let U := normUnit a b
  have hp : h * (U : A) ^ n = algebraMap R A p := by
    rw [normUnit_val, ← map_pow]
    exact IsLocalization.Away.sec_spec (quadNorm a b) h
  have hclear : (↑(U⁻¹) : A) ^ n * algebraMap R A p = h := by
    rw [← hp, mul_left_comm, ← mul_pow, Units.inv_mul, one_pow, mul_one]
  obtain ⟨w, hw, havg⟩ := average n p.totalDegree p le_rfl
  obtain ⟨support, c, hc⟩ := shiftedBlockSpan_has_expansion a b hw
  refine ⟨support, c, fun m hm ↦ ?_⟩
  have he := havg m hm
  change trace ((↑(U⁻¹) : A) ^ n * (algebraMap R A p * phase ^ m)) =
    trace ((↑(U⁻¹) : A) ^ n * (w * phase ^ m)) at he
  rwa [← mul_assoc, ← mul_assoc, hclear, hc] at he

/-- The sharp least positive cutoff for all polynomial multipliers of
degree at most D, given the block averaging proposition. -/
theorem draft_cutoff_polynomial {trace : A →ₗ[K] K} (h : GeometricTrace a b trace)
    (hM : 0 < M) (γ : ℕ) (hγ : M ≤ γ)
    (average : HasPolynomialBlockAverages a b trace (complementProductPhase a b M M γ))
    (D : ℕ) (hD : D < 2 * M) (m₀ : ℕ) (hm₀ : 0 < m₀) :
    (∀ p : R, p.totalDegree ≤ D → ∀ m ≥ m₀,
      trace (algebraMap R A p * complementProductPhase a b M M γ ^ m) = 0) ↔
      (2 * (D / 2) ⌈/⌉ M) + 1 ≤ m₀ := by
  let J := D / 2
  have hJ : J < M := by omega
  have ht := h.complementTrace.normalizer a b hM
  constructor
  · intro hall
    by_contra! he
    let m := 2 * J ⌈/⌉ M
    change m₀ < m + 1 at he
    have hm : m₀ ≤ m := by omega
    have hJ₀ : 0 < J := by
      by_contra! hh
      have hz : m = 0 := by simp [m, show J = 0 by omega, zero_ceilDiv]
      omega
    have hp : ((heightNorm b) ^ J).totalDegree ≤ D := by
      have hdeg := totalDegree_pow (heightNorm b) J
      have hbase := heightNorm_totalDegree_le b
      exact hdeg.trans ((Nat.mul_le_mul_left J hbase).trans (by omega))
    have hz := hall _ hp m hm
    rw [map_pow, heightNorm_image a b, mul_pow] at hz
    exact product_cutoff_sharp trace _ _ hM hγ ht.weight ht.height hJ₀ hJ hz
  · intro hc p hp m hm
    have hmpos : 0 < m := hm₀.trans_le hm
    obtain ⟨w, hw, havg⟩ := average D p hp
    rw [havg m hmpos]
    exact (product_cutoff_span trace _ _ hM hγ ht.weight ht.height hJ hm₀).2 hc w hw m hm

/-- The least detecting polynomial degree for the current product witness. -/
theorem draft_product_polynomial_detects_iff {trace : A →ₗ[K] K} (h : GeometricTrace a b trace)
    (hM : 0 < M) (γ : ℕ) (hγ : M ≤ γ)
    (average : HasPolynomialBlockAverages a b trace (complementProductPhase a b M M γ))
    (D : ℕ) :
    (∃ p : R, p.totalDegree ≤ D ∧ ¬ ∃ m₀ : ℕ, ∀ m ≥ m₀,
      trace (algebraMap R A p * complementProductPhase a b M M γ ^ m) = 0) ↔ 2 * M ≤ D := by
  constructor
  · rintro ⟨p, hp, hz⟩
    by_contra! hd
    let m₀ := (2 * (D / 2) ⌈/⌉ M) + 1
    exact hz ⟨m₀, (draft_cutoff_polynomial a b h hM γ hγ average D hd m₀ (by omega)).2 le_rfl p hp⟩
  · intro hd
    let p := quadNorm a b ^ M
    have hp : p.totalDegree ≤ D :=
      (totalDegree_pow _ _).trans ((Nat.mul_le_mul_left M (quadNorm_totalDegree_le a b)).trans (by simpa [mul_comm] using hd))
    refine ⟨p, hp, ?_⟩
    rintro ⟨m₀, hm⟩
    have hz := hm (m₀ + 1) (by omega)
    change trace (algebraMap R A (quadNorm a b ^ M) * complementProductPhase a b M M γ ^ (m₀ + 1)) = 0 at hz
    rw [map_pow, ← normUnit_val] at hz
    exact (draft_product_witness a b h hM γ hγ).1.2 _ (by omega) hz

/-- The least detecting polynomial degree for the perturbed witness. -/
theorem draft_perturbation_polynomial_detects_iff {trace : A →ₗ[K] K} (h : GeometricTrace a b trace)
    (hM : 4 ≤ M)
    (average : HasPolynomialBlockAverages a b trace (perturbedProductPhase (normUnit a b) (height a b) M))
    (D : ℕ) :
    (∃ p : R, p.totalDegree ≤ D ∧ ¬ ∃ m₀ : ℕ, ∀ m ≥ m₀,
      trace (algebraMap R A p * perturbedProductPhase (normUnit a b) (height a b) M ^ m) = 0) ↔
      2 * M ≤ D := by
  have ht := h.complementTrace.normalizer a b (by omega)
  constructor
  · rintro ⟨p, hp, hz⟩
    by_contra! hd
    obtain ⟨w, hw, havg⟩ := average D p hp
    refine hz ⟨M + D, fun m hm ↦ ?_⟩
    rw [havg m (by omega)]
    apply perturbation_cutoff_span trace _ _ hM ht.weight ht.height (show D / 2 < M by omega) hw
    have he := Nat.mul_le_mul_right m (show 1 ≤ M - 2 by omega)
    rw [one_mul] at he
    omega
  · intro hd
    let p := quadNorm a b ^ M
    have hp : p.totalDegree ≤ D :=
      (totalDegree_pow _ _).trans ((Nat.mul_le_mul_left M (quadNorm_totalDegree_le a b)).trans (by simpa [mul_comm] using hd))
    refine ⟨p, hp, ?_⟩
    rintro ⟨m₀, hm⟩
    have hz := hm (m₀ + 1) (by omega)
    change trace (algebraMap R A (quadNorm a b ^ M) *
      perturbedProductPhase (normUnit a b) (height a b) M ^ (m₀ + 1)) = 0 at hz
    rw [map_pow, ← normUnit_val] at hz
    exact (perturbation_allOrderWitness trace _ _ hM ht.weight ht.height).2 _ (by omega) hz

end MathieuQuadratic
