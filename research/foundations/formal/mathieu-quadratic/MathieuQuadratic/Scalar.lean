import MathieuQuadratic.Complement

/-!
# Scalar weights on the localized coordinate ring

In characteristic zero, invariance under the scalar two suffices to
annihilate every nonzero norm weight in the uniform moment calculation.
-/

open scoped BigOperators

namespace MathieuQuadratic

open MvPolynomial

variable {K : Type*} [Field K] {M : ℕ}
variable (a : Fin 2 → K) (b : Fin M → K)

private noncomputable def scalePolynomial (c : Kˣ) :
    QuadPolynomial K M →ₐ[K] QuadPolynomial K M :=
  aeval fun i ↦ C (c : K) * X i

private theorem scalePolynomial_plane (c : Kˣ) :
    scalePolynomial (M := M) c (planeNorm a) = C (c : K) ^ 2 * planeNorm a := by
  simp only [planeNorm, scalePolynomial, map_sum, map_mul, map_pow,
    aeval_C, algebraMap_eq, aeval_X, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  ring

private theorem scalePolynomial_height (c : Kˣ) :
    scalePolynomial c (heightNorm b) = C (c : K) ^ 2 * heightNorm b := by
  simp only [heightNorm, scalePolynomial, map_sum, map_mul, map_pow,
    aeval_C, algebraMap_eq, aeval_X, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  ring

private theorem scalePolynomial_norm (c : Kˣ) :
    scalePolynomial c (quadNorm a b) = C (c : K) ^ 2 * quadNorm a b := by
  change _ = C (c : K) ^ 2 * (planeNorm a + heightNorm b)
  simp only [quadNorm, map_add, scalePolynomial_plane, scalePolynomial_height, mul_add]

private theorem scalar_units (c : Kˣ) :
    IsUnit ((IsScalarTower.toAlgHom K (QuadPolynomial K M) (NormComplement a b)).comp
      (scalePolynomial c) (quadNorm a b)) := by
  simp only [AlgHom.comp_apply, scalePolynomial_norm, map_mul, map_pow, IsScalarTower.toAlgHom_apply]
  exact ((c.isUnit.map C).map (algebraMap _ _)).pow 2 |>.mul (IsLocalization.Away.algebraMap_isUnit _)

/-- Pullback by multiplication of every coordinate by the same unit. -/
noncomputable def scalarPullback (c : Kˣ) : NormComplement a b →ₐ[K] NormComplement a b :=
  IsLocalization.Away.liftAlgHom (quadNorm a b) (scalar_units a b c)

private theorem scalarPullback_polynomial (c : Kˣ) (p : QuadPolynomial K M) :
    scalarPullback a b c (algebraMap (QuadPolynomial K M) (NormComplement a b) p) =
      algebraMap (QuadPolynomial K M) (NormComplement a b) (scalePolynomial c p) :=
  IsLocalization.Away.lift_eq _ (scalar_units a b c) _

private theorem scalarPullback_unit (c : Kˣ) :
    Units.map (scalarPullback a b c).toMonoidHom (normUnit a b) =
      Units.map (algebraMap K (NormComplement a b)).toMonoidHom c ^ 2 * normUnit a b := by
  apply Units.ext
  change scalarPullback a b c (normUnit a b : NormComplement a b) =
    algebraMap K (NormComplement a b) (c : K) ^ 2 * (normUnit a b : NormComplement a b)
  rw [normUnit_val, scalarPullback_polynomial, scalePolynomial_norm, map_mul, map_pow,
    IsScalarTower.algebraMap_apply K (QuadPolynomial K M) (NormComplement a b)]
  rfl

private theorem scalarPullback_height (c : Kˣ) :
    scalarPullback a b c (height a b) = height a b := by
  let U := normUnit a b
  let C₀ := Units.map (algebraMap K (NormComplement a b)).toMonoidHom c
  have hi : scalarPullback a b c (↑(U⁻¹) : NormComplement a b) =
      (↑(U⁻¹ * (C₀ ^ 2)⁻¹) : NormComplement a b) := by
    change (↑(Units.map (scalarPullback a b c).toMonoidHom (U⁻¹)) : NormComplement a b) = _
    rw [map_inv, scalarPullback_unit, mul_inv_rev]
  unfold height
  rw [map_mul, scalarPullback_polynomial, scalePolynomial_height, map_mul, map_pow, hi]
  have cc : algebraMap (QuadPolynomial K M) (NormComplement a b) (C (c : K)) = (C₀ : NormComplement a b) := by
    change _ = algebraMap K _ (c : K)
    rw [IsScalarTower.algebraMap_apply K (QuadPolynomial K M) (NormComplement a b)]
    rfl
  rw [cc, Units.val_mul, ← Units.val_pow_eq_pow_val]
  calc
    _ = ((↑(C₀ ^ 2) : NormComplement a b) * ↑((C₀ ^ 2)⁻¹)) *
        (algebraMap (QuadPolynomial K M) (NormComplement a b) (heightNorm b) * ↑(U⁻¹)) := by ring
    _ = _ := by rw [Units.mul_inv, one_mul]

section CharZero

variable [CharZero K]

private theorem two_weight_ne_one (d : ℤ) (hd : d ≠ 0) :
    (↑((Units.mk0 (2 : K) two_ne_zero) ^ (2 * d)) : K) ≠ 1 := by
  rw [Units.val_zpow_eq_zpow_val, Units.val_mk0]
  intro he
  have hmap : algebraMap ℚ K ((2 : ℚ) ^ (2 * d)) = algebraMap ℚ K 1 := by simpa using he
  have hr := (RingHom.injective (algebraMap ℚ K)) hmap
  have hz := (zpow_right_strictMono₀ (by norm_num : (1 : ℚ) < 2)).injective
    (by simpa using hr : (2 : ℚ) ^ (2 * d) = 2 ^ (0 : ℤ))
  exact hd (by omega)

theorem trace_weight_of_scalar (trace : NormComplement a b →ₗ[K] K)
    (inv : ∀ f, trace (scalarPullback a b (Units.mk0 (2 : K) two_ne_zero) f) = trace f)
    (d : ℤ) (j : ℕ) (hd : d ≠ 0) :
    trace ((↑(normUnit a b ^ d) : NormComplement a b) * height a b ^ j) = 0 := by
  let c : Kˣ := Units.mk0 (2 : K) two_ne_zero
  let U := normUnit a b
  let V : NormComplement a b := (↑(U ^ d) : NormComplement a b) * height a b ^ j
  have he : scalarPullback a b c V = algebraMap K _ (↑(c ^ (2 * d)) : K) * V := by
    dsimp only [V]
    rw [map_mul, map_pow, scalarPullback_height]
    change (↑(Units.map (scalarPullback a b c).toMonoidHom (U ^ d)) : NormComplement a b) * _ = _
    rw [map_zpow, scalarPullback_unit, mul_zpow, ← zpow_natCast _ 2, ← zpow_mul,
      Units.val_mul, mul_assoc]
    have hc : (↑(Units.map (algebraMap K (NormComplement a b)).toMonoidHom c ^ (2 * d)) :
          NormComplement a b) = algebraMap K _ (↑(c ^ (2 * d)) : K) := by
      rw [← map_zpow]
      rfl
    exact congrArg (fun z ↦ z * V) hc
  have hs : trace (algebraMap K _ (↑(c ^ (2 * d)) : K) * V) = trace V := by
    rw [← he]
    exact inv V
  rw [← Algebra.smul_def, map_smul, smul_eq_mul] at hs
  have hz : ((↑(c ^ (2 * d)) : K) - 1) * trace V = 0 := by
    rw [sub_mul, one_mul, hs, sub_self]
  exact (mul_eq_zero.mp hz).resolve_left (sub_ne_zero.mpr (two_weight_ne_one d hd))

end CharZero

end MathieuQuadratic
