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

/-- A coordinate has scalar weight one. -/
theorem scalarPullback_X (c : Kˣ) (i : Fin 2 ⊕ Fin M) :
    scalarPullback a b c (algebraMap (QuadPolynomial K M) (NormComplement a b) (X i)) =
      algebraMap K (NormComplement a b) (c : K) *
        algebraMap (QuadPolynomial K M) (NormComplement a b) (X i) := by
  rw [scalarPullback_polynomial, scalePolynomial, aeval_X, map_mul,
    ← algebraMap_eq, ← IsScalarTower.algebraMap_apply]

/-- The plane norm has scalar weight two. -/
theorem scalarPullback_plane (c : Kˣ) :
    scalarPullback a b c (algebraMap (QuadPolynomial K M) (NormComplement a b) (planeNorm a)) =
      (algebraMap K (NormComplement a b) (c : K)) ^ 2 *
        algebraMap (QuadPolynomial K M) (NormComplement a b) (planeNorm a) := by
  rw [scalarPullback_polynomial, scalePolynomial_plane, map_mul, map_pow,
    ← algebraMap_eq, ← IsScalarTower.algebraMap_apply]

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

private theorem trace_of_scalar_eigenvector (trace : NormComplement a b →ₗ[K] K)
    (inv : ∀ f, trace (scalarPullback a b (Units.mk0 (2 : K) two_ne_zero) f) = trace f)
    (m d : ℕ) (hd : d ≠ 2 * m) (p : NormComplement a b)
    (hp : scalarPullback a b (Units.mk0 (2 : K) two_ne_zero) p =
      algebraMap K (NormComplement a b) ((2 : K) ^ d) * p) :
    trace ((↑((normUnit a b)⁻¹) : NormComplement a b) ^ m * p) = 0 := by
  let Cq := NormComplement a b
  let c : Kˣ := Units.mk0 2 two_ne_zero
  let U := normUnit a b
  let φ := algebraMap K Cq
  let s := scalarPullback a b c
  let v := (↑(U⁻¹) : Cq) ^ m * p
  have hU : s (U : Cq) = φ (c : K) ^ 2 * (U : Cq) :=
    congrArg (fun u : Cqˣ ↦ (u : Cq)) (scalarPullback_unit a b c)
  have hinv : φ (c : K) ^ 2 * s (↑(U⁻¹) : Cq) = (↑(U⁻¹) : Cq) := by
    have h : (φ (c : K) ^ 2 * (U : Cq)) * s (↑(U⁻¹) : Cq) = 1 := by
      rw [← hU, ← map_mul, Units.mul_inv, map_one]
    calc
      _ = (((↑(U⁻¹) : Cq) * (U : Cq)) * (φ (c : K) ^ 2 * s (↑(U⁻¹) : Cq))) := by
        rw [Units.inv_mul, one_mul]
      _ = (↑(U⁻¹) : Cq) * ((φ (c : K) ^ 2 * (U : Cq)) * s (↑(U⁻¹) : Cq)) := by ring
      _ = _ := by rw [h, mul_one]
  have hs : φ ((c : K) ^ (2 * m)) * s v = φ ((c : K) ^ d) * v := by
    change φ ((c : K) ^ (2 * m)) * s ((↑(U⁻¹) : Cq) ^ m * p) =
      φ ((c : K) ^ d) * ((↑(U⁻¹) : Cq) ^ m * p)
    rw [map_mul, map_pow, show s p = φ ((c : K) ^ d) * p from hp,
      map_pow, pow_mul]
    calc
      _ = (φ (c : K) ^ 2 * s (↑(U⁻¹) : Cq)) ^ m * (φ ((c : K) ^ d) * p) := by
        rw [mul_pow, mul_assoc]
      _ = _ := by rw [hinv]; ring
  have htrace := congrArg trace hs
  change trace (algebraMap K Cq ((c : K) ^ (2 * m)) * s v) =
    trace (algebraMap K Cq ((c : K) ^ d) * v) at htrace
  simp only [← Algebra.smul_def, map_smul, smul_eq_mul, show trace (s v) = trace v from inv v] at htrace
  have hne : (c : K) ^ (2 * m) ≠ (c : K) ^ d := by
    intro he
    have hA : algebraMap ℚ K ((2 : ℚ) ^ (2 * m)) = algebraMap ℚ K ((2 : ℚ) ^ d) := by
      simpa [c] using he
    have hi := (RingHom.injective (algebraMap ℚ K)) hA
    exact hd ((pow_right_strictMono₀ (by norm_num : (1 : ℚ) < 2)).injective hi).symm
  have hzero : ((c : K) ^ (2 * m) - (c : K) ^ d) * trace v = 0 := by
    rw [sub_mul, htrace, sub_self]
  exact (mul_eq_zero.mp hzero).resolve_left (sub_ne_zero.mpr hne)

/-- Scalar weights for monomials in a single coordinate and the plane norm. -/
theorem trace_coordinate_plane_weight (trace : NormComplement a b →ₗ[K] K)
    (inv : ∀ f, trace (scalarPullback a b (Units.mk0 (2 : K) two_ne_zero) f) = trace f)
    (i : Fin 2 ⊕ Fin M) (m n j : ℕ) (hn : n + 2 * j ≠ 2 * m) :
    trace ((↑((normUnit a b)⁻¹) : NormComplement a b) ^ m *
      algebraMap (QuadPolynomial K M) (NormComplement a b) (X i) ^ n *
      algebraMap (QuadPolynomial K M) (NormComplement a b) (planeNorm a) ^ j) = 0 := by
  rw [mul_assoc]
  apply trace_of_scalar_eigenvector a b trace inv m (n + 2 * j) hn
  rw [map_mul, map_pow, map_pow, scalarPullback_X, scalarPullback_plane,
    map_pow, pow_add, pow_mul, mul_pow, mul_pow, Units.val_mk0]
  ring

/-- A monomial divided by q^m has zero trace away from total degree 2m. -/
theorem trace_monomial_weight (trace : NormComplement a b →ₗ[K] K)
    (inv : ∀ f, trace (scalarPullback a b (Units.mk0 (2 : K) two_ne_zero) f) = trace f)
    (m : ℕ) (d : (Fin 2 ⊕ Fin M) →₀ ℕ) (r : K)
    (hd : d.sum (fun _ k ↦ k) ≠ 2 * m) :
    trace ((↑((normUnit a b)⁻¹) : NormComplement a b) ^ m *
      algebraMap (QuadPolynomial K M) (NormComplement a b) (monomial d r)) = 0 := by
  apply trace_of_scalar_eigenvector a b trace inv m _ hd
  rw [scalarPullback_polynomial, scalePolynomial, aeval_monomial, monomial_eq]
  simp only [Finsupp.prod, Finsupp.sum, mul_pow, Finset.prod_mul_distrib,
    Finset.prod_pow_eq_pow_sum, map_mul, map_pow, ← algebraMap_eq,
    ← IsScalarTower.algebraMap_apply]
  ac_rfl


/-- Scalar invariance and normalization recover the constant coefficient
of a polynomial. -/
theorem trace_polynomial_eq_constantCoeff (trace : NormComplement a b →ₗ[K] K)
    (normalized : trace 1 = 1)
    (inv : ∀ f, trace (scalarPullback a b (Units.mk0 (2 : K) two_ne_zero) f) = trace f)
    (p : QuadPolynomial K M) :
    trace (algebraMap (QuadPolynomial K M) (NormComplement a b) p) = coeff 0 p := by
  conv_lhs => rw [p.as_sum]
  rw [map_sum, map_sum, Finset.sum_eq_single 0]
  · rw [← C_apply, ← algebraMap_eq, ← IsScalarTower.algebraMap_apply, ← mul_one
      (algebraMap K (NormComplement a b) (coeff 0 p)), ← Algebra.smul_def,
      map_smul, normalized, smul_eq_mul, mul_one]
  · intro d _ hd
    have hdeg : d.sum (fun _ k ↦ k) ≠ 2 * 0 := by
      simpa using (Finsupp.sum_pos (g := fun _ k ↦ k)
        (fun _ hi ↦ Nat.pos_of_ne_zero (Finsupp.mem_support_iff.mp hi)) hd).ne'
    simpa using trace_monomial_weight a b trace inv 0 d (coeff d p) hdeg
  · intro h
    rw [notMem_support_iff.mp h, map_zero, map_zero, map_zero]

end CharZero

end MathieuQuadratic
