import MathieuQuadratic.InvariantConstants

/-!
# The action on regular functions

The orthogonal and scalar actions on the nonzero-norm locus induce
pullbacks on the actual localization. Over an algebraically closed
field, their fixed subspace is exactly the constants.
-/

open scoped BigOperators Matrix

namespace MathieuQuadratic

open MvPolynomial

/-- Polynomial pullback by a matrix, over any commutative coefficient ring. -/
noncomputable def matrixPolynomialPullback {R I J : Type*} [CommRing R] [Fintype J]
    (T : Matrix I J R) : MvPolynomial I R →ₐ[R] MvPolynomial J R :=
  aeval fun i ↦ ∑ j, C (T i j) * X j

@[simp] theorem matrixPolynomialPullback_X {R I J : Type*} [CommRing R] [Fintype J]
    (T : Matrix I J R) (i : I) :
    matrixPolynomialPullback T (X i) = ∑ j, C (T i j) * X j :=
  aeval_X _ _

variable {K J : Type*} [Field K] [Fintype J]

theorem eval_matrixPolynomialPullback {I : Type*} (T : Matrix I J K) (v : J → K) (p : MvPolynomial I K) :
    eval v (matrixPolynomialPullback T p) = eval (T *ᵥ v) p := by
  have he : (aeval v : MvPolynomial J K →ₐ[K] K).comp (matrixPolynomialPullback T) =
      aeval (T *ᵥ v) := by
    apply MvPolynomial.algHom_ext
    intro i
    simp [Matrix.mulVec, dotProduct]
  exact AlgHom.congr_fun he p

private noncomputable def linearPullback (T : Matrix J J K) (c : Kˣ) :
    MvPolynomial J K →ₐ[K] MvPolynomial J K :=
  aeval (fun i ↦ C (c : K) * ∑ j, C (T i j) * X j)

private theorem eval_linearPullback (T : Matrix J J K) (c : Kˣ) (v : J → K)
    (p : MvPolynomial J K) :
    eval v (linearPullback T c p) = eval (c • (T *ᵥ v)) p := by
  have he : (aeval v : MvPolynomial J K →ₐ[K] K).comp (linearPullback T c) =
      aeval (c • (T *ᵥ v)) := by
    apply MvPolynomial.algHom_ext
    intro i
    simp [linearPullback, Matrix.mulVec, dotProduct, Units.smul_def, smul_eq_mul]
  exact AlgHom.congr_fun he p

private theorem linearPullback_norm [Infinite K] (d : J → Kˣ) (T : Matrix J J K) (c : Kˣ)
    (hT : ∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v) :
    linearPullback T c (diagonalNormPolynomial d) =
      C ((c : K) ^ 2) * diagonalNormPolynomial d := by
  apply MvPolynomial.funext
  intro v
  rw [eval_linearPullback, eval_diagonalNormPolynomial, Units.smul_def, diagonalNorm_smul, hT,
    map_mul, eval_C, eval_diagonalNormPolynomial]

local notation "R" => MvPolynomial J K
variable (d : J → Kˣ) [Infinite K]

local notation "Cq" => Localization.Away (diagonalNormPolynomial d)

private theorem linearPullback_norm_unit (T : Matrix J J K) (c : Kˣ)
    (hT : ∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v) :
    IsUnit (algebraMap R Cq (linearPullback T c (diagonalNormPolynomial d))) := by
  rw [linearPullback_norm d T c hT, map_mul]
  apply IsUnit.mul
  · have he : algebraMap R Cq (C ((c : K) ^ 2)) = algebraMap K Cq ((c : K) ^ 2) :=
      (IsScalarTower.algebraMap_apply K R Cq _).symm
    rw [he]
    exact (c.isUnit.pow 2).map (algebraMap K Cq)
  · exact IsLocalization.Away.algebraMap_isUnit _

/-- The pullback by a scalar followed by a norm-preserving linear map. -/
noncomputable def geometricPullback (T : Matrix J J K) (c : Kˣ)
    (hT : ∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v) : Cq →ₐ[K] Cq :=
  IsLocalization.Away.liftAlgHom (diagonalNormPolynomial d)
    (f := (IsScalarTower.toAlgHom K R Cq).comp (linearPullback T c))
    (linearPullback_norm_unit d T c hT)

private theorem geometricPullback_polynomial (T : Matrix J J K) (c : Kˣ)
    (hT : ∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v) (p : R) :
    geometricPullback d T c hT (algebraMap R Cq p) = algebraMap R Cq (linearPullback T c p) :=
  IsLocalization.Away.lift_eq (diagonalNormPolynomial d)
    (g := (algebraMap R Cq).comp (linearPullback T c).toRingHom)
    (linearPullback_norm_unit d T c hT) p

omit [Infinite K] in
private theorem linearPullback_eq_matrixPullback (T : Matrix J J K) (c : Kˣ) :
    linearPullback T c = matrixPolynomialPullback ((c : K) • T) := by
  apply MvPolynomial.algHom_ext
  intro i
  simp only [linearPullback, aeval_X, matrixPolynomialPullback_X, Matrix.smul_apply,
    smul_eq_mul, map_mul, Finset.mul_sum, mul_assoc]

theorem scaledMatrixPullback_norm (T : Matrix J J K) (c : Kˣ)
    (hT : ∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v) :
    matrixPolynomialPullback ((c : K) • T) (diagonalNormPolynomial d) =
      C ((c : K) ^ 2) * diagonalNormPolynomial d := by
  rw [← linearPullback_eq_matrixPullback]
  exact linearPullback_norm d T c hT

theorem geometricPullback_polynomial_eq_matrixPullback (T : Matrix J J K) (c : Kˣ)
    (hT : ∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v) (p : R) :
    geometricPullback d T c hT (algebraMap R Cq p) =
      algebraMap R Cq (matrixPolynomialPullback ((c : K) • T) p) := by
  rw [geometricPullback_polynomial, linearPullback_eq_matrixPullback]

/-- After scalar-isometry pullback, clearing the same denominator
multiplies the pullback by the scalar weight of that denominator. -/
theorem geometricPullback_clear_denominator (T : Matrix J J K) (c : Kˣ)
    (hT : ∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v)
    (f : Cq) (p : R) (n : ℕ)
    (hp : f * algebraMap R Cq (diagonalNormPolynomial d ^ n) = algebraMap R Cq p) :
    ((c : K) ^ (2 * n) • geometricPullback d T c hT f) *
        algebraMap R Cq (diagonalNormPolynomial d ^ n) =
      algebraMap R Cq (matrixPolynomialPullback ((c : K) • T) p) := by
  let q := diagonalNormPolynomial d
  let ι := algebraMap R Cq
  let ψ := geometricPullback d T c hT
  have hψq : ψ (ι q) = algebraMap K Cq ((c : K) ^ 2) * ι q := by
    rw [show ψ (ι q) = ι (matrixPolynomialPullback ((c : K) • T) q) from
      geometricPullback_polynomial_eq_matrixPullback d T c hT q,
      scaledMatrixPullback_norm d T c hT, map_mul,
      ← algebraMap_eq, ← IsScalarTower.algebraMap_apply]
  have hs := congrArg ψ hp
  change ψ (f * ι (q ^ n)) = ψ (ι p) at hs
  rw [map_mul, map_pow, map_pow, hψq, mul_pow, ← map_pow, ← pow_mul, mul_left_comm,
    show ψ (ι p) = ι (matrixPolynomialPullback ((c : K) • T) p) from
      geometricPullback_polynomial_eq_matrixPullback d T c hT p] at hs
  change ((c : K) ^ (2 * n) • ψ f) * ι (q ^ n) = _
  rw [Algebra.smul_mul_assoc, Algebra.smul_def, ← hs]
  simp only [map_pow]

theorem eval_geometricPullback (T : Matrix J J K) (c : Kˣ)
    (hT : ∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v)
    (v : J → K) (hv : diagonalNorm d v ≠ 0) (hw : diagonalNorm d (c • (T *ᵥ v)) ≠ 0) (f : Cq) :
    evalOnComplement d v hv (geometricPullback d T c hT f) =
      evalOnComplement d (c • (T *ᵥ v)) hw f := by
  have h : (evalOnComplement d v hv).comp (geometricPullback d T c hT).toRingHom =
      evalOnComplement d (c • (T *ᵥ v)) hw := by
    apply IsLocalization.ringHom_ext (Submonoid.powers (diagonalNormPolynomial d))
    apply RingHom.ext
    intro p
    change evalOnComplement d v hv (geometricPullback d T c hT (algebraMap R Cq p)) =
      evalOnComplement d (c • (T *ᵥ v)) hw (algebraMap R Cq p)
    rw [geometricPullback_polynomial, evalOnComplement_polynomial, evalOnComplement_polynomial,
      eval_linearPullback]
  exact DFunLike.congr_fun h f

/-- Over an algebraic closure the regular functions fixed by the special
orthogonal and scalar pullbacks are exactly the constants. -/
theorem geometric_fixed_iff_constant [DecidableEq J] [Nonempty J] [CharZero K] [IsAlgClosed K]
    (f : Cq) :
    (∀ T : Matrix J J K, ∀ c : Kˣ, T.det = 1 →
      ∀ hT : ∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v,
        geometricPullback d T c hT f = f) ↔
      ∃! k : K, f = algebraMap K Cq k := by
  constructor
  · intro hf
    apply invariantRegularFunction_isConstant d f
    intro T c hd hT x y hx hy he
    have hy' : diagonalNorm d (c • (T *ᵥ x)) ≠ 0 := by
      rw [he]
      exact hy
    have h := eval_geometricPullback d T c hT x hx hy' f
    rw [hf T c hd hT] at h
    exact h.trans (by cases he; rfl)
  · rintro ⟨k, h, -⟩ T c _ hT
    rw [h, AlgHom.commutes]

end MathieuQuadratic
