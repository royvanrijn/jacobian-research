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

variable {K J : Type*} [Field K] [Fintype J]

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
