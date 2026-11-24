import MathieuQuadratic.Transitivity
import Mathlib.Algebra.MvPolynomial.Funext
import Mathlib.RingTheory.Localization.Away.Basic

/-!
# Constants on the geometric quadratic complement

Geometric transitivity in Section 2 of the draft implies that the
invariant regular functions are constants. Equality at the geometric
points is lifted to the localized coordinate ring by polynomial
evaluation over the infinite ground field.
-/

open scoped BigOperators Matrix

namespace MathieuQuadratic

open MvPolynomial

variable {K J : Type*} [Field K] [Fintype J]

/-- The norm polynomial of a nondegenerate diagonal form. -/
noncomputable def diagonalNormPolynomial (d : J → Kˣ) : MvPolynomial J K :=
  ∑ i, C (d i : K) * X i ^ 2

@[simp] theorem eval_diagonalNormPolynomial (d : J → Kˣ) (v : J → K) :
    eval v (diagonalNormPolynomial d) = diagonalNorm d v := by
  simp [diagonalNormPolynomial, diagonalNorm, diagonalPair, pow_two, mul_assoc]

variable [DecidableEq J]

private theorem diagonalNorm_single (d : J → Kˣ) (i : J) :
    diagonalNorm d (Pi.single i 1) = (d i : K) := by
  rw [diagonalNorm, diagonalPair_single_right, Pi.single_eq_same, mul_one]

private theorem diagonalNormPolynomial_ne_zero [Nonempty J] (d : J → Kˣ) :
    diagonalNormPolynomial d ≠ 0 := by
  obtain ⟨i⟩ := ‹Nonempty J›
  intro h
  have he := congrArg (eval (Pi.single i (1 : K))) h
  rw [eval_diagonalNormPolynomial, diagonalNorm_single, map_zero] at he
  exact Units.ne_zero (d i) he

/-- Evaluation at a point on the nonzero-norm locus. -/
noncomputable def evalOnComplement (d : J → Kˣ) (v : J → K)
    (hv : diagonalNorm d v ≠ 0) : Localization.Away (diagonalNormPolynomial d) →+* K :=
  IsLocalization.Away.lift (diagonalNormPolynomial d) (g := eval v)
    (isUnit_iff_ne_zero.mpr (by simpa only [eval_diagonalNormPolynomial] using hv))

omit [DecidableEq J] in
@[simp] theorem evalOnComplement_polynomial (d : J → Kˣ) (v : J → K)
    (hv : diagonalNorm d v ≠ 0) (p : MvPolynomial J K) :
    evalOnComplement d v hv
      (algebraMap (MvPolynomial J K) (Localization.Away (diagonalNormPolynomial d)) p) = eval v p := by
  unfold evalOnComplement
  rw [IsLocalization.Away.lift_eq]

omit [DecidableEq J] in
private theorem one_dimensional_scalar_transitive [Subsingleton J] [Nonempty J]
    (d : J → Kˣ) (x y : J → K) (hx : diagonalNorm d x ≠ 0) (hy : diagonalNorm d y ≠ 0) :
    ∃ c : Kˣ, c • x = y := by
  obtain ⟨i⟩ := ‹Nonempty J›
  have hx' : x i ≠ 0 := by
    intro h
    have hzero (j : J) : x j = 0 := (congrArg x (Subsingleton.elim j i)).trans h
    exact hx (by simp [diagonalNorm, diagonalPair, hzero])
  have hy' : y i ≠ 0 := by
    intro h
    have hzero (j : J) : y j = 0 := (congrArg y (Subsingleton.elim j i)).trans h
    exact hy (by simp [diagonalNorm, diagonalPair, hzero])
  refine ⟨Units.mk0 (y i / x i) (div_ne_zero hy' hx'), ?_⟩
  funext j
  rw [Subsingleton.elim j i]
  exact div_mul_cancel₀ _ hx'

private theorem eq_algebraMap_of_constant_evaluations [Infinite K] [Nonempty J]
    (d : J → Kˣ) (f : Localization.Away (diagonalNormPolynomial d)) (c : K)
    (hf : ∀ v (hv : diagonalNorm d v ≠ 0), evalOnComplement d v hv f = c) :
    f = algebraMap K _ c := by
  let q := diagonalNormPolynomial d
  let p := (IsLocalization.Away.sec q f).1
  let n := (IsLocalization.Away.sec q f).2
  let i : MvPolynomial J K →+* Localization.Away q := algebraMap _ _
  have hfrac : f * i (q ^ n) = i p := IsLocalization.Away.sec_spec q f
  have hp : p = C c * q ^ n := by
    have hz : q * (p - C c * q ^ n) = 0 := by
      apply MvPolynomial.funext
      intro v
      by_cases hv : diagonalNorm d v = 0
      · simp [q, hv]
      · have he := congrArg (evalOnComplement d v hv) hfrac
        change evalOnComplement d v hv (f * i (q ^ n)) = _ at he
        rw [map_mul, hf, show evalOnComplement d v hv (i (q ^ n)) = eval v (q ^ n) from
          evalOnComplement_polynomial d v hv _, show evalOnComplement d v hv (i p) = eval v p from
          evalOnComplement_polynomial d v hv _] at he
        simp [← he]
    exact sub_eq_zero.mp ((mul_eq_zero.mp hz).resolve_left (diagonalNormPolynomial_ne_zero d))
  refine ((IsLocalization.Away.algebraMap_isUnit (S := Localization.Away q) q).pow n).mul_right_cancel ?_
  rw [← map_pow, hfrac, hp, map_mul]
  congr 1

/-- Regular functions invariant under geometric scalar and special
orthogonal points are constant. In positive dimension every constant is
the value at a nonzero-norm point, so the scalar is unique. -/
theorem invariantRegularFunction_isConstant [CharZero K] [Infinite K] [Nonempty J] [IsAlgClosed K]
    (d : J → Kˣ) (f : Localization.Away (diagonalNormPolynomial d))
    (invariant : ∀ (T : Matrix J J K) (c : Kˣ), T.det = 1 →
      (∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v) →
      ∀ (x y : J → K) (hx : diagonalNorm d x ≠ 0) (hy : diagonalNorm d y ≠ 0),
        c • (T *ᵥ x) = y → evalOnComplement d x hx f = evalOnComplement d y hy f) :
    ∃! c : K, f = algebraMap K _ c := by
  obtain ⟨j⟩ := ‹Nonempty J›
  let y : J → K := Pi.single j (1 : K)
  have hy : diagonalNorm d y ≠ 0 := by
    rw [show y = Pi.single j (1 : K) from rfl, diagonalNorm_single]
    exact Units.ne_zero _
  let c := evalOnComplement d y hy f
  refine ⟨c, eq_algebraMap_of_constant_evaluations d f c ?_, ?_⟩
  swap
  · intro x hx
    dsimp only [c]
    rw [hx, IsScalarTower.algebraMap_apply K (MvPolynomial J K)
      (Localization.Away (diagonalNormPolynomial d)), evalOnComplement_polynomial]
    simp
  intro x hx
  by_cases hr : 1 < Fintype.card J
  · obtain ⟨T, c, hd, hnorm, he⟩ := exists_isometry_and_scalar d x y hr hx hy
    exact invariant T c hd hnorm x y hx hy he
  · letI : Subsingleton J := Fintype.card_le_one_iff_subsingleton.mp (not_lt.mp hr)
    obtain ⟨c, hc⟩ := one_dimensional_scalar_transitive d x y hx hy
    apply invariant (1 : Matrix J J K) c Matrix.det_one (by simp) x y hx hy
    simpa using hc

end MathieuQuadratic
