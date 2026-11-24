import MathieuQuadratic.GeometricProjection
import Mathlib.Algebra.MvPolynomial.Coeff

/-!
# Uniqueness of the geometric trace

Test a linear functional on powers of linear forms. The multinomial
coefficient of each parameter monomial recovers its value on the
corresponding coordinate monomial. On the weight-zero component of
the norm localization, the resulting parameter polynomial divided by
the corresponding power of the dual norm is invariant.
-/

open MvPolynomial
open scoped BigOperators Matrix

namespace MathieuQuadratic

variable {K J : Type*} [Field K] [CharZero K] [Fintype J] [DecidableEq J]

local notation "P" => MvPolynomial J K

private noncomputable def testPolynomial (L : P →ₗ[K] K) (p : MvPolynomial J P) : P :=
  AddMonoidAlgebra.ofCoeff ((AddMonoidAlgebra.coeff p).mapRange L (map_zero L))

omit [CharZero K] [Fintype J] [DecidableEq J] in
private theorem coeff_testPolynomial (L : P →ₗ[K] K) (p : MvPolynomial J P) (m : J →₀ ℕ) :
    coeff m (testPolynomial L p) = L (coeff m p) := rfl

omit [CharZero K] [Fintype J] [DecidableEq J] in
private theorem eval_testPolynomial (L : P →ₗ[K] K) (p : MvPolynomial J P) (x : J → K) :
    eval x (testPolynomial L p) = L (eval (fun i ↦ C (x i)) p) := by
  classical
  induction p using MvPolynomial.induction_on' with
  | add p r hp hr =>
    have h : testPolynomial L (p + r) = testPolynomial L p + testPolynomial L r := by
      apply MvPolynomial.ext
      intro m
      rw [coeff_testPolynomial, coeff_add, coeff_add, coeff_testPolynomial, coeff_testPolynomial, map_add]
    rw [h, map_add, hp, hr, map_add, map_add]
  | monomial m a =>
    have h : testPolynomial L (monomial m a) = monomial m (L a) := by
      apply MvPolynomial.ext
      intro r
      rw [coeff_testPolynomial]
      by_cases h : m = r <;> simp [coeff_monomial, h]
    rw [h, eval_monomial, eval_monomial]
    simp only [← map_pow, Finsupp.prod, ← map_prod,
      ← algebraMap_eq]
    rw [mul_comm a, ← Algebra.smul_def, map_smul, smul_eq_mul, mul_comm]

private noncomputable def powerTest (L : P →ₗ[K] K) (n : ℕ) : P :=
  testPolynomial L ((∑ i : J, (X i : P) • (X i : MvPolynomial J P)) ^ n)

private noncomputable def linearForm (x : J → K) : P := ∑ i, x i • X i

omit [CharZero K] [DecidableEq J] in
private theorem eval_powerTest (L : P →ₗ[K] K) (n : ℕ) (x : J → K) :
    eval x (powerTest L n) = L (linearForm x ^ n) := by
  simp only [powerTest, eval_testPolynomial, map_pow, map_sum, linearForm, smul_eq_C_mul,
    map_mul, eval_C, eval_X, mul_comm]

omit [CharZero K] [DecidableEq J] in
private theorem coeff_powerTest (L : P →ₗ[K] K) (n : ℕ) (m : J →₀ ℕ) (hm : m.degree = n) :
    coeff m (powerTest L n) = (m.multinomial : K) * L (monomial m 1) := by
  rw [powerTest, coeff_testPolynomial, coeff_linearCombination_X_pow_of_fintype,
    if_pos (show m.sum (fun _ z ↦ z) = n from hm), ← C_eq_coe_nat, ← smul_eq_C_mul,
    map_smul, smul_eq_mul]
  rw [monomial_eq, C_1, one_mul]

omit [DecidableEq J] in
private theorem powerTest_eq_on_monomial {L H : P →ₗ[K] K} {n : ℕ}
    (h : powerTest L n = powerTest H n) (m : J →₀ ℕ) (hm : m.degree = n) (a : K) :
    L (monomial m a) = H (monomial m a) := by
  have he := congrArg (coeff m) h
  rw [coeff_powerTest L n m hm, coeff_powerTest H n m hm] at he
  have hne : (m.multinomial : K) ≠ 0 := by
    exact_mod_cast Nat.ne_of_gt (Nat.multinomial_pos m.support m)
  have heq := mul_left_cancel₀ hne he
  rw [show monomial m a = a • monomial m 1 by rw [← LinearMap.map_smul, smul_eq_mul, mul_one],
    map_smul, map_smul, heq]

section Geometric

variable [Nonempty J] [IsAlgClosed K]

local notation "Q" => diagonalNormPolynomial
local notation "Loc[" d "]" => Localization.Away (diagonalNormPolynomial d)
local notation "incl[" d "]" => algebraMap P (Loc[d])

private noncomputable def normDenominator (d : J → Kˣ) : (Loc[d])ˣ :=
  (IsLocalization.Away.algebraMap_isUnit (S := Loc[d]) (Q d)).unit

omit [CharZero K] [DecidableEq J] [Nonempty J] [IsAlgClosed K] in
private theorem normDenominator_val (d : J → Kˣ) :
    (normDenominator d : Loc[d]) = incl[d] (Q d) := IsUnit.unit_spec _

private noncomputable def normFraction (d : J → Kˣ) (n : ℕ) (p : P) : Loc[d] :=
  incl[d] p * (↑((normDenominator d)⁻¹) : Loc[d]) ^ n

omit [CharZero K] [DecidableEq J] [Nonempty J] [IsAlgClosed K] in
private theorem normFraction_cancel (d : J → Kˣ) (n : ℕ) (p : P) :
    normFraction d n p * incl[d] (Q d ^ n) = incl[d] p := by
  rw [normFraction, map_pow, ← normDenominator_val, mul_assoc, ← mul_pow,
    Units.inv_mul, one_pow, mul_one]

omit [CharZero K] [DecidableEq J] [Nonempty J] [IsAlgClosed K] in
private theorem eq_normFraction (d : J → Kˣ) (f : Loc[d]) (p : P) (n : ℕ)
    (h : f * incl[d] (Q d ^ n) = incl[d] p) :
    f = normFraction d n p := by
  apply ((IsLocalization.Away.algebraMap_isUnit (S := Loc[d]) (Q d)).pow n).mul_right_cancel
  rw [← map_pow, h, normFraction_cancel]

private noncomputable def traceNumerator (d : J → Kˣ) (L : Loc[d] →ₗ[K] K) (n : ℕ) : P →ₗ[K] K :=
  L.comp ((LinearMap.mulRight K ((↑((normDenominator d)⁻¹) : Loc[d]) ^ n)).comp
    (IsScalarTower.toAlgHom K P (Loc[d])).toLinearMap)

omit [CharZero K] [DecidableEq J] [Nonempty J] [IsAlgClosed K] in
private theorem traceNumerator_apply (d : J → Kˣ) (L : Loc[d] →ₗ[K] K) (n : ℕ) (p : P) :
    traceNumerator d L n p = L (normFraction d n p) := rfl

omit [CharZero K] [DecidableEq J] [Nonempty J] [IsAlgClosed K] in
private theorem traceNumerator_norm_pow (d : J → Kˣ) (L : Loc[d] →ₗ[K] K) (n : ℕ) :
    traceNumerator d L n (Q d ^ n) = L 1 := by
  rw [traceNumerator_apply, ← eq_normFraction d 1 _ n (one_mul _)]

omit [CharZero K] [Nonempty J] in
private theorem traceNumerator_covariance (d : J → Kˣ) (L : Loc[d] →ₗ[K] K)
    (inv : ∀ g f, L (orthogonalScalarAction d g f) = L f)
    (T : Matrix J J K) (c : Kˣ) (hdet : T.det = 1)
    (hT : ∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v) (n : ℕ) (p : P) :
    traceNumerator d L n (matrixPolynomialPullback ((c : K) • T) p) =
      (c : K) ^ (2 * n) * traceNumerator d L n p := by
  have hs := geometricPullback_clear_denominator d T c hT (normFraction d n p) p n
    (normFraction_cancel d n p)
  have he := eq_normFraction d _ _ n hs
  simp only [traceNumerator_apply, ← he, map_smul, smul_eq_mul]
  rw [show L (geometricPullback d T c hT (normFraction d n p)) = L (normFraction d n p) from
    inv ⟨T, c, hdet, hT⟩ _]

omit [CharZero K] [Nonempty J] in
private theorem matrixPullback_one (p : P) : matrixPolynomialPullback (1 : Matrix J J K) p = p := by
  apply MvPolynomial.funext
  intro v
  rw [eval_matrixPolynomialPullback, Matrix.one_mulVec]

omit [Nonempty J] in
private theorem traceNumerator_monomial_unbalanced (d : J → Kˣ) (L : Loc[d] →ₗ[K] K)
    (inv : ∀ g f, L (orthogonalScalarAction d g f) = L f)
    (n : ℕ) (m : J →₀ ℕ) (a : K) (hm : m.degree ≠ 2 * n) :
    traceNumerator d L n (monomial m a) = 0 := by
  let c : Kˣ := Units.mk0 2 two_ne_zero
  have hv := traceNumerator_covariance d L inv (1 : Matrix J J K) c Matrix.det_one
    (by simp) n (monomial m a)
  rw [scaled_matrixPullback_monomial, matrixPullback_one, map_smul, smul_eq_mul] at hv
  have hne : (c : K) ^ (m.degree : ℕ) ≠ (c : K) ^ (2 * n) := by
    intro h
    have he : algebraMap ℚ K ((2 : ℚ) ^ (m.degree : ℕ)) = algebraMap ℚ K ((2 : ℚ) ^ (2 * n)) := by
      simpa [c] using h
    exact hm ((pow_right_strictMono₀ (by norm_num : (1 : ℚ) < 2)).injective
      ((RingHom.injective (algebraMap ℚ K)) he))
  have hz : (((c : K) ^ (m.degree : ℕ)) - (c : K) ^ (2 * n)) *
      traceNumerator d L n (monomial m a) = 0 := by
    rw [sub_mul, hv, sub_self]
  exact (mul_eq_zero.mp hz).resolve_left (sub_ne_zero.mpr hne)

omit [CharZero K] [DecidableEq J] [Nonempty J] [IsAlgClosed K] in
private theorem eval_linearForm (x y : J → K) : eval y (linearForm x) = x ⬝ᵥ y := by
  simp [linearForm, dotProduct, smul_eq_C_mul]

omit [CharZero K] [DecidableEq J] [Nonempty J] in
private theorem matrixPullback_linearForm (T : Matrix J J K) (c : K) (x : J → K) :
    matrixPolynomialPullback (c • T) (linearForm x) = linearForm (c • (Tᵀ *ᵥ x)) := by
  apply MvPolynomial.funext
  intro y
  rw [eval_matrixPolynomialPullback, eval_linearForm, eval_linearForm, Matrix.smul_mulVec,
    dotProduct_smul, smul_dotProduct, Matrix.dotProduct_mulVec,
    ← Matrix.vecMul_transpose Tᵀ x, Matrix.transpose_transpose]

omit [Nonempty J] [IsAlgClosed K] in
private theorem inverseTranspose_preserves (d : J → Kˣ) (T : Matrix J J K) (hdet : T.det = 1)
    (hT : ∀ x, diagonalNorm (fun i ↦ (d i)⁻¹) (T *ᵥ x) = diagonalNorm (fun i ↦ (d i)⁻¹) x) :
    ∀ x, diagonalNorm d (Tᵀ *ᵥ x) = diagonalNorm d x := by
  have h := inverseDiagonal_orthogonal (fun i ↦ (d i)⁻¹) T
    (by rw [hdet]; exact isUnit_one) (diagonalConjugate_of_preserves _ T hT)
  have hh : (Tᵀ)ᵀ * Matrix.diagonal (fun i ↦ (d i : K)) * Tᵀ = Matrix.diagonal (fun i ↦ (d i : K)) := by
    simpa only [inv_inv, Matrix.transpose_transpose] using h
  intro x
  have he := congrArg (eval x) (matrixPullback_norm_of_conjugate d Tᵀ hh)
  rwa [eval_matrixPolynomialPullback, eval_diagonalNormPolynomial, eval_diagonalNormPolynomial] at he

omit [CharZero K] [DecidableEq J] [Nonempty J] [IsAlgClosed K] in
private theorem eval_normFraction [Infinite K] (d : J → Kˣ) (p : P) (n : ℕ)
    (x : J → K) (hx : diagonalNorm d x ≠ 0) :
    evalOnComplement d x hx (normFraction d n p) = eval x p / diagonalNorm d x ^ n := by
  apply (eq_div_iff (pow_ne_zero _ hx)).2
  have h := congrArg (evalOnComplement d x hx) (normFraction_cancel d n p)
  simpa only [map_mul, evalOnComplement_polynomial, map_pow, eval_diagonalNormPolynomial] using h

private theorem powerTest_proportional_norm (d : J → Kˣ) (L : Loc[d] →ₗ[K] K)
    (inv : ∀ g f, L (orthogonalScalarAction d g f) = L f) (n : ℕ) :
    ∃ c : K, powerTest (traceNumerator d L n) (2 * n) =
      C c * (Q (fun j ↦ (d j)⁻¹)) ^ n := by
  let e : J → Kˣ := fun j ↦ (d j)⁻¹
  let B := powerTest (traceNumerator d L n) (2 * n)
  let f := normFraction e n B
  have hi (T : Matrix J J K) (c : Kˣ) (hdet : T.det = 1)
      (hT : ∀ v, diagonalNorm e (T *ᵥ v) = diagonalNorm e v)
      (x y : J → K) (hx : diagonalNorm e x ≠ 0) (hy : diagonalNorm e y ≠ 0)
      (hxy : c • (T *ᵥ x) = y) :
      evalOnComplement e x hx f = evalOnComplement e y hy f := by
    have ht := traceNumerator_covariance d L inv Tᵀ c
      (by rwa [Matrix.det_transpose]) (inverseTranspose_preserves d T hdet hT) n
      (linearForm x ^ (2 * n))
    rw [map_pow, matrixPullback_linearForm, Matrix.transpose_transpose, show (c : K) • (T *ᵥ x) = y from hxy] at ht
    rw [eval_normFraction e B n x hx, eval_normFraction e B n y hy]
    change eval x (powerTest _ _) / _ = eval y (powerTest _ _) / _
    rw [eval_powerTest, eval_powerTest, ht, ← hxy, Units.smul_def, diagonalNorm_smul, hT,
      mul_pow, ← pow_mul, mul_div_mul_left _ _ (pow_ne_zero _ c.ne_zero)]
  obtain ⟨c, hc, _⟩ := invariantRegularFunction_isConstant e f hi
  refine ⟨c, ?_⟩
  have h := congrArg (fun w : Loc[e] ↦ w * incl[e] (Q e ^ n)) hc
  change normFraction e n B * incl[e] (Q e ^ n) = _ at h
  rw [normFraction_cancel, IsScalarTower.algebraMap_apply K P (Loc[e]),
    algebraMap_eq, ← map_mul] at h
  have hinj : Function.Injective (incl[e]) := IsLocalization.injective (Loc[e])
    (Submonoid.powers_le.mpr (mem_nonZeroDivisors_iff_ne_zero.mpr (diagonalNormPolynomial_ne_zero e)))
  exact hinj h

omit [CharZero K] [DecidableEq J] [Nonempty J] in
private theorem powerTest_smul (L : P →ₗ[K] K) (c : K) (n : ℕ) :
    powerTest (c • L) n = C c * powerTest L n := by
  apply MvPolynomial.funext
  intro x
  rw [eval_powerTest, map_mul, eval_C, eval_powerTest, LinearMap.smul_apply, smul_eq_mul]

omit [DecidableEq J] [Nonempty J] [IsAlgClosed K] in
private theorem powerTest_eq_on_homogeneous (L H : P →ₗ[K] K) {n : ℕ}
    (h : powerTest L n = powerTest H n) (p : P) (hp : p.IsHomogeneous n) : L p = H p := by
  rw [p.as_sum, map_sum, map_sum]
  apply Finset.sum_congr rfl
  intro m hm
  apply powerTest_eq_on_monomial h m _ (coeff m p)
  rw [Finsupp.degree_eq_weight_one]
  exact hp (mem_support_iff.mp hm)

omit [DecidableEq J] [Nonempty J] [CharZero K] in
private theorem powerTest_zero (n : ℕ) : powerTest (0 : P →ₗ[K] K) n = 0 := by
  apply MvPolynomial.funext
  intro x
  rw [eval_powerTest, LinearMap.zero_apply, map_zero]

omit [CharZero K] [DecidableEq J] [Nonempty J] [IsAlgClosed K] in
private theorem norm_pow_homogeneous (d : J → Kˣ) (n : ℕ) : (Q d ^ n).IsHomogeneous (2 * n) := by
  have h : (Q d).IsHomogeneous 2 := by
    unfold diagonalNormPolynomial
    apply Submodule.sum_mem (homogeneousSubmodule J K 2)
    intro i _
    exact ((isHomogeneous_X K i).pow 2).C_mul (d i : K)
  exact h.pow n

private theorem traceNumerator_eq (d : J → Kˣ) (L H : Loc[d] →ₗ[K] K)
    (invL : ∀ g f, L (orthogonalScalarAction d g f) = L f)
    (invH : ∀ g f, H (orthogonalScalarAction d g f) = H f)
    (hL₁ : L 1 = 1) (hH₁ : H 1 = 1) (n : ℕ) :
    traceNumerator d L n = traceNumerator d H n := by
  let F := traceNumerator d L n
  let G := traceNumerator d H n
  let e : J → Kˣ := fun i ↦ (d i)⁻¹
  obtain ⟨c, hc⟩ := powerTest_proportional_norm d L invL n
  obtain ⟨a, ha⟩ := powerTest_proportional_norm d H invH n
  change powerTest F (2 * n) = C c * Q e ^ n at hc
  change powerTest G (2 * n) = C a * Q e ^ n at ha
  have ha0 : a ≠ 0 := by
    intro he
    have hz : powerTest G (2 * n) = powerTest (0 : P →ₗ[K] K) (2 * n) := by
      rw [ha, he, map_zero, zero_mul, powerTest_zero]
    have hp := powerTest_eq_on_homogeneous G 0 hz (Q d ^ n) (norm_pow_homogeneous d n)
    change traceNumerator d H n (Q d ^ n) = 0 at hp
    rw [traceNumerator_norm_pow, hH₁] at hp
    exact one_ne_zero hp
  have hmatch : powerTest (a • F) (2 * n) = powerTest (c • G) (2 * n) := by
    rw [powerTest_smul, powerTest_smul, hc, ha, mul_left_comm]
  have hac := powerTest_eq_on_homogeneous _ _ hmatch (Q d ^ n) (norm_pow_homogeneous d n)
  change a * traceNumerator d L n (Q d ^ n) = c * traceNumerator d H n (Q d ^ n) at hac
  rw [traceNumerator_norm_pow, traceNumerator_norm_pow, hL₁, hH₁, mul_one, mul_one] at hac
  rw [← hac] at hmatch
  apply LinearMap.ext
  intro p
  induction p using MvPolynomial.induction_on' with
  | add p r hp hr => rw [map_add, map_add, hp, hr]
  | monomial m r =>
    by_cases hm : m.degree = 2 * n
    · have h := powerTest_eq_on_monomial hmatch m hm r
      exact mul_left_cancel₀ ha0 h
    · rw [traceNumerator_monomial_unbalanced d L invL n m r hm,
        traceNumerator_monomial_unbalanced d H invH n m r hm]

/-- A normalized geometric invariant functional is unique. The proof
only uses transitivity and the coefficients of powers of linear forms. -/
theorem normalized_geometricPointTrace_unique (d : J → Kˣ) (L H : Loc[d] →ₗ[K] K)
    (invL : ∀ g f, L (orthogonalScalarAction d g f) = L f)
    (invH : ∀ g f, H (orthogonalScalarAction d g f) = H f)
    (hL₁ : L 1 = 1) (hH₁ : H 1 = 1) : L = H := by
  apply LinearMap.ext
  intro f
  let s := IsLocalization.Away.sec (Q d) f
  have hf := eq_normFraction d f s.1 s.2 (IsLocalization.Away.sec_spec (Q d) f)
  rw [hf]
  exact LinearMap.congr_fun (traceNumerator_eq d L H invL invH hL₁ hH₁ s.2) s.1

/-- The algebraic formula is the unique normalized trace on geometric
points over an algebraically closed ground field. -/
theorem existsUnique_geometricPointTrace (d : J → Kˣ) :
    ∃! trace : Loc[d] →ₗ[K] K, trace 1 = 1 ∧ ∀ g f,
      trace (orthogonalScalarAction d g f) = trace f := by
  have h := algebraicTrace_normalized_invariant d
  exact ⟨algebraicTrace d, h, fun trace ht ↦
    normalized_geometricPointTrace_unique d trace _ ht.2 h.2 ht.1 h.1⟩

end Geometric

end MathieuQuadratic
