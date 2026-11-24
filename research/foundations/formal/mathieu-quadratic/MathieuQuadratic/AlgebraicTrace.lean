import MathieuQuadratic.GeometricAction
import Mathlib.Algebra.MvPolynomial.PDeriv
import Mathlib.RingTheory.MvPolynomial.Homogeneous

/-!
# Algebraic moments for the normalized quadratic trace

The one-variable moments satisfy the Gaussian integration-by-parts
recursion over the coefficient ring. Their product gives a linear
functional on polynomials. To descend to the norm localization, keep
the part of degree twice the denominator exponent and normalize its
moment by the corresponding radial factor.
-/

open scoped BigOperators Matrix

namespace MathieuQuadratic

open MvPolynomial

section Polynomial

variable {R J : Type*} [CommRing R] [Fintype J]

/-- The formal moments of a coordinate whose inverse norm weight is e. -/
private def coordinateMoment (e : R) : ℕ → R
  | 0 => 1
  | 1 => 0
  | n + 2 => (n + 1 : ℕ) * e * coordinateMoment e n

private def monomialMoment (d : J → Rˣ) (α : J →₀ ℕ) : R :=
  ∏ i, coordinateMoment (↑((d i)⁻¹) : R) (α i)

/-- The product moment on polynomial monomials, extended linearly. -/
noncomputable def polynomialMoment (d : J → Rˣ) : MvPolynomial J R →ₗ[R] R :=
  (Finsupp.linearCombination R (monomialMoment d)).comp
    (AddMonoidAlgebra.coeffLinearEquiv R).toLinearMap

@[simp] theorem polynomialMoment_monomial (d : J → Rˣ) (α : J →₀ ℕ) (r : R) :
    polynomialMoment d (monomial α r) = r * monomialMoment d α := by
  change Finsupp.linearCombination R (monomialMoment d) (Finsupp.single α r) = _
  rw [Finsupp.linearCombination_single, smul_eq_mul]

@[simp] theorem polynomialMoment_one (d : J → Rˣ) :
    polynomialMoment d 1 = 1 := by
  rw [← C_1, C_apply, polynomialMoment_monomial]
  simp [monomialMoment, coordinateMoment]

private theorem monomialMoment_add_two [DecidableEq J] (d : J → Rˣ) (α : J →₀ ℕ) (i : J) :
    (d i : R) * monomialMoment d (α + Finsupp.single i 2) =
      (↑(α i) + 1 : R) * monomialMoment d α := by
  unfold monomialMoment
  rw [← Finset.mul_prod_erase Finset.univ _ (Finset.mem_univ i),
    ← Finset.mul_prod_erase Finset.univ (fun j ↦ coordinateMoment (↑((d j)⁻¹) : R) (α j))
      (Finset.mem_univ i)]
  have hprod : ∏ j ∈ Finset.univ.erase i,
      coordinateMoment (↑((d j)⁻¹) : R) ((α + Finsupp.single i 2 : J →₀ ℕ) j) =
      ∏ j ∈ Finset.univ.erase i, coordinateMoment (↑((d j)⁻¹) : R) (α j) := by
    apply Finset.prod_congr rfl
    intro j hj
    rw [Finsupp.add_apply, Finsupp.single_eq_of_ne (Finset.ne_of_mem_erase hj), add_zero]
  rw [hprod, Finsupp.add_apply, Finsupp.single_eq_same, coordinateMoment, Nat.cast_add, Nat.cast_one]
  calc
    _ = ((d i : R) * ↑((d i)⁻¹)) *
        ((↑(α i) + 1) * (coordinateMoment (↑((d i)⁻¹) : R) (α i) *
          ∏ j ∈ Finset.univ.erase i, coordinateMoment (↑((d j)⁻¹) : R) (α j))) := by ring
    _ = _ := by rw [Units.mul_inv, one_mul]

private theorem coordinateMoment_map {S : Type*} [CommRing S] (f : R →+* S) (e : R) (n : ℕ) :
    coordinateMoment (f e) n = f (coordinateMoment e n) := by
  induction n using Nat.strong_induction_on with
  | h n ih =>
    obtain _ | _ | n := n
    · simp [coordinateMoment]
    · simp [coordinateMoment]
    · rw [coordinateMoment, coordinateMoment, map_mul, map_mul, map_natCast, ih n (by omega)]

/-- Coefficient change commutes with taking polynomial moments. -/
theorem polynomialMoment_map {S : Type*} [CommRing S] (f : R →+* S)
    (d : J → Rˣ) (p : MvPolynomial J R) :
    polynomialMoment (Units.map f.toMonoidHom ∘ d) (MvPolynomial.map f p) =
      f (polynomialMoment d p) := by
  induction p using MvPolynomial.induction_on' with
  | add p q hp hq => rw [map_add, map_add, hp, hq, map_add, map_add]
  | monomial α r =>
    rw [map_monomial, polynomialMoment_monomial, polynomialMoment_monomial, map_mul]
    unfold monomialMoment
    rw [map_prod]
    congr 1
    apply Finset.prod_congr rfl
    intro i _
    change coordinateMoment (↑((Units.map f.toMonoidHom (d i))⁻¹) : S) (α i) = _
    rw [← map_inv (Units.map f.toMonoidHom)]
    exact coordinateMoment_map f _ _

variable [DecidableEq J]

/-- Algebraic integration by parts for the product moment. -/
theorem polynomialMoment_X_mul (d : J → Rˣ) (i : J) (p : MvPolynomial J R) :
    polynomialMoment d (X i * p) =
      (↑((d i)⁻¹) : R) * polynomialMoment d (pderiv i p) := by
  induction p using MvPolynomial.induction_on' with
  | add p q hp hq => rw [mul_add, map_add, hp, hq, map_add, map_add, mul_add]
  | monomial α r =>
    rw [X, monomial_mul, one_mul, pderiv_monomial, polynomialMoment_monomial,
      polynomialMoment_monomial, add_comm (Finsupp.single i 1) α]
    by_cases h : α i = 0
    · simp only [h, Nat.cast_zero, mul_zero, zero_mul]
      have hm : monomialMoment d (α + Finsupp.single i 1) = 0 := by
        unfold monomialMoment
        rw [← Finset.mul_prod_erase Finset.univ _ (Finset.mem_univ i),
          Finsupp.add_apply, Finsupp.single_eq_same, h, zero_add, coordinateMoment, zero_mul]
      rw [hm, mul_zero]
    · let β := α - Finsupp.single i 1
      have hle : Finsupp.single i 1 ≤ α := Finsupp.single_le_iff.mpr (Nat.pos_of_ne_zero h)
      have hb : β + Finsupp.single i 1 = α := tsub_add_cancel_of_le hle
      have hα : α + Finsupp.single i 1 = β + Finsupp.single i 2 := by
        calc
          _ = (β + Finsupp.single i 1) + Finsupp.single i 1 := by rw [hb]
          _ = _ := by rw [add_assoc, ← Finsupp.single_add]
      have hbi : (↑(β i) + 1 : R) = (↑(α i) : R) := by
        have he := congrArg (fun a : J →₀ ℕ ↦ a i) hb
        rw [Finsupp.add_apply, Finsupp.single_eq_same] at he
        rw [← Nat.cast_one, ← Nat.cast_add, he]
      rw [hα]
      have ht := monomialMoment_add_two d β i
      rw [hbi] at ht
      calc
        _ = r * ((↑((d i)⁻¹) : R) *
            ((d i : R) * monomialMoment d (β + Finsupp.single i 2))) := by
          rw [← mul_assoc (↑((d i)⁻¹) : R), Units.inv_mul, one_mul]
        _ = _ := by rw [ht]; ring

private theorem polynomialMoment_unique (d : J → Rˣ) (L : MvPolynomial J R →ₗ[R] R)
    (h₁ : L 1 = 1)
    (parts : ∀ (i : J) (p : MvPolynomial J R),
      L (X i * p) = (↑((d i)⁻¹) : R) * L (pderiv i p)) :
    L = polynomialMoment d := by
  have hm (n : ℕ) : ∀ (α : J →₀ ℕ), Finsupp.degree α = n →
      ∀ r : R, L (monomial α r) = polynomialMoment d (monomial α r) := by
    induction n using Nat.strong_induction_on with
    | h n ih =>
      intro α ha r
      by_cases hα : α = 0
      · subst α
        rw [← C_apply, ← algebraMap_eq, Algebra.algebraMap_eq_smul_one, map_smul, map_smul,
          h₁, polynomialMoment_one]
      · obtain ⟨i, hi⟩ := Finset.nonempty_iff_ne_empty.mpr (Finsupp.support_nonempty_iff.mpr hα).ne_empty
        let β := α - Finsupp.single i 1
        have hle : Finsupp.single i 1 ≤ α :=
          Finsupp.single_le_iff.mpr (Nat.pos_of_ne_zero (Finsupp.mem_support_iff.mp hi))
        have hb : β + Finsupp.single i 1 = α := tsub_add_cancel_of_le hle
        have hβ : Finsupp.degree β < n := by
          have he := congrArg Finsupp.degree hb
          rw [map_add, Finsupp.degree_single, ha] at he
          omega
        have hδ : Finsupp.degree (β - Finsupp.single i 1) < n :=
          lt_of_le_of_lt (Finsupp.degree_mono tsub_le_self) hβ
        have hmono : monomial α r = X i * monomial β r := by
          rw [X, monomial_mul, add_comm, hb, one_mul]
        rw [hmono, parts, polynomialMoment_X_mul, pderiv_monomial, ih _ hδ _ rfl]
  apply LinearMap.ext
  intro p
  induction p using MvPolynomial.induction_on' with
  | add p q hp hq => rw [map_add, map_add, hp, hq]
  | monomial α r => exact hm _ α rfl r

private theorem partial_matrixPolynomialPullback_X {I : Type*} (T : Matrix I J R) (i : I) (j : J) :
    pderiv j (matrixPolynomialPullback T (X i)) = C (T i j) := by
  rw [matrixPolynomialPullback_X, map_sum, Finset.sum_eq_single j]
  · rw [pderiv_C_mul, pderiv_X_self, mul_one]
  · intro k _ hk
    rw [pderiv_C_mul, pderiv_X_of_ne hk, mul_zero]
  · simp

/-- The chain rule for a linear substitution in the polynomial ring. -/
theorem partial_matrixPolynomialPullback {I : Type*} [Fintype I] [DecidableEq I]
    (T : Matrix I J R) (j : J) (p : MvPolynomial I R) :
    pderiv j (matrixPolynomialPullback T p) =
      ∑ i, C (T i j) * matrixPolynomialPullback T (pderiv i p) := by
  let L := matrixPolynomialPullback T
  change pderiv j (L p) = ∑ i, C (T i j) * L (pderiv i p)
  induction p using MvPolynomial.induction_on with
  | C r => simp [L]
  | add p q hp hq =>
    simp only [map_add, mul_add, Finset.sum_add_distrib]
    rw [hp, hq]
  | mul_X p k hp =>
    rw [map_mul, pderiv_mul, hp,
      show pderiv j (L (X k)) = C (T k j) from partial_matrixPolynomialPullback_X T k j]
    simp_rw [pderiv_mul, map_add, map_mul, mul_add, ← mul_assoc, Finset.sum_add_distrib,
      ← Finset.sum_mul]
    congr 1
    rw [Finset.sum_eq_single k]
    · rw [pderiv_X_self, map_one, mul_one, mul_comm]
    · intro i _ hi
      rw [pderiv_X_of_ne (Ne.symm hi), map_zero, mul_zero]
    · simp

private theorem polynomialMoment_mul_linear {I : Type*} (d : J → Rˣ) (T : Matrix I J R) (i : I)
    (p : MvPolynomial J R) :
    polynomialMoment d (matrixPolynomialPullback T (X i) * p) =
      ∑ j, (T i j * (↑((d j)⁻¹) : R)) * polynomialMoment d (pderiv j p) := by
  rw [matrixPolynomialPullback_X, Finset.sum_mul, map_sum]
  apply Finset.sum_congr rfl
  intro j _
  rw [mul_assoc, ← algebraMap_eq, ← Algebra.smul_def, map_smul, smul_eq_mul,
    polynomialMoment_X_mul, mul_assoc]

/-- A linear substitution with the transported inverse quadratic form
carries the polynomial moment to the original coordinates. -/
theorem polynomialMoment_matrixPullback_between {I : Type*} [Fintype I] [DecidableEq I]
    (d : I → Rˣ) (e : J → Rˣ) (T : Matrix I J R)
    (orthogonal : T * Matrix.diagonal (fun j ↦ (↑((e j)⁻¹) : R)) * Tᵀ =
      Matrix.diagonal (fun i ↦ (↑((d i)⁻¹) : R))) (p : MvPolynomial I R) :
    polynomialMoment e (matrixPolynomialPullback T p) = polynomialMoment d p := by
  let L := matrixPolynomialPullback T
  let moment := (polynomialMoment e).comp L.toLinearMap
  have hs (i k : I) : ∑ j, (T i j * (↑((e j)⁻¹) : R)) * T k j =
      Matrix.diagonal (fun j ↦ (↑((d j)⁻¹) : R)) i k := by
    have he := congrFun (congrFun orthogonal i) k
    simpa only [Matrix.mul_apply (M := T * Matrix.diagonal _), Matrix.mul_diagonal,
      Matrix.transpose_apply] using he
  have he : moment = polynomialMoment d := by
    apply polynomialMoment_unique d moment (by simp [moment])
    intro i g
    change polynomialMoment e (L (X i * g)) =
      (↑((d i)⁻¹) : R) * polynomialMoment e (L (pderiv i g))
    rw [map_mul, polynomialMoment_mul_linear]
    simp_rw [show ∀ j, pderiv j (L g) = ∑ k, C (T k j) * L (pderiv k g) from
      fun j ↦ partial_matrixPolynomialPullback T j g, map_sum, Finset.mul_sum]
    simp_rw [← algebraMap_eq, ← Algebra.smul_def, map_smul, smul_eq_mul]
    simp_rw [← mul_assoc]
    rw [Finset.sum_comm]
    simp_rw [← Finset.sum_mul, hs, Matrix.diagonal_apply, ite_mul, zero_mul,
      Finset.sum_ite_eq, Finset.mem_univ, if_true]
  exact LinearMap.congr_fun he p

/-- A matrix that preserves the inverse quadratic form preserves the
polynomial moment. The identity is valid over the coefficient ring. -/
theorem polynomialMoment_matrixPullback (d : J → Rˣ) (T : Matrix J J R)
    (orthogonal : T * Matrix.diagonal (fun j ↦ (↑((d j)⁻¹) : R)) * Tᵀ =
      Matrix.diagonal (fun j ↦ (↑((d j)⁻¹) : R))) (p : MvPolynomial J R) :
    polynomialMoment d (matrixPolynomialPullback T p) = polynomialMoment d p :=
  polynomialMoment_matrixPullback_between d d T orthogonal p

/-- The quadratic polynomial is invariant under a matrix preserving its
diagonal over the coefficient ring. -/
theorem matrixPullback_norm_of_conjugate (d : J → Rˣ) (T : Matrix J J R)
    (orthogonal : Tᵀ * Matrix.diagonal (fun j ↦ (d j : R)) * T =
      Matrix.diagonal (fun j ↦ (d j : R))) :
    matrixPolynomialPullback T (diagonalNormPolynomial d) = diagonalNormPolynomial d := by
  let A : Matrix J J (MvPolynomial J R) := T.map C
  let D : Matrix J J (MvPolynomial J R) := Matrix.diagonal fun j ↦ C (d j : R)
  let v : J → MvPolynomial J R := X
  have hD : (Matrix.diagonal (fun j ↦ (d j : R))).map (C : R →+* MvPolynomial J R) = D := by
    ext i j
    by_cases hij : i = j <;> simp [D, hij]
  have horth : Aᵀ * D * A = D := by
    have h := congrArg (fun B : Matrix J J R ↦ B.map (C : R →+* MvPolynomial J R)) orthogonal
    rw [Matrix.map_mul, Matrix.map_mul, hD] at h
    exact h
  have hquad (u : J → MvPolynomial J R) : u ⬝ᵥ D *ᵥ u = ∑ i, C (d i : R) * u i ^ 2 := by
    unfold dotProduct
    apply Finset.sum_congr rfl
    intro i _
    rw [Matrix.mulVec_diagonal, pow_two]
    ring
  have hA (i : J) : (A *ᵥ v) i = matrixPolynomialPullback T (X i) := by
    rw [matrixPolynomialPullback_X]
    rfl
  calc
    _ = (A *ᵥ v) ⬝ᵥ D *ᵥ (A *ᵥ v) := by
      rw [hquad]
      simp only [diagonalNormPolynomial, map_sum, map_mul, map_pow, hA,
        ← algebraMap_eq, AlgHom.commutes]
    _ = v ⬝ᵥ (Aᵀ * D * A) *ᵥ v := by
      rw [Matrix.dotProduct_mulVec, Matrix.dotProduct_mulVec, Matrix.dotProduct_mulVec,
        ← Matrix.vecMul_vecMul, ← Matrix.vecMul_vecMul, Matrix.vecMul_transpose]
    _ = _ := by rw [horth, hquad]; rfl

omit [DecidableEq J] in
private theorem matrixPullback_homogeneous {I : Type*} (T : Matrix I J R)
    {p : MvPolynomial I R} {m : ℕ}
    (hp : p.IsHomogeneous m) : (matrixPolynomialPullback T p).IsHomogeneous m := by
  have hh (i : I) : (∑ j, C (T i j) * X j : MvPolynomial J R).IsHomogeneous 1 := by
    apply Submodule.sum_mem (homogeneousSubmodule J R 1)
    intro j _
    exact (isHomogeneous_X _ _).C_mul _
  simpa [matrixPolynomialPullback] using hp.aeval _ hh

omit [DecidableEq J] in
/-- A linear substitution commutes with taking a homogeneous part. -/
theorem homogeneousComponent_matrixPullback {I : Type*} (T : Matrix I J R) (n : ℕ)
    (p : MvPolynomial I R) :
    homogeneousComponent n (matrixPolynomialPullback T p) =
      matrixPolynomialPullback T (homogeneousComponent n p) := by
  induction p using MvPolynomial.induction_on' with
  | add p q hp hq => rw [map_add, map_add, hp, hq, map_add, map_add]
  | monomial α r =>
    have hp := isHomogeneous_monomial r (rfl : Finsupp.degree α = Finsupp.degree α)
    rw [homogeneousComponent_of_mem (matrixPullback_homogeneous T hp), homogeneousComponent_of_mem hp]
    split_ifs <;> simp

end Polynomial

section Fractions

variable {K J : Type*} [Field K] [CharZero K] [Fintype J] [DecidableEq J] [Nonempty J]
variable (d : J → Kˣ)

omit [Nonempty J] in
private theorem polynomialMoment_norm_mul_monomial (α : J →₀ ℕ) (r : K) :
    polynomialMoment d (diagonalNormPolynomial d * monomial α r) =
      ((Fintype.card J : K) + ((Finsupp.degree α : ℕ) : K)) * polynomialMoment d (monomial α r) := by
  rw [diagonalNormPolynomial, Finset.sum_mul, map_sum, polynomialMoment_monomial]
  simp only [C_apply, X_pow_eq_monomial, monomial_mul, zero_add, mul_one, polynomialMoment_monomial]
  simp_rw [add_comm (Finsupp.single _ 2) α]
  have hf (i : J) :
      ((d i : K) * r) * monomialMoment d (α + Finsupp.single i 2) =
      (↑(α i) + 1) * (r * monomialMoment d α) := by
    calc
      _ = r * ((d i : K) * monomialMoment d (α + Finsupp.single i 2)) := by ring
      _ = _ := by rw [monomialMoment_add_two]; ring
  simp_rw [hf, ← Finset.sum_mul, Finset.sum_add_distrib, Finset.sum_const, Finset.card_univ,
    nsmul_eq_mul, mul_one, add_comm (∑ i : J, (α i : K))]
  congr 1
  rw [← Nat.cast_sum]
  norm_cast
  congr 1
  exact (Finsupp.sum_fintype α (fun _ n ↦ n) (fun _ ↦ rfl)).symm

private noncomputable def radialFactor (K : Type*) [Field K] (N n : ℕ) : K :=
  ∏ j ∈ Finset.range n, (↑(N + 2 * j) : K)

omit [DecidableEq J] in
private theorem radialFactor_ne_zero (n : ℕ) :
    radialFactor K (Fintype.card J) n ≠ 0 := by
  apply Finset.prod_ne_zero_iff.mpr
  intro j _
  exact_mod_cast Nat.ne_of_gt (by have := Fintype.card_pos (α := J); omega : 0 < Fintype.card J + 2 * j)

omit [CharZero K] in
private theorem radialFactor_succ (N n : ℕ) :
    radialFactor K N (n + 1) = radialFactor K N n * (↑(N + 2 * n) : K) :=
  Finset.prod_range_succ _ _

private noncomputable def normalizedMoment (n : ℕ) : MvPolynomial J K →ₗ[K] K :=
  (radialFactor K (Fintype.card J) n)⁻¹ •
    (polynomialMoment d).comp (homogeneousComponent (2 * n))

omit [CharZero K] [DecidableEq J] [Nonempty J] in
private theorem normalizedMoment_apply (n : ℕ) (p : MvPolynomial J K) :
    normalizedMoment d n p =
      polynomialMoment d (homogeneousComponent (2 * n) p) / radialFactor K (Fintype.card J) n := by
  simp only [normalizedMoment, LinearMap.smul_apply, LinearMap.comp_apply, smul_eq_mul]
  rw [div_eq_mul_inv, mul_comm]

omit [CharZero K] [DecidableEq J] [Nonempty J] in
private theorem norm_isHomogeneous : (diagonalNormPolynomial d).IsHomogeneous 2 := by
  rw [diagonalNormPolynomial]
  apply Submodule.sum_mem (homogeneousSubmodule J K 2)
  intro i _
  exact (isHomogeneous_X K i).pow 2 |>.C_mul (d i : K)

private theorem normalizedMoment_norm_mul (n : ℕ) (p : MvPolynomial J K) :
    normalizedMoment d (n + 1) (diagonalNormPolynomial d * p) = normalizedMoment d n p := by
  induction p using MvPolynomial.induction_on' with
  | add p q hp hq => rw [mul_add, map_add, hp, hq, map_add]
  | monomial α r =>
    have hp := isHomogeneous_monomial r (rfl : Finsupp.degree α = Finsupp.degree α)
    have hqp := (norm_isHomogeneous d).mul hp
    rw [normalizedMoment_apply, normalizedMoment_apply,
      homogeneousComponent_of_mem hqp, homogeneousComponent_of_mem hp]
    by_cases hn : 2 * n = Finsupp.degree α
    · have hsucc : 2 * (n + 1) = 2 + Finsupp.degree α := by omega
      rw [if_pos hsucc, if_pos hn, polynomialMoment_norm_mul_monomial d,
        radialFactor_succ]
      have hc : (↑(Fintype.card J + 2 * n) : K) ≠ 0 := by
        exact_mod_cast Nat.ne_of_gt (by have := Fintype.card_pos (α := J); omega :
          0 < Fintype.card J + 2 * n)
      rw [← hn, ← Nat.cast_add]
      field_simp
    · have hsucc : 2 * (n + 1) ≠ 2 + Finsupp.degree α := by omega
      simp [hn, hsucc]

private theorem normalizedMoment_norm_pow_mul (n k : ℕ) (p : MvPolynomial J K) :
    normalizedMoment d (n + k) (diagonalNormPolynomial d ^ k * p) = normalizedMoment d n p := by
  induction k with
  | zero => simp
  | succ k ih => rw [Nat.add_succ, pow_succ', mul_assoc, normalizedMoment_norm_mul, ih]

local notation "P" => MvPolynomial J K
local notation "q" => diagonalNormPolynomial d
local notation "A" => Localization.Away (diagonalNormPolynomial d)
local notation "ι" => algebraMap P A

private theorem same_fraction_moment (f : A) (p r : P) (n m : ℕ)
    (hp : f * ι (q ^ n) = ι p) (hr : f * ι (q ^ m) = ι r) :
    normalizedMoment d n p = normalizedMoment d m r := by
  have hi : Function.Injective ι :=
    IsLocalization.injective A (Submonoid.powers_le.mpr
      (mem_nonZeroDivisors_iff_ne_zero.mpr (diagonalNormPolynomial_ne_zero d)))
  have hcross : q ^ m * p = q ^ n * r := by
    apply hi
    rw [map_mul, map_mul, ← hp, ← hr]
    simp only [mul_left_comm, mul_comm]
  rw [← normalizedMoment_norm_pow_mul d n m p, ← normalizedMoment_norm_pow_mul d m n r,
    hcross, Nat.add_comm n m]

private noncomputable def fractionMoment (f : A) : K :=
  normalizedMoment d (IsLocalization.Away.sec q f).2 (IsLocalization.Away.sec q f).1

private theorem fractionMoment_eq (f : A) (p : P) (n : ℕ)
    (h : f * ι (q ^ n) = ι p) :
    fractionMoment d f = normalizedMoment d n p :=
  same_fraction_moment d f _ p _ n (IsLocalization.Away.sec_spec q f) h

/-- The normalized polynomial moment extended to the quadratic complement. -/
noncomputable def algebraicTrace : A →ₗ[K] K where
  toFun := fractionMoment d
  map_add' f g := by
    let n := (IsLocalization.Away.sec q f).2
    let m := (IsLocalization.Away.sec q g).2
    let p := (IsLocalization.Away.sec q f).1
    let r := (IsLocalization.Away.sec q g).1
    have hp : f * ι (q ^ n) = ι p := IsLocalization.Away.sec_spec _ _
    have hr : g * ι (q ^ m) = ι r := IsLocalization.Away.sec_spec _ _
    have hadd : (f + g) * ι (q ^ (n + m)) = ι (q ^ m * p + q ^ n * r) := by
      rw [add_mul, map_add, map_mul, map_mul, ← hp, ← hr, pow_add, map_mul]
      simp only [mul_left_comm, mul_comm]
    rw [fractionMoment_eq d _ _ _ hadd, map_add, normalizedMoment_norm_pow_mul,
      Nat.add_comm n m, normalizedMoment_norm_pow_mul]
    rfl
  map_smul' c f := by
    let n := (IsLocalization.Away.sec q f).2
    let p := (IsLocalization.Away.sec q f).1
    have hp : f * ι (q ^ n) = ι p := IsLocalization.Away.sec_spec _ _
    have hsmul : (c • f) * ι (q ^ n) = ι (c • p) := by
      rw [Algebra.smul_mul_assoc, hp, Algebra.smul_def, Algebra.smul_def, map_mul,
        ← IsScalarTower.algebraMap_apply]
    rw [fractionMoment_eq d _ _ _ hsmul, map_smul]
    rfl

theorem algebraicTrace_fraction (f : A) (p : P) (n : ℕ)
    (h : f * ι (q ^ n) = ι p) :
    algebraicTrace d f = polynomialMoment d (homogeneousComponent (2 * n) p) /
      radialFactor K (Fintype.card J) n :=
  (fractionMoment_eq d f p n h).trans (normalizedMoment_apply d n p)

@[simp] theorem algebraicTrace_one : algebraicTrace d 1 = 1 := by
  rw [algebraicTrace_fraction d 1 1 0 (by simp)]
  have h : homogeneousComponent 0 (1 : P) = 1 := by simp [homogeneousComponent_zero]
  simp [h, radialFactor]

omit [CharZero K] [DecidableEq J] [Nonempty J] in
theorem scaled_matrixPullback_monomial (T : Matrix J J K) (c : K)
    (α : J →₀ ℕ) (r : K) :
    matrixPolynomialPullback (c • T) (monomial α r) =
      (c ^ (Finsupp.degree α : ℕ)) • matrixPolynomialPullback T (monomial α r) := by
  simp only [matrixPolynomialPullback, aeval_monomial, Matrix.smul_apply, smul_eq_mul,
    map_mul, mul_assoc, ← Finset.mul_sum, ← algebraMap_eq, Algebra.smul_def,
    Finsupp.prod, Finsupp.degree_apply, mul_pow, Finset.prod_mul_distrib,
    Finset.prod_pow_eq_pow_sum, map_pow]
  ring

omit [CharZero K] [Nonempty J] in
private theorem normalizedMoment_scaled_pullback (T : Matrix J J K)
    (orthogonal : T * Matrix.diagonal (fun j ↦ (↑((d j)⁻¹) : K)) * Tᵀ =
      Matrix.diagonal (fun j ↦ (↑((d j)⁻¹) : K))) (c : K) (n : ℕ) (p : P) :
    normalizedMoment d n (matrixPolynomialPullback (c • T) p) =
      c ^ (2 * n) * normalizedMoment d n p := by
  induction p using MvPolynomial.induction_on' with
  | add p r hp hr => rw [map_add, map_add, hp, hr, map_add, mul_add]
  | monomial α r =>
    have hp := isHomogeneous_monomial r (rfl : Finsupp.degree α = Finsupp.degree α)
    have hL := matrixPullback_homogeneous (c • T) hp
    rw [normalizedMoment_apply, normalizedMoment_apply,
      homogeneousComponent_of_mem hL, homogeneousComponent_of_mem hp]
    by_cases hdeg : 2 * n = Finsupp.degree α
    · rw [if_pos hdeg, if_pos hdeg, scaled_matrixPullback_monomial, map_smul, smul_eq_mul,
        polynomialMoment_matrixPullback d T orthogonal, ← hdeg, mul_div_assoc]
    · simp [hdeg]

omit [Nonempty J] [CharZero K] in
/-- Congruence of two unit diagonal matrices implies the corresponding
transport law for their inverses. -/
theorem inverseDiagonal_conjugate_between {I R : Type*} [CommRing R]
    [Fintype I] [DecidableEq I] (d : I → Rˣ) (e : J → Rˣ)
    (T : Matrix I J R) (S : Matrix J I R) (hTS : T * S = 1)
    (hT : Tᵀ * Matrix.diagonal (fun i ↦ (d i : R)) * T = Matrix.diagonal (fun j ↦ (e j : R))) :
    T * Matrix.diagonal (fun j ↦ (↑((e j)⁻¹) : R)) * Tᵀ =
      Matrix.diagonal (fun i ↦ (↑((d i)⁻¹) : R)) := by
  let D := Matrix.diagonal fun i ↦ (d i : R)
  let E := Matrix.diagonal fun j ↦ (e j : R)
  let Dinv := Matrix.diagonal fun i ↦ (↑((d i)⁻¹) : R)
  let Einv := Matrix.diagonal fun j ↦ (↑((e j)⁻¹) : R)
  have he : Einv * E = 1 := by
    rw [Matrix.diagonal_mul_diagonal]
    simp [Units.inv_mul, Matrix.diagonal_one]
  have hd : D * Dinv = 1 := by
    rw [Matrix.diagonal_mul_diagonal]
    simp [Units.mul_inv, Matrix.diagonal_one]
  have hleft : Einv * Tᵀ * D = S := by
    calc
      _ = (Einv * Tᵀ * D) * (T * S) := by rw [hTS, Matrix.mul_one]
      _ = Einv * (Tᵀ * D * T) * S := by simp only [Matrix.mul_assoc]
      _ = _ := by rw [hT, he, Matrix.one_mul]
  change T * Einv * Tᵀ = Dinv
  calc
    _ = (T * Einv * Tᵀ) * (D * Dinv) := by rw [hd, mul_one]
    _ = T * (Einv * Tᵀ * D) * Dinv := by simp only [Matrix.mul_assoc]
    _ = _ := by rw [hleft, hTS, one_mul]

omit [Nonempty J] in
/-- Orthogonality of a matrix for a unit diagonal implies orthogonality
of the transposed matrix for the inverse diagonal. -/
theorem inverseDiagonal_orthogonal {R : Type*} [CommRing R] (e : J → Rˣ) (T : Matrix J J R)
    (hdet : IsUnit T.det)
    (hT : Tᵀ * Matrix.diagonal (fun j ↦ (e j : R)) * T = Matrix.diagonal (fun j ↦ (e j : R))) :
    T * Matrix.diagonal (fun j ↦ (↑((e j)⁻¹) : R)) * Tᵀ =
      Matrix.diagonal (fun j ↦ (↑((e j)⁻¹) : R)) :=
  inverseDiagonal_conjugate_between e e T T⁻¹ (Matrix.mul_nonsing_inv T hdet) hT

omit [Nonempty J] in
theorem diagonalConjugate_of_preserves_between {I : Type*} [Fintype I] [DecidableEq I]
    (d : I → Kˣ) (e : J → Kˣ) (T : Matrix I J K)
    (hT : ∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm e v) :
    Tᵀ * Matrix.diagonal (fun i ↦ (d i : K)) * T = Matrix.diagonal (fun j ↦ (e j : K)) := by
  have hp (v w : J → K) : diagonalPair d (T *ᵥ v) (T *ᵥ w) = diagonalPair e v w := by
    have he := hT (v + w)
    rw [Matrix.mulVec_add, diagonalNorm_add, diagonalNorm_add, hT, hT] at he
    linear_combination he / 2
  ext i j
  have he := hp (Pi.single i 1) (Pi.single j 1)
  rw [Matrix.mulVec_single_one, Matrix.mulVec_single_one, diagonalPair_single_right] at he
  change diagonalPair d (fun k ↦ T k i) (fun k ↦ T k j) = _ at he
  have hh : (Tᵀ * Matrix.diagonal (fun k ↦ (d k : K)) * T) i j =
      diagonalPair d (fun k ↦ T k i) (fun k ↦ T k j) := by
    rw [Matrix.mul_apply]
    apply Finset.sum_congr rfl
    intro k _
    rw [Matrix.mul_diagonal, Matrix.transpose_apply]
    ring
  rw [hh, he, Matrix.diagonal_apply, Pi.single_apply, mul_ite, mul_one, mul_zero]
  by_cases hij : i = j
  · subst j; simp
  · simp [hij, Ne.symm hij]

omit [Nonempty J] in
theorem diagonalConjugate_of_preserves (T : Matrix J J K)
    (hT : ∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v) :
    Tᵀ * Matrix.diagonal (fun j ↦ (d j : K)) * T = Matrix.diagonal (fun j ↦ (d j : K)) :=
  diagonalConjugate_of_preserves_between d d T hT

/-- The algebraic trace is invariant under the scalar and orthogonal point
actions on the norm complement. -/
theorem algebraicTrace_geometricPullback [Infinite K] (T : Matrix J J K) (c : Kˣ) (hdet : IsUnit T.det)
    (hT : ∀ v, diagonalNorm d (T *ᵥ v) = diagonalNorm d v) (f : A) :
    algebraicTrace d (geometricPullback d T c hT f) = algebraicTrace d f := by
  let ψ := geometricPullback d T c hT
  let p := (IsLocalization.Away.sec q f).1
  let n := (IsLocalization.Away.sec q f).2
  have hp : f * ι (q ^ n) = ι p := IsLocalization.Away.sec_spec _ _
  have hs' := geometricPullback_clear_denominator d T c hT f p n hp
  have ht := fractionMoment_eq d _ _ _ hs'
  change algebraicTrace d (((c : K) ^ (2 * n)) • ψ f) = _ at ht
  rw [map_smul, smul_eq_mul, normalizedMoment_scaled_pullback d T
    (inverseDiagonal_orthogonal d T hdet (diagonalConjugate_of_preserves d T hT)),
    ← fractionMoment_eq d f p n hp] at ht
  exact mul_left_cancel₀ (pow_ne_zero _ c.ne_zero) ht

omit [Fintype J] [DecidableEq J] [Nonempty J] in
theorem map_homogeneousComponent {R S : Type*} [CommRing R] [CommRing S]
    (f : R →+* S) (n : ℕ) (p : MvPolynomial J R) :
    homogeneousComponent n (MvPolynomial.map f p) =
      MvPolynomial.map f (homogeneousComponent n p) := by
  ext α
  rw [coeff_homogeneousComponent, coeff_map, coeff_map, coeff_homogeneousComponent]
  split_ifs <;> simp

/-- The trace value of an extended fraction is the extension of the
original trace value. -/
theorem algebraicTrace_map_fraction {L : Type*} [Field L] [CharZero L] [Algebra K L]
    (f : A) (g : Localization.Away
      (diagonalNormPolynomial (Units.map (algebraMap K L).toMonoidHom ∘ d)))
    (p : P) (n : ℕ) (hf : f * ι (q ^ n) = ι p)
    (hg : g * algebraMap (MvPolynomial J L) _
        (diagonalNormPolynomial (Units.map (algebraMap K L).toMonoidHom ∘ d) ^ n) =
      algebraMap (MvPolynomial J L) _ (MvPolynomial.map (algebraMap K L) p)) :
    algebraicTrace (Units.map (algebraMap K L).toMonoidHom ∘ d) g =
      algebraMap K L (algebraicTrace d f) := by
  rw [algebraicTrace_fraction _ g _ n hg, algebraicTrace_fraction d f p n hf,
    map_homogeneousComponent, polynomialMoment_map, map_div₀]
  congr 1
  simp only [radialFactor, map_prod, map_natCast]

/-- Fractions related by an isometric linear coordinate change have
the same trace. -/
theorem algebraicTrace_fraction_between {I : Type*} [Fintype I] [DecidableEq I] [Nonempty I]
    (e : I → Kˣ) (T : Matrix I J K) (hcard : Fintype.card I = Fintype.card J)
    (hT : T * Matrix.diagonal (fun j ↦ (↑((d j)⁻¹) : K)) * Tᵀ =
      Matrix.diagonal (fun i ↦ (↑((e i)⁻¹) : K)))
    (f : Localization.Away (diagonalNormPolynomial e)) (g : A)
    (p : MvPolynomial I K) (n : ℕ)
    (hf : f * algebraMap _ (Localization.Away (diagonalNormPolynomial e))
        (diagonalNormPolynomial e ^ n) = algebraMap _ _ p)
    (hg : g * ι (q ^ n) = ι (matrixPolynomialPullback T p)) :
    algebraicTrace d g = algebraicTrace e f := by
  rw [algebraicTrace_fraction d g _ n hg, algebraicTrace_fraction e f p n hf,
    homogeneousComponent_matrixPullback, polynomialMoment_matrixPullback_between e d T hT, hcard]

/-- The trace on a given presentation of the diagonal norm polynomial. -/
noncomputable def algebraicTraceAt (Q : P) (hQ : Q = q) : Localization.Away Q →ₗ[K] K := by
  subst Q
  exact algebraicTrace d

@[simp] theorem algebraicTraceAt_one (Q : P) (hQ : Q = q) :
    algebraicTraceAt d Q hQ 1 = 1 := by
  subst Q
  exact algebraicTrace_one d

/-- A pullback scaling every coordinate by the same unit preserves the
algebraic trace. -/
theorem algebraicTraceAt_scalar (Q : P) (hQ : Q = q) (c : Kˣ)
    (ψ : Localization.Away Q →ₐ[K] Localization.Away Q)
    (hψ : ∀ i, ψ (algebraMap P (Localization.Away Q) (X i)) =
      algebraMap K (Localization.Away Q) (c : K) * algebraMap P (Localization.Away Q) (X i))
    (f : Localization.Away Q) :
    algebraicTraceAt d Q hQ (ψ f) = algebraicTraceAt d Q hQ f := by
  subst Q
  letI : Infinite K := Infinite.of_injective (fun n : ℕ ↦ (n : K)) Nat.cast_injective
  have hT : ∀ x, diagonalNorm d ((1 : Matrix J J K) *ᵥ x) = diagonalNorm d x := by simp
  have h : ψ = geometricPullback d 1 c hT := by
    apply IsLocalization.algHom_ext (Submonoid.powers q)
    apply MvPolynomial.algHom_ext
    intro i
    change ψ (ι (X i)) = geometricPullback d 1 c hT (ι (X i))
    rw [hψ, geometricPullback_polynomial_eq_matrixPullback, matrixPolynomialPullback_X,
      Finset.sum_eq_single i]
    · rw [Matrix.smul_apply, Matrix.one_apply_eq, smul_eq_mul, mul_one,
        map_mul, ← algebraMap_eq, ← IsScalarTower.algebraMap_apply]
    · intro j _ hji
      rw [Matrix.smul_apply, Matrix.one_apply_ne hji.symm, smul_zero, map_zero, zero_mul]
    · simp
  rw [h]
  exact algebraicTrace_geometricPullback d 1 c (by rw [Matrix.det_one]; exact isUnit_one) hT f

end Fractions

end MathieuQuadratic
