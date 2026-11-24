import MathieuQuadratic.AlgebraicTrace
import Mathlib.Algebra.DualNumber

/-!
# First-order invariance of the algebraic trace

The moment formula is compatible with coefficient extension to dual
numbers. Applying orthogonal invariance over that ring and taking the
first-order coefficient gives the rotation identity on the localization.
-/

open scoped BigOperators DualNumber Matrix

namespace MathieuQuadratic

open MvPolynomial TrivSqZeroExt

variable {K J : Type*} [Field K] [Fintype J] [DecidableEq J]

/-- Extend both coefficients of a dual number. -/
def dualCoefficientMap {A : Type*} [CommRing A] [Algebra K A] : K[ε] →+* A[ε] where
  toFun x := ⟨algebraMap K A x.fst, algebraMap K A x.snd⟩
  map_zero' := by ext <;> simp
  map_one' := by ext <;> simp
  map_add' _ _ := by ext <;> simp
  map_mul' _ _ := by
    apply TrivSqZeroExt.ext
    · exact map_mul _ _ _
    · simp only [DualNumber.snd_mul, map_add, map_mul]
      rfl

@[simp] theorem dualCoefficientMap_inl {A : Type*} [CommRing A] [Algebra K A] (k : K) :
    dualCoefficientMap (A := A) (inl k) = inl (algebraMap K A k) := by
  ext <;> simp [dualCoefficientMap]

@[simp] theorem dualCoefficientMap_eps {A : Type*} [CommRing A] [Algebra K A] :
    dualCoefficientMap (A := A) (ε : K[ε]) = ε := by
  ext <;> simp [dualCoefficientMap]

private noncomputable abbrev embedDualPoly : MvPolynomial J K →+* MvPolynomial J K[ε] := MvPolynomial.map (TrivSqZeroExt.inlHom K K)

private abbrev embedDualWeight (d : J → Kˣ) :=
  Units.map (TrivSqZeroExt.inlHom K K).toMonoidHom ∘ d

private noncomputable def dualPolynomialFst (p : MvPolynomial J K[ε]) : MvPolynomial J K :=
  AddMonoidAlgebra.ofCoeff ((AddMonoidAlgebra.coeff p).mapRange (fun z : K[ε] ↦ z.fst) rfl)

private noncomputable def dualPolynomialSnd (p : MvPolynomial J K[ε]) : MvPolynomial J K :=
  AddMonoidAlgebra.ofCoeff ((AddMonoidAlgebra.coeff p).mapRange (fun z : K[ε] ↦ z.snd) rfl)

omit [Fintype J] [DecidableEq J] in
private theorem dualPolynomial_split (p : MvPolynomial J K[ε]) :
    p = embedDualPoly (dualPolynomialFst p) + C (ε : K[ε]) * embedDualPoly (dualPolynomialSnd p) := by
  apply MvPolynomial.ext
  intro α
  rw [coeff_add, coeff_C_mul, coeff_map, coeff_map]
  change coeff α p = inl (coeff α p).fst + ε * inl (coeff α p).snd
  apply TrivSqZeroExt.ext <;> simp

omit [DecidableEq J] in
private theorem snd_polynomialMoment (d : J → Kˣ) (p : MvPolynomial J K[ε]) (n : ℕ) :
    (polynomialMoment (embedDualWeight d) (homogeneousComponent n p)).snd =
      polynomialMoment d (homogeneousComponent n (dualPolynomialSnd p)) := by
  nth_rw 1 [dualPolynomial_split p]
  rw [map_add, homogeneousComponent_C_mul, map_homogeneousComponent, map_homogeneousComponent,
    map_add, ← algebraMap_eq, ← Algebra.smul_def, map_smul, smul_eq_mul,
    polynomialMoment_map, polynomialMoment_map]
  simp

private def evalDualPoly {A : Type*} [CommRing A] [Algebra K A] (x : J → A) :
    MvPolynomial J K[ε] →+* A[ε] :=
  eval₂Hom dualCoefficientMap (fun i ↦ inl (x i))

omit [Fintype J] [DecidableEq J] in
private theorem evalDualPoly_embed {A : Type*} [CommRing A] [Algebra K A] (x : J → A)
    (p : MvPolynomial J K) :
    evalDualPoly x (embedDualPoly p) = inl (eval₂ (algebraMap K A) x p) := by
  change ((evalDualPoly x).comp embedDualPoly) p =
    ((inlHom A A).comp (eval₂Hom (algebraMap K A) x)) p
  congr 1
  apply MvPolynomial.ringHom_ext
  · intro r
    simp [evalDualPoly, eval₂Hom_C]
  · intro i
    simp [evalDualPoly, eval₂Hom_X']

omit [Fintype J] [DecidableEq J] in
private theorem snd_evalDualPoly {A : Type*} [CommRing A] [Algebra K A] (x : J → A)
    (p : MvPolynomial J K[ε]) :
    (evalDualPoly x p).snd = eval₂ (algebraMap K A) x (dualPolynomialSnd p) := by
  nth_rw 1 [dualPolynomial_split p]
  rw [map_add, map_mul, evalDualPoly_embed, evalDualPoly_embed]
  simp [evalDualPoly, eval₂Hom_C]

variable [CharZero K] [Nonempty J] (d : J → Kˣ)

local notation "P" => MvPolynomial J K
local notation "q" => diagonalNormPolynomial d
local notation "A" => Localization.Away (diagonalNormPolynomial d)
local notation "ι" => algebraMap P A
local notation "v" => (fun i : J ↦ ι (X i))

/-- The first-order matrix acts on the coordinate generators by this
matrix-vector product. -/
def IsDualMatrixPoint (T : Matrix J J K[ε]) (ψ : A →ₐ[K] A[ε]) : Prop :=
  ∀ i, ψ (v i) = (T.map dualCoefficientMap *ᵥ (fun j ↦ inl (v j))) i

omit [DecidableEq J] [CharZero K] [Nonempty J] in
private theorem dualEval_matrix_X (T : Matrix J J K[ε]) (i : J) :
    evalDualPoly v (matrixPolynomialPullback T (embedDualPoly (X i))) =
      (T.map dualCoefficientMap *ᵥ (fun j ↦ inl (v j))) i := by
  rw [map_X, matrixPolynomialPullback_X]
  simp [evalDualPoly, Matrix.mulVec, dotProduct]

omit [DecidableEq J] [CharZero K] [Nonempty J] in
private theorem eval₂_localization_variables (p : P) :
    eval₂ (algebraMap K A) v p = ι p := by
  change (eval₂Hom (algebraMap K A) v) p = _
  congr 1
  apply MvPolynomial.ringHom_ext
  · intro r
    rw [eval₂Hom_C, ← algebraMap_eq, ← IsScalarTower.algebraMap_apply]
  · intro i
    exact eval₂Hom_X' _ _ _

omit [DecidableEq J] [CharZero K] [Nonempty J] in
private theorem dualMatrixPoint_polynomial (T : Matrix J J K[ε]) (ψ : A →ₐ[K] A[ε])
    (hψ : IsDualMatrixPoint d T ψ) (p : P) :
    ψ (ι p) = evalDualPoly v (matrixPolynomialPullback T (embedDualPoly p)) := by
  change (ψ.toRingHom.comp ι) p =
    (((evalDualPoly v).comp (matrixPolynomialPullback T).toRingHom).comp embedDualPoly) p
  congr 1
  apply MvPolynomial.ringHom_ext
  · intro r
    rw [RingHom.comp_apply, ← algebraMap_eq, ← IsScalarTower.algebraMap_apply,
      AlgHom.toRingHom_eq_coe, AlgHom.coe_toRingHom, AlgHom.commutes]
    change algebraMap K A[ε] r = evalDualPoly v (matrixPolynomialPullback T (embedDualPoly (C r)))
    rw [map_C, ← algebraMap_eq, AlgHom.commutes, algebraMap_eq, evalDualPoly, eval₂Hom_C]
    exact (dualCoefficientMap_inl r).symm
  · intro i
    exact (hψ i).trans (dualEval_matrix_X d T i).symm

/-- A determinant-one first-order isometry has zero first-order
coefficient against the algebraic trace. -/
theorem algebraicTrace_firstOrder_zero (T : Matrix J J K[ε]) (ψ : A →ₐ[K] A[ε])
    (hdet : T.det = 1)
    (hT : Tᵀ * Matrix.diagonal (fun j ↦ inl (d j : K)) * T =
      Matrix.diagonal (fun j ↦ inl (d j : K)))
    (hψ : IsDualMatrixPoint d T ψ) (f : A) :
    algebraicTrace d (ψ f).snd = 0 := by
  let p := (IsLocalization.Away.sec q f).1
  let n := (IsLocalization.Away.sec q f).2
  let u := matrixPolynomialPullback T (embedDualPoly p)
  have hp : f * ι (q ^ n) = ι p := IsLocalization.Away.sec_spec _ _
  have hnormMap : embedDualPoly q = diagonalNormPolynomial (embedDualWeight d) := by
    simp [diagonalNormPolynomial, embedDualWeight]
  have hnorm : ψ (ι q) = inl (ι q) := by
    rw [dualMatrixPoint_polynomial d T ψ hψ, hnormMap,
      matrixPullback_norm_of_conjugate (embedDualWeight d) T hT,
      ← hnormMap, evalDualPoly_embed, eval₂_localization_variables]
  have hpow : ψ (ι (q ^ n)) = inl (ι (q ^ n)) := by
    rw [map_pow, map_pow, hnorm]
    change (inlHom A A (ι q)) ^ n = inlHom A A (ι q ^ n)
    rw [map_pow]
  have hf := congrArg (fun g : A ↦ (ψ g).snd) hp
  change (ψ (f * ι (q ^ n))).snd = (ψ (ι p)).snd at hf
  rw [map_mul, hpow, DualNumber.snd_mul, snd_inl, fst_inl, mul_zero, zero_add,
    dualMatrixPoint_polynomial d T ψ hψ p, snd_evalDualPoly, eval₂_localization_variables] at hf
  have he := algebraicTrace_fraction d (ψ f).snd (dualPolynomialSnd u) n hf
  have hor := inverseDiagonal_orthogonal (embedDualWeight d) T
    (by rw [hdet]; exact isUnit_one) hT
  have hmoment : polynomialMoment d (homogeneousComponent (2 * n) (dualPolynomialSnd u)) = 0 := by
    rw [← snd_polynomialMoment]
    have hH : homogeneousComponent (2 * n) u =
        matrixPolynomialPullback T (homogeneousComponent (2 * n) (embedDualPoly p)) := by
      exact homogeneousComponent_matrixPullback T _ _
    rw [hH, polynomialMoment_matrixPullback _ T hor, map_homogeneousComponent, polynomialMoment_map]
    rfl
  rw [he, hmoment, zero_div]

theorem algebraicTraceAt_firstOrder_zero (Q : P) (hQ : Q = q) (T : Matrix J J K[ε])
    (ψ : Localization.Away Q →ₐ[K] (Localization.Away Q)[ε]) (hdet : T.det = 1)
    (hT : Tᵀ * Matrix.diagonal (fun j ↦ inl (d j : K)) * T =
      Matrix.diagonal (fun j ↦ inl (d j : K)))
    (hψ : ∀ i, ψ (algebraMap P (Localization.Away Q) (X i)) =
      (T.map dualCoefficientMap *ᵥ (fun j ↦ inl (algebraMap P (Localization.Away Q) (X j)))) i)
    (f : Localization.Away Q) :
    algebraicTraceAt d Q hQ (ψ f).snd = 0 := by
  subst Q
  exact algebraicTrace_firstOrder_zero d T ψ hdet hT hψ f

end MathieuQuadratic
