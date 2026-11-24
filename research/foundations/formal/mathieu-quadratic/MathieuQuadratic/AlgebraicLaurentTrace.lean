import MathieuQuadratic.Classification
import MathieuQuadratic.GeometricProjection
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure
import MathieuQuadratic.LaurentTrace
import MathieuQuadratic.Newton
import MathieuQuadratic.OneVariableLaurent

/-!
# Comparison of the algebraic trace with Laurent constant coefficient
-/

open scoped BigOperators Matrix TensorProduct

namespace MathieuQuadratic

open MvPolynomial

variable {K : Type*} [Field K] [CharZero K]

omit [CharZero K] in
private theorem unarySplit_eq_diagonal (a : Kˣ) :
    splitNorm a (fun _ : Fin 1 ↦ 2) = diagonalNormPolynomial (fun _ : Fin 1 ↦ a) := by
  simp [splitNorm, diagonalNormPolynomial]

/-- The algebraic trace in dimension one is exactly the constant coefficient
in the Laurent coordinate. -/
theorem algebraicTrace_unary_eq_laurentCT (a : Kˣ) :
    (algebraicTraceAt (fun _ : Fin 1 ↦ a) _ (unarySplit_eq_diagonal a)).comp
        (splitComplementEquivLaurent a (fun _ : Fin 1 ↦ 2) (by simp)).symm.toLinearMap =
      laurentCT := by
  apply splitUnaryTrace_eq_laurentCT a
  · exact algebraicTraceAt_one _ _ _
  · intro f
    exact algebraicTraceAt_scalar _ _ _ scalarTwo _ (splitDiagonalPullback_X _ _ _) f

private theorem algebraicTraceAt_isMathieu_iff {J : Type*} [Fintype J] [DecidableEq J] [Nonempty J]
    (d : J → Kˣ) (q : MvPolynomial J K) (hq : q = diagonalNormPolynomial d) :
    IsMathieu (algebraicTraceAt d q hq) ↔ IsMathieu (algebraicTrace d) := by
  subst q
  rfl

/-- The rank-one Mathieu conclusion uses only Laurent necessity. -/
theorem algebraicTrace_unary_isMathieu (a : Kˣ) (dvk : NewtonLaurentNecessity ℂ (Fin 1)) :
    IsMathieu (algebraicTrace (fun _ : Fin 1 ↦ a)) :=
  (algebraicTraceAt_isMathieu_iff _ _ (unarySplit_eq_diagonal a)).mp
    (isMathieu_of_laurent_equiv_trace_eq _ _ (algebraicTrace_unary_eq_laurentCT a)
      (laurentSeparation_of_complex (laurentSeparation_of_newton dvk)))

section Binary

variable (a b δ : Kˣ)

private abbrev binaryWeights : Fin 2 → Kˣ := ![a, b]

omit [CharZero K] in
private theorem binary_eq_diagonal :
    binaryNorm a b = diagonalNormPolynomial (binaryWeights a b) := by
  simp [binaryNorm, diagonalNormPolynomial, Fin.sum_univ_two, binaryWeights]

/-- The split chart of a binary diagonal quadratic form. -/
noncomputable def binarySplitChart (hδ : (δ : K) ^ 2 = -(b : K) / (a : K)) :
    SplitComplement a (fun _ : Fin 2 ↦ 1) ≃ₐ[K]
      Localization.Away (diagonalNormPolynomial (binaryWeights a b)) :=
  IsLocalization.algEquivOfAlgEquiv _ _ (binaryLinearCoordinates δ) (by
    rw [Submonoid.map_powers, binaryLinearCoordinates_norm a b δ hδ, binary_eq_diagonal])

private theorem binarySplitChart_polynomial (hδ : (δ : K) ^ 2 = -(b : K) / (a : K))
    (p : MvPolynomial (Fin 2) K) :
    binarySplitChart a b δ hδ (algebraMap _ _ p) =
      algebraMap (MvPolynomial (Fin 2) K) _ (binaryLinearCoordinates δ p) := by
  rw [binarySplitChart, IsLocalization.algEquivOfAlgEquiv_apply, IsLocalization.map_eq]
  rfl

private def torusDiagonalConjugate (c : Kˣ) : Matrix (Fin 2) (Fin 2) K :=
  let s := ((c : K) + ↑(c⁻¹)) / 2
  let t := ((c : K) - ↑(c⁻¹)) / 2
  !![s, (δ : K) * t; t / (δ : K), s]

private theorem torusDiagonalConjugate_det (c : Kˣ) :
    (torusDiagonalConjugate δ c).det = 1 := by
  rw [torusDiagonalConjugate, Matrix.det_fin_two_of]
  have hc := Units.mul_inv c
  field_simp [δ.ne_zero]
  linear_combination 4 * hc

private theorem torusDiagonalConjugate_preserves (hδ : (δ : K) ^ 2 = -(b : K) / (a : K))
    (c : Kˣ) (x : Fin 2 → K) :
    diagonalNorm (binaryWeights a b) (torusDiagonalConjugate δ c *ᵥ x) =
      diagonalNorm (binaryWeights a b) x := by
  let s := ((c : K) + ↑(c⁻¹)) / 2
  let t := ((c : K) - ↑(c⁻¹)) / 2
  have hs : s ^ 2 - t ^ 2 = 1 := by
    dsimp [s, t]
    linear_combination Units.mul_inv c
  have hb : (b : K) = -(a : K) * (δ : K) ^ 2 := by
    rw [hδ]
    field_simp [a.ne_zero]
  simp only [diagonalNorm, diagonalPair, Fin.sum_univ_two, binaryWeights, Matrix.cons_val_zero,
    Matrix.cons_val_one, torusDiagonalConjugate, Matrix.mulVec, dotProduct,
    Matrix.of_apply]
  change (a : K) * (s * x 0 + (δ : K) * t * x 1) * (s * x 0 + (δ : K) * t * x 1) +
      (b : K) * (t / (δ : K) * x 0 + s * x 1) * (t / (δ : K) * x 0 + s * x 1) =
      (a : K) * x 0 * x 0 + (b : K) * x 1 * x 1
  have he : (δ : K) * (t / (δ : K) * x 0 + s * x 1) =
      t * x 0 + (δ : K) * s * x 1 := by
    field_simp [δ.ne_zero]
  rw [hb]
  calc
    _ = (a : K) * ((s * x 0 + (δ : K) * t * x 1) ^ 2 -
        ((δ : K) * (t / (δ : K) * x 0 + s * x 1)) ^ 2) := by ring
    _ = (a : K) * ((x 0) ^ 2 - ((δ : K) * x 1) ^ 2) := by
      rw [he]
      linear_combination ((a : K) * ((x 0) ^ 2 - (δ : K) ^ 2 * (x 1) ^ 2)) * hs
    _ = _ := by ring

private theorem torusDiagonalConjugate_X [Infinite K] (c : Kˣ) (i : Fin 2) :
    matrixPolynomialPullback (torusDiagonalConjugate δ c) (binaryLinearCoordinates δ (X i)) =
      C (if i = 0 then (c : K) else (↑(c⁻¹) : K)) * binaryLinearCoordinates δ (X i) := by
  apply MvPolynomial.funext
  intro x
  fin_cases i <;>
    simp [binaryLinearCoordinates_X, matrixPolynomialPullback_X, Fin.sum_univ_two,
      torusDiagonalConjugate] <;>
    field_simp [δ.ne_zero] <;> ring

private theorem scalar_binary_X [Infinite K] (c : Kˣ) (i : Fin 2) :
    matrixPolynomialPullback ((c : K) • (1 : Matrix (Fin 2) (Fin 2) K))
        (binaryLinearCoordinates δ (X i)) =
      C (c : K) * binaryLinearCoordinates δ (X i) := by
  apply MvPolynomial.funext
  intro x
  fin_cases i <;>
    simp [binaryLinearCoordinates_X, matrixPolynomialPullback_X, Fin.sum_univ_two,
      Matrix.one_apply, Matrix.smul_apply, smul_eq_mul] <;> ring

private theorem binarySplitChart_intertwines [Infinite K] (hδ : (δ : K) ^ 2 = -(b : K) / (a : K))
    (T : Matrix (Fin 2) (Fin 2) K) (c : Kˣ)
    (hT : ∀ x, diagonalNorm (binaryWeights a b) (T *ᵥ x) = diagonalNorm (binaryWeights a b) x)
    (χ : Fin 2 → Kˣ)
    (hχ : ∀ i, matrixPolynomialPullback ((c : K) • T) (binaryLinearCoordinates δ (X i)) =
      C (χ i : K) * binaryLinearCoordinates δ (X i))
    (f : SplitComplement a (fun _ : Fin 2 ↦ 1)) :
    geometricPullback (binaryWeights a b) T c hT (binarySplitChart a b δ hδ f) =
      binarySplitChart a b δ hδ (splitDiagonalPullback a (fun _ : Fin 2 ↦ 1) χ f) := by
  let e := binarySplitChart a b δ hδ
  let g := geometricPullback (binaryWeights a b) T c hT
  let h := splitDiagonalPullback a (fun _ : Fin 2 ↦ 1) χ
  have he : g.comp e.toAlgHom = e.toAlgHom.comp h := by
    apply IsLocalization.algHom_ext (Submonoid.powers (splitNorm a (fun _ : Fin 2 ↦ 1)))
    apply MvPolynomial.algHom_ext
    intro i
    change g (e (algebraMap _ _ (X i))) = e (h (algebraMap _ _ (X i)))
    dsimp only [e, g, h]
    rw [binarySplitChart_polynomial, geometricPullback_polynomial_eq_matrixPullback, hχ,
      splitDiagonalPullback_X, map_mul, map_mul, AlgEquiv.commutes, binarySplitChart_polynomial, ← algebraMap_eq, ← IsScalarTower.algebraMap_apply]
  exact AlgHom.congr_fun he f

/-- After the binary linear change of coordinates, the algebraic trace
is the Laurent constant coefficient. -/
theorem algebraicTrace_binary_eq_laurentCT (hδ : (δ : K) ^ 2 = -(b : K) / (a : K)) :
    ((algebraicTrace (binaryWeights a b)).comp (binarySplitChart a b δ hδ).toLinearMap).comp
        (splitComplementEquivLaurent a (fun _ : Fin 2 ↦ 1) (by simp)).symm.toLinearMap = laurentCT := by
  letI : Infinite K := Infinite.of_injective (fun n : ℕ ↦ (n : K)) Nat.cast_injective
  let e := binarySplitChart a b δ hδ
  let τ := algebraicTrace (binaryWeights a b)
  apply splitBinaryTrace_eq_laurentCT a
  · simp
  · intro f
    have hI (x : Fin 2 → K) : diagonalNorm (binaryWeights a b) ((1 : Matrix (Fin 2) (Fin 2) K) *ᵥ x) =
        diagonalNorm (binaryWeights a b) x := by simp
    change τ (e (splitDiagonalPullback a _ _ f)) = τ (e f)
    rw [← binarySplitChart_intertwines a b δ hδ 1 scalarTwo hI _ (scalar_binary_X δ scalarTwo)]
    exact algebraicTrace_geometricPullback _ _ _ (by simp) hI _
  · intro f
    let T := torusDiagonalConjugate δ scalarTwo
    have hT := torusDiagonalConjugate_preserves a b δ hδ scalarTwo
    have hχ (i : Fin 2) : matrixPolynomialPullback ((↑(1 : Kˣ) : K) • T)
        (binaryLinearCoordinates δ (X i)) =
          C (↑(inverseTorusTwo (K := K) i) : K) * binaryLinearCoordinates δ (X i) := by
      rw [Units.val_one, one_smul]
      fin_cases i
      · simpa [inverseTorusTwo, T] using torusDiagonalConjugate_X δ scalarTwo 0
      · simpa [inverseTorusTwo, T] using torusDiagonalConjugate_X δ scalarTwo 1
    change τ (e (splitDiagonalPullback a _ _ f)) = τ (e f)
    rw [← binarySplitChart_intertwines a b δ hδ T 1 hT _ hχ]
    exact algebraicTrace_geometricPullback _ _ _
      (by rw [show T.det = 1 from torusDiagonalConjugate_det δ scalarTwo]; exact isUnit_one) hT _

/-- The rank-two Mathieu conclusion follows by scalar extension to the
split chart and the rank-two Laurent necessity theorem. -/
theorem algebraicTrace_binary_isMathieu (dvk : NewtonLaurentNecessity ℂ (Fin 2)) :
    IsMathieu (algebraicTrace (binaryWeights a b)) := by
  let L := AlgebraicClosure K
  let aL := Units.map (algebraMap K L).toMonoidHom a
  let bL := Units.map (algebraMap K L).toMonoidHom b
  obtain ⟨r, hr⟩ := IsAlgClosed.exists_pow_nat_eq (-(bL : L) / (aL : L)) (n := 2) (by decide)
  have hr0 : r ≠ 0 := by
    intro he
    exact (div_ne_zero (neg_ne_zero.mpr bL.ne_zero) aL.ne_zero) (by simpa [he] using hr.symm)
  let δL := Units.mk0 r hr0
  let e := binarySplitChart aL bL δL hr
  let z := splitComplementEquivLaurent aL (fun _ : Fin 2 ↦ 1) (by simp)
  have hw : extendedDiagonal (L := L) (binaryWeights a b) = binaryWeights aL bL := by
    funext i
    fin_cases i <;> rfl
  let ext : L ⊗[K] Localization.Away (diagonalNormPolynomial (binaryWeights a b)) ≃ₐ[L]
      Localization.Away (diagonalNormPolynomial (binaryWeights aL bL)) :=
    hw ▸ diagonalComplementBaseChange (L := L) (binaryWeights a b)
  apply isMathieu_of_split_laurent (L := L) _ (ext.trans (e.symm.trans z)) _
    (laurentSeparation_of_complex (laurentSeparation_of_newton dvk))
  intro f
  have hbaseAux (w : Fin 2 → Lˣ) (h : extendedDiagonal (L := L) (binaryWeights a b) = w) :
      (algebraicTrace w).comp (h ▸ diagonalComplementBaseChange (L := L) (binaryWeights a b)).toLinearMap =
        extendTrace L (algebraicTrace (binaryWeights a b)) := by
    subst w
    exact algebraicTrace_baseChange (L := L) (binaryWeights a b)
  have hbase : (algebraicTrace (binaryWeights aL bL)).comp ext.toLinearMap =
      extendTrace L (algebraicTrace (binaryWeights a b)) := hbaseAux _ hw
  have hvalue := LinearMap.congr_fun (algebraicTrace_binary_eq_laurentCT aL bL δL hr)
    (z (e.symm (ext f)))
  rw [← hbase, LinearMap.comp_apply, AlgEquiv.toLinearMap_apply]
  change algebraicTrace (binaryWeights aL bL) (ext f) = laurentCT (z (e.symm (ext f)))
  simpa only [LinearMap.comp_apply, AlgEquiv.toLinearMap_apply, AlgEquiv.symm_apply_apply,
    AlgEquiv.apply_symm_apply, e, z] using hvalue

end Binary

/-- For a diagonal form of positive dimension at most two, the only
remaining input to the Mathieu conclusion is Laurent necessity. -/
theorem algebraicTrace_lowRank_isMathieu (N : ℕ) (hN : 0 < N) (hN₂ : N ≤ 2)
    (d : Fin N → Kˣ) (dvk : LowRankLaurentTheorem) :
    letI : Nonempty (Fin N) := ⟨⟨0, hN⟩⟩
    IsMathieu (algebraicTrace d) := by
  interval_cases N
  · have h : d = fun _ ↦ d 0 := by funext i; fin_cases i; rfl
    exact (congrArg (fun w : Fin 1 → Kˣ ↦ IsMathieu (algebraicTrace w)) h).mpr
      (algebraicTrace_unary_isMathieu (d 0) (dvk 1 (Or.inl rfl)))
  · have h : d = binaryWeights (d 0) (d 1) := by funext i; fin_cases i <;> rfl
    exact (congrArg (fun w : Fin 2 → Kˣ ↦ IsMathieu (algebraicTrace w)) h).mpr
      (algebraicTrace_binary_isMathieu (d 0) (d 1) (dvk 2 (Or.inr rfl)))

/-- In characteristic zero a unary quadratic complement has a Mathieu
trace kernel, using the checked one-variable Laurent theorem. -/
theorem algebraicTrace_unary_mathieu (a : Kˣ) :
    IsMathieu (algebraicTrace (fun _ : Fin 1 ↦ a)) :=
  algebraicTrace_unary_isMathieu a newtonLaurentNecessity_fin_one

/-- The only remaining Laurent input for the dimension classification
is the binary case. -/
theorem lowRankLaurentTheorem_of_binary (binary : NewtonLaurentNecessity ℂ (Fin 2)) :
    LowRankLaurentTheorem := by
  intro n hn
  rcases hn with rfl | rfl
  · exact newtonLaurentNecessity_fin_one
  · exact binary

end MathieuQuadratic
