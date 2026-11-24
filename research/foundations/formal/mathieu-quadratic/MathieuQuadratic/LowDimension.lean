import MathieuQuadratic.SplitLocalization
import Mathlib.FieldTheory.IsAlgClosed.Basic

/-!
# The unary and split binary quadratic complements
-/

open scoped BigOperators

namespace MathieuQuadratic

open MvPolynomial

variable {K : Type*} [Field K] [CharZero K]

noncomputable def unaryNorm (a : Kˣ) : MvPolynomial (Fin 1) K :=
  C (a : K) * X 0 ^ 2

noncomputable def unaryLaurent (a : Kˣ) :
    Localization.Away (unaryNorm a) ≃ₐ[K] LaurentAlgebra K (Fin 1) := by
  have he : unaryNorm a = splitNorm a (fun _ : Fin 1 ↦ 2) := by
    simp [unaryNorm, splitNorm]
  exact he ▸ splitComplementEquivLaurent a (fun _ : Fin 1 ↦ 2) (by simp)

noncomputable def binaryNorm (a b : Kˣ) : MvPolynomial (Fin 2) K :=
  C (a : K) * X 0 ^ 2 + C (b : K) * X 1 ^ 2

private noncomputable def expandSplit (d : Kˣ) :
    MvPolynomial (Fin 2) K →ₐ[K] MvPolynomial (Fin 2) K :=
  aeval fun i ↦ if i = 0 then X 0 + C (d : K) * X 1 else X 0 - C (d : K) * X 1

private noncomputable def extractSplit (d : Kˣ) :
    MvPolynomial (Fin 2) K →ₐ[K] MvPolynomial (Fin 2) K :=
  aeval fun i ↦ if i = 0 then C ((2 : K)⁻¹) * (X 0 + X 1)
    else C ((2 * (d : K))⁻¹) * (X 0 - X 1)

private theorem C_two_inv :
    C ((2 : K)⁻¹) * 2 = (1 : MvPolynomial (Fin 2) K) := by
  have ht : C (2 : K) = (2 : MvPolynomial (Fin 2) K) := map_ofNat _ _
  rw [← ht, ← map_mul, inv_mul_cancel₀ two_ne_zero, map_one]

private theorem C_two_mul_inv (d : Kˣ) :
    C ((2 * (d : K))⁻¹) * (2 * C (d : K)) = (1 : MvPolynomial (Fin 2) K) := by
  have ht : C (2 : K) = (2 : MvPolynomial (Fin 2) K) := map_ofNat _ _
  rw [← ht, ← map_mul, ← map_mul, inv_mul_cancel₀
    (mul_ne_zero two_ne_zero (Units.ne_zero d)), map_one]

private theorem expand_extract (d : Kˣ) :
    (expandSplit d).comp (extractSplit d) = AlgHom.id K (MvPolynomial (Fin 2) K) := by
  apply MvPolynomial.algHom_ext
  intro i
  fin_cases i
  · change expandSplit d (extractSplit d (X (0 : Fin 2))) = X 0
    simp only [expandSplit, extractSplit, aeval_X, eq_self, if_true,
      if_neg (by decide : (1 : Fin 2) ≠ 0), map_mul, map_add, aeval_C, algebraMap_eq]
    trans (C ((2 : K)⁻¹) * 2) * (X 0 : MvPolynomial (Fin 2) K)
    · ring
    · rw [C_two_inv, one_mul]
  · change expandSplit d (extractSplit d (X (1 : Fin 2))) = X 1
    simp only [expandSplit, extractSplit, aeval_X, eq_self, if_true,
      if_neg (by decide : (1 : Fin 2) ≠ 0), map_mul, map_sub, aeval_C, algebraMap_eq]
    trans (C ((2 * (d : K))⁻¹) * (2 * C (d : K))) * (X 1 : MvPolynomial (Fin 2) K)
    · ring
    · rw [C_two_mul_inv, one_mul]

private theorem extract_expand (d : Kˣ) :
    (extractSplit d).comp (expandSplit d) = AlgHom.id K (MvPolynomial (Fin 2) K) := by
  have hc : C (d : K) * C ((2 * (d : K))⁻¹) =
      (C ((2 : K)⁻¹) : MvPolynomial (Fin 2) K) := by
    rw [← map_mul]
    congr 1
    rw [mul_inv, ← mul_assoc, mul_comm (d : K) ((2 : K)⁻¹), mul_assoc,
      mul_inv_cancel₀ (Units.ne_zero d), mul_one]
  apply MvPolynomial.algHom_ext
  intro i
  fin_cases i
  · change extractSplit d (expandSplit d (X (0 : Fin 2))) = X 0
    simp only [expandSplit, extractSplit, aeval_X, eq_self, if_true,
      if_neg (by decide : (1 : Fin 2) ≠ 0), map_mul, map_add, aeval_C, algebraMap_eq]
    rw [← mul_assoc, hc]
    trans (C ((2 : K)⁻¹) * 2) * (X 0 : MvPolynomial (Fin 2) K)
    · ring
    · rw [C_two_inv, one_mul]
  · change extractSplit d (expandSplit d (X (1 : Fin 2))) = X 1
    simp only [expandSplit, extractSplit, aeval_X, eq_self, if_true,
      if_neg (by decide : (1 : Fin 2) ≠ 0), map_mul, map_sub, aeval_C, algebraMap_eq]
    rw [← mul_assoc, hc]
    trans (C ((2 : K)⁻¹) * 2) * (X 1 : MvPolynomial (Fin 2) K)
    · ring
    · rw [C_two_inv, one_mul]

/-- Over a splitting field the two Laurent variables are x+δy and x−δy. -/
noncomputable def binaryLinearCoordinates (d : Kˣ) :
    MvPolynomial (Fin 2) K ≃ₐ[K] MvPolynomial (Fin 2) K :=
  AlgEquiv.ofAlgHom (expandSplit d) (extractSplit d) (expand_extract d) (extract_expand d)

theorem binaryLinearCoordinates_norm (a b d : Kˣ) (hd : (d : K) ^ 2 = -(b : K) / (a : K)) :
    binaryLinearCoordinates d (splitNorm a (fun _ : Fin 2 ↦ 1)) = binaryNorm a b := by
  have h : C (a : K) * C (d : K) ^ 2 = (- C (b : K) : MvPolynomial (Fin 2) K) := by
    rw [← map_pow, ← map_mul, hd, mul_comm, div_mul_cancel₀ _ (Units.ne_zero a), map_neg]
  change expandSplit d (splitNorm a (fun _ : Fin 2 ↦ 1)) = binaryNorm a b
  simp [splitNorm, binaryNorm, expandSplit, Fin.prod_univ_two, algebraMap_eq]
  trans C (a : K) * (X 0 : MvPolynomial (Fin 2) K) ^ 2 -
    (C (a : K) * C (d : K) ^ 2) * X 1 ^ 2
  · ring
  · rw [h]; ring

noncomputable def binaryLaurent (a b d : Kˣ) (hd : (d : K) ^ 2 = -(b : K) / (a : K)) :
    Localization.Away (binaryNorm a b) ≃ₐ[K] LaurentAlgebra K (Fin 2) :=
  let e : SplitComplement a (fun _ : Fin 2 ↦ 1) ≃ₐ[K] Localization.Away (binaryNorm a b) :=
    IsLocalization.algEquivOfAlgEquiv _ _ (binaryLinearCoordinates d) (by
      rw [Submonoid.map_powers, binaryLinearCoordinates_norm a b d hd])
  e.symm.trans (splitComplementEquivLaurent a (fun _ : Fin 2 ↦ 1) (by simp))

theorem exists_binaryLaurent [IsAlgClosed K] (a b : Kˣ) :
    Nonempty (Localization.Away (binaryNorm a b) ≃ₐ[K] LaurentAlgebra K (Fin 2)) := by
  obtain ⟨d, hd⟩ := IsAlgClosed.exists_pow_nat_eq (-(b : K) / (a : K)) (n := 2) (by decide)
  have hd0 : d ≠ 0 := by
    intro h
    exact (div_ne_zero (neg_ne_zero.mpr (Units.ne_zero b)) (Units.ne_zero a))
      (by simpa [h] using hd.symm)
  exact ⟨binaryLaurent a b (Units.mk0 d hd0) hd⟩

end MathieuQuadratic
