import MathieuQuadratic.LowDimension
import MathieuQuadratic.LocalizationBaseChange
import MathieuQuadratic.CoefficientField

/-!
# Trace identification in rank one and rank two

The normalized trace is constant coefficient if it is invariant under the
indicated diagonal torus actions on monomials in split coordinates.
-/

open scoped TensorProduct

namespace MathieuQuadratic

variable {K : Type*} [Field K] [CharZero K]

private theorem two_zpow_ne_one (n : ℤ) (hn : n ≠ 0) : (2 : K) ^ n ≠ 1 := by
  intro h
  have he : algebraMap ℚ K ((2 : ℚ) ^ n) = algebraMap ℚ K 1 := by simpa using h
  exact hn ((zpow_right_strictMono₀ (by norm_num : (1 : ℚ) < 2)).injective
    (by simpa using (RingHom.injective (algebraMap ℚ K)) he : (2 : ℚ) ^ n = 2 ^ (0 : ℤ)))

private theorem trace_monomial_eq_zero {J : Type*} (trace : LaurentAlgebra K J →ₗ[K] K)
    (g : J → ℤ) (n : ℤ) (hn : n ≠ 0)
    (h : trace ((2 : K) ^ n • AddMonoidAlgebra.single g 1) = trace (AddMonoidAlgebra.single g 1)) :
    trace (AddMonoidAlgebra.single g 1) = 0 := by
  rw [map_smul, smul_eq_mul] at h
  have hz : ((2 : K) ^ n - 1) * trace (AddMonoidAlgebra.single g 1) = 0 := by
    rw [sub_mul, one_mul, h, sub_self]
  exact (mul_eq_zero.mp hz).resolve_left (sub_ne_zero.mpr (two_zpow_ne_one n hn))

theorem unaryTrace_eq_laurentCT (trace : LaurentAlgebra K (Fin 1) →ₗ[K] K)
    (hn : trace 1 = 1)
    (scalar : ∀ g : Fin 1 → ℤ, trace ((2 : K) ^ g 0 • AddMonoidAlgebra.single g 1) =
      trace (AddMonoidAlgebra.single g 1)) :
    trace = laurentCT := by
  apply eq_laurentCT_of_monomials trace hn
  intro g hg
  have h0 : g 0 ≠ 0 := by
    intro he
    apply hg
    funext j
    have hj : j = 0 := Subsingleton.elim _ _
    simpa only [hj, Pi.zero_apply] using he
  exact trace_monomial_eq_zero trace g (g 0) h0 (scalar g)

theorem binaryTrace_eq_laurentCT (trace : LaurentAlgebra K (Fin 2) →ₗ[K] K)
    (hn : trace 1 = 1)
    (scalar : ∀ g : Fin 2 → ℤ, trace ((2 : K) ^ (g 0 + g 1) • AddMonoidAlgebra.single g 1) =
      trace (AddMonoidAlgebra.single g 1))
    (rotation : ∀ g : Fin 2 → ℤ, trace ((2 : K) ^ (g 0 - g 1) • AddMonoidAlgebra.single g 1) =
      trace (AddMonoidAlgebra.single g 1)) :
    trace = laurentCT := by
  apply eq_laurentCT_of_monomials trace hn
  intro g hg
  by_cases hsum : g 0 + g 1 = 0
  · have hdiff : g 0 - g 1 ≠ 0 := by
      intro he
      apply hg
      funext j
      fin_cases j
      · change g 0 = 0
        omega
      · change g 1 = 0
        omega
    exact trace_monomial_eq_zero trace g _ hdiff (rotation g)
  · exact trace_monomial_eq_zero trace g _ hsum (scalar g)

omit [CharZero K] in
open MvPolynomial in
theorem map_binaryNorm {L : Type*} [Field L] [Algebra K L] (a b : Kˣ) :
    MvPolynomial.map (algebraMap K L) (binaryNorm a b) =
      binaryNorm (Units.map (algebraMap K L).toMonoidHom a)
        (Units.map (algebraMap K L).toMonoidHom b) := by
  simp only [binaryNorm, map_add, map_mul, map_pow, MvPolynomial.map_C, MvPolynomial.map_X]
  rfl

/-- The binary split Laurent model after base change to an algebraically
closed extension. -/
noncomputable def binaryLaurentAfterBaseChange (L : Type*) [Field L] [Algebra K L]
    [CharZero L] [IsAlgClosed L] (a b : Kˣ) :
    L ⊗[K] Localization.Away (binaryNorm a b) ≃ₐ[L] LaurentAlgebra L (Fin 2) :=
  (localizationBaseChange (L := L) (binaryNorm a b)).trans <|
    (map_binaryNorm (L := L) a b).symm ▸ Classical.choice (exists_binaryLaurent
      (Units.map (algebraMap K L).toMonoidHom a) (Units.map (algebraMap K L).toMonoidHom b))

theorem unary_mathieu_of_character_weights (a : Kˣ)
    (trace : Localization.Away (unaryNorm a) →ₗ[K] K) (hn : trace 1 = 1)
    (dvk : LaurentSeparation ℂ (Fin 1 → ℤ)) :
    let trace' := trace.comp (unaryLaurent a).symm.toLinearMap
    (∀ g : Fin 1 → ℤ, trace' ((2 : K) ^ g 0 • AddMonoidAlgebra.single g 1) =
      trace' (AddMonoidAlgebra.single g 1)) →
    IsMathieu trace := by
  dsimp only
  intro scalar
  have hnorm : (trace.comp (unaryLaurent a).symm.toLinearMap) 1 = 1 := by
    change trace ((unaryLaurent a).symm 1) = 1
    rw [map_one, hn]
  have ht : trace.comp (unaryLaurent a).symm.toLinearMap = laurentCT :=
    @unaryTrace_eq_laurentCT K _ _ _ hnorm scalar
  exact isMathieu_of_laurent_equiv_trace_eq trace (unaryLaurent a) ht
    (laurentSeparation_of_complex dvk)

omit [CharZero K] in
theorem binary_mathieu_of_character_weights {L : Type*} [Field L] [Algebra K L]
    [CharZero L] [IsAlgClosed L] (a b : Kˣ)
    (trace : Localization.Away (binaryNorm a b) →ₗ[K] K) (hn : trace 1 = 1)
    (dvk : LaurentSeparation ℂ (Fin 2 → ℤ)) :
    let trace' := (extendTrace L trace).comp
      (binaryLaurentAfterBaseChange L a b).symm.toLinearMap
    (∀ g : Fin 2 → ℤ, trace' ((2 : L) ^ (g 0 + g 1) • AddMonoidAlgebra.single g 1) =
      trace' (AddMonoidAlgebra.single g 1)) →
    (∀ g : Fin 2 → ℤ, trace' ((2 : L) ^ (g 0 - g 1) • AddMonoidAlgebra.single g 1) =
      trace' (AddMonoidAlgebra.single g 1)) →
    IsMathieu trace := by
  dsimp only
  intro scalar rotation
  let e := binaryLaurentAfterBaseChange L a b
  have h1 : extendTrace L trace 1 = 1 := by
    change extendTrace L trace (1 ⊗ₜ[K] 1) = 1
    rw [extendTrace_tmul, hn, map_one, one_mul]
  have ht := binaryTrace_eq_laurentCT ((extendTrace L trace).comp e.symm.toLinearMap)
    (by simpa using h1) scalar rotation
  exact isMathieu_of_extended L trace
    (isMathieu_of_laurent_equiv_trace_eq (extendTrace L trace) e ht
      (laurentSeparation_of_complex dvk))

end MathieuQuadratic
