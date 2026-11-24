import MathieuQuadratic.Newton
import MathieuQuadratic.External.DvK.OneVariableDvK

/-!
# Laurent necessity in rank one

The verified one-variable root-sum argument supplies the necessity
direction for the rank-one quadratic complement.
-/

namespace MathieuQuadratic

open scoped BigOperators

private theorem laurentSeparation_int : LaurentSeparation ℂ ℤ := by
  intro f hp
  have hzero : (0 : ℤ) ∉ f.coeff.support := by
    apply Finsupp.notMem_support_iff.mpr
    simpa [laurentCT] using hp 1 Nat.one_pos
  by_cases hpos : ∃ i ∈ f.coeff.support, 0 < i
  · have hn : ¬ ∃ i ∈ f.coeff.support, i < 0 := by
      intro hneg
      obtain ⟨m, hm, hc⟩ := DvKProof.hasNonzeroConstantPower_of_dvkCoefficientExtraction
        DvKProof.dvkCoefficientExtraction f hneg hpos
      exact hc (hp m (by omega))
    refine ⟨AddMonoidHom.id ℤ, fun i hi ↦ ?_⟩
    have h0 : i ≠ 0 := fun h ↦ hzero (h ▸ hi)
    exact lt_of_le_of_ne (not_lt.mp fun h ↦ hn ⟨i, hi, h⟩) h0.symm
  · refine ⟨-(AddMonoidHom.id ℤ), fun i hi ↦ ?_⟩
    have h0 : i ≠ 0 := fun h ↦ hzero (h ▸ hi)
    have hn : i < 0 := lt_of_le_of_ne (not_lt.mp fun h ↦ hpos ⟨i, hi, h⟩) h0
    simpa only [AddMonoidHom.neg_apply, AddMonoidHom.id_apply] using neg_pos.mpr hn

private def unaryExponents : (Fin 1 → ℤ) ≃+ ℤ where
  toFun g := g 0
  invFun n := fun _ ↦ n
  left_inv g := by funext i; fin_cases i; rfl
  right_inv _ := rfl
  map_add' _ _ := rfl

/-- Every pure rank-one complex Laurent polynomial has strictly
one-sided support. -/
theorem laurentSeparation_fin_one : LaurentSeparation ℂ (Fin 1 → ℤ) := by
  intro f hp
  let E := AddMonoidAlgebra.domCongr ℂ ℂ unaryExponents
  have hp' : Pure laurentCT (E f) := by
    intro m hm
    change ((E f) ^ m).coeff 0 = 0
    rw [← map_pow, AddMonoidAlgebra.coeff_domCongr, map_zero]
    exact hp m hm
  obtain ⟨ℓ, hℓ⟩ := laurentSeparation_int (E f) hp'
  refine ⟨ℓ.comp unaryExponents.toAddMonoidHom, fun g hg ↦ hℓ _ ?_⟩
  rw [Finsupp.mem_support_iff]
  change (AddMonoidAlgebra.domCongr ℂ ℂ unaryExponents f).coeff (unaryExponents g) ≠ 0
  rw [AddMonoidAlgebra.coeff_domCongr, unaryExponents.symm_apply_apply]
  exact Finsupp.mem_support_iff.mp hg

/-- Duistermaat–van der Kallen necessity in rank one. -/
theorem newtonLaurentNecessity_fin_one : NewtonLaurentNecessity ℂ (Fin 1) := by
  intro f _ hp
  obtain ⟨ℓ, hℓ⟩ := laurentSeparation_fin_one f hp
  exact origin_not_mem_newton_of_positive_support f.coeff.support ℓ hℓ

end MathieuQuadratic
