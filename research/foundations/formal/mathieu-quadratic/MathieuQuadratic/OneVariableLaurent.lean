import MathieuQuadratic.CoefficientField
import MathieuQuadratic.Newton
import MathieuQuadratic.External.DvK.OneVariableDvK

/-!
# Laurent necessity in rank one

The verified one-variable root-sum argument supplies the necessity
direction for the rank-one quadratic complement.
-/

namespace MathieuQuadratic

open scoped BigOperators

theorem laurentSeparation_int : LaurentSeparation ℂ ℤ := by
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

/-- Strictly one-sided support, with the zero Laurent polynomial allowed. -/
def LaurentOneSided {K : Type*} [Field K] (f : AddMonoidAlgebra K ℤ) : Prop :=
  (∀ i ∈ f.coeff.support, 0 < i) ∨ (∀ i ∈ f.coeff.support, i < 0)

private theorem coeff_mul_at_minimum {K : Type*} [Field K] (w : ℤ →+ ℤ)
    (hw : Function.Injective w) (f g : AddMonoidAlgebra K ℤ) {e d : ℤ}
    (he : e ∈ f.coeff.support) (hd : d ∈ g.coeff.support)
    (hf : ∀ i ∈ f.coeff.support, w e ≤ w i)
    (hg : ∀ j ∈ g.coeff.support, w d ≤ w j) :
    (f * g).coeff (e + d) = f.coeff e * g.coeff d := by
  classical
  rw [AddMonoidAlgebra.coeff_mul, Finsupp.sum,
    Finset.sum_eq_single_of_mem e he, Finsupp.sum, Finset.sum_eq_single_of_mem d hd, if_pos rfl]
  · intro j hj hje
    rw [if_neg fun h ↦ hje (add_left_cancel h)]
  · intro i hi hie
    apply Finset.sum_eq_zero
    intro j hj
    change (if i + j = e + d then _ else (0 : K)) = 0
    rw [if_neg]
    intro h
    have hsum := congrArg w h
    simp only [map_add] at hsum
    have hwi := hf i hi
    have hwj := hg j hj
    exact hie (hw (by omega))

private theorem laurent_support_pos_of_pow {K : Type*} [Field K]
    (w : ℤ →+ ℤ) (hw : Function.Injective w) (f : AddMonoidAlgebra K ℤ)
    {n : ℕ} (h : ∀ i ∈ (f ^ n).coeff.support, 0 < w i) :
    ∀ i ∈ f.coeff.support, 0 < w i := by
  intro i hi
  obtain ⟨e, he, hmin⟩ := Finset.exists_min_image f.coeff.support w ⟨i, hi⟩
  have hbound (m : ℕ) :
      ∀ j ∈ (f ^ m).coeff.support, w ((m : ℤ) * e) ≤ w j := by
    induction m with
    | zero =>
      intro j hj
      have hj0 : j = 0 := by
        rw [pow_zero, AddMonoidAlgebra.one_def, AddMonoidAlgebra.coeff_single,
          Finsupp.support_single _ one_ne_zero, Finset.mem_singleton] at hj
        exact hj
      simp [hj0]
    | succ m ih =>
      rw [pow_succ, Nat.cast_succ, add_mul, one_mul, map_add]
      exact supportLowerBound_mul w (f ^ m) f _ _ ih hmin
  have hpower (m : ℕ) : ((m : ℤ) * e) ∈ (f ^ m).coeff.support := by
    induction m with
    | zero =>
      rw [pow_zero, AddMonoidAlgebra.one_def, AddMonoidAlgebra.coeff_single,
        Finsupp.support_single _ one_ne_zero, Finset.mem_singleton]
      exact Int.zero_mul e
    | succ m ih =>
      rw [pow_succ, Nat.cast_succ, add_mul, one_mul, Finsupp.mem_support_iff,
        coeff_mul_at_minimum w hw _ _ ih he (hbound m) hmin]
      exact mul_ne_zero (Finsupp.mem_support_iff.mp ih) (Finsupp.mem_support_iff.mp he)
  have hwe := h _ (hpower n)
  have hmap : w ((n : ℤ) * e) = (n : ℤ) * w e := w.map_zsmul (n : ℤ) e
  rw [hmap] at hwe
  exact (pos_of_mul_pos_right hwe (by exact_mod_cast Nat.zero_le n)).trans_le (hmin i hi)

/-- Eventual vanishing of the unmarked constant terms already forces a
rank-one Laurent polynomial to be zero or strictly one-sided. -/
theorem laurentOneSided_of_eventual_ct_zero {K : Type*} [Field K] [CharZero K]
    (f : AddMonoidAlgebra K ℤ)
    (eventual : ∃ cutoff : ℕ, ∀ n ≥ cutoff, laurentCT (f ^ n) = 0) :
    LaurentOneSided f := by
  obtain ⟨cutoff, hz⟩ := eventual
  let n := cutoff + 1
  have hp : Pure laurentCT (f ^ n) := by
    intro m hm
    rw [← pow_mul]
    exact hz _ ((Nat.le_succ _).trans (Nat.le_mul_of_pos_right n hm))
  obtain ⟨w, hw⟩ := laurentSeparation_of_complex laurentSeparation_int (f ^ n) hp
  have hform (i : ℤ) : w i = i * w 1 := by
    simpa only [Int.zsmul_eq_mul, mul_one] using w.map_zsmul i 1
  by_cases hf : f = 0
  · exact Or.inl (by simp [hf])
  have hnonzero : w 1 ≠ 0 := by
    intro h
    obtain ⟨i, hi⟩ := Finsupp.support_nonempty_iff.mpr
      (show (f ^ n).coeff ≠ 0 from fun hzero ↦
        pow_ne_zero n hf (AddMonoidAlgebra.coeff_injective (by simpa using hzero)))
    simpa only [hform i, h, mul_zero, lt_irrefl] using hw i hi
  have hinj : Function.Injective w := by
    intro i j hij
    rw [hform i, hform j] at hij
    exact mul_right_cancel₀ hnonzero hij
  have h := laurent_support_pos_of_pow w hinj f hw
  have hsign (i : ℤ) (hi : i ∈ f.coeff.support) := mul_pos_iff.mp (show 0 < i * w 1 by rw [← hform i]; exact h i hi)
  rcases lt_or_gt_of_ne hnonzero with hneg | hpos
  · exact Or.inr fun i hi ↦ ((hsign i hi).resolve_left (fun hp ↦ hneg.not_gt hp.2)).1
  · exact Or.inl fun i hi ↦ ((hsign i hi).resolve_right (fun hn ↦ hpos.not_gt hn.2)).1

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
