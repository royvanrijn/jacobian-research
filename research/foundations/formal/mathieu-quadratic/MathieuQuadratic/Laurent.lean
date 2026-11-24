import MathieuQuadratic.BaseChange
import Mathlib.Algebra.MonoidAlgebra.Basic
import Mathlib.Algebra.MonoidAlgebra.Support
import Mathlib.Tactic

/-!
# Laurent constant-term safety from separation

This is the marker argument in Section 3 of the dimension classification.
The separation hypothesis in the corollary is the Laurent theorem's
necessity direction.
-/

open scoped Pointwise TensorProduct

namespace MathieuQuadratic

variable {K : Type*} [Field K] {G : Type*} [AddCommGroup G]

/-- The constant coefficient of the Laurent polynomial. -/
noncomputable def laurentCT : AddMonoidAlgebra K G →ₗ[K] K where
  toFun f := f.coeff 0
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

/-- Normalization and vanishing on nonconstant monomials characterize the
Laurent constant coefficient functional. -/
theorem eq_laurentCT_of_monomials (trace : AddMonoidAlgebra K G →ₗ[K] K)
    (hn : trace 1 = 1) (hm : ∀ g : G, g ≠ 0 → trace (AddMonoidAlgebra.single g 1) = 0) :
    trace = laurentCT := by
  classical
  apply LinearMap.ext
  intro f
  induction f using AddMonoidAlgebra.induction_linear with
  | zero => simp
  | add f g hf hg => simp only [map_add, hf, hg]
  | single g c =>
      rw [show AddMonoidAlgebra.single g c = c • (AddMonoidAlgebra.single g 1) by simp,
        map_smul, map_smul]
      congr 1
      by_cases hg : g = 0
      · subst g
        simpa [laurentCT, AddMonoidAlgebra.one_def] using hn
      · rw [hm g hg]
        simp [laurentCT, AddMonoidAlgebra.coeff_single, hg]


theorem supportLowerBound_mul (ℓ : G →+ ℤ) (f h : AddMonoidAlgebra K G)
    (b c : ℤ) (hf : ∀ g ∈ f.coeff.support, b ≤ ℓ g)
    (hh : ∀ g ∈ h.coeff.support, c ≤ ℓ g) :
    ∀ g ∈ (f * h).coeff.support, b + c ≤ ℓ g := by
  classical
  intro g hg
  obtain ⟨g₁, hg₁, g₂, hg₂, rfl⟩ :=
    Finset.mem_add.mp (AddMonoidAlgebra.support_coeff_mul_subset f h hg)
  rw [map_add]
  exact add_le_add (hf _ hg₁) (hh _ hg₂)

private theorem supportLowerBound_pow (ℓ : G →+ ℤ) (f : AddMonoidAlgebra K G)
    (hf : ∀ g ∈ f.coeff.support, 1 ≤ ℓ g) (n : ℕ) :
    ∀ g ∈ (f ^ n).coeff.support, (n : ℤ) ≤ ℓ g := by
  classical
  induction n with
  | zero =>
      intro g hg
      have hg0 : g = 0 := by
        by_contra hn
        simp [AddMonoidAlgebra.one_def, AddMonoidAlgebra.coeff_single, hn] at hg
      simp [hg0]
  | succ n hn =>
      simpa only [pow_succ, Nat.cast_succ] using supportLowerBound_mul ℓ (f ^ n) f n 1 hn hf

/-- An integral covector positive on the support implies purity. -/
theorem pure_of_positive_support (f : AddMonoidAlgebra K G) (ℓ : G →+ ℤ)
    (h : ∀ g ∈ f.coeff.support, 0 < ℓ g) : Pure laurentCT f := by
  classical
  intro n hn
  by_contra hne
  have hz : (0 : G) ∈ (f ^ n).coeff.support := Finsupp.mem_support_iff.mpr hne
  have hb := supportLowerBound_pow ℓ f (fun g hg ↦ h g hg) n 0 hz
  rw [map_zero] at hb
  omega

/-- Every fixed Laurent marker eventually has zero mixed constant terms. -/
theorem safe_of_positive_support (f : AddMonoidAlgebra K G) (ℓ : G →+ ℤ)
    (hf : ∀ g ∈ f.coeff.support, 0 < ℓ g) : Safe laurentCT f := by
  classical
  intro marker
  let cutoff := (∑ g ∈ marker.coeff.support, (ℓ g).natAbs) + 1
  have hmarker (g : G) (hg : g ∈ marker.coeff.support) : 1 - (cutoff : ℤ) ≤ ℓ g := by
    have h := Finset.single_le_sum (f := fun g ↦ (ℓ g).natAbs) (fun _ _ ↦ Nat.zero_le _) hg
    dsimp [cutoff]
    omega
  refine ⟨cutoff, fun n hn ↦ ?_⟩
  by_contra hne
  have hz : (0 : G) ∈ (marker * f ^ n).coeff.support := Finsupp.mem_support_iff.mpr hne
  have hb := supportLowerBound_mul ℓ marker (f ^ n) (1 - (cutoff : ℤ)) n hmarker
    (supportLowerBound_pow ℓ f (fun g hg ↦ hf g hg) n) 0 hz
  rw [map_zero] at hb
  omega

/-- The necessity input asserts that a pure Laurent polynomial has support
strictly separated from zero by an integral covector. The empty support
also satisfies this condition. -/
def LaurentSeparation (K G : Type*) [Field K] [AddCommGroup G] : Prop :=
  ∀ f : AddMonoidAlgebra K G, Pure laurentCT f →
    ∃ ℓ : G →+ ℤ, ∀ g ∈ f.coeff.support, 0 < ℓ g

theorem isMathieu_laurentCT (necessity : LaurentSeparation K G) :
    IsMathieu (laurentCT (K := K) (G := G)) := by
  intro f hp
  obtain ⟨ℓ, h⟩ := necessity f hp
  exact safe_of_positive_support f ℓ h

/-- A constant-term Laurent model transports the Mathieu property. -/
theorem isMathieu_of_laurent_equiv {A : Type*} [CommRing A] [Algebra K A]
    (trace : A →ₗ[K] K) (coordinates : A ≃ₐ[K] AddMonoidAlgebra K G)
    (identifies : ∀ f, trace f = laurentCT (coordinates f))
    (necessity : LaurentSeparation K G) : IsMathieu trace := by
  intro f hp
  have hp' : Pure laurentCT (coordinates f) := by
    intro n hn
    rw [← map_pow, ← identifies]
    exact hp n hn
  intro g
  obtain ⟨cutoff, hc⟩ := isMathieu_laurentCT necessity _ hp' (coordinates g)
  refine ⟨cutoff, fun n hn ↦ ?_⟩
  rw [identifies, map_mul, map_pow]
  exact hc n hn

theorem isMathieu_of_laurent_equiv_trace_eq {A : Type*} [CommRing A] [Algebra K A]
    (trace : A →ₗ[K] K) (coordinates : A ≃ₐ[K] AddMonoidAlgebra K G)
    (identifies : trace.comp coordinates.symm.toLinearMap = laurentCT)
    (necessity : LaurentSeparation K G) : IsMathieu trace := by
  apply isMathieu_of_laurent_equiv trace coordinates _ necessity
  intro f
  rw [← identifies, LinearMap.comp_apply, AlgEquiv.toLinearMap_apply, AlgEquiv.symm_apply_apply]

/-- A Laurent model after a splitting extension gives the Mathieu
property over the original field if its trace is constant coefficient. -/
theorem isMathieu_of_split_laurent {A L : Type*} [CommRing A] [Algebra K A]
    [Field L] [Algebra K L] (trace : A →ₗ[K] K)
    (coordinates : L ⊗[K] A ≃ₐ[L] AddMonoidAlgebra L G)
    (identifies : ∀ f, extendTrace L trace f = laurentCT (coordinates f))
    (necessity : LaurentSeparation L G) : IsMathieu trace :=
  isMathieu_of_extended L trace (isMathieu_of_laurent_equiv _ coordinates identifies necessity)

end MathieuQuadratic
