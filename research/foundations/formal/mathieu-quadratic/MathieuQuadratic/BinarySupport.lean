import MathieuQuadratic.LaurentTrace

/-!
# Support symmetry and the nonsplit binary scalar weights
-/

namespace MathieuQuadratic

def binarySwapExponent (g : Fin 2 → ℤ) : Fin 2 → ℤ :=
  g ∘ Equiv.swap (0 : Fin 2) 1

private theorem add_swapExponent (g : Fin 2 → ℤ) :
    g + binarySwapExponent g = (g 0 + g 1) • (fun _ : Fin 2 ↦ (1 : ℤ)) := by
  funext i
  fin_cases i
  · change g 0 + binarySwapExponent g 0 = (g 0 + g 1) • (1 : ℤ)
    simp [binarySwapExponent]
  · change g 1 + binarySwapExponent g 1 = (g 0 + g 1) • (1 : ℤ)
    simp [binarySwapExponent, add_comm]

private def scalarWeight : (Fin 2 → ℤ) →+ ℤ where
  toFun g := g 0 + g 1
  map_zero' := by simp
  map_add' g h := by simp [Pi.add_apply]; ring

/-- A swap-invariant support admits a strictly positive integral covector
exactly when all its scalar weights have one strict sign. -/
theorem pure_binary_iff_scalar_sign {K : Type*} [Field K]
    (necessity : LaurentSeparation K (Fin 2 → ℤ))
    (f : LaurentAlgebra K (Fin 2))
    (symmetric : ∀ g ∈ f.coeff.support, binarySwapExponent g ∈ f.coeff.support) :
    Pure laurentCT f ↔
      (∀ g ∈ f.coeff.support, 0 < g 0 + g 1) ∨
      (∀ g ∈ f.coeff.support, g 0 + g 1 < 0) := by
  constructor
  · intro hp
    obtain ⟨ℓ, hℓ⟩ := necessity f hp
    have hpos (g : Fin 2 → ℤ) (hg : g ∈ f.coeff.support) :
        0 < (g 0 + g 1) * ℓ (fun _ ↦ 1) := by
      have h : 0 < ℓ (g + binarySwapExponent g) := by
        rw [map_add]
        exact add_pos (hℓ g hg) (hℓ _ (symmetric g hg))
      rwa [add_swapExponent, map_zsmul, zsmul_eq_mul] at h
    by_cases hd : 0 < ℓ (fun _ ↦ 1)
    · exact Or.inl fun g hg ↦ ((mul_pos_iff.mp (hpos g hg)).resolve_right
        (fun hn ↦ (not_lt_of_gt hd) hn.2)).1
    · exact Or.inr fun g hg ↦ ((mul_pos_iff.mp (hpos g hg)).resolve_left
        (fun hp' ↦ hd hp'.2)).1
  · rintro (h | h)
    · exact pure_of_positive_support f scalarWeight h
    · exact pure_of_positive_support f (-scalarWeight) (fun g hg ↦ neg_pos.mpr (h g hg))

end MathieuQuadratic
