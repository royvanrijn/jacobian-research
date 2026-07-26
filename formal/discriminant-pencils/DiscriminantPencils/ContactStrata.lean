/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import DiscriminantPencils.TangentPencil

/-!
# Contact-pattern factorization families

This file formalizes the factorization in equation (3.3) and the numerical
dimension calculation (3.4) for the three bad contact patterns.

The geometric theorem that the closure of the image of such a family has
dimension at most the source dimension is not present in mathlib; this file
does not postulate it.
-/

open Polynomial

namespace DiscriminantPencils

variable {𝕜 : Type*} [Field 𝕜]

/-- The product `∏ᵢ (X - rᵢ)^μᵢ` attached to a contact pattern. -/
noncomputable def contactFactor {ι : Type*} [Fintype ι]
    (r : ι → 𝕜) (μ : ι → ℕ) : 𝕜[X] :=
  ∏ i, (X - C (r i)) ^ μ i

/-- The polynomial family `M_μ Q + ℓ` in equation (3.3). -/
noncomputable def contactFamily {ι : Type*} [Fintype ι]
    (r : ι → 𝕜) (μ : ι → ℕ) (Q : 𝕜[X]) (a b : 𝕜) : 𝕜[X] :=
  contactFactor r μ * Q + C a * X + C b

@[simp] theorem contactFamily_sub_line {ι : Type*} [Fintype ι]
    (r : ι → 𝕜) (μ : ι → ℕ) (Q : 𝕜[X]) (a b : 𝕜) :
    contactFamily r μ Q a b - (C a * X + C b) = contactFactor r μ * Q := by
  simp [contactFamily]

/-- Every marked factor divides the graph-minus-line polynomial. -/
theorem marked_contact_dvd {ι : Type*} [Fintype ι]
    (r : ι → 𝕜) (μ : ι → ℕ) (Q : 𝕜[X]) (a b : 𝕜) (i : ι) :
    (X - C (r i)) ^ μ i ∣
      contactFamily r μ Q a b - (C a * X + C b) := by
  classical
  rw [contactFamily_sub_line]
  apply dvd_mul_of_dvd_left
  exact Finset.dvd_prod_of_mem (fun j ↦ (X - C (r j)) ^ μ j) (Finset.mem_univ i)

/-- The source dimension `k + (n - m + 1)` from equation (3.4). -/
def contactFamilyDimension (n : ℕ) (μ : List ℕ) : ℕ :=
  μ.length + (n - μ.sum + 1)

theorem contactFamilyDimension_four (n : ℕ) (hn : 4 ≤ n) :
    contactFamilyDimension n [4] = n - 2 := by
  simp [contactFamilyDimension]
  omega

theorem contactFamilyDimension_three_two (n : ℕ) (hn : 5 ≤ n) :
    contactFamilyDimension n [3, 2] = n - 2 := by
  simp [contactFamilyDimension]
  omega

theorem contactFamilyDimension_two_two_two (n : ℕ) (hn : 6 ≤ n) :
    contactFamilyDimension n [2, 2, 2] = n - 2 := by
  simp [contactFamilyDimension]
  omega

/-- Each possible bad-pattern family has one less parameter than the
`n - 1` dimensional coefficient space. -/
theorem bad_contact_family_codimension_one (n : ℕ)
    (μ : List ℕ)
    (hμ : (μ = [4] ∧ 4 ≤ n) ∨ (μ = [3, 2] ∧ 5 ≤ n) ∨
      (μ = [2, 2, 2] ∧ 6 ≤ n)) :
    contactFamilyDimension n μ + 1 = n - 1 := by
  rcases hμ with ⟨rfl, h⟩ | ⟨rfl, h⟩ | ⟨rfl, h⟩
  · rw [contactFamilyDimension_four n h]
    omega
  · rw [contactFamilyDimension_three_two n h]
    omega
  · rw [contactFamilyDimension_two_two_two n h]
    omega

end DiscriminantPencils
