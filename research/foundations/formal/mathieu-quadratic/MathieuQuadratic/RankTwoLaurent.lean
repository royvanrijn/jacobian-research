import MathieuQuadratic.LaurentFace
import MathieuQuadratic.OneVariableLaurent
import Mathlib.Analysis.Normed.Affine.AddTorsorBases

/-!
# The boundary case of planar Laurent necessity

A supporting line through the origin restricts the constant-term
calculation to its intersection with the exponent lattice. One
coordinate projects that lattice injectively to the integers, where
the checked one-variable Laurent theorem applies.
-/

open scoped BigOperators

namespace MathieuQuadratic

private abbrev Plane := Fin 2 → ℝ
private abbrev Lattice := Fin 2 → ℤ
private abbrev Laurent₂ := AddMonoidAlgebra ℂ Lattice

private def realExponentHom : Lattice →+ Plane where
  toFun := realExponent
  map_zero' := by ext i; simp [realExponent]
  map_add' a b := by ext i; simp [realExponent]

private def weight (F : Plane →ₗ[ℝ] ℝ) : Lattice →+ ℝ :=
  F.toAddMonoidHom.comp realExponentHom

private def projection (i : Fin 2) : Lattice →+ ℤ where
  toFun g := g i
  map_zero' := rfl
  map_add' _ _ := rfl

private theorem linearForm_coordinates (F : Plane →ₗ[ℝ] ℝ) (x : Plane) :
    F x = F (![1, 0] : Plane) * x 0 + F (![0, 1] : Plane) * x 1 := by
  have he : x = x 0 • (![1, 0] : Plane) + x 1 • (![0, 1] : Plane) := by
    funext j
    fin_cases j <;> simp
  calc
    F x = F (x 0 • (![1, 0] : Plane) + x 1 • (![0, 1] : Plane)) := congrArg F he
    _ = _ := by rw [map_add, map_smul, map_smul, smul_eq_mul, smul_eq_mul]; ring

private theorem exists_projection_injOn_kernel (F : Plane →ₗ[ℝ] ℝ) (hF : F ≠ 0) :
    ∃ i : Fin 2, Set.InjOn (projection i) {g | weight F g = 0} := by
  by_cases h₁ : F (![0, 1] : Plane) = 0
  · have h₀ : F (![1, 0] : Plane) ≠ 0 := by
      intro h₀
      apply hF
      apply LinearMap.ext
      intro x
      simp [linearForm_coordinates, h₀, h₁]
    refine ⟨1, ?_⟩
    intro x hx y hy hxy
    change x 1 = y 1 at hxy
    change F (realExponent x) = 0 at hx
    change F (realExponent y) = 0 at hy
    rw [linearForm_coordinates] at hx hy
    simp only [h₁, zero_mul, add_zero] at hx hy
    have h := (mul_eq_zero.mp hx).resolve_left h₀
    have h' := (mul_eq_zero.mp hy).resolve_left h₀
    funext j
    fin_cases j
    · change x 0 = y 0
      exact Int.cast_injective (h.trans h'.symm)
    · exact hxy
  · refine ⟨0, ?_⟩
    intro x hx y hy hxy
    change x 0 = y 0 at hxy
    change F (realExponent x) = 0 at hx
    change F (realExponent y) = 0 at hy
    rw [linearForm_coordinates] at hx hy
    have hxy' : realExponent x 0 = realExponent y 0 := by simp [realExponent, hxy]
    rw [hxy'] at hx
    have hv : F (![0, 1] : Plane) * realExponent x 1 =
        F (![0, 1] : Plane) * realExponent y 1 := by linarith
    have h := mul_left_cancel₀ h₁ hv
    funext j
    fin_cases j
    · exact hxy
    · change x 1 = y 1
      exact Int.cast_injective h

private theorem exists_separator_on_zeroFace (F : Plane →ₗ[ℝ] ℝ) (hF : F ≠ 0)
    (f : Laurent₂) (hp : Pure laurentCT f) (hn : ∀ g ∈ f.coeff.support, 0 ≤ F (realExponent g)) :
    ∃ ℓ : Lattice →+ ℤ, ∀ g ∈ f.coeff.support, F (realExponent g) = 0 → 0 < ℓ g := by
  obtain ⟨i, hi⟩ := exists_projection_injOn_kernel F hF
  let w := weight F
  let p := projection i
  let z := laurentZeroFace w f
  let h := AddMonoidAlgebra.mapDomainRingHom ℂ p z
  have hz : ∀ g ∈ z.coeff.support, w g = 0 :=
    fun g hg ↦ ((mem_support_laurentZeroFace w f g).mp hg).2
  have hp' : Pure laurentCT h := by
    intro m hm
    rw [laurentCT_pow_mapDomain_of_weight_zero w p hi z hz m,
      ← laurentCT_pow_zeroFace w f hn m]
    exact hp m hm
  obtain ⟨ℓ, hℓ⟩ := laurentSeparation_int h hp'
  refine ⟨ℓ.comp p, fun g hg hFg ↦ hℓ (p g) ?_⟩
  change w g = 0 at hFg
  have hgz : g ∈ z.coeff.support := (mem_support_laurentZeroFace w f g).mpr ⟨hg, hFg⟩
  apply Finsupp.mem_support_iff.mpr
  change (AddMonoidAlgebra.mapDomain p z).coeff (p g) ≠ 0
  rw [coeff_mapDomain_of_weight_zero w p hi z hz g hFg]
  exact Finsupp.mem_support_iff.mp hgz

private theorem zero_not_mem_newton_of_positive_form (s : Finset Lattice) (F : Plane →ₗ[ℝ] ℝ)
    (hF : ∀ g ∈ s, 0 < F (realExponent g)) : 0 ∉ newtonPolytope s := by
  have hs : newtonPolytope s ⊆ {x | 0 < F x} := by
    apply convexHull_min _ (convex_halfSpace_gt F.isLinear 0)
    rintro _ ⟨g, hg, rfl⟩
    exact hF g hg
  intro h
  have h' : 0 < F 0 := hs h
  simp at h'

/-- A pure planar Laurent polynomial with nonnegative support weights
for a nonzero linear functional has Newton polygon disjoint from zero.
The strict separator comes from a one-variable calculation on the zero
face and a finite perturbation. -/
theorem planarLaurent_pure_avoids_newton_of_supporting_form
    (F : Plane →ₗ[ℝ] ℝ) (hF : F ≠ 0) (f : Laurent₂) (hp : Pure laurentCT f)
    (hn : ∀ g ∈ f.coeff.support, 0 ≤ F (realExponent g)) :
    0 ∉ newtonPolytope f.coeff.support := by
  obtain ⟨ℓ, hℓ⟩ := exists_separator_on_zeroFace F hF f hp hn
  obtain ⟨E, hE⟩ := exists_realForm_of_integral ℓ
  let bound := ∑ g ∈ f.coeff.support, |E (realExponent g)| / F (realExponent g)
  have hb : 0 ≤ bound := Finset.sum_nonneg fun g hg ↦ div_nonneg (abs_nonneg _) (hn g hg)
  apply zero_not_mem_newton_of_positive_form f.coeff.support ((bound + 1) • F + E)
  intro g hg
  change 0 < (bound + 1) * F (realExponent g) + E (realExponent g)
  rcases (hn g hg).eq_or_lt with he | ht
  · have hs : 0 < E (realExponent g) := by rw [hE]; exact_mod_cast hℓ g hg he.symm
    simpa [← he] using hs
  · have hs : |E (realExponent g)| / F (realExponent g) ≤ bound :=
      Finset.single_le_sum (fun a ha ↦ div_nonneg (abs_nonneg (E (realExponent a))) (hn a ha)) hg
    have hb' : |E (realExponent g)| ≤ bound * F (realExponent g) := (div_le_iff₀ ht).mp hs
    have hp' : 0 < F (realExponent g) := ht
    have ha : - |E (realExponent g)| ≤ E (realExponent g) := neg_abs_le _
    nlinarith

private theorem exists_supporting_form_of_not_interior (s : Set Plane) (hs : Convex ℝ s)
    (h0 : 0 ∈ s) (hi : 0 ∉ interior s) :
    ∃ F : Plane →ₗ[ℝ] ℝ, F ≠ 0 ∧ ∀ x ∈ s, 0 ≤ F x := by
  by_cases hInt : (interior s).Nonempty
  · obtain ⟨F, hF, hFs⟩ := geometric_hahn_banach_of_nonempty_interior_point hs hi hInt
    refine ⟨(-F).toLinearMap, ?_, fun x hx ↦ ?_⟩
    · intro he
      apply hF
      ext x
      have h := LinearMap.congr_fun he x
      simpa using h
    · have h := hFs x hx
      simpa using h
  · let S := Submodule.span ℝ s
    have hspan : S ≠ ⊤ := by
      intro h
      apply hInt
      apply hs.interior_nonempty_iff_affineSpan_eq_top.mpr
      apply SetLike.coe_injective
      have h' : ((affineSpan ℝ s : AffineSubspace ℝ Plane) : Set Plane) = S := by
        simpa only [Set.insert_eq_of_mem h0, S] using affineSpan_insert_zero (k := ℝ) s
      rw [h', h]
      rfl
    obtain ⟨v, hv⟩ := SetLike.exists_of_lt (lt_top_iff_ne_top.mpr hspan)
    obtain ⟨F, u, hu, hS⟩ := geometric_hahn_banach_point_closed S.convex S.closed_of_finiteDimensional hv.2
    have hu0 : u < 0 := by simpa using hS (0 : Plane) S.zero_mem
    refine ⟨F.toLinearMap, ?_, fun x hx ↦ ?_⟩
    · intro hF
      have hz : F v = 0 := LinearMap.congr_fun hF v
      linarith
    · have hxS : x ∈ S := Submodule.subset_span hx
      by_cases hxF : F x = 0
      · simp [hxF]
      · exfalso
        have hh := hS (((u - 1) / F x) • x) (S.smul_mem _ hxS)
        rw [map_smul, smul_eq_mul, div_mul_cancel₀ _ hxF] at hh
        linarith

/-- The rank-two Laurent theorem at the boundary, including Newton
polygons of smaller affine dimension. -/
theorem planarLaurent_pure_avoids_newton_of_not_interior (f : Laurent₂)
    (hp : Pure laurentCT f) (hi : 0 ∉ interior (newtonPolytope f.coeff.support)) :
    0 ∉ newtonPolytope f.coeff.support := by
  intro hz
  obtain ⟨F, hF, h⟩ := exists_supporting_form_of_not_interior _
    (convex_convexHull ℝ _) hz hi
  apply planarLaurent_pure_avoids_newton_of_supporting_form F hF f hp _ hz
  intro g hg
  exact h _ (subset_convexHull ℝ _ ⟨g, hg, rfl⟩)

/-- The remaining planar Laurent necessity statement concerns phases
whose Newton polygon contains zero in its interior. -/
def PlanarInteriorLaurentNecessity : Prop :=
  ∀ f : Laurent₂, 0 ∈ interior (newtonPolytope f.coeff.support) → ¬ Pure laurentCT f

/-- The rank-one theorem and the planar supporting-face argument reduce
rank-two Laurent necessity to the interior case. -/
theorem newtonLaurentNecessity_fin_two_of_interior (inside : PlanarInteriorLaurentNecessity) :
    NewtonLaurentNecessity ℂ (Fin 2) := by
  intro f _ hp
  by_cases hi : 0 ∈ interior (newtonPolytope f.coeff.support)
  · exact (inside f hi hp).elim
  · exact planarLaurent_pure_avoids_newton_of_not_interior f hp hi

end MathieuQuadratic
