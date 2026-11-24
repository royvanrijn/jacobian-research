import MathieuQuadratic.CoefficientField
import Mathlib.Analysis.Convex.Topology
import Mathlib.Analysis.LocallyConvex.Separation

/-!
# Integral separation of a Laurent Newton polytope

Strictly separate the origin from the finite convex hull over the reals.
Scaling and rounding the finitely many coordinate coefficients gives an
integral covector, with a uniform error bound on the support.
-/

open scoped BigOperators

namespace MathieuQuadratic

def realExponent {J : Type*} (g : J → ℤ) : J → ℝ := fun i ↦ g i

/-- The real convex hull of a finite Laurent support. -/
def newtonPolytope {J : Type*} (s : Finset (J → ℤ)) : Set (J → ℝ) :=
  convexHull ℝ (realExponent '' (s : Set (J → ℤ)))

theorem exists_integral_separator {J : Type*} [Fintype J] (s : Finset (J → ℤ))
    (hz : 0 ∉ newtonPolytope s) :
    ∃ ℓ : (J → ℤ) →+ ℤ, ∀ g ∈ s, 0 < ℓ g := by
  classical
  have hfin : Set.Finite (realExponent '' (s : Set (J → ℤ))) := s.finite_toSet.image _
  obtain ⟨f, u, hu, hf⟩ := geometric_hahn_banach_point_closed
    (convex_convexHull ℝ _) ((hfin.isCompact_convexHull ℝ).isClosed) hz
  rw [map_zero] at hu
  let r (i : J) : ℝ := f fun j ↦ if i = j then 1 else 0
  have hsum (g : J → ℤ) : f (realExponent g) = ∑ i, r i * (g i : ℝ) := by
    have hs := f.toLinearMap.pi_apply_eq_sum_univ (realExponent g)
    change f (realExponent g) = ∑ i, realExponent g i • f (fun j ↦ if i = j then 1 else 0) at hs
    simpa only [r, realExponent, smul_eq_mul, mul_comm] using hs
  let B : ℝ := ∑ g ∈ s, ∑ i, |(g i : ℝ)|
  have hB : 0 ≤ B := Finset.sum_nonneg fun _ _ ↦ Finset.sum_nonneg fun _ _ ↦ abs_nonneg _
  obtain ⟨N, hN⟩ := exists_nat_gt (B / u)
  have hNpos : 0 < (N : ℝ) := (div_nonneg hB hu.le).trans_lt hN
  let m (i : J) : ℤ := ⌊(N : ℝ) * r i⌋
  let ℓ : (J → ℤ) →+ ℤ :=
    { toFun := fun g ↦ ∑ i, m i * g i
      map_zero' := by simp
      map_add' := fun g h ↦ by simp [Pi.add_apply, mul_add, Finset.sum_add_distrib] }
  refine ⟨ℓ, fun g hg ↦ ?_⟩
  have hg' := hf (realExponent g) (subset_convexHull ℝ _ ⟨g, hg, rfl⟩)
  have hNg : B < (N : ℝ) * f (realExponent g) :=
    ((div_lt_iff₀ hu).mp hN).trans (mul_lt_mul_of_pos_left hg' hNpos)
  have hdiff (i : J) : |(N : ℝ) * r i - (m i : ℝ)| ≤ 1 := by
    rw [abs_of_nonneg (sub_nonneg.mpr (Int.floor_le _))]
    change (N : ℝ) * r i - (↑⌊(N : ℝ) * r i⌋ : ℝ) ≤ 1
    linarith [Int.lt_floor_add_one ((N : ℝ) * r i)]
  have herr : (N : ℝ) * f (realExponent g) - (ℓ g : ℝ) ≤ B := by
    calc
      _ = ∑ i, ((N : ℝ) * r i - (m i : ℝ)) * (g i : ℝ) := by
        change (N : ℝ) * f (realExponent g) - (↑(∑ i, m i * g i) : ℝ) = _
        rw [hsum, Finset.mul_sum, Int.cast_sum, ← Finset.sum_sub_distrib]
        refine Finset.sum_congr rfl fun i _ ↦ ?_
        rw [Int.cast_mul]
        ring
      _ ≤ ∑ i, |(g i : ℝ)| := Finset.sum_le_sum fun i _ ↦ by
        apply (le_abs_self _).trans
        rw [abs_mul]
        exact mul_le_of_le_one_left (abs_nonneg _) (hdiff i)
      _ ≤ B := Finset.single_le_sum (f := fun x : J → ℤ ↦ ∑ i, |(x i : ℝ)|)
        (fun _ _ ↦ Finset.sum_nonneg fun _ _ ↦ abs_nonneg _) hg
  exact_mod_cast (show 0 < (ℓ g : ℝ) by linarith)

theorem exists_realForm_of_integral {J : Type*} [Fintype J] (ℓ : (J → ℤ) →+ ℤ) :
    ∃ F : (J → ℝ) →ₗ[ℝ] ℝ, ∀ g : J → ℤ, F (realExponent g) = (ℓ g : ℝ) := by
  classical
  let r (i : J) : ℤ := ℓ fun j ↦ if i = j then 1 else 0
  let F : (J → ℝ) →ₗ[ℝ] ℝ := ∑ i, (r i : ℝ) • LinearMap.proj i
  refine ⟨F, fun g ↦ ?_⟩
  have hi := ℓ.toIntLinearMap.pi_apply_eq_sum_univ g
  change ℓ g = ∑ i, g i • ℓ (fun j ↦ if i = j then 1 else 0) at hi
  simp only [zsmul_eq_mul, Int.cast_id] at hi
  simp only [F, LinearMap.sum_apply, LinearMap.smul_apply, LinearMap.proj_apply,
    realExponent, smul_eq_mul]
  rw [hi, Int.cast_sum]
  simp only [Int.cast_mul, r, mul_comm]

theorem origin_not_mem_newton_of_positive_support {J : Type*} [Fintype J]
    (s : Finset (J → ℤ)) (ℓ : (J → ℤ) →+ ℤ) (h : ∀ g ∈ s, 0 < ℓ g) :
    0 ∉ newtonPolytope s := by
  obtain ⟨F, hmap⟩ := exists_realForm_of_integral ℓ
  have hh : newtonPolytope s ⊆ {x | 0 < F x} := by
    apply convexHull_min _ (convex_halfSpace_gt F.isLinear 0)
    rintro _ ⟨g, hg, rfl⟩
    change 0 < F (realExponent g)
    rw [hmap]
    exact_mod_cast h g hg
  intro hzero
  have hc : 0 < F 0 := hh hzero
  simp at hc

/-- The exact Newton-polytope form of Laurent necessity. -/
def NewtonLaurentNecessity (K J : Type*) [Field K] : Prop :=
  ∀ f : AddMonoidAlgebra K (J → ℤ), f ≠ 0 → Pure laurentCT f → 0 ∉ newtonPolytope f.coeff.support

theorem laurentSeparation_of_newton {K J : Type*} [Field K] [Fintype J]
    (necessity : NewtonLaurentNecessity K J) : LaurentSeparation K (J → ℤ) := by
  intro f hp
  by_cases hz : f = 0
  · exact ⟨0, by simp [hz]⟩
  · exact exists_integral_separator f.coeff.support (necessity f hz hp)

theorem newtonLaurentNecessity_of_complex {K J : Type*} [Field K] [CharZero K] [Fintype J]
    (necessity : NewtonLaurentNecessity ℂ J) : NewtonLaurentNecessity K J := by
  intro f _ hp
  obtain ⟨ℓ, hℓ⟩ := laurentSeparation_of_complex (laurentSeparation_of_newton necessity) f hp
  exact origin_not_mem_newton_of_positive_support f.coeff.support ℓ hℓ

theorem pure_iff_origin_not_mem_newton {K J : Type*} [Field K] [CharZero K] [Fintype J]
    (necessity : NewtonLaurentNecessity ℂ J) (f : AddMonoidAlgebra K (J → ℤ)) :
    Pure laurentCT f ↔ 0 ∉ newtonPolytope f.coeff.support := by
  constructor
  · intro hf
    obtain ⟨ℓ, hℓ⟩ := laurentSeparation_of_complex (laurentSeparation_of_newton necessity) f hf
    exact origin_not_mem_newton_of_positive_support f.coeff.support ℓ hℓ
  · intro hz
    obtain ⟨ℓ, hℓ⟩ := exists_integral_separator f.coeff.support hz
    exact pure_of_positive_support f ℓ hℓ

end MathieuQuadratic
