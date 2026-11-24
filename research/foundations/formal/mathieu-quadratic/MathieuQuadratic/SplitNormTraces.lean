import MathieuQuadratic.NormOneChart
import Mathlib.LinearAlgebra.Dimension.Constructions

/-!
# Comparisons with other norm traces

The split quadratic-étale model has a vector and a covector.
A nonzero pairing splits the space into the selected line and the
kernel of the covector.
-/

namespace MathieuQuadratic

variable {K V : Type*} [Field K] [AddCommGroup V] [Module K V]

private noncomputable def vectorCovectorCoordinates (x : V) (y : V →ₗ[K] K) (hy : y x ≠ 0) :
    V ≃ₗ[K] K × LinearMap.ker y := by
  let c : V →ₗ[K] K := (y x)⁻¹ • y
  have hc (v : V) : c v * y x = y v := by
    change ((y x)⁻¹ * y v) * y x = _
    rw [mul_comm ((y x)⁻¹), mul_assoc, inv_mul_cancel₀ hy, mul_one]
  let line : K →ₗ[K] V := (LinearMap.id : K →ₗ[K] K).smulRight x
  let p : V →ₗ[K] V := LinearMap.id - line.comp c
  have hp (v : V) : p v ∈ LinearMap.ker y := by
    rw [LinearMap.mem_ker]
    change y (v - c v • x) = 0
    rw [map_sub, map_smul, smul_eq_mul, hc, sub_self]
  let p' : V →ₗ[K] LinearMap.ker y := p.codRestrict _ hp
  let out : V →ₗ[K] K × LinearMap.ker y := c.prod p'
  let back : K × LinearMap.ker y →ₗ[K] V :=
    line.coprod (LinearMap.ker y).subtype
  refine LinearEquiv.ofLinear out back ?_ ?_
  · apply LinearMap.ext
    rintro ⟨t, v⟩
    have hv : y v = 0 := v.2
    have ht : c (t • x + (v : V)) = t := by
      change (y x)⁻¹ * y (t • x + (v : V)) = _
      rw [map_add, map_smul, smul_eq_mul, hv, add_zero, mul_left_comm,
        inv_mul_cancel₀ hy, mul_one]
    apply Prod.ext
    · exact ht
    · apply Subtype.ext
      change t • x + (v : V) - c (t • x + (v : V)) • x = v
      rw [ht]
      abel
  · apply LinearMap.ext
    intro v
    change c v • x + (v - c v • x) = v
    abel

/-- Transitivity of the split general linear action on pairs with a
fixed nonzero vector-covector pairing. -/
theorem vector_covector_same_norm [FiniteDimensional K V]
    (x x' : V) (y y' : V →ₗ[K] K) (h : y x = y' x') (hn : y x ≠ 0) :
    ∃ e : V ≃ₗ[K] V, e x = x' ∧ ∀ v, y' (e v) = y v := by
  have hn' : y' x' ≠ 0 := h ▸ hn
  let e₁ := vectorCovectorCoordinates x y hn
  let e₂ := vectorCovectorCoordinates x' y' hn'
  have hdim : Module.finrank K (LinearMap.ker y) = Module.finrank K (LinearMap.ker y') := by
    have h₁ := LinearEquiv.finrank_eq e₁
    have h₂ := LinearEquiv.finrank_eq e₂
    rw [Module.finrank_prod, Module.finrank_self] at h₁ h₂
    omega
  let k := LinearEquiv.ofFinrankEq _ _ hdim
  let e := e₁.trans ((LinearEquiv.refl K K).prodCongr k) |>.trans e₂.symm
  have hfix₁ : e₁ x = (1, 0) := by
    have hx : (y x)⁻¹ * y x = 1 := inv_mul_cancel₀ hn
    apply Prod.ext
    · exact hx
    · apply Subtype.ext
      change x - ((y x)⁻¹ * y x) • x = 0
      rw [hx, one_smul, sub_self]
  have hfix₂ : e₂ x' = (1, 0) := by
    have hx : (y' x')⁻¹ * y' x' = 1 := inv_mul_cancel₀ hn'
    apply Prod.ext
    · exact hx
    · apply Subtype.ext
      change x' - ((y' x')⁻¹ * y' x') • x' = 0
      rw [hx, one_smul, sub_self]
  refine ⟨e, ?_, ?_⟩
  · change e₂.symm ((LinearEquiv.refl K K).prodCongr k (e₁ x)) = _
    rw [hfix₁, LinearEquiv.prodCongr_apply, LinearEquiv.refl_apply, map_zero,
      ← hfix₂, LinearEquiv.symm_apply_apply]
  · intro v
    have hker : y' (↑(k (e₁ v).2) : V) = 0 := (k (e₁ v).2).2
    change y' ((((y x)⁻¹ * y v) : K) • x' + ↑(k (e₁ v).2)) = y v
    rw [map_add, map_smul, smul_eq_mul, hker, add_zero, ← h, mul_comm ((y x)⁻¹), mul_assoc,
      inv_mul_cancel₀ hn, mul_one]

/-- A trace that detects every nonzero square gives a Mathieu–Zhao kernel. -/
theorem isMathieu_of_square_detection {A : Type*} [CommRing A] [Algebra K A]
    (trace : A →ₗ[K] K) (detect : ∀ f : A, f ≠ 0 → trace (f ^ 2) ≠ 0) : IsMathieu trace := by
  intro f hp
  have hf : f = 0 := by
    by_contra h
    exact detect f h (hp 2 (by decide))
  subst f
  intro g
  exact ⟨1, fun m hm ↦ by rw [zero_pow (by omega : m ≠ 0), mul_zero, map_zero]⟩

end MathieuQuadratic
