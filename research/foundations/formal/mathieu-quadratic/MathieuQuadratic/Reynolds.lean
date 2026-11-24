import MathieuQuadratic.Counterexample
import Mathlib.LinearAlgebra.Projection

/-!
# Projection from a Reynolds decomposition

This is the linear algebra of the geometric input in Section 2 of the
uniform construction. The direct sum of constants and an augmentation
subspace determines a unique normalized functional that kills that
subspace.
-/

namespace MathieuQuadratic

variable {K A : Type*} [Field K] [CommRing A] [Algebra K A]

/-- The normalized projection to the scalar coordinate. -/
noncomputable def reynoldsProjection (augmentation : Submodule K A)
    (hinj : Function.Injective (algebraMap K A))
    (h : IsCompl (LinearMap.range (Algebra.linearMap K A)) augmentation) : A →ₗ[K] K :=
  LinearMap.linearProjOfIsCompl augmentation (Algebra.linearMap K A) hinj h

theorem reynoldsProjection_normalized (augmentation : Submodule K A)
    (hinj : Function.Injective (algebraMap K A))
    (h : IsCompl (LinearMap.range (Algebra.linearMap K A)) augmentation) :
    reynoldsProjection augmentation hinj h 1 = 1 := by
  rw [← (algebraMap K A).map_one, reynoldsProjection]
  exact LinearMap.linearProjOfIsCompl_apply_left augmentation
    (Algebra.linearMap K A) hinj h 1

theorem reynoldsProjection_zero (augmentation : Submodule K A)
    (hinj : Function.Injective (algebraMap K A))
    (h : IsCompl (LinearMap.range (Algebra.linearMap K A)) augmentation) {f : A}
    (hf : f ∈ augmentation) :
    reynoldsProjection augmentation hinj h f = 0 :=
  LinearMap.linearProjOfIsCompl_apply_right' augmentation
    (Algebra.linearMap K A) hinj h f hf

theorem reynoldsProjection_unique (augmentation : Submodule K A)
    (hinj : Function.Injective (algebraMap K A))
    (h : IsCompl (LinearMap.range (Algebra.linearMap K A)) augmentation)
    (trace : A →ₗ[K] K) (hn : trace 1 = 1) (haug : ∀ f ∈ augmentation, trace f = 0) :
    trace = reynoldsProjection augmentation hinj h := by
  apply LinearMap.ext
  intro f
  have hsplit : f ∈ LinearMap.range (Algebra.linearMap K A) ⊔ augmentation :=
    h.sup_eq_top.symm ▸ Submodule.mem_top
  obtain ⟨g, ⟨c, rfl⟩, w, hw, rfl⟩ := Submodule.mem_sup.mp hsplit
  have hscalar : trace ((Algebra.linearMap K A) c) = c := by
    change trace (algebraMap K A c) = c
    rw [← mul_one (algebraMap K A c), ← Algebra.smul_def, map_smul, hn, smul_eq_mul, mul_one]
  rw [map_add, map_add, hscalar, haug w hw, reynoldsProjection, add_zero,
    LinearMap.linearProjOfIsCompl_apply_left, LinearMap.linearProjOfIsCompl_apply_right' _ _ _ _ _ hw,
    add_zero]

variable [CharZero K] {M : ℕ} (a : Fin 2 → K) (b : Fin M → K)

/-- The linear-algebraic Reynolds projection has the trace properties
used in the uniform construction once the augmentation subspace contains
the scalar differences and the rotation images. -/
theorem complementTrace_of_reynoldsDecomposition
    (augmentation : Submodule K (NormComplement a b))
    (hinj : Function.Injective (algebraMap K (NormComplement a b)))
    (h : IsCompl (LinearMap.range (Algebra.linearMap K (NormComplement a b))) augmentation)
    (scalar : ∀ f, scalarPullback a b (Units.mk0 (2 : K) two_ne_zero) f - f ∈ augmentation)
    (rotations : ∀ (i : Fin 2) (j : Fin M) (p : QuadPolynomial K M) (n : ℕ),
      algebraMap (QuadPolynomial K M) (NormComplement a b) (rotation a b i j p) *
        (↑((normUnit a b)⁻¹) : NormComplement a b) ^ n ∈ augmentation) :
    ComplementTrace a b (reynoldsProjection augmentation hinj h) := by
  refine ⟨reynoldsProjection_normalized augmentation hinj h, ?_, ?_⟩
  · intro f
    have hs := reynoldsProjection_zero augmentation hinj h (scalar f)
    rwa [map_sub, sub_eq_zero] at hs
  · intro i j p n
    exact reynoldsProjection_zero augmentation hinj h (rotations i j p n)

end MathieuQuadratic
