import MathieuQuadratic.BlockPolynomials
import MathieuQuadratic.GeometricAction
import MathieuQuadratic.Reynolds
import Mathlib.RepresentationTheory.Semisimple

/-!
# The normalized geometric trace from complete reducibility
-/

open Matrix
open scoped Matrix

namespace MathieuQuadratic

variable {K J : Type*} [Field K] [CharZero K] [Fintype J] [DecidableEq J] [Nonempty J]
  [IsAlgClosed K]
variable (d : J → Kˣ)

local notation "A" => Localization.Away (diagonalNormPolynomial d)

/-- A point of the product of the special orthogonal and scalar groups. -/
@[ext]
structure OrthogonalScalarPoint where
  matrix : Matrix J J K
  scalar : Kˣ
  determinant : matrix.det = 1
  preserves : ∀ v, diagonalNorm d (matrix *ᵥ v) = diagonalNorm d v

/-- Pullback reverses composition of the underlying point actions. -/
instance : Monoid (OrthogonalScalarPoint d) where
  one := ⟨1, 1, Matrix.det_one, fun v ↦ by rw [Matrix.one_mulVec]⟩
  mul g h := ⟨h.matrix * g.matrix, h.scalar * g.scalar, by
    rw [Matrix.det_mul, h.determinant, g.determinant, one_mul], fun v ↦ by
    rw [← Matrix.mulVec_mulVec, h.preserves, g.preserves]⟩
  one_mul g := by
    apply OrthogonalScalarPoint.ext
    · change g.matrix * 1 = g.matrix
      exact mul_one _
    · change g.scalar * 1 = g.scalar
      exact mul_one _
  mul_one g := by
    apply OrthogonalScalarPoint.ext
    · change 1 * g.matrix = g.matrix
      exact one_mul _
    · change 1 * g.scalar = g.scalar
      exact one_mul _
  mul_assoc g h i := by
    apply OrthogonalScalarPoint.ext
    · change i.matrix * (h.matrix * g.matrix) = (i.matrix * h.matrix) * g.matrix
      exact (mul_assoc _ _ _).symm
    · change i.scalar * (h.scalar * g.scalar) = (i.scalar * h.scalar) * g.scalar
      exact (mul_assoc _ _ _).symm

omit [CharZero K] [Nonempty J] [IsAlgClosed K] in
private theorem point_action_norm_ne_zero (g : OrthogonalScalarPoint d) (v : J → K)
    (hv : diagonalNorm d v ≠ 0) : diagonalNorm d (g.scalar • (g.matrix *ᵥ v)) ≠ 0 := by
  rw [Units.smul_def, diagonalNorm_smul, g.preserves]
  exact mul_ne_zero (pow_ne_zero _ g.scalar.ne_zero) hv

/-- The induced family of linear pullbacks on the coordinate ring. -/
noncomputable def orthogonalScalarAction (g : OrthogonalScalarPoint d) : A →ₗ[K] A :=
  (geometricPullback d g.matrix g.scalar g.preserves).toLinearMap

/-- The scalar-isometry point action is a linear representation. -/
noncomputable def orthogonalScalarRepresentation : Representation K (OrthogonalScalarPoint d) A where
  toFun := orthogonalScalarAction d
  map_one' := by
    ext f
    apply complement_funext d
    intro v hv
    let e : OrthogonalScalarPoint d := 1
    have hnorm : diagonalNorm d (e.scalar • (e.matrix *ᵥ v)) ≠ 0 :=
      point_action_norm_ne_zero d e v hv
    change evalOnComplement d v hv (geometricPullback d e.matrix e.scalar e.preserves f) =
      evalOnComplement d v hv f
    rw [eval_geometricPullback d e.matrix e.scalar e.preserves v hv hnorm]
    have he : e.scalar • (e.matrix *ᵥ v) = v := by
      change (1 : Kˣ) • ((1 : Matrix J J K) *ᵥ v) = v
      rw [Matrix.one_mulVec, one_smul]
    simp [he]
  map_mul' g h := by
    ext f
    apply complement_funext d
    intro v hv
    let vg := g.scalar • (g.matrix *ᵥ v)
    have hvg : diagonalNorm d vg ≠ 0 := point_action_norm_ne_zero d g v hv
    have hvgh := point_action_norm_ne_zero d (g * h) v hv
    have hvh := point_action_norm_ne_zero d h vg hvg
    change evalOnComplement d v hv (geometricPullback d (g * h).matrix (g * h).scalar (g * h).preserves f) =
      evalOnComplement d v hv (geometricPullback d g.matrix g.scalar g.preserves
        (geometricPullback d h.matrix h.scalar h.preserves f))
    rw [eval_geometricPullback d (g * h).matrix (g * h).scalar _ v hv hvgh,
      eval_geometricPullback d g.matrix g.scalar _ v hv hvg,
      eval_geometricPullback d h.matrix h.scalar _ vg hvg hvh]
    have he : (g * h).scalar • ((g * h).matrix *ᵥ v) = h.scalar • (h.matrix *ᵥ vg) := by
      change (h.scalar * g.scalar) • ((h.matrix * g.matrix) *ᵥ v) = _
      rw [← Matrix.mulVec_mulVec, mul_smul]
      unfold vg
      rw [Matrix.mulVec_smul]
    simp [he]

theorem orthogonalScalarAction_fixed :
    actionFixed (orthogonalScalarAction d) = LinearMap.range (Algebra.linearMap K A) := by
  ext f
  rw [mem_actionFixed_iff, LinearMap.mem_range]
  constructor
  · intro hf
    have hconst := (geometric_fixed_iff_constant d f).1
      (fun T c hdet hnorm ↦ hf ⟨T, c, hdet, hnorm⟩)
    obtain ⟨k, hk, _⟩ := hconst
    exact ⟨k, hk.symm⟩
  · rintro ⟨k, rfl⟩ g
    exact (geometricPullback d g.matrix g.scalar g.preserves).commutes k

private theorem geometric_algebraMap_injective : Function.Injective (algebraMap K A) := by
  intro x y hxy
  obtain ⟨_, _, h⟩ := (geometric_fixed_iff_constant d (algebraMap K A x)).1 (by
    intro T c hdet hT
    exact AlgHom.commutes _ _)
  exact (h x rfl).trans (h y hxy).symm

/-- The unique normalized trace invariant under the actual scalar and
special orthogonal pullbacks, given complete reducibility. -/
theorem existsUnique_geometric_point_trace
    (complete : ∀ W : Submodule K A, ActionStable (orthogonalScalarAction d) W →
      ∃ W', ActionStable (orthogonalScalarAction d) W' ∧ IsCompl W W') :
    ∃! trace : A →ₗ[K] K, trace 1 = 1 ∧ ∀ g f,
      trace (orthogonalScalarAction d g f) = trace f := by
  let D := blockDifferenceSpan (orthogonalScalarAction d)
  have hD : IsCompl (LinearMap.range (Algebra.linearMap K A)) D := by
    rw [← orthogonalScalarAction_fixed d]
    exact isCompl_actionFixed_differenceSpan _ complete
  let trace := reynoldsProjection D (geometric_algebraMap_injective d) hD
  refine ⟨trace, ⟨reynoldsProjection_normalized _ _ _, ?_⟩, ?_⟩
  · intro g f
    have hmem : orthogonalScalarAction d g f - f ∈ D :=
      Submodule.subset_span ⟨g, f, rfl⟩
    have hz : trace (orthogonalScalarAction d g f - f) = 0 :=
      reynoldsProjection_zero _ _ _ hmem
    rwa [map_sub, sub_eq_zero] at hz
  · intro F hF
    apply reynoldsProjection_unique D (geometric_algebraMap_injective d) hD F hF.1
    intro p hp
    induction hp using Submodule.span_induction with
    | mem x hx =>
      obtain ⟨g, f, rfl⟩ := hx
      rw [map_sub, hF.2, sub_self]
    | zero => simp
    | add x y _ _ hx hy => rw [map_add, hx, hy, add_zero]
    | smul c x _ hx => rw [map_smul, hx, smul_zero]

/-- The normalized geometric trace with the Mathlib semisimplicity
hypothesis on the actual coordinate-ring representation. -/
theorem existsUnique_geometric_point_trace_of_semisimple
    [Representation.IsSemisimpleRepresentation (orthogonalScalarRepresentation d)] :
    ∃! trace : A →ₗ[K] K, trace 1 = 1 ∧ ∀ g f,
      trace (orthogonalScalarAction d g f) = trace f := by
  apply existsUnique_geometric_point_trace
  intro W hW
  let U : Subrepresentation (orthogonalScalarRepresentation d) :=
    ⟨W, fun g v hv ↦ hW g v hv⟩
  obtain ⟨U', hc⟩ := exists_isCompl U
  refine ⟨U'.toSubmodule, fun g v hv ↦ U'.apply_mem_toSubmodule g hv, ?_⟩
  refine ⟨disjoint_iff.mpr ?_, codisjoint_iff.mpr ?_⟩
  · exact congrArg Subrepresentation.toSubmodule hc.inf_eq_bot
  · exact congrArg Subrepresentation.toSubmodule hc.sup_eq_top

end MathieuQuadratic
