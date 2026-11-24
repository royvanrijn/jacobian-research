import MathieuQuadratic.CoefficientField
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure
import Mathlib.RingTheory.TensorProduct.MonoidAlgebra

/-!
# Scalar extension of a split Laurent trace

An algebraic splitting field need not embed in the extension field.
It embeds in its algebraic closure. The Laurent constant-term model
then extends to that common field, and the Mathieu property descends to
the extension field.
-/

open scoped TensorProduct
open Algebra.TensorProduct

namespace MathieuQuadratic

section Tower

variable {K S E A : Type*} [Field K] [Field S] [Field E] [CommRing A]
  [Algebra K S] [Algebra S E] [Algebra K E] [IsScalarTower K S E] [Algebra K A]

theorem extendTrace_cancelBaseChange (trace : A →ₗ[K] K) (v : E ⊗[S] (S ⊗[K] A)) :
    extendTrace E trace (cancelBaseChange K S E E A v) =
      extendTrace E (extendTrace S trace) v := by
  induction v using TensorProduct.induction_on with
  | zero => simp
  | add x y hx hy => rw [map_add, map_add, map_add, hx, hy]
  | tmul e s =>
    induction s using TensorProduct.induction_on with
    | zero => simp
    | add x y hx hy => rw [TensorProduct.tmul_add, map_add, map_add, map_add, hx, hy]
    | tmul s a =>
      simp only [cancelBaseChange_tmul, extendTrace_tmul, Algebra.smul_def, map_mul]
      rw [IsScalarTower.algebraMap_apply K S E]
      ring

end Tower

section Laurent

variable {K S E J A : Type*} [Field K] [Field S] [Field E] [AddCommGroup J] [CommRing A]
  [Algebra K S] [Algebra S E] [Algebra K E] [IsScalarTower K S E] [Algebra K A]

private theorem scalarTensor_laurentCT (v : E ⊗[S] AddMonoidAlgebra S J) :
    laurentCT (AddMonoidAlgebra.scalarTensorEquiv S E v) =
      extendTrace E (laurentCT (K := S) (G := J)) v := by
  induction v using TensorProduct.induction_on with
  | zero => simp
  | add x y hx hy => rw [map_add, map_add, map_add, hx, hy]
  | tmul e p =>
    simp only [AddMonoidAlgebra.scalarTensorEquiv_tmul, extendTrace_tmul]
    change (e • AddMonoidAlgebra.mapAlgHom J (Algebra.ofId S E) p).coeff 0 =
      algebraMap S E (p.coeff 0) * e
    rw [AddMonoidAlgebra.coeff_smul_apply, AddMonoidAlgebra.coeff_mapAlgHom]
    exact mul_comm _ _

/-- Scalar extension carries the constant-term description of a split
trace to the constant-term description over the larger field. -/
theorem split_laurent_model_after_base_change (trace : A →ₗ[K] K)
    (coordinates : S ⊗[K] A ≃ₐ[S] AddMonoidAlgebra S J)
    (identifies : ∀ f, extendTrace S trace f = laurentCT (coordinates f)) :
    ∃ coordinates' : E ⊗[K] A ≃ₐ[E] AddMonoidAlgebra E J,
      ∀ f, extendTrace E trace f = laurentCT (coordinates' f) := by
  let c := cancelBaseChange K S E E A
  let e := Algebra.TensorProduct.congr (AlgEquiv.refl : E ≃ₐ[E] E) coordinates
  refine ⟨c.symm.trans (e.trans (AddMonoidAlgebra.scalarTensorEquiv S E)), ?_⟩
  intro f
  obtain ⟨v, rfl⟩ := c.surjective f
  simp only [AlgEquiv.trans_apply, AlgEquiv.symm_apply_apply]
  rw [extendTrace_cancelBaseChange, scalarTensor_laurentCT]
  induction v using TensorProduct.induction_on with
  | zero => simp
  | add x y hx hy => rw [map_add, map_add, map_add, hx, hy]
  | tmul r x =>
    simp [e, Algebra.TensorProduct.map_tmul, identifies]

end Laurent

/-- Mathieu ascends for a trace with an algebraic split Laurent model.
This uses the complex Laurent necessity theorem for the exponent group. -/
theorem isMathieu_after_extension_of_split_laurent
    {K S L A J : Type*} [Field K] [CharZero K] [Field S] [Algebra K S]
    [Algebra.IsAlgebraic K S] [Field L] [Algebra K L] [CommRing A] [Algebra K A]
    [AddCommGroup J] (trace : A →ₗ[K] K)
    (coordinates : S ⊗[K] A ≃ₐ[S] AddMonoidAlgebra S J)
    (identifies : ∀ f, extendTrace S trace f = laurentCT (coordinates f))
    (necessity : LaurentSeparation ℂ J) :
    IsMathieu (extendTrace L trace) := by
  let E := AlgebraicClosure L
  let φ : S →ₐ[K] E := IsAlgClosed.lift
  letI : Algebra S E := φ.toRingHom.toAlgebra
  letI : IsScalarTower K S E := IsScalarTower.of_algebraMap_eq (R := K) (S := S) (A := E) fun x ↦ (φ.commutes x).symm
  letI : CharZero E := charZero_of_injective_algebraMap (algebraMap K E).injective
  obtain ⟨model, hmodel⟩ := split_laurent_model_after_base_change (E := E) trace coordinates identifies
  let c := cancelBaseChange K L E E A
  exact isMathieu_of_split_laurent (extendTrace L trace) (c.trans model)
    (fun f ↦ (extendTrace_cancelBaseChange trace f).symm.trans (hmodel (c f)))
    (laurentSeparation_of_complex necessity)

end MathieuQuadratic
