import Mathlib.Algebra.MvPolynomial.Eval
import Mathlib.RingTheory.Localization.Away.Basic
import Mathlib.RingTheory.TensorProduct.Maps

/-!
# Scalar extension of the principal open coordinate ring
-/


open scoped TensorProduct

namespace MathieuQuadratic

open MvPolynomial Algebra.TensorProduct

variable {K L J : Type*} [Field K] [Field L] [Algebra K L]
variable (p : MvPolynomial J K)

local notation "R" => MvPolynomial J K
local notation "S" => MvPolynomial J L
local notation "A" => Localization.Away p

private noncomputable def coefficientChange : R →ₐ[K] S :=
  MvPolynomial.mapAlgHom (Algebra.ofId K L)

local notation "B" => Localization.Away (coefficientChange (L := L) p)

private noncomputable def toExtendedLocalization : A →ₐ[K] B :=
  IsLocalization.Away.liftAlgHom p
    (f := (IsScalarTower.toAlgHom K S B).comp coefficientChange)
    (IsLocalization.Away.algebraMap_isUnit (coefficientChange (L := L) p))

private theorem toExtendedLocalization_polynomial (q : R) :
    toExtendedLocalization (L := L) p (algebraMap R A q) = algebraMap S B (coefficientChange q) :=
  IsLocalization.Away.lift_eq
    (g := ((IsScalarTower.toAlgHom K S B).comp coefficientChange).toRingHom)
    p (IsLocalization.Away.algebraMap_isUnit (coefficientChange (L := L) p)) q

private noncomputable def polynomialToTensor : S →ₐ[L] L ⊗[K] A :=
  aeval fun i ↦ includeRight (algebraMap R A (X i))

private theorem polynomialToTensor_coefficientChange :
    ((polynomialToTensor (L := L) p).restrictScalars K).comp coefficientChange =
      includeRight.comp (IsScalarTower.toAlgHom K R A) := by
  apply MvPolynomial.algHom_ext
  intro i
  simp [polynomialToTensor, coefficientChange]

private theorem polynomialToTensor_on_base (q : R) :
    polynomialToTensor (L := L) p (coefficientChange q) = includeRight (algebraMap R A q) :=
  AlgHom.congr_fun (polynomialToTensor_coefficientChange p) q

private theorem polynomialToTensor_unit :
    IsUnit (polynomialToTensor (L := L) p (coefficientChange p)) := by
  rw [polynomialToTensor_on_base]
  have hu : IsUnit (algebraMap R A p) := IsLocalization.Away.algebraMap_isUnit p
  exact hu.map (includeRight : A →ₐ[K] L ⊗[K] A)

private noncomputable def extensionToTensor : B →ₐ[L] L ⊗[K] A :=
  IsLocalization.Away.liftAlgHom (coefficientChange p) (polynomialToTensor_unit (L := L) p)

private theorem extensionToTensor_apply (q : S) :
    extensionToTensor (L := L) p (algebraMap S B q) = polynomialToTensor p q :=
  IsLocalization.Away.lift_eq (coefficientChange p) (polynomialToTensor_unit p) q

private noncomputable def tensorToExtension : L ⊗[K] A →ₐ[L] B :=
  Algebra.TensorProduct.lift (Algebra.ofId L B) (toExtendedLocalization p) (fun _ _ ↦ Commute.all _ _)

private theorem tensorToExtension_includeRight (a : A) :
    tensorToExtension (L := L) p (includeRight a) = toExtendedLocalization p a := by
  simp [tensorToExtension, includeRight_apply]

private theorem tensorToExtension_extensionToTensor :
    (tensorToExtension (L := L) p).comp (extensionToTensor p) = AlgHom.id L B := by
  apply IsLocalization.algHom_ext (Submonoid.powers (coefficientChange (L := L) p))
  apply MvPolynomial.algHom_ext
  intro i
  change tensorToExtension (L := L) p (extensionToTensor (L := L) p (algebraMap S B (X i))) = algebraMap S B (X i)
  rw [extensionToTensor_apply, polynomialToTensor, aeval_X, tensorToExtension_includeRight,
    toExtendedLocalization_polynomial]
  simp [coefficientChange]

private theorem extensionToTensor_tensorToExtension :
    (extensionToTensor (L := L) p).comp (tensorToExtension p) = AlgHom.id L (L ⊗[K] A) := by
  apply Algebra.TensorProduct.ext
  · exact Subsingleton.elim _ _
  · apply IsLocalization.algHom_ext (Submonoid.powers p)
    apply AlgHom.ext
    intro q
    change extensionToTensor (L := L) p
        (tensorToExtension (L := L) p ((includeRight : A →ₐ[K] L ⊗[K] A) (algebraMap R A q))) =
      (includeRight : A →ₐ[K] L ⊗[K] A) (algebraMap R A q)
    rw [tensorToExtension_includeRight, toExtendedLocalization_polynomial,
      extensionToTensor_apply, polynomialToTensor_on_base]

/-- Extending coefficients of a principal open commutes with tensoring
its coordinate ring by the extension field. -/
noncomputable def localizationBaseChange :
    L ⊗[K] A ≃ₐ[L] Localization.Away (MvPolynomial.map (algebraMap K L) p) :=
  AlgEquiv.ofAlgHom (tensorToExtension p) (extensionToTensor p)
    (tensorToExtension_extensionToTensor p) (extensionToTensor_tensorToExtension p)

theorem localizationBaseChange_includeRight (q : R) :
    localizationBaseChange (L := L) p (includeRight (algebraMap R A q)) =
      algebraMap S (Localization.Away (MvPolynomial.map (algebraMap K L) p))
        (MvPolynomial.map (algebraMap K L) q) := by
  change tensorToExtension (L := L) p
    ((includeRight : A →ₐ[K] L ⊗[K] A) (algebraMap R A q)) = algebraMap S B (coefficientChange q)
  rw [tensorToExtension_includeRight, toExtendedLocalization_polynomial]

end MathieuQuadratic
