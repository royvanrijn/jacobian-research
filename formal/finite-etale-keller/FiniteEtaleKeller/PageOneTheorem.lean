/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeFunctionFieldComparison
import FiniteEtaleKeller.GeneralGaugeRealizationDegree
import FiniteEtaleKeller.GeneralGaugeFiberRank

/-!
# The page-one realization theorem

This module collects the determinant, geometric degree, represented literal
fiber, finite étaleness, rank, naturality, and effective degree bound in one
public certificate.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K] [CharZero K]

/-- The geometric degree of the actual determinant-one realization map.
The diagonal output normalization is a target automorphism, so its function
field degree is the degree computed by the comparison theorem for the
underlying gauge coordinates. -/
def automaticRealizationGeometricDegree
    (P : K[X]) (hdeg : 3 ≤ P.natDegree) : ℕ :=
  generalGaugeGeometricDegree
    (realizationSeed P (chosenAdmissibleTranslation P hdeg))
    (by
      simpa using
        chosenAdmissibleTranslation_linear_ne_zero P hdeg)
    (by
      simpa using
        chosenAdmissibleTranslation_cubic_ne_zero P hdeg)

/-- The determinant-one realization has geometric degree exactly
`P.natDegree`. -/
theorem automaticRealizationGeometricDegree_eq
    (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    automaticRealizationGeometricDegree P hdeg = P.natDegree := by
  let a := chosenAdmissibleTranslation P hdeg
  have h₁ : (realizationSeed P a).coeff 1 ≠ 0 := by
    simpa [a] using chosenAdmissibleTranslation_linear_ne_zero P hdeg
  have h₃ : (realizationSeed P a).coeff 3 ≠ 0 := by
    simpa [a] using chosenAdmissibleTranslation_cubic_ne_zero P hdeg
  have hseed :
      (realizationSeed P a).natDegree = P.natDegree :=
    realizationSeed_natDegree P a (by omega)
  have hseeddeg : 3 ≤ (realizationSeed P a).natDegree := by
    rwa [hseed]
  change generalGaugeGeometricDegree (realizationSeed P a) h₁ h₃ =
    P.natDegree
  rw [generalGaugeGeometricDegree_eq
    (realizationSeed P a) h₁ h₃ hseeddeg, hseed]

/-- The automatically chosen seed has nonzero linear coefficient. -/
theorem automaticRealizationSeed_linear_ne_zero
    (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    (realizationSeed P (chosenAdmissibleTranslation P hdeg)).coeff 1 ≠ 0 := by
  simpa using chosenAdmissibleTranslation_linear_ne_zero P hdeg

/-- The automatically chosen seed has nonzero cubic Hasse coefficient. -/
theorem automaticRealizationSeed_cubic_ne_zero
    (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    (realizationSeed P (chosenAdmissibleTranslation P hdeg)).coeff 3 ≠ 0 := by
  simpa using chosenAdmissibleTranslation_cubic_ne_zero P hdeg

/-- Translation and removal of the constant term preserve the degree of the
automatically chosen seed. -/
theorem automaticRealizationSeed_natDegree
    (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    (realizationSeed P (chosenAdmissibleTranslation P hdeg)).natDegree =
      P.natDegree :=
  realizationSeed_natDegree P (chosenAdmissibleTranslation P hdeg) (by omega)

/-- The canonical target presentation `K(Π,B)(C)` acts on the standard
three-variable source fraction field through the coordinate pullback of the
automatically chosen gauge map. -/
@[instance_reducible]
def automaticRealizationTargetFunctionFieldAlgebra
    (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    Algebra
      (RatFunc (FractionRing (MvPolynomial (Fin 2) K)))
      (FractionRing (MvPolynomial (Fin 3) K)) :=
  (generalGaugeFullyGenericTargetHom
    (realizationSeed P (chosenAdmissibleTranslation P hdeg))
    (automaticRealizationSeed_linear_ne_zero P hdeg)
    (automaticRealizationSeed_cubic_ne_zero P hdeg)).toRingHom.toAlgebra

/-- Explicit standard-object comparison between the inverse-root extension
and the source function field for the automatically chosen realization. -/
def automaticRealizationFunctionFieldComparison
    (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    letI : Algebra
        (RatFunc (FractionRing (MvPolynomial (Fin 2) K)))
        (FractionRing (MvPolynomial (Fin 3) K)) :=
      automaticRealizationTargetFunctionFieldAlgebra P hdeg
    AdjoinRoot
        (generalGaugeFullyGenericInversePolynomial
          (realizationSeed P (chosenAdmissibleTranslation P hdeg))) ≃ₐ[
      RatFunc (FractionRing (MvPolynomial (Fin 2) K))]
        FractionRing (MvPolynomial (Fin 3) K) := by
  letI : Algebra
      (RatFunc (FractionRing (MvPolynomial (Fin 2) K)))
      (FractionRing (MvPolynomial (Fin 3) K)) :=
    automaticRealizationTargetFunctionFieldAlgebra P hdeg
  exact generalGaugeFunctionFieldComparison
    (realizationSeed P (chosenAdmissibleTranslation P hdeg))
    (automaticRealizationSeed_linear_ne_zero P hdeg)
    (automaticRealizationSeed_cubic_ne_zero P hdeg)
    (by rw [automaticRealizationSeed_natDegree P hdeg]; exact hdeg)

/-- Direct function-field degree theorem in Mathlib's standard
`Module.finrank` and fraction-field objects.  The target algebra structure is
the pullback by the three displayed gauge coordinates; the determinant-one
output normalization is inverted by
`automaticRealizationMap_targetDenormalization`. -/
theorem automaticRealizationFunctionField_finrank
    (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    letI : Algebra
        (RatFunc (FractionRing (MvPolynomial (Fin 2) K)))
        (FractionRing (MvPolynomial (Fin 3) K)) :=
      automaticRealizationTargetFunctionFieldAlgebra P hdeg
    Module.finrank
      (RatFunc (FractionRing (MvPolynomial (Fin 2) K)))
      (FractionRing (MvPolynomial (Fin 3) K)) = P.natDegree := by
  letI : Algebra
      (RatFunc (FractionRing (MvPolynomial (Fin 2) K)))
      (FractionRing (MvPolynomial (Fin 3) K)) :=
    automaticRealizationTargetFunctionFieldAlgebra P hdeg
  rw [← (automaticRealizationFunctionFieldComparison P hdeg).toLinearEquiv.finrank_eq]
  rw [generalGaugeFullyGenericInverseAdjoinRoot_finrank
    (realizationSeed P (chosenAdmissibleTranslation P hdeg))
    (by rw [automaticRealizationSeed_natDegree P hdeg]; exact hdeg)]
  exact automaticRealizationSeed_natDegree P hdeg

/-- The determinant-one output normalization is explicitly inverted for
every seed by rescaling its second target coordinate by `-2`. -/
theorem generalGaugeJacobianOneMap_targetDenormalization
    (G : K[X]) :
    scaleOutput 1 (-2 : K) 1 (generalGaugeJacobianOneMap G) =
      generalGaugeMap G := by
  funext i
  fin_cases i <;>
    simp [generalGaugeJacobianOneMap, scaleOutput] <;>
    rw [← mul_assoc, ← MvPolynomial.C_mul] <;>
    norm_num

/-- The determinant-one output normalization is explicitly inverted by
rescaling its second target coordinate by `-2`. -/
theorem automaticRealizationMap_targetDenormalization
    (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    scaleOutput 1 (-2 : K) 1 (automaticRealizationMap P hdeg) =
      generalGaugeMap
        (realizationSeed P (chosenAdmissibleTranslation P hdeg)) := by
  exact generalGaugeJacobianOneMap_targetDenormalization
    (realizationSeed P (chosenAdmissibleTranslation P hdeg))

/-- A single record containing every assertion in the paper's page-one
polynomial-presentation theorem. -/
structure AutomaticPageOneCertificate
    (P : K[X]) (hP : Squarefree P) (hdeg : 3 ≤ P.natDegree) : Prop where
  targetNormalization :
    automaticRealizationMap P hdeg =
      scaleOutput 1 (-1 / 2 : K) 1
        (generalGaugeMap
          (realizationSeed P (chosenAdmissibleTranslation P hdeg)))
  targetDenormalization :
    scaleOutput 1 (-2 : K) 1 (automaticRealizationMap P hdeg) =
      generalGaugeMap
        (realizationSeed P (chosenAdmissibleTranslation P hdeg))
  jacobian :
    jacobianDet (automaticRealizationMap P hdeg) = 1
  geometricDegree :
    automaticRealizationGeometricDegree P hdeg = P.natDegree
  functionFieldFinrank :
    letI : Algebra
        (RatFunc (FractionRing (MvPolynomial (Fin 2) K)))
        (FractionRing (MvPolynomial (Fin 3) K)) :=
      automaticRealizationTargetFunctionFieldAlgebra P hdeg
    Module.finrank
      (RatFunc (FractionRing (MvPolynomial (Fin 2) K)))
      (FractionRing (MvPolynomial (Fin 3) K)) = P.natDegree
  fiber :
    ∀ (A : Type*) [CommRing A] [Algebra K A],
      Nonempty
        ((AdjoinRoot P →ₐ[K] A) ≃
          GeneralGaugeJacobianOneFiberPoint
            (realizationSeed P (chosenAdmissibleTranslation P hdeg)) 1
            (automaticRealizationTargetC P hdeg) A)
  fiber_natural :
    ∀ {A B : Type*} [CommRing A] [Algebra K A]
      [CommRing B] [Algebra K B]
      (f : A →ₐ[K] B) (φ : AdjoinRoot P →ₐ[K] A),
      GeneralGaugeJacobianOneFiberPoint.map f
          (automaticJacobianOneFiberRepresentingEquiv
            (A := A) P hP hdeg φ) =
        automaticJacobianOneFiberRepresentingEquiv
          (A := B) P hP hdeg (f.comp φ)
  finiteEtale : Algebra.Etale K (AdjoinRoot P)
  finite : Module.Finite K (AdjoinRoot P)
  rank : Module.finrank K (AdjoinRoot P) = P.natDegree
  degreeBound :
    ∀ i : Fin 3,
      (automaticRealizationMap P hdeg i).totalDegree ≤
        6 * P.natDegree + 2

/-- The final theorem: the actual displayed determinant-one map has geometric
degree `N`, its distinguished fiber is represented naturally by
`K[T]/(P)`, that quotient is finite étale of rank `N`, and every coordinate
has degree at most `6N+2`. -/
theorem automaticRealization_pageOne
    (P : K[X]) (hP : Squarefree P) (hdeg : 3 ≤ P.natDegree) :
    AutomaticPageOneCertificate P hP hdeg where
  targetNormalization := rfl
  targetDenormalization :=
    automaticRealizationMap_targetDenormalization P hdeg
  jacobian := automaticRealizationMap_jacobianDet P hdeg
  geometricDegree := automaticRealizationGeometricDegree_eq P hdeg
  functionFieldFinrank :=
    automaticRealizationFunctionField_finrank P hdeg
  fiber := fun A _ _ =>
    ⟨automaticJacobianOneFiberRepresentingEquiv (A := A) P hP hdeg⟩
  fiber_natural := fun f φ =>
    automaticJacobianOneFiberRepresentingEquiv_natural P hP hdeg f φ
  finiteEtale := automaticRepresentingAlgebra_etale P hP
  finite := automaticRepresentingAlgebra_finite P hP
  rank := automaticRealizationFiber_rank P hdeg
  degreeBound := automaticRealizationMap_totalDegree P hdeg

#print axioms generalGaugeFunctionFieldComparison
#print axioms generalGaugeGeometricDegree_eq
#print axioms automaticRealizationGeometricDegree_eq
#print axioms automaticRealizationFunctionFieldComparison
#print axioms automaticRealizationFunctionField_finrank
#print axioms generalGaugeJacobianOneMap_targetDenormalization
#print axioms automaticRealizationMap_targetDenormalization
#print axioms automaticRealization_pageOne

end FiniteEtaleKeller
