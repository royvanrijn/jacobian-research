/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeRealization

/-!
# Base change for the quadratic-gauge construction

The paper's scalar-extension statement keeps the supplied polynomial and
translation parameter fixed while applying a homomorphism of ground fields.
This module proves that coefficientwise base change commutes with the rooted
translation, the realization target, every coordinate of the general gauge,
and its determinant-one output normalization.

The automatically chosen translation is deliberately not used here:
`chosenAdmissibleTranslation` is defined by classical choice and is not
expected to be natural under field homomorphisms.
-/

noncomputable section

open Polynomial
open MvPolynomial
open scoped TensorProduct

namespace FiniteEtaleKeller

variable {K L : Type*} [Field K] [Field L]

section PolynomialTranslation

variable (f : K →+* L)

/-- Translation of the polynomial variable commutes with coefficientwise base
change. -/
@[simp]
theorem translatePolynomial_map (P : K[X]) (a : K) :
    (translatePolynomial P a).map f =
      translatePolynomial (P.map f) (f a) := by
  simp [translatePolynomial, Polynomial.map_comp]

/-- Removing the translated constant term also commutes with base change. -/
@[simp]
theorem rootedTranslate_map (P : K[X]) (a : K) :
    (rootedTranslate P a).map f =
      rootedTranslate (P.map f) (f a) := by
  simp [rootedTranslate]

/-- The supplied-translation seed used by the actual realization map commutes
with base change. -/
@[simp]
theorem realizationSeed_map (P : K[X]) (a : K) :
    (realizationSeed P a).map f =
      realizationSeed (P.map f) (f a) := by
  simp [realizationSeed]

/-- Every coefficient of the supplied realization seed is carried to the
corresponding coefficient after base change. -/
theorem realizationSeed_coeff_map (P : K[X]) (a : K) (n : ℕ) :
    (realizationSeed (P.map f) (f a)).coeff n =
      f ((realizationSeed P a).coeff n) := by
  rw [← realizationSeed_map f]
  exact Polynomial.coeff_map _ _

/-- In particular, the two nonvanishing seed coefficients required by the
quadratic gauge remain nonzero after extending the ground field. -/
theorem realizationSeed_admissible_map (P : K[X]) (a : K)
    (h₁ : (realizationSeed P a).coeff 1 ≠ 0)
    (h₃ : (realizationSeed P a).coeff 3 ≠ 0) :
    (realizationSeed (P.map f) (f a)).coeff 1 ≠ 0 ∧
      (realizationSeed (P.map f) (f a)).coeff 3 ≠ 0 := by
  rw [realizationSeed_coeff_map f, realizationSeed_coeff_map f]
  constructor
  · simpa using f.injective.ne h₁
  · simpa using f.injective.ne h₃

/-- Squarefreeness of the prescribed polynomial is preserved by a field
embedding. -/
theorem squarefree_map_fieldHom [CharZero K] [CharZero L]
    (P : K[X]) (hP : Squarefree P) :
    Squarefree (P.map f) :=
  (PerfectField.separable_iff_squarefree).1
    ((PerfectField.separable_iff_squarefree).2 hP).map

/-- The distinguished third target coordinate commutes with a homomorphism of
ground fields. -/
@[simp]
theorem realizationTargetC_map (P : K[X]) (a g₁ : K) :
    f (realizationTargetC P a g₁) =
      realizationTargetC (P.map f) (f a) (f g₁) := by
  simp [realizationTargetC]
  simp only [map_ofNat]

end PolynomialTranslation

section AdjoinRootBaseChange

/-- Adjoining a root commutes with extension of the ground field.

The temporary `Algebra K L` instance in the result is the one induced by the
given field homomorphism `f`. -/
noncomputable def adjoinRootBaseChangeEquiv (f : K →+* L) (P : K[X]) :
    letI : Algebra K L := f.toAlgebra
    L ⊗[K] AdjoinRoot P ≃ₐ[L] AdjoinRoot (P.map f) := by
  letI : Algebra K L := f.toAlgebra
  exact
    (AdjoinRoot.tensorAlgEquiv P
        (P.map Algebra.TensorProduct.includeRight.toRingHom) rfl).trans
      (AdjoinRoot.mapAlgEquiv (Algebra.TensorProduct.rid K L L)
        (P.map Algebra.TensorProduct.includeRight.toRingHom)
        (P.map f) (by
          apply Associated.of_eq
          rw [Polynomial.map_map]
          congr 1
          ext x
          simp [Algebra.smul_def, RingHom.algebraMap_toAlgebra]))

end AdjoinRootBaseChange

section GeneralGauge

variable (f : K →+* L)

/-- The recurrent polynomial `t = 1 + xy` is unchanged by coefficient base
change. -/
@[simp]
theorem generalGaugeT_map :
    MvPolynomial.map f (generalGaugeT : GaugePolynomial K) =
      (generalGaugeT : GaugePolynomial L) := by
  simp [generalGaugeT]

/-- The recurrent polynomial `q` commutes with coefficient base change. -/
@[simp]
theorem generalGaugeQ_map (G : K[X]) :
    MvPolynomial.map f (generalGaugeQ G) =
      generalGaugeQ (G.map f) := by
  simp [generalGaugeQ, Polynomial.coeff_map]
  simp only [map_ofNat]
  simp

/-- The first gauge coordinate commutes with coefficient base change. -/
@[simp]
theorem generalGaugePi_map (G : K[X]) :
    MvPolynomial.map f (generalGaugePi G) =
      generalGaugePi (G.map f) := by
  simp [generalGaugePi]

/-- The complete second gauge coordinate commutes with coefficient base
change, including its arbitrary finite tail. -/
@[simp]
theorem generalGaugeB_map (G : K[X]) :
    MvPolynomial.map f (generalGaugeB G) =
      generalGaugeB (G.map f) := by
  simp [generalGaugeB, Polynomial.natDegree_map f, Polynomial.coeff_map]
  simp only [map_ofNat]

/-- The complete third gauge coordinate commutes with coefficient base
change, including its arbitrary finite tail. -/
@[simp]
theorem generalGaugeC_map (G : K[X]) :
    MvPolynomial.map f (generalGaugeC G) =
      generalGaugeC (G.map f) := by
  simp [generalGaugeC, Polynomial.natDegree_map f, Polynomial.coeff_map]
  simp only [map_ofNat]

/-- The full three-coordinate general gauge commutes with coefficient base
change. -/
theorem generalGaugeMap_map (G : K[X]) :
    (fun i => MvPolynomial.map f (generalGaugeMap G i)) =
      generalGaugeMap (G.map f) := by
  funext i
  fin_cases i <;> simp [generalGaugeMap]

/-- The determinant-one output normalization commutes with coefficient base
change. -/
theorem generalGaugeJacobianOneMap_map (G : K[X]) :
    (fun i => MvPolynomial.map f (generalGaugeJacobianOneMap G i)) =
      generalGaugeJacobianOneMap (G.map f) := by
  funext i
  fin_cases i
  · simp [generalGaugeJacobianOneMap, scaleOutput, generalGaugeMap]
  · simp [generalGaugeJacobianOneMap, scaleOutput, generalGaugeMap]
    simp only [map_ofNat]
    simp
  · simp [generalGaugeJacobianOneMap, scaleOutput, generalGaugeMap]

end GeneralGauge

section SuppliedRealization

variable (f : K →+* L)

/-- The actual determinant-one realization map attached to a supplied
translation parameter commutes with extension of the ground field. -/
theorem realizationJacobianOneMap_map (P : K[X]) (a : K) :
    (fun i =>
      MvPolynomial.map f
        (generalGaugeJacobianOneMap (realizationSeed P a) i)) =
      generalGaugeJacobianOneMap
        (realizationSeed (P.map f) (f a)) := by
  rw [generalGaugeJacobianOneMap_map f, realizationSeed_map f]

/-- The complete supplied map-target pair used by the realization theorem is
compatible with extension of the ground field. -/
theorem realizationMapTarget_map (P : K[X]) (a : K) :
    (fun i =>
      MvPolynomial.map f
        (generalGaugeJacobianOneMap (realizationSeed P a) i)) =
        generalGaugeJacobianOneMap
          (realizationSeed (P.map f) (f a))
    ∧
      f (realizationTargetC P a (P.derivative.eval a)) =
        realizationTargetC (P.map f) (f a)
          ((P.map f).derivative.eval (f a)) := by
  constructor
  · exact realizationJacobianOneMap_map f P a
  · rw [realizationTargetC_map]
    simp

end SuppliedRealization

#print axioms translatePolynomial_map
#print axioms squarefree_map_fieldHom
#print axioms adjoinRootBaseChangeEquiv
#print axioms generalGaugeMap_map
#print axioms generalGaugeJacobianOneMap_map
#print axioms realizationMapTarget_map

end FiniteEtaleKeller
