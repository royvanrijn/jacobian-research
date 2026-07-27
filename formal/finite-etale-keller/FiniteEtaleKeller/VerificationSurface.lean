/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeRealizationDegree
import FiniteEtaleKeller.GeneralGaugeBaseChange
import FiniteEtaleKeller.AnnouncedCounterexample
import FiniteEtaleKeller.ExplicitFiber
import FiniteEtaleKeller.ExplicitThreeAdicPoint
import FiniteEtaleKeller.GeneralGaugeLocalizedFiber
import FiniteEtaleKeller.GenericInverseIrreducibility
import FiniteEtaleKeller.GeneralGaugeFunctionField
import FiniteEtaleKeller.GeneralGaugeFullGenericDegree
import FiniteEtaleKeller.PageOneTheorem

/-!
# Public verification surface

This module deliberately restates the load-bearing public certificates with
fully explicit types.  It is a compile-time regression guard: weakening or
silently changing a hypothesis, dropping naturality, or changing the actual
map/target represented by the final theorem must break this file.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K] [CharZero K]

/-- Signature guard for the raw determinant theorem.  Both nonvanishing
hypotheses are mathematically necessary. -/
example (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    jacobianDet (generalGaugeMap G) = MvPolynomial.C (-2) :=
  jacobianDet_generalGaugeMap G h₁ h₃

/-- Signature guard for the normalized determinant theorem. -/
example (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    jacobianDet (generalGaugeJacobianOneMap G) = 1 :=
  jacobianDet_generalGaugeJacobianOneMap G h₁ h₃

/-- Signature guard for the arbitrary-degree derivative factorization. -/
example (G : K[X]) (pi b c : K) (h₁ : G.coeff 1 ≠ 0) :
    (generalGaugeInversePolynomial G pi b c).derivative =
      C (G.coeff 1) * markedChartPolynomial pi b (generalGaugeBeta G pi) :=
  generalGaugeInversePolynomial_derivative G pi b c h₁

/-- For fixed `pi ≠ 0` and `b`, the one-parameter inverse equation over
`K(C)` is irreducible and has the full seed degree. -/
example (G : K[X]) (pi b : K)
    (h₁ : G.coeff 1 ≠ 0)
    (hdeg : 3 ≤ G.natDegree) (hpi : pi ≠ 0) :
    Irreducible (generalGaugeGenericInversePolynomial G pi b) ∧
      (generalGaugeGenericInversePolynomial G pi b).natDegree =
        G.natDegree :=
  generalGaugeGenericInversePolynomial_certificate
    G pi b h₁ hdeg hpi

/-- The quotient by the same fixed-parameter inverse equation has
`K(C)`-dimension `N`; together with the preceding irreducibility certificate,
this is the corresponding degree-`N` field extension. -/
example (G : K[X]) (pi b : K)
    (hdeg : 3 ≤ G.natDegree) (hpi : pi ≠ 0) :
    Module.finrank (RatFunc K)
        (AdjoinRoot (generalGaugeGenericInversePolynomial G pi b)) =
      G.natDegree :=
  generalGaugeGenericInverseAdjoinRoot_finrank
    G pi b hdeg hpi

/-- With `Pi` and `B` promoted to independent parameters, the full inverse
equation over the iterated presentation `K(Pi,B)(C)` is irreducible and has
the seed degree. -/
example (G : K[X]) (h₁ : G.coeff 1 ≠ 0)
    (hdeg : 3 ≤ G.natDegree) :
    Irreducible (generalGaugeFullyGenericInversePolynomial G) ∧
      (generalGaugeFullyGenericInversePolynomial G).natDegree =
        G.natDegree :=
  generalGaugeFullyGenericInversePolynomial_certificate G h₁ hdeg

/-- Its root quotient has the expected dimension over `K(Pi,B)(C)`. -/
example (G : K[X]) (hdeg : 3 ≤ G.natDegree) :
    Module.finrank (RatFunc (GaugeTargetParameterField K))
        (AdjoinRoot (generalGaugeFullyGenericInversePolynomial G)) =
      G.natDegree :=
  generalGaugeFullyGenericInverseAdjoinRoot_finrank G hdeg

/-- The displayed gauge coordinates are algebraically independent. -/
example (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    AlgebraicIndependent K (generalGaugeMap G) :=
  generalGaugeMap_algebraicIndependent G h₁ h₃

/-- The coordinate map induces an injective pullback on rational function
fields. -/
example (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0) :
    Function.Injective (generalGaugeFunctionFieldHom G h₁ h₃) :=
  generalGaugeFunctionFieldHom_injective G h₁ h₃

/-- The missing bridge: over the target field `K(Π,B)(C)`, the actual source
function field is explicitly the inverse-root extension. -/
example (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (hdeg : 3 ≤ G.natDegree) :
    letI : Algebra (RatFunc (GaugeTargetParameterField K))
        (GaugeFunctionField K) :=
      (generalGaugeFullyGenericTargetHom G h₁ h₃).toRingHom.toAlgebra
    GaugeFunctionField K ≃ₐ[RatFunc (GaugeTargetParameterField K)]
      AdjoinRoot (generalGaugeFullyGenericInversePolynomial G) :=
  generalGaugeSourceFunctionFieldComparison G h₁ h₃ hdeg

/-- Consequently the actual displayed map has geometric degree exactly the
seed degree. -/
example (G : K[X]) (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (hdeg : 3 ≤ G.natDegree) :
    generalGaugeGeometricDegree G h₁ h₃ = G.natDegree :=
  generalGaugeGeometricDegree_eq G h₁ h₃ hdeg

/-- The complete supplied-translation realization map and target commute
with coefficient extension.  The automatically chosen translation is
intentionally excluded because classical choice is not functorial. -/
example {L : Type*} [Field L] (f : K →+* L) (P : K[X]) (a : K) :
    (fun i =>
      MvPolynomial.map f
        (generalGaugeJacobianOneMap (realizationSeed P a) i)) =
        generalGaugeJacobianOneMap
          (realizationSeed (P.map f) (f a))
    ∧
      f (realizationTargetC P a (P.derivative.eval a)) =
        realizationTargetC (P.map f) (f a)
          ((P.map f).derivative.eval (f a)) :=
  realizationMapTarget_map f P a

/-- The localization away from the derivative represents roots at which the
derivative is invertible, without any separability hypothesis. -/
example (E : K[X]) (A : Type*) [CommRing A] [Algebra K A] :
    (LocalizedAdjoinRoot E →ₐ[K] A) ≃ LocalizedPolynomialRoot E A :=
  LocalizedPolynomialRoot.localizedAlgHomEquiv E A

/-- Full signature guard for the paper's localized-fiber theorem. -/
example (G : K[X]) (pi : Kˣ) (b c : K)
    (h₁ : G.coeff 1 ≠ 0) (h₃ : G.coeff 3 ≠ 0)
    (A : Type*) [CommRing A] [Algebra K A] :
    (LocalizedAdjoinRoot
        (generalGaugeInversePolynomial G (pi : K) b c) →ₐ[K] A) ≃
      GeneralGaugeRawFiberPoint G pi b c A :=
  generalGaugeLocalizedRawRepresentingEquiv G pi b c h₁ h₃ A

/-- The cubic seed must remain exactly, not merely linearly equivalent to, the
map in the original announcement. -/
example : generalGaugeMap announcedSeed = announcedCounterexampleMap :=
  generalGaugeMap_announcedSeed

/-- The represented algebra has exactly the polynomial degree claimed for the
special-fiber length. -/
example (P : K[X]) :
    Module.finrank K (AdjoinRoot P) = P.natDegree :=
  adjoinRoot_finrank_eq_natDegree P

/-- A squarefree represented quotient is finite étale over the base field. -/
example (P : K[X]) (hP : Squarefree P) :
    Algebra.Etale K (AdjoinRoot P) ∧ Module.Finite K (AdjoinRoot P) :=
  ⟨automaticRepresentingAlgebra_etale P hP,
    automaticRepresentingAlgebra_finite P hP⟩

/-- The exact denominator-free quintic map displayed in the paper has the
literal fiber represented by the explicit Berend--Bilu quotient. -/
example (A : Type*) [CommRing A] [Algebra ℚ A] :
    (AdjoinRoot ExplicitQuintic.p5 →ₐ[ℚ] A) ≃
      ExplicitQuintic.IntegralFiberPoint A :=
  ExplicitQuintic.integralFiberRepresentingEquiv

/-- The quotient representing the exact quintic fiber is finite étale. -/
example : Algebra.Etale ℚ (AdjoinRoot ExplicitQuintic.p5) :=
  ExplicitQuintic.p5_quotient_etale

/-- The same quotient is finite over `ℚ`. -/
example : Module.Finite ℚ (AdjoinRoot ExplicitQuintic.p5) :=
  ExplicitQuintic.p5_quotient_finite

/-- The exact displayed quintic fiber has no rational point. -/
example : IsEmpty (ExplicitQuintic.IntegralFiberPoint ℚ) :=
  ExplicitQuintic.integralFiberPoint_rat_isEmpty

/-- Its archimedean local fiber is nonempty. -/
example : Nonempty (ExplicitQuintic.IntegralFiberPoint ℝ) :=
  ExplicitQuintic.integralFiberPoint_real_nonempty

/-- A first nonarchimedean case is fully formal: the exact displayed fiber
has a point over the three-adic field. -/
example : Nonempty (ExplicitQuintic.IntegralFiberPoint ℚ_[3]) :=
  ExplicitQuintic.integralFiberPoint_threeAdic_nonempty

/-- Signature guard for the final construction-level certificate. -/
example (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    jacobianDet (automaticRealizationMap P hdeg) = 1 ∧
      ∀ i : Fin 3,
        (automaticRealizationMap P hdeg i).totalDegree ≤
          6 * P.natDegree + 2 :=
  automaticRealizationMap_certificate P hdeg

/-- Signature guard for the single page-one theorem. -/
example (P : K[X]) (hP : Squarefree P) (hdeg : 3 ≤ P.natDegree) :
    AutomaticPageOneCertificate P hP hdeg :=
  automaticRealization_pageOne P hP hdeg

#print axioms jacobianDet_generalGaugeMap
#print axioms jacobianDet_generalGaugeJacobianOneMap
#print axioms generalGaugeInversePolynomial_derivative
#print axioms generalGaugeGenericInversePolynomial_certificate
#print axioms generalGaugeGenericInverseAdjoinRoot_finrank
#print axioms generalGaugeFullyGenericInversePolynomial_certificate
#print axioms generalGaugeFullyGenericInverseAdjoinRoot_finrank
#print axioms generalGaugeMap_algebraicIndependent
#print axioms generalGaugeFunctionFieldHom_injective
#print axioms generalGaugeSourceFunctionFieldComparison
#print axioms generalGaugeGeometricDegree_eq
#print axioms generalGaugeJacobianOneMap_targetDenormalization
#print axioms automaticRealizationMap_targetDenormalization
#print axioms realizationMapTarget_map
#print axioms adjoinRootBaseChangeEquiv
#print axioms LocalizedPolynomialRoot.localizedAlgHomEquiv
#print axioms generalGaugeLocalizedRawRepresentingEquiv
#print axioms generalGaugeLocalizedRawRepresentingEquiv_natural
#print axioms generalGaugeLocalizedJacobianOneRepresentingEquiv_natural
#print axioms generalGaugeMap_announcedSeed
#print axioms adjoinRoot_finrank_eq_natDegree
#print axioms automaticRepresentingAlgebra_etale
#print axioms automaticRepresentingAlgebra_finite
#print axioms ExplicitQuintic.integralFiberRepresentingEquiv_natural
#print axioms ExplicitQuintic.p5_quotient_etale
#print axioms ExplicitQuintic.p5_quotient_finite
#print axioms ExplicitQuintic.integralFiberPoint_rat_isEmpty
#print axioms ExplicitQuintic.integralFiberPoint_real_nonempty
#print axioms ExplicitQuintic.integralFiberPoint_threeAdic_nonempty
#print axioms automaticRealizationMap_certificate
#print axioms automaticJacobianOneFiberRepresentingEquiv_natural
#print axioms automaticRealization_pageOne

end FiniteEtaleKeller
