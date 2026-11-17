/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeMap
import FiniteEtaleKeller.GeneralGaugeFunctionFieldComparison
import FiniteEtaleKeller.PageOneTheorem
import FiniteEtaleKeller.AbstractFiniteEtale
import FiniteEtaleKeller.CollisionFiber

/-!
# Certificate accompanying the finite-étale Keller-fiber paper

This is the publication-sized Lean artifact.  Its direct imports are limited
to the general gauge construction, its function-field comparison, the
polynomial page-one theorem, the abstract finite-étale realization, and the
presentation-independent collision-fiber interface.  Arithmetic examples and
unrelated companion results deliberately remain outside this module.

The paper states the polynomial theorem using `Polynomial.Separable`; the
internal construction uses the equivalent `Squarefree` hypothesis.  The
first declaration below is the exact interface between those statements.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

universe u

variable {K : Type*} [Field K] [CharZero K]

/-- The paper's characteristic-zero polynomial-presentation theorem, with
exactly its separability and degree hypotheses.  The returned certificate
contains the determinant-one map, geometric and function-field degree,
literal represented fiber and its naturality, finite étaleness, rank, and
coordinate-degree bound. -/
theorem paperPolynomialPresentation_pageOne
    (P : K[X]) (hP : P.Separable) (hdeg : 3 ≤ P.natDegree) :
    AutomaticPageOneCertificate P hP.squarefree hdeg :=
  automaticRealization_pageOne P hP.squarefree hdeg

/-- The paper's abstract finite-étale corollary: every finite étale algebra
of rank at least three over a characteristic-zero field has the complete
page-one realization certificate. -/
noncomputable def paperAbstractFiniteEtale_pageOne
    (A : Type u) [CommRing A] [Algebra K A] [Algebra.Etale K A]
    (hrank : 3 ≤ Module.finrank K A) :
    AbstractFiniteEtalePageOneCertificate (K := K) A hrank :=
  abstractFiniteEtale_pageOne (K := K) A hrank

/-- Compile-time guard for the exact headline theorem signature. -/
example (P : K[X]) (hP : P.Separable) (hdeg : 3 ≤ P.natDegree) :
    AutomaticPageOneCertificate P hP.squarefree hdeg :=
  paperPolynomialPresentation_pageOne P hP hdeg

/-- Compile-time guard for the exact abstract finite-étale corollary. -/
example (A : Type u) [CommRing A] [Algebra K A] [Algebra.Etale K A]
    (hrank : 3 ≤ Module.finrank K A) :
    AbstractFiniteEtalePageOneCertificate (K := K) A hrank :=
  paperAbstractFiniteEtale_pageOne A hrank

#print axioms generalGaugeFunctionFieldComparison
#print axioms generalGaugeGeometricDegree_eq
#print axioms automaticRealizationFunctionFieldComparison
#print axioms automaticRealizationFunctionField_finrank
#print axioms automaticRealization_pageOne
#print axioms paperPolynomialPresentation_pageOne
#print axioms finiteEtalePresentation
#print axioms abstractFiniteEtale_pageOne
#print axioms paperAbstractFiniteEtale_pageOne
#print axioms fiberCollisionDiagonal_surjective
#print axioms fiberCollisionPointPairsEquiv
#print axioms fiberCollisionDiagonal_not_injective
#print axioms finrank_fiberCollisionObstruction

end FiniteEtaleKeller
