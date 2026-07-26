/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeRealizationDegree
import FiniteEtaleKeller.AnnouncedCounterexample

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

/-- The cubic seed must remain exactly, not merely linearly equivalent to, the
map in the original announcement. -/
example : generalGaugeMap announcedSeed = announcedCounterexampleMap :=
  generalGaugeMap_announcedSeed

/-- Signature guard for the final construction-level certificate. -/
example (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    jacobianDet (automaticRealizationMap P hdeg) = 1 ∧
      ∀ i : Fin 3,
        (automaticRealizationMap P hdeg i).totalDegree ≤
          6 * P.natDegree + 2 :=
  automaticRealizationMap_certificate P hdeg

#print axioms jacobianDet_generalGaugeMap
#print axioms jacobianDet_generalGaugeJacobianOneMap
#print axioms generalGaugeInversePolynomial_derivative
#print axioms generalGaugeMap_announcedSeed
#print axioms automaticRealizationMap_certificate
#print axioms automaticJacobianOneFiberRepresentingEquiv_natural

end FiniteEtaleKeller
