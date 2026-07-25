/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeRealization
import Mathlib.RingTheory.AdjoinRoot

/-!
# Rank of the represented literal fiber

The final literal fiber is represented by `AdjoinRoot P = K[T]/(P)`.  This
module records explicitly that its `K`-dimension is the degree of `P`, so the
special-fiber length used by the paper is part of the formal certificate.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

/-- The quotient algebra representing the literal fiber has rank exactly the
polynomial degree. -/
theorem adjoinRoot_finrank_eq_natDegree (P : K[X]) :
    Module.finrank K (AdjoinRoot P) = P.natDegree := by
  change Module.finrank K (K[X] ⧸ Ideal.span {P}) = P.natDegree
  exact finrank_quotient_span_eq_natDegree

/-- In the range relevant to the realization theorem, the represented fiber
algebra is nonzero. -/
theorem automaticRepresentingAlgebra_nontrivial
    [CharZero K] (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    Nontrivial (AdjoinRoot P) := by
  apply AdjoinRoot.nontrivial
  rw [Polynomial.degree_eq_natDegree]
  · exact_mod_cast (show P.natDegree ≠ 0 by omega)
  · intro h
    simp [h] at hdeg

/-- The automatically realized literal fiber has length `P.natDegree`. -/
theorem automaticRealizationFiber_rank
    [CharZero K] (P : K[X]) (hdeg : 3 ≤ P.natDegree) :
    Module.finrank K (AdjoinRoot P) = P.natDegree :=
  adjoinRoot_finrank_eq_natDegree P

#print axioms adjoinRoot_finrank_eq_natDegree
#print axioms automaticRepresentingAlgebra_nontrivial
#print axioms automaticRealizationFiber_rank

end FiniteEtaleKeller
