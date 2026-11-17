/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.RealizationFiber

/-!
# Automatic finite-étale realization

The translated represented-fiber theorem previously required an explicitly
supplied admissible parameter `a`.  This module combines it with the existence
theorem in `Admissibility.lean`: for every squarefree polynomial of degree at
least three over a characteristic-zero field, Lean now chooses an admissible
translation and constructs the natural represented-fiber equivalence without
any remaining parameter hypothesis.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K] [CharZero K]

/-- The quadratic-gauge datum attached to a squarefree polynomial, using an
automatically chosen admissible translation parameter. -/
def automaticRealizationDatum (P : K[X]) (hP : Squarefree P)
    (hdeg : 3 ≤ P.natDegree) : QuadraticGaugeFiberDatum K :=
  realizationDatum P (chosenAdmissibleTranslation P hdeg) hP
    (chosenAdmissibleTranslation_linear_ne_zero P hdeg)
    (chosenAdmissibleTranslation_cubic_ne_zero P hdeg)

variable {A B : Type*}
variable [CommRing A] [Algebra K A]
variable [CommRing B] [Algebra K B]

/-- Every squarefree polynomial of degree at least three represents the full
abstract quadratic-gauge source fiber over every commutative test algebra. -/
def automaticFiberRepresentingEquiv (P : K[X]) (hP : Squarefree P)
    (hdeg : 3 ≤ P.natDegree) :
    (AdjoinRoot P →ₐ[K] A) ≃ (automaticRealizationDatum P hP hdeg).Point A := by
  unfold automaticRealizationDatum
  exact translatedFiberRepresentingEquiv (A := A) P
    (chosenAdmissibleTranslation P hdeg) hP
    (chosenAdmissibleTranslation_linear_ne_zero P hdeg)
    (chosenAdmissibleTranslation_cubic_ne_zero P hdeg)

/-- The automatic represented-fiber equivalence is natural under every
morphism of commutative test algebras. -/
theorem automaticFiberRepresentingEquiv_natural
    (P : K[X]) (hP : Squarefree P) (hdeg : 3 ≤ P.natDegree)
    (f : A →ₐ[K] B) (φ : AdjoinRoot P →ₐ[K] A) :
    GaugeFiberPoint.map f (automaticFiberRepresentingEquiv (A := A) P hP hdeg φ) =
      automaticFiberRepresentingEquiv (A := B) P hP hdeg (f.comp φ) := by
  exact translatedFiberRepresentingEquiv_natural P
    (chosenAdmissibleTranslation P hdeg) hP
    (chosenAdmissibleTranslation_linear_ne_zero P hdeg)
    (chosenAdmissibleTranslation_cubic_ne_zero P hdeg) f φ

#print axioms automaticFiberRepresentingEquiv_natural

end FiniteEtaleKeller
