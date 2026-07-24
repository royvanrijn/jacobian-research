/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.QuotientReconstruction
import Mathlib.FieldTheory.Perfect

/-!
# Reconstruction from separability

The paper assumes that the inverse polynomial is squarefree.  Over a
characteristic-zero field this is equivalent to separability, hence to a
Bézout identity with the derivative.  This module chooses those Bézout
coefficients internally and exposes the derivative class as a canonical unit
for reconstruction.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

/-- The first chosen coefficient in a Bézout identity for a separable
polynomial and its derivative. -/
def separableBezoutU (E : K[X]) (hE : E.Separable) : K[X] :=
  Classical.choose ((Polynomial.separable_def' E).mp hE)

/-- The second chosen coefficient in a Bézout identity for a separable
polynomial and its derivative. -/
def separableBezoutV (E : K[X]) (hE : E.Separable) : K[X] :=
  Classical.choose (Classical.choose_spec ((Polynomial.separable_def' E).mp hE))

/-- The chosen coefficients satisfy the required Bézout identity. -/
theorem separableBezout_identity (E : K[X]) (hE : E.Separable) :
    separableBezoutU E hE * E
      + separableBezoutV E hE * E.derivative = 1 :=
  Classical.choose_spec
    (Classical.choose_spec ((Polynomial.separable_def' E).mp hE))

/-- The derivative class of a separable polynomial, packaged as a unit in its
root quotient. -/
def derivativeUnitOfSeparable (E : K[X]) (hE : E.Separable) :
    (AdjoinRoot E)ˣ :=
  adjoinRootDerivativeUnit E (separableBezoutU E hE)
    (separableBezoutV E hE) (separableBezout_identity E hE)

@[simp]
theorem derivativeUnitOfSeparable_val (E : K[X]) (hE : E.Separable) :
    (derivativeUnitOfSeparable E hE : AdjoinRoot E) =
      AdjoinRoot.mk E E.derivative := rfl

/-- Separability alone supplies the quotient-ring reconstruction; no explicit
Bézout coefficients or localization are required in the theorem statement. -/
theorem separableQuotientReconstruction_identities
    (E : K[X]) (hE : E.Separable)
    (S Q pi a : AdjoinRoot E)
    (hD : (derivativeUnitOfSeparable E hE : AdjoinRoot E) =
      1 - S * Q + pi * S ^ 2) :
    let d := derivativeUnitOfSeparable E hE
    let t : AdjoinRoot E := ↑d⁻¹
    let x := S * t
    let y := Q - pi * S
    let q := pi * (d : AdjoinRoot E)
    let z := (d : AdjoinRoot E) ^ 2 *
      (q - a * y ^ 2 * (1 + 3 * t))
    (1 + x * y = t)
      ∧ (t ^ 2 * z + a * y ^ 2 * (1 + 3 * t) = q)
      ∧ (x * (d : AdjoinRoot E) = S)
      ∧ (y + x * q = Q)
      ∧ (t * q = pi)
      ∧ (t * (d : AdjoinRoot E) = 1) := by
  exact unitReconstruction_identities S Q pi a
    (derivativeUnitOfSeparable E hE) hD

/-- Over a characteristic-zero field, squarefreeness is enough to construct
the derivative unit. -/
def derivativeUnitOfSquarefree [CharZero K]
    (E : K[X]) (hE : Squarefree E) : (AdjoinRoot E)ˣ :=
  derivativeUnitOfSeparable E
    ((PerfectField.separable_iff_squarefree).2 hE)

@[simp]
theorem derivativeUnitOfSquarefree_val [CharZero K]
    (E : K[X]) (hE : Squarefree E) :
    (derivativeUnitOfSquarefree E hE : AdjoinRoot E) =
      AdjoinRoot.mk E E.derivative := rfl

/-- The reconstruction theorem in exactly the squarefree form used by the
paper. -/
theorem squarefreeQuotientReconstruction_identities [CharZero K]
    (E : K[X]) (hE : Squarefree E)
    (S Q pi a : AdjoinRoot E)
    (hD : (derivativeUnitOfSquarefree E hE : AdjoinRoot E) =
      1 - S * Q + pi * S ^ 2) :
    let d := derivativeUnitOfSquarefree E hE
    let t : AdjoinRoot E := ↑d⁻¹
    let x := S * t
    let y := Q - pi * S
    let q := pi * (d : AdjoinRoot E)
    let z := (d : AdjoinRoot E) ^ 2 *
      (q - a * y ^ 2 * (1 + 3 * t))
    (1 + x * y = t)
      ∧ (t ^ 2 * z + a * y ^ 2 * (1 + 3 * t) = q)
      ∧ (x * (d : AdjoinRoot E) = S)
      ∧ (y + x * q = Q)
      ∧ (t * q = pi)
      ∧ (t * (d : AdjoinRoot E) = 1) := by
  exact unitReconstruction_identities S Q pi a
    (derivativeUnitOfSquarefree E hE) hD

end FiniteEtaleKeller
