/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.Bezout
import FiniteEtaleKeller.Reconstruction

/-!
# Reconstruction inside the inverse-polynomial quotient

A Bézout identity makes the derivative class an explicit unit of `K[S]/(E)`.
This module feeds that unit directly into the universal reconstruction formulas,
so no localization is used anywhere in the construction of the source point.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [CommRing K]

/-- The derivative unit in the inverse-polynomial quotient. -/
def quotientDerivativeUnit (E U V : K[X])
    (hbez : U * E + V * E.derivative = 1) : (AdjoinRoot E)ˣ :=
  adjoinRootDerivativeUnit E U V hbez

@[simp]
theorem quotientDerivativeUnit_val (E U V : K[X])
    (hbez : U * E + V * E.derivative = 1) :
    (quotientDerivativeUnit E U V hbez : AdjoinRoot E) =
      AdjoinRoot.mk E E.derivative := rfl

@[simp]
theorem quotientDerivativeUnit_inv_val (E U V : K[X])
    (hbez : U * E + V * E.derivative = 1) :
    (↑((quotientDerivativeUnit E U V hbez)⁻¹) : AdjoinRoot E) =
      AdjoinRoot.mk E V := rfl

/-- The quotient-ring reconstruction theorem.  Its hypothesis is precisely the
marked-line identity `D = E'/g₁` after the harmless scalar normalization of the
derivative unit.  The conclusion constructs `t,x,y,q,z` as actual quotient-ring
elements and proves all source-chart identities simultaneously. -/
theorem quotientReconstruction_identities
    (E U V : K[X]) (hbez : U * E + V * E.derivative = 1)
    (S Q pi a : AdjoinRoot E)
    (hD : (quotientDerivativeUnit E U V hbez : AdjoinRoot E) =
      1 - S * Q + pi * S ^ 2) :
    let d := quotientDerivativeUnit E U V hbez
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
    (quotientDerivativeUnit E U V hbez) hD

/-- Specialization to the canonical root class `S mod E`. -/
theorem quotientRootReconstruction_identities
    (E U V : K[X]) (hbez : U * E + V * E.derivative = 1)
    (Q pi a : AdjoinRoot E)
    (hD : (quotientDerivativeUnit E U V hbez : AdjoinRoot E) =
      1 - AdjoinRoot.root E * Q + pi * AdjoinRoot.root E ^ 2) :
    let d := quotientDerivativeUnit E U V hbez
    let S := AdjoinRoot.root E
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
  exact quotientReconstruction_identities E U V hbez
    (AdjoinRoot.root E) Q pi a hD

end FiniteEtaleKeller
