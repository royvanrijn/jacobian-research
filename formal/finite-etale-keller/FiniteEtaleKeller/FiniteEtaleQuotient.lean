/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.RingTheory.Etale.StandardEtale

/-!
# Separable polynomial quotients are finite étale

This module packages a separable polynomial quotient as a Mathlib standard
étale presentation.  It handles nonmonic polynomials by unit normalization.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

/-- The standard étale presentation attached to a monic separable
polynomial.  The localization polynomial is `1`; separability supplies the
Bézout identity required by `StandardEtalePair`. -/
def standardEtalePairOfMonicSeparable
    (E : K[X]) (hmonic : E.Monic) (hsep : E.Separable) :
    StandardEtalePair K where
  f := E
  monic_f := hmonic
  g := 1
  cond := by
    obtain ⟨U, V, hbez⟩ := (Polynomial.separable_def' E).mp hsep
    refine ⟨V, U, 0, ?_⟩
    calc
      E.derivative * V + E * U =
          U * E + V * E.derivative := by ring
      _ = 1 := hbez
      _ = (1 : K[X]) ^ 0 := by simp

/-- The standard étale presentation at the unit localization is canonically
the usual polynomial quotient. -/
def standardEtalePairOfMonicSeparableEquiv
    (E : K[X]) (hmonic : E.Monic) (hsep : E.Separable) :
    (standardEtalePairOfMonicSeparable E hmonic hsep).Ring ≃ₐ[K]
      AdjoinRoot E := by
  let P := standardEtalePairOfMonicSeparable E hmonic hsep
  let eAway := P.equivAwayAdjoinRoot
  let eUnit :
      AdjoinRoot E ≃ₐ[AdjoinRoot E]
        Localization.Away (AdjoinRoot.mk E (1 : K[X])) :=
    IsLocalization.atUnit
      (AdjoinRoot E)
      (Localization.Away (AdjoinRoot.mk E (1 : K[X])))
      (AdjoinRoot.mk E (1 : K[X]))
      (by simp)
  exact eAway.trans (eUnit.symm.restrictScalars K)

/-- A monic separable polynomial presents an étale quotient algebra. -/
theorem adjoinRoot_etale_of_monic_separable
    (E : K[X]) (hmonic : E.Monic) (hsep : E.Separable) :
    Algebra.Etale K (AdjoinRoot E) := by
  let P := standardEtalePairOfMonicSeparable E hmonic hsep
  haveI : Algebra.Etale K P.Ring := inferInstance
  exact Algebra.Etale.of_equiv
    (standardEtalePairOfMonicSeparableEquiv E hmonic hsep)

/-- Over a field, every separable polynomial presents an étale quotient
algebra.  Monic normalization changes the polynomial only by a unit, hence
does not change its `AdjoinRoot`. -/
theorem adjoinRoot_etale_of_separable
    (E : K[X]) (hsep : E.Separable) :
    Algebra.Etale K (AdjoinRoot E) := by
  let Emonic := E * Polynomial.C E.leadingCoeff⁻¹
  have hE0 : E ≠ 0 := hsep.ne_zero
  have hunit : IsUnit (Polynomial.C E.leadingCoeff⁻¹) := by
    exact Polynomial.isUnit_C.mpr
      (isUnit_iff_ne_zero.mpr (inv_ne_zero (Polynomial.leadingCoeff_ne_zero.mpr hE0)))
  have hassoc : Associated E Emonic :=
    associated_mul_unit_right E (Polynomial.C E.leadingCoeff⁻¹) hunit
  have hmonic : Emonic.Monic :=
    Polynomial.monic_mul_leadingCoeff_inv hE0
  have hsepMonic : Emonic.Separable :=
    hassoc.separable hsep
  letI : Algebra.Etale K (AdjoinRoot Emonic) :=
    adjoinRoot_etale_of_monic_separable Emonic hmonic hsepMonic
  exact Algebra.Etale.of_equiv
    (AdjoinRoot.algEquivOfAssociated K Emonic E hassoc.symm)

theorem adjoinRoot_finite_of_separable
    (E : K[X]) (hsep : E.Separable) :
    Module.Finite K (AdjoinRoot E) := by
  let Emonic := E * Polynomial.C E.leadingCoeff⁻¹
  have hE0 : E ≠ 0 := hsep.ne_zero
  have hunit : IsUnit (Polynomial.C E.leadingCoeff⁻¹) := by
    exact Polynomial.isUnit_C.mpr
      (isUnit_iff_ne_zero.mpr (inv_ne_zero (Polynomial.leadingCoeff_ne_zero.mpr hE0)))
  have hassoc : Associated E Emonic :=
    associated_mul_unit_right E (Polynomial.C E.leadingCoeff⁻¹) hunit
  haveI : Module.Finite K (AdjoinRoot Emonic) :=
    (Polynomial.monic_mul_leadingCoeff_inv hE0).finite_adjoinRoot
  exact Module.Finite.equiv
    (AdjoinRoot.algEquivOfAssociated K Emonic E hassoc.symm).toLinearEquiv

#print axioms adjoinRoot_etale_of_separable
#print axioms adjoinRoot_finite_of_separable

end FiniteEtaleKeller
