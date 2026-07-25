/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeRealization
import FiniteEtaleKeller.GeneralGaugeDegree

/-!
# Effective degree of the final automatic realization map

Translation preserves polynomial degree, and removing the constant term does
not change it in positive degree.  Combining this with the all-degree map bound
gives the paper's `6N+2` estimate directly for the actual automatically chosen
map attached to the original polynomial `P`.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K] [CharZero K]

/-- Removing the constant term from a translated positive-degree polynomial
does not change its degree. -/
theorem rootedTranslate_natDegree
    (P : K[X]) (a : K) (hpos : 0 < P.natDegree) :
    (rootedTranslate P a).natDegree = P.natDegree := by
  have hupper : (rootedTranslate P a).natDegree ≤ P.natDegree := by
    change (Polynomial.taylor a P - C (P.eval a)).natDegree ≤ P.natDegree
    calc
      (Polynomial.taylor a P - C (P.eval a)).natDegree ≤
          max (Polynomial.taylor a P).natDegree (C (P.eval a)).natDegree :=
        Polynomial.natDegree_sub_le
      _ = P.natDegree := by simp
  have hP : P ≠ 0 := by
    intro h
    simp [h] at hpos
  have hcoeff : (rootedTranslate P a).coeff P.natDegree ≠ 0 := by
    change (Polynomial.taylor a P - C (P.eval a)).coeff P.natDegree ≠ 0
    rw [Polynomial.coeff_sub, Polynomial.coeff_taylor_natDegree]
    have hne : P.natDegree ≠ 0 := Nat.ne_of_gt hpos
    simp [hne, Polynomial.leadingCoeff_ne_zero.mpr hP]
  exact le_antisymm hupper (Polynomial.le_natDegree_of_ne_zero hcoeff)

@[simp]
theorem realizationSeed_natDegree
    (P : K[X]) (a : K) (hpos : 0 < P.natDegree) :
    (realizationSeed P a).natDegree = P.natDegree := by
  exact rootedTranslate_natDegree P a hpos

/-- Every coordinate of the final automatically chosen determinant-one map has
total degree at most `6 * P.natDegree + 2`. -/
theorem automaticRealizationMap_totalDegree
    (P : K[X]) (hdeg : 3 ≤ P.natDegree) (i : Fin 3) :
    (automaticRealizationMap P hdeg i).totalDegree ≤
      6 * P.natDegree + 2 := by
  let a := chosenAdmissibleTranslation P hdeg
  have hpos : 0 < P.natDegree := by omega
  have hseed : (realizationSeed P a).natDegree = P.natDegree :=
    realizationSeed_natDegree P a hpos
  have hseeddeg : 3 ≤ (realizationSeed P a).natDegree := by
    rw [hseed]
    exact hdeg
  change (generalGaugeJacobianOneMap (realizationSeed P a) i).totalDegree ≤
    6 * P.natDegree + 2
  have h := generalGaugeJacobianOneMap_totalDegree
    (realizationSeed P a) hseeddeg i
  rwa [hseed] at h

#print axioms rootedTranslate_natDegree
#print axioms automaticRealizationMap_totalDegree

end FiniteEtaleKeller
