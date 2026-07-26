/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeJacobian

/-!
# The announced three-dimensional counterexample

The universal quadratic gauge contains the map announced by Levent Alpöge as
its smallest seed.  Taking `G(S) = S + S^3` gives exactly the three displayed
coordinates of the announcement, not merely a linearly equivalent map.
-/

noncomputable section

open Matrix Function
open MvPolynomial
open Polynomial

namespace FiniteEtaleKeller

/-- The cubic seed underlying the announced counterexample. -/
def announcedSeed : ℚ[X] := X + X ^ 3

/-- The polynomial map in the original announcement. -/
def announcedCounterexampleMap : Fin 3 → GaugePolynomial ℚ :=
  let t : GaugePolynomial ℚ := 1 + MvPolynomial.X 0 * MvPolynomial.X 1
  ![
    t ^ 3 * MvPolynomial.X 2 +
      MvPolynomial.X 1 ^ 2 * t * (MvPolynomial.C 4 + MvPolynomial.C 3 *
        MvPolynomial.X 0 * MvPolynomial.X 1),
    MvPolynomial.X 1 +
      MvPolynomial.C 3 * MvPolynomial.X 0 * t ^ 2 * MvPolynomial.X 2 +
      MvPolynomial.C 3 * MvPolynomial.X 0 * MvPolynomial.X 1 ^ 2 *
        (MvPolynomial.C 4 + MvPolynomial.C 3 *
          MvPolynomial.X 0 * MvPolynomial.X 1),
    MvPolynomial.C 2 * MvPolynomial.X 0 -
      MvPolynomial.C 3 * MvPolynomial.X 0 ^ 2 * MvPolynomial.X 1 -
      MvPolynomial.X 0 ^ 3 * MvPolynomial.X 2]

private theorem announcedSeed_natDegree : announcedSeed.natDegree = 3 := by
  unfold announcedSeed
  compute_degree!

/-- Substituting `G(S)=S+S^3` into the all-degree construction recovers the
announced map coordinate for coordinate. -/
theorem generalGaugeMap_announcedSeed :
    generalGaugeMap announcedSeed = announcedCounterexampleMap := by
  funext i
  fin_cases i <;>
    simp [announcedCounterexampleMap, announcedSeed, generalGaugeMap,
      generalGaugePi, generalGaugeB, generalGaugeC, generalGaugeQ,
      generalGaugeT, announcedSeed_natDegree, Polynomial.coeff_X] <;>
    ring

/-- The constant Jacobian of the announced map is inherited from the general
theorem. -/
theorem jacobianDet_announcedCounterexampleMap :
    jacobianDet announcedCounterexampleMap = MvPolynomial.C (-2) := by
  rw [← generalGaugeMap_announcedSeed]
  apply jacobianDet_generalGaugeMap announcedSeed
  · norm_num [announcedSeed, Polynomial.coeff_X]
  · norm_num [announcedSeed, Polynomial.coeff_X]

#print axioms generalGaugeMap_announcedSeed
#print axioms jacobianDet_announcedCounterexampleMap

end FiniteEtaleKeller
