/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GVC.Definitions
import GVC.CuspIdentity
import GVC.EndpointCoefficients

/-!
# The concrete three-variable witness

This file fixes the manuscript's actual polynomials in `ℚ[x,y,t]`.  The
remaining analytic/algebraic transfer is represented by
`ConcreteCounterexampleBridge`; unlike an axiom, a value of this structure
has to be provided explicitly before the counterexample theorem can be
applied.
-/

namespace GVC

open MvPolynomial Polynomial

abbrev TernaryPolynomial := MvPolynomial (Fin 3) ℚ

noncomputable def gvcX : TernaryPolynomial := X 0
noncomputable def gvcY : TernaryPolynomial := X 1
noncomputable def gvcT : TernaryPolynomial := X 2

noncomputable def gvcRho : TernaryPolynomial := gvcT ^ 2 + gvcX * gvcY
noncomputable def gvcA : TernaryPolynomial := gvcRho + gvcX ^ 2
noncomputable def gvcC : TernaryPolynomial :=
  gvcY * gvcRho ^ 2 - 2 * gvcX * gvcT ^ 2 * gvcRho -
    gvcX ^ 3 * gvcT ^ 2
noncomputable def gvcP : TernaryPolynomial := gvcA * gvcC ^ 2
noncomputable def gvcDelta : TernaryPolynomial :=
  4 * gvcX * gvcY + gvcT ^ 2
noncomputable def gvcLambda : TernaryPolynomial := gvcDelta ^ 6
noncomputable def gvcQ : TernaryPolynomial := gvcX ^ 2

theorem gvc_cusp_identity :
    gvcX * gvcC = gvcRho ^ 3 - gvcT ^ 2 * gvcA ^ 2 := by
  exact cusp_identity gvcX gvcY gvcT

/-- The Reynolds factor for a degree-`12m` homogeneous polynomial. -/
noncomputable def reynoldsPureScale (m : ℕ) : ℚ :=
  (2 : ℚ) ^ (6 * m) * Nat.factorial (6 * m) *
    Nat.doubleFactorial (12 * m + 1)

/-- The Reynolds factor for a degree-`12m+2` homogeneous polynomial. -/
noncomputable def reynoldsMixedScale (m : ℕ) : ℚ :=
  (2 : ℚ) ^ (6 * m + 1) * Nat.factorial (6 * m + 1) *
    Nat.doubleFactorial (12 * m + 3)

theorem reynoldsMixedScale_pos (m : ℕ) : 0 < reynoldsMixedScale m := by
  unfold reynoldsMixedScale
  positivity

theorem reynolds_scale_mul_cuspMoment (m : ℕ) :
    reynoldsMixedScale m * cuspMoment m = mixedDerivativeValue m := by
  unfold reynoldsMixedScale cuspMoment mixedDerivativeValue
  ring

/-- Exact interface still needed for a complete proof of Theorem 8.1.

`pure_phase_eq` is the homogeneous Reynolds identity plus the quadric phase
extraction for the pure contraction.  `mixed_phase_eq` is the same bridge
after applying one additional `Δ`. -/
structure ConcreteCounterexampleBridge where
  pure_phase_eq : ∀ m, 0 < m →
    differentialAction (gvcLambda ^ m) (gvcP ^ m) =
      MvPolynomial.C (reynoldsPureScale m *
        (endpointKernel m (cuspMoment m)
          (endpointPrimitiveTail m)).coeff m)
  mixed_phase_eq : ∀ m, 0 < m →
    differentialAction gvcDelta
      (differentialAction (gvcLambda ^ m) (gvcQ * gvcP ^ m)) =
      MvPolynomial.C (reynoldsMixedScale m *
        (endpointKernel m (cuspMoment m)
          (endpointPrimitiveTail m)).coeff (m - 1))

/-- The manuscript's pure identity, conditional only on the explicitly
displayed Reynolds/phase bridge. -/
theorem gvc3_pure_identity
    (B : ConcreteCounterexampleBridge) (m : ℕ) (hm : 0 < m) :
    differentialAction (gvcLambda ^ m) (gvcP ^ m) = 0 := by
  rw [B.pure_phase_eq m hm, endpointKernel_coeff_pure m hm]
  simp

/-- The paper's exact scalar after one additional `Δ`. -/
theorem gvc3_exact_next_mixed
    (B : ConcreteCounterexampleBridge) (m : ℕ) (hm : 0 < m) :
    differentialAction gvcDelta
      (differentialAction (gvcLambda ^ m) (gvcQ * gvcP ^ m)) =
      MvPolynomial.C (mixedDerivativeValue m) := by
  rw [B.mixed_phase_eq m hm, endpointKernel_coeff_mixed m hm,
    reynolds_scale_mul_cuspMoment]

/-- The mixed GVC output cannot vanish: applying one further `Δ` gives the
nonzero exact scalar above. -/
theorem gvc3_mixed_ne_zero
    (B : ConcreteCounterexampleBridge) (m : ℕ) (hm : 0 < m) :
    differentialAction (gvcLambda ^ m) (gvcQ * gvcP ^ m) ≠ 0 := by
  intro hzero
  have hnext := gvc3_exact_next_mixed B m hm
  rw [hzero, differentialAction_zero_right] at hnext
  have hC : MvPolynomial.C (mixedDerivativeValue m) ≠
      (0 : TernaryPolynomial) := by
    simp [mixedDerivativeValue_ne_zero]
  exact hC hnext.symm

theorem gvc3_purePowersVanish
    (B : ConcreteCounterexampleBridge) :
    PurePowersVanish gvcLambda gvcP := by
  intro m hm
  exact gvc3_pure_identity B m hm

/-- A supplied bridge turns the displayed pair into a literal refutation of
GVC for this symbol and polynomial. -/
theorem gvc3_not_generalizedVanishingFor
    (B : ConcreteCounterexampleBridge) :
    ¬ GeneralizedVanishingFor gvcLambda gvcP := by
  intro hGVC
  have heventual := hGVC (gvc3_purePowersVanish B) gvcQ
  obtain ⟨M, hM⟩ := heventual
  let m := max M 1
  have hmM : M ≤ m := le_max_left _ _
  have hm : 0 < m := lt_of_lt_of_le Nat.zero_lt_one (le_max_right _ _)
  exact gvc3_mixed_ne_zero B m hm (hM m hmM)

end GVC
