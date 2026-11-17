/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GVC.Definitions
import GVC.CuspIdentity
import GVC.EndpointCoefficients
import GVC.TopContraction
import Mathlib.RingTheory.MvPolynomial.Homogeneous

/-!
# The concrete three-variable witness

This file fixes the manuscript's actual polynomials in `ℚ[x,y,t]`.  The
final transfer is factored through `ConcreteCounterexampleBridge`.  The
structure is not an axiom, and `GVC.PhaseKernel` constructs its value from
the Laurent restriction and coefficientwise formal integration.
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

theorem gvcRho_isHomogeneous : gvcRho.IsHomogeneous 2 := by
  simpa [gvcRho, gvcX, gvcY, gvcT] using
    ((MvPolynomial.isHomogeneous_X ℚ (2 : Fin 3)).pow 2).add
      ((MvPolynomial.isHomogeneous_X ℚ (0 : Fin 3)).mul
        (MvPolynomial.isHomogeneous_X ℚ (1 : Fin 3)))

theorem gvcA_isHomogeneous : gvcA.IsHomogeneous 2 := by
  simpa [gvcA, gvcX] using gvcRho_isHomogeneous.add
    ((MvPolynomial.isHomogeneous_X ℚ (0 : Fin 3)).pow 2)

theorem gvcC_isHomogeneous : gvcC.IsHomogeneous 5 := by
  have hx := MvPolynomial.isHomogeneous_X ℚ (0 : Fin 3)
  have hy := MvPolynomial.isHomogeneous_X ℚ (1 : Fin 3)
  have ht := MvPolynomial.isHomogeneous_X ℚ (2 : Fin 3)
  have hfirst : (gvcY * gvcRho ^ 2).IsHomogeneous 5 := by
    simpa [gvcY] using hy.mul (gvcRho_isHomogeneous.pow 2)
  have hmiddle :
      (2 * gvcX * gvcT ^ 2 * gvcRho).IsHomogeneous 5 := by
    have h := hx.mul ((ht.pow 2).mul gvcRho_isHomogeneous)
    have htwo : (2 : TernaryPolynomial) = MvPolynomial.C (2 : ℚ) :=
      (MvPolynomial.C_eq_coe_nat (R := ℚ) 2).symm
    rw [htwo]
    simpa [gvcX, gvcT, mul_assoc] using h.C_mul (2 : ℚ)
  have hlast : (gvcX ^ 3 * gvcT ^ 2).IsHomogeneous 5 := by
    simpa [gvcX, gvcT] using (hx.pow 3).mul (ht.pow 2)
  exact hfirst.sub hmiddle |>.sub hlast

theorem gvcP_isHomogeneous : gvcP.IsHomogeneous 12 := by
  simpa [gvcP] using gvcA_isHomogeneous.mul (gvcC_isHomogeneous.pow 2)

theorem gvcDelta_isHomogeneous : gvcDelta.IsHomogeneous 2 := by
  have hx := MvPolynomial.isHomogeneous_X ℚ (0 : Fin 3)
  have hy := MvPolynomial.isHomogeneous_X ℚ (1 : Fin 3)
  have ht := MvPolynomial.isHomogeneous_X ℚ (2 : Fin 3)
  have hxy : (4 * gvcX * gvcY).IsHomogeneous 2 := by
    have h := hx.mul hy
    have hfour : (4 : TernaryPolynomial) = MvPolynomial.C (4 : ℚ) :=
      (MvPolynomial.C_eq_coe_nat (R := ℚ) 4).symm
    rw [hfour]
    simpa [gvcX, gvcY, mul_assoc] using h.C_mul (4 : ℚ)
  simpa [gvcDelta, gvcT] using hxy.add (ht.pow 2)

theorem gvcLambda_isHomogeneous : gvcLambda.IsHomogeneous 12 := by
  simpa [gvcLambda] using gvcDelta_isHomogeneous.pow 6

theorem gvcQ_isHomogeneous : gvcQ.IsHomogeneous 2 := by
  simpa [gvcQ, gvcX] using
    (MvPolynomial.isHomogeneous_X ℚ (0 : Fin 3)).pow 2

theorem gvcP_eval_witness :
    eval ![(1 : ℚ), 0, 1] gvcP = 18 := by
  norm_num [gvcP, gvcA, gvcC, gvcRho, gvcX, gvcY, gvcT,
    Matrix.cons_val_two]

theorem gvcP_ne_zero : gvcP ≠ 0 := by
  intro hzero
  have h := congrArg (MvPolynomial.eval ![(1 : ℚ), 0, 1]) hzero
  rw [gvcP_eval_witness] at h
  norm_num at h

theorem gvcP_totalDegree : gvcP.totalDegree = 12 :=
  gvcP_isHomogeneous.totalDegree gvcP_ne_zero

theorem gvcLambda_eval_witness :
    eval ![(1 : ℚ), 1, 0] gvcLambda = 4096 := by
  norm_num [gvcLambda, gvcDelta, gvcX, gvcY, gvcT,
    Matrix.cons_val_two]

theorem gvcLambda_ne_zero : gvcLambda ≠ 0 := by
  intro hzero
  have h := congrArg (MvPolynomial.eval ![(1 : ℚ), 1, 0]) hzero
  rw [gvcLambda_eval_witness] at h
  norm_num at h

theorem gvcLambda_totalDegree : gvcLambda.totalDegree = 12 :=
  gvcLambda_isHomogeneous.totalDegree gvcLambda_ne_zero

/-- The universal Reynolds normalization in three variables. -/
noncomputable def reynoldsScale (k : ℕ) : ℚ :=
  (2 : ℚ) ^ k * Nat.factorial k * Nat.doubleFactorial (2 * k + 1)

theorem reynoldsScale_pos (k : ℕ) : 0 < reynoldsScale k := by
  unfold reynoldsScale
  positivity

theorem reynoldsScale_ne_zero (k : ℕ) : reynoldsScale k ≠ 0 :=
  ne_of_gt (reynoldsScale_pos k)

/-- The algebraic Reynolds functional, defined coefficientwise by the
normalized top contraction. -/
noncomputable def algebraicReynoldsMoment
    (k : ℕ) (p : TernaryPolynomial) : ℚ :=
  topContraction (gvcDelta ^ k) p / reynoldsScale k

/-- Algebraic Reynolds identity for every homogeneous ternary polynomial of
the matching degree.  The spherical integral in the paper is a convenient
real realization of this coefficientwise functional. -/
theorem differentialAction_delta_pow_eq_reynolds
    (k : ℕ) (p : TernaryPolynomial) (hp : p.IsHomogeneous (2 * k)) :
    differentialAction (gvcDelta ^ k) p =
      MvPolynomial.C (reynoldsScale k * algebraicReynoldsMoment k p) := by
  rw [differentialAction_eq_C_topContraction_of_homogeneous
    (gvcDelta_isHomogeneous.pow k) hp]
  congr 1
  rw [algebraicReynoldsMoment]
  exact (mul_div_cancel₀ _ (reynoldsScale_ne_zero k)).symm

theorem gvcLambda_pow_eq_delta_pow (m : ℕ) :
    gvcLambda ^ m = gvcDelta ^ (6 * m) := by
  rw [gvcLambda, ← pow_mul]

theorem gvcP_pow_isHomogeneous (m : ℕ) :
    (gvcP ^ m).IsHomogeneous (2 * (6 * m)) := by
  convert gvcP_isHomogeneous.pow m using 1
  omega

theorem gvcQ_mul_gvcP_pow_isHomogeneous (m : ℕ) :
    (gvcQ * gvcP ^ m).IsHomogeneous (2 * (6 * m + 1)) := by
  convert gvcQ_isHomogeneous.mul (gvcP_isHomogeneous.pow m) using 1
  omega

theorem gvc_next_action_eq_delta_pow (m : ℕ) :
    differentialAction gvcDelta
        (differentialAction (gvcLambda ^ m) (gvcQ * gvcP ^ m)) =
      differentialAction (gvcDelta ^ (6 * m + 1))
        (gvcQ * gvcP ^ m) := by
  rw [gvcLambda_pow_eq_delta_pow, ← differentialAction_mul_left,
    ← pow_succ']

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

theorem reynoldsScale_pure (m : ℕ) :
    reynoldsScale (6 * m) = reynoldsPureScale m := by
  unfold reynoldsScale reynoldsPureScale
  congr 2
  apply congrArg Nat.doubleFactorial
  omega

theorem reynoldsScale_mixed (m : ℕ) :
    reynoldsScale (6 * m + 1) = reynoldsMixedScale m := by
  unfold reynoldsScale reynoldsMixedScale
  congr 2
  apply congrArg Nat.doubleFactorial
  omega

/-- Factored interface for the two concrete quadric phase-extraction
formulas.  `GVC.PhaseKernel` proves both fields and constructs this record. -/
structure ConcreteCounterexampleBridge where
  pure_phase_eq : ∀ m, 0 < m →
    algebraicReynoldsMoment (6 * m) (gvcP ^ m) =
      (endpointKernel m (cuspMoment m)
        (endpointPrimitiveTail m)).coeff m
  mixed_phase_eq : ∀ m, 0 < m →
    algebraicReynoldsMoment (6 * m + 1) (gvcQ * gvcP ^ m) =
      (endpointKernel m (cuspMoment m)
        (endpointPrimitiveTail m)).coeff (m - 1)

/-- The manuscript's pure identity, conditional only on the explicitly
displayed quadric phase-extraction bridge. -/
theorem gvc3_pure_identity
    (B : ConcreteCounterexampleBridge) (m : ℕ) (hm : 0 < m) :
    differentialAction (gvcLambda ^ m) (gvcP ^ m) = 0 := by
  rw [gvcLambda_pow_eq_delta_pow,
    differentialAction_delta_pow_eq_reynolds _ _
      (gvcP_pow_isHomogeneous m),
    B.pure_phase_eq m hm, reynoldsScale_pure,
    endpointKernel_coeff_pure m hm]
  simp

/-- The paper's exact scalar after one additional `Δ`. -/
theorem gvc3_exact_next_mixed
    (B : ConcreteCounterexampleBridge) (m : ℕ) (hm : 0 < m) :
    differentialAction gvcDelta
      (differentialAction (gvcLambda ^ m) (gvcQ * gvcP ^ m)) =
      MvPolynomial.C (mixedDerivativeValue m) := by
  rw [gvc_next_action_eq_delta_pow,
    differentialAction_delta_pow_eq_reynolds _ _
      (gvcQ_mul_gvcP_pow_isHomogeneous m),
    B.mixed_phase_eq m hm, reynoldsScale_mixed,
    endpointKernel_coeff_mixed m hm,
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

/-- Any proof of the factored bridge turns the displayed pair into a literal
refutation of GVC.  `GVC.PhaseKernel` supplies such a proof. -/
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
