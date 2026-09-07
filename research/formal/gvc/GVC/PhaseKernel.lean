/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GVC.QuadricPhase

/-!
# The concrete even-phase kernel

This file converts the Laurent constant terms for the concrete ternary
witness into coefficients of the ordinary phase polynomial in `u = z²`.
The coefficients remain polynomials in the height variable.
-/

namespace GVC

open Polynomial
open scoped LaurentPolynomial

/-- Ordinary phase polynomials in `u`, with height-polynomial coefficients. -/
abbrev HeightPhasePolynomial := Polynomial (Polynomial ℚ)

/-- The polynomial before coefficientwise height integration:
`(1+u)^m (1-t²(1+u)²)^(2m)`. -/
noncomputable def heightPhaseKernel (m : ℕ) : HeightPhasePolynomial :=
  (1 + Polynomial.X) ^ m *
    (1 - Polynomial.C (Polynomial.X ^ 2) *
      (1 + Polynomial.X) ^ 2) ^ (2 * m)

/-- Substitute `u = z²` into an ordinary phase polynomial. -/
noncomputable def phaseToLaurent : HeightPhasePolynomial →+* PhaseLaurent :=
  Polynomial.eval₂RingHom LaurentPolynomial.C (LaurentPolynomial.T 2)

@[simp] theorem phaseToLaurent_X :
    phaseToLaurent Polynomial.X = LaurentPolynomial.T 2 := by
  simp [phaseToLaurent]

@[simp] theorem phaseToLaurent_C (a : Polynomial ℚ) :
    phaseToLaurent (Polynomial.C a) = LaurentPolynomial.C a := by
  simp [phaseToLaurent]

theorem phaseToLaurent_monomial (n : ℕ) (a : Polynomial ℚ) :
    phaseToLaurent (Polynomial.monomial n a) =
      LaurentPolynomial.C a * LaurentPolynomial.T (2 * n : ℤ) := by
  rw [← Polynomial.C_mul_X_pow_eq_monomial, map_mul, map_pow,
    phaseToLaurent_C, phaseToLaurent_X, LaurentPolynomial.T_pow]
  congr 2
  ring

/-- A negative even Laurent shift extracts the corresponding ordinary
`u`-coefficient. -/
theorem coeff_zero_T_neg_two_mul_phaseToLaurent
    (m : ℕ) (p : HeightPhasePolynomial) :
    ((LaurentPolynomial.T (-2 * (m : ℤ)) : PhaseLaurent) *
        phaseToLaurent p).coeff 0 = p.coeff m := by
  induction p using Polynomial.induction_on' with
  | add p q ihp ihq =>
      rw [map_add, mul_add]
      change
        ((LaurentPolynomial.T (-2 * (m : ℤ)) : PhaseLaurent) *
              phaseToLaurent p).coeff 0 +
            ((LaurentPolynomial.T (-2 * (m : ℤ)) : PhaseLaurent) *
              phaseToLaurent q).coeff 0 =
          p.coeff m + q.coeff m
      rw [ihp, ihq]
  | monomial n a =>
      rw [phaseToLaurent_monomial]
      have hlaurent :
          (LaurentPolynomial.T (-2 * (m : ℤ)) : PhaseLaurent) *
              (LaurentPolynomial.C a *
                LaurentPolynomial.T (2 * n : ℤ)) =
            LaurentPolynomial.C a *
              LaurentPolynomial.T (-2 * (m : ℤ) + 2 * (n : ℤ)) := by
        calc
          _ = LaurentPolynomial.C a *
                LaurentPolynomial.T (-2 * (m : ℤ)) *
                  LaurentPolynomial.T (2 * n : ℤ) := by ac_rfl
          _ = _ := by rw [LaurentPolynomial.mul_T_assoc]
      rw [hlaurent,
        ← LaurentPolynomial.single_eq_C_mul_T]
      change
        (Finsupp.single (-2 * (m : ℤ) + 2 * (n : ℤ)) a :
            ℤ →₀ Polynomial ℚ) 0 =
          (Polynomial.monomial n a).coeff m
      rw [Finsupp.single_apply, Polynomial.coeff_monomial]
      by_cases h : n = m
      · subst n
        simp
      · have hz : -(2 * (m : ℤ)) + 2 * (n : ℤ) ≠ 0 := by
          omega
        simp [hz, h]

/-- The adjacent positive shift extracts coefficient `m-1`. -/
theorem coeff_zero_T_two_mul_T_neg_two_mul_phaseToLaurent
    (m : ℕ) (p : HeightPhasePolynomial) (hm : 0 < m) :
    ((LaurentPolynomial.T 2 : PhaseLaurent) *
        (LaurentPolynomial.T (-2 * (m : ℤ)) * phaseToLaurent p)).coeff 0 =
      p.coeff (m - 1) := by
  induction p using Polynomial.induction_on' with
  | add p q ihp ihq =>
      rw [map_add, mul_add, mul_add]
      change
        ((LaurentPolynomial.T 2 : PhaseLaurent) *
              (LaurentPolynomial.T (-2 * (m : ℤ)) *
                phaseToLaurent p)).coeff 0 +
            ((LaurentPolynomial.T 2 : PhaseLaurent) *
              (LaurentPolynomial.T (-2 * (m : ℤ)) *
                phaseToLaurent q)).coeff 0 =
          p.coeff (m - 1) + q.coeff (m - 1)
      rw [ihp, ihq]
  | monomial n a =>
      rw [phaseToLaurent_monomial]
      have hlaurent :
          (LaurentPolynomial.T 2 : PhaseLaurent) *
              (LaurentPolynomial.T (-2 * (m : ℤ)) *
                (LaurentPolynomial.C a *
                  LaurentPolynomial.T (2 * n : ℤ))) =
            LaurentPolynomial.C a *
              LaurentPolynomial.T (2 - 2 * (m : ℤ) + 2 * (n : ℤ)) := by
        calc
          _ = LaurentPolynomial.C a * LaurentPolynomial.T 2 *
                LaurentPolynomial.T (-2 * (m : ℤ)) *
                  LaurentPolynomial.T (2 * n : ℤ) := by ac_rfl
          _ = LaurentPolynomial.C a *
                LaurentPolynomial.T
                  ((2 : ℤ) + (-2 * (m : ℤ) + 2 * (n : ℤ))) := by
              rw [LaurentPolynomial.mul_T_assoc,
                LaurentPolynomial.mul_T_assoc]
          _ = _ := by
            congr 2
            ring
      rw [hlaurent, ← LaurentPolynomial.single_eq_C_mul_T]
      change
        (Finsupp.single
            (2 - 2 * (m : ℤ) + 2 * (n : ℤ)) a :
          ℤ →₀ Polynomial ℚ) 0 =
            (Polynomial.monomial n a).coeff (m - 1)
      rw [Finsupp.single_apply, Polynomial.coeff_monomial]
      by_cases h : n = m - 1
      · subst n
        have hz :
            2 - 2 * (m : ℤ) + 2 * ((m - 1 : ℕ) : ℤ) = 0 := by
          omega
        simp [hz]
      · have hz :
            2 - 2 * (m : ℤ) + 2 * (n : ℤ) ≠ 0 := by
          intro hz
          apply h
          omega
        simp [hz, h]

/-- The Laurent restriction of `P^m` is the even phase kernel with its
overall winding shift exposed. -/
theorem quadricPhaseRestriction_gvcP_pow_eq_heightPhaseKernel (m : ℕ) :
    quadricPhaseRestriction (gvcP ^ m) =
      LaurentPolynomial.T (-2 * (m : ℤ)) *
        phaseToLaurent (heightPhaseKernel m) := by
  rw [quadricPhaseRestriction_gvcP_pow]
  have hmap :
      phaseToLaurent (heightPhaseKernel m) =
        (1 + LaurentPolynomial.T 2) ^ m *
          (1 - LaurentPolynomial.C (Polynomial.X ^ 2) *
            (1 + LaurentPolynomial.T 2) ^ 2) ^ (2 * m) := by
    simp only [heightPhaseKernel, map_mul, map_pow, map_add, map_one,
      map_sub, phaseToLaurent_X, phaseToLaurent_C]
  rw [hmap]
  rw [mul_pow, mul_pow, LaurentPolynomial.T_pow, ← pow_mul]
  have hexponent : (m : ℤ) * (-2) = -2 * (m : ℤ) := by ring
  rw [hexponent, mul_assoc]

theorem quadricPhaseConstant_gvcP_pow_eq_heightPhaseKernel_coeff (m : ℕ) :
    quadricPhaseConstant (gvcP ^ m) =
      (heightPhaseKernel m).coeff m := by
  rw [quadricPhaseConstant,
    quadricPhaseRestriction_gvcP_pow_eq_heightPhaseKernel,
    coeff_zero_T_neg_two_mul_phaseToLaurent]

theorem quadricPhaseConstant_gvcQ_mul_gvcP_pow_eq_heightPhaseKernel_coeff
    (m : ℕ) (hm : 0 < m) :
    quadricPhaseConstant (gvcQ * gvcP ^ m) =
      (heightPhaseKernel m).coeff (m - 1) := by
  rw [quadricPhaseConstant, map_mul, quadricPhaseRestriction_gvcQ,
    quadricPhaseRestriction_gvcP_pow_eq_heightPhaseKernel,
    coeff_zero_T_two_mul_T_neg_two_mul_phaseToLaurent m _ hm]

/-! ## Coefficientwise height integration -/

/-- Apply the formal integral on `[0,1]` independently to every height
coefficient of a phase polynomial. -/
noncomputable def integrateHeightCoefficients
    (p : HeightPhasePolynomial) : Polynomial ℚ :=
  p.sum fun n a ↦ Polynomial.monomial n (formalIntegral01 a)

@[simp] theorem integrateHeightCoefficients_monomial
    (n : ℕ) (a : Polynomial ℚ) :
    integrateHeightCoefficients (Polynomial.monomial n a) =
      Polynomial.monomial n (formalIntegral01 a) := by
  classical
  rw [integrateHeightCoefficients, Polynomial.sum_monomial_index]
  simp [formalIntegral01, polynomialPrimitive]

@[simp] theorem integrateHeightCoefficients_add
    (p q : HeightPhasePolynomial) :
    integrateHeightCoefficients (p + q) =
      integrateHeightCoefficients p + integrateHeightCoefficients q := by
  rw [integrateHeightCoefficients, integrateHeightCoefficients,
    integrateHeightCoefficients]
  apply Polynomial.sum_add_index
  · intro i
    simp [formalIntegral01, polynomialPrimitive]
  · intro i a b
    rw [formalIntegral01_add, map_add]

theorem integrateHeightCoefficients_coeff
    (p : HeightPhasePolynomial) (n : ℕ) :
    (integrateHeightCoefficients p).coeff n =
      formalIntegral01 (p.coeff n) := by
  induction p using Polynomial.induction_on' with
  | add p q ihp ihq =>
      rw [integrateHeightCoefficients_add, Polynomial.coeff_add,
        ihp, ihq, Polynomial.coeff_add, formalIntegral01_add]
  | monomial d a =>
      rw [integrateHeightCoefficients_monomial,
        Polynomial.coeff_monomial, Polynomial.coeff_monomial]
      split_ifs <;> simp [formalIntegral01, polynomialPrimitive]

/-- Phase polynomials with rational coefficients can be pulled through
coefficientwise height integration. -/
theorem integrateHeightCoefficients_map_mul_C
    (q a : Polynomial ℚ) :
    integrateHeightCoefficients
        (q.map Polynomial.C * Polynomial.C a) =
      q * Polynomial.C (formalIntegral01 a) := by
  induction q using Polynomial.induction_on' with
  | add p q ihp ihq =>
      rw [Polynomial.map_add, add_mul, integrateHeightCoefficients_add, ihp, ihq,
        add_mul]
  | monomial n c =>
      rw [Polynomial.map_monomial]
      have hproduct :
          Polynomial.monomial n (Polynomial.C c) * Polynomial.C a =
            Polynomial.monomial n (Polynomial.C c * a) := by
        change
          Polynomial.monomial n (Polynomial.C c) *
              Polynomial.monomial 0 a =
            Polynomial.monomial n (Polynomial.C c * a)
        rw [Polynomial.monomial_mul_monomial]
        simp
      rw [hproduct, integrateHeightCoefficients_monomial,
        formalIntegral01_C_mul]
      change
        Polynomial.monomial n (c * formalIntegral01 a) =
          Polynomial.monomial n c *
            Polynomial.monomial 0 (formalIntegral01 a)
      rw [Polynomial.monomial_mul_monomial]
      simp

theorem integrateHeightCoefficients_map_mul
    (q : Polynomial ℚ) (p : HeightPhasePolynomial) :
    integrateHeightCoefficients (q.map Polynomial.C * p) =
      q * integrateHeightCoefficients p := by
  induction p using Polynomial.induction_on' with
  | add p r ihp ihr =>
      rw [mul_add, integrateHeightCoefficients_add, ihp, ihr,
        integrateHeightCoefficients_add, mul_add]
  | monomial n a =>
      have hproduct :
          q.map Polynomial.C * Polynomial.monomial n a =
            (q * Polynomial.X ^ n).map Polynomial.C *
              Polynomial.C a := by
        rw [← Polynomial.C_mul_X_pow_eq_monomial,
          Polynomial.map_mul, Polynomial.map_pow]
        simp
        ring
      rw [hproduct, integrateHeightCoefficients_map_mul_C,
        integrateHeightCoefficients_monomial,
        ← Polynomial.C_mul_X_pow_eq_monomial]
      ring

/-! ## Algebraic change of height variable -/

/-- Substitute the bivariate height `t * B(u)` into a rational polynomial
in `t`, retaining `u` as the outer phase variable. -/
noncomputable def scaleHeightByPhase
    (B f : Polynomial ℚ) : HeightPhasePolynomial :=
  Polynomial.eval₂
    ((Polynomial.C : Polynomial ℚ →+* HeightPhasePolynomial).comp
      (Polynomial.C : ℚ →+* Polynomial ℚ))
    (Polynomial.C Polynomial.X * B.map Polynomial.C) f

@[simp] theorem scaleHeightByPhase_add
    (B f g : Polynomial ℚ) :
    scaleHeightByPhase B (f + g) =
      scaleHeightByPhase B f + scaleHeightByPhase B g := by
  simp [scaleHeightByPhase]

theorem scaleHeightByPhase_monomial
    (B : Polynomial ℚ) (n : ℕ) (a : ℚ) :
    scaleHeightByPhase B (Polynomial.monomial n a) =
      (B ^ n).map Polynomial.C *
        Polynomial.C (Polynomial.C a * Polynomial.X ^ n) := by
  rw [← Polynomial.C_mul_X_pow_eq_monomial]
  simp only [scaleHeightByPhase, Polynomial.eval₂_mul,
    Polynomial.eval₂_pow, Polynomial.eval₂_C, Polynomial.eval₂_X,
    RingHom.coe_comp, Function.comp_apply]
  rw [mul_pow, Polynomial.map_pow, map_mul, map_pow]
  ac_rfl

theorem polynomialPrimitive_monomial (n : ℕ) (a : ℚ) :
    polynomialPrimitive (Polynomial.monomial n a) =
      Polynomial.monomial (n + 1) (a / (n + 1)) := by
  classical
  rw [polynomialPrimitive, Polynomial.sum_monomial_index]
  simp

/-- Polynomial change of variables for the formal height integral:
`B * ∫ f(tB) dt` is the primitive of `f` evaluated at `B`. -/
theorem scaleHeightByPhase_integral
    (B f : Polynomial ℚ) :
    B * integrateHeightCoefficients (scaleHeightByPhase B f) =
      (polynomialPrimitive f).comp B := by
  induction f using Polynomial.induction_on' with
  | add f g ihf ihg =>
      rw [scaleHeightByPhase_add, integrateHeightCoefficients_add,
        mul_add, ihf, ihg, polynomialPrimitive_add, Polynomial.add_comp]
  | monomial n a =>
      rw [scaleHeightByPhase_monomial,
        integrateHeightCoefficients_map_mul_C,
        formalIntegral01_C_mul, formalIntegral01_X_pow,
        polynomialPrimitive_monomial,
        ← Polynomial.C_mul_X_pow_eq_monomial]
      rw [Polynomial.mul_comp, Polynomial.C_comp,
        Polynomial.X_pow_comp]
      rw [pow_succ]
      have hcoeff :
          a * (1 / ((n : ℚ) + 1)) = a / ((n : ℚ) + 1) := by
        simp only [div_eq_mul_inv, one_mul]
      rw [hcoeff]
      ac_rfl

theorem heightPhaseKernel_eq_scaleHeightByPhase (m : ℕ) :
    heightPhaseKernel m =
      (((1 + Polynomial.X) ^ m : Polynomial ℚ).map Polynomial.C) *
        scaleHeightByPhase (1 + Polynomial.X)
          ((1 - Polynomial.X ^ 2) ^ (2 * m)) := by
  simp only [heightPhaseKernel, scaleHeightByPhase,
    Polynomial.eval₂_pow, Polynomial.eval₂_sub,
    Polynomial.eval₂_one,
    Polynomial.eval₂_X, Polynomial.map_pow,
    Polynomial.map_add, Polynomial.map_one, Polynomial.map_X]
  have hbase :
      Polynomial.C (Polynomial.X ^ 2) *
          (1 + Polynomial.X : HeightPhasePolynomial) ^ 2 =
        (Polynomial.C Polynomial.X * (1 + Polynomial.X)) ^ 2 := by
    rw [map_pow, mul_pow]
  rw [hbase]

/-- The shifted endpoint primitive is literally the primitive of
`(1-t²)^(2m)` evaluated at `1+u`. -/
theorem polynomialPrimitive_one_sub_sq_comp_one_add_X (m : ℕ) :
    (polynomialPrimitive ((1 - Polynomial.X ^ 2) ^ (2 * m))).comp
        (1 + Polynomial.X) =
      shiftedEndpointPrimitive m := by
  let F : Polynomial ℚ :=
    (polynomialPrimitive ((1 - Polynomial.X ^ 2) ^ (2 * m))).comp
      (1 + Polynomial.X)
  have hderivative :
      derivative F = derivative (shiftedEndpointPrimitive m) := by
    dsimp [F]
    rw [Polynomial.derivative_comp, derivative_polynomialPrimitive,
      derivative_shiftedEndpointPrimitive,
      ← endpoint_derivative_factorization]
    simp
  have hconstant :
      (F - shiftedEndpointPrimitive m).coeff 0 = 0 := by
    rw [Polynomial.coeff_sub, Polynomial.coeff_zero_eq_eval_zero,
      Polynomial.coeff_zero_eq_eval_zero]
    dsimp [F]
    rw [Polynomial.eval_comp]
    simp only [eval_add, eval_one, eval_X, add_zero,
      shiftedEndpointPrimitive, eval_C]
    rw [cuspMoment_eq_formalIntegral, formalIntegral01]
    rw [← Polynomial.coeff_zero_eq_eval_zero,
      polynomialPrimitive_coeff_zero]
    ring
  have hderivative_zero :
      derivative (F - shiftedEndpointPrimitive m) = 0 := by
    rw [derivative_sub, hderivative, sub_self]
  have hconstant_poly :=
    Polynomial.eq_C_of_derivative_eq_zero hderivative_zero
  rw [hconstant, map_zero] at hconstant_poly
  exact sub_eq_zero.mp hconstant_poly

/-- Coefficientwise height integration of the concrete phase polynomial is
exactly the endpoint kernel already analyzed in `EndpointCoefficients`. -/
theorem integrateHeightCoefficients_heightPhaseKernel
    (m : ℕ) (hm : 0 < m) :
    integrateHeightCoefficients (heightPhaseKernel m) =
      endpointKernel m (cuspMoment m) (endpointPrimitiveTail m) := by
  rw [heightPhaseKernel_eq_scaleHeightByPhase,
    integrateHeightCoefficients_map_mul,
    endpointKernel_actual_eq]
  have hscale := scaleHeightByPhase_integral
    (1 + Polynomial.X) ((1 - Polynomial.X ^ 2) ^ (2 * m))
  rw [polynomialPrimitive_one_sub_sq_comp_one_add_X] at hscale
  have hmexp : m = (m - 1) + 1 := by omega
  have hpow :
      ((1 + Polynomial.X) ^ m : Polynomial ℚ) =
        (1 + Polynomial.X) ^ (m - 1) * (1 + Polynomial.X) := by
    calc
      ((1 + Polynomial.X) ^ m : Polynomial ℚ) =
          (1 + Polynomial.X) ^ ((m - 1) + 1) :=
        congrArg (fun e : ℕ ↦ ((1 + Polynomial.X) ^ e : Polynomial ℚ)) hmexp
      _ = (1 + Polynomial.X) ^ (m - 1) *
          (1 + Polynomial.X) := by rw [pow_succ]
  rw [hpow, mul_assoc, hscale]

/-- The concrete phase bridge is now inhabited entirely by kernel-checked
proofs. -/
theorem verifiedConcreteCounterexampleBridge :
    ConcreteCounterexampleBridge where
  pure_phase_eq m hm := by
    rw [algebraicReynoldsMoment_eq_formalIntegral_quadricPhaseConstant
        (6 * m) (gvcP_pow_isHomogeneous m),
      quadricPhaseConstant_gvcP_pow_eq_heightPhaseKernel_coeff,
      ← integrateHeightCoefficients_coeff,
      integrateHeightCoefficients_heightPhaseKernel m hm]
  mixed_phase_eq m hm := by
    rw [algebraicReynoldsMoment_eq_formalIntegral_quadricPhaseConstant
        (6 * m + 1) (gvcQ_mul_gvcP_pow_isHomogeneous m),
      quadricPhaseConstant_gvcQ_mul_gvcP_pow_eq_heightPhaseKernel_coeff m hm,
      ← integrateHeightCoefficients_coeff,
      integrateHeightCoefficients_heightPhaseKernel m hm]

/-- Unconditional rational form of Theorem 8.1's pure identity. -/
theorem verified_gvc3_pure_identity (m : ℕ) (hm : 0 < m) :
    differentialAction (gvcLambda ^ m) (gvcP ^ m) = 0 :=
  gvc3_pure_identity verifiedConcreteCounterexampleBridge m hm

/-- Unconditional rational form of the exact adjacent scalar. -/
theorem verified_gvc3_exact_next_mixed (m : ℕ) (hm : 0 < m) :
    differentialAction gvcDelta
      (differentialAction (gvcLambda ^ m) (gvcQ * gvcP ^ m)) =
      MvPolynomial.C (mixedDerivativeValue m) :=
  gvc3_exact_next_mixed verifiedConcreteCounterexampleBridge m hm

theorem verified_gvc3_mixed_ne_zero (m : ℕ) (hm : 0 < m) :
    differentialAction (gvcLambda ^ m) (gvcQ * gvcP ^ m) ≠ 0 :=
  gvc3_mixed_ne_zero verifiedConcreteCounterexampleBridge m hm

theorem verified_gvc3_purePowersVanish :
    PurePowersVanish gvcLambda gvcP :=
  gvc3_purePowersVanish verifiedConcreteCounterexampleBridge

/-- The displayed rational pair is an unconditional Lean-verified
counterexample to GVC in three variables. -/
theorem verified_gvc3_not_generalizedVanishingFor :
    ¬ GeneralizedVanishingFor gvcLambda gvcP :=
  gvc3_not_generalizedVanishingFor verifiedConcreteCounterexampleBridge

end GVC
