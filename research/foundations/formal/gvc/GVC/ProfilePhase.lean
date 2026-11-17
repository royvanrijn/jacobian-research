/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GVC.BaseChange
import GVC.PhaseKernel
import GVC.ProfileFamily

/-!
# The full cusp-profile phase bridge

This file discharges the phase-extraction interface left by
`GVC.ProfileFamily`.  The only additional hypothesis is the manuscript's
declared-degree condition on the profile polynomial.
-/

namespace GVC

open MvPolynomial Polynomial
open scoped LaurentPolynomial

/-- The height/phase polynomial before coefficientwise height integration
for the full winding--profile family. -/
noncomputable def profileHeightPhaseKernel
    (r m : ℕ) (S : ℚ[X]) : HeightPhasePolynomial :=
  (1 + Polynomial.X) ^ (r * m) *
    (1 - Polynomial.C (Polynomial.X ^ 2) *
      (1 + Polynomial.X) ^ 2) ^ (2 * r * m) *
    scaleHeightByPhase (1 + Polynomial.X)
      ((S.comp (Polynomial.X ^ 2)) ^ m)

/-- Any nonnegative even Laurent shift extracts the corresponding lower
ordinary phase coefficient. -/
theorem coeff_zero_T_two_mul_T_neg_two_mul_phaseToLaurent_general
    (n ell : ℕ) (p : HeightPhasePolynomial) (hell : ell ≤ n) :
    ((LaurentPolynomial.T (2 * (ell : ℤ)) : PhaseLaurent) *
        (LaurentPolynomial.T (-2 * (n : ℤ)) * phaseToLaurent p)).coeff 0 =
      p.coeff (n - ell) := by
  induction p using Polynomial.induction_on' with
  | add p q ihp ihq =>
      rw [map_add, mul_add, mul_add]
      change
        ((LaurentPolynomial.T (2 * (ell : ℤ)) : PhaseLaurent) *
              (LaurentPolynomial.T (-2 * (n : ℤ)) *
                phaseToLaurent p)).coeff 0 +
            ((LaurentPolynomial.T (2 * (ell : ℤ)) : PhaseLaurent) *
              (LaurentPolynomial.T (-2 * (n : ℤ)) *
                phaseToLaurent q)).coeff 0 =
          p.coeff (n - ell) + q.coeff (n - ell)
      rw [ihp, ihq]
  | monomial d a =>
      rw [phaseToLaurent_monomial]
      have hlaurent :
          (LaurentPolynomial.T (2 * (ell : ℤ)) : PhaseLaurent) *
              (LaurentPolynomial.T (-2 * (n : ℤ)) *
                (LaurentPolynomial.C a *
                  LaurentPolynomial.T (2 * d : ℤ))) =
            LaurentPolynomial.C a *
              LaurentPolynomial.T
                (2 * (ell : ℤ) - 2 * (n : ℤ) + 2 * (d : ℤ)) := by
        calc
          _ = LaurentPolynomial.C a *
                LaurentPolynomial.T (2 * (ell : ℤ)) *
                LaurentPolynomial.T (-2 * (n : ℤ)) *
                  LaurentPolynomial.T (2 * d : ℤ) := by ac_rfl
          _ = LaurentPolynomial.C a *
                LaurentPolynomial.T
                  (2 * (ell : ℤ) +
                    (-2 * (n : ℤ) + 2 * (d : ℤ))) := by
              rw [LaurentPolynomial.mul_T_assoc,
                LaurentPolynomial.mul_T_assoc]
          _ = _ := by
            congr 2
            ring
      rw [hlaurent, ← LaurentPolynomial.single_eq_C_mul_T]
      change
        (Finsupp.single
            (2 * (ell : ℤ) - 2 * (n : ℤ) + 2 * (d : ℤ)) a :
          ℤ →₀ Polynomial ℚ) 0 =
            (Polynomial.monomial d a).coeff (n - ell)
      rw [Finsupp.single_apply, Polynomial.coeff_monomial]
      by_cases h : d = n - ell
      · subst d
        have hz :
            2 * (ell : ℤ) - 2 * (n : ℤ) +
                2 * ((n - ell : ℕ) : ℤ) = 0 := by
          omega
        simp [hz]
      · have hz :
            2 * (ell : ℤ) - 2 * (n : ℤ) + 2 * (d : ℤ) ≠ 0 := by
          intro hz
          apply h
          omega
        simp [hz, h]

/-- On the split quadric, the declared homogeneous profile becomes the
literal univariate profile evaluated at `t²(1+z²)²`. -/
theorem quadricPhaseRestriction_gvcProfileHom
    (e : ℕ) (S : ℚ[X]) (hS : S.natDegree ≤ e) :
    quadricPhaseRestriction (gvcProfileHom e S) =
      Polynomial.eval₂
        ((LaurentPolynomial.C : Polynomial ℚ →+* PhaseLaurent).comp
          (Polynomial.C : ℚ →+* Polynomial ℚ))
        (LaurentPolynomial.C (Polynomial.X ^ 2) *
          (1 + LaurentPolynomial.T 2) ^ 2) S := by
  rw [gvcProfileHom, map_sum]
  rw [Polynomial.eval₂_eq_sum_range'
    ((LaurentPolynomial.C : Polynomial ℚ →+* PhaseLaurent).comp
      (Polynomial.C : ℚ →+* Polynomial ℚ))
      (p := S) (n := e + 1) (by omega)]
  apply Finset.sum_congr rfl
  intro j hj
  simp only [map_mul, map_pow,
    quadricPhaseRestriction_gvcT, quadricPhaseRestriction_gvcA,
    quadricPhaseRestriction_gvcRho, RingHom.coe_comp,
    Function.comp_apply, one_pow, mul_one]
  simp [quadricPhaseRestriction]

theorem phaseToLaurent_scaledProfile (S : ℚ[X]) :
    phaseToLaurent
        (scaleHeightByPhase (1 + Polynomial.X)
          (S.comp (Polynomial.X ^ 2))) =
      Polynomial.eval₂
        ((LaurentPolynomial.C : Polynomial ℚ →+* PhaseLaurent).comp
          (Polynomial.C : ℚ →+* Polynomial ℚ))
        (LaurentPolynomial.C (Polynomial.X ^ 2) *
          (1 + LaurentPolynomial.T 2) ^ 2) S := by
  rw [scaleHeightByPhase]
  rw [Polynomial.hom_eval₂]
  rw [Polynomial.eval₂_comp]
  simp only [Polynomial.map_add, Polynomial.map_one, Polynomial.map_X,
    map_mul, phaseToLaurent_C, map_add, map_one, phaseToLaurent_X,
    Polynomial.eval₂_X_pow, map_pow]
  have hcoeff :
      phaseToLaurent.comp (Polynomial.C.comp Polynomial.C) =
        (LaurentPolynomial.C : Polynomial ℚ →+* PhaseLaurent).comp
          (Polynomial.C : ℚ →+* Polynomial ℚ) := by
    ext a
    simp [phaseToLaurent]
  rw [hcoeff]
  have hbase :
      ((LaurentPolynomial.C Polynomial.X : PhaseLaurent) *
          (1 + LaurentPolynomial.T (2 : ℤ))) ^ 2 =
        LaurentPolynomial.C (Polynomial.X ^ 2) *
          (1 + LaurentPolynomial.T (2 : ℤ)) ^ 2 := by
    rw [mul_pow, map_pow]
  rw [hbase]
  have hmap :
      (LaurentPolynomial.C Polynomial.X : PhaseLaurent) ^ 2 =
        LaurentPolynomial.C (Polynomial.X ^ 2) := by
    rw [map_pow]
  rw [hmap]

theorem scaleHeightByPhase_pow
    (B f : Polynomial ℚ) (m : ℕ) :
    scaleHeightByPhase B (f ^ m) =
      scaleHeightByPhase B f ^ m := by
  change
    Polynomial.eval₂
        (Polynomial.C.comp Polynomial.C)
        (Polynomial.C Polynomial.X * B.map Polynomial.C) (f ^ m) =
      (Polynomial.eval₂
        (Polynomial.C.comp Polynomial.C)
        (Polynomial.C Polynomial.X * B.map Polynomial.C) f) ^ m
  rw [Polynomial.eval₂_pow]

theorem phaseToLaurent_profileHeightPhaseKernel
    (r m : ℕ) (S : ℚ[X]) :
    phaseToLaurent (profileHeightPhaseKernel r m S) =
      (1 + LaurentPolynomial.T 2) ^ (r * m) *
      (1 - LaurentPolynomial.C (Polynomial.X ^ 2) *
        (1 + LaurentPolynomial.T 2) ^ 2) ^ (2 * r * m) *
      (Polynomial.eval₂
        ((LaurentPolynomial.C : Polynomial ℚ →+* PhaseLaurent).comp
          (Polynomial.C : ℚ →+* Polynomial ℚ))
        (LaurentPolynomial.C (Polynomial.X ^ 2) *
          (1 + LaurentPolynomial.T 2) ^ 2) S) ^ m := by
  simp only [profileHeightPhaseKernel, map_mul, map_pow, map_add,
    map_one, map_sub, phaseToLaurent_X, phaseToLaurent_C,
    scaleHeightByPhase_pow]
  rw [phaseToLaurent_scaledProfile]
  rw [map_pow]

private theorem profilePowerRearrange
    {R : Type*} [CommMonoid R]
    (B E V T : R) (r m : ℕ) :
    (B ^ r * (T * E) ^ (2 * r) * V) ^ m =
      T ^ (2 * r * m) *
        (B ^ (r * m) * E ^ (2 * r * m) * V ^ m) := by
  simp only [mul_pow, ← pow_mul]
  ac_rfl

theorem quadricPhaseRestriction_gvcProfileP_pow
    (r e h m : ℕ) (S : ℚ[X]) (hS : S.natDegree ≤ e) :
    quadricPhaseRestriction (gvcProfileP r e h S ^ m) =
      LaurentPolynomial.T (-2 * ((r * m : ℕ) : ℤ)) *
        phaseToLaurent (profileHeightPhaseKernel r m S) := by
  rw [map_pow, gvcProfileP]
  simp only [map_mul, map_pow]
  rw [quadricPhaseRestriction_gvcRho,
    quadricPhaseRestriction_gvcA, quadricPhaseRestriction_gvcC,
    quadricPhaseRestriction_gvcProfileHom e S hS,
    phaseToLaurent_profileHeightPhaseKernel]
  simp only [one_pow, one_mul]
  rw [profilePowerRearrange]
  rw [LaurentPolynomial.T_pow]
  congr 1
  push_cast
  ring_nf

theorem quadricPhaseConstant_gvcProfileP_pow
    (r e h m : ℕ) (S : ℚ[X]) (hS : S.natDegree ≤ e) :
    quadricPhaseConstant (gvcProfileP r e h S ^ m) =
      (profileHeightPhaseKernel r m S).coeff (r * m) := by
  rw [quadricPhaseConstant,
    quadricPhaseRestriction_gvcProfileP_pow r e h m S hS,
    coeff_zero_T_neg_two_mul_phaseToLaurent]

theorem quadricPhaseConstant_gvcProfile_ladder
    (r e h m ell : ℕ) (S : ℚ[X])
    (hS : S.natDegree ≤ e) (hell : ell ≤ r * m) :
    quadricPhaseConstant
        (gvcX ^ (2 * ell) * gvcProfileP r e h S ^ m) =
      (profileHeightPhaseKernel r m S).coeff (r * m - ell) := by
  rw [quadricPhaseConstant, map_mul, map_pow,
    quadricPhaseRestriction_gvcX,
    quadricPhaseRestriction_gvcProfileP_pow r e h m S hS,
    LaurentPolynomial.T_pow]
  convert
    coeff_zero_T_two_mul_T_neg_two_mul_phaseToLaurent_general
      (r * m) ell (profileHeightPhaseKernel r m S) hell using 1
  all_goals norm_num

theorem scaleHeightByPhase_mul
    (B f g : Polynomial ℚ) :
    scaleHeightByPhase B (f * g) =
      scaleHeightByPhase B f * scaleHeightByPhase B g := by
  change
    Polynomial.eval₂
        (Polynomial.C.comp Polynomial.C)
        (Polynomial.C Polynomial.X * B.map Polynomial.C) (f * g) =
      Polynomial.eval₂
          (Polynomial.C.comp Polynomial.C)
          (Polynomial.C Polynomial.X * B.map Polynomial.C) f *
        Polynomial.eval₂
          (Polynomial.C.comp Polynomial.C)
          (Polynomial.C Polynomial.X * B.map Polynomial.C) g
  rw [Polynomial.eval₂_mul]

theorem scaleHeightByPhase_one_sub_sq_pow
    (B : Polynomial ℚ) (n : ℕ) :
    scaleHeightByPhase B ((1 - Polynomial.X ^ 2) ^ n) =
      (1 - Polynomial.C (Polynomial.X ^ 2) *
        (B.map Polynomial.C) ^ 2) ^ n := by
  rw [scaleHeightByPhase, Polynomial.eval₂_pow,
    Polynomial.eval₂_sub, Polynomial.eval₂_one,
    Polynomial.eval₂_pow, Polynomial.eval₂_X]
  rw [mul_pow, map_pow]

theorem profileHeightPhaseKernel_eq_scaleHeightByPhase
    (r m : ℕ) (S : ℚ[X]) :
    profileHeightPhaseKernel r m S =
      (((1 + Polynomial.X) ^ (r * m) : Polynomial ℚ).map
          Polynomial.C) *
        scaleHeightByPhase (1 + Polynomial.X)
          ((1 - Polynomial.X ^ 2) ^ (2 * r * m) *
            (S.comp (Polynomial.X ^ 2)) ^ m) := by
  rw [scaleHeightByPhase_mul, scaleHeightByPhase_one_sub_sq_pow,
    scaleHeightByPhase_pow]
  simp only [profileHeightPhaseKernel, Polynomial.map_pow,
    Polynomial.map_add, Polynomial.map_one, Polynomial.map_X]
  rw [scaleHeightByPhase_pow]
  ring

theorem gvcProfileEndpointDerivativeFactorization
    (r m : ℕ) (S : ℚ[X]) :
    (((1 - Polynomial.X ^ 2) ^ (2 * r * m) *
        (S.comp (Polynomial.X ^ 2)) ^ m).comp
      (1 + Polynomial.X)) =
      Polynomial.X ^ (2 * (r * m)) *
        gvcShiftedProfileFactor r m S := by
  rw [Polynomial.mul_comp, Polynomial.pow_comp, Polynomial.pow_comp,
    Polynomial.sub_comp, Polynomial.one_comp, Polynomial.pow_comp,
    Polynomial.X_comp]
  have hexp : 2 * r * m = 2 * (r * m) := Nat.mul_assoc 2 r m
  rw [hexp, endpoint_derivative_factorization (r * m),
    Polynomial.comp_assoc, Polynomial.pow_comp, Polynomial.X_comp]
  simp only [gvcShiftedProfileFactor, hexp]
  ring

theorem polynomialPrimitive_profile_comp_one_add_X
    (r m : ℕ) (S : ℚ[X]) :
    (polynomialPrimitive
        ((1 - Polynomial.X ^ 2) ^ (2 * r * m) *
          (S.comp (Polynomial.X ^ 2)) ^ m)).comp
        (1 + Polynomial.X) =
      shiftedProfilePrimitive (r * m) (gvcProfileMoment r m S)
        (gvcShiftedProfileFactor r m S) := by
  let F : Polynomial ℚ :=
    (polynomialPrimitive
      ((1 - Polynomial.X ^ 2) ^ (2 * r * m) *
        (S.comp (Polynomial.X ^ 2)) ^ m)).comp
      (1 + Polynomial.X)
  have hderivative :
      derivative F =
        derivative
          (shiftedProfilePrimitive (r * m) (gvcProfileMoment r m S)
            (gvcShiftedProfileFactor r m S)) := by
    dsimp [F]
    rw [Polynomial.derivative_comp, derivative_polynomialPrimitive,
      derivative_shiftedProfilePrimitive,
      gvcProfileEndpointDerivativeFactorization]
    simp
  have hconstant :
      (F - shiftedProfilePrimitive (r * m) (gvcProfileMoment r m S)
        (gvcShiftedProfileFactor r m S)).coeff 0 = 0 := by
    rw [Polynomial.coeff_sub, Polynomial.coeff_zero_eq_eval_zero,
      Polynomial.coeff_zero_eq_eval_zero]
    dsimp [F]
    rw [Polynomial.eval_comp]
    simp only [Polynomial.eval_add, Polynomial.eval_one,
      Polynomial.eval_X, add_zero, shiftedProfilePrimitive,
      Polynomial.eval_C]
    rw [gvcProfileMoment, formalIntegral01]
    rw [← Polynomial.coeff_zero_eq_eval_zero,
      polynomialPrimitive_coeff_zero]
    ring
  have hderivative_zero :
      derivative
          (F - shiftedProfilePrimitive (r * m) (gvcProfileMoment r m S)
            (gvcShiftedProfileFactor r m S)) = 0 := by
    rw [derivative_sub, hderivative, sub_self]
  have hconstant_poly :=
    Polynomial.eq_C_of_derivative_eq_zero hderivative_zero
  rw [hconstant, map_zero] at hconstant_poly
  exact sub_eq_zero.mp hconstant_poly

theorem integrateHeightCoefficients_profileHeightPhaseKernel
    (r m : ℕ) (S : ℚ[X]) (hr : 0 < r) (hm : 0 < m) :
    integrateHeightCoefficients (profileHeightPhaseKernel r m S) =
      gvcProfileEndpointKernel r m S := by
  rw [profileHeightPhaseKernel_eq_scaleHeightByPhase,
    integrateHeightCoefficients_map_mul]
  let core : Polynomial ℚ :=
    (1 - Polynomial.X ^ 2) ^ (2 * r * m) *
      (S.comp (Polynomial.X ^ 2)) ^ m
  have hscale := scaleHeightByPhase_integral
    (1 + Polynomial.X) core
  dsimp [core] at hscale
  rw [polynomialPrimitive_profile_comp_one_add_X] at hscale
  have hn : 0 < r * m := Nat.mul_pos hr hm
  have hnexp : r * m = (r * m - 1) + 1 := by omega
  have hpow :
      ((1 + Polynomial.X) ^ (r * m) : Polynomial ℚ) =
        (1 + Polynomial.X) ^ (r * m - 1) *
          (1 + Polynomial.X) := by
    calc
      ((1 + Polynomial.X) ^ (r * m) : Polynomial ℚ) =
          (1 + Polynomial.X) ^ ((r * m - 1) + 1) :=
        congrArg
          (fun d : ℕ => ((1 + Polynomial.X) ^ d : Polynomial ℚ))
          hnexp
      _ = (1 + Polynomial.X) ^ (r * m - 1) *
          (1 + Polynomial.X) := by rw [pow_succ]
  rw [hpow, mul_assoc, hscale, gvcProfileEndpointKernel,
    profileEndpointKernel]

/-- The phase bridge for Theorem 9.1 is constructed from the literal
multivariate family.  The degree hypothesis is exactly the manuscript's
declaration `S = sum_{j=0}^e s_j z^j`. -/
theorem verifiedProfileFamilyBridge
    (r e h : ℕ) (S : ℚ[X]) (hr : 0 < r)
    (hS : S.natDegree ≤ e) :
    ProfileFamilyBridge r e h S where
  pure_phase_eq m hm := by
    rw [algebraicReynoldsMoment_eq_formalIntegral_quadricPhaseConstant
        (gvcProfileOrder r e h * m)
        (gvcProfileP_pow_isHomogeneous r e h m S),
      quadricPhaseConstant_gvcProfileP_pow r e h m S hS,
      ← integrateHeightCoefficients_coeff,
      integrateHeightCoefficients_profileHeightPhaseKernel r m S hr hm]
  ladder_phase_eq m ell hm _hell helm := by
    rw [algebraicReynoldsMoment_eq_formalIntegral_quadricPhaseConstant
        (gvcProfileOrder r e h * m + ell)
        (gvcProfileMultiplier_isHomogeneous r e h m ell S),
      quadricPhaseConstant_gvcProfile_ladder r e h m ell S hS helm,
      ← integrateHeightCoefficients_coeff,
      integrateHeightCoefficients_profileHeightPhaseKernel r m S hr hm]

theorem verified_gvcProfile_pure_identity
    {r e h : ℕ} {S : ℚ[X]} (hr : 0 < r) (hS : S.natDegree ≤ e)
    (m : ℕ) (hm : 0 < m) :
    differentialAction (gvcProfileLambda r e h ^ m)
      (gvcProfileP r e h S ^ m) = 0 :=
  gvcProfile_pure_identity
    (verifiedProfileFamilyBridge r e h S hr hS) hr m hm

theorem verified_gvcProfile_exact_ladder
    {r e h : ℕ} {S : ℚ[X]} (hr : 0 < r) (hS : S.natDegree ≤ e)
    (m ell : ℕ) (hm : 0 < m) (hell : 1 ≤ ell)
    (helm : ell ≤ r * m) :
    differentialAction
        (gvcDelta ^ (gvcProfileOrder r e h * m + ell))
        (gvcX ^ (2 * ell) * gvcProfileP r e h S ^ m) =
      MvPolynomial.C
        (reynoldsScale (gvcProfileOrder r e h * m + ell) *
        (Nat.choose (r * m - 1) (ell - 1) : ℚ) *
        gvcProfileMoment r m S) :=
  gvcProfile_exact_ladder
    (verifiedProfileFamilyBridge r e h S hr hS) m ell hm hell helm

theorem verified_gvcProfile_not_generalizedVanishingFor
    {r e h : ℕ} {S : ℚ[X]} (hr : 0 < r) (hS : S.natDegree ≤ e)
    (hmoment : ∀ m, 0 < m → gvcProfileMoment r m S ≠ 0) :
    ¬ GeneralizedVanishingFor
      (gvcProfileLambda r e h) (gvcProfileP r e h S) :=
  gvcProfile_not_generalizedVanishingFor
    (verifiedProfileFamilyBridge r e h S hr hS) hr hmoment

/-- The full profile family after coefficient base change to an arbitrary
characteristic-zero field. -/
theorem verified_gvcProfile_charZero_not_generalizedVanishingFor
    {K : Type*} [Field K] [CharZero K]
    {r e h : ℕ} {S : ℚ[X]} (hr : 0 < r) (hS : S.natDegree ≤ e)
    (hmoment : ∀ m, 0 < m → gvcProfileMoment r m S ≠ 0) :
    ¬ GeneralizedVanishingFor
      (map (Rat.castHom K) (gvcProfileLambda r e h))
      (map (Rat.castHom K) (gvcProfileP r e h S)) := by
  let B := verifiedProfileFamilyBridge r e h S hr hS
  apply not_generalizedVanishingFor_map (Rat.castHom K)
    (RingHom.injective (Rat.castHom K))
    (gvcProfile_purePowersVanish B hr)
  · intro m hm
    exact gvcProfile_mixed_ne_zero B hr m hm (hmoment m hm)

end GVC
