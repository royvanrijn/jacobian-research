/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GVC.ConcreteWitness

/-!
# The winding--profile--radial family

This file defines the literal multivariate family from Theorem 9.1 and
checks its homogeneity and the claimed degree/order formula.  Its all-order
quadric phase extraction is constructed in `GVC.ProfilePhase`.
-/

namespace GVC

open MvPolynomial Polynomial

/-- The homogeneous lift of a univariate profile of declared degree at most
`e`.  Coefficients above `e` are intentionally ignored, matching the paper's
displayed finite sum. -/
noncomputable def gvcProfileHom (e : ℕ) (S : ℚ[X]) : TernaryPolynomial :=
  ∑ j ∈ Finset.range (e + 1),
    MvPolynomial.C (S.coeff j) * (gvcT ^ 2 * gvcA ^ 2) ^ j *
      gvcRho ^ (3 * (e - j))

theorem gvcProfileHom_isHomogeneous (e : ℕ) (S : ℚ[X]) :
    (gvcProfileHom e S).IsHomogeneous (6 * e) := by
  have ht : (gvcT ^ 2).IsHomogeneous 2 := by
    simpa [gvcT] using
      (MvPolynomial.isHomogeneous_X ℚ (2 : Fin 3)).pow 2
  have hbase : (gvcT ^ 2 * gvcA ^ 2).IsHomogeneous 6 := by
    convert ht.mul (gvcA_isHomogeneous.pow 2) using 1
  apply IsHomogeneous.sum
  intro j hj
  simp only [Finset.mem_range] at hj
  have hterm := (hbase.pow j).mul
    (gvcRho_isHomogeneous.pow (3 * (e - j)))
  have hterm' :
      ((gvcT ^ 2 * gvcA ^ 2) ^ j *
        gvcRho ^ (3 * (e - j))).IsHomogeneous (6 * e) := by
    convert hterm using 1
    omega
  simpa only [mul_assoc] using hterm'.C_mul (S.coeff j)

/-- The integer `N = 6r + 3e + h` controlling both polynomial degree and
Laplacian power in the family. -/
def gvcProfileOrder (r e h : ℕ) : ℕ := 6 * r + 3 * e + h

/-- The manuscript's literal family
`rho^h * A^r * C^(2r) * S^hom`. -/
noncomputable def gvcProfileP
    (r e h : ℕ) (S : ℚ[X]) : TernaryPolynomial :=
  ((gvcRho ^ h * gvcA ^ r) * gvcC ^ (2 * r)) * gvcProfileHom e S

/-- The corresponding homogeneous differential symbol `Delta^N`. -/
noncomputable def gvcProfileLambda
    (r e h : ℕ) : TernaryPolynomial :=
  gvcDelta ^ gvcProfileOrder r e h

theorem gvcProfileP_isHomogeneous (r e h : ℕ) (S : ℚ[X]) :
    (gvcProfileP r e h S).IsHomogeneous
      (2 * gvcProfileOrder r e h) := by
  have hpoly :=
    (((gvcRho_isHomogeneous.pow h).mul (gvcA_isHomogeneous.pow r)).mul
      (gvcC_isHomogeneous.pow (2 * r))).mul
      (gvcProfileHom_isHomogeneous e S)
  convert hpoly using 1
  · rfl
  · simp [gvcProfileOrder]
    omega

theorem gvcProfileLambda_isHomogeneous (r e h : ℕ) :
    (gvcProfileLambda r e h).IsHomogeneous
      (2 * gvcProfileOrder r e h) := by
  simpa [gvcProfileLambda] using
    gvcDelta_isHomogeneous.pow (gvcProfileOrder r e h)

@[simp] theorem gvcProfileHom_zero_one :
    gvcProfileHom 0 1 = 1 := by
  simp [gvcProfileHom]

/-- The concrete degree-twelve witness is the minimal profile
specialization `r=1, e=h=0, S=1`. -/
theorem gvcProfileP_minimal :
    gvcProfileP 1 0 0 1 = gvcP := by
  simp [gvcProfileP, gvcP]

/-- The paper's radial specialization is obtained by varying `h` while
keeping the minimal winding and constant profile. -/
theorem gvcProfileP_radial (h : ℕ) :
    gvcProfileP 1 0 h 1 = gvcRho ^ h * gvcP := by
  simp [gvcProfileP, gvcP, mul_assoc]

theorem gvcProfileLambda_radial (h : ℕ) :
    gvcProfileLambda 1 0 h = gvcDelta ^ (6 + h) := by
  simp [gvcProfileLambda, gvcProfileOrder]

/-- The paper's algebraic profile moment
`∫₀¹ (1-v²)^(2rm) S(v²)^m dv`. -/
noncomputable def gvcProfileMoment
    (r m : ℕ) (S : ℚ[X]) : ℚ :=
  formalIntegral01
    ((1 - Polynomial.X ^ 2) ^ (2 * r * m) *
      (S.comp (Polynomial.X ^ 2)) ^ m)

/-- The factor left after extracting the endpoint zero of order `2rm` from
the derivative of the shifted primitive. -/
noncomputable def gvcShiftedProfileFactor
    (r m : ℕ) (S : ℚ[X]) : ℚ[X] :=
  (2 + Polynomial.X) ^ (2 * r * m) *
    (S.comp ((1 + Polynomial.X) ^ 2)) ^ m

/-- The one-variable endpoint kernel appearing after phase extraction for
the full family. -/
noncomputable def gvcProfileEndpointKernel
    (r m : ℕ) (S : ℚ[X]) : ℚ[X] :=
  profileEndpointKernel (r * m) (gvcProfileMoment r m S)
    (gvcShiftedProfileFactor r m S)

theorem gvcProfileMoment_minimal (m : ℕ) :
    gvcProfileMoment 1 m 1 = cuspMoment m := by
  rw [gvcProfileMoment]
  simp [cuspMoment_eq_formalIntegral]

theorem gvcProfileEndpointKernel_minimal (m : ℕ) :
    gvcProfileEndpointKernel 1 m 1 =
      endpointKernel m (cuspMoment m) (endpointPrimitiveTail m) := by
  rw [gvcProfileEndpointKernel, gvcProfileMoment_minimal,
    profileEndpointKernel, gvcShiftedProfileFactor,
    endpointKernel_actual_eq]
  simp [shiftedProfilePrimitive, shiftedEndpointPrimitive]

theorem gvcProfileLambda_pow_eq_delta_pow
    (r e h m : ℕ) :
    gvcProfileLambda r e h ^ m =
      gvcDelta ^ (gvcProfileOrder r e h * m) := by
  rw [gvcProfileLambda, pow_mul]

theorem gvcProfileP_pow_isHomogeneous
    (r e h m : ℕ) (S : ℚ[X]) :
    (gvcProfileP r e h S ^ m).IsHomogeneous
      (2 * (gvcProfileOrder r e h * m)) := by
  convert (gvcProfileP_isHomogeneous r e h S).pow m using 1
  simp [mul_assoc]

theorem gvcProfileMultiplier_isHomogeneous
    (r e h m ell : ℕ) (S : ℚ[X]) :
    (gvcX ^ (2 * ell) * gvcProfileP r e h S ^ m).IsHomogeneous
      (2 * (gvcProfileOrder r e h * m + ell)) := by
  have hx : (gvcX ^ (2 * ell)).IsHomogeneous (2 * ell) := by
    simpa [gvcX] using
      (MvPolynomial.isHomogeneous_X ℚ (0 : Fin 3)).pow (2 * ell)
  convert hx.mul (gvcProfileP_pow_isHomogeneous r e h m S) using 1
  omega

/-- Exact interface for the pure and shifted quadric phase-extraction
formulas.  `GVC.ProfilePhase` constructs it from the literal family; the
endpoint ladder, Reynolds contraction, and all subsequent implications are
proved below. -/
structure ProfileFamilyBridge
    (r e h : ℕ) (S : ℚ[X]) where
  pure_phase_eq : ∀ m, 0 < m →
    algebraicReynoldsMoment (gvcProfileOrder r e h * m)
        (gvcProfileP r e h S ^ m) =
      (gvcProfileEndpointKernel r m S).coeff (r * m)
  ladder_phase_eq : ∀ m ell, 0 < m → 1 ≤ ell → ell ≤ r * m →
    algebraicReynoldsMoment (gvcProfileOrder r e h * m + ell)
        (gvcX ^ (2 * ell) * gvcProfileP r e h S ^ m) =
      (gvcProfileEndpointKernel r m S).coeff (r * m - ell)

theorem gvcProfile_pure_identity
    {r e h : ℕ} {S : ℚ[X]} (B : ProfileFamilyBridge r e h S)
    (hr : 0 < r) (m : ℕ) (hm : 0 < m) :
    differentialAction (gvcProfileLambda r e h ^ m)
      (gvcProfileP r e h S ^ m) = 0 := by
  rw [gvcProfileLambda_pow_eq_delta_pow,
    differentialAction_delta_pow_eq_reynolds _ _
      (gvcProfileP_pow_isHomogeneous r e h m S),
    B.pure_phase_eq m hm, gvcProfileEndpointKernel,
    profileEndpointKernel_coeff_pure (r * m) (by positivity)]
  simp

/-- The complete exact multiplier ladder of Theorem 9.1, conditional only
on its displayed phase extraction. -/
theorem gvcProfile_exact_ladder
    {r e h : ℕ} {S : ℚ[X]} (B : ProfileFamilyBridge r e h S)
    (m ell : ℕ) (hm : 0 < m) (hell : 1 ≤ ell) (helm : ell ≤ r * m) :
    differentialAction
        (gvcDelta ^ (gvcProfileOrder r e h * m + ell))
        (gvcX ^ (2 * ell) * gvcProfileP r e h S ^ m) =
      MvPolynomial.C
        (reynoldsScale (gvcProfileOrder r e h * m + ell) *
        (Nat.choose (r * m - 1) (ell - 1) : ℚ) *
        gvcProfileMoment r m S) := by
  rw [differentialAction_delta_pow_eq_reynolds _ _
      (gvcProfileMultiplier_isHomogeneous r e h m ell S),
    B.ladder_phase_eq m ell hm hell helm, gvcProfileEndpointKernel,
    profileEndpointKernel_coeff_ladder (r * m) ell hell helm]
  simp only [Nat.mul_comm]
  ring_nf

theorem gvcProfile_exact_ladder_ne_zero
    {r e h : ℕ} {S : ℚ[X]} (B : ProfileFamilyBridge r e h S)
    (m ell : ℕ) (hm : 0 < m) (hell : 1 ≤ ell) (helm : ell ≤ r * m)
    (hmoment : gvcProfileMoment r m S ≠ 0) :
    differentialAction
        (gvcDelta ^ (gvcProfileOrder r e h * m + ell))
        (gvcX ^ (2 * ell) * gvcProfileP r e h S ^ m) ≠ 0 := by
  rw [gvcProfile_exact_ladder B m ell hm hell helm,
    MvPolynomial.C_ne_zero]
  exact mul_ne_zero
    (mul_ne_zero (reynoldsScale_ne_zero _) (by
      exact_mod_cast Nat.choose_ne_zero (by omega))) hmoment

theorem gvcProfile_mixed_ne_zero
    {r e h : ℕ} {S : ℚ[X]} (B : ProfileFamilyBridge r e h S)
    (hr : 0 < r) (m : ℕ) (hm : 0 < m)
    (hmoment : gvcProfileMoment r m S ≠ 0) :
    differentialAction (gvcProfileLambda r e h ^ m)
      (gvcQ * gvcProfileP r e h S ^ m) ≠ 0 := by
  intro hzero
  have hnext :
      differentialAction gvcDelta
          (differentialAction (gvcProfileLambda r e h ^ m)
            (gvcQ * gvcProfileP r e h S ^ m)) =
        differentialAction
          (gvcDelta ^ (gvcProfileOrder r e h * m + 1))
          (gvcX ^ 2 * gvcProfileP r e h S ^ m) := by
    rw [gvcQ, gvcProfileLambda_pow_eq_delta_pow,
      ← differentialAction_mul_left, ← pow_succ']
  rw [hzero, differentialAction_zero_right] at hnext
  exact gvcProfile_exact_ladder_ne_zero B m 1 hm (by omega)
    (by have := Nat.mul_pos hr hm; omega) hmoment hnext.symm

theorem gvcProfile_purePowersVanish
    {r e h : ℕ} {S : ℚ[X]} (B : ProfileFamilyBridge r e h S)
    (hr : 0 < r) :
    PurePowersVanish (gvcProfileLambda r e h) (gvcProfileP r e h S) := by
  intro m hm
  exact gvcProfile_pure_identity B hr m hm

/-- Once the phase bridge and the paper's nonzero-moment hypothesis are
supplied, the literal family is a GVC counterexample.  The bridge is
constructed in `GVC.ProfilePhase`. -/
theorem gvcProfile_not_generalizedVanishingFor
    {r e h : ℕ} {S : ℚ[X]} (B : ProfileFamilyBridge r e h S)
    (hr : 0 < r) (hmoment : ∀ m, 0 < m → gvcProfileMoment r m S ≠ 0) :
    ¬ GeneralizedVanishingFor
      (gvcProfileLambda r e h) (gvcProfileP r e h S) := by
  intro hGVC
  obtain ⟨M, hM⟩ := hGVC (gvcProfile_purePowersVanish B hr) gvcQ
  let m := max M 1
  have hmM : M ≤ m := le_max_left _ _
  have hm : 0 < m := lt_of_lt_of_le Nat.zero_lt_one (le_max_right _ _)
  exact gvcProfile_mixed_ne_zero B hr m hm (hmoment m hm) (hM m hmM)

end GVC
