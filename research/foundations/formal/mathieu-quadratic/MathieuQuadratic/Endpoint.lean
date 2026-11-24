import Mathlib.Algebra.Algebra.Rat
import Mathlib.Analysis.Calculus.Deriv.Polynomial
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus
import Mathlib.RingTheory.PowerSeries.Binomial
import Mathlib.RingTheory.Polynomial.Pochhammer
import Mathlib.Tactic

/-!
# The polynomial endpoint identity

The beta antiderivative is rational. Its expansion at one agrees with its
constant value through the exponent of the beta density. The generalized binomial
power series allows an arbitrary power in the polynomial prefactor.
-/

open Polynomial
open scoped BigOperators Interval PowerSeries

namespace MathieuQuadratic

/-- The coefficients of the polynomial beta antiderivative after removing
its factor of degree M. -/
def powerBetaCoefficient (M k j : ℕ) : ℚ :=
  (M : ℚ) / (M + 2 * j : ℕ) * (k.choose j : ℚ) * (-1) ^ j

noncomputable def powerBetaQuotient (M k : ℕ) : ℚ[X] :=
  ∑ j ∈ Finset.range (k + 1), C (powerBetaCoefficient M k j) * X ^ (2 * j)

noncomputable def powerBetaAntiderivative (M k : ℕ) : ℚ[X] :=
  X ^ M * powerBetaQuotient M k

noncomputable def powerBetaMoment (M k : ℕ) : ℚ :=
  (powerBetaAntiderivative M k).eval 1

/-- The polynomial before the final Laurent coefficient extraction. -/
noncomputable def powerEndpointPolynomial (M v k : ℕ) : ℚ[X] :=
  (X + 1) ^ v * (powerBetaQuotient M k).comp (X + 1)

theorem powerBetaAntiderivative_derivative (M k : ℕ) (hM : 0 < M) :
    (powerBetaAntiderivative M k).derivative = C (M : ℚ) * X ^ (M - 1) * (1 - X ^ 2) ^ k := by
  rw [powerBetaAntiderivative, powerBetaQuotient, Finset.mul_sum, derivative_sum,
    sub_eq_add_neg, add_comm 1, add_pow, Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  rw [mul_left_comm (X ^ M), ← pow_add, derivative_C_mul_X_pow,
    one_pow, mul_one, neg_pow, ← pow_mul]
  have hn : (↑(M + 2 * j) : ℚ) ≠ 0 := by positivity
  have hcoeff : powerBetaCoefficient M k j * (↑(M + 2 * j) : ℚ) =
      M * (k.choose j : ℚ) * (-1) ^ j := by
    unfold powerBetaCoefficient
    field_simp
  rw [hcoeff, map_mul, map_mul, map_pow, map_neg, map_one, C_eq_natCast, C_eq_natCast]
  rw [show M + 2 * j - 1 = (M - 1) + 2 * j by omega, pow_add]
  ring

private def betaDensity (M n : ℕ) (x : ℝ) : ℝ :=
  (M : ℝ) * x ^ (M - 1) * (1 - x ^ 2) ^ n

private theorem betaDensity_continuous (M n : ℕ) : Continuous (betaDensity M n) := by
  unfold betaDensity
  fun_prop

theorem powerBetaMoment_cast_eq_integral (M k : ℕ) (hM : 0 < M) :
    (powerBetaMoment M k : ℝ) =
      ∫ x in (0 : ℝ)..1, (M : ℝ) * x ^ (M - 1) * (1 - x ^ 2) ^ k := by
  let p : ℝ[X] := (powerBetaAntiderivative M k).map (algebraMap ℚ ℝ)
  have hd (x : ℝ) : HasDerivAt (fun t ↦ p.eval t) (betaDensity M k x) x := by
    have he : p.derivative.eval x = betaDensity M k x := by
      simp [p, betaDensity, derivative_map, powerBetaAntiderivative_derivative _ _ hM]
    simpa only [he] using p.hasDerivAt x
  change _ = ∫ x in (0 : ℝ)..1, betaDensity M k x
  rw [intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ ↦ hd x)
    ((betaDensity_continuous M k).intervalIntegrable 0 1)]
  simp [p, eval_map, powerBetaMoment, powerBetaAntiderivative, hM.ne']

theorem powerBetaMoment_pos (M k : ℕ) (hM : 0 < M) : 0 < powerBetaMoment M k := by
  have hpos : 0 < ∫ x in (0 : ℝ)..1, betaDensity M k x := by
    apply intervalIntegral.intervalIntegral_pos_of_pos_on
      ((betaDensity_continuous M k).intervalIntegrable 0 1) _ zero_lt_one
    intro x hx
    have hx0 := hx.1
    have hx1 := hx.2
    dsimp [betaDensity]
    apply mul_pos (mul_pos (by exact_mod_cast hM) (pow_pos hx0 _))
    exact pow_pos (sub_pos.mpr (by nlinarith)) _
  change 0 < ∫ x in (0 : ℝ)..1, (M : ℝ) * x ^ (M - 1) * (1 - x ^ 2) ^ k at hpos
  rw [← powerBetaMoment_cast_eq_integral M k hM] at hpos
  exact_mod_cast hpos

private noncomputable def betaIntegral (M n : ℕ) : ℝ :=
  ∫ x in (0 : ℝ)..1, betaDensity M n x

private theorem betaIntegral_zero (M : ℕ) (hM : 0 < M) : betaIntegral M 0 = 1 := by
  have hd (x : ℝ) : HasDerivAt (fun t : ℝ ↦ t ^ M) (betaDensity M 0 x) x := by
    simpa [betaDensity, derivative_pow] using ((X : ℝ[X]) ^ M).hasDerivAt x
  rw [betaIntegral, intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ ↦ hd x)
    ((betaDensity_continuous M 0).intervalIntegrable 0 1)]
  simp [hM.ne']

private theorem betaIntegral_succ (M n : ℕ) (hM : 0 < M) :
    ((M : ℝ) + 2 * (n + 1 : ℕ)) * betaIntegral M (n + 1) =
      (2 * (n + 1 : ℕ) : ℝ) * betaIntegral M n := by
  let p : ℝ[X] := C (M : ℝ) * X ^ M * (1 - X ^ 2) ^ (n + 1)
  let c : ℝ := (M : ℝ) + 2 * (n + 1 : ℕ)
  let d : ℝ := 2 * (n + 1 : ℕ)
  have hd (x : ℝ) : HasDerivAt (fun t ↦ p.eval t)
      (c * betaDensity M (n + 1) x - d * betaDensity M n x) x := by
    have he : p.derivative.eval x = c * betaDensity M (n + 1) x - d * betaDensity M n x := by
      simp only [p, c, d, betaDensity, derivative_mul, derivative_C, derivative_pow, derivative_sub, derivative_one, derivative_X, zero_sub,
        eval_add, eval_sub, eval_mul, eval_neg, eval_C, eval_pow, eval_X, eval_one, eval_zero,
        show 2 - 1 = 1 by decide, pow_one, mul_one, Nat.add_sub_cancel]
      have powM : x ^ M = x ^ (M - 1) * x := by rw [← pow_succ, Nat.sub_add_cancel hM]
      rw [powM, pow_succ (1 - x ^ 2)]
      ring
    simpa only [he] using p.hasDerivAt x
  have hc := ((betaDensity_continuous M (n + 1)).const_mul c).intervalIntegrable
    (μ := MeasureTheory.volume) 0 1
  have hi := ((betaDensity_continuous M n).const_mul d).intervalIntegrable
    (μ := MeasureTheory.volume) 0 1
  have he := intervalIntegral.integral_eq_sub_of_hasDerivAt (fun x _ ↦ hd x) (hc.sub hi)
  rw [intervalIntegral.integral_sub hc hi, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_const_mul] at he
  have hv : p.eval 1 - p.eval 0 = 0 := by
    simp [p, hM.ne']
  rw [hv, sub_eq_zero] at he
  exact he

open Polynomial in
/-- The factorial and rising-factorial evaluation of the beta moment. -/
theorem powerBetaMoment_closed_form (M k : ℕ) (hM : 0 < M) :
    powerBetaMoment M k =
      (k.factorial : ℚ) / (ascPochhammer ℚ k).eval (1 + (M : ℚ) / 2) := by
  have hc (n : ℕ) : betaIntegral M n =
      (n.factorial : ℝ) / (ascPochhammer ℝ n).eval (1 + (M : ℝ) / 2) := by
    induction n with
    | zero => simp [betaIntegral_zero M hM]
    | succ n ih =>
        have hpos : 0 < (ascPochhammer ℝ n).eval (1 + (M : ℝ) / 2) :=
          ascPochhammer_pos _ _ (by positivity)
        have h := betaIntegral_succ M n hM
        rw [ih] at h
        rw [Nat.factorial_succ, Nat.cast_mul, ascPochhammer_succ_eval]
        apply (eq_div_iff (mul_ne_zero hpos.ne' (by positivity))).2
        simp only [Nat.cast_succ] at h ⊢
        let D := (ascPochhammer ℝ n).eval (1 + (M : ℝ) / 2)
        change betaIntegral M (n + 1) * (D * (1 + (M : ℝ) / 2 + n)) =
          ((n : ℝ) + 1) * n.factorial
        calc
          _ = (((M : ℝ) + 2 * ((n : ℝ) + 1)) * betaIntegral M (n + 1)) * D / 2 := by ring
          _ = (2 * ((n : ℝ) + 1) * n.factorial / D) * D / 2 := by rw [h]; ring
          _ = _ := by rw [div_mul_cancel₀ _ hpos.ne']; ring
  have he := hc k
  have hb : (powerBetaMoment M k : ℝ) = betaIntegral M k := powerBetaMoment_cast_eq_integral M k hM
  rw [← hb] at he
  apply (RingHom.injective (algebraMap ℚ ℝ))
  change (powerBetaMoment M k : ℝ) = (algebraMap ℚ ℝ) _
  rw [he, map_div₀, map_natCast, ← Polynomial.eval₂_at_apply, ← ascPochhammer_eval₂]
  simp

private theorem X_pow_succ_dvd_of_derivative {p : ℚ[X]} {n : ℕ} (h0 : p.coeff 0 = 0)
    (h : X ^ n ∣ p.derivative) : X ^ (n + 1) ∣ p := by
  rw [Polynomial.X_pow_dvd_iff] at h ⊢
  intro d hd
  cases d with
  | zero => exact h0
  | succ j =>
      have hj := h j (Nat.lt_of_succ_lt_succ hd)
      rw [coeff_derivative] at hj
      exact (mul_eq_zero.mp hj).resolve_right (by positivity)

theorem powerBetaAntiderivative_comp_sub_dvd (M k : ℕ) (hM : 0 < M) :
    X ^ (k + 1) ∣ (powerBetaAntiderivative M k).comp (X + 1) - C (powerBetaMoment M k) := by
  apply X_pow_succ_dvd_of_derivative
  · rw [coeff_zero_eq_eval_zero]
    simp [powerBetaMoment, eval_comp]
  · rw [derivative_sub, derivative_C, sub_zero, derivative_comp,
      show derivative (X + 1 : ℚ[X]) = 1 by simp, one_mul,
      powerBetaAntiderivative_derivative _ _ hM]
    simp only [mul_comp, pow_comp, C_comp, X_comp, sub_comp, one_comp]
    rw [show (1 : ℚ[X]) - (X + 1) ^ 2 = X * (-(X + 2)) by ring, mul_pow, ← mul_assoc]
    exact dvd_mul_left _ _ |>.trans (dvd_mul_right _ _)

private theorem powerEndpointPolynomial_series (M v k : ℕ) :
    (powerEndpointPolynomial M v k : ℚ⟦X⟧) =
      PowerSeries.binomialSeries ℚ ((v : ℤ) - M) *
        ((powerBetaAntiderivative M k).comp (X + 1) : ℚ⟦X⟧) := by
  rw [powerEndpointPolynomial, powerBetaAntiderivative, mul_comp, pow_comp, X_comp]
  push_cast
  have hb (n : ℕ) : (PowerSeries.X + 1 : ℚ⟦X⟧) ^ n =
      PowerSeries.binomialSeries ℚ (n : ℤ) := by
    rw [PowerSeries.binomialSeries_nat, add_comm]
  simp only [← mul_assoc, hb, ← PowerSeries.binomialSeries_add, sub_add_cancel]

/-- Coefficient agreement through degree k with the generalized binomial. -/
theorem powerEndpointPolynomial_coeff (M v k d : ℕ) (hM : 0 < M) (hd : d ≤ k) :
    (powerEndpointPolynomial M v k).coeff d =
      powerBetaMoment M k * ((Ring.choose ((v : ℤ) - M) d : ℤ) : ℚ) := by
  obtain ⟨error, he⟩ := powerBetaAntiderivative_comp_sub_dvd M k hM
  have hs : (powerBetaAntiderivative M k).comp (X + 1) =
      C (powerBetaMoment M k) + X ^ (k + 1) * error := by
    rw [← he]; ring
  rw [← Polynomial.coeff_coe, powerEndpointPolynomial_series, hs]
  push_cast
  rw [mul_add, map_add, mul_left_comm _ (PowerSeries.X ^ (k + 1)), PowerSeries.coeff_X_pow_mul']
  rw [if_neg (by omega : ¬ k + 1 ≤ d), add_zero,
    PowerSeries.coeff_mul_C, PowerSeries.binomialSeries_coeff, zsmul_eq_mul, mul_one, mul_comm]

private theorem powerBetaAntiderivative_first_error (M k : ℕ) (hM : 0 < M) :
    ((powerBetaAntiderivative M k).comp (X + 1) - C (powerBetaMoment M k)).coeff (k + 1) =
      (M : ℚ) * (-2) ^ k / (k + 1 : ℕ) := by
  let p := (powerBetaAntiderivative M k).comp (X + 1) - C (powerBetaMoment M k)
  have hp : derivative p = X ^ k * (C (M : ℚ) * (X + 1) ^ (M - 1) * (-(X + 2)) ^ k) := by
    dsimp only [p]
    rw [derivative_sub, derivative_C, sub_zero, derivative_comp,
      show derivative (X + 1 : ℚ[X]) = 1 by simp, one_mul,
      powerBetaAntiderivative_derivative _ _ hM]
    simp only [mul_comp, pow_comp, C_comp, X_comp, sub_comp, one_comp]
    rw [show (1 : ℚ[X]) - (X + 1) ^ 2 = X * (-(X + 2)) by ring, mul_pow]
    ring
  have hc : p.derivative.coeff k = (M : ℚ) * (-2) ^ k := by
    rw [hp, coeff_X_pow_mul', if_pos le_rfl, Nat.sub_self, coeff_zero_eq_eval_zero]
    simp
  rw [coeff_derivative] at hc
  exact (eq_div_iff (by positivity)).2 (by simpa only [p, Nat.cast_succ] using hc)

/-- The first coefficient beyond the range of generalized binomial agreement. -/
theorem powerEndpointPolynomial_coeff_succ (M v k : ℕ) (hM : 0 < M) :
    (powerEndpointPolynomial M v k).coeff (k + 1) =
      powerBetaMoment M k * ((Ring.choose ((v : ℤ) - M) (k + 1) : ℤ) : ℚ) +
        (M : ℚ) * (-2) ^ k / (k + 1 : ℕ) := by
  obtain ⟨error, he⟩ := powerBetaAntiderivative_comp_sub_dvd M k hM
  have herr := powerBetaAntiderivative_first_error M k hM
  rw [he, coeff_X_pow_mul', if_pos le_rfl, Nat.sub_self] at herr
  have hs : (powerBetaAntiderivative M k).comp (X + 1) =
      C (powerBetaMoment M k) + X ^ (k + 1) * error := by
    rw [← he]; ring
  rw [← Polynomial.coeff_coe, powerEndpointPolynomial_series, hs]
  push_cast
  rw [mul_add, map_add, mul_left_comm _ (PowerSeries.X ^ (k + 1)), PowerSeries.coeff_X_pow_mul',
    if_pos le_rfl, Nat.sub_self, PowerSeries.coeff_mul_C, PowerSeries.binomialSeries_coeff,
    zsmul_eq_mul, mul_one, mul_comm _ (powerBetaMoment M k), PowerSeries.coeff_zero_eq_constantCoeff,
    map_mul, PowerSeries.binomialSeries_constantCoeff, one_mul,
    ← PowerSeries.coeff_zero_eq_constantCoeff, Polynomial.coeff_coe, herr]
  simp only [Nat.cast_succ]

/-! The squared normalizer uses exponent 2k and prefactor degree k. -/

abbrev betaCoefficient (M k j : ℕ) : ℚ := powerBetaCoefficient M (2 * k) j
noncomputable abbrev betaQuotient (M k : ℕ) : ℚ[X] := powerBetaQuotient M (2 * k)
noncomputable abbrev betaAntiderivative (M k : ℕ) : ℚ[X] := powerBetaAntiderivative M (2 * k)
noncomputable abbrev betaMoment (M k : ℕ) : ℚ := powerBetaMoment M (2 * k)
noncomputable abbrev endpointPolynomial (M k : ℕ) : ℚ[X] := powerEndpointPolynomial M k (2 * k)

theorem betaAntiderivative_derivative (M k : ℕ) (hM : 0 < M) :
    (betaAntiderivative M k).derivative = C (M : ℚ) * X ^ (M - 1) * (1 - X ^ 2) ^ (2 * k) :=
  powerBetaAntiderivative_derivative M (2 * k) hM

theorem betaAntiderivative_comp_sub_dvd (M k : ℕ) (hM : 0 < M) :
    X ^ (2 * k + 1) ∣ (betaAntiderivative M k).comp (X + 1) - C (betaMoment M k) :=
  powerBetaAntiderivative_comp_sub_dvd M (2 * k) hM

theorem betaMoment_cast_eq_integral (M k : ℕ) (hM : 0 < M) :
    (betaMoment M k : ℝ) =
      ∫ x in (0 : ℝ)..1, (M : ℝ) * x ^ (M - 1) * (1 - x ^ 2) ^ (2 * k) :=
  powerBetaMoment_cast_eq_integral M (2 * k) hM

theorem betaMoment_pos (M k : ℕ) (hM : 0 < M) : 0 < betaMoment M k :=
  powerBetaMoment_pos M (2 * k) hM

theorem betaMoment_closed_form (M k : ℕ) (hM : 0 < M) :
    betaMoment M k =
      ((2 * k).factorial : ℚ) / (ascPochhammer ℚ (2 * k)).eval (1 + (M : ℚ) / 2) :=
  powerBetaMoment_closed_form M (2 * k) hM

theorem endpointPolynomial_coeff (M k d : ℕ) (hM : 0 < M) (hd : d ≤ 2 * k) :
    (endpointPolynomial M k).coeff d =
      betaMoment M k * ((Ring.choose ((k : ℤ) - M) d : ℤ) : ℚ) :=
  powerEndpointPolynomial_coeff M k (2 * k) d hM hd

end MathieuQuadratic
