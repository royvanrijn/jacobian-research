import QuadraticNorms.PointwiseResidueEnergy
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

open scoped BigOperators

namespace QuadraticNorms

private theorem pointwise_phase_mul (θ : ℝ) (a b : ℤ) :
    pointwiseAdditivePhase θ a * pointwiseAdditivePhase θ b =
      pointwiseAdditivePhase θ (a + b) := by
  unfold pointwiseAdditivePhase
  rw [← Complex.exp_add, Int.cast_add]
  congr 1
  ring

private theorem pointwise_phase_neg_frequency (θ : ℝ) (a : ℤ) :
    pointwiseAdditivePhase (-θ) a = pointwiseAdditivePhase θ (-a) := by
  unfold pointwiseAdditivePhase
  rw [Complex.ofReal_neg, Int.cast_neg]
  congr 1
  ring

private theorem pointwise_phase_mul_frequency (θ : ℝ) (a b : ℤ) :
    pointwiseAdditivePhase ((a : ℝ) * θ) b = pointwiseAdditivePhase θ (a * b) := by
  unfold pointwiseAdditivePhase
  rw [Complex.ofReal_mul, Complex.ofReal_intCast, Int.cast_mul]
  congr 1
  ring

/-- Integrating e(θn) over one period detects the zero integer frequency. -/
theorem pointwise_integral_additive_phase (n : ℤ) :
    (∫ θ in (0 : ℝ)..1, pointwiseAdditivePhase θ n) = if n = 0 then 1 else 0 := by
  by_cases hn : n = 0
  · subst n
    simp [pointwiseAdditivePhase, intervalIntegral.integral_const]
  · let c : ℂ := (n : ℂ) * (2 * Real.pi * Complex.I)
    have hc : c ≠ 0 := by
      exact mul_ne_zero (by exact_mod_cast hn)
        (mul_ne_zero (mul_ne_zero (by norm_num) (by exact_mod_cast Real.pi_ne_zero)) Complex.I_ne_zero)
    have hphase (θ : ℝ) : pointwiseAdditivePhase θ n = Complex.exp (c * (θ : ℂ)) := by
      unfold pointwiseAdditivePhase
      congr 1
      dsimp [c]
      ring
    simp only [hphase, if_neg hn]
    rw [integral_exp_mul_complex hc]
    simp [c, Complex.exp_int_mul_two_pi_mul_I]

/-- The common Fourier reduction detects the signed shift n + kD = m.
The finite weighted family W represents the coefficients of the product
multiplier and can include repeated product labels. -/
theorem pointwise_fourier_comparison {ι κ ν : Type*}
    (s : Finset ι) (t : Finset κ) (u : Finset ν)
    (n : ι → ℤ) (m : κ → ℤ) (D : ν → ℤ)
    (F : ι → ℂ) (G : κ → ℂ) (W : ν → ℂ) (k : ℤ) :
    (∫ θ in (0 : ℝ)..1, pointwiseFourierSum s n F θ *
      pointwiseFourierSum t m G (-θ) * pointwiseFourierSum u D W ((k : ℝ) * θ)) =
      ∑ i ∈ s, ∑ j ∈ t, ∑ l ∈ u,
        if n i + k * D l = m j then F i * G j * W l else 0 := by
  let A (i : ι) (j : κ) (l : ν) (θ : ℝ) :=
    F i * G j * W l * pointwiseAdditivePhase θ (n i + k * D l - m j)
  have hcont (i : ι) (j : κ) (l : ν) : Continuous (A i j l) :=
    continuous_const.mul (pointwise_additive_phase_continuous _)
  have hproduct (θ : ℝ) : pointwiseFourierSum s n F θ *
      pointwiseFourierSum t m G (-θ) * pointwiseFourierSum u D W ((k : ℝ) * θ) =
        ∑ i ∈ s, ∑ j ∈ t, ∑ l ∈ u, A i j l θ := by
    unfold pointwiseFourierSum
    simp only [Finset.sum_mul]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Finset.mul_sum (s := t), Finset.sum_mul]
    simp only [Finset.mul_sum]
    refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun l _ => ?_
    rw [pointwise_phase_neg_frequency, pointwise_phase_mul_frequency]
    dsimp [A]
    have hs : (n i + -m j) + k * D l = n i + k * D l - m j := by ring
    calc
      _ = F i * G j * W l * ((pointwiseAdditivePhase θ (n i) *
          pointwiseAdditivePhase θ (-m j)) * pointwiseAdditivePhase θ (k * D l)) := by ring
      _ = _ := by rw [pointwise_phase_mul, pointwise_phase_mul, hs]
  simp only [hproduct]
  rw [intervalIntegral.integral_finsetSum fun i _ =>
    (continuous_finsetSum _ fun j _ =>
      continuous_finsetSum _ fun l _ => hcont i j l).intervalIntegrable _ _]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [intervalIntegral.integral_finsetSum fun j _ =>
    (continuous_finsetSum _ fun l _ => hcont i j l).intervalIntegrable _ _]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [intervalIntegral.integral_finsetSum fun l _ => (hcont i j l).intervalIntegrable _ _]
  refine Finset.sum_congr rfl fun l _ => ?_
  dsimp [A]
  rw [intervalIntegral.integral_const_mul, pointwise_integral_additive_phase]
  by_cases h : n i + k * D l = m j
  · simp [h]
  · simp [sub_eq_zero, h]

private theorem pointwise_fourier_sum_conj {ι : Type*} (s : Finset ι) (n : ι → ℤ)
    (F : ι → ℂ) (θ : ℝ) :
    pointwiseFourierSum s n (fun i => star (F i)) (-θ) = star (pointwiseFourierSum s n F θ) := by
  unfold pointwiseFourierSum
  simp only [← starRingEnd_apply, map_sum, map_mul]
  refine Finset.sum_congr rfl fun i _ => ?_
  unfold pointwiseAdditivePhase
  rw [← Complex.exp_conj]
  simp only [map_mul, map_ofNat, Complex.conj_ofReal, Complex.conj_I, map_intCast,
    Complex.ofReal_neg]
  apply congrArg (fun z => (starRingEnd ℂ) (F i) * Complex.exp z)
  ring

/-- Continuous Parseval on a full period for an endpoint with distinct integer labels. -/
theorem pointwise_fourier_period_parseval {ι : Type*} (s : Finset ι) (n : ι → ℤ)
    (F : ι → ℂ) (hn : Set.InjOn n (s : Set ι)) :
    (∫ θ in (0 : ℝ)..1, ‖pointwiseFourierSum s n F θ‖ ^ 2) = ∑ i ∈ s, ‖F i‖ ^ 2 := by
  classical
  have hm (θ : ℝ) :
      pointwiseFourierSum ({()} : Finset Unit) (fun _ => 0) (fun _ => 1) θ = 1 := by
    simp [pointwiseFourierSum, pointwiseAdditivePhase]
  have hcont : Continuous (fun θ => ‖pointwiseFourierSum s n F θ‖ ^ 2) :=
    (pointwise_fourier_sum_continuous s n F).norm.pow 2
  apply Complex.ofReal_injective
  change Complex.ofRealCLM (∫ θ in (0 : ℝ)..1, ‖pointwiseFourierSum s n F θ‖ ^ 2) = _
  rw [← Complex.ofRealCLM.intervalIntegral_comp_comm (hcont.intervalIntegrable _ _)]
  have hsq (z : ℂ) : Complex.ofRealCLM (‖z‖ ^ 2) = z * star z := by
    change ((‖z‖ ^ 2 : ℝ) : ℂ) = _
    rw [← Complex.normSq_eq_norm_sq, ← Complex.mul_conj, starRingEnd_apply]
  simp only [hsq]
  have h := pointwise_fourier_comparison s s ({()} : Finset Unit) n n (fun _ => 0) F
    (fun i => star (F i)) (fun _ => 1) 0
  simp only [pointwise_fourier_sum_conj, hm, mul_one, Finset.sum_singleton, zero_mul, add_zero] at h
  rw [h, Complex.ofReal_sum]
  refine Finset.sum_congr rfl fun i hi => ?_
  rw [Finset.sum_eq_single_of_mem i hi]
  · simp [← starRingEnd_apply, Complex.mul_conj, Complex.normSq_eq_norm_sq]
  · intro j hj hji
    rw [if_neg (fun he => hji (hn hi hj he).symm)]

end QuadraticNorms
