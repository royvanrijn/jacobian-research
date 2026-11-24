import QuadraticNorms.MinusElevenCuspBounds
import QuadraticNorms.MinusElevenDivisorShellBounds
import QuadraticNorms.MinusElevenRepresentationCounts
import QuadraticNorms.OrbitAveraging

open scoped BigOperators

namespace QuadraticNorms

/-- Finite operator data for the filtered shell at an auxiliary level
prime to the target. The return identity connects the operator to the
actual signed-coordinate bad count. -/
structure MinusElevenHeckeOrbit (Y : Type*) [Fintype Y] [DecidableEq Y]
    (n p : ℕ) (D : ℝ) where
  measure : Y → ℝ
  orbit : Finset Y
  base : Y
  character : Y → ℝ
  operator : (Y → ℝ) →ₗ[ℝ] (Y → ℝ)
  characterEigenvalue : ℝ
  measure_nonneg : ∀ y, 0 ≤ measure y
  probability : ∑ y, measure y = 1
  character_orthogonal : weightedPair measure 1 character = 0
  character_normalized : weightedPair measure character character = 1
  constant_eigen : operator 1 = minusElevenFilteredSigma n • 1
  character_eigen : operator character = characterEigenvalue • character
  self_adjoint : ∀ f g, weightedPair measure (operator f) g = weightedPair measure f (operator g)
  orbit_constant : ∀ y ∈ orbit, operator (orbitIndicator orbit) y = operator (orbitIndicator orbit) base
  complementary_bound :
    let f := orbitIndicator orbit
    let h := f - (∑ y ∈ orbit, measure y) • (1 : Y → ℝ) - weightedPair measure f character • character
    weightedPair measure (operator h) (operator h) ≤ D ^ 2 * weightedPair measure h h
  return_count : ((minusElevenBadShell n p).card : ℝ) = 4 * operator (orbitIndicator orbit) base

variable {Y : Type*} [Fintype Y] [DecidableEq Y]

noncomputable def MinusElevenHeckeOrbit.mass {n p : ℕ} {D : ℝ}
    (O : MinusElevenHeckeOrbit Y n p D) : ℝ :=
  ∑ y ∈ O.orbit, O.measure y

noncomputable def MinusElevenHeckeOrbit.characterMean {n p : ℕ} {D : ℝ}
    (O : MinusElevenHeckeOrbit Y n p D) : ℝ :=
  weightedPair O.measure (orbitIndicator O.orbit) O.character

private theorem MinusElevenHeckeOrbit.deviation_le {n p : ℕ} {D : ℝ}
    (O : MinusElevenHeckeOrbit Y n p D) (hD : 0 ≤ D) :
    O.mass * O.operator (orbitIndicator O.orbit) O.base -
      (minusElevenFilteredSigma n * O.mass ^ 2 + O.characterEigenvalue * O.characterMean ^ 2) ≤
        D * (O.mass - O.mass ^ 2 - O.characterMean ^ 2) := by
  exact (abs_le.mp (spectral_orbit_two_eigenfunctions O.measure_nonneg O.probability O.orbit O.base
    O.operator O.character hD O.constant_eigen O.character_eigen O.character_orthogonal
    O.character_normalized O.self_adjoint O.orbit_constant O.complementary_bound)).2

/-- At an inert unit level the character mean vanishes and the
principal-fiber mass gives the stated absolute bad-count bound. -/
theorem minusEleven_inert_bad_count_of_hecke_orbit {n p : ℕ} {D : ℝ}
    (hp : MinusElevenInertPrime p) (hD : 0 ≤ D) (O : MinusElevenHeckeOrbit Y n p D)
    (hMass : O.mass = (3 / 5) * minusElevenLocalMass p) (hMean : O.characterMean = 0) :
    ((minusElevenBadShell n p).card : ℝ) ≤
      (12 / 5) * minusElevenFilteredSigma n * minusElevenLocalMass p + 4 * D := by
  have hpR : 2 < (p : ℝ) := by exact_mod_cast hp.2.1
  have hp₀ : 0 < (p : ℝ) := by linarith
  have hp₁ : 0 < (p : ℝ) - 1 := by linarith
  have hμ : 0 < O.mass := by
    rw [hMass, minusElevenLocalMass]
    positivity
  have h := O.deviation_le hD
  rw [hMean] at h
  have hpoint : O.operator (orbitIndicator O.orbit) O.base ≤
      minusElevenFilteredSigma n * O.mass + D * (1 - O.mass) := by
    apply (mul_le_mul_iff_of_pos_left hμ).mp
    nlinarith only [h]
  rw [O.return_count]
  calc
    _ ≤ 4 * (minusElevenFilteredSigma n * O.mass + D * (1 - O.mass)) :=
      mul_le_mul_of_nonneg_left hpoint (by norm_num)
    _ = (12 / 5) * minusElevenFilteredSigma n * minusElevenLocalMass p + 4 * D - 4 * D * O.mass := by
      conv_lhs => rw [hMass]
      conv_rhs => rw [hMass]
      ring
    _ ≤ _ := sub_le_self _ (by positivity)

/-- At eleven the base-point indicator has mass and character mean
both equal to one tenth. -/
theorem minusEleven_ramified_bad_count_of_hecke_orbit {n : ℕ} {D : ℝ}
    (hD : 0 ≤ D) (O : MinusElevenHeckeOrbit Y n 11 D)
    (hMass : O.mass = 1 / 10) (hMean : O.characterMean = 1 / 10)
    (hEigen : O.characterEigenvalue ≤ minusElevenFilteredSigma n) :
    ((minusElevenBadShell n 11).card : ℝ) ≤
      (4 / 5) * minusElevenFilteredSigma n + (16 / 5) * D := by
  have h := O.deviation_le hD
  rw [hMass, hMean] at h
  rw [O.return_count]
  linarith only [h, hEigen]

private theorem minusEleven_relative_bad_bound {t m p : ℕ} {c : ℝ}
    (ht : 0 < t) (hc : 0 ≤ c)
    (hShell : (12 / 5 : ℝ) * minusElevenFilteredSigma (2 ^ t * m) * (1 - minusElevenFilteredError t m) ≤
      (minusElevenOddShell (2 ^ t * m)).card)
    (hBad : ((minusElevenBadShell (2 ^ t * m) p).card : ℝ) ≤
      (12 / 5) * minusElevenFilteredSigma (2 ^ t * m) * c) :
    (1 - minusElevenFilteredError t m) * (minusElevenBadShell (2 ^ t * m) p).card ≤
      c * (minusElevenOddShell (2 ^ t * m)).card := by
  have hscale : 0 ≤ 1 - minusElevenFilteredError t m :=
    sub_nonneg.mpr (minusEleven_filtered_error_lt_one ht).le
  calc
    _ ≤ _ := mul_le_mul_of_nonneg_left hBad hscale
    _ = c * ((12 / 5) * minusElevenFilteredSigma (2 ^ t * m) * (1 - minusElevenFilteredError t m)) := by
      ring
    _ ≤ _ := mul_le_mul_of_nonneg_left hShell hc

/-- The unrestricted −11 row follows from the theta identity,
unramified Hecke coefficient data, and finite operator data at eleven
and the inert levels prime to the reduced target. -/
theorem P11_sum_of_theta_and_hecke_orbits
    (a : ℕ → ℝ) (hTheta : MinusElevenThetaFormula a) (ha : MinusElevenHeckeCoefficients a)
    (spaces : ℕ → Type*) [∀ p, Fintype (spaces p)] [∀ p, DecidableEq (spaces p)]
    (hRamified : ∀ t m : ℕ, 0 < t → Odd m → 2 ^ 25 < 2 ^ t * m → ¬ 11 ∣ 2 ^ t * m →
      IsReduced C11 (2 ^ t * m) →
      ∃ O : MinusElevenHeckeOrbit (spaces 11) (2 ^ t * m) 11 (minusElevenCuspErrorBound t m),
        O.mass = 1 / 10 ∧ O.characterMean = 1 / 10 ∧
          O.characterEigenvalue ≤ minusElevenFilteredSigma (2 ^ t * m))
    (hInert : ∀ t m p : ℕ, 0 < t → Odd m → 2 ^ 25 < 2 ^ t * m → ¬ 11 ∣ 2 ^ t * m →
      IsReduced C11 (2 ^ t * m) → MinusElevenInertPrime p → ¬ p ∣ 2 ^ t * m →
      ∃ O : MinusElevenHeckeOrbit (spaces p) (2 ^ t * m) p (minusElevenCuspErrorBound t m),
        O.mass = (3 / 5) * minusElevenLocalMass p ∧ O.characterMean = 0)
    (n : ℕ) :
    InSumset P11 P11 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 121 ∣ n := by
  apply P11_sum_of_spectral_shell_bounds ?_ n
  intro t m ht hm hcut h₁₁ hred
  let v := 2 ^ t * m
  let E := minusElevenFilteredError t m
  let S := minusElevenFilteredSigma v
  let D := minusElevenCuspErrorBound t m
  have hE : 0 ≤ E := minusEleven_filtered_error_nonneg t m
  have hD : 0 ≤ D := by unfold D minusElevenCuspErrorBound; positivity
  have hSE : (3 * S) * E = 2 * D := minusEleven_filtered_error_mul ht hm
  have hCusp := minusEleven_filtered_cusp_bound ha ht hm h₁₁
  have hShell := (minusEleven_odd_shell_bounds_of_theta a hTheta ht hm h₁₁ hCusp).1
  refine ⟨minusEleven_shell_lower_of_theta a hTheta ht hm h₁₁ hCusp
    (minusEleven_filtered_error_lt_one ht).le, ?_, fun p hp => ?_⟩
  · obtain ⟨O, hM, hB, hEig⟩ := hRamified t m ht hm hcut h₁₁ hred
    apply minusEleven_relative_bad_bound ht (by positivity : 0 ≤ 1 / 3 + 2 * E) hShell
    have hb := minusEleven_ramified_bad_count_of_hecke_orbit hD O hM hB hEig
    change ((minusElevenBadShell v 11).card : ℝ) ≤ 4 / 5 * S + 16 / 5 * D at hb
    change ((minusElevenBadShell v 11).card : ℝ) ≤ 12 / 5 * S * (1 / 3 + 2 * E)
    nlinarith only [hb, hSE]
  · by_cases hpv : p ∣ v
    · exact minusEleven_inert_divisor_small_bound a hTheta
        (fun _ _ ht' hm' h₁₁' => minusEleven_filtered_cusp_bound ha ht' hm' h₁₁') ht hm h₁₁ hp hpv
    · obtain ⟨O, hM, hB⟩ := hInert t m p ht hm hcut h₁₁ hred hp hpv
      have hpR : 2 < (p : ℝ) := by exact_mod_cast hp.2.1
      have hm₀ : 0 ≤ minusElevenLocalMass p := by
        unfold minusElevenLocalMass
        have h₀ : 0 < (p : ℝ) := by linarith
        have h₁ : 0 < (p : ℝ) - 1 := by linarith
        positivity
      apply minusEleven_relative_bad_bound ht (by positivity : 0 ≤ minusElevenLocalMass p + (5 / 2) * E) hShell
      have hb := minusEleven_inert_bad_count_of_hecke_orbit hp hD O hM hB
      change ((minusElevenBadShell v p).card : ℝ) ≤ 12 / 5 * S * minusElevenLocalMass p + 4 * D at hb
      change ((minusElevenBadShell v p).card : ℝ) ≤ 12 / 5 * S * (minusElevenLocalMass p + (5 / 2) * E)
      nlinarith only [hb, hSE]

end QuadraticNorms
