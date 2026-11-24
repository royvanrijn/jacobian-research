import QuadraticNorms.HeckeCoefficientBounds
import QuadraticNorms.MinusEightCompletion
import QuadraticNorms.OrbitAveraging

open scoped BigOperators ArithmeticFunction.sigma

namespace QuadraticNorms

/-- The odd level-32 coefficient bound follows from the unramified
prime bounds and Hecke recurrences. -/
theorem minusEight_cusp_coefficient_bound {a : ℕ → ℝ}
    (ha : WeightTwoHeckeCoefficients a (· ≠ 2)) {m : ℕ} (hm : Odd m) :
    |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m := by
  apply weight_two_hecke_coefficient_bound ha (by obtain ⟨j, hj⟩ := hm; omega)
  intro p _ hpd hp₂
  exact hm.not_two_dvd_nat (hp₂ ▸ hpd)

/-- The finite data at an auxiliary inert unit level. The odd operator
has the unramified bound off the two global norm characters. The
dyadic operator is a nonnegative reversible kernel of degree c_e. -/
structure MinusEightHeckeOrbit (Y : Type*) [Fintype Y] [DecidableEq Y]
    (e m p : ℕ) where
  measure : Y → ℝ
  orbit : Finset Y
  base : Y
  character : Y → ℝ
  oddOperator : (Y → ℝ) →ₗ[ℝ] (Y → ℝ)
  dyadicOperator : (Y → ℝ) →ₗ[ℝ] (Y → ℝ)
  dyadicKernel : Y → Y → ℝ
  measure_nonneg : ∀ x, 0 ≤ measure x
  probability : ∑ x, measure x = 1
  orbit_mass : ∑ x ∈ orbit, measure x = minusEightLocalMass p
  character_orthogonal : weightedPair measure 1 character = 0
  character_mean_zero : weightedPair measure (orbitIndicator orbit) character = 0
  odd_degree : oddOperator 1 = (σ 1 m : ℝ) • 1
  odd_self_adjoint : ∀ f g, weightedPair measure (oddOperator f) g = weightedPair measure f (oddOperator g)
  odd_complementary_bound : ∀ f, weightedPair measure f 1 = 0 →
    weightedPair measure f character = 0 →
    weightedPair measure (oddOperator f) (oddOperator f) ≤
      ((σ 0 m : ℝ) * Real.sqrt m) ^ 2 * weightedPair measure f f
  dyadic_apply : ∀ f x, dyadicOperator f x = ∑ y, dyadicKernel x y * f y
  dyadic_nonneg : ∀ x y, 0 ≤ dyadicKernel x y
  dyadic_balance : ∀ x y, measure x * dyadicKernel x y = measure y * dyadicKernel y x
  dyadic_row : ∀ x, ∑ y, dyadicKernel x y = minusEightShellFactor e
  commute : ∀ f, oddOperator (dyadicOperator f) = dyadicOperator (oddOperator f)
  orbit_constant : ∀ x ∈ orbit,
    dyadicOperator (oddOperator (orbitIndicator orbit)) x =
      dyadicOperator (oddOperator (orbitIndicator orbit)) base
  return_count : ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) =
    4 * dyadicOperator (oddOperator (orbitIndicator orbit)) base

/-- The dyadic kernel and the odd-index spectral estimate give the
sharp local bad-count bound at auxiliary primes coprime to the target. -/
theorem minusEight_unit_level_bad_count_of_hecke_orbit {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {Y : Type*} [Fintype Y] [DecidableEq Y] {e m p : ℕ}
    (he : e ≤ 3) (hm : Odd m) (hp : MinusEightInertPrime p)
    (O : MinusEightHeckeOrbit Y e m p) :
    ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
      (minusEightLocalMass p + (1 - minusEightLocalMass p) * minusEightShellError m) *
        (minusEightFullShell (2 ^ e * m)).card := by
  let c : ℝ := minusEightShellFactor e
  let σm : ℝ := σ 1 m
  let D := (σ 0 m : ℝ) * Real.sqrt m
  let μ := ∑ x ∈ O.orbit, O.measure x
  let ϕ := orbitIndicator O.orbit
  let U := O.dyadicOperator.comp O.oddOperator
  have hm₀ : 0 < m := by obtain ⟨k, hk⟩ := hm; omega
  have hσ : 0 < σm := Nat.cast_pos.mpr (ArithmeticFunction.sigma_pos _ _ hm₀.ne')
  have hcN : 0 < minusEightShellFactor e := by unfold minusEightShellFactor; split_ifs <;> decide
  have hc : 0 < c := Nat.cast_pos.mpr hcN
  have hμ : 0 < μ := by
    have hp₂ : 1 < p := hp.1.one_lt
    have hpR : (1 : ℝ) < p := by exact_mod_cast hp₂
    dsimp [μ]
    rw [O.orbit_mass]
    unfold minusEightLocalMass
    exact div_pos (by norm_num) (mul_pos (by linarith) (sub_pos.mpr hpR))
  have hA : O.dyadicOperator 1 = c • 1 := by
    ext x
    rw [O.dyadic_apply]
    simp only [Pi.one_apply, mul_one, Pi.smul_apply, smul_eq_mul]
    exact O.dyadic_row x
  have hU : U 1 = (c * σm) • 1 := by
    dsimp [U]
    rw [O.odd_degree, map_smul, hA, smul_smul, mul_comm]
  have hselfA (f g : Y → ℝ) :
      weightedPair O.measure (O.dyadicOperator f) g = weightedPair O.measure f (O.dyadicOperator g) := by
    simpa only [← O.dyadic_apply] using
      weightedPair_kernel_self_adjoint O.measure f g O.dyadicKernel O.dyadic_balance
  have hselfU (f g : Y → ℝ) :
      weightedPair O.measure (U f) g = weightedPair O.measure f (U g) := by
    dsimp [U]
    rw [hselfA, O.odd_self_adjoint, O.commute]
  have hspec :
      let h := ϕ - μ • (1 : Y → ℝ)
      weightedPair O.measure (U h) (U h) ≤ (c * D) ^ 2 * weightedPair O.measure h h := by
    intro h
    have h₁ : weightedPair O.measure h 1 = 0 := by
      have hone : weightedPair O.measure (1 : Y → ℝ) 1 = 1 := by simpa [weightedPair] using O.probability
      change weightedPair O.measure (ϕ - μ • (1 : Y → ℝ)) 1 = 0
      rw [weightedPair_sub_left, weightedPair_smul_left, hone, weightedPair_indicator_left]
      simp [μ]
    have hχ : weightedPair O.measure h O.character = 0 := by
      change weightedPair O.measure (ϕ - μ • (1 : Y → ℝ)) O.character = 0
      rw [weightedPair_sub_left, weightedPair_smul_left, O.character_orthogonal, O.character_mean_zero]
      ring
    have hT := O.odd_complementary_bound h h₁ hχ
    have hNorm := weightedPair_nonnegative_kernel_le O.measure (O.oddOperator h) O.dyadicKernel
      O.measure_nonneg O.dyadic_nonneg O.dyadic_balance O.dyadic_row
    have hA' : weightedPair O.measure (U h) (U h) ≤
        c ^ 2 * weightedPair O.measure (O.oddOperator h) (O.oddOperator h) := by
      simpa only [U, LinearMap.comp_apply, ← O.dyadic_apply] using hNorm
    exact hA'.trans ((mul_le_mul_of_nonneg_left hT (sq_nonneg c)).trans_eq (by dsimp [D]; ring))
  have hdev := spectral_orbit_deviation O.measure_nonneg O.probability O.orbit O.base U
    (mul_pos hc hσ) (by dsimp [D]; positivity) hμ hU hselfU O.orbit_constant hspec
  have hquot : c * D / (c * σm) = minusEightShellError m := by
    rw [minusEightShellError]
    change c * D / (c * σm) = D / σm
    field_simp [hc.ne']
  have hrel := (abs_le.mp hdev).2
  change U ϕ O.base / (c * σm) - μ ≤ c * D / (c * σm) * (1 - μ) at hrel
  rw [hquot] at hrel
  have hrel' : U ϕ O.base / (c * σm) ≤ μ + (1 - μ) * minusEightShellError m := by
    rw [mul_comm (1 - μ)]
    linarith only [hrel]
  have hb := (div_le_iff₀ (mul_pos hc hσ)).mp hrel'
  have hr : ((minusEightFullShell (2 ^ e * m)).card : ℝ) = 4 * (c * σm) := by
    rw [minusEight_full_shell_dyadic hθ he hm]
    dsimp [c, σm]
    ring
  rw [hr, O.return_count]
  have hmass : μ = minusEightLocalMass p := O.orbit_mass
  rw [← hmass]
  calc
    _ ≤ 4 * ((μ + (1 - μ) * minusEightShellError m) * (c * σm)) :=
      mul_le_mul_of_nonneg_left hb (by norm_num)
    _ = _ := by ring

/-- The unit-level sieve estimate, in the form used for classification. -/
theorem minusEight_unit_level_bad_bound {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    {Y : Type*} [Fintype Y] [DecidableEq Y] {e m p : ℕ} (he : e ≤ 3) (hm : Odd m)
    (hp : MinusEightInertPrime p) (O : MinusEightHeckeOrbit Y e m p) :
    ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
      (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card := by
  refine (minusEight_unit_level_bad_count_of_hecke_orbit hθ he hm hp O).trans ?_
  have h := mul_nonneg (minusEight_local_mass_nonneg p) (minusEight_shell_error_nonneg m)
  apply mul_le_mul_of_nonneg_right _ (by positivity)
  nlinarith only [h]

/-- The exact −8 sumsets and quantitative count follow from the
theta identities, weight-two coefficient data, and finite unit-level
Hecke orbit data. The sieve and its full interval certificate are checked. -/
theorem minusEight_sumset_theorem_of_hecke_data {a : ℕ → ℝ} (hθ : MinusEightThetaFormula a)
    (ha : WeightTwoHeckeCoefficients a (· ≠ 2))
    (spaces : ℕ → Type*) [∀ p, Fintype (spaces p)] [∀ p, DecidableEq (spaces p)]
    (hOrbits : ∀ e m p : ℕ, e ≤ 3 → Odd m → p ≤ 2 ^ e * m → MinusEightInertPrime p →
      ¬ p ∣ 2 ^ e * m → Nonempty (MinusEightHeckeOrbit (spaces p) e m p)) :
    (∀ n, (InSumset C8 C8 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 8 ∣ n) ∧
      (InSumset C8 (Dilate 2 C8) n ↔ 3 ≤ n ∧ n % 2 = 1) ∧
      (InSumset P8 P8 n ↔ 2 ≤ n ∧ ¬ 16 ∣ n)) ∧
    (∀ n : ℕ, 2 ^ 78 ≤ n → ¬ 16 ∣ n →
      (n : ℝ) ^ (2 / 3 : ℝ) / 4800 ≤ (minusEightOrderedPairCount n : ℝ)) := by
  have hsmall (e m p : ℕ) (he : e ≤ 3) (hm : Odd m) (hpn : p ≤ 2 ^ e * m)
      (hp : MinusEightInertPrime p) (hunit : ¬ p ∣ 2 ^ e * m) :
      ((minusEightEitherDivisible (2 ^ e * m) p).card : ℝ) ≤
        (minusEightLocalMass p + minusEightShellError m) * (minusEightFullShell (2 ^ e * m)).card := by
    obtain ⟨O⟩ := hOrbits e m p he hm hpn hp hunit
    exact minusEight_unit_level_bad_bound hθ he hm hp O
  exact ⟨minusEight_exact_sumsets_of_theta_and_unit hθ (fun _ => minusEight_cusp_coefficient_bound ha) hsmall,
    fun _ => minusEight_large_ordered_pair_count hθ (fun _ => minusEight_cusp_coefficient_bound ha) hsmall⟩

end QuadraticNorms
