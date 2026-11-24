import MathieuQuadratic.Endpoint
import MathieuQuadratic.Witness

/-!
# From the weight and height identities to the mixed moments

The trace is on a commutative K-algebra with an invertible norm U and a
height S. The hypotheses here are the scalar weight calculation and the
rank-two Lie recurrence in the uniform quadratic-complement proof.
-/

open Polynomial
open scoped BigOperators Interval

namespace MathieuQuadratic

variable {K A : Type*} [Field K] [CharZero K] [CommRing A] [Algebra K A]

/-- The circuit (11), with division only by the norm unit. -/
def normalizerPhase (U : Aˣ) (S : A) : A :=
  (↑(U⁻¹) : A) * (1 + (U : A)) * (1 - (1 + (U : A)) ^ 2 * S) ^ 2

theorem normalizerPhase_eq_regular (U : Aˣ) (B : A) :
    normalizerPhase U (B * (↑(U⁻¹) : A)) =
      ((U : A) + 1) * ((U : A) - ((U : A) + 1) ^ 2 * B) ^ 2 * (↑(U⁻¹) : A) ^ 3 := by
  have hunit : (U : A) * (↑(U⁻¹) : A) = 1 := Units.mul_inv _
  have hs : 1 - (1 + (U : A)) ^ 2 * (B * (↑(U⁻¹) : A)) =
      ((U : A) - ((U : A) + 1) ^ 2 * B) * (↑(U⁻¹) : A) := by
    rw [sub_mul, hunit]
    ring
  rw [normalizerPhase, hs]
  ring

/-- The solution of the height recurrence (9). -/
theorem heightMoment_of_recurrence (trace : A →ₗ[K] K) (S : A) (M : ℕ) (hM : 0 < M)
    (norm : trace 1 = 1)
    (recurrence : ∀ j : ℕ, (↑(M + 2 * (j + 1)) : K) * trace (S ^ (j + 1)) =
      (↑(M + 2 * j) : K) * trace (S ^ j)) :
    ∀ j : ℕ, trace (S ^ j) = (M : K) / (↑(M + 2 * j) : K) := by
  have hn (j : ℕ) : (↑(M + 2 * j) : K) ≠ 0 := by
    exact_mod_cast (show M + 2 * j ≠ 0 by omega)
  intro j
  induction j with
  | zero => simp [norm, div_self (by exact_mod_cast (Nat.ne_of_gt hM) : (M : K) ≠ 0)]
  | succ j ih =>
      apply (eq_div_iff (hn (j + 1))).2
      rw [mul_comm, recurrence, ih, mul_div_cancel₀ _ (hn j)]

private def evalAtUnit (K : Type*) [Field K] [CharZero K] [Algebra K A] (U : Aˣ) :
    ℚ[X] →+* A :=
  eval₂RingHom ((algebraMap K A).comp (algebraMap ℚ K)) (U : A)

private theorem trace_scalar (trace : A →ₗ[K] K) (r : ℚ) (y : A) :
    trace (algebraMap K A (algebraMap ℚ K r) * y) = algebraMap ℚ K r * trace y := by
  rw [← Algebra.smul_def, map_smul, smul_eq_mul]

private theorem trace_evalAtUnit (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
    (weight : ∀ (a : ℤ) (j : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ j) = 0)
    (d : ℤ) (j : ℕ) (p : ℚ[X]) :
    trace ((↑(U ^ (-d)) : A) * evalAtUnit K U p * S ^ j) =
      if 0 ≤ d then algebraMap ℚ K (p.coeff d.toNat) * trace (S ^ j) else 0 := by
  induction p using Polynomial.induction_on' with
  | add p q hp hq =>
      rw [map_add, mul_add, add_mul, map_add, hp, hq, coeff_add, map_add]
      split_ifs <;> ring
  | monomial n r =>
      have he : evalAtUnit K U (monomial n r) =
          algebraMap K A (algebraMap ℚ K r) * (U : A) ^ n :=
        eval₂_monomial _ _
      rw [he, mul_left_comm (↑(U ^ (-d)) : A), mul_assoc, trace_scalar,
        ← Units.val_pow_eq_pow_val, ← zpow_natCast, ← Units.val_mul, ← zpow_add]
      by_cases hn : (n : ℤ) = d
      · subst d
        simp
      · rw [weight _ _ (by omega), mul_zero, coeff_monomial]
        split_ifs with hd hnd
        · exact (hn (by omega)).elim
        · simp
        · rfl

/-- The product phase on the norm complement. -/
def productPhase (U : Aˣ) (S : A) (α β γ : ℕ) : A :=
  (↑(U⁻¹) : A) ^ α * (1 + (U : A)) ^ β * (1 - (1 + (U : A)) ^ 2 * S) ^ γ

def perturbedProductPhase (U : Aˣ) (S : A) (M : ℕ) : A :=
  productPhase U S M M M * (1 + ((U : A) * S) * ((M - 1 : ℕ) + (M : A) * (U : A)))

private noncomputable def productTerm (v k i : ℕ) : ℚ[X] :=
  C ((k.choose i : ℚ) * (-1) ^ i) * (X + 1) ^ (v + 2 * i)

private theorem productPhase_pow (U : Aˣ) (S : A) (α β γ m : ℕ) :
    productPhase U S α β γ ^ m =
      (↑(U ^ (-(α * m : ℕ) : ℤ)) : A) *
        ∑ i ∈ Finset.range (γ * m + 1), evalAtUnit K U (productTerm (β * m) (γ * m) i) * S ^ i := by
  have expand : (1 - (1 + (U : A)) ^ 2 * S) ^ (γ * m) =
      ∑ i ∈ Finset.range (γ * m + 1),
        (-((1 + (U : A)) ^ 2 * S)) ^ i * ((γ * m).choose i : A) := by
    rw [sub_eq_add_neg, add_comm 1, add_pow]
    simp only [one_pow, mul_one]
  rw [productPhase, mul_pow, mul_pow, ← pow_mul, ← pow_mul, ← pow_mul, zpow_neg,
    zpow_natCast, ← inv_pow, Units.val_pow_eq_pow_val, mul_assoc, expand, Finset.mul_sum]
  congr 1
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  change _ = (productTerm (β * m) (γ * m) i).eval₂
    ((algebraMap K A).comp (algebraMap ℚ K)) (U : A) * S ^ i
  rw [productTerm, eval₂_mul, eval₂_C, eval₂_pow, eval₂_add, eval₂_X, eval₂_one]
  simp only [map_mul, map_natCast, map_pow, map_neg, map_one]
  rw [neg_pow, mul_pow ((1 + (U : A)) ^ 2) S, ← pow_mul, pow_add]
  ring

private theorem powerEndpointPolynomial_eq_sum (M v k : ℕ) :
    powerEndpointPolynomial M v k =
      ∑ i ∈ Finset.range (k + 1),
        C ((M : ℚ) / (↑(M + 2 * i) : ℚ)) * productTerm v k i := by
  rw [powerEndpointPolynomial, powerBetaQuotient, sum_comp, Finset.mul_sum]
  refine Finset.sum_congr rfl fun i _ ↦ ?_
  simp only [mul_comp, C_comp, pow_comp, X_comp,
    productTerm, powerBetaCoefficient, map_mul, pow_add]
  ring

/-- The rational beta coefficient of the height moment. -/
noncomputable def mixedBetaMoment (M k j : ℕ) : ℚ :=
  ((M : ℚ) / (M + 2 * j : ℕ)) * powerBetaMoment (M + 2 * j) k

/-- Laurent coefficient extraction for all indices, before Taylor truncation. -/
theorem productMoment_of_weights (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
    (M α β γ m j : ℕ) (hM : 0 < M) (ℓ : ℤ)
    (weight : ∀ (a : ℤ) (i : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ i) = 0)
    (height : ∀ i : ℕ, trace (S ^ i) = (M : K) / (↑(M + 2 * i) : K)) :
    trace ((↑(U ^ ℓ) : A) * S ^ j * productPhase U S α β γ ^ m) =
      if ℓ ≤ (α * m : ℕ) then
        algebraMap ℚ K (((M : ℚ) / (M + 2 * j : ℕ)) *
          (powerEndpointPolynomial (M + 2 * j) (β * m) (γ * m)).coeff
            ((α * m : ℕ) - ℓ : ℤ).toNat)
      else 0 := by
  rw [productPhase_pow (K := K), ← mul_assoc, mul_right_comm _ (S ^ j),
    ← Units.val_mul, ← zpow_add,
    show ℓ + -(α * m : ℕ) = -((α * m : ℕ) - ℓ : ℤ) by omega, mul_assoc]
  simp_rw [Finset.mul_sum, map_sum,
    show ∀ p : A, ∀ i : ℕ, (↑(U ^ (-((α * m : ℕ) - ℓ : ℤ))) : A) *
        (S ^ j * (p * S ^ i)) =
        (↑(U ^ (-((α * m : ℕ) - ℓ : ℤ))) : A) * p * S ^ (j + i) by
      intros; rw [pow_add]; ring,
    trace_evalAtUnit trace U S weight, height, sub_nonneg]
  by_cases hℓ : ℓ ≤ (α * m : ℕ)
  · simp only [hℓ, if_true, powerEndpointPolynomial_eq_sum, finsetSum_coeff, map_sum,
      Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ ↦ ?_
    rw [coeff_C_mul, show M + 2 * (j + i) = M + 2 * j + 2 * i by omega,
      map_mul, map_mul, map_div₀, map_div₀]
    simp only [map_natCast]
    have hne : (↑(M + 2 * j) : K) ≠ 0 := by exact_mod_cast (show M + 2 * j ≠ 0 by omega)
    field_simp
  · simp only [hℓ, if_false, Finset.sum_const_zero]

theorem mixedBetaMoment_pos (M k j : ℕ) (hM : 0 < M) :
    0 < mixedBetaMoment M k j :=
  mul_pos (by positivity) (powerBetaMoment_pos _ _ (by omega))

theorem mixedBetaMoment_cast_eq_integral (M k j : ℕ) (hM : 0 < M) :
    (mixedBetaMoment M k j : ℝ) =
      (M : ℝ) * ∫ w in (0 : ℝ)..1, w ^ (M + 2 * j - 1) * (1 - w ^ 2) ^ k := by
  unfold mixedBetaMoment
  push_cast
  rw [powerBetaMoment_cast_eq_integral _ _ (by omega)]
  simp_rw [mul_assoc]
  rw [intervalIntegral.integral_const_mul]
  push_cast
  have hne : (M : ℝ) + 2 * j ≠ 0 := by positivity
  rw [← mul_assoc, div_mul_cancel₀ _ hne]

/-- The beta recurrence in the height index used in the tail test. -/
theorem mixedBetaMoment_succ_height (M k j : ℕ) (hM : 0 < M) :
    (↑(M + 2 * j) : ℚ) * mixedBetaMoment M k j =
      (↑(M + 2 * j + 2 * k + 2) : ℚ) * mixedBetaMoment M k (j + 1) := by
  let r := M + 2 * j
  have hr : 0 < r := by omega
  rw [mixedBetaMoment, mixedBetaMoment, show M + 2 * (j + 1) = r + 2 by omega]
  change (r : ℚ) * ((M : ℚ) / r * powerBetaMoment r k) =
    (↑(r + 2 * k + 2) : ℚ) * ((M : ℚ) / (↑(r + 2) : ℚ) * powerBetaMoment (r + 2) k)
  rw [powerBetaMoment_closed_form _ _ hr, powerBetaMoment_closed_form _ _ (by omega)]
  let x : ℚ := 1 + (r : ℚ) / 2
  have hx : 0 < x := by dsimp [x]; positivity
  have he : (ascPochhammer ℚ k).eval x * (x + k) =
      x * (ascPochhammer ℚ k).eval (x + 1) := by
    rw [← ascPochhammer_succ_eval, ascPochhammer_succ_left, eval_mul, eval_X, eval_comp,
      eval_add, eval_X, eval_one]
  have hx₀ : 0 < (ascPochhammer ℚ k).eval x := ascPochhammer_pos _ _ hx
  have hx₁ : 0 < (ascPochhammer ℚ k).eval (x + 1) := ascPochhammer_pos _ _ (by positivity)
  rw [show (1 : ℚ) + (r : ℚ) / 2 = x by rfl,
    show (1 : ℚ) + (↑(r + 2) : ℚ) / 2 = x + 1 by push_cast; dsimp [x]; ring]
  have hr₀ : (r : ℚ) ≠ 0 := by positivity
  have hr₂ : (↑(r + 2) : ℚ) ≠ 0 := by positivity
  have hmul : (↑(r + 2) : ℚ) * (ascPochhammer ℚ k).eval (x + 1) =
      (↑(r + 2 * k + 2) : ℚ) * (ascPochhammer ℚ k).eval x := by
    rw [show (↑(r + 2) : ℚ) = 2 * x by push_cast; dsimp [x]; ring,
      show (↑(r + 2 * k + 2) : ℚ) = 2 * (x + k) by push_cast; dsimp [x]; ring,
      mul_assoc, ← he]
    ring
  field_simp [hr₀, hr₂, hx₀.ne', hx₁.ne']
  linear_combination hmul

/-- The height-dependent factor in the rational multiplier test. -/
theorem mixedBetaMoment_height_factor (M k j : ℕ) (hM : 0 < M) :
    mixedBetaMoment M k j = mixedBetaMoment M k 0 *
      ((ascPochhammer ℚ j).eval ((M : ℚ) / 2) /
        (ascPochhammer ℚ j).eval ((M : ℚ) / 2 + k + 1)) := by
  induction j with
  | zero => simp
  | succ j ih =>
    have hc := mixedBetaMoment_succ_height M k j hM
    have hv : 0 < (ascPochhammer ℚ j).eval ((M : ℚ) / 2 + k + 1) :=
      ascPochhammer_pos _ _ (by positivity)
    have hW : (↑(M + 2 * j + 2 * k + 2) : ℚ) ≠ 0 := by positivity
    apply (mul_right_inj' hW).mp
    rw [← hc, ih, ascPochhammer_succ_eval, ascPochhammer_succ_eval]
    push_cast
    have hsum : (M : ℚ) / 2 + k + 1 + j ≠ 0 := by positivity
    field_simp [hv.ne', hsum]
    ring

/-- The mixed-moment formula of the draft, including the negative-index vanishing. -/
theorem product_mixedMoment (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
    (M α β γ m j : ℕ) (hM : 0 < M) (ℓ : ℤ) (hℓ : (α : ℤ) * m - (γ : ℤ) * m ≤ ℓ)
    (weight : ∀ (a : ℤ) (i : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ i) = 0)
    (height : ∀ i : ℕ, trace (S ^ i) = (M : K) / (↑(M + 2 * i) : K)) :
    trace ((↑(U ^ ℓ) : A) * S ^ j * productPhase U S α β γ ^ m) =
      if ℓ ≤ (α * m : ℕ) then algebraMap ℚ K
        (mixedBetaMoment M (γ * m) j *
          ((Ring.choose ((β * m : ℕ) - (M + 2 * j : ℕ) : ℤ)
            ((α * m : ℕ) - ℓ : ℤ).toNat : ℤ) : ℚ))
      else 0 := by
  rw [productMoment_of_weights trace U S M α β γ m j hM ℓ weight height]
  split_ifs with h
  · rw [powerEndpointPolynomial_coeff _ _ _ _ (by omega) (by omega), mixedBetaMoment, mul_assoc]
  · rfl

/-- The exact discrepancy at one step below the threshold in the draft. -/
theorem product_mixedMoment_first_error (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
    (M α β γ m j : ℕ) (hM : 0 < M)
    (weight : ∀ (a : ℤ) (i : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ i) = 0)
    (height : ∀ i : ℕ, trace (S ^ i) = (M : K) / (↑(M + 2 * i) : K)) :
    trace ((↑(U ^ ((α : ℤ) * m - (γ : ℤ) * m - 1)) : A) * S ^ j *
        productPhase U S α β γ ^ m) =
      algebraMap ℚ K (mixedBetaMoment M (γ * m) j *
        ((Ring.choose ((β * m : ℕ) - (M + 2 * j : ℕ) : ℤ) (γ * m + 1) : ℤ) : ℚ) +
          (M : ℚ) * (-2) ^ (γ * m) / (γ * m + 1 : ℕ)) := by
  rw [productMoment_of_weights trace U S M α β γ m j hM _ weight height,
    if_pos (by have := mul_nonneg (Int.natCast_nonneg γ) (Int.natCast_nonneg m); omega : (α : ℤ) * m - (γ : ℤ) * m - 1 ≤ (α * m : ℕ)),
    show ((α * m : ℕ) - ((α : ℤ) * m - (γ : ℤ) * m - 1) : ℤ).toNat = γ * m + 1 by omega,
    powerEndpointPolynomial_coeff_succ _ _ _ (by omega), mul_add, mixedBetaMoment, mul_assoc]
  congr 2
  have hn : (↑(M + 2 * j) : ℚ) ≠ 0 := by exact_mod_cast (show M + 2 * j ≠ 0 by omega)
  field_simp

/-- All positive pure moments vanish, and the norm-power multiplier
has exactly the constant beta moment. This includes the ternary case. -/
theorem product_uniform_moments (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
    (M γ m : ℕ) (hM : 0 < M) (hγ : M ≤ γ) (hm : 0 < m)
    (weight : ∀ (a : ℤ) (i : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ i) = 0)
    (height : ∀ i : ℕ, trace (S ^ i) = (M : K) / (↑(M + 2 * i) : K)) :
    trace (productPhase U S M M γ ^ m) = 0 ∧
      trace ((U : A) ^ M * productPhase U S M M γ ^ m) =
        algebraMap ℚ K (mixedBetaMoment M (γ * m) 0) := by
  have hγm := Nat.mul_le_mul_right m hγ
  have hMm := Nat.le_mul_of_pos_right M hm
  have hu : ((M * m : ℕ) - (M + 2 * 0 : ℕ) : ℤ) = ↑(M * m - M : ℕ) := by
    simp only [mul_zero, add_zero]
    exact (Int.ofNat_sub hMm).symm
  have hpow : trace (productPhase U S M M γ ^ m) =
      trace ((↑(U ^ (0 : ℤ)) : A) * S ^ 0 * productPhase U S M M γ ^ m) := by simp
  constructor
  · rw [hpow, product_mixedMoment trace U S M M M γ m 0 hM 0 (by omega) weight height,
      if_pos (by positivity : (0 : ℤ) ≤ (M * m : ℕ)), hu, sub_zero, Int.toNat_natCast,
      Ring.choose_natCast, Nat.choose_eq_zero_of_lt (show M * m - M < M * m by omega)]
    simp
  · have hz : (U : A) ^ M = (↑(U ^ (M : ℤ)) : A) := by simp
    rw [hz, ← mul_one (↑(U ^ (M : ℤ)) : A), ← pow_zero S,
      product_mixedMoment trace U S M M M γ m 0 hM M (by omega) weight height,
      if_pos (by exact_mod_cast hMm), hu, ← Int.ofNat_sub hMm, Int.toNat_natCast,
      Ring.choose_natCast, Nat.choose_self]
    simp

/-- The squared normalizer is one instance of the general extraction. -/
theorem mixedMoment_of_weights (trace : A →ₗ[K] K) (U : Aˣ) (S : A) (M k : ℕ) (ℓ : ℤ)
    (weight : ∀ (a : ℤ) (j : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ j) = 0)
    (height : ∀ j : ℕ, trace (S ^ j) = (M : K) / (↑(M + 2 * j) : K)) :
    trace ((↑(U ^ ℓ) : A) * normalizerPhase U S ^ k) =
      if ℓ ≤ (k : ℤ) then
        algebraMap ℚ K ((endpointPolynomial M k).coeff ((k : ℤ) - ℓ).toNat)
      else 0 := by
  by_cases hM : 0 < M
  · simpa [productPhase, normalizerPhase, endpointPolynomial, Nat.ne_of_gt hM] using
      productMoment_of_weights trace U S M 1 1 2 k 0 hM ℓ weight height
  · have hall (j : ℕ) : trace (S ^ j) = 0 := by simp [height, show M = 0 by omega]
    have he := productPhase_pow (K := K) U S 1 1 2 k
    have hphase : normalizerPhase U S = productPhase U S 1 1 2 := by
      simp [normalizerPhase, productPhase]
    rw [hphase, he, ← mul_assoc, ← Units.val_mul, ← zpow_add,
      show ℓ + -(1 * k : ℕ) = -((k : ℤ) - ℓ) by push_cast; ring]
    simp_rw [Finset.mul_sum, map_sum, ← mul_assoc, trace_evalAtUnit trace U S weight, hall,
      mul_zero, ite_self, Finset.sum_const_zero]
    have hpoly : endpointPolynomial M k = 0 := by
      simp [endpointPolynomial, powerEndpointPolynomial_eq_sum, productTerm, show M = 0 by omega]
    simp [hpoly]

/-- The generalized binomial of equation (3), including the independent
vanishing for ℓ > k. -/
theorem normalizer_mixedMoment (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
    (M k : ℕ) (hM : 0 < M) (ℓ : ℤ) (hℓ : -(k : ℤ) ≤ ℓ)
    (weight : ∀ (a : ℤ) (j : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ j) = 0)
    (height : ∀ j : ℕ, trace (S ^ j) = (M : K) / (↑(M + 2 * j) : K)) :
    trace ((↑(U ^ ℓ) : A) * normalizerPhase U S ^ k) =
      if ℓ ≤ (k : ℤ) then algebraMap ℚ K
        (betaMoment M k * ((Ring.choose ((k : ℤ) - M) ((k : ℤ) - ℓ).toNat : ℤ) : ℚ))
      else 0 := by
  rw [mixedMoment_of_weights trace U S M k ℓ weight height]
  split_ifs
  · rw [endpointPolynomial_coeff M k _ hM (by omega)]
  · rfl

end MathieuQuadratic
