import MathieuQuadratic.Endpoint
import MathieuQuadratic.Witness

/-!
# From the weight and height identities to the mixed moments

The trace is on a commutative K-algebra with an invertible norm U and a
height S. The hypotheses here are the scalar weight calculation and the
rank-two Lie recurrence in the uniform quadratic-complement proof.
-/

open Polynomial
open scoped BigOperators

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

private noncomputable def normalizerTerm (k j : ℕ) : ℚ[X] :=
  C (((2 * k).choose j : ℚ) * (-1) ^ j) * (X + 1) ^ (k + 2 * j)

private theorem normalizerPhase_pow (U : Aˣ) (S : A) (k : ℕ) :
    normalizerPhase U S ^ k =
      (↑(U ^ (-(k : ℤ))) : A) *
        ∑ j ∈ Finset.range (2 * k + 1), evalAtUnit K U (normalizerTerm k j) * S ^ j := by
  have expand : (1 - (1 + (U : A)) ^ 2 * S) ^ (2 * k) =
      ∑ j ∈ Finset.range (2 * k + 1),
        (-((1 + (U : A)) ^ 2 * S)) ^ j * ((2 * k).choose j : A) := by
    rw [sub_eq_add_neg, add_comm 1, add_pow]
    simp only [one_pow, mul_one]
  rw [normalizerPhase, mul_pow, mul_pow, ← pow_mul (1 - _), zpow_neg, zpow_natCast,
    ← inv_pow, Units.val_pow_eq_pow_val, mul_assoc, expand, Finset.mul_sum]
  congr 1
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  change _ = (normalizerTerm k j).eval₂
    ((algebraMap K A).comp (algebraMap ℚ K)) (U : A) * S ^ j
  rw [normalizerTerm, eval₂_mul, eval₂_C, eval₂_pow, eval₂_add, eval₂_X, eval₂_one]
  simp only [map_mul, map_natCast, map_pow, map_neg, map_one]
  rw [neg_pow, mul_pow ((1 + (U : A)) ^ 2) S, ← pow_mul, pow_add]
  ring

private theorem endpointPolynomial_eq_sum (M k : ℕ) :
    endpointPolynomial M k =
      ∑ j ∈ Finset.range (2 * k + 1),
        C ((M : ℚ) / (↑(M + 2 * j) : ℚ)) * normalizerTerm k j := by
  rw [endpointPolynomial, betaQuotient, sum_comp, Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  simp only [mul_comp, C_comp, pow_comp, X_comp,
    normalizerTerm, betaCoefficient, map_mul, pow_add]
  ring

/-- Coefficient extraction before endpoint truncation. The negative-index
case gives the independent scalar-weight vanishing for ℓ > k. -/
theorem mixedMoment_of_weights (trace : A →ₗ[K] K) (U : Aˣ) (S : A) (M k : ℕ) (ℓ : ℤ)
    (weight : ∀ (a : ℤ) (j : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ j) = 0)
    (height : ∀ j : ℕ, trace (S ^ j) = (M : K) / (↑(M + 2 * j) : K)) :
    trace ((↑(U ^ ℓ) : A) * normalizerPhase U S ^ k) =
      if ℓ ≤ (k : ℤ) then
        algebraMap ℚ K ((endpointPolynomial M k).coeff ((k : ℤ) - ℓ).toNat)
      else 0 := by
  rw [normalizerPhase_pow (K := K), ← mul_assoc, ← Units.val_mul, ← zpow_add,
    show ℓ + -(k : ℤ) = -((k : ℤ) - ℓ) by omega]
  simp_rw [Finset.mul_sum, map_sum, ← mul_assoc, trace_evalAtUnit trace U S weight,
    height, sub_nonneg]
  by_cases hℓ : ℓ ≤ (k : ℤ)
  · simp only [hℓ, if_true, endpointPolynomial_eq_sum, finsetSum_coeff, map_sum]
    refine Finset.sum_congr rfl fun j _ ↦ ?_
    rw [coeff_C_mul, map_mul, map_div₀, map_natCast, map_natCast, mul_comm]
  · simp [hℓ]

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
