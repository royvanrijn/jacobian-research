/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.Algebra.Polynomial.Coeff
import Mathlib.Algebra.Polynomial.Derivative
import Mathlib.Algebra.Polynomial.Div
import Mathlib.Data.Nat.Factorial.DoubleFactorial
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.Positivity

/-!
# The adjacent-coefficient mechanism

The three-variable counterexample reduces to the following all-order
polynomial fact.  Endpoint contact gives

`J_m(1 + u) = c_m + u^(2m+1) R_m(u)`.

After multiplication by `(1+u)^(m-1)`, the coefficient of `u^m` is zero
while the neighboring coefficient of `u^(m-1)` is `c_m`.
-/

namespace GVC

open Polynomial

/-- The coefficientwise primitive normalized to have constant term zero. -/
noncomputable def polynomialPrimitive (p : ℚ[X]) : ℚ[X] :=
  p.sum fun n a ↦ monomial (n + 1) (a / (n + 1))

theorem polynomialPrimitive_coeff_succ (p : ℚ[X]) (n : ℕ) :
    (polynomialPrimitive p).coeff (n + 1) = p.coeff n / (n + 1) := by
  classical
  rw [polynomialPrimitive, coeff_sum, sum_def, Finset.sum_eq_single n]
  · rw [coeff_monomial, if_pos rfl]
  · intro i _hi hin
    rw [coeff_monomial]
    simp [hin]
  · intro h
    simp [Polynomial.notMem_support_iff.mp h]

@[simp] theorem polynomialPrimitive_coeff_zero (p : ℚ[X]) :
    (polynomialPrimitive p).coeff 0 = 0 := by
  classical
  simp [polynomialPrimitive, sum_def, coeff_monomial]

theorem derivative_polynomialPrimitive (p : ℚ[X]) :
    derivative (polynomialPrimitive p) = p := by
  classical
  ext n
  rw [coeff_derivative, polynomialPrimitive_coeff_succ]
  field_simp

/-- Integrating a polynomial divisible by `X^n` raises its order of
vanishing by one. -/
theorem X_pow_succ_dvd_polynomialPrimitive_X_pow_mul
    (n : ℕ) (p : ℚ[X]) :
    X ^ (n + 1) ∣ polynomialPrimitive (X ^ n * p) := by
  rw [X_pow_dvd_iff]
  intro d hd
  cases d with
  | zero => exact polynomialPrimitive_coeff_zero _
  | succ d =>
      rw [polynomialPrimitive_coeff_succ, coeff_X_pow_mul']
      have hdn : ¬n ≤ d := by omega
      simp [hdn]

/-- The actual tail obtained by integrating the endpoint derivative
`u^(2m) (2+u)^(2m)`. -/
noncomputable def endpointPrimitiveTail (m : ℕ) : ℚ[X] :=
  Classical.choose
    (X_pow_succ_dvd_polynomialPrimitive_X_pow_mul
      (2 * m) ((2 + X) ^ (2 * m)))

theorem polynomialPrimitive_endpoint_eq (m : ℕ) :
    polynomialPrimitive (X ^ (2 * m) * (2 + X) ^ (2 * m)) =
      X ^ (2 * m + 1) * endpointPrimitiveTail m :=
  Classical.choose_spec
    (X_pow_succ_dvd_polynomialPrimitive_X_pow_mul
      (2 * m) ((2 + X) ^ (2 * m)))

/-- The tail for an arbitrary rational profile whose endpoint derivative has
an explicit factor `X^(2n)`.  This is the local algebra needed by the full
cusp-profile family. -/
noncomputable def profilePrimitiveTail (n : ℕ) (profile : ℚ[X]) : ℚ[X] :=
  Classical.choose
    (X_pow_succ_dvd_polynomialPrimitive_X_pow_mul (2 * n) profile)

theorem polynomialPrimitive_profile_eq (n : ℕ) (profile : ℚ[X]) :
    polynomialPrimitive (X ^ (2 * n) * profile) =
      X ^ (2 * n + 1) * profilePrimitiveTail n profile :=
  Classical.choose_spec
    (X_pow_succ_dvd_polynomialPrimitive_X_pow_mul (2 * n) profile)

/-- A shifted primitive with arbitrary constant moment and arbitrary endpoint
profile. -/
noncomputable def shiftedProfilePrimitive
    (n : ℕ) (c : ℚ) (profile : ℚ[X]) : ℚ[X] :=
  C c + polynomialPrimitive (X ^ (2 * n) * profile)

theorem shiftedProfilePrimitive_flat
    (n : ℕ) (c : ℚ) (profile : ℚ[X]) :
    shiftedProfilePrimitive n c profile =
      C c + X ^ (2 * n + 1) * profilePrimitiveTail n profile := by
  rw [shiftedProfilePrimitive, polynomialPrimitive_profile_eq]

theorem derivative_shiftedProfilePrimitive
    (n : ℕ) (c : ℚ) (profile : ℚ[X]) :
    derivative (shiftedProfilePrimitive n c profile) =
      X ^ (2 * n) * profile := by
  rw [shiftedProfilePrimitive, map_add, derivative_C,
    zero_add, derivative_polynomialPrimitive]

/-- The explicit nonzero height coefficient appearing in the minimal
counterexample. -/
noncomputable def cuspMoment (m : ℕ) : ℚ :=
  ((2 : ℚ) ^ (2 * m) * (Nat.factorial (2 * m) : ℚ)) /
    (Nat.doubleFactorial (4 * m + 1) : ℚ)

theorem cuspMoment_pos (m : ℕ) : 0 < cuspMoment m := by
  unfold cuspMoment
  positivity

theorem cuspMoment_ne_zero (m : ℕ) : cuspMoment m ≠ 0 :=
  ne_of_gt (cuspMoment_pos m)

/-- The shifted primitive `J_m(1+u)`, normalized by its value `c_m` at
`u=0`. -/
noncomputable def shiftedEndpointPrimitive (m : ℕ) : ℚ[X] :=
  C (cuspMoment m) +
    polynomialPrimitive (X ^ (2 * m) * (2 + X) ^ (2 * m))

theorem shiftedEndpointPrimitive_flat (m : ℕ) :
    shiftedEndpointPrimitive m =
      C (cuspMoment m) + X ^ (2 * m + 1) * endpointPrimitiveTail m := by
  rw [shiftedEndpointPrimitive, polynomialPrimitive_endpoint_eq]

theorem derivative_shiftedEndpointPrimitive (m : ℕ) :
    derivative (shiftedEndpointPrimitive m) =
      X ^ (2 * m) * (2 + X) ^ (2 * m) := by
  rw [shiftedEndpointPrimitive, map_add, derivative_C,
    zero_add, derivative_polynomialPrimitive]

theorem endpoint_derivative_factorization (m : ℕ) :
    (1 - (1 + X) ^ 2 : ℚ[X]) ^ (2 * m) =
      X ^ (2 * m) * (2 + X) ^ (2 * m) := by
  have hbase : (1 - (1 + X) ^ 2 : ℚ[X]) = -(X * (2 + X)) := by
    ring
  rw [hbase]
  calc
    (-(X * (2 + X))) ^ (2 * m) =
        ((-(X * (2 + X))) ^ 2) ^ m := by rw [pow_mul]
    _ = ((X * (2 + X)) ^ 2) ^ m := by
      congr 1
      ring
    _ = (X ^ 2 * (2 + X) ^ 2) ^ m := by
      congr 1
      ring
    _ = X ^ (2 * m) * (2 + X) ^ (2 * m) := by
      rw [mul_pow, ← pow_mul, ← pow_mul]

/-- The kernel after the endpoint-flat primitive has been written as a
constant plus its order-`2m+1` tail. -/
noncomputable def endpointKernel
    {K : Type*} [Field K] (m : ℕ) (c : K) (tail : K[X]) : K[X] :=
  (1 + X) ^ (m - 1) * (C c + X ^ (2 * m + 1) * tail)

theorem endpointKernel_actual_eq (m : ℕ) :
    endpointKernel m (cuspMoment m) (endpointPrimitiveTail m) =
      (1 + X) ^ (m - 1) * shiftedEndpointPrimitive m := by
  rw [endpointKernel, shiftedEndpointPrimitive_flat]

/-- The endpoint kernel generated by an arbitrary rational profile. -/
noncomputable def profileEndpointKernel
    (n : ℕ) (c : ℚ) (profile : ℚ[X]) : ℚ[X] :=
  (1 + X) ^ (n - 1) * shiftedProfilePrimitive n c profile

theorem profileEndpointKernel_eq
    (n : ℕ) (c : ℚ) (profile : ℚ[X]) :
    profileEndpointKernel n c profile =
      endpointKernel n c (profilePrimitiveTail n profile) := by
  rw [profileEndpointKernel, shiftedProfilePrimitive_flat, endpointKernel]

private theorem coeff_endpoint_error
    {K : Type*} [Field K] (m d : ℕ) (tail : K[X])
    (hd : d < 2 * m + 1) :
    (((1 + X) ^ (m - 1)) * (X ^ (2 * m + 1) * tail)).coeff d = 0 := by
  rw [← mul_assoc, mul_comm ((1 + X) ^ (m - 1)) (X ^ (2 * m + 1)),
    mul_assoc, coeff_X_pow_mul']
  simp [Nat.not_le.mpr hd]

/-- The pure phase coefficient is just beyond the binomial endpoint. -/
theorem endpointKernel_coeff_pure
    {K : Type*} [Field K] (m : ℕ) (hm : 0 < m)
    (c : K) (tail : K[X]) :
    (endpointKernel m c tail).coeff m = 0 := by
  have hlt : m - 1 < m := by omega
  rw [endpointKernel, mul_add, coeff_add]
  rw [coeff_endpoint_error m m tail (by omega)]
  rw [coeff_mul_C, Polynomial.coeff_one_add_X_pow,
    Nat.choose_eq_zero_of_lt hlt]
  simp

/-- Multiplication by the phase variable shifts extraction back to the last
nonzero binomial coefficient. -/
theorem endpointKernel_coeff_mixed
    {K : Type*} [Field K] (m : ℕ) (hm : 0 < m)
    (c : K) (tail : K[X]) :
    (endpointKernel m c tail).coeff (m - 1) = c := by
  rw [endpointKernel, mul_add, coeff_add]
  rw [coeff_endpoint_error m (m - 1) tail (by omega)]
  simp [Polynomial.coeff_one_add_X_pow]

/-- The complete neighboring-coefficient ladder.  In the notation of the
cusp-profile family, instantiate `m` here with `r * M`: multiplying by the
phase variable to the `ell`-th power changes extraction from coefficient
`m` to coefficient `m - ell`, producing the binomial factor in the paper. -/
theorem endpointKernel_coeff_ladder
    {K : Type*} [Field K] (m ell : ℕ)
    (hell : 1 ≤ ell) (helm : ell ≤ m)
    (c : K) (tail : K[X]) :
    (endpointKernel m c tail).coeff (m - ell) =
      (Nat.choose (m - 1) (ell - 1) : K) * c := by
  rw [endpointKernel, mul_add, coeff_add]
  rw [coeff_endpoint_error m (m - ell) tail (by omega)]
  rw [coeff_mul_C, Polynomial.coeff_one_add_X_pow]
  have hsplit : m - 1 = (m - ell) + (ell - 1) := by omega
  rw [Nat.choose_symm_of_eq_add hsplit]
  simp

/-- The full ladder for the explicit endpoint primitive used in the concrete
counterexample. -/
theorem endpointKernel_actual_coeff_ladder
    (m ell : ℕ) (hell : 1 ≤ ell) (helm : ell ≤ m) :
    (endpointKernel m (cuspMoment m) (endpointPrimitiveTail m)).coeff
        (m - ell) =
      (Nat.choose (m - 1) (ell - 1) : ℚ) * cuspMoment m :=
  endpointKernel_coeff_ladder m ell hell helm _ _

theorem endpointKernel_actual_coeff_ladder_ne_zero
    (m ell : ℕ) (hell : 1 ≤ ell) (helm : ell ≤ m) :
    (endpointKernel m (cuspMoment m) (endpointPrimitiveTail m)).coeff
        (m - ell) ≠ 0 := by
  rw [endpointKernel_actual_coeff_ladder m ell hell helm]
  exact mul_ne_zero (by exact_mod_cast Nat.choose_ne_zero (by omega))
    (cuspMoment_ne_zero m)

theorem profileEndpointKernel_coeff_pure
    (n : ℕ) (hn : 0 < n) (c : ℚ) (profile : ℚ[X]) :
    (profileEndpointKernel n c profile).coeff n = 0 := by
  rw [profileEndpointKernel_eq]
  exact endpointKernel_coeff_pure n hn c _

/-- Complete local coefficient calculation for the paper's cusp-profile
family.  Taking `n = r * m` gives precisely its binomial ladder. -/
theorem profileEndpointKernel_coeff_ladder
    (n ell : ℕ) (hell : 1 ≤ ell) (heln : ell ≤ n)
    (c : ℚ) (profile : ℚ[X]) :
    (profileEndpointKernel n c profile).coeff (n - ell) =
      (Nat.choose (n - 1) (ell - 1) : ℚ) * c := by
  rw [profileEndpointKernel_eq]
  exact endpointKernel_coeff_ladder n ell hell heln c _

theorem profileEndpointKernel_coeff_ladder_ne_zero
    (n ell : ℕ) (hell : 1 ≤ ell) (heln : ell ≤ n)
    (c : ℚ) (hc : c ≠ 0) (profile : ℚ[X]) :
    (profileEndpointKernel n c profile).coeff (n - ell) ≠ 0 := by
  rw [profileEndpointKernel_coeff_ladder n ell hell heln c profile]
  exact mul_ne_zero (by exact_mod_cast Nat.choose_ne_zero (by omega)) hc

/-- The exact scalar obtained after one additional Laplacian. -/
noncomputable def mixedDerivativeValue (m : ℕ) : ℚ :=
  (2 : ℚ) ^ (8 * m + 1) * Nat.factorial (6 * m + 1) *
    Nat.factorial (2 * m) * Nat.doubleFactorial (12 * m + 3) /
      Nat.doubleFactorial (4 * m + 1)

theorem mixedDerivativeValue_pos (m : ℕ) :
    0 < mixedDerivativeValue m := by
  unfold mixedDerivativeValue
  positivity

theorem mixedDerivativeValue_ne_zero (m : ℕ) :
    mixedDerivativeValue m ≠ 0 :=
  ne_of_gt (mixedDerivativeValue_pos m)

end GVC
