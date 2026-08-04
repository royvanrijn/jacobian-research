/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.Algebra.Polynomial.Coeff
import Mathlib.Data.Nat.Factorial.DoubleFactorial
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

/-- The kernel after the endpoint-flat primitive has been written as a
constant plus its order-`2m+1` tail. -/
noncomputable def endpointKernel
    {K : Type*} [Field K] (m : ℕ) (c : K) (tail : K[X]) : K[X] :=
  (1 + X) ^ (m - 1) * (C c + X ^ (2 * m + 1) * tail)

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
