import QuadraticNorms.MinusEightThetaNormalization
import Mathlib.Algebra.Polynomial.Expand

open scoped BigOperators ArithmeticFunction.sigma
open Polynomial

namespace QuadraticNorms

private noncomputable def minusEightEtaProduct (k : ℕ) : ℤ[X] :=
  ∏ j ∈ Finset.range k, (1 - X ^ (j + 1)) ^ 2 * (1 - X ^ (2 * (j + 1))) ^ 2

private theorem minusEight_eta_product_coeff_stable {n k : ℕ} (h : n ≤ k) :
    (minusEightEtaProduct k).coeff n = (minusEightEtaProduct n).coeff n := by
  induction k, h using Nat.le_induction with
  | base => rfl
  | succ k hk ih =>
    let y : ℤ[X] := X ^ (k + 1)
    let v : ℤ[X] := -2 - y + 4 * y ^ 2 - y ^ 3 - 2 * y ^ 4 + y ^ 5
    have hfactor : (1 - X ^ (k + 1)) ^ 2 * (1 - X ^ (2 * (k + 1))) ^ 2 =
        1 + y * v := by
      rw [mul_comm 2 (k + 1), pow_mul]
      change (1 - y) ^ 2 * (1 - y ^ 2) ^ 2 = 1 + y * v
      dsimp [v]
      ring
    have hprod : minusEightEtaProduct (k + 1) =
        minusEightEtaProduct k + (minusEightEtaProduct k * v) * X ^ (k + 1) := by
      unfold minusEightEtaProduct
      rw [Finset.prod_range_succ, hfactor]
      dsimp [y]
      ring
    rw [hprod, coeff_add, coeff_mul_X_pow', if_neg (by omega), add_zero, ih]

/-- The coefficient of qⁿ in the level-32 eta product, obtained from
the stable finite Euler product after writing the nonconstant factors in q⁴. -/
noncomputable def minusEightEtaCoefficient : ℕ → ℤ
  | 0 => 0
  | n + 1 => if 4 ∣ n then (minusEightEtaProduct (n / 4)).coeff (n / 4) else 0

/-- Every cutoff at least as large as the requested degree recovers the exact
coefficient in the shifted level-32 product. -/
theorem minusEight_eta_coefficient_eq_truncation {n k : ℕ} (h : n ≤ k) :
    minusEightEtaCoefficient n = (X * expand ℤ 4 (minusEightEtaProduct k)).coeff n := by
  cases n with
  | zero => simp [minusEightEtaCoefficient]
  | succ n =>
    rw [minusEightEtaCoefficient, coeff_X_mul, coeff_expand (by decide)]
    rw [minusEight_eta_product_coeff_stable ((Nat.div_le_self n 4).trans ((Nat.le_succ n).trans h))]

theorem minusEight_eta_even_zero {n : ℕ} (hn : 2 ∣ n) : minusEightEtaCoefficient n = 0 := by
  cases n with
  | zero => rfl
  | succ n =>
    rw [minusEightEtaCoefficient, if_neg]
    omega

private theorem minusEight_eta_five : minusEightEtaCoefficient 5 = -2 := by
  have h : (((1 - (X : ℤ[X])) ^ 2) * (1 - X ^ 2) ^ 2) =
      1 - 2 * X - X ^ 2 + 4 * X ^ 3 - X ^ 4 - 2 * X ^ 5 + X ^ 6 := by ring
  change ((minusEightEtaProduct 1).coeff 1) = -2
  rw [minusEightEtaProduct, Finset.prod_range_one]
  norm_num only [Nat.zero_add, pow_one, Nat.mul_one]
  rw [h]
  norm_num [coeff_add, coeff_sub, coeff_X, coeff_X_pow, coeff_one, ← C_ofNat, coeff_C_mul]

private theorem minusEight_eta_small_table :
    ∀ n : Fin 9, minusEightEtaCoefficient n.val = if n.val = 1 then 1 else if n.val = 5 then -2 else 0 := by
  intro n
  have hn : n.val ≤ 8 := by omega
  generalize hv : n.val = j at *
  interval_cases j <;> first
    | exact minusEight_eta_five
    | simp [minusEightEtaCoefficient, minusEightEtaProduct]

/-- The two quaternary theta shell counts through the respective
Sturm cutoffs, including the constant terms. -/
theorem minusEight_theta_shell_small_tables :
    (∀ n : Fin 3, (minusEightFullShell n.val).card = [1, 4, 8][n.val]!) ∧
    (∀ n : Fin 9, (minusEightWeightedShell 4 n.val).card =
      [1, 2, 2, 4, 4, 4, 8, 8, 8][n.val]!) := by
  decide +kernel

private def sigmaQuotientInt (n d : ℕ) : ℤ :=
  if d ∣ n then (σ 1 (n / d) : ℤ) else 0

/-- The proposed level-8 and level-32 right-hand sides agree with
the shell counts through the entire coefficient checks in the manuscript. -/
theorem minusEight_theta_coefficients_through_sturm :
    (∀ n : Fin 3, 0 < n.val → ((minusEightFullShell n.val).card : ℤ) =
      4 * sigmaQuotientInt n.val 1 - 4 * sigmaQuotientInt n.val 2 +
      8 * sigmaQuotientInt n.val 4 - 32 * sigmaQuotientInt n.val 8) ∧
    (∀ n : Fin 9, 0 < n.val → ((minusEightWeightedShell 4 n.val).card : ℤ) =
      sigmaQuotientInt n.val 1 - sigmaQuotientInt n.val 2 +
      8 * sigmaQuotientInt n.val 16 - 32 * sigmaQuotientInt n.val 32 + minusEightEtaCoefficient n.val) := by
  constructor
  · decide +kernel
  · intro n hn
    rw [minusEight_eta_small_table n]
    exact (by decide +kernel : ∀ n : Fin 9, 0 < n.val →
      ((minusEightWeightedShell 4 n.val).card : ℤ) =
        sigmaQuotientInt n.val 1 - sigmaQuotientInt n.val 2 +
        8 * sigmaQuotientInt n.val 16 - 32 * sigmaQuotientInt n.val 32 +
        (if n.val = 1 then 1 else if n.val = 5 then -2 else 0)) n hn

end QuadraticNorms
