/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GVC.ConcreteWitness

/-!
# Explicit ternary Reynolds expansion

This file expands the algebraic Reynolds functional into the finite diagonal
coefficient sum selected by `(4XY + T²)^k`.  It is the first monomialwise
part of the remaining quadric phase extraction.
-/

namespace GVC

open MvPolynomial

/-- The exponent `(j,j,2(k-j))` selected by the `j`-th binomial term of
`Delta^k`. -/
noncomputable def reynoldsIndex (k j : ℕ) : Fin 3 →₀ ℕ :=
  Finsupp.single 0 j + Finsupp.single 1 j +
    Finsupp.single 2 (2 * (k - j))

@[simp] theorem reynoldsIndex_zero (k j : ℕ) :
    reynoldsIndex k j 0 = j := by
  simp [reynoldsIndex]

@[simp] theorem reynoldsIndex_one (k j : ℕ) :
    reynoldsIndex k j 1 = j := by
  simp [reynoldsIndex]

@[simp] theorem reynoldsIndex_two (k j : ℕ) :
    reynoldsIndex k j 2 = 2 * (k - j) := by
  simp [reynoldsIndex]

theorem gvcDelta_pow_expansion (k : ℕ) :
    gvcDelta ^ k =
      ∑ j ∈ Finset.range (k + 1),
        monomial (reynoldsIndex k j)
          ((4 : ℚ) ^ j * Nat.choose k j) := by
  have hxy : 4 * gvcX * gvcY =
      monomial (Finsupp.single 0 1 + Finsupp.single 1 1) (4 : ℚ) := by
    have hfour : (4 : TernaryPolynomial) = C (4 : ℚ) :=
      (C_eq_coe_nat (R := ℚ) 4).symm
    rw [hfour, gvcX, gvcY, C_mul_X_eq_monomial,
      ← pow_one (X (1 : Fin 3)), ← monomial_add_single]
  have ht : gvcT ^ 2 =
      monomial (Finsupp.single 2 2) (1 : ℚ) := by
    simp [gvcT, X_pow_eq_monomial]
  rw [gvcDelta, hxy, ht, add_pow]
  apply Finset.sum_congr rfl
  intro j hj
  simp only [Finset.mem_range] at hj
  rw [monomial_pow, monomial_pow, monomial_mul,
    ← C_eq_coe_nat (R := ℚ), mul_comm, C_mul_monomial]
  apply congrArg₂ (fun α a ↦ (monomial α a : TernaryPolynomial))
  · ext i
    fin_cases i
    · simp [reynoldsIndex]
    · simp [reynoldsIndex]
    · simp [reynoldsIndex]
      omega
  · ring_nf

theorem multiDescFactorial_reynoldsIndex_self (k j : ℕ) :
    multiDescFactorial (reynoldsIndex k j) (reynoldsIndex k j) =
      Nat.factorial j * Nat.factorial j *
        Nat.factorial (2 * (k - j)) := by
  rw [multiDescFactorial_self]
  simp [Fin.prod_univ_succ, mul_assoc]

/-- The unnormalized top contraction is the explicit finite diagonal
coefficient sum. -/
theorem topContraction_delta_pow_eq_sum (k : ℕ) (p : TernaryPolynomial) :
    topContraction (gvcDelta ^ k) p =
      ∑ j ∈ Finset.range (k + 1),
        ((4 : ℚ) ^ j * Nat.choose k j) *
          p.coeff (reynoldsIndex k j) *
          (Nat.factorial j * Nat.factorial j *
            Nat.factorial (2 * (k - j))) := by
  classical
  rw [gvcDelta_pow_expansion]
  generalize Finset.range (k + 1) = s
  induction s using Finset.induction_on with
  | empty => simp
  | @insert j s hj ih =>
      rw [Finset.sum_insert hj, topContraction_add_left, ih,
        Finset.sum_insert hj, topContraction_monomial_left,
        multiDescFactorial_reynoldsIndex_self]
      push_cast
      ring

theorem algebraicReynoldsMoment_eq_diagonal_sum
    (k : ℕ) (p : TernaryPolynomial) :
    algebraicReynoldsMoment k p =
      (∑ j ∈ Finset.range (k + 1),
        ((4 : ℚ) ^ j * Nat.choose k j) *
          p.coeff (reynoldsIndex k j) *
          (Nat.factorial j * Nat.factorial j *
            Nat.factorial (2 * (k - j)))) /
        reynoldsScale k := by
  rw [algebraicReynoldsMoment, topContraction_delta_pow_eq_sum]

section PolynomialIntegral

open Polynomial

/-- Algebraic integration of a monomial on `[0,1]`. -/
theorem formalIntegral01_X_pow (n : ℕ) :
    formalIntegral01 (Polynomial.X ^ n) = 1 / (n + 1 : ℚ) := by
  have h := formalIntegral01_derivative (Polynomial.X ^ (n + 1))
  rw [derivative_X_pow_succ, formalIntegral01_C_mul] at h
  simp at h
  have hn : (n + 1 : ℚ) ≠ 0 := by positivity
  apply (eq_div_iff hn).2
  nlinarith

/-- Integration by parts for the two-parameter beta monomials. -/
theorem formalIntegral01_beta_recurrence (j r : ℕ) :
    ((2 * r + 1 : ℕ) : ℚ) *
        formalIntegral01
          ((1 - Polynomial.X ^ 2) ^ (j + 1) *
            Polynomial.X ^ (2 * r)) =
      ((2 * (j + 1) : ℕ) : ℚ) *
        formalIntegral01
          ((1 - Polynomial.X ^ 2) ^ j *
            Polynomial.X ^ (2 * (r + 1))) := by
  let f : ℚ[X] := 1 - Polynomial.X ^ 2
  have hjcoeff :
      Polynomial.C ((2 * (j + 1) : ℕ) : ℚ) =
        Polynomial.C ((j + 1 : ℕ) : ℚ) * Polynomial.C 2 := by
    rw [← Polynomial.C_mul]
    congr 1
    push_cast
    ring
  have hderivative :
      derivative (Polynomial.X ^ (2 * r + 1) * f ^ (j + 1)) =
        Polynomial.C ((2 * r + 1 : ℕ) : ℚ) *
            (f ^ (j + 1) * Polynomial.X ^ (2 * r)) -
          Polynomial.C ((2 * (j + 1) : ℕ) : ℚ) *
            (f ^ j * Polynomial.X ^ (2 * (r + 1))) := by
    dsimp [f]
    rw [derivative_mul, derivative_X_pow_succ, derivative_pow_succ,
      derivative_sub, derivative_one, derivative_X_sq, zero_sub, hjcoeff]
    push_cast
    ring
  have hboundary :
      formalIntegral01
        (derivative (Polynomial.X ^ (2 * r + 1) * f ^ (j + 1))) = 0 := by
    rw [formalIntegral01_derivative]
    simp [f]
  rw [hderivative, formalIntegral01_sub, formalIntegral01_C_mul,
    formalIntegral01_C_mul] at hboundary
  push_cast at hboundary ⊢
  nlinarith

/-- Closed beta evaluation for the diagonal monomials selected by the
ternary Reynolds functional. -/
theorem formalIntegral01_beta_monomial (j r : ℕ) :
    formalIntegral01
        ((1 - Polynomial.X ^ 2) ^ j * Polynomial.X ^ (2 * r)) =
      ((2 : ℚ) ^ j * Nat.factorial j *
          Nat.doubleFactorial (2 * r - 1)) /
        Nat.doubleFactorial (2 * (r + j) + 1) := by
  induction j generalizing r with
  | zero =>
      rw [pow_zero, one_mul, formalIntegral01_X_pow]
      have hdf := Nat.doubleFactorial_add_one (2 * r)
      have hdenom :
          Nat.doubleFactorial (2 * r + 1) =
            (2 * r + 1) * Nat.doubleFactorial (2 * r - 1) := by
        simpa using hdf
      rw [show 2 * (r + 0) + 1 = 2 * r + 1 by omega, hdenom]
      push_cast
      have hodd : (2 * (r : ℚ) + 1) ≠ 0 := by positivity
      have hdfne : (Nat.doubleFactorial (2 * r - 1) : ℚ) ≠ 0 := by
        positivity
      field_simp [hodd, hdfne]
      norm_num
  | succ j ih =>
      have hrec := formalIntegral01_beta_recurrence j r
      rw [ih (r + 1)] at hrec
      have hdf :
          Nat.doubleFactorial (2 * (r + 1) - 1) =
            (2 * r + 1) * Nat.doubleFactorial (2 * r - 1) := by
        rw [show 2 * (r + 1) - 1 = 2 * r + 1 by omega,
          Nat.doubleFactorial_add_one]
      have hdenom :
          2 * ((r + 1) + j) + 1 = 2 * (r + (j + 1)) + 1 := by
        omega
      rw [hdf, hdenom] at hrec
      have hpow : (2 : ℚ) ^ (j + 1) = 2 ^ j * 2 := by
        rw [pow_succ]
      rw [hpow, Nat.factorial_succ]
      have hscale : ((2 * r + 1 : ℕ) : ℚ) ≠ 0 := by positivity
      have hdenom_ne :
          (Nat.doubleFactorial (2 * (r + (j + 1)) + 1) : ℚ) ≠ 0 := by
        positivity
      push_cast at hrec ⊢
      field_simp [hdenom_ne] at hrec ⊢
      nlinarith

theorem factorial_two_mul_eq_pow_mul_doubleFactorial (r : ℕ) :
    Nat.factorial (2 * r) =
      2 ^ r * Nat.factorial r * Nat.doubleFactorial (2 * r - 1) := by
  induction r with
  | zero => simp
  | succ r ih =>
      rw [show 2 * (r + 1) = (2 * r + 1) + 1 by omega,
        Nat.factorial_succ, Nat.factorial_succ, ih, pow_succ,
        Nat.factorial_succ r,
        show 2 * r + 1 + 1 - 1 = 2 * r + 1 by omega,
        Nat.doubleFactorial_add_one]
      ring

/-- Each normalized diagonal contraction weight is exactly the formal beta
integral of its phase-height monomial. -/
theorem reynolds_diagonal_weight_eq_beta
    (k j : ℕ) (hj : j ≤ k) :
    (((4 : ℚ) ^ j * Nat.choose k j) *
          (Nat.factorial j * Nat.factorial j *
            Nat.factorial (2 * (k - j)))) /
        reynoldsScale k =
      formalIntegral01
        ((1 - Polynomial.X ^ 2) ^ j *
          Polynomial.X ^ (2 * (k - j))) := by
  rw [formalIntegral01_beta_monomial, reynoldsScale]
  have hindex : 2 * ((k - j) + j) + 1 = 2 * k + 1 := by omega
  rw [hindex]
  have hchoose := Nat.choose_mul_factorial_mul_factorial hj
  have heven := factorial_two_mul_eq_pow_mul_doubleFactorial (k - j)
  have hfour : 4 ^ j = 2 ^ (2 * j) := by
    rw [show 4 = 2 ^ 2 by norm_num, ← pow_mul]
  have hpowers :
      2 ^ (2 * j) * 2 ^ (k - j) = 2 ^ k * 2 ^ j := by
    rw [← pow_add, ← pow_add]
    congr 1
    omega
  have hnat :
      4 ^ j * Nat.choose k j *
          Nat.factorial j * Nat.factorial (2 * (k - j)) =
        (2 ^ k * Nat.factorial k) * 2 ^ j *
          Nat.doubleFactorial (2 * (k - j) - 1) := by
    rw [hfour, heven]
    calc
      2 ^ (2 * j) * Nat.choose k j *
          Nat.factorial j *
            (2 ^ (k - j) * Nat.factorial (k - j) *
              Nat.doubleFactorial (2 * (k - j) - 1)) =
          (2 ^ (2 * j) * 2 ^ (k - j)) *
            (Nat.choose k j * Nat.factorial j *
              Nat.factorial (k - j)) *
              Nat.doubleFactorial (2 * (k - j) - 1) := by ring
      _ = (2 ^ k * 2 ^ j) * Nat.factorial k *
              Nat.doubleFactorial (2 * (k - j) - 1) := by
          rw [hpowers, hchoose]
      _ = (2 ^ k * Nat.factorial k) *
            2 ^ j * Nat.doubleFactorial (2 * (k - j) - 1) := by ring
  have hscale :
      (2 : ℚ) ^ k * Nat.factorial k *
        Nat.doubleFactorial (2 * k + 1) ≠ 0 := by positivity
  have hdouble :
      (Nat.doubleFactorial (2 * k + 1) : ℚ) ≠ 0 := by positivity
  field_simp [hscale, hdouble]
  exact_mod_cast hnat

theorem formalIntegral01_finset_sum
    {ι : Type*} (s : Finset ι) (f : ι → ℚ[X]) :
    formalIntegral01 (∑ i ∈ s, f i) =
      ∑ i ∈ s, formalIntegral01 (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp [formalIntegral01, polynomialPrimitive]
  | @insert i s hi ih =>
      simp [hi, ih, formalIntegral01_add]

/-- The phase-height polynomial obtained by extracting equal `x`/`y`
exponents after imposing `xy = 1-t²`. -/
noncomputable def reynoldsPhasePolynomial
    (k : ℕ) (p : TernaryPolynomial) : ℚ[X] :=
  ∑ j ∈ Finset.range (k + 1),
    Polynomial.C (p.coeff (reynoldsIndex k j)) *
      ((1 - Polynomial.X ^ 2) ^ j *
        Polynomial.X ^ (2 * (k - j)))

/-- Fully algebraic Reynolds--phase identity: normalized contraction by
`Delta^k` is formal integration of the diagonal phase-height polynomial. -/
theorem algebraicReynoldsMoment_eq_formalIntegral_phase
    (k : ℕ) (p : TernaryPolynomial) :
    algebraicReynoldsMoment k p =
      formalIntegral01 (reynoldsPhasePolynomial k p) := by
  rw [algebraicReynoldsMoment_eq_diagonal_sum,
    reynoldsPhasePolynomial, formalIntegral01_finset_sum]
  simp_rw [formalIntegral01_C_mul]
  rw [div_eq_mul_inv, Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro j hj
  have hjk : j ≤ k := Nat.le_of_lt_succ (Finset.mem_range.mp hj)
  rw [← reynolds_diagonal_weight_eq_beta k j hjk]
  ring

end PolynomialIntegral

end GVC
