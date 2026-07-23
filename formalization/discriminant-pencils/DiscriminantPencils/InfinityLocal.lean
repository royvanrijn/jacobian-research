/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import DiscriminantPencils.ProjectiveParametrization
import Mathlib.RingTheory.PowerSeries.Inverse

/-!
# Local coordinates at infinity

This file formalizes the coefficient calculation in equation (2.2).  On the
source chart `R = 1`, the local parameter is `q = U/R`.  The homogeneous
coordinate forms become reversed coefficient polynomials in `q`.
-/

open Polynomial

namespace DiscriminantPencils

variable {𝕜 : Type*} [Field 𝕜]

/-- The fixed-degree reversal
`p∞(q) = ∑_{m=0}^n p_{n-m} q^m`. -/
noncomputable def infinityPolynomial (p : 𝕜[X]) (n : ℕ) : 𝕜[X] :=
  ∑ m ∈ Finset.range (n + 1), monomial m (p.coeff (n - m))

@[simp] theorem coeff_infinityPolynomial (p : 𝕜[X]) (n m : ℕ) :
    (infinityPolynomial p n).coeff m =
      if m ≤ n then p.coeff (n - m) else 0 := by
  classical
  simp [infinityPolynomial, coeff_monomial]

/-- Fixed-degree reversal is evaluation of the homogenization on `[1:q]`. -/
theorem infinityPolynomial_eq_aeval_homogenize (p : 𝕜[X]) (n : ℕ) :
    infinityPolynomial p n =
      MvPolynomial.aeval ![1, X] (p.homogenize n) := by
  ext m
  classical
  simp [infinityPolynomial, Polynomial.homogenize,
    MvPolynomial.aeval_monomial, coeff_monomial]
  rw [Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk]
  by_cases hmn : m ≤ n
  · rw [if_pos hmn]
    rw [Finset.sum_eq_single (n - m)]
    · simp [Nat.sub_sub_self hmn]
    · intro b hb hne
      have hb' : b ≤ n := Nat.lt_succ_iff.mp (Finset.mem_range.mp hb)
      simp only
      split_ifs with h
      · omega
      · rfl
    · simp [hmn]
  · rw [if_neg hmn]
    symm
    apply Finset.sum_eq_zero
    intro b hb
    have hb' : b ≤ n := Nat.lt_succ_iff.mp (Finset.mem_range.mp hb)
    simp only
    split_ifs with h
    · omega
    · rfl

/-- The three local numerator polynomials in the `R = 1` source chart. -/
noncomputable def infinityCoordinates (H : 𝕜[X]) (n : ℕ) : Fin 3 → 𝕜[X]
  | 0 => infinityPolynomial H.derivative n
  | 1 => infinityPolynomial (X * H.derivative - H) n
  | 2 => X ^ n

theorem infinityCoordinates_eq_projectiveCoordinates (H : 𝕜[X]) (n : ℕ)
    (i : Fin 3) :
    infinityCoordinates H n i =
      MvPolynomial.aeval ![1, X] (projectiveCoordinates H n i) := by
  fin_cases i
  · exact infinityPolynomial_eq_aeval_homogenize _ _
  · exact infinityPolynomial_eq_aeval_homogenize _ _
  · simp [infinityCoordinates, projectiveCoordinates]

/-- The `S` numerator vanishes at infinity. -/
@[simp] theorem infinity_s_coeff_zero (H : 𝕜[X]) (n : ℕ) :
    (infinityCoordinates H n 0).coeff 0 = H.derivative.coeff n := by
  simp [infinityCoordinates]

/-- Its linear coefficient is `n h_n`. -/
theorem infinity_s_coeff_one (H : 𝕜[X]) (n : ℕ) (hn : 1 ≤ n) :
    (infinityCoordinates H n 0).coeff 1 = (n : 𝕜) * H.coeff n := by
  simp [infinityCoordinates, coeff_derivative, hn]
  ring

/-- The `T` numerator has constant coefficient `(n-1)h_n`. -/
theorem infinity_t_coeff_zero (H : 𝕜[X]) (n : ℕ) (hn : 1 ≤ n) :
    (infinityCoordinates H n 1).coeff 0 =
      (n - 1 : ℕ) * H.coeff n := by
  simp [infinityCoordinates, coeff_X_mul_derivative_sub H n hn]

/-- For exact degree `n ≥ 2`, `T` is a unit in the formal power-series
local ring at `q = 0`. -/
theorem infinity_t_isUnit [CharZero 𝕜] (H : 𝕜[X]) (n : ℕ) (hn : 2 ≤ n)
    (hdeg : H.natDegree = n) :
    IsUnit ((infinityCoordinates H n 1 : 𝕜[X]) : PowerSeries 𝕜) := by
  rw [PowerSeries.isUnit_iff_constantCoeff]
  change IsUnit ((infinityCoordinates H n 1).coeff 0)
  rw [infinity_t_coeff_zero H n (by omega)]
  apply isUnit_iff_ne_zero.mpr
  have hH : H ≠ 0 := by
    intro h
    subst H
    simp at hdeg
    omega
  have hlc : H.coeff n ≠ 0 := by
    simpa [leadingCoeff, hdeg] using leadingCoeff_ne_zero.mpr hH
  exact mul_ne_zero (by exact_mod_cast (show n - 1 ≠ 0 by omega)) hlc

end DiscriminantPencils
