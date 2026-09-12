/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import DiscriminantPencils.ProjectiveParametrization
import Mathlib.RingTheory.PowerSeries.Inverse
import Mathlib.RingTheory.PowerSeries.Order
import Mathlib.RingTheory.PowerSeries.Substitution

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
    · simp
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

/-- The first affine coordinate `S/T` in the target chart `T ≠ 0`. -/
noncomputable def infinitySOverT (H : 𝕜[X]) (n : ℕ) : PowerSeries 𝕜 :=
  (infinityCoordinates H n 0 : PowerSeries 𝕜) *
    (infinityCoordinates H n 1 : PowerSeries 𝕜)⁻¹

/-- The second affine coordinate `Z/T` in the target chart `T ≠ 0`. -/
noncomputable def infinityZOverT (H : 𝕜[X]) (n : ℕ) : PowerSeries 𝕜 :=
  (PowerSeries.X : PowerSeries 𝕜) ^ n *
    (infinityCoordinates H n 1 : PowerSeries 𝕜)⁻¹

theorem exactDegree_coeff_ne_zero (H : 𝕜[X]) (n : ℕ)
    (hdeg : H.natDegree = n) (hn : 0 < n) :
    H.coeff n ≠ 0 := by
  have hH : H ≠ 0 := by
    intro h
    subst H
    simp at hdeg
    omega
  simpa [leadingCoeff, hdeg] using leadingCoeff_ne_zero.mpr hH

/-- Exact degree `n` makes the constant coefficient of the `S` numerator
zero. -/
theorem infinity_s_coeff_zero_of_exactDegree (H : 𝕜[X]) (n : ℕ)
    (hdeg : H.natDegree = n) :
    (infinityCoordinates H n 0).coeff 0 = 0 := by
  rw [infinity_s_coeff_zero]
  rw [coeff_derivative]
  have hcoeff : H.coeff (n + 1) = 0 := by
    apply coeff_eq_zero_of_natDegree_lt
    omega
  simp [hcoeff]

/-- The chart coordinate `S/T` vanishes at the point at infinity. -/
theorem infinitySOverT_coeff_zero (H : 𝕜[X]) (n : ℕ)
    (hdeg : H.natDegree = n) :
    PowerSeries.coeff 0 (infinitySOverT H n) = 0 := by
  simp [infinitySOverT, infinity_s_coeff_zero_of_exactDegree H n hdeg]

/-- The linear coefficient of `S/T` is `n/(n-1)`, the first expansion in
equation (2.2). -/
theorem infinitySOverT_coeff_one [CharZero 𝕜] (H : 𝕜[X]) (n : ℕ)
    (hn : 2 ≤ n) (hdeg : H.natDegree = n) :
    PowerSeries.coeff 1 (infinitySOverT H n) =
      (n : 𝕜) / (n - 1 : ℕ) := by
  rw [infinitySOverT, PowerSeries.coeff_mul]
  rw [Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk]
  simp [Finset.sum_range_succ, infinity_s_coeff_zero_of_exactDegree H n hdeg,
    infinity_s_coeff_one H n (by omega), PowerSeries.constantCoeff_inv,
    infinity_t_coeff_zero H n (by omega)]
  have hlc := exactDegree_coeff_ne_zero H n hdeg (by omega)
  have hn1 : ((n - 1 : ℕ) : 𝕜) ≠ 0 := by
    exact_mod_cast (show n - 1 ≠ 0 by omega)
  field_simp

theorem infinitySOverT_coeff_one_ne_zero [CharZero 𝕜] (H : 𝕜[X]) (n : ℕ)
    (hn : 2 ≤ n) (hdeg : H.natDegree = n) :
    PowerSeries.coeff 1 (infinitySOverT H n) ≠ 0 := by
  rw [infinitySOverT_coeff_one H n hn hdeg]
  apply div_ne_zero
  · exact_mod_cast (show n ≠ 0 by omega)
  · exact_mod_cast (show n - 1 ≠ 0 by omega)

/-- Consequently `S/T` is a uniformizing parameter: it has formal order
exactly one. -/
theorem infinitySOverT_order [CharZero 𝕜] (H : 𝕜[X]) (n : ℕ)
    (hn : 2 ≤ n) (hdeg : H.natDegree = n) :
    PowerSeries.order (infinitySOverT H n) = 1 := by
  change PowerSeries.order (infinitySOverT H n) = (↑(1 : ℕ) : ℕ∞)
  apply PowerSeries.order_eq_nat.mpr
  refine ⟨infinitySOverT_coeff_one_ne_zero H n hn hdeg, ?_⟩
  intro i hi
  have hi0 : i = 0 := by omega
  subst i
  exact infinitySOverT_coeff_zero H n hdeg

theorem infinitySOverT_constantCoeff (H : 𝕜[X]) (n : ℕ)
    (hdeg : H.natDegree = n) :
    (infinitySOverT H n).constantCoeff = 0 := by
  rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]
  exact infinitySOverT_coeff_zero H n hdeg

/-- `Z/T` has no terms below degree `n`. -/
theorem infinityZOverT_coeff_of_lt (H : 𝕜[X]) (n m : ℕ) (hm : m < n) :
    PowerSeries.coeff m (infinityZOverT H n) = 0 := by
  simp [infinityZOverT, PowerSeries.coeff_X_pow_mul', hm.not_ge]

/-- The coefficient of `q^n` in `Z/T` is
`1 / ((n-1)h_n)`, the second expansion in equation (2.2). -/
theorem infinityZOverT_coeff_degree [CharZero 𝕜] (H : 𝕜[X]) (n : ℕ)
    (hn : 2 ≤ n) (hdeg : H.natDegree = n) :
    PowerSeries.coeff n (infinityZOverT H n) =
      ((n - 1 : ℕ) * H.coeff n)⁻¹ := by
  rw [infinityZOverT, PowerSeries.coeff_X_pow_mul']
  simp [PowerSeries.constantCoeff_inv, infinity_t_coeff_zero H n (by omega)]

/-- The second chart coordinate has formal order exactly `n`. -/
theorem infinityZOverT_order [CharZero 𝕜] (H : 𝕜[X]) (n : ℕ)
    (hn : 2 ≤ n) (hdeg : H.natDegree = n) :
    PowerSeries.order (infinityZOverT H n) = n := by
  rw [PowerSeries.order_eq_nat]
  constructor
  · rw [infinityZOverT_coeff_degree H n hn hdeg]
    apply inv_ne_zero
    exact mul_ne_zero
      (by exact_mod_cast (show n - 1 ≠ 0 by omega))
      (exactDegree_coeff_ne_zero H n hdeg (by omega))
  · intro i hi
    exact infinityZOverT_coeff_of_lt H n i hi

/-- The order-one chart coordinate admits a two-sided compositional inverse.
This is the formal inverse lemma used implicitly in the paper's smoothness
argument at infinity. -/
theorem infinitySOverT_has_compositional_inverse [CharZero 𝕜]
    (H : 𝕜[X]) (n : ℕ) (hn : 2 ≤ n) (hdeg : H.natDegree = n) :
    ∃ Q : PowerSeries 𝕜,
      Q.constantCoeff = 0 ∧
      PowerSeries.HasSubst Q ∧
      (infinitySOverT H n).subst Q = PowerSeries.X ∧
      Q.subst (infinitySOverT H n) = PowerSeries.X := by
  let P := infinitySOverT H n
  have hP : P.constantCoeff = 0 :=
    infinitySOverT_constantCoeff H n hdeg
  have hunit : IsUnit (P.coeff 1) :=
    isUnit_iff_ne_zero.mpr (infinitySOverT_coeff_one_ne_zero H n hn hdeg)
  let Q := P.substInvOfIsUnit hunit
  refine ⟨Q, ?_, ?_, ?_, ?_⟩
  · exact PowerSeries.constantCoeff_substInvOfIsUnit P hunit
  · exact PowerSeries.HasSubst.substInvOfIsUnit P hunit
  · exact PowerSeries.subst_substInvOfIsUnit_right P hP hunit
  · exact PowerSeries.subst_substInvOfIsUnit_left P hP hunit

/-- In the completed local chart, `Z/T` is a formal graph over the
uniformizing coordinate `S/T`. -/
theorem infinity_is_formal_graph [CharZero 𝕜]
    (H : 𝕜[X]) (n : ℕ) (hn : 2 ≤ n) (hdeg : H.natDegree = n) :
    ∃ G : PowerSeries 𝕜,
      G.constantCoeff = 0 ∧
      G.subst (infinitySOverT H n) = infinityZOverT H n := by
  obtain ⟨Q, hQ0, hQ, _, hQP⟩ :=
    infinitySOverT_has_compositional_inverse H n hn hdeg
  let G : PowerSeries 𝕜 := (infinityZOverT H n).subst Q
  refine ⟨G, ?_, ?_⟩
  · change ((infinityZOverT H n).subst Q).constantCoeff = 0
    apply PowerSeries.constantCoeff_subst_eq_zero hQ0
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]
    exact infinityZOverT_coeff_of_lt H n 0 (by omega)
  · change PowerSeries.subst (R := 𝕜) (S := 𝕜) (τ := Unit)
      (infinitySOverT H n)
      (PowerSeries.subst (R := 𝕜) (S := 𝕜) (τ := Unit)
        Q (infinityZOverT H n)) = infinityZOverT H n
    rw [PowerSeries.subst_comp_subst_apply hQ
      (PowerSeries.HasSubst.of_constantCoeff_zero
        (infinitySOverT_constantCoeff H n hdeg))]
    rw [hQP]
    exact PowerSeries.X_subst _

end DiscriminantPencils
