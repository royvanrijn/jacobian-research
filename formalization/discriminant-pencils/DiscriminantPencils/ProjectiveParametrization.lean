/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import DiscriminantPencils.TangentPencil
import Mathlib.Algebra.Polynomial.Homogenize

/-!
# The homogenized tangent-map coordinates

This file formalizes the homogeneous coordinate forms in equation (2.1) and
their behavior on the affine chart and at infinity.
-/

open Polynomial

namespace DiscriminantPencils

variable {𝕜 : Type*} [Field 𝕜]

/-- The three degree-`n` homogeneous forms
`[S_H : T_H : Z]` from equation (2.1). -/
noncomputable def projectiveCoordinates (H : 𝕜[X]) (n : ℕ) :
    Fin 3 → MvPolynomial (Fin 2) 𝕜
  | 0 => H.derivative.homogenize n
  | 1 => (X * H.derivative - H).homogenize n
  | 2 => MvPolynomial.X 1 ^ n

@[simp] theorem projectiveCoordinates_zero (H : 𝕜[X]) (n : ℕ) :
    projectiveCoordinates H n 0 = H.derivative.homogenize n := rfl

@[simp] theorem projectiveCoordinates_one (H : 𝕜[X]) (n : ℕ) :
    projectiveCoordinates H n 1 =
      (X * H.derivative - H).homogenize n := rfl

@[simp] theorem projectiveCoordinates_two (H : 𝕜[X]) (n : ℕ) :
    projectiveCoordinates H n 2 = MvPolynomial.X 1 ^ n := rfl

/-- Evaluation of a degree-`n` homogenization at `[R : 0]` keeps only its
degree-`n` coefficient. -/
theorem eval_homogenize_at_infinity (p : 𝕜[X]) (n : ℕ) (R : 𝕜) :
    MvPolynomial.eval ![R, 0] (p.homogenize n) = p.coeff n * R ^ n := by
  induction p using Polynomial.induction_on' with
  | add p q ihp ihq =>
      simp [ihp, ihq, add_mul]
  | monomial m a =>
      rcases le_or_gt m n with hmn | hnm
      · rw [homogenize_monomial hmn]
        by_cases h : m = n
        · subst m
          simp [MvPolynomial.eval_monomial]
        · have hlt : m < n := lt_of_le_of_ne hmn h
          have hsub : n - m ≠ 0 := by omega
          simp [MvPolynomial.eval_monomial, coeff_monomial, h,
            hsub]
      · rw [homogenize_monomial_of_lt hnm]
        have hne : m ≠ n := by omega
        simp [coeff_monomial, hne]

/-- The coefficient of degree `n` in `X H' - H`. -/
theorem coeff_X_mul_derivative_sub (H : 𝕜[X]) (n : ℕ) (hn : 1 ≤ n) :
    (X * H.derivative - H).coeff n = (n - 1 : ℕ) * H.coeff n := by
  obtain ⟨m, rfl⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : n ≠ 0)
  simp [coeff_derivative, Nat.cast_succ]
  ring

/-- On the affine chart `U = 1`, the homogeneous coordinates recover the
affine tangent map and `Z = 1`. -/
theorem projectiveCoordinates_affine (H : 𝕜[X]) (n : ℕ) (hn : 1 ≤ n)
    (hdeg : H.natDegree ≤ n) (r : 𝕜) :
    (MvPolynomial.eval ![r, 1] (projectiveCoordinates H n 0),
      MvPolynomial.eval ![r, 1] (projectiveCoordinates H n 1),
      MvPolynomial.eval ![r, 1] (projectiveCoordinates H n 2)) =
      (H.derivative.eval r, r * H.derivative.eval r - H.eval r, 1) := by
  have hd : H.derivative.natDegree ≤ n :=
    (natDegree_derivative_le H).trans ((Nat.sub_le _ _).trans hdeg)
  have ht : (X * H.derivative - H).natDegree ≤ n := by
    apply (natDegree_sub_le _ _).trans
    apply max_le
    · calc
        (X * H.derivative).natDegree ≤ X.natDegree + H.derivative.natDegree :=
          natDegree_mul_le
        _ ≤ 1 + (H.natDegree - 1) := by
          simpa using Nat.add_le_add_left (natDegree_derivative_le H) 1
        _ ≤ n := by omega
    · exact hdeg
  have hS :
      MvPolynomial.eval ![r, 1] (projectiveCoordinates H n 0) =
        H.derivative.eval r := by
    simpa [projectiveCoordinates] using
      (Polynomial.eval_homogenize hd ![r, 1] (by simp))
  have hT :
      MvPolynomial.eval ![r, 1] (projectiveCoordinates H n 1) =
        (X * H.derivative - H).eval r := by
    simpa [projectiveCoordinates] using
      (Polynomial.eval_homogenize ht ![r, 1] (by simp))
  rw [hS, hT]
  simp

/-- The `Z` coordinate vanishes exactly when `U = 0` (for positive degree). -/
theorem projective_z_eq_zero_iff {n : ℕ} (hn : 0 < n) (R U : 𝕜) :
    MvPolynomial.eval ![R, U] (projectiveCoordinates (0 : 𝕜[X]) n 2) = 0 ↔
      U = 0 := by
  constructor
  · simpa [projectiveCoordinates] using
      (show U ^ n = 0 → U = 0 from eq_zero_of_pow_eq_zero)
  · rintro rfl
    simp [projectiveCoordinates, hn.ne']

/-- For an exact-degree polynomial in characteristic zero, the homogeneous
coordinate forms have no common nontrivial zero. -/
theorem projectiveCoordinates_no_basepoint [CharZero 𝕜] (H : 𝕜[X]) (n : ℕ)
    (hn : 2 ≤ n) (hdeg : H.natDegree = n) {R U : 𝕜}
    (hRU : R ≠ 0 ∨ U ≠ 0) :
    ¬(MvPolynomial.eval ![R, U] (projectiveCoordinates H n 0) = 0 ∧
      MvPolynomial.eval ![R, U] (projectiveCoordinates H n 1) = 0 ∧
      MvPolynomial.eval ![R, U] (projectiveCoordinates H n 2) = 0) := by
  rintro ⟨_, hT, hZ⟩
  have hU : U = 0 := by
    have : U ^ n = 0 := by simpa [projectiveCoordinates] using hZ
    exact eq_zero_of_pow_eq_zero this
  have hR : R ≠ 0 := by simpa [hU] using hRU
  change MvPolynomial.eval ![R, U]
    ((X * H.derivative - H).homogenize n) = 0 at hT
  rw [hU, eval_homogenize_at_infinity] at hT
  rw [coeff_X_mul_derivative_sub H n (by omega)] at hT
  have hH : H ≠ 0 := by
    intro h
    subst H
    simp at hdeg
    omega
  have hlc : H.coeff n ≠ 0 := by
    simpa [leadingCoeff, hdeg] using leadingCoeff_ne_zero.mpr hH
  exact mul_ne_zero (mul_ne_zero (by exact_mod_cast (show n - 1 ≠ 0 by omega)) hlc)
    (pow_ne_zero n hR) hT

/-- The unique source point over the line at infinity is `[1 : 0]`, in the
set-theoretic homogeneous-coordinate sense. -/
theorem unique_source_point_at_infinity {n : ℕ} (hn : 0 < n) {R U : 𝕜}
    (hRU : R ≠ 0 ∨ U ≠ 0)
    (hZ : MvPolynomial.eval ![R, U]
      (projectiveCoordinates (0 : 𝕜[X]) n 2) = 0) :
    U = 0 ∧ R ≠ 0 := by
  have hU := (projective_z_eq_zero_iff hn R U).1 hZ
  exact ⟨hU, by simpa [hU] using hRU⟩

end DiscriminantPencils
