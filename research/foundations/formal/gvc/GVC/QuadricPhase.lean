/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GVC.ReynoldsExpansion
import Mathlib.Algebra.Polynomial.Laurent

/-!
# Laurent restriction to the split quadric

The ring map in this file substitutes

`x ↦ z`, `y ↦ (1-t²)z⁻¹`, `t ↦ t`.

It therefore imposes `t²+xy=1` while retaining the phase exponent as a
Laurent degree.  This is the multiplicative algebraic version of the
quadric parameterization used in the manuscript.
-/

namespace GVC

open MvPolynomial
open scoped LaurentPolynomial

abbrev PhaseLaurent := LaurentPolynomial (Polynomial ℚ)

noncomputable def quadricPhaseRestriction :
    TernaryPolynomial →+* PhaseLaurent :=
  MvPolynomial.eval₂Hom
    (LaurentPolynomial.C.comp Polynomial.C)
    ![LaurentPolynomial.T 1,
      LaurentPolynomial.C (1 - Polynomial.X ^ 2) * LaurentPolynomial.T (-1),
      LaurentPolynomial.C Polynomial.X]

@[simp] theorem quadricPhaseRestriction_gvcX :
    quadricPhaseRestriction gvcX = LaurentPolynomial.T 1 := by
  simp [quadricPhaseRestriction, gvcX]

@[simp] theorem quadricPhaseRestriction_gvcY :
    quadricPhaseRestriction gvcY =
      LaurentPolynomial.C (1 - Polynomial.X ^ 2) *
        LaurentPolynomial.T (-1) := by
  simp [quadricPhaseRestriction, gvcY, Matrix.cons_val_one]

@[simp] theorem quadricPhaseRestriction_gvcT :
    quadricPhaseRestriction gvcT =
      LaurentPolynomial.C Polynomial.X := by
  simp [quadricPhaseRestriction, gvcT, Matrix.cons_val_two]

theorem quadricPhaseRestriction_gvcRho :
    quadricPhaseRestriction gvcRho = 1 := by
  rw [gvcRho, map_add, map_pow, map_mul,
    quadricPhaseRestriction_gvcT, quadricPhaseRestriction_gvcX,
    quadricPhaseRestriction_gvcY]
  rw [LaurentPolynomial.T_mul, LaurentPolynomial.mul_T_assoc]
  norm_num

theorem quadricPhaseRestriction_gvcA :
    quadricPhaseRestriction gvcA =
      1 + LaurentPolynomial.T 2 := by
  rw [gvcA, map_add, quadricPhaseRestriction_gvcRho,
    map_pow, quadricPhaseRestriction_gvcX,
    LaurentPolynomial.T_pow]
  norm_num

theorem quadricPhaseRestriction_gvcC :
    quadricPhaseRestriction gvcC =
      LaurentPolynomial.T (-1) *
        (1 - LaurentPolynomial.C (Polynomial.X ^ 2) *
          (1 + LaurentPolynomial.T 2) ^ 2) := by
  have hcusp₀ := congrArg quadricPhaseRestriction gvc_cusp_identity
  simp only [map_mul, map_sub, map_pow] at hcusp₀
  rw [quadricPhaseRestriction_gvcX, quadricPhaseRestriction_gvcRho,
    quadricPhaseRestriction_gvcT, quadricPhaseRestriction_gvcA] at hcusp₀
  have hcusp :
      LaurentPolynomial.T 1 * quadricPhaseRestriction gvcC =
        1 - LaurentPolynomial.C Polynomial.X ^ 2 *
          (1 + LaurentPolynomial.T 2) ^ 2 := by
    simpa using hcusp₀
  calc
    quadricPhaseRestriction gvcC =
        LaurentPolynomial.T (-1) *
          (LaurentPolynomial.T 1 * quadricPhaseRestriction gvcC) := by
      rw [← mul_assoc, ← LaurentPolynomial.T_add]
      norm_num
    _ = LaurentPolynomial.T (-1) *
          (1 - LaurentPolynomial.C Polynomial.X ^ 2 *
            (1 + LaurentPolynomial.T 2) ^ 2) := by rw [hcusp]
    _ = LaurentPolynomial.T (-1) *
          (1 - LaurentPolynomial.C (Polynomial.X ^ 2) *
            (1 + LaurentPolynomial.T 2) ^ 2) := by rw [map_pow]

/-- The paper's Laurent formula for the concrete polynomial on `rho=1`. -/
theorem quadricPhaseRestriction_gvcP :
    quadricPhaseRestriction gvcP =
      LaurentPolynomial.T (-2) * (1 + LaurentPolynomial.T 2) *
        (1 - LaurentPolynomial.C (Polynomial.X ^ 2) *
          (1 + LaurentPolynomial.T 2) ^ 2) ^ 2 := by
  rw [gvcP, map_mul, map_pow, quadricPhaseRestriction_gvcA,
    quadricPhaseRestriction_gvcC, mul_pow,
    LaurentPolynomial.T_pow]
  norm_num
  ac_rfl

set_option linter.flexible false in
/-- Restriction of one ternary monomial: the Laurent exponent is the phase
difference `α₀-α₁`, and its coefficient is the height polynomial obtained
from `xy=1-t²`. -/
theorem quadricPhaseRestriction_monomial
    (α : Fin 3 →₀ ℕ) (a : ℚ) :
    quadricPhaseRestriction (MvPolynomial.monomial α a) =
      LaurentPolynomial.C
          (Polynomial.C a * (1 - Polynomial.X ^ 2) ^ (α 1) *
            Polynomial.X ^ (α 2)) *
        LaurentPolynomial.T ((α 0 : ℤ) - (α 1 : ℤ)) := by
  rw [quadricPhaseRestriction, MvPolynomial.eval₂Hom_monomial,
    Finsupp.prod_fintype]
  · simp [Fin.prod_univ_succ, mul_pow, LaurentPolynomial.T_pow]
    have hT :
        (LaurentPolynomial.T (- (α 1 : ℤ)) : PhaseLaurent) *
            LaurentPolynomial.T (α 0 : ℤ) =
          LaurentPolynomial.T ((α 0 : ℤ) - (α 1 : ℤ)) := by
      rw [← LaurentPolynomial.T_add]
      congr 1
      ring
    rw [← hT]
    ring
  · intro i
    simp

/-- The Laurent constant term of the split-quadric restriction. -/
noncomputable def quadricPhaseConstant
    (p : TernaryPolynomial) : Polynomial ℚ :=
  (quadricPhaseRestriction p).coeff 0

@[simp] theorem quadricPhaseConstant_zero :
    quadricPhaseConstant (0 : TernaryPolynomial) = 0 := by
  simp [quadricPhaseConstant]

@[simp] theorem quadricPhaseConstant_add
    (p q : TernaryPolynomial) :
    quadricPhaseConstant (p + q) =
      quadricPhaseConstant p + quadricPhaseConstant q := by
  simp [quadricPhaseConstant]

theorem quadricPhaseConstant_monomial
    (α : Fin 3 →₀ ℕ) (a : ℚ) :
    quadricPhaseConstant (MvPolynomial.monomial α a) =
      if α 0 = α 1 then
        Polynomial.C a * (1 - Polynomial.X ^ 2) ^ (α 1) *
          Polynomial.X ^ (α 2)
      else 0 := by
  rw [quadricPhaseConstant, quadricPhaseRestriction_monomial]
  change
    (LaurentPolynomial.C
        (Polynomial.C a * (1 - Polynomial.X ^ 2) ^ (α 1) *
          Polynomial.X ^ (α 2)) *
      LaurentPolynomial.T ((α 0 : ℤ) - (α 1 : ℤ))).coeff 0 = _
  rw [← LaurentPolynomial.single_eq_C_mul_T]
  change
    (Finsupp.single ((α 0 : ℤ) - (α 1 : ℤ))
      (Polynomial.C a * (1 - Polynomial.X ^ 2) ^ (α 1) *
        Polynomial.X ^ (α 2)) : ℤ →₀ Polynomial ℚ) 0 = _
  simp only [Finsupp.single_apply]
  have hz :
      ((α 0 : ℤ) - (α 1 : ℤ) = 0) ↔ α 0 = α 1 := by
    omega
  simp only [hz]

@[simp] theorem reynoldsPhasePolynomial_zero (k : ℕ) :
    reynoldsPhasePolynomial k (0 : TernaryPolynomial) = 0 := by
  simp [reynoldsPhasePolynomial]

@[simp] theorem reynoldsPhasePolynomial_add
    (k : ℕ) (p q : TernaryPolynomial) :
    reynoldsPhasePolynomial k (p + q) =
      reynoldsPhasePolynomial k p + reynoldsPhasePolynomial k q := by
  simp only [reynoldsPhasePolynomial, coeff_add, map_add, add_mul,
    Finset.sum_add_distrib]

private theorem ternary_degree_eq_sum (α : Fin 3 →₀ ℕ) :
    α.degree = α 0 + α 1 + α 2 := by
  rw [Finsupp.degree_eq_sum]
  simp [Fin.sum_univ_succ, add_assoc]

theorem quadricPhaseConstant_monomial_eq_reynoldsPhasePolynomial
    (k : ℕ) (α : Fin 3 →₀ ℕ) (a : ℚ)
    (hdegree : α.degree = 2 * k) :
    quadricPhaseConstant (MvPolynomial.monomial α a) =
      reynoldsPhasePolynomial k (MvPolynomial.monomial α a) := by
  by_cases hbalanced : α 0 = α 1
  · have hsum : α 0 + α 1 + α 2 = 2 * k := by
      rw [← ternary_degree_eq_sum, hdegree]
    have hjle : α 0 ≤ k := by omega
    have htwo : α 2 = 2 * (k - α 0) := by omega
    have halpha : α = reynoldsIndex k (α 0) := by
      ext i
      fin_cases i
      · simp
      · simpa using hbalanced.symm
      · simpa using htwo
    rw [quadricPhaseConstant_monomial, if_pos hbalanced,
      reynoldsPhasePolynomial]
    apply Eq.symm
    calc
      (∑ j ∈ Finset.range (k + 1),
          Polynomial.C
              ((MvPolynomial.monomial α a).coeff (reynoldsIndex k j)) *
            ((1 - Polynomial.X ^ 2) ^ j *
              Polynomial.X ^ (2 * (k - j)))) =
          Polynomial.C
              ((MvPolynomial.monomial α a).coeff
                (reynoldsIndex k (α 0))) *
            ((1 - Polynomial.X ^ 2) ^ (α 0) *
              Polynomial.X ^ (2 * (k - α 0))) := by
        apply Finset.sum_eq_single (α 0)
        · intro j hj hne
          have hindex : α ≠ reynoldsIndex k j := by
            intro h
            have hzero := congrArg (fun β : Fin 3 →₀ ℕ ↦ β 0) h
            simp only [reynoldsIndex_zero] at hzero
            exact hne hzero.symm
          simp [MvPolynomial.coeff_monomial, hindex]
        · intro hnot
          exact (hnot (Finset.mem_range.mpr (Nat.lt_succ_of_le hjle))).elim
      _ = Polynomial.C a * (1 - Polynomial.X ^ 2) ^ α 1 *
            Polynomial.X ^ α 2 := by
        have hcoeff :
            (MvPolynomial.monomial α a).coeff
                (reynoldsIndex k (α 0)) = a := by
          rw [← halpha]
          simp
        rw [hcoeff, ← hbalanced, htwo]
        ring
  · rw [quadricPhaseConstant_monomial, if_neg hbalanced,
      reynoldsPhasePolynomial]
    apply Eq.symm
    apply Finset.sum_eq_zero
    intro j hj
    have hindex : α ≠ reynoldsIndex k j := by
      intro h
      apply hbalanced
      calc
        α 0 = reynoldsIndex k j 0 := congrArg (fun β : Fin 3 →₀ ℕ ↦ β 0) h
        _ = reynoldsIndex k j 1 := by simp
        _ = α 1 := (congrArg (fun β : Fin 3 →₀ ℕ ↦ β 1) h).symm
    simp [MvPolynomial.coeff_monomial, hindex]

/-- For a homogeneous ternary polynomial of degree `2k`, the Laurent
constant term on the split quadric is exactly the diagonal phase-height
polynomial used by the Reynolds expansion. -/
theorem quadricPhaseConstant_eq_reynoldsPhasePolynomial
    (k : ℕ) {p : TernaryPolynomial} (hp : p.IsHomogeneous (2 * k)) :
    quadricPhaseConstant p = reynoldsPhasePolynomial k p := by
  induction hp using IsWeightedHomogeneous.induction_on with
  | zero => simp
  | add p q hp hq ihp ihq =>
      rw [quadricPhaseConstant_add, reynoldsPhasePolynomial_add, ihp, ihq]
  | monomial α a hα =>
      apply quadricPhaseConstant_monomial_eq_reynoldsPhasePolynomial
      calc
        α.degree = (Finsupp.weight (1 : Fin 3 → ℕ)) α :=
          congrArg (fun f ↦ f α) Finsupp.degree_eq_weight_one
        _ = 2 * k := hα

/-- Complete algebraic Reynolds--phase bridge for homogeneous inputs:
normalized apolar contraction equals formal integration of the Laurent
constant term after restriction to `t²+xy=1`. -/
theorem algebraicReynoldsMoment_eq_formalIntegral_quadricPhaseConstant
    (k : ℕ) {p : TernaryPolynomial} (hp : p.IsHomogeneous (2 * k)) :
    algebraicReynoldsMoment k p =
      formalIntegral01 (quadricPhaseConstant p) := by
  rw [algebraicReynoldsMoment_eq_formalIntegral_phase,
    quadricPhaseConstant_eq_reynoldsPhasePolynomial k hp]

@[simp] theorem quadricPhaseRestriction_gvcQ :
    quadricPhaseRestriction gvcQ = LaurentPolynomial.T 2 := by
  rw [gvcQ, map_pow, quadricPhaseRestriction_gvcX,
    LaurentPolynomial.T_pow]
  norm_num

theorem quadricPhaseRestriction_gvcP_pow (m : ℕ) :
    quadricPhaseRestriction (gvcP ^ m) =
      (LaurentPolynomial.T (-2) * (1 + LaurentPolynomial.T 2) *
        (1 - LaurentPolynomial.C (Polynomial.X ^ 2) *
          (1 + LaurentPolynomial.T 2) ^ 2) ^ 2) ^ m := by
  rw [map_pow, quadricPhaseRestriction_gvcP]

theorem quadricPhaseRestriction_gvcQ_mul_gvcP_pow (m : ℕ) :
    quadricPhaseRestriction (gvcQ * gvcP ^ m) =
      LaurentPolynomial.T 2 *
        (LaurentPolynomial.T (-2) * (1 + LaurentPolynomial.T 2) *
          (1 - LaurentPolynomial.C (Polynomial.X ^ 2) *
            (1 + LaurentPolynomial.T 2) ^ 2) ^ 2) ^ m := by
  rw [map_mul, quadricPhaseRestriction_gvcQ,
    quadricPhaseRestriction_gvcP_pow]

/-- The pure moment has now been reduced, without an analytic sphere model,
to formal height integration of the literal Laurent constant term appearing
in the manuscript. -/
theorem algebraicReynoldsMoment_gvcP_pow_eq_laurent (m : ℕ) :
    algebraicReynoldsMoment (6 * m) (gvcP ^ m) =
      formalIntegral01
        ((((LaurentPolynomial.T (-2) : PhaseLaurent) *
          (1 + LaurentPolynomial.T 2) *
          (1 - LaurentPolynomial.C (Polynomial.X ^ 2) *
            (1 + LaurentPolynomial.T 2) ^ 2) ^ 2) ^ m).coeff 0) := by
  rw [algebraicReynoldsMoment_eq_formalIntegral_quadricPhaseConstant
      (6 * m) (gvcP_pow_isHomogeneous m),
    quadricPhaseConstant, quadricPhaseRestriction_gvcP_pow]

/-- The adjacent mixed moment is the corresponding Laurent constant term
after the extra phase shift `x²`. -/
theorem algebraicReynoldsMoment_gvcQ_mul_gvcP_pow_eq_laurent (m : ℕ) :
    algebraicReynoldsMoment (6 * m + 1) (gvcQ * gvcP ^ m) =
      formalIntegral01
        (((LaurentPolynomial.T 2 : PhaseLaurent) *
          (LaurentPolynomial.T (-2) * (1 + LaurentPolynomial.T 2) *
            (1 - LaurentPolynomial.C (Polynomial.X ^ 2) *
              (1 + LaurentPolynomial.T 2) ^ 2) ^ 2) ^ m).coeff 0) := by
  rw [algebraicReynoldsMoment_eq_formalIntegral_quadricPhaseConstant
      (6 * m + 1) (gvcQ_mul_gvcP_pow_isHomogeneous m),
    quadricPhaseConstant, quadricPhaseRestriction_gvcQ_mul_gvcP_pow]

end GVC
