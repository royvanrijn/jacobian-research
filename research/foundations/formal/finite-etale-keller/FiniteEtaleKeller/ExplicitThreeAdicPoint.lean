/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.NumberTheory.Padics.Hensel
import FiniteEtaleKeller.ExplicitFiber

/-!
# A three-adic point on the explicit quintic

The cubic factor `X³ - 19` has a root in `ℤ_[3]`.  Indeed, `-2` is a
sufficiently accurate approximate root for the strong form of Hensel's lemma:
the polynomial value is `-27` and the derivative value is `12`.
-/

noncomputable section

namespace FiniteEtaleKeller.ExplicitQuintic

/-- The integral cubic factor of the explicit quintic. -/
def cubic19Int : Polynomial ℤ :=
  Polynomial.X ^ 3 - 19

private theorem cubic19Int_at_neg_two :
    Polynomial.aeval (-2 : ℤ_[3]) cubic19Int = -27 := by
  norm_num [cubic19Int, Polynomial.aeval_def]

private theorem cubic19Int_derivative_at_neg_two :
    Polynomial.aeval (-2 : ℤ_[3]) cubic19Int.derivative = 12 := by
  norm_num [cubic19Int, Polynomial.aeval_def]

private theorem norm_neg_twenty_seven :
    ‖(-27 : ℤ_[3])‖ = (1 / 27 : ℝ) := by
  have h : (-27 : ℤ_[3]) = -(3 : ℤ_[3]) ^ 3 := by norm_num
  have h3 : ‖(3 : ℤ_[3])‖ = (3 : ℝ)⁻¹ := by
    change ‖((3 : ℕ) : ℤ_[3])‖ = (3 : ℝ)⁻¹
    exact PadicInt.norm_p
  calc
    ‖(-27 : ℤ_[3])‖ = ‖(3 : ℤ_[3])‖ ^ 3 := by rw [h, norm_neg, norm_pow]
    _ = ((3 : ℝ)⁻¹) ^ 3 := congrArg (· ^ 3) h3
    _ = 1 / 27 := by norm_num

private theorem norm_twelve :
    ‖(12 : ℤ_[3])‖ = (1 / 3 : ℝ) := by
  have h : (12 : ℤ_[3]) = (3 : ℤ_[3]) * 4 := by norm_num
  have h3 : ‖(3 : ℤ_[3])‖ = (3 : ℝ)⁻¹ := by
    change ‖((3 : ℕ) : ℤ_[3])‖ = (3 : ℝ)⁻¹
    exact PadicInt.norm_p
  have h4 : ‖(4 : ℤ_[3])‖ = 1 := by
    change ‖((4 : ℕ) : ℤ_[3])‖ = 1
    exact PadicInt.norm_natCast_eq_one_iff.mpr (by norm_num)
  calc
    ‖(12 : ℤ_[3])‖ = ‖(3 : ℤ_[3])‖ * ‖(4 : ℤ_[3])‖ := by rw [h, norm_mul]
    _ = (3 : ℝ)⁻¹ * ‖(4 : ℤ_[3])‖ := congrArg (· * ‖(4 : ℤ_[3])‖) h3
    _ = (3 : ℝ)⁻¹ * 1 := congrArg ((3 : ℝ)⁻¹ * ·) h4
    _ = 1 / 3 := by norm_num

/-- The cubic factor `X³ - 19` has a root in the three-adic integers. -/
theorem cubic19_has_threeAdicInt_root :
    ∃ z : ℤ_[3], Polynomial.aeval z cubic19Int = 0 := by
  have hnorm :
      ‖Polynomial.aeval (-2 : ℤ_[3]) cubic19Int‖ <
        ‖Polynomial.aeval (-2 : ℤ_[3]) cubic19Int.derivative‖ ^ 2 := by
    rw [cubic19Int_at_neg_two, cubic19Int_derivative_at_neg_two,
      norm_neg_twenty_seven, norm_twelve]
    norm_num
  obtain ⟨z, hz, -⟩ := hensels_lemma hnorm
  exact ⟨z, hz⟩

/-- The displayed quintic has a root in the three-adic field. -/
theorem p5_has_threeAdic_root :
    ∃ z : ℚ_[3], Polynomial.aeval z p5 = 0 := by
  obtain ⟨z, hz⟩ := cubic19_has_threeAdicInt_root
  refine ⟨(z : ℚ_[3]), ?_⟩
  have hzExpr : z ^ 3 - 19 = 0 := by
    simpa [cubic19Int, Polynomial.aeval_def] using hz
  have hzExpr' : (z : ℚ_[3]) ^ 3 - 19 = 0 := by
    have hmap := congrArg (algebraMap ℤ_[3] ℚ_[3]) hzExpr
    simp only [map_sub, map_pow, map_zero, PadicInt.algebraMap_apply] at hmap
    have h19 : ((19 : ℤ_[3]) : ℚ_[3]) = (19 : ℚ_[3]) :=
      PadicInt.coe_natCast 19
    rw [h19] at hmap
    exact hmap
  simp [p5, Polynomial.aeval_def, hzExpr']

/-- The literal displayed fiber over `(1, 0, -38)` has a three-adic point. -/
theorem integralFiberPoint_threeAdic_nonempty :
    Nonempty (IntegralFiberPoint ℚ_[3]) := by
  obtain ⟨z, hz⟩ := p5_has_threeAdic_root
  let root : PolynomialRoot p5 ℚ_[3] := ⟨z, hz⟩
  exact ⟨integralFiberRepresentingEquiv
    ((PolynomialRoot.algHomEquiv p5 ℚ_[3]).symm root)⟩

#print axioms cubic19_has_threeAdicInt_root
#print axioms p5_has_threeAdic_root
#print axioms integralFiberPoint_threeAdic_nonempty

end FiniteEtaleKeller.ExplicitQuintic
