/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib

/-!
# Stable-separation certificates for the gauge families

This file formalizes the exact lattice arithmetic used after stable
normalization has identified the intrinsic Fitting divisor or the complete
boundary.  It proves translation and unimodular invariance of normalized
quadrilateral area, evaluates the power-shift polygon, and proves that both
the area and cubic boundary-count certificates recover the family parameter.

The geometric statement that stable polynomial left--right equivalence
transports the selected normalized Fitting divisor and the complete boundary
is deliberately kept as an input to the final implication theorems.
-/

namespace FiniteEtaleKeller

/-- The determinant pairing on the rank-two lattice. -/
def latticeCross (u v : ℤ × ℤ) : ℤ :=
  u.1 * v.2 - u.2 * v.1

/-- Twice the oriented area of an ordered lattice quadrilateral. -/
def latticeQuadrilateralDoubleArea
    (p₀ p₁ p₂ p₃ : ℤ × ℤ) : ℤ :=
  latticeCross p₀ p₁ + latticeCross p₁ p₂ +
    latticeCross p₂ p₃ + latticeCross p₃ p₀

/-- The normalized (unoriented) lattice area. -/
def normalizedLatticeQuadrilateralArea
    (p₀ p₁ p₂ p₃ : ℤ × ℤ) : ℕ :=
  (latticeQuadrilateralDoubleArea p₀ p₁ p₂ p₃).natAbs

/-- Translation of a lattice point. -/
def translateLatticePoint (u p : ℤ × ℤ) : ℤ × ℤ :=
  (p.1 + u.1, p.2 + u.2)

/-- An integral two-by-two linear transformation. -/
def transformLatticePoint
    (a b c d : ℤ) (p : ℤ × ℤ) : ℤ × ℤ :=
  (a * p.1 + b * p.2, c * p.1 + d * p.2)

/-- Translation, corresponding to multiplication of a Laurent polynomial by
a unit, preserves normalized Newton area. -/
theorem normalizedLatticeQuadrilateralArea_translate
    (u p₀ p₁ p₂ p₃ : ℤ × ℤ) :
    normalizedLatticeQuadrilateralArea
        (translateLatticePoint u p₀) (translateLatticePoint u p₁)
        (translateLatticePoint u p₂) (translateLatticePoint u p₃) =
      normalizedLatticeQuadrilateralArea p₀ p₁ p₂ p₃ := by
  simp only [normalizedLatticeQuadrilateralArea,
    latticeQuadrilateralDoubleArea, latticeCross,
    translateLatticePoint]
  congr 1
  ring

/-- A lattice linear transformation scales the determinant pairing by its
determinant. -/
theorem latticeCross_transform
    (a b c d : ℤ) (p q : ℤ × ℤ) :
    latticeCross (transformLatticePoint a b c d p)
        (transformLatticePoint a b c d q) =
      (a * d - b * c) * latticeCross p q := by
  simp [latticeCross, transformLatticePoint]
  ring

/-- A lattice linear transformation scales oriented quadrilateral area by
its determinant. -/
theorem latticeQuadrilateralDoubleArea_transform
    (a b c d : ℤ) (p₀ p₁ p₂ p₃ : ℤ × ℤ) :
    latticeQuadrilateralDoubleArea
        (transformLatticePoint a b c d p₀)
        (transformLatticePoint a b c d p₁)
        (transformLatticePoint a b c d p₂)
        (transformLatticePoint a b c d p₃) =
      (a * d - b * c) *
        latticeQuadrilateralDoubleArea p₀ p₁ p₂ p₃ := by
  simp only [latticeQuadrilateralDoubleArea, latticeCross_transform]
  ring

/-- Every `GL₂(ℤ)` coordinate change preserves normalized Newton area. -/
theorem normalizedLatticeQuadrilateralArea_unimodular
    (a b c d : ℤ) (p₀ p₁ p₂ p₃ : ℤ × ℤ)
    (hdet : a * d - b * c = 1 ∨ a * d - b * c = -1) :
    normalizedLatticeQuadrilateralArea
        (transformLatticePoint a b c d p₀)
        (transformLatticePoint a b c d p₁)
        (transformLatticePoint a b c d p₂)
        (transformLatticePoint a b c d p₃) =
      normalizedLatticeQuadrilateralArea p₀ p₁ p₂ p₃ := by
  unfold normalizedLatticeQuadrilateralArea
  rw [latticeQuadrilateralDoubleArea_transform]
  rcases hdet with hdet | hdet
  · rw [hdet, one_mul]
  · rw [hdet, neg_one_mul, Int.natAbs_neg]

/-- The four ordered vertices of the power-shift Fitting Newton polygon.
For `N = 4`, the middle two vertices coincide. -/
def powerShiftFittingNewtonVertices (N m : ℕ) : Fin 4 → ℤ × ℤ :=
  ![
    (0, 0),
    ((4 + m : ℕ), (3 : ℕ)),
    ((N + m : ℕ), (N - 1 : ℕ)),
    ((1 : ℕ), (2 : ℕ))
  ]

/-- The exact shoelace evaluation before removing orientation. -/
theorem powerShiftFittingNewtonDoubleArea
    (N m : ℕ) (hN : 4 ≤ N) :
    latticeQuadrilateralDoubleArea
        (powerShiftFittingNewtonVertices N m 0)
        (powerShiftFittingNewtonVertices N m 1)
        (powerShiftFittingNewtonVertices N m 2)
        (powerShiftFittingNewtonVertices N m 3) =
      ((2 * N - 3 + (N - 2) * m : ℕ) : ℤ) := by
  have hN1 : 1 ≤ N := by omega
  have hN2 : 2 ≤ N := by omega
  have hN3 : 3 ≤ 2 * N := by omega
  simp only [powerShiftFittingNewtonVertices, Matrix.cons_val_zero,
    Matrix.cons_val_one, Matrix.cons_val,
    latticeQuadrilateralDoubleArea, latticeCross]
  push_cast [Nat.cast_sub hN1, Nat.cast_sub hN2, Nat.cast_sub hN3]
  ring

/-- The normalized Fitting Newton area certificate. -/
def powerShiftFittingNewtonArea (N m : ℕ) : ℕ :=
  2 * N - 3 + (N - 2) * m

/-- The ordered polygon has the advertised normalized area. -/
theorem normalized_powerShiftFittingNewtonArea
    (N m : ℕ) (hN : 4 ≤ N) :
    normalizedLatticeQuadrilateralArea
        (powerShiftFittingNewtonVertices N m 0)
        (powerShiftFittingNewtonVertices N m 1)
        (powerShiftFittingNewtonVertices N m 2)
        (powerShiftFittingNewtonVertices N m 3) =
      powerShiftFittingNewtonArea N m := by
  unfold normalizedLatticeQuadrilateralArea
  rw [powerShiftFittingNewtonDoubleArea N m hN]
  exact Int.natAbs_natCast _

/-- For fixed degree at least four, the Fitting Newton area strictly
increases with the common shift. -/
theorem powerShiftFittingNewtonArea_strictMono
    (N : ℕ) (hN : 4 ≤ N) :
    StrictMono (powerShiftFittingNewtonArea N) := by
  intro m m' hmm'
  rw [powerShiftFittingNewtonArea, powerShiftFittingNewtonArea]
  exact Nat.add_lt_add_left
    ((Nat.mul_lt_mul_left (by omega : 0 < N - 2)).2 hmm') _

/-- Equality of the intrinsic Fitting-area certificates forces equality of
the power-shift parameters. -/
theorem powerShift_eq_of_equalFittingNewtonArea
    (N m m' : ℕ) (hN : 4 ≤ N)
    (harea :
      powerShiftFittingNewtonArea N m =
        powerShiftFittingNewtonArea N m') :
    m = m' :=
  (powerShiftFittingNewtonArea_strictMono N hN).injective harea

/-- Abstract stable-separation bridge for the power-shift family: the only
geometric input needed here is preservation of the certified area. -/
theorem powerShift_eq_of_stable_preservesFittingNewtonArea
    {Map : Type*} (family : ℕ → Map)
    (StableEquivalent : Map → Map → Prop)
    (N m m' : ℕ) (hN : 4 ≤ N)
    (hpreserves :
      StableEquivalent (family m) (family m') →
        powerShiftFittingNewtonArea N m =
          powerShiftFittingNewtonArea N m')
    (hstable : StableEquivalent (family m) (family m')) :
    m = m' :=
  powerShift_eq_of_equalFittingNewtonArea N m m' hN (hpreserves hstable)

/-- The complete cubic boundary count: one repeated-root discriminant
component and `n-1` degree-drop components. -/
def cubicLiftBoundaryTargetComponentCount (n : ℕ) : ℕ :=
  1 + (n - 1)

/-- In the cubic range, the complete boundary count is exactly the lift
parameter. -/
theorem cubicLiftBoundaryTargetComponentCount_eq
    (n : ℕ) (hn : 4 ≤ n) :
    cubicLiftBoundaryTargetComponentCount n = n := by
  simp [cubicLiftBoundaryTargetComponentCount, Nat.add_sub_of_le (by omega : 1 ≤ n)]

/-- Equality of complete cubic boundary counts forces equality of the lift
parameters. -/
theorem cubicLift_eq_of_equalBoundaryTargetComponentCount
    (n n' : ℕ) (hn : 4 ≤ n) (hn' : 4 ≤ n')
    (hcount :
      cubicLiftBoundaryTargetComponentCount n =
        cubicLiftBoundaryTargetComponentCount n') :
    n = n' := by
  rw [cubicLiftBoundaryTargetComponentCount_eq n hn,
    cubicLiftBoundaryTargetComponentCount_eq n' hn'] at hcount
  exact hcount

/-- Abstract stable-separation bridge for cubic lifts: the geometric input is
preservation of the complete canonical boundary-component count. -/
theorem cubicLift_eq_of_stable_preservesBoundaryTargetComponentCount
    {Map : Type*} (family : ℕ → Map)
    (StableEquivalent : Map → Map → Prop)
    (n n' : ℕ) (hn : 4 ≤ n) (hn' : 4 ≤ n')
    (hpreserves :
      StableEquivalent (family n) (family n') →
        cubicLiftBoundaryTargetComponentCount n =
          cubicLiftBoundaryTargetComponentCount n')
    (hstable : StableEquivalent (family n) (family n')) :
    n = n' :=
  cubicLift_eq_of_equalBoundaryTargetComponentCount
    n n' hn hn' (hpreserves hstable)

#print axioms normalizedLatticeQuadrilateralArea_translate
#print axioms normalizedLatticeQuadrilateralArea_unimodular
#print axioms powerShiftFittingNewtonDoubleArea
#print axioms normalized_powerShiftFittingNewtonArea
#print axioms powerShift_eq_of_equalFittingNewtonArea
#print axioms powerShift_eq_of_stable_preservesFittingNewtonArea
#print axioms cubicLiftBoundaryTargetComponentCount_eq
#print axioms cubicLift_eq_of_stable_preservesBoundaryTargetComponentCount

end FiniteEtaleKeller
