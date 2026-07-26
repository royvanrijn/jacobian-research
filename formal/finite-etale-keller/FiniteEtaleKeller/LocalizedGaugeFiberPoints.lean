/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.LocalizedFiberPoints
import FiniteEtaleKeller.FiberNaturality

/-!
# Localized roots and quadratic-gauge source fibers

The separable fiber theorem reconstructs a source point from every root because
`E'` is then automatically invertible modulo `E`.  This module removes that
hypothesis: a root reconstructs a source point exactly when its derivative is a
unit.  Combining this with the localization universal property gives the
scheme represented by `(K[S]/E)[1/E']`.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K A B : Type*} [Field K]
variable [CommRing A] [Algebra K A]
variable [CommRing B] [Algebra K B]

namespace LocalizedPolynomialRoot

variable {E β : K[X]} {pi b : K}

/-- Divide the derivative unit by the nonzero scalar `g₁`. -/
def normalizedDerivativeUnit
    (g₁ : Kˣ) (s : LocalizedPolynomialRoot E A) : Aˣ :=
  Units.map (algebraMap K A) g₁⁻¹ * s.derivativeUnit

@[simp]
theorem normalizedDerivativeUnit_val
    (g₁ : Kˣ) (s : LocalizedPolynomialRoot E A) :
    (s.normalizedDerivativeUnit g₁ : A) =
      algebraMap K A (↑g₁⁻¹ : K) * Polynomial.aeval s.val E.derivative := by
  simp [normalizedDerivativeUnit, s.derivativeUnit_eq]

/-- The normalized derivative is the unit appearing in the marked source
chart. -/
theorem normalizedDerivativeUnit_eq_chartFactor
    (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β)
    (s : LocalizedPolynomialRoot E A) :
    (s.normalizedDerivativeUnit g₁ : A) =
      1 - s.val * (algebraMap K A b - Polynomial.aeval s.val β) +
        algebraMap K A pi * s.val ^ 2 := by
  rw [normalizedDerivativeUnit_val, hderiv]
  simp only [map_mul, Polynomial.aeval_C]
  have hg :
      algebraMap K A (↑g₁⁻¹ : K) * algebraMap K A (g₁ : K) = 1 := by
    rw [← map_mul]
    simp
  rw [← mul_assoc, hg, one_mul]
  simp [markedChartPolynomial]
  ring

/-- A root with invertible derivative determines the unique marked chart
point. -/
def toMarkedFiberPoint
    (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β)
    (s : LocalizedPolynomialRoot E A) : MarkedFiberPoint E β pi b A where
  chart :=
    { S := s.val
      Q := algebraMap K A b - Polynomial.aeval s.val β
      d := s.normalizedDerivativeUnit g₁
      chart_eq := s.normalizedDerivativeUnit_eq_chartFactor g₁ hderiv }
  root_eq := s.root_eq
  marked_eq := by ring

end LocalizedPolynomialRoot

namespace MarkedFiberPoint

variable {E β : K[X]} {pi b : K}

/-- A marked source point carries a canonical derivative unit. -/
def toLocalizedPolynomialRoot
    (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β)
    (p : MarkedFiberPoint E β pi b A) : LocalizedPolynomialRoot E A where
  val := p.chart.S
  root_eq := p.root_eq
  derivativeUnit := Units.map (algebraMap K A) g₁ * p.chart.d
  derivativeUnit_eq := by
    simp only [Units.coe_mul, Units.coe_map]
    rw [hderiv]
    simp only [map_mul, Polynomial.aeval_C]
    congr 1
    simp only [markedChartPolynomial, Polynomial.aeval_sub,
      Polynomial.aeval_one, Polynomial.aeval_mul, Polynomial.aeval_C,
      Polynomial.aeval_X, Polynomial.aeval_add, Polynomial.aeval_pow]
    rw [p.chart.chart_eq]
    linear_combination -p.chart.S * p.marked_eq

@[simp]
theorem toLocalizedPolynomialRoot_val
    (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β)
    (p : MarkedFiberPoint E β pi b A) :
    (p.toLocalizedPolynomialRoot g₁ hderiv).val = p.chart.S := rfl

end MarkedFiberPoint

/-- Roots with invertible derivative and marked chart points are equivalent
over every commutative test algebra. -/
def localizedRootEquivMarkedFiberPoint
    {E β : K[X]} {pi b : K}
    (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β) :
    LocalizedPolynomialRoot E A ≃ MarkedFiberPoint E β pi b A where
  toFun := LocalizedPolynomialRoot.toMarkedFiberPoint g₁ hderiv
  invFun := MarkedFiberPoint.toLocalizedPolynomialRoot g₁ hderiv
  left_inv := by
    intro s
    apply LocalizedPolynomialRoot.ext
    rfl
  right_inv := by
    intro p
    have hQ :
        algebraMap K A b - Polynomial.aeval p.chart.S β = p.chart.Q := by
      calc
        algebraMap K A b - Polynomial.aeval p.chart.S β =
            (p.chart.Q + Polynomial.aeval p.chart.S β) -
              Polynomial.aeval p.chart.S β := by rw [p.marked_eq]
        _ = p.chart.Q := by ring
    apply MarkedFiberPoint.ext
    apply GaugeChart.ext
    · rfl
    · exact hQ
    · apply Units.ext
      simp [LocalizedPolynomialRoot.toMarkedFiberPoint,
        MarkedFiberPoint.toLocalizedPolynomialRoot,
        LocalizedPolynomialRoot.normalizedDerivativeUnit]

/-- Localized roots reconstruct full source-fiber points without a separability
hypothesis. -/
def localizedRootEquivGaugeFiberPoint
    {E β : K[X]} {pi b : K} (a : K)
    (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β) :
    LocalizedPolynomialRoot E A ≃ GaugeFiberPoint E β pi b a A :=
  (localizedRootEquivMarkedFiberPoint (A := A) g₁ hderiv).trans
    (markedFiberPointEquivGaugeFiberPoint a)

/-- The localized-root/source equivalence is natural in the commutative test
algebra. -/
theorem localizedRootEquivGaugeFiberPoint_natural
    {E β : K[X]} {pi b : K} (a : K)
    (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β)
    (f : A →ₐ[K] B) (s : LocalizedPolynomialRoot E A) :
    GaugeFiberPoint.map f
        (localizedRootEquivGaugeFiberPoint (A := A) a g₁ hderiv s) =
      localizedRootEquivGaugeFiberPoint (A := B) a g₁ hderiv (s.map f) := by
  apply (localizedRootEquivGaugeFiberPoint (A := B) a g₁ hderiv).symm.injective
  apply LocalizedPolynomialRoot.ext
  rfl

#print axioms localizedRootEquivGaugeFiberPoint
#print axioms localizedRootEquivGaugeFiberPoint_natural

end FiniteEtaleKeller
