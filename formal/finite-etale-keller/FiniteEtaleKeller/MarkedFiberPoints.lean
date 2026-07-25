/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.RootPoints
import FiniteEtaleKeller.SourceEquivalence

/-!
# Roots and marked fiber points

This module isolates the broad algebraic core of the fiber theorem.  Suppose
`E'` factors as `g₁` times the marked chart polynomial.  Then, in every
commutative `K`-algebra, roots of `E` are equivalent to marked chart points.
No field structure on the test algebra and no localization are used.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K A : Type*} [Field K]
variable [CommRing A] [Algebra K A]

/-- The polynomial whose value is the marked chart factor after imposing
`Q + β(S) = b`. -/
def markedChartPolynomial (pi b : K) (β : K[X]) : K[X] :=
  1 - C b * X + X * β + C pi * X ^ 2

/-- A marked point over a test algebra: a chart point whose marked coordinate
is `b` and whose `S` coordinate is a root of `E`. -/
@[ext]
structure MarkedFiberPoint (E β : K[X]) (pi b : K)
    (A : Type*) [CommRing A] [Algebra K A] where
  chart : GaugeChart A (algebraMap K A pi)
  root_eq : Polynomial.aeval chart.S E = 0
  marked_eq : chart.Q + Polynomial.aeval chart.S β = algebraMap K A b

namespace PolynomialRoot

variable {E β : K[X]} {pi b : K}

/-- Evaluating the derivative factorization at a root identifies the
normalized derivative unit with the chart factor. -/
theorem normalizedDerivativeUnit_eq_chartFactor
    (hE : E.Separable) (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β)
    (s : PolynomialRoot E A) :
    (s.normalizedDerivativeUnit E hE g₁ : A) =
      1 - s.1 * (algebraMap K A b - Polynomial.aeval s.1 β)
        + algebraMap K A pi * s.1 ^ 2 := by
  rw [normalizedDerivativeUnit_val, hderiv]
  simp only [map_mul, Polynomial.aeval_C]
  have hg :
      algebraMap K A (↑g₁⁻¹ : K) * algebraMap K A (g₁ : K) = 1 := by
    rw [← map_mul]
    simp
  rw [← mul_assoc, hg, one_mul]
  simp [markedChartPolynomial]
  ring

/-- A root determines the unique marked chart point. -/
def toMarkedFiberPoint
    (hE : E.Separable) (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β)
    (s : PolynomialRoot E A) : MarkedFiberPoint E β pi b A where
  chart :=
    { S := s.1
      Q := algebraMap K A b - Polynomial.aeval s.1 β
      d := s.normalizedDerivativeUnit E hE g₁
      chart_eq := s.normalizedDerivativeUnit_eq_chartFactor hE g₁ hderiv }
  root_eq := s.2
  marked_eq := by ring

end PolynomialRoot

namespace MarkedFiberPoint

variable {E β : K[X]} {pi b : K}

/-- Forgetting the marked coordinates leaves the underlying root. -/
def toRoot (p : MarkedFiberPoint E β pi b A) : PolynomialRoot E A :=
  ⟨p.chart.S, p.root_eq⟩

@[simp]
theorem toRoot_val (p : MarkedFiberPoint E β pi b A) : p.toRoot.1 = p.chart.S := rfl

end MarkedFiberPoint

/-- Roots and marked fiber points are equivalent over every commutative test
algebra. -/
def rootEquivMarkedFiberPoint
    {E β : K[X]} {pi b : K}
    (hE : E.Separable) (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β) :
    PolynomialRoot E A ≃ MarkedFiberPoint E β pi b A where
  toFun := PolynomialRoot.toMarkedFiberPoint hE g₁ hderiv
  invFun := MarkedFiberPoint.toRoot
  left_inv := by
    intro s
    apply PolynomialRoot.ext
    rfl
  right_inv := by
    intro p
    have hQ :
        algebraMap K A b - Polynomial.aeval p.chart.S β = p.chart.Q := by
      linear_combination p.marked_eq
    apply MarkedFiberPoint.ext
    apply GaugeChart.ext
    · rfl
    · exact hQ
    · apply Units.ext
      change
        (p.toRoot.normalizedDerivativeUnit E hE g₁ : A) =
          (p.chart.d : A)
      rw [PolynomialRoot.normalizedDerivativeUnit_eq_chartFactor hE g₁ hderiv]
      rw [hQ, ← p.chart.chart_eq]

/-- Characteristic-zero squarefreeness gives the paper's public form of the
same equivalence. -/
def squarefreeRootEquivMarkedFiberPoint [CharZero K]
    {E β : K[X]} {pi b : K}
    (hE : Squarefree E) (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β) :
    PolynomialRoot E A ≃ MarkedFiberPoint E β pi b A :=
  rootEquivMarkedFiberPoint
    ((PerfectField.separable_iff_squarefree).2 hE) g₁ hderiv

end FiniteEtaleKeller
