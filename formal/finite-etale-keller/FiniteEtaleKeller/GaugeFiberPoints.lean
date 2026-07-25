/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.MarkedFiberPoints

/-!
# Roots and source fiber points

The marked chart is equivalent to the original source coordinates over every
commutative test ring.  Combining that equivalence with the root description
produces the functor-of-points form of the quadratic-gauge fiber theorem.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K A : Type*} [Field K]
variable [CommRing A] [Algebra K A]

namespace GaugeChart

variable {R : Type*} [CommRing R] {pi : R}

/-- Reconstruction preserves the global marked coordinate `S`. -/
@[simp]
theorem toSource_S (p : GaugeChart R pi) (a : R) :
    (p.toSource a).S = p.S := by
  have h := congrArg GaugeChart.S (p.toSource_toChart a)
  simpa using h

/-- Reconstruction preserves the global marked coordinate `Q`. -/
@[simp]
theorem toSource_Q (p : GaugeChart R pi) (a : R) :
    (p.toSource a).Q = p.Q := by
  have h := congrArg GaugeChart.Q (p.toSource_toChart a)
  simpa using h

end GaugeChart

/-- A full source point on the abstract quadratic-gauge fiber over a test
algebra.  The two source equations are stored by `GaugeSource`; the remaining
two marked equations are the root equation and the fixed `B` coordinate. -/
@[ext]
structure GaugeFiberPoint (E β : K[X]) (pi b a : K)
    (A : Type*) [CommRing A] [Algebra K A] where
  source : GaugeSource A (algebraMap K A pi) (algebraMap K A a)
  root_eq : Polynomial.aeval source.S E = 0
  marked_eq : source.Q + Polynomial.aeval source.S β = algebraMap K A b

namespace MarkedFiberPoint

variable {E β : K[X]} {pi b : K}

/-- Reconstruct a full source fiber point from a marked chart point. -/
def toGaugeFiberPoint
    (p : MarkedFiberPoint E β pi b A) (a : K) :
    GaugeFiberPoint E β pi b a A where
  source := p.chart.toSource (algebraMap K A a)
  root_eq := by
    rw [GaugeChart.toSource_S]
    exact p.root_eq
  marked_eq := by
    rw [GaugeChart.toSource_Q]
    exact p.marked_eq

end MarkedFiberPoint

namespace GaugeFiberPoint

variable {E β : K[X]} {pi b a : K}

/-- Pass from source coordinates to the global marked chart. -/
def toMarkedFiberPoint
    (p : GaugeFiberPoint E β pi b a A) :
    MarkedFiberPoint E β pi b A where
  chart := p.source.toChart
  root_eq := by simpa using p.root_eq
  marked_eq := by simpa using p.marked_eq

end GaugeFiberPoint

/-- Marked and source fiber points are equivalent over every commutative test
algebra. -/
def markedFiberPointEquivGaugeFiberPoint
    {E β : K[X]} {pi b : K} (a : K) :
    MarkedFiberPoint E β pi b A ≃ GaugeFiberPoint E β pi b a A where
  toFun := fun p => p.toGaugeFiberPoint a
  invFun := GaugeFiberPoint.toMarkedFiberPoint
  left_inv := by
    intro p
    apply MarkedFiberPoint.ext
    exact p.chart.toSource_toChart (algebraMap K A a)
  right_inv := by
    intro p
    apply GaugeFiberPoint.ext
    exact p.source.toChart_toSource

/-- The broad functor-of-points theorem: roots of `E` are equivalent to full
quadratic-gauge source fiber points in every commutative `K`-algebra. -/
def rootEquivGaugeFiberPoint
    {E β : K[X]} {pi b : K} (a : K)
    (hE : E.Separable) (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β) :
    PolynomialRoot E A ≃ GaugeFiberPoint E β pi b a A :=
  (rootEquivMarkedFiberPoint hE g₁ hderiv).trans
    (markedFiberPointEquivGaugeFiberPoint a)

/-- Characteristic-zero squarefreeness gives exactly the public hypothesis of
the paper. -/
def squarefreeRootEquivGaugeFiberPoint [CharZero K]
    {E β : K[X]} {pi b : K} (a : K)
    (hE : Squarefree E) (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β) :
    PolynomialRoot E A ≃ GaugeFiberPoint E β pi b a A :=
  rootEquivGaugeFiberPoint a
    ((PerfectField.separable_iff_squarefree).2 hE) g₁ hderiv

end FiniteEtaleKeller
