/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GaugeAssembly
import FiniteEtaleKeller.Jacobian

/-!
# The all-degree quadratic-gauge polynomial map

This module packages the displayed construction from the paper as one
`MvPolynomial (Fin 3) K` map for an arbitrary seed polynomial `G`.

The preceding `GaugeAssembly` module proves the coefficientwise identities in
an arbitrary commutative ring.  Here the finite coefficient sums are assembled
into actual multivariate polynomials and their evaluations are exposed in a
form suitable for the represented-fiber bridge and later Jacobian and degree
certificates.
-/

noncomputable section

open Polynomial
open MvPolynomial

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

/-- Three-variable polynomial ring used by the universal quadratic gauge. -/
abbrev GaugePolynomial (K : Type*) [CommSemiring K] :=
  MvPolynomial (Fin 3) K

/-- The recurrent source polynomial `t = 1 + x*y`. -/
def generalGaugeT : GaugePolynomial K :=
  1 + MvPolynomial.X 0 * MvPolynomial.X 1

/-- The recurrent source polynomial
`q = t^2*z + (g₁/g₃)*y^2*(1+3*t)`. -/
def generalGaugeQ (G : K[X]) : GaugePolynomial K :=
  generalGaugeT ^ 2 * MvPolynomial.X 2 +
    MvPolynomial.C (G.coeff 1 / G.coeff 3) * MvPolynomial.X 1 ^ 2 *
      (1 + MvPolynomial.C 3 * generalGaugeT)

/-- The first displayed coordinate `Π = t*q`. -/
def generalGaugePi (G : K[X]) : GaugePolynomial K :=
  generalGaugeT * generalGaugeQ G

/-- The complete displayed second coordinate. -/
def generalGaugeB (G : K[X]) : GaugePolynomial K :=
  MvPolynomial.X 1 +
    MvPolynomial.C (3 * (G.coeff 3 / G.coeff 1)) *
      MvPolynomial.X 0 * generalGaugeQ G +
    MvPolynomial.C (2 * (G.coeff 2 / G.coeff 1)) *
      generalGaugeT * generalGaugeQ G +
    ∑ k ∈ Finset.Icc 4 G.natDegree,
      MvPolynomial.C ((k : K) * (G.coeff k / G.coeff 1)) *
        generalGaugeT ^ 2 * MvPolynomial.X 0 ^ (k - 2) *
          generalGaugeQ G ^ k

/-- The complete displayed third coordinate. -/
def generalGaugeC (G : K[X]) : GaugePolynomial K :=
  MvPolynomial.X 0 * (MvPolynomial.C 5 - MvPolynomial.C 3 * generalGaugeT) -
    MvPolynomial.C (G.coeff 3 / G.coeff 1) *
      MvPolynomial.X 0 ^ 3 * MvPolynomial.X 2 -
    ∑ k ∈ Finset.Icc 4 G.natDegree,
      MvPolynomial.C (((k - 2 : ℕ) : K) * (G.coeff k / G.coeff 1)) *
        (MvPolynomial.X 0 * generalGaugeQ G) ^ k

/-- The all-degree determinant-`-2` quadratic-gauge map. -/
def generalGaugeMap (G : K[X]) : Fin 3 → GaugePolynomial K :=
  ![generalGaugePi G, generalGaugeB G, generalGaugeC G]

/-- The target-preserving determinant-one output normalization. -/
def generalGaugeJacobianOneMap (G : K[X]) : Fin 3 → GaugePolynomial K :=
  scaleOutput 1 (-1 / 2 : K) 1 (generalGaugeMap G)

section Evaluation

variable {A : Type*} [CommRing A] [Algebra K A]

/-- Evaluate a polynomial map after extending its coefficients to a test
`K`-algebra. -/
def eval₂Map (F : Fin 3 → GaugePolynomial K) (p : Fin 3 → A) : Fin 3 → A :=
  fun i => MvPolynomial.eval₂ (algebraMap K A) p (F i)

@[simp]
theorem eval₂_generalGaugeT (p : Fin 3 → A) :
    MvPolynomial.eval₂ (algebraMap K A) p generalGaugeT =
      1 + p 0 * p 1 := by
  simp [generalGaugeT]

@[simp]
theorem eval₂_generalGaugeQ (G : K[X]) (p : Fin 3 → A) :
    MvPolynomial.eval₂ (algebraMap K A) p (generalGaugeQ G) =
      (1 + p 0 * p 1) ^ 2 * p 2 +
        algebraMap K A (G.coeff 1 / G.coeff 3) * p 1 ^ 2 *
          (1 + 3 * (1 + p 0 * p 1)) := by
  simp [generalGaugeQ, generalGaugeT]
  norm_num

@[simp]
theorem eval₂_generalGaugePi (G : K[X]) (p : Fin 3 → A) :
    MvPolynomial.eval₂ (algebraMap K A) p (generalGaugePi G) =
      (1 + p 0 * p 1) *
        ((1 + p 0 * p 1) ^ 2 * p 2 +
          algebraMap K A (G.coeff 1 / G.coeff 3) * p 1 ^ 2 *
            (1 + 3 * (1 + p 0 * p 1))) := by
  simp [generalGaugePi]

/-- Evaluation of the complete second coordinate, including the arbitrary
finite high-degree sum. -/
theorem eval₂_generalGaugeB (G : K[X]) (p : Fin 3 → A) :
    MvPolynomial.eval₂ (algebraMap K A) p (generalGaugeB G) =
      p 1 +
        algebraMap K A (3 * (G.coeff 3 / G.coeff 1)) * p 0 *
          MvPolynomial.eval₂ (algebraMap K A) p (generalGaugeQ G) +
        algebraMap K A (2 * (G.coeff 2 / G.coeff 1)) *
          (1 + p 0 * p 1) *
          MvPolynomial.eval₂ (algebraMap K A) p (generalGaugeQ G) +
        ∑ k ∈ Finset.Icc 4 G.natDegree,
          algebraMap K A ((k : K) * (G.coeff k / G.coeff 1)) *
            (1 + p 0 * p 1) ^ 2 * p 0 ^ (k - 2) *
              MvPolynomial.eval₂ (algebraMap K A) p (generalGaugeQ G) ^ k := by
  simp [generalGaugeB, generalGaugeT]

/-- Evaluation of the complete third coordinate, including the arbitrary
finite high-degree sum. -/
theorem eval₂_generalGaugeC (G : K[X]) (p : Fin 3 → A) :
    MvPolynomial.eval₂ (algebraMap K A) p (generalGaugeC G) =
      p 0 * (5 - 3 * (1 + p 0 * p 1)) -
        algebraMap K A (G.coeff 3 / G.coeff 1) * p 0 ^ 3 * p 2 -
        ∑ k ∈ Finset.Icc 4 G.natDegree,
          algebraMap K A (((k - 2 : ℕ) : K) * (G.coeff k / G.coeff 1)) *
            (p 0 * MvPolynomial.eval₂ (algebraMap K A) p (generalGaugeQ G)) ^ k := by
  simp [generalGaugeC, generalGaugeT]
  norm_num

@[simp]
theorem eval₂Map_generalGaugeMap_zero (G : K[X]) (p : Fin 3 → A) :
    eval₂Map (generalGaugeMap G) p 0 =
      MvPolynomial.eval₂ (algebraMap K A) p (generalGaugePi G) := by
  rfl

@[simp]
theorem eval₂Map_generalGaugeMap_one (G : K[X]) (p : Fin 3 → A) :
    eval₂Map (generalGaugeMap G) p 1 =
      MvPolynomial.eval₂ (algebraMap K A) p (generalGaugeB G) := by
  rfl

@[simp]
theorem eval₂Map_generalGaugeMap_two (G : K[X]) (p : Fin 3 → A) :
    eval₂Map (generalGaugeMap G) p 2 =
      MvPolynomial.eval₂ (algebraMap K A) p (generalGaugeC G) := by
  rfl

end Evaluation

#print axioms eval₂_generalGaugeB
#print axioms eval₂_generalGaugeC

end FiniteEtaleKeller
