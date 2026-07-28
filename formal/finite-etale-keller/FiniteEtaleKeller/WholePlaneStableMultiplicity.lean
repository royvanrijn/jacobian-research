/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.StableGaugeFiber
import Mathlib.RingTheory.Polynomial.Resultant.Basic
import Mathlib.RingTheory.Localization.Away.AdjoinRoot

/-!
# Relative whole-plane stable multiplicity

This file upgrades the pointwise stable-fiber equivalences to a relative
functor-of-points statement.  A `GaugePlaneRestrictionPoint` is a point of
the common source divisor `Π = 1`.  Evaluating the last two coordinates of a
map gives its morphism from this divisor to the `(B,C)` target plane.

For every test algebra, all common power shifts induce the same target
morphism, as do all cubic lifts.  The equality is natural in the test
algebra.  It is also restricted to the principal open cut out by the
discriminant of the universal inverse polynomial, yielding one common
discriminant-open relative family rather than separate equivalences on
individual scalar fibers.
-/

noncomputable section

open Polynomial
open MvPolynomial

namespace FiniteEtaleKeller

variable {K A B : Type*} [Field K] [CharZero K]
variable [CommRing A] [Algebra K A]
variable [CommRing B] [Algebra K B]

/-- The coordinate ring of the `(B,C)` target plane. -/
abbrev GaugeTargetPlanePolynomial (K : Type*) [CommSemiring K] :=
  MvPolynomial (Fin 2) K

/-- The inverse equation over the whole target plane `Π = 1`. -/
def generalGaugePlaneInversePolynomial (G : K[X]) :
    (GaugeTargetPlanePolynomial K)[X] :=
  G.map (MvPolynomial.C : K →+* GaugeTargetPlanePolynomial K) -
    Polynomial.C (MvPolynomial.C (G.coeff 1 / 2)) *
      (Polynomial.C (MvPolynomial.X 0) * Polynomial.X ^ 2 +
        Polynomial.C (MvPolynomial.X 1))

/-- The discriminant cutting out the squarefree inverse locus in the target
plane. -/
def generalGaugePlaneDiscriminant (G : K[X]) :
    GaugeTargetPlanePolynomial K :=
  (generalGaugePlaneInversePolynomial G).discr

/-- The coordinate ring of the squarefree target-plane open. -/
abbrev GaugeTargetPlaneDiscriminantOpenRing (G : K[X]) :=
  Localization.Away (generalGaugePlaneDiscriminant G)

/-- The universal inverse equation after restriction to the discriminant
open. -/
def generalGaugePlaneOpenInversePolynomial (G : K[X]) :
    (GaugeTargetPlaneDiscriminantOpenRing G)[X] :=
  (generalGaugePlaneInversePolynomial G).map
    (algebraMap (GaugeTargetPlanePolynomial K)
      (GaugeTargetPlaneDiscriminantOpenRing G))

/-- The affine coordinate algebra of the universal root cover over the
discriminant-open target plane. -/
abbrev GeneralGaugePlaneRootCoverAlgebra (G : K[X]) :=
  AdjoinRoot (generalGaugePlaneOpenInversePolynomial G)

/-- The universal discriminant is a unit on its principal open. -/
theorem generalGaugePlaneDiscriminant_isUnit :
    IsUnit
      (algebraMap (GaugeTargetPlanePolynomial K)
        (GaugeTargetPlaneDiscriminantOpenRing G)
        (generalGaugePlaneDiscriminant G)) :=
  IsLocalization.Away.algebraMap_isUnit
    (S := GaugeTargetPlaneDiscriminantOpenRing G)
    (generalGaugePlaneDiscriminant G)

/-- Evaluation of the universal target-plane coefficients at a scalar target
`(b,c)`. -/
def gaugeTargetPlaneSpecialization (b c : K) :
    GaugeTargetPlanePolynomial K →+* K :=
  MvPolynomial.eval₂Hom (RingHom.id K) ![b, c]

/-- Specializing the universal inverse equation at `(b,c)` gives the scalar
inverse polynomial used by the fiberwise reconstruction theorem. -/
theorem generalGaugePlaneInversePolynomial_specialize
    (G : K[X]) (hzero : G.coeff 0 = 0) (b c : K) :
    (generalGaugePlaneInversePolynomial G).map
        (gaugeTargetPlaneSpecialization b c) =
      generalGaugeInversePolynomial G 1 b c := by
  rw [generalGaugeInversePolynomial,
    generalGaugeSeedPolynomial_one_eq G hzero]
  simp [generalGaugePlaneInversePolynomial,
    gaugeTargetPlaneSpecialization, Polynomial.map_map]

/-- A point of the common source divisor `Π = 1`.  The type deliberately
depends only on `G`, not on a shift parameter: all maps in either stable
family have literally the same first coordinate. -/
@[ext]
structure GaugePlaneRestrictionPoint
    (G : K[X]) (A : Type*) [CommRing A] [Algebra K A] where
  point : Fin 3 → A
  first_eq :
    MvPolynomial.eval₂ (algebraMap K A) point (generalGaugePi G) = 1

private theorem eval₂_map_algHom_relative
    {σ : Type*} (f : A →ₐ[K] B) (point : σ → A)
    (P : MvPolynomial σ K) :
    MvPolynomial.eval₂ (algebraMap K B) (fun i => f (point i)) P =
      f (MvPolynomial.eval₂ (algebraMap K A) point P) := by
  have hcomp : f.toRingHom.comp (algebraMap K A) = algebraMap K B := by
    ext r
    exact f.commutes r
  calc
    MvPolynomial.eval₂ (algebraMap K B) (fun i => f (point i)) P =
        MvPolynomial.eval₂ (f.toRingHom.comp (algebraMap K A))
          (fun i => f (point i)) P := by rw [hcomp]
    _ = f (MvPolynomial.eval₂ (algebraMap K A) point P) :=
      (MvPolynomial.hom_eval₂ P (algebraMap K A) f.toRingHom point).symm

namespace GaugePlaneRestrictionPoint

variable {G : K[X]}

/-- The relative `(B,C)` target of a restricted source point under `F`. -/
def target (F : Fin 3 → GaugePolynomial K)
    (p : GaugePlaneRestrictionPoint G A) : Fin 2 → A :=
  ![
    MvPolynomial.eval₂ (algebraMap K A) p.point (F 1),
    MvPolynomial.eval₂ (algebraMap K A) p.point (F 2)
  ]

/-- The common source divisor is functorial in the test algebra. -/
def map (f : A →ₐ[K] B) (p : GaugePlaneRestrictionPoint G A) :
    GaugePlaneRestrictionPoint G B where
  point := fun i => f (p.point i)
  first_eq := by
    rw [eval₂_map_algHom_relative f p.point, p.first_eq, map_one]

@[simp]
theorem map_point (f : A →ₐ[K] B)
    (p : GaugePlaneRestrictionPoint G A) (i : Fin 3) :
    (p.map f).point i = f (p.point i) := rfl

/-- Every relative target map is natural in the test algebra. -/
theorem target_map (F : Fin 3 → GaugePolynomial K)
    (f : A →ₐ[K] B) (p : GaugePlaneRestrictionPoint G A) :
    (p.map f).target F = fun i => f (p.target F i) := by
  funext i
  fin_cases i <;>
    simp [target, map, eval₂_map_algHom_relative]

end GaugePlaneRestrictionPoint

/-- All common power shifts induce the same morphism from the complete source
divisor `Π = 1` to the `(B,C)` target plane. -/
theorem powerShiftedGaugeWholePlane_target
    (G : K[X]) (m m' : ℕ)
    (p : GaugePlaneRestrictionPoint G A) :
    p.target (powerShiftedGaugeMap G m) =
      p.target (powerShiftedGaugeMap G m') := by
  have hm := eval₂_powerShiftedGaugeMap_eq G m p.point p.first_eq
  have hm' := eval₂_powerShiftedGaugeMap_eq G m' p.point p.first_eq
  funext i
  fin_cases i
  · change eval₂Map (powerShiftedGaugeMap G m) p.point 1 =
      eval₂Map (powerShiftedGaugeMap G m') p.point 1
    exact (congrFun hm 1).trans (congrFun hm' 1).symm
  · change eval₂Map (powerShiftedGaugeMap G m) p.point 2 =
      eval₂Map (powerShiftedGaugeMap G m') p.point 2
    exact (congrFun hm 2).trans (congrFun hm' 2).symm

/-- The whole-plane power-shift equality commutes with change of test
algebra. -/
theorem powerShiftedGaugeWholePlane_target_natural
    (G : K[X]) (m m' : ℕ) (f : A →ₐ[K] B)
    (p : GaugePlaneRestrictionPoint G A) :
    (GaugePlaneRestrictionPoint.map f p).target
        (powerShiftedGaugeMap G m) =
      (GaugePlaneRestrictionPoint.map f p).target
        (powerShiftedGaugeMap G m') :=
  powerShiftedGaugeWholePlane_target G m m'
    (GaugePlaneRestrictionPoint.map f p)

/-- The relative source functor equivalence for two power shifts.  It is
literally the identity because the common source divisor is unchanged. -/
def powerShiftedGaugeWholePlaneEquiv
    (G : K[X]) (_m _m' : ℕ) :
    GaugePlaneRestrictionPoint G A ≃ GaugePlaneRestrictionPoint G A :=
  Equiv.refl _

/-- The identity source equivalence lies over the common target morphism. -/
theorem powerShiftedGaugeWholePlaneEquiv_overTarget
    (G : K[X]) (m m' : ℕ)
    (p : GaugePlaneRestrictionPoint G A) :
    (powerShiftedGaugeWholePlaneEquiv (A := A) G m m' p).target
        (powerShiftedGaugeMap G m') =
      p.target (powerShiftedGaugeMap G m) := by
  exact (powerShiftedGaugeWholePlane_target G m m' p).symm

/-- The relative power-shift equivalence is natural in the test algebra. -/
theorem powerShiftedGaugeWholePlaneEquiv_natural
    (G : K[X]) (m m' : ℕ) (f : A →ₐ[K] B)
    (p : GaugePlaneRestrictionPoint G A) :
    GaugePlaneRestrictionPoint.map f
        (powerShiftedGaugeWholePlaneEquiv (A := A) G m m' p) =
      powerShiftedGaugeWholePlaneEquiv
        (A := B) G m m' (GaugePlaneRestrictionPoint.map f p) :=
  rfl

/-- All cubic lifts induce the same morphism from the complete source divisor
`Π = 1` to the `(B,C)` target plane. -/
theorem cubicLiftGaugeWholePlane_target
    (G : K[X]) (n n' : ℕ) (hn : 4 ≤ n) (hn' : 4 ≤ n')
    (p : GaugePlaneRestrictionPoint G A) :
    p.target (cubicLiftGaugeMap G n) =
      p.target (cubicLiftGaugeMap G n') := by
  have hnEval := eval₂_cubicLiftGaugeMap_eq G n hn p.point p.first_eq
  have hnEval' := eval₂_cubicLiftGaugeMap_eq G n' hn' p.point p.first_eq
  funext i
  fin_cases i
  · change eval₂Map (cubicLiftGaugeMap G n) p.point 1 =
      eval₂Map (cubicLiftGaugeMap G n') p.point 1
    exact (congrFun hnEval 1).trans (congrFun hnEval' 1).symm
  · change eval₂Map (cubicLiftGaugeMap G n) p.point 2 =
      eval₂Map (cubicLiftGaugeMap G n') p.point 2
    exact (congrFun hnEval 2).trans (congrFun hnEval' 2).symm

/-- The whole-plane cubic equality commutes with change of test algebra. -/
theorem cubicLiftGaugeWholePlane_target_natural
    (G : K[X]) (n n' : ℕ) (hn : 4 ≤ n) (hn' : 4 ≤ n')
    (f : A →ₐ[K] B) (p : GaugePlaneRestrictionPoint G A) :
    (GaugePlaneRestrictionPoint.map f p).target
        (cubicLiftGaugeMap G n) =
      (GaugePlaneRestrictionPoint.map f p).target
        (cubicLiftGaugeMap G n') :=
  cubicLiftGaugeWholePlane_target G n n' hn hn'
    (GaugePlaneRestrictionPoint.map f p)

/-- The relative source functor equivalence for two cubic lifts. -/
def cubicLiftGaugeWholePlaneEquiv
    (G : K[X]) (_n _n' : ℕ) :
    GaugePlaneRestrictionPoint G A ≃ GaugePlaneRestrictionPoint G A :=
  Equiv.refl _

/-- The cubic identity source equivalence lies over the common target
morphism. -/
theorem cubicLiftGaugeWholePlaneEquiv_overTarget
    (G : K[X]) (n n' : ℕ) (hn : 4 ≤ n) (hn' : 4 ≤ n')
    (p : GaugePlaneRestrictionPoint G A) :
    (cubicLiftGaugeWholePlaneEquiv (A := A) G n n' p).target
        (cubicLiftGaugeMap G n') =
      p.target (cubicLiftGaugeMap G n) := by
  exact (cubicLiftGaugeWholePlane_target G n n' hn hn' p).symm

/-- The relative cubic equivalence is natural in the test algebra. -/
theorem cubicLiftGaugeWholePlaneEquiv_natural
    (G : K[X]) (n n' : ℕ) (f : A →ₐ[K] B)
    (p : GaugePlaneRestrictionPoint G A) :
    GaugePlaneRestrictionPoint.map f
        (cubicLiftGaugeWholePlaneEquiv (A := A) G n n' p) =
      cubicLiftGaugeWholePlaneEquiv
        (A := B) G n n' (GaugePlaneRestrictionPoint.map f p) :=
  rfl

/-- Evaluation of a target-plane polynomial at a relative `(B,C)` point. -/
def evalGaugeTarget
    (u : Fin 2 → A) (D : GaugeTargetPlanePolynomial K) : A :=
  MvPolynomial.eval₂ (algebraMap K A) u D

/-- The common restriction over the discriminant principal open.  The
undeformed target is used in the definition; the whole-plane theorems prove
that every stable deformation gives exactly the same unit condition. -/
def GaugePlaneDiscriminantOpenPoint
    (G : K[X]) (A : Type*) [CommRing A] [Algebra K A] :=
  { p : GaugePlaneRestrictionPoint G A //
    IsUnit (evalGaugeTarget (p.target (generalGaugeMap G))
      (generalGaugePlaneDiscriminant G)) }

namespace GaugePlaneDiscriminantOpenPoint

variable {G : K[X]}

/-- The discriminant-open source restriction is functorial. -/
def map (f : A →ₐ[K] B)
    (p : GaugePlaneDiscriminantOpenPoint G A) :
    GaugePlaneDiscriminantOpenPoint G B := by
  refine ⟨GaugePlaneRestrictionPoint.map f p.1, ?_⟩
  rw [GaugePlaneRestrictionPoint.target_map]
  change IsUnit
    (MvPolynomial.eval₂ (algebraMap K B)
      (fun i => f (p.1.target (generalGaugeMap G) i))
      (generalGaugePlaneDiscriminant G))
  rw [eval₂_map_algHom_relative f
    (p.1.target (generalGaugeMap G)) (generalGaugePlaneDiscriminant G)]
  exact p.2.map f

end GaugePlaneDiscriminantOpenPoint

/-- Every power shift has a unit discriminant at every point of the common
discriminant-open restriction. -/
theorem powerShiftedGaugeDiscriminantOpen_isUnit
    (G : K[X]) (m : ℕ)
    (p : GaugePlaneDiscriminantOpenPoint G A) :
    IsUnit (evalGaugeTarget
      (p.1.target (powerShiftedGaugeMap G m))
      (generalGaugePlaneDiscriminant G)) := by
  rw [powerShiftedGaugeWholePlane_target G m 0 p.1]
  have hzero :
      p.1.target (powerShiftedGaugeMap G 0) =
        p.1.target (generalGaugeMap G) := by
    have heval :=
      eval₂_powerShiftedGaugeMap_eq G 0 p.1.point p.1.first_eq
    funext i
    fin_cases i
    · exact congrFun heval 1
    · exact congrFun heval 2
  rw [hzero]
  exact p.2

/-- The common discriminant-open relative source equivalence for power
shifts. -/
def powerShiftedGaugeDiscriminantOpenEquiv
    (G : K[X]) (_m _m' : ℕ) :
    GaugePlaneDiscriminantOpenPoint G A ≃
      GaugePlaneDiscriminantOpenPoint G A :=
  Equiv.refl _

/-- The discriminant-open power-shift equivalence lies over the common target
morphism. -/
theorem powerShiftedGaugeDiscriminantOpenEquiv_overTarget
    (G : K[X]) (m m' : ℕ)
    (p : GaugePlaneDiscriminantOpenPoint G A) :
    (powerShiftedGaugeDiscriminantOpenEquiv (A := A) G m m' p).1.target
        (powerShiftedGaugeMap G m') =
      p.1.target (powerShiftedGaugeMap G m) := by
  exact (powerShiftedGaugeWholePlane_target G m m' p.1).symm

/-- The discriminant-open power-shift equivalence is natural. -/
theorem powerShiftedGaugeDiscriminantOpenEquiv_natural
    (G : K[X]) (m m' : ℕ) (f : A →ₐ[K] B)
    (p : GaugePlaneDiscriminantOpenPoint G A) :
    GaugePlaneDiscriminantOpenPoint.map f
        (powerShiftedGaugeDiscriminantOpenEquiv (A := A) G m m' p) =
      powerShiftedGaugeDiscriminantOpenEquiv
        (A := B) G m m' (GaugePlaneDiscriminantOpenPoint.map f p) :=
  rfl

/-- Every cubic lift has a unit discriminant at every point of the common
discriminant-open restriction. -/
theorem cubicLiftGaugeDiscriminantOpen_isUnit
    (G : K[X]) (n : ℕ) (hn : 4 ≤ n)
    (p : GaugePlaneDiscriminantOpenPoint G A) :
    IsUnit (evalGaugeTarget
      (p.1.target (cubicLiftGaugeMap G n))
      (generalGaugePlaneDiscriminant G)) := by
  have htarget :
      p.1.target (cubicLiftGaugeMap G n) =
        p.1.target (generalGaugeMap G) := by
    have heval :=
      eval₂_cubicLiftGaugeMap_eq G n hn p.1.point p.1.first_eq
    funext i
    fin_cases i
    · exact congrFun heval 1
    · exact congrFun heval 2
  rw [htarget]
  exact p.2

/-- The common discriminant-open relative source equivalence for cubic
lifts. -/
def cubicLiftGaugeDiscriminantOpenEquiv
    (G : K[X]) (_n _n' : ℕ) :
    GaugePlaneDiscriminantOpenPoint G A ≃
      GaugePlaneDiscriminantOpenPoint G A :=
  Equiv.refl _

/-- The discriminant-open cubic equivalence lies over the common target
morphism. -/
theorem cubicLiftGaugeDiscriminantOpenEquiv_overTarget
    (G : K[X]) (n n' : ℕ) (hn : 4 ≤ n) (hn' : 4 ≤ n')
    (p : GaugePlaneDiscriminantOpenPoint G A) :
    (cubicLiftGaugeDiscriminantOpenEquiv (A := A) G n n' p).1.target
        (cubicLiftGaugeMap G n') =
      p.1.target (cubicLiftGaugeMap G n) := by
  exact (cubicLiftGaugeWholePlane_target G n n' hn hn' p.1).symm

/-- The discriminant-open cubic equivalence is natural. -/
theorem cubicLiftGaugeDiscriminantOpenEquiv_natural
    (G : K[X]) (n n' : ℕ) (f : A →ₐ[K] B)
    (p : GaugePlaneDiscriminantOpenPoint G A) :
    GaugePlaneDiscriminantOpenPoint.map f
        (cubicLiftGaugeDiscriminantOpenEquiv (A := A) G n n' p) =
      cubicLiftGaugeDiscriminantOpenEquiv
        (A := B) G n n' (GaugePlaneDiscriminantOpenPoint.map f p) :=
  rfl

#print axioms powerShiftedGaugeWholePlane_target
#print axioms generalGaugePlaneInversePolynomial_specialize
#print axioms generalGaugePlaneDiscriminant_isUnit
#print axioms powerShiftedGaugeWholePlaneEquiv_overTarget
#print axioms powerShiftedGaugeWholePlaneEquiv_natural
#print axioms powerShiftedGaugeDiscriminantOpenEquiv_overTarget
#print axioms powerShiftedGaugeDiscriminantOpenEquiv_natural
#print axioms cubicLiftGaugeWholePlane_target
#print axioms cubicLiftGaugeWholePlaneEquiv_overTarget
#print axioms cubicLiftGaugeWholePlaneEquiv_natural
#print axioms cubicLiftGaugeDiscriminantOpenEquiv_overTarget
#print axioms cubicLiftGaugeDiscriminantOpenEquiv_natural

end FiniteEtaleKeller
