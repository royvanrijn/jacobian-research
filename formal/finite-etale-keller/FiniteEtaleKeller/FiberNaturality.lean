/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GaugeFiberPoints

/-!
# Naturality of the fiber equivalence

The root/source-point equivalence commutes with morphisms of commutative test
algebras.  This is the precise functor-of-points statement needed to pass from
the explicit formulas to the represented scheme isomorphism.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K A B : Type*} [Field K]
variable [CommRing A] [Algebra K A]
variable [CommRing B] [Algebra K B]

namespace GaugeChart

variable {pi : K}

/-- Marked charts are functorial in the test algebra. -/
def map (f : A →ₐ[K] B)
    (p : GaugeChart A (algebraMap K A pi)) :
    GaugeChart B (algebraMap K B pi) where
  S := f p.S
  Q := f p.Q
  d := Units.map f.toRingHom p.d
  chart_eq := by
    change f (p.d : A) =
      1 - f p.S * f p.Q + algebraMap K B pi * (f p.S) ^ 2
    simpa using congrArg f p.chart_eq

@[simp]
theorem map_S (f : A →ₐ[K] B)
    (p : GaugeChart A (algebraMap K A pi)) : (p.map f).S = f p.S := rfl

@[simp]
theorem map_Q (f : A →ₐ[K] B)
    (p : GaugeChart A (algebraMap K A pi)) : (p.map f).Q = f p.Q := rfl

@[simp]
theorem map_d (f : A →ₐ[K] B)
    (p : GaugeChart A (algebraMap K A pi)) :
    (p.map f).d = Units.map f.toRingHom p.d := rfl

end GaugeChart

namespace GaugeSource

variable {pi a : K}

/-- Source-chart data are functorial in the test algebra. -/
def map (f : A →ₐ[K] B)
    (p : GaugeSource A (algebraMap K A pi) (algebraMap K A a)) :
    GaugeSource B (algebraMap K B pi) (algebraMap K B a) where
  t := Units.map f.toRingHom p.t
  x := f p.x
  y := f p.y
  z := f p.z
  t_eq := by
    change f (p.t : A) = 1 + f p.x * f p.y
    simpa using congrArg f p.t_eq
  pi_eq := by
    change f (p.t : A) *
        (f (p.t : A) ^ 2 * f p.z
          + algebraMap K B a * (f p.y) ^ 2 * (1 + 3 * f (p.t : A))) =
      algebraMap K B pi
    simpa only [map_mul, map_add, map_pow, map_one, map_ofNat,
      AlgHom.commutes] using congrArg f p.pi_eq

@[simp]
theorem map_t (f : A →ₐ[K] B)
    (p : GaugeSource A (algebraMap K A pi) (algebraMap K A a)) :
    (p.map f).t = Units.map f.toRingHom p.t := rfl

@[simp]
theorem map_x (f : A →ₐ[K] B)
    (p : GaugeSource A (algebraMap K A pi) (algebraMap K A a)) :
    (p.map f).x = f p.x := rfl

@[simp]
theorem map_y (f : A →ₐ[K] B)
    (p : GaugeSource A (algebraMap K A pi) (algebraMap K A a)) :
    (p.map f).y = f p.y := rfl

@[simp]
theorem map_z (f : A →ₐ[K] B)
    (p : GaugeSource A (algebraMap K A pi) (algebraMap K A a)) :
    (p.map f).z = f p.z := rfl

@[simp]
theorem map_q (f : A →ₐ[K] B)
    (p : GaugeSource A (algebraMap K A pi) (algebraMap K A a)) :
    (p.map f).q = f p.q := by
  simp [GaugeSource.map, GaugeSource.q, map_ofNat]

@[simp]
theorem map_S (f : A →ₐ[K] B)
    (p : GaugeSource A (algebraMap K A pi) (algebraMap K A a)) :
    (p.map f).S = f p.S := by
  simp [GaugeSource.map, GaugeSource.S]

@[simp]
theorem map_Q (f : A →ₐ[K] B)
    (p : GaugeSource A (algebraMap K A pi) (algebraMap K A a)) :
    (p.map f).Q = f p.Q := by
  simp [GaugeSource.Q]

end GaugeSource

namespace GaugeFiberPoint

variable {E β : K[X]} {pi b a : K}

/-- Forgetting a full source fiber point leaves its polynomial root. -/
def toRoot (p : GaugeFiberPoint E β pi b a A) : PolynomialRoot E A :=
  ⟨p.source.S, p.root_eq⟩

@[simp]
theorem toRoot_val (p : GaugeFiberPoint E β pi b a A) :
    p.toRoot.1 = p.source.S := rfl

/-- Full source fiber points are functorial in the test algebra. -/
def map (f : A →ₐ[K] B) (p : GaugeFiberPoint E β pi b a A) :
    GaugeFiberPoint E β pi b a B where
  source := p.source.map f
  root_eq := by
    rw [GaugeSource.map_S, Polynomial.aeval_algHom_apply f]
    rw [p.root_eq, map_zero]
  marked_eq := by
    rw [GaugeSource.map_Q, GaugeSource.map_S,
      Polynomial.aeval_algHom_apply f]
    calc
      f p.source.Q + f (Polynomial.aeval p.source.S β) =
          f (p.source.Q + Polynomial.aeval p.source.S β) := by rw [map_add]
      _ = f (algebraMap K A b) := congrArg f p.marked_eq
      _ = algebraMap K B b := f.commutes b

@[simp]
theorem map_toRoot (f : A →ₐ[K] B)
    (p : GaugeFiberPoint E β pi b a A) :
    (p.map f).toRoot = p.toRoot.map f := by
  apply PolynomialRoot.ext
  change (p.source.map f).S = f p.source.S
  exact GaugeSource.map_S f p.source

end GaugeFiberPoint

@[simp]
theorem rootEquivGaugeFiberPoint_symm_apply
    {E β : K[X]} {pi b : K} (a : K)
    (hE : E.Separable) (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β)
    (p : GaugeFiberPoint E β pi b a A) :
    (rootEquivGaugeFiberPoint (A := A) a hE g₁ hderiv).symm p = p.toRoot := rfl

@[simp]
theorem rootEquivGaugeFiberPoint_toRoot_apply
    {E β : K[X]} {pi b : K} (a : K)
    (hE : E.Separable) (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β)
    (s : PolynomialRoot E A) :
    (rootEquivGaugeFiberPoint (A := A) a hE g₁ hderiv s).toRoot = s := by
  have h :=
    (rootEquivGaugeFiberPoint (A := A) a hE g₁ hderiv).symm_apply_apply s
  simpa using h

/-- The root/source-point equivalence is natural in every commutative test
`K`-algebra. -/
theorem rootEquivGaugeFiberPoint_natural
    {E β : K[X]} {pi b : K} (a : K)
    (hE : E.Separable) (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β)
    (f : A →ₐ[K] B) (s : PolynomialRoot E A) :
    GaugeFiberPoint.map f
        (rootEquivGaugeFiberPoint (A := A) a hE g₁ hderiv s) =
      rootEquivGaugeFiberPoint (A := B) a hE g₁ hderiv (s.map f) := by
  apply (rootEquivGaugeFiberPoint (A := B) a hE g₁ hderiv).symm.injective
  rw [rootEquivGaugeFiberPoint_symm_apply,
    GaugeFiberPoint.map_toRoot,
    rootEquivGaugeFiberPoint_toRoot_apply,
    rootEquivGaugeFiberPoint_symm_apply,
    rootEquivGaugeFiberPoint_toRoot_apply]

/-- The naturality theorem in the characteristic-zero squarefree form used by
the paper. -/
theorem squarefreeRootEquivGaugeFiberPoint_natural [CharZero K]
    {E β : K[X]} {pi b : K} (a : K)
    (hE : Squarefree E) (g₁ : Kˣ)
    (hderiv : E.derivative = C (g₁ : K) * markedChartPolynomial pi b β)
    (f : A →ₐ[K] B) (s : PolynomialRoot E A) :
    GaugeFiberPoint.map f
        (squarefreeRootEquivGaugeFiberPoint (A := A) a hE g₁ hderiv s) =
      squarefreeRootEquivGaugeFiberPoint (A := B) a hE g₁ hderiv (s.map f) := by
  exact rootEquivGaugeFiberPoint_natural a
    ((PerfectField.separable_iff_squarefree).2 hE) g₁ hderiv f s

end FiniteEtaleKeller
