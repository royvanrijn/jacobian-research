/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.SeparableReconstruction

/-!
# Roots as points of a polynomial quotient

For every commutative `K`-algebra `A`, roots of `E` in `A` are equivalent to
`K`-algebra homomorphisms `K[S]/(E) → A`.  This is the functor-of-points form
of the quotient universal property used by the fiber theorem.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K A B : Type*} [Field K]
variable [CommRing A] [Algebra K A]
variable [CommRing B] [Algebra K B]

/-- A root of `E` in a commutative `K`-algebra. -/
def PolynomialRoot (E : K[X]) (A : Type*) [CommRing A] [Algebra K A] :=
  { s : A // Polynomial.aeval s E = 0 }

namespace PolynomialRoot

variable {E : K[X]}

instance : Coe (PolynomialRoot E A) A := ⟨Subtype.val⟩

@[ext]
theorem ext {s t : PolynomialRoot E A} (h : (s : A) = t) : s = t :=
  Subtype.ext h

/-- A root defines the corresponding algebra homomorphism out of the quotient. -/
def liftAlgHom (s : PolynomialRoot E A) : AdjoinRoot E →ₐ[K] A :=
  AdjoinRoot.liftAlgHom E (Algebra.ofId K A) (s : A) s.property

@[simp]
theorem liftAlgHom_root (s : PolynomialRoot E A) :
    s.liftAlgHom (AdjoinRoot.root E) = (s : A) := by
  simp [liftAlgHom]

/-- An algebra homomorphism out of the quotient is determined by the image of
the distinguished root. -/
def ofAlgHom (f : AdjoinRoot E →ₐ[K] A) : PolynomialRoot E A :=
  ⟨f (AdjoinRoot.root E), AdjoinRoot.aeval_algHom_eq_zero E f⟩

@[simp]
theorem ofAlgHom_val (f : AdjoinRoot E →ₐ[K] A) :
    ((ofAlgHom f : PolynomialRoot E A) : A) = f (AdjoinRoot.root E) := rfl

/-- The universal property of `K[S]/(E)`, as an explicit equivalence of
points. -/
def algHomEquiv (E : K[X]) (A : Type*) [CommRing A] [Algebra K A] :
    (AdjoinRoot E →ₐ[K] A) ≃ PolynomialRoot E A where
  toFun := ofAlgHom
  invFun := liftAlgHom
  left_inv := by
    intro f
    apply AdjoinRoot.algHom_ext
    simp [ofAlgHom, liftAlgHom]
  right_inv := by
    intro s
    apply PolynomialRoot.ext
    simp [ofAlgHom, liftAlgHom]

@[simp]
theorem algHomEquiv_apply (f : AdjoinRoot E →ₐ[K] A) :
    algHomEquiv E A f = ofAlgHom f := rfl

@[simp]
theorem algHomEquiv_symm_apply (s : PolynomialRoot E A) :
    (algHomEquiv E A).symm s = s.liftAlgHom := rfl

/-- Roots are functorial in the test algebra. -/
def map (f : A →ₐ[K] B) (s : PolynomialRoot E A) : PolynomialRoot E B where
  val := f s
  property := by
    rw [Polynomial.aeval_algHom_apply f]
    simp [s.property]

@[simp]
theorem map_val (f : A →ₐ[K] B) (s : PolynomialRoot E A) :
    ((s.map f : PolynomialRoot E B) : B) = f s := rfl

@[simp]
theorem map_id (s : PolynomialRoot E A) :
    s.map (AlgHom.id K A) = s := by
  ext
  rfl

@[simp]
theorem map_comp {C : Type*} [CommRing C] [Algebra K C]
    (f : A →ₐ[K] B) (g : B →ₐ[K] C) (s : PolynomialRoot E A) :
    (s.map f).map g = s.map (g.comp f) := by
  ext
  rfl

/-- Evaluate the canonical normalized derivative unit at a root. -/
def normalizedDerivativeUnit
    (E : K[X]) (hE : E.Separable) (g₁ : Kˣ)
    (s : PolynomialRoot E A) : Aˣ :=
  Units.map s.liftAlgHom.toRingHom
    (normalizedDerivativeUnitOfSeparable E hE g₁)

@[simp]
theorem normalizedDerivativeUnit_val
    (E : K[X]) (hE : E.Separable) (g₁ : Kˣ)
    (s : PolynomialRoot E A) :
    (s.normalizedDerivativeUnit E hE g₁ : A) =
      algebraMap K A (↑g₁⁻¹ : K) * Polynomial.aeval (s : A) E.derivative := by
  rw [normalizedDerivativeUnit, Units.coe_map,
    normalizedDerivativeUnitOfSeparable_val]
  simp [liftAlgHom, Polynomial.aeval_def]

end PolynomialRoot

end FiniteEtaleKeller
