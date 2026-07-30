/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas
import Mathlib.RingTheory.TensorProduct.Finite
import Mathlib.RingTheory.TensorProduct.Maps

/-!
# Collision algebras of finite étale Keller fibers

For a finite fiber represented by a commutative `K`-algebra `A`, its ordered
self-collision fiber is represented by `A ⊗[K] A`.  The diagonal is induced
by multiplication, and its kernel is the collision obstruction ideal.

This module is independent of the quadratic-gauge presentation.  It applies
to every literal fiber equivalence produced by the finite-étale Keller
realization theorem.
-/

noncomputable section

namespace FiniteEtaleKeller

open scoped TensorProduct

universe u v w

/-- The coordinate algebra of the ordered self-collision fiber of
`Spec A → Spec K`. -/
abbrev FiberCollisionAlgebra
    (K : Type u) (A : Type v) [CommRing K] [CommRing A] [Algebra K A] :=
  A ⊗[K] A

/-- Multiplication restricts an ordered collision to the diagonal. -/
def fiberCollisionDiagonal
    (K : Type u) (A : Type v) [CommRing K] [CommRing A] [Algebra K A] :
    FiberCollisionAlgebra K A →ₐ[K] A :=
  Algebra.TensorProduct.productMap (AlgHom.id K A) (AlgHom.id K A)

@[simp]
theorem fiberCollisionDiagonal_tmul
    (K : Type u) (A : Type v) [CommRing K] [CommRing A] [Algebra K A]
    (a b : A) :
    fiberCollisionDiagonal K A (a ⊗ₜ b) = a * b :=
  rfl

/-- The scheme-theoretic obstruction to the collision fiber being only its
diagonal. -/
def fiberCollisionObstruction
    (K : Type u) (A : Type v) [CommRing K] [CommRing A] [Algebra K A] :
    Ideal (FiberCollisionAlgebra K A) :=
  RingHom.ker (fiberCollisionDiagonal K A).toRingHom

/-- The diagonal map is always onto: `a` is the image of `a ⊗ 1`. -/
theorem fiberCollisionDiagonal_surjective
    (K : Type u) (A : Type v) [CommRing K] [CommRing A] [Algebra K A] :
    Function.Surjective (fiberCollisionDiagonal K A) := by
  intro a
  exact ⟨a ⊗ₜ (1 : A), by simp⟩

/-- Maps from the collision algebra are exactly ordered pairs of maps from
the original fiber algebra.  This is the affine self-fiber-product universal
property. -/
def fiberCollisionPointPairsEquiv
    (K : Type u) (A : Type v) (R : Type w)
    [CommRing K] [CommRing A] [CommRing R] [Algebra K A] [Algebra K R] :
    (FiberCollisionAlgebra K A →ₐ[K] R) ≃
      (A →ₐ[K] R) × (A →ₐ[K] R) where
  toFun f :=
    (f.comp Algebra.TensorProduct.includeLeft,
      f.comp Algebra.TensorProduct.includeRight)
  invFun f :=
    Algebra.TensorProduct.lift f.1 f.2 (fun _ _ ↦ .all _ _)
  left_inv f := by
    ext <;> simp
  right_inv f := by
    ext <;> simp

/-- A fiber of rank greater than one has a nonzero collision obstruction.
Equivalently, multiplication on its tensor square cannot be injective. -/
theorem fiberCollisionDiagonal_not_injective
    (K : Type u) (A : Type v) [Field K] [CommRing A] [Algebra K A]
    [Module.Finite K A] (hdeg : 1 < Module.finrank K A) :
    ¬ Function.Injective (fiberCollisionDiagonal K A) := by
  intro hinj
  have hle :=
    (fiberCollisionDiagonal K A).toLinearMap.finrank_le_finrank_of_injective hinj
  rw [Module.finrank_tensorProduct] at hle
  nlinarith

/-- The obstruction has the expected dimension `N²-N`: the collision algebra
has dimension `N²`, and the diagonal multiplication map is surjective onto a
rank-`N` algebra. -/
theorem finrank_fiberCollisionObstruction
    (K : Type u) (A : Type v) [Field K] [CommRing A] [Nontrivial A]
    [Algebra K A] [Module.Finite K A] :
    Module.finrank K (fiberCollisionObstruction K A) =
      Module.finrank K A ^ 2 - Module.finrank K A := by
  let f := (fiberCollisionDiagonal K A).toLinearMap
  change Module.finrank K (LinearMap.ker f) =
    Module.finrank K A ^ 2 - Module.finrank K A
  have hrange : LinearMap.range f = ⊤ :=
    LinearMap.range_eq_top.mpr (fiberCollisionDiagonal_surjective K A)
  have h := f.finrank_range_add_finrank_ker
  rw [hrange, finrank_top, Module.finrank_tensorProduct] at h
  simp only [pow_two]
  omega

#print axioms fiberCollisionDiagonal_surjective
#print axioms fiberCollisionPointPairsEquiv
#print axioms fiberCollisionDiagonal_not_injective
#print axioms finrank_fiberCollisionObstruction

end FiniteEtaleKeller
