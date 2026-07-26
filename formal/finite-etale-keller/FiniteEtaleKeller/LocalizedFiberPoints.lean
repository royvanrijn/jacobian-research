/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeRawFiber
import Mathlib.RingTheory.Localization.Away.AdjoinRoot

/-!
# Localized inverse-polynomial fibers

For an arbitrary inverse polynomial `E`, a point of an étale polynomial-map
fiber is a root of `E` at which `E'` is invertible.  This module packages that
condition and proves that it is represented by

`(K[S] / (E))[1 / E'] = Localization.Away (AdjoinRoot.mk E E.derivative)`.

No separability assumption is used here.  When `E` is separable, the derivative
is already a unit modulo `E`, so the localization reduces to `AdjoinRoot E` and
the earlier squarefree fiber theorem is recovered.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller

variable {K A B : Type*} [Field K] [CharZero K]
variable [CommRing A] [Algebra K A]
variable [CommRing B] [Algebra K B]

/-- A root of `E` together with the unit represented by its derivative.  The
unit is bundled because the reconstruction formulas need its inverse. -/
structure LocalizedPolynomialRoot (E : K[X])
    (A : Type*) [CommRing A] [Algebra K A] where
  val : A
  root_eq : Polynomial.aeval val E = 0
  derivativeUnit : Aˣ
  derivativeUnit_eq : (derivativeUnit : A) = Polynomial.aeval val E.derivative

namespace LocalizedPolynomialRoot

variable {E : K[X]}

@[ext]
theorem ext {s t : LocalizedPolynomialRoot E A} (h : s.val = t.val) : s = t := by
  cases s with
  | mk sval sroot sunit sequnit =>
      cases t with
      | mk tval troot tunit tequnit =>
          cases h
          have hu : sunit = tunit := by
            apply Units.ext
            rw [sequnit, tequnit]
          cases hu
          rfl

/-- Localized roots are functorial in the test algebra. -/
def map (f : A →ₐ[K] B) (s : LocalizedPolynomialRoot E A) :
    LocalizedPolynomialRoot E B where
  val := f s.val
  root_eq := by
    rw [Polynomial.aeval_algHom_apply f s.val E, s.root_eq, map_zero]
  derivativeUnit := Units.map f.toRingHom s.derivativeUnit
  derivativeUnit_eq := by
    simp only [Units.coe_map]
    rw [s.derivativeUnit_eq, Polynomial.aeval_algHom_apply]
    rfl

@[simp]
theorem map_val (f : A →ₐ[K] B) (s : LocalizedPolynomialRoot E A) :
    (s.map f).val = f s.val := rfl

/-- Forgetting the derivative unit leaves the underlying root. -/
def toPolynomialRoot (s : LocalizedPolynomialRoot E A) : PolynomialRoot E A :=
  ⟨s.val, s.root_eq⟩

@[simp]
theorem toPolynomialRoot_val (s : LocalizedPolynomialRoot E A) :
    s.toPolynomialRoot.1 = s.val := rfl

end LocalizedPolynomialRoot

/-- The class of `E'` in `K[S]/(E)`. -/
def derivativeClass (E : K[X]) : AdjoinRoot E :=
  AdjoinRoot.mk E E.derivative

/-- The quotient by `E`, localized away from the derivative class. -/
abbrev LocalizedAdjoinRoot (E : K[X]) :=
  Localization.Away (derivativeClass E)

private def localizedBaseAlgHom {E : K[X]}
    (φ : LocalizedAdjoinRoot E →ₐ[K] A) : AdjoinRoot E →ₐ[K] A :=
  φ.comp (Algebra.algHom K (AdjoinRoot E) (LocalizedAdjoinRoot E))

namespace LocalizedPolynomialRoot

variable {E : K[X]}

/-- A map out of the localized quotient gives a root whose derivative is a
unit. -/
def ofLocalizedAlgHom
    (φ : LocalizedAdjoinRoot E →ₐ[K] A) : LocalizedPolynomialRoot E A where
  val := φ (algebraMap (AdjoinRoot E) (LocalizedAdjoinRoot E) (AdjoinRoot.root E))
  root_eq := AdjoinRoot.aeval_algHom_eq_zero E (localizedBaseAlgHom φ)
  derivativeUnit :=
    Units.map φ.toRingHom
      (IsLocalization.Away.algebraMap_isUnit (derivativeClass E)).unit
  derivativeUnit_eq := by
    simp only [Units.coe_map, IsUnit.unit_spec]
    change
      (localizedBaseAlgHom φ) (AdjoinRoot.mk E E.derivative) =
        Polynomial.aeval
          ((localizedBaseAlgHom φ) (AdjoinRoot.root E)) E.derivative
    rw [← AdjoinRoot.aeval_eq, Polynomial.aeval_algHom_apply]

/-- A root with invertible derivative extends uniquely to the localization. -/
def toLocalizedAlgHom
    (s : LocalizedPolynomialRoot E A) : LocalizedAdjoinRoot E →ₐ[K] A :=
  IsLocalization.Away.liftAlgHom (derivativeClass E)
    (f := s.toPolynomialRoot.liftAlgHom)
    (by
      change IsUnit
        (s.toPolynomialRoot.liftAlgHom (AdjoinRoot.mk E E.derivative))
      rw [PolynomialRoot.liftAlgHom_mk]
      rw [← s.derivativeUnit_eq]
      exact s.derivativeUnit.isUnit)

/-- The universal property of `(K[S]/E)[1/E']`, in functor-of-points form. -/
def localizedAlgHomEquiv (E : K[X])
    (A : Type*) [CommRing A] [Algebra K A] :
    (LocalizedAdjoinRoot E →ₐ[K] A) ≃ LocalizedPolynomialRoot E A where
  toFun := ofLocalizedAlgHom
  invFun := toLocalizedAlgHom
  left_inv := by
    intro φ
    apply Localization.algHom_ext (Submonoid.powers (derivativeClass E))
    apply AdjoinRoot.algHom_ext
    simp [toLocalizedAlgHom, ofLocalizedAlgHom, Algebra.algHom]
  right_inv := by
    intro s
    apply LocalizedPolynomialRoot.ext
    simp [toLocalizedAlgHom, ofLocalizedAlgHom, localizedBaseAlgHom]

/-- The localized quotient universal property is natural under
postcomposition of test-algebra maps. -/
theorem localizedAlgHomEquiv_natural
    (f : A →ₐ[K] B) (φ : LocalizedAdjoinRoot E →ₐ[K] A) :
    (localizedAlgHomEquiv E B) (f.comp φ) =
      ((localizedAlgHomEquiv E A) φ).map f := by
  apply LocalizedPolynomialRoot.ext
  rfl

end LocalizedPolynomialRoot

#print axioms LocalizedPolynomialRoot.localizedAlgHomEquiv
#print axioms LocalizedPolynomialRoot.localizedAlgHomEquiv_natural

end FiniteEtaleKeller
