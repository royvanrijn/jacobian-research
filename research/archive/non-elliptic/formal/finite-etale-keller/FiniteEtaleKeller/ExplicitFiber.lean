/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.ExplicitMap
import FiniteEtaleKeller.FiniteEtaleQuotient
import FiniteEtaleKeller.GeneralGaugeFiberRank
import FiniteEtaleKeller.GeneralGaugeRealization

/-!
# The literal fiber of the explicit optimal quintic map

The paper displays a denominator-free map with Jacobian determinant `-722`
and target `(1, 0, -38)`.  This module completes the concrete certificate:
over every commutative rational algebra, maps out of the Berend--Bilu quintic
quotient are naturally equivalent to literal points of that displayed map
fiber.
-/

noncomputable section

open Polynomial
open MvPolynomial

namespace FiniteEtaleKeller.ExplicitQuintic

variable {A B : Type*}
variable [CommRing A] [Algebra ℚ A]
variable [CommRing B] [Algebra ℚ B]

/-- The Bézout identity already checked for `p5` is a separability
certificate. -/
theorem p5_separable : p5.Separable :=
  (Polynomial.separable_def' p5).2 ⟨bezoutU, bezoutV, p5_bezout⟩

/-- The displayed Berend--Bilu quotient is an étale rational algebra. -/
theorem p5_quotient_etale :
    Algebra.Etale ℚ (AdjoinRoot p5) :=
  adjoinRoot_etale_of_separable p5 p5_separable

/-- The displayed quotient is finite as a rational module. -/
theorem p5_quotient_finite :
    Module.Finite ℚ (AdjoinRoot p5) :=
  adjoinRoot_finite_of_separable p5 p5_separable

/-- The displayed polynomial has degree five. -/
theorem p5_natDegree : p5.natDegree = 5 := by
  rw [p5_expanded]
  compute_degree!

/-- The represented algebra in the explicit example has rank five. -/
theorem p5_quotient_rank :
    Module.finrank ℚ (AdjoinRoot p5) = 5 := by
  rw [adjoinRoot_finrank_eq_natDegree, p5_natDegree]

/-- At the normalized target `(1, 0, -2)`, the general inverse polynomial is
literally the Berend--Bilu quintic. -/
theorem generalGaugeInversePolynomial_explicitTarget :
    generalGaugeInversePolynomial g5 1 0 (-2) = p5 := by
  rw [generalGaugeInversePolynomial]
  rw [generalGaugeSeedPolynomial_one_eq g5]
  · have hcoeff : g5.coeff 1 = -19 := by
      norm_num [g5, Polynomial.coeff_X]
    rw [hcoeff]
    norm_num
    rw [← Polynomial.C_mul]
    norm_num
    exact inversePolynomial_eq_p5
  · norm_num [g5, Polynomial.coeff_X]

private theorem g5_coeff_one_ne_zero : g5.coeff 1 ≠ 0 := by
  norm_num [g5, Polynomial.coeff_X]

private theorem g5_coeff_three_ne_zero : g5.coeff 3 ≠ 0 := by
  norm_num [g5, Polynomial.coeff_X]

/-- Separability in the exact syntactic form consumed by the general literal
fiber theorem. -/
theorem generalGaugeInversePolynomial_explicitTarget_separable :
    (generalGaugeInversePolynomial g5 ((1 : ℚˣ) : ℚ) 0 (-2)).Separable := by
  have hInv :
      generalGaugeInversePolynomial g5 ((1 : ℚˣ) : ℚ) 0 (-2) = p5 := by
    simpa using generalGaugeInversePolynomial_explicitTarget
  rw [hInv]
  exact p5_separable

/-- A root of the displayed quintic is exactly a root of the generic inverse
equation at the explicit target.  Transporting roots instead of quotient types
keeps this closed computation lightweight. -/
def p5RootEquivExplicitInverseRoot :
    PolynomialRoot p5 A ≃
      PolynomialRoot
        (generalGaugeInversePolynomial g5 ((1 : ℚˣ) : ℚ) 0 (-2)) A where
  toFun := fun s =>
    ⟨s.1, by
        have hInv :
            generalGaugeInversePolynomial g5 ((1 : ℚˣ) : ℚ) 0 (-2) = p5 := by
          simpa using generalGaugeInversePolynomial_explicitTarget
        rw [hInv]
        exact s.2⟩
  invFun := fun s =>
    ⟨s.1, by
        have hInv :
            generalGaugeInversePolynomial g5 ((1 : ℚˣ) : ℚ) 0 (-2) = p5 := by
          simpa using generalGaugeInversePolynomial_explicitTarget
        rw [← hInv]
        exact s.2⟩
  left_inv := by
    intro s
    apply PolynomialRoot.ext
    rfl
  right_inv := by
    intro s
    apply PolynomialRoot.ext
    rfl

/-- The concrete root identification commutes with morphisms of test
algebras. -/
theorem p5RootEquivExplicitInverseRoot_natural
    (f : A →ₐ[ℚ] B) (s : PolynomialRoot p5 A) :
    (p5RootEquivExplicitInverseRoot s).map f =
      p5RootEquivExplicitInverseRoot (s.map f) := by
  apply PolynomialRoot.ext
  rfl

/-- A literal point of the displayed denominator-free map fiber at the target
`(1, 0, -38)`. -/
@[ext]
structure IntegralFiberPoint
    (A : Type*) [CommRing A] [Algebra ℚ A] where
  point : Fin 3 → A
  first_eq :
    MvPolynomial.eval₂ (algebraMap ℚ A) point (integralMap 0) = 1
  second_eq :
    MvPolynomial.eval₂ (algebraMap ℚ A) point (integralMap 1) = 0
  third_eq :
    MvPolynomial.eval₂ (algebraMap ℚ A) point (integralMap 2) =
      algebraMap ℚ A (-38)

@[simp]
theorem eval₂_integralMap_zero (point : Fin 3 → A) :
    MvPolynomial.eval₂ (algebraMap ℚ A) point (integralMap 0) =
      MvPolynomial.eval₂ (algebraMap ℚ A) point (generalGaugePi g5) := by
  rw [integralMap_eq_scaled_normalized, normalizedMap_eq_generalGaugeMap]
  simp [scaleOutput, generalGaugeMap]

@[simp]
theorem eval₂_integralMap_one (point : Fin 3 → A) :
    MvPolynomial.eval₂ (algebraMap ℚ A) point (integralMap 1) =
      algebraMap ℚ A 19 *
        MvPolynomial.eval₂ (algebraMap ℚ A) point (generalGaugeB g5) := by
  rw [integralMap_eq_scaled_normalized, normalizedMap_eq_generalGaugeMap]
  simp [scaleOutput, generalGaugeMap]

@[simp]
theorem eval₂_integralMap_two (point : Fin 3 → A) :
    MvPolynomial.eval₂ (algebraMap ℚ A) point (integralMap 2) =
      algebraMap ℚ A 19 *
        MvPolynomial.eval₂ (algebraMap ℚ A) point (generalGaugeC g5) := by
  rw [integralMap_eq_scaled_normalized, normalizedMap_eq_generalGaugeMap]
  simp [scaleOutput, generalGaugeMap]

namespace GeneralGaugeRawFiberPoint

/-- Scale a point of the normalized `(1,0,-2)` fiber to the displayed
`(1,0,-38)` fiber.  The source point itself is unchanged. -/
def toIntegral
    (p : GeneralGaugeRawFiberPoint g5 1 0 (-2) A) :
    IntegralFiberPoint A where
  point := p.point
  first_eq := by
    simpa using p.pi_eq
  second_eq := by
    rw [eval₂_integralMap_one, p.b_eq]
    simp
  third_eq := by
    rw [eval₂_integralMap_two, p.c_eq, ← map_mul]
    norm_num

end GeneralGaugeRawFiberPoint

namespace IntegralFiberPoint

/-- Undo the output scaling on the displayed fiber. -/
def toRaw (p : IntegralFiberPoint A) :
    GeneralGaugeRawFiberPoint g5 1 0 (-2) A where
  point := p.point
  pi_eq := by
    simpa using p.first_eq
  b_eq := by
    have hb :
        algebraMap ℚ A 19 *
          MvPolynomial.eval₂ (algebraMap ℚ A) p.point (generalGaugeB g5) = 0 := by
      simpa only [eval₂_integralMap_one] using p.second_eq
    have hscale : algebraMap ℚ A (1 / 19) * algebraMap ℚ A 19 = 1 := by
      rw [← map_mul]
      norm_num
    calc
      MvPolynomial.eval₂ (algebraMap ℚ A) p.point (generalGaugeB g5) =
          1 * MvPolynomial.eval₂ (algebraMap ℚ A) p.point
            (generalGaugeB g5) := by rw [one_mul]
      _ = (algebraMap ℚ A (1 / 19) * algebraMap ℚ A 19) *
          MvPolynomial.eval₂ (algebraMap ℚ A) p.point
            (generalGaugeB g5) := by rw [hscale]
      _ = algebraMap ℚ A (1 / 19) *
          (algebraMap ℚ A 19 *
            MvPolynomial.eval₂ (algebraMap ℚ A) p.point
              (generalGaugeB g5)) := by rw [mul_assoc]
      _ = 0 := by rw [hb, mul_zero]
      _ = algebraMap ℚ A 0 := by simp
  c_eq := by
    have hc :
        algebraMap ℚ A 19 *
          MvPolynomial.eval₂ (algebraMap ℚ A) p.point (generalGaugeC g5) =
            algebraMap ℚ A (-38) := by
      simpa only [eval₂_integralMap_two] using p.third_eq
    have hscale : algebraMap ℚ A (1 / 19) * algebraMap ℚ A 19 = 1 := by
      rw [← map_mul]
      norm_num
    calc
      MvPolynomial.eval₂ (algebraMap ℚ A) p.point (generalGaugeC g5) =
          1 * MvPolynomial.eval₂ (algebraMap ℚ A) p.point
            (generalGaugeC g5) := by rw [one_mul]
      _ = (algebraMap ℚ A (1 / 19) * algebraMap ℚ A 19) *
          MvPolynomial.eval₂ (algebraMap ℚ A) p.point
            (generalGaugeC g5) := by rw [hscale]
      _ = algebraMap ℚ A (1 / 19) *
          (algebraMap ℚ A 19 *
            MvPolynomial.eval₂ (algebraMap ℚ A) p.point
              (generalGaugeC g5)) := by rw [mul_assoc]
      _ = algebraMap ℚ A (1 / 19) * algebraMap ℚ A (-38) := by rw [hc]
      _ = algebraMap ℚ A ((1 / 19) * (-38)) := by rw [map_mul]
      _ = algebraMap ℚ A (-2) := by norm_num

private theorem eval₂_map_algHom
    (f : A →ₐ[ℚ] B) (point : Fin 3 → A) (P : M) :
    MvPolynomial.eval₂ (algebraMap ℚ B) (fun i => f (point i)) P =
      f (MvPolynomial.eval₂ (algebraMap ℚ A) point P) := by
  have hcomp : f.toRingHom.comp (algebraMap ℚ A) = algebraMap ℚ B := by
    ext r
    exact f.commutes r
  calc
    MvPolynomial.eval₂ (algebraMap ℚ B) (fun i => f (point i)) P =
        MvPolynomial.eval₂ (f.toRingHom.comp (algebraMap ℚ A))
          (fun i => f (point i)) P := by rw [hcomp]
    _ = f (MvPolynomial.eval₂ (algebraMap ℚ A) point P) :=
      (MvPolynomial.hom_eval₂ P (algebraMap ℚ A) f.toRingHom point).symm

/-- Literal displayed fibers are functorial in the test algebra. -/
def map (f : A →ₐ[ℚ] B) (p : IntegralFiberPoint A) :
    IntegralFiberPoint B where
  point := fun i => f (p.point i)
  first_eq := by
    rw [eval₂_map_algHom f p.point (integralMap 0), p.first_eq, map_one]
  second_eq := by
    rw [eval₂_map_algHom f p.point (integralMap 1), p.second_eq, map_zero]
  third_eq := by
    rw [eval₂_map_algHom f p.point (integralMap 2), p.third_eq]
    exact f.commutes (-38)

end IntegralFiberPoint

/-- Output scaling identifies the normalized quintic fiber with the literal
fiber of the denominator-free map. -/
def rawFiberEquivIntegral :
    GeneralGaugeRawFiberPoint g5 1 0 (-2) A ≃ IntegralFiberPoint A where
  toFun := GeneralGaugeRawFiberPoint.toIntegral
  invFun := IntegralFiberPoint.toRaw
  left_inv := by
    intro p
    apply GeneralGaugeRawFiberPoint.ext
    rfl
  right_inv := by
    intro p
    apply IntegralFiberPoint.ext
    rfl

/-- The output-scaling equivalence commutes with morphisms of test
algebras. -/
theorem rawFiberEquivIntegral_natural
    (f : A →ₐ[ℚ] B) (p : GeneralGaugeRawFiberPoint g5 1 0 (-2) A) :
    IntegralFiberPoint.map f (rawFiberEquivIntegral p) =
      rawFiberEquivIntegral (GeneralGaugeRawFiberPoint.map f p) := by
  apply IntegralFiberPoint.ext
  rfl

/-- Roots of the explicit generic inverse equation reconstruct literal points
of the normalized gauge fiber. -/
def explicitInverseRootEquivRawFiber :
    PolynomialRoot
        (generalGaugeInversePolynomial g5 ((1 : ℚˣ) : ℚ) 0 (-2)) A ≃
      GeneralGaugeRawFiberPoint g5 1 0 (-2) A :=
  (PolynomialRoot.algHomEquiv
      (generalGaugeInversePolynomial g5 ((1 : ℚˣ) : ℚ) 0 (-2)) A).symm.trans
    (generalGaugeRawRepresentingEquiv g5 1 0 (-2)
      g5_coeff_one_ne_zero g5_coeff_three_ne_zero
      generalGaugeInversePolynomial_explicitTarget_separable A)

/-- Root reconstruction for the explicit inverse equation is natural in the
test algebra. -/
theorem explicitInverseRootEquivRawFiber_natural
    (f : A →ₐ[ℚ] B)
    (s : PolynomialRoot
      (generalGaugeInversePolynomial g5 ((1 : ℚˣ) : ℚ) 0 (-2)) A) :
    GeneralGaugeRawFiberPoint.map f (explicitInverseRootEquivRawFiber s) =
      explicitInverseRootEquivRawFiber (s.map f) := by
  change
    GeneralGaugeRawFiberPoint.map f
        (generalGaugeRawRepresentingEquiv g5 1 0 (-2)
          g5_coeff_one_ne_zero g5_coeff_three_ne_zero
          generalGaugeInversePolynomial_explicitTarget_separable A
          (PolynomialRoot.liftAlgHom s)) =
      generalGaugeRawRepresentingEquiv g5 1 0 (-2)
        g5_coeff_one_ne_zero g5_coeff_three_ne_zero
        generalGaugeInversePolynomial_explicitTarget_separable B
        (PolynomialRoot.liftAlgHom (s.map f))
  rw [generalGaugeRawRepresentingEquiv_natural]
  congr 1
  apply AdjoinRoot.algHom_ext
  rw [AlgHom.comp_apply, PolynomialRoot.liftAlgHom_root,
    PolynomialRoot.liftAlgHom_root, PolynomialRoot.map_val]

/-- The explicit quotient represents the literal displayed map fiber over
every commutative rational algebra. -/
def integralFiberRepresentingEquiv :
    (AdjoinRoot p5 →ₐ[ℚ] A) ≃ IntegralFiberPoint A :=
  (PolynomialRoot.algHomEquiv p5 A).trans
    (p5RootEquivExplicitInverseRoot.trans
      (explicitInverseRootEquivRawFiber.trans rawFiberEquivIntegral))

/-- The concrete quotient/fiber equivalence is natural in every commutative
rational test algebra. -/
theorem integralFiberRepresentingEquiv_natural
    (f : A →ₐ[ℚ] B) (φ : AdjoinRoot p5 →ₐ[ℚ] A) :
    IntegralFiberPoint.map f (integralFiberRepresentingEquiv φ) =
      integralFiberRepresentingEquiv (f.comp φ) := by
  change
    IntegralFiberPoint.map f
        (rawFiberEquivIntegral
          (explicitInverseRootEquivRawFiber
            (p5RootEquivExplicitInverseRoot
              (PolynomialRoot.ofAlgHom φ)))) =
      rawFiberEquivIntegral
        (explicitInverseRootEquivRawFiber
          (p5RootEquivExplicitInverseRoot
            (PolynomialRoot.ofAlgHom (f.comp φ))))
  rw [rawFiberEquivIntegral_natural]
  rw [explicitInverseRootEquivRawFiber_natural]
  rw [p5RootEquivExplicitInverseRoot_natural]
  congr 2

/-- The explicit quintic quotient admits no rational point. -/
theorem p5_adjoinRoot_no_rat_algHom :
    IsEmpty (AdjoinRoot p5 →ₐ[ℚ] ℚ) := by
  constructor
  intro φ
  let root : PolynomialRoot p5 ℚ :=
    PolynomialRoot.algHomEquiv p5 ℚ φ
  exact p5_no_rational_root root.1 root.2

/-- The cubic factor changes sign between `2` and `3`, so the quintic has a
real root. -/
theorem p5_has_real_root : Nonempty (PolynomialRoot p5 ℝ) := by
  have hcont : Continuous (fun x : ℝ => x ^ 3 - 19) := by
    fun_prop
  have hzero :
      (0 : ℝ) ∈ Set.Icc ((2 : ℝ) ^ 3 - 19) ((3 : ℝ) ^ 3 - 19) := by
    norm_num
  obtain ⟨r, -, hr⟩ := (Set.mem_image ..).mp
    (intermediate_value_Icc (by norm_num : (2 : ℝ) ≤ 3)
      hcont.continuousOn hzero)
  refine ⟨⟨r, ?_⟩⟩
  simp [p5, Polynomial.aeval_def, hr]

/-- Consequently the literal fiber of the displayed map over
`(1, 0, -38)` has no rational point. -/
theorem integralFiberPoint_rat_isEmpty :
    IsEmpty (IntegralFiberPoint ℚ) := by
  constructor
  intro point
  exact p5_adjoinRoot_no_rat_algHom.false
    (integralFiberRepresentingEquiv.symm point)

/-- The same literal fiber has a real point. -/
theorem integralFiberPoint_real_nonempty :
    Nonempty (IntegralFiberPoint ℝ) := by
  obtain ⟨root⟩ := p5_has_real_root
  exact ⟨integralFiberRepresentingEquiv
    ((PolynomialRoot.algHomEquiv p5 ℝ).symm root)⟩

#print axioms p5_separable
#print axioms p5_quotient_etale
#print axioms p5_quotient_finite
#print axioms p5_quotient_rank
#print axioms generalGaugeInversePolynomial_explicitTarget
#print axioms integralFiberRepresentingEquiv
#print axioms integralFiberRepresentingEquiv_natural
#print axioms integralFiberPoint_rat_isEmpty
#print axioms integralFiberPoint_real_nonempty

end FiniteEtaleKeller.ExplicitQuintic
