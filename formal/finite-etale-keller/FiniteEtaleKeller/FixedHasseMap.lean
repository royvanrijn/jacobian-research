/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.FixedHasseFamily

/-!
# The exact Jacobian-one normalization of the fixed Hasse map

The fixed-map paper does not use the universal output normalization: it
halves the first source coordinate and negates the first target coordinate.
This keeps the moving target

`(-1, 32a/9, (8a+1)/3)`

at its smallest natural height.  This module formalizes that exact affine
normalization and proves that its Jacobian determinant is one.
-/

noncomputable section

open Matrix Function
open MvPolynomial

namespace FiniteEtaleKeller

variable {R : Type*} [CommRing R]

/-- Diagonal substitution in the three source variables. -/
def diagonalSourceSubstitution (a b c : R) :
    Fin 3 → MvPolynomial (Fin 3) R :=
  ![C a * X 0, C b * X 1, C c * X 2]

/-- Precompose a three-variable polynomial map with a diagonal source
scaling. -/
def scaleInput (a b c : R)
    (F : Fin 3 → MvPolynomial (Fin 3) R) :
    Fin 3 → MvPolynomial (Fin 3) R :=
  fun i => bind₁ (diagonalSourceSubstitution a b c) (F i)

private theorem pderiv_diagonalSourceSubstitution
    (a b c : R) (i j : Fin 3) :
    pderiv j (diagonalSourceSubstitution a b c i) =
      if i = j then C (![a, b, c] j) else 0 := by
  fin_cases i <;> fin_cases j <;>
    simp [diagonalSourceSubstitution]

private theorem pderiv_bind_diagonal
    (a b c : R) (P : MvPolynomial (Fin 3) R) (j : Fin 3) :
    pderiv j
        (bind₁ (diagonalSourceSubstitution a b c) P) =
      C (![a, b, c] j) *
        bind₁ (diagonalSourceSubstitution a b c) (pderiv j P) := by
  classical
  induction P using MvPolynomial.induction_on with
  | C r =>
      simp [diagonalSourceSubstitution]
  | add P Q hP hQ =>
      simp only [map_add, hP, hQ]
      ring
  | mul_X P i hP =>
      simp only [map_mul, bind₁_X_right, pderiv_mul, hP,
        pderiv_X, Pi.single_apply,
        pderiv_diagonalSourceSubstitution]
      by_cases hij : i = j
      · subst i
        simp
        ring
      · simp [hij]
        ring

/-- Diagonal source scaling multiplies the Jacobian by the product of the
three scaling factors and substitutes the scaled variables into the old
Jacobian. -/
theorem jacobianDet_scaleInput (a b c : R)
    (F : Fin 3 → MvPolynomial (Fin 3) R) :
    jacobianDet (scaleInput a b c F) =
      C (a * b * c) *
        bind₁ (diagonalSourceSubstitution a b c) (jacobianDet F) := by
  simp only [jacobianDet, jacobianMatrix, det_fin_three, of_apply,
    scaleInput, pderiv_bind_diagonal, map_add, map_sub, map_mul]
  simp only [diagonalSourceSubstitution, cons_val_zero, cons_val_one,
    cons_val_two, head_cons, tail_cons]
  ring

section Evaluation

variable {A : Type*} [CommRing A] [Algebra R A]

/-- Evaluation after source scaling is evaluation of the original map at the
scaled point. -/
theorem eval₂_scaleInput (a b c : R)
    (F : Fin 3 → MvPolynomial (Fin 3) R)
    (point : Fin 3 → A) (i : Fin 3) :
    MvPolynomial.eval₂ (algebraMap R A) point
        (scaleInput a b c F i) =
      MvPolynomial.eval₂ (algebraMap R A)
        (fun j => algebraMap R A (![a, b, c] j) * point j)
        (F i) := by
  change
    eval₂Hom (algebraMap R A) point
        (bind₁ (diagonalSourceSubstitution a b c) (F i)) =
      eval₂Hom (algebraMap R A)
        (fun j => algebraMap R A (![a, b, c] j) * point j) (F i)
  rw [eval₂Hom_bind₁]
  congr 2
  funext j
  fin_cases j <;> simp [diagonalSourceSubstitution]

end Evaluation

/-- Successive diagonal source substitutions multiply their scaling
factors.  This is the algebraic identity which makes every nonzero
diagonal source scaling a polynomial automorphism. -/
theorem scaleInput_scaleInput (a b c d e f : R)
    (F : Fin 3 → MvPolynomial (Fin 3) R) :
    scaleInput a b c (scaleInput d e f F) =
      scaleInput (d * a) (e * b) (f * c) F := by
  funext i
  simp only [scaleInput, MvPolynomial.bind₁_bind₁]
  congr 1
  apply MvPolynomial.algHom_ext
  intro j
  fin_cases j <;> simp [diagonalSourceSubstitution, mul_assoc]

/-- Scaling every source coordinate by one is the identity. -/
@[simp]
theorem scaleInput_one
    (F : Fin 3 → MvPolynomial (Fin 3) R) :
    scaleInput 1 1 1 F = F := by
  funext i
  change
    MvPolynomial.bind₁ (diagonalSourceSubstitution 1 1 1) (F i) = F i
  have h :
      diagonalSourceSubstitution (1 : R) 1 1 =
        (MvPolynomial.X : Fin 3 → MvPolynomial (Fin 3) R) := by
    funext j
    fin_cases j <;> simp [diagonalSourceSubstitution]
  rw [h]
  exact DFunLike.congr_fun MvPolynomial.bind₁_X_left (F i)

/-- Source and output diagonal scalings commute. -/
theorem scaleInput_scaleOutput (a b c d e f : R)
    (F : Fin 3 → MvPolynomial (Fin 3) R) :
    scaleInput a b c (scaleOutput d e f F) =
      scaleOutput d e f (scaleInput a b c F) := by
  funext i
  fin_cases i <;> simp [scaleInput, scaleOutput]

namespace FixedHasseFamily

/-- The exact fixed map `Φ = L ∘ F₀ ∘ A` from the paper, where
`A(x,y,z)=(x/2,y,z)` and `L(Π,B,C)=(-Π,B,C)`. -/
def paperMap : Fin 3 → GaugePolynomial ℚ :=
  scaleOutput (-1) 1 1 (scaleInput (1 / 2) 1 1 baseMap)

/-- The paper's fixed map has Jacobian determinant one. -/
theorem jacobianDet_paperMap :
    jacobianDet paperMap = 1 := by
  rw [paperMap, jacobianDet_scaleOutput, jacobianDet_scaleInput,
    jacobianDet_baseMap]
  simp [diagonalSourceSubstitution]
  norm_num [← MvPolynomial.C_mul]

/-- Applying the displayed inverse source and target scalings recovers the
base map exactly.  Thus `paperMap` and `baseMap` are related by polynomial
automorphisms on both sides, rather than merely having matching Jacobians. -/
theorem paperMap_normalization_inverse :
    scaleInput 2 1 1 (scaleOutput (-1) 1 1 paperMap) = baseMap := by
  rw [paperMap, scaleInput_scaleOutput, scaleInput_scaleOutput,
    scaleInput_scaleInput]
  simp only [one_div, ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true,
    inv_mul_cancel₀, mul_one, scaleInput_one]
  funext i
  fin_cases i <;> simp [scaleOutput]

/-- The geometric degree of the exact displayed paper map, transported
through the certified source and target polynomial automorphisms to the
canonical function-field presentation of `baseMap`. -/
def paperMapGeometricDegree : ℕ :=
  generalGaugeGeometricDegree seed
    seed_coeff_one_ne_zero seed_coeff_three_ne_zero

/-- The exact displayed Jacobian-one paper map has geometric degree five. -/
theorem paperMap_geometricDegree :
    paperMapGeometricDegree = 5 :=
  baseMap_geometricDegree

variable {A : Type*} [CommRing A] [Algebra ℚ A]

/-- The source point after applying `A(x,y,z)=(x/2,y,z)`. -/
def halfFirstPoint (point : Fin 3 → A) : Fin 3 → A :=
  ![algebraMap ℚ A (1 / 2) * point 0, point 1, point 2]

/-- Evaluation of the first paper-map coordinate. -/
@[simp]
theorem eval₂_paperMap_zero (point : Fin 3 → A) :
    MvPolynomial.eval₂ (algebraMap ℚ A) point (paperMap 0) =
      -MvPolynomial.eval₂ (algebraMap ℚ A) (halfFirstPoint point)
        (baseMap 0) := by
  change
    MvPolynomial.eval₂ (algebraMap ℚ A) point
        (C (-1) * scaleInput (1 / 2) 1 1 baseMap 0) =
      -MvPolynomial.eval₂ (algebraMap ℚ A) (halfFirstPoint point)
        (baseMap 0)
  rw [eval₂_mul, eval₂_C]
  rw [eval₂_scaleInput]
  have hp :
      (fun j => algebraMap ℚ A (![(1 / 2 : ℚ), 1, 1] j) * point j) =
        halfFirstPoint point := by
    funext j
    fin_cases j <;> simp [halfFirstPoint]
  rw [hp]
  simp

/-- Evaluation of the second paper-map coordinate. -/
@[simp]
theorem eval₂_paperMap_one (point : Fin 3 → A) :
    MvPolynomial.eval₂ (algebraMap ℚ A) point (paperMap 1) =
      MvPolynomial.eval₂ (algebraMap ℚ A) (halfFirstPoint point)
        (baseMap 1) := by
  change
    MvPolynomial.eval₂ (algebraMap ℚ A) point
        (C 1 * scaleInput (1 / 2) 1 1 baseMap 1) =
      MvPolynomial.eval₂ (algebraMap ℚ A) (halfFirstPoint point)
        (baseMap 1)
  rw [eval₂_mul, eval₂_C]
  rw [eval₂_scaleInput]
  have hp :
      (fun j => algebraMap ℚ A (![(1 / 2 : ℚ), 1, 1] j) * point j) =
        halfFirstPoint point := by
    funext j
    fin_cases j <;> simp [halfFirstPoint]
  rw [hp]
  simp

/-- Evaluation of the third paper-map coordinate. -/
@[simp]
theorem eval₂_paperMap_two (point : Fin 3 → A) :
    MvPolynomial.eval₂ (algebraMap ℚ A) point (paperMap 2) =
      MvPolynomial.eval₂ (algebraMap ℚ A) (halfFirstPoint point)
        (baseMap 2) := by
  change
    MvPolynomial.eval₂ (algebraMap ℚ A) point
        (C 1 * scaleInput (1 / 2) 1 1 baseMap 2) =
      MvPolynomial.eval₂ (algebraMap ℚ A) (halfFirstPoint point)
        (baseMap 2)
  rw [eval₂_mul, eval₂_C]
  rw [eval₂_scaleInput]
  have hp :
      (fun j => algebraMap ℚ A (![(1 / 2 : ℚ), 1, 1] j) * point j) =
        halfFirstPoint point := by
    funext j
    fin_cases j <;> simp [halfFirstPoint]
  rw [hp]
  simp

/-- A literal point of the paper's Jacobian-one map over its moving target
`(-1, 32a/9, (8a+1)/3)`. -/
@[ext]
structure PaperFiberPoint (a : ℚ)
    (A : Type*) [CommRing A] [Algebra ℚ A] where
  point : Fin 3 → A
  first_eq :
    MvPolynomial.eval₂ (algebraMap ℚ A) point (paperMap 0) =
      algebraMap ℚ A (-1)
  second_eq :
    MvPolynomial.eval₂ (algebraMap ℚ A) point (paperMap 1) =
      algebraMap ℚ A (targetB a)
  third_eq :
    MvPolynomial.eval₂ (algebraMap ℚ A) point (paperMap 2) =
      algebraMap ℚ A (targetC a)

/-- The inverse source scaling `A⁻¹(x,y,z)=(2x,y,z)`. -/
def doubleFirstPoint (point : Fin 3 → A) : Fin 3 → A :=
  ![algebraMap ℚ A 2 * point 0, point 1, point 2]

@[simp]
theorem halfFirstPoint_doubleFirstPoint (point : Fin 3 → A) :
    halfFirstPoint (doubleFirstPoint point) = point := by
  funext i
  fin_cases i
  · change
      algebraMap ℚ A (1 / 2) * (algebraMap ℚ A 2 * point 0) =
        point 0
    rw [← mul_assoc, ← map_mul]
    norm_num
  · rfl
  · rfl

@[simp]
theorem doubleFirstPoint_halfFirstPoint (point : Fin 3 → A) :
    doubleFirstPoint (halfFirstPoint point) = point := by
  funext i
  fin_cases i
  · change
      algebraMap ℚ A 2 * (algebraMap ℚ A (1 / 2) * point 0) =
        point 0
    rw [← mul_assoc, ← map_mul]
    norm_num
  · rfl
  · rfl

namespace PaperFiberPoint

variable {a : ℚ}

/-- A paper-map point gives a point of the determinant-`-2` base fiber after
applying the source scaling. -/
def toRaw (p : PaperFiberPoint a A) :
    GeneralGaugeRawFiberPoint seed 1 (targetB a) (targetC a) A where
  point := halfFirstPoint p.point
  pi_eq := by
    have h := p.first_eq
    rw [eval₂_paperMap_zero] at h
    have hneg := congrArg Neg.neg h
    change
      MvPolynomial.eval₂ (algebraMap ℚ A)
          (halfFirstPoint p.point) (baseMap 0) =
        algebraMap ℚ A (1 : ℚ)
    simpa only [map_neg, map_one, neg_neg] using hneg
  b_eq := by
    have h := p.second_eq
    rw [eval₂_paperMap_one] at h
    change
      MvPolynomial.eval₂ (algebraMap ℚ A)
          (halfFirstPoint p.point) (baseMap 1) =
        algebraMap ℚ A (targetB a)
    exact h
  c_eq := by
    have h := p.third_eq
    rw [eval₂_paperMap_two] at h
    change
      MvPolynomial.eval₂ (algebraMap ℚ A)
          (halfFirstPoint p.point) (baseMap 2) =
        algebraMap ℚ A (targetC a)
    exact h

end PaperFiberPoint

namespace GeneralGaugeRawFiberPoint

variable {a : ℚ}

/-- Undo the source scaling to obtain a literal point of the paper map. -/
def toPaper
    (p : GeneralGaugeRawFiberPoint seed 1 (targetB a) (targetC a) A) :
    PaperFiberPoint a A where
  point := doubleFirstPoint p.point
  first_eq := by
    rw [eval₂_paperMap_zero, halfFirstPoint_doubleFirstPoint]
    change
      -MvPolynomial.eval₂ (algebraMap ℚ A) p.point
          (generalGaugePi seed) = algebraMap ℚ A (-1)
    rw [p.pi_eq]
    simp
  second_eq := by
    rw [eval₂_paperMap_one, halfFirstPoint_doubleFirstPoint]
    change
      MvPolynomial.eval₂ (algebraMap ℚ A) p.point
          (generalGaugeB seed) = algebraMap ℚ A (targetB a)
    exact p.b_eq
  third_eq := by
    rw [eval₂_paperMap_two, halfFirstPoint_doubleFirstPoint]
    change
      MvPolynomial.eval₂ (algebraMap ℚ A) p.point
          (generalGaugeC seed) = algebraMap ℚ A (targetC a)
    exact p.c_eq

end GeneralGaugeRawFiberPoint

/-- The affine source and target normalizations identify the complete base
fiber with the literal paper-map fiber. -/
def rawFiberEquivPaper (a : ℚ) :
    GeneralGaugeRawFiberPoint seed 1 (targetB a) (targetC a) A ≃
      PaperFiberPoint a A where
  toFun := GeneralGaugeRawFiberPoint.toPaper
  invFun := PaperFiberPoint.toRaw
  left_inv := by
    intro p
    apply GeneralGaugeRawFiberPoint.ext
    exact halfFirstPoint_doubleFirstPoint p.point
  right_inv := by
    intro p
    apply PaperFiberPoint.ext
    exact doubleFirstPoint_halfFirstPoint p.point

namespace PaperFiberPoint

variable {a : ℚ}
variable {B : Type*} [CommRing B] [Algebra ℚ B]

/-- Literal fibers of the exact paper map are functorial in the commutative
test algebra. -/
def map (f : A →ₐ[ℚ] B) (p : PaperFiberPoint a A) :
    PaperFiberPoint a B :=
  rawFiberEquivPaper (A := B) a
    (GeneralGaugeRawFiberPoint.map f
      ((rawFiberEquivPaper (A := A) a).symm p))

end PaperFiberPoint

/-- The centered inverse quotient naturally represents the literal fiber of
the exact Jacobian-one paper map. -/
def centeredPaperFiberRepresentingEquiv
    (a : ℚ) (ha0 : a ≠ 0) (ha1 : a ≠ 1) :
    (AdjoinRoot
        (generalGaugeInversePolynomial
          seed 1 (targetB a) (targetC a)) →ₐ[ℚ] A) ≃
      PaperFiberPoint a A :=
  (centeredRawFiberRepresentingEquiv (A := A) a ha0 ha1).trans
    (rawFiberEquivPaper (A := A) a)

/-- Precomposition with the translation equivalence identifies maps from the
uncentered arithmetic quotient with maps from the centered inverse
quotient. -/
def inverseQuotientHomEquiv (a : ℚ) :
    (AdjoinRoot (polynomial a) →ₐ[ℚ] A) ≃
      (AdjoinRoot
        (generalGaugeInversePolynomial
          seed 1 (targetB a) (targetC a)) →ₐ[ℚ] A) := by
  rw [inversePolynomial_eq_translate]
  exact translationQuotientHomEquiv (A := A) (polynomial a) (-1 / 2)

/-- The uncentered arithmetic quotient naturally represents the literal
fiber of the exact Jacobian-one paper map. -/
def paperFiberRepresentingEquiv
    (a : ℚ) (ha0 : a ≠ 0) (ha1 : a ≠ 1) :
    (AdjoinRoot (polynomial a) →ₐ[ℚ] A) ≃ PaperFiberPoint a A :=
  (inverseQuotientHomEquiv (A := A) a).trans
    (centeredPaperFiberRepresentingEquiv (A := A) a ha0 ha1)

/-- If `a` is not a rational cube, the literal paper fiber has no rational
point. -/
theorem paperFiberPoint_rat_isEmpty
    (a : ℚ) (ha0 : a ≠ 0) (ha1 : a ≠ 1)
    (hnoncube : ¬∃ r : ℚ, r ^ 3 = a) :
    IsEmpty (PaperFiberPoint a ℚ) := by
  constructor
  intro point
  let φ : AdjoinRoot (polynomial a) →ₐ[ℚ] ℚ :=
    (paperFiberRepresentingEquiv
      (A := ℚ) a ha0 ha1).symm point
  let root : PolynomialRoot (polynomial a) ℚ :=
    PolynomialRoot.algHomEquiv (polynomial a) ℚ φ
  exact polynomial_no_rational_root a hnoncube root.1 root.2

/-- Every literal paper fiber in the family has a real point. -/
theorem paperFiberPoint_real_nonempty
    (a : ℚ) (ha0 : a ≠ 0) (ha1 : a ≠ 1) :
    Nonempty (PaperFiberPoint a ℝ) := by
  obtain ⟨root⟩ := polynomial_has_real_root a
  exact ⟨paperFiberRepresentingEquiv (A := ℝ) a ha0 ha1
    ((PolynomialRoot.algHomEquiv (polynomial a) ℝ).symm root)⟩

#print axioms jacobianDet_scaleInput
#print axioms jacobianDet_paperMap
#print axioms paperMap_normalization_inverse
#print axioms paperMap_geometricDegree
#print axioms centeredPaperFiberRepresentingEquiv
#print axioms paperFiberRepresentingEquiv
#print axioms paperFiberPoint_rat_isEmpty
#print axioms paperFiberPoint_real_nonempty

end FixedHasseFamily
end FiniteEtaleKeller
