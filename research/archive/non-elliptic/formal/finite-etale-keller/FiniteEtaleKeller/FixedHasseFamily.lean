/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeFunctionFieldComparison
import FiniteEtaleKeller.GeneralGaugeRawFiber
import FiniteEtaleKeller.GeneralGaugeRealization
import FiniteEtaleKeller.GeneralGaugeFiberRank
import FiniteEtaleKeller.TranslationQuotient

/-!
# The fixed-map Hasse family

This module formalizes the algebraic spine of the fixed-map paper.  A single
quintic quadratic-gauge seed is fixed once and for all, while the parameter
`a` enters only through the target.  The corresponding inverse polynomial is
proved to be the translate of

`(X^3 - a) * (X^2 + X + 1)`.

For `a ≠ 0, 1`, explicit Bézout identities prove separability, so the
quotient is finite étale of rank five and naturally represents the complete
literal fiber of the fixed determinant-`-2` map.  The elementary rational and
real point assertions are also recorded.

The analytic counting theorem and the uniform local-solubility theorem for
the paper's restricted integer parameter set are deliberately separate from
this algebraic module.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller.FixedHasseFamily

/-- The fixed seed left after centering the family at `S = X + 1/2`. -/
def seed : ℚ[X] :=
  X ^ 5 - C (3 / 2 : ℚ) * X ^ 4 + C (3 / 2 : ℚ) * X ^ 3
    - C (5 / 4 : ℚ) * X ^ 2 + C (9 / 16 : ℚ) * X

/-- The variable cubic factor. -/
def cubic (a : ℚ) : ℚ[X] :=
  X ^ 3 - C a

/-- The fixed cyclotomic quadratic factor. -/
def quadratic : ℚ[X] :=
  X ^ 2 + X + 1

/-- The uncentered arithmetic quintic. -/
def polynomial (a : ℚ) : ℚ[X] :=
  cubic a * quadratic

/-- The same quintic in the marked inverse coordinate `S = X + 1/2`. -/
def centeredPolynomial (a : ℚ) : ℚ[X] :=
  ((X - C (1 / 2 : ℚ)) ^ 3 - C a) *
    (X ^ 2 + C (3 / 4 : ℚ))

/-- The moving second target coordinate of the determinant-`-2` gauge. -/
def targetB (a : ℚ) : ℚ :=
  32 * a / 9

/-- The moving third target coordinate of the determinant-`-2` gauge. -/
def targetC (a : ℚ) : ℚ :=
  (8 * a + 1) / 3

/-- The fixed determinant-`-2` polynomial map before the paper's affine
Jacobian-one normalization. -/
def baseMap : Fin 3 → GaugePolynomial ℚ :=
  generalGaugeMap seed

@[simp]
theorem seed_coeff_zero : seed.coeff 0 = 0 := by
  norm_num [seed, Polynomial.coeff_X]

@[simp]
theorem seed_coeff_one : seed.coeff 1 = 9 / 16 := by
  norm_num [seed, Polynomial.coeff_X]

@[simp]
theorem seed_coeff_three : seed.coeff 3 = 3 / 2 := by
  norm_num [seed, Polynomial.coeff_X]

theorem seed_coeff_one_ne_zero : seed.coeff 1 ≠ 0 := by
  norm_num

theorem seed_coeff_three_ne_zero : seed.coeff 3 ≠ 0 := by
  norm_num

@[simp]
theorem seed_natDegree : seed.natDegree = 5 := by
  unfold seed
  compute_degree!

/-- Centering by `-1/2` gives the displayed factorization in the inverse
coordinate. -/
theorem translate_polynomial_eq_centered (a : ℚ) :
    translatePolynomial (polynomial a) (-1 / 2) =
      centeredPolynomial a := by
  have hshift :
      X + C (-1 / 2 : ℚ) = X - C (1 / 2 : ℚ) := by
    have h : (-1 / 2 : ℚ) = -(1 / 2 : ℚ) := by norm_num
    rw [h, Polynomial.C_neg]
    ring
  have hconstant :
      C (1 / 2 : ℚ) ^ 2 - C (1 / 2 : ℚ) + 1 =
        C (3 / 4 : ℚ) := by
    have h : (1 / 2 : ℚ) ^ 2 - 1 / 2 + 1 = 3 / 4 := by norm_num
    simpa only [Polynomial.C_pow, Polynomial.C_sub, Polynomial.C_add,
      map_one] using congrArg (C : ℚ → ℚ[X]) h
  have hlinear :
      1 - C (1 / 2 : ℚ) * 2 = (0 : ℚ[X]) := by
    have h : (1 : ℚ) - (1 / 2) * 2 = 0 := by norm_num
    simpa only [Polynomial.C_sub, Polynomial.C_mul, map_one, map_ofNat,
      map_zero] using congrArg (C : ℚ → ℚ[X]) h
  simp only [translatePolynomial, polynomial, cubic, quadratic,
    Polynomial.mul_comp, Polynomial.sub_comp, Polynomial.pow_comp,
    Polynomial.X_comp, Polynomial.C_comp, Polynomial.add_comp,
    Polynomial.one_comp]
  rw [hshift]
  unfold centeredPolynomial
  have hquad :
    (X - C (1 / 2 : ℚ)) ^ 2 + (X - C (1 / 2 : ℚ)) + 1 =
        X ^ 2 + (1 - C (1 / 2 : ℚ) * 2) * X +
          (C (1 / 2 : ℚ) ^ 2 - C (1 / 2 : ℚ) + 1) := by ring
  rw [hquad, hlinear, hconstant]
  ring

/-- The elementary centered identity from which the fixed seed is
discovered. -/
theorem centered_identity :
    (X - C (1 / 2 : ℚ)) ^ 3 * (X ^ 2 + C (3 / 4 : ℚ)) =
      seed - C (3 / 32 : ℚ) := by
  unfold seed
  have hX :
      C (1 / 2 : ℚ) ^ 2 * C (3 / 4 : ℚ) * 3 =
        C (9 / 16 : ℚ) := by
    have h : (1 / 2 : ℚ) ^ 2 * (3 / 4) * 3 = 9 / 16 := by norm_num
    simpa only [Polynomial.C_pow, Polynomial.C_mul, map_ofNat] using
      congrArg (C : ℚ → ℚ[X]) h
  have hX2 :
      C (1 / 2 : ℚ) * C (3 / 4 : ℚ) * 3 +
          C (1 / 2 : ℚ) ^ 3 =
        C (5 / 4 : ℚ) := by
    have h :
        (1 / 2 : ℚ) * (3 / 4) * 3 + (1 / 2) ^ 3 = 5 / 4 := by
      norm_num
    simpa only [Polynomial.C_pow, Polynomial.C_mul, Polynomial.C_add,
      map_ofNat] using congrArg (C : ℚ → ℚ[X]) h
  have hX3 :
      C (1 / 2 : ℚ) ^ 2 * 3 + C (3 / 4 : ℚ) =
        C (3 / 2 : ℚ) := by
    have h : (1 / 2 : ℚ) ^ 2 * 3 + 3 / 4 = 3 / 2 := by norm_num
    simpa only [Polynomial.C_pow, Polynomial.C_mul, Polynomial.C_add,
      map_ofNat] using congrArg (C : ℚ → ℚ[X]) h
  have hX4 :
      C (1 / 2 : ℚ) * 3 = C (3 / 2 : ℚ) := by
    have h : (1 / 2 : ℚ) * 3 = 3 / 2 := by norm_num
    simpa only [Polynomial.C_mul, map_ofNat] using
      congrArg (C : ℚ → ℚ[X]) h
  have hC :
      C (1 / 2 : ℚ) ^ 3 * C (3 / 4 : ℚ) =
        C (3 / 32 : ℚ) := by
    have h : (1 / 2 : ℚ) ^ 3 * (3 / 4) = 3 / 32 := by norm_num
    simpa only [Polynomial.C_pow, Polynomial.C_mul] using
      congrArg (C : ℚ → ℚ[X]) h
  ring_nf
  linear_combination
    X * hX - X ^ 2 * hX2 + X ^ 3 * hX3 - X ^ 4 * hX4 - hC

/-- At the moving target, the inverse equation of the one fixed map is the
centered arithmetic quintic. -/
theorem inversePolynomial_eq_centered (a : ℚ) :
    generalGaugeInversePolynomial seed 1 (targetB a) (targetC a) =
      centeredPolynomial a := by
  rw [generalGaugeInversePolynomial]
  rw [generalGaugeSeedPolynomial_one_eq seed seed_coeff_zero]
  simp only [seed_coeff_one]
  have hB :
      C (9 / 32 : ℚ) * C (targetB a) = C a := by
    rw [← Polynomial.C_mul]
    unfold targetB
    congr 1
    ring
  have hC :
      C (9 / 32 : ℚ) * C (targetC a) =
        C (3 / 32 : ℚ) + C a * C (3 / 4 : ℚ) := by
    rw [← Polynomial.C_mul, ← Polynomial.C_mul, ← Polynomial.C_add]
    unfold targetC
    congr 1
    ring
  have hcenter :
      centeredPolynomial a =
        (X - C (1 / 2 : ℚ)) ^ 3 *
            (X ^ 2 + C (3 / 4 : ℚ)) -
          C a * (X ^ 2 + C (3 / 4 : ℚ)) := by
    unfold centeredPolynomial
    ring
  have hhalf :
      C ((9 / 16 : ℚ) / 2) = C (9 / 32 : ℚ) := by
    congr 1
    norm_num
  rw [hcenter, centered_identity]
  rw [hhalf]
  calc
    seed -
          C (9 / 32 : ℚ) *
            (C (targetB a) * X ^ 2 + C (targetC a)) =
        seed -
          (C (9 / 32 : ℚ) * C (targetB a)) * X ^ 2 -
          C (9 / 32 : ℚ) * C (targetC a) := by ring
    _ = seed - C a * X ^ 2 -
          (C (3 / 32 : ℚ) + C a * C (3 / 4 : ℚ)) := by
      rw [hB, hC]
    _ = seed - C (3 / 32 : ℚ) -
          C a * (X ^ 2 + C (3 / 4 : ℚ)) := by ring

/-- Equivalently, the inverse equation is the translation by `-1/2` of
`(X^3-a)(X^2+X+1)`. -/
theorem inversePolynomial_eq_translate (a : ℚ) :
    generalGaugeInversePolynomial seed 1 (targetB a) (targetC a) =
      translatePolynomial (polynomial a) (-1 / 2) := by
  rw [inversePolynomial_eq_centered, translate_polynomial_eq_centered]

/-- An explicit Bézout certificate for the cubic and its derivative. -/
theorem cubic_separable (a : ℚ) (ha : a ≠ 0) :
    (cubic a).Separable := by
  apply (Polynomial.separable_def' (cubic a)).2
  refine ⟨C (-1 / a), C (1 / (3 * a)) * X, ?_⟩
  have hcoeff :
      C (-1 / a) + C (1 / (3 * a)) * C (3 : ℚ) = 0 := by
    have h : (-1 / a) + (1 / (3 * a)) * 3 = (0 : ℚ) := by
      field_simp [ha]
      ring
    simpa only [Polynomial.C_add, Polynomial.C_mul, map_zero] using
      congrArg (C : ℚ → ℚ[X]) h
  have hconstant :
      -(C (-1 / a) * C a) = (1 : ℚ[X]) := by
    have h : -((-1 / a) * a) = (1 : ℚ) := by
      field_simp [ha]
    simpa only [Polynomial.C_neg, Polynomial.C_mul, map_one] using
      congrArg (C : ℚ → ℚ[X]) h
  simp only [cubic, Polynomial.derivative_sub, Polynomial.derivative_pow,
    Polynomial.derivative_X, Polynomial.derivative_C, mul_one, sub_zero]
  calc
    C (-1 / a) * (X ^ 3 - C a) +
          C (1 / (3 * a)) * X * (C 3 * X ^ 2) =
        (C (-1 / a) + C (1 / (3 * a)) * C 3) * X ^ 3 -
          C (-1 / a) * C a := by ring
    _ = 0 * X ^ 3 + (-(C (-1 / a) * C a)) := by rw [hcoeff]; ring
    _ = 1 := by rw [hconstant]; simp

/-- An explicit Bézout certificate for the cyclotomic quadratic. -/
theorem quadratic_separable : quadratic.Separable := by
  apply (Polynomial.separable_def' quadratic).2
  refine ⟨C (4 / 3 : ℚ), -(C (1 / 3 : ℚ) * (2 * X + 1)), ?_⟩
  have hlinear :
      C (4 / 3 : ℚ) - C (1 / 3 : ℚ) * C 4 = 0 := by
    rw [← Polynomial.C_mul, ← Polynomial.C_sub]
    norm_num
  have hconstant :
      C (4 / 3 : ℚ) - C (1 / 3 : ℚ) = 1 := by
    rw [← Polynomial.C_sub]
    norm_num
  have htwoAdd :
      C (2 : ℚ) + C 2 = C 4 := by
    rw [← Polynomial.C_add]
    norm_num
  have htwoMul :
      C (2 : ℚ) * C 2 = C 4 := by
    rw [← Polynomial.C_mul]
    norm_num
  have htwoCast : (2 : ℚ[X]) = C 2 := by
    exact (Polynomial.C_ofNat (R := ℚ) 2).symm
  simp only [quadratic, Polynomial.derivative_add,
    Polynomial.derivative_pow, Polynomial.derivative_X,
    Polynomial.derivative_one, mul_one, add_zero]
  rw [htwoCast]
  ring_nf
  linear_combination X * hlinear + X ^ 2 * hlinear + hconstant -
    X * C (1 / 3 : ℚ) * htwoAdd -
    X ^ 2 * C (1 / 3 : ℚ) * htwoMul

/-- The two factors are coprime away from the unique collision `a = 1`.
The certificate is the identity

`-(X^3-a) + (X-1)(X^2+X+1) = a-1`.
-/
theorem cubic_isCoprime_quadratic (a : ℚ) (ha : a ≠ 1) :
    IsCoprime (cubic a) quadratic := by
  refine ⟨C (-1 / (a - 1)), C (1 / (a - 1)) * (X - 1), ?_⟩
  have hconstant :
      C (-1 / (a - 1)) * (-C a) -
          C (1 / (a - 1)) = (1 : ℚ[X]) := by
    have h :
        (-1 / (a - 1)) * (-a) - 1 / (a - 1) = (1 : ℚ) := by
      field_simp [sub_ne_zero.mpr ha]
    simpa only [Polynomial.C_neg, Polynomial.C_mul, Polynomial.C_sub,
      map_one] using congrArg (C : ℚ → ℚ[X]) h
  have hcoeff :
      C (-1 / (a - 1)) + C (1 / (a - 1)) = 0 := by
    have h : -1 / (a - 1) + 1 / (a - 1) = (0 : ℚ) := by ring
    simpa only [Polynomial.C_add, map_zero] using
      congrArg (C : ℚ → ℚ[X]) h
  unfold cubic quadratic
  calc
    C (-1 / (a - 1)) * (X ^ 3 - C a) +
          (C (1 / (a - 1)) * (X - 1)) * (X ^ 2 + X + 1) =
        (C (-1 / (a - 1)) + C (1 / (a - 1))) * X ^ 3 +
          (C (-1 / (a - 1)) * (-C a) -
            C (1 / (a - 1))) := by ring
    _ = 1 := by rw [hcoeff, hconstant]; ring

/-- The family polynomial is separable exactly on the range used in the
paper.  The hypotheses also make visible the two possible degenerations:
the cubic ramifies at `a = 0`, and the factors meet at `a = 1`. -/
theorem polynomial_separable (a : ℚ) (ha0 : a ≠ 0) (ha1 : a ≠ 1) :
    (polynomial a).Separable := by
  exact (cubic_separable a ha0).mul quadratic_separable
    (cubic_isCoprime_quadratic a ha1)

/-- Separability survives the centering used by the fixed inverse equation. -/
theorem centeredPolynomial_separable
    (a : ℚ) (ha0 : a ≠ 0) (ha1 : a ≠ 1) :
    (centeredPolynomial a).Separable := by
  rw [← translate_polynomial_eq_centered]
  exact translatePolynomial_separable (polynomial a) (-1 / 2)
    (polynomial_separable a ha0 ha1)

/-- Separability in the exact syntactic form consumed by the literal-fiber
theorem. -/
theorem inversePolynomial_separable
    (a : ℚ) (ha0 : a ≠ 0) (ha1 : a ≠ 1) :
    (generalGaugeInversePolynomial
      seed ((1 : ℚˣ) : ℚ) (targetB a) (targetC a)).Separable := by
  change
    (generalGaugeInversePolynomial
      seed 1 (targetB a) (targetC a)).Separable
  rw [inversePolynomial_eq_centered]
  exact centeredPolynomial_separable a ha0 ha1

@[simp]
theorem polynomial_natDegree (a : ℚ) : (polynomial a).natDegree = 5 := by
  unfold polynomial cubic quadratic
  compute_degree!

/-- Every reduced member of the family represents a rank-five algebra. -/
theorem quotient_rank (a : ℚ) :
    Module.finrank ℚ (AdjoinRoot (polynomial a)) = 5 := by
  rw [adjoinRoot_finrank_eq_natDegree, polynomial_natDegree]

/-- Every reduced member of the family is finite étale. -/
theorem quotient_etale (a : ℚ) (ha0 : a ≠ 0) (ha1 : a ≠ 1) :
    Algebra.Etale ℚ (AdjoinRoot (polynomial a)) :=
  adjoinRoot_etale_of_separable (polynomial a)
    (polynomial_separable a ha0 ha1)

/-- Every reduced member of the family is finite as a rational module. -/
theorem quotient_finite (a : ℚ) (ha0 : a ≠ 0) (ha1 : a ≠ 1) :
    Module.Finite ℚ (AdjoinRoot (polynomial a)) :=
  adjoinRoot_finite_of_separable (polynomial a)
    (polynomial_separable a ha0 ha1)

/-- The fixed base map has the universal quadratic-gauge determinant `-2`. -/
theorem jacobianDet_baseMap :
    jacobianDet baseMap = MvPolynomial.C (-2) := by
  exact jacobianDet_generalGaugeMap seed
    seed_coeff_one_ne_zero seed_coeff_three_ne_zero

/-- The fixed base map has geometric degree five. -/
theorem baseMap_geometricDegree :
    generalGaugeGeometricDegree seed
      seed_coeff_one_ne_zero seed_coeff_three_ne_zero = 5 := by
  rw [generalGaugeGeometricDegree_eq seed
    seed_coeff_one_ne_zero seed_coeff_three_ne_zero (by simp)]
  exact seed_natDegree

/-- Translation identifies the centered inverse quotient with the
uncentered arithmetic quotient. -/
def inverseQuotientAlgEquiv (a : ℚ) :
    AdjoinRoot
        (generalGaugeInversePolynomial
          seed 1 (targetB a) (targetC a)) ≃ₐ[ℚ]
      AdjoinRoot (polynomial a) := by
  rw [inversePolynomial_eq_translate]
  exact translationQuotientEquiv (polynomial a) (-1 / 2)

section Fiber

variable {A B : Type*}
variable [CommRing A] [Algebra ℚ A]
variable [CommRing B] [Algebra ℚ B]

/-- The complete literal fiber of the one fixed determinant-`-2` map is
naturally represented by its centered inverse quotient.  The preceding
`inverseQuotientAlgEquiv` identifies this quotient with the uncentered
arithmetic quotient. -/
def centeredRawFiberRepresentingEquiv
    (a : ℚ) (ha0 : a ≠ 0) (ha1 : a ≠ 1) :
    (AdjoinRoot
        (generalGaugeInversePolynomial
          seed 1 (targetB a) (targetC a)) →ₐ[ℚ] A) ≃
      GeneralGaugeRawFiberPoint
        seed 1 (targetB a) (targetC a) A :=
  generalGaugeRawRepresentingEquiv
    seed 1 (targetB a) (targetC a)
    seed_coeff_one_ne_zero seed_coeff_three_ne_zero
    (inversePolynomial_separable a ha0 ha1) A

/-- Naturality of the complete fixed-map fiber identification. -/
theorem centeredRawFiberRepresentingEquiv_natural
    (a : ℚ) (ha0 : a ≠ 0) (ha1 : a ≠ 1)
    (f : A →ₐ[ℚ] B)
    (φ : AdjoinRoot
      (generalGaugeInversePolynomial
        seed 1 (targetB a) (targetC a)) →ₐ[ℚ] A) :
    GeneralGaugeRawFiberPoint.map f
        (centeredRawFiberRepresentingEquiv
          (A := A) a ha0 ha1 φ) =
      centeredRawFiberRepresentingEquiv
        (A := B) a ha0 ha1 (f.comp φ) :=
  generalGaugeRawRepresentingEquiv_natural
    seed 1 (targetB a) (targetC a)
    seed_coeff_one_ne_zero seed_coeff_three_ne_zero
    (inversePolynomial_separable a ha0 ha1) f φ

end Fiber

/-- The cyclotomic quadratic has no rational root. -/
theorem quadratic_no_rational_root (r : ℚ) :
    ¬quadratic.IsRoot r := by
  intro hr
  have h : r ^ 2 + r + 1 = 0 := by
    simpa [quadratic, Polynomial.IsRoot] using hr
  nlinarith [sq_nonneg (r + 1 / 2)]

/-- If `a` is not a rational cube, the family polynomial has no rational
root. -/
theorem polynomial_no_rational_root
    (a : ℚ) (ha : ¬∃ r : ℚ, r ^ 3 = a) (r : ℚ) :
    ¬(polynomial a).IsRoot r := by
  intro hr
  have hprod :
      (r ^ 3 - a) * (r ^ 2 + r + 1) = 0 := by
    simpa [polynomial, cubic, quadratic, Polynomial.IsRoot] using hr
  rcases mul_eq_zero.mp hprod with hcubic | hquadratic
  · apply ha
    exact ⟨r, sub_eq_zero.mp hcubic⟩
  · exact quadratic_no_rational_root r (by
      simpa [quadratic, Polynomial.IsRoot] using hquadratic)

/-- The family polynomial always has a real root, supplied by the cubic
factor. -/
theorem polynomial_has_real_root (a : ℚ) :
    Nonempty (PolynomialRoot (polynomial a) ℝ) := by
  let A : ℝ := |(a : ℝ)| + 1
  have hA1 : 1 ≤ A := by
    dsimp [A]
    exact le_add_of_nonneg_left (abs_nonneg _)
  have hAleCube : A ≤ A ^ 3 := by
    nlinarith [sq_nonneg A]
  have haUpper : (a : ℝ) ≤ A := by
    dsimp [A]
    linarith [le_abs_self (a : ℝ)]
  have haLower : -A ≤ (a : ℝ) := by
    dsimp [A]
    linarith [neg_abs_le (a : ℝ)]
  have hcont : Continuous (fun x : ℝ => x ^ 3 - (a : ℝ)) := by
    fun_prop
  have hzero :
      (0 : ℝ) ∈
        Set.Icc ((-A) ^ 3 - (a : ℝ)) (A ^ 3 - (a : ℝ)) := by
    constructor <;> nlinarith
  obtain ⟨r, -, hr⟩ := (Set.mem_image ..).mp
    (intermediate_value_Icc (by linarith : -A ≤ A)
      hcont.continuousOn hzero)
  refine ⟨⟨r, ?_⟩⟩
  simp [polynomial, cubic, quadratic, Polynomial.aeval_def, hr]

#print axioms inversePolynomial_eq_centered
#print axioms polynomial_separable
#print axioms jacobianDet_baseMap
#print axioms baseMap_geometricDegree
#print axioms centeredRawFiberRepresentingEquiv_natural
#print axioms polynomial_no_rational_root
#print axioms polynomial_has_real_root

end FiniteEtaleKeller.FixedHasseFamily
