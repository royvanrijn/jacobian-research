/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeInverse
import Mathlib.Algebra.MvPolynomial.Equiv
import Mathlib.FieldTheory.RatFunc.Basic
import Mathlib.RingTheory.Polynomial.GaussLemma

/-!
# Irreducibility and degree of the generic inverse equation

This module formalizes the paper's Gauss-lemma argument for
`E(S) = H(S) - λ C`. A polynomial-variable swap turns this equation into a
degree-one polynomial in `C`, whose leading coefficient is a unit. Transport
through the swap proves irreducibility over `K[C]`; Gauss's lemma then proves
irreducibility over `K(C)`.

The final certificate applies this argument to the all-degree quadratic gauge.
For a seed of degree `N ≥ 3`, fixed parameters `pi ≠ 0` and `b`, and a
nonzero linear seed coefficient, the one-parameter inverse polynomial
`E_(pi,b,C)` is irreducible of degree `N` over `K(C)`. Specializing `C`
recovers the inverse polynomial already used by the literal-fiber theorems.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller.GenericInverse

variable (K : Type*) [Field K]

/-- Swap the two variables of an iterated univariate polynomial ring. -/
def polynomialVariableSwap : K[X][X] ≃+* K[X][X] :=
  (Polynomial.mapEquiv
      (MvPolynomial.uniqueAlgEquiv K Unit).symm.toRingEquiv).trans
    (MvPolynomial.optionEquivLeft K Unit).symm.toRingEquiv |>.trans
    (MvPolynomial.optionEquivRight K Unit).toRingEquiv |>.trans
    (MvPolynomial.uniqueAlgEquiv K[X] Unit).toRingEquiv

@[simp]
theorem polynomialVariableSwap_X :
    polynomialVariableSwap K (X : K[X][X]) = C (X : K[X]) := by
  simp [polynomialVariableSwap]

@[simp]
theorem polynomialVariableSwap_C_X :
    polynomialVariableSwap K (C (X : K[X])) = (X : K[X][X]) := by
  simp [polynomialVariableSwap]

@[simp]
theorem polynomialVariableSwap_C_C (r : K) :
    polynomialVariableSwap K (C (C r)) = C (C r) := by
  simp [polynomialVariableSwap]

theorem polynomialVariableSwap_map_C (H : K[X]) :
    polynomialVariableSwap K (H.map (C : K →+* K[X])) = C H := by
  induction H using Polynomial.induction_on' with
  | add p q hp hq =>
      rw [Polynomial.map_add, map_add, hp, hq, map_add]
  | monomial n r =>
      rw [← Polynomial.C_mul_X_pow_eq_monomial]
      simp

/-- The bivariate equation `H(S) - λ C`, with `C` in the coefficient ring. -/
def linearTargetPolynomial (H : K[X]) (lambda : K) : K[X][X] :=
  H.map (C : K →+* K[X]) - C (C lambda * X)

theorem polynomialVariableSwap_linearTargetPolynomial
    (H : K[X]) (lambda : K) :
    polynomialVariableSwap K (linearTargetPolynomial K H lambda) =
      C (-(C lambda : K[X])) * X + C H := by
  rw [linearTargetPolynomial, map_sub,
    polynomialVariableSwap_map_C]
  simp
  ring

/-- A nontrivial equation `H(S) - λ C` is irreducible over `K[C]`. -/
theorem linearTargetPolynomial_irreducible
    (H : K[X]) (lambda : K) (hlambda : lambda ≠ 0) :
    Irreducible (linearTargetPolynomial K H lambda) := by
  have hunit : IsUnit (-(C lambda : K[X])) :=
    (Polynomial.isUnit_C.mpr (isUnit_iff_ne_zero.mpr hlambda)).neg
  have hirr :
      Irreducible (C (-(C lambda : K[X])) * X + C H) :=
    irreducible_C_mul_X_add_C hunit.ne_zero hunit.isRelPrime_left
  rw [← polynomialVariableSwap_linearTargetPolynomial K H lambda] at hirr
  exact Irreducible.of_map hirr

/-- Subtracting the parameter term does not change a positive source degree. -/
theorem linearTargetPolynomial_natDegree
    (H : K[X]) (lambda : K) (hH : H.natDegree ≠ 0) :
    (linearTargetPolynomial K H lambda).natDegree = H.natDegree := by
  rw [linearTargetPolynomial, natDegree_sub_eq_left_of_natDegree_lt]
  · exact natDegree_map_eq_of_injective C_injective H
  · rw [natDegree_C, natDegree_map_eq_of_injective C_injective]
    exact Nat.pos_of_ne_zero hH

/-- Gauss's lemma transports the irreducibility certificate to `K(C)[S]`. -/
theorem linearTargetPolynomial_ratFunc_irreducible
    (H : K[X]) (lambda : K) (hlambda : lambda ≠ 0)
    (hH : H.natDegree ≠ 0) :
    Irreducible
      ((linearTargetPolynomial K H lambda).map
        (algebraMap K[X] (RatFunc K))) := by
  have hirr := linearTargetPolynomial_irreducible K H lambda hlambda
  have hprimitive :
      (linearTargetPolynomial K H lambda).IsPrimitive :=
    hirr.isPrimitive (by rw [linearTargetPolynomial_natDegree K H lambda hH]; exact hH)
  exact hprimitive.irreducible_iff_irreducible_map_fraction_map.mp hirr

end FiniteEtaleKeller.GenericInverse

namespace FiniteEtaleKeller

variable {K : Type*} [Field K]

/-- The source-variable part of the inverse equation before subtracting the
third target parameter. -/
def generalGaugeTargetFreeInversePolynomial
    (G : K[X]) (pi b : K) : K[X] :=
  generalGaugeSeedPolynomial G pi -
    C (G.coeff 1 / 2) * C b * X ^ 2

private theorem generalGaugeSeedPolynomial_natDegree_le
    (G : K[X]) (pi : K) (hdeg : 3 ≤ G.natDegree) :
    (generalGaugeSeedPolynomial G pi).natDegree ≤ G.natDegree := by
  rw [natDegree_le_iff_coeff_eq_zero]
  intro n hn
  rw [generalGaugeSeedPolynomial, coeff_add, coeff_add]
  have hn1 : n ≠ 1 := by omega
  have hn2 : n ≠ 2 := by omega
  have hn3 : n ≠ 3 := by omega
  simp only [coeff_C_mul, coeff_X_pow, coeff_add]
  have htail :
      (∑ k ∈ Finset.Icc 4 G.natDegree,
        C (G.coeff k * pi ^ k) * X ^ k).coeff n = 0 := by
    rw [Polynomial.finsetSum_coeff]
    apply Finset.sum_eq_zero
    intro k hk
    have hkN : k ≤ G.natDegree := (Finset.mem_Icc.mp hk).2
    have hnk : n ≠ k := by omega
    rw [coeff_C_mul_X_pow, if_neg hnk]
  rw [htail]
  simp [Polynomial.coeff_X, Ne.symm hn1, hn2, hn3]

private theorem generalGaugeSeedPolynomial_coeff_natDegree_ne_zero
    (G : K[X]) (pi : K) (hdeg : 3 ≤ G.natDegree)
    (hpi : pi ≠ 0) :
    (generalGaugeSeedPolynomial G pi).coeff G.natDegree ≠ 0 := by
  have hG : G ≠ 0 := by
    intro h
    simp [h] at hdeg
  have hlead : G.coeff G.natDegree ≠ 0 := by
    rw [coeff_natDegree]
    exact leadingCoeff_ne_zero.mpr hG
  rcases hdeg.eq_or_lt with hN | hN
  · have hN' : G.natDegree = 3 := hN.symm
    have hlead3 : G.coeff 3 ≠ 0 := by
      rw [← hN']
      exact hlead
    rw [generalGaugeSeedPolynomial, hN']
    simp [Polynomial.coeff_X_pow, hlead3, hpi]
  · have hN4 : 4 ≤ G.natDegree := by omega
    have hN1 : G.natDegree ≠ 1 := by omega
    have hN2 : G.natDegree ≠ 2 := by omega
    have hN3 : G.natDegree ≠ 3 := by omega
    have hlow1 :
        (C (G.coeff 1) * X).coeff G.natDegree = 0 := by
      rw [coeff_C_mul]
      simp [Polynomial.coeff_X, Ne.symm hN1]
    have hlow2 :
        (C pi * (C (G.coeff 2) * X ^ 2 +
          C (G.coeff 3) * X ^ 3)).coeff G.natDegree = 0 := by
      rw [coeff_C_mul]
      simp [hN2, hN3]
    have htail :
        (∑ k ∈ Finset.Icc 4 G.natDegree,
          C (G.coeff k * pi ^ k) * X ^ k).coeff G.natDegree =
            G.coeff G.natDegree * pi ^ G.natDegree := by
      rw [Polynomial.finsetSum_coeff, Finset.sum_eq_single G.natDegree]
      · rw [coeff_C_mul_X_pow, if_pos rfl]
      · intro k hk hkne
        rw [coeff_C_mul_X_pow, if_neg hkne.symm]
      · simp [hN4]
    rw [generalGaugeSeedPolynomial, coeff_add, coeff_add, hlow1, hlow2,
      htail, zero_add]
    simpa using mul_ne_zero hlead (pow_ne_zero G.natDegree hpi)

/-- A nonzero chart parameter preserves the seed degree in the gauge
substitution. -/
theorem generalGaugeSeedPolynomial_natDegree
    (G : K[X]) (pi : K) (hdeg : 3 ≤ G.natDegree)
    (hpi : pi ≠ 0) :
    (generalGaugeSeedPolynomial G pi).natDegree = G.natDegree :=
  natDegree_eq_of_le_of_coeff_ne_zero
    (generalGaugeSeedPolynomial_natDegree_le G pi hdeg)
    (generalGaugeSeedPolynomial_coeff_natDegree_ne_zero G pi hdeg hpi)

/-- The quadratic `B` term does not change the degree of a seed of degree at
least three. -/
theorem generalGaugeTargetFreeInversePolynomial_natDegree
    (G : K[X]) (pi b : K) (hdeg : 3 ≤ G.natDegree)
    (hpi : pi ≠ 0) :
    (generalGaugeTargetFreeInversePolynomial G pi b).natDegree =
      G.natDegree := by
  rw [generalGaugeTargetFreeInversePolynomial,
    natDegree_sub_eq_left_of_natDegree_lt]
  · exact generalGaugeSeedPolynomial_natDegree G pi hdeg hpi
  · calc
      (C (G.coeff 1 / 2) * C b * X ^ 2).natDegree ≤ 2 := by
        rw [← C_mul]
        exact natDegree_C_mul_X_pow_le ((G.coeff 1 / 2) * b) 2
      _ < (generalGaugeSeedPolynomial G pi).natDegree := by
        rw [generalGaugeSeedPolynomial_natDegree G pi hdeg hpi]
        omega

/-- The exact generic inverse equation over `K[C]`. -/
def generalGaugeTargetPolynomial
    (G : K[X]) (pi b : K) : K[X][X] :=
  GenericInverse.linearTargetPolynomial K
    (generalGaugeTargetFreeInversePolynomial G pi b)
    (G.coeff 1 / 2)

/-- Evaluating the formal target variable at `c` recovers the inverse
polynomial used by the literal-fiber construction. -/
theorem generalGaugeTargetPolynomial_specialize
    (G : K[X]) (pi b c : K) :
    (generalGaugeTargetPolynomial G pi b).map
        (Polynomial.evalRingHom c) =
      generalGaugeInversePolynomial G pi b c := by
  have hmap (P : K[X]) :
      (P.map (C : K →+* K[X])).map (Polynomial.evalRingHom c) = P := by
    ext n
    simp
  simp [generalGaugeTargetPolynomial,
    GenericInverse.linearTargetPolynomial,
    generalGaugeTargetFreeInversePolynomial,
    generalGaugeInversePolynomial, hmap]
  ring

/-- The exact gauge target polynomial is irreducible over `K[C]`. -/
theorem generalGaugeTargetPolynomial_irreducible
    [CharZero K] (G : K[X]) (pi b : K)
    (h₁ : G.coeff 1 ≠ 0) :
    Irreducible (generalGaugeTargetPolynomial G pi b) := by
  apply GenericInverse.linearTargetPolynomial_irreducible
  exact div_ne_zero h₁ (by norm_num)

/-- The generic inverse equation over the rational-function field `K(C)`. -/
def generalGaugeGenericInversePolynomial
    (G : K[X]) (pi b : K) : (RatFunc K)[X] :=
  (generalGaugeTargetPolynomial G pi b).map
    (algebraMap K[X] (RatFunc K))

/-- The paper's generic inverse equation is irreducible over `K(C)`. -/
theorem generalGaugeGenericInversePolynomial_irreducible
    [CharZero K] (G : K[X]) (pi b : K)
    (h₁ : G.coeff 1 ≠ 0)
    (hdeg : 3 ≤ G.natDegree) (hpi : pi ≠ 0) :
    Irreducible (generalGaugeGenericInversePolynomial G pi b) := by
  exact GenericInverse.linearTargetPolynomial_ratFunc_irreducible K
    (generalGaugeTargetFreeInversePolynomial G pi b)
    (G.coeff 1 / 2) (div_ne_zero h₁ (by norm_num))
    (by
      rw [generalGaugeTargetFreeInversePolynomial_natDegree G pi b hdeg hpi]
      omega)

/-- The generic inverse equation has exactly the seed degree. -/
theorem generalGaugeGenericInversePolynomial_natDegree
    (G : K[X]) (pi b : K)
    (hdeg : 3 ≤ G.natDegree) (hpi : pi ≠ 0) :
    (generalGaugeGenericInversePolynomial G pi b).natDegree =
      G.natDegree := by
  calc
    (generalGaugeGenericInversePolynomial G pi b).natDegree =
        (generalGaugeTargetPolynomial G pi b).natDegree := by
      exact natDegree_map_eq_of_injective
        (IsFractionRing.injective K[X] (RatFunc K)) _
    _ = (generalGaugeTargetFreeInversePolynomial G pi b).natDegree := by
      exact GenericInverse.linearTargetPolynomial_natDegree K
        (generalGaugeTargetFreeInversePolynomial G pi b)
        (G.coeff 1 / 2)
        (by
          rw [generalGaugeTargetFreeInversePolynomial_natDegree G pi b hdeg hpi]
          omega)
    _ = G.natDegree :=
      generalGaugeTargetFreeInversePolynomial_natDegree G pi b hdeg hpi

/-- Combined irreducibility-and-degree certificate for the paper's
fixed-`pi,b`, one-parameter inverse equation. -/
theorem generalGaugeGenericInversePolynomial_certificate
    [CharZero K] (G : K[X]) (pi b : K)
    (h₁ : G.coeff 1 ≠ 0)
    (hdeg : 3 ≤ G.natDegree) (hpi : pi ≠ 0) :
    Irreducible (generalGaugeGenericInversePolynomial G pi b) ∧
      (generalGaugeGenericInversePolynomial G pi b).natDegree =
        G.natDegree :=
  ⟨generalGaugeGenericInversePolynomial_irreducible
      G pi b h₁ hdeg hpi,
    generalGaugeGenericInversePolynomial_natDegree G pi b hdeg hpi⟩

/-- The quotient generated by a root of the fixed-parameter inverse equation
has extension dimension exactly the seed degree. Together with irreducibility,
Mathlib's `AdjoinRoot` field instance makes this the corresponding
degree-`N` field extension. -/
theorem generalGaugeGenericInverseAdjoinRoot_finrank
    (G : K[X]) (pi b : K)
    (hdeg : 3 ≤ G.natDegree) (hpi : pi ≠ 0) :
    Module.finrank (RatFunc K)
        (AdjoinRoot (generalGaugeGenericInversePolynomial G pi b)) =
      G.natDegree := by
  change Module.finrank (RatFunc K)
      ((RatFunc K)[X] ⧸
        Ideal.span {generalGaugeGenericInversePolynomial G pi b}) =
    G.natDegree
  rw [finrank_quotient_span_eq_natDegree,
    generalGaugeGenericInversePolynomial_natDegree G pi b hdeg hpi]

#print axioms generalGaugeTargetPolynomial_specialize
#print axioms generalGaugeTargetPolynomial_irreducible
#print axioms generalGaugeGenericInversePolynomial_certificate
#print axioms generalGaugeGenericInverseAdjoinRoot_finrank

end FiniteEtaleKeller
