/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.NumberTheory.Padics.Hensel
import Mathlib.RingTheory.ZMod.UnitsCyclic
import FiniteEtaleKeller.FixedHasseMap

/-!
# Uniform local points in the fixed Hasse family

This module formalizes the prime-by-prime local argument for

`(X^3-a)(X^2+X+1)`.

The hypotheses are stated in the exact elementary form used by the paper:
`a` is odd, `a ≡ 1 (mod 9)`, and every prime divisor of `a` is `1 (mod 3)`.
The oddness hypothesis is separated out to keep the local theorem independent
of a small auxiliary divisor argument.  The resulting theorem supplies a root
over `ℚ_[p]` for every prime `p`.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller.FixedHasseFamily

/-- Integral cubic factor for a natural parameter. -/
def cubicInt (a : ℕ) : ℤ[X] :=
  X ^ 3 - C (a : ℤ)

/-- Integral cyclotomic quadratic factor. -/
def quadraticInt : ℤ[X] :=
  X ^ 2 + X + 1

/-- If every prime divisor of `a` is `1 (mod 3)`, then `a` is odd.  This
discharges the auxiliary parity hypothesis from the paper-facing theorem. -/
theorem odd_of_prime_support
    (a : ℕ)
    (hsupport :
      ∀ q : ℕ, q.Prime → q ∣ a → q % 3 = 1) :
    Odd a := by
  rcases Nat.even_or_odd a with heven | hodd
  · have hmod := hsupport 2 Nat.prime_two heven.two_dvd
    norm_num at hmod
  · exact hodd

/-- The simple-root form of Hensel's lemma used throughout the local table. -/
private theorem hensel_of_simple_zmod_root
    (p : ℕ) [Fact p.Prime] (F : ℤ[X]) (a : ZMod p)
    (hroot : Polynomial.aeval a F = 0)
    (hderiv : Polynomial.aeval a F.derivative ≠ 0) :
    ∃ z : ℤ_[p], Polynomial.aeval z F = 0 := by
  let a₀ : ℤ_[p] := (a.val : ℕ)
  have ha₀ : PadicInt.toZMod a₀ = a := by
    dsimp [a₀]
    rw [map_natCast, ZMod.natCast_zmod_val]
  have hcomp :
      PadicInt.toZMod.comp (algebraMap ℤ ℤ_[p]) =
        algebraMap ℤ (ZMod p) := by
    ext n
    simp
  have hroot₀ :
      PadicInt.toZMod (Polynomial.aeval a₀ F) = 0 := by
    rw [Polynomial.aeval_def,
      Polynomial.hom_eval₂ F (algebraMap ℤ ℤ_[p]) PadicInt.toZMod,
      hcomp, ha₀]
    exact hroot
  have hderiv₀ :
      PadicInt.toZMod (Polynomial.aeval a₀ F.derivative) ≠ 0 := by
    rw [Polynomial.aeval_def,
      Polynomial.hom_eval₂ F.derivative (algebraMap ℤ ℤ_[p])
        PadicInt.toZMod, hcomp, ha₀]
    exact hderiv
  have hvalue_lt_one :
      ‖Polynomial.aeval a₀ F‖ < 1 := by
    rw [PadicInt.norm_lt_one_iff_dvd]
    rw [← Ideal.mem_span_singleton, ← PadicInt.maximalIdeal_eq_span_p,
      ← PadicInt.ker_toZMod, RingHom.mem_ker]
    exact hroot₀
  have hderiv_norm :
      ‖Polynomial.aeval a₀ F.derivative‖ = 1 := by
    apply le_antisymm (PadicInt.norm_le_one _)
    by_contra hnot
    have hlt : ‖Polynomial.aeval a₀ F.derivative‖ < 1 :=
      lt_of_not_ge hnot
    rw [PadicInt.norm_lt_one_iff_dvd] at hlt
    have hker :
        PadicInt.toZMod (Polynomial.aeval a₀ F.derivative) = 0 := by
      rw [← RingHom.mem_ker, PadicInt.ker_toZMod,
        PadicInt.maximalIdeal_eq_span_p, Ideal.mem_span_singleton]
      exact hlt
    exact hderiv₀ hker
  have hnorm :
      ‖Polynomial.aeval a₀ F‖ <
        ‖Polynomial.aeval a₀ F.derivative‖ ^ 2 := by
    rw [hderiv_norm, one_pow]
    exact hvalue_lt_one
  obtain ⟨z, hz, -⟩ := hensels_lemma hnorm
  exact ⟨z, hz⟩

/-- At primes `p ≡ 1 (mod 3)`, the cyclotomic quadratic has a simple
residue-field root. -/
private theorem quadratic_simple_zmod_root_of_mod_three_eq_one
    (p : ℕ) (hp : p.Prime) (hmod : p % 3 = 1) :
    ∃ r : ZMod p,
      Polynomial.aeval r quadraticInt = 0 ∧
        Polynomial.aeval r quadraticInt.derivative ≠ 0 := by
  letI : Fact p.Prime := ⟨hp⟩
  have hdiv : 3 ∣ Fintype.card (ZMod p)ˣ := by
    rw [ZMod.card_units, Nat.dvd_iff_mod_eq_zero]
    omega
  obtain ⟨r, hrorder⟩ :=
    exists_prime_orderOf_dvd_card (G := (ZMod p)ˣ) 3 hdiv
  have hrpowUnits : r ^ 3 = 1 := by
    simpa [hrorder] using pow_orderOf_eq_one r
  have hrpow : (r : ZMod p) ^ 3 = 1 :=
    congrArg Units.val hrpowUnits
  have hrne : (r : ZMod p) ≠ 1 := by
    intro hr
    have hru : r = 1 := Units.ext hr
    rw [hru, orderOf_one] at hrorder
    norm_num at hrorder
  have hfactor :
      ((r : ZMod p) - 1) *
          ((r : ZMod p) ^ 2 + (r : ZMod p) + 1) = 0 := by
    calc
      ((r : ZMod p) - 1) *
          ((r : ZMod p) ^ 2 + (r : ZMod p) + 1) =
          (r : ZMod p) ^ 3 - 1 := by ring
      _ = 0 := sub_eq_zero.mpr hrpow
  have hquad :
      (r : ZMod p) ^ 2 + (r : ZMod p) + 1 = 0 :=
    (mul_eq_zero.mp hfactor).resolve_left (sub_ne_zero.mpr hrne)
  refine ⟨r, ?_, ?_⟩
  · simpa [quadraticInt, Polynomial.aeval_def] using hquad
  · have hp3 : p ≠ 3 := by omega
    have h3 : (3 : ZMod p) ≠ 0 := by
      apply (CharP.cast_eq_zero_iff (ZMod p) p 3).not.mpr
      intro hdiv3
      have hor : p = 1 ∨ p = 3 :=
        (Nat.dvd_prime Nat.prime_three).mp hdiv3
      exact hor.elim hp.ne_one hp3
    have hlinear : (2 : ZMod p) * (r : ZMod p) + 1 ≠ 0 := by
      intro hz
      have hid :
          (4 : ZMod p) *
                ((r : ZMod p) ^ 2 + (r : ZMod p) + 1) -
              ((2 : ZMod p) * (r : ZMod p) + 1) ^ 2 = 3 := by
        ring
      rw [hquad, hz] at hid
      norm_num at hid
      exact h3 hid.symm
    simpa [quadraticInt, Polynomial.aeval_def] using hlinear

/-- At primes `p ≡ 2 (mod 3)` not dividing `a`, cubing is bijective in the
residue field, so the cubic has a simple root. -/
private theorem cubic_simple_zmod_root_of_mod_three_eq_two
    (a p : ℕ) (hp : p.Prime) (hpa : ¬p ∣ a)
    (hmod : p % 3 = 2) :
    ∃ r : ZMod p,
      Polynomial.aeval r (cubicInt a) = 0 ∧
        Polynomial.aeval r (cubicInt a).derivative ≠ 0 := by
  letI : Fact p.Prime := ⟨hp⟩
  have ha : (a : ZMod p) ≠ 0 := by
    exact (CharP.cast_eq_zero_iff (ZMod p) p a).not.mpr hpa
  let u : (ZMod p)ˣ := Units.mk0 a ha
  have hnot : ¬3 ∣ p - 1 := by
    intro hdiv
    rw [Nat.dvd_iff_mod_eq_zero] at hdiv
    omega
  have hcoprime : (Nat.card (ZMod p)ˣ).Coprime 3 := by
    rw [Nat.card_eq_fintype_card, ZMod.card_units]
    exact ((Nat.prime_three.coprime_iff_not_dvd.mpr hnot).symm)
  obtain ⟨r, hr⟩ :=
    (Nat.Coprime.pow_left_bijective (G := (ZMod p)ˣ) hcoprime).2 u
  refine ⟨r, ?_, ?_⟩
  · have hrcoe : (r : ZMod p) ^ 3 = a :=
      congrArg Units.val hr
    simp [cubicInt, Polynomial.aeval_def, hrcoe]
  · have hp3 : p ≠ 3 := by omega
    have h3 : (3 : ZMod p) ≠ 0 := by
      apply (CharP.cast_eq_zero_iff (ZMod p) p 3).not.mpr
      intro hdiv3
      have hor : p = 1 ∨ p = 3 :=
        (Nat.dvd_prime Nat.prime_three).mp hdiv3
      exact hor.elim hp.ne_one hp3
    have hr0 : (r : ZMod p) ≠ 0 := Units.ne_zero r
    simp [cubicInt, Polynomial.aeval_def, h3, hr0]

/-- An odd parameter gives the simple residue root `1` at `p=2`. -/
theorem cubicInt_has_twoAdicInt_root
    (a : ℕ) (ha : Odd a) :
    ∃ z : ℤ_[2], Polynomial.aeval z (cubicInt a) = 0 := by
  have haZMod : (a : ZMod 2) = 1 := ha.natCast_zmod_two
  apply hensel_of_simple_zmod_root 2 (cubicInt a) (1 : ZMod 2)
  · simp [cubicInt, Polynomial.aeval_def, haZMod]
  · norm_num [cubicInt, Polynomial.aeval_def]
    decide

/-- The natural quotient appearing after writing `a = 1 + 9m`. -/
def parameterQuotient (a : ℕ) : ℕ :=
  a / 9

/-- The transformed simple-root polynomial used at `p=3`. -/
def threeAuxiliaryInt (a : ℕ) : ℤ[X] :=
  X + 3 * X ^ 2 + 3 * X ^ 3 - C (parameterQuotient a : ℤ)

/-- If `a ≡ 1 (mod 9)`, the cubic has a `3`-adic integral root. -/
theorem cubicInt_has_threeAdicInt_root
    (a : ℕ) (ha : a % 9 = 1) :
    ∃ z : ℤ_[3], Polynomial.aeval z (cubicInt a) = 0 := by
  let r : ZMod 3 := parameterQuotient a
  have hthree : (3 : ZMod 3) = 0 := by decide
  have hroot :
      Polynomial.aeval r (threeAuxiliaryInt a) = 0 := by
    simp [threeAuxiliaryInt, Polynomial.aeval_def, r, hthree]
  have hderiv :
      Polynomial.aeval r (threeAuxiliaryInt a).derivative ≠ 0 := by
    simp [threeAuxiliaryInt, Polynomial.aeval_def, r, hthree]
  obtain ⟨t, ht⟩ :=
    hensel_of_simple_zmod_root 3 (threeAuxiliaryInt a) r hroot hderiv
  refine ⟨1 + 3 * t, ?_⟩
  have haDecomp : a = 1 + 9 * parameterQuotient a := by
    have h := Nat.mod_add_div a 9
    simp only [parameterQuotient]
    omega
  have htExpr :
      t + 3 * t ^ 2 + 3 * t ^ 3 -
          (parameterQuotient a : ℤ_[3]) = 0 := by
    simpa [threeAuxiliaryInt, Polynomial.aeval_def] using ht
  have haCast :
      (a : ℤ_[3]) =
        1 + 9 * (parameterQuotient a : ℤ_[3]) := by
    simpa using congrArg (fun n : ℕ => (n : ℤ_[3])) haDecomp
  have hgoal : (1 + 3 * t) ^ 3 - (a : ℤ_[3]) = 0 := by
    rw [haCast]
    linear_combination 9 * htExpr
  simpa [cubicInt, Polynomial.aeval_def] using hgoal

private theorem quadratic_has_padicInt_root_of_mod_three_eq_one
    (p : ℕ) [Fact p.Prime] (hmod : p % 3 = 1) :
    ∃ z : ℤ_[p], Polynomial.aeval z quadraticInt = 0 := by
  obtain ⟨r, hr, hdr⟩ :=
    quadratic_simple_zmod_root_of_mod_three_eq_one p Fact.out hmod
  exact hensel_of_simple_zmod_root p quadraticInt r hr hdr

private theorem cubic_has_padicInt_root_of_mod_three_eq_two
    (a p : ℕ) [Fact p.Prime] (hpa : ¬p ∣ a)
    (hmod : p % 3 = 2) :
    ∃ z : ℤ_[p], Polynomial.aeval z (cubicInt a) = 0 := by
  obtain ⟨r, hr, hdr⟩ :=
    cubic_simple_zmod_root_of_mod_three_eq_two a p Fact.out hpa hmod
  exact hensel_of_simple_zmod_root p (cubicInt a) r hr hdr

private theorem polynomial_has_padic_root_of_cubicInt_root
    (a p : ℕ) [Fact p.Prime] (z : ℤ_[p])
    (hz : Polynomial.aeval z (cubicInt a) = 0) :
    ∃ x : ℚ_[p],
      Polynomial.aeval x (polynomial (a : ℚ)) = 0 := by
  refine ⟨(z : ℚ_[p]), ?_⟩
  have hzExpr : z ^ 3 - a = 0 := by
    simpa [cubicInt, Polynomial.aeval_def] using hz
  have hzExpr' : (z : ℚ_[p]) ^ 3 - a = 0 := by
    have hmap := congrArg (algebraMap ℤ_[p] ℚ_[p]) hzExpr
    simp only [map_sub, map_pow, map_natCast, map_zero,
      PadicInt.algebraMap_apply] at hmap
    exact hmap
  simp [polynomial, cubic, quadratic, Polynomial.aeval_def, hzExpr']

private theorem polynomial_has_padic_root_of_quadraticInt_root
    (a p : ℕ) [Fact p.Prime] (z : ℤ_[p])
    (hz : Polynomial.aeval z quadraticInt = 0) :
    ∃ x : ℚ_[p],
      Polynomial.aeval x (polynomial (a : ℚ)) = 0 := by
  refine ⟨(z : ℚ_[p]), ?_⟩
  have hzExpr : z ^ 2 + z + 1 = 0 := by
    simpa [quadraticInt, Polynomial.aeval_def] using hz
  have hzExpr' : (z : ℚ_[p]) ^ 2 + (z : ℚ_[p]) + 1 = 0 := by
    have hmap := congrArg (algebraMap ℤ_[p] ℚ_[p]) hzExpr
    simp only [map_add, map_pow, map_one, map_zero,
      PadicInt.algebraMap_apply] at hmap
    exact hmap
  simp [polynomial, cubic, quadratic, Polynomial.aeval_def, hzExpr']

/-- Uniform finite-place local solubility for the paper's parameter
conditions. -/
theorem polynomial_has_padic_root
    (a : ℕ) (hodd : Odd a) (hmod9 : a % 9 = 1)
    (hsupport :
      ∀ q : ℕ, q.Prime → q ∣ a → q % 3 = 1)
    (p : ℕ) [Fact p.Prime] :
    ∃ z : ℚ_[p],
      Polynomial.aeval z (polynomial (a : ℚ)) = 0 := by
  have hp : p.Prime := Fact.out
  by_cases hp2 : p = 2
  · subst p
    obtain ⟨z, hz⟩ := cubicInt_has_twoAdicInt_root a hodd
    exact polynomial_has_padic_root_of_cubicInt_root a 2 z hz
  by_cases hp3 : p = 3
  · subst p
    obtain ⟨z, hz⟩ := cubicInt_has_threeAdicInt_root a hmod9
    exact polynomial_has_padic_root_of_cubicInt_root a 3 z hz
  have hmod :
      p % 3 = 1 ∨ p % 3 = 2 := by
    have hlt : p % 3 < 3 := Nat.mod_lt p (by norm_num)
    have hnezero : p % 3 ≠ 0 := by
      intro hz
      have hdiv : 3 ∣ p := Nat.dvd_of_mod_eq_zero hz
      have hor : 3 = 1 ∨ 3 = p := (Nat.dvd_prime hp).mp hdiv
      omega
    omega
  rcases hmod with hmod1 | hmod2
  · obtain ⟨z, hz⟩ :=
      quadratic_has_padicInt_root_of_mod_three_eq_one p hmod1
    exact polynomial_has_padic_root_of_quadraticInt_root a p z hz
  · have hpa : ¬p ∣ a := by
      intro hdiv
      have := hsupport p hp hdiv
      omega
    obtain ⟨z, hz⟩ :=
      cubic_has_padicInt_root_of_mod_three_eq_two a p hpa hmod2
    exact polynomial_has_padic_root_of_cubicInt_root a p z hz

/-- Under the uniform arithmetic hypotheses, every nonarchimedean completion
contains a point of the literal fiber of the exact paper map. -/
theorem paperFiberPoint_padic_nonempty
    (a : ℕ) (hodd : Odd a) (hmod9 : a % 9 = 1)
    (hsupport :
      ∀ q : ℕ, q.Prime → q ∣ a → q % 3 = 1)
    (ha1 : a ≠ 1)
    (p : ℕ) [Fact p.Prime] :
    Nonempty (PaperFiberPoint (a : ℚ) ℚ_[p]) := by
  have ha0 : a ≠ 0 := by
    intro ha
    subst a
    norm_num at hmod9
  have ha0Q : (a : ℚ) ≠ 0 := by exact_mod_cast ha0
  have ha1Q : (a : ℚ) ≠ 1 := by exact_mod_cast ha1
  obtain ⟨z, hz⟩ :=
    polynomial_has_padic_root a hodd hmod9 hsupport p
  let root : PolynomialRoot (polynomial (a : ℚ)) ℚ_[p] := ⟨z, hz⟩
  exact ⟨paperFiberRepresentingEquiv
    (A := ℚ_[p]) (a : ℚ) ha0Q ha1Q
    ((PolynomialRoot.algHomEquiv
      (polynomial (a : ℚ)) ℚ_[p]).symm root)⟩

/-- End-to-end Hasse-principle certificate with parity stated explicitly.
The public paper-facing theorem below derives parity from prime support. -/
theorem paperFiberPoint_hasse_certificate_of_odd
    (a : ℕ) (hodd : Odd a) (hmod9 : a % 9 = 1)
    (hsupport :
      ∀ q : ℕ, q.Prime → q ∣ a → q % 3 = 1)
    (hnoncube : ¬∃ r : ℚ, r ^ 3 = (a : ℚ)) :
    IsEmpty (PaperFiberPoint (a : ℚ) ℚ) ∧
      Nonempty (PaperFiberPoint (a : ℚ) ℝ) ∧
        ∀ (p : ℕ) [Fact p.Prime],
          Nonempty (PaperFiberPoint (a : ℚ) ℚ_[p]) := by
  have ha0 : a ≠ 0 := by
    intro ha
    subst a
    norm_num at hmod9
  have ha1 : a ≠ 1 := by
    intro ha
    subst a
    apply hnoncube
    exact ⟨1, by norm_num⟩
  have ha0Q : (a : ℚ) ≠ 0 := by exact_mod_cast ha0
  have ha1Q : (a : ℚ) ≠ 1 := by exact_mod_cast ha1
  exact
    ⟨paperFiberPoint_rat_isEmpty (a : ℚ) ha0Q ha1Q hnoncube,
      paperFiberPoint_real_nonempty (a : ℚ) ha0Q ha1Q,
      paperFiberPoint_padic_nonempty a hodd hmod9 hsupport ha1⟩

/-- End-to-end Hasse-principle certificate under exactly the elementary
parameter conditions displayed in the paper. -/
theorem paperFiberPoint_hasse_certificate
    (a : ℕ) (_ha : 1 < a) (hmod9 : a % 9 = 1)
    (hsupport :
      ∀ q : ℕ, q.Prime → q ∣ a → q % 3 = 1)
    (hnoncube : ¬∃ r : ℚ, r ^ 3 = (a : ℚ)) :
    IsEmpty (PaperFiberPoint (a : ℚ) ℚ) ∧
      Nonempty (PaperFiberPoint (a : ℚ) ℝ) ∧
        ∀ (p : ℕ) [Fact p.Prime],
          Nonempty (PaperFiberPoint (a : ℚ) ℚ_[p]) :=
  paperFiberPoint_hasse_certificate_of_odd
    a (odd_of_prime_support a hsupport) hmod9 hsupport hnoncube

#print axioms odd_of_prime_support
#print axioms cubicInt_has_twoAdicInt_root
#print axioms cubicInt_has_threeAdicInt_root
#print axioms polynomial_has_padic_root
#print axioms paperFiberPoint_padic_nonempty
#print axioms paperFiberPoint_hasse_certificate_of_odd
#print axioms paperFiberPoint_hasse_certificate

end FiniteEtaleKeller.FixedHasseFamily
