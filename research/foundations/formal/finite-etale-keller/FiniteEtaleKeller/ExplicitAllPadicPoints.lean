/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.NumberTheory.Padics.Hensel
import Mathlib.RingTheory.ZMod.UnitsCyclic
import FiniteEtaleKeller.ExplicitThreeAdicPoint

/-!
# Points at every finite completion for the explicit quintic

This module formalizes the complete prime-by-prime table for
`(X³ - 19)(X² + X + 1)`.  The primes `2`, `3`, and `19` have explicit Hensel
witnesses; all other primes are divided by their residue modulo three.
-/

noncomputable section

namespace FiniteEtaleKeller.ExplicitQuintic

private instance : Fact (Nat.Prime 19) := ⟨by norm_num⟩

private theorem hensel_of_simple_zmod_root
    (p : ℕ) [Fact p.Prime] (F : Polynomial ℤ) (a : ZMod p)
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
      Polynomial.hom_eval₂ F.derivative (algebraMap ℤ ℤ_[p]) PadicInt.toZMod,
      hcomp, ha₀]
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

/-- The integral cyclotomic quadratic factor. -/
def quadraticThreeInt : Polynomial ℤ :=
  Polynomial.X ^ 2 + Polynomial.X + 1

private theorem cubic19_simple_zmod_root_of_mod_three_eq_two
    (p : ℕ) (hp : p.Prime) (hp19 : p ≠ 19)
    (hmod : p % 3 = 2) :
    ∃ a : ZMod p,
      Polynomial.aeval a cubic19Int = 0 ∧
        Polynomial.aeval a cubic19Int.derivative ≠ 0 := by
  letI : Fact p.Prime := ⟨hp⟩
  have h19 : (19 : ZMod p) ≠ 0 := by
    apply (CharP.cast_eq_zero_iff (ZMod p) p 19).not.mpr
    intro hdiv
    have hor : p = 1 ∨ p = 19 :=
      (Nat.dvd_prime (by norm_num : Nat.Prime 19)).mp hdiv
    exact hor.elim hp.ne_one hp19
  let u : (ZMod p)ˣ := Units.mk0 19 h19
  have hnot : ¬3 ∣ p - 1 := by
    intro hdiv
    rw [Nat.dvd_iff_mod_eq_zero] at hdiv
    omega
  have hcoprime : (Nat.card (ZMod p)ˣ).Coprime 3 := by
    rw [Nat.card_eq_fintype_card, ZMod.card_units]
    exact ((Nat.prime_three.coprime_iff_not_dvd.mpr hnot).symm)
  obtain ⟨a, ha⟩ :=
    (Nat.Coprime.pow_left_bijective (G := (ZMod p)ˣ) hcoprime).2 u
  refine ⟨a, ?_, ?_⟩
  · have hacoe : (a : ZMod p) ^ 3 = 19 :=
      congrArg Units.val ha
    simp [cubic19Int, Polynomial.aeval_def, hacoe]
  · have h3 : (3 : ZMod p) ≠ 0 := by
      apply (CharP.cast_eq_zero_iff (ZMod p) p 3).not.mpr
      intro hdiv
      have hor : p = 1 ∨ p = 3 :=
        (Nat.dvd_prime Nat.prime_three).mp hdiv
      exact hor.elim hp.ne_one (by omega)
    have ha0 : (a : ZMod p) ≠ 0 := Units.ne_zero a
    simp [cubic19Int, Polynomial.aeval_def, h3, ha0]

private theorem quadratic_simple_zmod_root_of_mod_three_eq_one
    (p : ℕ) (hp : p.Prime) (hmod : p % 3 = 1) :
    ∃ a : ZMod p,
      Polynomial.aeval a quadraticThreeInt = 0 ∧
        Polynomial.aeval a quadraticThreeInt.derivative ≠ 0 := by
  letI : Fact p.Prime := ⟨hp⟩
  have hdiv : 3 ∣ Fintype.card (ZMod p)ˣ := by
    rw [ZMod.card_units, Nat.dvd_iff_mod_eq_zero]
    omega
  obtain ⟨a, haorder⟩ :=
    exists_prime_orderOf_dvd_card (G := (ZMod p)ˣ) 3 hdiv
  have hapowUnits : a ^ 3 = 1 := by
    simpa [haorder] using pow_orderOf_eq_one a
  have hapow : (a : ZMod p) ^ 3 = 1 :=
    congrArg Units.val hapowUnits
  have hane : (a : ZMod p) ≠ 1 := by
    intro ha
    have hau : a = 1 := Units.ext ha
    rw [hau, orderOf_one] at haorder
    norm_num at haorder
  have hfactor :
      ((a : ZMod p) - 1) *
          ((a : ZMod p) ^ 2 + (a : ZMod p) + 1) = 0 := by
    calc
      ((a : ZMod p) - 1) *
          ((a : ZMod p) ^ 2 + (a : ZMod p) + 1) =
          (a : ZMod p) ^ 3 - 1 := by ring
      _ = 0 := sub_eq_zero.mpr hapow
  have hquad :
      (a : ZMod p) ^ 2 + (a : ZMod p) + 1 = 0 :=
    (mul_eq_zero.mp hfactor).resolve_left (sub_ne_zero.mpr hane)
  refine ⟨a, ?_, ?_⟩
  · simpa [quadraticThreeInt, Polynomial.aeval_def] using hquad
  · have hp3 : p ≠ 3 := by omega
    have h3 : (3 : ZMod p) ≠ 0 := by
      apply (CharP.cast_eq_zero_iff (ZMod p) p 3).not.mpr
      intro hdiv
      have hor : p = 1 ∨ p = 3 :=
        (Nat.dvd_prime Nat.prime_three).mp hdiv
      exact hor.elim hp.ne_one hp3
    have hlinear : (2 : ZMod p) * (a : ZMod p) + 1 ≠ 0 := by
      intro hz
      have hid :
          (4 : ZMod p) *
                ((a : ZMod p) ^ 2 + (a : ZMod p) + 1) -
              ((2 : ZMod p) * (a : ZMod p) + 1) ^ 2 = 3 := by
        ring
      rw [hquad, hz] at hid
      norm_num at hid
      exact h3 hid.symm
    simpa [quadraticThreeInt, Polynomial.aeval_def] using hlinear

/-- The exceptional prime `2`: the residue class `1` is a simple root of
`X³ - 19`. -/
theorem cubic19_has_twoAdicInt_root :
    ∃ z : ℤ_[2], Polynomial.aeval z cubic19Int = 0 := by
  apply hensel_of_simple_zmod_root 2 cubic19Int (1 : ZMod 2)
  · norm_num [cubic19Int, Polynomial.aeval_def]
    decide
  · norm_num [cubic19Int, Polynomial.aeval_def]
    decide

/-- The exceptional prime `19`: the residue class `7` is a simple root of
`X² + X + 1`. -/
theorem quadraticThree_has_nineteenAdicInt_root :
    ∃ z : ℤ_[19], Polynomial.aeval z quadraticThreeInt = 0 := by
  apply hensel_of_simple_zmod_root 19 quadraticThreeInt (7 : ZMod 19)
  · norm_num [quadraticThreeInt, Polynomial.aeval_def]
    decide
  · norm_num [quadraticThreeInt, Polynomial.aeval_def]
    decide

private theorem cubic19_has_padicInt_root_of_mod_three_eq_two
    (p : ℕ) [Fact p.Prime] (hp19 : p ≠ 19)
    (hmod : p % 3 = 2) :
    ∃ z : ℤ_[p], Polynomial.aeval z cubic19Int = 0 := by
  obtain ⟨a, ha, hda⟩ :=
    cubic19_simple_zmod_root_of_mod_three_eq_two p Fact.out hp19 hmod
  exact hensel_of_simple_zmod_root p cubic19Int a ha hda

private theorem quadraticThree_has_padicInt_root_of_mod_three_eq_one
    (p : ℕ) [Fact p.Prime] (hmod : p % 3 = 1) :
    ∃ z : ℤ_[p], Polynomial.aeval z quadraticThreeInt = 0 := by
  obtain ⟨a, ha, hda⟩ :=
    quadratic_simple_zmod_root_of_mod_three_eq_one p Fact.out hmod
  exact hensel_of_simple_zmod_root p quadraticThreeInt a ha hda

private theorem p5_has_padic_root_of_cubic19Int_root
    (p : ℕ) [Fact p.Prime] (z : ℤ_[p])
    (hz : Polynomial.aeval z cubic19Int = 0) :
    ∃ x : ℚ_[p], Polynomial.aeval x p5 = 0 := by
  refine ⟨(z : ℚ_[p]), ?_⟩
  have hzExpr : z ^ 3 - 19 = 0 := by
    simpa [cubic19Int, Polynomial.aeval_def] using hz
  have hzExpr' : (z : ℚ_[p]) ^ 3 - 19 = 0 := by
    have hmap := congrArg (algebraMap ℤ_[p] ℚ_[p]) hzExpr
    simp only [map_sub, map_pow, map_zero, PadicInt.algebraMap_apply] at hmap
    have h19 : ((19 : ℤ_[p]) : ℚ_[p]) = (19 : ℚ_[p]) :=
      PadicInt.coe_natCast 19
    rw [h19] at hmap
    exact hmap
  simp [p5, Polynomial.aeval_def, hzExpr']

private theorem p5_has_padic_root_of_quadraticThreeInt_root
    (p : ℕ) [Fact p.Prime] (z : ℤ_[p])
    (hz : Polynomial.aeval z quadraticThreeInt = 0) :
    ∃ x : ℚ_[p], Polynomial.aeval x p5 = 0 := by
  refine ⟨(z : ℚ_[p]), ?_⟩
  have hzExpr : z ^ 2 + z + 1 = 0 := by
    simpa [quadraticThreeInt, Polynomial.aeval_def] using hz
  have hzExpr' : (z : ℚ_[p]) ^ 2 + (z : ℚ_[p]) + 1 = 0 := by
    have hmap := congrArg (algebraMap ℤ_[p] ℚ_[p]) hzExpr
    simp only [map_add, map_pow, map_one, map_zero,
      PadicInt.algebraMap_apply] at hmap
    exact hmap
  simp [p5, Polynomial.aeval_def, hzExpr']

/-- The direct exceptional witness at `p = 2` gives a root of the displayed
quintic over `ℚ_[2]`. -/
theorem p5_has_twoAdic_root :
    ∃ z : ℚ_[2], Polynomial.aeval z p5 = 0 := by
  obtain ⟨z, hz⟩ := cubic19_has_twoAdicInt_root
  exact p5_has_padic_root_of_cubic19Int_root 2 z hz

/-- The direct exceptional witness at `p = 19` gives a root of the displayed
quintic over `ℚ_[19]`. -/
theorem p5_has_nineteenAdic_root :
    ∃ z : ℚ_[19], Polynomial.aeval z p5 = 0 := by
  obtain ⟨z, hz⟩ := quadraticThree_has_nineteenAdicInt_root
  exact p5_has_padic_root_of_quadraticThreeInt_root 19 z hz

/-- The Berend--Bilu quintic has a root over the `p`-adic field for every
prime `p`.  The proof follows the paper's five-row table exactly. -/
theorem p5_has_padic_root (p : ℕ) [Fact p.Prime] :
    ∃ z : ℚ_[p], Polynomial.aeval z p5 = 0 := by
  have hp : p.Prime := Fact.out
  by_cases hp2 : p = 2
  · subst p
    exact p5_has_twoAdic_root
  by_cases hp3 : p = 3
  · subst p
    exact p5_has_threeAdic_root
  by_cases hp19 : p = 19
  · subst p
    exact p5_has_nineteenAdic_root
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
      quadraticThree_has_padicInt_root_of_mod_three_eq_one p hmod1
    exact p5_has_padic_root_of_quadraticThreeInt_root p z hz
  · obtain ⟨z, hz⟩ :=
      cubic19_has_padicInt_root_of_mod_three_eq_two p hp19 hmod2
    exact p5_has_padic_root_of_cubic19Int_root p z hz

/-- The literal displayed Keller fiber has a point over `ℚ_[p]` for every
prime `p`. -/
theorem integralFiberPoint_padic_nonempty (p : ℕ) [Fact p.Prime] :
    Nonempty (IntegralFiberPoint ℚ_[p]) := by
  obtain ⟨z, hz⟩ := p5_has_padic_root p
  let root : PolynomialRoot p5 ℚ_[p] := ⟨z, hz⟩
  exact ⟨integralFiberRepresentingEquiv
    ((PolynomialRoot.algHomEquiv p5 ℚ_[p]).symm root)⟩

/-- End-to-end Hasse certificate for the explicit literal fiber: it has no
rational point, has a real point, and has a point over every `p`-adic
completion. -/
theorem integralFiberPoint_hasse_certificate :
    IsEmpty (IntegralFiberPoint ℚ) ∧
      Nonempty (IntegralFiberPoint ℝ) ∧
        ∀ (p : ℕ) [Fact p.Prime], Nonempty (IntegralFiberPoint ℚ_[p]) :=
  ⟨integralFiberPoint_rat_isEmpty, integralFiberPoint_real_nonempty,
    integralFiberPoint_padic_nonempty⟩

#print axioms p5_has_twoAdic_root
#print axioms p5_has_nineteenAdic_root
#print axioms p5_has_padic_root
#print axioms integralFiberPoint_padic_nonempty
#print axioms integralFiberPoint_hasse_certificate

end FiniteEtaleKeller.ExplicitQuintic
