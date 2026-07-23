import GMC2.GaussianModel
import GMC2.DuistermaatVanDerKallen
import Mathlib.Algebra.MonoidAlgebra.Degree
import Mathlib.Algebra.MonoidAlgebra.Support

/-!
# One-sided angular support

This file proves the final support argument of GMC(2): if every angular
weight of `P` has one strict sign, then every fixed mixed moment eventually
vanishes.
-/

namespace GMC2

open LaurentPolynomial

noncomputable def lowerWeight
    {R : Type*} [CommRing R] (P : CircularPolynomial R) : WithTop ℤ :=
  P.coeff.support.inf fun k ↦ (k : WithTop ℤ)

theorem lowerWeight_le_of_mem
    {R : Type*} [CommRing R] {P : CircularPolynomial R} {k : ℤ}
    (hk : k ∈ angularSupport P) :
    lowerWeight P ≤ (k : WithTop ℤ) :=
  Finset.inf_le hk

theorem one_le_lowerWeight
    {R : Type*} [CommRing R] {P : CircularPolynomial R}
    (hP : PositiveAngularSupport P) :
    (1 : WithTop ℤ) ≤ lowerWeight P := by
  apply Finset.le_inf
  intro k hk
  exact_mod_cast hP k hk

theorem lowerWeight_pow
    {R : Type*} [CommRing R] (P : CircularPolynomial R) (m : ℕ) :
    m • lowerWeight P ≤ lowerWeight (P ^ m) := by
  exact AddMonoidAlgebra.le_inf_support_pow (by simp) (by simp) m P

theorem pow_weight_ge
    {R : Type*} [CommRing R] {P : CircularPolynomial R}
    (hP : PositiveAngularSupport P) {m : ℕ} {k : ℤ}
    (hk : k ∈ angularSupport (P ^ m)) :
    (m : ℤ) ≤ k := by
  have h₁ := nsmul_le_nsmul_right (one_le_lowerWeight hP) m
  have h₂ := lowerWeight_pow P m
  have h₃ := lowerWeight_le_of_mem hk
  have h : (m : WithTop ℤ) ≤ (k : WithTop ℤ) := by
    simpa using h₁.trans (h₂.trans h₃)
  exact WithTop.coe_le_coe.mp h

theorem constantTerm_mul_pow_eq_zero_of_positive
    {R : Type*} [CommRing R] (Q P : CircularPolynomial R)
    (hP : PositiveAngularSupport P) :
    ∃ M : ℕ, ∀ m ≥ M, constantTerm (Q * P ^ m) = 0 := by
  classical
  let B : ℕ := (angularSupport Q).sup fun q ↦ q.natAbs
  refine ⟨B + 1, ?_⟩
  intro m hm
  apply Finsupp.notMem_support_iff.mp
  intro hzero
  have hsub := AddMonoidAlgebra.support_coeff_mul_subset Q (P ^ m) hzero
  obtain ⟨q, hq, k, hk, hqk⟩ := Finset.mem_add.mp hsub
  have hkge : (m : ℤ) ≤ k := pow_weight_ge hP hk
  have hqB : q.natAbs ≤ B := Finset.le_sup hq
  have hqneg : -(B : ℤ) ≤ q := by
    have habsZ : (q.natAbs : ℤ) ≤ (B : ℤ) := by exact_mod_cast hqB
    have hqbound : -(q.natAbs : ℤ) ≤ q := by
      rcases le_total 0 q with h | h
      · exact le_trans (neg_nonpos.mpr (Int.natCast_nonneg _)) h
      · exact (Int.eq_neg_natAbs_of_nonpos h).ge
    exact le_trans (neg_le_neg habsZ) hqbound
  have hmB : (B : ℤ) < m := by
    exact_mod_cast (lt_of_lt_of_le (Nat.lt_succ_self B) hm)
  omega

theorem eventuallyMixedMomentsVanish_of_positive
    {R : Type*} [CommRing R] (P : CircularPolynomial R)
    (hP : PositiveAngularSupport P) :
    EventuallyMixedMomentsVanish P := by
  intro Q
  obtain ⟨M, hM⟩ := constantTerm_mul_pow_eq_zero_of_positive Q P hP
  refine ⟨M, fun m hm ↦ ?_⟩
  simp [mixedMoment, circularExpectation, hM m hm, factorialFunctional]

noncomputable def upperWeight
    {R : Type*} [CommRing R] (P : CircularPolynomial R) : WithBot ℤ :=
  P.coeff.support.sup fun k ↦ (k : WithBot ℤ)

theorem upperWeight_ge_of_mem
    {R : Type*} [CommRing R] {P : CircularPolynomial R} {k : ℤ}
    (hk : k ∈ angularSupport P) :
    (k : WithBot ℤ) ≤ upperWeight P :=
  Finset.le_sup hk

theorem upperWeight_le_neg_one
    {R : Type*} [CommRing R] {P : CircularPolynomial R}
    (hP : NegativeAngularSupport P) :
    upperWeight P ≤ ((-1 : ℤ) : WithBot ℤ) := by
  apply Finset.sup_le
  intro k hk
  have hk' : k ≤ (-1 : ℤ) := by
    have := hP k hk
    omega
  exact_mod_cast hk'

theorem upperWeight_pow
    {R : Type*} [CommRing R] (P : CircularPolynomial R) (m : ℕ) :
    upperWeight (P ^ m) ≤ m • upperWeight P := by
  exact AddMonoidAlgebra.sup_support_pow_le (by simp) (by simp) m P

theorem pow_weight_le_neg
    {R : Type*} [CommRing R] {P : CircularPolynomial R}
    (hP : NegativeAngularSupport P) {m : ℕ} {k : ℤ}
    (hk : k ∈ angularSupport (P ^ m)) :
    k ≤ -(m : ℤ) := by
  have h₁ := nsmul_le_nsmul_right (upperWeight_le_neg_one hP) m
  have h₂ := upperWeight_pow P m
  have h₃ := upperWeight_ge_of_mem hk
  have h : (k : WithBot ℤ) ≤ ((-(m : ℤ) : ℤ) : WithBot ℤ) := by
    calc
      (k : WithBot ℤ) ≤ m • ((-1 : ℤ) : WithBot ℤ) := h₃.trans (h₂.trans h₁)
      _ = ((-(m : ℤ) : ℤ) : WithBot ℤ) := by
        rw [← WithBot.coe_nsmul]
        congr 1
        simp
  exact WithBot.coe_le_coe.mp h

theorem constantTerm_mul_pow_eq_zero_of_negative
    {R : Type*} [CommRing R] (Q P : CircularPolynomial R)
    (hP : NegativeAngularSupport P) :
    ∃ M : ℕ, ∀ m ≥ M, constantTerm (Q * P ^ m) = 0 := by
  classical
  let B : ℕ := (angularSupport Q).sup fun q ↦ q.natAbs
  refine ⟨B + 1, ?_⟩
  intro m hm
  apply Finsupp.notMem_support_iff.mp
  intro hzero
  have hsub := AddMonoidAlgebra.support_coeff_mul_subset Q (P ^ m) hzero
  obtain ⟨q, hq, k, hk, hqk⟩ := Finset.mem_add.mp hsub
  have hkle : k ≤ -(m : ℤ) := pow_weight_le_neg hP hk
  have hqB : q.natAbs ≤ B := Finset.le_sup hq
  have hqle : q ≤ (B : ℤ) := by
    have habsZ : (q.natAbs : ℤ) ≤ (B : ℤ) := by exact_mod_cast hqB
    exact le_trans Int.le_natAbs habsZ
  have hmB : (B : ℤ) < m := by
    exact_mod_cast (lt_of_lt_of_le (Nat.lt_succ_self B) hm)
  omega

theorem eventuallyMixedMomentsVanish_of_negative
    {R : Type*} [CommRing R] (P : CircularPolynomial R)
    (hP : NegativeAngularSupport P) :
    EventuallyMixedMomentsVanish P := by
  intro Q
  obtain ⟨M, hM⟩ := constantTerm_mul_pow_eq_zero_of_negative Q P hP
  refine ⟨M, fun m hm ↦ ?_⟩
  simp [mixedMoment, circularExpectation, hM m hm, factorialFunctional]

theorem eventuallyMixedMomentsVanish_of_oneSided
    {R : Type*} [CommRing R] (P : CircularPolynomial R)
    (hP : OneSidedAngularSupport P) :
    EventuallyMixedMomentsVanish P := by
  rcases hP with hpos | hneg
  · exact eventuallyMixedMomentsVanish_of_positive P hpos
  · exact eventuallyMixedMomentsVanish_of_negative P hneg

theorem oneSidedAngularSupport_of_not_straddles
    {R : Type*} [CommRing R] (P : CircularPolynomial R)
    (hP : ¬ SupportStraddlesZero P) :
    OneSidedAngularSupport P := by
  classical
  rw [SupportStraddlesZero, not_and_or] at hP
  rcases hP with hnonpos | hnonneg
  · left
    intro k hk
    exact lt_of_not_ge fun hk0 ↦ hnonpos ⟨k, hk, hk0⟩
  · right
    intro k hk
    exact lt_of_not_ge fun hk0 ↦ hnonneg ⟨k, hk, hk0⟩

theorem eventuallyMixedMomentsVanish_of_not_straddles
    {R : Type*} [CommRing R] (P : CircularPolynomial R)
    (hP : ¬ SupportStraddlesZero P) :
    EventuallyMixedMomentsVanish P :=
  eventuallyMixedMomentsVanish_of_oneSided P
    (oneSidedAngularSupport_of_not_straddles P hP)

end GMC2
