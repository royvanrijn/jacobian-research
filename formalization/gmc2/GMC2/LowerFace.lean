/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GMC2.ConstantTerm
import GMC2.DuistermaatVanDerKallen
import GMC2.GaussianModel
import GMC2.OneSided
import GMC2.PrimeIsolation
import GMC2.Specialization
import GMC2.SupportingFace

/-!
# The lower-face route to GMC(2)

This file records the proof's dependency graph and exposes the final
contradiction as a small theorem.  The geometric construction supplies a
lowest radial polynomial `F_r`; specialization supplies a good prime; the
theorem below is then exactly the prime-isolation paragraph of the paper.
-/

namespace GMC2

open Polynomial

/-- Once the lower-face polynomial and a good characteristic-`p`
specialization have been constructed, vanishing of the normalized Gaussian
moment contradicts survival of its exposed coefficient. -/
theorem lowerFace_prime_contradiction
    {K : Type*} [Field K] {p n : ℕ} [CharP K p]
    (hp : p.Prime) (f : K[X])
    (hlow : ∀ j < n, f.coeff j = 0)
    (hc : f.coeff n ≠ 0)
    (hmoment : primeDilatedMoment p n f = 0) :
    False :=
  hc (lowestCoeff_eq_zero_of_primeDilatedMoment_eq_zero hp f hlow hmoment)

/-- Strictly one-sided integral weights eventually miss every fixed bounded
set of weights.  This is the final, purely order-theoretic GMC step. -/
theorem eventually_no_zero_weight
    (S Q : Finset ℤ) (hS : ∀ k ∈ S, 0 < k) :
    ∃ M : ℕ, ∀ m ≥ M, ∀ k ∈ S, ∀ q ∈ Q, m * k + q ≠ 0 := by
  classical
  let B : ℕ := Q.sup (fun q ↦ q.natAbs)
  refine ⟨B + 1, ?_⟩
  intro m hm k hk q hq hzero
  have hk1 : (1 : ℤ) ≤ k := hS k hk
  have hmk : (m : ℤ) ≤ (m : ℤ) * k := by
    nlinarith
  have hqB : q.natAbs ≤ B := Finset.le_sup hq
  have hqneg : -(B : ℤ) ≤ q := by
    have habsZ : (q.natAbs : ℤ) ≤ (B : ℤ) := by exact_mod_cast hqB
    have hqbound : -(q.natAbs : ℤ) ≤ q := by
      rcases le_total 0 q with h | h
      · exact le_trans (neg_nonpos.mpr (Int.natCast_nonneg _)) h
      · exact (Int.eq_neg_natAbs_of_nonpos h).ge
    exact le_trans (neg_le_neg habsZ) hqbound
  have hmB : (B : ℤ) < m := by exact_mod_cast (lt_of_lt_of_le (Nat.lt_succ_self B) hm)
  nlinarith

/-- The negative-weight counterpart of `eventually_no_zero_weight`. -/
theorem eventually_no_zero_weight_neg
    (S Q : Finset ℤ) (hS : ∀ k ∈ S, k < 0) :
    ∃ M : ℕ, ∀ m ≥ M, ∀ k ∈ S, ∀ q ∈ Q, m * k + q ≠ 0 := by
  classical
  let B : ℕ := Q.sup (fun q ↦ q.natAbs)
  refine ⟨B + 1, ?_⟩
  intro m hm k hk q hq hzero
  have hkneg : k ≤ (-1 : ℤ) := by
    have := hS k hk
    omega
  have hmk : (m : ℤ) * k ≤ -(m : ℤ) := by
    have hmnonneg : (0 : ℤ) ≤ m := Int.natCast_nonneg m
    nlinarith
  have hqB : q.natAbs ≤ B := Finset.le_sup hq
  have habsZ : (q.natAbs : ℤ) ≤ (B : ℤ) := by exact_mod_cast hqB
  have hqle : q ≤ (B : ℤ) := le_trans Int.le_natAbs habsZ
  have hmB : (B : ℤ) < m := by exact_mod_cast (lt_of_lt_of_le (Nat.lt_succ_self B) hm)
  nlinarith

/-- Either strict side of the angular support gives the eventual
constant-term vanishing required by GMC(2). -/
theorem eventually_no_zero_weight_of_oneSided
    (S Q : Finset ℤ)
    (hS : (∀ k ∈ S, 0 < k) ∨ (∀ k ∈ S, k < 0)) :
    ∃ M : ℕ, ∀ m ≥ M, ∀ k ∈ S, ∀ q ∈ Q, m * k + q ≠ 0 := by
  rcases hS with hpos | hneg
  · exact eventually_no_zero_weight S Q hpos
  · exact eventually_no_zero_weight_neg S Q hneg

end GMC2
