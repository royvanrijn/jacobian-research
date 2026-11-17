/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.NumberTheory.Padics.PadicVal.Basic
import Mathlib.Data.Nat.Prime.Factorial

/-!
# Factorial valuations in the prime-dilation argument

This file proves the rational-prime arithmetic core of the paper's
factorial-valuation lemma.  The transfer from the rational `p`-adic
valuation to a valuation at a prime ideal of an unramified number field is
kept outside this module.
-/

namespace GVC

/-- The integral quotient `j! / n!`, used only with `n ≤ j`. -/
def factorialQuotient (n j : ℕ) : ℕ :=
  Nat.factorial j / Nat.factorial n

theorem factorialQuotient_mul {n j : ℕ} (hnj : n ≤ j) :
    Nat.factorial n * factorialQuotient n j = Nat.factorial j := by
  exact Nat.mul_div_cancel' (Nat.factorial_dvd_factorial hnj)

theorem factorialQuotient_ne_zero {n j : ℕ} (hnj : n ≤ j) :
    factorialQuotient n j ≠ 0 := by
  intro hzero
  have h := factorialQuotient_mul hnj
  rw [hzero, mul_zero] at h
  exact Nat.factorial_ne_zero j h.symm

private theorem padicValNat_le_of_dvd
    {p a b : ℕ} [Fact p.Prime]
    (ha : a ≠ 0) (hb : b ≠ 0) (hab : a ∣ b) :
    padicValNat p a ≤ padicValNat p b := by
  obtain ⟨k, rfl⟩ := hab
  have hk : k ≠ 0 := by
    intro hk
    simp [hk] at hb
  rw [padicValNat.mul ha hk]
  omega

private theorem padicValNat_factorial_eq_zero_of_lt
    {p n : ℕ} [hp : Fact p.Prime] (hn : n < p) :
    padicValNat p (Nat.factorial n) = 0 := by
  have hpprime : p.Prime := Fact.out
  apply padicValNat.eq_zero_of_not_dvd
  rw [Nat.Prime.dvd_factorial hpprime, not_le]
  exact hn

/-- One-coordinate lower bound from the factorial-valuation lemma. -/
theorem padicVal_factorialQuotient_lower
    {p : ℕ} [hp : Fact p.Prime] (gamma eta : ℕ) :
    gamma / p ≤
      padicValNat p (factorialQuotient (p * eta) (gamma + p * eta)) := by
  have hdenom_le : p * eta ≤ gamma + p * eta := by omega
  have htop :
      padicValNat p (Nat.factorial (gamma + p * eta)) =
        padicValNat p (Nat.factorial (gamma / p + eta)) +
          (gamma / p + eta) := by
    calc
      padicValNat p (Nat.factorial (gamma + p * eta)) =
          padicValNat p
            (Nat.factorial (p * ((gamma + p * eta) / p))) := by
              rw [padicValNat_mul_div_factorial]
      _ = padicValNat p (Nat.factorial ((gamma + p * eta) / p)) +
          ((gamma + p * eta) / p) :=
            padicValNat_factorial_mul _
      _ = padicValNat p (Nat.factorial (gamma / p + eta)) +
          (gamma / p + eta) := by
            have hp_pos : 0 < p := (show p.Prime from Fact.out).pos
            rw [Nat.add_mul_div_left gamma eta hp_pos]
  have hdenom :
      padicValNat p (Nat.factorial (p * eta)) =
        padicValNat p (Nat.factorial eta) + eta :=
    padicValNat_factorial_mul eta
  have heta_le : eta ≤ gamma / p + eta := Nat.le_add_left _ _
  have hvaluation_le :
      padicValNat p (Nat.factorial eta) ≤
        padicValNat p (Nat.factorial (gamma / p + eta)) :=
    padicValNat_le_of_dvd
      (Nat.factorial_ne_zero eta)
      (Nat.factorial_ne_zero (gamma / p + eta))
      (Nat.factorial_dvd_factorial heta_le)
  rw [factorialQuotient,
    padicValNat.div_of_dvd (Nat.factorial_dvd_factorial hdenom_le),
    htop, hdenom]
  apply Nat.le_sub_of_add_le
  omega

/-- One-coordinate exact valuation when a whole prime block is added and
the block remains below `p`. -/
theorem padicVal_factorialQuotient_prime_block
    {p : ℕ} [hp : Fact p.Prime] (beta eta : ℕ)
    (hsmall : beta + eta < p) :
    padicValNat p
        (factorialQuotient (p * eta) (p * (beta + eta))) = beta := by
  have heta : eta < p := by omega
  have hdenom_le : p * eta ≤ p * (beta + eta) := by
    exact Nat.mul_le_mul_left p (by omega)
  rw [factorialQuotient,
    padicValNat.div_of_dvd (Nat.factorial_dvd_factorial hdenom_le),
    padicValNat_factorial_mul, padicValNat_factorial_mul,
    padicValNat_factorial_eq_zero_of_lt hsmall,
    padicValNat_factorial_eq_zero_of_lt heta]
  omega

/-- The two-coordinate factorial product `F_p(γ,η)` from the paper. -/
def factorialValuationProduct
    (p gammaX gammaY etaX etaY : ℕ) : ℕ :=
  factorialQuotient (p * etaX) (gammaX + p * etaX) *
    factorialQuotient (p * etaY) (gammaY + p * etaY)

theorem padicVal_factorialValuationProduct_lower
    {p : ℕ} [hp : Fact p.Prime]
    (gammaX gammaY etaX etaY : ℕ) :
    gammaX / p + gammaY / p ≤
      padicValNat p
        (factorialValuationProduct p gammaX gammaY etaX etaY) := by
  have hxne :
      factorialQuotient (p * etaX) (gammaX + p * etaX) ≠ 0 :=
    factorialQuotient_ne_zero (by omega)
  have hyne :
      factorialQuotient (p * etaY) (gammaY + p * etaY) ≠ 0 :=
    factorialQuotient_ne_zero (by omega)
  rw [factorialValuationProduct, padicValNat.mul hxne hyne]
  exact Nat.add_le_add
    (padicVal_factorialQuotient_lower gammaX etaX)
    (padicVal_factorialQuotient_lower gammaY etaY)

/-- The floor estimate used for non-Frobenius exponents in the shifted-ray
proof.  It is the integral form of the observation that two remainders are
strictly smaller than `2p`. -/
theorem floor_sum_ge_sub_one
    {p gammaX gammaY s : ℕ} (hp : 0 < p)
    (htotal : p * s ≤ gammaX + gammaY) :
    s - 1 ≤ gammaX / p + gammaY / p := by
  have hx := Nat.div_add_mod gammaX p
  have hy := Nat.div_add_mod gammaY p
  have hrx := Nat.mod_lt gammaX hp
  have hry := Nat.mod_lt gammaY hp
  have hxlt : gammaX < p * (gammaX / p) + p := by
    calc
      gammaX = p * (gammaX / p) + gammaX % p := hx.symm
      _ < p * (gammaX / p) + p := Nat.add_lt_add_left hrx _
  have hylt : gammaY < p * (gammaY / p) + p := by
    calc
      gammaY = p * (gammaY / p) + gammaY % p := hy.symm
      _ < p * (gammaY / p) + p := Nat.add_lt_add_left hry _
  have hupper :
      gammaX + gammaY < p * (gammaX / p + gammaY / p + 2) := by
    calc
      gammaX + gammaY <
          (p * (gammaX / p) + p) + (p * (gammaY / p) + p) :=
        Nat.add_lt_add hxlt hylt
      _ = p * (gammaX / p + gammaY / p + 2) := by ring
  by_contra hbound
  have hs : 0 < s := by
    by_contra hs0
    have : s = 0 := Nat.eq_zero_of_not_pos hs0
    simp [this] at hbound
  have hq : gammaX / p + gammaY / p + 2 ≤ s := by omega
  have hmul :
      p * (gammaX / p + gammaY / p + 2) ≤ p * s :=
    Nat.mul_le_mul_left p hq
  omega

theorem padicVal_factorialValuationProduct_ge_sub_one_of_total
    {p : ℕ} [hp : Fact p.Prime]
    (gammaX gammaY etaX etaY s : ℕ)
    (htotal : p * s ≤ gammaX + gammaY) :
    s - 1 ≤
      padicValNat p
        (factorialValuationProduct p gammaX gammaY etaX etaY) :=
  (floor_sum_ge_sub_one (show 0 < p from hp.out.pos) htotal).trans
    (padicVal_factorialValuationProduct_lower gammaX gammaY etaX etaY)

theorem padicVal_factorialValuationProduct_exact
    {p : ℕ} [hp : Fact p.Prime]
    (betaX betaY etaX etaY : ℕ)
    (hxsmall : betaX + etaX < p)
    (hysmall : betaY + etaY < p) :
    padicValNat p
        (factorialValuationProduct p (p * betaX) (p * betaY) etaX etaY) =
      betaX + betaY := by
  have hxne :
      factorialQuotient (p * etaX) (p * betaX + p * etaX) ≠ 0 :=
    factorialQuotient_ne_zero (by omega)
  have hyne :
      factorialQuotient (p * etaY) (p * betaY + p * etaY) ≠ 0 :=
    factorialQuotient_ne_zero (by omega)
  rw [factorialValuationProduct, padicValNat.mul hxne hyne]
  simpa [Nat.mul_add] using congrArg₂ (fun a b : ℕ ↦ a + b)
    (padicVal_factorialQuotient_prime_block betaX etaX hxsmall)
    (padicVal_factorialQuotient_prime_block betaY etaY hysmall)

end GVC
