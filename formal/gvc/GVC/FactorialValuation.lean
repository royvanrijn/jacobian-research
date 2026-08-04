/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib.NumberTheory.Padics.PadicVal.Basic

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
  apply padicValNat.eq_zero_of_not_dvd
  rw [hp.out.dvd_factorial, not_le]
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
            simp [Nat.mul_comm, Nat.add_mul_div_left]
  have hdenom :
      padicValNat p (Nat.factorial (p * eta)) =
        padicValNat p (Nat.factorial eta) + eta :=
    padicValNat_factorial_mul eta
  have heta_le : eta ≤ gamma / p + eta := by omega
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
  convert congrArg₂ (.+.)
    (padicVal_factorialQuotient_prime_block betaX etaX hxsmall)
    (padicVal_factorialQuotient_prime_block betaY etaY hysmall) using 1 <;>
    ring

end GVC
