/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.FixedHasseHeight

/-!
# The multiplicative arithmetic family for the fixed Hasse map

This module isolates the elementary arithmetic conditions behind the
counting family.  The congruence

`a ≡ 1 (mod 9)`

together with support on primes `1 (mod 3)` is closed under multiplication
and natural powers.  The noncube condition is deliberately kept separate:
it is needed for failure of the Hasse principle, but it is not itself
multiplicatively closed.

The module also formalizes the prime line.  A prime `ℓ ≡ 1 (mod 9)`
automatically satisfies the support condition and cannot be a cube in `ℚ`.
Consequently every such prime produces the complete paper parameter
certificate.
-/

noncomputable section

open Polynomial

namespace FiniteEtaleKeller.FixedHasseFamily

/-- The multiplicatively closed part of the paper's parameter conditions. -/
structure HasseCoreCondition (a : ℕ) : Prop where
  modNine : a % 9 = 1
  primeSupport :
    ∀ q : ℕ, q.Prime → q ∣ a → q % 3 = 1

namespace HasseCoreCondition

/-- The unit parameter satisfies the multiplicative core conditions. -/
theorem one : HasseCoreCondition 1 := by
  refine ⟨by norm_num, ?_⟩
  intro q hq hqone
  have : q = 1 := Nat.dvd_one.mp hqone
  exact (hq.ne_one this).elim

/-- The congruence and prime-support conditions are closed under
multiplication. -/
theorem mul {a b : ℕ}
    (ha : HasseCoreCondition a)
    (hb : HasseCoreCondition b) :
    HasseCoreCondition (a * b) := by
  refine ⟨?_, ?_⟩
  · rw [Nat.mul_mod, ha.modNine, hb.modNine]
  · intro q hq hqab
    rcases hq.dvd_mul.mp hqab with hqa | hqb
    · exact ha.primeSupport q hq hqa
    · exact hb.primeSupport q hq hqb

/-- The multiplicative core conditions are closed under natural powers. -/
theorem pow {a : ℕ}
    (ha : HasseCoreCondition a) (n : ℕ) :
    HasseCoreCondition (a ^ n) := by
  induction n with
  | zero =>
      simpa using one
  | succ n ih =>
      simpa [pow_succ] using mul ih ha

end HasseCoreCondition

/-- The exact nonanalytic parameter predicate used by the fixed-map
theorem. -/
structure AdmissibleHasseParameter (a : ℕ) : Prop
    extends HasseCoreCondition a where
  one_lt : 1 < a
  notRationalCube : ¬∃ r : ℚ, r ^ 3 = (a : ℚ)

private theorem cubicInt_monic (a : ℕ) :
    (cubicInt a).Monic := by
  simpa [cubicInt] using
    (monic_X_pow_sub_C (a : ℤ) (by norm_num : (3 : ℕ) ≠ 0))

/-- A natural prime is not a cube in the rational numbers. -/
theorem prime_not_rational_cube
    (ℓ : ℕ) (hprime : ℓ.Prime) :
    ¬∃ r : ℚ, r ^ 3 = (ℓ : ℚ) := by
  rintro ⟨r, hr⟩
  have hroot : Polynomial.aeval r (cubicInt ℓ) = 0 := by
    simpa [cubicInt, Polynomial.aeval_def] using sub_eq_zero.mpr hr
  obtain ⟨z, hz, -⟩ :=
    exists_integer_of_is_root_of_monic (cubicInt_monic ℓ) hroot
  subst r
  have hzpow : z ^ 3 = (ℓ : ℤ) := by
    have hcast : ((z ^ 3 : ℤ) : ℚ) = (((ℓ : ℕ) : ℤ) : ℚ) := by
      simpa using hr
    exact Int.cast_injective hcast
  have hznonneg : 0 ≤ z := by
    by_contra hneg
    have hzlt : z < 0 := lt_of_not_ge hneg
    have hzsqpos : 0 < z ^ 2 := sq_pos_of_ne_zero (ne_of_lt hzlt)
    have hzcube : z ^ 3 < 0 := by
      nlinarith
    have hℓpos : (0 : ℤ) < ℓ := by exact_mod_cast hprime.pos
    omega
  obtain ⟨n, rfl⟩ := Int.eq_ofNat_of_zero_le hznonneg
  have hnpow : n ^ 3 = ℓ := by exact_mod_cast hzpow
  exact Nat.Prime.not_prime_pow (x := n) (n := 3) (by norm_num)
    (hnpow ▸ hprime)

/-- A prime in the progression `1 (mod 9)` satisfies the multiplicative
core conditions. -/
theorem prime_hasseCoreCondition
    (ℓ : ℕ) (hprime : ℓ.Prime) (hmod9 : ℓ % 9 = 1) :
    HasseCoreCondition ℓ := by
  refine ⟨hmod9, ?_⟩
  intro q hq hqℓ
  rcases (Nat.dvd_prime hprime).mp hqℓ with hq1 | hqeq
  · exact (hq.ne_one hq1).elim
  · subst q
    omega

/-- Every prime `ℓ ≡ 1 (mod 9)` is an admissible Hasse parameter. -/
theorem prime_admissibleHasseParameter
    (ℓ : ℕ) (hprime : ℓ.Prime) (hmod9 : ℓ % 9 = 1) :
    AdmissibleHasseParameter ℓ :=
  { toHasseCoreCondition :=
      prime_hasseCoreCondition ℓ hprime hmod9
    one_lt := hprime.one_lt
    notRationalCube := prime_not_rational_cube ℓ hprime }

/-- Every admissible arithmetic parameter supplies the complete fixed-map
certificate. -/
theorem AdmissibleHasseParameter.certificate
    {a : ℕ} (ha : AdmissibleHasseParameter a) :
    PaperParameterCertificate a :=
  paperParameter_certificate a ha.one_lt ha.modNine
    ha.primeSupport ha.notRationalCube

/-- Paper-facing endpoint for the prime progression `ℓ ≡ 1 (mod 9)`. -/
theorem primeParameter_certificate
    (ℓ : ℕ) (hprime : ℓ.Prime) (hmod9 : ℓ % 9 = 1) :
    PaperParameterCertificate ℓ :=
  (prime_admissibleHasseParameter ℓ hprime hmod9).certificate

#print axioms HasseCoreCondition.mul
#print axioms HasseCoreCondition.pow
#print axioms prime_not_rational_cube
#print axioms prime_admissibleHasseParameter
#print axioms AdmissibleHasseParameter.certificate
#print axioms primeParameter_certificate

end FiniteEtaleKeller.FixedHasseFamily
