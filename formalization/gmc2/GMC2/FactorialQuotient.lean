import Mathlib.Data.Nat.Factorial.Basic
import Mathlib.Data.Nat.Prime.Basic

/-!
# Factorial quotients

The elementary divisibility step in the lower-face proof.  If `j > n`,
then `(jp)! / (np)!` contains the factor `(n+1)p`.
-/

namespace GMC2

/-- The integral quotient `j! / n!`.  It is only used when `n ≤ j`. -/
def factorialQuotient (n j : ℕ) : ℕ := j ! / n !

theorem factorialQuotient_mul (h : n ≤ j) :
    n ! * factorialQuotient n j = j ! := by
  exact Nat.mul_div_cancel' (Nat.factorial_dvd_factorial h)

theorem factor_dvd_factorialQuotient
    {a n j : ℕ} (hnj : n ≤ j) (hna : n < a) (haj : a ≤ j) :
    a ∣ factorialQuotient n j := by
  rw [factorialQuotient, Nat.dvd_div_iff_mul_dvd
    (Nat.factorial_dvd_factorial hnj)]
  exact (Nat.mul_dvd_mul_left a (Nat.factorial_dvd_factorial hna)).trans
    (Nat.factorial_dvd_factorial haj)

/-- The precise divisibility used after prime dilation. -/
theorem prime_dvd_factorialQuotient_mul
    {p n j : ℕ} (hp : 0 < p) (hnj : n < j) :
    p ∣ factorialQuotient (n * p) (j * p) := by
  have h₁ : n * p < (n + 1) * p := Nat.mul_lt_mul_of_pos_right (Nat.lt_succ_self n) hp
  have h₂ : (n + 1) * p ≤ j * p :=
    Nat.mul_le_mul_right p (Nat.succ_le_iff.mpr hnj)
  exact (dvd_mul_left p (n + 1)).trans
    (factor_dvd_factorialQuotient h₁.le h₁ h₂)

end GMC2

