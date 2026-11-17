import Mathlib

/-!
# Factorial quotients

The elementary divisibility step in the lower-face proof.  If `j > n`,
then `(jp)! / (np)!` contains the factor `(n+1)p`.
-/

namespace GMC2

/-- The integral quotient `j! / n!`.  It is only used when `n ≤ j`. -/
def factorialQuotient (n j : ℕ) : ℕ :=
  Nat.factorial j / Nat.factorial n

theorem factorialQuotient_mul (h : n ≤ j) :
    Nat.factorial n * factorialQuotient n j = Nat.factorial j := by
  exact Nat.mul_div_cancel' (Nat.factorial_dvd_factorial h)

theorem factor_dvd_ascFactorial
    {a start len : ℕ} (h₁ : start ≤ a) (h₂ : a < start + len) :
    a ∣ start.ascFactorial len := by
  induction len with
  | zero => omega
  | succ len ih =>
      rw [Nat.ascFactorial_succ]
      by_cases ha : a = start + len
      · exact ha ▸ dvd_mul_right _ _
      · exact dvd_mul_of_dvd_right (ih (by omega)) _

theorem factor_dvd_factorialQuotient
    {a n j : ℕ} (hnj : n ≤ j) (hna : n < a) (haj : a ≤ j) :
    a ∣ factorialQuotient n j := by
  have hj : n + (j - n) = j := Nat.add_sub_of_le hnj
  rw [factorialQuotient, ← hj, ← Nat.ascFactorial_eq_div]
  exact factor_dvd_ascFactorial (by omega) (by omega)

/-- The precise divisibility used after prime dilation. -/
theorem prime_dvd_factorialQuotient_mul
    {p n j : ℕ} (hp : 0 < p) (hnj : n < j) :
    p ∣ factorialQuotient (n * p) (j * p) := by
  have h₁ : n * p < (n + 1) * p := Nat.mul_lt_mul_of_pos_right (Nat.lt_succ_self n) hp
  have h₂ : (n + 1) * p ≤ j * p :=
    Nat.mul_le_mul_right p (Nat.succ_le_iff.mpr hnj)
  have hdiv : p ∣ (n + 1) * p := ⟨n + 1, Nat.mul_comm _ _⟩
  exact hdiv.trans
    (factor_dvd_factorialQuotient (h₁.le.trans h₂) h₁ h₂)

end GMC2
