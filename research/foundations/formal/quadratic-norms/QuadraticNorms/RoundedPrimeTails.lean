import Mathlib.Algebra.Order.Field.Rat
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.NormNum

namespace QuadraticNorms

def inverseSquareFixedPointScale : ℕ := 2 ^ 64

/-- The upward-rounded fixed-point weight for an inverse square. -/
def roundedInverseSquareWeight (p : ℕ) : ℕ :=
  (inverseSquareFixedPointScale + p ^ 2 - 1) / p ^ 2

theorem inverse_square_le_rounded {p : ℕ} (hp : 0 < p) :
    (1 : ℚ) / (p : ℚ) ^ 2 ≤ (roundedInverseSquareWeight p : ℕ) / (inverseSquareFixedPointScale : ℚ) := by
  unfold roundedInverseSquareWeight
  have hpow : 0 < p ^ 2 := pow_pos hp _
  have hscale : 0 < inverseSquareFixedPointScale := by decide
  have hdiv := Nat.div_add_mod (inverseSquareFixedPointScale + p ^ 2 - 1) (p ^ 2)
  have hmod := Nat.mod_lt (inverseSquareFixedPointScale + p ^ 2 - 1) hpow
  have hceil : inverseSquareFixedPointScale ≤ ((inverseSquareFixedPointScale + p ^ 2 - 1) / p ^ 2) * p ^ 2 := by
    rw [Nat.mul_comm] at hdiv
    omega
  have hpQ : (0 : ℚ) < p := by exact_mod_cast hp
  have hsQ : (0 : ℚ) < inverseSquareFixedPointScale := by exact_mod_cast hscale
  apply (div_le_div_iff₀ (pow_pos hpQ _) hsQ).mpr
  simpa only [one_mul, Nat.cast_mul, Nat.cast_pow] using
    (show (inverseSquareFixedPointScale : ℚ) ≤ (((inverseSquareFixedPointScale + p ^ 2 - 1) / p ^ 2 : ℕ) : ℚ) * (p : ℚ) ^ 2 from by
      exact_mod_cast hceil)

theorem inverse_square_sum_le_rounded (primes : List ℕ) (hpos : ∀ p ∈ primes, 0 < p) :
    (primes.map (fun (p : ℕ) => (1 : ℚ) / (p : ℚ) ^ 2)).sum ≤
      ((primes.map roundedInverseSquareWeight).sum : ℕ) / (inverseSquareFixedPointScale : ℚ) := by
  induction primes with
  | nil => simp
  | cons p primes ih =>
    have hp := inverse_square_le_rounded (hpos p (by simp))
    have hrest := ih (fun q hq => hpos q (by simp [hq]))
    simpa only [List.map_cons, List.sum_cons, Nat.cast_add, add_div] using add_le_add hp hrest

end QuadraticNorms
