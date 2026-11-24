import Mathlib.NumberTheory.ArithmeticFunction.Misc
import Mathlib.Analysis.Real.Sqrt

open scoped BigOperators ArithmeticFunction.sigma

namespace QuadraticNorms

/-- The mean of paired positive divisors d and n/d bounds √n. -/
theorem sigma_zero_mul_sqrt_le_sigma_one (n : ℕ) :
    (σ 0 n : ℝ) * Real.sqrt n ≤ σ 1 n := by
  have hpoint (d : ℕ) (hd : d ∈ n.divisors) :
      2 * Real.sqrt n ≤ (d : ℝ) + ((n / d : ℕ) : ℝ) := by
    have hmul : (d : ℝ) * ((n / d : ℕ) : ℝ) = n := by
      exact_mod_cast Nat.mul_div_cancel' (Nat.mem_divisors.mp hd).1
    have hsq : Real.sqrt (n : ℝ) ^ 2 = n := Real.sq_sqrt (by positivity)
    apply (sq_le_sq₀ (by positivity) (by positivity)).mp
    nlinarith only [hsq, hmul, sq_nonneg ((d : ℝ) - ((n / d : ℕ) : ℝ))]
  have hsum := Finset.sum_le_sum (s := n.divisors) hpoint
  rw [Finset.sum_const, nsmul_eq_mul, Finset.sum_add_distrib,
    Nat.sum_div_divisors n (fun d => (d : ℝ))] at hsum
  rw [ArithmeticFunction.sigma_zero_apply, ArithmeticFunction.sigma_one_apply]
  push_cast
  linarith only [hsum]

/-- Scaling divisors by a embeds them in the divisors of a*b. -/
theorem pow_mul_sigma_le {a b k : ℕ} (ha : 0 < a) (hb : 0 < b) :
    a ^ k * σ k b ≤ σ k (a * b) := by
  let f : ℕ ↪ ℕ := ⟨fun d => a * d, fun _ _ h => Nat.eq_of_mul_eq_mul_left ha h⟩
  have hsub : b.divisors.map f ⊆ (a * b).divisors := by
    intro d hd
    obtain ⟨e, he, rfl⟩ := Finset.mem_map.mp hd
    exact Nat.mem_divisors.mpr ⟨Nat.mul_dvd_mul_left a (Nat.mem_divisors.mp he).1, (Nat.mul_pos ha hb).ne'⟩
  simp only [ArithmeticFunction.sigma_apply]
  calc
    _ = ∑ d ∈ b.divisors.map f, d ^ k := by
      rw [Finset.sum_map, Finset.mul_sum]
      exact Finset.sum_congr rfl (fun _ _ => (mul_pow ..).symm)
    _ ≤ _ := Finset.sum_le_sum_of_subset_of_nonneg hsub (fun _ _ _ => Nat.zero_le _)

end QuadraticNorms
