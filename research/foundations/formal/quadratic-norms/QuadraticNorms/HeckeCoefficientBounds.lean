import Mathlib.NumberTheory.ArithmeticFunction.Misc
import Mathlib.Analysis.Complex.Norm
import Mathlib.Analysis.Complex.Basic
import Mathlib.Analysis.Real.Sqrt
import Mathlib.Tactic.LinearCombination

open scoped BigOperators ArithmeticFunction.sigma

namespace QuadraticNorms

/-- Normalization, multiplicativity, and the weight-two recurrence at
the primes where the coefficient bound is available. -/
structure WeightTwoHeckeCoefficients (a : ℕ → ℝ) (allowed : ℕ → Prop) : Prop where
  normalized : a 1 = 1
  coprime_mul : ∀ {m n : ℕ}, m.Coprime n → a (m * n) = a m * a n
  prime_bound : ∀ p : ℕ, p.Prime → allowed p → |a p| ≤ 2 * Real.sqrt p
  prime_power : ∀ p : ℕ, p.Prime → allowed p → ∀ k : ℕ,
    a (p ^ (k + 2)) = a p * a (p ^ (k + 1)) - (p : ℝ) * a (p ^ k)

def heckePowerSum (u v : ℂ) (k : ℕ) : ℂ :=
  ∑ j ∈ Finset.range (k + 1), u ^ j * v ^ (k - j)

private theorem hecke_power_sum_succ (u v : ℂ) (k : ℕ) :
    heckePowerSum u v (k + 1) = v * heckePowerSum u v k + u ^ (k + 1) := by
  unfold heckePowerSum
  rw [Finset.sum_range_succ, Finset.mul_sum, Nat.sub_self, pow_zero, mul_one]
  congr 1
  refine Finset.sum_congr rfl (fun j hj => ?_)
  have hj' := Finset.mem_range.mp hj
  rw [show k + 1 - j = k - j + 1 by omega, pow_succ]
  ring

private theorem hecke_power_sum_succ' (u v : ℂ) (k : ℕ) :
    heckePowerSum u v (k + 1) = u * heckePowerSum u v k + v ^ (k + 1) := by
  unfold heckePowerSum
  rw [Finset.sum_range_succ', Finset.mul_sum, Nat.sub_zero, pow_zero, one_mul]
  congr 1
  refine Finset.sum_congr rfl (fun j _ => ?_)
  rw [Nat.add_sub_add_right, pow_succ]
  ring

private theorem hecke_power_sum_recurrence (u v : ℂ) (k : ℕ) :
    heckePowerSum u v (k + 2) =
      (u + v) * heckePowerSum u v (k + 1) - (u * v) * heckePowerSum u v k := by
  rw [hecke_power_sum_succ u v (k + 1)]
  have h := hecke_power_sum_succ u v k
  rw [pow_succ u (k + 1)]
  linear_combination -u * h

theorem hecke_power_sum_two_endpoints (u v : ℂ) (k : ℕ) :
    heckePowerSum u v (k + 2) =
      u ^ (k + 2) + v ^ (k + 2) + (u * v) * heckePowerSum u v k := by
  rw [hecke_power_sum_succ u v (k + 1), hecke_power_sum_succ' u v k, pow_succ v (k + 1)]
  ring

theorem hecke_prime_power_roots {a : ℕ → ℝ} {allowed : ℕ → Prop}
    (ha : WeightTwoHeckeCoefficients a allowed) {p : ℕ} (hp : p.Prime) (haP : allowed p) :
    ∃ u v : ℂ, ‖u‖ = Real.sqrt p ∧ ‖v‖ = Real.sqrt p ∧ u * v = (p : ℂ) ∧
      ∀ k : ℕ, (a (p ^ k) : ℂ) = heckePowerSum u v k := by
  let x := a p / 2
  let y := Real.sqrt ((p : ℝ) - x ^ 2)
  let u : ℂ := ⟨x, y⟩
  let v : ℂ := ⟨x, -y⟩
  have hr2 : Real.sqrt (p : ℝ) ^ 2 = p := Real.sq_sqrt (by positivity)
  have hx : x ^ 2 ≤ (p : ℝ) := by
    have h := (sq_le_sq₀ (abs_nonneg _) (by positivity)).mpr (ha.prime_bound p hp haP)
    rw [sq_abs, mul_pow, hr2] at h
    dsimp [x]
    nlinarith only [h]
  have hy2 : y ^ 2 = (p : ℝ) - x ^ 2 := Real.sq_sqrt (sub_nonneg.mpr hx)
  have hu : ‖u‖ = Real.sqrt p := by
    rw [Complex.norm_def]
    congr 1
    rw [Complex.normSq_apply]
    change x * x + y * y = p
    nlinarith only [hy2]
  have hv : ‖v‖ = Real.sqrt p := by
    rw [Complex.norm_def]
    congr 1
    rw [Complex.normSq_apply]
    change x * x + (-y) * (-y) = p
    nlinarith only [hy2]
  have huv : u * v = (p : ℂ) := by
    apply Complex.ext
    · rw [Complex.mul_re, Complex.natCast_re]
      change x * x - y * (-y) = p
      nlinarith only [hy2]
    · rw [Complex.mul_im, Complex.natCast_im]
      change x * (-y) + y * x = 0
      ring
  have hadd : u + v = (a p : ℂ) := by
    apply Complex.ext
    · change x + x = a p
      dsimp [x]
      ring
    · change y + (-y) = 0
      ring
  refine ⟨u, v, hu, hv, huv, fun k => ?_⟩
  induction k using Nat.twoStepInduction with
  | zero => simp [heckePowerSum, ha.normalized]
  | one =>
    rw [hecke_power_sum_succ]
    simpa [heckePowerSum, add_comm] using hadd.symm
  | more k ih₀ ih₁ =>
    rw [ha.prime_power p hp haP k, Complex.ofReal_sub, Complex.ofReal_mul, Complex.ofReal_mul,
      Complex.ofReal_natCast, ih₀, ih₁, hecke_power_sum_recurrence, hadd, huv]

theorem norm_hecke_power_sum_le {u v : ℂ} {r : ℝ}
    (hu : ‖u‖ = r) (hv : ‖v‖ = r) (k : ℕ) :
    ‖heckePowerSum u v k‖ ≤ ((k : ℝ) + 1) * r ^ k := by
  unfold heckePowerSum
  refine (norm_sum_le _ _).trans ?_
  calc
    _ = ∑ _j ∈ Finset.range (k + 1), r ^ k := by
      refine Finset.sum_congr rfl (fun j hj => ?_)
      rw [norm_mul, norm_pow, norm_pow, hu, hv, ← pow_add]
      congr 1
      have := Finset.mem_range.mp hj
      omega
    _ ≤ _ := by rw [Finset.sum_const, nsmul_eq_mul, Finset.card_range, Nat.cast_add, Nat.cast_one]

theorem sqrt_nat_power {x : ℝ} (hx : 0 ≤ x) (k : ℕ) :
    Real.sqrt (x ^ k) = (Real.sqrt x) ^ k := by
  induction k with
  | zero => simp
  | succ k ih => rw [pow_succ x, Real.sqrt_mul (pow_nonneg hx _), ih, pow_succ]

/-- The prime-power recurrence and coprime multiplicativity give the
divisor bound at every index supported on the allowed primes. -/
theorem weight_two_hecke_coefficient_bound {a : ℕ → ℝ} {allowed : ℕ → Prop}
    (ha : WeightTwoHeckeCoefficients a allowed)
    {m : ℕ} (hm : 0 < m) (hP : ∀ p, p.Prime → p ∣ m → allowed p) :
    |a m| ≤ (σ 0 m : ℝ) * Real.sqrt m := by
  induction m using Nat.recOnPosPrimePosCoprime with
  | prime_pow p k hp hk =>
    obtain ⟨u, v, hu, hv, _, hpow⟩ := hecke_prime_power_roots ha hp (hP p hp (dvd_pow_self p hk.ne'))
    rw [← Real.norm_eq_abs, ← Complex.norm_real, hpow, ArithmeticFunction.sigma_zero_apply_prime_pow hp,
      Nat.cast_add, Nat.cast_one, Nat.cast_pow, sqrt_nat_power (by positivity)]
    exact norm_hecke_power_sum_le hu hv k
  | zero => omega
  | one => simp [ha.normalized]
  | coprime x y hx hy hc ihx ihy =>
    have hxpos : 0 < x := by omega
    have hypos : 0 < y := by omega
    have hPx (p : ℕ) (hp : p.Prime) (hd : p ∣ x) := hP p hp (hd.trans (dvd_mul_right ..))
    have hPy (p : ℕ) (hp : p.Prime) (hd : p ∣ y) := hP p hp (hd.trans (dvd_mul_left ..))
    rw [ha.coprime_mul hc, abs_mul, ArithmeticFunction.isMultiplicative_sigma.map_mul_of_coprime hc,
      Nat.cast_mul, Nat.cast_mul, Real.sqrt_mul (by positivity)]
    calc
      _ ≤ ((σ 0 x : ℝ) * Real.sqrt x) * ((σ 0 y : ℝ) * Real.sqrt y) :=
        mul_le_mul (ihx hxpos hPx) (ihy hypos hPy) (abs_nonneg _) (by positivity)
      _ = _ := by ring

end QuadraticNorms
