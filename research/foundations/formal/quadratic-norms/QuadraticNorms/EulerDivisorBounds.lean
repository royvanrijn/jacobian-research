import Mathlib.NumberTheory.ArithmeticFunction.Misc
import Mathlib.Data.Finset.Lattice.Fold
import Mathlib.Algebra.Order.Field.Rat
import Mathlib.Algebra.Order.Field.Basic
import Mathlib.Tactic.FieldSimp

open scoped BigOperators ArithmeticFunction.sigma

namespace QuadraticNorms

private def finiteFactorMax (factor : ℕ → ℚ) (last : ℕ) : ℚ :=
  (Finset.range (last + 1)).sup' ⟨0, by simp⟩ factor

theorem factor_le_finite_max (f : ℕ → ℚ) {last : ℕ}
    (hstep : ∀ a ≥ last, f (a + 1) ≤ f a) (a : ℕ) :
    f a ≤ finiteFactorMax f last := by
  induction a with
  | zero => exact Finset.le_sup' f (Finset.mem_range.mpr (by omega : 0 < last + 1))
  | succ a ih =>
    by_cases ha : a < last
    · exact Finset.le_sup' f (Finset.mem_range.mpr (by omega : a + 1 < last + 1))
    · exact (hstep a (Nat.le_of_not_gt ha)).trans ih

def boundedEulerCap (factor : ℕ → ℕ → ℚ) (last p : ℕ) : ℚ :=
  if p < 64 then finiteFactorMax (factor p) last else 1

theorem one_le_euler_cap (factor : ℕ → ℕ → ℚ) (last p : ℕ)
    (hzero : factor p 0 = 1) : 1 ≤ boundedEulerCap factor last p := by
  unfold boundedEulerCap
  split
  · rw [← hzero]
    exact Finset.le_sup' (factor p) (Finset.mem_range.mpr (by omega : 0 < last + 1))
  · exact le_rfl

theorem euler_product_le_small_primes (factor : ℕ → ℕ → ℚ) (last : ℕ)
    {support small : Finset ℕ} (exponent : ℕ → ℕ)
    (hzero : ∀ p, factor p 0 = 1)
    (hsmall : ∀ p ∈ support, p < 64 → p ∈ small)
    (hbound : ∀ p ∈ support, factor p (exponent p) ≤ boundedEulerCap factor last p)
    (hnonneg : ∀ p ∈ support, 0 ≤ factor p (exponent p)) :
    (∏ p ∈ support, factor p (exponent p)) ≤ ∏ p ∈ small, boundedEulerCap factor last p := by
  refine (Finset.prod_le_prod hnonneg hbound).trans ?_
  rw [← Finset.prod_filter_of_ne (p := (· < 64))
    (fun p _ hp => by by_contra! h; exact hp (if_neg (Nat.not_lt_of_ge h)))]
  exact Finset.prod_le_prod_of_subset_of_one_le
    (fun p hp => let ⟨hp', hp₆₄⟩ := Finset.mem_filter.mp hp; hsmall p hp' hp₆₄)
    (fun p _ => (one_le_euler_cap factor last p (hzero p)).trans' (by norm_num))
    (fun p _ _ => one_le_euler_cap factor last p (hzero p))

private theorem divisor_ratio_large {k p : ℕ} (hk : k ≤ 6) (hp : 64 ≤ p) (a : ℕ) :
    ((a : ℚ) + 1) ^ k ≤ (p : ℚ) ^ a := by
  have ha : (a : ℚ) + 1 ≤ (2 : ℚ) ^ a := by
    exact_mod_cast Nat.succ_le_of_lt (Nat.lt_two_pow_self (n := a))
  have hstep := pow_le_pow_left₀ (show 0 ≤ (a : ℚ) + 1 by positivity) ha k
  rw [← pow_mul, mul_comm a, pow_mul] at hstep
  have hbase : (2 : ℚ) ^ k ≤ p := (pow_le_pow_right₀ (by norm_num) hk).trans
    (by exact_mod_cast hp : (2 : ℚ) ^ 6 ≤ p)
  exact hstep.trans (pow_le_pow_left₀ (by positivity) hbase _)

private theorem next_exponent_bound {k p a : ℕ} {A : ℚ}
    (ha : (a : ℚ) + 2 ≤ A * ((a : ℚ) + 1)) (hA : A ^ k ≤ p) :
    ((a : ℚ) + 2) ^ k ≤ (p : ℚ) * ((a : ℚ) + 1) ^ k := by
  have hp := pow_le_pow_left₀ (show 0 ≤ (a : ℚ) + 2 by positivity) ha k
  rw [mul_pow] at hp
  exact hp.trans (mul_le_mul_of_nonneg_right hA (by positivity))

def divisorFactor (k p a : ℕ) : ℚ := ((a : ℚ) + 1) ^ k / (p : ℚ) ^ a

private theorem divisor_factor_step {k p a : ℕ} (hk : k ≤ 6) (hp : 2 ≤ p) (ha : 10 ≤ a) :
    divisorFactor k p (a + 1) ≤ divisorFactor k p a := by
  have hpQ : (2 : ℚ) ≤ p := by exact_mod_cast hp
  have haQ : (10 : ℚ) ≤ a := by exact_mod_cast ha
  have hp0 : (0 : ℚ) < p := by linarith
  have hratio : (12 / 11 : ℚ) ^ k ≤ p :=
    (pow_le_pow_right₀ (by norm_num) hk).trans ((by norm_num : (12 / 11 : ℚ) ^ 6 ≤ 2).trans hpQ)
  have hh := next_exponent_bound (p := p) (a := a) (A := 12 / 11) (by linarith) hratio
  unfold divisorFactor
  push_cast
  rw [pow_succ (p : ℚ), div_le_div_iff₀ (mul_pos (pow_pos hp0 _) hp0) (pow_pos hp0 _)]
  calc
    ((a : ℚ) + 1 + 1) ^ k * (p : ℚ) ^ a ≤
        ((p : ℚ) * ((a : ℚ) + 1) ^ k) * (p : ℚ) ^ a := by
      apply mul_le_mul_of_nonneg_right _ (pow_nonneg hp0.le _)
      simpa only [add_assoc, one_add_one_eq_two] using hh
    _ = _ := by ring

theorem divisor_factor_le_cap {k p : ℕ} (hk : k ≤ 6) (hp : 2 ≤ p) (a : ℕ) :
    divisorFactor k p a ≤ boundedEulerCap (divisorFactor k) 10 p := by
  unfold boundedEulerCap
  by_cases h : p < 64
  · rw [if_pos h]
    exact factor_le_finite_max _ (fun _ ha => divisor_factor_step hk hp ha) a
  · rw [if_neg h]
    have hpN : 0 < p := by omega
    have hpQ : (0 : ℚ) < p := by exact_mod_cast hpN
    exact (div_le_one (pow_pos hpQ _)).mpr (divisor_ratio_large hk (Nat.le_of_not_gt h) a)

def sigmaRatioFactor (p a : ℕ) : ℚ :=
  ((a : ℚ) + 1) ^ 6 * (p : ℚ) ^ (5 * a) / (σ 1 (p ^ a) : ℚ) ^ 6

theorem sigma_one_prime_pow_succ {p a : ℕ} (hp : p.Prime) :
    σ 1 (p ^ (a + 1)) = p * σ 1 (p ^ a) + 1 := by
  rw [ArithmeticFunction.sigma_one_apply_prime_pow hp, ArithmeticFunction.sigma_one_apply_prime_pow hp,
    Finset.sum_range_succ']
  simp only [pow_succ, pow_zero, ← Finset.mul_sum, mul_comm]

private theorem sigma_power_ge_power {p a : ℕ} (hp : p.Prime) :
    p ^ a ≤ σ 1 (p ^ a) := by
  rw [ArithmeticFunction.sigma_one_apply_prime_pow hp, Finset.sum_range_succ]
  exact Nat.le_add_left _ _

theorem sigma_ratio_factor_step {p a : ℕ} (hp : p.Prime)
    {A : ℚ} (ha : (a : ℚ) + 2 ≤ A * ((a : ℚ) + 1)) (hA : A ^ 6 ≤ p) :
    sigmaRatioFactor p (a + 1) ≤ sigmaRatioFactor p a := by
  have hp0 : (0 : ℚ) < p := by exact_mod_cast hp.pos
  let s : ℚ := σ 1 (p ^ a)
  have hs0 : 0 < s := by
    have h : 0 < σ 1 (p ^ a) := ArithmeticFunction.sigma_pos _ _ (pow_ne_zero _ hp.ne_zero)
    dsimp [s]
    exact_mod_cast h
  have hh := next_exponent_bound ha hA
  have hS : (p : ℚ) * s ≤ (σ 1 (p ^ (a + 1)) : ℚ) := by
    have hn : p * σ 1 (p ^ a) ≤ σ 1 (p ^ (a + 1)) := by rw [sigma_one_prime_pow_succ hp]; omega
    dsimp [s]
    exact_mod_cast hn
  calc
    sigmaRatioFactor p (a + 1) ≤
        ((a : ℚ) + 2) ^ 6 * (p : ℚ) ^ (5 * (a + 1)) / ((p : ℚ) * s) ^ 6 := by
      unfold sigmaRatioFactor
      push_cast
      rw [show (a : ℚ) + 1 + 1 = (a : ℚ) + 2 by ring]
      exact div_le_div_of_nonneg_left (by positivity) (by positivity)
        (pow_le_pow_left₀ (by positivity) hS _)
    _ = (((a : ℚ) + 2) ^ 6 / p) * ((p : ℚ) ^ (5 * a) / s ^ 6) := by
      rw [Nat.mul_add, pow_add, mul_pow]
      norm_num only [Nat.mul_one]
      field_simp [hp0.ne', hs0.ne']
    _ ≤ ((a : ℚ) + 1) ^ 6 * ((p : ℚ) ^ (5 * a) / s ^ 6) := by
      apply mul_le_mul_of_nonneg_right _ (div_nonneg (pow_nonneg hp0.le _) (pow_nonneg hs0.le _))
      exact (div_le_iff₀ hp0).mpr (by simpa only [mul_comm] using hh)
    _ = sigmaRatioFactor p a := by unfold sigmaRatioFactor s; ring

theorem sigma_ratio_factor_large {p a : ℕ} (hp : p.Prime) (hle : 64 ≤ p) :
    sigmaRatioFactor p a ≤ 1 := by
  have hp0 : (0 : ℚ) < p := by exact_mod_cast hp.pos
  have hs : (p : ℚ) ^ a ≤ (σ 1 (p ^ a) : ℚ) := by exact_mod_cast sigma_power_ge_power (a := a) hp
  unfold sigmaRatioFactor
  have hs0 := (pow_pos hp0 a).trans_le hs
  apply (div_le_one (pow_pos hs0 _)).mpr
  calc
    _ ≤ (p : ℚ) ^ a * (p : ℚ) ^ (5 * a) :=
      mul_le_mul_of_nonneg_right (divisor_ratio_large (by decide : 6 ≤ 6) hle a) (by positivity)
    _ = ((p : ℚ) ^ a) ^ 6 := by rw [← pow_add, ← pow_mul]; congr 1; omega
    _ ≤ _ := pow_le_pow_left₀ (by positivity) hs _

theorem sigma_ratio_factor_product_eq {m : ℕ} (hm : 0 < m) :
    (∏ p ∈ m.primeFactors, sigmaRatioFactor p (m.factorization p)) =
      (σ 0 m : ℚ) ^ 6 * (m : ℚ) ^ 5 / (σ 1 m : ℚ) ^ 6 := by
  let s := m.primeFactors
  let e := m.factorization
  have hτ : (σ 0 m : ℚ) = ∏ p ∈ s, (((e p : ℕ) : ℚ) + 1) := by
    have h := ArithmeticFunction.isMultiplicative_sigma.multiplicative_factorization (σ 0) hm.ne'
    rw [Finsupp.prod, Nat.support_factorization] at h
    have he : σ 0 m = ∏ p ∈ s, (e p + 1) :=
      h.trans (Finset.prod_congr rfl (fun p hp =>
        ArithmeticFunction.sigma_zero_apply_prime_pow (Nat.prime_of_mem_primeFactors hp)))
    exact_mod_cast he
  have hσ : (σ 1 m : ℚ) = ∏ p ∈ s, (σ 1 (p ^ e p) : ℚ) := by
    have h := ArithmeticFunction.isMultiplicative_sigma.multiplicative_factorization (σ 1) hm.ne'
    rw [Finsupp.prod, Nat.support_factorization] at h
    exact_mod_cast h
  have hmprod : (m : ℚ) = ∏ p ∈ s, (p : ℚ) ^ e p := by
    have h := Nat.prod_factorization_pow_eq_self hm.ne'
    rw [Finsupp.prod, Nat.support_factorization] at h
    exact_mod_cast h.symm
  unfold sigmaRatioFactor
  simp only [Finset.prod_div_distrib, Finset.prod_mul_distrib, pow_mul' (m := 5), Finset.prod_pow]
  rw [← hτ, ← hσ, ← hmprod]

/-- Bound the divisor product on any selected prime factors using the
finite local maxima for primes below 64. -/
theorem divisor_product_power_le {n k : ℕ} (hn : 0 < n) (hk : k ≤ 6)
    {s small : Finset ℕ} (hsub : s ⊆ n.primeFactors)
    (hsmall : ∀ p ∈ s, p < 64 → p ∈ small) :
    (∏ p ∈ s, ((n.factorization p : ℚ) + 1)) ^ k ≤
      (∏ p ∈ small, boundedEulerCap (divisorFactor k) 10 p) * n := by
  let e := n.factorization
  have h := euler_product_le_small_primes (divisorFactor k) 10 e
    (by intro p; simp [divisorFactor]) hsmall
    (fun p hp => divisor_factor_le_cap hk (Nat.prime_of_mem_primeFactors (hsub hp)).two_le _)
    (fun _ _ => by unfold divisorFactor; positivity)
  unfold divisorFactor at h
  rw [Finset.prod_div_distrib, Finset.prod_pow] at h
  have hpos : 0 < ∏ p ∈ s, (p : ℚ) ^ e p := by
    exact Finset.prod_pos fun p hp => pow_pos (Nat.cast_pos.mpr (Nat.prime_of_mem_primeFactors (hsub hp)).pos) _
  have hle : ∏ p ∈ s, (p : ℚ) ^ e p ≤ n := by
    have hall := Nat.prod_primeFactors_pow_factorization hn.ne'
    have hnat : ∏ p ∈ s, p ^ e p ≤ n := by
      conv_rhs => rw [hall]
      exact Finset.prod_le_prod_of_subset_of_one_le' hsub
        (fun p hp _ => Nat.one_le_pow _ _ (Nat.prime_of_mem_primeFactors hp).pos)
    exact_mod_cast hnat
  exact ((div_le_iff₀ hpos).mp h).trans (mul_le_mul_of_nonneg_left hle
    (Finset.prod_nonneg fun p _ => (one_le_euler_cap (divisorFactor k) 10 p (by simp [divisorFactor])).trans'
      (by norm_num)))

end QuadraticNorms
