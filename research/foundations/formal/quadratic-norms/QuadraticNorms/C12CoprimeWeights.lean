import QuadraticNorms.C12EulerProduct
import QuadraticNorms.PointwiseSingularSeries
import Mathlib.NumberTheory.ArithmeticFunction.Misc

open scoped BigOperators ArithmeticFunction.Moebius Classical

namespace QuadraticNorms

/-- Allowed prime divisors of the target. -/
def c12AllowedPrimeFactors (n : ℕ) : Finset ℕ :=
  n.primeFactors.filter (fun p => p % 12 = 1)

noncomputable def c12AllowedDivisors (n : ℕ) : Finset ℕ := n.divisors.filter C12

/-- The coprime factor P₁₂(N) in the pointwise main term. -/
noncomputable def c12CoprimeFactor (n : ℕ) : ℝ :=
  ∏ p ∈ c12AllowedPrimeFactors n, (1 - (p : ℝ)⁻¹)

noncomputable def c12AbsoluteDivisorMass (n : ℕ) : ℝ :=
  ∑ d ∈ c12AllowedDivisors n, |(μ d : ℝ)| / d

noncomputable def c12DivisorLogMoment (n : ℕ) : ℝ :=
  ∑ d ∈ c12AllowedDivisors n, |(μ d : ℝ)| * Real.log d / d

private theorem allowed_prime {n p : ℕ} (hp : p ∈ c12AllowedPrimeFactors n) : p.Prime :=
  Nat.prime_of_mem_primeFactors (Finset.mem_filter.mp hp).1

private theorem c12_allowed_prime_prod_squarefree (n : ℕ) :
    Squarefree (∏ p ∈ c12AllowedPrimeFactors n, p) := by
  apply Finset.squarefree_prod_of_pairwise_isCoprime
  · intro p hp q hq hne
    exact Nat.coprime_iff_isRelPrime.mp ((Nat.coprime_primes (allowed_prime hp) (allowed_prime hq)).mpr hne)
  · exact fun p hp => (allowed_prime hp).squarefree

private theorem prime_subset_dvd {s t : Finset ℕ} (hst : s ⊆ t) :
    (∏ p ∈ s, p) ∣ ∏ p ∈ t, p :=
  Finset.prod_dvd_prod_of_subset s t id hst

private theorem c12_allowed_squarefree_divisors (n : ℕ) (hn : n ≠ 0) :
    (c12AllowedDivisors n).filter Squarefree = (∏ p ∈ c12AllowedPrimeFactors n, p).divisors := by
  let D := ∏ p ∈ c12AllowedPrimeFactors n, p
  have hDsf := c12_allowed_prime_prod_squarefree n
  have hD : D ∣ n := (prime_subset_dvd (Finset.filter_subset ..)).trans (Nat.prod_primeFactors_dvd n)
  have hpf : D.primeFactors = c12AllowedPrimeFactors n := Nat.primeFactors_prod (fun _ hp => allowed_prime hp)
  ext d
  simp only [c12AllowedDivisors, Finset.mem_filter, Nat.mem_divisors]
  change ((d ∣ n ∧ n ≠ 0) ∧ C12 d) ∧ Squarefree d ↔ d ∣ D ∧ D ≠ 0
  constructor
  · rintro ⟨⟨⟨hd, _⟩, hc⟩, hsf⟩
    refine ⟨?_, hDsf.ne_zero⟩
    rw [← Nat.prod_primeFactors_of_squarefree hsf]
    apply prime_subset_dvd
    intro p hp
    have hp' := Nat.mem_primeFactors.mp hp
    exact Finset.mem_filter.mpr ⟨Nat.primeFactors_mono hd hn hp, hc.2 p hp'.1 hp'.2.1⟩
  · rintro ⟨hd, _⟩
    have hsf := hDsf.squarefree_of_dvd hd
    have hdpos := Nat.pos_of_ne_zero hsf.ne_zero
    refine ⟨⟨⟨hd.trans hD, hn⟩, ⟨hdpos, fun p hp hpd => ?_⟩⟩, hsf⟩
    have hmem : p ∈ D.primeFactors := Nat.mem_primeFactors.mpr ⟨hp, hpd.trans hd, hDsf.ne_zero⟩
    exact (Finset.mem_filter.mp (hpf ▸ hmem)).2

private noncomputable def reciprocalArithmetic : ArithmeticFunction ℝ :=
  ⟨fun n => (n : ℝ)⁻¹, by simp⟩

private theorem reciprocal_multiplicative : reciprocalArithmetic.IsMultiplicative := by
  refine ⟨by simp [reciprocalArithmetic], ?_⟩
  intro m n _
  simp [reciprocalArithmetic, Nat.cast_mul, mul_comm]

/-- The signed Möbius weight in the coprime main term is exactly the
Euler product over allowed prime divisors. -/
theorem c12_signed_divisor_sum {n : ℕ} (hn : n ≠ 0) :
    (∑ d ∈ c12AllowedDivisors n, (μ d : ℝ) / d) = c12CoprimeFactor n := by
  let D := ∏ p ∈ c12AllowedPrimeFactors n, p
  have hsf := c12_allowed_prime_prod_squarefree n
  have h := ArithmeticFunction.IsMultiplicative.prodPrimeFactors_one_sub_of_squarefree
    reciprocalArithmetic reciprocal_multiplicative hsf
  change (∏ p ∈ D.primeFactors, (1 - reciprocalArithmetic p)) =
    ∑ d ∈ D.divisors, (μ d : ℝ) * reciprocalArithmetic d at h
  rw [Nat.primeFactors_prod (fun _ hp => allowed_prime hp),
    ← c12_allowed_squarefree_divisors n hn] at h
  unfold c12CoprimeFactor
  simp only [reciprocalArithmetic, ArithmeticFunction.coe_mk, ← div_eq_mul_inv] at h
  rw [h]
  symm
  exact Finset.sum_subset (Finset.filter_subset ..) (fun d hd hnot => by
    have hnsf : ¬ Squarefree d := fun hsf' => hnot (Finset.mem_filter.mpr ⟨hd, hsf'⟩)
    simp [ArithmeticFunction.moebius_eq_zero_of_not_squarefree hnsf])

/-- The coprime Euler factor is positive for every target. -/
theorem c12_coprime_factor_pos (n : ℕ) : 0 < c12CoprimeFactor n := by
  apply Finset.prod_pos
  intro p hp
  have : (1 : ℝ) < p := by exact_mod_cast (allowed_prime hp).one_lt
  exact sub_pos.mpr (inv_lt_one_of_one_lt₀ this)

private theorem c12_squarefree_divisors_sum {n : ℕ} (hn : n ≠ 0) (f : ℕ → ℝ) :
    (∑ d ∈ (c12AllowedDivisors n).filter Squarefree, f d) =
      ∑ t ∈ (c12AllowedPrimeFactors n).powerset, f (∏ p ∈ t, p) := by
  rw [c12_allowed_squarefree_divisors n hn]
  let s := c12AllowedPrimeFactors n
  let D := ∏ p ∈ s, p
  have hsf : Squarefree D := c12_allowed_prime_prod_squarefree n
  have hEq : D.primeFactors = s := Nat.primeFactors_prod (fun _ hp => allowed_prime hp)
  change (∑ d ∈ D.divisors, f d) = ∑ t ∈ s.powerset, f (∏ p ∈ t, p)
  symm
  refine Finset.sum_bij' (fun t _ => ∏ p ∈ t, p) (fun d _ => d.primeFactors) ?_ ?_ ?_ ?_ ?_
  · intro t ht
    exact Nat.mem_divisors.mpr ⟨prime_subset_dvd (Finset.mem_powerset.mp ht), hsf.ne_zero⟩
  · intro d hd
    exact Finset.mem_powerset.mpr (hEq ▸ Nat.primeFactors_mono (Nat.dvd_of_mem_divisors hd) hsf.ne_zero)
  · intro t ht
    exact Nat.primeFactors_prod (fun p hp => allowed_prime (Finset.mem_powerset.mp ht hp))
  · intro d hd
    exact Nat.prod_primeFactors_of_squarefree (hsf.squarefree_of_dvd (Nat.dvd_of_mem_divisors hd))
  · intros
    rfl

private theorem c12_weighted_divisor_sum {n : ℕ} (hn : n ≠ 0) (f : ℕ → ℝ) :
    (∑ d ∈ c12AllowedDivisors n, |(μ d : ℝ)| * f d / d) =
      ∑ t ∈ (c12AllowedPrimeFactors n).powerset,
        (∏ p ∈ t, (p : ℝ)⁻¹) * f (∏ p ∈ t, p) := by
  calc
    _ = ∑ d ∈ (c12AllowedDivisors n).filter Squarefree, |(μ d : ℝ)| * f d / d := by
      symm
      exact Finset.sum_subset (Finset.filter_subset ..) (fun d hd hnot => by
        have hnsf : ¬ Squarefree d := fun hs => hnot (Finset.mem_filter.mpr ⟨hd, hs⟩)
        simp [ArithmeticFunction.moebius_eq_zero_of_not_squarefree hnsf])
    _ = ∑ d ∈ (c12AllowedDivisors n).filter Squarefree, f d / d := by
      apply Finset.sum_congr rfl
      intro d hd
      have h : |(μ d : ℝ)| = 1 := by
        exact_mod_cast ArithmeticFunction.abs_moebius_eq_one_of_squarefree (Finset.mem_filter.mp hd).2
      rw [h, one_mul]
    _ = ∑ t ∈ (c12AllowedPrimeFactors n).powerset, f (∏ p ∈ t, p) / (∏ p ∈ t, p : ℕ) :=
      c12_squarefree_divisors_sum hn (fun d => f d / d)
    _ = _ := by
      apply Finset.sum_congr rfl
      intro t _
      rw [div_eq_mul_inv, Nat.cast_prod, ← Finset.prod_inv_distrib, mul_comm]

/-- The exact Euler product for the absolute divisor mass S₀. -/
theorem c12_absolute_divisor_mass_product {n : ℕ} (hn : n ≠ 0) :
    c12AbsoluteDivisorMass n = ∏ p ∈ c12AllowedPrimeFactors n, (1 + (p : ℝ)⁻¹) := by
  have h := c12_weighted_divisor_sum hn (fun _ => 1)
  simpa only [mul_one, one_mul, c12AbsoluteDivisorMass, Finset.prod_one_add] using h

private theorem powerset_weighted_moment {ι : Type*} [DecidableEq ι] (s : Finset ι)
    (a w : ι → ℝ) (ha : ∀ i ∈ s, 1 + a i ≠ 0) :
    (∑ t ∈ s.powerset, (∏ i ∈ t, a i) * ∑ i ∈ t, w i) =
      (∏ i ∈ s, (1 + a i)) * ∑ i ∈ s, a i * w i / (1 + a i) := by
  induction s using Finset.induction_on with
  | empty => simp
  | insert x s hx ih =>
    have hs := ih (fun i hi => ha i (Finset.mem_insert_of_mem hi))
    have hnz := ha x (Finset.mem_insert_self ..)
    rw [Finset.sum_powerset_insert hx, hs, Finset.prod_insert hx, Finset.sum_insert hx]
    have hterm (t : Finset ι) (ht : t ∈ s.powerset) :
        (∏ i ∈ insert x t, a i) * (∑ i ∈ insert x t, w i) =
          (a x * w x) * (∏ i ∈ t, a i) + a x * ((∏ i ∈ t, a i) * ∑ i ∈ t, w i) := by
      have hxt : x ∉ t := fun h => hx (Finset.mem_powerset.mp ht h)
      rw [Finset.prod_insert hxt, Finset.sum_insert hxt]
      ring
    rw [Finset.sum_congr rfl hterm, Finset.sum_add_distrib, ← Finset.mul_sum,
      ← Finset.mul_sum, ← Finset.prod_one_add, hs]
    field_simp
    ring

/-- The exact logarithmic moment S₁, including the p+1 normalization
from deleting that prime from the product over the other divisors. -/
theorem c12_divisor_log_moment_product {n : ℕ} (hn : n ≠ 0) :
    c12DivisorLogMoment n = c12AbsoluteDivisorMass n *
      ∑ p ∈ c12AllowedPrimeFactors n, Real.log p / ((p : ℝ) + 1) := by
  let s := c12AllowedPrimeFactors n
  rw [c12_absolute_divisor_mass_product hn]
  unfold c12DivisorLogMoment
  rw [c12_weighted_divisor_sum hn]
  have hlog (t : Finset ℕ) (ht : t ∈ s.powerset) :
      Real.log ((∏ p ∈ t, p : ℕ) : ℝ) = ∑ p ∈ t, Real.log (p : ℝ) := by
    rw [Nat.cast_prod]
    exact Real.log_prod (fun p hp => by exact_mod_cast (allowed_prime (Finset.mem_powerset.mp ht hp)).ne_zero)
  have hsum : (∑ t ∈ s.powerset, (∏ p ∈ t, (p : ℝ)⁻¹) * Real.log ((∏ p ∈ t, p : ℕ) : ℝ)) =
      ∑ t ∈ s.powerset, (∏ p ∈ t, (p : ℝ)⁻¹) * ∑ p ∈ t, Real.log (p : ℝ) := by
    exact Finset.sum_congr rfl (fun t ht => by rw [hlog t ht])
  rw [hsum]
  rw [powerset_weighted_moment s (fun p => (p : ℝ)⁻¹) (fun p => Real.log p)
    (fun p _ => ne_of_gt (by positivity))]
  congr 1
  apply Finset.sum_congr rfl
  intro p hp
  have hp0 : (p : ℝ) ≠ 0 := by exact_mod_cast (allowed_prime hp).ne_zero
  have hp1 : (p : ℝ) + 1 ≠ 0 := by positivity
  field_simp [hp0, hp1]

/-- Division by an allowed common divisor preserves the admissible
residue and every exceptional prime divisor of the target. -/
theorem c12_target_factor_div {n d : ℕ} (hn : 0 < n) (hd : d ∣ n) (hc : C12 d) :
    (n / d) % 12 = n % 12 ∧
      pointwiseTargetFactor c12ExceptionalPrime (n / d) = pointwiseTargetFactor c12ExceptionalPrime n := by
  have hq : 0 < n / d := Nat.div_pos (Nat.le_of_dvd hn hd) hc.1
  have hm : d * (n / d) = n := Nat.mul_div_cancel' hd
  constructor
  · conv_rhs => rw [← hm, Nat.mul_mod, C12_mod_twelve hc, one_mul, Nat.mod_mod]
  · unfold pointwiseTargetFactor
    have heq : (n / d).primeFactors.filter c12ExceptionalPrime = n.primeFactors.filter c12ExceptionalPrime := by
      ext p
      simp only [Finset.mem_filter, Nat.mem_primeFactors]
      constructor
      · rintro ⟨⟨hp, hpq, _⟩, hb⟩
        exact ⟨⟨hp, hpq.trans (Nat.div_dvd_of_dvd hd), hn.ne'⟩, hb⟩
      · rintro ⟨⟨hp, hpn, _⟩, hb⟩
        have hgood : ¬p ∣ d := fun h => hb.2.2.2 (hc.2 p hp h)
        have hpq := ((hp.dvd_mul.mp (hm.symm ▸ hpn))).resolve_left hgood
        exact ⟨⟨hp, hpq, hq.ne'⟩, hb⟩
    rw [heq]

end QuadraticNorms
