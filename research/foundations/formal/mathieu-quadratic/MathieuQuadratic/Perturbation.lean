import MathieuQuadratic.SharpMultiplier

/-!
# The product perturbation

The coefficients are obtained by expanding the two binomial factors
before applying the mixed-moment formula.
-/

open scoped BigOperators

namespace MathieuQuadratic

variable {K A : Type*} [Field K] [CharZero K] [CommRing A] [Algebra K A]

def perturbationCoefficient (M m j r : ℕ) : ℕ :=
  m.choose j * j.choose r * (M - 1) ^ (j - r) * M ^ r

private theorem perturbation_factor_pow (U : Aˣ) (S : A) (M m : ℕ) :
    (1 + ((U : A) * S) * ((M - 1 : ℕ) + (M : A) * (U : A))) ^ m =
      ∑ j ∈ Finset.range (m + 1), ∑ r ∈ Finset.range (j + 1),
        (perturbationCoefficient M m j r : A) * ((U : A) ^ (j + r) * S ^ j) := by
  rw [add_comm 1, add_pow]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  rw [one_pow, mul_one, mul_pow, add_comm (↑(M - 1) : A), add_pow,
    Finset.mul_sum, Finset.sum_mul]
  refine Finset.sum_congr rfl fun r _ ↦ ?_
  rw [perturbationCoefficient, Nat.cast_mul, Nat.cast_mul, Nat.cast_mul, Nat.cast_pow,
    Nat.cast_pow, mul_pow ((M : A)), mul_pow (U : A), pow_add]
  ring

omit [CharZero K] in
/-- The double sum in the perturbation proposition of the draft. -/
theorem perturbation_moment_expansion (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
    (M m ℓ j₀ : ℕ) :
    trace ((U : A) ^ ℓ * S ^ j₀ * perturbedProductPhase U S M ^ m) =
      ∑ j ∈ Finset.range (m + 1), ∑ r ∈ Finset.range (j + 1),
        (perturbationCoefficient M m j r : K) *
          trace ((U : A) ^ (ℓ + j + r) * S ^ (j₀ + j) * productPhase U S M M M ^ m) := by
  rw [perturbedProductPhase, mul_pow, perturbation_factor_pow, ← mul_assoc, Finset.mul_sum, map_sum]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  rw [Finset.mul_sum, map_sum]
  refine Finset.sum_congr rfl fun r _ ↦ ?_
  rw [show (↑U ^ ℓ * S ^ j₀ * productPhase U S M M M ^ m) *
      ((perturbationCoefficient M m j r : A) * ((U : A) ^ (j + r) * S ^ j)) =
        (perturbationCoefficient M m j r : A) *
          ((U : A) ^ (ℓ + j + r) * S ^ (j₀ + j) * productPhase U S M M M ^ m) by
      simp only [pow_add]
      ring,
    ← map_natCast (algebraMap K A), ← Algebra.smul_def, map_smul, smul_eq_mul]

private theorem choose_neg_two (n : ℕ) :
    ((Ring.choose (-2 : ℤ) n : ℤ) : ℚ) = (-1 : ℚ) ^ n * (n + 1 : ℕ) := by
  rw [Ring.choose_neg, show (2 : ℤ) + n - 1 = (↑(n + 1) : ℤ) by omega, Ring.choose_natCast]
  rw [show (n + 1).choose n = n + 1 by
    rw [← Nat.choose_symm (Nat.le_succ n), show n + 1 - n = 1 by omega, Nat.choose_one_right]]
  rw [Units.smul_def, zsmul_eq_mul, Int.cast_mul, Int.coe_negOnePow_natCast,
    Int.cast_pow, Int.cast_neg, Int.cast_one, Int.cast_natCast]
  push_cast
  rfl

/-- The first-order cancellation in the perturbation proposition. -/
theorem perturbation_first_cancellation (M : ℕ) (hM : 4 ≤ M) :
    (↑(M - 1) : ℚ) * ((Ring.choose (-2 : ℤ) (M - 1) : ℤ) : ℚ) +
      M * ((Ring.choose (-2 : ℤ) (M - 2) : ℤ) : ℚ) = 0 := by
  rw [choose_neg_two, choose_neg_two, show M - 1 = (M - 2) + 1 by omega, pow_succ,
    show M - 2 + 1 + 1 = M by omega]
  ring

variable (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
variable {M : ℕ} (hM : 4 ≤ M)
variable (weight : ∀ (a : ℤ) (i : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ i) = 0)
variable (height : ∀ i : ℕ, trace (S ^ i) = (M : K) / (↑(M + 2 * i) : K))

include hM weight height in
/-- Every below-threshold block monomial eventually gives a zero marked
moment for the perturbed phase. -/
theorem perturbation_monomial_eventually_zero {ℓ j₀ m : ℕ}
    (hℓ : ℓ < M) (hm : M + 2 * j₀ ≤ (M - 2) * m) :
    trace ((U : A) ^ ℓ * S ^ j₀ * perturbedProductPhase U S M ^ m) = 0 := by
  rw [perturbation_moment_expansion]
  apply Finset.sum_eq_zero
  intro j hj
  apply Finset.sum_eq_zero
  intro r hr
  have hjm : j ≤ m := by have := Finset.mem_range.mp hj; omega
  have hrj : r ≤ j := by have := Finset.mem_range.mp hr; omega
  have hMm : M * m = (M - 2) * m + 2 * m := by
    rw [← Nat.add_mul, Nat.sub_add_cancel (by omega)]
  have hz : (U : A) ^ (ℓ + j + r) = (↑(U ^ ((ℓ + j + r : ℕ) : ℤ)) : A) := by
    rw [zpow_natCast, Units.val_pow_eq_pow_val]
  rw [hz, product_mixedMoment trace U S M M M M m (j₀ + j) (by omega) _ (by omega) weight height]
  split_ifs with hd
  · let n := M * m - (M + 2 * (j₀ + j))
    have hn : M + 2 * (j₀ + j) ≤ M * m := by omega
    have hcast : (n : ℤ) = ((M * m : ℕ) - (M + 2 * (j₀ + j) : ℕ) : ℤ) := Int.ofNat_sub hn
    rw [← hcast, Ring.choose_natCast, Nat.choose_eq_zero_of_lt (by omega)]
    simp
  · simp

include hM weight height in
/-- The positive marked sum for the perturbed witness at orders m ≥ 2. -/
theorem perturbation_marked_of_two_le {m : ℕ} (hm : 2 ≤ m) :
    trace ((U : A) ^ M * perturbedProductPhase U S M ^ m) =
      algebraMap ℚ K (∑ j ∈ Finset.range (m + 1),
        ((m.choose j : ℚ) * (M : ℚ) ^ j) * mixedBetaMoment M (M * m) j) := by
  have hupper : M ≤ (M - 2) * m := by
    have h := Nat.mul_le_mul_left (M - 2) hm
    omega
  have hMm : M * m = (M - 2) * m + 2 * m := by
    rw [← Nat.add_mul, Nat.sub_add_cancel (by omega)]
  have he : (U : A) ^ M * perturbedProductPhase U S M ^ m =
      (U : A) ^ M * S ^ 0 * perturbedProductPhase U S M ^ m := by simp
  rw [he, perturbation_moment_expansion, map_sum]
  refine Finset.sum_congr rfl fun j hj ↦ ?_
  have hjm : j ≤ m := by have := Finset.mem_range.mp hj; omega
  have ht : j ∈ Finset.range (j + 1) := Finset.mem_range.mpr (by omega)
  have hj0 : M + 2 * j ≤ M * m := by omega
  let n := M * m - (M + 2 * j)
  have hn : (n : ℤ) = ((M * m : ℕ) - (M + 2 * j : ℕ) : ℤ) := Int.ofNat_sub hj0
  rw [Finset.sum_eq_single_of_mem j ht]
  · simp only [zero_add]
    have hz : (U : A) ^ (M + j + j) = (↑(U ^ ((M + j + j : ℕ) : ℤ)) : A) := by
      rw [zpow_natCast, Units.val_pow_eq_pow_val]
    rw [hz, product_mixedMoment trace U S M M M M m j (by omega) _ (by omega) weight height,
      if_pos (by omega), ← hn,
      show ((M * m : ℕ) - (M + j + j : ℕ) : ℤ).toNat = n by omega,
      Ring.choose_natCast, Nat.choose_self]
    simp [perturbationCoefficient]
  · intro r hr hrne
    have hrj : r < j := by have := Finset.mem_range.mp hr; omega
    simp only [zero_add]
    have hz : (U : A) ^ (M + j + r) = (↑(U ^ ((M + j + r : ℕ) : ℤ)) : A) := by
      rw [zpow_natCast, Units.val_pow_eq_pow_val]
    rw [hz, product_mixedMoment trace U S M M M M m j (by omega) _ (by omega) weight height,
      if_pos (by omega), ← hn, Ring.choose_natCast, Nat.choose_eq_zero_of_lt (by omega)]
    simp

include hM weight height in
/-- The perturbed phase is pure at every order at least two. -/
theorem perturbation_unmarked_of_two_le {m : ℕ} (hm : 2 ≤ m) :
    trace (perturbedProductPhase U S M ^ m) = 0 := by
  have hupper : M ≤ (M - 2) * m := by
    have h := Nat.mul_le_mul_left (M - 2) hm
    omega
  simpa using perturbation_monomial_eventually_zero trace U S hM weight height
    (ℓ := 0) (j₀ := 0) (by omega) (by simpa using hupper)

omit [CharZero K] in
private theorem perturbation_first_expanded (ℓ : ℕ) :
    trace ((U : A) ^ ℓ * perturbedProductPhase U S M) =
      trace ((U : A) ^ ℓ * productPhase U S M M M) +
      (↑(M - 1) : K) * trace ((U : A) ^ (ℓ + 1) * S * productPhase U S M M M) +
      M * trace ((U : A) ^ (ℓ + 2) * S * productPhase U S M M M) := by
  have he := perturbation_moment_expansion trace U S M 1 ℓ 0
  norm_num [Finset.sum_range_succ, perturbationCoefficient, Nat.choose] at he
  simpa [Nat.add_assoc, add_assoc] using he

include hM weight height in
/-- The pure and marked first moments for the perturbed phase. -/
theorem perturbation_first_moments :
    trace (perturbedProductPhase U S M) = 0 ∧
      trace ((U : A) ^ M * perturbedProductPhase U S M) =
        algebraMap ℚ K (mixedBetaMoment M M 0) := by
  have huniform := product_uniform_moments trace U S M M 1 (by omega) le_rfl Nat.one_pos weight height
  have hp0 : trace (productPhase U S M M M) = 0 := by simpa using huniform.1
  have hmark : trace ((U : A) ^ M * productPhase U S M M M) =
      algebraMap ℚ K (mixedBetaMoment M M 0) := by simpa using huniform.2
  have hterm (ℓ : ℕ) (hℓ : 1 ≤ ℓ) :
      trace ((U : A) ^ ℓ * S * productPhase U S M M M) =
        if ℓ ≤ M then algebraMap ℚ K
          (mixedBetaMoment M M 1 * ((Ring.choose (-2 : ℤ) (M - ℓ) : ℤ) : ℚ))
        else 0 := by
    have he : (U : A) ^ ℓ * S * productPhase U S M M M =
      (↑(U ^ (ℓ : ℤ)) : A) * S ^ 1 * productPhase U S M M M ^ 1 := by simp
    rw [he, product_mixedMoment trace U S M M M M 1 1 (by omega) ℓ (by omega) weight height]
    simp only [mul_one]
    have huv : ((M : ℤ) - (M + 2 * 1 : ℕ)) = -2 := by omega
    rw [huv]
    split_ifs with h₁ h₂
    · rw [← Int.ofNat_sub h₂, Int.toNat_natCast]
    · omega
    · omega
    · rfl
  constructor
  · have he : trace (perturbedProductPhase U S M) =
        trace ((U : A) ^ 0 * perturbedProductPhase U S M) := by simp
    rw [he, perturbation_first_expanded, pow_zero, one_mul, hp0, zero_add,
      hterm 1 (by omega), hterm 2 (by omega), if_pos (by omega), if_pos (by omega)]
    rw [← map_natCast (algebraMap ℚ K), ← map_natCast (algebraMap ℚ K), ← map_mul, ← map_mul,
      ← map_add]
    have hcancel := perturbation_first_cancellation M hM
    calc
      _ = algebraMap ℚ K (mixedBetaMoment M M 1 *
          ((↑(M - 1) : ℚ) * ((Ring.choose (-2 : ℤ) (M - 1) : ℤ) : ℚ) +
            M * ((Ring.choose (-2 : ℤ) (M - 2) : ℤ) : ℚ))) := by congr 1; ring
      _ = _ := by rw [hcancel, mul_zero, map_zero]
  · rw [perturbation_first_expanded, hmark, hterm (M + 1) (by omega),
      hterm (M + 2) (by omega), if_neg (by omega), if_neg (by omega)]
    ring

include hM weight height in
/-- All positive powers of the perturbation have zero pure moment, and
the prescribed norm-power multiplier detects every order. -/
theorem perturbation_allOrderWitness :
    AllOrderWitness trace (perturbedProductPhase U S M) ((U : A) ^ M) := by
  have hfirst := perturbation_first_moments trace U S hM weight height
  have hbeta (r : ℚ) (hr : 0 < r) : algebraMap ℚ K r ≠ 0 :=
    (map_ne_zero_iff (algebraMap ℚ K) (RingHom.injective (algebraMap ℚ K))).mpr hr.ne'
  constructor
  · intro m hm
    by_cases he : m = 1
    · simpa [he] using hfirst.1
    · exact perturbation_unmarked_of_two_le trace U S hM weight height (by omega)
  · intro m hm
    by_cases he : m = 1
    · rw [he, pow_one, hfirst.2]
      exact hbeta _ (mixedBetaMoment_pos M M 0 (by omega))
    · rw [perturbation_marked_of_two_le trace U S hM weight height (by omega)]
      apply hbeta
      apply Finset.sum_pos
      · intro j hj
        have hjm : j ≤ m := by have := Finset.mem_range.mp hj; omega
        exact mul_pos (mul_pos (by exact_mod_cast Nat.choose_pos hjm) (by positivity))
          (mixedBetaMoment_pos M _ _ (by omega))
      · exact ⟨0, Finset.mem_range.mpr (by omega)⟩

include hM weight height in
/-- Vanishing below the minimum detecting degree on the averaged block span. -/
theorem perturbation_cutoff_span {J : ℕ} (hJ : J < M) {h : A}
    (hspan : h ∈ Submodule.span K {v : A | ∃ ℓ j : ℕ, j ≤ ℓ ∧ ℓ ≤ J ∧ v = (U : A) ^ ℓ * S ^ j})
    (m : ℕ) (hm : M + 2 * J ≤ (M - 2) * m) :
    trace (h * perturbedProductPhase U S M ^ m) = 0 := by
  induction hspan using Submodule.span_induction with
  | mem x hx =>
      obtain ⟨ℓ, j, hj, hℓ, rfl⟩ := hx
      exact perturbation_monomial_eventually_zero trace U S hM weight height (by omega) (by omega)
  | zero => simp
  | add x y _ _ hx hy => rw [add_mul, map_add, hx, hy, add_zero]
  | smul c x _ hx => rw [smul_mul_assoc, map_smul, hx, smul_zero]

end MathieuQuadratic
