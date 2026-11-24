import MathieuQuadratic.Perturbation

/-!
# Powers of the norm needed to clear the pole

The stated orders are detected by restricting the numerator to the
norm-zero locus with nonzero height norm.
-/

open scoped BigOperators TensorProduct

namespace MathieuQuadratic

section Denominators

variable {R : Type*} [CommRing R] [IsDomain R]

/-- A norm power clears the denominator if the result lies in the polynomial ring. -/
def ClearsNormDenominators (q : R) (f : Localization.Away q) (n : ℕ) : Prop :=
  (algebraMap R (Localization.Away q) q) ^ n * f ∈
    Set.range (algebraMap R (Localization.Away q))

theorem clearsNormDenominators_iff {q p : R} (hq : q ≠ 0) (hp : ¬ q ∣ p)
    (f : Localization.Away q) (r n : ℕ)
    (hf : (algebraMap R (Localization.Away q) q) ^ r * f =
      algebraMap R (Localization.Away q) p) :
    ClearsNormDenominators q f n ↔ r ≤ n := by
  let φ := algebraMap R (Localization.Away q)
  constructor
  · rintro ⟨g, hg⟩
    by_contra! hn
    have he : p = q ^ (r - n) * g := by
      apply IsLocalization.injective (Localization.Away q) (powers_le_nonZeroDivisors_of_noZeroDivisors hq)
      change φ p = φ (q ^ (r - n) * g)
      rw [map_mul, map_pow, hg, ← mul_assoc, ← pow_add, Nat.sub_add_cancel hn.le]
      exact hf.symm
    exact hp ⟨q ^ (r - n - 1) * g, by
      calc
        p = q ^ (r - n) * g := he
        _ = q ^ (1 + (r - n - 1)) * g := by congr 2; omega
        _ = _ := by rw [pow_add, pow_one, mul_assoc]⟩
  · intro hn
    refine ⟨q ^ (n - r) * p, ?_⟩
    change φ (q ^ (n - r) * p) = φ q ^ n * f
    rw [map_mul, map_pow, ← hf, ← mul_assoc, ← pow_add, Nat.sub_add_cancel hn]

end Denominators

variable {K : Type*} [Field K] [CharZero K] {M : ℕ}
variable (a : Fin 2 → K) (b : Fin M → K)

/-- Polynomial numerator of the product phase. -/
noncomputable def productNumerator (β γ : ℕ) : QuadPolynomial K M :=
  (1 + quadNorm a b) ^ β * (quadNorm a b - (1 + quadNorm a b) ^ 2 * heightNorm b) ^ γ

omit [CharZero K] in
/-- The displayed numerator and denominator of the product phase. -/
theorem productPhase_clear_order (α β γ : ℕ) :
    (normUnit a b : NormComplement a b) ^ (α + γ) *
      productPhase (normUnit a b) (height a b) α β γ =
        algebraMap (QuadPolynomial K M) (NormComplement a b) (productNumerator a b β γ) := by
  let U := normUnit a b
  let V := 1 + (U : NormComplement a b)
  let B := algebraMap (QuadPolynomial K M) (NormComplement a b) (heightNorm b)
  let W := 1 - V ^ 2 * height a b
  have hw : (U : NormComplement a b) * W = (U : NormComplement a b) - V ^ 2 * B := by
    change (U : NormComplement a b) * (1 - V ^ 2 * (B * (↑(U⁻¹) : NormComplement a b))) = _
    calc
      _ = (U : NormComplement a b) - V ^ 2 * B * ((U : NormComplement a b) * (↑(U⁻¹) : NormComplement a b)) := by ring
      _ = _ := by rw [Units.mul_inv, mul_one]
  have hu : (U : NormComplement a b) ^ α * (↑(U⁻¹) : NormComplement a b) ^ α = 1 := by
    rw [← mul_pow, Units.mul_inv, one_pow]
  rw [productNumerator, map_mul, map_pow, map_pow, map_sub, map_mul, map_pow, map_add, map_one, ← normUnit_val]
  change (U : NormComplement a b) ^ (α + γ) *
      ((↑(U⁻¹) : NormComplement a b) ^ α * V ^ β * W ^ γ) =
        V ^ β * ((U : NormComplement a b) - V ^ 2 * B) ^ γ
  rw [← hw, mul_pow, pow_add]
  calc
    _ = ((U : NormComplement a b) ^ α * (↑(U⁻¹) : NormComplement a b) ^ α) *
      (V ^ β * ((U : NormComplement a b) ^ γ * W ^ γ)) := by ring
    _ = _ := by rw [hu, one_mul]

omit [CharZero K] in
private theorem exists_norm_height_values (ha : a 0 ≠ 0) (i : Fin M) (hb : b i ≠ 0)
    (u v : AlgebraicClosure K) :
    ∃ f : QuadPolynomial K M →ₐ[K] AlgebraicClosure K,
      f (quadNorm a b) = u ∧ f (heightNorm b) = v := by
  let L := AlgebraicClosure K
  let a₀ := algebraMap K L (a 0)
  let bᵢ := algebraMap K L (b i)
  have ha₀ : a₀ ≠ 0 := by
    change algebraMap K L (a 0) ≠ 0
    exact (map_ne_zero_iff (algebraMap K L) (RingHom.injective (algebraMap K L))).mpr ha
  have hbᵢ : bᵢ ≠ 0 := by
    change algebraMap K L (b i) ≠ 0
    exact (map_ne_zero_iff (algebraMap K L) (RingHom.injective (algebraMap K L))).mpr hb
  obtain ⟨x, hx⟩ := IsAlgClosed.exists_pow_nat_eq ((u - v) / a₀) (n := 2) (by decide)
  obtain ⟨z, hz⟩ := IsAlgClosed.exists_pow_nat_eq (v / bᵢ) (n := 2) (by decide)
  let f : QuadPolynomial K M →ₐ[K] L := MvPolynomial.aeval
    (Sum.elim (fun j : Fin 2 ↦ if j = 0 then x else 0) (fun j : Fin M ↦ if j = i then z else 0))
  have hA : f (planeNorm a) = u - v := by
    rw [planeNorm, map_sum, Finset.sum_eq_single (0 : Fin 2)]
    · simp [f, hx, a₀, mul_div_cancel₀ _ ha₀]
    · intro j _ hj
      simp [f, hj]
    · simp
  have hB : f (heightNorm b) = v := by
    rw [heightNorm, map_sum, Finset.sum_eq_single i]
    · simp [f, hz, bᵢ, mul_div_cancel₀ _ hbᵢ]
    · intro j _ hj
      simp [f, hj]
    · simp
  exact ⟨f, by rw [quadNorm, map_add, hA, hB, sub_add_cancel], hB⟩

/-- The product phase has exact denominator order α + γ. -/
theorem productPhase_pole_order (ha : a 0 ≠ 0) (i : Fin M) (hb : b i ≠ 0)
    (α β γ n : ℕ) :
    ClearsNormDenominators (quadNorm a b)
      (productPhase (normUnit a b) (height a b) α β γ) n ↔ α + γ ≤ n := by
  obtain ⟨g, hg, _⟩ := exists_norm_height_values a b ha i hb 1 0
  have hq : quadNorm a b ≠ 0 := by
    intro h
    rw [h, map_zero] at hg
    exact zero_ne_one hg
  obtain ⟨g, hg, hB⟩ := exists_norm_height_values a b ha i hb 0 1
  have hp : ¬ quadNorm a b ∣ productNumerator a b β γ := by
    rintro ⟨p, hp⟩
    have he := congrArg g hp
    simp only [productNumerator, map_mul, map_pow, map_add, map_one, hg, zero_mul, map_sub,
      hB, add_zero, one_pow, one_mul, zero_sub] at he
    exact pow_ne_zero _ (neg_ne_zero.mpr one_ne_zero) he
  apply clearsNormDenominators_iff hq hp
  rw [← normUnit_val]
  exact productPhase_clear_order a b α β γ

/-- The perturbing polynomial, whose value is M at q = 0 and B = 1. -/
noncomputable def perturbationPolynomial : QuadPolynomial K M :=
  1 + heightNorm b * ((M - 1 : ℕ) + (M : QuadPolynomial K M) * quadNorm a b)

/-- The perturbation still has exact denominator order 2M. -/
theorem perturbedProductPhase_pole_order
    (ha : a 0 ≠ 0) (i : Fin M) (hb : b i ≠ 0) (n : ℕ) :
    ClearsNormDenominators (quadNorm a b)
      (perturbedProductPhase (normUnit a b) (height a b) M) n ↔ 2 * M ≤ n := by
  obtain ⟨g, hg, _⟩ := exists_norm_height_values a b ha i hb 1 0
  have hq : quadNorm a b ≠ 0 := by
    intro h
    rw [h, map_zero] at hg
    exact zero_ne_one hg
  obtain ⟨g, hg, hB⟩ := exists_norm_height_values a b ha i hb 0 1
  have hp : ¬ quadNorm a b ∣ productNumerator a b M M * perturbationPolynomial a b := by
    rintro ⟨p, hp⟩
    have he := congrArg g hp
    simp only [productNumerator, perturbationPolynomial, map_mul, map_pow, map_add, map_one,
      map_natCast, hg, zero_mul, mul_zero, map_sub, hB, add_zero, one_pow, one_mul, zero_sub] at he
    have hnon : (1 + (↑(M - 1) : AlgebraicClosure K)) ≠ 0 := by
      norm_cast
      omega
    exact mul_ne_zero (pow_ne_zero _ (neg_ne_zero.mpr one_ne_zero)) hnon he
  rw [two_mul]
  apply clearsNormDenominators_iff hq hp
  rw [← normUnit_val, perturbedProductPhase, ← mul_assoc, productPhase_clear_order, map_mul]
  congr 1
  have hB' : (normUnit a b : NormComplement a b) * height a b =
      algebraMap (QuadPolynomial K M) (NormComplement a b) (heightNorm b) := by
    rw [height, mul_left_comm, Units.mul_inv, mul_one]
  simp only [perturbationPolynomial, map_add, map_mul, map_one, map_natCast]
  rw [← normUnit_val, hB']

/-- At the minimum possible pole exponent in height rank at least two,
the only pure product phase is a proper power. -/
theorem productFamily_pole_lower_bound
    {α β γ : ℕ} (hM : 2 ≤ M) (hβ : M ≤ β) (hβα : β ≤ α) (hαγ : α ≤ γ)
    (hn : NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := AlgebraicClosure K)
      (productPhase (normUnit a b) (height a b) α β γ))) :
    2 * M + 1 ≤ α + γ := by
  by_contra! h
  have heα : α = M := by omega
  have heβ : β = M := by omega
  have heγ : γ = M := by omega
  subst α
  subst β
  subst γ
  let U := normUnit a b
  let V := (↑(U⁻¹) : NormComplement a b) * (1 + (U : NormComplement a b)) *
    (1 - (1 + (U : NormComplement a b)) ^ 2 * height a b)
  have he : productPhase U (height a b) M M M = V ^ M := by
    simp [V, productPhase, mul_pow]
  apply hn 1 (Algebra.TensorProduct.includeRight V) M hM
  simp [he]

/-- The perturbation in the draft on the localized diagonal chart. -/
theorem draft_perturbation (ha : a 0 ≠ 0) (i : Fin M) (hb : b i ≠ 0) (hM : 4 ≤ M)
    {trace : NormComplement a b →ₗ[K] K} (h : GeometricTrace a b trace) :
    let F := perturbedProductPhase (normUnit a b) (height a b) M
    NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := AlgebraicClosure K) F) ∧
      AllOrderWitness trace F ((normUnit a b : NormComplement a b) ^ M) ∧
      ∀ n : ℕ, ClearsNormDenominators (quadNorm a b) F n ↔ 2 * M ≤ n := by
  have ht := h.complementTrace.normalizer a b (by omega)
  exact ⟨perturbedProductPhase_absolute_nonpower a b ha i hb hM,
    perturbation_allOrderWitness trace _ _ hM ht.weight ht.height,
    perturbedProductPhase_pole_order a b ha i hb⟩

end MathieuQuadratic
