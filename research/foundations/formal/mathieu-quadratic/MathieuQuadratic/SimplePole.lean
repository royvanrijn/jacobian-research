import MathieuQuadratic.LinearPhaseEndpoint
import MathieuQuadratic.PoleOrder

/-!
# The simple-pole phase with a linear multiplier
-/

open Polynomial
open scoped BigOperators

namespace MathieuQuadratic

variable {K A : Type*} [Field K] [CharZero K] [CommRing A] [Algebra K A]

/-- The normalized moment for x^(2k) times the (m−k)-th power
of the rank-two squared norm. -/
noncomputable def simplePoleGaussianMoment (N m k : ℕ) : ℚ :=
  (ascPochhammer ℚ k).eval (1 / 2) * ((m - k).factorial : ℚ) /
    (ascPochhammer ℚ m).eval ((N : ℚ) / 2)

private def evalAtX (K : Type*) [Field K] [CharZero K] [Algebra K A] (x : A) : ℚ[X] →+* A :=
  eval₂RingHom ((algebraMap K A).comp (algebraMap ℚ K)) x

private theorem simplePole_trace_eval (trace : A →ₗ[K] K) (U : Aˣ) (x B : A)
    (N m k : ℕ) (P : ℚ[X])
    (weight : ∀ i : ℕ, i ≠ 2 * k →
      trace ((↑(U⁻¹) : A) ^ m * x ^ i * B ^ (m - k)) = 0)
    (height : trace ((↑(U⁻¹) : A) ^ m * x ^ (2 * k) * B ^ (m - k)) =
      algebraMap ℚ K (simplePoleGaussianMoment N m k)) :
    trace ((↑(U⁻¹) : A) ^ m * evalAtX K x P * B ^ (m - k)) =
      algebraMap ℚ K (P.coeff (2 * k) * simplePoleGaussianMoment N m k) := by
  induction P using Polynomial.induction_on' with
  | add p q hp hq =>
      rw [map_add, mul_add, add_mul, map_add, hp, hq, coeff_add, add_mul, map_add]
  | monomial i r =>
      have he : evalAtX K x (monomial i r) =
          algebraMap K A (algebraMap ℚ K r) * x ^ i := eval₂_monomial _ _
      rw [he, mul_left_comm _ (algebraMap K A (algebraMap ℚ K r)), mul_assoc,
        ← Algebra.smul_def, map_smul, smul_eq_mul, coeff_monomial]
      by_cases hi : i = 2 * k
      · subst i
        rw [height, if_pos rfl, map_mul]
      · rw [weight i hi, if_neg hi, mul_zero, zero_mul, map_zero]

private theorem simplePole_term_factor (N m k : ℕ) (hN : 0 < N) (hkm : k ≤ m) (c : ℚ) :
    ((m.choose k : ℚ) * (-1) ^ (m - k)) * (c * simplePoleGaussianMoment N m k) =
      ((-1 : ℚ) ^ m * m.factorial / (ascPochhammer ℚ m).eval ((N : ℚ) / 2)) *
        (negativeHalfBinomial k * c) := by
  have hs : (-1 : ℚ) ^ (m - k) = (-1 : ℚ) ^ m * (-1 : ℚ) ^ k := by
    have h : (-1 : ℚ) ^ k * (-1 : ℚ) ^ k = 1 := by rw [← mul_pow]; norm_num
    calc
      _ = (-1 : ℚ) ^ (m - k) * ((-1 : ℚ) ^ k * (-1 : ℚ) ^ k) := by rw [h, mul_one]
      _ = _ := by rw [← mul_assoc, ← pow_add, Nat.sub_add_cancel hkm]
  have hf : (m.choose k : ℚ) * k.factorial * (m - k).factorial = m.factorial := by
    exact_mod_cast Nat.choose_mul_factorial_mul_factorial hkm
  have hd : 0 < (ascPochhammer ℚ m).eval ((N : ℚ) / 2) := ascPochhammer_pos _ _ (by positivity)
  rw [hs, simplePoleGaussianMoment, negativeHalfBinomial]
  rw [← hf]
  field_simp

/-- Coefficient extraction for the simple-pole phase, from the spherical
mixed moments and the scalar weights in the draft. -/
theorem simplePole_moment_coeff (trace : A →ₗ[K] K) (U : Aˣ) (x B : A)
    (N m ℓ : ℕ) (hN : 0 < N)
    (weight : ∀ k i : ℕ, k ≤ m → i ≠ 2 * k →
      trace ((↑(U⁻¹) : A) ^ m * x ^ i * B ^ (m - k)) = 0)
    (height : ∀ k : ℕ, k ≤ m →
      trace ((↑(U⁻¹) : A) ^ m * x ^ (2 * k) * B ^ (m - k)) =
        algebraMap ℚ K (simplePoleGaussianMoment N m k)) :
    trace (x ^ ℓ * simplePolePhase U x B ^ m) =
      algebraMap ℚ K (((-1 : ℚ) ^ m * m.factorial / (ascPochhammer ℚ m).eval ((N : ℚ) / 2)) *
        (linearPhaseEndpoint m).coeff (m + ℓ)) := by
  rw [simplePolePhase, mul_pow, mul_pow, sub_eq_add_neg, add_pow (1 + 2 * x) (-B), ← mul_assoc,
    Finset.mul_sum, map_sum, linearPhaseEndpoint_coeff_eq_sum, Finset.mul_sum, map_sum]
  refine Finset.sum_congr rfl fun k hk ↦ ?_
  have hk' : k ≤ m := by have := Finset.mem_range.mp hk; omega
  let P : ℚ[X] := (1 + X) ^ m * (1 + 2 * X) ^ k * X ^ ℓ
  have he : evalAtX K x P = (1 + x) ^ m * (1 + 2 * x) ^ k * x ^ ℓ := by
    change P.eval₂ ((algebraMap K A).comp (algebraMap ℚ K)) x = _
    norm_num [P, eval₂_mul, eval₂_pow, eval₂_add, ← C_ofNat]
  have hc : x ^ ℓ * ((↑(U⁻¹) : A) ^ m * (1 + x) ^ m) *
      ((1 + 2 * x) ^ k * (-B) ^ (m - k) * (m.choose k : A)) =
        algebraMap K A (algebraMap ℚ K ((m.choose k : ℚ) * (-1) ^ (m - k))) *
          ((↑(U⁻¹) : A) ^ m * evalAtX K x P * B ^ (m - k)) := by
    rw [he, neg_pow]
    simp only [map_mul, map_natCast, map_pow, map_neg, map_one]
    ring
  rw [hc, ← Algebra.smul_def, map_smul, smul_eq_mul,
    simplePole_trace_eval trace U x B N m k P (fun i hi ↦ weight k i hk' hi) (height k hk'),
    ← map_mul, simplePole_term_factor N m k hN hk']

private theorem simplePole_marked_factor (N m : ℕ) :
    ((-1 : ℚ) ^ m * m.factorial / (ascPochhammer ℚ m).eval ((N : ℚ) / 2)) *
      ((2 * (m : ℚ) + 1) * negativeHalfBinomial m * 2 ^ m / (m + 1 : ℕ)) =
      2 ^ m * (ascPochhammer ℚ m).eval (3 / 2) /
        (((m + 1 : ℕ) : ℚ) * (ascPochhammer ℚ m).eval ((N : ℚ) / 2)) := by
  have hp : (2 * (m : ℚ) + 1) * (ascPochhammer ℚ m).eval (1 / 2) =
      (ascPochhammer ℚ m).eval (3 / 2) := by
    have he : (ascPochhammer ℚ (m + 1)).eval (1 / 2) =
        (1 / 2) * (ascPochhammer ℚ m).eval (3 / 2) := by
      rw [ascPochhammer_succ_left, eval_mul, eval_X, eval_comp, eval_add, eval_X, eval_one]
      norm_num
    rw [ascPochhammer_succ_eval] at he
    linear_combination 2 * he
  have hs : (-1 : ℚ) ^ m * (-1 : ℚ) ^ m = 1 := by rw [← mul_pow]; norm_num
  have hf : (m.factorial : ℚ) ≠ 0 := by exact_mod_cast Nat.factorial_ne_zero m
  rw [negativeHalfBinomial, ← hp]
  calc
    _ = (((-1 : ℚ) ^ m * (-1 : ℚ) ^ m) * ((m.factorial : ℚ) / m.factorial)) *
      (2 ^ m * ((2 * (m : ℚ) + 1) * (ascPochhammer ℚ m).eval (1 / 2)) /
        (((m + 1 : ℕ) : ℚ) * (ascPochhammer ℚ m).eval ((N : ℚ) / 2))) := by ring
    _ = _ := by rw [hs, div_self hf, one_mul, one_mul]

/-- The opening phase and linear marker, from the spherical moment
formula in the draft. -/
theorem simplePole_moments (trace : A →ₗ[K] K) (U : Aˣ) (x B : A)
    (N m : ℕ) (hN : 0 < N) (hm : 0 < m)
    (weight : ∀ k i : ℕ, k ≤ m → i ≠ 2 * k →
      trace ((↑(U⁻¹) : A) ^ m * x ^ i * B ^ (m - k)) = 0)
    (height : ∀ k : ℕ, k ≤ m →
      trace ((↑(U⁻¹) : A) ^ m * x ^ (2 * k) * B ^ (m - k)) =
        algebraMap ℚ K (simplePoleGaussianMoment N m k)) :
    trace (simplePolePhase U x B ^ m) = 0 ∧
      trace (x * simplePolePhase U x B ^ m) =
        algebraMap ℚ K (2 ^ m * (ascPochhammer ℚ m).eval (3 / 2) /
          (((m + 1 : ℕ) : ℚ) * (ascPochhammer ℚ m).eval ((N : ℚ) / 2))) := by
  have he := linearPhaseEndpoint_coefficients m hm
  have hzero := simplePole_moment_coeff trace U x B N m 0 hN weight height
  have hone := simplePole_moment_coeff trace U x B N m 1 hN weight height
  rw [pow_zero, one_mul, add_zero, he.1, mul_zero, map_zero] at hzero
  rw [pow_one, he.2, simplePole_marked_factor] at hone
  exact ⟨hzero, hone⟩

private theorem trace_beta_height (trace : A →ₗ[K] K) (S : A) (M k j : ℕ) (hM : 0 < M)
    (height : ∀ r : ℕ, trace (S ^ r) = (M : K) / (↑(M + 2 * r) : K)) :
    trace (S ^ j * (1 - S) ^ k) = algebraMap ℚ K (mixedBetaMoment M k j) := by
  rw [sub_eq_add_neg, add_comm 1, add_pow, Finset.mul_sum, map_sum,
    mixedBetaMoment, powerBetaMoment, powerBetaAntiderivative, eval_mul, eval_pow, eval_X,
    one_pow, one_mul, powerBetaQuotient, eval_finsetSum, Finset.mul_sum, map_sum]
  refine Finset.sum_congr rfl fun r _ ↦ ?_
  have he : S ^ j * ((-S) ^ r * 1 ^ (k - r) * (k.choose r : A)) =
      algebraMap K A (algebraMap ℚ K ((k.choose r : ℚ) * (-1) ^ r)) * S ^ (j + r) := by
    rw [one_pow, mul_one, neg_pow, pow_add S]
    simp only [map_mul, map_natCast, map_neg, map_one, map_pow]
    ring
  rw [he, ← Algebra.smul_def, map_smul, smul_eq_mul, height, eval_mul, eval_C, eval_pow, eval_X,
    one_pow, mul_one, powerBetaCoefficient]
  simp only [map_mul, map_div₀, map_natCast, map_pow, map_neg, map_one]
  have hn : (↑(M + 2 * j) : K) ≠ 0 := by exact_mod_cast (show M + 2 * j ≠ 0 by omega)
  rw [show M + 2 * (j + r) = M + 2 * j + 2 * r by omega]
  field_simp

section Coordinates

open MvPolynomial

variable {M : ℕ} (a : Fin 2 → K) (b : Fin M → K)
local notation "D" => NormComplement a b

/-- The joint spherical moments are forced by the rank-two rotations.
The proof uses the beta moment of the plane norm as its initial value. -/
theorem trace_coordinate_plane_joint (i : Fin M) (hb : b i = 1)
    {trace : D →ₗ[K] K} (h : GeometricTrace a b trace) (k r : ℕ) :
    let U := normUnit a b
    let x := algebraMap (QuadPolynomial K M) D (X (Sum.inr i))
    let B := algebraMap (QuadPolynomial K M) D (planeNorm a)
    trace ((↑(U⁻¹) : D) ^ (k + r) * x ^ (2 * k) * B ^ r) =
      algebraMap ℚ K ((ascPochhammer ℚ k).eval (1 / 2) * (r.factorial : ℚ) /
        (ascPochhammer ℚ (k + r)).eval ((↑(M + 2 : ℕ) : ℚ) / 2)) := by
  intro U x B
  let v : D := ↑(U⁻¹)
  let S := height a b
  have hM : 0 < M := by have := i.isLt; omega
  have ht := h.complementTrace.normalizer a b hM
  have hB : B * v = 1 - S := by
    have he : (B + algebraMap (QuadPolynomial K M) D (heightNorm b)) * v = 1 := by
      change (algebraMap (QuadPolynomial K M) D (planeNorm a) +
        algebraMap (QuadPolynomial K M) D (heightNorm b)) *
          (↑((normUnit a b)⁻¹) : D) = 1
      rw [← map_add]
      change (algebraMap (QuadPolynomial K M) D (quadNorm a b)) * _ = _
      rw [← normUnit_val, Units.mul_inv]
    rw [add_mul] at he
    exact eq_sub_of_add_eq he
  induction k generalizing r with
  | zero =>
    rw [zero_add, mul_zero, pow_zero, mul_one, mul_comm, ← mul_pow, hB]
    have he := trace_beta_height trace S M r 0 hM ht.height
    rw [pow_zero, one_mul, mixedBetaMoment, mul_zero, add_zero] at he
    have hn : (M : ℚ) ≠ 0 := by positivity
    rw [div_self hn, one_mul, powerBetaMoment_closed_form M r hM] at he
    simpa [show (1 : ℚ) + (M : ℚ) / 2 = (↑(M + 2 : ℕ) : ℚ) / 2 by push_cast; ring] using he
  | succ k ih =>
    have hc := coordinate_plane_recurrence_of_rotations a b trace h.complementTrace.rotations
      i hb k r (k + 1 + r)
    change (2 * (r : K) + 2) * trace (v ^ (k + 1 + r) * x ^ (2 * k + 2) * B ^ r) =
      (2 * (k : K) + 1) * trace (v ^ (k + 1 + r) * x ^ (2 * k) * B ^ (r + 1)) at hc
    have he : k + 1 + r = k + (r + 1) := by omega
    rw [he, ih (r + 1), ← he, show r + 1 = r.succ by rfl, Nat.factorial_succ] at hc
    have hne : (2 * (r : K) + 2) ≠ 0 := by exact_mod_cast (show 2 * r + 2 ≠ 0 by omega)
    change trace (v ^ (k + 1 + r) * x ^ (2 * (k + 1)) * B ^ r) = _
    apply (mul_right_inj' hne).mp
    rw [show 2 * (k + 1) = 2 * k + 2 by omega, hc, ascPochhammer_succ_eval]
    simp only [Nat.cast_mul, Nat.cast_succ]
    have h₁ : (2 * (k : K) + 1) = algebraMap ℚ K (2 * (k : ℚ) + 1) := by simp
    have h₂ : (2 * (r : K) + 2) = algebraMap ℚ K (2 * (r : ℚ) + 2) := by simp
    rw [h₁, h₂, ← map_mul, ← map_mul]
    congr 1
    ring

/-- The new simple-pole proposition on the normalized diagonal chart. -/
theorem draft_simplePole_moments (i : Fin M) (hb : b i = 1)
    {trace : D →ₗ[K] K} (h : GeometricTrace a b trace) (m : ℕ) (hm : 0 < m) :
    let U := normUnit a b
    let x := algebraMap (QuadPolynomial K M) D (X (Sum.inr i))
    let B := algebraMap (QuadPolynomial K M) D (planeNorm a)
    trace (simplePolePhase U x B ^ m) = 0 ∧
      trace (x * simplePolePhase U x B ^ m) =
        algebraMap ℚ K (2 ^ m * (ascPochhammer ℚ m).eval (3 / 2) /
          (((m + 1 : ℕ) : ℚ) * (ascPochhammer ℚ m).eval ((↑(M + 2 : ℕ) : ℚ) / 2))) ∧
      trace (x * simplePolePhase U x B ^ m) ≠ 0 := by
  intro U x B
  have hweight (k n : ℕ) (hk : k ≤ m) (hi : n ≠ 2 * k) :
      trace ((↑(U⁻¹) : D) ^ m * x ^ n * B ^ (m - k)) = 0 :=
    trace_coordinate_plane_weight a b trace h.complementTrace.scalar _ m n (m - k) (by omega)
  have hgauss (k : ℕ) (hk : k ≤ m) :
      trace ((↑(U⁻¹) : D) ^ m * x ^ (2 * k) * B ^ (m - k)) =
        algebraMap ℚ K (simplePoleGaussianMoment (M + 2) m k) := by
    have he := trace_coordinate_plane_joint a b i hb h k (m - k)
    simpa only [Nat.add_sub_of_le hk, simplePoleGaussianMoment] using he
  obtain ⟨h₁, h₂⟩ := simplePole_moments trace U x B (M + 2) m (by omega) hm hweight hgauss
  refine ⟨h₁, h₂, h₂.trans_ne ?_⟩
  apply (map_ne_zero_iff (algebraMap ℚ K) (RingHom.injective (algebraMap ℚ K))).mpr
  exact (div_pos (mul_pos (by positivity) (ascPochhammer_pos _ _ (by norm_num)))
    (mul_pos (by positivity) (ascPochhammer_pos _ _ (by positivity)))).ne'

/-- The simple-pole proposition on the chosen diagonal chart,
including the fixed linear marker and the geometric nonpower property. -/
theorem draft_simplePole_witness (i : Fin M) (hb : b i = 1)
    {trace : D →ₗ[K] K} (h : GeometricTrace a b trace) :
    let U := normUnit a b
    let x := algebraMap (QuadPolynomial K M) D (X (Sum.inr i))
    let B := algebraMap (QuadPolynomial K M) D (planeNorm a)
    let p := simplePolePhase U x B
    NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := AlgebraicClosure K) p) ∧
      AllOrderWitness trace p x ∧
      ∀ n : ℕ, ClearsNormDenominators (quadNorm a b) p n ↔ 1 ≤ n := by
  intro U x B p
  have hb' : b i ≠ 0 := by rw [hb]; exact one_ne_zero
  exact ⟨simplePolePhase_absolute_nonpower a b i hb',
    ⟨fun m hm ↦ (draft_simplePole_moments a b i hb h m hm).1,
      fun m hm ↦ (draft_simplePole_moments a b i hb h m hm).2.2⟩,
    simplePolePhase_pole_order a b i hb'⟩

open MvPolynomial

omit [CharZero K] in
private theorem coeff_mul_vanishes_below (f g : QuadPolynomial K M) (r s : ℕ)
    (hf : ∀ d : (Fin 2 ⊕ Fin M) →₀ ℕ, d.sum (fun _ k ↦ k) < r → coeff d f = 0)
    (hg : ∀ d : (Fin 2 ⊕ Fin M) →₀ ℕ, d.sum (fun _ k ↦ k) < s → coeff d g = 0)
    (d : (Fin 2 ⊕ Fin M) →₀ ℕ) (hd : d.sum (fun _ k ↦ k) < r + s) :
    coeff d (f * g) = 0 := by
  rw [MvPolynomial.coeff_mul]
  apply Finset.sum_eq_zero
  intro e he
  have hs := Finset.HasAntidiagonal.mem_antidiagonal.mp he
  rw [← hs, Finsupp.sum_add_index' (fun _ ↦ rfl) (fun _ _ _ ↦ rfl)] at hd
  by_cases h : e.1.sum (fun _ k ↦ k) < r
  · rw [hf _ h, zero_mul]
  · rw [hg _ (by omega), mul_zero]

omit [CharZero K] in
private theorem coeff_pow_vanishes_below (f : QuadPolynomial K M) (hf : coeff 0 f = 0)
    (m : ℕ) (d : (Fin 2 ⊕ Fin M) →₀ ℕ) (hd : d.sum (fun _ k ↦ k) < m) :
    coeff d (f ^ m) = 0 := by
  induction m generalizing d with
  | zero => omega
  | succ m ih =>
    rw [pow_succ]
    apply coeff_mul_vanishes_below _ _ m 1 ih _ d hd
    intro e he
    have h : e = 0 := by
      by_contra hn
      have hpos := Finsupp.sum_pos (g := fun _ k ↦ k)
        (fun _ hi ↦ Nat.pos_of_ne_zero (Finsupp.mem_support_iff.mp hi)) hn
      omega
    simpa only [h] using hf

/-- For a pure polynomial phase, a multiplier with denominator q^s
has zero marked moment at every order m > 2s. -/
theorem trace_regular_mul_polynomial_phase (trace : D →ₗ[K] K)
    (normalized : trace 1 = 1)
    (inv : ∀ f, trace (scalarPullback a b (Units.mk0 (2 : K) two_ne_zero) f) = trace f)
    (f g : QuadPolynomial K M) (hf : Pure trace (algebraMap _ D f))
    (s m : ℕ) (hm : 2 * s < m) :
    trace (((↑((normUnit a b)⁻¹) : D) ^ s * algebraMap _ D g) *
      (algebraMap _ D f) ^ m) = 0 := by
  have hf₀ : coeff 0 f = 0 := by
    rw [← trace_polynomial_eq_constantCoeff a b trace normalized inv]
    simpa using hf 1 Nat.one_pos
  have hp (d : (Fin 2 ⊕ Fin M) →₀ ℕ) (hd : d.sum (fun _ k ↦ k) < m) :
      coeff d (g * f ^ m) = 0 :=
    coeff_mul_vanishes_below g (f ^ m) 0 m (fun _ h ↦ by omega)
      (coeff_pow_vanishes_below f hf₀ m) d (by simpa using hd)
  rw [mul_assoc, ← map_pow, ← map_mul, (g * f ^ m).as_sum,
    map_sum, Finset.mul_sum, map_sum]
  apply Finset.sum_eq_zero
  intro d hd
  apply trace_monomial_weight a b trace inv
  intro h
  exact notMem_support_iff.mpr (hp d (h ▸ hm)) hd

/-- A normalized scalar-invariant trace makes every pure polynomial phase safe. -/
theorem safe_polynomial_phase (trace : D →ₗ[K] K) (normalized : trace 1 = 1)
    (inv : ∀ f, trace (scalarPullback a b (Units.mk0 (2 : K) two_ne_zero) f) = trace f)
    (f : QuadPolynomial K M) (hf : Pure trace (algebraMap _ D f)) :
    Safe trace (algebraMap _ D f) := by
  intro h
  let g := (IsLocalization.Away.sec (quadNorm a b) h).1
  let s := (IsLocalization.Away.sec (quadNorm a b) h).2
  let U := normUnit a b
  have hp : h * (U : D) ^ s = algebraMap _ D g := by
    rw [normUnit_val, ← map_pow]
    simpa using IsLocalization.Away.sec_spec (quadNorm a b) h
  have hclear : (↑(U⁻¹) : D) ^ s * algebraMap _ D g = h := by
    rw [← hp, mul_left_comm, ← mul_pow, Units.inv_mul, one_pow, mul_one]
  refine ⟨2 * s + 1, fun m hm ↦ ?_⟩
  rw [← hclear]
  exact trace_regular_mul_polynomial_phase a b trace normalized inv f g hf s m (by omega)

/-- Every all-order witness requires a positive norm denominator order
and a polynomial marker of positive degree whenever the marker is polynomial. -/
theorem draft_simplePole_lower_bounds {trace : D →ₗ[K] K}
    (h : GeometricTrace a b trace) {p g : D} (hw : AllOrderWitness trace p g) :
    (∀ s : ℕ, ClearsNormDenominators (quadNorm a b) p s → 1 ≤ s) ∧
      ∀ r : QuadPolynomial K M, algebraMap _ D r = g → 1 ≤ r.totalDegree := by
  constructor
  · intro s hs
    by_contra hn
    have hs₀ : s = 0 := by omega
    obtain ⟨f, hf⟩ : p ∈ Set.range (algebraMap (QuadPolynomial K M) D) := by
      simpa [hs₀, ClearsNormDenominators] using hs
    obtain ⟨n, hn⟩ := safe_polynomial_phase a b trace h.complementTrace.normalized
      h.complementTrace.scalar f (hf ▸ hw.1) g
    exact hw.2 (n + 1) (by omega) (hf ▸ hn (n + 1) (by omega))
  · intro r hr
    by_contra hn
    have hz : r.totalDegree = 0 := by omega
    have hc := totalDegree_eq_zero_iff_eq_C.mp hz
    have hmark := hw.2 1 Nat.one_pos
    rw [← hr, hc, ← algebraMap_eq, ← IsScalarTower.algebraMap_apply, ← Algebra.smul_def,
      map_smul, hw.1 1 Nat.one_pos, smul_zero] at hmark
    exact hmark rfl

end Coordinates

open MvPolynomial in
/-- The current opening display for a ternary sum of squares. -/
theorem draft_simplePole_ternary_opening (a : Fin 2 → K)
    {trace : NormComplement a (fun _ : Fin 1 ↦ (1 : K)) →ₗ[K] K}
    (h : GeometricTrace a (fun _ : Fin 1 ↦ (1 : K)) trace) (m : ℕ) (hm : 0 < m) :
    let b : Fin 1 → K := fun _ ↦ 1
    let U := normUnit a b
    let x := algebraMap (QuadPolynomial K 1) (NormComplement a b) (X (Sum.inr 0))
    let B := algebraMap (QuadPolynomial K 1) (NormComplement a b) (planeNorm a)
    trace (simplePolePhase U x B ^ m) = 0 ∧
      trace (x * simplePolePhase U x B ^ m) = algebraMap ℚ K (2 ^ m / (m + 1 : ℕ)) ∧
      trace (x * simplePolePhase U x B ^ m) ≠ 0 := by
  intro b U x B
  obtain ⟨h₁, h₂, h₃⟩ := draft_simplePole_moments a b 0 rfl h m hm
  refine ⟨h₁, h₂.trans ?_, h₃⟩
  have he : (ascPochhammer ℚ m).eval (3 / 2) ≠ 0 :=
    (ascPochhammer_pos _ _ (by norm_num)).ne'
  congr 1
  norm_num
  rw [mul_div_mul_right _ _ he]

end MathieuQuadratic
