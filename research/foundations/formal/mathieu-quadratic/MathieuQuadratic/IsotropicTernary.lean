import MathieuQuadratic.Draft
import Mathlib.Algebra.MvPolynomial.PDeriv

/-!
# The linear marker in the isotropic ternary example

The computation uses a Laurent model for the two coordinates x and y/q.
The comparison uses only their scalar-rotation weights and the recurrence
for their equal-exponent moments.
-/

open Polynomial
open scoped BigOperators Interval

namespace MathieuQuadratic

variable {K A : Type*} [Field K] [CharZero K] [CommRing A] [Algebra K A] [Algebra ℚ A]

/-- The draft's ternary phase, with u standing for y/q. -/
def isotropicTernaryPhase (x u : A) : A :=
  (1 + x) ^ 3 * u - (1 + x) * (2 + x) * algebraMap ℚ A (1 / 2)

private abbrev Model := AddMonoidAlgebra ℚ[X] ℤ

private noncomputable def laurentUnitAt (n : ℤ) : Modelˣ where
  val := AddMonoidAlgebra.single n 1
  inv := AddMonoidAlgebra.single (-n) 1
  val_inv := by
    rw [AddMonoidAlgebra.single_mul_single, add_neg_cancel, mul_one, AddMonoidAlgebra.one_def]
  inv_val := by
    rw [AddMonoidAlgebra.single_mul_single, neg_add_cancel, mul_one, AddMonoidAlgebra.one_def]

private noncomputable def laurentUnitMap : Multiplicative ℤ →* Modelˣ where
  toFun n := laurentUnitAt n.toAdd
  map_one' := by
    apply Units.ext
    rfl
  map_mul' a b := by
    apply Units.ext
    change AddMonoidAlgebra.single (a.toAdd + b.toAdd) (1 : ℚ[X]) =
      AddMonoidAlgebra.single a.toAdd 1 * AddMonoidAlgebra.single b.toAdd 1
    rw [AddMonoidAlgebra.single_mul_single, mul_one]

private noncomputable def xUnit : Modelˣ := laurentUnitAt 1

private theorem xUnit_zpow (a : ℤ) :
    (↑(xUnit ^ a) : Model) = AddMonoidAlgebra.single a 1 := by
  have he : (Multiplicative.ofAdd (1 : ℤ)) ^ a = Multiplicative.ofAdd a := by
    change a • (1 : ℤ) = a
    simp
  change (↑((laurentUnitMap (Multiplicative.ofAdd (1 : ℤ))) ^ a) : Model) = _
  rw [← map_zpow, he]
  rfl

private noncomputable def betaFunctional (M : ℕ) : ℚ[X] →ₗ[ℚ] ℚ :=
  Polynomial.lsum (fun j ↦ LinearMap.mulRight ℚ ((M : ℚ) / (↑(M + 2 * j) : ℚ)))

private noncomputable def modelTrace (M : ℕ) : Model →ₗ[ℚ] ℚ where
  toFun p := betaFunctional M (p.coeff 0)
  map_add' p q := by change betaFunctional M (p.coeff 0 + q.coeff 0) = _; rw [map_add]
  map_smul' c p := by change betaFunctional M (c • p.coeff 0) = _; rw [map_smul]; rfl

private noncomputable def embedPoly : ℚ[X] →+* Model := AddMonoidAlgebra.singleZeroRingHom

private noncomputable def modelHeight : Model := embedPoly X

private theorem modelTrace_coeff (M : ℕ) (a : ℤ) (g : ℚ[X]) :
    modelTrace M ((↑(xUnit ^ a) : Model) * embedPoly g) =
      if a = 0 then betaFunctional M g else 0 := by
  rw [xUnit_zpow]
  change betaFunctional M ((AddMonoidAlgebra.single a (1 : ℚ[X]) *
      AddMonoidAlgebra.single 0 g).coeff 0) = _
  rw [AddMonoidAlgebra.single_mul_single, add_zero, one_mul]
  by_cases h : a = 0
  · simp only [h, AddMonoidAlgebra.coeff_single, if_true]
    rw [Finsupp.single_eq_same]
  · rw [AddMonoidAlgebra.coeff_single, Finsupp.single_eq_of_ne (Ne.symm h), if_neg h, map_zero]

private theorem betaFunctional_power (M j : ℕ) :
    betaFunctional M (X ^ j) = (M : ℚ) / (↑(M + 2 * j) : ℚ) := by
  rw [← one_mul (X ^ j), ← C_1, C_mul_X_pow_eq_monomial, betaFunctional,
    lsum_apply, sum_monomial_index]
  · norm_num [LinearMap.mulRight_apply]
  · exact map_zero _

private theorem model_weight (M : ℕ) (a : ℤ) (j : ℕ) (ha : a ≠ 0) :
    modelTrace M ((↑(xUnit ^ a) : Model) * modelHeight ^ j) = 0 := by
  rw [modelHeight, ← map_pow, modelTrace_coeff, if_neg ha]

private theorem model_height (M j : ℕ) :
    modelTrace M (modelHeight ^ j) = (M : ℚ) / (↑(M + 2 * j) : ℚ) := by
  simpa [modelHeight, map_pow] using (modelTrace_coeff M 0 (X ^ j)).trans
    (by simpa using betaFunctional_power M j)

private theorem betaFunctional_beta (M n : ℕ) :
    betaFunctional M ((1 - X) ^ n) = powerBetaMoment M n := by
  rw [sub_eq_add_neg, add_comm, add_pow, map_sum, powerBetaMoment, powerBetaAntiderivative,
    powerBetaQuotient, eval_mul, eval_pow, eval_X, one_pow, one_mul, eval_finsetSum]
  refine Finset.sum_congr rfl fun j _ ↦ ?_
  rw [one_pow, mul_one, neg_pow, mul_comm, ← C_eq_natCast, ← C_1, ← C_neg, ← map_pow,
    ← mul_assoc, ← C_mul, ← smul_eq_C_mul, map_smul, betaFunctional_power,
    eval_mul, eval_C, eval_pow, eval_X, one_pow, mul_one, powerBetaCoefficient]
  simp only [smul_eq_mul]
  ring

private noncomputable def modelU : Model :=
  (↑(xUnit⁻¹) : Model) * embedPoly (1 - X) * algebraMap ℚ Model (1 / 2)

/-- The equal-exponent trace of x and y/q in the isotropic ternary calculation. -/
noncomputable def hyperbolicDiagonalMoment (M n : ℕ) : ℚ :=
  (1 / 2 : ℚ) ^ n * powerBetaMoment M n

noncomputable abbrev isotropicDiagonalMoment (n : ℕ) : ℚ :=
  hyperbolicDiagonalMoment 1 n

private theorem model_pair_moment (M i j : ℕ) :
    modelTrace M ((xUnit : Model) ^ i * modelU ^ j) =
      if i = j then hyperbolicDiagonalMoment M j else 0 := by
  have he : (xUnit : Model) ^ i * modelU ^ j =
      (1 / 2 : ℚ) ^ j • ((↑(xUnit ^ ((i : ℤ) - j)) : Model) * embedPoly ((1 - X) ^ j)) := by
    rw [modelU, mul_pow, mul_pow, ← map_pow, ← map_pow, Algebra.smul_def,
      zpow_sub, zpow_natCast, zpow_natCast, div_eq_mul_inv, ← inv_pow,
      Units.val_mul, Units.val_pow_eq_pow_val, Units.val_pow_eq_pow_val]
    ring
  rw [he, map_smul, smul_eq_mul, modelTrace_coeff, betaFunctional_beta]
  simp [show ((i : ℤ) - j = 0) ↔ i = j by omega, mul_ite, hyperbolicDiagonalMoment]

private theorem phase_in_model :
    isotropicTernaryPhase (xUnit : Model) modelU =
      (1 / 2 : ℚ) • productPhase xUnit modelHeight 1 1 1 := by
  let x : Model := xUnit
  let y : Model := ↑(xUnit⁻¹)
  let Z := modelHeight
  let r : Model := algebraMap ℚ Model (1 / 2)
  have hxy : x * y = 1 := Units.mul_inv _
  have hZ : embedPoly (1 - X) = 1 - modelHeight := by
    rw [map_sub, map_one]
    rfl
  simp only [isotropicTernaryPhase, productPhase, modelU, hZ, pow_one, Algebra.smul_def]
  change (1 + x) ^ 3 * (y * (1 - Z) * r) - (1 + x) * (2 + x) * r =
    r * (y * (1 + x) * (1 - (1 + x) ^ 2 * Z))
  linear_combination (r * (1 + x) * (2 + x)) * hxy

private theorem model_phase_moments (m : ℕ) (hm : 0 < m) :
    modelTrace 1 (isotropicTernaryPhase (xUnit : Model) modelU ^ m) = 0 ∧
      modelTrace 1 ((xUnit : Model) * isotropicTernaryPhase (xUnit : Model) modelU ^ m) =
        isotropicDiagonalMoment m := by
  obtain ⟨hp, ht⟩ := product_uniform_moments (modelTrace 1) xUnit modelHeight 1 1 m
    Nat.one_pos le_rfl hm (model_weight 1) (model_height 1)
  have hbeta : mixedBetaMoment 1 (1 * m) 0 = powerBetaMoment 1 m := by simp [mixedBetaMoment]
  rw [phase_in_model, _root_.smul_pow, map_smul, hp, smul_zero, mul_smul_comm, map_smul, smul_eq_mul]
  simp only [pow_one, hbeta] at ht
  exact ⟨rfl, by rw [ht]; rfl⟩

/-- The equal-exponent recurrence derived from the hyperbolic rotation. -/
theorem isotropicDiagonalMoment_succ (n : ℕ) :
    (2 * (n + 1 : ℕ) + 1 : ℚ) * isotropicDiagonalMoment (n + 1) =
      (n + 1 : ℕ) * isotropicDiagonalMoment n := by
  have h := powerBetaMoment_succ 1 n Nat.one_pos
  simp only [isotropicDiagonalMoment, hyperbolicDiagonalMoment, pow_succ]
  push_cast at *
  linear_combination ((1 / 2 : ℚ) ^ n / 2) * h

theorem isotropicDiagonalMoment_cast_eq_integral (n : ℕ) :
    (isotropicDiagonalMoment n : ℝ) =
      (1 / 2 : ℝ) ^ n * ∫ w in (0 : ℝ)..1, (1 - w ^ 2) ^ n := by
  simpa [isotropicDiagonalMoment, hyperbolicDiagonalMoment] using
    congrArg ((1 / 2 : ℝ) ^ n * ·) (powerBetaMoment_cast_eq_integral 1 n Nat.one_pos)

private theorem trace_aeval_two_monomials (trace : A →ₗ[K] K) (x u : A) (M : ℕ)
    (compat : ∀ c : ℚ, algebraMap ℚ A c = algebraMap K A (algebraMap ℚ K c))
    (moments : ∀ i j : ℕ, trace (x ^ i * u ^ j) = if i = j then
      algebraMap ℚ K (hyperbolicDiagonalMoment M j) else 0)
    (p : MvPolynomial (Fin 2) ℚ) :
    trace (MvPolynomial.aeval ![x, u] p) =
      algebraMap ℚ K (modelTrace M (MvPolynomial.aeval ![(xUnit : Model), modelU] p)) := by
  induction p using MvPolynomial.induction_on' with
  | add p q hp hq =>
      simp only [map_add, hp, hq]
  | monomial d c =>
      simp only [MvPolynomial.aeval_monomial, Finsupp.prod_fintype _ (fun _ k ↦ _ ^ k) (fun _ ↦ pow_zero _),
        Fin.prod_univ_two, Matrix.cons_val_zero, Matrix.cons_val_one, compat,
        ← Algebra.smul_def, map_smul, smul_eq_mul]
      rw [moments, model_pair_moment M]
      split_ifs <;> simp [Algebra.smul_def]

/-- The two moments in the isotropic example, once the geometric trace
satisfies the two-coordinate weight and rotation identities. -/
theorem isotropicTernary_moments (trace : A →ₗ[K] K) (x u : A)
    (compat : ∀ c : ℚ, algebraMap ℚ A c = algebraMap K A (algebraMap ℚ K c))
    (normalized : trace 1 = 1)
    (weight : ∀ i j : ℕ, i ≠ j → trace (x ^ i * u ^ j) = 0)
    (rotation : ∀ n : ℕ, (2 * (n + 1 : ℕ) + 1 : K) * trace ((x * u) ^ (n + 1)) =
      (n + 1 : ℕ) * trace ((x * u) ^ n))
    (m : ℕ) (hm : 0 < m) :
    trace (isotropicTernaryPhase x u ^ m) = 0 ∧
      trace (x * isotropicTernaryPhase x u ^ m) =
        algebraMap ℚ K (isotropicDiagonalMoment m) ∧
      trace (x * isotropicTernaryPhase x u ^ m) ≠ 0 := by
  have he (n : ℕ) : trace ((x * u) ^ n) = algebraMap ℚ K (isotropicDiagonalMoment n) := by
    induction n with
    | zero => simp [normalized, isotropicDiagonalMoment, hyperbolicDiagonalMoment, powerBetaMoment, powerBetaAntiderivative,
        powerBetaQuotient, powerBetaCoefficient]
    | succ n ih =>
      have hn : (2 * (n + 1 : ℕ) + 1 : K) ≠ 0 := by exact_mod_cast (show 2 * (n + 1) + 1 ≠ 0 by omega)
      apply (mul_right_inj' hn).mp
      rw [rotation, ih]
      have hc := congrArg (algebraMap ℚ K) (isotropicDiagonalMoment_succ n)
      simpa using hc.symm
  have hpair (i j : ℕ) : trace (x ^ i * u ^ j) =
      if i = j then algebraMap ℚ K (isotropicDiagonalMoment j) else 0 := by
    split_ifs with h
    · subst i
      rw [← mul_pow, he]
    · exact weight i j h
  let p : MvPolynomial (Fin 2) ℚ :=
    isotropicTernaryPhase (MvPolynomial.X 0) (MvPolynomial.X 1)
  have hp : MvPolynomial.aeval ![x, u] p = isotropicTernaryPhase x u := by
    simp only [p, isotropicTernaryPhase, map_sub, map_mul, map_pow, map_add, map_one,
      map_ofNat, MvPolynomial.aeval_X, Matrix.cons_val_zero, Matrix.cons_val_one, AlgHom.commutes]
  have hp' : MvPolynomial.aeval ![(xUnit : Model), modelU] p =
      isotropicTernaryPhase (xUnit : Model) modelU := by
    simp only [p, isotropicTernaryPhase, map_sub, map_mul, map_pow, map_add, map_one,
      map_ofNat, MvPolynomial.aeval_X, Matrix.cons_val_zero, Matrix.cons_val_one, AlgHom.commutes]
  have h₁ := trace_aeval_two_monomials trace x u 1 compat hpair (p ^ m)
  have h₂ := trace_aeval_two_monomials trace x u 1 compat hpair (MvPolynomial.X 0 * p ^ m)
  simp only [map_pow, hp, hp', map_mul, MvPolynomial.aeval_X, Matrix.cons_val_zero] at h₁ h₂
  obtain ⟨h₁', h₂'⟩ := model_phase_moments m hm
  rw [h₁', map_zero] at h₁
  rw [h₂'] at h₂
  refine ⟨h₁, h₂, h₂.trans_ne ?_⟩
  apply (map_ne_zero_iff (algebraMap ℚ K) (RingHom.injective (algebraMap ℚ K))).mpr
  exact mul_ne_zero (pow_ne_zero _ (by norm_num)) (powerBetaMoment_pos 1 m Nat.one_pos).ne'

/-- The norm-one polynomial phase on a hyperbolic plane. -/
def normOnePhase (x y : A) (M : ℕ) : A :=
  x * (1 + x) ^ M * (2 * (1 + x) ^ 2 * y - 2 - x) ^ (M + 1)

private theorem normOnePhase_in_model (M : ℕ) :
    normOnePhase (xUnit : Model) modelU M = productPhase xUnit modelHeight M M (M + 1) := by
  let x : Model := xUnit
  let z : Model := ↑(xUnit⁻¹)
  let S := modelHeight
  let c : Model := algebraMap ℚ Model (1 / 2)
  let T := 1 + x
  have hxz : x * z = 1 := Units.mul_inv _
  have hc : 2 * c = 1 := by
    change 2 * algebraMap ℚ Model (1 / 2) = _
    rw [← map_ofNat (algebraMap ℚ Model) 2, ← map_mul]
    norm_num
  have hxy : x * modelU = (1 - S) * c := by
    rw [modelU, map_sub, map_one]
    change x * (z * (1 - S) * c) = _
    rw [← mul_assoc x, ← mul_assoc x z, hxz, one_mul]
  have hf : 2 * T ^ 2 * modelU - 2 - x = z * (1 - T ^ 2 * S) := by
    apply xUnit.isUnit.mul_left_cancel
    change x * (2 * T ^ 2 * modelU - 2 - x) = x * (z * (1 - T ^ 2 * S))
    rw [← mul_assoc x z, hxz, one_mul]
    have he : x * (2 * T ^ 2 * modelU - 2 - x) = 2 * T ^ 2 * (x * modelU) - x * (2 + x) := by ring
    rw [he, hxy]
    calc
      _ = (2 * c) * (T ^ 2 * (1 - S)) - x * (2 + x) := by ring
      _ = _ := by rw [hc, one_mul]; dsimp [T]; ring
  change x * T ^ M * (2 * T ^ 2 * modelU - 2 - x) ^ (M + 1) =
    z ^ M * T ^ M * (1 - T ^ 2 * S) ^ (M + 1)
  rw [hf, mul_pow, pow_succ z]
  calc
    _ = (x * z) * (z ^ M * T ^ M * (1 - T ^ 2 * S) ^ (M + 1)) := by ring
    _ = _ := by rw [hxz, one_mul]

private theorem trace_eval_of_height (trace : A →ₗ[K] K) (S : A) (M : ℕ)
    (compat : ∀ c : ℚ, algebraMap ℚ A c = algebraMap K A (algebraMap ℚ K c))
    (height : ∀ j : ℕ, trace (S ^ j) = (M : K) / (↑(M + 2 * j) : K))
    (p : ℚ[X]) :
    trace (Polynomial.aeval S p) = algebraMap ℚ K (betaFunctional M p) := by
  induction p using Polynomial.induction_on' with
  | add p q hp hq => simp only [map_add, hp, hq]
  | monomial d c =>
      rw [Polynomial.aeval_monomial, compat, ← Algebra.smul_def, map_smul, smul_eq_mul, height,
        ← C_mul_X_pow_eq_monomial, ← smul_eq_C_mul, map_smul, smul_eq_mul, betaFunctional_power,
        map_mul, map_div₀, map_natCast, map_natCast]

/-- The norm-one marked moments, from the quadratic relation and the
opposite weights and height moments of the geometric trace. -/
theorem normOnePhase_moments (trace : A →ₗ[K] K) (x y S : A) (M : ℕ) (hM : 0 < M)
    (compat : ∀ c : ℚ, algebraMap ℚ A c = algebraMap K A (algebraMap ℚ K c))
    (relation : x * y = algebraMap ℚ A (1 / 2) * (1 - S))
    (weight : ∀ i j : ℕ, i ≠ j → trace (x ^ i * y ^ j) = 0)
    (height : ∀ j : ℕ, trace (S ^ j) = (M : K) / (↑(M + 2 * j) : K))
    (m : ℕ) (hm : 0 < m) :
    trace (normOnePhase x y M ^ m) = 0 ∧
      trace (x ^ M * normOnePhase x y M ^ m) = algebraMap ℚ K (mixedBetaMoment M ((M + 1) * m) 0) ∧
      trace (x ^ M * normOnePhase x y M ^ m) ≠ 0 := by
  have he (j : ℕ) : trace ((x * y) ^ j) = algebraMap ℚ K (hyperbolicDiagonalMoment M j) := by
    rw [relation, mul_pow, ← map_pow, compat, ← Algebra.smul_def, map_smul, smul_eq_mul]
    have hg := trace_eval_of_height trace S M compat height ((1 - X) ^ j)
    rw [betaFunctional_beta] at hg
    simp only [map_pow, map_sub, map_one, aeval_X] at hg
    rw [hg, hyperbolicDiagonalMoment, map_mul]
  have hpair (i j : ℕ) :
      trace (x ^ i * y ^ j) = if i = j then algebraMap ℚ K (hyperbolicDiagonalMoment M j) else 0 := by
    split_ifs with h
    · subst i
      rw [← mul_pow, he]
    · exact weight i j h
  let p : MvPolynomial (Fin 2) ℚ := normOnePhase (MvPolynomial.X 0) (MvPolynomial.X 1) M
  have hp : MvPolynomial.aeval ![x, y] p = normOnePhase x y M := by
    simp [p, normOnePhase]
  have hp' : MvPolynomial.aeval ![(xUnit : Model), modelU] p =
      productPhase xUnit modelHeight M M (M + 1) := by
    rw [show MvPolynomial.aeval ![(xUnit : Model), modelU] p =
      normOnePhase (xUnit : Model) modelU M by simp [p, normOnePhase], normOnePhase_in_model]
  have h₁ := trace_aeval_two_monomials trace x y M compat hpair (p ^ m)
  have h₂ := trace_aeval_two_monomials trace x y M compat hpair (MvPolynomial.X 0 ^ M * p ^ m)
  simp only [map_pow, hp, hp', map_mul, MvPolynomial.aeval_X, Matrix.cons_val_zero] at h₁ h₂
  obtain ⟨h₁', h₂'⟩ := product_uniform_moments (modelTrace M) xUnit modelHeight M (M + 1) m
    hM (by omega) hm (model_weight M) (model_height M)
  rw [h₁', map_zero] at h₁
  rw [h₂'] at h₂
  have hn : algebraMap ℚ K (mixedBetaMoment M ((M + 1) * m) 0) ≠ 0 :=
    (map_ne_zero_iff (algebraMap ℚ K) (RingHom.injective (algebraMap ℚ K))).mpr
      (mixedBetaMoment_pos M _ _ hM).ne'
  exact ⟨h₁, h₂, h₂.trans_ne hn⟩

section HyperbolicCoordinates

open MvPolynomial
open scoped DualNumber TensorProduct

/-- The ternary norm in hyperbolic coordinates. -/
noncomputable def hyperbolicNorm (a : K) : MvPolynomial (Fin 3) K :=
  2 * (MvPolynomial.X (0 : Fin 3) * MvPolynomial.X (1 : Fin 3)) + MvPolynomial.C a * MvPolynomial.X (2 : Fin 3) ^ 2

abbrev HyperbolicComplement (a : K) := Localization.Away (hyperbolicNorm a)

noncomputable def hyperbolicNormUnit (a : K) : (HyperbolicComplement a)ˣ :=
  (IsLocalization.Away.algebraMap_isUnit (S := HyperbolicComplement a) (hyperbolicNorm a)).unit

omit [CharZero K] in
private theorem hyperbolicNormUnit_val (a : K) :
    (hyperbolicNormUnit a : HyperbolicComplement a) = algebraMap _ _ (hyperbolicNorm a) :=
  IsUnit.unit_spec _

noncomputable def hyperbolicX (a : K) : HyperbolicComplement a :=
  algebraMap (MvPolynomial (Fin 3) K) _ (MvPolynomial.X (0 : Fin 3))

noncomputable def hyperbolicYOverNorm (a : K) : HyperbolicComplement a :=
  algebraMap (MvPolynomial (Fin 3) K) _ (MvPolynomial.X (1 : Fin 3)) * (↑((hyperbolicNormUnit a)⁻¹) : HyperbolicComplement a)

/-- The infinitesimal rotation used in the hyperbolic moment recurrence. -/
noncomputable def hyperbolicRotation (a : K) :
    Derivation K (MvPolynomial (Fin 3) K) (MvPolynomial (Fin 3) K) :=
  (MvPolynomial.C a * MvPolynomial.X (2 : Fin 3)) • pderiv (R := K) (1 : Fin 3) - MvPolynomial.X (R := K) (0 : Fin 3) • pderiv (R := K) (2 : Fin 3)

omit [CharZero K] in
private theorem hyperbolicRotation_monomial (a : K) (n : ℕ) :
    hyperbolicRotation a (MvPolynomial.X (0 : Fin 3) ^ n * MvPolynomial.X (1 : Fin 3) ^ (n + 1) * MvPolynomial.X (2 : Fin 3)) =
      ((n + 1 : ℕ) : MvPolynomial (Fin 3) K) * (MvPolynomial.X (0 : Fin 3) * MvPolynomial.X (1 : Fin 3)) ^ n * hyperbolicNorm a -
        (2 * (n + 1 : ℕ) + 1 : MvPolynomial (Fin 3) K) * (MvPolynomial.X (0 : Fin 3) * MvPolynomial.X (1 : Fin 3)) ^ (n + 1) := by
  simp [hyperbolicRotation, hyperbolicNorm, smul_eq_mul, Derivation.leibniz,
    Derivation.leibniz_pow, pderiv_X_of_ne (by decide : (0 : Fin 3) ≠ 1),
    pderiv_X_of_ne (by decide : (0 : Fin 3) ≠ 2),
    pderiv_X_of_ne (by decide : (1 : Fin 3) ≠ 2),
    pderiv_X_of_ne (by decide : (2 : Fin 3) ≠ 1)]
  rw [mul_pow, mul_pow, pow_succ (MvPolynomial.X (0 : Fin 3) : MvPolynomial (Fin 3) K) n,
    pow_succ (MvPolynomial.X (1 : Fin 3) : MvPolynomial (Fin 3) K) n]
  ring

omit [CharZero K] in
private theorem hyperbolicRotation_recurrence (a : K)
    (trace : HyperbolicComplement a →ₗ[K] K)
    (invariance : ∀ (p : MvPolynomial (Fin 3) K) (r : ℕ),
      trace (algebraMap (MvPolynomial (Fin 3) K) _ (hyperbolicRotation a p) *
        (↑((hyperbolicNormUnit a)⁻¹) : HyperbolicComplement a) ^ r) = 0) (n : ℕ) :
    (2 * (n + 1 : ℕ) + 1 : K) * trace ((hyperbolicX a * hyperbolicYOverNorm a) ^ (n + 1)) =
      (n + 1 : ℕ) * trace ((hyperbolicX a * hyperbolicYOverNorm a) ^ n) := by
  let R := MvPolynomial (Fin 3) K
  let U := hyperbolicNormUnit a
  let φ := algebraMap R (HyperbolicComplement a)
  have h := invariance ((MvPolynomial.X (0 : Fin 3) : R) ^ n * MvPolynomial.X (1 : Fin 3) ^ (n + 1) * MvPolynomial.X (2 : Fin 3)) (n + 1)
  rw [hyperbolicRotation_monomial] at h
  simp only [map_sub, map_mul, map_pow, map_add, map_natCast, map_one, map_ofNat] at h
  change trace ((((n + 1 : ℕ) : HyperbolicComplement a) * (φ (MvPolynomial.X (0 : Fin 3)) * φ (MvPolynomial.X (1 : Fin 3))) ^ n *
      φ (hyperbolicNorm a) -
        (2 * ((n + 1 : ℕ) : HyperbolicComplement a) + 1) * (φ (MvPolynomial.X (0 : Fin 3)) * φ (MvPolynomial.X (1 : Fin 3))) ^ (n + 1)) *
      (↑(U⁻¹) : HyperbolicComplement a) ^ (n + 1)) = 0 at h
  have hcancel : φ (hyperbolicNorm a) * (↑(U⁻¹) : HyperbolicComplement a) ^ (n + 1) =
      (↑(U⁻¹) : HyperbolicComplement a) ^ n := by
    rw [← hyperbolicNormUnit_val, pow_succ', ← mul_assoc, Units.mul_inv, one_mul]
  rw [sub_mul, mul_assoc _ (φ (hyperbolicNorm a)), hcancel, map_sub, sub_eq_zero] at h
  let x := hyperbolicX a
  let u := hyperbolicYOverNorm a
  have hpow (r : ℕ) : (φ (MvPolynomial.X (0 : Fin 3)) * φ (MvPolynomial.X (1 : Fin 3))) ^ r * (↑(U⁻¹) : HyperbolicComplement a) ^ r =
      (x * u) ^ r := by
    change _ = (φ (MvPolynomial.X (0 : Fin 3)) * (φ (MvPolynomial.X (1 : Fin 3)) * (↑(U⁻¹) : HyperbolicComplement a))) ^ r
    rw [← mul_pow, mul_assoc]
  rw [mul_assoc, hpow, mul_assoc, hpow] at h
  have hc (c : K) (p : HyperbolicComplement a) :
      trace (algebraMap K _ c * p) = c * trace p := by
    rw [← Algebra.smul_def, map_smul, smul_eq_mul]
  rw [← map_natCast (algebraMap K (HyperbolicComplement a)),
    ← map_one (algebraMap K (HyperbolicComplement a)),
    ← map_ofNat (algebraMap K (HyperbolicComplement a)) 2,
    ← map_mul, ← map_add, hc, hc] at h
  exact h.symm

/-- The hyperbolic chart phase with the rational coefficient supplied by K. -/
noncomputable def hyperbolicPhase (a : K) : HyperbolicComplement a :=
  letI : Algebra ℚ (HyperbolicComplement a) :=
    Algebra.compHom _ (algebraMap ℚ K)
  isotropicTernaryPhase (hyperbolicX a) (hyperbolicYOverNorm a)

/-- The draft's linear marker with its exact moments on the hyperbolic
chart. The inputs are the torus character and rotation trace identities. -/
theorem draft_isotropic_linear_moments (a : K)
    (trace : HyperbolicComplement a →ₗ[K] K) (normalized : trace 1 = 1)
    (torus : ∀ i j : ℕ, trace ((2 : K) ^ ((i : ℤ) - j) •
      (hyperbolicX a ^ i * hyperbolicYOverNorm a ^ j)) =
        trace (hyperbolicX a ^ i * hyperbolicYOverNorm a ^ j))
    (rotation : ∀ (p : MvPolynomial (Fin 3) K) (r : ℕ),
      trace (algebraMap (MvPolynomial (Fin 3) K) _ (hyperbolicRotation a p) *
        (↑((hyperbolicNormUnit a)⁻¹) : HyperbolicComplement a) ^ r) = 0) (m : ℕ) (hm : 0 < m) :
    trace (hyperbolicPhase a ^ m) = 0 ∧
      trace (hyperbolicX a * hyperbolicPhase a ^ m) =
        algebraMap ℚ K (isotropicDiagonalMoment m) ∧
      trace (hyperbolicX a * hyperbolicPhase a ^ m) ≠ 0 := by
  letI : Algebra ℚ (HyperbolicComplement a) := Algebra.compHom _ (algebraMap ℚ K)
  refine isotropicTernary_moments trace _ _ (fun _ ↦ rfl) normalized ?_
    (hyperbolicRotation_recurrence a trace rotation) m hm
  intro i j hij
  have h := torus i j
  rw [map_smul, smul_eq_mul] at h
  have he : (2 : K) ^ ((i : ℤ) - j) ≠ 1 := by
    intro hz
    have hh : algebraMap ℚ K ((2 : ℚ) ^ ((i : ℤ) - j)) = algebraMap ℚ K 1 := by simpa using hz
    have hi := (RingHom.injective (algebraMap ℚ K)) hh
    have hz₀ := (zpow_right_strictMono₀ (by norm_num : (1 : ℚ) < 2)).injective
      (show (2 : ℚ) ^ ((i : ℤ) - j) = (2 : ℚ) ^ (0 : ℤ) by simpa using hi)
    omega
  have hz : ((2 : K) ^ ((i : ℤ) - j) - 1) *
      trace (hyperbolicX a ^ i * hyperbolicYOverNorm a ^ j) = 0 := by
    rw [sub_mul, one_mul, h, sub_self]
  exact (mul_eq_zero.mp hz).resolve_left (sub_ne_zero.mpr he)

/-- The restriction to y = 1 and z = 0 has a simple zero at x = −1. -/
theorem hyperbolicPhase_absolute_nonpower (a : K) :
    NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := AlgebraicClosure K)
      (hyperbolicPhase a)) := by
  let L := AlgebraicClosure K
  let R := MvPolynomial (Fin 3) K
  let D := HyperbolicComplement a
  letI : Algebra ℚ D := Algebra.compHom _ (algebraMap ℚ K)
  let x₀ : L[ε] := TrivSqZeroExt.inl (-1) + TrivSqZeroExt.inr 1
  let e : R →ₐ[K] L[ε] := MvPolynomial.aeval ![x₀, 1, 0]
  have hEval (i : Fin 3) : e (MvPolynomial.X i) = (![x₀, 1, 0] : Fin 3 → L[ε]) i :=
    MvPolynomial.aeval_X _ _
  have hq : e (hyperbolicNorm a) = 2 * x₀ := by
    simp [hyperbolicNorm, hEval]
    rw [map_ofNat]
  have htwo : (2 : L[ε]).fst = 2 :=
    map_ofNat (TrivSqZeroExt.fstHom L L L) 2
  have hunit : IsUnit (e (hyperbolicNorm a)) := by
    rw [hq, TrivSqZeroExt.isUnit_iff_isUnit_fst]
    apply isUnit_iff_ne_zero.2
    simp [x₀, TrivSqZeroExt.fst_mul, htwo]
  let jet : D →ₐ[K] L[ε] := IsLocalization.Away.liftAlgHom (hyperbolicNorm a) hunit
  have hmap (p : R) : jet (algebraMap R D p) = e p := IsLocalization.Away.lift_eq _ hunit _
  have hx : jet (hyperbolicX a) = x₀ := by
    rw [hyperbolicX, hmap]
    simp [hEval]
  have hy : jet (hyperbolicYOverNorm a) *
      (2 * x₀) = 1 := by
    let U := hyperbolicNormUnit a
    have he : jet (U : D) = 2 * x₀ := by rw [hyperbolicNormUnit_val, hmap, hq]
    have hyy : e (MvPolynomial.X (1 : Fin 3)) = 1 := by simp [hEval]
    rw [← he, hyperbolicYOverNorm, map_mul, hmap, hyy, one_mul, ← map_mul, Units.inv_mul, map_one]
  let c := jet (algebraMap ℚ D (1 / 2))
  have hc : 2 * c = 1 := by
    change 2 * jet (algebraMap ℚ D (1 / 2)) = _
    rw [← map_ofNat jet 2, ← map_mul, ← map_ofNat (algebraMap ℚ D) 2, ← map_mul]
    norm_num
  have he : 2 * x₀ * jet (hyperbolicPhase a) = 1 + x₀ := by
    change 2 * x₀ * jet (isotropicTernaryPhase (hyperbolicX a) (hyperbolicYOverNorm a)) = _
    rw [isotropicTernaryPhase, map_sub, map_mul, map_mul, map_mul, map_pow, map_add, map_add,
      map_one, map_ofNat, hx]
    change 2 * x₀ * ((1 + x₀) ^ 3 * jet (hyperbolicYOverNorm a) -
      (1 + x₀) * (2 + x₀) * c) = _
    linear_combination (1 + x₀) ^ 3 * hy - x₀ * (1 + x₀) * (2 + x₀) * hc
  have hfst : (jet (hyperbolicPhase a)).fst = 0 := by
    have h := congrArg (fun w : L[ε] ↦ w.fst) he
    simpa [x₀, TrivSqZeroExt.fst_mul, htwo] using h
  have hsnd : (jet (hyperbolicPhase a)).snd ≠ 0 := by
    intro hz
    have h := congrArg (fun w : L[ε] ↦ w.snd) he
    simp [x₀, hfst, hz] at h
  let extJet : L ⊗[K] D →ₐ[L] L[ε] :=
    Algebra.TensorProduct.lift (Algebra.ofId L L[ε]) jet (fun _ _ ↦ Commute.all _ _)
  have hex : extJet (Algebra.TensorProduct.includeRight (hyperbolicPhase a)) = jet (hyperbolicPhase a) := by
    dsimp only [extJet]
    rw [Algebra.TensorProduct.includeRight_apply, Algebra.TensorProduct.lift_tmul, map_one, one_mul]
  apply nonpower_of_dualNumber extJet.toRingHom
  · change (extJet (Algebra.TensorProduct.includeRight _)).fst = 0
    rwa [hex]
  · change (extJet (Algebra.TensorProduct.includeRight _)).snd ≠ 0
    rwa [hex]

/-- The isotropic ternary example has least detecting degree one,
given its all-order linear witness on the hyperbolic chart. -/
theorem draft_isotropic_least_degree (a : K) (trace : HyperbolicComplement a →ₗ[K] K)
    (hw : AllOrderWitness trace (hyperbolicPhase a) (hyperbolicX a)) (D : ℕ) :
    (∃ p : MvPolynomial (Fin 3) K, p.totalDegree ≤ D ∧ ¬ ∃ m₀ : ℕ, ∀ m ≥ m₀,
      trace (algebraMap (MvPolynomial (Fin 3) K) _ p * hyperbolicPhase a ^ m) = 0) ↔ 1 ≤ D := by
  constructor
  · rintro ⟨p, hp, h⟩
    by_contra! hd
    have hconst := MvPolynomial.totalDegree_eq_zero_iff_eq_C.mp (by omega : p.totalDegree = 0)
    apply h
    refine ⟨1, fun m hm ↦ ?_⟩
    rw [hconst, ← MvPolynomial.algebraMap_eq, ← IsScalarTower.algebraMap_apply K
      (MvPolynomial (Fin 3) K) (HyperbolicComplement a),
      ← Algebra.smul_def, map_smul, hw.1 m hm, smul_zero]
  · intro hD
    refine ⟨MvPolynomial.X 0, by simpa [MvPolynomial.totalDegree_X] using hD, ?_⟩
    rintro ⟨m₀, h⟩
    exact hw.2 (m₀ + 1) (by omega) (h (m₀ + 1) (by omega))

end HyperbolicCoordinates

end MathieuQuadratic
