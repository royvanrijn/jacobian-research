import MathieuQuadratic.OneVariableLaurent
import MathieuQuadratic.Uniform
import Mathlib.Algebra.Polynomial.Derivative

/-!
# Angular endpoints of eventually vanishing phases

The angular integrations in the draft are coefficient identities over
the original characteristic-zero field. The Laurent variable is the
norm, and the polynomial variable is the height ratio.
-/

open Polynomial

namespace MathieuQuadratic

variable {K : Type*} [Field K] [CharZero K]

local notation "L" => AddMonoidAlgebra K ℤ

private noncomputable def angularIntegral (w : ℕ → K) : L[X] →ₗ[K] L :=
  Polynomial.lsum fun n ↦ w n • LinearMap.id

omit [CharZero K] in
@[simp] private theorem angularIntegral_monomial (w : ℕ → K) (n : ℕ) (c : L) :
    angularIntegral w (monomial n c) = w n • c := by
  rw [angularIntegral, lsum_apply, sum_monomial_index]
  · rfl
  · simp

private noncomputable def angularWeight (M : ℕ) (j : ℕ) : K :=
  (M : K) / (↑(M + 2 * j) : K)

/-- Termwise angular integration followed by the Laurent constant term. -/
noncomputable def angularTrace (M : ℕ) : L[X] →ₗ[K] K :=
  laurentCT.comp (angularIntegral (angularWeight M))

private theorem angularIntegral_upper_endpoint (M : ℕ) (hM : 0 < M) (p : L[X]) :
    (M : K) • p.eval 1 =
      (M : K) • angularIntegral (angularWeight M) p +
      (2 : K) • angularIntegral (angularWeight M) (X * p.derivative) := by
  induction p using Polynomial.induction_on' with
  | add p q hp hq =>
    simp only [eval_add, smul_add, map_add, mul_add]
    rw [hp, hq]
    abel
  | monomial n c =>
    have hd : X * (monomial n c).derivative = monomial n ((n : K) • c) := by
      cases n with
      | zero => simp only [monomial_zero_left, derivative_C, mul_zero, Nat.cast_zero, zero_smul, map_zero]
      | succ n =>
        rw [derivative_monomial_succ, ← monomial_one_one_eq_X, monomial_mul_monomial,
          Algebra.smul_def, map_natCast, Nat.cast_succ]
        congr 1 <;> ring
    rw [hd, eval_monomial, one_pow, mul_one, angularIntegral_monomial,
      angularIntegral_monomial, smul_comm (angularWeight M n), ← smul_assoc, ← smul_assoc,
      ← smul_assoc, ← add_smul]
    congr 1
    have hne : (↑(M + 2 * n) : K) ≠ 0 := by exact_mod_cast (by omega : M + 2 * n ≠ 0)
    simp only [angularWeight, smul_eq_mul]
    field_simp
    push_cast
    ring

omit [CharZero K] in
private theorem eventual_iterate_derivative_moments (M : ℕ) (f : L[X])
    (safe : Safe (angularTrace M) f) (k : ℕ) (g : L[X]) :
    ∃ cutoff : ℕ, ∀ m ≥ cutoff, angularTrace M (derivative^[k] (g * f ^ m)) = 0 := by
  induction k generalizing g with
  | zero => exact safe g
  | succ k ih =>
    obtain ⟨c, hc⟩ := ih g.derivative
    obtain ⟨d, hd⟩ := ih (g * f.derivative)
    refine ⟨max c (d + 1), fun m hm ↦ ?_⟩
    have hc' : c ≤ m := (le_max_left _ _).trans hm
    have hd' : d ≤ m - 1 := by have := (le_max_right c (d + 1)).trans hm; omega
    have hder :
        derivative (g * f ^ m) = g.derivative * f ^ m + (m : K) • ((g * f.derivative) * f ^ (m - 1)) := by
      simp only [derivative_mul, derivative_pow, Algebra.smul_def, map_natCast]
      ring
    rw [Function.iterate_succ_apply, hder, ← Module.End.pow_apply, map_add, map_add]
    simp only [Module.End.pow_apply, iterate_derivative_smul, map_smul, hc _ hc', hd _ hd',
      smul_zero, add_zero]

/-- Eventual marked vanishing forces the upper angular endpoint to be
zero or to have Laurent exponents of a single strict sign. -/
theorem angular_upper_endpoint_oneSided (M : ℕ) (hM : 0 < M) (f : L[X])
    (safe : Safe (angularTrace M) f) : LaurentOneSided (f.eval 1) := by
  apply laurentOneSided_of_eventual_ct_zero
  obtain ⟨a, ha⟩ := safe 1
  obtain ⟨b, hb⟩ := safe (X * f.derivative)
  refine ⟨max a (b + 1), fun m hm ↦ ?_⟩
  have hma : a ≤ m := (le_max_left _ _).trans hm
  have hmb : b ≤ m - 1 := by have := (le_max_right a (b + 1)).trans hm; omega
  have h := congrArg laurentCT (angularIntegral_upper_endpoint M hM (f ^ m))
  have hder :
      X * derivative (f ^ m) = (m : K) • ((X * f.derivative) * f ^ (m - 1)) := by
    simp only [derivative_pow, Algebra.smul_def, map_natCast]
    ring
  rw [hder] at h
  simp only [eval_pow, map_smul, map_add, smul_eq_mul] at h
  change (M : K) * laurentCT ((f.eval 1) ^ m) =
    (M : K) * angularTrace M (f ^ m) +
    2 * ((m : K) * angularTrace M ((X * f.derivative) * f ^ (m - 1))) at h
  have ha' := ha m hma
  rw [one_mul] at ha'
  rw [ha', hb _ hmb, mul_zero, mul_zero, mul_zero, add_zero] at h
  exact (mul_eq_zero.mp h).resolve_left (by exact_mod_cast hM.ne')

private noncomputable def unitIntegral : L[X] →ₗ[K] L :=
  angularIntegral fun n ↦ (↑(n + 1) : K)⁻¹

private theorem unitIntegral_derivative (p : L[X]) :
    unitIntegral p.derivative = p.eval 1 - p.eval 0 := by
  induction p using Polynomial.induction_on' with
  | add p q hp hq =>
    simp only [map_add, eval_add, hp, hq]
    abel
  | monomial n c =>
    cases n with
    | zero => simp [monomial_zero_left]
    | succ n =>
      rw [derivative_monomial_succ, unitIntegral, angularIntegral_monomial]
      have hc : c * ((n : L) + 1) = (↑(n + 1) : K) • c := by
        rw [Algebra.smul_def, map_natCast, Nat.cast_succ, mul_comm c]
      rw [hc, ← smul_assoc, smul_eq_mul, inv_mul_cancel₀ (by exact_mod_cast Nat.succ_ne_zero n)]
      simp [eval_monomial]

private theorem angularIntegral_eq_unitIntegral (r : ℕ) (p : L[X]) :
    angularIntegral (angularWeight (2 * (r + 1))) p =
      (↑(r + 1) : K) • unitIntegral (X ^ r * p) := by
  induction p using Polynomial.induction_on' with
  | add p q hp hq => rw [map_add, mul_add (X ^ r) p q, map_add, smul_add, hp, hq]
  | monomial n c =>
    rw [X_pow_mul_monomial, unitIntegral, angularIntegral_monomial, angularIntegral_monomial,
      ← smul_assoc, smul_eq_mul, angularWeight]
    congr 1
    have hz : (↑(n + r + 1) : K) ≠ 0 := by exact_mod_cast Nat.succ_ne_zero (n + r)
    have he : 2 * (r + 1) + 2 * n = 2 * (n + r + 1) := by omega
    rw [he, Nat.cast_mul, div_eq_mul_inv, Nat.cast_mul, mul_inv_rev]
    field_simp

private theorem unitIntegral_iterate_derivative (r : ℕ) (p : L[X])
    (upper : ∀ j ≤ r, (derivative^[j] p).eval 1 = 0) :
    unitIntegral (X ^ r * derivative^[r + 1] p) =
      ((-1 : K) ^ (r + 1) * (r.factorial : K)) • p.eval 0 := by
  induction r with
  | zero =>
    have h := upper 0 le_rfl
    change p.eval 1 = 0 at h
    simp [unitIntegral_derivative, h]
  | succ r ih =>
    have ih' := ih fun j hj ↦ upper j (hj.trans (Nat.le_succ _))
    have hp := unitIntegral_derivative (X ^ (r + 1) * derivative^[r + 1] p)
    have hu := upper (r + 1) le_rfl
    have hd : derivative (X ^ (r + 1) * derivative^[r + 1] p) =
        (↑(r + 1) : K) • (X ^ r * derivative^[r + 1] p) + X ^ (r + 1) * derivative^[r + 1 + 1] p := by
      rw [derivative_mul, derivative_X_pow_succ, mul_assoc,
        Function.iterate_succ_apply' (⇑(derivative : L[X] →ₗ[L] L[X])) (r + 1) p,
        Algebra.smul_def, map_natCast, ← C_eq_natCast, Nat.cast_succ]
    rw [eval_mul, eval_mul, eval_pow, eval_pow, eval_X, eval_X, one_pow,
      hu, mul_zero, zero_pow (Nat.succ_ne_zero r), zero_mul, sub_zero, hd,
      map_add, map_smul, ih', ← smul_assoc, smul_eq_mul] at hp
    rw [eq_neg_of_add_eq_zero_right hp, ← neg_smul]
    congr 1
    rw [Nat.factorial_succ, Nat.cast_mul, pow_succ]
    ring

private theorem angularIntegral_lower_endpoint (r : ℕ) (p : L[X]) :
    angularIntegral (angularWeight (2 * (r + 1))) (derivative^[r + 1] ((1 - X) ^ (r + 1) * p)) =
      ((-1 : K) ^ (r + 1) * ((r + 1).factorial : K)) • p.eval 0 := by
  rw [angularIntegral_eq_unitIntegral, unitIntegral_iterate_derivative]
  · simp only [eval_mul, eval_pow, eval_sub, eval_one, eval_X, sub_zero,
      one_pow, one_mul, ← smul_assoc, smul_eq_mul]
    congr 1
    rw [Nat.factorial_succ, Nat.cast_mul]
    ring
  · intro j hj
    obtain ⟨g, hg⟩ := pow_sub_dvd_iterate_derivative_of_pow_dvd j
      (dvd_mul_right ((1 - X) ^ (r + 1)) p)
    rw [hg, eval_mul, eval_pow, eval_sub, eval_one, eval_X, sub_self,
      zero_pow (by omega : r + 1 - j ≠ 0), zero_mul]

/-- In even dimension, eventual marked vanishing also forces the lower
angular endpoint to be zero or strictly one-sided. -/
theorem angular_lower_endpoint_oneSided (M : ℕ) (hM : 0 < M) (even : Even M)
    (f : L[X]) (safe : Safe (angularTrace M) f) : LaurentOneSided (f.eval 0) := by
  obtain ⟨k, hk⟩ := even
  obtain ⟨r, hr⟩ := Nat.exists_eq_succ_of_ne_zero (by omega : k ≠ 0)
  have hMeq : M = 2 * (r + 1) := by omega
  rw [hMeq] at safe
  apply laurentOneSided_of_eventual_ct_zero
  obtain ⟨cutoff, hc⟩ := eventual_iterate_derivative_moments (2 * (r + 1)) f safe (r + 1) ((1 - X) ^ (r + 1))
  refine ⟨cutoff, fun m hm ↦ ?_⟩
  have h := congrArg laurentCT (angularIntegral_lower_endpoint r (f ^ m))
  change angularTrace (2 * (r + 1)) (derivative^[r + 1] ((1 - X) ^ (r + 1) * f ^ m)) = _ at h
  rw [hc m hm, map_smul, eval_pow, smul_eq_mul] at h
  have hcoeff : (-1 : K) ^ (r + 1) * ((r + 1).factorial : K) ≠ 0 :=
    mul_ne_zero (pow_ne_zero _ (neg_ne_zero.mpr one_ne_zero))
      (by exact_mod_cast Nat.factorial_ne_zero (r + 1))
  exact (mul_eq_zero.mp h.symm).resolve_left hcoeff

section Realization

variable {A : Type*} [CommRing A] [Algebra K A]

/-- Evaluate the block-invariant ring by sending the Laurent variable
to the norm unit and the angular variable to the height ratio. -/
noncomputable def angularPhase (U : Aˣ) (S : A) : L[X] →ₐ[K] A :=
  let χ : Multiplicative ℤ →* A :=
    { toFun := fun i ↦ (↑(U ^ i.toAdd) : A)
      map_one' := by simp
      map_mul' := fun i j ↦ by
        change (↑(U ^ (i.toAdd + j.toAdd)) : A) = _
        rw [zpow_add, Units.val_mul] }
  Polynomial.eval₂AlgHom (AddMonoidAlgebra.lift K A ℤ χ) S (fun _ ↦ Commute.all _ _)

omit [CharZero K] in
theorem angularTrace_eq_realization (M : ℕ) (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
    (h : NormalizerTrace M trace U S) (p : L[X]) :
    trace (angularPhase U S p) = angularTrace M p := by
  induction p using Polynomial.induction_on' with
  | add p q hp hq => simp [hp, hq]
  | monomial n c =>
    induction c using AddMonoidAlgebra.induction_linear with
    | zero => simp
    | add c d hc hd => simp only [map_add, hc, hd]
    | single i c =>
      simp only [angularTrace, LinearMap.comp_apply, angularIntegral_monomial, map_smul]
      change trace (Polynomial.eval₂ _ S (monomial n (AddMonoidAlgebra.single i c))) = _
      rw [Polynomial.eval₂_monomial]
      change trace ((((AddMonoidAlgebra.lift K A ℤ) _) (AddMonoidAlgebra.single i c)) * S ^ n) = _
      rw [AddMonoidAlgebra.lift_single, smul_mul_assoc, map_smul]
      change c • trace ((↑(U ^ i) : A) * S ^ n) = angularWeight M n •
        laurentCT (AddMonoidAlgebra.single i c)
      by_cases hi : i = 0
      · subst i
        rw [zpow_zero, Units.val_one, one_mul, h.height,
          show laurentCT (AddMonoidAlgebra.single (0 : ℤ) c) = c by
            change (Finsupp.single (0 : ℤ) c) 0 = c
            exact Finsupp.single_eq_same]
        exact mul_comm _ _
      · rw [h.weight i n hi]
        change c • (0 : K) = angularWeight M n • (Finsupp.single i c) 0
        rw [Finsupp.single_eq_of_ne' hi, smul_zero, smul_zero]

omit [CharZero K] in
private theorem safe_angularTrace_of_safe_realization (M : ℕ) (trace : A →ₗ[K] K)
    (U : Aˣ) (S : A) (h : NormalizerTrace M trace U S) (f : L[X])
    (safe : Safe trace (angularPhase U S f)) : Safe (angularTrace M) f := by
  intro g
  obtain ⟨n, hn⟩ := safe (angularPhase U S g)
  refine ⟨n, fun m hm ↦ ?_⟩
  rw [← angularTrace_eq_realization M trace U S h, map_mul, map_pow]
  exact hn m hm

/-- The upper-endpoint proposition for the norm and height chart. -/
theorem draft_angular_upper_endpoint (M : ℕ) (hM : 0 < M) (trace : A →ₗ[K] K)
    (U : Aˣ) (S : A) (h : NormalizerTrace M trace U S) (f : L[X])
    (safe : Safe trace (angularPhase U S f)) : LaurentOneSided (f.eval 1) :=
  angular_upper_endpoint_oneSided M hM f (safe_angularTrace_of_safe_realization M trace U S h f safe)

/-- The draft's angular endpoint proposition on the block chart. -/
theorem draft_angular_endpoints (M : ℕ) (hM : 0 < M) (trace : A →ₗ[K] K)
    (U : Aˣ) (S : A) (h : NormalizerTrace M trace U S) (f : L[X])
    (safe : Safe trace (angularPhase U S f)) :
    LaurentOneSided (f.eval 1) ∧ (Even M → LaurentOneSided (f.eval 0)) :=
  ⟨draft_angular_upper_endpoint M hM trace U S h f safe,
    fun even ↦ angular_lower_endpoint_oneSided M hM even f
      (safe_angularTrace_of_safe_realization M trace U S h f safe)⟩

end Realization

end MathieuQuadratic
