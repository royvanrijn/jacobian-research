import MathieuQuadratic.Complement
import Mathlib.Algebra.DualNumber
import Mathlib.FieldTheory.IsAlgClosed.Basic
import Mathlib.RingTheory.TensorProduct.Maps

/-!
# Absolute unit-times-power exclusion

A simple zero gives a map to dual numbers with zero residue and nonzero
linear coefficient. Such an image cannot be a unit times any proper power.
The jets are defined at norm −1/2 for the uniform phase and norm −1 for
the direct ternary phase.
-/

open scoped DualNumber TensorProduct BigOperators
open TrivSqZeroExt

namespace MathieuQuadratic

variable {K L : Type*} [Field K] [CharZero K] [Field L] [Algebra K L] [CharZero L]
variable {M : ℕ}

omit [CharZero L] in
theorem nonpower_of_dualNumber {D : Type*} [CommRing D] (f : D →+* L[ε]) (v : D)
    (h0 : (f v).fst = 0) (h1 : (f v).snd ≠ 0) : NonpowerUpToUnits v := by
  intro u g e he h
  have hprod := congrArg (fun x ↦ (f x).fst) h
  simp only [h0, map_mul, map_pow, fst_mul, fst_pow] at hprod
  have hu : (f (u : D)).fst ≠ 0 :=
    (((u.isUnit.map f).map (TrivSqZeroExt.fstHom L L L)).ne_zero)
  have hg : (f g).fst = 0 := by
    by_contra hn
    exact (mul_ne_zero hu (pow_ne_zero e hn)) hprod.symm
  have hz : f g ^ 2 = 0 := by
    apply TrivSqZeroExt.ext
    · simp [hg]
    · simp [pow_two, hg]
  have heq := congrArg f h
  rw [map_mul, map_pow, show e = 2 + (e - 2) by omega, pow_add, hz, zero_mul, mul_zero] at heq
  exact h1 (by rw [heq, snd_zero])

private def jetCoordinates (r : L) :
    Fin 2 ⊕ Fin M → L[ε] :=
  Sum.elim (fun i ↦ if i = 0 then inl r + inr 1 else 0) (fun _ ↦ 0)

private noncomputable def jetPolynomial (r : L) :
    QuadPolynomial K M →ₐ[K] L[ε] :=
  MvPolynomial.aeval (jetCoordinates r)

omit [CharZero K] [CharZero L] in
private theorem jetPolynomial_height (r : L) (b : Fin M → K) :
    jetPolynomial (K := K) r (heightNorm b) = 0 := by
  simp [jetPolynomial, heightNorm, jetCoordinates]

omit [CharZero K] [CharZero L] in
private theorem jetPolynomial_norm (a : Fin 2 → K) (b : Fin M → K) (r : L) :
    jetPolynomial r (quadNorm a b) =
      inl (algebraMap K L (a 0)) * (inl r + inr 1) ^ 2 := by
  rw [quadNorm, map_add, jetPolynomial_height, add_zero, planeNorm, map_sum,
    Finset.sum_eq_single (0 : Fin 2)]
  · simp [jetPolynomial, jetCoordinates, IsScalarTower.algebraMap_apply K L L[ε]]
    rfl
  · intro i _ hi
    simp [jetPolynomial, jetCoordinates, hi]
  · simp

omit [CharZero K] in
private theorem exists_complement_jet [IsAlgClosed L] (a : Fin 2 → K) (b : Fin M → K)
    (ha : a 0 ≠ 0) (c : L) (hc : c ≠ 0) :
    ∃ jet : NormComplement a b →ₐ[K] L[ε],
      (jet (normUnit a b : NormComplement a b)).fst = c ∧
      (jet (normUnit a b : NormComplement a b)).snd ≠ 0 ∧
      jet (height a b) = 0 := by
  have haL : algebraMap K L (a 0) ≠ 0 := by
    simpa using (RingHom.injective (algebraMap K L)).ne ha
  obtain ⟨r, hr⟩ := IsAlgClosed.exists_pow_nat_eq (c / algebraMap K L (a 0)) (n := 2) (by decide)
  have hr0 : r ≠ 0 := by
    intro h
    exact (div_ne_zero hc haL) (by simpa [h] using hr.symm)
  have hfst : (jetPolynomial r (quadNorm a b)).fst = c := by
    rw [jetPolynomial_norm, fst_mul, fst_inl, fst_pow, fst_add, fst_inl, fst_inr, add_zero,
      hr, mul_div_cancel₀ _ haL]
  have hsnd : (jetPolynomial r (quadNorm a b)).snd ≠ 0 := by
    rw [jetPolynomial_norm, DualNumber.snd_mul, fst_inl, snd_inl, zero_mul, add_zero,
      pow_two, DualNumber.snd_mul]
    simp only [fst_add, fst_inl, fst_inr, add_zero, snd_add, snd_inl, snd_inr, zero_add,
      mul_one, one_mul]
    exact mul_ne_zero haL (by rwa [← two_mul, mul_ne_zero_iff, and_iff_right two_ne_zero])
  have hunit : IsUnit (jetPolynomial r (quadNorm a b)) :=
    TrivSqZeroExt.isUnit_iff_isUnit_fst.mpr (hfst ▸ isUnit_iff_ne_zero.mpr hc)
  let jet : NormComplement a b →ₐ[K] L[ε] :=
    IsLocalization.Away.liftAlgHom (quadNorm a b) hunit
  have hmap (p : QuadPolynomial K M) :
      jet (algebraMap (QuadPolynomial K M) (NormComplement a b) p) = jetPolynomial r p :=
    IsLocalization.Away.lift_eq _ hunit _
  exact ⟨jet, by rwa [normUnit_val, hmap], by rwa [normUnit_val, hmap],
    by rw [height, map_mul, hmap, jetPolynomial_height, zero_mul]⟩

omit [CharZero K] [CharZero L] in
private theorem normalizer_jet {A : Type*} [CommRing A] [Algebra K A]
    (jet : A →ₐ[K] L[ε]) (U : Aˣ) (S : A) (c : L) (hc : c ≠ 0)
    (hfst : (jet (U : A)).fst = c) (hsnd : (jet (U : A)).snd ≠ 0) (hS : jet S = 0) :
    (jet (normalizerPhase U S)).fst = 1 + c⁻¹ ∧ (jet (normalizerPhase U S)).snd ≠ 0 := by
  have hu : jet (↑(U⁻¹) : A) * jet (U : A) = 1 := by
    rw [← map_mul, Units.inv_mul, map_one]
  have hf := congrArg (fun v : L[ε] ↦ v.fst) hu
  simp only [fst_mul, fst_one, hfst] at hf
  have hf' : (jet (↑(U⁻¹) : A)).fst = c⁻¹ := (eq_div_iff hc).2 hf |>.trans (one_div c)
  have hn : (jet (↑(U⁻¹) : A)).snd ≠ 0 := by
    intro h
    have he := congrArg (fun v : L[ε] ↦ v.snd) hu
    simp only [DualNumber.snd_mul, snd_one, h, zero_mul, add_zero, hf'] at he
    exact (mul_ne_zero (inv_ne_zero hc) hsnd) he
  have value : jet (normalizerPhase U S) = 1 + jet (↑(U⁻¹) : A) := by
    simp only [normalizerPhase, map_mul, map_add, map_one, map_pow, map_sub, hS, mul_zero,
      sub_zero, one_pow, mul_one, mul_add, hu]
    rw [add_comm]
  rw [value]
  exact ⟨by simp [hf'], by simpa using hn⟩

private noncomputable def extendJet {A : Type*} [CommRing A] [Algebra K A]
    (jet : A →ₐ[K] L[ε]) : L ⊗[K] A →ₐ[L] L[ε] :=
  Algebra.TensorProduct.lift (Algebra.ofId L L[ε]) jet (fun _ _ ↦ Commute.all _ _)

omit [CharZero K] [CharZero L] in
private theorem extendJet_includeRight {A : Type*} [CommRing A] [Algebra K A]
    (jet : A →ₐ[K] L[ε]) (x : A) :
    extendJet jet (Algebra.TensorProduct.includeRight (R := K) (A := L) x) = jet x := by
  rw [Algebra.TensorProduct.includeRight_apply, extendJet, Algebra.TensorProduct.lift_tmul,
    map_one, one_mul]

omit [CharZero K] in
/-- Equation (16) gives the absolute non-power condition over any
algebraically closed extension. -/
theorem uniformPhase_absolute_nonpower [IsAlgClosed L] (a : Fin 2 → K) (b : Fin M → K)
    (ha : a 0 ≠ 0) :
    NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := L)
      (uniformPhase M (normalizerPhase (normUnit a b) (height a b)))) := by
  obtain ⟨jet, hfst, hsnd, hS⟩ :=
    exists_complement_jet a b ha (-1 / 2 : L) (by norm_num)
  obtain ⟨hp, hpn⟩ := normalizer_jet jet (normUnit a b) (height a b) _ (by norm_num) hfst hsnd hS
  have hpf : (jet (normalizerPhase (normUnit a b) (height a b))).fst = -1 := by
    rw [hp]
    norm_num
  apply nonpower_of_dualNumber (extendJet jet).toRingHom
  · change (extendJet jet (Algebra.TensorProduct.includeRight _)).fst = _
    rw [extendJet_includeRight]
    simp [uniformPhase, hpf]
  · change (extendJet jet (Algebra.TensorProduct.includeRight _)).snd ≠ _
    rw [extendJet_includeRight]
    simp only [uniformPhase, map_mul, map_pow, map_add, map_one, DualNumber.snd_mul, fst_pow,
      hpf, fst_add, fst_one, snd_add, snd_one, zero_add, add_neg_cancel, mul_zero, add_zero]
    exact mul_ne_zero (pow_ne_zero _ (neg_ne_zero.mpr one_ne_zero)) hpn

omit [CharZero K] in
theorem directPhase_absolute_nonpower [IsAlgClosed L] (a : Fin 2 → K) (b : Fin M → K)
    (ha : a 0 ≠ 0) :
    NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := L)
      (normalizerPhase (normUnit a b) (height a b))) := by
  obtain ⟨jet, hfst, hsnd, hS⟩ :=
    exists_complement_jet a b ha (-1 : L) (by norm_num)
  obtain ⟨hp, hpn⟩ := normalizer_jet jet (normUnit a b) (height a b) _ (by norm_num) hfst hsnd hS
  apply nonpower_of_dualNumber (extendJet jet).toRingHom
  · change (extendJet jet (Algebra.TensorProduct.includeRight _)).fst = _
    rw [extendJet_includeRight, hp]
    norm_num
  · change (extendJet jet (Algebra.TensorProduct.includeRight _)).snd ≠ _
    rw [extendJet_includeRight]
    exact hpn

end MathieuQuadratic
