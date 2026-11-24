import MathieuQuadratic.Complement
import Mathlib.Algebra.DualNumber
import Mathlib.FieldTheory.IsAlgClosed.Basic
import Mathlib.RingTheory.TensorProduct.Maps
import Mathlib.RingTheory.PowerSeries.Order
import Mathlib.RingTheory.PowerSeries.Inverse

/-!
# Absolute unit-times-power exclusion

A simple zero gives a map to dual numbers with zero residue and nonzero
linear coefficient. Such an image cannot be a unit times any proper power.
For the product phase, the two first-order restrictions give orders β
and γ, whose greatest common divisor controls the power assertion.
The perturbation uses a simple zero on the norm-one level.
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

private noncomputable def extendJet {A D : Type*} [CommRing A] [Algebra K A]
    [CommRing D] [Algebra K D] [Algebra L D] [IsScalarTower K L D]
    (jet : A →ₐ[K] D) : L ⊗[K] A →ₐ[L] D :=
  Algebra.TensorProduct.lift (Algebra.ofId L D) jet (fun _ _ ↦ Commute.all _ _)

omit [CharZero K] [CharZero L] in
private theorem extendJet_includeRight {A D : Type*} [CommRing A] [Algebra K A]
    [CommRing D] [Algebra K D] [Algebra L D] [IsScalarTower K L D]
    (jet : A →ₐ[K] D) (x : A) :
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

open scoped PowerSeries ENat
open PowerSeries

/-- The constant and linear coefficients of an arc, viewed as dual numbers. -/
private noncomputable def seriesJet : L⟦X⟧ →ₐ[L] L[ε] where
  toFun f := ⟨constantCoeff f, coeff 1 f⟩
  map_one' := by ext <;> simp [coeff_one]
  map_zero' := by ext <;> simp
  map_add' f g := by ext <;> simp
  map_mul' f g := by
    ext
    · exact map_mul constantCoeff f g
    · change coeff 1 (f * g) = _
      rw [coeff_one_mul, DualNumber.snd_mul]
      change coeff 1 f * constantCoeff g + coeff 1 g * constantCoeff f =
        constantCoeff f * coeff 1 g + coeff 1 f * constantCoeff g
      ring
  commutes' r := by
    change (⟨constantCoeff (C r), coeff 1 (C r)⟩ : L[ε]) = inl r
    apply TrivSqZeroExt.ext <;> simp

omit [CharZero L] in
private theorem seriesJet_lift (z : L[ε]) :
    seriesJet (C z.fst + X * C z.snd) = z := by
  apply TrivSqZeroExt.ext
  · change constantCoeff (C z.fst + X * C z.snd) = z.fst
    simp
  · change coeff 1 (C z.fst + X * C z.snd) = z.snd
    simp [coeff_mul_C]

omit [CharZero L] in
private theorem order_one_of_jet {f : L⟦X⟧}
    (hf : (seriesJet f).fst = 0) (hg : (seriesJet f).snd ≠ 0) : f.order = 1 := by
  apply (order_eq_nat).2
  refine ⟨hg, ?_⟩
  intro n hn
  change constantCoeff f = 0 at hf
  simpa [show n = 0 by omega, coeff_zero_eq_constantCoeff] using hf

omit [CharZero L] in
private theorem order_zero_of_jet {f : L⟦X⟧} (hf : (seriesJet f).fst ≠ 0) :
    f.order = 0 :=
  order_zero_of_unit (isUnit_iff_constantCoeff.2 (isUnit_iff_ne_zero.2 hf))

omit [CharZero K] [CharZero L] in
private theorem restriction_dvd_order {A : Type*} [CommRing A]
    (f : A →+* L⟦X⟧) {u : Aˣ} {g v : A} {e n : ℕ}
    (hv : v = u * g ^ e) (hn : (f v).order = n) : e ∣ n := by
  rw [hv, map_mul, map_pow, order_mul, order_zero_of_unit (u.isUnit.map f), zero_add, order_pow] at hn
  refine ⟨(f g).order.toNat, ?_⟩
  apply congrArg ENat.toNat at hn
  simpa [nsmul_eq_mul] using hn.symm

omit [CharZero K] [CharZero L] in
/-- An arc through the given dual-number point, defined by lifting each
coordinate to its affine-linear series. -/
private theorem exists_series_arc (a : Fin 2 → K) (b : Fin M → K)
    (z : Fin 2 ⊕ Fin M → L[ε])
    (hn : ((MvPolynomial.aeval z) (quadNorm a b)).fst ≠ 0) :
    ∃ f : NormComplement a b →ₐ[K] L⟦X⟧,
      seriesJet (f (normUnit a b : NormComplement a b)) = MvPolynomial.aeval z (quadNorm a b) ∧
      seriesJet (f (algebraMap (QuadPolynomial K M) (NormComplement a b) (heightNorm b))) =
        MvPolynomial.aeval z (heightNorm b) := by
  let eval : QuadPolynomial K M →ₐ[K] L⟦X⟧ :=
    MvPolynomial.aeval (fun i ↦ C (z i).fst + X * C (z i).snd)
  have he (p : QuadPolynomial K M) : seriesJet (eval p) = MvPolynomial.aeval z p := by
    change ((seriesJet (L := L)).restrictScalars K) (eval p) = _
    simp [eval, MvPolynomial.comp_aeval_apply, seriesJet_lift]
  have hu : IsUnit (eval (quadNorm a b)) :=
    isUnit_iff_constantCoeff.2 (isUnit_iff_ne_zero.2 (by
      change (seriesJet (eval (quadNorm a b))).fst ≠ 0
      rwa [he]))
  let f : NormComplement a b →ₐ[K] L⟦X⟧ := IsLocalization.Away.liftAlgHom (quadNorm a b) hu
  have hv (p : QuadPolynomial K M) :
      f (algebraMap (QuadPolynomial K M) (NormComplement a b) p) = eval p :=
    IsLocalization.Away.lift_eq _ hu _
  exact ⟨f, by rw [normUnit_val, hv, he], by rw [hv, he]⟩

omit [CharZero K] in
/-- Given the jets of q and B, their elementary relations give the two
orders used to rule out any common exponent. -/
private theorem product_order_at_jet (a : Fin 2 → K) (b : Fin M → K)
    (z : Fin 2 ⊕ Fin M → L[ε]) (α β γ : ℕ)
    (mode : Bool)
    (hn : ((MvPolynomial.aeval z) (quadNorm a b)).fst ≠ 0)
    (hs : if mode then ((MvPolynomial.aeval z) (quadNorm a b)).fst = -1 ∧
        ((MvPolynomial.aeval z) (quadNorm a b)).snd ≠ 0
      else ((MvPolynomial.aeval z) (quadNorm a b)).fst = 1 ∧
        ((MvPolynomial.aeval z) (heightNorm b)).fst = 1 / 4 ∧
        ((MvPolynomial.aeval z) (heightNorm b)).snd =
          ((MvPolynomial.aeval z) (quadNorm a b)).snd ∧
        ((MvPolynomial.aeval z) (heightNorm b)).snd ≠ 0) :
    ∃ f : NormComplement a b →ₐ[K] L⟦X⟧,
      (f (productPhase (normUnit a b) (height a b) α β γ)).order =
        ((if mode then β else γ) : ℕ) := by
  obtain ⟨f, hf, hb⟩ := exists_series_arc a b z hn
  let U := normUnit a b
  let W := 1 - (1 + (U : NormComplement a b)) ^ 2 * height a b
  let V := 1 + (U : NormComplement a b)
  let H := algebraMap (QuadPolynomial K M) (NormComplement a b) (heightNorm b)
  let q := MvPolynomial.aeval z (quadNorm a b)
  let B := MvPolynomial.aeval z (heightNorm b)
  let jet := (seriesJet (L := L)).restrictScalars K |>.comp f
  change jet (U : NormComplement a b) = q at hf
  change jet H = B at hb
  have hWU : W * (U : NormComplement a b) = (U : NormComplement a b) - V ^ 2 * H := by
    change (1 - V ^ 2 * (H * (↑(U⁻¹) : NormComplement a b))) * (U : NormComplement a b) = _
    rw [sub_mul, one_mul, mul_assoc (V ^ 2), mul_assoc H, Units.inv_mul, mul_one]
  have hw := congrArg jet hWU
  have hV : jet V = 1 + q := by
    simp [V, hf]
  simp only [map_mul, map_sub, map_pow, hf, hV, hb] at hw
  have hi : (f (↑(U⁻¹) : NormComplement a b)).order = 0 :=
    order_zero_of_unit ((U⁻¹).isUnit.map f.toRingHom)
  refine ⟨f, ?_⟩
  change (f ((↑(U⁻¹) : NormComplement a b) ^ α * V ^ β * W ^ γ)).order = _
  rw [map_mul, map_mul, map_pow, map_pow, map_pow, order_mul, order_mul,
    order_pow, order_pow, order_pow, hi, nsmul_zero, zero_add]
  cases mode with
  | true =>
    obtain ⟨hq, hq'⟩ := hs
    change q.fst = -1 at hq
    change q.snd ≠ 0 at hq'
    have hv0 : (jet V).fst = 0 := by simp [hV, hq]
    have hv1 : (jet V).snd ≠ 0 := by simpa [hV] using hq'
    have hw0 : (jet W).fst ≠ 0 := by
      have he := congrArg (fun x : L[ε] ↦ x.fst) hw
      simp [hq] at he
      simp [he]
    have hv : (f V).order = 1 := order_one_of_jet hv0 hv1
    have hw' : (f W).order = 0 := order_zero_of_jet hw0
    rw [hv, hw']
    simp
  | false =>
    obtain ⟨hq, hB, hBq, hB'⟩ := hs
    change q.fst = 1 at hq
    change B.fst = 1 / 4 at hB
    change B.snd = q.snd at hBq
    change B.snd ≠ 0 at hB'
    have hv0 : (jet V).fst ≠ 0 := by norm_num [hV, hq]
    have hw0 : (jet W).fst = 0 := by
      have he := congrArg (fun x : L[ε] ↦ x.fst) hw
      norm_num [hq, hB] at he
      exact he
    have hw1 : (jet W).snd ≠ 0 := by
      have he := congrArg (fun x : L[ε] ↦ x.snd) hw
      simp only [pow_two, DualNumber.snd_mul, fst_mul, fst_add, fst_one, hq, snd_sub, hw0,
        snd_add, snd_one, zero_add, hB, ← hBq, mul_one] at he
      have hval : (jet W).snd = -4 * B.snd := by
        convert he using 1 <;> ring
      rw [hval]
      exact mul_ne_zero (by norm_num) hB'
    have hv : (f V).order = 0 := order_zero_of_jet hv0
    have hw' : (f W).order = 1 := order_one_of_jet hw0 hw1
    rw [hv, hw']
    simp

private def blockJetCoordinates (i : Fin M) (x z : L[ε]) :
    Fin 2 ⊕ Fin M → L[ε] :=
  Sum.elim (fun j ↦ if j = 0 then x else 0) (fun j ↦ if j = i then z else 0)

omit [CharZero K] [CharZero L] in
private theorem blockJet_height (b : Fin M → K) (i : Fin M) (x z : L[ε]) :
    (MvPolynomial.aeval (blockJetCoordinates i x z)) (heightNorm b) =
      inl (algebraMap K L (b i)) * z ^ 2 := by
  rw [heightNorm, map_sum, Finset.sum_eq_single i]
  · simp [blockJetCoordinates, IsScalarTower.algebraMap_apply K L L[ε]]
    rfl
  · intro j _ hj
    simp [blockJetCoordinates, hj]
  · simp

omit [CharZero K] [CharZero L] in
private theorem blockJet_norm (a : Fin 2 → K) (b : Fin M → K) (i : Fin M) (x z : L[ε]) :
    (MvPolynomial.aeval (blockJetCoordinates i x z)) (quadNorm a b) =
      inl (algebraMap K L (a 0)) * x ^ 2 + inl (algebraMap K L (b i)) * z ^ 2 := by
  rw [quadNorm, map_add, blockJet_height, planeNorm, map_sum, Finset.sum_eq_single (0 : Fin 2)]
  · simp [blockJetCoordinates, IsScalarTower.algebraMap_apply K L L[ε]]
    rfl
  · intro j _ hj
    simp [blockJetCoordinates, hj]
  · simp

omit [CharZero K] in
/-- Every power exponent divides the two orders β and γ of the draft. -/
theorem productPhase_absolute_nonpower [IsAlgClosed L] (a : Fin 2 → K) (b : Fin M → K)
    (ha : a 0 ≠ 0) (i : Fin M) (hb : b i ≠ 0)
    (α β γ : ℕ) (hcoprime : β.Coprime γ) :
    NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := L)
      (productPhase (normUnit a b) (height a b) α β γ)) := by
  let a₀ := algebraMap K L (a 0)
  let bᵢ := algebraMap K L (b i)
  have ha₀ : a₀ ≠ 0 := by
    change algebraMap K L (a 0) ≠ 0
    exact (map_ne_zero_iff (algebraMap K L) (RingHom.injective (algebraMap K L))).mpr ha
  have hbᵢ : bᵢ ≠ 0 := by
    change algebraMap K L (b i) ≠ 0
    exact (map_ne_zero_iff (algebraMap K L) (RingHom.injective (algebraMap K L))).mpr hb
  obtain ⟨r, hr⟩ := IsAlgClosed.exists_pow_nat_eq ((-1 : L) / a₀) (n := 2) (by decide)
  obtain ⟨s, hs⟩ := IsAlgClosed.exists_pow_nat_eq ((3 / 4 : L) / a₀) (n := 2) (by decide)
  obtain ⟨t, ht⟩ := IsAlgClosed.exists_pow_nat_eq ((1 / 4 : L) / bᵢ) (n := 2) (by decide)
  have hr₀ : r ≠ 0 := by
    intro h
    exact (div_ne_zero (neg_ne_zero.mpr one_ne_zero) ha₀) (by rw [← hr, h]; simp)
  have ht₀ : t ≠ 0 := by
    intro h
    exact (div_ne_zero (by norm_num : (1 / 4 : L) ≠ 0) hbᵢ) (by rw [← ht, h]; simp)
  have hfst₁ : (inl a₀ * (inl r + inr 1) ^ 2 : L[ε]).fst = -1 := by rw [fst_mul, fst_inl, fst_pow, fst_add, fst_inl, fst_inr, add_zero, hr, mul_div_cancel₀ _ ha₀]
  have hfst₂ : (inl a₀ * (inl s : L[ε]) ^ 2).fst = 3 / 4 := by rw [fst_mul, fst_inl, fst_pow, fst_inl, hs, mul_div_cancel₀ _ ha₀]
  have hfst₃ : (inl bᵢ * (inl t + inr 1) ^ 2 : L[ε]).fst = 1 / 4 := by rw [fst_mul, fst_inl, fst_pow, fst_add, fst_inl, fst_inr, add_zero, ht, mul_div_cancel₀ _ hbᵢ]
  have hlin (c v : L) :
      (inl c * (inl v + inr 1) ^ 2 : L[ε]).snd = c * (2 * v) := by
    simp only [pow_two, DualNumber.snd_mul, fst_inl, snd_inl, fst_add, fst_inr, add_zero,
      snd_add, snd_inr, zero_add, one_mul, mul_one]
    ring
  have haxis := product_order_at_jet a b (blockJetCoordinates i (inl r + inr 1) 0) α β γ true
  simp only [blockJet_norm, zero_pow (by decide : (2 : ℕ) ≠ 0), mul_zero, add_zero, if_true] at haxis
  change (inl a₀ * (inl r + inr 1) ^ 2 : L[ε]).fst ≠ 0 →
    ((inl a₀ * (inl r + inr 1) ^ 2 : L[ε]).fst = -1 ∧
      (inl a₀ * (inl r + inr 1) ^ 2 : L[ε]).snd ≠ 0) → _ at haxis
  obtain ⟨f, hf⟩ := haxis (by rw [hfst₁]; norm_num)
    ⟨hfst₁, by rw [hlin]; exact mul_ne_zero ha₀ (mul_ne_zero two_ne_zero hr₀)⟩
  have hwall := product_order_at_jet a b (blockJetCoordinates i (inl s) (inl t + inr 1)) α β γ false
  simp only [blockJet_norm, blockJet_height, Bool.false_eq_true, if_false] at hwall
  have hq₁ : (inl a₀ * (inl s : L[ε]) ^ 2 + inl bᵢ * (inl t + inr 1) ^ 2).fst = 1 := by
    rw [fst_add, hfst₂, hfst₃]
    norm_num
  have hq₂ : (inl a₀ * (inl s : L[ε]) ^ 2 + inl bᵢ * (inl t + inr 1) ^ 2).snd =
      (inl bᵢ * (inl t + inr 1) ^ 2 : L[ε]).snd := by
    simp only [snd_add, pow_two, DualNumber.snd_mul, fst_inl, snd_inl, mul_zero,
      zero_mul, add_zero, zero_add]
  obtain ⟨g, hg⟩ := hwall (by change (inl a₀ * (inl s : L[ε]) ^ 2 + inl bᵢ * (inl t + inr 1) ^ 2).fst ≠ 0; rw [hq₁]; exact one_ne_zero)
    ⟨hq₁, hfst₃, hq₂.symm, by rw [hlin]; exact mul_ne_zero hbᵢ (mul_ne_zero two_ne_zero ht₀)⟩
  intro u p e he h
  have hdβ : e ∣ β := restriction_dvd_order (extendJet (L := L) f).toRingHom h (by
    change ((extendJet (L := L) f) (Algebra.TensorProduct.includeRight _)).order = _
    rw [extendJet_includeRight]
    exact hf)
  have hdγ : e ∣ γ := restriction_dvd_order (extendJet (L := L) g).toRingHom h (by
    change ((extendJet (L := L) g) (Algebra.TensorProduct.includeRight _)).order = _
    rw [extendJet_includeRight]
    exact hg)
  have hdiv := Nat.dvd_gcd hdβ hdγ
  rw [show Nat.gcd β γ = 1 from hcoprime] at hdiv
  have he₁ := Nat.dvd_one.mp hdiv
  omega

omit [CharZero L] in
private theorem dual_scaled_square (c r s : L) :
    inl c * (inl r + inr s : L[ε]) ^ 2 = inl (c * r ^ 2) + inr (2 * c * r * s) := by
  apply TrivSqZeroExt.ext
  · simp
  · simp only [pow_two, DualNumber.snd_mul, fst_inl, snd_inl, fst_add, fst_inr, add_zero,
      snd_add, snd_inr, zero_add]
    ring

omit [CharZero K] in
/-- The added factor has a simple zero on the norm-one level, where
the base product phase is nonzero. -/
theorem perturbedProductPhase_absolute_nonpower [IsAlgClosed L]
    (a : Fin 2 → K) (b : Fin M → K) (ha : a 0 ≠ 0) (i : Fin M) (hb : b i ≠ 0)
    (hM : 4 ≤ M) :
    NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := L)
      (perturbedProductPhase (normUnit a b) (height a b) M)) := by
  let a₀ := algebraMap K L (a 0)
  let bᵢ := algebraMap K L (b i)
  have ha₀ : a₀ ≠ 0 := by
    change algebraMap K L (a 0) ≠ 0
    exact (map_ne_zero_iff (algebraMap K L) (RingHom.injective (algebraMap K L))).mpr ha
  have hbᵢ : bᵢ ≠ 0 := by
    change algebraMap K L (b i) ≠ 0
    exact (map_ne_zero_iff (algebraMap K L) (RingHom.injective (algebraMap K L))).mpr hb
  let c : L := (M - 1 : ℕ) + M
  let d : L := -1 / c
  have hc : c ≠ 0 := by
    dsimp [c]
    norm_cast
    omega
  have hc₁ : c + 1 ≠ 0 := by
    dsimp [c]
    norm_cast
  have hc₄ : c + 4 ≠ 0 := by
    dsimp [c]
    norm_cast
  have hdc : d * c = -1 := div_mul_cancel₀ _ hc
  have hd : d ≠ 0 := div_ne_zero (by norm_num) hc
  have hd₁ : 1 - d ≠ 0 := by
    intro hz
    have he := congrArg (fun k : L ↦ k * c) hz
    simp only [sub_mul, one_mul, hdc, zero_mul, sub_neg_eq_add] at he
    exact hc₁ he
  have hw : 1 - 4 * d ≠ 0 := by
    intro hz
    have he := congrArg (fun k : L ↦ k * c) hz
    rw [sub_mul, one_mul, mul_assoc, hdc, zero_mul] at he
    apply hc₄
    linear_combination he
  obtain ⟨r, hr⟩ := IsAlgClosed.exists_pow_nat_eq ((1 - d) / a₀) (n := 2) (by decide)
  obtain ⟨s, hs⟩ := IsAlgClosed.exists_pow_nat_eq (d / bᵢ) (n := 2) (by decide)
  have hr₀ : r ≠ 0 := by
    intro h
    exact (div_ne_zero hd₁ ha₀) (by rw [← hr, h]; simp)
  have hs₀ : s ≠ 0 := by
    intro h
    exact (div_ne_zero hd hbᵢ) (by rw [← hs, h]; simp)
  let x : L[ε] := inl r + inr ((-1 : L) / (2 * a₀ * r))
  let z : L[ε] := inl s + inr ((1 : L) / (2 * bᵢ * s))
  let v := blockJetCoordinates i x z
  have hx : inl a₀ * x ^ 2 = inl (1 - d) + inr (-1) := by
    dsimp only [x]
    rw [dual_scaled_square, hr, mul_div_cancel₀ _ ha₀, mul_div_cancel₀ _ (mul_ne_zero (mul_ne_zero two_ne_zero ha₀) hr₀)]
  have hz : inl bᵢ * z ^ 2 = inl d + inr 1 := by
    dsimp only [z]
    rw [dual_scaled_square, hs, mul_div_cancel₀ _ hbᵢ, mul_div_cancel₀ _ (mul_ne_zero (mul_ne_zero two_ne_zero hbᵢ) hs₀)]
  have hq : MvPolynomial.aeval v (quadNorm a b) = inl 1 := by
    dsimp only [v]
    rw [blockJet_norm, show inl (algebraMap K L (a 0)) * x ^ 2 = _ from hx,
      show inl (algebraMap K L (b i)) * z ^ 2 = _ from hz]
    apply TrivSqZeroExt.ext <;> simp
  have hB : MvPolynomial.aeval v (heightNorm b) = inl d + inr 1 := by
    dsimp only [v]
    rw [blockJet_height]
    exact hz
  obtain ⟨f, hfu, hfB⟩ := exists_series_arc a b v (by rw [hq, fst_inl]; exact one_ne_zero)
  let jet := (seriesJet (L := L)).restrictScalars K |>.comp f
  let U := normUnit a b
  have hU : jet (U : NormComplement a b) = 1 := hfu.trans (hq.trans (by simp))
  have hi : jet (↑(U⁻¹) : NormComplement a b) = 1 := by
    have he := congrArg jet (Units.inv_mul U)
    simpa only [map_mul, hU, mul_one, map_one] using he
  have hS : jet (height a b) = inl d + inr 1 := by
    rw [height, map_mul, hi, mul_one]
    exact hfB.trans hB
  let factor := 1 + ((U : NormComplement a b) * height a b) * ((M - 1 : ℕ) + (M : NormComplement a b) * U)
  have hfactor : jet factor = inr c := by
    dsimp only [factor]
    simp only [map_add, map_mul, map_one, map_natCast, hU, hS, one_mul, mul_one]
    rw [← inl_natCast (M - 1), ← inl_natCast M]
    apply TrivSqZeroExt.ext
    · simp only [fst_add, fst_one, fst_mul, fst_inl, fst_inr, add_zero]
      change 1 + d * c = 0
      rw [hdc, add_neg_cancel]
    · simp only [snd_add, snd_one, DualNumber.snd_mul, fst_add, fst_inl, fst_inr, snd_inl,
        snd_inr, add_zero, zero_add, mul_zero, one_mul]
      rfl
  have hbase : (jet (productPhase U (height a b) M M M)).fst ≠ 0 := by
    simp only [productPhase, map_mul, map_pow, map_add, map_sub, map_one, hU, hi, hS,
      fst_mul, fst_pow, fst_sub, fst_add, fst_one, fst_inl, fst_inr, add_zero]
    apply mul_ne_zero (mul_ne_zero (pow_ne_zero _ one_ne_zero)
      (pow_ne_zero _ (by norm_num)))
    rw [show (1 + 1 : L) ^ 2 = 4 by norm_num]
    exact pow_ne_zero M hw
  apply nonpower_of_dualNumber (extendJet (L := L) jet).toRingHom
  · change ((extendJet jet) (Algebra.TensorProduct.includeRight _)).fst = 0
    rw [extendJet_includeRight]
    change (jet (productPhase U (height a b) M M M * factor)).fst = 0
    simp [hfactor]
  · change ((extendJet jet) (Algebra.TensorProduct.includeRight _)).snd ≠ 0
    rw [extendJet_includeRight]
    change (jet (productPhase U (height a b) M M M * factor)).snd ≠ 0
    simpa only [map_mul, hfactor, DualNumber.snd_mul, fst_inr, snd_inr, mul_zero, add_zero] using
      mul_ne_zero hbase hc

omit [CharZero K] in
/-- The axis restriction of the simple-pole witness has a simple
zero at x = −1. -/
theorem simplePolePhase_absolute_nonpower (a : Fin 2 → K) (b : Fin M → K)
    (i : Fin M) (hb : b i ≠ 0) :
    let x := algebraMap (QuadPolynomial K M) (NormComplement a b) (MvPolynomial.X (Sum.inr i))
    let B := algebraMap (QuadPolynomial K M) (NormComplement a b) (planeNorm a)
    NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := L)
      (simplePolePhase (normUnit a b) x B)) := by
  intro x B
  let D := NormComplement a b
  let R := QuadPolynomial K M
  let U := normUnit a b
  let t : L[ε] := inl (-1) + inr 1
  let f : R →ₐ[K] L[ε] := MvPolynomial.aeval (blockJetCoordinates i 0 t)
  have hq : (f (quadNorm a b)).fst = algebraMap K L (b i) := by
    change ((MvPolynomial.aeval (blockJetCoordinates i 0 t)) (quadNorm a b)).fst = _
    rw [blockJet_norm]
    simp [t]
  have hunit : IsUnit (f (quadNorm a b)) := by
    rw [TrivSqZeroExt.isUnit_iff_isUnit_fst, hq]
    exact (isUnit_iff_ne_zero.2 hb).map (algebraMap K L)
  let jet : D →ₐ[K] L[ε] := IsLocalization.Away.liftAlgHom (quadNorm a b) hunit
  have hmap (p : R) : jet (algebraMap R D p) = f p := IsLocalization.Away.lift_eq _ hunit _
  have hx : jet x = t := by
    change jet (algebraMap R D (MvPolynomial.X (Sum.inr i))) = t
    rw [hmap]
    exact (MvPolynomial.aeval_X _ _).trans (by simp [blockJetCoordinates])
  have hB : jet B = 0 := by
    change jet (algebraMap R D (planeNorm a)) = 0
    rw [hmap, planeNorm, map_sum]
    apply Finset.sum_eq_zero
    intro j _
    have hz : f (MvPolynomial.X (Sum.inl j)) = 0 :=
      (MvPolynomial.aeval_X _ _).trans (by simp [blockJetCoordinates])
    rw [map_mul, map_pow, hz, zero_pow (by decide), mul_zero]
  have hTwo : (2 : L[ε]).fst = 2 := map_ofNat (TrivSqZeroExt.fstHom L L L) 2
  have hinv : (jet (↑(U⁻¹) : D)).fst ≠ 0 :=
    (((U⁻¹).isUnit.map jet).map (TrivSqZeroExt.fstHom L L L)).ne_zero
  have he : jet (simplePolePhase U x B) =
      jet (↑(U⁻¹) : D) * (1 + t) * (1 + 2 * t) := by
    simp only [simplePolePhase, map_mul, map_add, map_sub, map_ofNat, map_one, hx, hB, sub_zero]
  have hfactor : (1 + 2 * t : L[ε]).fst = -1 := by
    rw [fst_add, fst_one, fst_mul, hTwo]
    simp only [t, fst_add, fst_inl, fst_inr, add_zero]
    norm_num
  apply nonpower_of_dualNumber (extendJet (L := L) jet).toRingHom
  · change (extendJet jet (Algebra.TensorProduct.includeRight _)).fst = 0
    rw [extendJet_includeRight, he, fst_mul, fst_mul]
    simp [t]
  · change (extendJet jet (Algebra.TensorProduct.includeRight _)).snd ≠ 0
    rw [extendJet_includeRight, he, DualNumber.snd_mul, DualNumber.snd_mul, fst_mul, hfactor]
    simp only [t, fst_add, fst_one, fst_inl, fst_inr, add_zero, add_neg_cancel, mul_zero,
      snd_add, snd_inl, snd_inr, snd_one, zero_mul, zero_add, mul_one, add_zero]
    exact mul_ne_zero hinv (neg_ne_zero.mpr one_ne_zero)

end MathieuQuadratic
