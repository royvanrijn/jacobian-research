import MathieuQuadratic.IsotropicTernary

/-!
# The norm-one witness on the hyperbolic chart
-/

open MvPolynomial TrivSqZeroExt
open scoped BigOperators DualNumber TensorProduct

namespace MathieuQuadratic

variable {K : Type*} [Field K] [CharZero K] {M : ℕ}
variable (b : Fin M → K)

noncomputable def normOneRelation : QuadPolynomial K M :=
  2 * (X (Sum.inl 0) * X (Sum.inl 1)) + heightNorm b - 1

noncomputable def normOneIdeal : Ideal (QuadPolynomial K M) :=
  Ideal.span {normOneRelation b}

abbrev NormOneChart := QuadPolynomial K M ⧸ normOneIdeal b

noncomputable def normOneX : NormOneChart b :=
  algebraMap (QuadPolynomial K M) _ (X (Sum.inl 0))

noncomputable def normOneY : NormOneChart b :=
  algebraMap (QuadPolynomial K M) _ (X (Sum.inl 1))

noncomputable def normOneHeight : NormOneChart b :=
  algebraMap (QuadPolynomial K M) _ (heightNorm b)

omit [CharZero K] in
theorem normOneChart_relation :
    2 * (normOneX b * normOneY b) + normOneHeight b = 1 := by
  have h : algebraMap (QuadPolynomial K M) (NormOneChart b) (normOneRelation b) = 0 :=
    Ideal.Quotient.eq_zero_iff_mem.2 (Ideal.subset_span (Set.mem_singleton _))
  rw [normOneRelation, map_sub, map_add, map_mul, map_mul, map_ofNat, map_one, sub_eq_zero] at h
  exact h

private noncomputable def normOneTorusPolynomial (c : Kˣ) :
    QuadPolynomial K M →ₐ[K] QuadPolynomial K M :=
  aeval (Sum.elim ![C (c : K) * X (Sum.inl 0), C (↑(c⁻¹) : K) * X (Sum.inl 1)]
    (fun i ↦ X (Sum.inr i)))

omit [CharZero K] in
private theorem normOneTorusPolynomial_relation (c : Kˣ) :
    normOneTorusPolynomial (M := M) c (normOneRelation b) = normOneRelation b := by
  have hc : C (c : K) * C (↑(c⁻¹) : K) = (1 : QuadPolynomial K M) := by
    rw [← map_mul, Units.mul_inv, map_one]
  simp only [normOneRelation, heightNorm, map_sub, map_add, map_mul, map_sum, map_pow,
    map_one, map_ofNat, normOneTorusPolynomial, aeval_X, aeval_C, algebraMap_eq,
    Sum.elim_inl, Sum.elim_inr, Matrix.cons_val_zero, Matrix.cons_val_one]
  rw [mul_mul_mul_comm, hc, one_mul]

/-- Pullback by the hyperbolic torus on the norm-one quotient. -/
noncomputable def normOneTorusPullback (c : Kˣ) : NormOneChart b →ₐ[K] NormOneChart b :=
  Ideal.Quotient.liftₐ _ ((IsScalarTower.toAlgHom K (QuadPolynomial K M) (NormOneChart b)).comp
    (normOneTorusPolynomial c)) (by
      intro p hp
      obtain ⟨w, hw⟩ := Ideal.mem_span_singleton.mp hp
      have hzero : algebraMap (QuadPolynomial K M) (NormOneChart b) (normOneRelation b) = 0 :=
        Ideal.Quotient.eq_zero_iff_mem.mpr (Ideal.subset_span (Set.mem_singleton _))
      rw [hw, AlgHom.comp_apply, map_mul, normOneTorusPolynomial_relation,
        IsScalarTower.toAlgHom_apply, map_mul, hzero, zero_mul])

omit [CharZero K] in
private theorem normOneTorusPullback_X (c : Kˣ) :
    normOneTorusPullback b c (normOneX b) =
      algebraMap K (NormOneChart b) (c : K) * normOneX b := by
  change algebraMap (QuadPolynomial K M) (NormOneChart b)
    (normOneTorusPolynomial c (X (Sum.inl 0))) = _
  simp only [normOneTorusPolynomial, aeval_X, Sum.elim_inl, Matrix.cons_val_zero, map_mul,
    ← algebraMap_eq, ← IsScalarTower.algebraMap_apply, normOneX]

omit [CharZero K] in
private theorem normOneTorusPullback_Y (c : Kˣ) :
    normOneTorusPullback b c (normOneY b) =
      algebraMap K (NormOneChart b) (↑(c⁻¹) : K) * normOneY b := by
  change algebraMap (QuadPolynomial K M) (NormOneChart b)
    (normOneTorusPolynomial c (X (Sum.inl 1))) = _
  simp only [normOneTorusPolynomial, aeval_X, Sum.elim_inl, Matrix.cons_val_one, Matrix.cons_val_zero, map_mul,
    ← algebraMap_eq, ← IsScalarTower.algebraMap_apply, normOneY]

/-- Invariance under one hyperbolic torus element kills unequal exponents. -/
theorem normOne_weight_of_torus (trace : NormOneChart b →ₗ[K] K)
    (invariant : ∀ g, trace (normOneTorusPullback b (Units.mk0 (2 : K) two_ne_zero) g) = trace g)
    (i j : ℕ) (hij : i ≠ j) : trace (normOneX b ^ i * normOneY b ^ j) = 0 := by
  let c : Kˣ := Units.mk0 2 two_ne_zero
  let g := normOneX b ^ i * normOneY b ^ j
  have he : normOneTorusPullback b c g =
      algebraMap K (NormOneChart b) (↑(c ^ i * (c⁻¹) ^ j) : K) * g := by
    dsimp only [g]
    rw [map_mul, map_pow, map_pow, normOneTorusPullback_X, normOneTorusPullback_Y,
      Units.val_mul, Units.val_pow_eq_pow_val, Units.val_pow_eq_pow_val,
      map_mul, map_pow, map_pow, mul_pow, mul_pow, mul_mul_mul_comm]
  have ht := invariant g
  rw [he, ← Algebra.smul_def, map_smul, smul_eq_mul] at ht
  have hd : (↑(c ^ i * (c⁻¹) ^ j) : K) ≠ 1 := by
    intro hf
    have hu : c ^ i * (c⁻¹) ^ j = 1 := Units.ext hf
    have hp : (c ^ i : Kˣ) = c ^ j := by
      rw [inv_pow] at hu
      exact eq_of_mul_inv_eq_one hu
    have hK : (2 : K) ^ i = (2 : K) ^ j := by
      simpa only [Units.val_pow_eq_pow_val, Units.val_mk0, c] using congrArg Units.val hp
    have hQ : algebraMap ℚ K ((2 : ℚ) ^ i) = algebraMap ℚ K ((2 : ℚ) ^ j) := by simpa using hK
    exact hij ((pow_right_strictMono₀ (by norm_num : (1 : ℚ) < 2)).injective
      ((RingHom.injective (algebraMap ℚ K)) hQ))
  have hz : ((↑(c ^ i * (c⁻¹) ^ j) : K) - 1) * trace g = 0 := by
    rw [sub_mul, one_mul, ht, sub_self]
  exact (mul_eq_zero.mp hz).resolve_left (sub_ne_zero.mpr hd)

/-- The plane in squared coordinates u = x + y and v = x - y. -/
private noncomputable def normOnePlane : Fin 2 → K :=
  ![(1 : K) / 2, -(1 / 2)]

private noncomputable def normOneDiagonalEval : QuadPolynomial K M →ₐ[K] NormOneChart b :=
  aeval (Sum.elim ![normOneX b + normOneY b, normOneX b - normOneY b]
    (fun i ↦ algebraMap (QuadPolynomial K M) (NormOneChart b) (X (Sum.inr i))))

private theorem normOneDiagonalEval_plane :
    normOneDiagonalEval b (planeNorm (M := M) (normOnePlane (K := K))) =
      1 - normOneHeight b := by
  have hhalf : algebraMap K (NormOneChart b) (1 / 2) * 2 = 1 := by
    rw [← map_ofNat (algebraMap K (NormOneChart b)) 2, ← map_mul]
    norm_num
  rw [planeNorm, map_sum, Fin.sum_univ_two]
  simp only [map_mul, map_pow, normOneDiagonalEval, aeval_C, aeval_X, Sum.elim_inl,
    normOnePlane, Matrix.cons_val_zero, Matrix.cons_val_one, map_neg]
  calc
    _ = (algebraMap K (NormOneChart b) (1 / 2) * 2) * (2 * (normOneX b * normOneY b)) := by ring
    _ = _ := by rw [hhalf, one_mul]; exact eq_sub_of_add_eq (normOneChart_relation b)

omit [CharZero K] in
private theorem normOneDiagonalEval_height :
    normOneDiagonalEval b (heightNorm b) = normOneHeight b := by
  simp only [normOneHeight, heightNorm, map_sum, map_mul, map_pow, normOneDiagonalEval,
    aeval_X, Sum.elim_inr, ← algebraMap_eq, ← IsScalarTower.algebraMap_apply, AlgHom.commutes]

/-- The normalized trace annihilates the restrictions of infinitesimal
rotations from the squared coordinates. -/
def NormOneRotationInvariant (trace : NormOneChart b →ₗ[K] K) : Prop :=
  ∀ (i : Fin 2) (j : Fin M) (f : QuadPolynomial K M),
    trace (normOneDiagonalEval b (rotation normOnePlane b i j f)) = 0

/-- The norm-one height formula follows from the same summed rotation
identity as for the complement. -/
theorem normOne_height_of_rotations (trace : NormOneChart b →ₗ[K] K)
    (normalized : trace 1 = 1) (rotations : NormOneRotationInvariant b trace)
    (hM : 0 < M) (j : ℕ) :
    trace (normOneHeight b ^ j) = (M : K) / (↑(M + 2 * j) : K) := by
  apply heightMoment_of_recurrence trace _ M hM normalized _ j
  intro n
  let f := normOneDiagonalEval b
  let S := normOneHeight b
  have hzero : trace (f (∑ i : Fin 2, ∑ j : Fin M,
      rotation normOnePlane b i j (heightNorm b ^ n * (X (Sum.inl i) * X (Sum.inr j))))) = 0 := by
    simp only [map_sum]
    exact Finset.sum_eq_zero fun i _ ↦ Finset.sum_eq_zero fun j _ ↦ rotations i j _
  have hsum := sum_rotation_numerators (a := normOnePlane (K := K)) b n
  change (∑ i : Fin 2, ∑ j : Fin M,
    rotation normOnePlane b i j (heightNorm b ^ n * (X (Sum.inl i) * X (Sum.inr j)))) = _ at hsum
  rw [hsum] at hzero
  simp only [map_sub, map_mul, map_pow, map_add, map_natCast, map_ofNat,
    show f (planeNorm normOnePlane) = 1 - S from normOneDiagonalEval_plane b,
    show f (heightNorm b) = S from normOneDiagonalEval_height b] at hzero
  have hfor : (2 : NormOneChart b) * S ^ (n + 1) -
      ((M : NormOneChart b) + 2 * (n : NormOneChart b)) * (1 - S) * S ^ n =
      algebraMap K (NormOneChart b) (↑(M + 2 * (n + 1))) * S ^ (n + 1) -
        algebraMap K (NormOneChart b) (↑(M + 2 * n)) * S ^ n := by
    simp only [Nat.cast_add, Nat.cast_mul, Nat.cast_succ, map_add, map_mul, map_one,
      map_natCast, pow_succ]
    ring
  rw [← map_sub] at hzero
  rwa [hfor, map_sub, ← Algebra.smul_def, ← Algebra.smul_def,
    map_smul, map_smul, smul_eq_mul, smul_eq_mul, sub_eq_zero] at hzero

/-- The marked moments for the norm-one corollary in the draft. -/
theorem draft_normOne_moments (trace : NormOneChart b →ₗ[K] K) (hM : 0 < M)
    (weight : ∀ i j : ℕ, i ≠ j → trace (normOneX b ^ i * normOneY b ^ j) = 0)
    (height : ∀ j : ℕ, trace (normOneHeight b ^ j) = (M : K) / (↑(M + 2 * j) : K))
    (m : ℕ) (hm : 0 < m) :
    trace (normOnePhase (normOneX b) (normOneY b) M ^ m) = 0 ∧
      trace (normOneX b ^ M * normOnePhase (normOneX b) (normOneY b) M ^ m) =
        algebraMap ℚ K (mixedBetaMoment M ((M + 1) * m) 0) ∧
      trace (normOneX b ^ M * normOnePhase (normOneX b) (normOneY b) M ^ m) ≠ 0 := by
  letI : Algebra ℚ (NormOneChart b) := Algebra.compHom _ (algebraMap ℚ K)
  have hc : algebraMap ℚ (NormOneChart b) (1 / 2) * 2 = 1 := by
    rw [← map_ofNat (algebraMap ℚ (NormOneChart b)) 2, ← map_mul]
    norm_num
  have hr : normOneX b * normOneY b =
      algebraMap ℚ (NormOneChart b) (1 / 2) * (1 - normOneHeight b) := by
    calc
      _ = (algebraMap ℚ (NormOneChart b) (1 / 2) * 2) * (normOneX b * normOneY b) := by
        rw [hc, one_mul]
      _ = _ := by
        rw [mul_assoc, show 2 * (normOneX b * normOneY b) = 1 - normOneHeight b from
          eq_sub_of_add_eq (normOneChart_relation b)]
  exact normOnePhase_moments trace _ _ _ M hM (fun _ ↦ rfl) hr weight height m hm

/-- On the norm-one level, setting y = 0 and fixing B = 1 gives a
simple zero in the x direction. -/
theorem normOnePhase_absolute_nonpower (i : Fin M) (hi : b i ≠ 0) :
    NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := AlgebraicClosure K)
      (normOnePhase (normOneX b) (normOneY b) M)) := by
  classical
  let D := NormOneChart b
  let R := QuadPolynomial K M
  let c := algebraMap K (AlgebraicClosure K) (b i)
  have hc : c ≠ 0 := (map_ne_zero_iff (algebraMap K (AlgebraicClosure K)) (RingHom.injective _)).mpr hi
  obtain ⟨r, hr⟩ := IsAlgClosed.exists_pow_nat_eq (1 / c) (n := 2) (by decide)
  let v : Fin 2 ⊕ Fin M → (AlgebraicClosure K)[ε] :=
    Sum.elim ![ε, 0] (fun j ↦ if j = i then inl r else 0)
  let e : R →ₐ[K] (AlgebraicClosure K)[ε] := aeval v
  have heX (j : Fin 2 ⊕ Fin M) : e (X j) = v j := aeval_X _ _
  have hx : e (X (Sum.inl 0)) = ε := by rw [heX]; rfl
  have hy : e (X (Sum.inl 1)) = 0 := by rw [heX]; rfl
  have hB : e (heightNorm b) = 1 := by
    rw [heightNorm, map_sum, Finset.sum_eq_single i]
    · rw [map_mul, map_pow, heX]
      change e (C (b i)) * (if i = i then inl r else 0) ^ 2 = 1
      rw [if_pos rfl, show e (C (b i)) = algebraMap K (AlgebraicClosure K)[ε] (b i) from aeval_C _ _,
        IsScalarTower.algebraMap_apply K (AlgebraicClosure K) (AlgebraicClosure K)[ε]]
      change (inl c : (AlgebraicClosure K)[ε]) * (inl r) ^ 2 = 1
      apply TrivSqZeroExt.ext
      · rw [fst_mul, fst_inl, fst_pow, fst_inl, hr, mul_div_cancel₀ _ hc, fst_one]
      · simp [pow_two]
    · intro j _ hj
      rw [map_mul, map_pow, heX]
      change e (C (b j)) * (if j = i then inl r else 0) ^ 2 = 0
      rw [if_neg hj, zero_pow (by decide), mul_zero]
    · simp
  have hz : e (normOneRelation b) = 0 := by
    rw [normOneRelation, map_sub, map_add, map_mul, map_mul, map_ofNat, hx, hy, hB,
      mul_zero, mul_zero, zero_add, map_one, sub_self]
  let jet : D →ₐ[K] (AlgebraicClosure K)[ε] := Ideal.Quotient.liftₐ _ e (by
    intro p hp
    obtain ⟨w, hw⟩ := Ideal.mem_span_singleton.mp hp
    rw [hw, map_mul, hz, zero_mul])
  have hJ (p : R) : jet (algebraMap R D p) = e p := rfl
  have hx' : jet (normOneX b) = ε := (hJ _).trans hx
  have hy' : jet (normOneY b) = 0 := (hJ _).trans hy
  have hp : jet (normOnePhase (normOneX b) (normOneY b) M) =
      (ε : (AlgebraicClosure K)[ε]) * (1 + ε) ^ M * (-2 - ε) ^ (M + 1) := by
    simp only [normOnePhase, map_mul, map_pow, map_add, map_one, map_sub, map_ofNat,
      hx', hy', mul_zero, zero_sub]
  let extJet : (AlgebraicClosure K) ⊗[K] D →ₐ[(AlgebraicClosure K)] (AlgebraicClosure K)[ε] :=
    Algebra.TensorProduct.lift (Algebra.ofId (AlgebraicClosure K) (AlgebraicClosure K)[ε]) jet (fun _ _ ↦ Commute.all _ _)
  have hex : extJet (Algebra.TensorProduct.includeRight (normOnePhase (normOneX b) (normOneY b) M)) =
      jet (normOnePhase (normOneX b) (normOneY b) M) := by
    dsimp only [extJet]
    rw [Algebra.TensorProduct.includeRight_apply, Algebra.TensorProduct.lift_tmul, map_one, one_mul]
  have h2 : (2 : (AlgebraicClosure K)[ε]).fst = 2 := map_ofNat (TrivSqZeroExt.fstHom (AlgebraicClosure K) (AlgebraicClosure K) (AlgebraicClosure K)) 2
  apply nonpower_of_dualNumber extJet.toRingHom
  · change (extJet (Algebra.TensorProduct.includeRight _)).fst = 0
    rw [hex, hp]
    simp
  · change (extJet (Algebra.TensorProduct.includeRight _)).snd ≠ 0
    rw [hex, hp, DualNumber.snd_mul, DualNumber.snd_mul, fst_mul, fst_pow, fst_pow,
      fst_sub, fst_neg, h2]
    simp only [DualNumber.fst_eps, DualNumber.snd_eps, fst_add, fst_one, add_zero, one_pow,
      zero_mul, add_zero, zero_add, one_mul, sub_zero]
    exact pow_ne_zero _ (neg_ne_zero.mpr two_ne_zero)

/-- The norm-one witness with torus and infinitesimal invariance. -/
theorem draft_normOne_witness (trace : NormOneChart b →ₗ[K] K) (hM : 0 < M)
    (i : Fin M) (hi : b i ≠ 0)
    (invariant : ∀ g, trace (normOneTorusPullback b (Units.mk0 (2 : K) two_ne_zero) g) = trace g)
    (normalized : trace 1 = 1) (rotations : NormOneRotationInvariant b trace) :
    let P := normOnePhase (normOneX b) (normOneY b) M
    NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := AlgebraicClosure K) P) ∧
      AllOrderWitness trace P (normOneX b ^ M) := by
  have hmoment := draft_normOne_moments b trace hM (normOne_weight_of_torus b trace invariant)
    (normOne_height_of_rotations b trace normalized rotations hM)
  exact ⟨normOnePhase_absolute_nonpower b i hi,
    ⟨fun m hm ↦ (hmoment m hm).1, fun m hm ↦ (hmoment m hm).2.2⟩⟩

end MathieuQuadratic
