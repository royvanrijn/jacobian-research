import MathieuQuadratic.DraftBounds

/-!
# Weight bounds and eventual marked vanishing
-/

namespace MathieuQuadratic

variable {K A : Type*} [Field K] [CharZero K] [CommRing A] [Algebra K A]

omit [CharZero K] in
/-- The least-weight estimate in the instability lemma. -/
theorem trace_vanishes_of_positive_weight_bound (trace : A →ₗ[K] K)
    (W : ℤ → Submodule K A) (hone : (1 : A) ∈ W 0)
    (hmul : ∀ {d e : ℤ} {f h : A}, f ∈ W d → h ∈ W e → f * h ∈ W (d + e))
    (htrace : ∀ {d : ℤ} {f : A}, 0 < d → f ∈ W d → trace f = 0)
    {f h : A} {δ b : ℤ} (hf : f ∈ W δ) (hh : h ∈ W b) (m : ℕ)
    (hm : 0 < (m : ℤ) * δ + b) :
    trace (h * f ^ m) = 0 := by
  have hp (n : ℕ) : f ^ n ∈ W ((n : ℤ) * δ) := by
    induction n with
    | zero => simpa using hone
    | succ n ih =>
      rw [pow_succ]
      simpa [Int.natCast_succ, add_mul] using hmul ih hf
  have hpos : 0 < b + (m : ℤ) * δ := by omega
  exact htrace hpos (hmul hh (hp m))

omit [CharZero K] in
/-- A positive weight bound forces both purity and eventual vanishing
against every fixed multiplier with a lower weight bound. -/
theorem pure_eventual_of_weight_bounds (trace : A →ₗ[K] K)
    (W : ℤ → Submodule K A) (hone : (1 : A) ∈ W 0)
    (hmul : ∀ {d e : ℤ} {f h : A}, f ∈ W d → h ∈ W e → f * h ∈ W (d + e))
    (htrace : ∀ {d : ℤ} {f : A}, 0 < d → f ∈ W d → trace f = 0)
    (bounded : ∀ h : A, ∃ b : ℤ, h ∈ W b)
    {f : A} {δ : ℤ} (hf : f ∈ W δ) (hδ : 0 < δ) :
    Pure trace f ∧ ∀ h : A, ∃ m₀ : ℕ, ∀ m ≥ m₀, trace (h * f ^ m) = 0 := by
  constructor
  · intro m hm
    simpa using trace_vanishes_of_positive_weight_bound trace W hone hmul htrace hf hone m
      (by simpa using (mul_pos (by exact_mod_cast hm : 0 < (m : ℤ)) hδ))
  · intro h
    obtain ⟨b, hb⟩ := bounded h
    refine ⟨b.natAbs + 1, fun m hm ↦ trace_vanishes_of_positive_weight_bound trace W hone hmul htrace hf hb m ?_⟩
    have hmul := mul_le_mul_of_nonneg_left (show (1 : ℤ) ≤ δ by omega) (Int.natCast_nonneg m)
    rw [mul_one] at hmul
    omega

/-- The product phase used in the draft to separate eventual marked
vanishing from geometric instability. -/
def nullconeSeparationPhase (U : Aˣ) (S : A) (M : ℕ) : A :=
  productPhase U S (M + 1) M (2 * M + 1)

/-- The pure but eventually vanishing part of the strict-inclusion
example, using the regular block averaging input. -/
theorem nullconeSeparationPhase_pure_eventual (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
    (M : ℕ) (hM : 0 < M)
    (weight : ∀ (a : ℤ) (i : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ i) = 0)
    (height : ∀ i : ℕ, trace (S ^ i) = (M : K) / (↑(M + 2 * i) : K))
    (average : ∀ h : A, ∃ (support : Finset (ℤ × ℕ)) (c : ℤ × ℕ → K),
      ∀ m : ℕ, 0 < m → trace (h * nullconeSeparationPhase U S M ^ m) =
        trace ((∑ p ∈ support, c p • ((↑(U ^ p.1) : A) * S ^ p.2)) *
          nullconeSeparationPhase U S M ^ m)) :
    Pure trace (nullconeSeparationPhase U S M) ∧
      ∀ h : A, ∃ m₀ : ℕ, ∀ m ≥ m₀, trace (h * nullconeSeparationPhase U S M ^ m) = 0 := by
  constructor
  · exact (product_pure_iff trace U S hM weight height (by omega : M + 1 ≤ 2 * M + 1)).2
      ⟨le_rfl, by omega⟩
  · intro h
    by_contra hn
    have he := (product_family_detects_iff trace U S hM weight height le_rfl
      (by omega : M ≤ M + 1) (by omega : M + 1 ≤ 2 * M + 1) average).1 ⟨h, hn⟩
    omega

/-- The two nonzero extreme moments used in the separating invariant. -/
theorem nullconeSeparationPhase_extreme_moments (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
    (M : ℕ) (hM : 0 < M)
    (weight : ∀ (a : ℤ) (i : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ i) = 0)
    (height : ∀ i : ℕ, trace (S ^ i) = (M : K) / (↑(M + 2 * i) : K)) :
    let f := nullconeSeparationPhase U S M
    trace ((↑(U ^ (↑(M + 1 : ℕ) : ℤ)) : A) * f) =
      algebraMap ℚ K (mixedBetaMoment M (2 * M + 1) 0) ∧
      trace ((↑(U ^ (-(↑(4 * M + 1 : ℕ) : ℤ))) : A) * f) =
        -((M : K) / (↑(M + 2 * (2 * M + 1)) : K)) := by
  intro f
  constructor
  · have he := product_mixedMoment trace U S M (M + 1) M (2 * M + 1) 1 0 hM
      (↑(M + 1 : ℕ) : ℤ) (by omega) weight height
    simp only [mul_one, pow_one, pow_zero, mul_zero, add_zero, le_rfl, if_true,
      sub_self, Int.toNat_zero, Ring.choose_zero_right, Int.cast_one] at he
    simpa [f, nullconeSeparationPhase] using he
  · have he := product_highest_moment trace U S M (M + 1) M (2 * M + 1) hM weight height
    rw [show ((↑(M + 1 : ℕ) : ℤ) - (M + 2 * (2 * M + 1) : ℕ)) =
      -(↑(4 * M + 1 : ℕ) : ℤ) by omega] at he
    rw [show (-1 : ℚ) ^ (2 * M + 1) = -1 by rw [pow_succ, pow_mul]; norm_num] at he
    simpa [f, nullconeSeparationPhase] using he

/-- The polynomial in the two extreme trace functionals used to
separate the draft's phase from the nullcone. -/
noncomputable def extremeMomentInvariant (trace : A →ₗ[K] K) (U : Aˣ) (M : ℕ) (f : A) : K :=
  trace ((↑(U ^ ((M + 1 : ℕ) : ℤ)) : A) * f) ^ (4 * M + 1) *
    trace ((↑(U ^ (-(↑(4 * M + 1 : ℕ) : ℤ))) : A) * f) ^ (M + 1)

omit [CharZero K] in
/-- The opposite characters cancel in the product of extreme moments. -/
theorem extremeMomentInvariant_eq_of_covariance (trace : A →ₗ[K] K) (U : Aˣ) (M : ℕ)
    (f g : A) (c : Kˣ)
    (h : ∀ ℓ : ℤ, trace ((↑(U ^ ℓ) : A) * g) =
      (↑(c ^ ℓ) : K) * trace ((↑(U ^ ℓ) : A) * f)) :
    extremeMomentInvariant trace U M g = extremeMomentInvariant trace U M f := by
  unfold extremeMomentInvariant
  rw [h, h, mul_pow, mul_pow]
  have hc : (↑(c ^ ((M + 1 : ℕ) : ℤ)) : K) ^ (4 * M + 1) *
      (↑(c ^ (-(↑(4 * M + 1 : ℕ) : ℤ))) : K) ^ (M + 1) = 1 := by
    rw [← Units.val_pow_eq_pow_val, ← Units.val_pow_eq_pow_val, ← Units.val_mul,
      ← zpow_natCast _ (4 * M + 1), ← zpow_natCast _ (M + 1), ← zpow_mul, ← zpow_mul, ← zpow_add]
    rw [show (↑(M + 1 : ℕ) : ℤ) * ↑(4 * M + 1 : ℕ) +
      -(↑(4 * M + 1 : ℕ) : ℤ) * ↑(M + 1 : ℕ) = 0 by ring]
    simp
  calc
    _ = ((↑(c ^ ((M + 1 : ℕ) : ℤ)) : K) ^ (4 * M + 1) *
        (↑(c ^ (-(↑(4 * M + 1 : ℕ) : ℤ))) : K) ^ (M + 1)) *
        extremeMomentInvariant trace U M f := by unfold extremeMomentInvariant; ring
    _ = _ := by rw [hc, one_mul]; rfl

/-- The separating invariant has a nonzero value on the explicit phase. -/
theorem extremeMomentInvariant_ne_zero (trace : A →ₗ[K] K) (U : Aˣ) (S : A)
    (M : ℕ) (hM : 0 < M)
    (weight : ∀ (a : ℤ) (i : ℕ), a ≠ 0 → trace ((↑(U ^ a) : A) * S ^ i) = 0)
    (height : ∀ i : ℕ, trace (S ^ i) = (M : K) / (↑(M + 2 * i) : K)) :
    extremeMomentInvariant trace U M (nullconeSeparationPhase U S M) ≠ 0 := by
  obtain ⟨h₁, h₂⟩ := nullconeSeparationPhase_extreme_moments trace U S M hM weight height
  unfold extremeMomentInvariant
  rw [h₁, h₂]
  have hn : algebraMap ℚ K (mixedBetaMoment M (2 * M + 1) 0) ≠ 0 :=
    (map_ne_zero_iff (algebraMap ℚ K) (RingHom.injective (algebraMap ℚ K))).mpr
      (mixedBetaMoment_pos M _ _ hM).ne'
  have hm : (M : K) ≠ 0 := by exact_mod_cast hM.ne'
  have hd : (↑(M + 2 * (2 * M + 1)) : K) ≠ 0 := by exact_mod_cast (show M + 2 * (2 * M + 1) ≠ 0 by omega)
  exact mul_ne_zero (pow_ne_zero _ hn) (pow_ne_zero _ (neg_ne_zero.mpr (div_ne_zero hm hd)))

omit [CharZero K] in
/-- An invariant with a closed nonzero fiber excludes zero from the
closure of any orbit contained in that fiber. -/
theorem zero_not_mem_orbit_closure_of_extreme_moments [TopologicalSpace A]
    (trace : A →ₗ[K] K) (U : Aˣ) (M : ℕ) (f : A) (orbit : Set A)
    (hInv : ∀ g ∈ orbit, extremeMomentInvariant trace U M g = extremeMomentInvariant trace U M f)
    (hclosed : IsClosed {g | extremeMomentInvariant trace U M g = extremeMomentInvariant trace U M f})
    (hne : extremeMomentInvariant trace U M f ≠ 0) :
    (0 : A) ∉ closure orbit := by
  intro hz
  have he := closure_minimal hInv hclosed hz
  apply hne
  have hi : extremeMomentInvariant trace U M 0 = 0 := by
    simp [extremeMomentInvariant, zero_pow (show 4 * M + 1 ≠ 0 by omega)]
  exact he.symm.trans hi

/-- The explicit separator of the geometric nullcone is geometrically
nonpower by the two coprime restriction orders. -/
theorem nullconeSeparationPhase_absolute_nonpower {M : ℕ} (a : Fin 2 → K) (b : Fin M → K)
    (ha : a 0 ≠ 0) (j : Fin M) (hb : b j ≠ 0) :
    NonpowerUpToUnits (Algebra.TensorProduct.includeRight (R := K) (A := AlgebraicClosure K)
      (nullconeSeparationPhase (normUnit a b) (height a b) M)) :=
  productPhase_absolute_nonpower a b ha j hb (M + 1) M (2 * M + 1)
    (by simp [two_mul, add_assoc])

end MathieuQuadratic
