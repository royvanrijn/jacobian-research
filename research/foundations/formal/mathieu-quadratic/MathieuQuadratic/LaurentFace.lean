import MathieuQuadratic.Laurent

/-!
# Constant terms on a supporting face

When all support weights are nonnegative, a product contributing to
weight zero uses only factors of weight zero. This preserves constant
terms of all powers on restriction to a supporting face through zero.
-/

open scoped Pointwise

namespace MathieuQuadratic

variable {K G W : Type*} [Field K] [AddCommGroup G] [AddCommGroup W] [LinearOrder W]
  [IsOrderedAddMonoid W]

/-- Keep the coefficients whose exponent has weight zero. -/
noncomputable def laurentZeroFace (w : G →+ W) (f : AddMonoidAlgebra K G) : AddMonoidAlgebra K G :=
  AddMonoidAlgebra.ofCoeff (f.coeff.filter (fun g ↦ w g = 0))

omit [IsOrderedAddMonoid W] in
theorem coeff_laurentZeroFace (w : G →+ W) (f : AddMonoidAlgebra K G) (g : G) :
    (laurentZeroFace w f).coeff g = if w g = 0 then f.coeff g else 0 := by
  classical
  exact Finsupp.filter_apply _ _ _

omit [IsOrderedAddMonoid W] in
theorem mem_support_laurentZeroFace (w : G →+ W) (f : AddMonoidAlgebra K G) (g : G) :
    g ∈ (laurentZeroFace w f).coeff.support ↔ g ∈ f.coeff.support ∧ w g = 0 := by
  classical
  rw [Finsupp.mem_support_iff, Finsupp.mem_support_iff, coeff_laurentZeroFace]
  by_cases h : w g = 0 <;> simp [h]

omit [LinearOrder W] [IsOrderedAddMonoid W] in
private theorem support_mul_property (w : G →+ W) (f h : AddMonoidAlgebra K G)
    (P Q R : W → Prop)
    (hf : ∀ g ∈ f.coeff.support, P (w g))
    (hh : ∀ g ∈ h.coeff.support, Q (w g))
    (hadd : ∀ x y, P x → Q y → R (x + y)) :
    ∀ g ∈ (f * h).coeff.support, R (w g) := by
  classical
  intro g hg
  obtain ⟨a, ha, b, hb, rfl⟩ :=
    Finset.mem_add.mp (AddMonoidAlgebra.support_coeff_mul_subset f h hg)
  rw [map_add]
  exact hadd _ _ (hf _ ha) (hh _ hb)

omit [LinearOrder W] [IsOrderedAddMonoid W] in
private theorem support_pow_property (w : G →+ W) (f : AddMonoidAlgebra K G)
    (P : W → Prop) (h0 : P 0) (hadd : ∀ x y, P x → P y → P (x + y))
    (hf : ∀ g ∈ f.coeff.support, P (w g)) (n : ℕ) :
    ∀ g ∈ (f ^ n).coeff.support, P (w g) := by
  classical
  induction n with
  | zero =>
    intro g hg
    have hg0 : g = 0 := by
      by_contra hn
      simp [AddMonoidAlgebra.one_def, AddMonoidAlgebra.coeff_single, hn] at hg
    rwa [hg0, map_zero]
  | succ n hn =>
    rw [pow_succ]
    exact support_mul_property w _ _ P P P hn hf hadd

omit [LinearOrder W] [IsOrderedAddMonoid W] in
/-- Powers of a phase supported on the zero-weight subgroup stay there. -/
theorem support_pow_of_weight_zero (w : G →+ W) (f : AddMonoidAlgebra K G)
    (hf : ∀ g ∈ f.coeff.support, w g = 0) (n : ℕ) :
    ∀ g ∈ (f ^ n).coeff.support, w g = 0 :=
  support_pow_property w f (fun x ↦ x = 0) rfl (fun _ _ hx hy ↦ by rw [hx, hy, add_zero]) hf n

omit [IsOrderedAddMonoid W] in
private theorem support_sub_zeroFace (w : G →+ W) (f : AddMonoidAlgebra K G)
    (hf : ∀ g ∈ f.coeff.support, 0 ≤ w g) :
    ∀ g ∈ (f - laurentZeroFace w f).coeff.support, 0 < w g := by
  classical
  intro g hg
  rw [Finsupp.mem_support_iff, AddMonoidAlgebra.coeff_sub, Finsupp.sub_apply,
    coeff_laurentZeroFace] at hg
  by_cases h : w g = 0
  · simp [h] at hg
  · rw [if_neg h, sub_zero] at hg
    exact lt_of_le_of_ne (hf g (Finsupp.mem_support_iff.mpr hg)) (Ne.symm h)

/-- Restriction to the zero face preserves each pure-power constant
coefficient when the original support has nonnegative weights. -/
theorem laurentCT_pow_zeroFace (w : G →+ W) (f : AddMonoidAlgebra K G)
    (hf : ∀ g ∈ f.coeff.support, 0 ≤ w g) (n : ℕ) :
    laurentCT (f ^ n) = laurentCT (laurentZeroFace w f ^ n) := by
  classical
  let z := laurentZeroFace w f
  have hz (m : ℕ) : ∀ g ∈ (z ^ m).coeff.support, w g = 0 :=
    support_pow_of_weight_zero w z (fun g hg ↦ (mem_support_laurentZeroFace w f g |>.mp hg).2) m
  have hd : ∀ m : ℕ, ∀ g ∈ (f ^ m - z ^ m).coeff.support, 0 < w g := by
    intro m
    induction m with
    | zero => simp
    | succ m ih =>
      have he : f ^ (m + 1) - z ^ (m + 1) =
          (f ^ m - z ^ m) * f + z ^ m * (f - z) := by
        rw [pow_succ, pow_succ]
        ring
      rw [he]
      intro g hg
      have hs₁ := support_mul_property w (f ^ m - z ^ m) f (0 < ·) (0 ≤ ·) (0 < ·)
        ih hf (fun _ _ hx hy ↦ add_pos_of_pos_of_nonneg hx hy)
      have hs₂ := support_mul_property w (z ^ m) (f - z) (· = 0) (0 < ·) (0 < ·)
        (hz m) (support_sub_zeroFace w f hf) (fun _ _ hx hy ↦ by simpa [hx] using hy)
      rcases Finset.mem_union.mp (Finsupp.support_add hg) with h | h
      · exact hs₁ g h
      · exact hs₂ g h
  have hv : laurentCT (f ^ n - z ^ n) = 0 := by
    change (f ^ n - z ^ n).coeff 0 = 0
    apply Finsupp.notMem_support_iff.mp
    intro h
    simpa using hd n 0 h
  exact sub_eq_zero.mp (by simpa only [map_sub] using hv)

omit [LinearOrder W] [IsOrderedAddMonoid W] in
/-- A coordinate projection injective on the zero-weight subgroup
preserves its Laurent constant term after taking powers. -/
theorem coeff_mapDomain_of_weight_zero {H : Type*} [AddCommGroup H]
    (w : G →+ W) (π : G →+ H) (hinj : Set.InjOn π {g | w g = 0})
    (f : AddMonoidAlgebra K G) (hf : ∀ g ∈ f.coeff.support, w g = 0)
    (x : G) (hx : w x = 0) :
    (AddMonoidAlgebra.mapDomain π f).coeff (π x) = f.coeff x := by
  exact Finsupp.mapDomain_apply' {g | w g = 0} f.coeff (fun g hg ↦ hf g hg) hinj hx

omit [LinearOrder W] [IsOrderedAddMonoid W] in
theorem laurentCT_pow_mapDomain_of_weight_zero {H : Type*} [AddCommGroup H]
    (w : G →+ W) (π : G →+ H) (hinj : Set.InjOn π {g | w g = 0})
    (f : AddMonoidAlgebra K G) (hf : ∀ g ∈ f.coeff.support, w g = 0) (n : ℕ) :
    laurentCT ((AddMonoidAlgebra.mapDomainRingHom K π f) ^ n) = laurentCT (f ^ n) := by
  rw [← map_pow]
  change (AddMonoidAlgebra.mapDomain π (f ^ n)).coeff 0 = (f ^ n).coeff 0
  rw [← π.map_zero]
  exact coeff_mapDomain_of_weight_zero w π hinj (f ^ n)
    (support_pow_of_weight_zero w f hf n) 0 w.map_zero

end MathieuQuadratic
