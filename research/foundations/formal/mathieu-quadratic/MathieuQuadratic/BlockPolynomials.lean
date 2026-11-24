import MathieuQuadratic.QuadraticCoordinates
import MathieuQuadratic.Transitivity
import Mathlib.Algebra.MvPolynomial.Expand
import Mathlib.FieldTheory.IsAlgClosed.AlgebraicClosure

/-!
# Block-invariant polynomials and the degree filtration
-/

open MvPolynomial
open scoped BigOperators

namespace MathieuQuadratic

variable {K I : Type*} [Field K] [CharZero K]

noncomputable def blockCoordinateSignChange (i : I) : MvPolynomial I K →ₐ[K] MvPolynomial I K := by
  classical
  exact aeval fun j ↦ if j = i then -X j else X j

omit [CharZero K] in
private theorem signChange_monomial (i : I) (d : I →₀ ℕ) (c : K) :
    blockCoordinateSignChange i (monomial d c) = monomial d (((-1 : K) ^ d i) * c) := by
  classical
  rw [blockCoordinateSignChange, aeval_monomial, algebraMap_eq, monomial_eq]
  by_cases hi : i ∈ d.support
  · rw [Finsupp.prod, Finsupp.prod, ← Finset.mul_prod_erase _ _ hi,
      ← Finset.mul_prod_erase _ _ hi, if_pos rfl, neg_pow]
    have hprod :
        (∏ x ∈ d.support.erase i,
          (if x = i then -X x else X x : MvPolynomial I K) ^ d x) =
        ∏ x ∈ d.support.erase i, X x ^ d x := by
      refine Finset.prod_congr rfl fun j hj ↦ ?_
      rw [if_neg (Finset.ne_of_mem_erase hj)]
    rw [hprod, map_mul, map_pow, map_neg, map_one]
    ring
  · have hdi : d i = 0 := Finsupp.notMem_support_iff.1 hi
    rw [hdi, pow_zero, one_mul]
    simp only [Finsupp.prod]
    congr 1
    refine Finset.prod_congr rfl fun j hj ↦ ?_
    have hne : j ≠ i := by
      intro he
      subst j
      exact hi hj
    rw [if_neg hne]

omit [CharZero K] in
private theorem signChange_coeff (i : I) (d : I →₀ ℕ) (p : MvPolynomial I K) :
    (blockCoordinateSignChange i p).coeff d = (-1 : K) ^ d i * p.coeff d := by
  classical
  induction p using MvPolynomial.induction_on' with
  | add p q hp hq => rw [map_add, coeff_add, hp, hq, coeff_add, mul_add]
  | monomial e c =>
      rw [signChange_monomial, coeff_monomial, coeff_monomial]
      by_cases he : e = d <;> simp [he]

private theorem even_exponent_of_sign_invariant (p : MvPolynomial I K)
    (hs : ∀ i, blockCoordinateSignChange i p = p) (d : I →₀ ℕ) (hd : d ∈ p.support) (i : I) :
    2 ∣ d i := by
  classical
  have hc : p.coeff d ≠ 0 := mem_support_iff.1 hd
  have he := signChange_coeff i d p
  rw [hs i] at he
  have hv : (-1 : K) ^ d i = 1 := mul_right_cancel₀ hc (by simpa using he.symm)
  rw [← even_iff_two_dvd]
  by_contra h
  have hodd := Nat.not_even_iff_odd.mp h
  rw [hodd.neg_one_pow] at hv
  exact one_ne_zero (CharZero.neg_eq_self_iff.mp hv)

/-- Separately even polynomials are precisely polynomial expressions in
the squares of the coordinates. The degree of the reduced polynomial
is at most half the degree of the original polynomial. -/
theorem exists_expand_two_of_sign_invariant (p : MvPolynomial I K)
    (hs : ∀ i, blockCoordinateSignChange i p = p) :
    ∃ r : MvPolynomial I K, expand 2 r = p ∧ 2 * r.totalDegree = p.totalDegree := by
  classical
  let halve (d : I →₀ ℕ) := d.mapRange (· / 2) (by decide)
  let r := ∑ d ∈ p.support, monomial (halve d) (p.coeff d)
  have he (d : I →₀ ℕ) (hd : d ∈ p.support) : 2 • halve d = d := by
    ext i
    simp only [halve, Finsupp.mapRange_apply, Finsupp.smul_apply, smul_eq_mul]
    exact Nat.mul_div_cancel' (even_exponent_of_sign_invariant p hs d hd i)
  have hr : expand 2 r = p := by
    simp only [r, map_sum, expand_monomial]
    calc
      _ = ∑ d ∈ p.support, monomial d (p.coeff d) :=
        Finset.sum_congr rfl fun d hd ↦ by rw [he d hd]
      _ = p := by
        exact MvPolynomial.support_sum_monomial_coeff p
  exact ⟨r, hr, by rw [← hr, totalDegree_expand, mul_comm]⟩

section TwoBlocks

variable {M : ℕ} (a : Fin 2 → K) (b : Fin M → K)

private def axes (i : Fin M) {D : Type*} [Zero D] (v : Fin 2 → D) :
    Fin 2 ⊕ Fin M → D :=
  Sum.elim (fun j ↦ if j = 0 then v 0 else 0) (fun j ↦ if j = i then v 1 else 0)

private noncomputable def restrictAxes (i : Fin M) :
    QuadPolynomial K M →ₐ[K] MvPolynomial (Fin 2) K :=
  aeval (axes i X)

omit [CharZero K] in
private theorem aeval_restrictAxes {L : Type*} [CommRing L] [Algebra K L] (i : Fin M)
    (v : Fin 2 → L) (p : QuadPolynomial K M) :
    aeval v (restrictAxes i p) = aeval (axes i v) p := by
  rw [restrictAxes, comp_aeval_apply]
  congr with j
  cases j with
  | inl j => by_cases h : j = 0 <;> simp [axes, h]
  | inr j => by_cases h : j = i <;> simp [axes, h]

omit [CharZero K] in
private theorem aeval_axes_plane {L : Type*} [CommRing L] [Algebra K L] (i : Fin M)
    (v : Fin 2 → L) :
    aeval (axes i v) (planeNorm a) = algebraMap K L (a 0) * v 0 ^ 2 := by
  rw [planeNorm, map_sum, Finset.sum_eq_single (0 : Fin 2)]
  · simp [axes]
  · intro j _ hj
    simp [axes, hj]
  · simp

omit [CharZero K] in
private theorem aeval_axes_height {L : Type*} [CommRing L] [Algebra K L] (i : Fin M)
    (v : Fin 2 → L) :
    aeval (axes i v) (heightNorm b) = algebraMap K L (b i) * v 1 ^ 2 := by
  rw [heightNorm, map_sum, Finset.sum_eq_single i]
  · simp [axes]
  · intro j _ hj
    simp [axes, hj]
  · simp

omit [CharZero K] in
private theorem restrictAxes_totalDegree (i : Fin M) (p : QuadPolynomial K M) :
    (restrictAxes i p).totalDegree ≤ p.totalDegree := by
  classical
  have hp : restrictAxes i p =
      ∑ d ∈ p.support, restrictAxes i (monomial d (p.coeff d)) := by
    conv_lhs => rw [p.as_sum]
    rw [map_sum]
  rw [hp]
  apply totalDegree_finsetSum_le
  intro d hd
  rw [restrictAxes, aeval_monomial, Finsupp.prod, algebraMap_eq]
  refine (totalDegree_mul _ _).trans ?_
  rw [totalDegree_C, zero_add]
  refine (totalDegree_finsetProd _ _).trans ?_
  refine (Finset.sum_le_sum fun j (_ : j ∈ d.support) ↦ ?_).trans (le_totalDegree hd)
  refine (totalDegree_pow _ _).trans ?_
  cases j with
  | inl j => by_cases hj : j = 0 <;> simp [axes, hj]
  | inr j => by_cases hj : j = i <;> simp [axes, hj]

omit [CharZero K] in
private theorem aeval_axisSign {L : Type*} [CommRing L] [Algebra K L]
    (i : Fin 2) (v : Fin 2 → L) (p : MvPolynomial (Fin 2) K) :
    aeval v (blockCoordinateSignChange i p) = aeval (fun j ↦ if j = i then -v j else v j) p := by
  rw [blockCoordinateSignChange, comp_aeval_apply]
  congr with j
  by_cases h : j = i <;> simp [h]

/-- Invariance on the nonzero block norm levels over the algebraic closure. -/
def ConstantOnBlockNormLevels (p : QuadPolynomial K M) : Prop :=
  ∀ v w : Fin 2 ⊕ Fin M → AlgebraicClosure K,
    aeval v (planeNorm a) = aeval w (planeNorm a) →
    aeval v (heightNorm b) = aeval w (heightNorm b) →
    aeval v (planeNorm a) ≠ 0 → aeval v (heightNorm b) ≠ 0 →
    aeval v p = aeval w p

/-- The polynomial invariant in the block averaging argument is a
polynomial over K in A and B. Its degree in these two arguments is at
most half the degree in the original coordinates. -/
theorem blockPolynomial_of_constant_norm_levels (p : QuadPolynomial K M)
    (hi : ConstantOnBlockNormLevels a b p) (ha : a 0 ≠ 0) (i : Fin M) (hb : b i ≠ 0) :
    ∃ r : MvPolynomial (Fin 2) K,
      p = aeval ![C (a 0)⁻¹ * planeNorm a, C (b i)⁻¹ * heightNorm b] r ∧
      2 * r.totalDegree ≤ p.totalDegree := by
  let L := AlgebraicClosure K
  let f := algebraMap K L
  have hf : Function.Injective f := RingHom.injective f
  have ha' : f (a 0) ≠ 0 := (map_ne_zero_iff f hf).mpr ha
  have hb' : f (b i) ≠ 0 := (map_ne_zero_iff f hf).mpr hb
  let P := restrictAxes i p
  have hsign (k : Fin 2) : blockCoordinateSignChange k P = P := by
    apply MvPolynomial.map_injective f hf
    apply MvPolynomial.funext_set (s := fun _ ↦ {x : L | x ≠ 0})
    · intro j
      exact (Set.infinite_univ.sdiff (Set.finite_singleton (0 : L))).mono
        (fun x hx hzero ↦ hx.2 (Set.mem_singleton_iff.mpr hzero))
    · intro v hv
      simp only [eval_map]
      change aeval v (blockCoordinateSignChange k P) = aeval v P
      rw [aeval_axisSign, aeval_restrictAxes, aeval_restrictAxes]
      have hv₀ : v 0 ≠ 0 := Set.mem_pi.1 hv 0 (Set.mem_univ _)
      have hv₁ : v 1 ≠ 0 := Set.mem_pi.1 hv 1 (Set.mem_univ _)
      apply hi
      · rw [aeval_axes_plane, aeval_axes_plane]
        split_ifs <;> simp
      · rw [aeval_axes_height, aeval_axes_height]
        split_ifs <;> simp
      · rw [aeval_axes_plane]
        exact mul_ne_zero ha' (by split_ifs <;> simp [hv₀])
      · rw [aeval_axes_height]
        exact mul_ne_zero hb' (by split_ifs <;> simp [hv₁])
  obtain ⟨r, hr, hd⟩ := exists_expand_two_of_sign_invariant P hsign
  refine ⟨r, ?_, hd.le.trans (restrictAxes_totalDegree i p)⟩
  let target := aeval ![C (a 0)⁻¹ * planeNorm a, C (b i)⁻¹ * heightNorm b] r
  have hA : planeNorm (M := M) a ≠ 0 := by
    intro h
    have hz := aeval_axes_plane (L := K) a i (fun _ ↦ 1)
    rw [h, map_zero] at hz
    exact ha (by simpa using hz.symm)
  have hB : heightNorm b ≠ 0 := by
    intro h
    have hz := aeval_axes_height (L := K) b i (fun _ ↦ 1)
    rw [h, map_zero] at hz
    exact hb (by simpa using hz.symm)
  have hz : planeNorm a * heightNorm b * (p - target) = 0 := by
    apply MvPolynomial.map_injective f hf
    apply MvPolynomial.funext
    intro v
    simp only [eval_map]
    change aeval v (planeNorm a * heightNorm b * (p - target)) = aeval v (0 : QuadPolynomial K M)
    rw [map_mul, map_mul, map_sub, map_zero]
    by_cases hAv : aeval v (planeNorm a) = 0
    · simp [hAv]
    by_cases hBv : aeval v (heightNorm b) = 0
    · simp [hBv]
    obtain ⟨x, hx⟩ := IsAlgClosed.exists_pow_nat_eq
      (aeval v (planeNorm a) / f (a 0)) (n := 2) (by decide)
    obtain ⟨z, hz⟩ := IsAlgClosed.exists_pow_nat_eq
      (aeval v (heightNorm b) / f (b i)) (n := 2) (by decide)
    let w : Fin 2 → L := ![x, z]
    have hA' : aeval v (planeNorm a) = aeval (axes i w) (planeNorm a) := by
      rw [aeval_axes_plane]
      change _ = f (a 0) * x ^ 2
      rw [hx, mul_div_cancel₀ _ ha']
    have hB' : aeval v (heightNorm b) = aeval (axes i w) (heightNorm b) := by
      rw [aeval_axes_height]
      change _ = f (b i) * z ^ 2
      rw [hz, mul_div_cancel₀ _ hb']
    have hpv := hi v (axes i w) hA' hB' hAv hBv
    rw [← aeval_restrictAxes] at hpv
    change aeval v p = aeval w P at hpv
    rw [← hr, expand, ← AlgHom.comp_apply, aeval_comp_bind₁] at hpv
    have htarget : aeval v target = aeval (fun j ↦ w j ^ 2) r := by
      dsimp only [target]
      rw [comp_aeval_apply]
      congr with j
      fin_cases j
      · change aeval v (C (a 0)⁻¹ * planeNorm a) = x ^ 2
        rw [map_mul, aeval_C, map_inv₀, hx, div_eq_mul_inv, mul_comm]
      · change aeval v (C (b i)⁻¹ * heightNorm b) = z ^ 2
        rw [map_mul, aeval_C, map_inv₀, hz, div_eq_mul_inv, mul_comm]
    rw [htarget] at *
    have hh : aeval v p = aeval (fun j ↦ w j ^ 2) r := by simpa using hpv
    rw [hh, sub_self, mul_zero]
  exact sub_eq_zero.mp ((mul_eq_zero.mp hz).resolve_left (mul_ne_zero hA hB))

end TwoBlocks

section BlockIsometries

open Matrix
open scoped Matrix

variable {M : ℕ}
local notation "L" => AlgebraicClosure K

private noncomputable def liftDiagonal {J : Type*} (d : J → Kˣ) : J → Lˣ :=
  Units.map (algebraMap K L).toMonoidHom ∘ d

/-- Polynomial invariance under the two orthogonal blocks, tested on
geometric points. -/
def BlockIsometryInvariant (a : Fin 2 → Kˣ) (b : Fin M → Kˣ)
    (p : QuadPolynomial K M) : Prop :=
  ∀ (T : Matrix (Fin 2) (Fin 2) L) (V : Matrix (Fin M) (Fin M) L),
    (∀ v, diagonalNorm (liftDiagonal a) (T *ᵥ v) = diagonalNorm (liftDiagonal a) v) →
    (∀ v, diagonalNorm (liftDiagonal b) (V *ᵥ v) = diagonalNorm (liftDiagonal b) v) →
    ∀ v : Fin 2 ⊕ Fin M → L,
      aeval (Sum.elim (T *ᵥ (v ∘ Sum.inl)) (V *ᵥ (v ∘ Sum.inr))) p = aeval v p

/-- Full orthogonal transitivity on the two blocks gives the polynomial
normal form and degree bound used in the draft. -/
theorem blockPolynomial_of_blockIsometryInvariant (a : Fin 2 → Kˣ) (b : Fin M → Kˣ)
    (p : QuadPolynomial K M) (hi : BlockIsometryInvariant a b p) (i : Fin M) :
    ∃ r : MvPolynomial (Fin 2) K,
      p = aeval ![C ((a 0 : K)⁻¹) * planeNorm (fun j ↦ (a j : K)),
        C ((b i : K)⁻¹) * heightNorm (fun j ↦ (b j : K))] r ∧
      2 * r.totalDegree ≤ p.totalDegree := by
  apply blockPolynomial_of_constant_norm_levels (fun j ↦ (a j : K)) (fun j ↦ (b j : K)) p _
    (a 0).ne_zero i (b i).ne_zero
  intro v w hA hB hA₀ hB₀
  have heA (f : Fin 2 ⊕ Fin M → L) :
      aeval f (planeNorm (fun j ↦ (a j : K))) =
        diagonalNorm (liftDiagonal a) (f ∘ Sum.inl) := by
    simp [planeNorm, diagonalNorm, diagonalPair, liftDiagonal, pow_two, mul_assoc]
  have heB (f : Fin 2 ⊕ Fin M → L) :
      aeval f (heightNorm (fun j ↦ (b j : K))) =
        diagonalNorm (liftDiagonal b) (f ∘ Sum.inr) := by
    simp [heightNorm, diagonalNorm, diagonalPair, liftDiagonal, pow_two, mul_assoc]
  rw [heA, heA] at hA
  rw [heB, heB] at hB
  rw [heA, hA] at hA₀
  rw [heB, hB] at hB₀
  obtain ⟨T, hT, heT⟩ := exists_orthogonal_same_norm (liftDiagonal a) _ _ hA hA₀
  obtain ⟨V, hV, heV⟩ := exists_orthogonal_same_norm (liftDiagonal b) _ _ hB hB₀
  have hh : Sum.elim (w ∘ Sum.inl) (w ∘ Sum.inr) = w := by
    funext j
    cases j <;> rfl
  rw [← hi T V hT hV v, heT, heV, hh]

end BlockIsometries

section Projection

variable {V A : Type*} [AddCommGroup V] [Module K V] [CommRing A] [Algebra K A]

/-- The augmentation span for block averaging on a stable vector space. -/
def blockDifferenceSpan {J : Type*} (action : J → V →ₗ[K] V) : Submodule K V :=
  Submodule.span K {v | ∃ j x, v = action j x - x}

/-- The elements fixed by every map in an action. -/
def actionFixed {J : Type*} (action : J → V →ₗ[K] V) : Submodule K V :=
  ⨅ j, LinearMap.ker (action j - LinearMap.id)

omit [CharZero K] in
theorem mem_actionFixed_iff {J : Type*} (action : J → V →ₗ[K] V) (v : V) :
    v ∈ actionFixed action ↔ ∀ j, action j v = v := by
  simp [actionFixed, Submodule.mem_iInf, LinearMap.mem_ker, sub_eq_zero]

/-- Stability under the action is the hypothesis for an invariant complement. -/
def ActionStable {J : Type*} (action : J → V →ₗ[K] V) (W : Submodule K V) : Prop :=
  ∀ j v, v ∈ W → action j v ∈ W

omit [CharZero K] in
/-- Complete reducibility gives the direct sum of the fixed subspace
and the augmentation span. -/
theorem isCompl_actionFixed_differenceSpan {J : Type*} (action : J → V →ₗ[K] V)
    (complete : ∀ W : Submodule K V, ActionStable action W →
      ∃ W' : Submodule K V, ActionStable action W' ∧ IsCompl W W') :
    IsCompl (actionFixed action) (blockDifferenceSpan action) := by
  let D := blockDifferenceSpan action
  let E := actionFixed action
  have hd (j : J) (v : V) : action j v - v ∈ D :=
    Submodule.subset_span ⟨j, v, rfl⟩
  have hD : ActionStable action D := by
    intro j v hv
    simpa only [sub_add_cancel] using D.add_mem (hd j v) hv
  have hE : ActionStable action E := by
    intro j v hv
    rwa [(mem_actionFixed_iff action v).1 hv j]
  obtain ⟨D', hD', hcD⟩ := complete D hD
  obtain ⟨E', hE', hcE⟩ := complete E hE
  have hDE : D' ≤ E := by
    intro v hv
    apply (mem_actionFixed_iff action v).2
    intro j
    have he := D'.sub_mem (hD' j v hv) hv
    exact sub_eq_zero.mp (Submodule.disjoint_def.mp hcD.disjoint _ (hd j v) he)
  have hED : D ≤ E' := by
    apply Submodule.span_le.2
    rintro _ ⟨j, v, rfl⟩
    have hv : v ∈ E ⊔ E' := hcE.sup_eq_top.symm ▸ Submodule.mem_top
    obtain ⟨x, hx, y, hy, rfl⟩ := Submodule.mem_sup.mp hv
    rw [map_add, (mem_actionFixed_iff action x).1 hx j, add_sub_add_left_eq_sub]
    exact E'.sub_mem (hE' j y hy) hy
  exact ⟨hcE.disjoint.mono_right hED, (hcD.codisjoint.mono_right hDE).symm⟩

omit [CharZero K] in
/-- A block-invariant functional vanishes on the augmentation span. Once
the fixed subspace complements it, one representative works for all
marked moments. -/
theorem exists_block_average_of_complement {J : Type*} (action : J → V →ₗ[K] V)
    (fixed : Submodule K V) (decomp : IsCompl fixed (blockDifferenceSpan action))
    (incl : V →ₗ[K] A) (trace : A →ₗ[K] K) (phase : A)
    (inv : ∀ j v m, 0 < m → trace (incl (action j v) * phase ^ m) = trace (incl v * phase ^ m))
    (v : V) :
    ∃ w ∈ fixed, ∀ m : ℕ, 0 < m →
      trace (incl v * phase ^ m) = trace (incl w * phase ^ m) := by
  have hv : v ∈ fixed ⊔ blockDifferenceSpan action := decomp.sup_eq_top.symm ▸ Submodule.mem_top
  obtain ⟨w, hw, u, hu, rfl⟩ := Submodule.mem_sup.mp hv
  clear hv
  refine ⟨w, hw, fun m hm ↦ ?_⟩
  have hu' : trace (incl u * phase ^ m) = 0 := by
    induction hu using Submodule.span_induction with
    | mem x hx =>
        obtain ⟨j, y, rfl⟩ := hx
        rw [map_sub, sub_mul, map_sub, inv j y m hm, sub_self]
    | zero => simp
    | add x y _ _ hx hy => rw [map_add, add_mul, map_add, hx, hy, add_zero]
    | smul c x _ hx => rw [map_smul, smul_mul_assoc, map_smul, hx, smul_zero]
  rw [map_add, add_mul, map_add, hu', add_zero]

end Projection

section BlockSpan

variable {A : Type*} [CommRing A] [Algebra K A]

/-- The span of the draft's normalized block monomials at degree at most 2J. -/
def blockMonomialSpan (U : Aˣ) (S : A) (J : ℕ) : Submodule K A :=
  Submodule.span K {v | ∃ ℓ j : ℕ, j ≤ ℓ ∧ ℓ ≤ J ∧ v = (U : A) ^ ℓ * S ^ j}

omit [CharZero K] in
private theorem blockPower_mem_span (U : Aˣ) (S : A) {J k j : ℕ} (h : k + j ≤ J) :
    ((U : A) * (1 - S)) ^ k * ((U : A) * S) ^ j ∈ blockMonomialSpan (K := K) U S J := by
  rw [mul_pow, mul_pow]
  have he : (U : A) ^ k * (1 - S) ^ k * ((U : A) ^ j * S ^ j) =
      (1 - S) ^ k * ((U : A) ^ (k + j) * S ^ j) := by
    rw [pow_add]
    ring
  rw [he, sub_eq_add_neg, add_comm 1, add_pow, Finset.sum_mul]
  apply Submodule.sum_mem
  intro r hr
  have hrk : r ≤ k := by have := Finset.mem_range.mp hr; omega
  rw [one_pow, mul_one]
  have hterm : (-S) ^ r * (k.choose r : A) * ((U : A) ^ (k + j) * S ^ j) =
      (((-1 : K) ^ r * (k.choose r : K)) • ((U : A) ^ (k + j) * S ^ (j + r))) := by
    rw [Algebra.smul_def, map_mul, map_pow, map_neg, map_one, map_natCast, neg_pow, pow_add S]
    ring
  rw [hterm]
  exact Submodule.smul_mem _ _ (Submodule.subset_span ⟨k + j, j + r, by omega, h, rfl⟩)

omit [CharZero K] in
/-- Replacing A by q − B expresses the block average in the normalized
monomials with the same degree bound. -/
theorem blockPolynomial_mem_monomialSpan (U : Aˣ) (S : A) (c₀ c₁ : K)
    (r : MvPolynomial (Fin 2) K) {J : ℕ} (hJ : r.totalDegree ≤ J) :
    aeval ![algebraMap K A c₀ * ((U : A) * (1 - S)),
      algebraMap K A c₁ * ((U : A) * S)] r ∈ blockMonomialSpan (K := K) U S J := by
  classical
  nth_rw 1 [r.as_sum]
  rw [map_sum]
  apply Submodule.sum_mem
  intro d hd
  rw [aeval_monomial, Finsupp.prod_fintype _ _ (fun _ ↦ pow_zero _), Fin.prod_univ_two]
  change algebraMap K A (r.coeff d) *
      ((algebraMap K A c₀ * ((U : A) * (1 - S))) ^ d 0 *
        (algebraMap K A c₁ * ((U : A) * S)) ^ d 1) ∈ _
  rw [mul_pow (algebraMap K A c₀), mul_pow (algebraMap K A c₁)]
  have he : d 0 + d 1 ≤ J := by
    have h := (le_totalDegree hd).trans hJ
    rwa [Finsupp.sum_fintype _ _ (fun _ ↦ rfl), Fin.sum_univ_two] at h
  rw [← map_pow, ← map_pow]
  rw [show algebraMap K A (r.coeff d) *
      (algebraMap K A (c₀ ^ d 0) * ((U : A) * (1 - S)) ^ d 0 *
        (algebraMap K A (c₁ ^ d 1) * ((U : A) * S) ^ d 1)) =
      (r.coeff d * c₀ ^ d 0 * c₁ ^ d 1) •
        (((U : A) * (1 - S)) ^ d 0 * ((U : A) * S) ^ d 1) by
      rw [Algebra.smul_def, map_mul, map_mul]
      ring]
  exact Submodule.smul_mem _ _ (blockPower_mem_span U S he)

end BlockSpan

section PolynomialAveraging

variable {M : ℕ} (a : Fin 2 → Kˣ) (b : Fin M → Kˣ)

local notation "aK" => (fun i ↦ (a i : K))
local notation "bK" => (fun i ↦ (b i : K))
local notation "R" => QuadPolynomial K M
local notation "A" => NormComplement aK bK

/-- The polynomial filtration used for averaging a fixed multiplier. -/
def polynomialDegreeAtMost (D : ℕ) : Submodule K R where
  carrier := {p | p.totalDegree ≤ D}
  zero_mem' := by simp
  add_mem' hp hq := (totalDegree_add _ _).trans (max_le hp hq)
  smul_mem' _ _ hp := (totalDegree_smul_le _ _).trans hp

/-- An invariant polynomial of degree D has an expansion in the block
monomials with norm exponent at most ⌊D/2⌋. -/
theorem blockInvariant_mem_monomialSpan (i : Fin M) (p : R)
    (hi : BlockIsometryInvariant a b p) (D : ℕ) (hp : p.totalDegree ≤ D) :
    algebraMap R A p ∈ blockMonomialSpan (K := K) (normUnit aK bK) (height aK bK) (D / 2) := by
  obtain ⟨r, hr, hd⟩ := blockPolynomial_of_blockIsometryInvariant a b p hi i
  have hbound : r.totalDegree ≤ D / 2 := by omega
  let U := normUnit aK bK
  let S := height aK bK
  let φ := IsScalarTower.toAlgHom K R A
  have hB : (U : A) * S = φ (heightNorm bK) := by
    change (U : A) * (φ (heightNorm bK) * (↑(U⁻¹) : A)) = _
    rw [mul_left_comm, Units.mul_inv, mul_one]
  have hA : (U : A) * (1 - S) = φ (planeNorm aK) := by
    rw [mul_sub, mul_one, hB]
    change (↑(normUnit aK bK) : A) - φ (heightNorm bK) = _
    rw [normUnit_val]
    change φ (quadNorm aK bK) - φ (heightNorm bK) = _
    change φ (planeNorm aK + heightNorm bK) - φ (heightNorm bK) = _
    rw [map_add, add_sub_cancel_right]
  change φ p ∈ _
  rw [hr, comp_aeval_apply]
  have he : (fun j ↦ φ ((![C ((a 0 : K)⁻¹) * planeNorm aK,
      C ((b i : K)⁻¹) * heightNorm bK] : Fin 2 → R) j)) =
      ![algebraMap K A (↑(a 0) : K)⁻¹ * ((U : A) * (1 - S)),
        algebraMap K A (↑(b i) : K)⁻¹ * ((U : A) * S)] := by
    funext j
    fin_cases j
    · change φ (C ((a 0 : K)⁻¹) * planeNorm aK) =
        algebraMap K A (↑(a 0) : K)⁻¹ * ((U : A) * (1 - S))
      rw [hA, map_mul, ← algebraMap_eq, AlgHom.commutes]
    · change φ (C ((b i : K)⁻¹) * heightNorm bK) =
        algebraMap K A (↑(b i) : K)⁻¹ * ((U : A) * S)
      rw [hB, map_mul, ← algebraMap_eq, AlgHom.commutes]
  rw [he]
  exact blockPolynomial_mem_monomialSpan U S _ _ r hbound

/-- The degree-bounded block average preserves all marked moments.
Complete reducibility and geometric compatibility are stated on the
finite polynomial representation used for the average. -/
theorem draft_block_average_polynomial {J : Type*} (i : Fin M) (D : ℕ)
    (action : J → polynomialDegreeAtMost (K := K) (M := M) D →ₗ[K]
      polynomialDegreeAtMost (K := K) (M := M) D)
    (complete : ∀ W : Submodule K (polynomialDegreeAtMost (K := K) (M := M) D),
      ActionStable action W → ∃ W', ActionStable action W' ∧ IsCompl W W')
    (fixed_geom : ∀ p ∈ actionFixed action, BlockIsometryInvariant a b (p : R))
    (phase : A) (trace : A →ₗ[K] K)
    (inv : ∀ j (p : polynomialDegreeAtMost (K := K) (M := M) D) m, 0 < m →
      trace (algebraMap R A (action j p : R) * phase ^ m) =
        trace (algebraMap R A (p : R) * phase ^ m))
    (p : polynomialDegreeAtMost (K := K) (M := M) D) :
    ∃ w ∈ blockMonomialSpan (K := K) (normUnit aK bK) (height aK bK) (D / 2),
      ∀ m : ℕ, 0 < m → trace (algebraMap R A (p : R) * phase ^ m) = trace (w * phase ^ m) := by
  let incl := (IsScalarTower.toAlgHom K R A).toLinearMap.comp
    (polynomialDegreeAtMost (K := K) (M := M) D).subtype
  obtain ⟨w, hw, havg⟩ := exists_block_average_of_complement action (actionFixed action)
    (isCompl_actionFixed_differenceSpan action complete) incl trace phase inv p
  exact ⟨algebraMap R A (w : R), blockInvariant_mem_monomialSpan a b i _ (fixed_geom w hw) D w.2, havg⟩

end PolynomialAveraging

end MathieuQuadratic
