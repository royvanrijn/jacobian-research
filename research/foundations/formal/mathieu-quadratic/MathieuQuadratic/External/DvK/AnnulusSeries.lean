/-
Adapted from GMC2/AnnulusSeries.lean at MurrellGroup/GMC-2-lean commit
1782de7ff6c97eb1d98e63e7ff34df18b9cd322e. MIT license in this directory.
-/
import MathieuQuadratic.External.DvK.Compat
import Mathlib.Topology.Algebra.InfiniteSum.Nonarchimedean
import Mathlib.Topology.Algebra.InfiniteSum.NatInt
import Mathlib.Topology.Algebra.Valued.WithZeroMulInt
import Mathlib.Algebra.Polynomial.Laurent

open scoped WithZero
open scoped Topology

set_option maxHeartbeats 1500000

namespace MathieuQuadratic.DvKProof

local instance coeLaurentAnnulusSeries {R : Type*} [Semiring R] : CoeFun (LaurentPolynomial R) (fun _ ↦ ℤ → R) :=
  ⟨fun p ↦ p.coeff⟩


open Filter

noncomputable section

/-- The scalar submodule of summable integer-indexed coefficient families. -/
def annulusSeriesSubmodule
    (K : Type*) [Field K] [TopologicalSpace K]
    [IsTopologicalRing K] :
    Submodule K (ℤ → K) where
  carrier := Summable
  zero_mem' := summable_zero
  add_mem' ha hb := ha.add hb
  smul_mem' c a ha := by
    change Summable (fun i ↦ c * a i)
    exact ha.mul_left c

/--
A bilateral Laurent series whose coefficients tend to zero away from every
finite set. Over a complete nonarchimedean field this is equivalently a
summable family of coefficients.
-/
abbrev AnnulusSeries
    (K : Type*) [Field K] [TopologicalSpace K]
    [IsTopologicalRing K] :=
  annulusSeriesSubmodule K

namespace AnnulusSeries

variable {K : Type*} [Field K] [Valued K ℤᵐ⁰]
  [NonarchimedeanRing K] [CompleteSpace K]

instance : CoeFun (AnnulusSeries K) fun _ ↦ ℤ → K :=
  ⟨fun a ↦ a.1⟩

omit [CompleteSpace K] in
@[ext]
theorem ext {a b : AnnulusSeries K} (h : ∀ n, a n = b n) : a = b :=
  Subtype.ext (funext h)

/-- The coefficient of a bilateral series at an integer exponent. -/
def coeff (n : ℤ) : AnnulusSeries K →+ K where
  toFun a := a n
  map_zero' := rfl
  map_add' _ _ := rfl

omit [CompleteSpace K] in
@[simp]
theorem coeff_zero (n : ℤ) : (0 : AnnulusSeries K) n = 0 :=
  rfl

omit [CompleteSpace K] in
@[simp]
theorem coeff_add (a b : AnnulusSeries K) (n : ℤ) :
    (a + b) n = a n + b n :=
  rfl

omit [CompleteSpace K] in
@[simp]
theorem coeff_neg (a : AnnulusSeries K) (n : ℤ) :
    (-a) n = -a n :=
  rfl

omit [CompleteSpace K] in
@[simp]
theorem coeff_sub (a b : AnnulusSeries K) (n : ℤ) :
    (a - b) n = a n - b n := by
  rfl

/-- The bilateral monomial supported at one integer exponent. -/
def single (n : ℤ) (c : K) : AnnulusSeries K :=
  ⟨Pi.single n c, summable_of_hasFiniteSupport <|
    (Set.finite_singleton n).subset (by
      intro m hm
      simp only [Set.mem_singleton_iff]
      by_contra hmn
      exact hm (by simp [hmn]))⟩

omit [CompleteSpace K] in
@[simp]
theorem single_apply (n m : ℤ) (c : K) :
    single n c m = if n = m then c else 0 := by
  simp [single, Pi.single_apply, eq_comm]

/--
The coefficient sequence of the annulus expansion of `(x-a)⁻¹`. When
`v(a) < 1` it is supported at negative exponents; otherwise it is supported
at nonnegative exponents.
-/
def linearInverseCoeff (a : K) : ℤ → K
  | .ofNat n =>
      if 1 < Valued.v a then -(a⁻¹) ^ (n + 1) else 0
  | .negSucc n =>
      if Valued.v a < 1 then a ^ n else 0

omit [NonarchimedeanRing K] [CompleteSpace K] in
@[simp]
theorem linearInverseCoeff_ofNat (a : K) (n : ℕ) :
    linearInverseCoeff a (n : ℤ) =
      if 1 < Valued.v a then -(a⁻¹) ^ (n + 1) else 0 :=
  rfl

omit [NonarchimedeanRing K] [CompleteSpace K] in
@[simp]
theorem linearInverseCoeff_negSucc (a : K) (n : ℕ) :
    linearInverseCoeff a (.negSucc n) =
      if Valued.v a < 1 then a ^ n else 0 :=
  rfl

theorem summable_pow_of_valuation_lt_one
    (a : K) (ha : Valued.v a < 1) :
    Summable fun n : ℕ ↦ a ^ n := by
  rw [NonarchimedeanAddGroup.summable_iff_tendsto_cofinite_zero,
    Nat.cofinite_eq_atTop]
  exact Valued.tendsto_zero_pow_of_v_lt_one ha

theorem summable_linearInverseCoeff
    (a : K) (ha : Valued.v a ≠ 1) :
    Summable (linearInverseCoeff a) := by
  rcases lt_or_gt_of_ne ha with haSmall | haLarge
  · apply Summable.of_nat_of_neg_add_one
    · simp [linearInverseCoeff, not_lt.mpr haSmall.le]
    · have hp := summable_pow_of_valuation_lt_one a haSmall
      refine hp.congr fun n ↦ ?_
      rw [show -((n : ℤ) + 1) = Int.negSucc n by omega]
      simp [linearInverseCoeff, haSmall]
  · have ha0 : a ≠ 0 := by
      intro hzero
      simp [hzero] at haLarge
    have haInv : Valued.v a⁻¹ < 1 := by
      rw [map_inv₀]
      exact (inv_lt_one₀ (lt_trans (by simp) haLarge)).2 haLarge
    apply Summable.of_nat_of_neg_add_one
    · simpa [linearInverseCoeff, haLarge, pow_succ] using
        (summable_pow_of_valuation_lt_one a⁻¹ haInv).mul_right (-a⁻¹)
    · have hz : Summable fun _ : ℕ ↦ (0 : K) := summable_zero
      refine hz.congr fun n ↦ ?_
      rw [show -((n : ℤ) + 1) = Int.negSucc n by omega]
      simp [linearInverseCoeff, not_lt.mpr haLarge.le]

/-- The annulus expansion of `(x-a)⁻¹`. -/
def linearInverse (a : K) (ha : Valued.v a ≠ 1) :
    AnnulusSeries K :=
  ⟨linearInverseCoeff a, summable_linearInverseCoeff a ha⟩

@[simp]
theorem linearInverse_ofNat
    (a : K) (ha : Valued.v a ≠ 1) (n : ℕ) :
    linearInverse a ha (n : ℤ) =
      if 1 < Valued.v a then -(a⁻¹) ^ (n + 1) else 0 :=
  rfl

@[simp]
theorem linearInverse_negSucc
    (a : K) (ha : Valued.v a ≠ 1) (n : ℕ) :
    linearInverse a ha (.negSucc n) =
      if Valued.v a < 1 then a ^ n else 0 :=
  rfl

/-- The constant coefficient in the small-root expansion is zero. -/
theorem linearInverse_coeff_zero_of_lt
    (a : K) (ha : Valued.v a < 1) :
    linearInverse a (ne_of_lt ha) 0 = 0 := by
  rw [show (0 : ℤ) = (0 : ℕ) by rfl, linearInverse_ofNat]
  simp [not_lt.mpr ha.le]

/-- The constant coefficient in the large-root expansion is `-a⁻¹`. -/
theorem linearInverse_coeff_zero_of_gt
    (a : K) (ha : 1 < Valued.v a) :
    linearInverse a (ne_of_gt ha) 0 = -a⁻¹ := by
  rw [show (0 : ℤ) = (0 : ℕ) by rfl, linearInverse_ofNat]
  simp [ha]

/-- Coefficientwise scalar multiplication of an annulus series. -/
def scale (c : K) (u : AnnulusSeries K) : AnnulusSeries K :=
  ⟨fun n ↦ c * u n, u.2.mul_left c⟩

omit [CompleteSpace K] in
@[simp]
theorem scale_apply (c : K) (u : AnnulusSeries K) (n : ℤ) :
    scale c u n = c * u n :=
  rfl

/-- A finite linear combination of geometric linear-factor inverses. -/
def linearInverseSum
    (roots : Finset K) (c : K → K)
    (hunit : ∀ a ∈ roots, Valued.v a ≠ 1) :
    AnnulusSeries K :=
  ∑ a ∈ roots.attach,
    scale (c a.1) (linearInverse a.1 (hunit a.1 a.2))

/--
The constant coefficient of a finite partial-fraction expansion is the sum
over exactly those roots whose valuation is greater than one.
-/
theorem linearInverseSum_coeff_zero
    (roots : Finset K) (c : K → K)
    (hunit : ∀ a ∈ roots, Valued.v a ≠ 1) :
    linearInverseSum roots c hunit 0 =
      ∑ a ∈ roots.filter fun a ↦ 1 < Valued.v a,
        c a * (-a⁻¹) := by
  classical
  rw [linearInverseSum]
  change
    coeff 0
        (∑ a ∈ roots.attach,
          scale (c a.1) (linearInverse a.1 (hunit a.1 a.2))) =
      _
  rw [map_sum, Finset.sum_filter,
    ← roots.sum_attach (fun a ↦
      if 1 < Valued.v a then c a * -a⁻¹ else 0)]
  apply Finset.sum_congr rfl
  intro a ha
  by_cases hlarge : 1 < Valued.v a.1
  · simp [coeff, hlarge, linearInverse_coeff_zero_of_gt]
  · have hsmall : Valued.v a.1 < 1 :=
      lt_of_le_of_ne (le_of_not_gt hlarge) (hunit a.1 a.2)
    change
      c a.1 * linearInverse a.1 (hunit a.1 a.2) 0 =
        (if 1 < Valued.v a.1 then c a.1 * -a.1⁻¹ else 0)
    rw [if_neg hlarge, linearInverse_coeff_zero_of_lt a.1 hsmall,
      mul_zero]

/-- A finite linear-factor expansion with roots indexed by another type. -/
def indexedLinearInverseSum
    {ι : Type*} (indices : Finset ι) (root c : ι → K)
    (hunit : ∀ i ∈ indices, Valued.v (root i) ≠ 1) :
    AnnulusSeries K :=
  ∑ i ∈ indices.attach,
    scale (c i.1) (linearInverse (root i.1) (hunit i.1 i.2))

/--
Constant-coefficient extraction for an indexed finite partial-fraction
expansion.
-/
theorem indexedLinearInverseSum_coeff_zero
    {ι : Type*} (indices : Finset ι) (root c : ι → K)
    (hunit : ∀ i ∈ indices, Valued.v (root i) ≠ 1) :
    indexedLinearInverseSum indices root c hunit 0 =
      ∑ i ∈ indices.filter fun i ↦ 1 < Valued.v (root i),
        c i * (-(root i)⁻¹) := by
  classical
  rw [indexedLinearInverseSum]
  change
    coeff 0
        (∑ i ∈ indices.attach,
          scale (c i.1)
            (linearInverse (root i.1) (hunit i.1 i.2))) =
      _
  rw [map_sum, Finset.sum_filter,
    ← indices.sum_attach (fun i ↦
      if 1 < Valued.v (root i) then c i * -(root i)⁻¹ else 0)]
  apply Finset.sum_congr rfl
  intro i hi
  by_cases hlarge : 1 < Valued.v (root i.1)
  · simp [coeff, hlarge, linearInverse_coeff_zero_of_gt]
  · have hsmall : Valued.v (root i.1) < 1 :=
      lt_of_le_of_ne (le_of_not_gt hlarge) (hunit i.1 i.2)
    change
      c i.1 * linearInverse (root i.1) (hunit i.1 i.2) 0 =
        (if 1 < Valued.v (root i.1) then
          c i.1 * -(root i.1)⁻¹ else 0)
    rw [if_neg hlarge,
      linearInverse_coeff_zero_of_lt (root i.1) hsmall, mul_zero]

/-- Translation of exponents by an integer. -/
def shift (k : ℤ) : AnnulusSeries K →ₗ[K] AnnulusSeries K where
  toFun u :=
    ⟨fun n ↦ u (n - k), by
      change Summable (fun n ↦ u (n - k))
      exact (((Equiv.addRight (-k)).summable_iff).mpr u.2).congr
        (fun n ↦ by simp [sub_eq_add_neg])⟩
  map_add' _ _ := rfl
  map_smul' _ _ := rfl

omit [CompleteSpace K] in
@[simp]
theorem shift_apply (k n : ℤ) (u : AnnulusSeries K) :
    shift k u n = u (n - k) :=
  rfl

omit [CompleteSpace K] in
theorem shift_injective (k : ℤ) :
    Function.Injective (shift k : AnnulusSeries K → AnnulusSeries K) := by
  intro u v huv
  ext n
  have h := congrArg (fun w : AnnulusSeries K ↦ w (n + k)) huv
  simpa [shift_apply] using h

omit [CompleteSpace K] in
theorem shift_scale
    (k : ℤ) (c : K) (u : AnnulusSeries K) :
    shift k (scale c u) = scale c (shift k u) := by
  ext n
  rfl

/-- Integer shifts as a multiplicative family of linear endomorphisms. -/
def shiftMonoidHom :
    Multiplicative ℤ →* Module.End K (AnnulusSeries K) where
  toFun k := shift k.toAdd
  map_one' := by
    ext u n
    simp [shift]
  map_mul' k l := by
    ext u n
    simp [shift, sub_sub]

/--
The representation of Laurent polynomials by finite convolution on annulus
series.
-/
def laurentActionHom :
    LaurentPolynomial K →+* Module.End K (AnnulusSeries K) :=
  AddMonoidAlgebra.liftNCRingHom
    (algebraMap K (Module.End K (AnnulusSeries K)))
    shiftMonoidHom
    (by
      intro c k
      apply LinearMap.ext
      intro u
      rfl)

/-- Apply a Laurent polynomial to an annulus series by finite convolution. -/
def laurentAction
    (q : LaurentPolynomial K) (u : AnnulusSeries K) :
    AnnulusSeries K :=
  laurentActionHom q u

omit [CompleteSpace K] in
theorem laurentAction_mul
    (q r : LaurentPolynomial K) (u : AnnulusSeries K) :
    laurentAction (q * r) u =
      laurentAction q (laurentAction r u) := by
  rw [laurentAction, map_mul]
  rfl

omit [CompleteSpace K] in
theorem laurentAction_add
    (q : LaurentPolynomial K) (u v : AnnulusSeries K) :
    laurentAction q (u + v) =
      laurentAction q u + laurentAction q v :=
  (laurentActionHom q).map_add u v

omit [CompleteSpace K] in
theorem laurentAction_one (u : AnnulusSeries K) :
    laurentAction 1 u = u := by
  change (laurentActionHom 1) u = u
  rw [map_one]
  rfl

omit [CompleteSpace K] in
theorem laurentAction_sub_left
    (q r : LaurentPolynomial K) (u : AnnulusSeries K) :
    laurentAction (q - r) u =
      laurentAction q u - laurentAction r u := by
  change (laurentActionHom (q - r)) u = _
  rw [map_sub, LinearMap.sub_apply]
  rfl

omit [CompleteSpace K] in
theorem laurentAction_smul
    (q : LaurentPolynomial K) (c : K) (u : AnnulusSeries K) :
    laurentAction q (c • u) =
      c • laurentAction q u :=
  (laurentActionHom q).map_smul c u

omit [CompleteSpace K] in
theorem laurentAction_T
    (k : ℤ) (u : AnnulusSeries K) :
    laurentAction (LaurentPolynomial.T k) u = shift k u := by
  rw [laurentAction, laurentActionHom, LaurentPolynomial.T,
    AddMonoidAlgebra.liftNCRingHom_single]
  rw [map_one, one_mul]
  rfl

omit [CompleteSpace K] in
theorem laurentAction_C
    (c : K) (u : AnnulusSeries K) :
    laurentAction (LaurentPolynomial.C c) u = c • u := by
  change laurentAction (AddMonoidAlgebra.single (0 : ℤ) c) u = c • u
  rw [laurentAction, laurentActionHom,
    AddMonoidAlgebra.liftNCRingHom_single]
  ext n
  change c * u (n - 0) = c * u n
  rw [sub_zero]

omit [CompleteSpace K] in
theorem laurentAction_apply
    (q : LaurentPolynomial K) (u : AnnulusSeries K) (n : ℤ) :
    laurentAction q u n =
      ∑ k ∈ q.support, q k * u (n - k) := by
  classical
  change
    (((Finsupp.liftAddHom fun k : ℤ ↦
        (AddMonoidHom.mulRight
          (shift k : Module.End K (AnnulusSeries K))).comp
            (algebraMap K (Module.End K (AnnulusSeries K))).toAddMonoidHom)
      q.coeff) u) n =
      ∑ k ∈ q.support, q k * u (n - k)
  rw [Finsupp.liftAddHom_apply]
  change
    ((∑ k ∈ q.support, (q k) • shift k) u) n =
      ∑ k ∈ q.support, q k * u (n - k)
  simp [shift]

omit [CompleteSpace K] in
theorem laurentAction_single_zero_apply
    (q : LaurentPolynomial K) (c : K) (n : ℤ) :
    laurentAction q (single 0 c) n = q n * c := by
  classical
  rw [laurentAction_apply]
  by_cases hn : n ∈ q.support
  · rw [Finset.sum_eq_single n]
    · simp [single_apply]
    · intro k hk hkn
      have hzero : (0 : ℤ) ≠ n - k := by omega
      simp [single_apply, hzero]
    · intro hn'
      exact (hn' hn).elim
  · rw [Finset.sum_eq_zero]
    · rw [Finsupp.notMem_support_iff.mp hn, zero_mul]
    · intro k hk
      have hkn : k ≠ n := by
        intro h
        exact hn (h ▸ hk)
      have hzero : (0 : ℤ) ≠ n - k := by omega
      simp [single_apply, hzero]

/--
A nonzero annulus series has a coefficient of maximal valuation.
-/
theorem exists_valuation_max
    (u : AnnulusSeries K) (hu : u ≠ 0) :
    ∃ n : ℤ, u n ≠ 0 ∧ ∀ m : ℤ, Valued.v (u m) ≤ Valued.v (u n) := by
  classical
  have hnonzero : ∃ n : ℤ, u n ≠ 0 := by
    by_contra h
    push Not at h
    apply hu
    ext n
    exact h n
  obtain ⟨n₀, hn₀⟩ := hnonzero
  let v : Valuation K ℤᵐ⁰ := Valued.v
  have hrestrict₀ : v.restrict (u n₀) ≠ 0 := by
    intro hzero
    apply (v.ne_zero_iff.mpr hn₀)
    exact (Valuation.restrict_eq_zero_iff v).mp hzero
  let γ : (MonoidWithZeroHom.ValueGroup₀ (.ofClass v))ˣ :=
    Units.mk0 (v.restrict (u n₀)) hrestrict₀
  have htend :
      Tendsto (fun n : ℤ ↦ u n) cofinite (𝓝 0) :=
    (NonarchimedeanAddGroup.summable_iff_tendsto_cofinite_zero _).mp u.2
  rw [(Valued.hasBasis_nhds_zero K ℤᵐ⁰).tendsto_right_iff] at htend
  have hevent :
      ∀ᶠ n : ℤ in cofinite, v.restrict (u n) < γ.1 :=
    htend γ trivial
  let exceptional : Finset ℤ :=
    (mem_cofinite.mp hevent).toFinset
  let candidates : Finset ℤ := insert n₀ exceptional
  have hcandidates : candidates.Nonempty :=
    ⟨n₀, Finset.mem_insert_self n₀ exceptional⟩
  obtain ⟨n, hn, hmax⟩ :=
    candidates.exists_max_image (fun m ↦ Valued.v (u m)) hcandidates
  have hn₀le : Valued.v (u n₀) ≤ Valued.v (u n) :=
    hmax n₀ (Finset.mem_insert_self n₀ exceptional)
  have hn0 : u n ≠ 0 := by
    exact Valued.v.pos_iff.mp
      ((Valued.v.pos_iff.mpr hn₀).trans_le hn₀le)
  refine ⟨n, hn0, fun m ↦ ?_⟩
  by_cases hm : m ∈ candidates
  · exact hmax m hm
  · have hmExceptional : m ∉ exceptional := by
      intro hm'
      exact hm (Finset.mem_insert_of_mem hm')
    have hmEvent : v.restrict (u m) < γ.1 := by
      have hmSet :
          m ∈ {j | v.restrict (u j) < γ.1} := by
        have hmNotCompl :
          m ∉ {j | v.restrict (u j) < γ.1}ᶜ := by
          simpa [exceptional] using hmExceptional
        exact (Set.notMem_compl_iff).mp hmNotCompl
      exact hmSet
    have hmLt : Valued.v (u m) < Valued.v (u n₀) := by
      have hmLt' :=
        (Valuation.restrict_lt_iff_lt_embedding v).mp hmEvent
      change v (u m) < v (u n₀)
      rw [← Valuation.embedding_restrict v (u n₀)]
      exact hmLt'
    exact hmLt.le.trans hn₀le

/--
If every coefficient of `q` has valuation less than one, multiplication by
`1-q` is injective on annulus series.
-/
theorem laurentAction_one_sub_injective
    (q : LaurentPolynomial K)
    (hq : ∀ k ∈ q.support, Valued.v (q k) < 1) :
    Function.Injective (laurentAction (1 - q)) := by
  intro u w huw
  rw [← sub_eq_zero]
  let d : AnnulusSeries K := u - w
  have hzero :
      laurentAction (1 - q) d = 0 := by
    calc
      laurentAction (1 - q) d =
          laurentAction (1 - q) u -
            laurentAction (1 - q) w := by
        change
          (laurentActionHom (1 - q)) (u - w) =
            (laurentActionHom (1 - q)) u -
              (laurentActionHom (1 - q)) w
        exact (laurentActionHom (1 - q)).map_sub u w
      _ = 0 := sub_eq_zero.mpr huw
  have hfixed : d = laurentAction q d := by
    rw [laurentAction_sub_left, laurentAction_one] at hzero
    exact sub_eq_zero.mp hzero
  by_contra hd
  obtain ⟨n, hn0, hmax⟩ := exists_valuation_max d hd
  have hnpos : 0 < Valued.v (d n) :=
    Valued.v.pos_iff.mpr hn0
  have hcoeff :
      d n = ∑ k ∈ q.support, q k * d (n - k) := by
    calc
      d n = laurentAction q d n :=
        congrArg (fun x : AnnulusSeries K ↦ x n) hfixed
      _ = _ := laurentAction_apply q d n
  have hterm :
      ∀ k ∈ q.support,
        Valued.v (q k * d (n - k)) < Valued.v (d n) := by
    intro k hk
    rw [map_mul]
    by_cases hdk : d (n - k) = 0
    · simp [hdk, hnpos]
    · exact
        (mul_lt_of_lt_one_left
          (Valued.v.pos_iff.mpr hdk) (hq k hk)).trans_le
            (hmax (n - k))
  have hsum :
      Valued.v
          (∑ k ∈ q.support, q k * d (n - k)) <
        Valued.v (d n) :=
    Valued.v.map_sum_lt hnpos.ne' hterm
  rw [← hcoeff] at hsum
  exact (lt_irrefl _) hsum

omit [CompleteSpace K] in
/--
A coefficientwise valuation bound for finite Laurent convolution.
-/
theorem valuation_laurentAction_le
    (q : LaurentPolynomial K) (u : AnnulusSeries K)
    (ρ M : ℤᵐ⁰)
    (hq : ∀ k ∈ q.support, Valued.v (q k) ≤ ρ)
    (hu : ∀ n : ℤ, Valued.v (u n) ≤ M)
    (n : ℤ) :
    Valued.v (laurentAction q u n) ≤ ρ * M := by
  rw [laurentAction_apply]
  apply Valued.v.map_sum_le
  intro k hk
  rw [map_mul]
  exact mul_le_mul' (hq k hk) (hu (n - k))

omit [CompleteSpace K] in
/--
Iterating a Laurent convolution multiplies its uniform valuation bound.
-/
theorem valuation_laurentAction_pow_le
    (q : LaurentPolynomial K) (u : AnnulusSeries K)
    (ρ M : ℤᵐ⁰)
    (hq : ∀ k ∈ q.support, Valued.v (q k) ≤ ρ)
    (hu : ∀ n : ℤ, Valued.v (u n) ≤ M) :
    ∀ (m : ℕ) (n : ℤ),
      Valued.v (laurentAction (q ^ m) u n) ≤ ρ ^ m * M := by
  intro m
  induction m with
  | zero =>
      intro n
      simpa [laurentAction_one] using hu n
  | succ m ih =>
      intro n
      rw [pow_succ']
      rw [laurentAction_mul]
      have hstep :=
        valuation_laurentAction_le q
          (laurentAction (q ^ m) u) ρ (ρ ^ m * M)
          hq ih
          n
      simpa [pow_succ', mul_assoc] using hstep

/--
Finite geometric iteration extracts the constant coefficient of a strict
Laurent-contraction fixed point.
-/
theorem fixedPoint_coeff_zero
    (q : LaurentPolynomial K) (u : AnnulusSeries K)
    (z : K) (ρ : ℤᵐ⁰)
    (hfixed :
      u = scale z (single 0 1) + laurentAction q u)
    (hρ : ρ < 1)
    (hq : ∀ k ∈ q.support, Valued.v (q k) ≤ ρ)
    (hpowzero : ∀ m : ℕ, 1 ≤ m → (q ^ m) 0 = 0) :
    u 0 = z := by
  by_cases hu0 : u = 0
  · have hcoeff :=
      congrArg (fun w : AnnulusSeries K ↦ w 0) hfixed
    simp [hu0, laurentAction, scale_apply, single_apply] at hcoeff
    simpa [hu0] using hcoeff
  have hiterate :
      ∀ t : ℕ,
        u 0 =
          z + laurentAction (q ^ (t + 1)) u 0 := by
    intro t
    induction t with
    | zero =>
        have hcoeff :=
          congrArg (fun w : AnnulusSeries K ↦ w 0) hfixed
        simpa [pow_one, scale_apply, single_apply] using hcoeff
    | succ t ih =>
        let N := t + 1
        have hNpos : 1 ≤ N := by omega
        have hvanish :
            laurentAction (q ^ N)
                (scale z (single 0 1)) 0 = 0 := by
          change
            laurentAction (q ^ N)
                (z • single 0 1) 0 = 0
          rw [laurentAction_smul]
          change
            z * laurentAction (q ^ N) (single 0 1) 0 = 0
          rw [laurentAction_single_zero_apply,
            hpowzero N hNpos, zero_mul, mul_zero]
        have hremainder :
            laurentAction (q ^ N) u 0 =
              laurentAction (q ^ N)
                  (scale z (single 0 1)) 0 +
                laurentAction (q ^ (N + 1)) u 0 := by
          calc
            laurentAction (q ^ N) u 0 =
                laurentAction (q ^ N)
                  (scale z (single 0 1) +
                    laurentAction q u) 0 := by rw [← hfixed]
            _ = laurentAction (q ^ N)
                    (scale z (single 0 1)) 0 +
                  laurentAction (q ^ N)
                    (laurentAction q u) 0 := by
              rw [laurentAction_add]
              rfl
            _ = laurentAction (q ^ N)
                    (scale z (single 0 1)) 0 +
                  laurentAction (q ^ (N + 1)) u 0 := by
              rw [← laurentAction_mul, ← pow_succ]
        change
          u 0 = z + laurentAction (q ^ (N + 1)) u 0
        rw [ih, hremainder, hvanish, zero_add]
  by_contra hne
  have hdelta0 : u 0 - z ≠ 0 :=
    sub_ne_zero.mpr hne
  obtain ⟨n, hn0, hmax⟩ :=
    exists_valuation_max u hu0
  let D : ℤᵐ⁰ := Valued.v (u 0 - z)
  let M : ℤᵐ⁰ := Valued.v (u n)
  have hDpos : 0 < D := by
    exact Valued.v.pos_iff.mpr hdelta0
  have hMpos : 0 < M := by
    exact Valued.v.pos_iff.mpr hn0
  have hratio0 : D * M⁻¹ ≠ 0 :=
    mul_ne_zero hDpos.ne' (inv_ne_zero hMpos.ne')
  let γ : ℤᵐ⁰ˣ := Units.mk0 (D * M⁻¹) hratio0
  obtain ⟨N, hN⟩ := exists_pow_lt₀ hρ γ
  have hpow :
      ρ ^ (N + 1) < D * M⁻¹ := by
    exact
      (pow_le_pow_right_of_le_one' hρ.le
        (Nat.le_succ N)).trans_lt hN
  have hsmall : ρ ^ (N + 1) * M < D := by
    have hmul :=
      mul_lt_mul_of_pos_right hpow hMpos
    simpa [mul_assoc, hMpos.ne'] using hmul
  have hdelta :
      u 0 - z =
        laurentAction (q ^ (N + 1)) u 0 := by
    rw [hiterate N]
    ring
  have hbound :=
    valuation_laurentAction_pow_le
      q u ρ M hq (by
        intro m
        exact hmax m) (N + 1) 0
  have hle :
      D ≤ ρ ^ (N + 1) * M := by
    change
      Valued.v (u 0 - z) ≤ ρ ^ (N + 1) * M
    rw [hdelta]
    exact hbound
  exact (not_lt_of_ge hle) hsmall

omit [CompleteSpace K] in
theorem shift_single_zero
    (k : ℤ) (c : K) :
    shift k (single 0 c) = single k c := by
  ext n
  simp only [shift_apply, single_apply]
  by_cases hkn : k = n
  · subst n
    simp
  · have hzero : (0 : ℤ) ≠ n - k := by omega
    rw [if_neg hzero, if_neg hkn]

omit [CompleteSpace K] in
theorem laurentAction_T_single_zero
    (k : ℤ) (c : K) :
    laurentAction (LaurentPolynomial.T k) (single 0 c) =
      single k c := by
  rw [laurentAction_T, shift_single_zero]

/-- Coefficientwise multiplication of a sequence by `x-a`. -/
def linearFactorCoeff (a : K) (u : ℤ → K) : ℤ → K
  | .ofNat 0 => u (.negSucc 0) - a * u 0
  | .ofNat (n + 1) => u n - a * u (n + 1)
  | .negSucc n => u (.negSucc (n + 1)) - a * u (.negSucc n)

omit [Valued K ℤᵐ⁰] [NonarchimedeanRing K] [CompleteSpace K] in
theorem linearFactorCoeff_eq
    (a : K) (u : ℤ → K) (n : ℤ) :
    linearFactorCoeff a u n = u (n - 1) - a * u n := by
  cases n with
  | ofNat n =>
      cases n <;> simp [linearFactorCoeff]
  | negSucc n =>
      simp [linearFactorCoeff]

/-- Coefficientwise multiplication of a bilateral series by `x-a`. -/
def linearFactorApply (a : K) (u : AnnulusSeries K) :
    AnnulusSeries K :=
  ⟨linearFactorCoeff a u, by
    have hshift : Summable fun n : ℤ ↦ u (n - 1) := by
      exact (((Equiv.addRight (-1 : ℤ)).summable_iff).mpr u.2).congr
        (fun n ↦ by simp [sub_eq_add_neg])
    have hmain : Summable fun n : ℤ ↦ u (n - 1) - a * u n :=
      hshift.sub (u.2.mul_left a)
    exact hmain.congr (fun n ↦ (linearFactorCoeff_eq a u n).symm)⟩

omit [CompleteSpace K] in
@[simp]
theorem linearFactorApply_apply
    (a : K) (u : AnnulusSeries K) (n : ℤ) :
    linearFactorApply a u n = u (n - 1) - a * u n :=
  linearFactorCoeff_eq a u n

omit [CompleteSpace K] in
/--
The Laurent-polynomial action of `T-a` is the coefficient recurrence
implemented by `linearFactorApply`.
-/
theorem laurentAction_T_sub_C
    (a : K) (u : AnnulusSeries K) :
    laurentAction
        (LaurentPolynomial.T 1 - LaurentPolynomial.C a) u =
      linearFactorApply a u := by
  change
    (laurentActionHom
      (LaurentPolynomial.T 1 - LaurentPolynomial.C a)) u =
      linearFactorApply a u
  rw [map_sub, LinearMap.sub_apply]
  change
    laurentAction (LaurentPolynomial.T 1) u -
        laurentAction (LaurentPolynomial.C a) u =
      linearFactorApply a u
  rw [laurentAction_T, laurentAction_C]
  ext n
  simp [linearFactorApply_apply, shift_apply]

/--
Monsky's geometric annulus expansion is a coefficientwise inverse of the
linear factor `x-a`.
-/
theorem linearFactorApply_linearInverse
    (a : K) (ha : Valued.v a ≠ 1) :
    linearFactorApply a (linearInverse a ha) = single 0 1 := by
  ext n
  rcases lt_or_gt_of_ne ha with haSmall | haLarge
  · cases n with
    | ofNat n =>
        cases n with
        | zero =>
            simp [linearFactorApply, linearFactorCoeff, linearInverse, linearInverseCoeff,
              haSmall, not_lt.mpr haSmall.le]
        | succ n =>
            change
              linearFactorCoeff a (linearInverseCoeff a) (.ofNat (n + 1)) =
                if (0 : ℤ) = .ofNat (n + 1) then (1 : K) else 0
            change
              linearInverseCoeff a (n : ℤ) -
                    a * linearInverseCoeff a ((n + 1 : ℕ) : ℤ) =
                if (0 : ℤ) = ((n + 1 : ℕ) : ℤ) then (1 : K) else 0
            rw [linearInverseCoeff_ofNat, linearInverseCoeff_ofNat]
            simp [not_lt.mpr haSmall.le]
            omega
    | negSucc n =>
        simp [linearFactorApply, linearFactorCoeff, linearInverse, linearInverseCoeff,
          haSmall, pow_succ']
  · have ha0 : a ≠ 0 := by
      intro hzero
      simp [hzero] at haLarge
    cases n with
    | ofNat n =>
        cases n with
        | zero =>
            simp [linearFactorApply, linearFactorCoeff, linearInverse, linearInverseCoeff,
              haLarge, not_lt.mpr haLarge.le, ha0]
        | succ n =>
            change
              linearFactorCoeff a (linearInverseCoeff a) (.ofNat (n + 1)) =
                if (0 : ℤ) = .ofNat (n + 1) then (1 : K) else 0
            change
              linearInverseCoeff a (n : ℤ) -
                    a * linearInverseCoeff a ((n + 1 : ℕ) : ℤ) =
                if (0 : ℤ) = ((n + 1 : ℕ) : ℤ) then (1 : K) else 0
            rw [linearInverseCoeff_ofNat, linearInverseCoeff_ofNat]
            simp only [haLarge, ↓reduceIte]
            have hn : (0 : ℤ) ≠ ((n + 1 : ℕ) : ℤ) := by omega
            rw [if_neg hn]
            change
              -(a⁻¹) ^ (n + 1) -
                  a * (-(a⁻¹) ^ (n + 2)) = 0
            calc
              _ = -(a⁻¹) ^ (n + 1) +
                    (a⁻¹) ^ (n + 1) * (a * a⁻¹) := by
                      rw [show n + 2 = (n + 1) + 1 by omega, pow_succ]
                      ring
              _ = 0 := by rw [mul_inv_cancel₀ ha0]; ring
    | negSucc n =>
        simp [linearFactorApply, linearFactorCoeff, linearInverse, linearInverseCoeff,
          not_lt.mpr haLarge.le]

theorem laurentAction_T_sub_C_linearInverse
    (a : K) (ha : Valued.v a ≠ 1) :
    laurentAction
        (LaurentPolynomial.T 1 - LaurentPolynomial.C a)
        (linearInverse a ha) =
      single 0 1 := by
  rw [laurentAction_T_sub_C, linearFactorApply_linearInverse]

end AnnulusSeries

end

end MathieuQuadratic.DvKProof
