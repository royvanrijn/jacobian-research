/-
Adapted from GMC2/OneVariableDvK.lean at MurrellGroup/GMC-2-lean commit
1782de7ff6c97eb1d98e63e7ff34df18b9cd322e. MIT license in this directory.
-/
import MathieuQuadratic.External.DvK.OneVariableConstantTerm
import MathieuQuadratic.External.DvK.SimpleRootPartialFractions
import MathieuQuadratic.External.DvK.AnnulusPartialFractions
import Mathlib.FieldTheory.SplittingField.Construction
import Mathlib.NumberTheory.FunctionField

open scoped WithZero

set_option maxHeartbeats 1500000

namespace MathieuQuadratic.DvKProof

local instance coeLaurentOneVariableDvK {R : Type*} [Semiring R] : CoeFun (LaurentPolynomial R) (fun _ ↦ ℤ → R) :=
  ⟨fun p ↦ p.coeff⟩


open LaurentPolynomial

noncomputable section

local instance : DecidableEq (RatFunc ℂ) := Classical.decEq _

/-- The ordinary polynomial whose roots solve `f(a) = z⁻¹`. -/
def laurentRootPolynomial
    {L : Type*} [Field L]
    (g : ℂ →+* L) (z : L) (s : ℕ) (p : Polynomial ℂ) :
    Polynomial L :=
  Polynomial.X ^ s - Polynomial.C z * p.map g

theorem laurentRootPolynomial_eval
    {L : Type*} [Field L]
    (g : ℂ →+* L) (z a : L) (s : ℕ) (p : Polynomial ℂ) :
    (laurentRootPolynomial g z s p).eval a =
      a ^ s - z * p.eval₂ g a := by
  simp [laurentRootPolynomial, Polynomial.eval_map]

/--
At the valuation at infinity, every partial-fraction summand arising from a
root of `X^s - z p(X)` has a nonzero derivative denominator and valuation
strictly below that of `z`.
-/
theorem valuation_partialFraction_root_data
    {L Γ : Type*} [Field L] [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation L Γ) (g : ℂ →+* L)
    (hunit : ∀ c : ℂ, c ≠ 0 → v (g c) = 1)
    (f : LaurentPolynomial ℂ)
    (hnegative : ∃ i ∈ f.support, i < 0)
    (hpositive : ∃ i ∈ f.support, 0 < i)
    (s : ℕ) (p : Polynomial ℂ)
    (hspos : 0 < s)
    (hshift : Polynomial.toLaurent p = f * LaurentPolynomial.T s)
    (hpzero : p.coeff 0 ≠ 0)
    (z a : L)
    (hz : 1 < v z)
    (ha : (laurentRootPolynomial g z s p).IsRoot a) :
    (laurentRootPolynomial g z s p).derivative.eval a ≠ 0 ∧
      v (-z * a ^ (s - 1) /
        (laurentRootPolynomial g z s p).derivative.eval a) < v z := by
  have hz0 : z ≠ 0 := by
    intro hzero
    simp [hzero] at hz
  have ha0 : a ≠ 0 := by
    intro hzero
    rw [Polynomial.IsRoot, laurentRootPolynomial_eval, hzero,
      zero_pow hspos.ne'] at ha
    have hpg0 : p.eval₂ g 0 ≠ 0 := by
      simpa [Polynomial.eval₂_at_zero] using
        (map_ne_zero g).mpr hpzero
    have hzeroEval : z * p.eval₂ g 0 = 0 := by
      simpa using ha
    exact hpg0 ((mul_eq_zero.mp hzeroEval).resolve_left hz0)
  let au : Lˣ := Units.mk0 a ha0
  have hrootEq :
      p.eval₂ g a = z⁻¹ * a ^ s := by
    rw [Polynomial.IsRoot, laurentRootPolynomial_eval,
      sub_eq_zero] at ha
    exact (eq_inv_mul_iff_mul_eq₀ hz0).2 ha.symm
  have hlaurentEq :
      LaurentPolynomial.eval₂ g au f = z⁻¹ := by
    have hEvalShift :=
      congrArg (LaurentPolynomial.eval₂ g au) hshift
    rw [LaurentPolynomial.eval₂_toLaurent,
      map_mul, LaurentPolynomial.eval₂_T] at hEvalShift
    change p.eval₂ g a =
      LaurentPolynomial.eval₂ g au f * a ^ s at hEvalShift
    rw [hrootEq] at hEvalShift
    exact (mul_right_cancel₀ (pow_ne_zero s ha0)) hEvalShift.symm
  have hzinv : v z⁻¹ < 1 := by
    rw [map_inv₀]
    exact (inv_lt_one₀ (lt_trans (by simp) hz)).2 hz
  have hava : v a = 1 := by
    apply valuation_laurent_eval₂_eq_one_of_two_sided
      v g au f hnegative hpositive
    · intro i hi
      exact hunit (f i) (Finsupp.mem_support_iff.mp hi)
    · simpa [hlaurentEq] using hzinv
  have hp : p ≠ 0 := by
    intro hp
    exact hpzero (by simp [hp])
  have hpEvalNe : p.eval₂ g a ≠ 0 := by
    rw [hrootEq]
    exact mul_ne_zero (inv_ne_zero hz0) (pow_ne_zero s ha0)
  have hpEvalVal :
      v (p.eval₂ g a) = v z⁻¹ := by
    rw [hrootEq, map_mul, map_pow, hava, one_pow, mul_one]
  have hpEvalLt : v (p.eval₂ g a) < 1 := by
    rw [hpEvalVal]
    exact hzinv
  have hpDeriv :
      v (p.eval₂ g a) < v (p.derivative.eval₂ g a) :=
    valuation_polynomial_derivative_gt_of_eval_lt_one
      v g hunit p hp a hava hpEvalNe hpEvalLt
  have hzDeriv :
      1 < v (z * p.derivative.eval₂ g a) := by
    rw [map_mul]
    have hzinvEq : v z⁻¹ = (v z)⁻¹ := map_inv₀ v z
    rw [hpEvalVal, hzinvEq] at hpDeriv
    rw [mul_comm]
    exact (inv_lt_iff_one_lt_mul₀
      (lt_trans (by simp) hz)).mp hpDeriv
  have hfirst :
      v ((s : L) * a ^ (s - 1)) ≤ 1 := by
    rw [map_mul, map_pow, hava, one_pow, mul_one]
    exact valuation_natCast_le_one v s
  have hderivEval :
      (laurentRootPolynomial g z s p).derivative.eval a =
        (s : L) * a ^ (s - 1) -
          z * p.derivative.eval₂ g a := by
    simp [laurentRootPolynomial, Polynomial.derivative_sub,
      Polynomial.derivative_pow, Polynomial.derivative_map,
      Polynomial.eval_map]
  have hUderiv :
      v ((laurentRootPolynomial g z s p).derivative.eval a) =
        v (z * p.derivative.eval₂ g a) := by
    rw [hderivEval]
    exact v.map_sub_eq_of_lt_right
      (lt_of_le_of_lt hfirst hzDeriv)
  constructor
  · intro hzero
    have hpositive :
        0 < v (z * p.derivative.eval₂ g a) :=
      lt_trans (by simp) hzDeriv
    rw [← hUderiv, hzero, map_zero] at hpositive
    exact (lt_irrefl 0) hpositive
  · rw [v.map_div, v.map_mul, v.map_neg, v.map_pow, hava,
      one_pow, mul_one, hUderiv]
    apply (div_lt_iff₀ (lt_trans (by simp) hzDeriv)).2
    exact lt_mul_of_one_lt_right
      (lt_trans (by simp) hz) hzDeriv

/--
The valuation bound from `valuation_partialFraction_root_data`.
-/
theorem valuation_partialFraction_summand_lt
    {L Γ : Type*} [Field L] [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation L Γ) (g : ℂ →+* L)
    (hunit : ∀ c : ℂ, c ≠ 0 → v (g c) = 1)
    (f : LaurentPolynomial ℂ)
    (hnegative : ∃ i ∈ f.support, i < 0)
    (hpositive : ∃ i ∈ f.support, 0 < i)
    (s : ℕ) (p : Polynomial ℂ)
    (hspos : 0 < s)
    (hshift : Polynomial.toLaurent p = f * LaurentPolynomial.T s)
    (hpzero : p.coeff 0 ≠ 0)
    (z a : L)
    (hz : 1 < v z)
    (ha : (laurentRootPolynomial g z s p).IsRoot a) :
    v (-z * a ^ (s - 1) /
      (laurentRootPolynomial g z s p).derivative.eval a) < v z :=
  (valuation_partialFraction_root_data
    v g hunit f hnegative hpositive s p hspos hshift hpzero
    z a hz ha).2

/-- Root polynomial over the rational function field `ℂ(z)`. -/
def dvkRootPolynomial (s : ℕ) (p : Polynomial ℂ) :
    Polynomial (RatFunc ℂ) :=
  laurentRootPolynomial RatFunc.C RatFunc.X s p

/--
The DvK root polynomial has the same degree as the shifted polynomial when
the latter has a term above the shift.
-/
theorem dvkRootPolynomial_natDegree
    (s : ℕ) (p : Polynomial ℂ)
    (hhigh : ∃ j ∈ p.support, s < j) :
    (dvkRootPolynomial s p).natDegree = p.natDegree := by
  obtain ⟨j, hj, hsj⟩ := hhigh
  have hsp : s < p.natDegree :=
    lt_of_lt_of_le hsj
      (Polynomial.le_natDegree_of_mem_supp j hj)
  have hright :
      (Polynomial.C RatFunc.X * p.map RatFunc.C).natDegree =
        p.natDegree := by
    rw [Polynomial.natDegree_C_mul RatFunc.X_ne_zero,
      Polynomial.natDegree_map]
  rw [dvkRootPolynomial, laurentRootPolynomial,
    Polynomial.natDegree_sub_eq_right_of_natDegree_lt]
  · exact hright
  · rw [Polynomial.natDegree_X_pow, hright]
    exact hsp

/-- A fixed splitting field for the DvK root polynomial. -/
abbrev DvKSplittingField (s : ℕ) (p : Polynomial ℂ) :=
  (dvkRootPolynomial s p).SplittingField

local instance (s : ℕ) (p : Polynomial ℂ) :
    DecidableEq (DvKSplittingField s p) :=
  Classical.decEq _

/-- The DvK polynomial after mapping to its fixed splitting field. -/
def dvkSplittingPolynomial (s : ℕ) (p : Polynomial ℂ) :
    Polynomial (DvKSplittingField s p) :=
  (dvkRootPolynomial s p).map
    (algebraMap (RatFunc ℂ) (DvKSplittingField s p))

/-- Extension to the splitting field of the valuation at `z = ∞`. -/
noncomputable def dvkInfinityValuation (s : ℕ) (p : Polynomial ℂ) :
    Valuation (DvKSplittingField s p) ℤᵐ⁰ :=
  extendedDiscreteValuation
    (L := DvKSplittingField s p)
    (RatFunc.inftyValuation ℂ)

instance dvkInfinityValuation_hasExtension (s : ℕ) (p : Polynomial ℂ) :
    (RatFunc.inftyValuation ℂ).HasExtension
      (dvkInfinityValuation s p) :=
  extendedDiscreteValuation_hasExtension
    (L := DvKSplittingField s p)
    (RatFunc.inftyValuation ℂ)

/-- Extension to the splitting field of the valuation at `z = 0`. -/
noncomputable def dvkZeroValuation (s : ℕ) (p : Polynomial ℂ) :
    Valuation (DvKSplittingField s p) ℤᵐ⁰ :=
  extendedDiscreteValuation
    (L := DvKSplittingField s p)
    ((Polynomial.idealX ℂ).valuation (RatFunc ℂ))

instance dvkZeroValuation_hasExtension (s : ℕ) (p : Polynomial ℂ) :
    ((Polynomial.idealX ℂ).valuation (RatFunc ℂ)).HasExtension
      (dvkZeroValuation s p) :=
  extendedDiscreteValuation_hasExtension
    (L := DvKSplittingField s p)
    ((Polynomial.idealX ℂ).valuation (RatFunc ℂ))

/-- The zero-adic topology on the DvK splitting field. -/
noncomputable instance dvkZeroValued (s : ℕ) (p : Polynomial ℂ) :
    Valued (DvKSplittingField s p) ℤᵐ⁰ :=
  Valued.mk' (dvkZeroValuation s p)

/-- The zero-adic topology is nonarchimedean. -/
instance dvkZeroNonarchimedean (s : ℕ) (p : Polynomial ℂ) :
    NonarchimedeanRing (DvKSplittingField s p) :=
  (dvkZeroValuation s p).subgroups_basis.nonarchimedean

/-- Completion of the DvK splitting field at the extended `z`-adic valuation. -/
abbrev DvKZeroCompletion (s : ℕ) (p : Polynomial ℂ) :=
  UniformSpace.Completion (DvKSplittingField s p)

/-- Nonzero complex constants are zero-adic valuation units. -/
theorem dvk_zeroValuation_constant
    (s : ℕ) (p : Polynomial ℂ) (c : ℂ) (hc : c ≠ 0) :
    dvkZeroValuation s p
        ((algebraMap (RatFunc ℂ) (DvKSplittingField s p))
          (RatFunc.C c)) = 1 := by
  let K := RatFunc ℂ
  let L := DvKSplittingField s p
  let vK := (Polynomial.idealX ℂ).valuation K
  let vL := dvkZeroValuation s p
  letI : vK.HasExtension vL := by
    simpa [vK, vL] using
      (dvkZeroValuation_hasExtension s p)
  apply
    (Valuation.HasExtension.val_map_eq_one_iff vK vL
      (RatFunc.C c)).mpr
  change
    (Polynomial.idealX ℂ).valuation K (RatFunc.C c) = 1
  rw [← RatFunc.algebraMap_C,
    IsDedekindDomain.HeightOneSpectrum.valuation_of_algebraMap,
    IsDedekindDomain.HeightOneSpectrum.intValuation_eq_one_iff_mem_primeCompl]
  simp [Ideal.primeCompl, Polynomial.idealX_span,
    Ideal.mem_span_singleton, Polynomial.X_dvd_iff, hc]

/-- The parameter `z` has zero-adic valuation strictly below one. -/
theorem dvk_zeroValuation_X_lt_one
    (s : ℕ) (p : Polynomial ℂ) :
    dvkZeroValuation s p
        ((algebraMap (RatFunc ℂ) (DvKSplittingField s p))
          RatFunc.X) < 1 := by
  let K := RatFunc ℂ
  let L := DvKSplittingField s p
  let vK := (Polynomial.idealX ℂ).valuation K
  let vL := dvkZeroValuation s p
  letI : vK.HasExtension vL := by
    simpa [vK, vL] using
      (dvkZeroValuation_hasExtension s p)
  have hbase : vK RatFunc.X < vK (1 : K) := by
    have hX :
        vK RatFunc.X = WithZero.exp (-1 : ℤ) := by
      change
        (Polynomial.idealX ℂ).valuation K RatFunc.X =
          WithZero.exp (-1 : ℤ)
      exact Polynomial.valuation_X_eq_neg_one ℂ
    rw [hX, show vK (1 : K) = 1 by simp,
      ← WithZero.exp_zero, WithZero.exp_lt_exp]
    omega
  simpa [K, L, vK, vL] using
    (Valuation.HasExtension.val_map_lt_iff vK vL
      RatFunc.X (1 : K)).mpr hbase

/--
No root of the DvK polynomial lies on the unit circle for the extended
`z`-adic valuation.
-/
theorem dvk_root_zeroValuation_ne_one
    (s : ℕ) (p : Polynomial ℂ)
    (a : DvKSplittingField s p)
    (ha : (dvkSplittingPolynomial s p).IsRoot a) :
    dvkZeroValuation s p a ≠ 1 := by
  let K := RatFunc ℂ
  let L := DvKSplittingField s p
  let vK := (Polynomial.idealX ℂ).valuation K
  let vL := dvkZeroValuation s p
  letI : vK.HasExtension vL := by
    simpa [vK, vL] using
      (dvkZeroValuation_hasExtension s p)
  let g : ℂ →+* L :=
    (algebraMap K L).comp RatFunc.C
  let z : L := algebraMap K L RatFunc.X
  have hunit : ∀ c : ℂ, c ≠ 0 → vL (g c) = 1 := by
    intro c hc
    apply
      (Valuation.HasExtension.val_map_eq_one_iff vK vL
        (RatFunc.C c)).mpr
    change
      (Polynomial.idealX ℂ).valuation K (RatFunc.C c) = 1
    rw [← RatFunc.algebraMap_C,
      IsDedekindDomain.HeightOneSpectrum.valuation_of_algebraMap,
      IsDedekindDomain.HeightOneSpectrum.intValuation_eq_one_iff_mem_primeCompl]
    simp [Ideal.primeCompl, Polynomial.idealX_span,
      Ideal.mem_span_singleton, Polynomial.X_dvd_iff, hc]
  have hz : vL z < 1 := by
    have hbase : vK RatFunc.X < vK (1 : K) := by
      have hX :
          vK RatFunc.X = WithZero.exp (-1 : ℤ) := by
        change
          (Polynomial.idealX ℂ).valuation K RatFunc.X =
            WithZero.exp (-1 : ℤ)
        exact Polynomial.valuation_X_eq_neg_one ℂ
      rw [hX, show vK (1 : K) = 1 by simp,
        ← WithZero.exp_zero, WithZero.exp_lt_exp]
      omega
    simpa [z, vK, vL] using
      (Valuation.HasExtension.val_map_lt_iff vK vL
        RatFunc.X (1 : K)).mpr hbase
  have hpoly :
      dvkSplittingPolynomial s p =
        laurentRootPolynomial g z s p := by
    ext n
    simp [dvkSplittingPolynomial, dvkRootPolynomial,
      laurentRootPolynomial, g, z, K, L]
  have ha' :
      (laurentRootPolynomial g z s p).IsRoot a := by
    rw [← hpoly]
    exact ha
  intro haunit
  have hpEval :
      vL (p.eval₂ g a) ≤ 1 := by
    apply valuation_polynomial_eval₂_le_one
      vL g p a
    · intro i hi
      exact (hunit (p.coeff i)
        (Polynomial.mem_support_iff.mp hi)).le
    · exact haunit.le
  have hright :
      vL (z * p.eval₂ g a) < 1 := by
    rw [map_mul]
    exact mul_lt_one_of_nonneg_of_lt_one_left
      zero_le hz hpEval
  have hroot :
      a ^ s = z * p.eval₂ g a := by
    rw [Polynomial.IsRoot, laurentRootPolynomial_eval,
      sub_eq_zero] at ha'
    exact ha'
  have hleft : vL (a ^ s) = 1 := by
    rw [map_pow, haunit, one_pow]
  rw [hroot] at hleft
  rw [hleft] at hright
  exact (lt_irrefl 1) hright

/-- Roots selected by the `z`-adic annulus expansion. -/
noncomputable def dvkLargeRoots (s : ℕ) (p : Polynomial ℂ) :
    Finset (DvKSplittingField s p) :=
  (dvkSplittingPolynomial s p).roots.toFinset.filter fun a =>
    1 < dvkZeroValuation s p a

/-- The root sum selected by the `z`-adic annulus expansion. -/
noncomputable def dvkSelectedRootSum
    (s : ℕ) (p : Polynomial ℂ) :
    DvKSplittingField s p :=
  let z :=
    algebraMap (RatFunc ℂ) (DvKSplittingField s p) RatFunc.X
  ∑ a ∈ dvkLargeRoots s p,
    -z * a ^ (s - 1) /
      (dvkSplittingPolynomial s p).derivative.eval a

/--
The annulus-series expansion of the DvK partial fractions. Its coefficients
live in the completion of the splitting field at `z = 0`.
-/
noncomputable def dvkAnnulusPartialFractionSeries
    (s : ℕ) (p : Polynomial ℂ) :
    AnnulusSeries (DvKZeroCompletion s p) := by
  let L := DvKSplittingField s p
  let C := DvKZeroCompletion s p
  let U := dvkSplittingPolynomial s p
  let roots := U.roots.toFinset
  let z := algebraMap (RatFunc ℂ) L RatFunc.X
  apply AnnulusSeries.indexedLinearInverseSum roots
    (fun a : L ↦ (a : C))
    (fun a : L ↦
      ((z * a ^ s / U.derivative.eval a : L) : C))
  intro a ha
  rw [Valued.valuedCompletion_apply]
  change dvkZeroValuation s p a ≠ 1
  apply dvk_root_zeroValuation_ne_one s p a
  apply Polynomial.isRoot_of_mem_roots
  exact Multiset.mem_toFinset.mp ha

/--
The constant coefficient of the annulus partial-fraction series is the
selected DvK root sum, after embedding into the zero-adic completion.
-/
theorem dvkAnnulusPartialFractionSeries_coeff_zero
    (s : ℕ) (p : Polynomial ℂ) (hspos : 0 < s) :
    dvkAnnulusPartialFractionSeries s p 0 =
      (dvkSelectedRootSum s p : DvKZeroCompletion s p) := by
  let L := DvKSplittingField s p
  let C := DvKZeroCompletion s p
  let U := dvkSplittingPolynomial s p
  let roots := U.roots.toFinset
  let z := algebraMap (RatFunc ℂ) L RatFunc.X
  rw [dvkAnnulusPartialFractionSeries,
    AnnulusSeries.indexedLinearInverseSum_coeff_zero]
  change
    (∑ a ∈ roots.filter fun a : L ↦
      1 < Valued.v (a : C),
      ((z * a ^ s / U.derivative.eval a : L) : C) *
        (-((a : C)⁻¹))) =
      ((∑ a ∈ dvkLargeRoots s p,
        -z * a ^ (s - 1) / U.derivative.eval a : L) : C)
  change
    (∑ a ∈ roots.filter fun a : L ↦
      1 < Valued.v (a : C),
      ((z * a ^ s / U.derivative.eval a : L) : C) *
        (-((a : C)⁻¹))) =
      UniformSpace.Completion.coeRingHom
        (∑ a ∈ dvkLargeRoots s p,
          -z * a ^ (s - 1) / U.derivative.eval a : L)
  rw [map_sum]
  simp only [map_div₀, map_mul, map_neg, map_pow]
  rw [show roots.filter (fun a : L ↦ 1 < Valued.v (a : C)) =
      dvkLargeRoots s p by
    ext a
    simp only [Finset.mem_filter, dvkLargeRoots]
    change
      (a ∈ U.roots.toFinset ∧ 1 < Valued.v (a : C)) ↔
        (a ∈ U.roots.toFinset ∧
          1 < dvkZeroValuation s p a)
    constructor
    · rintro ⟨ha, hval⟩
      rw [Valued.valuedCompletion_apply] at hval
      change 1 < dvkZeroValuation s p a at hval
      exact ⟨ha, hval⟩
    · rintro ⟨ha, hval⟩
      refine ⟨ha, ?_⟩
      rw [Valued.valuedCompletion_apply]
      change 1 < dvkZeroValuation s p a
      exact hval]
  apply Finset.sum_congr rfl
  intro a ha
  have haRoot : U.IsRoot a := by
    apply Polynomial.isRoot_of_mem_roots
    exact Multiset.mem_toFinset.mp
      (Finset.mem_filter.mp ha).1
  have haValNe : dvkZeroValuation s p a ≠ 1 := by
    apply dvk_root_zeroValuation_ne_one s p a
    simpa [U, L, dvkSplittingPolynomial] using haRoot
  have ha0 : a ≠ 0 := by
    intro hzero
    have hlarge := (Finset.mem_filter.mp ha).2
    simp [hzero] at hlarge
  have hs : s - 1 + 1 = s := by omega
  have ha0C : (a : C) ≠ 0 := by
    simpa using ha0
  have hmap :
      ((z * a ^ s / U.derivative.eval a : L) : C) =
        (z : C) * (a : C) ^ s /
          (U.derivative.eval a : L) := by
    change
      UniformSpace.Completion.coeRingHom
          (z * a ^ s / U.derivative.eval a) =
        UniformSpace.Completion.coeRingHom z *
            UniformSpace.Completion.coeRingHom a ^ s /
          UniformSpace.Completion.coeRingHom
            (U.derivative.eval a)
    rw [map_div₀, map_mul, map_pow]
  rw [hmap]
  change
    ((z : C) * (a : C) ^ s /
          ((U.derivative.eval a : L) : C)) *
        (-((a : C)⁻¹)) =
      (-(z : C)) * (a : C) ^ (s - 1) /
        ((U.derivative.eval a : L) : C)
  have hpow :
      (a : C) ^ s = (a : C) ^ (s - 1) * (a : C) :=
    (congrArg (fun n : ℕ ↦ (a : C) ^ n) hs.symm).trans
      (pow_succ _ _)
  rw [hpow, div_eq_mul_inv, div_eq_mul_inv]
  calc
    _ = -((z : C) * (a : C) ^ (s - 1) *
          ((U.derivative.eval a : L) : C)⁻¹) *
        ((a : C) * (a : C)⁻¹) := by
          rw [mul_neg, neg_mul]
          congr 1
          simp only [mul_assoc]
          rw [mul_left_comm
            (a : C)
            (((U.derivative.eval a : L) : C)⁻¹)
            ((a : C)⁻¹)]
    _ = _ := by
      rw [mul_inv_cancel₀ ha0C]
      rw [mul_one, neg_mul, neg_mul]

/--
The annulus partial-fraction series satisfies the DvK polynomial equation
in the zero-adic completion.
-/
theorem dvkAnnulusPartialFractionSeries_polynomialEquation
    (s : ℕ) (p : Polynomial ℂ)
    (hspos : 0 < s)
    (hhigh : ∃ j ∈ p.support, s < j)
    (hsimple :
      (dvkSplittingPolynomial s p).roots.Nodup) :
    let L := DvKSplittingField s p
    let C := DvKZeroCompletion s p
    let U := dvkSplittingPolynomial s p
    let φ : L →+* C := UniformSpace.Completion.coeRingHom
    let z := algebraMap (RatFunc ℂ) L RatFunc.X
    AnnulusSeries.laurentAction
        (Polynomial.toLaurent (U.map φ))
        (dvkAnnulusPartialFractionSeries s p) =
      AnnulusSeries.scale (φ z)
        (AnnulusSeries.single s 1) := by
  classical
  let L := DvKSplittingField s p
  let C := DvKZeroCompletion s p
  let U := dvkSplittingPolynomial s p
  let roots := U.roots.toFinset
  let φ : L →+* C := UniformSpace.Completion.coeRingHom
  let z := algebraMap (RatFunc ℂ) L RatFunc.X
  let Q : L → Polynomial C := fun a ↦
    (U /ₘ (Polynomial.X - Polynomial.C a)).map φ
  have hbaseDegree :
      (dvkRootPolynomial s p).natDegree = p.natDegree :=
    dvkRootPolynomial_natDegree s p hhigh
  have hUdegree : U.natDegree = p.natDegree := by
    rw [show U =
      (dvkRootPolynomial s p).map
        (algebraMap (RatFunc ℂ) L) by
          rfl,
      Polynomial.natDegree_map, hbaseDegree]
  have hdegree : s < U.natDegree := by
    rw [hUdegree]
    obtain ⟨j, hj, hsj⟩ := hhigh
    exact lt_of_lt_of_le hsj
      (Polynomial.le_natDegree_of_mem_supp j hj)
  have hU : U ≠ 0 := by
    intro hzero
    have hpositiveDegree : 0 < U.natDegree :=
      lt_trans hspos hdegree
    rw [hzero, Polynomial.natDegree_zero] at hpositiveDegree
    omega
  have hsplit : U.Splits := by
    simpa [U, L, dvkSplittingPolynomial] using
      (Polynomial.SplittingField.splits
        (dvkRootPolynomial s p))
  have hsimpleU : U.roots.Nodup := by
    simpa [U] using hsimple
  have hpoly :
      Polynomial.X ^ s =
        ∑ a ∈ roots,
          Polynomial.C (a ^ s / U.derivative.eval a) *
            (U /ₘ (Polynomial.X - Polynomial.C a)) :=
    polynomial_eq_sum_divByMonic_over_simple_roots
      U hU hsplit hsimpleU s hdegree
  have hunit :
      ∀ a ∈ roots, Valued.v (φ a) ≠ 1 := by
    intro a ha
    change Valued.v (a : C) ≠ 1
    rw [Valued.valuedCompletion_apply]
    change dvkZeroValuation s p a ≠ 1
    apply dvk_root_zeroValuation_ne_one s p a
    apply Polynomial.isRoot_of_mem_roots
    exact Multiset.mem_toFinset.mp ha
  have hfactor :
      ∀ a ∈ roots,
        U.map φ =
          Q a * (Polynomial.X - Polynomial.C (φ a)) := by
    intro a ha
    have haRoot : U.IsRoot a := by
      apply Polynomial.isRoot_of_mem_roots
      exact Multiset.mem_toFinset.mp ha
    have hfactorL :
        U =
          (U /ₘ (Polynomial.X - Polynomial.C a)) *
            (Polynomial.X - Polynomial.C a) := by
      rw [mul_comm,
        Polynomial.mul_divByMonic_eq_iff_isRoot.mpr haRoot]
    calc
      U.map φ =
          ((U /ₘ (Polynomial.X - Polynomial.C a)) *
            (Polynomial.X - Polynomial.C a)).map φ := by
        exact congrArg (Polynomial.map φ) hfactorL
      _ = Q a * (Polynomial.X - Polynomial.C (φ a)) := by
        simp [Q]
  have hidentity :
      Polynomial.C (φ z) * Polynomial.X ^ s =
        ∑ a ∈ roots,
          Polynomial.C
              (φ (z * a ^ s / U.derivative.eval a)) *
            Q a := by
    have hmapped :=
      congrArg (Polynomial.map φ) hpoly
    simp only [Polynomial.map_pow, Polynomial.map_X,
      Polynomial.map_sum, Polynomial.map_mul,
      Polynomial.map_C] at hmapped
    calc
      Polynomial.C (φ z) * Polynomial.X ^ s =
          Polynomial.C (φ z) *
            ∑ a ∈ roots,
              Polynomial.C
                  (φ (a ^ s / U.derivative.eval a)) *
                Q a := by
        rw [hmapped]
      _ = ∑ a ∈ roots,
          Polynomial.C
              (φ (z * a ^ s / U.derivative.eval a)) *
            Q a := by
        rw [Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro a ha
        simp only [← mul_assoc, ← Polynomial.C_mul]
        congr 2
        simp only [map_div₀, map_mul, map_pow]
        ring
  have haction :=
    AnnulusSeries.laurentAction_indexedLinearInverseSum_of_identity
      roots
      (fun a : L ↦ φ a)
      (fun a : L ↦ φ (z * a ^ s / U.derivative.eval a))
      hunit
      (U.map φ)
      (Polynomial.C (φ z) * Polynomial.X ^ s)
      Q hfactor hidentity
  change
    AnnulusSeries.laurentAction
        (Polynomial.toLaurent (U.map φ))
        (AnnulusSeries.indexedLinearInverseSum roots
          (fun a : L ↦ φ a)
          (fun a : L ↦ φ
            (z * a ^ s / U.derivative.eval a))
          hunit) =
      AnnulusSeries.scale (φ z)
        (AnnulusSeries.single s 1)
  rw [haction]
  rw [map_mul, Polynomial.toLaurent_C,
    Polynomial.toLaurent_X_pow,
    AnnulusSeries.laurentAction_mul,
    AnnulusSeries.laurentAction_T_single_zero,
    AnnulusSeries.laurentAction_C]
  rfl

/--
After cancelling the shift `T^s`, the DvK annulus series is a fixed point
of multiplication by `zf`.
-/
theorem dvkAnnulusPartialFractionSeries_fixedPoint
    (f : LaurentPolynomial ℂ)
    (s : ℕ) (p : Polynomial ℂ)
    (hspos : 0 < s)
    (hshift : Polynomial.toLaurent p = f * LaurentPolynomial.T s)
    (hhigh : ∃ j ∈ p.support, s < j)
    (hsimple :
      (dvkSplittingPolynomial s p).roots.Nodup) :
    let L := DvKSplittingField s p
    let C := DvKZeroCompletion s p
    let φ : L →+* C := UniformSpace.Completion.coeRingHom
    let gL : ℂ →+* L :=
      (algebraMap (RatFunc ℂ) L).comp RatFunc.C
    let g : ℂ →+* C := φ.comp gL
    let z := algebraMap (RatFunc ℂ) L RatFunc.X
    let fC : LaurentPolynomial C :=
      AddMonoidAlgebra.mapRingHom ℤ g f
    let q := LaurentPolynomial.C (φ z) * fC
    AnnulusSeries.laurentAction (1 - q)
        (dvkAnnulusPartialFractionSeries s p) =
      AnnulusSeries.scale (φ z)
        (AnnulusSeries.single 0 1) := by
  classical
  let L := DvKSplittingField s p
  let C := DvKZeroCompletion s p
  let U := dvkSplittingPolynomial s p
  let φ : L →+* C := UniformSpace.Completion.coeRingHom
  let gL : ℂ →+* L :=
    (algebraMap (RatFunc ℂ) L).comp RatFunc.C
  let g : ℂ →+* C := φ.comp gL
  let z := algebraMap (RatFunc ℂ) L RatFunc.X
  let fC : LaurentPolynomial C :=
    AddMonoidAlgebra.mapRingHom ℤ g f
  let q := LaurentPolynomial.C (φ z) * fC
  have htoLaurent :
      AddMonoidAlgebra.mapRingHom ℤ g
          (Polynomial.toLaurent p) =
        Polynomial.toLaurent (p.map g) := by
    ext n
    cases n with
    | ofNat n =>
        rw [AddMonoidAlgebra.coeff_mapRingHom]
        change g (Polynomial.toLaurent p (n : ℤ)) = Polynomial.toLaurent (p.map g) (n : ℤ)
        rw [polynomial_toLaurent_apply_nat, polynomial_toLaurent_apply_nat, Polynomial.coeff_map]
    | negSucc n =>
        have hpneg :
            Int.negSucc n ∉ (Polynomial.toLaurent p).coeff.support := by
          rw [LaurentPolynomial.support_coeff_toLaurent]
          simp
        have hmapneg :
            Int.negSucc n ∉
              (Polynomial.toLaurent (p.map g)).coeff.support := by
          rw [LaurentPolynomial.support_coeff_toLaurent]
          simp
        rw [AddMonoidAlgebra.coeff_mapRingHom,
          Finsupp.notMem_support_iff.mp hpneg,
          Finsupp.notMem_support_iff.mp hmapneg,
          map_zero]
  have hshiftC :
      Polynomial.toLaurent (p.map g) =
        fC * LaurentPolynomial.T s := by
    have hT :
        AddMonoidAlgebra.mapRingHom ℤ g
            (LaurentPolynomial.T s) =
          LaurentPolynomial.T s := by
      change
        AddMonoidAlgebra.mapRingHom ℤ g
            (AddMonoidAlgebra.single (s : ℤ) 1) =
          AddMonoidAlgebra.single (s : ℤ) 1
      rw [AddMonoidAlgebra.mapRingHom_single, map_one]
    calc
      Polynomial.toLaurent (p.map g) =
          AddMonoidAlgebra.mapRingHom ℤ g
            (Polynomial.toLaurent p) :=
        htoLaurent.symm
      _ = AddMonoidAlgebra.mapRingHom ℤ g
          (f * LaurentPolynomial.T s) := by rw [hshift]
      _ = AddMonoidAlgebra.mapRingHom ℤ g f *
          AddMonoidAlgebra.mapRingHom ℤ g
            (LaurentPolynomial.T s) := by rw [map_mul]
      _ = fC * LaurentPolynomial.T s := by rw [hT]
  have hUmap :
      U.map φ =
        Polynomial.X ^ s -
          Polynomial.C (φ z) * p.map g := by
    ext n
    simp [U, dvkSplittingPolynomial, dvkRootPolynomial,
      laurentRootPolynomial, φ, g, gL, z, L]
  have hfactor :
      Polynomial.toLaurent (U.map φ) =
        LaurentPolynomial.T s * (1 - q) := by
    calc
      Polynomial.toLaurent (U.map φ) =
          Polynomial.toLaurent
            (Polynomial.X ^ s -
              Polynomial.C (φ z) * p.map g) := by rw [hUmap]
      _ = LaurentPolynomial.T s -
          LaurentPolynomial.C (φ z) *
            Polynomial.toLaurent (p.map g) := by simp
      _ = LaurentPolynomial.T s -
          LaurentPolynomial.C (φ z) *
            (fC * LaurentPolynomial.T s) := by rw [hshiftC]
      _ = LaurentPolynomial.T s * (1 - q) := by
        simp only [q]
        ring
  have hpolyEq :=
    dvkAnnulusPartialFractionSeries_polynomialEquation
      s p hspos hhigh hsimple
  change
    AnnulusSeries.laurentAction
        (Polynomial.toLaurent (U.map φ))
        (dvkAnnulusPartialFractionSeries s p) =
      AnnulusSeries.scale (φ z)
        (AnnulusSeries.single s 1) at hpolyEq
  rw [hfactor, AnnulusSeries.laurentAction_mul,
    AnnulusSeries.laurentAction_T] at hpolyEq
  have hrhs :
      AnnulusSeries.scale (φ z)
          (AnnulusSeries.single s 1) =
        AnnulusSeries.shift s
          (AnnulusSeries.scale (φ z)
            (AnnulusSeries.single 0 1)) := by
    rw [AnnulusSeries.shift_scale,
      AnnulusSeries.shift_single_zero]
  rw [hrhs] at hpolyEq
  exact AnnulusSeries.shift_injective s hpolyEq

/--
The splitting-field specialization of the derivative and valuation data at
infinity.
-/
theorem dvk_root_data_at_infinity
    (f : LaurentPolynomial ℂ)
    (hnegative : ∃ i ∈ f.support, i < 0)
    (hpositive : ∃ i ∈ f.support, 0 < i)
    (s : ℕ) (p : Polynomial ℂ)
    (hspos : 0 < s)
    (hshift : Polynomial.toLaurent p = f * LaurentPolynomial.T s)
    (hpzero : p.coeff 0 ≠ 0)
    (a : DvKSplittingField s p)
    (ha :
      ((dvkRootPolynomial s p).map
        (algebraMap (RatFunc ℂ) (DvKSplittingField s p))).IsRoot a) :
    let v := dvkInfinityValuation s p
    let z := algebraMap (RatFunc ℂ) (DvKSplittingField s p) RatFunc.X
    ((dvkRootPolynomial s p).map
        (algebraMap (RatFunc ℂ) (DvKSplittingField s p))).derivative.eval a ≠ 0 ∧
      v (-z * a ^ (s - 1) /
        ((dvkRootPolynomial s p).map
          (algebraMap (RatFunc ℂ) (DvKSplittingField s p))).derivative.eval a) <
        v z := by
  let K := RatFunc ℂ
  let L := DvKSplittingField s p
  let vK := RatFunc.inftyValuation ℂ
  let vL := dvkInfinityValuation s p
  letI : vK.HasExtension vL := by
    simpa [vK, vL] using
      (dvkInfinityValuation_hasExtension s p)
  let g : ℂ →+* L :=
    (algebraMap K L).comp RatFunc.C
  let z : L := algebraMap K L RatFunc.X
  have hunit : ∀ c : ℂ, c ≠ 0 → vL (g c) = 1 := by
    intro c hc
    exact
      (Valuation.HasExtension.val_map_eq_one_iff vK vL
        (RatFunc.C c)).mpr
        (RatFunc.inftyValuation.C ℂ hc)
  have hz : 1 < vL z := by
    have hbase : vK (1 : K) < vK RatFunc.X := by
      rw [show vK (1 : K) = 1 by simp,
        show vK RatFunc.X = WithZero.exp (1 : ℤ) by
          change
            RatFunc.inftyValuation ℂ RatFunc.X =
              WithZero.exp (1 : ℤ)
          exact RatFunc.inftyValuation.X ℂ]
      rw [← WithZero.exp_zero, WithZero.exp_lt_exp]
      omega
    simpa [z, vK, vL] using
      (Valuation.HasExtension.val_map_lt_iff vK vL
        (1 : K) RatFunc.X).mpr hbase
  have hpoly :
      (dvkRootPolynomial s p).map (algebraMap K L) =
        laurentRootPolynomial g z s p := by
    ext n
    simp [dvkRootPolynomial, laurentRootPolynomial, g, z, K, L]
  have ha' :
      (laurentRootPolynomial g z s p).IsRoot a := by
    rw [← hpoly]
    exact ha
  have hdata :=
    valuation_partialFraction_root_data
      vL g hunit f hnegative hpositive s p hspos hshift hpzero
      z a hz ha'
  have hderiv :
      ((dvkRootPolynomial s p).map (algebraMap K L)).derivative.eval a =
        (laurentRootPolynomial g z s p).derivative.eval a := by
    rw [hpoly]
  constructor
  · rw [hderiv]
    exact hdata.1
  · change
      vL (-z * a ^ (s - 1) /
        ((dvkRootPolynomial s p).map
          (algebraMap K L)).derivative.eval a) < vL z
    rw [hderiv]
    exact hdata.2

/--
Every root of the DvK polynomial in its splitting field is simple.
-/
theorem dvk_root_derivative_ne_zero_at_infinity
    (f : LaurentPolynomial ℂ)
    (hnegative : ∃ i ∈ f.support, i < 0)
    (hpositive : ∃ i ∈ f.support, 0 < i)
    (s : ℕ) (p : Polynomial ℂ)
    (hspos : 0 < s)
    (hshift : Polynomial.toLaurent p = f * LaurentPolynomial.T s)
    (hpzero : p.coeff 0 ≠ 0)
    (a : DvKSplittingField s p)
    (ha :
      ((dvkRootPolynomial s p).map
        (algebraMap (RatFunc ℂ) (DvKSplittingField s p))).IsRoot a) :
    ((dvkRootPolynomial s p).map
      (algebraMap (RatFunc ℂ) (DvKSplittingField s p))).derivative.eval a ≠ 0 :=
  (dvk_root_data_at_infinity
    f hnegative hpositive s p hspos hshift hpzero a ha).1

/--
The splitting-field specialization of the per-root estimate at infinity.
-/
theorem dvk_summand_lt_at_infinity
    (f : LaurentPolynomial ℂ)
    (hnegative : ∃ i ∈ f.support, i < 0)
    (hpositive : ∃ i ∈ f.support, 0 < i)
    (s : ℕ) (p : Polynomial ℂ)
    (hspos : 0 < s)
    (hshift : Polynomial.toLaurent p = f * LaurentPolynomial.T s)
    (hpzero : p.coeff 0 ≠ 0)
    (a : DvKSplittingField s p)
    (ha :
      ((dvkRootPolynomial s p).map
        (algebraMap (RatFunc ℂ) (DvKSplittingField s p))).IsRoot a) :
    let v := dvkInfinityValuation s p
    let z := algebraMap (RatFunc ℂ) (DvKSplittingField s p) RatFunc.X
    v (-z * a ^ (s - 1) /
      ((dvkRootPolynomial s p).map
        (algebraMap (RatFunc ℂ) (DvKSplittingField s p))).derivative.eval a) <
      v z :=
  (dvk_root_data_at_infinity
    f hnegative hpositive s p hspos hshift hpzero a ha).2

/--
The complete root sum selected at `z = 0` has valuation strictly below that
of `z` at infinity.
-/
theorem dvk_selectedRootSum_lt_at_infinity
    (f : LaurentPolynomial ℂ)
    (hnegative : ∃ i ∈ f.support, i < 0)
    (hpositive : ∃ i ∈ f.support, 0 < i)
    (s : ℕ) (p : Polynomial ℂ)
    (hspos : 0 < s)
    (hshift : Polynomial.toLaurent p = f * LaurentPolynomial.T s)
    (hpzero : p.coeff 0 ≠ 0) :
    let v := dvkInfinityValuation s p
    let z :=
      algebraMap (RatFunc ℂ) (DvKSplittingField s p) RatFunc.X
    v (dvkSelectedRootSum s p) < v z := by
  let L := DvKSplittingField s p
  let v := dvkInfinityValuation s p
  let z := algebraMap (RatFunc ℂ) L RatFunc.X
  have hz0 : z ≠ 0 :=
    (map_ne_zero (algebraMap (RatFunc ℂ) L)).mpr
      RatFunc.X_ne_zero
  change
    v (∑ a ∈ dvkLargeRoots s p,
      -z * a ^ (s - 1) /
        (dvkSplittingPolynomial s p).derivative.eval a) < v z
  apply v.map_sum_lt
  · exact (v.ne_zero_iff).mpr hz0
  · intro a ha
    have haroots :
        a ∈ (dvkSplittingPolynomial s p).roots := by
      exact Multiset.mem_toFinset.mp
        (Finset.mem_filter.mp ha).1
    have haRoot :
        (dvkSplittingPolynomial s p).IsRoot a :=
      Polynomial.isRoot_of_mem_roots haroots
    simpa [z, L, dvkSplittingPolynomial] using
      (dvk_summand_lt_at_infinity
        f hnegative hpositive s p hspos hshift hpzero
        a (by simpa [dvkSplittingPolynomial] using haRoot))

/--
The coefficient-extraction statement in the one-variable DvK argument. It
identifies the selected root sum when all positive-power constant terms vanish.
-/
def DvKCoefficientExtraction : Prop :=
  ∀ (f : LaurentPolynomial ℂ)
    (_hnegative : ∃ i ∈ f.support, i < 0)
    (_hpositive : ∃ i ∈ f.support, 0 < i)
    (s : ℕ) (p : Polynomial ℂ)
    (_hspos : 0 < s)
    (_hshift :
      Polynomial.toLaurent p = f * LaurentPolynomial.T s)
    (_hpzero : p.coeff 0 ≠ 0)
    (_hhigh : ∃ j ∈ p.support, s < j),
    (∀ m : ℕ, 1 ≤ m → (f ^ m) 0 = 0) →
      dvkSelectedRootSum s p =
        algebraMap (RatFunc ℂ) (DvKSplittingField s p) RatFunc.X

/--
Coefficient extraction for the one-variable DvK annulus expansion.
-/
theorem dvkCoefficientExtraction : DvKCoefficientExtraction := by
  intro f hnegative hpositive s p hspos hshift hpzero hhigh hall
  classical
  let L := DvKSplittingField s p
  let C := DvKZeroCompletion s p
  let U := dvkSplittingPolynomial s p
  let φ : L →+* C := UniformSpace.Completion.coeRingHom
  let gL : ℂ →+* L :=
    (algebraMap (RatFunc ℂ) L).comp RatFunc.C
  let g : ℂ →+* C := φ.comp gL
  let z := algebraMap (RatFunc ℂ) L RatFunc.X
  let fC : LaurentPolynomial C :=
    AddMonoidAlgebra.mapRingHom ℤ g f
  let q : LaurentPolynomial C :=
    LaurentPolynomial.C (φ z) * fC
  let S := dvkAnnulusPartialFractionSeries s p
  let ρ : ℤᵐ⁰ := Valued.v (φ z)
  have hbaseDegree :
      (dvkRootPolynomial s p).natDegree = p.natDegree :=
    dvkRootPolynomial_natDegree s p hhigh
  have hUdegree : U.natDegree = p.natDegree := by
    rw [show U =
      (dvkRootPolynomial s p).map
        (algebraMap (RatFunc ℂ) L) by
          rfl,
      Polynomial.natDegree_map, hbaseDegree]
  have hdegree : s < U.natDegree := by
    rw [hUdegree]
    obtain ⟨j, hj, hsj⟩ := hhigh
    exact lt_of_lt_of_le hsj
      (Polynomial.le_natDegree_of_mem_supp j hj)
  have hU : U ≠ 0 := by
    intro hzero
    have hpositiveDegree : 0 < U.natDegree :=
      lt_trans hspos hdegree
    rw [hzero, Polynomial.natDegree_zero] at hpositiveDegree
    omega
  have hsimple : U.roots.Nodup := by
    apply nodup_roots_of_derivative_ne_zero_at_roots U hU
    intro a ha
    apply dvk_root_derivative_ne_zero_at_infinity
      f hnegative hpositive s p hspos hshift hpzero a
    simpa [U, L, dvkSplittingPolynomial] using ha
  have hfixedEquation :=
    dvkAnnulusPartialFractionSeries_fixedPoint
      f s p hspos hshift hhigh
        (by simpa [U] using hsimple)
  change
    AnnulusSeries.laurentAction (1 - q) S =
      AnnulusSeries.scale (φ z)
        (AnnulusSeries.single 0 1) at hfixedEquation
  have hfixed :
      S =
        AnnulusSeries.scale (φ z)
            (AnnulusSeries.single 0 1) +
          AnnulusSeries.laurentAction q S := by
    rw [AnnulusSeries.laurentAction_sub_left,
      AnnulusSeries.laurentAction_one] at hfixedEquation
    exact sub_eq_iff_eq_add.mp hfixedEquation
  have hρ : ρ < 1 := by
    change Valued.v (z : C) < 1
    rw [Valued.valuedCompletion_apply]
    exact dvk_zeroValuation_X_lt_one s p
  have hg_le_one :
      ∀ c : ℂ, Valued.v (g c) ≤ 1 := by
    intro c
    by_cases hc : c = 0
    · simp [hc]
    · have heq : Valued.v (g c) = 1 := by
        change Valued.v (gL c : C) = 1
        rw [Valued.valuedCompletion_apply]
        exact dvk_zeroValuation_constant s p c hc
      rw [heq]
  have hqcoeff :
      ∀ k : ℤ, q k = φ z * g (f k) := by
    intro k
    change
      ((AddMonoidAlgebra.single (0 : ℤ) (φ z) :
          LaurentPolynomial C) * fC) k =
        φ z * g (f k)
    rw [AddMonoidAlgebra.coeff_single_mul_apply]
    simp only [neg_zero, zero_add]
    rw [AddMonoidAlgebra.coeff_mapRingHom]
  have hq :
      ∀ k ∈ q.support, Valued.v (q k) ≤ ρ := by
    intro k hk
    rw [hqcoeff, map_mul]
    calc
      Valued.v (φ z) * Valued.v (g (f k)) ≤
          Valued.v (φ z) * 1 :=
        mul_le_mul_right (hg_le_one (f k)) _
      _ = ρ := by simp [ρ]
  have hpowzero :
      ∀ m : ℕ, 1 ≤ m → (q ^ m) 0 = 0 := by
    intro m hm
    have hfCpow :
        fC ^ m =
          AddMonoidAlgebra.mapRingHom ℤ g (f ^ m) := by
      rw [map_pow]
    have hCpow :
        (LaurentPolynomial.C (φ z)) ^ m =
          LaurentPolynomial.C ((φ z) ^ m) := by
      rw [map_pow]
    change
      ((LaurentPolynomial.C (φ z) * fC) ^ m) 0 = 0
    rw [mul_pow, hCpow, hfCpow]
    change
      ((AddMonoidAlgebra.single (0 : ℤ) ((φ z) ^ m) :
          LaurentPolynomial C) *
        AddMonoidAlgebra.mapRingHom ℤ g (f ^ m)) 0 = 0
    rw [AddMonoidAlgebra.coeff_single_mul_apply]
    simp only [neg_zero, zero_add,
      AddMonoidAlgebra.coeff_mapRingHom, hall m hm, map_zero, mul_zero]
  have hcoeff :
      S 0 = φ z :=
    AnnulusSeries.fixedPoint_coeff_zero
      q S (φ z) ρ hfixed hρ hq hpowzero
  have hselected :
      φ (dvkSelectedRootSum s p) = φ z := by
    calc
      φ (dvkSelectedRootSum s p) = S 0 := by
        symm
        exact
          dvkAnnulusPartialFractionSeries_coeff_zero
            s p hspos
      _ = φ z := hcoeff
  exact UniformSpace.Completion.coe_injective L hselected

/--
The coefficient-extraction statement implies the one-variable constant-term
theorem.
-/
theorem hasNonzeroConstantPower_of_dvkCoefficientExtraction
    (hExtraction : DvKCoefficientExtraction)
    (f : LaurentPolynomial ℂ)
    (hnegative : ∃ i ∈ f.support, i < 0)
    (hpositive : ∃ i ∈ f.support, 0 < i) :
    (∃ m : ℕ, 1 ≤ m ∧ (f ^ m) 0 ≠ 0) := by
  obtain ⟨s, p, hspos, hshift, hpzero, hhigh⟩ :=
    exists_exact_laurent_shift f hnegative hpositive
  by_contra hnone
  have hall : ∀ m : ℕ, 1 ≤ m → (f ^ m) 0 = 0 := by
    intro m hm
    by_contra hcoeff
    exact hnone ⟨m, hm, hcoeff⟩
  have hsum :
      dvkSelectedRootSum s p =
        algebraMap (RatFunc ℂ) (DvKSplittingField s p) RatFunc.X :=
    hExtraction f hnegative hpositive s p hspos hshift
      hpzero hhigh hall
  have hlt :=
    dvk_selectedRootSum_lt_at_infinity
      f hnegative hpositive s p hspos hshift hpzero
  rw [hsum] at hlt
  exact (lt_irrefl _) hlt

end

end MathieuQuadratic.DvKProof
