/-
Adapted from GMC2/OneVariableConstantTerm.lean at MurrellGroup/GMC-2-lean commit
1782de7ff6c97eb1d98e63e7ff34df18b9cd322e. MIT license in this directory.
-/
import MathieuQuadratic.External.DvK.Compat
import Mathlib.Algebra.Polynomial.Laurent
import Mathlib.Analysis.Complex.Polynomial.Basic
import MathieuQuadratic.External.DvK.ValuationExtension
import Mathlib.FieldTheory.IsAlgClosed.Basic
import Mathlib.RingTheory.LaurentSeries
import Mathlib.RingTheory.LocalRing.ResidueField.Basic

open scoped WithZero

set_option maxHeartbeats 1500000

namespace MathieuQuadratic.DvKProof

local instance coeLaurentOneVariableConstantTerm {R : Type*} [Semiring R] : CoeFun (LaurentPolynomial R) (fun _ ↦ ℤ → R) :=
  ⟨fun p ↦ p.coeff⟩


open LaurentPolynomial

noncomputable section

/-- Truncation preserves every nonnegative Laurent coefficient. -/
@[simp] theorem laurent_trunc_coeff
    {R : Type*} [Semiring R]
    (f : LaurentPolynomial R) (n : ℕ) :
    (LaurentPolynomial.trunc f).coeff n = f n := by
  rfl

/-- Embedding a polynomial as a Laurent polynomial preserves its coefficients. -/
@[simp] theorem polynomial_toLaurent_apply_nat
    {R : Type*} [Semiring R]
    (p : Polynomial R) (n : ℕ) :
    Polynomial.toLaurent p (n : ℤ) = p.coeff n := by
  rw [Polynomial.toLaurent_apply]
  change p.toFinsupp.mapDomain (fun m : ℕ => (m : ℤ)) (n : ℤ) =
    p.coeff n
  exact Finsupp.mapDomain_apply Int.ofNat_injective p.toFinsupp.coeff n

/-- Truncation is inverse to `toLaurent` when no negative term is present. -/
theorem polynomial_toLaurent_trunc_of_nonnegative_support
    {R : Type*} [Semiring R]
    (f : LaurentPolynomial R)
    (hnonnegative : ∀ i ∈ f.support, 0 ≤ i) :
    Polynomial.toLaurent (LaurentPolynomial.trunc f) = f := by
  ext i
  by_cases hi : 0 ≤ i
  · lift i to ℕ using hi
    exact polynomial_toLaurent_apply_nat (trunc f) i
  · have hif : f i = 0 := by
      apply Finsupp.notMem_support_iff.mp
      intro himem
      exact hi (hnonnegative i himem)
    rw [hif]
    apply Finsupp.notMem_support_iff.mp
    rw [LaurentPolynomial.support_coeff_toLaurent]
    intro himem
    obtain ⟨n, hn, hni⟩ := Finset.mem_map.mp himem
    apply hi
    rw [← hni]
    exact Int.natCast_nonneg n

/--
Shifting a two-sided Laurent polynomial by its minimal exponent gives an
ordinary polynomial with nonzero constant coefficient and a term above the
shift.
-/
theorem exists_exact_laurent_shift
    (f : LaurentPolynomial ℂ)
    (hnegative : ∃ i ∈ f.support, i < 0)
    (hpositive : ∃ i ∈ f.support, 0 < i) :
    ∃ (s : ℕ) (p : Polynomial ℂ),
      0 < s ∧
      Polynomial.toLaurent p = f * LaurentPolynomial.T s ∧
      p.coeff 0 ≠ 0 ∧
      ∃ j ∈ p.support, s < j := by
  classical
  have hsupport : f.support.Nonempty :=
    ⟨hnegative.choose, hnegative.choose_spec.1⟩
  let lo := f.support.min' hsupport
  have hloMem : lo ∈ f.support :=
    Finset.min'_mem f.support hsupport
  have hloNeg : lo < 0 := by
    exact lt_of_le_of_lt
      (Finset.min'_le f.support hnegative.choose
        hnegative.choose_spec.1)
      hnegative.choose_spec.2
  let s := lo.natAbs
  have hspos : 0 < s := by
    simpa [s] using Int.natAbs_pos.mpr (ne_of_lt hloNeg)
  have hslo : (s : ℤ) = -lo := by
    simpa [s] using Int.ofNat_natAbs_of_nonpos hloNeg.le
  let shifted := f * LaurentPolynomial.T s
  have hshiftedSupport :
      ∀ i ∈ shifted.support, 0 ≤ i := by
    intro i hi
    change i ∈ (f * LaurentPolynomial.T (s : ℤ)).coeff.support at hi
    rw [LaurentPolynomial.T,
      AddMonoidAlgebra.support_coeff_mul_single f 1 (by simp) (s : ℤ)] at hi
    obtain ⟨e, he, hei⟩ := Finset.mem_map.mp hi
    have hlole : lo ≤ e :=
      Finset.min'_le f.support e he
    change e + (s : ℤ) = i at hei
    subst i
    omega
  let p := LaurentPolynomial.trunc shifted
  have hpLaurent : Polynomial.toLaurent p = shifted :=
    polynomial_toLaurent_trunc_of_nonnegative_support
      shifted hshiftedSupport
  have hshiftedZero : shifted 0 = f lo := by
    change
      ((f * LaurentPolynomial.T (s : ℤ) :
        LaurentPolynomial ℂ) 0) = f lo
    rw [LaurentPolynomial.T, AddMonoidAlgebra.coeff_mul_single_apply]
    simp [hslo]
  have hpzero : p.coeff 0 ≠ 0 := by
    rw [show p.coeff 0 = shifted 0 by rfl, hshiftedZero]
    exact Finsupp.mem_support_iff.mp hloMem
  let pos := hpositive.choose
  have hposMem : pos ∈ f.support := hpositive.choose_spec.1
  have hposPos : 0 < pos := hpositive.choose_spec.2
  let j := Int.toNat (pos + s)
  have hposShift : 0 ≤ pos + (s : ℤ) := by omega
  have hjCast : (j : ℤ) = pos + s := by
    simp [j, Int.toNat_of_nonneg hposShift]
  have hshiftedPos : shifted (j : ℤ) ≠ 0 := by
    change
      ((f * LaurentPolynomial.T (s : ℤ) :
        LaurentPolynomial ℂ) (j : ℤ)) ≠ 0
    rw [LaurentPolynomial.T, AddMonoidAlgebra.coeff_mul_single_apply]
    simp [hjCast]
    exact Finsupp.mem_support_iff.mp hposMem
  have hjMem : j ∈ p.support := by
    apply Polynomial.mem_support_iff.mpr
    change shifted (j : ℤ) ≠ 0
    exact hshiftedPos
  have hsj : s < j := by
    exact_mod_cast (show (s : ℤ) < (j : ℤ) by omega)
  exact ⟨s, p, hspos, by simpa [shifted] using hpLaurent,
    hpzero, j, hjMem, hsj⟩

/-- Laurent evaluation written as an additive homomorphism over the finite support. -/
def laurentEvalAddHom
    {R S : Type*} [CommRing R] [CommRing S]
    (g : R →+* S) (x : Sˣ) :
    LaurentPolynomial R →+ S :=
  (Finsupp.liftAddHom fun n : ℤ =>
    (AddMonoidHom.mulRight (x ^ n).val).comp g.toAddMonoidHom).comp
    AddMonoidAlgebra.coeffAddEquiv.toAddMonoidHom

theorem laurentEvalAddHom_apply
    {R S : Type*} [CommRing R] [CommRing S]
    (g : R →+* S) (x : Sˣ) (f : LaurentPolynomial R) :
    laurentEvalAddHom g x f =
      f.coeff.sum fun n c => g c * (x ^ n).val := by
  rfl

theorem laurentEvalAddHom_eq_eval₂
    {R S : Type*} [CommRing R] [CommRing S]
    (g : R →+* S) (x : Sˣ) :
    laurentEvalAddHom g x = (LaurentPolynomial.eval₂ g x).toAddMonoidHom := by
  apply AddMonoidHom.ext
  intro f
  induction f using LaurentPolynomial.induction_on' with
  | add p q hp hq => rw [map_add, map_add, hp, hq]
  | C_mul_T n c =>
    rw [← LaurentPolynomial.single_eq_C_mul_T c n]
    change (Finsupp.liftAddHom fun n : ℤ =>
      (AddMonoidHom.mulRight (x ^ n).val).comp g.toAddMonoidHom)
        (Finsupp.single n c) = LaurentPolynomial.eval₂ g x (AddMonoidAlgebra.single n c)
    rw [Finsupp.liftAddHom_apply_single, LaurentPolynomial.single_eq_C_mul_T,
      LaurentPolynomial.eval₂_C_mul_T]
    rfl

theorem laurent_eval₂_eq_sum
    {R S : Type*} [CommRing R] [CommRing S]
    (g : R →+* S) (x : Sˣ) (f : LaurentPolynomial R) :
    LaurentPolynomial.eval₂ g x f =
      ∑ n ∈ f.support, g (f n) * (x ^ n).val := by
  change (LaurentPolynomial.eval₂ g x).toAddMonoidHom f = _
  rw [← laurentEvalAddHom_eq_eval₂]
  rfl

/--
A polynomial evaluated at an element of the valuation ring stays in the
valuation ring when all of its coefficients do.
-/
theorem valuation_polynomial_eval₂_le_one
    {R S Γ : Type*} [CommRing R] [Field S]
    [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation S Γ) (g : R →+* S) (p : Polynomial R) (x : S)
    (hcoeff : ∀ i ∈ p.support, v (g (p.coeff i)) ≤ 1)
    (hx : v x ≤ 1) :
    v (p.eval₂ g x) ≤ 1 := by
  rw [Polynomial.eval₂_eq_sum]
  apply v.map_sum_le
  intro i hi
  rw [map_mul, map_pow]
  exact mul_le_one₀ (hcoeff i hi) zero_le
    (pow_le_one₀ zero_le hx)

/--
If the constant coefficient vanishes, evaluating at an element of valuation
strictly below one gives valuation strictly below one.
-/
theorem valuation_polynomial_eval₂_lt_one_of_coeff_zero
    {R S Γ : Type*} [CommRing R] [Field S]
    [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation S Γ) (g : R →+* S) (p : Polynomial R) (x : S)
    (hcoeff : ∀ i ∈ p.support, v (g (p.coeff i)) ≤ 1)
    (hzero : p.coeff 0 = 0)
    (hx : v x < 1) :
    v (p.eval₂ g x) < 1 := by
  rw [Polynomial.eval₂_eq_sum]
  apply v.map_sum_lt one_ne_zero
  intro i hi
  have hi0 : i ≠ 0 := by
    intro hi0
    subst i
    exact (Polynomial.mem_support_iff.mp hi) hzero
  rw [map_mul, map_pow]
  exact mul_lt_one_of_nonneg_of_lt_one_right
    (hcoeff i hi) zero_le (pow_lt_one₀ zero_le hx hi0)

theorem valuation_natCast_le_one
    {S Γ : Type*} [Field S] [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation S Γ) (n : ℕ) :
    v (n : S) ≤ 1 := by
  induction n with
  | zero => simp
  | succ n ih =>
      rw [Nat.cast_succ]
      exact v.map_add_le ih (by simp)

theorem valuation_hasseDeriv_coeff_le_one
    {S Γ : Type*} [Field S] [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation S Γ) (p : Polynomial S)
    (hcoeff : ∀ i ∈ p.support, v (p.coeff i) ≤ 1)
    (k i : ℕ) :
    v ((Polynomial.hasseDeriv k p).coeff i) ≤ 1 := by
  rw [Polynomial.hasseDeriv_coeff, map_mul]
  apply mul_le_one₀ (valuation_natCast_le_one v _) zero_le
  by_cases hzero : p.coeff (i + k) = 0
  · simp [hzero]
  · exact hcoeff (i + k) (Polynomial.mem_support_iff.mpr hzero)

theorem valuation_taylor_coeff_le_one
    {S Γ : Type*} [Field S] [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation S Γ) (p : Polynomial S) (c : S)
    (hcoeff : ∀ i ∈ p.support, v (p.coeff i) ≤ 1)
    (hc : v c ≤ 1) (i : ℕ) :
    v ((Polynomial.taylor c p).coeff i) ≤ 1 := by
  rw [Polynomial.taylor_coeff]
  rw [← Polynomial.eval₂_id]
  apply valuation_polynomial_eval₂_le_one v (RingHom.id S)
  · intro j hj
    simpa using valuation_hasseDeriv_coeff_le_one v p hcoeff i j
  · exact hc

/--
Evaluation of an integral polynomial at two residue-equivalent elements has
the same valuation whenever one of the evaluations is a valuation unit.
-/
theorem valuation_polynomial_eval_eq_one_of_sub_lt_one
    {S Γ : Type*} [Field S] [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation S Γ) (p : Polynomial S) (c x : S)
    (hcoeff : ∀ i ∈ p.support, v (p.coeff i) ≤ 1)
    (hc : v c ≤ 1)
    (hpc : v (p.eval c) = 1)
    (hxc : v (x - c) < 1) :
    v (p.eval x) = 1 := by
  let r := Polynomial.taylor c p - Polynomial.C (p.eval c)
  have hrcoeff :
      ∀ i ∈ r.support, v (r.coeff i) ≤ 1 := by
    intro i hi
    rw [show r = Polynomial.taylor c p - Polynomial.C (p.eval c) by rfl,
      Polynomial.coeff_sub]
    apply v.map_sub_le
    · exact valuation_taylor_coeff_le_one v p c hcoeff hc i
    · by_cases hi0 : i = 0
      · subst i
        simpa using hpc.le
      · rw [Polynomial.coeff_C_of_ne_zero hi0]
        simp
  have hrzero : r.coeff 0 = 0 := by
    simp [r]
  have hr :
      v (r.eval (x - c)) < 1 := by
    rw [← Polynomial.eval₂_id]
    exact valuation_polynomial_eval₂_lt_one_of_coeff_zero
      v (RingHom.id S) r (x - c) hrcoeff hrzero hxc
  have heval :
      r.eval (x - c) = p.eval x - p.eval c := by
    simp [r, Polynomial.taylor_eval_sub]
  rw [heval] at hr
  calc
    v (p.eval x) = v (p.eval c) := v.map_eq_of_sub_lt (by simpa [hpc] using hr)
    _ = 1 := hpc

/--
An element integral over an algebraically closed field already belongs to
that field. This elementwise form avoids assuming the whole extension is
integral.
-/
theorem exists_algebraMap_eq_of_isIntegral_isAlgClosed
    {k K : Type*} [Field k] [Ring K] [IsDomain K]
    [IsAlgClosed k] [Algebra k K]
    (x : K) (hx : IsIntegral k x) :
    ∃ y : k, algebraMap k K y = x := by
  refine ⟨-(minpoly k x).coeff 0, ?_⟩
  have hq : (minpoly k x).leadingCoeff = 1 :=
    minpoly.monic hx
  have hdegree : (minpoly k x).degree = 1 :=
    IsAlgClosed.degree_eq_one_of_irreducible k
      (minpoly.irreducible hx)
  have hroot : Polynomial.aeval x (minpoly k x) = 0 :=
    minpoly.aeval k x
  rw [Polynomial.eq_X_add_C_of_degree_eq_one hdegree, hq,
    Polynomial.C_1, one_mul, map_add, Polynomial.aeval_X,
    Polynomial.aeval_C,
    add_eq_zero_iff_eq_neg] at hroot
  simpa using hroot.symm

/-- Restrict a coefficient map to the valuation subring. -/
noncomputable def coefficientToValuationSubring
    {R L Γ : Type*} [CommRing R] [Field L]
    [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation L Γ) (g : R →+* L)
    (hcoeff : ∀ c, v (g c) ≤ 1) :
    R →+* v.valuationSubring :=
  g.codRestrict v.valuationSubring hcoeff

/--
If the residue of a valuation-integral element satisfies a monic polynomial
over `ℂ`, then the element is congruent modulo the maximal ideal to a complex
constant.
-/
theorem exists_complex_approximation_of_monic_residue_root
    {L Γ : Type*} [Field L] [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation L Γ) (g : ℂ →+* L)
    (hcoeff : ∀ c, v (g c) ≤ 1)
    (x : L) (hx : v x ≤ 1)
    (q : Polynomial ℂ) (hq : q.Monic)
    (hroot :
      q.eval₂
        ((IsLocalRing.residue v.valuationSubring).comp
          (coefficientToValuationSubring v g hcoeff))
        (IsLocalRing.residue v.valuationSubring ⟨x, hx⟩) = 0) :
    ∃ c : ℂ, v (x - g c) < 1 := by
  let O := v.valuationSubring
  let gO : ℂ →+* O := coefficientToValuationSubring v g hcoeff
  letI : Algebra ℂ O := gO.toAlgebra
  let E := IsLocalRing.ResidueField O
  letI : Algebra ℂ E := ((IsLocalRing.residue O).comp gO).toAlgebra
  let xO : O := ⟨x, hx⟩
  let xbar : E := IsLocalRing.residue O xO
  have hroot' :
      q.eval₂ (algebraMap ℂ E) xbar = 0 := by
    exact hroot
  have hxbarIntegral : IsIntegral ℂ xbar :=
    ⟨q, hq, hroot'⟩
  obtain ⟨c, hc⟩ :=
    exists_algebraMap_eq_of_isIntegral_isAlgClosed xbar hxbarIntegral
  let y : O := xO - algebraMap ℂ O c
  have hyres : IsLocalRing.residue O y = 0 := by
    rw [map_sub]
    change xbar - algebraMap ℂ E c = 0
    rw [hc]
    exact sub_self xbar
  have hynonunit : ¬IsUnit y := by
    rw [← IsLocalRing.residue_ne_zero_iff_isUnit]
    exact not_ne_iff.mpr hyres
  have hyne : v (y : L) ≠ 1 := by
    intro hy
    apply hynonunit
    exact
      ((Valuation.valuationSubring.integers v).isUnit_iff_valuation_eq_one).mpr
        (by simpa [O] using hy)
  have hylt : v (y : L) < 1 :=
    lt_of_le_of_ne y.property hyne
  refine ⟨c, ?_⟩
  exact hylt

/--
A valuation-unit argument at which a nonzero complex polynomial has value
below one is residue-equivalent to an actual complex root of the reduced
polynomial.
-/
theorem exists_complex_approximation_of_polynomial_eval_lt_one
    {L Γ : Type*} [Field L] [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation L Γ) (g : ℂ →+* L)
    (hunit : ∀ c : ℂ, c ≠ 0 → v (g c) = 1)
    (p : Polynomial ℂ) (hp : p ≠ 0)
    (x : L) (hx : v x = 1)
    (heval : v (p.eval₂ g x) < 1) :
    ∃ c : ℂ, v (x - g c) < 1 := by
  have hcoeff : ∀ c : ℂ, v (g c) ≤ 1 := by
    intro c
    by_cases hc : c = 0
    · simp [hc]
    · exact (hunit c hc).le
  let O := v.valuationSubring
  let gO : ℂ →+* O := coefficientToValuationSubring v g hcoeff
  let xO : O := ⟨x, hx.le⟩
  let zO : O := p.eval₂ gO xO
  have hzcoe : (zO : L) = p.eval₂ g x := by
    change
      ((O.subtype : O →+* L) (p.eval₂ gO xO)) =
        p.eval₂ g x
    rw [Polynomial.hom_eval₂]
    rfl
  have hznonunit : ¬IsUnit zO := by
    intro hzunit
    have hzval : v (zO : L) = 1 :=
      ((Valuation.valuationSubring.integers v).isUnit_iff_valuation_eq_one).mp
        (by simpa [O] using hzunit)
    rw [hzcoe] at hzval
    exact (ne_of_lt heval) hzval
  have hzres : IsLocalRing.residue O zO = 0 := by
    rw [IsLocalRing.residue_eq_zero_iff, IsLocalRing.mem_maximalIdeal]
    exact hznonunit
  let q := p * Polynomial.C p.leadingCoeff⁻¹
  have hq : q.Monic :=
    Polynomial.monic_mul_leadingCoeff_inv hp
  have hpbar :
      p.eval₂
        ((IsLocalRing.residue O).comp gO)
        (IsLocalRing.residue O xO) = 0 := by
    rw [← Polynomial.hom_eval₂]
    exact hzres
  have hqbar :
      q.eval₂
        ((IsLocalRing.residue O).comp gO)
        (IsLocalRing.residue O xO) = 0 := by
    rw [show q = p * Polynomial.C p.leadingCoeff⁻¹ by rfl,
      Polynomial.eval₂_mul, hpbar, zero_mul]
  apply exists_complex_approximation_of_monic_residue_root
    v g hcoeff x hx.le q hq
  simpa [O, gO, xO] using hqbar

/--
An integral polynomial cannot have value below one at a point strictly
residue-equivalent to a nonroot.
-/
theorem complex_isRoot_of_close_polynomial_eval_lt_one
    {L Γ : Type*} [Field L] [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation L Γ) (g : ℂ →+* L)
    (hunit : ∀ z : ℂ, z ≠ 0 → v (g z) = 1)
    (p : Polynomial ℂ) (x : L) (c : ℂ)
    (hclose : v (x - g c) < 1)
    (heval : v (p.eval₂ g x) < 1) :
    p.IsRoot c := by
  rw [Polynomial.IsRoot]
  by_contra hpc
  have hcoeff :
      ∀ i ∈ (p.map g).support,
        v ((p.map g).coeff i) ≤ 1 := by
    intro i hi
    have hmapne : g (p.coeff i) ≠ 0 := by
      simpa using Polynomial.mem_support_iff.mp hi
    have hpne : p.coeff i ≠ 0 := by
      intro hpzero
      exact hmapne (by simp [hpzero])
    simpa using (hunit (p.coeff i) hpne).le
  have hc : v (g c) ≤ 1 := by
    by_cases hc0 : c = 0
    · simp [hc0]
    · exact (hunit c hc0).le
  have hcenter :
      v ((p.map g).eval (g c)) = 1 := by
    rw [Polynomial.eval_map_apply]
    exact hunit (p.eval c) hpc
  have hxval :=
    valuation_polynomial_eval_eq_one_of_sub_lt_one
      v (p.map g) (g c) x hcoeff hc hcenter hclose
  rw [Polynomial.eval_map] at hxval
  exact (ne_of_lt heval) hxval

/--
Near `c`, the valuation of a nonzero complex polynomial is exactly the
valuation of `x-c` raised to the root multiplicity at `c`.
-/
theorem valuation_polynomial_eval₂_eq_pow_rootMultiplicity
    {L Γ : Type*} [Field L] [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation L Γ) (g : ℂ →+* L)
    (hunit : ∀ z : ℂ, z ≠ 0 → v (g z) = 1)
    (p : Polynomial ℂ) (hp : p ≠ 0)
    (x : L) (c : ℂ)
    (hclose : v (x - g c) < 1) :
    v (p.eval₂ g x) =
      v (x - g c) ^ p.rootMultiplicity c := by
  let q :=
    p /ₘ
      (Polynomial.X - Polynomial.C c) ^ p.rootMultiplicity c
  have hqnonzero : q.eval c ≠ 0 := by
    simpa [q] using
      Polynomial.eval_divByMonic_pow_rootMultiplicity_ne_zero c hp
  have hqcoeff :
      ∀ i ∈ (q.map g).support,
        v ((q.map g).coeff i) ≤ 1 := by
    intro i hi
    have hmapne : g (q.coeff i) ≠ 0 := by
      simpa using Polynomial.mem_support_iff.mp hi
    have hcoeffne : q.coeff i ≠ 0 := by
      intro hzero
      exact hmapne (by simp [hzero])
    simpa using (hunit (q.coeff i) hcoeffne).le
  have hc : v (g c) ≤ 1 := by
    by_cases hc0 : c = 0
    · simp [hc0]
    · exact (hunit c hc0).le
  have hqcenter :
      v ((q.map g).eval (g c)) = 1 := by
    rw [Polynomial.eval_map_apply]
    exact hunit (q.eval c) hqnonzero
  have hqval :
      v (q.eval₂ g x) = 1 := by
    have h :=
      valuation_polynomial_eval_eq_one_of_sub_lt_one
        v (q.map g) (g c) x hqcoeff hc hqcenter hclose
    simpa [Polynomial.eval_map] using h
  have hfactor :
      (Polynomial.X - Polynomial.C c) ^ p.rootMultiplicity c * q = p := by
    simpa [q] using
      Polynomial.pow_mul_divByMonic_rootMultiplicity_eq p c
  calc
    v (p.eval₂ g x) =
        v (((Polynomial.X - Polynomial.C c) ^
          p.rootMultiplicity c * q).eval₂ g x) := by rw [hfactor]
    _ = v (x - g c) ^ p.rootMultiplicity c := by
      rw [Polynomial.eval₂_mul, map_mul, Polynomial.eval₂_pow, map_pow,
        Polynomial.eval₂_sub, Polynomial.eval₂_X, Polynomial.eval₂_C,
        hqval, mul_one]

/--
At a nonzero subunit value of a complex polynomial on a valuation unit, the
derivative has strictly larger valuation.
-/
theorem valuation_polynomial_derivative_gt_of_eval_lt_one
    {L Γ : Type*} [Field L] [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation L Γ) (g : ℂ →+* L)
    (hunit : ∀ z : ℂ, z ≠ 0 → v (g z) = 1)
    (p : Polynomial ℂ) (hp : p ≠ 0)
    (x : L) (hx : v x = 1)
    (hevalne : p.eval₂ g x ≠ 0)
    (heval : v (p.eval₂ g x) < 1) :
    v (p.eval₂ g x) < v (p.derivative.eval₂ g x) := by
  obtain ⟨c, hclose⟩ :=
    exists_complex_approximation_of_polynomial_eval_lt_one
      v g hunit p hp x hx heval
  have hroot :
      p.IsRoot c :=
    complex_isRoot_of_close_polynomial_eval_lt_one
      v g hunit p x c hclose heval
  have hderiv : p.derivative ≠ 0 := by
    intro hd
    have peq :
        p = Polynomial.C (p.coeff 0) :=
      Polynomial.eq_C_of_derivative_eq_zero hd
    have hc0 : p.coeff 0 = 0 := by
      rw [Polynomial.IsRoot, peq, Polynomial.eval_C] at hroot
      exact hroot
    exact hp (by rw [peq, hc0, Polynomial.C_0])
  have hpval :=
    valuation_polynomial_eval₂_eq_pow_rootMultiplicity
      v g hunit p hp x c hclose
  have hdval :=
    valuation_polynomial_eval₂_eq_pow_rootMultiplicity
      v g hunit p.derivative hderiv x c hclose
  have hdeltaNe : x - g c ≠ 0 := by
    intro hdelta
    have hxgc : x = g c := sub_eq_zero.mp hdelta
    apply hevalne
    rw [hxgc, Polynomial.eval₂_at_apply, hroot.eq_zero, map_zero]
  have hdeltapos : 0 < v (x - g c) :=
    (v.pos_iff).2 hdeltaNe
  have hmpos : 0 < p.rootMultiplicity c :=
    (Polynomial.rootMultiplicity_pos hp).2 hroot
  rw [hpval, hdval,
    Polynomial.derivative_rootMultiplicity_of_root hroot]
  exact pow_lt_pow_right_of_lt_one₀ hdeltapos hclose
    (Nat.sub_lt hmpos zero_lt_one)

/--
If one Laurent monomial has strictly larger valuation than every other
supported monomial, it determines the valuation of the evaluation.
-/
theorem valuation_laurent_eval₂_eq_of_unique_max
    {R S Γ : Type*} [CommRing R] [Field S]
    [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation S Γ) (g : R →+* S) (x : Sˣ)
    (f : LaurentPolynomial R) (j : ℤ)
    (hj : j ∈ f.support)
    (hmax : ∀ i ∈ f.support \ {j},
      v (g (f i) * (x ^ i).val) <
        v (g (f j) * (x ^ j).val)) :
    v (LaurentPolynomial.eval₂ g x f) =
      v (g (f j) * (x ^ j).val) := by
  rw [laurent_eval₂_eq_sum]
  exact v.map_sum_eq_of_lt hj hmax

theorem valuation_laurent_monomial
    {R S Γ : Type*} [CommRing R] [Field S]
    [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation S Γ) (g : R →+* S) (x : Sˣ)
    (f : LaurentPolynomial R) (i : ℤ)
    (hcoeff : v (g (f i)) = 1) :
    v (g (f i) * (x ^ i).val) = v x.val ^ i := by
  rw [map_mul, hcoeff, one_mul]
  rw [Units.val_zpow_eq_zpow_val]
  exact map_zpow₀ v x.val i

/--
For a two-sided Laurent polynomial whose supported coefficients are valuation
units, an evaluation of value below one can only occur at an argument of value
one.
-/
theorem valuation_laurent_eval₂_eq_one_of_two_sided
    {R S Γ : Type*} [CommRing R] [Field S]
    [LinearOrderedCommGroupWithZero Γ]
    (v : Valuation S Γ) (g : R →+* S) (x : Sˣ)
    (f : LaurentPolynomial R)
    (hnegative : ∃ i ∈ f.support, i < 0)
    (hpositive : ∃ i ∈ f.support, 0 < i)
    (hcoeff : ∀ i ∈ f.support, v (g (f i)) = 1)
    (heval : v (LaurentPolynomial.eval₂ g x f) < 1) :
    v x.val = 1 := by
  have hsupport : f.support.Nonempty :=
    ⟨hnegative.choose, hnegative.choose_spec.1⟩
  have hxne : x.val ≠ 0 := Units.ne_zero x
  have hxvalpos : 0 < v x.val := (v.pos_iff).2 hxne
  apply le_antisymm
  · apply not_lt.mp
    intro hxgt
    let hi := f.support.max' hsupport
    have hiMem : hi ∈ f.support :=
      Finset.max'_mem f.support hsupport
    have hiPos : 0 < hi := by
      exact lt_of_lt_of_le hpositive.choose_spec.2
        (Finset.le_max' f.support hpositive.choose
          hpositive.choose_spec.1)
    have hmax : ∀ i ∈ f.support \ {hi},
        v (g (f i) * (x ^ i).val) <
          v (g (f hi) * (x ^ hi).val) := by
      intro i hiDiff
      have hiSupport : i ∈ f.support :=
        (Finset.mem_sdiff.mp hiDiff).1
      have hine : i ≠ hi := by
        simpa using (Finset.mem_sdiff.mp hiDiff).2
      have hiLt : i < hi :=
        lt_of_le_of_ne
          (Finset.le_max' f.support i hiSupport) hine
      rw [valuation_laurent_monomial v g x f i
          (hcoeff i hiSupport),
        valuation_laurent_monomial v g x f hi
          (hcoeff hi hiMem)]
      exact zpow_lt_zpow_right₀ hxgt hiLt
    have hevalEq :=
      valuation_laurent_eval₂_eq_of_unique_max
        v g x f hi hiMem hmax
    have hiVal :
        1 < v (g (f hi) * (x ^ hi).val) := by
      rw [valuation_laurent_monomial v g x f hi
        (hcoeff hi hiMem)]
      simpa only [zpow_zero] using
        zpow_lt_zpow_right₀ hxgt hiPos
    rw [hevalEq] at heval
    exact (not_lt_of_ge hiVal.le) heval
  · apply not_lt.mp
    intro hxlt
    let lo := f.support.min' hsupport
    have hloMem : lo ∈ f.support :=
      Finset.min'_mem f.support hsupport
    have hloNeg : lo < 0 := by
      exact lt_of_le_of_lt
        (Finset.min'_le f.support hnegative.choose
          hnegative.choose_spec.1)
        hnegative.choose_spec.2
    have hmax : ∀ i ∈ f.support \ {lo},
        v (g (f i) * (x ^ i).val) <
          v (g (f lo) * (x ^ lo).val) := by
      intro i hi
      have hiMem : i ∈ f.support := (Finset.mem_sdiff.mp hi).1
      have hine : i ≠ lo := by
        simpa using (Finset.mem_sdiff.mp hi).2
      have hloLt : lo < i :=
        lt_of_le_of_ne
          (Finset.min'_le f.support i hiMem) (Ne.symm hine)
      rw [valuation_laurent_monomial v g x f i
          (hcoeff i hiMem),
        valuation_laurent_monomial v g x f lo
          (hcoeff lo hloMem)]
      exact zpow_lt_zpow_right_of_lt_one₀ hxvalpos hxlt hloLt
    have hevalEq :=
      valuation_laurent_eval₂_eq_of_unique_max
        v g x f lo hloMem hmax
    have hloVal :
        1 < v (g (f lo) * (x ^ lo).val) := by
      rw [valuation_laurent_monomial v g x f lo
        (hcoeff lo hloMem)]
      simpa only [zpow_zero] using
        zpow_lt_zpow_right_of_lt_one₀ hxvalpos hxlt hloNeg
    rw [hevalEq] at heval
    exact (not_lt_of_ge hloVal.le) heval

theorem complex_laurentSeries_valuation_algebraMap_eq_one
    {c : ℂ} (hc : c ≠ 0) :
    (Valued.v : Valuation (LaurentSeries ℂ) ℤᵐ⁰)
        (algebraMap ℂ (LaurentSeries ℂ) c) = 1 := by
  rw [show
    algebraMap ℂ (LaurentSeries ℂ) c =
      algebraMap (PowerSeries ℂ) (LaurentSeries ℂ)
        (algebraMap ℂ (PowerSeries ℂ) c) by
      rw [LaurentSeries.algebraMap_apply]
      simp]
  rw [LaurentSeries.valuation_def,
    IsDedekindDomain.HeightOneSpectrum.valuation_of_algebraMap,
    IsDedekindDomain.HeightOneSpectrum.intValuation_eq_one_iff]
  intro hmem
  rw [PowerSeries.idealX, Ideal.mem_span_singleton,
    PowerSeries.X_dvd_iff] at hmem
  exact hc (by simpa using hmem)

/--
Specialization of the two-sided valuation lemma to the `t`-adic valuation on
complex Laurent series and any finite-extension valuation lying above it.
-/
theorem complex_laurentSeries_root_valuation_eq_one
    {L : Type*} [Field L] [Algebra (LaurentSeries ℂ) L]
    (vL : Valuation L ℤᵐ⁰)
    [(Valued.v : Valuation (LaurentSeries ℂ) ℤᵐ⁰).HasExtension vL]
    (f : LaurentPolynomial ℂ) (x : Lˣ)
    (hnegative : ∃ i ∈ f.support, i < 0)
    (hpositive : ∃ i ∈ f.support, 0 < i)
    (heval :
      vL (LaurentPolynomial.eval₂
        ((algebraMap (LaurentSeries ℂ) L).comp
          (algebraMap ℂ (LaurentSeries ℂ))) x f) < 1) :
    vL x.val = 1 := by
  apply valuation_laurent_eval₂_eq_one_of_two_sided
    vL ((algebraMap (LaurentSeries ℂ) L).comp
      (algebraMap ℂ (LaurentSeries ℂ))) x f hnegative hpositive
  · intro i hi
    have hcoeffne : f i ≠ 0 :=
      Finsupp.mem_support_iff.mp hi
    have hbase :=
      complex_laurentSeries_valuation_algebraMap_eq_one hcoeffne
    have hext :
        vL (algebraMap (LaurentSeries ℂ) L
          (algebraMap ℂ (LaurentSeries ℂ) (f i))) = 1 :=
      (Valuation.HasExtension.val_map_eq_one_iff
        (Valued.v : Valuation (LaurentSeries ℂ) ℤᵐ⁰) vL _).mpr hbase
    exact hext
  · exact heval

end

end MathieuQuadratic.DvKProof
