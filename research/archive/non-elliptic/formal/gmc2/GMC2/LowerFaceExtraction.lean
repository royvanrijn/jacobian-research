import GMC2.CircularFlatten
import GMC2.DuistermaatVanDerKallen
import GMC2.SupportingFace
import Mathlib.Algebra.Polynomial.Degree.TrailingDegree
import Mathlib.Tactic

/-!
# Extraction of the exposed lowest radial coefficient

This file connects the rational supporting-face certificate to the
Duistermaat--van der Kallen theorem.  The weighted associated-graded
calculation is supplied by `WeightedInitial`.
-/

namespace GMC2

open Polynomial LaurentPolynomial

/-- Radial order of the coefficient at angular weight `k`. -/
noncomputable def radialOrder {K : Type*} [CommRing K]
    (P : CircularPolynomial K) (k : ℤ) : ℕ :=
  (P.coeff k).natTrailingDegree

private theorem certificate_lower_weight
    {K : Type*} [CommRing K]
    (P : CircularPolynomial K)
    (C : LowerFaceCertificate (angularSupport P) (radialOrder P)) :
    HasLowerBidegreeWeight C.slope C.intercept
      (circularBigradedEquiv K P) := by
  intro x hx
  have hcoeff :
      (P.coeff x.1).coeff x.2 ≠ 0 := by
    rw [← circularBigradedEquiv_coeff P x.1 x.2]
    simpa only [Prod.eta] using (Finsupp.mem_support_iff.mp hx)
  have hkpoly : P.coeff x.1 ≠ 0 := by
    intro h
    rw [h] at hcoeff
    exact hcoeff (by simp)
  have hk : x.1 ∈ angularSupport P := by
    exact Finsupp.mem_support_iff.mpr hkpoly
  have hface := C.lower_bound x.1 hk
  have hradial : radialOrder P x.1 ≤ x.2 := by
    by_contra h
    have hj : x.2 < radialOrder P x.1 := Nat.lt_of_not_ge h
    exact hcoeff (Polynomial.coeff_eq_zero_of_lt_natTrailingDegree hj)
  have hradialQ : (radialOrder P x.1 : ℚ) ≤ (x.2 : ℚ) := by
    exact_mod_cast hradial
  dsimp [bidegreeWeight]
  nlinarith

private noncomputable def faceInitial
    {K : Type*} [CommRing K]
    (P : CircularPolynomial K)
    (C : LowerFaceCertificate (angularSupport P) (radialOrder P)) :
    Bigraded K :=
  weightedInitial C.slope C.intercept (circularBigradedEquiv K P)

/-- The angular lower-face polynomial obtained by forgetting radial degree
in the exact initial form. -/
noncomputable def lowerFacePolynomial
    {K : Type*} [CommRing K]
    (P : CircularPolynomial K)
    (C : LowerFaceCertificate (angularSupport P) (radialOrder P)) :
    LaurentPolynomial K :=
  forgetRadial K (faceInitial P C)

private theorem faceInitial_coeff_ne_zero
    {K : Type*} [CommRing K]
    (P : CircularPolynomial K)
    (C : LowerFaceCertificate (angularSupport P) (radialOrder P))
    {k : ℤ} (hk : k ∈ C.face) :
    (faceInitial P C).coeff (k, radialOrder P k) ≠ 0 := by
  have hkS := C.face_subset hk
  have hpoly : P.coeff k ≠ 0 :=
    Finsupp.mem_support_iff.mp hkS
  have htrail : (P.coeff k).coeff (radialOrder P k) ≠ 0 := by
    exact Polynomial.coeff_natTrailingDegree_ne_zero.mpr hpoly
  rw [faceInitial, weightedInitial_coeff, if_pos]
  · simpa using htrail
  · have heq := C.eq_on_face k hk
    simp [bidegreeWeight, radialOrder] at heq ⊢
    linarith

private theorem faceInitial_unique_radial
    {K : Type*} [CommRing K]
    (P : CircularPolynomial K)
    (C : LowerFaceCertificate (angularSupport P) (radialOrder P))
    {k : ℤ} {x : ℤ × ℕ}
    (hx : x ∈ (faceInitial P C).coeff.support) (hxk : x.1 = k) :
    x = (k, radialOrder P k) := by
  have hinit := Finsupp.mem_support_iff.mp hx
  have hweight : bidegreeWeight C.slope x = C.intercept := by
    by_contra h
    rw [faceInitial, weightedInitial_coeff, if_neg h] at hinit
    exact hinit rfl
  have hflat :
      (circularBigradedEquiv K P).coeff x ≠ 0 := by
    rw [faceInitial, weightedInitial_coeff, if_pos hweight] at hinit
    exact hinit
  have hkpoly : P.coeff k ≠ 0 := by
    intro h
    apply hflat
    rw [← Prod.eta x, hxk, circularBigradedEquiv_coeff, h]
    simp
  have hkS : k ∈ angularSupport P :=
    Finsupp.mem_support_iff.mpr hkpoly
  have hlower := C.lower_bound k hkS
  have hcoeff :
      (P.coeff k).coeff x.2 ≠ 0 := by
    rw [← Prod.eta x, hxk] at hflat
    rw [circularBigradedEquiv_coeff] at hflat
    exact hflat
  have hradial : radialOrder P k ≤ x.2 := by
    by_contra h
    exact hcoeff (Polynomial.coeff_eq_zero_of_lt_natTrailingDegree
      (Nat.lt_of_not_ge h))
  have hradialQ : (radialOrder P k : ℚ) ≤ (x.2 : ℚ) := by
    exact_mod_cast hradial
  have hxweight : C.intercept + C.slope * (k : ℚ) = (x.2 : ℚ) := by
    simp [bidegreeWeight, hxk] at hweight
    linarith
  have hdegree : x.2 = radialOrder P k := by
    exact_mod_cast (le_antisymm (by linarith) hradialQ)
  exact Prod.ext hxk hdegree

theorem lowerFacePolynomial_coeff_ne_zero
    {K : Type*} [CommRing K]
    (P : CircularPolynomial K)
    (C : LowerFaceCertificate (angularSupport P) (radialOrder P))
    {k : ℤ} (hk : k ∈ C.face) :
    (lowerFacePolynomial P C).coeff k ≠ 0 := by
  classical
  let x : ℤ × ℕ := (k, radialOrder P k)
  have hx : x ∈ (faceInitial P C).coeff.support :=
    Finsupp.mem_support_iff.mpr (faceInitial_coeff_ne_zero P C hk)
  rw [lowerFacePolynomial, forgetRadial_coeff]
  rw [Finset.sum_eq_single x]
  · exact faceInitial_coeff_ne_zero P C hk
  · intro y hy hyne
    have hySupport : y ∈ (faceInitial P C).coeff.support :=
      (Finset.mem_filter.mp hy).1
    have hyk : y.1 = k := (Finset.mem_filter.mp hy).2
    exact (hyne (faceInitial_unique_radial P C hySupport hyk)).elim
  · intro h
    exact (h (by simpa [x] using hx)).elim

theorem lowerFacePolynomial_straddles
    {K : Type*} [CommRing K]
    (P : CircularPolynomial K)
    (C : LowerFaceCertificate (angularSupport P) (radialOrder P)) :
    SupportStraddlesZero (lowerFacePolynomial P C) := by
  rcases C.has_nonpos with ⟨a, ha, ha0⟩
  rcases C.has_nonneg with ⟨b, hb, hb0⟩
  exact
    ⟨⟨a, Finsupp.mem_support_iff.mpr
        (lowerFacePolynomial_coeff_ne_zero P C ha), ha0⟩,
      ⟨b, Finsupp.mem_support_iff.mpr
        (lowerFacePolynomial_coeff_ne_zero P C hb), hb0⟩⟩

/-- Data extracted from a straddling support: a positive power whose
angular constant term has a genuine lowest radial coefficient. -/
structure ExposedLowestTerm
    {K : Type*} [CommRing K] (P : CircularPolynomial K) where
  r : ℕ
  r_pos : 0 < r
  n : ℕ
  coeff_ne_zero : (constantTerm (P ^ r)).coeff n ≠ 0
  coeff_below_eq_zero :
    ∀ j < n, (constantTerm (P ^ r)).coeff j = 0
  scaled_coeff_below_eq_zero :
    ∀ q j, j < n * q → (constantTerm (P ^ (r * q))).coeff j = 0

/-- The lower-face lemma of the paper, including the DvdK step and the
associated-graded identification of the exposed coefficient. -/
noncomputable def exposedLowestTerm_of_straddles
    {K : Type*} [CommRing K] [IsDomain K] [CharZero K]
    (P : CircularPolynomial K) (hconv : SupportStraddlesZero P) :
    ExposedLowestTerm P := by
  classical
  let C := rational_supportingFace (angularSupport P) (radialOrder P)
    hconv.1 hconv.2
  let P₀ := lowerFacePolynomial P C
  have hP₀conv : SupportStraddlesZero P₀ :=
    lowerFacePolynomial_straddles P C
  have hdvk := duistermaat_van_der_kallen_domain P₀ hP₀conv
  let r := Classical.choose hdvk
  have hrSpec := Classical.choose_spec hdvk
  have hr : 0 < r := hrSpec.1
  have hct : (P₀ ^ r).coeff 0 ≠ 0 := hrSpec.2
  let f := circularBigradedEquiv K P
  let g := weightedInitial C.slope ((r : ℚ) * C.intercept) (f ^ r)
  have hpow :
      g = faceInitial P C ^ r := by
    exact weightedInitial_pow C.slope C.intercept f
      (certificate_lower_weight P C) r
  have hproject : (forgetRadial K g).coeff 0 ≠ 0 := by
    rw [hpow, map_pow]
    exact hct
  have hsum :
      (∑ x ∈ g.coeff.support with x.1 = 0, g.coeff x) ≠ 0 := by
    simpa [forgetRadial_coeff] using hproject
  have hxExists := Finset.exists_ne_zero_of_sum_ne_zero hsum
  let x := Classical.choose hxExists
  have hxSpec := Classical.choose_spec hxExists
  have hx : x ∈ g.coeff.support.filter (fun x ↦ x.1 = 0) := hxSpec.1
  have hxcoeff : g.coeff x ≠ 0 := hxSpec.2
  have hxSupport : x ∈ g.coeff.support := (Finset.mem_filter.mp hx).1
  have hxzero : x.1 = 0 := (Finset.mem_filter.mp hx).2
  let n := x.2
  have hxweight :
      bidegreeWeight C.slope x = (r : ℚ) * C.intercept := by
    have := Finsupp.mem_support_iff.mp hxSupport
    by_contra h
    change
      (weightedInitial C.slope ((r : ℚ) * C.intercept) (f ^ r)).coeff x ≠ 0
      at this
    rw [weightedInitial_coeff, if_neg h] at this
    exact this rfl
  have hnweight : (n : ℚ) = (r : ℚ) * C.intercept := by
    simpa [n, bidegreeWeight, hxzero] using hxweight
  have hprojectCoeff : (forgetRadial K g).coeff 0 = g.coeff (0, n) := by
    rw [forgetRadial_coeff]
    rw [Finset.sum_eq_single x]
    · congr 1
      exact Prod.ext hxzero (by rfl)
    · intro y hy hyne
      have hySupport : y ∈ g.coeff.support := (Finset.mem_filter.mp hy).1
      have hyzero : y.1 = 0 := (Finset.mem_filter.mp hy).2
      have hyweight :
          bidegreeWeight C.slope y = (r : ℚ) * C.intercept := by
        have := Finsupp.mem_support_iff.mp hySupport
        by_contra h
        change
          (weightedInitial C.slope ((r : ℚ) * C.intercept) (f ^ r)).coeff y ≠ 0
          at this
        rw [weightedInitial_coeff, if_neg h] at this
        exact this rfl
      have hyn : y.2 = n := by
        have : (y.2 : ℚ) = (n : ℚ) := by
          simp [bidegreeWeight, hyzero] at hyweight
          linarith
        exact_mod_cast this
      exact (hyne (Prod.ext (hyzero.trans hxzero.symm)
        (by simpa [n] using hyn))).elim
    · intro h
      exact (h (by simpa [n, hxzero] using hxSupport)).elim
  have hgnonzero : g.coeff (0, n) ≠ 0 := by
    rw [← hprojectCoeff]
    exact hproject
  have hweight0 :
      bidegreeWeight C.slope (0, n) = (r : ℚ) * C.intercept := by
    simp [bidegreeWeight, hnweight]
  have hcoeff :
      (constantTerm (P ^ r)).coeff n ≠ 0 := by
    have hgflat :
        g.coeff (0, n) = (f ^ r).coeff (0, n) := by
      change
        (weightedInitial C.slope ((r : ℚ) * C.intercept) (f ^ r)).coeff (0, n) =
          (f ^ r).coeff (0, n)
      rw [weightedInitial_coeff, if_pos hweight0]
    rw [hgflat] at hgnonzero
    have hmap : f ^ r = circularBigradedEquiv K (P ^ r) := by
      simp [f]
    rw [hmap, circularBigradedEquiv_coeff] at hgnonzero
    simpa using hgnonzero
  refine
    { r := r
      r_pos := hr
      n := n
      coeff_ne_zero := hcoeff
      coeff_below_eq_zero := ?_
      scaled_coeff_below_eq_zero := ?_ }
  · intro j hj
    by_contra hjcoeff
    have hsupport :
        (0, j) ∈ (f ^ r).coeff.support := by
      apply Finsupp.mem_support_iff.mpr
      have hmap : f ^ r = circularBigradedEquiv K (P ^ r) := by
        simp [f]
      rw [hmap, circularBigradedEquiv_coeff]
      simpa using hjcoeff
    have hlower :=
      hasLowerBidegreeWeight_pow C.slope C.intercept f
        (certificate_lower_weight P C) r (0, j) hsupport
    have hlower' : (r : ℚ) * C.intercept ≤ (j : ℚ) := by
      simpa only [bidegreeWeight, Int.cast_zero, mul_zero, sub_zero] using hlower
    have hjQ : (j : ℚ) < (n : ℚ) := by exact_mod_cast hj
    linarith
  · intro q j hj
    by_contra hjcoeff
    have hsupport :
        (0, j) ∈ (f ^ (r * q)).coeff.support := by
      apply Finsupp.mem_support_iff.mpr
      have hmap :
          f ^ (r * q) = circularBigradedEquiv K (P ^ (r * q)) := by
        simp [f]
      rw [hmap, circularBigradedEquiv_coeff]
      simpa using hjcoeff
    have hlower :=
      hasLowerBidegreeWeight_pow C.slope C.intercept f
        (certificate_lower_weight P C) (r * q) (0, j) hsupport
    have hlower' :
        ((r * q : ℕ) : ℚ) * C.intercept ≤ (j : ℚ) := by
      simpa only [bidegreeWeight, Int.cast_zero, mul_zero, sub_zero] using hlower
    have hjQ : (j : ℚ) < ((n * q : ℕ) : ℚ) := by exact_mod_cast hj
    push_cast at hlower' hjQ
    nlinarith

end GMC2
