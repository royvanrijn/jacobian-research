/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GVC.Definitions
import Mathlib.Combinatorics.Hall.Basic
import Mathlib.Tactic.FieldSimp
import Mathlib.Tactic.FinCases
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Order
import Mathlib.Tactic.Push
import Mathlib.Tactic.Ring

/-!
# The finite Hall-localization core

This file formalizes the combinatorial heart of Lemma 3.1 in the GVC
manuscript.  Derivative copies are indexed by `Fin r`, polynomial factors by
`Fin d`, and `adj i j` says that derivative direction `i` acts nontrivially
on polynomial factor `j`.

The translated Duistermaat--van der Kallen argument is responsible for the
absence of a matching.  From that exact input, Hall's theorem supplies a
deficient set.  In two variables, any two nonparallel derivative directions
together see every nonzero linear factor.  Since `d > r`, a deficient set
therefore lies in one direction class.  Counting the complementary factors
gives the exponent `d - e + 1` in the Hall normal form.
-/

namespace GVC

/-- Polynomial factors on which one derivative copy acts nontrivially. -/
def hallNeighbors {r d : ℕ} (adj : Fin r → Fin d → Prop)
    [DecidableRel adj] (i : Fin r) : Finset (Fin d) :=
  Finset.univ.filter fun j => adj i j

/-- The union of the neighbors of a finite set of derivative copies. -/
def hallNeighborhood {r d : ℕ} (adj : Fin r → Fin d → Prop)
    [DecidableRel adj] (S : Finset (Fin r)) : Finset (Fin d) :=
  S.biUnion (hallNeighbors adj)

/-- The full multiplicity class of one derivative direction. -/
def hallDirectionClass {r : ℕ} (parallel : Fin r → Fin r → Prop)
    [DecidableRel parallel] (i : Fin r) : Finset (Fin r) :=
  Finset.univ.filter fun k => parallel i k

/-- Polynomial factors annihilated by one derivative direction. -/
def hallAnnihilators {r d : ℕ} (adj : Fin r → Fin d → Prop)
    [DecidableRel adj] (i : Fin r) : Finset (Fin d) :=
  Finset.univ.filter fun j => ¬adj i j

/-- The exact finite Hall-localization statement used in Lemma 3.1.

`nonparallel_cover` is the specifically binary input: two nonparallel
directions cannot both be annihilated by a nonzero binary linear form.  The
conclusion chooses a full direction class of multiplicity `e`; at least
`d - e + 1` polynomial factors annihilate that direction. -/
theorem binaryHall_localization_core
    {r d : ℕ}
    (adj : Fin r → Fin d → Prop) [DecidableRel adj]
    (parallel : Fin r → Fin r → Prop) [DecidableRel parallel]
    (parallel_refl : ∀ i, parallel i i)
    (nonparallel_cover :
      ∀ i k, ¬parallel i k → ∀ j, adj i j ∨ adj k j)
    (hdegree : r < d)
    (hnoMatching : ¬∃ f : Fin r → Fin d,
      Function.Injective f ∧ ∀ i, adj i (f i)) :
    ∃ i : Fin r,
      1 ≤ (hallDirectionClass parallel i).card ∧
      (hallDirectionClass parallel i).card ≤ r ∧
      d - (hallDirectionClass parallel i).card + 1 ≤
        (hallAnnihilators adj i).card := by
  classical
  have hnotHall :
      ¬∀ S : Finset (Fin r), S.card ≤ (hallNeighborhood adj S).card := by
    intro hHall
    apply hnoMatching
    obtain ⟨f, hfInjective, hfMem⟩ :=
      (Finset.all_card_le_biUnion_card_iff_exists_injective
        (hallNeighbors adj)).mp (by
          intro S
          simpa [hallNeighborhood] using hHall S)
    refine ⟨f, hfInjective, ?_⟩
    intro i
    simpa [hallNeighbors] using hfMem i
  push Not at hnotHall
  obtain ⟨S, hdeficient⟩ := hnotHall
  have hSNonempty : S.Nonempty := by
    rw [Finset.nonempty_iff_ne_empty]
    intro hS
    subst S
    simp [hallNeighborhood] at hdeficient
  obtain ⟨i, hiS⟩ := hSNonempty
  have hparallel : ∀ k ∈ S, parallel i k := by
    intro k hkS
    by_contra hik
    have hfull : hallNeighborhood adj S = Finset.univ := by
      apply Finset.eq_univ_of_forall
      intro j
      rcases nonparallel_cover i k hik j with hij | hkj
      · exact Finset.mem_biUnion.mpr
          ⟨i, hiS, by simp [hallNeighbors, hij]⟩
      · exact Finset.mem_biUnion.mpr
          ⟨k, hkS, by simp [hallNeighbors, hkj]⟩
    have hScard : S.card ≤ r := by simpa using Finset.card_le_univ S
    have hNcard : (hallNeighborhood adj S).card = d := by
      rw [hfull]
      simp
    omega
  let E := hallDirectionClass parallel i
  let N := hallNeighbors adj i
  let A := hallAnnihilators adj i
  have hiE : i ∈ E := by
    simp [E, hallDirectionClass, parallel_refl]
  have hSE : S ⊆ E := by
    intro k hkS
    simp [E, hallDirectionClass, hparallel k hkS]
  have hNE : N ⊆ hallNeighborhood adj S := by
    intro j hjN
    exact Finset.mem_biUnion.mpr ⟨i, hiS, hjN⟩
  have hNltS : N.card < S.card :=
    lt_of_le_of_lt (Finset.card_le_card hNE) hdeficient
  have hSleE : S.card ≤ E.card := Finset.card_le_card hSE
  have hEle : E.card ≤ r := by simpa [E] using Finset.card_le_univ E
  have hpartition : N.card + A.card = d := by
    simpa [N, A, hallNeighbors, hallAnnihilators] using
      (Finset.card_filter_add_card_filter_not
        (s := (Finset.univ : Finset (Fin d))) (fun j => adj i j))
  refine ⟨i, ?_, ?_, ?_⟩
  · exact Finset.card_pos.mpr ⟨i, hiE⟩
  · simpa [E] using hEle
  · change d - E.card + 1 ≤ A.card
    omega

/-- Evaluation of a binary linear form, represented by its two
coefficients, on a binary direction. -/
def binaryLinearEvaluation {K : Type*} [Semiring K]
    (ell v : Fin 2 → K) : K :=
  ell 0 * v 0 + ell 1 * v 1

/-- The determinant detecting parallel binary directions. -/
def binaryDirectionDet {K : Type*} [CommRing K]
    (u v : Fin 2 → K) : K :=
  u 0 * v 1 - u 1 * v 0

def BinaryDirectionsParallel {K : Type*} [CommRing K]
    (u v : Fin 2 → K) : Prop :=
  binaryDirectionDet u v = 0

/-- A binary linear form as a polynomial in the two coordinate variables. -/
noncomputable def binaryLinearPolynomial {K : Type*} [CommSemiring K]
    (a : Fin 2 → K) : MvPolynomial (Fin 2) K :=
  MvPolynomial.C (a 0) * MvPolynomial.X 0 +
    MvPolynomial.C (a 1) * MvPolynomial.X 1

@[simp] theorem binaryDirectionsParallel_refl
    {K : Type*} [CommRing K] (v : Fin 2 → K) :
    BinaryDirectionsParallel v v := by
  simp only [BinaryDirectionsParallel, binaryDirectionDet]
  ring

/-- A nonzero binary linear form cannot annihilate two nonparallel
directions.  This is the two-dimensional linear-algebra input in the Hall
localization proof. -/
theorem nonparallel_binaryDirections_cover
    {K : Type*} [Field K]
    (ell u v : Fin 2 → K) (hell : ell ≠ 0)
    (huv : ¬BinaryDirectionsParallel u v) :
    binaryLinearEvaluation ell u ≠ 0 ∨
      binaryLinearEvaluation ell v ≠ 0 := by
  by_contra hvanish
  push Not at hvanish
  rcases hvanish with ⟨hu, hv⟩
  have hdet : binaryDirectionDet u v ≠ 0 := by
    simpa [BinaryDirectionsParallel] using huv
  have hell0mul : ell 0 * binaryDirectionDet u v = 0 := by
    dsimp [binaryLinearEvaluation] at hu hv
    dsimp [binaryDirectionDet]
    linear_combination (v 1) * hu - (u 1) * hv
  have hell1mul : ell 1 * binaryDirectionDet u v = 0 := by
    dsimp [binaryLinearEvaluation] at hu hv
    dsimp [binaryDirectionDet]
    linear_combination (u 0) * hv - (v 0) * hu
  have hell0 : ell 0 = 0 := (mul_eq_zero.mp hell0mul).resolve_right hdet
  have hell1 : ell 1 = 0 := (mul_eq_zero.mp hell1mul).resolve_right hdet
  apply hell
  funext j
  fin_cases j
  · exact hell0
  · exact hell1

/-- A fixed linear form annihilating `u`; after taking `u` as the first
coordinate it becomes the second coordinate form. -/
def binaryPerpendicular {K : Type*} [CommRing K]
    (u : Fin 2 → K) : Fin 2 → K :=
  fun i => if i = 0 then -u 1 else u 0

theorem binaryPerpendicular_ne_zero
    {K : Type*} [CommRing K] {u : Fin 2 → K} (hu : u ≠ 0) :
    binaryPerpendicular u ≠ 0 := by
  intro hperp
  apply hu
  funext i
  fin_cases i
  · have h := congrFun hperp 1
    simpa [binaryPerpendicular] using h
  · have h := congrFun hperp 0
    have hneg : -u 1 = 0 := by simpa [binaryPerpendicular] using h
    simpa using hneg

theorem binaryDirectionDet_perpendicular
    {K : Type*} [CommRing K] (ell u : Fin 2 → K) :
    binaryDirectionDet (binaryPerpendicular u) ell =
      -binaryLinearEvaluation ell u := by
  simp [binaryDirectionDet, binaryPerpendicular, binaryLinearEvaluation]
  ring

/-- Over a field, determinant-zero binary directions are scalar
multiples. -/
theorem eq_smul_of_binaryDirectionsParallel
    {K : Type*} [Field K] {u v : Fin 2 → K}
    (hu : u ≠ 0) (huv : BinaryDirectionsParallel u v) :
    ∃ c : K, v = c • u := by
  have hdet : u 0 * v 1 - u 1 * v 0 = 0 := by
    simpa [BinaryDirectionsParallel, binaryDirectionDet] using huv
  by_cases hu0 : u 0 = 0
  · have hu1 : u 1 ≠ 0 := by
      intro hu1
      apply hu
      funext i
      fin_cases i
      · exact hu0
      · exact hu1
    have hv0 : v 0 = 0 := by
      rw [hu0, zero_mul, zero_sub] at hdet
      have hmul : u 1 * v 0 = 0 := neg_eq_zero.mp hdet
      exact (mul_eq_zero.mp hmul).resolve_left hu1
    refine ⟨v 1 / u 1, ?_⟩
    funext i
    fin_cases i
    · simp [hu0, hv0]
    · simp [hu1]
  · have heq : u 0 * v 1 = u 1 * v 0 := sub_eq_zero.mp hdet
    refine ⟨v 0 / u 0, ?_⟩
    funext i
    fin_cases i
    · simp [hu0]
    · simp only [Pi.smul_apply, smul_eq_mul]
      field_simp
      simpa [mul_comm] using heq

/-- Every binary linear form annihilating a nonzero direction is a scalar
multiple of the fixed perpendicular form. -/
theorem eq_smul_binaryPerpendicular_of_evaluation_eq_zero
    {K : Type*} [Field K] {ell u : Fin 2 → K}
    (hu : u ≠ 0) (hvanish : binaryLinearEvaluation ell u = 0) :
    ∃ c : K, ell = c • binaryPerpendicular u := by
  apply eq_smul_of_binaryDirectionsParallel (binaryPerpendicular_ne_zero hu)
  rw [BinaryDirectionsParallel, binaryDirectionDet_perpendicular, hvanish]
  simp

theorem binaryLinearPolynomial_smul
    {K : Type*} [CommRing K] (c : K) (u : Fin 2 → K) :
    binaryLinearPolynomial (c • u) =
      MvPolynomial.C c * binaryLinearPolynomial u := by
  simp [binaryLinearPolynomial]
  ring

theorem binaryLinearPolynomial_dvd_of_parallel
    {K : Type*} [Field K] {u v : Fin 2 → K}
    (hu : u ≠ 0) (huv : BinaryDirectionsParallel u v) :
    binaryLinearPolynomial u ∣ binaryLinearPolynomial v := by
  obtain ⟨c, rfl⟩ := eq_smul_of_binaryDirectionsParallel hu huv
  refine ⟨MvPolynomial.C c, ?_⟩
  rw [binaryLinearPolynomial_smul]
  ring

theorem binaryLinearPolynomial_perpendicular_dvd_of_evaluation_eq_zero
    {K : Type*} [Field K] {ell u : Fin 2 → K}
    (hu : u ≠ 0) (hvanish : binaryLinearEvaluation ell u = 0) :
    binaryLinearPolynomial (binaryPerpendicular u) ∣
      binaryLinearPolynomial ell := by
  obtain ⟨c, rfl⟩ :=
    eq_smul_binaryPerpendicular_of_evaluation_eq_zero hu hvanish
  refine ⟨MvPolynomial.C c, ?_⟩
  rw [binaryLinearPolynomial_smul]
  ring

/-- If one element divides every factor in a finite product, the
corresponding cardinal power divides the product. -/
theorem pow_card_dvd_finset_prod_of_dvd
    {ι M : Type*} [CommMonoid M]
    (S : Finset ι) (a : M) (f : ι → M)
    (h : ∀ i ∈ S, a ∣ f i) :
    a ^ S.card ∣ ∏ i ∈ S, f i := by
  rw [← Finset.prod_const]
  exact Finset.prod_dvd_prod_of_dvd (fun _ => a) f h

/-- Split-factor incidence data after scalar extension.  The
Duistermaat--van der Kallen/polarization part of Lemma 3.1 supplies
`no_matching`; all subsequent Hall localization is proved below. -/
structure SplitBinaryHallSystem
    (K : Type*) [Field K] (r d : ℕ) where
  symbol : MvPolynomial (Fin 2) K
  polynomial : MvPolynomial (Fin 2) K
  symbolScalar : K
  polynomialScalar : K
  symbolScalar_ne_zero : symbolScalar ≠ 0
  polynomialScalar_ne_zero : polynomialScalar ≠ 0
  derivativeDirection : Fin r → Fin 2 → K
  polynomialFactor : Fin d → Fin 2 → K
  derivative_ne_zero : ∀ i, derivativeDirection i ≠ 0
  polynomialFactor_ne_zero : ∀ j, polynomialFactor j ≠ 0
  symbol_factorization :
    symbol = MvPolynomial.C symbolScalar *
      ∏ i, binaryLinearPolynomial (derivativeDirection i)
  polynomial_factorization :
    polynomial = MvPolynomial.C polynomialScalar *
      ∏ j, binaryLinearPolynomial (polynomialFactor j)
  no_matching : ¬∃ f : Fin r → Fin d,
    Function.Injective f ∧ ∀ i,
      binaryLinearEvaluation (polynomialFactor (f i))
        (derivativeDirection i) ≠ 0

noncomputable def SplitBinaryHallSystem.directionMultiplicity
    {K : Type*} [Field K] {r d : ℕ}
    (H : SplitBinaryHallSystem K r d) (i : Fin r) : ℕ := by
  classical
  exact (Finset.univ.filter fun k =>
    BinaryDirectionsParallel (H.derivativeDirection i)
      (H.derivativeDirection k)).card

noncomputable def SplitBinaryHallSystem.annihilatorCount
    {K : Type*} [Field K] {r d : ℕ}
    (H : SplitBinaryHallSystem K r d) (i : Fin r) : ℕ := by
  classical
  exact (Finset.univ.filter fun j =>
    binaryLinearEvaluation (H.polynomialFactor j)
      (H.derivativeDirection i) = 0).card

/-- The coordinate-free operator half of the Hall normal form: the full
parallel-direction power divides the split symbol. -/
theorem SplitBinaryHallSystem.directionPower_dvd_symbol
    {K : Type*} [Field K] {r d : ℕ}
    (H : SplitBinaryHallSystem K r d) (i : Fin r) :
    binaryLinearPolynomial (H.derivativeDirection i) ^
        H.directionMultiplicity i ∣ H.symbol := by
  classical
  let S := Finset.univ.filter fun k =>
    BinaryDirectionsParallel (H.derivativeDirection i)
      (H.derivativeDirection k)
  have hfactor : ∀ k ∈ S,
      binaryLinearPolynomial (H.derivativeDirection i) ∣
        binaryLinearPolynomial (H.derivativeDirection k) := by
    intro k hk
    exact binaryLinearPolynomial_dvd_of_parallel
      (H.derivative_ne_zero i) (Finset.mem_filter.mp hk).2
  have hselected :
      binaryLinearPolynomial (H.derivativeDirection i) ^ S.card ∣
        ∏ k ∈ S, binaryLinearPolynomial (H.derivativeDirection k) :=
    pow_card_dvd_finset_prod_of_dvd S
      (binaryLinearPolynomial (H.derivativeDirection i))
      (fun k => binaryLinearPolynomial (H.derivativeDirection k)) hfactor
  have hallFactors :
      binaryLinearPolynomial (H.derivativeDirection i) ^ S.card ∣
        ∏ k, binaryLinearPolynomial (H.derivativeDirection k) :=
    hselected.trans (Finset.prod_dvd_prod_of_subset S Finset.univ _
      (Finset.subset_univ S))
  rw [H.symbol_factorization]
  simpa [SplitBinaryHallSystem.directionMultiplicity, S] using
    dvd_mul_of_dvd_right hallFactors (MvPolynomial.C H.symbolScalar)

/-- The coordinate-free polynomial half of the Hall normal form: the power
of the fixed perpendicular form indexed by all annihilated factors divides
the split polynomial. -/
theorem SplitBinaryHallSystem.annihilatorPower_dvd_polynomial
    {K : Type*} [Field K] {r d : ℕ}
    (H : SplitBinaryHallSystem K r d) (i : Fin r) :
    binaryLinearPolynomial (binaryPerpendicular (H.derivativeDirection i)) ^
        H.annihilatorCount i ∣ H.polynomial := by
  classical
  let A := Finset.univ.filter fun j =>
    binaryLinearEvaluation (H.polynomialFactor j)
      (H.derivativeDirection i) = 0
  have hfactor : ∀ j ∈ A,
      binaryLinearPolynomial
          (binaryPerpendicular (H.derivativeDirection i)) ∣
        binaryLinearPolynomial (H.polynomialFactor j) := by
    intro j hj
    exact binaryLinearPolynomial_perpendicular_dvd_of_evaluation_eq_zero
      (H.derivative_ne_zero i) (Finset.mem_filter.mp hj).2
  have hselected :
      binaryLinearPolynomial
          (binaryPerpendicular (H.derivativeDirection i)) ^ A.card ∣
        ∏ j ∈ A, binaryLinearPolynomial (H.polynomialFactor j) :=
    pow_card_dvd_finset_prod_of_dvd A
      (binaryLinearPolynomial
        (binaryPerpendicular (H.derivativeDirection i)))
      (fun j => binaryLinearPolynomial (H.polynomialFactor j)) hfactor
  have hallFactors :
      binaryLinearPolynomial
          (binaryPerpendicular (H.derivativeDirection i)) ^ A.card ∣
        ∏ j, binaryLinearPolynomial (H.polynomialFactor j) :=
    hselected.trans (Finset.prod_dvd_prod_of_subset A Finset.univ _
      (Finset.subset_univ A))
  rw [H.polynomial_factorization]
  simpa [SplitBinaryHallSystem.annihilatorCount, A] using
    dvd_mul_of_dvd_right hallFactors (MvPolynomial.C H.polynomialScalar)

/-- The factor-count form of binary Hall localization.  A direction occurs
with full multiplicity `e` between one and `r`, and at least `d-e+1`
polynomial linear factors annihilate it.  Choosing this direction as the
first coordinate gives exactly the exponents in Lemma 3.1. -/
theorem SplitBinaryHallSystem.exists_localized_direction
    {K : Type*} [Field K] {r d : ℕ}
    (H : SplitBinaryHallSystem K r d) (hdegree : r < d) :
    ∃ i : Fin r,
      let e := H.directionMultiplicity i
      let c := H.annihilatorCount i
      1 ≤ e ∧ e ≤ r ∧ d - e + 1 ≤ c := by
  classical
  let adj : Fin r → Fin d → Prop := fun i j =>
    binaryLinearEvaluation (H.polynomialFactor j)
      (H.derivativeDirection i) ≠ 0
  let parallel : Fin r → Fin r → Prop := fun i k =>
    BinaryDirectionsParallel (H.derivativeDirection i)
      (H.derivativeDirection k)
  have hcover : ∀ i k, ¬parallel i k → ∀ j, adj i j ∨ adj k j := by
    intro i k hik j
    exact nonparallel_binaryDirections_cover
      (H.polynomialFactor j) (H.derivativeDirection i)
      (H.derivativeDirection k) (H.polynomialFactor_ne_zero j) hik
  obtain ⟨i, hiPos, hiLe, hiCount⟩ :=
    binaryHall_localization_core adj parallel
      (fun i => binaryDirectionsParallel_refl (H.derivativeDirection i))
      hcover hdegree H.no_matching
  refine ⟨i, ?_, ?_, ?_⟩
  · simpa only [SplitBinaryHallSystem.directionMultiplicity, parallel,
      hallDirectionClass] using hiPos
  · simpa only [SplitBinaryHallSystem.directionMultiplicity, parallel,
      hallDirectionClass] using hiLe
  · simpa only [SplitBinaryHallSystem.directionMultiplicity,
      SplitBinaryHallSystem.annihilatorCount, adj, parallel,
      hallDirectionClass, hallAnnihilators, not_ne_iff] using hiCount

/-- The factor-count conclusion together with the two coordinate-free
divisibility statements equivalent to the displayed Hall normal forms after
choosing the localized direction as the first coordinate. -/
theorem SplitBinaryHallSystem.exists_localized_normal_form
    {K : Type*} [Field K] {r d : ℕ}
    (H : SplitBinaryHallSystem K r d) (hdegree : r < d) :
    ∃ i : Fin r,
      let e := H.directionMultiplicity i
      let c := H.annihilatorCount i
      1 ≤ e ∧ e ≤ r ∧ d - e + 1 ≤ c ∧
        binaryLinearPolynomial (H.derivativeDirection i) ^ e ∣ H.symbol ∧
        binaryLinearPolynomial
          (binaryPerpendicular (H.derivativeDirection i)) ^ c ∣
            H.polynomial := by
  obtain ⟨i, hiPos, hiLe, hiCount⟩ := H.exists_localized_direction hdegree
  exact ⟨i, hiPos, hiLe, hiCount,
    H.directionPower_dvd_symbol i,
    H.annihilatorPower_dvd_polynomial i⟩

end GVC
