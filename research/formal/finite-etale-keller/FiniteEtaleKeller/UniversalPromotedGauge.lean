/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.GeneralGaugeJacobian
import FiniteEtaleKeller.UniversalPromotedMap
import Mathlib.RingTheory.Localization.FractionRing

/-!
# The literal universal promoted quadratic-gauge map

For a rank bound `N`, the parameter type is the finite interval
`{j // 4 ≤ j ≤ N}`.  Its coordinate ring supplies the coefficients
`u₄,...,u_N`.  This module first defines the normalized quadratic-gauge map
over an arbitrary coefficient ring without using coefficient division.  It
then specializes the coefficients to the parameter variables and promotes
them to unchanged coordinates.

The determinant proof is transferred to the fraction field of the parameter
ring, where the existing all-degree quadratic-gauge theorem applies to the
seed

`S + S³ + ∑ u_j S^j`.

Injectivity of the coefficient map returns the identity to the polynomial
parameter ring, and `UniversalPromotedMap` supplies the literal full
Jacobian.
-/

noncomputable section

open Matrix
open MvPolynomial
open Polynomial

namespace FiniteEtaleKeller

/-- The parameter indices `4,...,N`. -/
abbrev UniversalHighParameter (N : ℕ) :=
  {j : ℕ // j ∈ Finset.Icc 4 N}

/-- Polynomial coordinate ring of the high-degree parameters. -/
abbrev UniversalParameterRing (K : Type*) [CommRing K] (N : ℕ) :=
  MvPolynomial (UniversalHighParameter N) K

/-- The coefficient `u_j`, extended by zero away from `4 ≤ j ≤ N`. -/
def universalHighCoefficient
    {K : Type*} [CommRing K] (N j : ℕ) :
    UniversalParameterRing K N :=
  if h : j ∈ Finset.Icc 4 N then X ⟨j, h⟩ else 0

@[simp]
theorem universalHighCoefficient_of_mem
    {K : Type*} [CommRing K] (N j : ℕ)
    (h : j ∈ Finset.Icc 4 N) :
    universalHighCoefficient (K := K) N j = X ⟨j, h⟩ := by
  simp [universalHighCoefficient, h]

theorem universalHighCoefficient_of_not_mem
    {K : Type*} [CommRing K] (N j : ℕ)
    (h : j ∉ Finset.Icc 4 N) :
    universalHighCoefficient (K := K) N j = 0 := by
  simp [universalHighCoefficient, h]

section NormalizedGauge

variable {A : Type*} [CommRing A]

/-- The normalized recurrent polynomial `t=1+xy`. -/
def normalizedGaugeT : MvPolynomial (Fin 3) A :=
  1 + MvPolynomial.X 0 * MvPolynomial.X 1

/-- The normalized recurrent polynomial
`q=t²z+y²(1+3t)`. -/
def normalizedGaugeQ : MvPolynomial (Fin 3) A :=
  normalizedGaugeT ^ 2 * MvPolynomial.X 2 +
    MvPolynomial.X 1 ^ 2 *
      (1 + MvPolynomial.C 3 * normalizedGaugeT)

/-- The first vertical coordinate `Π=tq`. -/
def normalizedGaugePi : MvPolynomial (Fin 3) A :=
  normalizedGaugeT * normalizedGaugeQ

/-- The unscaled second vertical coordinate for coefficient family `d`. -/
def normalizedGaugeB (d : ℕ → A) (N : ℕ) :
    MvPolynomial (Fin 3) A :=
  MvPolynomial.X 1 +
    MvPolynomial.C 3 * MvPolynomial.X 0 * normalizedGaugeQ +
    ∑ k ∈ Finset.Icc 4 N,
      MvPolynomial.C ((k : A) * d k) * normalizedGaugeT ^ 2 *
        MvPolynomial.X 0 ^ (k - 2) * normalizedGaugeQ ^ k

/-- The third vertical coordinate for coefficient family `d`. -/
def normalizedGaugeC (d : ℕ → A) (N : ℕ) :
    MvPolynomial (Fin 3) A :=
  MvPolynomial.X 0 *
      (MvPolynomial.C 5 - MvPolynomial.C 3 * normalizedGaugeT) -
    MvPolynomial.X 0 ^ 3 * MvPolynomial.X 2 -
    ∑ k ∈ Finset.Icc 4 N,
      MvPolynomial.C (((k - 2 : ℕ) : A) * d k) *
        (MvPolynomial.X 0 * normalizedGaugeQ) ^ k

/-- The normalized determinant-`-2` vertical map. -/
def normalizedGaugeMap (d : ℕ → A) (N : ℕ) :
    Fin 3 → MvPolynomial (Fin 3) A :=
  ![normalizedGaugePi, normalizedGaugeB d N, normalizedGaugeC d N]

/-- Scale the second output by `s`; taking `s=-1/2` gives determinant one. -/
def normalizedGaugeScaledMap (d : ℕ → A) (N : ℕ) (s : A) :
    Fin 3 → MvPolynomial (Fin 3) A :=
  scaleOutput 1 s 1 (normalizedGaugeMap d N)

end NormalizedGauge

section Seed

variable {L : Type*} [Field L]

/-- Field-valued seed used to invoke the existing all-degree theorem. -/
def normalizedGaugeSeed (d : ℕ → L) (N : ℕ) : L[X] :=
  Polynomial.X + Polynomial.X ^ 3 +
    ∑ k ∈ Finset.Icc 4 N,
      Polynomial.C (d k) * Polynomial.X ^ k

theorem normalizedGaugeSeed_coeff
    (d : ℕ → L) (N j : ℕ) :
    (normalizedGaugeSeed d N).coeff j =
      (if j = 1 then 1 else 0) +
      (if j = 3 then 1 else 0) +
      (if j ∈ Finset.Icc 4 N then d j else 0) := by
  classical
  simp [normalizedGaugeSeed, Polynomial.coeff_X,
    Polynomial.coeff_X_pow, eq_comm]

@[simp]
theorem normalizedGaugeSeed_coeff_one (d : ℕ → L) (N : ℕ) :
    (normalizedGaugeSeed d N).coeff 1 = 1 := by
  rw [normalizedGaugeSeed_coeff]
  simp

@[simp]
theorem normalizedGaugeSeed_coeff_two (d : ℕ → L) (N : ℕ) :
    (normalizedGaugeSeed d N).coeff 2 = 0 := by
  rw [normalizedGaugeSeed_coeff]
  simp [Finset.mem_Icc]

@[simp]
theorem normalizedGaugeSeed_coeff_three (d : ℕ → L) (N : ℕ) :
    (normalizedGaugeSeed d N).coeff 3 = 1 := by
  rw [normalizedGaugeSeed_coeff]
  simp [Finset.mem_Icc]

@[simp]
theorem normalizedGaugeSeed_coeff_of_mem
    (d : ℕ → L) (N j : ℕ) (hj : j ∈ Finset.Icc 4 N) :
    (normalizedGaugeSeed d N).coeff j = d j := by
  rw [normalizedGaugeSeed_coeff]
  have hj1 : j ≠ 1 := by
    have := (Finset.mem_Icc.mp hj).1
    omega
  have hj3 : j ≠ 3 := by
    have := (Finset.mem_Icc.mp hj).1
    omega
  simp [hj1, hj3, hj]

theorem normalizedGaugeSeed_natDegree_le
    (d : ℕ → L) (N : ℕ) (hN : 3 ≤ N) :
    (normalizedGaugeSeed d N).natDegree ≤ N := by
  rw [Polynomial.natDegree_le_iff_coeff_eq_zero]
  intro j hj
  rw [normalizedGaugeSeed_coeff]
  have hj1 : j ≠ 1 := by omega
  have hj3 : j ≠ 3 := by omega
  have hjIcc : j ∉ Finset.Icc 4 N := by
    simp [Finset.mem_Icc]
    omega
  simp [hj1, hj3, hjIcc]

theorem normalizedGaugeSeed_natDegree_eq_three
    (d : ℕ → L) :
    (normalizedGaugeSeed d 3).natDegree = 3 := by
  apply Polynomial.natDegree_eq_of_le_of_coeff_ne_zero
  · exact normalizedGaugeSeed_natDegree_le d 3 (by omega)
  · simp

theorem normalizedGaugeSeed_natDegree_eq
    (d : ℕ → L) (N : ℕ) (hN : 4 ≤ N) (hdN : d N ≠ 0) :
    (normalizedGaugeSeed d N).natDegree = N := by
  apply Polynomial.natDegree_eq_of_le_of_coeff_ne_zero
  · exact normalizedGaugeSeed_natDegree_le d N (by omega)
  · rw [normalizedGaugeSeed_coeff]
    have hN1 : N ≠ 1 := by omega
    have hN3 : N ≠ 3 := by omega
    simp only [hN1, hN3, ↓reduceIte]
    simpa [Finset.mem_Icc, hN] using hdN

end Seed

section Comparison

variable {A B : Type*} [CommRing A] [CommRing B]

/-- The division-free normalized formula is natural under coefficient-ring
maps. -/
theorem normalizedGaugeScaledMap_map
    (φ : A →+* B) (d : ℕ → A) (N : ℕ) (s : A) :
    (fun i => MvPolynomial.map φ (normalizedGaugeScaledMap d N s i)) =
      normalizedGaugeScaledMap (fun k => φ (d k)) N (φ s) := by
  funext i
  fin_cases i <;>
    simp [normalizedGaugeScaledMap, scaleOutput, normalizedGaugeMap,
      normalizedGaugePi, normalizedGaugeB, normalizedGaugeC,
      normalizedGaugeQ, normalizedGaugeT, map_ofNat]

variable {L : Type*} [Field L]

/-- When the seed has the expected top degree, the division-free formula is
definitionally the existing all-degree gauge construction. -/
theorem normalizedGaugeMap_eq_generalGaugeMap
    (d : ℕ → L) (N : ℕ)
    (hdeg : (normalizedGaugeSeed d N).natDegree = N) :
    normalizedGaugeMap d N =
      generalGaugeMap (normalizedGaugeSeed d N) := by
  funext i
  fin_cases i
  · simp [normalizedGaugeMap, normalizedGaugePi, normalizedGaugeB,
      normalizedGaugeC, normalizedGaugeQ, normalizedGaugeT,
      generalGaugeMap, generalGaugePi, generalGaugeB, generalGaugeC,
      generalGaugeQ, generalGaugeT, hdeg]
  · simp [normalizedGaugeMap, normalizedGaugePi, normalizedGaugeB,
      normalizedGaugeC, normalizedGaugeQ, normalizedGaugeT,
      generalGaugeMap, generalGaugePi, generalGaugeB, generalGaugeC,
      generalGaugeQ, generalGaugeT, hdeg]
    apply Finset.sum_congr rfl
    intro k hk
    rw [normalizedGaugeSeed_coeff_of_mem d N k hk]
  · simp [normalizedGaugeMap, normalizedGaugePi, normalizedGaugeB,
      normalizedGaugeC, normalizedGaugeQ, normalizedGaugeT,
      generalGaugeMap, generalGaugePi, generalGaugeB, generalGaugeC,
      generalGaugeQ, generalGaugeT, hdeg]
    apply Finset.sum_congr rfl
    intro k hk
    rw [normalizedGaugeSeed_coeff_of_mem d N k hk]

end Comparison

section UniversalMap

variable (K : Type*) [Field K] [CharZero K]

/-- The three vertical formulas over the polynomial parameter ring. -/
def universalNestedGaugeMap (N : ℕ) :
    Fin 3 → MvPolynomial (Fin 3) (UniversalParameterRing K N) :=
  normalizedGaugeScaledMap
    (universalHighCoefficient (K := K) N) N
    (MvPolynomial.C (-1 / 2 : K))

/-- The literal promoted map: the three gauge coordinates followed by all
unchanged high-degree parameter coordinates. -/
def universalPromotedGaugeMap (N : ℕ) :
    (Fin 3 ⊕ UniversalHighParameter N) →
      MvPolynomial (Fin 3 ⊕ UniversalHighParameter N) K :=
  unchangedParameterPromotion (universalNestedGaugeMap K N)

/-- The nested three-variable family already has determinant one over its
polynomial parameter ring.  The proof checks the identity after the injective
map to the fraction field, where it is the established general gauge map. -/
theorem jacobianDet_universalNestedGaugeMap
    (N : ℕ) (hN : 3 ≤ N) :
    jacobianDet (universalNestedGaugeMap K N) = 1 := by
  let A := UniversalParameterRing K N
  let L := FractionRing A
  let φ : A →+* L := algebraMap A L
  let d : ℕ → A := universalHighCoefficient (K := K) N
  let dL : ℕ → L := fun k => φ (d k)
  let G : L[X] := normalizedGaugeSeed dL N
  have hφ : Function.Injective φ := IsFractionRing.injective A L
  have hdeg : G.natDegree = N := by
    rcases hN.eq_or_lt with rfl | hN'
    · exact normalizedGaugeSeed_natDegree_eq_three dL
    · apply normalizedGaugeSeed_natDegree_eq dL N (by omega)
      change φ (d N) ≠ 0
      rw [← map_zero φ]
      apply hφ.ne
      have hmem : N ∈ Finset.Icc 4 N := by
        simp [Finset.mem_Icc]
        omega
      change universalHighCoefficient (K := K) N N ≠ 0
      rw [universalHighCoefficient_of_mem N N hmem]
      exact MvPolynomial.X_ne_zero
        (⟨N, hmem⟩ : UniversalHighParameter N)
  have hbase :
      normalizedGaugeScaledMap dL N (-1 / 2 : L) =
        generalGaugeJacobianOneMap G := by
    rw [normalizedGaugeScaledMap, normalizedGaugeMap_eq_generalGaugeMap dL N hdeg]
    rfl
  have hmaps :
      (fun i => MvPolynomial.map φ (universalNestedGaugeMap K N i)) =
        generalGaugeJacobianOneMap G := by
    rw [universalNestedGaugeMap,
      normalizedGaugeScaledMap_map φ d N (MvPolynomial.C (-1 / 2 : K))]
    have hs : φ (MvPolynomial.C (-1 / 2 : K)) = (-1 / 2 : L) := by
      calc
        φ (MvPolynomial.C (-1 / 2 : K)) =
            algebraMap A L (algebraMap K A (-1 / 2 : K)) := by
              rw [MvPolynomial.algebraMap_eq]
        _ = algebraMap K L (-1 / 2 : K) := by
              rw [IsScalarTower.algebraMap_apply K A L]
        _ = (-1 / 2 : L) := by
              simp only [map_div₀, map_neg, map_one, map_ofNat]
    rw [hs, hbase]
  apply MvPolynomial.map_injective φ hφ
  rw [← jacobianDet_map φ, hmaps]
  simpa using
    jacobianDet_generalGaugeJacobianOneMap G (by simp [G]) (by simp [G])

/-- The actual full Jacobian of the literal unchanged-parameter promotion is
one. -/
theorem jacobianDet_universalPromotedGaugeMap
    (N : ℕ) (hN : 3 ≤ N) :
    jacobianDet (universalPromotedGaugeMap K N) = 1 :=
  jacobianDet_unchangedParameterPromotion
    (universalNestedGaugeMap K N)
    (jacobianDet_universalNestedGaugeMap K N hN)

/-- There are exactly `N-3` promoted parameter coordinates when `3 ≤ N`. -/
theorem card_universalHighParameter
    (N : ℕ) (hN : 3 ≤ N) :
    Fintype.card (UniversalHighParameter N) = N - 3 := by
  rw [Fintype.card_coe, Nat.card_Icc]
  omega

/-- The literal promoted map is a self-map on an `N`-element coordinate
type. -/
theorem card_universalPromotedCoordinates
    (N : ℕ) (hN : 3 ≤ N) :
    Fintype.card (Fin 3 ⊕ UniversalHighParameter N) = N := by
  rw [Fintype.card_sum, Fintype.card_fin,
    card_universalHighParameter N hN]
  omega

/-- A chosen reindexing of the literal promoted coordinates by `Fin N`. -/
def universalPromotedCoordinateEquiv
    (N : ℕ) (hN : 3 ≤ N) :
    (Fin 3 ⊕ UniversalHighParameter N) ≃ Fin N :=
  Fintype.equivOfCardEq (by
    rw [card_universalPromotedCoordinates N hN, Fintype.card_fin])

#print axioms normalizedGaugeScaledMap_map
#print axioms normalizedGaugeMap_eq_generalGaugeMap
#print axioms jacobianDet_universalNestedGaugeMap
#print axioms jacobianDet_universalPromotedGaugeMap
#print axioms card_universalPromotedCoordinates

end UniversalMap

end FiniteEtaleKeller
