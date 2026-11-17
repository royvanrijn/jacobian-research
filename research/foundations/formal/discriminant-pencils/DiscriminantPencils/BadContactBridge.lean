/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import DiscriminantPencils.ContactStrata
import DiscriminantPencils.CuspImages
import Mathlib.Algebra.Polynomial.FieldDivision

/-!
# From equal images to the `(3,2)` contact pattern

This file proves the algebraic bridge used when Section 4 asserts that cusp
images are distinct.  If a zero of `H''` and a distinct parameter have the
same tangent-map image, their common tangent line has contact multiplicities
at least three and two respectively.
-/

open Polynomial

namespace DiscriminantPencils

variable {𝕜 : Type*} [Field 𝕜]

/-- The graph-minus-tangent-line polynomial at `r`. -/
noncomputable def tangentDeviation (H : 𝕜[X]) (r : 𝕜) : 𝕜[X] :=
  pencil H (H.derivative.eval r)
    (r * H.derivative.eval r - H.eval r)

theorem tangentDeviation_ne_zero_of_exactDegree (H : 𝕜[X]) (n : ℕ)
    (hn : 2 ≤ n) (hdeg : H.natDegree = n) (r : 𝕜) :
    tangentDeviation H r ≠ 0 := by
  intro hp
  have hcoeff : (tangentDeviation H r).coeff n = H.coeff n := by
    have hn0 : n ≠ 0 := by omega
    have hn1 : n ≠ 1 := by omega
    simp only [tangentDeviation, pencil, coeff_add, coeff_sub, coeff_C_mul]
    rw [coeff_X_of_ne_one hn1]
    simp [coeff_C, hn0]
  have hH : H ≠ 0 := by
    intro h
    subst H
    simp at hdeg
    omega
  have hlc : H.coeff n ≠ 0 := by
    simpa [leadingCoeff, hdeg] using leadingCoeff_ne_zero.mpr hH
  apply hlc
  rw [← hcoeff, hp]
  simp

theorem tangentDeviation_repeated_at (H : 𝕜[X]) (r : 𝕜) :
    IsRepeatedRoot (tangentDeviation H r) r := by
  simpa [tangentDeviation, tangentMap] using repeatedRoot_at_tangentMap H r

theorem tangentDeviation_repeated_of_sameImage (H : 𝕜[X]) {r u : 𝕜}
    (himage : tangentMap H r = tangentMap H u) :
    IsRepeatedRoot (tangentDeviation H r) u := by
  rw [tangentDeviation]
  apply (isRepeatedRoot_pencil_iff H _ _ u).2
  simpa [tangentMap] using himage.symm

/-- A zero of `H''` gives tangent contact at least three. -/
theorem cube_X_sub_C_dvd_tangentDeviation [CharZero 𝕜]
    (H : 𝕜[X]) (n : ℕ) (hn : 2 ≤ n) (hdeg : H.natDegree = n)
    {r : 𝕜} (hr : H.derivative.derivative.eval r = 0) :
    (X - C r) ^ 3 ∣ tangentDeviation H r := by
  let p := tangentDeviation H r
  have hp : p ≠ 0 := tangentDeviation_ne_zero_of_exactDegree H n hn hdeg r
  have hrep := tangentDeviation_repeated_at H r
  have hroot0 : p.IsRoot r := by
    exact hrep.1
  have hroot1 : p.derivative.IsRoot r := by
    exact hrep.2
  have hroot2 : p.derivative.derivative.IsRoot r := by
    change p.derivative.derivative.eval r = 0
    simpa [p, tangentDeviation, derivative_pencil] using hr
  have hmult : 2 < p.rootMultiplicity r := by
    apply lt_rootMultiplicity_of_isRoot_iterate_derivative hp
    intro m hm
    interval_cases m
    · simpa using hroot0
    · simpa [Function.iterate_succ_apply'] using hroot1
    · simpa [Function.iterate_succ_apply'] using hroot2
  exact (pow_dvd_pow (X - C r) (by omega)).trans
    (pow_rootMultiplicity_dvd p r)

/-- A repeated root gives tangent contact at least two. -/
theorem sq_X_sub_C_dvd_of_repeated [CharZero 𝕜]
    {p : 𝕜[X]} {u : 𝕜} (hp : p ≠ 0) (hu : IsRepeatedRoot p u) :
    (X - C u) ^ 2 ∣ p := by
  have hmult : 1 < p.rootMultiplicity u :=
    (one_lt_rootMultiplicity_iff_isRoot hp).2 hu
  exact (pow_dvd_pow (X - C u) (by omega)).trans
    (pow_rootMultiplicity_dvd p u)

/-- A ramified equal-image pair produces the bad `(3,2)` factor. -/
theorem three_two_contact_of_ramified_sameImage [CharZero 𝕜]
    (H : 𝕜[X]) (n : ℕ) (hn : 2 ≤ n) (hdeg : H.natDegree = n)
    {r u : 𝕜} (hru : r ≠ u)
    (hr : H.derivative.derivative.eval r = 0)
    (himage : tangentMap H r = tangentMap H u) :
    (X - C r) ^ 3 * (X - C u) ^ 2 ∣ tangentDeviation H r := by
  let p := tangentDeviation H r
  have hp : p ≠ 0 := tangentDeviation_ne_zero_of_exactDegree H n hn hdeg r
  apply IsCoprime.mul_dvd
  · exact (isCoprime_X_sub_C_of_isUnit_sub
      (sub_ne_zero.mpr hru).isUnit).pow_left.pow_right
  · exact cube_X_sub_C_dvd_tangentDeviation H n hn hdeg hr
  · exact sq_X_sub_C_dvd_of_repeated hp
      (tangentDeviation_repeated_of_sameImage H himage)

/-- Excluding `(3,2)` factorizations implies the no-collision hypothesis
used by `card_cuspImageFinset`. -/
theorem noRamifiedImageCollision_of_noThreeTwo [CharZero 𝕜]
    (H : 𝕜[X]) (n : ℕ) (hn : 2 ≤ n) (hdeg : H.natDegree = n)
    (hno : ∀ ⦃r u : 𝕜⦄, r ≠ u →
      ¬((X - C r) ^ 3 * (X - C u) ^ 2 ∣ tangentDeviation H r)) :
    NoRamifiedImageCollision H := by
  intro r u hr himage
  by_contra hru
  exact hno hru (three_two_contact_of_ramified_sameImage
    H n hn hdeg hru hr himage)

/-- Simultaneous vanishing of `H''` and `H'''` gives tangent contact at
least four. -/
theorem fourthPower_X_sub_C_dvd_tangentDeviation [CharZero 𝕜]
    (H : 𝕜[X]) (n : ℕ) (hn : 2 ≤ n) (hdeg : H.natDegree = n)
    {r : 𝕜} (hr₂ : H.derivative.derivative.eval r = 0)
    (hr₃ : H.derivative.derivative.derivative.eval r = 0) :
    (X - C r) ^ 4 ∣ tangentDeviation H r := by
  let p := tangentDeviation H r
  have hp : p ≠ 0 := tangentDeviation_ne_zero_of_exactDegree H n hn hdeg r
  have hrep := tangentDeviation_repeated_at H r
  have hroot₀ : p.IsRoot r := hrep.1
  have hroot₁ : p.derivative.IsRoot r := hrep.2
  have hroot₂ : p.derivative.derivative.IsRoot r := by
    change p.derivative.derivative.eval r = 0
    simpa [p, tangentDeviation, derivative_pencil] using hr₂
  have hroot₃ : p.derivative.derivative.derivative.IsRoot r := by
    change p.derivative.derivative.derivative.eval r = 0
    simpa [p, tangentDeviation, derivative_pencil] using hr₃
  have hmult : 3 < p.rootMultiplicity r := by
    apply lt_rootMultiplicity_of_isRoot_iterate_derivative hp
    intro m hm
    interval_cases m
    · simpa using hroot₀
    · simpa [Function.iterate_succ_apply'] using hroot₁
    · simpa [Function.iterate_succ_apply'] using hroot₂
    · simpa [Function.iterate_succ_apply'] using hroot₃
  exact (pow_dvd_pow (X - C r) (by omega)).trans
    (pow_rootMultiplicity_dvd p r)

/-- Excluding fourth-order tangent contact makes every ramification
parameter simple. -/
theorem onlySimpleRamification_of_noFour [CharZero 𝕜]
    (H : 𝕜[X]) (n : ℕ) (hn : 2 ≤ n) (hdeg : H.natDegree = n)
    (hno : ∀ r : 𝕜, ¬((X - C r) ^ 4 ∣ tangentDeviation H r)) :
    ∀ ⦃r : 𝕜⦄,
      H.derivative.derivative.eval r = 0 →
      H.derivative.derivative.derivative.eval r ≠ 0 := by
  intro r hr₂ hr₃
  exact hno r
    (fourthPower_X_sub_C_dvd_tangentDeviation H n hn hdeg hr₂ hr₃)

/-- Three distinct parameters in one tangent-map fiber produce the bad
`(2,2,2)` contact factor. -/
theorem two_two_two_contact_of_sameImage [CharZero 𝕜]
    (H : 𝕜[X]) (n : ℕ) (hn : 2 ≤ n) (hdeg : H.natDegree = n)
    {r u v : 𝕜} (hru : r ≠ u) (hrv : r ≠ v) (huv : u ≠ v)
    (hruImage : tangentMap H r = tangentMap H u)
    (hrvImage : tangentMap H r = tangentMap H v) :
    (X - C r) ^ 2 * (X - C u) ^ 2 * (X - C v) ^ 2 ∣
      tangentDeviation H r := by
  let p := tangentDeviation H r
  have hp : p ≠ 0 := tangentDeviation_ne_zero_of_exactDegree H n hn hdeg r
  have hr₂ : (X - C r) ^ 2 ∣ p :=
    sq_X_sub_C_dvd_of_repeated hp (tangentDeviation_repeated_at H r)
  have hu₂ : (X - C u) ^ 2 ∣ p :=
    sq_X_sub_C_dvd_of_repeated hp
      (tangentDeviation_repeated_of_sameImage H hruImage)
  have hv₂ : (X - C v) ^ 2 ∣ p :=
    sq_X_sub_C_dvd_of_repeated hp
      (tangentDeviation_repeated_of_sameImage H hrvImage)
  have hcop_ru : IsCoprime ((X - C r) ^ 2) ((X - C u) ^ 2) :=
    (isCoprime_X_sub_C_of_isUnit_sub
      (sub_ne_zero.mpr hru).isUnit).pow_left.pow_right
  have hcop_rv : IsCoprime ((X - C r) ^ 2) ((X - C v) ^ 2) :=
    (isCoprime_X_sub_C_of_isUnit_sub
      (sub_ne_zero.mpr hrv).isUnit).pow_left.pow_right
  have hcop_uv : IsCoprime ((X - C u) ^ 2) ((X - C v) ^ 2) :=
    (isCoprime_X_sub_C_of_isUnit_sub
      (sub_ne_zero.mpr huv).isUnit).pow_left.pow_right
  have hru₂ : (X - C r) ^ 2 * (X - C u) ^ 2 ∣ p :=
    hcop_ru.mul_dvd hr₂ hu₂
  exact (hcop_rv.mul_left hcop_uv).mul_dvd hru₂ hv₂

/-- Excluding `(2,2,2)` contact factorizations rules out three-point
tangent-map fibers. -/
theorem noTripleImageFiber_of_noTwoTwoTwo [CharZero 𝕜]
    (H : 𝕜[X]) (n : ℕ) (hn : 2 ≤ n) (hdeg : H.natDegree = n)
    (hno : ∀ ⦃r u v : 𝕜⦄, r ≠ u → r ≠ v → u ≠ v →
      ¬((X - C r) ^ 2 * (X - C u) ^ 2 * (X - C v) ^ 2 ∣
        tangentDeviation H r)) :
    ∀ ⦃r u v : 𝕜⦄,
      tangentMap H r = tangentMap H u →
      tangentMap H r = tangentMap H v →
      r = u ∨ r = v ∨ u = v := by
  intro r u v hruImage hrvImage
  by_contra hdistinct
  push Not at hdistinct
  exact hno hdistinct.1 hdistinct.2.1 hdistinct.2.2
    (two_two_two_contact_of_sameImage H n hn hdeg
      hdistinct.1 hdistinct.2.1 hdistinct.2.2 hruImage hrvImage)

end DiscriminantPencils
