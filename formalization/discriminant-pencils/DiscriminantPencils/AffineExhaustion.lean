/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import DiscriminantPencils.BadContactBridge

/-!
# Exhaustion of affine normalization fibers

This file gives a precise, entirely algebraic version of the exhaustion
argument in Section 3 of *Generic Discriminants of Polynomial Tangent
Pencils*.

After excluding higher ramification, collisions involving a ramified
parameter, and fibers with three distinct parameters, every fiber of the
tangent map is exactly one of the following:

* one unramified parameter;
* one ordinary-cusp parameter;
* two distinct unramified parameters with transverse velocities.

Connecting these alternatives to the scheme-theoretic predicates "smooth
point", "ordinary cusp", and "ordinary node" is deliberately left to the
future local-ring layer.  The normalization-fiber classification itself is
proved here without axioms.
-/

open Polynomial

namespace DiscriminantPencils

variable {𝕜 : Type*} [Field 𝕜]

/-- The set-theoretic fiber of the affine tangent map over the image of
`r`. -/
def tangentFiber (H : 𝕜[X]) (r : 𝕜) : Set 𝕜 :=
  {u | tangentMap H u = tangentMap H r}

@[simp] theorem mem_tangentFiber (H : 𝕜[X]) (r u : 𝕜) :
    u ∈ tangentFiber H r ↔ tangentMap H u = tangentMap H r :=
  Iff.rfl

@[simp] theorem self_mem_tangentFiber (H : 𝕜[X]) (r : 𝕜) :
    r ∈ tangentFiber H r := by
  rfl

/-- No affine tangent-map fiber contains three distinct parameters.  This
is the fiber-level consequence of excluding the `(2,2,2)` contact
stratum. -/
def NoTripleImageFiber (H : 𝕜[X]) : Prop :=
  ∀ ⦃r u v : 𝕜⦄,
    tangentMap H r = tangentMap H u →
    tangentMap H r = tangentMap H v →
    r = u ∨ r = v ∨ u = v

/-- Every ramification parameter is simple.  Equivalently for the present
parametrization, every ramified branch has the ordinary-cusp jet.  This is
the parameter-level consequence of excluding the `(4)` contact stratum. -/
def OnlySimpleRamification (H : 𝕜[X]) : Prop :=
  ∀ ⦃r : 𝕜⦄,
    H.derivative.derivative.eval r = 0 →
    H.derivative.derivative.derivative.eval r ≠ 0

/-- The three algebraic conditions needed for the affine exhaustion. -/
structure AffineAdmissible (H : 𝕜[X]) : Prop where
  simpleRamification : OnlySimpleRamification H
  noRamifiedCollision : NoRamifiedImageCollision H
  noTripleFiber : NoTripleImageFiber H

/-- Direct divisibility-form exclusion of the three bad contact patterns
from the paper. -/
structure ContactExclusions (H : 𝕜[X]) : Prop where
  noFour :
    ∀ r : 𝕜, ¬((X - C r) ^ 4 ∣ tangentDeviation H r)
  noThreeTwo :
    ∀ ⦃r u : 𝕜⦄, r ≠ u →
      ¬((X - C r) ^ 3 * (X - C u) ^ 2 ∣ tangentDeviation H r)
  noTwoTwoTwo :
    ∀ ⦃r u v : 𝕜⦄, r ≠ u → r ≠ v → u ≠ v →
      ¬((X - C r) ^ 2 * (X - C u) ^ 2 * (X - C v) ^ 2 ∣
        tangentDeviation H r)

/-- The three contact-factor exclusions imply all hypotheses used by the
affine normalization-fiber exhaustion. -/
theorem affineAdmissible_of_contactExclusions [CharZero 𝕜]
    (H : 𝕜[X]) (n : ℕ) (hn : 2 ≤ n) (hdeg : H.natDegree = n)
    (h : ContactExclusions H) :
    AffineAdmissible H where
  simpleRamification :=
    onlySimpleRamification_of_noFour H n hn hdeg h.noFour
  noRamifiedCollision :=
    noRamifiedImageCollision_of_noThreeTwo H n hn hdeg h.noThreeTwo
  noTripleFiber :=
    noTripleImageFiber_of_noTwoTwoTwo H n hn hdeg h.noTwoTwoTwo

/-- Over an algebraically closed field, simple ramification is equivalent
to separability of the second derivative in the direction needed here. -/
theorem separable_secondDerivative_of_onlySimpleRamification
    [IsAlgClosed 𝕜] (H : 𝕜[X]) (h : OnlySimpleRamification H) :
    H.derivative.derivative.Separable := by
  rw [separable_def,
    Polynomial.isCoprime_iff_aeval_ne_zero_of_isAlgClosed
      (k := 𝕜) 𝕜]
  intro r
  by_cases hr : H.derivative.derivative.eval r = 0
  · right
    simpa [aeval_def] using h hr
  · left
    simpa [aeval_def] using hr

/-- Contact admissibility supplies the squarefreeness hypothesis used in
the cusp count. -/
theorem squarefree_secondDerivative_of_contactExclusions
    [CharZero 𝕜] [IsAlgClosed 𝕜]
    (H : 𝕜[X]) (n : ℕ) (hn : 2 ≤ n) (hdeg : H.natDegree = n)
    (h : ContactExclusions H) :
    Squarefree H.derivative.derivative := by
  rw [← PerfectField.separable_iff_squarefree]
  exact separable_secondDerivative_of_onlySimpleRamification H
    (affineAdmissible_of_contactExclusions H n hn hdeg h).simpleRamification

/-- A singleton unramified normalization fiber. -/
def IsRegularFiber (H : 𝕜[X]) (r : 𝕜) : Prop :=
  IsUnramified H r ∧ tangentFiber H r = {r}

/-- A singleton normalization fiber whose unique branch has the
ordinary-cusp jet. -/
def IsCuspFiber (H : 𝕜[X]) (r : 𝕜) : Prop :=
  IsOrdinaryCuspParameter H r ∧ tangentFiber H r = {r}

/-- A two-point normalization fiber whose two branches are unramified.
Their transversality follows formally from `r ≠ u`. -/
def IsNodeFiber (H : 𝕜[X]) (r u : 𝕜) : Prop :=
  r ≠ u ∧
    tangentMap H r = tangentMap H u ∧
    IsUnramified H r ∧
    IsUnramified H u ∧
    tangentFiber H r = {r, u}

theorem tangentFiber_eq_of_sameImage (H : 𝕜[X]) {r u : 𝕜}
    (h : tangentMap H r = tangentMap H u) :
    tangentFiber H r = tangentFiber H u := by
  ext v
  simp only [mem_tangentFiber]
  rw [h]

theorem tangentFiber_eq_singleton_of_ramified (H : 𝕜[X])
    (hcollision : NoRamifiedImageCollision H) {r : 𝕜}
    (hr : H.derivative.derivative.eval r = 0) :
    tangentFiber H r = {r} := by
  ext u
  simp only [mem_tangentFiber, Set.mem_singleton_iff]
  constructor
  · intro hu
    exact (hcollision (r := r) (u := u) hr hu.symm).symm
  · rintro rfl
    rfl

theorem unramified_of_distinct_sameImage (H : 𝕜[X])
    (hcollision : NoRamifiedImageCollision H) {r u : 𝕜}
    (hru : r ≠ u) (himage : tangentMap H r = tangentMap H u) :
    IsUnramified H r ∧ IsUnramified H u := by
  constructor
  · intro hr
    exact hru (hcollision hr himage)
  · intro hu
    exact hru.symm (hcollision hu himage.symm)

theorem tangentFiber_eq_pair (H : 𝕜[X])
    (htriple : NoTripleImageFiber H) {r u : 𝕜}
    (hru : r ≠ u) (himage : tangentMap H r = tangentMap H u) :
    tangentFiber H r = {r, u} := by
  ext v
  simp only [mem_tangentFiber, Set.mem_insert_iff, Set.mem_singleton_iff]
  constructor
  · intro hv
    rcases htriple himage hv.symm with h | h | h
    · exact (hru h).elim
    · exact Or.inl h.symm
    · exact Or.inr h.symm
  · rintro (rfl | rfl)
    · rfl
    · exact himage.symm

/-- The two velocities in a node fiber are transverse. -/
theorem IsNodeFiber.transverse [CharZero 𝕜] {H : 𝕜[X]} {r u : 𝕜}
    (h : IsNodeFiber H r u) :
    let ar := H.derivative.derivative.eval r
    let au := H.derivative.derivative.eval u
    ar * (u * au) - (r * ar) * au ≠ 0 :=
  distinct_unramified_branches_transverse H h.1 h.2.2.1 h.2.2.2.1

/-- Complete affine normalization-fiber exhaustion.

This is the set-theoretic and jet-theoretic content of the paper's claim
that, after the three bad contact patterns are removed, every affine image
point is regular, an ordinary cusp, or an ordinary transverse node.
-/
theorem affine_fiber_exhaustion [CharZero 𝕜] (H : 𝕜[X])
    (hgood : AffineAdmissible H) (r : 𝕜) :
    IsRegularFiber H r ∨
      IsCuspFiber H r ∨
      ∃ u : 𝕜, IsNodeFiber H r u := by
  by_cases hr : H.derivative.derivative.eval r = 0
  · right
    left
    refine ⟨⟨hr, hgood.simpleRamification hr⟩, ?_⟩
    exact tangentFiber_eq_singleton_of_ramified H
      hgood.noRamifiedCollision hr
  · by_cases hsingleton : tangentFiber H r = {r}
    · exact Or.inl ⟨hr, hsingleton⟩
    · right
      right
      have hexists : ∃ u ∈ tangentFiber H r, u ≠ r := by
        by_contra h
        push Not at h
        apply hsingleton
        ext u
        simp only [Set.mem_singleton_iff]
        constructor
        · exact fun hu ↦ h u hu
        · rintro rfl
          rfl
      obtain ⟨u, hu, hur⟩ := hexists
      have hru : r ≠ u := hur.symm
      have himage : tangentMap H r = tangentMap H u := hu.symm
      obtain ⟨hr', hu'⟩ :=
        unramified_of_distinct_sameImage H hgood.noRamifiedCollision hru himage
      exact ⟨u, hru, himage, hr', hu',
        tangentFiber_eq_pair H hgood.noTripleFiber hru himage⟩

/-- The complete checked affine conclusion supplied by the three contact
exclusions. -/
structure AffineClassification (H : 𝕜[X]) (n : ℕ) : Prop where
  fiberExhaustion :
    ∀ r : 𝕜,
      IsRegularFiber H r ∨
        IsCuspFiber H r ∨
        ∃ u : 𝕜, IsNodeFiber H r u
  cuspImageCount :
    (cuspImageFinset H).card = n - 2

/-- Combined affine theorem: the direct exclusion of `(4)`, `(3,2)`, and
`(2,2,2)` contact factors yields the exhaustive regular/cusp/node
classification and exactly `n - 2` distinct cusp images.

No scheme-theoretic singularity vocabulary or global genus formula is used
in this theorem.
-/
theorem affineClassification_of_contactExclusions
    [CharZero 𝕜] [IsAlgClosed 𝕜]
    (H : 𝕜[X]) (n : ℕ) (hn : 2 ≤ n) (hdeg : H.natDegree = n)
    (h : ContactExclusions H) :
    AffineClassification H n := by
  let hgood := affineAdmissible_of_contactExclusions H n hn hdeg h
  refine ⟨?_, ?_⟩
  · intro r
    exact affine_fiber_exhaustion H hgood r
  · exact card_cuspImageFinset H n hdeg
      (squarefree_secondDerivative_of_contactExclusions H n hn hdeg h)
      hgood.noRamifiedCollision

end DiscriminantPencils
