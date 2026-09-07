/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import DiscriminantPencils.CuspCount

/-!
# Distinct cusp images

The roots of `H''` count cusp parameters.  This file makes precise the extra
step in Section 4 saying that their images are distinct: it follows from the
absence of a ramified equal-image pair, the parameter-level content of
excluding the `(3,2)` contact stratum.
-/

open Polynomial

namespace DiscriminantPencils

variable {𝕜 : Type*} [Field 𝕜]

/-- No root of `H''` shares its tangent-map image with a distinct parameter.
This is the equal-image consequence of avoiding the `(3,2)` contact locus. -/
def NoRamifiedImageCollision (H : 𝕜[X]) : Prop :=
  ∀ ⦃r u : 𝕜⦄,
    H.derivative.derivative.eval r = 0 →
    tangentMap H r = tangentMap H u →
    r = u

/-- The finite set underlying the multiset of roots of `H''`. -/
noncomputable def cuspParameterFinset (H : 𝕜[X]) : Finset 𝕜 :=
  @Multiset.toFinset 𝕜 (Classical.decEq 𝕜) (cuspParameters H)

/-- The finite set of tangent-map images of cusp parameters. -/
noncomputable def cuspImageFinset (H : 𝕜[X]) : Finset (𝕜 × 𝕜) :=
  letI : DecidableEq 𝕜 := Classical.decEq 𝕜
  Finset.image (tangentMap H) (cuspParameterFinset H)

theorem card_cuspParameterFinset [CharZero 𝕜] [IsAlgClosed 𝕜]
    (H : 𝕜[X]) (n : ℕ) (hdeg : H.natDegree = n)
    (hsq : Squarefree H.derivative.derivative) :
    (cuspParameterFinset H).card = n - 2 := by
  classical
  rw [cuspParameterFinset]
  rw [Multiset.toFinset_card_of_nodup]
  · exact card_cuspParameters H n hdeg
  · exact cuspParameters_nodup H
      (PerfectField.separable_iff_squarefree.mpr hsq)

theorem tangentMap_injOn_cuspParameterFinset (H : 𝕜[X])
    (hcollision : NoRamifiedImageCollision H) :
    Set.InjOn (tangentMap H) (cuspParameterFinset H) := by
  classical
  intro r hr u hu hru
  apply hcollision
  · have hr' : r ∈ cuspParameters H := by
      simpa [cuspParameterFinset] using hr
    exact (isRoot_of_mem_roots hr').eq_zero
  · exact hru

/-- Under squarefreeness and exclusion of ramified collisions, the
discriminant has exactly `n - 2` distinct cusp image points. -/
theorem card_cuspImageFinset [CharZero 𝕜] [IsAlgClosed 𝕜]
    (H : 𝕜[X]) (n : ℕ) (hdeg : H.natDegree = n)
    (hsq : Squarefree H.derivative.derivative)
    (hcollision : NoRamifiedImageCollision H) :
    (cuspImageFinset H).card = n - 2 := by
  classical
  calc
    (cuspImageFinset H).card = (cuspParameterFinset H).card := by
      unfold cuspImageFinset
      exact Finset.card_image_of_injOn
        (tangentMap_injOn_cuspParameterFinset H hcollision)
    _ = n - 2 := card_cuspParameterFinset H n hdeg hsq

end DiscriminantPencils
