/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.FixedHasseLocal

/-!
# Primitive coordinates and height of the fixed Hasse target line

This module formalizes the elementary arithmetic passage

`(-1, 32a/9, (8a+1)/3) ↦ [9 : -9 : 32a : 24a+3]`.

For `a ≡ 1 (mod 9)`, these integer coordinates are primitive.  For
`a ≥ 1`, their maximum absolute value is exactly `32a`.  Thus the standard
projective height calculation used by the counting argument is now inside
the Lean certificate; only the analytic counting theorem remains external.
-/

noncomputable section

open Function

namespace FiniteEtaleKeller.FixedHasseFamily

/-- The moving affine target of the exact Jacobian-one paper map. -/
def rationalTarget (a : ℕ) : Fin 3 → ℚ :=
  ![-1, targetB (a : ℚ), targetC (a : ℚ)]

/-- Integral homogeneous coordinates for `[1 : rationalTarget a]`. -/
def targetProjectiveCoordinates (a : ℕ) : Fin 4 → ℤ :=
  ![9, -9, ((32 * a : ℕ) : ℤ), ((24 * a + 3 : ℕ) : ℤ)]

/-- Absolute values of the displayed homogeneous coordinates. -/
def targetProjectiveAbsCoordinates (a : ℕ) : Fin 4 → ℕ :=
  ![9, 9, 32 * a, 24 * a + 3]

/-- The gcd of the four displayed homogeneous coordinates. -/
def targetProjectiveContent (a : ℕ) : ℕ :=
  Nat.gcd 9 (Nat.gcd 9 (Nat.gcd (32 * a) (24 * a + 3)))

/-- The maximum absolute value of the four displayed homogeneous
coordinates. -/
def targetProjectiveHeight (a : ℕ) : ℕ :=
  max 9 (max 9 (max (32 * a) (24 * a + 3)))

/-- Dehomogenizing the displayed integral coordinates recovers the moving
affine target exactly. -/
theorem targetProjectiveCoordinates_dehomogenize
    (a : ℕ) (i : Fin 3) :
    (targetProjectiveCoordinates a i.succ : ℚ) /
        (targetProjectiveCoordinates a 0 : ℚ) =
      rationalTarget a i := by
  fin_cases i <;>
    simp [targetProjectiveCoordinates, rationalTarget, targetB, targetC]
  all_goals ring

/-- The natural absolute-coordinate tuple is the coordinatewise integer
absolute value of the signed tuple. -/
theorem targetProjectiveCoordinates_natAbs
    (a : ℕ) (i : Fin 4) :
    (targetProjectiveCoordinates a i).natAbs =
      targetProjectiveAbsCoordinates a i := by
  fin_cases i
  · norm_num [targetProjectiveCoordinates,
      targetProjectiveAbsCoordinates]
  · norm_num [targetProjectiveCoordinates,
      targetProjectiveAbsCoordinates]
  · change Int.natAbs ((32 * a : ℕ) : ℤ) = 32 * a
    exact Int.natAbs_natCast _
  · change Int.natAbs ((24 * a + 3 : ℕ) : ℤ) = 24 * a + 3
    exact Int.natAbs_natCast _

/-- The congruence `a ≡ 1 (mod 9)` makes `9` coprime to `32a`. -/
theorem nine_coprime_thirtyTwo_mul
    (a : ℕ) (hmod9 : a % 9 = 1) :
    Nat.Coprime 9 (32 * a) := by
  have h9a : Nat.Coprime 9 a := by
    rw [Nat.coprime_iff_gcd_eq_one, Nat.gcd_rec, hmod9]
    norm_num
  rw [Nat.coprime_mul_iff_right]
  exact ⟨by norm_num, h9a⟩

/-- The displayed projective coordinates are primitive. -/
theorem targetProjectiveContent_eq_one
    (a : ℕ) (hmod9 : a % 9 = 1) :
    targetProjectiveContent a = 1 := by
  have h9inner :
      Nat.Coprime 9 (Nat.gcd (32 * a) (24 * a + 3)) :=
    Nat.Coprime.coprime_dvd_right
      (Nat.gcd_dvd_left (32 * a) (24 * a + 3))
      (nine_coprime_thirtyTwo_mul a hmod9)
  unfold targetProjectiveContent
  rw [h9inner.gcd_eq_one]
  norm_num

/-- Every displayed coordinate is bounded by `32a` once `a ≥ 1`. -/
theorem targetProjectiveAbsCoordinates_le
    (a : ℕ) (ha : 1 ≤ a) (i : Fin 4) :
    targetProjectiveAbsCoordinates a i ≤ 32 * a := by
  fin_cases i <;>
    simp [targetProjectiveAbsCoordinates] <;>
    omega

/-- The projective height of the moving target is exactly `32a`. -/
theorem targetProjectiveHeight_eq
    (a : ℕ) (ha : 1 ≤ a) :
    targetProjectiveHeight a = 32 * a := by
  have h9 : 9 ≤ 32 * a := by omega
  have hlast : 24 * a + 3 ≤ 32 * a := by omega
  simp [targetProjectiveHeight, max_eq_left hlast, max_eq_right h9]

/-- A height cutoff on the target is exactly the corresponding parameter
cutoff. -/
theorem targetProjectiveHeight_le_iff
    (a B : ℕ) (ha : 1 ≤ a) :
    targetProjectiveHeight a ≤ B ↔ a ≤ B / 32 := by
  rw [targetProjectiveHeight_eq a ha]
  omega

/-- Distinct natural parameters give distinct rational targets. -/
theorem rationalTarget_injective :
    Injective rationalTarget := by
  intro a b hab
  have hcoord := congrFun hab (1 : Fin 3)
  simp [rationalTarget, targetB] at hcoord
  exact_mod_cast hcoord

/-- Passing a finite parameter set to the corresponding targets does not
change its cardinality. -/
theorem card_image_rationalTarget (s : Finset ℕ) :
    (s.image rationalTarget).card = s.card :=
  Finset.card_image_iff.mpr rationalTarget_injective.injOn

/-- One bundled paper-facing certificate for a single admissible parameter. -/
structure PaperParameterCertificate (a : ℕ) : Prop where
  jacobianOne : jacobianDet paperMap = 1
  geometricDegree : paperMapGeometricDegree = 5
  normalizationInverse :
    scaleInput 2 1 1 (scaleOutput (-1) 1 1 paperMap) = baseMap
  fiberRepresented :
    ∀ (A : Type) [CommRing A] [Algebra ℚ A],
      Nonempty
        ((AdjoinRoot (polynomial (a : ℚ)) →ₐ[ℚ] A) ≃
          PaperFiberPoint (a : ℚ) A)
  finiteEtale :
    Algebra.Etale ℚ (AdjoinRoot (polynomial (a : ℚ)))
  rankFive :
    Module.finrank ℚ (AdjoinRoot (polynomial (a : ℚ))) = 5
  noRationalPoint : IsEmpty (PaperFiberPoint (a : ℚ) ℚ)
  realPoint : Nonempty (PaperFiberPoint (a : ℚ) ℝ)
  padicPoint :
    ∀ (p : ℕ) [Fact p.Prime],
      Nonempty (PaperFiberPoint (a : ℚ) ℚ_[p])
  primitiveTarget : targetProjectiveContent a = 1
  targetHeight : targetProjectiveHeight a = 32 * a

/-- The exact parameter conditions from the paper produce the complete
algebraic, local-global, and height certificate. -/
theorem paperParameter_certificate
    (a : ℕ) (ha : 1 < a) (hmod9 : a % 9 = 1)
    (hsupport :
      ∀ q : ℕ, q.Prime → q ∣ a → q % 3 = 1)
    (hnoncube : ¬∃ r : ℚ, r ^ 3 = (a : ℚ)) :
    PaperParameterCertificate a := by
  have ha0 : a ≠ 0 := by omega
  have ha1 : a ≠ 1 := by omega
  have ha0Q : (a : ℚ) ≠ 0 := by exact_mod_cast ha0
  have ha1Q : (a : ℚ) ≠ 1 := by exact_mod_cast ha1
  obtain ⟨hRat, hReal, hPadic⟩ :=
    paperFiberPoint_hasse_certificate
      a ha hmod9 hsupport hnoncube
  exact
    { jacobianOne := jacobianDet_paperMap
      geometricDegree := paperMap_geometricDegree
      normalizationInverse := paperMap_normalization_inverse
      fiberRepresented := fun A => ⟨paperFiberRepresentingEquiv
        (A := A) (a : ℚ) ha0Q ha1Q⟩
      finiteEtale := quotient_etale (a : ℚ) ha0Q ha1Q
      rankFive := quotient_rank (a : ℚ)
      noRationalPoint := hRat
      realPoint := hReal
      padicPoint := hPadic
      primitiveTarget := targetProjectiveContent_eq_one a hmod9
      targetHeight := targetProjectiveHeight_eq a (by omega) }

#print axioms targetProjectiveCoordinates_dehomogenize
#print axioms targetProjectiveContent_eq_one
#print axioms targetProjectiveHeight_eq
#print axioms rationalTarget_injective
#print axioms paperParameter_certificate

end FiniteEtaleKeller.FixedHasseFamily
