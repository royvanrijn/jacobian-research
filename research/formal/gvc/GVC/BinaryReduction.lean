/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import GVC.TopContraction

/-!
# Binary GVC reduction surface

This file records the exact logical interface between Sections 4--6 of the
manuscript.  It does not assume the binary theorem as one opaque proposition:
the remaining arithmetic/geometric obligations are separately named as
envelope closure and ordering of the equal common face.  The complete
common-threshold cutoff is proved here from those finite data.
-/

namespace GVC

open MvPolynomial

abbrev BinaryPolynomial (K : Type*) [CommSemiring K] :=
  MvPolynomial (Fin 2) K

def binaryWeight (w : Fin 2 → ℕ) (α : Fin 2 →₀ ℕ) : ℕ :=
  ∑ i, w i * α i

@[simp] theorem binaryWeight_zero (w : Fin 2 → ℕ) :
    binaryWeight w 0 = 0 := by
  simp [binaryWeight]

theorem binaryWeight_add
    (w : Fin 2 → ℕ) (α β : Fin 2 →₀ ℕ) :
    binaryWeight w (α + β) = binaryWeight w α + binaryWeight w β := by
  simp [binaryWeight, mul_add, Finset.sum_add_distrib]

theorem binaryWeight_mono
    (w : Fin 2 → ℕ) {α β : Fin 2 →₀ ℕ} (hαβ : α ≤ β) :
    binaryWeight w α ≤ binaryWeight w β := by
  apply Finset.sum_le_sum
  intro i _hi
  exact Nat.mul_le_mul_left (w i) (hαβ i)

def PositiveUnequalWeight (w : Fin 2 → ℕ) : Prop :=
  (∀ i, 0 < w i) ∧ w 0 ≠ w 1

/-- The support inequalities at the first common envelope threshold. -/
def HasCommonThreshold
    {K : Type*} [CommSemiring K]
    (symbol p : BinaryPolynomial K) (w : Fin 2 → ℕ) (W : ℕ) : Prop :=
  PositiveUnequalWeight w ∧
    (∀ α ∈ symbol.support, W ≤ binaryWeight w α) ∧
    (∀ β ∈ p.support, binaryWeight w β ≤ W)

/-- A positive integral weight which puts every symbol monomial strictly
above one integer threshold and every polynomial monomial at or below it. -/
def HasStrictWeightThreshold
    {K : Type*} [CommSemiring K]
    (symbol p : BinaryPolynomial K) (w : Fin 2 → ℕ) (W : ℕ) : Prop :=
  (∀ i, 0 < w i) ∧
    (∀ α ∈ symbol.support, W < binaryWeight w α) ∧
    (∀ β ∈ p.support, binaryWeight w β ≤ W)

/-- The equal-weight faces are separated by one integral coordinate cut.
This is the finite cut derived below from the `delta = 0` shifted-ray
ordering; no differential or asymptotic conclusion is built into it. -/
def HasSeparatedCommonFace
    {K : Type*} [CommSemiring K]
    (symbol p : BinaryPolynomial K) (w : Fin 2 → ℕ) (W : ℕ) : Prop :=
  ∃ i c,
    (∀ α ∈ symbol.support, binaryWeight w α = W → c < α i) ∧
    (∀ β ∈ p.support, binaryWeight w β = W → β i ≤ c)

/-- One coordinate orders every symbol exponent on the common face strictly
to the right of every polynomial exponent.  Unlike
`HasSeparatedCommonFace`, this is the direct Newton-segment conclusion and
does not package a preselected integral cut. -/
def HasOrderedCommonFace
    {K : Type*} [CommSemiring K]
    (symbol p : BinaryPolynomial K) (w : Fin 2 → ℕ) (W : ℕ) : Prop :=
  ∃ i, ∀ α ∈ symbol.support, binaryWeight w α = W →
    ∀ β ∈ p.support, binaryWeight w β = W → β i < α i

/-- Finite ordered common faces admit an integral coordinate cut.  The
maximum coordinate on the polynomial face is the required cut. -/
theorem separatedCommonFace_of_orderedCommonFace
    {K : Type*} [CommSemiring K]
    {symbol p : BinaryPolynomial K} {w : Fin 2 → ℕ} {W : ℕ}
    (hordered : HasOrderedCommonFace symbol p w W)
    (hpFace : ∃ β ∈ p.support, binaryWeight w β = W) :
    HasSeparatedCommonFace symbol p w W := by
  rcases hordered with ⟨i, horder⟩
  let face := p.support.filter fun β => binaryWeight w β = W
  have hface : face.Nonempty := by
    obtain ⟨β, hβ, hweight⟩ := hpFace
    exact ⟨β, Finset.mem_filter.mpr ⟨hβ, hweight⟩⟩
  let coordinates := face.image fun β => β i
  have hcoordinates : coordinates.Nonempty := hface.image fun β => β i
  let c := coordinates.max' hcoordinates
  refine ⟨i, c, ?_, ?_⟩
  · intro α hα hαweight
    have hcMem : c ∈ coordinates := Finset.max'_mem coordinates hcoordinates
    obtain ⟨β, hβface, hβcoordinate⟩ := Finset.mem_image.mp hcMem
    have hβData := Finset.mem_filter.mp hβface
    rw [← hβcoordinate]
    exact horder α hα hαweight β hβData.1 hβData.2
  · intro β hβ hβweight
    exact Finset.le_max' coordinates (β i)
      (Finset.mem_image.mpr
        ⟨β, Finset.mem_filter.mpr ⟨hβ, hβweight⟩, rfl⟩)

def binaryCoordinateSupportSum
    {K : Type*} [CommSemiring K]
    (symbol p : BinaryPolynomial K) (i : Fin 2) : ℕ :=
  (∑ α ∈ symbol.support, α i) + ∑ β ∈ p.support, β i

theorem symbolCoordinate_le_binaryCoordinateSupportSum
    {K : Type*} [CommSemiring K]
    (symbol p : BinaryPolynomial K) (i : Fin 2)
    {α : Fin 2 →₀ ℕ} (hα : α ∈ symbol.support) :
    α i ≤ binaryCoordinateSupportSum symbol p i := by
  rw [binaryCoordinateSupportSum]
  have hle : α i ≤ ∑ γ ∈ symbol.support, γ i := by
    exact Finset.single_le_sum (fun γ _ => Nat.zero_le (γ i)) hα
  exact hle.trans (Nat.le_add_right _ _)

theorem polynomialCoordinate_le_binaryCoordinateSupportSum
    {K : Type*} [CommSemiring K]
    (symbol p : BinaryPolynomial K) (i : Fin 2)
    {β : Fin 2 →₀ ℕ} (hβ : β ∈ p.support) :
    β i ≤ binaryCoordinateSupportSum symbol p i := by
  rw [binaryCoordinateSupportSum]
  have hle : β i ≤ ∑ γ ∈ p.support, γ i := by
    exact Finset.single_le_sum (fun γ _ => Nat.zero_le (γ i)) hβ
  exact hle.trans (Nat.le_add_left _ _)

def refinedBinaryWeight
    (w : Fin 2 → ℕ) (i : Fin 2) (scale : ℕ) : Fin 2 → ℕ :=
  fun j => scale * w j + if j = i then 1 else 0

theorem binaryWeight_refinedBinaryWeight
    (w : Fin 2 → ℕ) (i : Fin 2) (scale : ℕ) (α : Fin 2 →₀ ℕ) :
    binaryWeight (refinedBinaryWeight w i scale) α =
      scale * binaryWeight w α + α i := by
  fin_cases i <;>
    simp [binaryWeight, refinedBinaryWeight, Fin.sum_univ_two] <;>
    ring

/-- A separated equal-weight face refines a common threshold to one strict
global integral-weight threshold.  The scale is chosen larger than every
coordinate occurring in either finite support, so off-face weight defects
cannot cross the coordinate cut. -/
theorem strictWeightThreshold_of_separatedCommonFace
    {K : Type*} [CommSemiring K]
    {symbol p : BinaryPolynomial K} {w : Fin 2 → ℕ} {W : ℕ}
    (hcommon : HasCommonThreshold symbol p w W)
    (hseparated : HasSeparatedCommonFace symbol p w W) :
    ∃ v V, HasStrictWeightThreshold symbol p v V := by
  rcases hcommon with ⟨hpositive, hsymbol, hp⟩
  rcases hseparated with ⟨i, c, hsymbolFace, hpFace⟩
  let D := binaryCoordinateSupportSum symbol p i
  let scale := c + D + 1
  let v := refinedBinaryWeight w i scale
  let V := scale * W + c
  have hscalePos : 0 < scale := by
    dsimp [scale]
    omega
  have hcscale : c < scale := by
    dsimp [scale]
    omega
  have hDscale : D < scale := by
    dsimp [scale]
    omega
  refine ⟨v, V, ?_, ?_, ?_⟩
  · intro j
    dsimp [v, refinedBinaryWeight]
    exact Nat.add_pos_left (Nat.mul_pos hscalePos (hpositive.1 j)) _
  · intro α hα
    dsimp [V, v]
    rw [binaryWeight_refinedBinaryWeight]
    by_cases heq : binaryWeight w α = W
    · have hcut := hsymbolFace α hα heq
      rw [heq]
      exact Nat.add_lt_add_left hcut (scale * W)
    · have hweight : W + 1 ≤ binaryWeight w α := by
        exact Nat.succ_le_iff.mpr
          (lt_of_le_of_ne (hsymbol α hα) (Ne.symm heq))
      have hscaleWeight :
          scale * (W + 1) ≤ scale * binaryWeight w α :=
        Nat.mul_le_mul_left scale hweight
      have hbefore : scale * W + c < scale * (W + 1) := by
        rw [Nat.mul_succ]
        omega
      exact lt_of_lt_of_le hbefore
        (le_trans hscaleWeight (Nat.le_add_right _ _))
  · intro β hβ
    dsimp [V, v]
    rw [binaryWeight_refinedBinaryWeight]
    by_cases heq : binaryWeight w β = W
    · have hcut := hpFace β hβ heq
      rw [heq]
      exact Nat.add_le_add_left hcut (scale * W)
    · have hweight : binaryWeight w β + 1 ≤ W := by
        exact Nat.succ_le_iff.mpr (lt_of_le_of_ne (hp β hβ) heq)
      have hcoord : β i ≤ D := by
        exact polynomialCoordinate_le_binaryCoordinateSupportSum
          symbol p i hβ
      have hcoordScale : β i < scale := lt_of_le_of_lt hcoord hDscale
      have hstep :
          scale * binaryWeight w β + β i <
            scale * (binaryWeight w β + 1) := by
        rw [Nat.mul_succ]
        omega
      have hscaled :
          scale * (binaryWeight w β + 1) ≤ scale * W :=
        Nat.mul_le_mul_left scale hweight
      exact le_trans (le_of_lt hstep)
        (le_trans hscaled (Nat.le_add_right _ _))

theorem support_pow_binaryWeight_lower
    {K : Type*} [CommSemiring K] [Nontrivial K]
    (p : BinaryPolynomial K) (w : Fin 2 → ℕ) (W m : ℕ)
    (hlower : ∀ α ∈ p.support, W ≤ binaryWeight w α) :
    ∀ α ∈ (p ^ m).support, m * W ≤ binaryWeight w α := by
  induction m with
  | zero =>
      intro α hα
      simp only [pow_zero] at hα
      have hzero : α = 0 := by
        rw [MvPolynomial.support_one] at hα
        simpa using hα
      subst α
      simp
  | succ m ih =>
      intro α hα
      rw [pow_succ] at hα
      have hadd := MvPolynomial.support_mul (p ^ m) p hα
      obtain ⟨β, hβ, γ, hγ, hsum⟩ := Finset.mem_add.mp hadd
      subst α
      rw [binaryWeight_add]
      have hm := ih β hβ
      have hγW := hlower γ hγ
      simpa [Nat.succ_mul] using Nat.add_le_add hm hγW

theorem support_pow_binaryWeight_upper
    {K : Type*} [CommSemiring K] [Nontrivial K]
    (p : BinaryPolynomial K) (w : Fin 2 → ℕ) (W m : ℕ)
    (hupper : ∀ α ∈ p.support, binaryWeight w α ≤ W) :
    ∀ α ∈ (p ^ m).support, binaryWeight w α ≤ m * W := by
  induction m with
  | zero =>
      intro α hα
      simp only [pow_zero] at hα
      have hzero : α = 0 := by
        rw [MvPolynomial.support_one] at hα
        simpa using hα
      subst α
      simp
  | succ m ih =>
      intro α hα
      rw [pow_succ] at hα
      have hadd := MvPolynomial.support_mul (p ^ m) p hα
      obtain ⟨β, hβ, γ, hγ, hsum⟩ := Finset.mem_add.mp hadd
      subst α
      rw [binaryWeight_add]
      have hm := ih β hβ
      have hγW := hupper γ hγ
      simpa [Nat.succ_mul] using Nat.add_le_add hm hγW

/-- A deliberately simple finite upper bound for the weights in a support.
Using a sum avoids choosing a maximum and also handles the zero polynomial. -/
def binarySupportWeightSum
    {K : Type*} [CommSemiring K]
    (p : BinaryPolynomial K) (w : Fin 2 → ℕ) : ℕ :=
  ∑ α ∈ p.support, binaryWeight w α

theorem binaryWeight_le_supportWeightSum
    {K : Type*} [CommSemiring K]
    (p : BinaryPolynomial K) (w : Fin 2 → ℕ)
    {α : Fin 2 →₀ ℕ} (hα : α ∈ p.support) :
    binaryWeight w α ≤ binarySupportWeightSum p w := by
  rw [binarySupportWeightSum]
  exact Finset.single_le_sum (fun _ _ => Nat.zero_le _) hα

/-- If no supported differential monomial is coordinatewise bounded by a
supported input monomial, then the coefficientwise differential action is
zero. -/
theorem differentialAction_eq_zero_of_support_separated
    {K : Type*} [CommSemiring K]
    (symbol p : BinaryPolynomial K)
    (hsep : ∀ α ∈ symbol.support, ∀ β ∈ p.support, ¬α ≤ β) :
    differentialAction symbol p = 0 := by
  rw [differentialAction, MvPolynomial.sum_def]
  apply Finset.sum_eq_zero
  intro α hα
  rw [MvPolynomial.sum_def]
  apply Finset.sum_eq_zero
  intro β hβ
  rw [multiDescFactorial_eq_zero_of_not_le β α (hsep α hα β hβ)]
  simp

/-- A strict global integral-weight threshold is terminal for GVC.  This is
the finite-support part of Proposition 5.1: powers amplify the unit weight
gap linearly, while a fixed multiplier contributes only a fixed defect. -/
theorem eventuallyMixedPowersVanish_of_strictWeightThreshold
    {K : Type*} [CommSemiring K] [Nontrivial K]
    {symbol p : BinaryPolynomial K} {w : Fin 2 → ℕ} {W : ℕ}
    (hthreshold : HasStrictWeightThreshold symbol p w W) :
    EventuallyMixedPowersVanish symbol p := by
  rcases hthreshold with ⟨_hpositive, hsymbol, hp⟩
  intro q
  let Q := binarySupportWeightSum q w
  refine ⟨Q + 1, ?_⟩
  intro m hm
  apply differentialAction_eq_zero_of_support_separated
  intro α hα β hβ
  have hsymbolLower :
      ∀ γ ∈ symbol.support, W + 1 ≤ binaryWeight w γ := by
    intro γ hγ
    exact Nat.succ_le_iff.mpr (hsymbol γ hγ)
  have hαweight : m * (W + 1) ≤ binaryWeight w α :=
    support_pow_binaryWeight_lower symbol w (W + 1) m hsymbolLower α hα
  have hadd := MvPolynomial.support_mul q (p ^ m) hβ
  obtain ⟨γ, hγ, δ, hδ, hsum⟩ := Finset.mem_add.mp hadd
  subst β
  have hγweight : binaryWeight w γ ≤ Q :=
    binaryWeight_le_supportWeightSum q w hγ
  have hδweight : binaryWeight w δ ≤ m * W :=
    support_pow_binaryWeight_upper p w W m hp δ hδ
  have hβweight :
      binaryWeight w (γ + δ) ≤ Q + m * W := by
    rw [binaryWeight_add]
    exact Nat.add_le_add hγweight hδweight
  intro hαβ
  have hmonotone := binaryWeight_mono w hαβ
  have hQm : Q < m := by omega
  have hmiddle : Q + m * W < m * W + m := by omega
  have hstrict : binaryWeight w (γ + δ) < binaryWeight w α := by
    calc
      binaryWeight w (γ + δ) ≤ Q + m * W := hβweight
      _ < m * W + m := hmiddle
      _ = m * (W + 1) := by rw [Nat.mul_succ]
      _ ≤ binaryWeight w α := hαweight
  exact (not_lt_of_ge hmonotone) hstrict

/-- The complete finite-support deduction in Proposition 5.1.  Once the
equal common faces have one coordinate cut, the refined-weight construction
and the growing power gap prove eventual mixed vanishing. -/
theorem eventuallyMixedPowersVanish_of_separatedCommonFace
    {K : Type*} [CommSemiring K] [Nontrivial K]
    {symbol p : BinaryPolynomial K} {w : Fin 2 → ℕ} {W : ℕ}
    (hcommon : HasCommonThreshold symbol p w W)
    (hseparated : HasSeparatedCommonFace symbol p w W) :
    EventuallyMixedPowersVanish symbol p := by
  obtain ⟨v, V, hstrict⟩ :=
    strictWeightThreshold_of_separatedCommonFace hcommon hseparated
  exact eventuallyMixedPowersVanish_of_strictWeightThreshold hstrict

/-- The common-threshold cutoff stated using only the ordered-face
conclusion of shifted-ray separation.  Empty equality faces are terminal
without any coordinate cut; otherwise finite maximization constructs the
cut used by `eventuallyMixedPowersVanish_of_separatedCommonFace`. -/
theorem eventuallyMixedPowersVanish_of_orderedCommonFace
    {K : Type*} [CommSemiring K] [Nontrivial K]
    {symbol p : BinaryPolynomial K} {w : Fin 2 → ℕ} {W : ℕ}
    (hcommon : HasCommonThreshold symbol p w W)
    (hordered : HasOrderedCommonFace symbol p w W) :
    EventuallyMixedPowersVanish symbol p := by
  rcases hcommon with ⟨hpositive, hsymbol, hp⟩
  by_cases hsymbolFace :
      ∃ α ∈ symbol.support, binaryWeight w α = W
  · by_cases hpFace : ∃ β ∈ p.support, binaryWeight w β = W
    · exact eventuallyMixedPowersVanish_of_separatedCommonFace
        ⟨hpositive, hsymbol, hp⟩
        (separatedCommonFace_of_orderedCommonFace hordered hpFace)
    · by_cases hW : W = 0
      · subst W
        have hpSupport : p.support = ∅ := by
          apply Finset.not_nonempty_iff_eq_empty.mp
          intro hnonempty
          obtain ⟨β, hβ⟩ := hnonempty
          apply hpFace
          have hweight := hp β hβ
          exact ⟨β, hβ, by omega⟩
        have hpzero : p = 0 := MvPolynomial.support_eq_empty.mp hpSupport
        subst p
        intro q
        refine ⟨1, ?_⟩
        intro m hm
        have hmzero : m ≠ 0 := by omega
        simp [hmzero]
      · apply eventuallyMixedPowersVanish_of_strictWeightThreshold
          (w := w) (W := W - 1)
        refine ⟨hpositive.1, ?_, ?_⟩
        · intro α hα
          have := hsymbol α hα
          omega
        · intro β hβ
          have hle := hp β hβ
          have hne : binaryWeight w β ≠ W := by
            intro heq
            exact hpFace ⟨β, hβ, heq⟩
          omega
  · apply eventuallyMixedPowersVanish_of_strictWeightThreshold
      (w := w) (W := W)
    refine ⟨hpositive.1, ?_, hp⟩
    intro α hα
    have hle := hsymbol α hα
    have hne : binaryWeight w α ≠ W := by
      intro heq
      exact hsymbolFace ⟨α, hα, heq⟩
    omega

/-- The two irreducible proof obligations left by the current Lean
development of the binary theorem. -/
structure BinaryEnvelopeBridge
    (K : Type*) [Field K] [CharZero K] where
  /-- Sections 3, 4, and 6: Hall localization plus shifted-ray separation
  force the moving envelopes to reach a common threshold. -/
  envelope_closure : ∀ symbol p : BinaryPolynomial K,
    PurePowersVanish symbol p →
    ∃ w W, HasCommonThreshold symbol p w W
  /-- The `delta = 0` case of shifted-ray separation orders the two complete
  equal-weight faces in one coordinate.  The integral cut itself is
  constructed above from their finite supports. -/
  common_face_ordering : ∀ symbol p : BinaryPolynomial K,
    PurePowersVanish symbol p →
    ∀ w W, HasCommonThreshold symbol p w W →
      HasOrderedCommonFace symbol p w W

/-- The binary theorem follows from envelope closure and equal-face
ordering; common-threshold terminality is now a proved theorem. -/
theorem binary_gvc_of_envelope_bridge
    {K : Type*} [Field K] [CharZero K]
    (B : BinaryEnvelopeBridge K) :
    GeneralizedVanishingConjecture (Fin 2) K := by
  intro symbol p hpure
  obtain ⟨w, W, hcommon⟩ := B.envelope_closure symbol p hpure
  exact eventuallyMixedPowersVanish_of_orderedCommonFace hcommon
    (B.common_face_ordering symbol p hpure w W hcommon)

end GVC
