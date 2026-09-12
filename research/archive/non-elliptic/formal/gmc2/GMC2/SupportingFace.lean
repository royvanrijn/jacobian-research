import Mathlib.Data.Rat.Lemmas
import Mathlib.Data.Finset.Basic
import Mathlib.Tactic

/-!
# Rational supporting-face certificates

The geometric content needed downstream is packaged as a finite certificate.
`ρ` and `θ` define a rational supporting line, and `face` records the
nonempty zero-straddling set on which equality holds.
-/

namespace GMC2

/-- A rational lower supporting face for radial orders indexed by integral
angular weights. -/
structure LowerFaceCertificate (S : Finset ℤ) (ν : ℤ → ℕ) where
  intercept : ℚ
  slope : ℚ
  face : Finset ℤ
  face_subset : face ⊆ S
  lower_bound : ∀ k, k ∈ S →
    intercept + slope * (k : ℚ) ≤ (ν k : ℚ)
  eq_on_face : ∀ k, k ∈ face →
    intercept + slope * (k : ℚ) = (ν k : ℚ)
  has_nonpos : ∃ k ∈ face, k ≤ 0
  has_nonneg : ∃ k ∈ face, 0 ≤ k

/-- The height at the vertical axis of the chord joining two weighted
points.  It is only used when `a < b`. -/
private noncomputable def chordHeight (ν : ℤ → ℕ) (a b : ℤ) : ℚ :=
  ((b : ℚ) * (ν a : ℚ) - (a : ℚ) * (ν b : ℚ)) /
    ((b : ℚ) - (a : ℚ))

/-- The slope of the chord joining two weighted points. -/
private noncomputable def chordSlope (ν : ℤ → ℕ) (a b : ℤ) : ℚ :=
  ((ν b : ℚ) - (ν a : ℚ)) / ((b : ℚ) - (a : ℚ))

private noncomputable def raySlope (ν : ℤ → ℕ) (k : ℤ) : ℚ :=
  ((ν k : ℚ) - (ν 0 : ℚ)) / (k : ℚ)

private def crossingPairs (S : Finset ℤ) : Finset (ℤ × ℤ) :=
  (S ×ˢ S).filter fun ab ↦ ab.1 < 0 ∧ 0 < ab.2

/-- Rational supporting-face extraction for a finite weighted set.

The assumptions are the one-dimensional form of `0 ∈ conv(S)`.
The proof constructs the lower convex hull directly: when points occur
strictly on both sides, minimize the height at zero among all crossing
chords.  If zero is an endpoint of the convex hull, an explicit steep line
supports the singleton face at zero.
-/
noncomputable def rational_supportingFace
    (S : Finset ℤ) (ν : ℤ → ℕ)
    (hneg : ∃ k ∈ S, k ≤ 0) (hpos : ∃ k ∈ S, 0 ≤ k) :
    LowerFaceCertificate S ν := by
  classical
  by_cases hstrictNeg : ∃ k ∈ S, k < 0
  · by_cases hstrictPos : ∃ k ∈ S, 0 < k
    · let an := Classical.choose hstrictNeg
      have hanSpec := Classical.choose_spec hstrictNeg
      have hanS : an ∈ S := hanSpec.1
      have han : an < 0 := hanSpec.2
      let bp := Classical.choose hstrictPos
      have hbpSpec := Classical.choose_spec hstrictPos
      have hbpS : bp ∈ S := hbpSpec.1
      have hbp : 0 < bp := hbpSpec.2
      have hcp : (crossingPairs S).Nonempty := by
        refine ⟨(an, bp), ?_⟩
        simp [crossingPairs, hanS, hbpS, han, hbp]
      have hminExists :=
        Finset.exists_min_image (crossingPairs S)
          (fun ab ↦ chordHeight ν ab.1 ab.2) hcp
      let ab := Classical.choose hminExists
      have habSpec := Classical.choose_spec hminExists
      have habC : ab ∈ crossingPairs S := habSpec.1
      have habMin : ∀ ab' ∈ crossingPairs S,
          chordHeight ν ab.1 ab.2 ≤ chordHeight ν ab'.1 ab'.2 :=
        habSpec.2
      let a := ab.1
      let b := ab.2
      have haS : a ∈ S := by
        have hp := (Finset.mem_filter.mp habC).1
        exact (Finset.mem_product.mp hp).1
      have hbS : b ∈ S := by
        have hp := (Finset.mem_filter.mp habC).1
        exact (Finset.mem_product.mp hp).2
      have habPred : a < 0 ∧ 0 < b := by
        simpa [a, b, crossingPairs] using (Finset.mem_filter.mp habC).2
      have ha : a < 0 := habPred.1
      have hb : 0 < b := habPred.2
      let ρ := chordHeight ν a b
      let θ := chordSlope ν a b
      let F : Finset ℤ := S.filter fun k ↦ ρ + θ * (k : ℚ) = (ν k : ℚ)
      have habQ : (0 : ℚ) < (b : ℚ) - (a : ℚ) := by
        have hab : a < b := lt_trans ha hb
        exact sub_pos.mpr (by exact_mod_cast hab)
      have hlineA : ρ + θ * (a : ℚ) = (ν a : ℚ) := by
        dsimp [ρ, θ, chordHeight, chordSlope]
        field_simp
        ring
      have hlineB : ρ + θ * (b : ℚ) = (ν b : ℚ) := by
        dsimp [ρ, θ, chordHeight, chordSlope]
        field_simp
        ring
      by_cases hzeroBelow : 0 ∈ S ∧ (ν 0 : ℚ) < ρ
      · let Pos : Finset ℤ := S.filter fun k ↦ 0 < k
        have hPos : Pos.Nonempty := by
          exact ⟨bp, by simp [Pos, hbpS, hbp]⟩
        have hslopeExists :=
          Finset.exists_min_image Pos (raySlope ν) hPos
        let t := Classical.choose hslopeExists
        have htSpec := Classical.choose_spec hslopeExists
        have htPos : t ∈ Pos := htSpec.1
        have htMin : ∀ k ∈ Pos, raySlope ν t ≤ raySlope ν k := htSpec.2
        have htS : t ∈ S := (Finset.mem_filter.mp htPos).1
        have ht : 0 < t := (Finset.mem_filter.mp htPos).2
        let θ₀ := raySlope ν t
        refine
          { intercept := (ν 0 : ℚ)
            slope := θ₀
            face := {0}
            face_subset := by
              intro k hk
              simpa only [Finset.mem_singleton.mp hk] using hzeroBelow.1
            lower_bound := by
              intro k hkS
              rcases lt_trichotomy k 0 with hk | hk | hk
              · have hpair : (k, t) ∈ crossingPairs S := by
                  simp [crossingPairs, hkS, htS, hk, ht]
                have hρmin : ρ ≤ chordHeight ν k t := by
                  simpa [ρ, a, b] using habMin (k, t) hpair
                have hchord : (ν 0 : ℚ) < chordHeight ν k t :=
                  lt_of_lt_of_le hzeroBelow.2 hρmin
                have hkt : k < t := lt_trans hk ht
                have hden : (0 : ℚ) < (t : ℚ) - (k : ℚ) :=
                  sub_pos.mpr (by exact_mod_cast hkt)
                have hcross :
                    (ν 0 : ℚ) * ((t : ℚ) - (k : ℚ)) <
                      (t : ℚ) * (ν k : ℚ) - (k : ℚ) * (ν t : ℚ) := by
                  exact (lt_div_iff₀ hden).mp (by
                    simpa [chordHeight] using hchord)
                have htQ : (0 : ℚ) < (t : ℚ) := by exact_mod_cast ht
                dsimp [θ₀, raySlope]
                rw [div_mul_eq_mul_div,
                  add_div' _ _ _ (ne_of_gt htQ), div_le_iff₀ htQ]
                nlinarith
              · subst k
                simp
              · have hminSlope : raySlope ν t ≤ raySlope ν k :=
                  htMin k (by simp [Pos, hkS, hk])
                have htQ : (0 : ℚ) < (t : ℚ) := by exact_mod_cast ht
                have hkQ : (0 : ℚ) < (k : ℚ) := by exact_mod_cast hk
                have hcross :
                    ((ν t : ℚ) - (ν 0 : ℚ)) * (k : ℚ) ≤
                      ((ν k : ℚ) - (ν 0 : ℚ)) * (t : ℚ) := by
                  exact (div_le_div_iff₀ htQ hkQ).mp (by
                    simpa [raySlope] using hminSlope)
                dsimp [θ₀, raySlope]
                rw [div_mul_eq_mul_div,
                  add_div' _ _ _ (ne_of_gt htQ), div_le_iff₀ htQ]
                nlinarith
            eq_on_face := by
              intro k hk
              have hk0 : k = 0 := Finset.mem_singleton.mp hk
              subst k
              simp
            has_nonpos := ⟨0, by simp, le_rfl⟩
            has_nonneg := ⟨0, by simp, le_rfl⟩ }
      · have hlower : ∀ k, k ∈ S → ρ + θ * (k : ℚ) ≤ (ν k : ℚ) := by
          intro k hkS
          rcases lt_trichotomy k 0 with hk | hk | hk
          · have hkb : k < b := lt_trans hk hb
            have hpair : (k, b) ∈ crossingPairs S := by
              simp [crossingPairs, hkS, hbS, hk, hb]
            have hmin := habMin (k, b) hpair
            have hdenKB : (0 : ℚ) < (b : ℚ) - (k : ℚ) :=
              sub_pos.mpr (by exact_mod_cast hkb)
            have hcross :
                (((b : ℚ) * (ν a : ℚ) - (a : ℚ) * (ν b : ℚ)) *
                    ((b : ℚ) - (k : ℚ))) ≤
                  (((b : ℚ) * (ν k : ℚ) - (k : ℚ) * (ν b : ℚ)) *
                    ((b : ℚ) - (a : ℚ))) := by
              exact (div_le_div_iff₀ habQ hdenKB).mp (by
                simpa [ρ, a, b, chordHeight] using hmin)
            dsimp [ρ, θ, chordHeight, chordSlope]
            rw [div_mul_eq_mul_div, ← add_div, div_le_iff₀ habQ]
            have hbQ : (0 : ℚ) < (b : ℚ) := by exact_mod_cast hb
            nlinarith
          · subst k
            dsimp [θ]
            simp only [mul_zero, add_zero]
            exact le_of_not_gt fun h ↦ hzeroBelow ⟨hkS, h⟩
          · have hak : a < k := lt_trans ha hk
            have hpair : (a, k) ∈ crossingPairs S := by
              simp [crossingPairs, haS, hkS, ha, hk]
            have hmin := habMin (a, k) hpair
            have hdenAK : (0 : ℚ) < (k : ℚ) - (a : ℚ) :=
              sub_pos.mpr (by exact_mod_cast hak)
            have hcross :
                (((b : ℚ) * (ν a : ℚ) - (a : ℚ) * (ν b : ℚ)) *
                    ((k : ℚ) - (a : ℚ))) ≤
                  (((k : ℚ) * (ν a : ℚ) - (a : ℚ) * (ν k : ℚ)) *
                    ((b : ℚ) - (a : ℚ))) := by
              exact (div_le_div_iff₀ habQ hdenAK).mp (by
                simpa [ρ, a, b, chordHeight] using hmin)
            dsimp [ρ, θ, chordHeight, chordSlope]
            rw [div_mul_eq_mul_div, ← add_div, div_le_iff₀ habQ]
            have haQ : (a : ℚ) < 0 := by exact_mod_cast ha
            nlinarith
        refine
          { intercept := ρ
            slope := θ
            face := F
            face_subset := by
              intro k hk
              exact (Finset.mem_filter.mp hk).1
            lower_bound := hlower
            eq_on_face := by
              intro k hk
              exact (Finset.mem_filter.mp hk).2
            has_nonpos := ⟨a, ?_, ha.le⟩
            has_nonneg := ⟨b, ?_, hb.le⟩ }
        · simp [F, haS, hlineA]
        · simp [F, hbS, hlineB]
    · have hzeroS : 0 ∈ S := by
        let k := Classical.choose hpos
        have hkSpec := Classical.choose_spec hpos
        have hkS : k ∈ S := hkSpec.1
        have hk : 0 ≤ k := hkSpec.2
        have hk0 : k = 0 := le_antisymm (le_of_not_gt fun h ↦ hstrictPos ⟨k, hkS, h⟩) hk
        simpa [hk0] using hkS
      let ρ : ℚ := ν 0
      let θ : ℚ := ν 0
      let F : Finset ℤ := S.filter fun k ↦ ρ + θ * (k : ℚ) = (ν k : ℚ)
      refine
        { intercept := ρ
          slope := θ
          face := F
          face_subset := by
            intro k hk
            exact (Finset.mem_filter.mp hk).1
          lower_bound := by
            intro k hkS
            have hk : k ≤ 0 := le_of_not_gt fun h ↦ hstrictPos ⟨k, hkS, h⟩
            by_cases hk0 : k = 0
            · subst k
              simp [ρ, θ]
            · have hk1 : k ≤ -1 := by omega
              dsimp [ρ, θ]
              have hν : (0 : ℚ) ≤ (ν 0 : ℚ) := by positivity
              have hνk : (0 : ℚ) ≤ (ν k : ℚ) := by positivity
              have hkQ : (k : ℚ) ≤ -1 := by exact_mod_cast hk1
              nlinarith
          eq_on_face := by
            intro k hk
            exact (Finset.mem_filter.mp hk).2
          has_nonpos := ⟨0, by simp [F, hzeroS, ρ, θ], le_rfl⟩
          has_nonneg := ⟨0, by simp [F, hzeroS, ρ, θ], le_rfl⟩ }
  · have hzeroS : 0 ∈ S := by
      let k := Classical.choose hneg
      have hkSpec := Classical.choose_spec hneg
      have hkS : k ∈ S := hkSpec.1
      have hk : k ≤ 0 := hkSpec.2
      have hk0 : k = 0 := le_antisymm hk (le_of_not_gt fun h ↦ hstrictNeg ⟨k, hkS, h⟩)
      simpa [hk0] using hkS
    let ρ : ℚ := ν 0
    let θ : ℚ := -(ν 0 : ℚ)
    let F : Finset ℤ := S.filter fun k ↦ ρ + θ * (k : ℚ) = (ν k : ℚ)
    refine
      { intercept := ρ
        slope := θ
        face := F
        face_subset := by
          intro k hk
          exact (Finset.mem_filter.mp hk).1
        lower_bound := by
          intro k hkS
          have hk : 0 ≤ k := le_of_not_gt fun h ↦ hstrictNeg ⟨k, hkS, h⟩
          by_cases hk0 : k = 0
          · subst k
            simp [ρ, θ]
          · have hk1 : 1 ≤ k := by omega
            dsimp [ρ, θ]
            have hν : (0 : ℚ) ≤ (ν 0 : ℚ) := by positivity
            have hνk : (0 : ℚ) ≤ (ν k : ℚ) := by positivity
            have hkQ : (1 : ℚ) ≤ (k : ℚ) := by exact_mod_cast hk1
            nlinarith
        eq_on_face := by
          intro k hk
          exact (Finset.mem_filter.mp hk).2
        has_nonpos := ⟨0, by simp [F, hzeroS, ρ, θ], le_rfl⟩
        has_nonneg := ⟨0, by simp [F, hzeroS, ρ, θ], le_rfl⟩ }

end GMC2
