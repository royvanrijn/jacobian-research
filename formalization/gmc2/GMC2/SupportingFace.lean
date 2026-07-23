import Mathlib.Data.Rat.Lemmas
import Mathlib.Data.Finset.Basic

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

/-- Rational supporting-face extraction for a finite weighted set.

The assumptions are the one-dimensional form of `0 ∈ conv(S)`.
This finite rational linear-programming lemma is isolated so its eventual
replacement by a general polyhedral-duality theorem is local.
-/
axiom rational_supportingFace
    (S : Finset ℤ) (ν : ℤ → ℕ)
    (hneg : ∃ k ∈ S, k ≤ 0) (hpos : ∃ k ∈ S, 0 ≤ k) :
    LowerFaceCertificate S ν

end GMC2
