import QuadraticNorms.MinusElevenLattice
import QuadraticNorms.NormRepresentationBounds

namespace QuadraticNorms

/-- The signed representations whose coordinate gcd is one. -/
def minusElevenPrimitiveNormFiber (n : ℕ) : Finset MinusElevenVector :=
  (minusElevenNormFiber n).filter fun v => Int.gcd v.1 v.2 = 1

/-- The roots of the defining polynomial at a positive norm. -/
noncomputable def minusElevenRoots (n : ℕ) [NeZero n] : Finset (ZMod n) :=
  normRootResidues 1 3 n

/-- A root of the defining polynomial identifies a signed pair of
opposite vectors of a given positive norm. -/
theorem minusEleven_eq_or_neg_of_root_congruence {n : ℕ} (hn : 0 < n)
    {v w : MinusElevenVector} {t : ZMod n}
    (hv : minusElevenNorm v = n) (hw : minusElevenNorm w = n)
    (hvt : (v.1 : ZMod n) + (v.2 : ZMod n) * t = 0)
    (hwt : (w.1 : ZMod n) + (w.2 : ZMod n) * t = 0) :
    w = v ∨ w = (-v.1, -v.2) :=
  binary_form_eq_or_neg_of_root_congruence (by decide : 4 < 4 * (3 : ℤ) - 1 ^ 2) hn
    (by rw [← minusEleven_norm_cast, hv]) (by rw [← minusEleven_norm_cast, hw]) hvt hwt

/-- Each polynomial root accounts for at most two primitive signed
representations. -/
theorem minusEleven_primitive_norm_fiber_card_le {n : ℕ} (hn : 0 < n) :
    letI : NeZero n := ⟨hn.ne'⟩
    (minusElevenPrimitiveNormFiber n).card ≤ 2 * (minusElevenRoots n).card := by
  letI : NeZero n := ⟨hn.ne'⟩
  rw [minusElevenRoots, ← norm_root_count_apply]
  apply binary_primitive_norm_card_le (by decide) hn
  intro v hv
  obtain ⟨h, hg⟩ := Finset.mem_filter.mp hv
  exact ⟨hg, by rw [← minusEleven_norm_cast, mem_minusEleven_norm_fiber.mp h]⟩

end QuadraticNorms
