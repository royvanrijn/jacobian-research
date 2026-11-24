import QuadraticNorms.EisensteinCertificate
import QuadraticNorms.CoreDescent

namespace QuadraticNorms

/-! Replay of both complete finite scans in the Eisenstein proof.
The mask specification and its soundness are proved in EisensteinCertificate.
The two closed Boolean evaluations in this module use Lean's native compiler. -/

theorem eisensteinOrdinaryCheck_passes : eisensteinOrdinaryCheck = true := by native_decide

theorem eisensteinMixedCheck_passes : eisensteinMixedCheck = true := by native_decide

theorem eisensteinOrdinary_finite {n : ℕ} (hbound : n < 100000000) (hclass : n % 6 = 2) :
    InSumset C3 C3 n := by
  let k := n / 6
  have heq : n = 6 * k + 2 := by dsimp [k]; omega
  obtain ⟨i, _, hle, hi, hk⟩ := eisensteinOrdinaryCertificate_sound
    eisensteinOrdinaryCheck_passes (k := k) (by omega)
  exact ⟨6 * i + 1, 6 * (k - i) + 1, hi, hk, by omega⟩

theorem eisensteinMixed_finite {n : ℕ} (hbound : n < 100000000) (hclass : n % 6 = 4) :
    InSumset C3 (Dilate 3 C3) n := by
  let k := n / 6
  have heq : n = 6 * k + 4 := by dsimp [k]; omega
  obtain ⟨j, _, hle, hj, hk⟩ := eisensteinMixedCertificate_sound
    eisensteinMixedCheck_passes (k := k) (by omega)
  exact ⟨6 * (k - 3 * j) + 1, 3 * (6 * j + 1), hk, ⟨6 * j + 1, hj, rfl⟩, by omega⟩

/-- The finite Eisenstein scans meet the two large reduced-target
obligations at exactly 10⁸. Split-prime descent then covers all
positive targets in the two admissible classes. -/
theorem eisenstein_core_sumsets_of_large_reduced
    (hregular : ∀ n ≥ 100000000, n % 6 = 2 → IsReduced C3 n → InSumset C3 C3 n)
    (hmixed : ∀ n ≥ 100000000, n % 6 = 4 → IsReduced C3 n → InSumset C3 (Dilate 3 C3) n) :
    (∀ n, n % 6 = 2 → InSumset C3 C3 n) ∧
      (∀ n, n % 6 = 4 → InSumset C3 (Dilate 3 C3) n) := by
  have hOne {n : ℕ} : InSumset C3 (Dilate 1 C3) n ↔ InSumset C3 C3 n := by
    simp [InSumset, Dilate]
  constructor
  · intro n hn
    apply hOne.mp
    refine C3_core_coverage (k := 1) ?_ (by omega) hn
    intro m _ hm hred
    apply hOne.mpr
    by_cases hlt : m < 100000000
    · exact eisensteinOrdinary_finite hlt hm
    · exact hregular m (by omega) hm hred
  · intro n hn
    refine C3_core_coverage (k := 3) ?_ (by omega) hn
    intro m _ hm hred
    by_cases hlt : m < 100000000
    · exact eisensteinMixed_finite hlt hm
    · exact hmixed m (by omega) hm hred

end QuadraticNorms
