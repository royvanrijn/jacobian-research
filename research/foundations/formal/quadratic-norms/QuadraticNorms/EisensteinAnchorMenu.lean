import QuadraticNorms.PrimeCores
import Mathlib.NumberTheory.LSeries.PrimesInAP

namespace QuadraticNorms

private theorem progression_avoiding_anchors (c : ℕ) (hc : c = 1 ∨ c = 3) (anchors : List ℕ) :
    ∃ a m : ℕ, 0 < m ∧ 6 ∣ m ∧ a % 6 = c + 1 ∧
      ∀ f ∈ anchors, ∃ p, p.Prime ∧ p % 6 = 5 ∧ p ∣ m ∧ a ≡ c * f [MOD p] := by
  induction anchors with
  | nil => exact ⟨c + 1, 6, by decide, dvd_rfl, by omega, by simp⟩
  | cons f rest ih =>
    obtain ⟨a, m, hm, h₆, ha, hrest⟩ := ih
    obtain ⟨p, hpm, hp, hp₆⟩ := Nat.forall_exists_prime_gt_and_modEq m
      (q := 6) (a := 5) (by decide) (by decide)
    have hcop : m.Coprime p := hp.coprime_iff_not_dvd.mpr
      (fun h => not_lt_of_ge (Nat.le_of_dvd hm h) hpm) |>.symm
    let witness := Nat.chineseRemainder hcop a (c * f)
    have haw : (witness : ℕ) ≡ a [MOD m] := witness.2.1
    refine ⟨witness, m * p, Nat.mul_pos hm hp.pos, h₆.trans (dvd_mul_right ..), ?_, ?_⟩
    · exact (show (witness : ℕ) % 6 = a % 6 from haw.of_dvd h₆).trans ha
    · intro g hg
      rcases List.mem_cons.mp hg with rfl | h
      · exact ⟨p, hp, hp₆, dvd_mul_left .., witness.2.2⟩
      · obtain ⟨q, hq, hq₆, hqm, hqa⟩ := hrest g h
        exact ⟨q, hq, hq₆, hqm.trans (dvd_mul_right ..), (haw.of_dvd hqm).trans hqa⟩

/-- No finite anchor list covers an entire Eisenstein target class,
whether the anchors enter with coefficient one or three. -/
theorem no_finite_Eisenstein_anchor_menu (c : ℕ) (hc : c = 1 ∨ c = 3)
    (anchors : List ℕ) (B : ℕ) :
    ∃ n > B, n % 6 = c + 1 ∧ ¬ InSumset C3 (Dilate c (fun x => x ∈ anchors)) n := by
  obtain ⟨a, m, hm, h₆, ha, hcover⟩ := progression_avoiding_anchors c hc anchors
  let n := a + m * (B + 1)
  have hnm : n ≡ a [MOD m] := by simp [n, Nat.ModEq]
  have hlarge : B < n := by
    have hmul := Nat.mul_le_mul_right (B + 1) (show 1 ≤ m by omega)
    simp only [one_mul] at hmul
    dsimp [n]
    omega
  refine ⟨n, hlarge, (show n % 6 = a % 6 from hnm.of_dvd h₆).trans ha, ?_⟩
  rintro ⟨x, _, hx, ⟨f, hf, rfl⟩, heq⟩
  obtain ⟨p, hp, hp₆, hpm, hcong⟩ := hcover f hf
  have hnat := (hnm.of_dvd hpm).trans hcong
  rw [← heq, ← zero_add (c * f)] at hnat
  have hd := Nat.modEq_zero_iff_dvd.mp (Nat.ModEq.add_right_cancel' (c * f) (by simpa using hnat))
  have := hx.2 p hp hd
  omega

end QuadraticNorms
