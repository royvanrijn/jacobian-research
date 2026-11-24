import QuadraticNorms.FiniteValues
import QuadraticNorms.PrimeCores

namespace QuadraticNorms

/-- The two inclusions needed to restore a ramified factor after
finding the split-prime core representation. -/
structure PrimitiveRestoration (C P : ℕ → Prop) (r : ℕ) : Prop where
  core : ∀ n, C n → P n
  ramified : ∀ n, C n → P (r * n)

private theorem PrimitiveRestoration.lift {C P : ℕ → Prop} {r n : ℕ}
    (h : PrimitiveRestoration C P r) (hn : InSumset C C n) : InSumset P P n := by
  obtain ⟨a, b, ha, hb, hab⟩ := hn
  exact ⟨a, b, h.core a ha, h.core b hb, hab⟩

private theorem PrimitiveRestoration.liftMixed {C P : ℕ → Prop} {r n : ℕ}
    (h : PrimitiveRestoration C P r) (hn : InSumset C (Dilate r C) n) : InSumset P P n := by
  obtain ⟨a, _, ha, ⟨b, hb, rfl⟩, hab⟩ := hn
  exact ⟨a, r * b, h.core a ha, h.ramified b hb, hab⟩

private theorem PrimitiveRestoration.liftDilation {C P : ℕ → Prop} {r n : ℕ}
    (h : PrimitiveRestoration C P r) (hn : InSumset C C n) : InSumset P P (r * n) := by
  obtain ⟨a, b, ha, hb, rfl⟩ := hn
  exact ⟨r * a, r * b, h.ramified a ha, h.ramified b hb, (Nat.mul_add ..).symm⟩

/-- The Eisenstein row of the exact overview follows from the two
core coverage statements and the primitive witnesses for core values. -/
theorem P3_sum_classification_up_to (B : ℕ)
    (hrestore : PrimitiveRestoration C3 P3 3)
    (hordinary : ∀ n, n ≤ B → n % 6 = 2 → InSumset C3 C3 n)
    (hmixed : ∀ n, n ≤ B → n % 6 = 4 → InSumset C3 (Dilate 3 C3) n) (n : ℕ) (hbound : n ≤ B) :
    InSumset P3 P3 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ n % 18 ≠ 0 ∧ n % 18 ≠ 12 := by
  refine ⟨P3_sum_necessary, ?_⟩
  rintro ⟨hn, heven, hzero, htwelve⟩
  by_cases h₂ : n % 6 = 2
  · exact hrestore.lift (hordinary n hbound h₂)
  by_cases h₄ : n % 6 = 4
  · exact hrestore.liftMixed (hmixed n hbound h₄)
  have h₆ : n % 18 = 6 := by omega
  have hquot : (n / 3) % 6 = 2 := by omega
  obtain hn' := hrestore.liftDilation (hordinary (n / 3) ((Nat.div_le_self ..).trans hbound) hquot)
  have heq : 3 * (n / 3) = n := by omega
  rwa [heq] at hn'

theorem P3_sum_classification
    (hrestore : PrimitiveRestoration C3 P3 3)
    (hordinary : ∀ n, n % 6 = 2 → InSumset C3 C3 n)
    (hmixed : ∀ n, n % 6 = 4 → InSumset C3 (Dilate 3 C3) n) (n : ℕ) :
    InSumset P3 P3 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ n % 18 ≠ 0 ∧ n % 18 ≠ 12 :=
  P3_sum_classification_up_to n hrestore (fun m _ => hordinary m)
    (fun m _ => hmixed m) n le_rfl

/-- The −7 row, including the two patches where the core alone fails. -/
theorem P7_sum_classification_up_to (B : ℕ)
    (hrestore : PrimitiveRestoration C7 P7 7)
    (hcore : ∀ n, n ≤ B → 2 ≤ n → ¬ 7 ∣ n → n ≠ 11 → n ≠ 29 → InSumset C7 C7 n)
    (n : ℕ) (hbound : n ≤ B) :
    InSumset P7 P7 n ↔ 2 ≤ n ∧ ¬ 49 ∣ n ∧ n ≠ 7 ∧ n ≠ 77 ∧ n ≠ 203 := by
  refine ⟨P7_sum_necessary, ?_⟩
  rintro ⟨hn, h49, h7, h77, h203⟩
  by_cases hd : 7 ∣ n
  · obtain ⟨m, rfl⟩ := hd
    refine hrestore.liftDilation (hcore m (by omega) (by omega) ?_ (by omega) (by omega))
    intro h
    obtain ⟨k, rfl⟩ := h
    exact h49 ⟨k, by omega⟩
  by_cases h11 : n = 11
  · subst n
    exact ⟨4, 7, ⟨by decide, 1, 1, by decide, by decide⟩,
      ⟨by decide, -1, 2, by decide, by decide⟩, by decide⟩
  by_cases h29 : n = 29
  · subst n
    exact ⟨1, 28, ⟨by decide, 1, 0, by decide, by decide⟩,
      ⟨by decide, 2, 3, by decide, by decide⟩, by decide⟩
  exact hrestore.lift (hcore n hbound hn hd h11 h29)

theorem P7_sum_classification
    (hrestore : PrimitiveRestoration C7 P7 7)
    (hcore : ∀ n, 2 ≤ n → ¬ 7 ∣ n → n ≠ 11 → n ≠ 29 → InSumset C7 C7 n)
    (n : ℕ) :
    InSumset P7 P7 n ↔ 2 ≤ n ∧ ¬ 49 ∣ n ∧ n ≠ 7 ∧ n ≠ 77 ∧ n ≠ 203 :=
  P7_sum_classification_up_to n hrestore (fun m _ => hcore m) n le_rfl

/-- The −8 row follows from the even and mixed-core decompositions. -/
theorem P8_sum_classification_up_to (B : ℕ)
    (hrestore : PrimitiveRestoration C8 P8 2)
    (hordinary : ∀ n, n ≤ B → 2 ≤ n → n % 2 = 0 → ¬ 8 ∣ n → InSumset C8 C8 n)
    (hmixed : ∀ n, n ≤ B → 3 ≤ n → n % 2 = 1 → InSumset C8 (Dilate 2 C8) n)
    (n : ℕ) (hbound : n ≤ B) :
    InSumset P8 P8 n ↔ 2 ≤ n ∧ ¬ 16 ∣ n := by
  refine ⟨P8_sum_necessary, ?_⟩
  rintro ⟨hn, h16⟩
  by_cases hodd : n % 2 = 1
  · exact hrestore.liftMixed (hmixed n hbound (by omega) hodd)
  by_cases h8 : 8 ∣ n
  · obtain ⟨m, rfl⟩ := h8
    have h := hrestore.liftDilation (hordinary (4 * m) (by omega) (by omega) (by omega) (by omega))
    simpa [← Nat.mul_assoc] using h
  · exact hrestore.lift (hordinary n hbound hn (by omega) h8)

theorem P8_sum_classification
    (hrestore : PrimitiveRestoration C8 P8 2)
    (hordinary : ∀ n, 2 ≤ n → n % 2 = 0 → ¬ 8 ∣ n → InSumset C8 C8 n)
    (hmixed : ∀ n, 3 ≤ n → n % 2 = 1 → InSumset C8 (Dilate 2 C8) n)
    (n : ℕ) :
    InSumset P8 P8 n ↔ 2 ≤ n ∧ ¬ 16 ∣ n :=
  P8_sum_classification_up_to n hrestore (fun m _ => hordinary m)
    (fun m _ => hmixed m) n le_rfl

/-- The −11 row follows by restoring at most one ramified factor. -/
theorem P11_sum_classification_up_to (B : ℕ)
    (hrestore : PrimitiveRestoration C11 P11 11)
    (hcore : ∀ n, n ≤ B → 2 ≤ n → n % 2 = 0 → ¬ 11 ∣ n → InSumset C11 C11 n)
    (n : ℕ) (hbound : n ≤ B) :
    InSumset P11 P11 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 121 ∣ n := by
  refine ⟨P11_sum_necessary, ?_⟩
  rintro ⟨hn, heven, h121⟩
  by_cases hd : 11 ∣ n
  · obtain ⟨k, rfl⟩ := hd
    refine hrestore.liftDilation (hcore k (by omega) (by omega) (by omega) ?_)
    intro h
    obtain ⟨m, rfl⟩ := h
    exact h121 ⟨m, by omega⟩
  · exact hrestore.lift (hcore n hbound hn heven hd)

theorem P11_sum_classification
    (hrestore : PrimitiveRestoration C11 P11 11)
    (hcore : ∀ n, 2 ≤ n → n % 2 = 0 → ¬ 11 ∣ n → InSumset C11 C11 n)
    (n : ℕ) :
    InSumset P11 P11 n ↔ 2 ≤ n ∧ n % 2 = 0 ∧ ¬ 121 ∣ n :=
  P11_sum_classification_up_to n hrestore (fun m _ => hcore m) n le_rfl

end QuadraticNorms
