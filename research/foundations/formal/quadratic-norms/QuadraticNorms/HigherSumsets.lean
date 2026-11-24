import QuadraticNorms.PointwiseCounts
import QuadraticNorms.PrimitiveClassification

namespace QuadraticNorms

/-- Exactly k positive numerical summands, with repetitions allowed. -/
def LengthSum (P : ℕ → Prop) (k n : ℕ) : Prop :=
  ∃ xs : List ℕ, xs.length = k ∧ (∀ x ∈ xs, P x) ∧ xs.sum = n

/-- At most k positive numerical summands. -/
def AtMostSum (P : ℕ → Prop) (k n : ℕ) : Prop :=
  ∃ j ≤ k, LengthSum P j n

theorem LengthSum.append {P : ℕ → Prop} {k j n m : ℕ}
    (hn : LengthSum P k n) (hm : LengthSum P j m) : LengthSum P (k + j) (n + m) := by
  obtain ⟨xs, hx, hxP, hsumx⟩ := hn
  obtain ⟨ys, hy, hyP, hsumy⟩ := hm
  refine ⟨xs ++ ys, by simp [hx, hy], ?_, by simp [hsumx, hsumy]⟩
  intro x hxmem
  exact (List.mem_append.mp hxmem).elim (hxP x) (hyP x)

theorem LengthSum.of_binary {P : ℕ → Prop} {n : ℕ}
    (h : InSumset P P n) : LengthSum P 2 n := by
  obtain ⟨a, b, ha, hb, rfl⟩ := h
  exact ⟨[a, b], rfl, by simpa using And.intro ha hb, by simp⟩

theorem LengthSum.singleton {P : ℕ → Prop} {n : ℕ} (h : P n) : LengthSum P 1 n :=
  ⟨[n], rfl, by simpa, by simp⟩

theorem LengthSum.ones {P : ℕ → Prop} (h : P 1) (k : ℕ) : LengthSum P k k := by
  exact ⟨List.replicate k 1, by simp, fun x hx => List.eq_of_mem_replicate hx ▸ h, by simp⟩

theorem LengthSum.sub_ones {P : ℕ → Prop} {k j n : ℕ}
    (hone : P 1) (hj : j ≤ k) (hn : j ≤ n)
    (h : LengthSum P (k - j) (n - j)) : LengthSum P k n := by
  simpa only [Nat.sub_add_cancel hj, Nat.sub_add_cancel hn] using h.append (LengthSum.ones hone j)

theorem LengthSum.bounds {P : ℕ → Prop} {k n m : ℕ}
    (hpos : ∀ a, P a → 0 < a) (hmod : ∀ a, P a → a ≡ 1 [MOD m])
    (h : LengthSum P k n) : k ≤ n ∧ n ≡ k [MOD m] := by
  obtain ⟨xs, rfl, hxP, rfl⟩ := h
  induction xs with
  | nil => exact ⟨by simp, rfl⟩
  | cons a xs ih =>
    have ha := hxP a (by simp)
    obtain ⟨hle, hxm⟩ := ih (fun x hx => hxP x (List.mem_cons_of_mem _ hx))
    refine ⟨?_, ?_⟩
    · simp only [List.length_cons, List.sum_cons]
      have := hpos a ha
      omega
    · simpa [Nat.add_comm] using (hmod a ha).add hxm

theorem C3_higher_sumsets
    (hbinary : ∀ n, n % 6 = 2 → InSumset C3 C3 n)
    {k n : ℕ} (hk : 2 ≤ k) :
    LengthSum C3 k n ↔ k ≤ n ∧ n ≡ k [MOD 6] := by
  refine ⟨LengthSum.bounds (fun _ h => h.1) (fun _ h => C3_mod_six h), ?_⟩
  rintro ⟨hn, hclass⟩
  have hres : (n - (k - 2)) % 6 = 2 := by
    change n % 6 = k % 6 at hclass
    omega
  apply LengthSum.sub_ones (primeCore_one _) (j := k - 2) (by omega) (by omega)
  have hlen : k - (k - 2) = 2 := by omega
  rw [hlen]
  exact LengthSum.of_binary (hbinary _ hres)

private theorem P3_three_summands
    (hbinary : ∀ n, n % 6 = 2 → InSumset C3 C3 n)
    (hmixed : ∀ n, n % 6 = 4 → InSumset C3 (Dilate 3 C3) n)
    {n : ℕ} (h₃ : 3 ≤ n) (hodd : n % 2 = 1) : LengthSum P3 3 n := by
  have h₁ : P3 1 := restoration3.core _ (primeCore_one _)
  have hp₃ : P3 3 := restoration3.ramified _ (primeCore_one _)
  have h6 : n % 6 = 1 ∨ n % 6 = 3 ∨ n % 6 = 5 := by omega
  suffices ∃ a, InSumset P3 P3 (n - a) ∧ P3 a ∧ a ≤ n by
    obtain ⟨a, hab, ha, hal⟩ := this
    simpa only [Nat.sub_add_cancel hal] using
      (LengthSum.of_binary hab).append (LengthSum.singleton ha)
  rcases h6 with h | h | h
  · refine ⟨3, ?_, hp₃, h₃⟩
    obtain ⟨a, _, ha, ⟨b, hb, rfl⟩, hab⟩ := hmixed (n - 3) (by omega)
    exact ⟨a, 3 * b, restoration3.core a ha, restoration3.ramified b hb, hab⟩
  · refine ⟨1, ?_, h₁, by omega⟩
    obtain ⟨a, b, ha, hb, hab⟩ := hbinary (n - 1) (by omega)
    exact ⟨a, b, restoration3.core a ha, restoration3.core b hb, hab⟩
  · refine ⟨3, ?_, hp₃, h₃⟩
    obtain ⟨a, b, ha, hb, hab⟩ := hbinary (n - 3) (by omega)
    exact ⟨a, b, restoration3.core a ha, restoration3.core b hb, hab⟩

theorem P3_higher_sumsets
    (hbinary : ∀ n, n % 6 = 2 → InSumset C3 C3 n)
    (hmixed : ∀ n, n % 6 = 4 → InSumset C3 (Dilate 3 C3) n)
    {k n : ℕ} (hk : 3 ≤ k) :
    LengthSum P3 k n ↔ k ≤ n ∧ n ≡ k [MOD 2] := by
  refine ⟨LengthSum.bounds (fun _ h => h.1) (fun _ h => ?_), ?_⟩
  · have := P3_residues h
    change _ % 2 = 1
    omega
  · rintro ⟨hn, hclass⟩
    change n % 2 = k % 2 at hclass
    apply LengthSum.sub_ones (restoration3.core _ (primeCore_one _)) (j := k - 3)
      (by omega) (by omega)
    have hlen : k - (k - 3) = 3 := by omega
    rw [hlen]
    exact P3_three_summands hbinary hmixed (by omega) (by omega)

/-- Exact order thirteen, conditional on eventual C₁₂ binary
coverage of the single admissible congruence class. -/
theorem C12_asymptotic_basis_order_thirteen
    (hbinary : ∃ B, ∀ n ≥ B, n % 12 = 2 → InSumset C12 C12 n) :
    (∃ B, ∀ n ≥ B, AtMostSum C12 13 n) ∧
      (∀ B, ∃ n ≥ B, ¬ AtMostSum C12 12 n) := by
  constructor
  · obtain ⟨B, hB⟩ := hbinary
    refine ⟨B + 11, fun n hn => ?_⟩
    let j := (n + 10) % 12
    have hj : j < 12 := Nat.mod_lt _ (by decide)
    have hnjeq : (n - j) % 12 = 2 := by dsimp [j]; omega
    refine ⟨2 + j, by omega, ?_⟩
    have hpair := LengthSum.of_binary (hB (n - j) (by omega) hnjeq)
    simpa [Nat.sub_add_cancel (by omega : j ≤ n)] using
      hpair.append (LengthSum.ones (P := C12) (primeCore_one _) j)
  · intro B
    refine ⟨25 + 60 * B, by omega, ?_⟩
    rintro ⟨j, hj, hsum⟩
    have heq := (LengthSum.bounds (fun _ h => h.1) (fun _ h => C12_mod_twelve h) hsum).2
    change (25 + 60 * B) % 12 = j % 12 at heq
    have hj₁ : j = 1 := by omega
    obtain ⟨xs, hlen, hP, hsum⟩ := hsum
    cases xs with
    | nil => simp [hj₁] at hlen
    | cons a xs =>
      have hxs : xs = [] := List.length_eq_zero_iff.mp (by simpa [hj₁] using hlen)
      subst xs
      have ha : a = 25 + 60 * B := by simpa using hsum
      have hcore := hP a (by simp)
      exact (by decide : ¬ 5 % 12 = 1) (hcore.2 5 (by decide) (by omega))

end QuadraticNorms
