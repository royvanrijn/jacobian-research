import QuadraticNorms.Definitions
import Mathlib.Data.Int.ModEq
import Mathlib.Data.Nat.Prime.Int
import Mathlib.Tactic.NormNum
import Mathlib.Tactic.LinearCombination

namespace QuadraticNorms

private theorem form_modEq {m b c x y x' y' : ℤ}
    (hx : x ≡ x' [ZMOD m]) (hy : y ≡ y' [ZMOD m]) :
    form b c x y ≡ form b c x' y' [ZMOD m] := by
  unfold form
  exact ((hx.pow 2).add (((Int.ModEq.refl b).mul hx).mul hy)).add
    ((Int.ModEq.refl c).mul (hy.pow 2))

private def formResidue {m : ℕ} (b c : ℤ) (x y : Fin m) : ℤ :=
  form b c x.val y.val % m

private def primitiveAt {m : ℕ} (p : ℕ) (x y : Fin m) : Prop :=
  ¬ (p ∣ x.val ∧ p ∣ y.val)

private instance {m : ℕ} (p : ℕ) (x y : Fin m) : Decidable (primitiveAt p x y) :=
  inferInstanceAs (Decidable (¬ (p ∣ x.val ∧ p ∣ y.val)))

private theorem residue_of_primitive {m n : ℕ} {B C x y : ℤ}
    (hm : 0 < m) (hxy : Int.gcd x y = 1) (hn : form B C x y = n) :
    ∃ u v : Fin m, (∀ p : ℕ, 1 < p → p ∣ m → primitiveAt p u v) ∧
      formResidue B C u v = (n : ℤ) % m := by
  have hmz : (0 : ℤ) < m := by exact_mod_cast hm
  have hxlt := Int.emod_lt_of_pos x hmz
  have hylt := Int.emod_lt_of_pos y hmz
  let u : Fin m := ⟨(x % m).toNat, by omega⟩
  let v : Fin m := ⟨(y % m).toNat, by omega⟩
  have hu : (u.val : ℤ) = x % m :=
    Int.toNat_of_nonneg (Int.emod_nonneg _ (by omega))
  have hv : (v.val : ℤ) = y % m :=
    Int.toNat_of_nonneg (Int.emod_nonneg _ (by omega))
  have hx : (u.val : ℤ) ≡ x [ZMOD (m : ℤ)] := by rw [hu]; exact Int.mod_modEq ..
  have hy : (v.val : ℤ) ≡ y [ZMOD (m : ℤ)] := by rw [hv]; exact Int.mod_modEq ..
  refine ⟨u, v, ?_, ?_⟩
  · rintro p hp hpm ⟨hpu, hpv⟩
    have hd : (p : ℤ) ∣ (m : ℤ) := by exact_mod_cast hpm
    have hpx : (p : ℤ) ∣ x :=
      (hx.of_dvd hd).dvd_iff.mp (by exact_mod_cast hpu)
    have hpy : (p : ℤ) ∣ y :=
      (hy.of_dvd hd).dvd_iff.mp (by exact_mod_cast hpv)
    have := Int.dvd_gcd hpx hpy
    rw [hxy] at this
    have := Nat.dvd_one.mp this
    omega
  · simpa only [formResidue, Int.ModEq, hn] using form_modEq hx hy (b := B) (c := C)

/-- The possible residues of a primitive Eisenstein value modulo eighteen. -/
theorem P3_residues {n : ℕ} (h : P3 n) :
    n % 18 = 1 ∨ n % 18 = 3 ∨ n % 18 = 7 ∨ n % 18 = 13 := by
  obtain ⟨_, x, y, hxy, hn⟩ := h
  obtain ⟨u, v, hprim, hres⟩ := residue_of_primitive (m := 18) (by decide) hxy hn
  have table : ∀ u v : Fin 18, primitiveAt 2 u v → primitiveAt 3 u v →
      formResidue (-1) 1 u v = 1 ∨ formResidue (-1) 1 u v = 3 ∨
      formResidue (-1) 1 u v = 7 ∨ formResidue (-1) 1 u v = 13 := by decide
  have := table u v (hprim 2 (by decide) (by decide)) (hprim 3 (by decide) (by decide))
  omega

/-- The possible residues of a primitive x² + 2y² value modulo sixteen. -/
theorem P8_residues {n : ℕ} (h : P8 n) :
    n % 16 = 1 ∨ n % 16 = 2 ∨ n % 16 = 3 ∨ n % 16 = 6 ∨ n % 16 = 9 ∨ n % 16 = 11 := by
  obtain ⟨_, x, y, hxy, hn⟩ := h
  obtain ⟨u, v, hprim, hres⟩ := residue_of_primitive (m := 16) (by decide) hxy hn
  have table : ∀ u v : Fin 16, primitiveAt 2 u v →
      formResidue 0 2 u v = 1 ∨ formResidue 0 2 u v = 2 ∨ formResidue 0 2 u v = 3 ∨
      formResidue 0 2 u v = 6 ∨ formResidue 0 2 u v = 9 ∨ formResidue 0 2 u v = 11 := by decide
  have := table u v (hprim 2 (by decide) (by decide))
  omega

theorem P11_odd {n : ℕ} (h : P11 n) : n % 2 = 1 := by
  obtain ⟨_, x, y, hxy, hn⟩ := h
  obtain ⟨u, v, hprim, hres⟩ := residue_of_primitive (m := 2) (by decide) hxy hn
  have table : ∀ u v : Fin 2, primitiveAt 2 u v → formResidue 1 3 u v = 1 := by decide
  have := table u v (hprim 2 (by decide) (by decide))
  omega

theorem P7_mod_seven {n : ℕ} (h : P7 n) :
    n % 7 = 0 ∨ n % 7 = 1 ∨ n % 7 = 2 ∨ n % 7 = 4 := by
  obtain ⟨_, x, y, hxy, hn⟩ := h
  obtain ⟨u, v, _, hres⟩ := residue_of_primitive (m := 7) (by decide) hxy hn
  have table : ∀ u v : Fin 7, formResidue 1 2 u v = 0 ∨ formResidue 1 2 u v = 1 ∨
      formResidue 1 2 u v = 2 ∨ formResidue 1 2 u v = 4 := by decide
  have := table u v
  omega

theorem P11_mod_eleven {n : ℕ} (h : P11 n) :
    n % 11 = 0 ∨ n % 11 = 1 ∨ n % 11 = 3 ∨ n % 11 = 4 ∨
      n % 11 = 5 ∨ n % 11 = 9 := by
  obtain ⟨_, x, y, hxy, hn⟩ := h
  obtain ⟨u, v, _, hres⟩ := residue_of_primitive (m := 11) (by decide) hxy hn
  have table : ∀ u v : Fin 11, formResidue 1 3 u v = 0 ∨ formResidue 1 3 u v = 1 ∨
      formResidue 1 3 u v = 3 ∨ formResidue 1 3 u v = 4 ∨
      formResidue 1 3 u v = 5 ∨ formResidue 1 3 u v = 9 := by decide
  have := table u v
  omega

/-- If the sum is divisible by the ramified prime, both norm values
are divisible by that prime. -/
theorem P7_ramified_pair {a b : ℕ} (ha : P7 a) (hb : P7 b) (h : 7 ∣ a + b) :
    7 ∣ a ∧ 7 ∣ b := by
  have := P7_mod_seven ha
  have := P7_mod_seven hb
  omega

theorem P11_ramified_pair {a b : ℕ} (ha : P11 a) (hb : P11 b) (h : 11 ∣ a + b) :
    11 ∣ a ∧ 11 ∣ b := by
  have := P11_mod_eleven ha
  have := P11_mod_eleven hb
  omega

/-- The ramified-square obstruction uses the quotient of the sum of
two norms after their linear terms have been divided by the ramified prime. -/
private theorem ramified_square_obstruction {p : ℕ} {c x y z w : ℤ}
    (hp : Nat.Prime p) (hp₂ : p ≠ 2) (hc : 4 * c = (p : ℤ) + 1)
    (hAnisotropic : ∀ a b : ℤ, (p : ℤ) ∣ a ^ 2 + b ^ 2 → (p : ℤ) ∣ a ∧ (p : ℤ) ∣ b)
    (hx : Int.gcd x y = 1) :
    ¬ ((p : ℤ) ^ 2 ∣ form 1 c x y + form 1 c z w) := by
  have heq : 4 * (form 1 c x y + form 1 c z w) =
      (2 * x + y) ^ 2 + (2 * z + w) ^ 2 + (p : ℤ) * (y ^ 2 + w ^ 2) := by
    unfold form
    linear_combination (y ^ 2 + w ^ 2) * hc
  intro hd
  obtain ⟨k, hk⟩ := hd
  have hdlin : (p : ℤ) ∣ (2 * x + y) ^ 2 + (2 * z + w) ^ 2 := by
    use 4 * p * k - (y ^ 2 + w ^ 2)
    rw [hk] at heq
    linear_combination -heq
  obtain ⟨⟨r, hr⟩, ⟨s, hs⟩⟩ := hAnisotropic _ _ hdlin
  rw [hk, hr, hs] at heq
  have hpne : (p : ℤ) ≠ 0 := by exact_mod_cast hp.ne_zero
  have hquot : (p : ℤ) * (4 * k - r ^ 2 - s ^ 2) = y ^ 2 + w ^ 2 := by
    apply mul_left_cancel₀ hpne
    linear_combination heq
  have hpy := (hAnisotropic y w ⟨_, hquot.symm⟩).1
  have hpx2 : (p : ℤ) ∣ 2 * x := by
    simpa using Int.dvd_sub (show (p : ℤ) ∣ 2 * x + y from ⟨r, hr⟩) hpy
  have hxg : p ∣ Int.gcd x y := by
    apply Int.dvd_gcd _ hpy
    have hprime : Prime (p : ℤ) := Nat.prime_iff_prime_int.mp hp
    refine (hprime.dvd_mul.mp hpx2).resolve_left ?_
    intro hp2
    exact hp₂ (Nat.dvd_prime Nat.prime_two |>.mp (by exact_mod_cast hp2) |>.resolve_left hp.ne_one)
  rw [hx] at hxg
  exact hp.ne_one (Nat.dvd_one.mp hxg)

private theorem anisotropic_of_residues (p : ℕ) (hp : 0 < p)
    (h : ∀ a b : Fin p, ((a.val : ℤ) ^ 2 + (b.val : ℤ) ^ 2) % p = 0 →
      a.val = 0 ∧ b.val = 0) :
    ∀ a b : ℤ, (p : ℤ) ∣ a ^ 2 + b ^ 2 → (p : ℤ) ∣ a ∧ (p : ℤ) ∣ b := by
  intro a b hd
  have hpz : (0 : ℤ) < p := by exact_mod_cast hp
  have halt := Int.emod_lt_of_pos a hpz
  have hblt := Int.emod_lt_of_pos b hpz
  let u : Fin p := ⟨(a % p).toNat, by omega⟩
  let v : Fin p := ⟨(b % p).toNat, by omega⟩
  have hu : (u.val : ℤ) = a % p := Int.toNat_of_nonneg (Int.emod_nonneg _ (by omega))
  have hv : (v.val : ℤ) = b % p := Int.toNat_of_nonneg (Int.emod_nonneg _ (by omega))
  have hmod := ((Int.mod_modEq a p).pow 2).add ((Int.mod_modEq b p).pow 2)
  have huv : ((u.val : ℤ) ^ 2 + (v.val : ℤ) ^ 2) % p = 0 := by
    rw [hu, hv]
    exact (show _ % (p : ℤ) = _ % (p : ℤ) from hmod).trans (Int.emod_eq_zero_of_dvd hd)
  obtain ⟨ha, hb⟩ := h u v huv
  exact ⟨Int.dvd_of_emod_eq_zero (by omega : a % p = 0),
    Int.dvd_of_emod_eq_zero (by omega : b % p = 0)⟩

theorem P7_sum_not_dvd_49 {a b : ℕ} (ha : P7 a) (hb : P7 b) : ¬ 49 ∣ a + b := by
  obtain ⟨_, x, y, hxy, hqa⟩ := ha
  obtain ⟨_, z, w, _, hqb⟩ := hb
  have anisotropic := anisotropic_of_residues 7 (by decide) (by decide)
  intro hd
  apply ramified_square_obstruction (c := 2) (z := z) (w := w)
    (by decide : Nat.Prime 7) (by decide) (by decide) anisotropic hxy
  rw [hqa, hqb]
  exact_mod_cast hd

theorem P11_sum_not_dvd_121 {a b : ℕ} (ha : P11 a) (hb : P11 b) : ¬ 121 ∣ a + b := by
  obtain ⟨_, x, y, hxy, hqa⟩ := ha
  obtain ⟨_, z, w, _, hqb⟩ := hb
  have anisotropic := anisotropic_of_residues 11 (by decide) (by decide)
  intro hd
  apply ramified_square_obstruction (c := 3) (z := z) (w := w)
    (by decide : Nat.Prime 11) (by decide) (by decide) anisotropic hxy
  rw [hqa, hqb]
  exact_mod_cast hd

theorem P3_sum_necessary {n : ℕ} (h : InSumset P3 P3 n) :
    2 ≤ n ∧ n % 2 = 0 ∧ n % 18 ≠ 0 ∧ n % 18 ≠ 12 := by
  obtain ⟨a, b, ha, hb, rfl⟩ := h
  have := ha.1
  have := hb.1
  have := P3_residues ha
  have := P3_residues hb
  omega

theorem P8_sum_necessary {n : ℕ} (h : InSumset P8 P8 n) :
    2 ≤ n ∧ ¬ 16 ∣ n := by
  obtain ⟨a, b, ha, hb, rfl⟩ := h
  have := ha.1
  have := hb.1
  have := P8_residues ha
  have := P8_residues hb
  omega

theorem P11_sum_necessary {n : ℕ} (h : InSumset P11 P11 n) :
    2 ≤ n ∧ n % 2 = 0 ∧ ¬ 121 ∣ n := by
  obtain ⟨a, b, ha, hb, rfl⟩ := h
  have := ha.1
  have := hb.1
  have := P11_odd ha
  have := P11_odd hb
  exact ⟨by omega, by omega, P11_sum_not_dvd_121 ha hb⟩

end QuadraticNorms
