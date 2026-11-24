import QuadraticNorms.PrimeCores
import Mathlib.NumberTheory.EulerProduct.DirichletLSeries
import Mathlib.NumberTheory.LegendreSymbol.QuadraticChar.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Summable

open scoped BigOperators Topology LSeries.notation Classical

namespace QuadraticNorms

private instance : Fact (Nat.Prime 3) := ⟨by decide⟩

/-- The quadratic Dirichlet characters at the primitive levels of the
two imaginary quadratic factors in the C₁₂ Euler identity. -/
def c12ChiNegFour : DirichletCharacter ℂ 4 := ZMod.χ₄.ringHomComp (Int.castRingHom ℂ)
noncomputable def c12ChiNegThree : DirichletCharacter ℂ 3 :=
  (quadraticChar (ZMod 3)).ringHomComp (Int.castRingHom ℂ)

noncomputable def c12ChiTwelve : DirichletCharacter ℂ 12 :=
  DirichletCharacter.changeLevel (by decide : 4 ∣ 12) c12ChiNegFour *
    DirichletCharacter.changeLevel (by decide : 3 ∣ 12) c12ChiNegThree

private theorem chiFour_nat (n : ℕ) :
    c12ChiNegFour n = if n % 2 = 0 then 0 else if n % 4 = 1 then 1 else -1 := by
  unfold c12ChiNegFour
  rw [MulChar.ringHomComp_apply, ZMod.χ₄_nat_eq_if_mod_four]
  split_ifs <;> simp

private theorem chiThree_nat (n : ℕ) :
    c12ChiNegThree n = if n % 3 = 0 then 0 else if n % 3 = 1 then 1 else -1 := by
  have h₂ : quadraticChar (ZMod 3) 2 = -1 :=
    quadraticChar_neg_one_iff_not_isSquare.mpr (by decide)
  unfold c12ChiNegThree
  rw [MulChar.ringHomComp_apply, ← ZMod.natCast_mod n 3]
  have hn := Nat.mod_lt n (by decide : 0 < 3)
  have h : n % 3 = 0 ∨ n % 3 = 1 ∨ n % 3 = 2 := by omega
  rcases h with he | he | he
  · simp [he]
  · simp [he, map_one]
  · simp [he, h₂]

private theorem chiTwelve_nat_of_coprime {n : ℕ} (hn : n.Coprime 12) :
    c12ChiTwelve n = c12ChiNegFour n * c12ChiNegThree n := by
  let u := ZMod.unitOfCoprime n hn
  have hu : ((u : (ZMod 12)ˣ) : ZMod 12) = (n : ZMod 12) := ZMod.coe_unitOfCoprime n hn
  unfold c12ChiTwelve
  rw [← hu, MulChar.mul_apply, DirichletCharacter.changeLevel_eq_cast_of_dvd,
    DirichletCharacter.changeLevel_eq_cast_of_dvd, hu]
  simp [ZMod.cast_natCast, (by decide : 3 ∣ 12), (by decide : 4 ∣ 12)]

private theorem chiTwelve_nat_of_not_coprime {n : ℕ} (hn : ¬ n.Coprime 12) :
    c12ChiTwelve n = 0 :=
  c12ChiTwelve.map_nonunit (fun h => hn ((ZMod.isUnit_iff_coprime ..).mp h))

private theorem coprime_twelve_of_prime {p : ℕ} (hp : p.Prime) (h₂ : p ≠ 2) (h₃ : p ≠ 3) :
    p.Coprime 12 := by
  apply (hp.coprime_iff_not_dvd).mpr
  rw [show (12 : ℕ) = 2 ^ 2 * 3 by decide]
  intro hd
  rcases hp.dvd_mul.mp hd with h | h
  · exact h₂ (((Nat.prime_dvd_prime_iff_eq hp Nat.prime_two).mp (hp.dvd_of_dvd_pow h)))
  · exact h₃ ((Nat.prime_dvd_prime_iff_eq hp (by decide)).mp h)

/-- Good primes excluded by the C₁₂ condition. -/
def c12ExceptionalPrime (p : ℕ) : Prop := p.Prime ∧ p ≠ 2 ∧ p ≠ 3 ∧ p % 12 ≠ 1
instance : DecidablePred c12ExceptionalPrime := fun p => by unfold c12ExceptionalPrime; infer_instance

/-- The complete correction at p in the Euler product for F₁₂(s)^4. -/
noncomputable def c12EulerCorrectionFactor (p : ℕ) (z : ℂ) : ℂ :=
  (if p = 2 then 1 - z ^ 2 else 1) * (if p = 3 then 1 - z ^ 2 else 1) *
    (if c12ExceptionalPrime p then (1 - z ^ 2) ^ 2 else 1)

private theorem c12_local_polynomial {p : ℕ} (hp : p.Prime) (z : ℂ) :
    (1 - z) * (1 - c12ChiNegFour p * z) * (1 - c12ChiNegThree p * z) *
      (1 - c12ChiTwelve p * z) =
      (1 - (if p % 12 = 1 then z else 0)) ^ 4 * c12EulerCorrectionFactor p z := by
  by_cases h₂ : p = 2
  · subst p
    have h : ¬ Nat.Coprime 2 12 := by decide
    rw [chiFour_nat, chiThree_nat, chiTwelve_nat_of_not_coprime h]
    norm_num [c12EulerCorrectionFactor, c12ExceptionalPrime]
    ring
  by_cases h₃ : p = 3
  · subst p
    have h : ¬ Nat.Coprime 3 12 := by decide
    rw [chiFour_nat, chiThree_nat, chiTwelve_nat_of_not_coprime h]
    norm_num [c12EulerCorrectionFactor, c12ExceptionalPrime]
    ring
  have hc := coprime_twelve_of_prime hp h₂ h₃
  have h4 : p % 4 = 1 ∨ p % 4 = 3 := by
    have hodd := hp.mod_two_eq_one_iff_ne_two.mpr h₂
    omega
  have h3 : p % 3 ≠ 0 := by
    intro h
    obtain h | h := (Nat.dvd_prime hp).mp (Nat.dvd_iff_mod_eq_zero.mpr h)
    · omega
    · exact h₃ h.symm
  have hmod : p % 12 = 1 ∨ p % 12 = 5 ∨ p % 12 = 7 ∨ p % 12 = 11 := by omega
  rw [chiTwelve_nat_of_coprime hc, chiFour_nat, chiThree_nat]
  rcases hmod with h | h | h | h
  all_goals
    have hmod2 : p % 2 ≠ 0 := by omega
    have hmod4 : p % 4 = (p % 12) % 4 := by omega
    have hmod3 : p % 3 = (p % 12) % 3 := by omega
    simp [c12EulerCorrectionFactor, c12ExceptionalPrime, hp, h₂, h₃, h, hmod2, hmod3, hmod4]
    ring

/-- The C₁₂ Dirichlet series in its right half-plane of absolute convergence. -/
noncomputable def c12DirichletSeries (s : ℂ) : ℂ := by
  classical
  exact ∑' n : ℕ, if C12 n then (n : ℂ) ^ (-s) else 0

private theorem c12_mul_iff (n m : ℕ) : C12 (n * m) ↔ C12 n ∧ C12 m :=
  ⟨fun h => ⟨primeCore_of_dvd h (dvd_mul_right ..), primeCore_of_dvd h (dvd_mul_left ..)⟩,
    fun h => primeCore_mul h.1 h.2⟩

private noncomputable def c12SummandHom {s : ℂ} (hs : s ≠ 0) : ℕ →*₀ ℂ := by
  classical
  let z := riemannZetaSummandHom hs
  exact
    { toFun := fun n => if C12 n then z n else 0
      map_one' := by
        have h : C12 1 := primeCore_one _
        rw [if_pos h, map_one]
      map_zero' := by
        rw [if_neg (fun h => (Nat.lt_irrefl _ h.1))]
      map_mul' := fun n m => by
        rw [c12_mul_iff]
        by_cases hn : C12 n <;> by_cases hm : C12 m <;> simp [hn, hm, map_mul] }

private theorem c12_summand_formula {s : ℂ} (hs : s ≠ 0) (n : ℕ) :
    c12SummandHom hs n = if C12 n then (n : ℂ) ^ (-s) else 0 := by
  classical
  unfold c12SummandHom riemannZetaSummandHom
  rfl

private theorem c12_s_ne_zero {s : ℂ} (hs : 1 < s.re) : s ≠ 0 := by
  rintro rfl
  norm_num at hs

private theorem c12_summable {s : ℂ} (hs : 1 < s.re) :
    Summable (fun n => ‖c12SummandHom (c12_s_ne_zero hs) n‖) := by
  classical
  apply Summable.of_nonneg_of_le (fun _ => norm_nonneg _) _ (summable_riemannZetaSummand hs)
  intro n
  unfold c12SummandHom
  change ‖if C12 n then riemannZetaSummandHom (c12_s_ne_zero hs) n else 0‖ ≤ _
  split <;> simp

private theorem c12_euler_product_hasProd {s : ℂ} (hs : 1 < s.re) :
    HasProd (fun p : Nat.Primes => (1 - (if (p : ℕ) % 12 = 1 then (p : ℂ) ^ (-s) else 0))⁻¹)
      (c12DirichletSeries s) := by
  have h := EulerProduct.eulerProduct_completely_multiplicative_hasProd (c12_summable hs)
  simp_rw [c12_summand_formula] at h
  unfold c12DirichletSeries
  exact HasProd.congr_fun h (fun p => by
    have hp : C12 (p : ℕ) ↔ (p : ℕ) % 12 = 1 := primeCore_prime p.2
    simp only [hp])

private theorem norm_prime_power_lt_one {s : ℂ} (hs : 1 < s.re) (p : Nat.Primes) :
    ‖(p : ℂ) ^ (-s)‖ < 1 := by
  let q := riemannZetaSummandHom (c12_s_ne_zero hs)
  change ‖q p.val‖ < 1
  exact Summable.norm_lt_one (f := q.toMonoidHom)
    (summable_riemannZetaSummand hs).of_norm p.2.one_lt

private theorem two_and_three_are_prime : Nat.Prime 2 ∧ Nat.Prime 3 := by decide

/-- The Euler identity for C₁₂, on the full half-plane Re(s)>1. -/
theorem c12_euler_product_fourth_power {s : ℂ} (hs : 1 < s.re) :
    c12DirichletSeries s ^ 4 =
      riemannZeta s * (L ↗c12ChiNegFour s) * (L ↗c12ChiNegThree s) * (L ↗c12ChiTwelve s) *
        (1 - ((2 : ℂ) ^ (-s)) ^ 2) * (1 - ((3 : ℂ) ^ (-s)) ^ 2) *
        ∏' p : Nat.Primes, (if c12ExceptionalPrime p then (1 - ((p : ℂ) ^ (-s)) ^ 2) ^ 2 else 1) := by
  -- Keep both exceptional factors and the convergent good-prime correction explicit.
  let z (p : Nat.Primes) := (p : ℂ) ^ (-s)
  let d (p : Nat.Primes) := if c12ExceptionalPrime p then -(z p ^ 2) else 0
  have hsum : Summable (fun p : Nat.Primes => ‖d p‖) := by
    refine Summable.of_nonneg_of_le (fun _ => norm_nonneg _) ?_
      ((summable_riemannZetaSummand hs).subtype Nat.Prime)
    intro p
    change ‖d p‖ ≤ ‖z p‖
    unfold d
    split
    · rw [norm_neg, norm_pow]
      rw [pow_two]
      exact mul_le_of_le_one_right (norm_nonneg _) (norm_prime_power_lt_one hs p).le
    · simp
  have hcorrection : Multipliable (fun p : Nat.Primes =>
      if c12ExceptionalPrime p then (1 - (z p) ^ 2) ^ 2 else 1) := by
    have h := (multipliable_one_add_of_summable hsum).hasProd.pow 2
    have h' : HasProd (fun p : Nat.Primes =>
        if c12ExceptionalPrime p then (1 - (z p) ^ 2) ^ 2 else 1) _ :=
      HasProd.congr_fun h (fun p => by by_cases hp : c12ExceptionalPrime p <;> simp [d, hp, sub_eq_add_neg])
    exact h'.multipliable
  let p₂ : Nat.Primes := ⟨2, two_and_three_are_prime.1⟩
  let p₃ : Nat.Primes := ⟨3, two_and_three_are_prime.2⟩
  have hsingle (q : Nat.Primes) : HasProd (fun p : Nat.Primes => if p = q then 1 - z p ^ 2 else 1)
      (1 - z q ^ 2) := by
    let f (p : Nat.Primes) : ℂ := if p = q then 1 - z p ^ 2 else 1
    have h : HasProd f (∏ p ∈ ({q} : Finset Nat.Primes), f p) :=
      hasProd_prod_of_ne_finset_one (fun p hp => if_neg (by simpa using hp))
    simpa [f] using h
  have hζ := riemannZeta_eulerProduct_hasProd hs
  have h₄ := c12ChiNegFour.LSeries_eulerProduct_hasProd hs
  have h₃ := c12ChiNegThree.LSeries_eulerProduct_hasProd hs
  have h₁₂ := c12ChiTwelve.LSeries_eulerProduct_hasProd hs
  have hR := (((((hζ.mul h₄).mul h₃).mul h₁₂).mul (hsingle p₂)).mul (hsingle p₃)).mul hcorrection.hasProd
  have heq (p : Nat.Primes) :
      (1 - (if (p : ℕ) % 12 = 1 then z p else 0))⁻¹ ^ 4 =
        (1 - z p)⁻¹ * (1 - c12ChiNegFour p * z p)⁻¹ *
          (1 - c12ChiNegThree p * z p)⁻¹ * (1 - c12ChiTwelve p * z p)⁻¹ *
          c12EulerCorrectionFactor p (z p) := by
    have hz := norm_prime_power_lt_one hs p
    have hunit (a : ℂ) (ha : ‖a‖ ≤ 1) : (1 - a * z p) ≠ 0 := by
      intro he
      have hone : a * z p = 1 := (sub_eq_zero.mp he).symm
      have hlt := (mul_le_of_le_one_left (norm_nonneg (z p)) ha).trans_lt hz
      rw [← norm_mul, hone, norm_one] at hlt
      exact (lt_irrefl _ hlt)
    have hz1 : (1 - z p) ≠ 0 := by simpa using hunit 1 (by simp)
    have hA := hunit (c12ChiNegFour p) (c12ChiNegFour.norm_le_one p)
    have hB := hunit (c12ChiNegThree p) (c12ChiNegThree.norm_le_one p)
    have hC := hunit (c12ChiTwelve p) (c12ChiTwelve.norm_le_one p)
    have hnorm : (1 - (if (p : ℕ) % 12 = 1 then z p else 0)) ≠ 0 := by split <;> simp [hz1]
    rw [inv_pow, ← mul_inv, ← mul_inv, ← mul_inv, inv_mul_eq_div]
    apply (eq_div_iff (mul_ne_zero (mul_ne_zero (mul_ne_zero hz1 hA) hB) hC)).mpr
    rw [c12_local_polynomial p.2 (z p), ← mul_assoc, inv_mul_cancel₀ (pow_ne_zero _ hnorm), one_mul]
  have hL := (c12_euler_product_hasProd hs).pow 4
  apply hL.unique
  refine HasProd.congr_fun hR fun p => ?_
  rw [heq p]
  have htwo : p = p₂ ↔ (p : ℕ) = 2 := Subtype.ext_iff
  have hthree : p = p₃ ↔ (p : ℕ) = 3 := Subtype.ext_iff
  have hz (q : Nat.Primes) : (q : ℂ) ^ (-s) = z q := rfl
  simp only [hz, htwo, hthree, c12EulerCorrectionFactor]
  ring

end QuadraticNorms
