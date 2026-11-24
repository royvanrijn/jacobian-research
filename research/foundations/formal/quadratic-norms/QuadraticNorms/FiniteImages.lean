import QuadraticNorms.PrimeCores
import Mathlib.Data.Nat.ChineseRemainder
import Mathlib.NumberTheory.LSeries.PrimesInAP
import Mathlib.GroupTheory.OrderOfElement

namespace QuadraticNorms

/-- Membership in the reduction of a subgroup of units, viewed as a predicate on integers. -/
def HasSubgroupResidue {L : ℕ} (H : Subgroup (ZMod L)ˣ) (n : ℕ) : Prop :=
  ∃ u ∈ H, ((u : (ZMod L)ˣ) : ZMod L) = (n : ZMod L)

theorem core_hasSubgroupResidue {L : ℕ} (H : Subgroup (ZMod L)ˣ) {n : ℕ}
    (hn : PrimeCore (HasSubgroupResidue H) n) : HasSubgroupResidue H n := by
  refine primeCore_induction ⟨1, H.one_mem, by simp⟩ (fun _ _ h => h) ?_ hn
  rintro a b ⟨u, hu, hau⟩ ⟨v, hv, hbv⟩
  exact ⟨u * v, H.mul_mem hu hv, by simp [hau, hbv, Nat.cast_mul]⟩

private theorem coprime_of_modEq {a b m : ℕ} (h : a ≡ b [MOD m]) (hb : b.Coprime m) :
    a.Coprime m := by
  apply (ZMod.isUnit_iff_coprime a m).mp
  rw [(ZMod.natCast_eq_natCast_iff a b m).mpr h]
  exact (ZMod.isUnit_iff_coprime b m).mpr hb

private theorem coprime_lcm_of_coprime {a m n : ℕ}
    (ha : a.Coprime m) (hb : a.Coprime n) : a.Coprime (m.lcm n) :=
  (ha.mul_right hb).coprime_dvd_right (Nat.lcm_dvd (dvd_mul_right ..) (dvd_mul_left ..))

/-- Membership in an allowed set of unit classes, viewed as a predicate
on the integer representing the class. -/
def HasGeneratorResidue {L : ℕ} (S : Set (ZMod L)ˣ) (n : ℕ) : Prop :=
  ∃ u ∈ S, (u : ZMod L) = (n : ZMod L)

private theorem generatorCore_realizes_coprime {L M : ℕ} (hL : 0 < L) (hM : 0 < M)
    (S : Set (ZMod L)ˣ) {u : (ZMod L)ˣ} (hu : u ∈ Subgroup.closure S) :
    ∃ x, PrimeCore (HasGeneratorResidue S) x ∧ x.Coprime M ∧ (x : ZMod L) = (u : ZMod L) := by
  letI : NeZero L := ⟨hL.ne'⟩
  let reachable : Submonoid (ZMod L)ˣ :=
    { carrier := {v | ∃ x, PrimeCore (HasGeneratorResidue S) x ∧ x.Coprime M ∧
        (x : ZMod L) = (v : ZMod L)}
      one_mem' := ⟨1, ⟨by decide, fun p hp hd => (hp.not_dvd_one hd).elim⟩, by simp, by simp⟩
      mul_mem' := by
        rintro a b ⟨x, hx, hxM, ha⟩ ⟨y, hy, hyM, hb⟩
        exact ⟨x * y, primeCore_mul hx hy, hxM.mul_left hyM, by simp [ha, hb]⟩ }
  have hS : S ⊆ reachable := by
    intro v hv
    obtain ⟨p, hpM, hp, hmod⟩ := Nat.forall_exists_prime_gt_and_modEq M hL.ne'
      (ZMod.val_coe_unit_coprime v)
    have hpM' : p.Coprime M := (Nat.Prime.coprime_iff_not_dvd hp).mpr
      (fun hd => (not_le_of_gt hpM) (Nat.le_of_dvd hM hd))
    have hcast := (ZMod.natCast_eq_natCast_iff ..).mpr hmod
    rw [ZMod.natCast_zmod_val] at hcast
    exact ⟨p, (primeCore_prime hp).mpr ⟨v, hv, hcast.symm⟩, hpM', hcast⟩
  have hu' : u ∈ Submonoid.closure S := by
    rw [← Subgroup.closure_toSubmonoid_of_finite]
    exact hu
  exact Submonoid.closure_le.mpr hS hu'

private theorem exists_prime_with_two_residues {M L r : ℕ} (hM : 0 < M) (hL : 0 < L)
    (hr : r.Coprime M) (v : (ZMod L)ˣ)
    (hcompat : r ≡ (v : ZMod L).val [MOD M.gcd L]) :
    ∃ p, p.Prime ∧ p ≡ r [MOD M] ∧ (p : ZMod L) = (v : ZMod L) := by
  letI : NeZero L := ⟨hL.ne'⟩
  let q := Nat.chineseRemainder' hcompat
  have hqM := coprime_of_modEq q.2.1 hr
  have hqL := coprime_of_modEq q.2.2 (ZMod.val_coe_unit_coprime v)
  obtain ⟨p, _, hp, hpq⟩ := Nat.forall_exists_prime_gt_and_modEq 0
    (Nat.lcm_pos hM hL).ne' (coprime_lcm_of_coprime hqM hqL)
  have hpL := (hpq.of_dvd (Nat.dvd_lcm_right M L)).trans q.2.2
  exact ⟨p, hp, (hpq.of_dvd (Nat.dvd_lcm_left M L)).trans q.2.1,
    ((ZMod.natCast_eq_natCast_iff ..).mpr hpL).trans (ZMod.natCast_zmod_val _)⟩

/-- The coset finite-image lemma for arbitrary nonempty allowed sets of
unit classes. Prime support restricts the prefix, while a product of
allowed primes realizes the generated subgroup state of the tail. -/
theorem generator_finite_image {L M r : ℕ} (hL : 0 < L) (hM : 0 < M)
    (S : Set (ZMod L)ˣ) (hne : S.Nonempty) {u : (ZMod L)ˣ} (hu : u ∈ Subgroup.closure S) :
    (∃ x, PrimeCore (HasGeneratorResidue S) x ∧
      (x : ZMod L) = (u : ZMod L) ∧ x ≡ r [MOD M]) ↔
    r ≡ (u : ZMod L).val [MOD M.gcd L] ∧
      ∀ p, p.Prime → p ∣ M → ¬ p ∣ L → ¬ HasGeneratorResidue S p → ¬ p ∣ r := by
  letI : NeZero L := ⟨hL.ne'⟩
  constructor
  · rintro ⟨x, hx, hxu, hxr⟩
    have hxL : x ≡ (u : ZMod L).val [MOD L] := by
      rw [← ZMod.natCast_eq_natCast_iff, ZMod.natCast_zmod_val]
      exact hxu
    refine ⟨(hxr.of_dvd (Nat.gcd_dvd_left ..)).symm.trans
      (hxL.of_dvd (Nat.gcd_dvd_right ..)), ?_⟩
    exact fun p hp hpM _ hforbid hpr => hforbid (hx.2 p hp ((hxr.dvd_iff hpM).mpr hpr))
  · rintro ⟨hr, havoid⟩
    let D := M.gcd r
    have hDpos : 0 < D := Nat.gcd_pos_of_pos_left _ hM
    have hDM : D ∣ M := Nat.gcd_dvd_left ..
    have hDr : D ∣ r := Nat.gcd_dvd_right ..
    have hDC : PrimeCore (HasGeneratorResidue S) D := by
      refine ⟨hDpos, fun p hp hpD => ?_⟩
      by_contra! hn
      apply havoid p hp (hpD.trans hDM) _ hn (hpD.trans hDr)
      intro hpL
      have hpval := (hr.dvd_iff (Nat.dvd_gcd (hpD.trans hDM) hpL)).mp (hpD.trans hDr)
      exact hp.not_dvd_one ((Nat.dvd_gcd hpval hpL).trans (dvd_of_eq (ZMod.val_coe_unit_coprime u)))
    have hDH : PrimeCore (HasSubgroupResidue (Subgroup.closure S)) D := by
      refine ⟨hDpos, fun p hp hpd => ?_⟩
      obtain ⟨w, hw, he⟩ := hDC.2 p hp hpd
      exact ⟨w, Subgroup.subset_closure hw, he⟩
    obtain ⟨d, hd, hdD⟩ := core_hasSubgroupResidue (Subgroup.closure S) hDH
    let M' := M / D
    let r' := r / D
    letI : NeZero M' := ⟨(Nat.div_pos (Nat.le_of_dvd hM hDM) hDpos).ne'⟩
    have hM'eq : D * M' = M := Nat.mul_div_cancel' hDM
    have hr'eq : D * r' = r := Nat.mul_div_cancel' hDr
    have hr'M' : r'.Coprime M' := (Nat.coprime_div_gcd_div_gcd hDpos).symm
    obtain ⟨a, ha⟩ := hne
    let t := (d * a)⁻¹ * u
    have ht : t ∈ Subgroup.closure S :=
      Subgroup.mul_mem _ (Subgroup.inv_mem _ (Subgroup.mul_mem _ hd (Subgroup.subset_closure ha))) hu
    obtain ⟨x, hx, hxM, hxt⟩ := generatorCore_realizes_coprime hL (NeZero.pos M') S ht
    let y := ((ZMod.unitOfCoprime x hxM)⁻¹ : (ZMod M')ˣ)
    let z := (y : ZMod M').val * r'
    have hz : z.Coprime M' := (ZMod.val_coe_unit_coprime y).mul_left hr'M'
    have hxz : x * z ≡ r' [MOD M'] := by
      rw [← ZMod.natCast_eq_natCast_iff, Nat.cast_mul]
      dsimp [z]
      rw [Nat.cast_mul, ZMod.natCast_zmod_val, ← mul_assoc, ← ZMod.coe_unitOfCoprime x hxM,
        ← Units.val_mul]
      simp [y]
    have hcast : (D * x : ℕ) ≡ (d * t : (ZMod L)ˣ).val.val [MOD L] := by
      rw [← ZMod.natCast_eq_natCast_iff, Nat.cast_mul, ZMod.natCast_zmod_val, ← hdD, hxt]
      exact (Units.val_mul d t).symm
    have hDXL : (D * x).Coprime L := coprime_of_modEq hcast (ZMod.val_coe_unit_coprime (d * t))
    have hDa : (D * x) * (a : ZMod L).val ≡ (u : ZMod L).val [MOD L] := by
      rw [← ZMod.natCast_eq_natCast_iff, Nat.cast_mul, Nat.cast_mul, ZMod.natCast_zmod_val,
        ZMod.natCast_zmod_val, ← hdD, hxt, ← Units.val_mul, ← Units.val_mul]
      congr 1
      simp [t, mul_comm, mul_left_comm]
    have hgd : M'.gcd L ∣ M.gcd L :=
      Nat.dvd_gcd ((Nat.gcd_dvd_left M' L).trans (Nat.div_dvd_of_dvd hDM)) (Nat.gcd_dvd_right ..)
    have hzcompat : z ≡ (a : ZMod L).val [MOD M'.gcd L] := by
      apply Nat.ModEq.cancel_left_of_coprime
        ((hDXL.coprime_dvd_right (Nat.gcd_dvd_right M' L)).symm.gcd_eq_one)
      have heq := (hxz.of_dvd (Nat.gcd_dvd_left M' L)).mul_left D
      rw [← Nat.mul_assoc, hr'eq] at heq
      exact heq.trans ((hr.of_dvd hgd).trans (hDa.of_dvd (Nat.gcd_dvd_right M' L)).symm)
    obtain ⟨p, hp, hpz, hpa⟩ := exists_prime_with_two_residues (NeZero.pos M') hL hz a hzcompat
    refine ⟨D * (x * p), primeCore_mul hDC (primeCore_mul hx
      ((primeCore_prime hp).mpr ⟨a, ha, hpa.symm⟩)), ?_, ?_⟩
    · push_cast
      have hfin := (ZMod.natCast_eq_natCast_iff ..).mpr hDa
      push_cast at hfin
      rw [ZMod.natCast_zmod_val, ZMod.natCast_zmod_val] at hfin
      simpa [hpa, mul_assoc] using hfin
    · have hprod := ((hpz.mul_left x).trans hxz).mul_left' D
      simpa [hM'eq, hr'eq] using hprod

/-- The fixed-state form of the finite-images lemma for subgroup semigroups.
The gcd extracts the entire prime-power prefix, including saturated valuations
and the case of the zero residue. -/
theorem subgroup_finite_image {L M r : ℕ} (hL : 0 < L) (hM : 0 < M)
    (H : Subgroup (ZMod L)ˣ) {u : (ZMod L)ˣ} (hu : u ∈ H) :
    (∃ x, PrimeCore (HasSubgroupResidue H) x ∧
      (x : ZMod L) = (u : ZMod L) ∧ x ≡ r [MOD M]) ↔
    r ≡ (u : ZMod L).val [MOD M.gcd L] ∧
      ∀ p, p.Prime → p ∣ M → ¬ p ∣ L → ¬ HasSubgroupResidue H p → ¬ p ∣ r := by
  have hu' : u ∈ Subgroup.closure (H : Set (ZMod L)ˣ) := by simpa only [Subgroup.closure_eq] using hu
  exact generator_finite_image hL hM H ⟨1, H.one_mem⟩ hu'

private theorem exists_avoiding_two_affine_roots {p : ℕ} (hp : p.Prime) (hbig : 2 < p)
    (x y a b : ℤ) (ha : ¬ (p : ℤ) ∣ a) (hb : ¬ (p : ℤ) ∣ b) :
    ∃ t : ℕ, ¬ (p : ℤ) ∣ x + a * t ∧ ¬ (p : ℤ) ∣ y - b * t := by
  classical
  letI : Fact p.Prime := ⟨hp⟩
  let u : ZMod p := -(x : ZMod p) / (a : ZMod p)
  let v : ZMod p := (y : ZMod p) / (b : ZMod p)
  obtain ⟨z, hz⟩ : ∃ z : ZMod p, z ∉ ({u, v} : Finset (ZMod p)) := by
    by_contra! h
    have hcard := Finset.card_le_card (show Finset.univ ⊆ ({u, v} : Finset (ZMod p))
      from fun z _ => h z)
    rw [Finset.card_univ, ZMod.card] at hcard
    exact (not_le_of_gt hbig) (hcard.trans Finset.card_le_two)
  have haz : (a : ZMod p) ≠ 0 := fun h => ha ((ZMod.intCast_zmod_eq_zero_iff_dvd a p).mp h)
  have hbz : (b : ZMod p) ≠ 0 := fun h => hb ((ZMod.intCast_zmod_eq_zero_iff_dvd b p).mp h)
  have hcast : (z.val : ZMod p) = z := ZMod.natCast_zmod_val z
  refine ⟨z.val, ?_, ?_⟩
  · intro hd
    have he := (ZMod.intCast_zmod_eq_zero_iff_dvd _ p).mpr hd
    push_cast at he
    rw [hcast] at he
    have hzu : z = u := by
      apply (eq_div_iff haz).mpr
      rwa [mul_comm, eq_neg_iff_add_eq_zero, add_comm]
    exact hz (by simp [hzu])
  · intro hd
    have he := (ZMod.intCast_zmod_eq_zero_iff_dvd _ p).mpr hd
    push_cast at he
    rw [hcast] at he
    have hzv : z = v := by
      apply (eq_div_iff hbz).mpr
      rw [mul_comm]
      exact (sub_eq_zero.mp he).symm
    exact hz (by simp [hzv])

private theorem choose_line_parameter {M L : ℕ} (hM : 0 < M) (heven : 2 ∣ L)
    (x y a b : ℤ)
    (hgood : ∀ p, p.Prime → ¬ p ∣ L → ¬ (p : ℤ) ∣ a ∧ ¬ (p : ℤ) ∣ b) :
    ∃ t : ℕ, ∀ p, p.Prime → p ∣ M → ¬ p ∣ L →
      ¬ (p : ℤ) ∣ x + a * t ∧ ¬ (p : ℤ) ∣ y - b * t := by
  classical
  let good := M.primeFactors.filter (fun p => ¬ p ∣ L)
  have hs : ∀ p ∈ good, ∃ t : ℕ, ¬ (p : ℤ) ∣ x + a * t ∧ ¬ (p : ℤ) ∣ y - b * t := by
    intro p hp
    obtain ⟨hpM, hnL⟩ := Finset.mem_filter.mp hp
    have hprime := Nat.prime_of_mem_primeFactors hpM
    have hbig : 2 < p := by
      have hne : p ≠ 2 := fun h => hnL (h.symm ▸ heven)
      have := hprime.two_le
      omega
    exact exists_avoiding_two_affine_roots hprime hbig x y a b
      (hgood p hprime hnL).1 (hgood p hprime hnL).2
  let f (p : ℕ) := if h : p ∈ good then (hs p h).choose else 0
  have hf {p : ℕ} (hp : p ∈ good) :
      ¬ (p : ℤ) ∣ x + a * f p ∧ ¬ (p : ℤ) ∣ y - b * f p := by
    simpa [f, hp] using (hs p hp).choose_spec
  have hne : ∀ p ∈ good, p ≠ 0 := fun p hp =>
    (Nat.prime_of_mem_primeFactors (Finset.mem_filter.mp hp).1).ne_zero
  have hcop : Set.Pairwise (↑good : Set ℕ) Nat.Coprime := by
    intro p hp q hq hne
    exact (Nat.coprime_primes (Nat.prime_of_mem_primeFactors (Finset.mem_filter.mp hp).1)
      (Nat.prime_of_mem_primeFactors (Finset.mem_filter.mp hq).1)).mpr hne
  let t := Nat.chineseRemainderOfFinset f id good hne hcop
  refine ⟨t, fun p hp hd hnL => ?_⟩
  have hmem : p ∈ good := Finset.mem_filter.mpr
    ⟨Nat.mem_primeFactors.mpr ⟨hp, hd, hM.ne'⟩, hnL⟩
  have hc : (t : ℕ) ≡ f p [MOD p] := t.2 p hmem
  have hcZ := Int.natCast_modEq_iff.mpr hc
  exact ⟨fun h => (hf hmem).1 (((Int.ModEq.refl x).add (hcZ.mul_left a)).dvd_iff.mp h),
    fun h => (hf hmem).2 (((Int.ModEq.refl y).sub (hcZ.mul_left b)).dvd_iff.mp h)⟩

private theorem repr_with_reduced_residue {K : ℕ} (hK : 0 < K) (z : ℤ) :
    ∃ r : ℕ, (r : ℤ) ≡ z [ZMOD (K : ℤ)] := by
  refine ⟨(z % (K : ℤ)).toNat, ?_⟩
  rw [Int.toNat_of_nonneg (Int.emod_nonneg _ (by exact_mod_cast hK.ne'))]
  exact Int.mod_modEq z K

/-- Lift an admissible pair of states along its integer solution line, then
choose the line parameter to avoid both zero coordinates at every good prime. -/
private theorem lift_residue_pair {L M A B N u v : ℕ} (hL : 0 < L) (hM : 0 < M)
    (heven : 2 ∣ L) (hcop : A.Coprime B)
    (hA : ∀ p, p.Prime → p ∣ A → p ∣ L) (hB : ∀ p, p.Prime → p ∣ B → p ∣ L)
    (hroot : A * u + B * v ≡ N [MOD L]) :
    ∃ r s : ℕ, r ≡ u [MOD L] ∧ s ≡ v [MOD L] ∧ A * r + B * s ≡ N [MOD M] ∧
      ∀ p, p.Prime → p ∣ M → ¬ p ∣ L → ¬ p ∣ r ∧ ¬ p ∣ s := by
  have hz : (A : ℤ) * u + (B : ℤ) * v ≡ (N : ℤ) [ZMOD (L : ℤ)] := by
    exact_mod_cast (Int.natCast_modEq_iff.mpr hroot)
  obtain ⟨k, hk⟩ := hz.dvd
  obtain ⟨i, j, hij⟩ : IsCoprime (A : ℤ) (B : ℤ) := hcop.isCoprime
  let x : ℤ := u + L * k * i
  let y : ℤ := v + L * k * j
  have hx : x ≡ (u : ℤ) [ZMOD (L : ℤ)] := Int.ModEq.symm
    (Int.modEq_iff_add_fac.mpr ⟨k * i, by dsimp [x]; ring⟩)
  have hy : y ≡ (v : ℤ) [ZMOD (L : ℤ)] := Int.ModEq.symm
    (Int.modEq_iff_add_fac.mpr ⟨k * j, by dsimp [y]; ring⟩)
  have hxy : (A : ℤ) * x + (B : ℤ) * y = N := by
    dsimp [x, y]
    linear_combination (L : ℤ) * k * hij - hk
  have hgood (p : ℕ) (hp : p.Prime) (hl : ¬ p ∣ L) :
      ¬ (p : ℤ) ∣ (B : ℤ) * L ∧ ¬ (p : ℤ) ∣ (A : ℤ) * L := by
    constructor
    · intro hd
      exact (hp.dvd_mul.mp (by exact_mod_cast hd)).elim (fun h => hl (hB p hp h)) hl
    · intro hd
      exact (hp.dvd_mul.mp (by exact_mod_cast hd)).elim (fun h => hl (hA p hp h)) hl
  obtain ⟨t, ht⟩ := choose_line_parameter hM heven x y ((B : ℤ) * L) ((A : ℤ) * L) hgood
  let X := x + (B : ℤ) * L * t
  let Y := y - (A : ℤ) * L * t
  let K := M.lcm L
  obtain ⟨r, hr⟩ := repr_with_reduced_residue (Nat.lcm_pos hM hL) X
  obtain ⟨s, hs⟩ := repr_with_reduced_residue (Nat.lcm_pos hM hL) Y
  have hKL : (L : ℤ) ∣ (K : ℤ) := by exact_mod_cast Nat.dvd_lcm_right M L
  have hKM : (M : ℤ) ∣ (K : ℤ) := by exact_mod_cast Nat.dvd_lcm_left M L
  have hXL : X ≡ x [ZMOD (L : ℤ)] := Int.ModEq.symm
    (Int.modEq_iff_add_fac.mpr ⟨(B : ℤ) * t, by dsimp [X]; ring⟩)
  have hYL : Y ≡ y [ZMOD (L : ℤ)] := Int.ModEq.symm
    (Int.modEq_iff_add_fac.mpr ⟨-(A : ℤ) * t, by dsimp [Y]; ring⟩)
  have hsum : (A : ℤ) * X + (B : ℤ) * Y = N := by
    dsimp [X, Y]
    linear_combination hxy
  refine ⟨r, s, Int.natCast_modEq_iff.mp ((hr.of_dvd hKL).trans (hXL.trans hx)),
    Int.natCast_modEq_iff.mp ((hs.of_dvd hKL).trans (hYL.trans hy)), ?_, ?_⟩
  · apply Int.natCast_modEq_iff.mp
    push_cast
    exact (((hr.of_dvd hKM).mul_left _).add ((hs.of_dvd hKM).mul_left _)).trans
      (hsum ▸ Int.ModEq.refl _)
  · intro p hp hpM hnL
    have hpK : (p : ℤ) ∣ (K : ℤ) := by exact_mod_cast hpM.trans (Nat.dvd_lcm_left M L)
    exact ⟨fun hd => (ht p hp hpM hnL).1 ((hr.of_dvd hpK).dvd_iff.mp (by exact_mod_cast hd)),
      fun hd => (ht p hp hpM hnL).2 ((hs.of_dvd hpK).dvd_iff.mp (by exact_mod_cast hd))⟩

/-- Local solubility for arbitrary nonempty sets of allowed unit classes, without any asymptotic
hypothesis. The left side is the existence of a pair counted by the joint
residue factor. -/
theorem generator_sum_local_iff {L A B N : ℕ} (hL : 0 < L) (heven : 2 ∣ L)
    (hcop : A.Coprime B) (hA : ∀ p, p.Prime → p ∣ A → p ∣ L)
    (hB : ∀ p, p.Prime → p ∣ B → p ∣ L)
    (S₁ S₂ U V : Set (ZMod L)ˣ) (h₁ : S₁.Nonempty) (h₂ : S₂.Nonempty)
    (hU : U ⊆ Subgroup.closure S₁) (hV : V ⊆ Subgroup.closure S₂) :
    (∃ u ∈ U, ∃ v ∈ V, (A : ZMod L) * (u : ZMod L) + (B : ZMod L) * (v : ZMod L) = N) ↔
    ∀ M > 0, ∃ x y, PrimeCore (HasGeneratorResidue S₁) x ∧ PrimeCore (HasGeneratorResidue S₂) y ∧
      (∃ u ∈ U, (x : ZMod L) = (u : ZMod L)) ∧ (∃ v ∈ V, (y : ZMod L) = (v : ZMod L)) ∧
      A * x + B * y ≡ N [MOD M] := by
  letI : NeZero L := ⟨hL.ne'⟩
  constructor
  · rintro ⟨u, hu, v, hv, huv⟩ M hM
    have hroot : A * (u : ZMod L).val + B * (v : ZMod L).val ≡ N [MOD L] := by
      rw [← ZMod.natCast_eq_natCast_iff]
      simpa only [Nat.cast_add, Nat.cast_mul, ZMod.natCast_zmod_val] using huv
    obtain ⟨r, s, huL, hvL, hrs, hgood⟩ :=
      lift_residue_pair hL hM heven hcop hA hB hroot
    obtain ⟨x, hx, hxu, hxr⟩ := (generator_finite_image hL hM S₁ h₁ (hU hu)).mpr
      ⟨huL.of_dvd (Nat.gcd_dvd_right ..), fun p hp hpd hnL _ => (hgood p hp hpd hnL).1⟩
    obtain ⟨y, hy, hyv, hys⟩ := (generator_finite_image hL hM S₂ h₂ (hV hv)).mpr
      ⟨hvL.of_dvd (Nat.gcd_dvd_right ..), fun p hp hpd hnL _ => (hgood p hp hpd hnL).2⟩
    exact ⟨x, y, hx, hy, ⟨u, hu, hxu⟩, ⟨v, hv, hyv⟩, ((hxr.mul_left A).add (hys.mul_left B)).trans hrs⟩
  · intro h
    obtain ⟨x, y, _, _, ⟨u, hu, hxu⟩, ⟨v, hv, hyv⟩, hxy⟩ := h L hL
    refine ⟨u, hu, v, hv, ?_⟩
    rw [← hxu, ← hyv, ← Nat.cast_mul, ← Nat.cast_mul, ← Nat.cast_add]
    exact (ZMod.natCast_eq_natCast_iff ..).mpr hxy

/-- The local equivalence of the subgroup-sum theorem, without any asymptotic
hypothesis. The left side is the existence of a pair counted by the joint
residue factor. -/
theorem subgroup_sum_local_iff {L A B N : ℕ} (hL : 0 < L) (heven : 2 ∣ L)
    (hcop : A.Coprime B) (hA : ∀ p, p.Prime → p ∣ A → p ∣ L)
    (hB : ∀ p, p.Prime → p ∣ B → p ∣ L)
    (H₁ H₂ : Subgroup (ZMod L)ˣ) (U V : Set (ZMod L)ˣ)
    (hU : U ⊆ H₁) (hV : V ⊆ H₂) :
    (∃ u ∈ U, ∃ v ∈ V, (A : ZMod L) * (u : ZMod L) + (B : ZMod L) * (v : ZMod L) = N) ↔
    ∀ M > 0, ∃ x y, PrimeCore (HasSubgroupResidue H₁) x ∧ PrimeCore (HasSubgroupResidue H₂) y ∧
      (∃ u ∈ U, (x : ZMod L) = (u : ZMod L)) ∧ (∃ v ∈ V, (y : ZMod L) = (v : ZMod L)) ∧
      A * x + B * y ≡ N [MOD M] := by
  exact generator_sum_local_iff hL heven hcop hA hB H₁ H₂ U V
    ⟨1, H₁.one_mem⟩ ⟨1, H₂.one_mem⟩
    (by simpa only [Subgroup.closure_eq] using hU)
    (by simpa only [Subgroup.closure_eq] using hV)

end QuadraticNorms
