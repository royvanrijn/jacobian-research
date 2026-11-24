import Mathlib.NumberTheory.DirichletCharacter.Bounds
import Mathlib.NumberTheory.MulChar.Lemmas

open scoped BigOperators

namespace QuadraticNorms

/-- Characters are equivalent on a progression if their quotient
descends to its modulus. Both characters are lifted to a common level. -/
def BankCharacterEquivalent (L : ℕ) {M : ℕ} (χ ψ : DirichletCharacter ℂ M) : Prop :=
  (χ * ψ⁻¹).conductor ∣ L

/-- The conductor criterion is equality on the kernel of the unit
reduction map. -/
theorem bank_character_equivalent_iff {L M : ℕ} [NeZero M] (h : L ∣ M)
    (χ ψ : DirichletCharacter ℂ M) :
    BankCharacterEquivalent L χ ψ ↔
      ∀ u : (ZMod M)ˣ, ZMod.unitsMap h u = 1 → χ u = ψ u := by
  rw [BankCharacterEquivalent, ← DirichletCharacter.mem_conductorSet_iff_conductor_dvd _ h,
    DirichletCharacter.mem_conductorSet_iff, DirichletCharacter.factorsThrough_iff_ker_unitsMap h]
  constructor
  · intro hker u hu
    have hc : (χ * ψ⁻¹).toUnitHom u = 1 := hker hu
    have he := congrArg (fun z : ℂˣ => (z : ℂ)) hc
    have hp : ψ u ≠ 0 := MulChar.apply_ne_zero_iff.mpr u.isUnit
    rw [MulChar.coe_toUnitHom, MulChar.mul_apply, MulChar.inv_apply_eq_inv', Units.val_one,
      ← div_eq_mul_inv, div_eq_iff hp, one_mul] at he
    exact he
  · intro he u hu
    change (χ * ψ⁻¹).toUnitHom u = 1
    apply Units.ext
    rw [MulChar.coe_toUnitHom, MulChar.mul_apply, MulChar.inv_apply_eq_inv', he u hu, Units.val_one]
    exact mul_inv_cancel₀ (MulChar.apply_ne_zero_iff.mpr u.isUnit)

theorem bank_character_equivalence {L M : ℕ} [NeZero M] (h : L ∣ M) :
    Equivalence (BankCharacterEquivalent L : DirichletCharacter ℂ M → _ → Prop) where
  refl χ := (bank_character_equivalent_iff h χ χ).mpr (fun _ _ => rfl)
  symm {χ ψ} he := (bank_character_equivalent_iff h ψ χ).mpr (fun u hu =>
    ((bank_character_equivalent_iff h χ ψ).mp he u hu).symm)
  trans {χ ψ θ} h₁ h₂ := (bank_character_equivalent_iff h χ θ).mpr (fun u hu =>
    ((bank_character_equivalent_iff h χ ψ).mp h₁ u hu).trans
      ((bank_character_equivalent_iff h ψ θ).mp h₂ u hu))

/-- Equivalent characters differ by a character of the bank modulus,
evaluated at the reduction of the unit. -/
theorem bank_character_equivalent_phase {L M : ℕ} [NeZero M] (h : L ∣ M)
    {χ ψ : DirichletCharacter ℂ M} (he : BankCharacterEquivalent L χ ψ) :
    ∃ θ : DirichletCharacter ℂ L, ∀ u : (ZMod M)ˣ, χ u = θ (ZMod.unitsMap h u) * ψ u := by
  have hfactor : (χ * ψ⁻¹).FactorsThrough L :=
    (DirichletCharacter.mem_conductorSet_iff_conductor_dvd _ h).mpr he
  obtain ⟨h', θ, hθ⟩ := hfactor
  refine ⟨θ, fun u => ?_⟩
  have heq := congrArg (fun η : DirichletCharacter ℂ M => η u) hθ
  change (χ * ψ⁻¹) u = DirichletCharacter.changeLevel h' θ u at heq
  rw [MulChar.mul_apply, MulChar.inv_apply_eq_inv',
    DirichletCharacter.changeLevel_eq_cast_of_dvd, ← ZMod.unitsMap_val h u, ← div_eq_mul_inv] at heq
  exact (div_eq_iff (MulChar.apply_ne_zero_iff.mpr u.isUnit)).mp heq

/-- For a product of labels all in the same progression, the phase is
the phase of one label raised to the number of slots. -/
theorem bank_character_product_phase {L M : ℕ} [NeZero M] (h : L ∣ M)
    {χ ψ : DirichletCharacter ℂ M} (he : BankCharacterEquivalent L χ ψ)
    (c : (ZMod L)ˣ) :
    ∃ ζ : ℂ, ‖ζ‖ = 1 ∧ ∀ {ι : Type*} (s : Finset ι) (u : ι → (ZMod M)ˣ),
      (∀ i ∈ s, ZMod.unitsMap h (u i) = c) →
      χ (↑(∏ i ∈ s, u i) : ZMod M) = ζ ^ s.card * ψ (↑(∏ i ∈ s, u i) : ZMod M) := by
  obtain ⟨θ, hθ⟩ := bank_character_equivalent_phase h he
  refine ⟨θ c, θ.unit_norm_eq_one c, fun s u hu => ?_⟩
  rw [hθ (∏ i ∈ s, u i), map_prod]
  have hmap : (∏ i ∈ s, ZMod.unitsMap h (u i)) = c ^ s.card := by
    rw [← Finset.prod_const]
    exact Finset.prod_congr rfl hu
  rw [hmap, Units.val_pow_eq_pow_val, map_pow]

/-- The quotient of two inequivalent characters has zero sum on
every unit residue class of the bank modulus. -/
theorem bank_inequivalent_character_fiber_sum_zero {L M : ℕ} [NeZero M] (h : L ∣ M)
    {χ ψ : DirichletCharacter ℂ M} (he : ¬ BankCharacterEquivalent L χ ψ)
    (c : (ZMod L)ˣ) :
    (∑ u ∈ Finset.univ.filter (fun u => ZMod.unitsMap h u = c), χ u * star (ψ u)) = 0 := by
  classical
  let ξ := χ * ψ⁻¹
  have hex : ∃ a : (ZMod M)ˣ, ZMod.unitsMap h a = 1 ∧ ξ a ≠ 1 := by
    have hker : ¬ (ZMod.unitsMap h).ker ≤ ξ.toUnitHom.ker := fun hc => he <|
      (DirichletCharacter.mem_conductorSet_iff_conductor_dvd _ h).mp <|
        (DirichletCharacter.factorsThrough_iff_ker_unitsMap h).mpr hc
    by_contra! hnone
    exact hker (fun a ha => Units.ext (by
      rw [MulChar.coe_toUnitHom, Units.val_one]
      exact hnone a ha))
  obtain ⟨a, ha, ha1⟩ := hex
  let s := Finset.univ.filter (fun u => ZMod.unitsMap h u = c)
  have ha' : ZMod.unitsMap h a⁻¹ = 1 := by rw [map_inv, ha, inv_one]
  have hperm : (∑ u ∈ s, ξ (↑(a * u) : ZMod M)) = ∑ u ∈ s, ξ u := by
    refine Finset.sum_bij' (fun u _ => a * u) (fun v _ => a⁻¹ * v) ?_ ?_ ?_ ?_ ?_
    · intro u hu
      have huc := (Finset.mem_filter.mp hu).2
      exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, by rw [map_mul, ha, huc, one_mul]⟩
    · intro v hv
      have hvc := (Finset.mem_filter.mp hv).2
      exact Finset.mem_filter.mpr ⟨Finset.mem_univ _, by rw [map_mul, ha', hvc, one_mul]⟩
    · intro u _
      simp
    · intro u _
      simp
    · intro u _
      rfl
  have hmul : ξ a * (∑ u ∈ s, ξ u) = ∑ u ∈ s, ξ u := by
    rw [Finset.mul_sum, ← hperm]
    exact Finset.sum_congr rfl (fun u _ => by simp only [Units.val_mul, map_mul])
  have hzero := eq_zero_of_mul_eq_self_left ha1 hmul
  simpa only [ξ, s, MulChar.mul_apply, ← MulChar.star_apply'] using hzero

end QuadraticNorms
