/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import Mathlib

/-!
# The degree-four fixed-point lemma

This file formalizes the finite-group lemma used in the paper's Chebotarev
proof of the degree-four local--global theorem:

> If a finite group acts on at most four points and every group element fixes
> a point, then the whole group fixes a point.

The proof uses orbit decomposition and Burnside's lemma.  If there is no
global fixed point, every orbit has at least two elements, so there are at
most two orbits.  In the one-orbit case, the identity makes the average
number of fixed points too large.  In the two-orbit case, both orbits have
size two; every nonempty fixed-point set therefore has at least two elements,
and the identity again makes the Burnside average too large.

The arithmetic passage from local points at almost all finite places to
fixed points of individual Galois elements remains the separate Chebotarev
input.
-/

noncomputable section

namespace FiniteEtaleKeller

open scoped BigOperators

/-- **Degree-four fixed-point lemma.** If a finite group acts on at most four
points and every group element fixes some point, then the whole group fixes
one point. -/
theorem degreeFour_fixedPoint
    (G Ω : Type*) [Group G] [Finite G] [Fintype Ω] [MulAction G Ω]
    (hcard : Fintype.card Ω ≤ 4)
    (hlocal : ∀ g : G, Set.Nonempty (MulAction.fixedBy Ω g)) :
    Set.Nonempty (MulAction.fixedPoints G Ω) := by
  classical
  letI : Fintype G := Fintype.ofFinite G
  by_contra hglobal
  obtain ⟨x₀, hx₀⟩ := hlocal (1 : G)
  let Q := MulAction.orbitRel.Quotient G Ω
  letI : Fintype Q := Fintype.ofFinite Q
  letI (x : Ω) : Fintype (MulAction.orbit G x) := Fintype.ofFinite _
  letI (q : Q) : Fintype q.orbit := Fintype.ofFinite _
  letI (g : G) : Fintype (MulAction.fixedBy Ω g) := Fintype.ofFinite _
  have hnotfixed : ∀ x : Ω, x ∉ MulAction.fixedPoints G Ω := by
    intro x hx
    exact hglobal ⟨x, hx⟩
  have horbit2 : ∀ x : Ω, 2 ≤ Fintype.card (MulAction.orbit G x) := by
    intro x
    have hne : Fintype.card (MulAction.orbit G x) ≠ 1 := by
      intro h
      exact hnotfixed x
        (MulAction.mem_fixedPoints_iff_card_orbit_eq_one.mpr h)
    have hpos : 0 < Fintype.card (MulAction.orbit G x) :=
      Fintype.card_pos_iff.mpr ⟨⟨x, MulAction.mem_orbit_self x⟩⟩
    omega
  have hdecomp :
      Fintype.card Ω = ∑ q : Q, Fintype.card q.orbit := by
    rw [← Fintype.card_sigma]
    exact Fintype.card_congr (MulAction.selfEquivSigmaOrbits' G Ω)
  have hqorbit2 : ∀ q : Q, 2 ≤ Fintype.card q.orbit := by
    intro q
    have hc :
        Fintype.card q.orbit =
          Fintype.card (MulAction.orbit G q.out) :=
      Fintype.card_congr (Equiv.setCongr
        (MulAction.orbitRel.Quotient.orbit_eq_orbit_out q Quotient.out_eq'))
    rw [hc]
    exact horbit2 q.out
  have htwor : 2 * Fintype.card Q ≤ Fintype.card Ω := by
    rw [hdecomp]
    calc
      2 * Fintype.card Q = ∑ _q : Q, 2 := by simp [Nat.mul_comm]
      _ ≤ ∑ q : Q, Fintype.card q.orbit :=
        Finset.sum_le_sum fun q _ => hqorbit2 q
  have hQpos : 0 < Fintype.card Q := by
    exact Fintype.card_pos_iff.mpr
      ⟨Quotient.mk (MulAction.orbitRel G Ω) x₀⟩
  have hQ : Fintype.card Q = 1 ∨ Fintype.card Q = 2 := by
    omega
  have hburn :=
    MulAction.sum_card_fixedBy_eq_card_orbits_mul_card_group G Ω
  rcases hQ with hQ | hQ
  · have hfixpos : ∀ g : G, 1 ≤ Fintype.card (MulAction.fixedBy Ω g) :=
      fun g => Fintype.card_pos_iff.mpr (hlocal g).to_subtype
    have hn2 : 2 ≤ Fintype.card Ω := by
      have hne : Fintype.card Ω ≠ 1 := by
        intro h
        haveI : Subsingleton Ω :=
          Fintype.card_le_one_iff_subsingleton.mp (by omega)
        exact hglobal ⟨x₀,
          MulAction.mem_fixedPoints.mpr fun g => Subsingleton.elim _ _⟩
      have hp : 0 < Fintype.card Ω := Fintype.card_pos_iff.mpr ⟨x₀⟩
      omega
    have hstrict : Fintype.card G <
        ∑ g : G, Fintype.card (MulAction.fixedBy Ω g) := by
      calc
        Fintype.card G = ∑ _g : G, 1 := by simp
        _ < ∑ g : G, Fintype.card (MulAction.fixedBy Ω g) := by
          apply Finset.sum_lt_sum
          · intro g _
            exact hfixpos g
          · exact ⟨1, Finset.mem_univ _, by
              simpa using (show 1 < Fintype.card Ω by omega)⟩
    rw [hburn, hQ, one_mul] at hstrict
    exact (Nat.lt_irrefl _ hstrict)
  · have hΩ4 : Fintype.card Ω = 4 := by omega
    have hqcard2 : ∀ q : Q, Fintype.card q.orbit = 2 := by
      intro q
      have hrest :
          2 * (Finset.univ.erase q).card ≤
            ∑ r ∈ Finset.univ.erase q, Fintype.card r.orbit := by
        calc
          2 * (Finset.univ.erase q).card =
              ∑ _r ∈ Finset.univ.erase q, 2 := by
                simp [Nat.mul_comm]
          _ ≤ ∑ r ∈ Finset.univ.erase q, Fintype.card r.orbit :=
            Finset.sum_le_sum fun r _ => hqorbit2 r
      have herase : (Finset.univ.erase q).card = 1 := by
        simp [hQ]
      have hsplit :
          Fintype.card Ω =
            Fintype.card q.orbit +
              ∑ r ∈ Finset.univ.erase q, Fintype.card r.orbit := by
        rw [hdecomp, ← Finset.add_sum_erase _ _ (Finset.mem_univ q)]
      have hle : Fintype.card q.orbit ≤ 2 := by omega
      exact Nat.le_antisymm hle (hqorbit2 q)
    have hfix2 : ∀ g : G, 2 ≤ Fintype.card (MulAction.fixedBy Ω g) := by
      intro g
      obtain ⟨x, hx⟩ := hlocal g
      let q : Q := Quotient.mk'' x
      have horb2 : Fintype.card (MulAction.orbit G x) = 2 := by
        calc
          Fintype.card (MulAction.orbit G x) = Fintype.card q.orbit :=
            Fintype.card_congr (Equiv.setCongr (by rfl))
          _ = 2 := hqcard2 q
      let xs : MulAction.orbit G x := ⟨x, MulAction.mem_orbit_self x⟩
      have hnat : Nat.card (MulAction.orbit G x) = 2 := by
        rw [Nat.card_eq_fintype_card]
        exact horb2
      obtain ⟨ys, hysne, hysunique⟩ :=
        (Nat.card_eq_two_iff' xs).mp hnat
      have hgys_mem : g • (ys : Ω) ∈ MulAction.orbit G x :=
        MulAction.mapsTo_smul_orbit g x ys.property
      let gys : MulAction.orbit G x := ⟨g • (ys : Ω), hgys_mem⟩
      have hgysne : gys ≠ xs := by
        intro heq
        apply hysne
        apply Subtype.ext
        apply smul_left_cancel g
        calc
          g • (ys : Ω) = x := congrArg Subtype.val heq
          _ = g • x := (MulAction.mem_fixedBy.mp hx).symm
      have hgys_eq : gys = ys := hysunique gys hgysne
      have hysfix : g • (ys : Ω) = (ys : Ω) :=
        congrArg Subtype.val hgys_eq
      have hone : 1 < Fintype.card (MulAction.fixedBy Ω g) := by
        rw [Fintype.one_lt_card_iff]
        exact ⟨⟨x, hx⟩, ⟨ys, MulAction.mem_fixedBy.mpr hysfix⟩,
          fun h => hysne (Subtype.ext (congrArg Subtype.val h).symm)⟩
      omega
    have hstrict : 2 * Fintype.card G <
        ∑ g : G, Fintype.card (MulAction.fixedBy Ω g) := by
      calc
        2 * Fintype.card G = ∑ _g : G, 2 := by simp [Nat.mul_comm]
        _ < ∑ g : G, Fintype.card (MulAction.fixedBy Ω g) := by
          apply Finset.sum_lt_sum
          · intro g _
            exact hfix2 g
          · exact ⟨1, Finset.mem_univ _, by
              simpa using (show 2 < Fintype.card Ω by omega)⟩
    rw [hburn, hQ] at hstrict
    exact Nat.lt_irrefl _ hstrict

#print axioms degreeFour_fixedPoint

end FiniteEtaleKeller
