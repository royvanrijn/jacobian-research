/-
Copyright (c) 2026 Roy van Rijn. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Roy van Rijn
-/
import FiniteEtaleKeller.PageOneTheorem
import Mathlib.RingTheory.Etale.Field
import Mathlib.RingTheory.Ideal.Quotient.Operations
import Mathlib.RingTheory.Trace.Basic

/-!
# Abstract finite étale algebras as Keller fibers

The polynomial-presentation theorem realizes `AdjoinRoot P` as the literal
fiber of a determinant-one polynomial map.  This file supplies the missing
monogenicity bridge: over a characteristic-zero field, every finite étale
algebra has a power basis and hence a squarefree polynomial presentation.

For a product of finite separable fields, primitive elements in the factors
are translated so that their traces are pairwise distinct.  Their minimal
polynomials are then pairwise coprime, and the Chinese remainder theorem shows
that the resulting tuple generates the whole product.
-/

noncomputable section

open Polynomial Function

namespace FiniteEtaleKeller

universe u v

section ProductMonogenicity

variable {K : Type*} [Field K] [CharZero K]
variable {I : Type u} [Finite I]
variable (L : I → Type v)
variable [∀ i, Field (L i)] [∀ i, Algebra K (L i)]
variable [∀ i, Module.Finite K (L i)]
variable [∀ i, Algebra.IsSeparable K (L i)]

private noncomputable def factorPowerBasis (i : I) : PowerBasis K (L i) :=
  Field.powerBasisOfFiniteOfSeparable K (L i)

private noncomputable def traceLabel (i : I) : K := by
  letI := Fintype.ofFinite I
  exact ((Fintype.equivFin I i : ℕ) : K)

private theorem traceLabel_injective :
    Function.Injective (traceLabel (K := K) (I := I)) := by
  classical
  letI := Fintype.ofFinite I
  intro i j hij
  change
    (Nat.castEmbedding (R := K) (Fintype.equivFin I i) =
      Nat.castEmbedding (R := K) (Fintype.equivFin I j)) at hij
  apply (Fintype.equivFin I).injective
  exact Fin.ext ((Nat.castEmbedding (R := K)).injective hij)

private noncomputable def factorShift (i : I) : K :=
  (traceLabel (K := K) (I := I) i -
      Algebra.trace K (L i) (factorPowerBasis (K := K) L i).gen) /
    (Module.finrank K (L i) : K)

private noncomputable def productGenerator : ∀ i, L i :=
  fun i =>
    (factorPowerBasis (K := K) L i).gen +
      algebraMap K (L i) (factorShift (K := K) L i)

omit [Finite I] [∀ i, Algebra.IsSeparable K (L i)] in
private theorem finrank_cast_ne_zero (i : I) :
    (Module.finrank K (L i) : K) ≠ 0 := by
  exact_mod_cast (Module.finrank_pos (R := K) (M := L i)).ne'

private theorem trace_productGenerator (i : I) :
    Algebra.trace K (L i) (productGenerator (K := K) L i) =
      traceLabel (K := K) (I := I) i := by
  rw [productGenerator, map_add, Algebra.trace_algebraMap]
  simp only [nsmul_eq_mul]
  dsimp [factorShift]
  field_simp [finrank_cast_ne_zero (K := K) L i]
  ring

omit [CharZero K] in
private theorem factorGenerator_adjoin_eq_top (i : I) :
    Algebra.adjoin K {productGenerator (K := K) L i} = ⊤ := by
  apply (factorPowerBasis (K := K) L i).adjoin_eq_top_of_gen_mem_adjoin
  have hx :
      productGenerator (K := K) L i ∈
        Algebra.adjoin K {productGenerator (K := K) L i} :=
    Algebra.subset_adjoin (Set.mem_singleton _)
  have hc :
      algebraMap K (L i) (factorShift (K := K) L i) ∈
        Algebra.adjoin K {productGenerator (K := K) L i} :=
    (Algebra.adjoin K {productGenerator (K := K) L i}).algebraMap_mem _
  simpa [productGenerator] using
    (Algebra.adjoin K {productGenerator (K := K) L i}).sub_mem hx hc

private noncomputable def shiftedFactorPowerBasis (i : I) : PowerBasis K (L i) :=
  PowerBasis.ofAdjoinEqTop'
    (Algebra.IsIntegral.isIntegral
      (R := K) (productGenerator (K := K) L i))
    (factorGenerator_adjoin_eq_top (K := K) L i)

omit [CharZero K] in
private theorem shiftedFactorPowerBasis_gen (i : I) :
    (shiftedFactorPowerBasis (K := K) L i).gen =
      productGenerator (K := K) L i := by
  simp [shiftedFactorPowerBasis]

private noncomputable def factorPolynomial (i : I) : K[X] :=
  minpoly K (productGenerator (K := K) L i)

omit [CharZero K] in
private theorem factorPolynomial_monic (i : I) :
    (factorPolynomial (K := K) L i).Monic :=
  minpoly.monic (Algebra.IsIntegral.isIntegral _)

omit [CharZero K] in
private theorem factorPolynomial_irreducible (i : I) :
    Irreducible (factorPolynomial (K := K) L i) :=
  minpoly.irreducible (Algebra.IsIntegral.isIntegral _)

private theorem factorPolynomial_pairwise_ne :
    Pairwise fun i j =>
      factorPolynomial (K := K) L i ≠ factorPolynomial (K := K) L j := by
  intro i j hij hp
  apply hij
  apply traceLabel_injective (K := K) (I := I)
  calc
    traceLabel (K := K) (I := I) i =
        Algebra.trace K (L i) (productGenerator (K := K) L i) :=
      (trace_productGenerator (K := K) L i).symm
    _ = -(factorPolynomial (K := K) L i).nextCoeff := by
      simpa [factorPolynomial, shiftedFactorPowerBasis_gen (K := K) L i] using
        (shiftedFactorPowerBasis (K := K) L i).trace_gen_eq_nextCoeff_minpoly
    _ = -(factorPolynomial (K := K) L j).nextCoeff := by rw [hp]
    _ = Algebra.trace K (L j) (productGenerator (K := K) L j) := by
      simpa [factorPolynomial, shiftedFactorPowerBasis_gen (K := K) L j] using
        (shiftedFactorPowerBasis (K := K) L j).trace_gen_eq_nextCoeff_minpoly.symm
    _ = traceLabel (K := K) (I := I) j :=
      trace_productGenerator (K := K) L j

private theorem factorPolynomial_pairwise_coprime :
    Pairwise (IsCoprime on factorPolynomial (K := K) L) := by
  intro i j hij
  change IsCoprime
    (factorPolynomial (K := K) L i)
    (factorPolynomial (K := K) L j)
  rw [(factorPolynomial_irreducible (K := K) L i).coprime_iff_not_dvd]
  intro hdvd
  have hassoc :
      Associated (factorPolynomial (K := K) L i)
        (factorPolynomial (K := K) L j) :=
    (factorPolynomial_irreducible (K := K) L i).associated_of_dvd
      (factorPolynomial_irreducible (K := K) L j) hdvd
  exact factorPolynomial_pairwise_ne (K := K) L hij
    (eq_of_monic_of_associated
      (factorPolynomial_monic (K := K) L i)
      (factorPolynomial_monic (K := K) L j) hassoc)

private theorem productGenerator_surjective :
    Function.Surjective
      (aeval (productGenerator (K := K) L) :
        K[X] →ₐ[K] (∀ i, L i)) := by
  classical
  letI := Fintype.ofFinite I
  intro y
  choose q hq using fun i =>
    (shiftedFactorPowerBasis (K := K) L i).exists_eq_aeval' (y i)
  let J : I → Ideal K[X] :=
    fun i => Ideal.span {factorPolynomial (K := K) L i}
  have hJ : Pairwise (IsCoprime on J) := by
    intro i j hij
    change IsCoprime (J i) (J j)
    rw [Ideal.isCoprime_span_singleton_iff]
    exact factorPolynomial_pairwise_coprime (K := K) L hij
  obtain ⟨qAll, hqAll⟩ :=
    Ideal.pi_quotient_surjective hJ
      (fun i => Ideal.Quotient.mk (J i) (q i))
  refine ⟨qAll, funext fun i => ?_⟩
  have hmem : qAll - q i ∈ J i :=
    (Submodule.Quotient.eq (J i)).mp (hqAll i)
  have hker :
      J i =
        RingHom.ker
          (aeval (R := K)
            (productGenerator (K := K) L i)).toRingHom := by
    dsimp [J, factorPolynomial]
    exact
      (minpoly.ker_aeval_eq_span_minpoly
        (A := K) (x := productGenerator (K := K) L i)).symm
  have heval :
      aeval (productGenerator (K := K) L i) qAll =
        aeval (productGenerator (K := K) L i) (q i) := by
    have hmem' :
        qAll - q i ∈
          RingHom.ker
            (aeval (R := K)
              (productGenerator (K := K) L i)).toRingHom := by
      rw [← hker]
      exact hmem
    rw [← sub_eq_zero, ← map_sub]
    exact hmem'
  rw [show
      aeval (productGenerator (K := K) L) qAll i =
        aeval (productGenerator (K := K) L i) qAll by
      simpa using
        (Polynomial.aeval_algHom_apply
          (Pi.evalAlgHom K L i) (productGenerator (K := K) L) qAll).symm]
  rw [heval, ← shiftedFactorPowerBasis_gen (K := K) L i, ← hq i]

private theorem productGenerator_adjoin_eq_top :
    Algebra.adjoin K {productGenerator (K := K) L} = ⊤ := by
  rw [Algebra.adjoin_singleton_eq_range_aeval, AlgHom.range_eq_top]
  exact productGenerator_surjective (K := K) L

/-- A finite product of finite separable field extensions over a
characteristic-zero field has a power basis. -/
noncomputable def productPowerBasis : PowerBasis K (∀ i, L i) := by
  classical
  letI := Fintype.ofFinite I
  exact
    (IsAdjoinRootMonic.mkOfAdjoinEqTop'
      (productGenerator_adjoin_eq_top (K := K) L)).powerBasis

private theorem productPowerBasis_gen :
    (productPowerBasis (K := K) L).gen =
      productGenerator (K := K) L := by
  simp [productPowerBasis]

private noncomputable def productPolynomial : K[X] := by
  classical
  letI := Fintype.ofFinite I
  exact ∏ i, factorPolynomial (K := K) L i

omit [CharZero K] in
private theorem productPolynomial_aeval :
    aeval (productGenerator (K := K) L)
        (productPolynomial (K := K) L) = 0 := by
  classical
  letI := Fintype.ofFinite I
  ext i
  rw [show
      aeval (productGenerator (K := K) L)
          (productPolynomial (K := K) L) i =
        aeval (productGenerator (K := K) L i)
          (productPolynomial (K := K) L) by
      simpa using
        (Polynomial.aeval_algHom_apply
          (Pi.evalAlgHom K L i) (productGenerator (K := K) L)
          (productPolynomial (K := K) L)).symm]
  simp only [productPolynomial, map_prod]
  exact Finset.prod_eq_zero (Finset.mem_univ i)
    (minpoly.aeval K (productGenerator (K := K) L i))

private theorem productPowerBasis_minpoly_eq :
    minpoly K (productPowerBasis (K := K) L).gen =
      productPolynomial (K := K) L := by
  classical
  letI := Fintype.ofFinite I
  let Q := productPolynomial (K := K) L
  have hQmonic : Q.Monic :=
    Polynomial.monic_prod_of_monic Finset.univ
      (factorPolynomial (K := K) L) fun i _ =>
        factorPolynomial_monic (K := K) L i
  have hroot :
      aeval (productPowerBasis (K := K) L).gen Q = 0 := by
    simpa [Q, productPowerBasis_gen (K := K) L] using
      productPolynomial_aeval (K := K) L
  have hdvd :
      minpoly K (productPowerBasis (K := K) L).gen ∣ Q :=
    minpoly.dvd K _ hroot
  symm
  apply eq_of_monic_of_dvd_of_natDegree_le
    (minpoly.monic (Algebra.IsIntegral.isIntegral _))
    hQmonic
    hdvd
  rw [(productPowerBasis (K := K) L).natDegree_minpoly,
    ← (productPowerBasis (K := K) L).finrank,
    Module.finrank_pi_fintype]
  dsimp [Q, productPolynomial]
  rw [Polynomial.natDegree_prod_of_monic]
  · apply Finset.sum_le_sum
    intro i _
    rw [factorPolynomial,
      ← shiftedFactorPowerBasis_gen (K := K) L i,
      (shiftedFactorPowerBasis (K := K) L i).natDegree_minpoly,
      ← (shiftedFactorPowerBasis (K := K) L i).finrank]
  · intro i _
    exact factorPolynomial_monic (K := K) L i

/-- The minimal polynomial of the chosen product generator is separable. -/
theorem productPowerBasis_minpoly_separable :
    (minpoly K (productPowerBasis (K := K) L).gen).Separable := by
  classical
  letI := Fintype.ofFinite I
  rw [productPowerBasis_minpoly_eq (K := K) L]
  exact Polynomial.separable_prod
    (factorPolynomial_pairwise_coprime (K := K) L)
    (fun i => (factorPolynomial_irreducible (K := K) L i).separable)

end ProductMonogenicity

section AbstractMonogenicity

variable {K : Type*} [Field K] [CharZero K]
variable (A : Type u) [CommRing A] [Algebra K A] [Algebra.Etale K A]

private theorem exists_finiteEtalePowerBasis :
    ∃ pb : PowerBasis K A, (minpoly K pb.gen).Separable := by
  classical
  obtain ⟨I, _, L, _, _, e, hL⟩ :=
    (Algebra.Etale.iff_exists_algEquiv_prod K A).mp inferInstance
  letI (i : I) : Module.Finite K (L i) := (hL i).1
  letI (i : I) : Algebra.IsSeparable K (L i) := (hL i).2
  let pb := productPowerBasis (K := K) L
  let pbA := pb.map e.symm
  refine ⟨pbA, ?_⟩
  have hmin :
      minpoly K pbA.gen = minpoly K pb.gen := by
    rw [← pbA.minpolyGen_eq, ← pb.minpolyGen_eq]
    exact PowerBasis.minpolyGen_map pb e.symm
  rw [hmin]
  exact productPowerBasis_minpoly_separable (K := K) L

/-- Every finite étale algebra over a characteristic-zero field is
monogenic, expressed as the existence of a power basis. -/
noncomputable def finiteEtalePowerBasis : PowerBasis K A :=
  (exists_finiteEtalePowerBasis (K := K) A).choose

/-- The polynomial selected from the noncanonical finite-étale power basis. -/
noncomputable def finiteEtalePolynomial : K[X] :=
  minpoly K (finiteEtalePowerBasis (K := K) A).gen

/-- The selected polynomial is separable. -/
theorem finiteEtalePolynomial_separable :
    (finiteEtalePolynomial (K := K) A).Separable :=
  (exists_finiteEtalePowerBasis (K := K) A).choose_spec

/-- The selected polynomial is squarefree. -/
theorem finiteEtalePolynomial_squarefree :
    Squarefree (finiteEtalePolynomial (K := K) A) :=
  (finiteEtalePolynomial_separable (K := K) A).squarefree

/-- The selected polynomial presentation of an abstract finite étale
algebra. -/
noncomputable def finiteEtalePresentation :
    AdjoinRoot (finiteEtalePolynomial (K := K) A) ≃ₐ[K] A :=
  AdjoinRoot.equiv'
    (finiteEtalePolynomial (K := K) A)
    (finiteEtalePowerBasis (K := K) A)
    (by
      rw [AdjoinRoot.aeval_eq]
      exact AdjoinRoot.mk_self)
    (minpoly.aeval K (finiteEtalePowerBasis (K := K) A).gen)

/-- The degree of the selected polynomial is the rank of the abstract
finite étale algebra. -/
theorem finiteEtalePolynomial_natDegree :
    (finiteEtalePolynomial (K := K) A).natDegree =
      Module.finrank K A := by
  rw [finiteEtalePolynomial,
    (finiteEtalePowerBasis (K := K) A).natDegree_minpoly,
    ← (finiteEtalePowerBasis (K := K) A).finrank]

variable {R S : Type*}
variable [CommRing R] [Algebra K R]
variable [CommRing S] [Algebra K S]

/-- Precomposition with the selected polynomial presentation identifies maps
out of the abstract finite étale algebra with maps out of its quotient
presentation. -/
noncomputable def finiteEtaleHomEquiv :
    (A →ₐ[K] R) ≃
      (AdjoinRoot (finiteEtalePolynomial (K := K) A) →ₐ[K] R) where
  toFun φ := φ.comp (finiteEtalePresentation (K := K) A).toAlgHom
  invFun ψ := ψ.comp (finiteEtalePresentation (K := K) A).symm.toAlgHom
  left_inv φ := by
    ext x
    simp
  right_inv ψ := by
    ext
    simp

/-- Naturality of the change from the abstract algebra to its selected
polynomial presentation. -/
theorem finiteEtaleHomEquiv_natural
    (f : R →ₐ[K] S) (φ : A →ₐ[K] R) :
    finiteEtaleHomEquiv (K := K) (A := A) (R := S) (f.comp φ) =
      f.comp (finiteEtaleHomEquiv (K := K) (A := A) (R := R) φ) := by
  rfl

end AbstractMonogenicity

section AbstractRealization

variable {K : Type*} [Field K] [CharZero K]
variable (A : Type u) [CommRing A] [Algebra K A] [Algebra.Etale K A]

/-- The three target coordinates of the automatically selected realization. -/
def abstractFiniteEtaleTarget
    (P : K[X]) (hdeg : 3 ≤ P.natDegree) : Fin 3 → K
  | 0 => 1
  | 1 => 0
  | 2 => automaticRealizationTargetC P hdeg

set_option maxHeartbeats 800000 in
-- Elaborating the universe-polymorphic naturality fields exceeds the default.
/-- The public data package for realizing an abstract finite étale algebra as
the full fiber of a Keller map.  It explicitly exposes the polynomial
presentation, map, target, natural represented-fiber equivalence, geometric
degree, and coordinate degree bound. -/
structure AbstractFiniteEtalePageOneCertificate
    (hrank : 3 ≤ Module.finrank K A) where
  P : K[X]
  squarefree : Squarefree P
  degreeAtLeastThree : 3 ≤ P.natDegree
  presentation : AdjoinRoot P ≃ₐ[K] A
  map : Fin 3 → GaugePolynomial K
  target : Fin 3 → K
  map_eq : map = automaticRealizationMap P degreeAtLeastThree
  target_eq : target = abstractFiniteEtaleTarget (K := K) P degreeAtLeastThree
  jacobian : jacobianDet map = 1
  geometricDegree :
    automaticRealizationGeometricDegree P degreeAtLeastThree =
      Module.finrank K A
  fiberEquiv :
    ∀ (R : Type v) [CommRing R] [Algebra K R],
      (A →ₐ[K] R) ≃
        GeneralGaugeJacobianOneFiberPoint
          (realizationSeed P
            (chosenAdmissibleTranslation P degreeAtLeastThree))
          1 (target 2) R
  fiber_natural :
    ∀ {R S : Type v} [CommRing R] [Algebra K R]
      [CommRing S] [Algebra K S]
      (f : R →ₐ[K] S) (φ : A →ₐ[K] R),
      GeneralGaugeJacobianOneFiberPoint.map f
          (fiberEquiv R φ) =
        fiberEquiv S (f.comp φ)
  rank : P.natDegree = Module.finrank K A
  degreeBound :
    ∀ i : Fin 3, (map i).totalDegree ≤ 6 * Module.finrank K A + 2

variable {R S : Type*}
variable [CommRing R] [Algebra K R]
variable [CommRing S] [Algebra K S]

/-- The natural equivalence between maps out of an abstract finite étale
algebra and points of the distinguished literal Keller fiber. -/
noncomputable def abstractFiniteEtaleFiberRepresentingEquiv
    (hrank : 3 ≤ Module.finrank K A) :
    (A →ₐ[K] R) ≃
      GeneralGaugeJacobianOneFiberPoint
        (realizationSeed (finiteEtalePolynomial (K := K) A)
          (chosenAdmissibleTranslation
            (finiteEtalePolynomial (K := K) A)
            (finiteEtalePolynomial_natDegree (K := K) A ▸ hrank)))
        1
        (automaticRealizationTargetC
          (finiteEtalePolynomial (K := K) A)
          (finiteEtalePolynomial_natDegree (K := K) A ▸ hrank)) R :=
  (finiteEtaleHomEquiv (K := K) (A := A) (R := R)).trans
    (automaticJacobianOneFiberRepresentingEquiv
      (A := R)
      (finiteEtalePolynomial (K := K) A)
      (finiteEtalePolynomial_squarefree (K := K) A)
      (finiteEtalePolynomial_natDegree (K := K) A ▸ hrank))

/-- The abstract represented-fiber equivalence commutes with every morphism
of commutative test algebras. -/
theorem abstractFiniteEtaleFiberRepresentingEquiv_natural
    (hrank : 3 ≤ Module.finrank K A)
    (f : R →ₐ[K] S) (φ : A →ₐ[K] R) :
    GeneralGaugeJacobianOneFiberPoint.map f
        (abstractFiniteEtaleFiberRepresentingEquiv
          (K := K) (A := A) (R := R) hrank φ) =
      abstractFiniteEtaleFiberRepresentingEquiv
        (K := K) (A := A) (R := S) hrank (f.comp φ) := by
  exact automaticJacobianOneFiberRepresentingEquiv_natural
    (finiteEtalePolynomial (K := K) A)
    (finiteEtalePolynomial_squarefree (K := K) A)
    (finiteEtalePolynomial_natDegree (K := K) A ▸ hrank)
    f
    (finiteEtaleHomEquiv (K := K) (A := A) (R := R) φ)

/-- Every abstract finite étale algebra of rank at least three over a
characteristic-zero field is naturally represented by the full distinguished
fiber of an explicitly selected determinant-one polynomial map.  The choice
of primitive element, and hence this construction, is noncomputable. -/
noncomputable def abstractFiniteEtale_pageOne
    (hrank : 3 ≤ Module.finrank K A) :
    AbstractFiniteEtalePageOneCertificate (K := K) A hrank := by
  let P := finiteEtalePolynomial (K := K) A
  have hP : Squarefree P :=
    finiteEtalePolynomial_squarefree (K := K) A
  have hrankP : P.natDegree = Module.finrank K A :=
    finiteEtalePolynomial_natDegree (K := K) A
  have hdeg : 3 ≤ P.natDegree := hrankP.symm ▸ hrank
  refine
    { P := P
      squarefree := hP
      degreeAtLeastThree := hdeg
      presentation := finiteEtalePresentation (K := K) A
      map := automaticRealizationMap P hdeg
      target := abstractFiniteEtaleTarget (K := K) P hdeg
      map_eq := rfl
      target_eq := rfl
      jacobian := automaticRealizationMap_jacobianDet P hdeg
      geometricDegree := ?_
      fiberEquiv := fun R _ _ =>
        abstractFiniteEtaleFiberRepresentingEquiv
          (K := K) (A := A) (R := R) hrank
      fiber_natural := ?_
      rank := hrankP
      degreeBound := ?_ }
  · exact (automaticRealizationGeometricDegree_eq P hdeg).trans hrankP
  · intro R S _ _ _ _ f φ
    exact abstractFiniteEtaleFiberRepresentingEquiv_natural
      (K := K) (A := A) hrank f φ
  · intro i
    simpa [hrankP] using automaticRealizationMap_totalDegree P hdeg i

#print axioms productPowerBasis_minpoly_separable
#print axioms finiteEtalePowerBasis
#print axioms finiteEtalePolynomial_squarefree
#print axioms finiteEtalePresentation
#print axioms abstractFiniteEtaleFiberRepresentingEquiv_natural
#print axioms abstractFiniteEtale_pageOne

end AbstractRealization

end FiniteEtaleKeller
