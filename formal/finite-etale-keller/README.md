# Lean formalization: finite étale Keller fibers

This project formalizes the polynomial-presentation construction in
*Prescribed Finite Étale Algebras as Full Fibers of Keller Maps with Symmetric
Monodromy*. It uses Lean `v4.33.0-rc1` and Mathlib at the matching release
candidate.

## Proof status

| Stage | Scope | Status |
|---|---|---|
| 1 | Explicit quintic map, output scaling, and Bézout inverse | implemented |
| 2 | Universal marked-line identities and Jacobian cancellation core | implemented |
| 3 | Two-sided source/chart reconstruction over arbitrary commutative rings | implemented |
| 4 | Roots versus full abstract source-fiber points, including naturality | implemented |
| 5 | Representation by `K[S]/(E)` and transport to `K[T]/(P)` | implemented |
| 6 | Existence and automatic choice of an admissible translation | implemented |
| 7 | Complete finite-sum all-degree gauge assembly identities | implemented |
| 8 | One arbitrary-degree `MvPolynomial (Fin 3) K` map and exact evaluations | implemented |
| 9 | Generic inverse polynomial, explicit `β`, and derivative factorization | implemented |
| 10 | Actual displayed equations versus represented source equations | implemented |
| 11 | Literal raw map fiber and naturality over arbitrary test algebras | implemented |
| 12 | Determinant-one normalization of the literal fiber | implemented |
| 13 | General arbitrary-degree Jacobian `-2`, normalized Jacobian `1` | implemented |
| 14 | General and final `6N+2` coordinate-degree bounds | implemented |
| 15 | Final automatic actual-map realization from squarefree `P` | implemented |
| 16 | Exact represented-fiber rank and nontriviality | implemented |
| 17 | Exact displayed quintic literal fiber, naturality, and rank five | implemented |
| 18 | Fixed-`π,b` inverse irreducibility and exact degree over `K(C)` | implemented |
| 19 | Full independent-parameter inverse irreducibility and root-extension finrank over `K(Π,B)(C)` | implemented |
| 20 | Explicit quintic: no rational point, a real point, and a three-adic point | implemented |
| 21 | Finite étaleness of every squarefree quotient and of the explicit quintic quotient | implemented |
| 22 | Coordinate algebraic independence and injective function-field pullback | implemented |
| 23 | Supplied-parameter map and target compatibility under scalar extension | implemented |
| 24 | Tensor-product base change of the representing polynomial quotient | implemented |
| 25 | Comparison with the actual pullback-field extension and geometric degree | implemented |
| 26 | Monogenicity and the passage from arbitrary finite étale algebras | paper proof; not yet Lean |
| 27 | Historical degree-two Galois exclusion | external theorem; not yet Lean |

## Final polynomial-presentation theorem

Let `K` be a characteristic-zero field, let `P : K[X]` be squarefree, and
assume `3 ≤ P.natDegree`. Lean now chooses an admissible translation parameter
internally and defines

```text
automaticRealizationMap P hdeg : Fin 3 → MvPolynomial (Fin 3) K
```

together with its distinguished target. The formal development proves:

```text
jacobianDet (automaticRealizationMap P hdeg) = 1
```

and, for every coordinate `i`,

```text
(automaticRealizationMap P hdeg i).totalDegree ≤ 6 * P.natDegree + 2.
```

The new function-field comparison proves

```text
automaticRealizationGeometricDegree P hdeg = P.natDegree
```

and `automaticRealization_pageOne` bundles this equality with determinant
one and all of the fiber assertions below.

For every commutative test `K`-algebra `A`, it constructs an equivalence

```text
(AdjoinRoot P →ₐ[K] A) ≃
  GeneralGaugeJacobianOneFiberPoint ... A
```

where the right side is the literal fiber of the three actual
`MvPolynomial` coordinates at the chosen target. The equivalence commutes with
every algebra homomorphism `A →ₐ[K] B`.

Lean also records

```text
Module.finrank K (AdjoinRoot P) = P.natDegree
Algebra.Etale K (AdjoinRoot P)
Module.Finite K (AdjoinRoot P)
```

and specializes the literal-fiber theorem to the exact denominator-free
quintic map printed in the paper. At target `(1,0,-38)`, the resulting
equivalence is natural from maps out of
`ℚ[T]/((T^3-19)(T^2+T+1))`, whose rank is proved to be five.
The explicit quotient is proved finite étale. Lean also proves that this
literal fiber has no rational point and has both a real point and a
three-adic point.

The principal final declarations are:

```text
automaticRealizationMap_certificate
automaticRealization_pageOne
automaticRealizationGeometricDegree_eq
generalGaugeJacobianOneMap_targetDenormalization
automaticJacobianOneFiberRepresentingEquiv
automaticJacobianOneFiberRepresentingEquiv_natural
ExplicitQuintic.integralFiberRepresentingEquiv
ExplicitQuintic.integralFiberRepresentingEquiv_natural
automaticRepresentingAlgebra_etale
automaticRepresentingAlgebra_finite
ExplicitQuintic.p5_quotient_etale
ExplicitQuintic.p5_quotient_finite
ExplicitQuintic.p5_quotient_rank
ExplicitQuintic.p5_no_rational_root
ExplicitQuintic.integralFiberPoint_rat_isEmpty
ExplicitQuintic.integralFiberPoint_real_nonempty
ExplicitQuintic.integralFiberPoint_threeAdic_nonempty
generalGaugeGenericInversePolynomial_certificate
generalGaugeGenericInverseAdjoinRoot_finrank
generalGaugeFullyGenericInversePolynomial_certificate
generalGaugeFullyGenericInverseAdjoinRoot_finrank
generalGaugeMap_algebraicIndependent
generalGaugeFunctionFieldHom_injective
generalGaugeSourceFunctionFieldComparison
generalGaugeGeometricDegree_eq
realizationMapTarget_map
adjoinRootBaseChangeEquiv
```

No translation parameter, coefficient nonvanishing proof, chart unit, or
auxiliary abstract-fiber hypothesis remains as an external input.

## Formal chain

`GaugeAssembly.lean` and `GaugeInverseAssembly.lean` prove the low-degree and
complete finite high-degree coefficient identities over arbitrary commutative
rings. `GeneralGaugeMap.lean` packages them as one actual arbitrary-degree
three-variable polynomial map.

`GeneralGaugeJacobian.lean` expands the complete three-by-three Jacobian,
including all finite coefficient sums, and proves the determinant is exactly
`-2`. The fixed output scaling then has determinant `1`.

`GeneralGaugeInverse.lean` defines `G_π`, `β(π,S)`, and `E_{π,b,c}` and proves
the exact derivative factorization required by the represented-fiber theorem.
`GenericInverseIrreducibility.lean` constructs an exact polynomial-variable
swap, proves `H(S)-λC` irreducible over `K[C]`, applies Gauss's lemma over
`K(C)`, and verifies that, for each fixed `π ≠ 0` and `b`, the inverse equation
`E_{π,b,C}` is irreducible of exactly the seed degree. It also proves that
specializing the formal target variable recovers
`generalGaugeInversePolynomial` and that adjoining a root of this
fixed-parameter inverse equation gives an extension of exactly that degree.
`GeneralGaugeFullGenericDegree.lean` then takes the base field to be the
two-parameter rational function field `K(Π,B)`. It proves that the full
independent-parameter inverse equation over `K(Π,B)(C)` is irreducible of the
seed degree and that its root quotient has exactly that finrank.
`GeneralGaugeFunctionField.lean` proves coordinate-substitution injectivity
from the nonzero Jacobian, deduces algebraic independence of the displayed
coordinates, and constructs the induced injective pullback on rational
function fields. `GeneralGaugeFunctionFieldComparison.lean` identifies the
actual source function field over the canonical target presentation
`K(Π,B)(C)` with the inverse-root quotient and transfers its finrank to prove
that the displayed map has geometric degree exactly the seed degree.
`GeneralGaugeDisplayedFiber.lean` proves that evaluating the actual `B` and
`C` coordinates gives precisely the marked equations.

`GeneralGaugeRawFiber.lean` starts from a literal triple satisfying the three
polynomial equations. From `t*q = π` with unit `π`, it constructs the source
chart unit over an arbitrary commutative test ring and proves that the literal
fiber is naturally represented by the generic inverse quotient.
`GeneralGaugeNormalization.lean` transports this theorem to the determinant-one
output normalization at the zero second target coordinate.

`GeneralGaugeRealization.lean` specializes to
`G(S) = P(a+S)-P(a)`, proves that the chosen inverse polynomial is exactly
`P(a+S)`, translates the quotient back to `K[T]/(P)`, and then removes the
supplied parameter through `chosenAdmissibleTranslation`.
`GeneralGaugeRealizationDegree.lean` transports both the Jacobian and degree
certificates to the final automatically chosen map in terms of the original
input degree.
`GeneralGaugeBaseChange.lean` proves that polynomial translation, the full
general gauge, determinant-one normalization, the supplied-translation
realization map, and its distinguished target commute with a homomorphism of
ground fields. It also proves the tensor-product equivalence
`L ⊗[K] AdjoinRoot P ≃ₐ[L] AdjoinRoot (P.map f)`. The automatic classical
choice of a translation is deliberately excluded because it is not
functorial.

`FiniteEtaleQuotient.lean` packages every separable polynomial quotient as a
standard étale algebra and proves it finite. `GeneralGaugeFiberRank.lean`
records both finite étaleness and the exact dimension of the squarefree
representing quotient. `PageOneTheorem.lean` combines determinant one,
geometric degree, literal fiber representation, naturality, finite étaleness,
rank, and the effective degree bound in one theorem. `ExplicitFiber.lean`
specializes the reconstruction to
the exact denominator-free quintic map and target displayed in the paper,
including finite étaleness, naturality, and rank five.
`ExplicitPolynomial.lean` applies the rational-root theorem to certify that
the quintic has no rational root; `ExplicitFiber.lean` transfers this to
emptiness of the literal rational fiber and uses the intermediate value
theorem to construct a real point. `ExplicitThreeAdicPoint.lean` applies the
strong form of Hensel's lemma to `X^3-19` at `-2`, lifts the resulting root to
`ℚ_[3]`, and constructs a point on the literal displayed fiber.

## Foundations and axioms

The project contains no `sorry` and introduces no project-specific axioms.
The final functor-of-points theorem uses only Lean's standard foundations
`propext`, `Classical.choice`, and `Quot.sound`. The determinant and degree
certificates are algebraic theorems with no additional axioms.

The reconstruction works over arbitrary commutative test algebras. Units are
carried explicitly, so there is no hidden localization, omitted component, or
reduction to field-valued points. Separability makes the derivative class
invertible by Bézout; the unit first target coordinate makes the chart global
on the entire fiber functor.

## Remaining formal boundary

The actual map, determinant, geometric degree, effective degree, literal
fiber, finite étaleness, quotient translation, naturality, coordinate
algebraic independence, and explicit source-over-target function-field
comparison are formalized. The paper now proves exact reduced nonproperness
and boundary-sheet accounting separately. The remaining formalization tasks
for paper-level theorems are:

1. formalize the exact nonproperness locus, boundary-sheet accounting, and
   discriminant-order statements;
2. formalize the symmetric-monodromy and Hilbertian-specialization arguments;
3. formalize monogenicity of arbitrary finite étale algebras over an infinite
   field and compose it with the polynomial-presentation theorem;
4. formalize, or explicitly isolate as a classical theorem interface, the
   Campbell--Razar--Wright degree-two Galois case;
5. formalize the nonarchimedean local-number-theoretic and prime-distribution
   inputs in the Hasse-principle applications if those corollaries are to be
   machine-checked end to end. The rational obstruction, archimedean local
   point, and three-adic local point for the explicit quintic are already
   formalized.

The current certificate therefore proves the complete constructive,
scheme-theoretic, and geometric-degree polynomial-presentation layer, while
keeping the separately proved nonproperness theorem, monodromy, monogenicity,
the classical rank-two obstruction, and the remaining nonarchimedean
arithmetic inputs explicitly outside the Lean certificate.

## Build

```bash
cd formal/finite-etale-keller
lake build
```

Repository CI builds this project independently of the external Lean certificate
for the foundational three-dimensional map.
