# Lean formalization: finite étale Keller fibers

This project formalizes the polynomial-presentation construction in
*Every Finite Étale Algebra of Rank at Least Three Is a Full Keller Fiber*.
It uses Lean `v4.33.0-rc1` and Mathlib at the matching release
candidate.

The paper's supplied-presentation theorem is stated over fields of
characteristic different from two, using nonvanishing of the first and third
Hasse coefficients at the supplied translation. The current end-to-end Lean
development retains the stronger assumption `CharZero K`; it formalizes the
complete characteristic-zero specialization and the automatic abstract
finite-étale corollary.

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
| 20 | Explicit quintic: no rational point, a real point, and points over every `p`-adic field | implemented |
| 21 | Finite étaleness of every squarefree quotient and of the explicit quintic quotient | implemented |
| 22 | Coordinate algebraic independence and injective function-field pullback | implemented |
| 23 | Supplied-parameter map and target compatibility under scalar extension | implemented |
| 24 | Tensor-product base change of the representing polynomial quotient | implemented |
| 25 | Comparison with the actual pullback-field extension and geometric degree | implemented |
| 26 | Full algebraic degree-four barrier: finite-étale decomposition, local-sheet bound, tensor surplus, and positive-moment contradiction | implemented |
| 27 | Finite-group fixed-point lemma for actions on at most four points | implemented |
| 28 | Dedekind-zeta first prime moment (Euler-coefficient extraction) | ordinary proof; not yet Lean |
| 29 | Monogenicity and the passage from arbitrary finite étale algebras in characteristic zero | implemented |
| 30 | Historical degree-two Galois exclusion | external theorem; not yet Lean |
| 31 | Common power-shift and cubic-lift Jacobians, selected-fiber invariance, and quotient representation | implemented |

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

It also exposes the underlying standard Mathlib statement directly:

```text
letI := automaticRealizationTargetFunctionFieldAlgebra P hdeg
Module.finrank
  (RatFunc (FractionRing (MvPolynomial (Fin 2) K)))
  (FractionRing (MvPolynomial (Fin 3) K)) = P.natDegree
```

The theorem type itself therefore uses only standard Mathlib fraction-field
objects. Its algebra structure is induced by the three target-coordinate
pullbacks. The explicit comparison
`automaticRealizationFunctionFieldComparison` is an `AlgEquiv` between this
source fraction field and the corresponding `AdjoinRoot` extension over
`K(Π,B)(C)`. Thus the custom geometric-degree abbreviation is not the
public endpoint on which the degree claim rests.

`automaticRealization_pageOne` bundles the direct finrank equality as
`functionFieldFinrank`, together with the convenience equality above,
determinant one, and all of the fiber assertions below.

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

and separately specializes the literal-fiber theorem to an exact
denominator-free quintic certificate. At target `(1,0,-38)`, the resulting
equivalence is natural from maps out of
`ℚ[T]/((T^3-19)(T^2+T+1))`, whose rank is proved to be five.
The explicit quotient is proved finite étale. Lean also proves that this
literal fiber has no rational point, has a real point, and has a point over
`ℚ_[p]` for every prime `p`.

`StableGaugeFiber.lean` formalizes the deformation layer used by universal
multiplicity.  For every common power shift, and for every cubic lift exponent
`n ≥ 4`, it proves the normalized Jacobian determinant is one.  On `Π=1` the
deformed coordinates evaluate identically to the undeformed coordinates, so
the original quotient algebra represents the literal stable-map fiber over
every commutative test algebra, naturally under algebra homomorphisms.  The
two generated stable arithmetic modules
instantiate these statements for the ramified quintic at `m=2` and
`T^3-T-1` at `n=7`.  The separate boundary invariants proving pairwise stable
inequivalence are not formalized here.

The principal final declarations are:

```text
automaticRealizationMap_certificate
automaticRealization_pageOne
finiteEtalePowerBasis
finiteEtalePolynomial_squarefree
finiteEtalePresentation
finiteEtalePolynomial_natDegree
abstractFiniteEtaleFiberRepresentingEquiv
abstractFiniteEtaleFiberRepresentingEquiv_natural
abstractFiniteEtale_pageOne
automaticRealizationGeometricDegree_eq
automaticRealizationFunctionFieldComparison
automaticRealizationFunctionField_finrank
generalGaugeJacobianOneMap_targetDenormalization
automaticJacobianOneFiberRepresentingEquiv
automaticJacobianOneFiberRepresentingEquiv_natural
ExplicitQuintic.integralFiberRepresentingEquiv
ExplicitQuintic.integralFiberRepresentingEquiv_natural
ExplicitQuintic.integralFiberPoint_padic_nonempty
ExplicitQuintic.integralFiberPoint_hasse_certificate
automaticRepresentingAlgebra_etale
automaticRepresentingAlgebra_finite
powerShiftedGaugeRealizationFiberRepresentingEquiv
cubicLiftGaugeRealizationFiberRepresentingEquiv
powerShiftedGaugeRealizationFiberRepresentingEquiv_natural
cubicLiftGaugeRealizationFiberRepresentingEquiv_natural
jacobianDet_powerShiftedGaugeJacobianOneMap
jacobianDet_cubicLiftGaugeJacobianOneMap
ExplicitQuintic.p5_quotient_etale
ExplicitQuintic.p5_quotient_finite
ExplicitQuintic.p5_quotient_rank
ExplicitQuintic.p5_no_rational_root
ExplicitQuintic.integralFiberPoint_rat_isEmpty
ExplicitQuintic.integralFiberPoint_real_nonempty
ExplicitQuintic.integralFiberPoint_threeAdic_nonempty
localPointCount_tensor_self
PositiveNormalizedMean.second_moment_eq_sq_of_bounds
PositiveNormalizedMean.contradiction_of_component_surplus
second_moment_eq_sq_of_dirichletPrimeMean
no_rank_le_four_hasse_failure_of_rationalPrimeMomentStatement
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

## Abstract finite étale theorem

Let `A` be any finite étale algebra over a characteristic-zero field `K`, and
assume `3 ≤ Module.finrank K A`. Lean noncomputably selects a power basis and
its squarefree minimal polynomial

```text
finiteEtalePolynomial K A : K[X]
```

and constructs a presentation

```text
finiteEtalePresentation K A :
  AdjoinRoot (finiteEtalePolynomial K A) ≃ₐ[K] A.
```

The selected polynomial has degree `Module.finrank K A`. Composing this
presentation with the polynomial-presentation theorem gives the public
certificate

```text
abstractFiniteEtale_pageOne
```

It supplies a determinant-one map, its target, geometric degree equal to the
rank of `A`, the `6 * rank + 2` coordinate-degree bound, and, naturally in
every commutative test `K`-algebra `R`, an equivalence

```text
(A →ₐ[K] R) ≃ GeneralGaugeJacobianOneFiberPoint ... R.
```

The primitive element is noncanonical, so the construction is intentionally
noncomputable.

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
rank, and the effective degree bound in one theorem.
`AbstractFiniteEtale.lean` proves monogenicity for every finite étale algebra
over a characteristic-zero field: after decomposing it into a finite product
of finite separable field extensions, it translates primitive generators to
give them pairwise distinct traces and applies the Chinese remainder theorem.
It then composes the resulting squarefree polynomial presentation with
`PageOneTheorem.lean`, including the natural equivalence on all test
algebras. `ExplicitFiber.lean` separately specializes the reconstruction to
an exact denominator-free quintic map and target,
including finite étaleness, naturality, and rank five.
`ExplicitPolynomial.lean` applies the rational-root theorem to certify that
the quintic has no rational root; `ExplicitFiber.lean` transfers this to
emptiness of the literal rational fiber and uses the intermediate value
theorem to construct a real point. `ExplicitThreeAdicPoint.lean` applies the
strong form of Hensel's lemma to `X^3-19` at `-2`, lifts the resulting root to
`ℚ_[3]`, and constructs a point on the literal displayed fiber.
`ExplicitAllPadicPoints.lean` proves a generic simple-root Hensel interface,
formalizes the direct witnesses at `2` and `19`, handles the two remaining
prime classes through the unit group of `ZMod p`, and packages the rational,
real, and all-`p`-adic assertions into one Hasse certificate.

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

`DegreeFourMomentBarrier.lean` formalizes the rank-minimality proof after its
single analytic input has been supplied.  It proves that algebra maps from a
tensor product are pairs, hence
`localPointCount K (A ⊗[K] A) L = localPointCount K A L ^ 2`.  It then proves
Mathlib's finite-étale field-product decomposition and identifies component
count with the number of factors.  For rank at most four it proves directly
that every local point supplies at least that many local sheets (including
the two-quadratic case), and proves that every nontrivial diagonal tensor
block contributes an extra connected component.  Lean now also proves
absolute convergence, linearity, and order preservation for the actual
normalized Dirichlet prime sums on bounded functions.  The public theorem
`no_rank_le_four_hasse_failure_of_rationalPrimeMomentStatement` therefore
derives the contradiction directly from
`RationalFiniteEtalePrimeMomentStatement`; no abstract mean functional
remains on the critical path.  This route uses neither a Galois action nor a
monogenic presentation.

`DegreeFourFixedPoint.lean` independently formalizes the exact finite-group
lemma used by the shorter arithmetic proof.  Its endpoint
`degreeFour_fixedPoint` says that an action on at most four points has a
global fixed point whenever every group element has a fixed point.  The proof
uses orbit decomposition and Mathlib's Burnside lemma.  The arithmetic
Chebotarev passage from local points to the elementwise fixed-point hypothesis
is not part of this declaration.

## Remaining formal boundary

The actual map, determinant, geometric degree, effective degree, literal
fiber, finite étaleness, quotient translation, naturality, coordinate
algebraic independence, explicit source-over-target function-field
comparison, base change for supplied data, monogenicity, and the abstract
finite-étale corollary are formalized in characteristic zero. The remaining
formal boundary of the focused paper is:

1. generalize the supplied-seed, reconstruction, function-field, and
   base-change layers from `CharZero K` to the exact characteristic-not-two
   hypotheses used by the paper;
2. formalize, or explicitly isolate as a classical theorem interface, the
   Campbell--Razar--Wright degree-two Galois case;
3. for the separate arithmetic development, formalize either the Chebotarev
   passage or the first-prime-moment theorem from the Dedekind-zeta Euler
   product if rank-minimality is to be machine-checked end to end.  The finite-group
   fixed-point lemma used after Chebotarev is now formalized as
   `degreeFour_fixedPoint`.  Mathlib's
   nonzero simple-pole theorem is exposed as
   `dedekindZeta_simplePole_input`; the remaining step is the Euler-product
   coefficient extraction.  Absolute convergence, linearity and positivity
   of the normalized prime mean, the finite-étale component adapter, low-rank
   local-sheet bound, tensor surplus, moment contradiction, explicit
   quintic's rational obstruction, archimedean point, and points over every
   nonarchimedean completion are formalized.

The current certificates therefore prove the complete constructive,
scheme-theoretic, and geometric-degree layers for both polynomial
presentations and abstract finite étale algebras in characteristic zero.
The classical rank-two obstruction and the analytic first-prime-moment
extraction remain explicitly outside the Lean certificate. Nonproperness,
monodromy, and stable atomicity are companion results outside the focused
paper.

## Build

```bash
cd formal/finite-etale-keller
lake build
```

Repository CI builds this project independently of the external Lean certificate
for the foundational three-dimensional map.
