# Verification matrix

This file records the proof layer for every load-bearing claim in
*Every Finite Étale Algebra of Rank at Least Three Is a Full Keller Fiber*.
“Lean” means a theorem in the pinned
`formal/finite-etale-keller` project with no `sorry` and no project-specific
axiom.  Symbolic checks are exact independent audits, not substitutes for the
corresponding mathematical argument.

| Claim | Paper proof | Lean certificate | Independent audit / external input |
|---|---|---|---|
| Existence and automatic choice of an admissible translation | Polynomial realization theorem and its effective corollary | `Admissibility.lean`, `AutomaticRealization.lean` | Concrete translated seeds |
| Two-sided source/chart equivalence over arbitrary commutative rings | General localized quadratic-gauge fiber theorem | `SourceEquivalence.lean` | Universal source-chart identity |
| For a separable inverse polynomial, its derivative is a unit at every test-algebra root | Derivative-units lemma | `Bezout.lean`, `SeparableReconstruction.lean` | Quotient-ring regressions |
| Roots represent the derivative-localized fiber, naturally in the test algebra and without a separability hypothesis | General localized quadratic-gauge fiber theorem | `LocalizedFiberPoints.lean`, `LocalizedGaugeFiberPoints.lean`, `GeneralGaugeLocalizedFiber.lean` | Exact two-sided reconstruction |
| Squarefree inverse equations give finite étale quotient algebras representing the literal fibers | Finite-étale (squarefree) fiber corollary | `FiniteEtaleQuotient.lean`, `GeneralGaugeFiberRank.lean`, `GeneralGaugeRawFiber.lean` | Quotient-ring reconstruction |
| Uniform finite-sum assembly of the displayed coordinates | Compact design in Section 2; coefficientwise pullback in Appendix A | `GaugeAssembly.lean`, `GaugeInverseAssembly.lean` | Structural SymPy certificate |
| One arbitrary-degree `MvPolynomial (Fin 3) K` map and exact coordinate evaluation | Polynomial-pullback subsection and displayed gauge map | `GeneralGaugeMap.lean` | Structural checker and six-coefficient bridge |
| Generic inverse polynomial, explicit `β`, and derivative factorization | Compact inverse-design equations | `GeneralGaugeInverse.lean` | Concrete inverse-polynomial regressions |
| Actual displayed coordinate equations equal the represented source equations | Compact design in Section 2 and coefficientwise verification in Appendix A | `GeneralGaugeDisplayedFiber.lean` | SymPy and Singular expansions |
| Literal three-coordinate map fiber is represented, naturally over every test algebra | General localized quadratic-gauge fiber theorem | `GeneralGaugeRawFiber.lean`, `GeneralGaugeLocalizedFiber.lean` | Quotient-ring reconstruction |
| Zero-second-coordinate fiber is preserved by determinant-one output normalization | Jacobian-one effective-form corollary | `GeneralGaugeNormalization.lean` | Exact scalar-normalization checks |
| Translation `K[S]/(P(a+S)) ≃ K[T]/(P)` | Polynomial realization theorem | `TranslationQuotient.lean`, `GeneralGaugeRealization.lean` | Concrete translated examples |
| Final automatic page-one certificate from squarefree `P` | Polynomial realization theorem and its effective corollary | `PageOneTheorem.lean`; final proposition `automaticRealization_pageOne` simultaneously contains target normalization, determinant, geometric degree, literal fiber representation, naturality, finite étaleness, finiteness, rank, and degree bound | Full Lean build and axiom reports |
| Represented special-fiber length is `deg P` | Polynomial realization theorem and finite-étale fiber corollary | `GeneralGaugeFiberRank.lean`; theorem `automaticRealizationFiber_rank` | Standard polynomial-quotient dimension theorem |
| Constant determinant `-2` and normalized determinant `1` for the arbitrary-degree actual map | Compact Jacobian factorization, general localized theorem, and effective-form corollary | `GeneralGaugeJacobian.lean` | Structural SymPy audit; independent generic degree-six Singular audit; concrete degrees 3–5 |
| Final coordinate-degree bound `6N+2` in terms of `N = deg P` | Polynomial realization theorem and effective normalization corollary | `GeneralGaugeDegree.lean`, `GeneralGaugeRealizationDegree.lean` | Structural termwise audit; Singular degree-six profile |
| Irreducibility and degree `N` of the fully independent inverse equation over the iterated target field `K(Π,B)(C)` | General localized quadratic-gauge fiber theorem, primitive linear-parameter argument | `GenericInverseIrreducibility.lean` proves the fixed-`π,b` engine; `GeneralGaugeFullGenericDegree.lean` promotes `Π,B` to independent parameters and proves `generalGaugeFullyGenericInversePolynomial_certificate` and `generalGaugeFullyGenericInverseAdjoinRoot_finrank` | Exact polynomial-variable swap and Mathlib Gauss lemma; concrete degree 3–5 regressions |
| Function-field reconstruction and geometric degree `N` | General localized quadratic-gauge fiber theorem: explicit equality `K(x,y,z)=K(Π,B,C)(S)` followed by inverse-polynomial irreducibility | Complete for the polynomial-presentation map: `GeneralGaugeFunctionField.lean` proves algebraic independence and injective pullback; `GeneralGaugeFunctionFieldComparison.lean` constructs `generalGaugeSourceFunctionFieldComparison : K(x,y,z) ≃ K(Π,B)(C)[S]/(E)` over the actual target embedding and proves `generalGaugeGeometricDegree_eq`; `PageOneTheorem.lean` transports the result to the determinant-one realization | Independent rational reconstruction in the paper |
| Every finite étale algebra over an infinite field is monogenic | Monogenicity lemma, used only in the finite-étale realization corollary | `AbstractFiniteEtale.lean` proves the characteristic-zero case needed here, constructs `finiteEtalePresentation`, and composes it with the polynomial certificate in `abstractFiniteEtale_pageOne` | Paper's discriminant/Vandermonde proof; Lean uses translated primitive elements with distinct traces and the Chinese remainder theorem |
| No characteristic-zero Keller map has generic degree two | Rank-classification corollary, explicitly separated from the constructive theorem | Not formalized | Campbell's unnumbered theorem on p. 244 (normal complex function-field extension), scalar-invariance of generic rank under extension, quadratic separability/normality, and faithfully flat descent; Razar and Wright are cited as later algebraic treatments |
| Compatibility of the realization with extension of the ground field | Base-change proposition in Section 4 | `GeneralGaugeBaseChange.lean` proves coefficientwise compatibility of translation, the full supplied-parameter map, normalization, admissibility, squarefreeness, and distinguished target, together with `L ⊗[K] AdjoinRoot P ≃ₐ[L] AdjoinRoot (P.map f)` | Coefficientwise paper proof |
| Explicit Berend--Bilu quintic Hasse fiber | Explicit arithmetic example: displayed map, target, quotient, and local/global root audit | `ExplicitFiber.lean` and `ExplicitAllPadicPoints.lean`; endpoint `integralFiberPoint_hasse_certificate` | `verify_finite_etale_keller_fibers.py` |

## Independent exact audits rerun

The following independent commands were rerun successfully on 27 July 2026:

- `verify_universal_quadratic_gauge.py`: source-chart reciprocal identities,
  marked-line Jacobian cancellation, generic coefficient assembly, the
  degree-six bridge, the `6N+2` bound, and determinant-one normalization;
- `verify_root_engineered_quadratic_gauge.py`: coefficient engineering through
  degree six, the cubic and quartic regressions, the discriminant differential
  `dC/dB = -r²`, quotient reconstruction identities, and translated seeds
  through degree twelve;
- `verify_finite_etale_keller_fibers.py`: degrees three through five,
  with quotient-ring reconstruction in both directions;
- `verify_universal_quadratic_gauge.sing`: a fresh expansion of the generic
  degree-six Jacobian over a rational function field in six independent
  coefficients, together with the `(7,38,36)` degree profile.

These are exact symbolic computations from implementations independent of the
Lean development. They audit the displayed algebra but do not replace the
ordinary proof of the classical degree-two interface.

## Formal theorem now obtained

For every characteristic-zero field `K`, every squarefree polynomial
`P : K[X]`, and every proof that `3 ≤ P.natDegree`, Lean now:

1. chooses an admissible translation parameter internally;
2. constructs the actual arbitrary-degree map
   `automaticRealizationMap P hdeg : Fin 3 → MvPolynomial (Fin 3) K`;
3. proves its Jacobian determinant is `1`;
4. proves the inverse target scaling from the determinant-one map back to the
   raw gauge;
5. constructs the source-field equivalence with the independent inverse-root
   quotient and proves geometric degree `P.natDegree`;
6. proves every coordinate has total degree at most `6 * P.natDegree + 2`;
7. constructs the literal target fiber of those three polynomial coordinates;
8. gives a natural equivalence from `AdjoinRoot P` to that literal fiber over
   every commutative test `K`-algebra;
9. proves that the represented quotient algebra has dimension
   `P.natDegree`; and
10. proves that the quotient algebra is finite étale over `K`.

The single combined proposition is `automaticRealization_pageOne`.  Its
geometric-degree field is supplied by
`automaticRealizationGeometricDegree_eq`; the explicit bridge is
`generalGaugeSourceFunctionFieldComparison`, and its degree theorem is
`generalGaugeGeometricDegree_eq`.  The construction, fiber, and finiteness
layers remain separately available as
`automaticRealizationMap_certificate`,
`automaticJacobianOneFiberRepresentingEquiv_natural`,
`automaticRepresentingAlgebra_etale`, and
`automaticRepresentingAlgebra_finite`.

For an admissible seed `G`, Lean also promotes `Π` and `B` to independent
parameters, proves the resulting inverse equation over `K(Π,B)(C)`
irreducible with degree and root-quotient finrank `G.natDegree`, proves the
three displayed coordinates algebraically independent, constructs the
injective pullback on rational function fields, and proves that the resulting
extension is the actual source-over-target extension. The principal
declarations are
`generalGaugeFullyGenericInversePolynomial_certificate`,
`generalGaugeFullyGenericInverseAdjoinRoot_finrank`,
`generalGaugeMap_algebraicIndependent`, and
`generalGaugeFunctionFieldHom_injective`, followed by
`generalGaugeSourceFunctionFieldComparison` and
`generalGaugeGeometricDegree_eq`. For a supplied translation
parameter, `realizationMapTarget_map` additionally certifies coefficientwise
compatibility of the complete map-target pair under scalar extension, while
`adjoinRootBaseChangeEquiv` certifies the tensor-product base change of the
representing quotient.

## Remaining formal boundary

The polynomial realization theorem is formalized end to end, including
actual geometric degree. Monogenicity and its composition with that theorem
are also formalized in characteristic zero, yielding a single certificate
starting from an abstract finite étale algebra. The sole remaining formal
boundary of the focused paper is:

1. either formalize the Campbell--Razar--Wright Galois case or keep it as a
   clearly isolated classical theorem interface;

Symmetric monodromy, stable atomicity, Hilbertian specialization, and exact
nonproperness are maintained as companion results outside this paper and
outside this verification matrix.
