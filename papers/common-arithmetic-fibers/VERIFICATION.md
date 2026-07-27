# Verification matrix

This file records the proof layer for every load-bearing claim and every
separately stated appendix theorem in *Prescribed Finite Étale Algebras as
Full Fibers of Keller Maps with Symmetric Monodromy*.  “Lean” means a theorem in the pinned
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
| Geometric and arithmetic generic monodromy are `S_N` for every admissible seed | Morse-slice lemma, the standalone ordered-root-cover specialization lemma, and full symmetric generic monodromy theorem | Not formalized | Exact `dC/dB=-r^2` audit in `verify_root_engineered_quadratic_gauge.py`; Serre's Morse-polynomial theorem |
| Infinitely many connected full Keller fibers with splitting-field group `S_N` over a Hilbertian field | Hilbertian specialization corollary, including the arithmetic/geometric equality and regularity argument | Not formalized | Hilbert irreducibility |
| Every finite étale algebra over an infinite field is monogenic | Monogenicity lemma, used only in the finite-étale realization corollary | Not yet formalized | Discriminant/Vandermonde proof |
| No characteristic-zero Keller map has generic degree two | Degree-two lemma; full descent proof in Appendix B | Not formalized | Campbell's unnumbered theorem on p. 244 (normal complex function-field extension), scalar-invariance of generic rank on a finite-locally-free open, quadratic separability/normality, and faithfully flat descent; Razar and Wright are cited as later algebraic treatments |
| Degree-four local–global theorem over every number field | Degree-four fixed-point lemma plus the stated Chebotarev density input | `DegreeFourFixedPoint.lean` proves the finite-group lemma as `degreeFour_fixedPoint`; the Chebotarev passage is not formalized | Chebotarev |
| Exact displayed quintic full Keller fiber at `(1,0,-38)` is represented by the finite étale algebra `ℚ[T]/((T³-19)(T²+T+1))`, naturally in test algebras, and has rank five | Minimal finite-étale Hasse theorem and optimal Keller Hasse failure corollary | `ExplicitPolynomial.lean`, `ExplicitMap.lean`, `ExplicitFiber.lean`; theorems `integralFiberRepresentingEquiv_natural`, `p5_quotient_etale`, `p5_quotient_finite`, and `p5_quotient_rank` | Exact quintic checker |
| Explicit quintic has no rational point and has points over `ℝ` and every `ℚ_p` | Minimal finite-étale Hasse theorem: exhaustive prime table and direct Hensel data at `2,3,19` | `integralFiberPoint_hasse_certificate`, assembled from the rational obstruction, real point, direct Hensel witnesses at `2,3,19`, and the two generic residue-class theorems | Rational-root theorem, intermediate value theorem, cyclicity and coprime-power bijectivity in `(ZMod p)ˣ`, and Hensel's lemma |
| Optimal rank-five Hasse failure | First proved for finite étale schemes by the fixed-point/Chebotarev theorem, then transferred by the realization theorem | The exact polynomial, map, literal quotient fiber, naturality, determinant, rank, rational obstruction, and local points at every completion are formalized. `DegreeFourMomentBarrier.lean` follows the alternative zeta-moment route: it proves the finite-étale component adapter, rank-at-most-four local-sheet inequality, tensor-square identity, strict tensor surplus, absolute convergence and positivity of the normalized prime sums, and the concrete first/second-moment contradiction. The combined endpoint is `no_rank_le_four_hasse_failure_of_rationalPrimeMomentStatement`. Only the Dedekind-zeta first-prime-moment extraction is not yet formalized | Exact exceptional-prime assertions and residue-cover regression in `verify_minimal_hasse_keller_fiber.py`; paper proof uses Chebotarev |
| Compatibility of the realization with extension of the ground field | Base-change proposition in Section 4 | `GeneralGaugeBaseChange.lean` proves coefficientwise compatibility of translation, the full supplied-parameter map, normalization, admissibility, squarefreeness, and distinguished target, together with `L ⊗[K] AdjoinRoot P ≃ₐ[L] AdjoinRoot (P.map f)` | Coefficientwise paper proof |
| Exact reduced nonproperness locus `S_F = V_red(Disc_S E)` over every algebraically closed characteristic-zero field | Appendix D, exact reduced quadratic-gauge nonproperness theorem: graph-boundary definition, fiber-cardinality criterion, localized fibers, and complete `Π=0` table | Not formalized; explicitly listed in the remaining formal boundary | Jelonek Proposition 6 is stated over `ℂ`; the paper supplies descent to a finitely generated field, embedding into `ℂ`, and graph-boundary/fiber base-change invariance |
| Complete `Π=0` fiber table and exact global discriminant factor `Π^(N²-3N-2)` for `N ≥ 4` | Appendix D: direct `q=0` and `t=0` source charts; standalone exact-order lemma with all three residual polynomials, a root-pair valuation ledger, and the first nonzero coefficient in `K[B,C][Π]` | Not formalized | `verify_quadratic_gauge_nonproperness.py` checks the charts, exact orders and saturated-slice coefficients through degree ten, and the Newton ledger through degree 64 |

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
  quotient-ring reconstruction in both directions, the exact quintic scaling,
  determinant, target, and inverse polynomial;
- `verify_quadratic_gauge_nonproperness.py`: the cubic and quartic
  discriminant specializations, direct `q=0` and `t=0` charts, exact
  discriminant orders and saturated-slice coefficients through degree ten,
  and the Newton ledger through degree 64; and
- `verify_universal_quadratic_gauge.sing`: a fresh expansion of the generic
  degree-six Jacobian over a rational function field in six independent
  coefficients, together with the `(7,38,36)` degree profile.

These are exact symbolic computations from implementations independent of the
Lean development. They audit the displayed algebra but do not replace the
ordinary proofs of nonproperness, monodromy, Hilbert
specialization, degree two, or the remaining local arithmetic.

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

For the paper's explicit quintic map, Lean additionally proves that the exact
denominator-free coordinates at target `(1,0,-38)` are naturally represented
by the displayed Berend--Bilu quotient.  The declarations are
`integralFiberRepresentingEquiv`, `integralFiberRepresentingEquiv_natural`,
`p5_quotient_etale`, `p5_quotient_finite`, and `p5_quotient_rank`. Lean also
proves the global obstruction and the archimedean local point in
`integralFiberPoint_rat_isEmpty` and `integralFiberPoint_real_nonempty`.
The Hensel certificate `integralFiberPoint_threeAdic_nonempty` additionally
constructs a point on the literal fiber over `ℚ_[3]`.

## Remaining formal boundary

The polynomial realization theorem is now formalized end to end, including
actual geometric degree.  The
remaining steps needed for a single theorem starting from an arbitrary finite
étale algebra are:

1. formalize monogenicity of finite étale products over an infinite field;
2. compose monogenicity with the polynomial-presentation certificate;
3. either formalize the Campbell--Razar--Wright Galois case or keep it as a
   clearly isolated classical theorem interface;
4. formalize symmetric monodromy and Hilbertian specialization;
5. formalize the Chebotarev passage used in the paper (the finite-`G`-set
   endpoint `degreeFour_fixedPoint` is now formalized), or
   complete the alternative Lean route by formalizing the Dedekind-zeta
   first-prime-moment extraction.  The finite-étale
   component adapter, rank-at-most-four local-sheet inequality, tensor-square
   identity, strict tensor surplus, positive-moment contradiction, explicit
   quintic's rational obstruction, and its points over the real and every
   `p`-adic completion are already formalized.  Mathlib's nonzero simple-pole
   theorem is already connected; the missing step is the Euler-product
   coefficient extraction.

The repository keeps these boundaries explicit so that a compiled certificate
is never described as proving a stronger statement than it does.
