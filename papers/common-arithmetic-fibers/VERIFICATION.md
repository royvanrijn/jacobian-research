# Verification matrix

This file records the proof layer for every load-bearing claim in *Prescribed
Finite Étale Algebras as Full Fibers of Keller Maps with Symmetric
Monodromy*.  “Lean” means a theorem in the pinned
`formal/finite-etale-keller` project with no `sorry` and no project-specific
axiom.  Symbolic checks are exact independent audits, not substitutes for the
corresponding mathematical argument.

| Claim | Paper proof | Lean certificate | Independent audit / external input |
|---|---|---|---|
| Existence and automatic choice of an admissible translation | Prescribed-algebra realization theorem (Theorem 2.3), final paragraph | `Admissibility.lean`, `AutomaticRealization.lean` | Concrete translated seeds |
| Two-sided source/chart equivalence over arbitrary commutative rings | General localized quadratic-gauge fiber theorem (Theorem 3.2) | `SourceEquivalence.lean` | Universal source-chart identity |
| For a separable inverse polynomial, its derivative is a unit at every test-algebra root | Lemma 3.1 | `Bezout.lean`, `SeparableReconstruction.lean` | Quotient-ring regressions |
| Roots represent the derivative-localized fiber, naturally in the test algebra and without a separability hypothesis | General localized quadratic-gauge fiber theorem (Theorem 3.2) | `LocalizedFiberPoints.lean`, `LocalizedGaugeFiberPoints.lean`, `GeneralGaugeLocalizedFiber.lean` | Exact two-sided reconstruction |
| Squarefree inverse equations give finite étale quotient algebras representing the literal fibers | Finite-étale (squarefree) fiber corollary (Corollary 3.3) | `FiniteEtaleQuotient.lean`, `GeneralGaugeFiberRank.lean`, `GeneralGaugeRawFiber.lean` | Quotient-ring reconstruction |
| Uniform finite-sum assembly of the displayed coordinates | Subsection 3.1 | `GaugeAssembly.lean`, `GaugeInverseAssembly.lean` | Structural SymPy certificate |
| One arbitrary-degree `MvPolynomial (Fin 3) K` map and exact coordinate evaluation | Equations (3.1)–(3.3) | `GeneralGaugeMap.lean` | Structural checker and six-coefficient bridge |
| Generic inverse polynomial, explicit `β`, and derivative factorization | Equations (3.4)–(3.7) | `GeneralGaugeInverse.lean` | Concrete inverse-polynomial regressions |
| Actual displayed coordinate equations equal the represented source equations | General localized quadratic-gauge fiber theorem (Theorem 3.2) | `GeneralGaugeDisplayedFiber.lean` | SymPy and Singular expansions |
| Literal three-coordinate map fiber is represented, naturally over every test algebra | General localized quadratic-gauge fiber theorem (Theorem 3.2) | `GeneralGaugeRawFiber.lean`, `GeneralGaugeLocalizedFiber.lean` | Quotient-ring reconstruction |
| Zero-second-coordinate fiber is preserved by determinant-one output normalization | Effective normalization corollary (Corollary 3.7) | `GeneralGaugeNormalization.lean` | Exact scalar-normalization checks |
| Translation `K[S]/(P(a+S)) ≃ K[T]/(P)` | Prescribed-algebra realization theorem (Theorem 2.3) | `TranslationQuotient.lean`, `GeneralGaugeRealization.lean` | Concrete translated examples |
| Final automatic literal-fiber realization from squarefree `P` | Prescribed-algebra realization theorem (Theorem 2.3) | `GeneralGaugeRealization.lean`; final theorem `automaticJacobianOneFiberRepresentingEquiv_natural` | Full Lean build and axiom reports |
| Represented special-fiber length is `deg P` | Prescribed-algebra realization theorem (Theorem 2.3) and finite-étale fiber corollary (Corollary 3.3) | `GeneralGaugeFiberRank.lean`; theorem `automaticRealizationFiber_rank` | Standard polynomial-quotient dimension theorem |
| Constant determinant `-2` and normalized determinant `1` for the arbitrary-degree actual map | General localized theorem (Theorem 3.2) and effective normalization corollary (Corollary 3.7) | `GeneralGaugeJacobian.lean` | Structural SymPy audit; independent generic degree-six Singular audit; concrete degrees 3–5 |
| Final coordinate-degree bound `6N+2` in terms of `N = deg P` | Effective normalization corollary (Corollary 3.7) and prescribed-algebra realization theorem (Theorem 2.3) | `GeneralGaugeDegree.lean`, `GeneralGaugeRealizationDegree.lean` | Structural termwise audit; Singular degree-six profile |
| Irreducibility and degree `N` of the fully independent inverse equation over the iterated target field `K(Π,B)(C)` | General localized quadratic-gauge fiber theorem (Theorem 3.2), primitive linear-parameter argument | `GenericInverseIrreducibility.lean` proves the fixed-`π,b` engine; `GeneralGaugeFullGenericDegree.lean` promotes `Π,B` to independent parameters and proves `generalGaugeFullyGenericInversePolynomial_certificate` and `generalGaugeFullyGenericInverseAdjoinRoot_finrank` | Exact polynomial-variable swap and Mathlib Gauss lemma; concrete degree 3–5 regressions |
| Function-field reconstruction and geometric degree `N` | General localized quadratic-gauge fiber theorem (Theorem 3.2), generic reconstruction | Partial: `GeneralGaugeFunctionField.lean` proves coordinate substitution is injective, the displayed coordinates are algebraically independent, and the induced function-field pullback is injective; `GeneralGaugeFullGenericDegree.lean` proves the independent inverse-root extension has finrank `N`; the explicit comparison with the actual pullback-field extension and generic reconstruction remain | Paper's reconstruction argument |
| Exact nonproperness set `S_(F_G)=V(Disc_S(E))` | Exact quadratic-gauge nonproperness theorem | Not yet formalized | Independent SymPy/Singular discriminants; direct `t=0`, `q=0` fibers; normalization boundary ledger |
| Exact sheet loss over every inverse-root partition | Exact nonproperness theorem | Partial: the localized fiber equivalence for `Pi != 0` excludes derivative-zero roots, but the multiplicity-by-multiplicity missing-sheet count is not formalized | Repeated-root reconstruction and quartic image regressions |
| `Pi=0` fibers and `N-3` generic missing sheets | Exact nonproperness theorem and boundary-arc proposition | Not yet formalized | Newton polygon through degree 64; direct source-divisor calculation |
| Discriminant order `N^2-3N-2` and saturated slice `B^2(1-BC)` | Exact nonproperness theorem | Not yet formalized | Dedicated nonproperness checker through degree 10; independent root-valuation derivation |
| Distinguished target lies in the maximal finite-étale locus | Intrinsic-fullness corollary | Partial: inverse specialization to `P(a+S)` and finite étaleness of its quotient are formalized; identification of the maximal locus via nonproperness is not | Nonzero discriminant of squarefree `P` |
| Geometric and arithmetic generic monodromy are `S_N` for every admissible seed | Full symmetric generic monodromy theorem (Theorem 3.4) | Not formalized | Exact `dC/dB=-r^2` audit in `verify_root_engineered_quadratic_gauge.py`; Serre's Morse-polynomial theorem |
| Infinitely many connected full `S_N`-fibers over a Hilbertian field | Hilbertian specialization corollary (Corollary 3.5) | Not formalized | Hilbert irreducibility |
| Every finite étale algebra over an infinite field is monogenic | Lemma 4.1 | Not yet formalized | Discriminant/Vandermonde proof |
| No characteristic-zero Keller map has generic degree two | Lemma 2.2 | Not formalized | Campbell–Razar–Wright plus faithfully flat descent |
| Exact displayed quintic fiber at `(1,0,-38)` is represented by the finite étale algebra `ℚ[T]/((T³-19)(T²+T+1))`, naturally in test algebras, and has rank five | Theorem 5.2, displayed map and fiber | `ExplicitPolynomial.lean`, `ExplicitMap.lean`, `ExplicitFiber.lean`; theorems `integralFiberRepresentingEquiv_natural`, `p5_quotient_etale`, `p5_quotient_finite`, and `p5_quotient_rank` | Exact quintic checker |
| Explicit quintic has no rational point and has real and three-adic points | Theorem 5.2 | `p5_no_rational_root`, `integralFiberPoint_rat_isEmpty`, `integralFiberPoint_real_nonempty`, `integralFiberPoint_threeAdic_nonempty` | Rational-root theorem, intermediate value theorem, and Hensel's lemma at `-2` |
| Optimal rank-five Hasse failure | Theorem 5.2 | The exact polynomial, map, literal quotient fiber, naturality, determinant, rank, rational obstruction, archimedean point, and `ℚ_3` point are formalized; the remaining nonarchimedean local points and optimality are not | Local-solubility proof and degree-four barrier |
| One fixed map has infinitely many Hasse-failing fibers | Theorem 6.1 | Not yet formalized | Exact inverse-family checker plus local and prime-counting proof |
| Compatibility of the realization with extension of the ground field | Base-change proposition in Section 4 | `GeneralGaugeBaseChange.lean` proves coefficientwise compatibility of translation, the full supplied-parameter map, normalization, admissibility, squarefreeness, and distinguished target, together with `L ⊗[K] AdjoinRoot P ≃ₐ[L] AdjoinRoot (P.map f)` | Coefficientwise paper proof |

## Formal theorem now obtained

For every characteristic-zero field `K`, every squarefree polynomial
`P : K[X]`, and every proof that `3 ≤ P.natDegree`, Lean now:

1. chooses an admissible translation parameter internally;
2. constructs the actual arbitrary-degree map
   `automaticRealizationMap P hdeg : Fin 3 → MvPolynomial (Fin 3) K`;
3. proves its Jacobian determinant is `1`;
4. proves every coordinate has total degree at most `6 * P.natDegree + 2`;
5. constructs the literal target fiber of those three polynomial coordinates;
6. gives a natural equivalence from `AdjoinRoot P` to that literal fiber over
   every commutative test `K`-algebra;
7. proves that the represented quotient algebra has dimension
   `P.natDegree`; and
8. proves that the quotient algebra is finite étale over `K`.

The final construction statements are
`automaticRealizationMap_certificate`,
`automaticJacobianOneFiberRepresentingEquiv_natural`,
`automaticRepresentingAlgebra_etale`, and
`automaticRepresentingAlgebra_finite`.

For an admissible seed `G`, Lean also promotes `Π` and `B` to independent
parameters, proves the resulting inverse equation over `K(Π,B)(C)`
irreducible with degree and root-quotient finrank `G.natDegree`, proves the
three displayed coordinates algebraically independent, and constructs the
injective pullback on rational function fields. The principal declarations are
`generalGaugeFullyGenericInversePolynomial_certificate`,
`generalGaugeFullyGenericInverseAdjoinRoot_finrank`,
`generalGaugeMap_algebraicIndependent`, and
`generalGaugeFunctionFieldHom_injective`. For a supplied translation
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

The polynomial-presentation part of the prescribed-algebra realization theorem
is now formalized end to end.  The remaining steps needed for a single theorem
starting from an arbitrary finite étale algebra are:

1. identify the actual source-over-target pullback-field extension with the
   now-formal degree-`N` inverse-root extension over the iterated presentation
   `K(Π,B)(C)`, and connect that comparison to the definition of geometric
   degree;
2. formalize monogenicity of finite étale products over an infinite field;
3. compose monogenicity with the polynomial-presentation certificate;
4. either formalize the Campbell--Razar--Wright Galois case or keep it as a
   clearly isolated classical theorem interface;
5. formalize exact nonproperness, boundary-sheet accounting, and discriminant
   orders;
6. formalize symmetric monodromy and Hilbertian specialization;
7. formalize the nonarchimedean local-solubility, Chebotarev, Dirichlet, and
   prime-counting inputs used in the arithmetic corollaries if complete
   machine verification of those applications is desired. The explicit
   quintic's rational obstruction, real point, and three-adic point are already
   formalized.

The repository keeps these boundaries explicit so that a compiled certificate
is never described as proving a stronger statement than it does.
