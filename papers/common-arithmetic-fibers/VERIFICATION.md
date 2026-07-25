# Verification matrix

This file records the proof layer for every load-bearing claim in *Every
Nonzero Finite Étale Algebra Except Rank Two Is a Keller Fiber*.  “Lean” means
a theorem in the pinned `formal/finite-etale-keller` project with no `sorry`
and no project-specific axiom.  Symbolic checks are exact independent audits,
not substitutes for the corresponding mathematical argument.

| Claim | Paper proof | Lean certificate | Independent audit / external input |
|---|---|---|---|
| Existence and automatic choice of an admissible translation | Theorem 1 and its final paragraph | `Admissibility.lean`, `AutomaticRealization.lean` | Concrete translated seeds |
| Two-sided source/chart equivalence over arbitrary commutative rings | Proposition 2.1 | `SourceEquivalence.lean` | Universal source-chart identity |
| Derivative is a unit at every test-algebra root | Lemma 2.2 | `Bezout.lean`, `SeparableReconstruction.lean` | Quotient-ring regressions |
| Roots represent the complete abstract source fiber, naturally in the test algebra | Proposition 2.1 | `GaugeFiberPoints.lean`, `FiberNaturality.lean`, `UniversalFiber.lean` | Exact two-sided reconstruction |
| Uniform finite-sum assembly of the displayed coordinates | Subsection 2.2 | `GaugeAssembly.lean`, `GaugeInverseAssembly.lean` | Structural SymPy certificate |
| One arbitrary-degree `MvPolynomial (Fin 3) K` map and exact coordinate evaluation | Equations (2.1)–(2.3) | `GeneralGaugeMap.lean` | Structural checker and six-coefficient bridge |
| Generic inverse polynomial, explicit `β`, and derivative factorization | Equations (2.4)–(2.7) | `GeneralGaugeInverse.lean` | Concrete inverse-polynomial regressions |
| Actual displayed coordinate equations equal the represented source equations | Proposition 2.1 | `GeneralGaugeDisplayedFiber.lean` | SymPy and Singular expansions |
| Literal three-coordinate map fiber is represented, naturally over every test algebra | Proposition 2.1 | `GeneralGaugeRawFiber.lean` | Quotient-ring reconstruction |
| Zero-second-coordinate fiber is preserved by determinant-one output normalization | Corollary 2.4 | `GeneralGaugeNormalization.lean` | Exact scalar-normalization checks |
| Translation `K[S]/(P(a+S)) ≃ K[T]/(P)` | Theorem 1 | `TranslationQuotient.lean`, `GeneralGaugeRealization.lean` | Concrete translated examples |
| Final automatic literal-fiber realization from squarefree `P` | Theorem 1 | `GeneralGaugeRealization.lean`; final theorem `automaticJacobianOneFiberRepresentingEquiv_natural` | Full Lean build and axiom reports |
| Constant determinant `-2` and normalized determinant `1` for the arbitrary-degree actual map | Proposition 2.1 and Corollary 2.4 | `GeneralGaugeJacobian.lean` | Structural SymPy audit; independent generic degree-six Singular audit; concrete degrees 3–5 |
| Final coordinate-degree bound `6N+2` in terms of `N = deg P` | Corollary 2.4 and Theorem 1 | `GeneralGaugeDegree.lean`, `GeneralGaugeRealizationDegree.lean` | Structural termwise audit; Singular degree-six profile |
| Generic inverse irreducibility and geometric degree `N` | Proposition 2.1, primitive linear-parameter argument | Not yet formalized | Gauss-lemma proof; concrete degree 3–5 regressions |
| Every finite étale algebra over an infinite field is monogenic | Lemma 3.1 | Not yet formalized | Discriminant/Vandermonde proof |
| No characteristic-zero Keller map has generic degree two | Lemma 1.2 | Not formalized | Campbell–Razar–Wright plus faithfully flat descent |
| Optimal rank-five Hasse failure | Theorem 4.2 | Explicit polynomial, Bézout, and map certificates | Exact quintic checker plus local-solubility proof |
| One fixed map has infinitely many Hasse-failing fibers | Theorem 5.1 | Not yet formalized | Exact inverse-family checker plus local and prime-counting proof |

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
   every commutative test `K`-algebra.

The final construction statements are
`automaticRealizationMap_certificate` and
`automaticJacobianOneFiberRepresentingEquiv_natural`.

## Remaining formal boundary

The polynomial-presentation realization theorem is now formalized end to end.
The remaining steps needed for a single theorem starting from an arbitrary
finite étale algebra are:

1. formalize the primitive linear-in-`C` proof that the generic inverse
   polynomial is irreducible of degree `N`, and connect it to the function-field
   definition of geometric degree;
2. formalize monogenicity of finite étale products over an infinite field;
3. compose monogenicity with the polynomial-presentation certificate;
4. either formalize the Campbell--Razar--Wright Galois case or keep it as a
   clearly isolated classical theorem interface;
5. formalize the local-solubility, Chebotarev, Dirichlet, and prime-counting
   inputs used in the arithmetic corollaries if complete machine verification
   of those applications is desired.

The repository keeps these boundaries explicit so that a compiled certificate
is never described as proving a stronger statement than it does.
