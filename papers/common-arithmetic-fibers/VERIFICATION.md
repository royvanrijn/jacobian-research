# Verification matrix

This file records the proof layer for every load-bearing claim in *Every Finite
Étale Algebra Except Rank Two Is a Keller Fiber*.  “Lean” means a theorem in the
pinned `formal/finite-etale-keller` project with no `sorry` and no
project-specific axiom.  Symbolic checks are exact independent audits, not
substitutes for the corresponding mathematical argument.

| Claim | Paper proof | Lean certificate | Independent audit / external input |
|---|---|---|---|
| Existence of an admissible translation | Theorem 1 and its final paragraph | `Admissibility.lean`, `AutomaticRealization.lean` | Concrete translated seeds |
| Source/chart equivalence over arbitrary commutative rings | Proposition 2.1 | `SourceEquivalence.lean` | Universal source-chart identity |
| Derivative is a unit at every test-algebra root | Lemma 2.2 | `Bezout.lean`, `SeparableReconstruction.lean` | Quotient-ring regressions |
| Roots represent the complete source fiber, naturally in the test algebra | Proposition 2.1 and Theorem 1 | `GaugeFiberPoints.lean`, `FiberNaturality.lean`, `UniversalFiber.lean` | Exact two-sided reconstruction |
| Translation `K[S]/(P(a+S)) ≃ K[T]/(P)` | Theorem 1 | `TranslationQuotient.lean`, `RealizationFiber.lean` | Concrete translated examples |
| Automatic polynomial-level represented-fiber theorem | Theorem 1 | `AutomaticRealization.lean` | Full Lean build and axiom report |
| Uniform finite-sum assembly of the displayed coordinates | Subsection 2.2 | `GaugeAssembly.lean` | `verify_universal_quadratic_gauge.py` |
| One all-degree `MvPolynomial (Fin 3) K` map and exact coordinate evaluation | Equations (2.1)–(2.3) | `GeneralGaugeMap.lean` | Structural checker and degree-six bridge |
| Constant determinant and determinant-one normalization | Proposition 2.1 and Corollary 2.4 | Universal cancellation core; explicit quintic map | Structural Jacobian audit and concrete map regressions |
| Coordinate-degree bound `6N+2` | Corollary 2.4 | Not yet packaged as a `totalDegree` theorem | Structural termwise audit |
| Generic degree `N` | Proposition 2.1, resultant paragraph | Not yet formalized | Paper resultant argument; concrete degree 3–5 regressions |
| Every finite étale algebra over an infinite field is monogenic | Lemma 3.1 | Not yet formalized | Discriminant/Vandermonde proof |
| No characteristic-zero Keller map has generic degree two | Lemma 1.2 | Not formalized | Campbell–Razar–Wright plus descent |
| Optimal rank-five Hasse failure | Theorem 4.2 | Explicit polynomial and map certificate | Exact quintic checker plus local-solubility proof |
| One fixed map has infinitely many Hasse-failing fibers | Theorem 5.1 | Not yet formalized | Exact inverse-family checker plus local and prime-counting proof |

## Remaining formal boundary

The intended end state is one theorem which starts with a finite étale
`K`-algebra of rank `N ≥ 3` and produces the displayed Jacobian-one map, its
geometric degree, coordinate-degree bound, and the natural scheme-fiber
isomorphism.  The main remaining modules are:

1. map-level Jacobian and `totalDegree` theorems for `GeneralGaugeMap.lean`;
2. the nonzero-resultant/generic-degree theorem;
3. monogenicity for finite étale products over an infinite field;
4. the final composition theorem;
5. formalization or an explicitly isolated axiom/theorem interface for the
   classical degree-two Galois result and the number-theoretic applications.

The repository keeps these boundaries explicit so that a compiled lower-level
certificate is never described as proving a stronger statement than it does.
