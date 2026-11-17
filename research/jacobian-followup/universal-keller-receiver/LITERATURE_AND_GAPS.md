# Receiver literature, completed gates and remaining limits

14 November 2026. This continues the [trace-quadric literature transfer](../universal-keller-quadrics/LITERATURE_AND_GAPS.md) with a categorical construction and affine quotient obstructions.

| Primary source | Use here |
|---|---|
| [Arpin--Bozlee--Herr--Smith, Scheme of Monogenic Generators II, Sections3--4](https://arxiv.org/html/2205.04620v2) | Generator schemes, configuration spaces and scalar-twisted generators supply the established framework. Our H_A is the specific denominator-open chart for the Keller compiler. |
| [Stacks Project, quotient stacks](https://stacks.math.columbia.edu/tag/04UI) | Quotients by commuting torsor actions and the interpretation of objects as torsors with equivariant maps. |
| [Stacks Project, inertia and quotient stacks](https://stacks.math.columbia.edu/tag/06PA) | Stabilizers are the automorphism groups of quotient objects. A section must induce a split map on those groups. |
| [Stacks Project, etale morphisms of stacks](https://stacks.math.columbia.edu/tag/0CIK) | Etaleness of our representable quotient morphism is checked after the smooth target presentation. Relative dimension zero alone would not justify it. |
| [Alper, Good moduli spaces for Artin stacks](https://arxiv.org/abs/0804.2242) and [On the local quotient structure of Artin stacks, Section6](https://sites.math.washington.edu/~jarod/papers/quotient-structure.pdf) | Affine invariant quotients are not geometric orbit spaces. Etale stack morphisms do not automatically give etale maps of those quotients. Our explicit contracted curve verifies the failure directly. |
| [Brassil--Reichstein, Hermite--Joubert problem over p-fields, Corollary8.1](https://www.math.uni-bielefeld.de/LAG/man/560.pdf) | The rank-six two-trace receiver does not remove its negative binary arithmetic test in characteristic zero. We do not infer universality from its dimension or smoothness. |

## Gates completed

The previously proposed identification of the full trace-root quotient receiver is now proved, with its exact inertia. In fact the unconstrained receiver works for every N>=3. The trace receiver includes rank five without requiring the old d4*d5 denominator chart; two traces can be retained at rank six before a polynomial coefficient section exists.

The resulting no-section statement is stronger in scope than the earlier algebraic-space inertia observation: these explicit stacks and all enhancements mapping to them still have too small an image on algebra automorphisms. This does not obstruct every conceivable stack receiver.

The global coarse-quotient and C=1 normalization ideas now have written, uniform failure certificates. The former contracts Gamma; the latter retains a nonconstant unit. Neither requires a high-degree search.

## Candidate-screening rules

1. State whether the receiver retains a marked fibre, its algebra identification, and a generator line. A fixed unmarked universal map is a different object.
2. Before seeking a section over BS_N, verify that its automorphism image can contain the full S_N. A kernel of irrelevant gauge automorphisms cannot substitute for this image.
3. Before interpreting a quotient dimension as a Keller target dimension, identify the actual source and target schemes and prove a polynomial extension on the whole affine spaces.
4. Check quasi-finiteness of the coarse map. Generic degree N does not prevent an exceptional positive-dimensional fibre.
5. Check units and singularities of the actual chosen charts. Smoothness, rationality or a quotient presentation does not imply affine space.
6. Keep the finite target chart separate from the full affine quotient: removing a contracted locus or filling an omitted divisor changes the problem.

## Remaining work

- UAFROB1 gives effective polynomial incidence over B_AGL1(F_q) on an affine source bundle without an origin, including BS3 in characteristic three. Characteristic-zero and general S_N receivers forgetting the good generator line remain open.
- The later UABLEVEL1 [coordinate-level theorem](../universal-keller-levels/PROOF.md) classifies B=constant for independent promoted seeds in every odd characteristic, including its precise characteristic-three affine exception. Relations among coefficients and other completions are separate source problems.
- Filling the divisors in the even-rank charts requires a separate construction. Existing polynomial seed-only sections show why a unit on an open chart is not itself an obstruction to every possible completion.
- Low-rank A3 universality and arbitrary-base dimension minima remain governed by both arithmetic and affine geometry.
- Independent mathematical review of the new proofs is pending; the checker is not a stack formalization.

## A separate necessary condition after forgetting the line

UAEQUIV1 proves directly that a polynomial equivariant Keller map restricts to a Keller map on each subgroup-fixed linear space, with the exact fixed-label fibre. It requires no classification of finite group representations. For the permutation bundle alone the full fixed space has dimension one and the natural N-point fibre has no full-group-fixed point, a contradiction. Adding one trivial line makes that restriction a plane Keller counterexample; two extra trivial lines are the first such architecture not excluded here. This condition concerns natural marked representation-bundle lifts and supplies no global lower bound for arbitrary universal maps.

## From common power shifts to arbitrary profiles

The retained UPGM1 theorem already handles a common nonnegative extra power of Pi. UAPROFILE1 proves the general Jacobian-preserving substitution rule and classifies every rational coefficient profile that is polynomial in the original source, first excluding finite nonzero poles. Its characteristic-zero monodromy proof uses the quadratic-tilt envelope, transposition generation and [Riemann--Hurwitz](https://stacks.math.columbia.edu/tag/0C1B). Negative profiles have a precise characteristic-p exception, so the old Pi-zero fibre theorem must not be reused on those exceptional extensions.
