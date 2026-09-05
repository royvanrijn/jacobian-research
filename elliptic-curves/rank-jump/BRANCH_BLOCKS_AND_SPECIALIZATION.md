# Branch characters and exact specialization collapse

The 37 already-split bisections on five historic fibres do not form a
three-cover genus-zero or genus-one block. Their branch quadratics are
irreducible, distinct and pairwise coprime. **Every retained triple has
genus five.** More significantly for interpreting simultaneous solubility,
the 25 soluble covers at the lower-gain control specialize to only four
independent rational quotient directions.

The branch conclusion is a retrospective replay of existing
[Theorems F2–F4](../../elkies-k3/RANK_MUTATION_AND_LIFT_THEOREMS.md);
those theorems already establish disjoint branching for the full atlas.
It is not a new exclusion of the general low-genus-carrier mechanism.
The additional calculation here solves the stored exact group relations
over the rationals and records the specialization kernels and proportional
direction blocks. It does not identify global relations from finite
Kummer fingerprints alone.

## Frozen experiment and evidence

The [protocol](BRANCH_BLOCK_PROTOCOL.json) selects every hit in the existing
[five-fibre census](../notes/ELKIES_BISECTION_SPECIALIZATION_CONTROLS.md).
There are 37 retained records, 319 pairs and 2,321 triples. No atlas was
rebuilt, parameter tested for a new hit, or new point sought.

The [compact inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_branch_block_inputs_v1.json)
pin the census bytes, full-atlas hash and exact selected equations.
The [certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_branch_blocks_v1.json)
retains branch factors, rational quotient coordinates, integer kernel
vectors, all pair ranks and triple genus counts. Per-fibre checkpoints
are under the ignored local-artifact directory.

    python3 elliptic-curves/rank-jump/branch_blocks.py check
    sage -python elliptic-curves/rank-jump/branch_blocks.py verify

The portable replay checks all 37 discriminant identities
`b^2-4ac=h^2 q`, nonzero square values and roots, branch character ranks,
and exact rational elimination of the census relation matrices. Sage
independently factors the quadratics, checks every pair gcd and repeats
the relation-matrix inversion. Exact elliptic group-law verification of
the input relations is inherited from the pinned census; it is not
rerun by these commands.

## Geometry versus specialized directions

Let `k` be the number of retained covers and `C_k` their connected,
smooth projective fibre-product normalization. All five cases have
geometric and arithmetic character rank `k`, with `2k` distinct branch
points. The generic rank bounds below are over **Q(C_k)**; they are not
new generic ranks over Q(t).

| Parameter | Known quotient over MW17 | Soluble covers k | Branch points | Genus C_k | Generic rank lower bound over Q(C_k) | Exact quotient rank of these lifts | Kernel modulo specialized M_Q |
|---|---:|---:|---:|---:|---:|---:|---:|
| -2/377 | 8 | 6 | 12 | 129 | 23 | 5 | 1 |
| -308/251 | 9 | 3 | 6 | 5 | 20 | 3 | 0 |
| 2456/135 | 10 | 2 | 4 | 1 | 19 | 2 | 0 |
| -9529/5471 | 11 | 1 | 2 | 0 | 18 | 1 | 0 |
| 3/8 | 4 | 25 | 50 | 385,875,969 | 42 | 4 | 21 |

Every normalization has a rational point above the retained parameter:
all its square roots are supplied and nonzero, so this point is away
from the branch divisor. High genus here is compatible with a known
rational point. It neither proves rank nor supplies a probability of
simultaneous solubility.

For `n` independent geometric characters and `b` branch points, every
nontrivial inertia group has order two. Riemann–Hurwitz gives

    2g-2 = -2*2^n + b*2^(n-1),
    g = 1 + 2^(n-2)*(b-4).

This is the elementary-abelian case of the branch-matrix construction
in [Frediani, Section 2](https://link.springer.com/article/10.1007/s00229-024-01556-0).
Here each quadratic supplies two private branch points, hence `n=k`,
`b=2k`, and `g=1+2^(k-1)*(k-2)`.
More generally, a rational unbranched simultaneous lift rules out a
nonsquare constant relation: evaluating `product(q_i)=c*h^2` at that
point forces `c` to be a rational square.

The atlas gives a non-torsion anti-invariant section
`T_i=P_i-sigma_i(P_i)=2P_i-tau_i` of height 12 on each double cover.
Theorem F4 places them in distinct character spaces, with height block

    12*2^(k-1) * I_k

on `C_k`. Thus their generic independence is established. At a chosen
rational point of `C_k`, specialization need not preserve that independence
modulo the specialized generic subgroup. The table measures this second
map exactly. The full curve ranks and the full generic rank over `C_k`
remain UNKNOWN beyond their stated lower bounds.

## The collapse blocks are rational, not just mod two

At `-2/377`, write `P_label` for the positive branch used in the census.
The unique relation among the six quotient directions is

    -P_1cb25 + P_0cff7 - P_1ea09 + P_0d4ca = 0 modulo M tensor Q.

Every pair remains independent. Thus the dependency is a four-direction
specialization event, although its four branch characters remain
geometrically independent. The other two lifts supply separate directions.

At `3/8`, the 25 lifts have the following exact rational line groups in
the pinned public complement `(Q1,Q2,Q3,Q4)`:

| Rational line | Number of lifts |
|---|---:|
| Q1 | 8 |
| Q2 | 5 |
| Q1+Q2 | 3 |
| Q3-Q4 | 2 |
| Q2-Q4 | 1 |
| Q2+Q3-Q4 | 1 |
| Q3 | 1 |
| Zero modulo M_Q | 4 |

The four zero directions are `orbit-0be21`, `orbit-10aaa`,
`orbit-06f04`, and `orbit-126e6`. Here zero means membership in the
rational span of M, not an integral membership or torsion assertion.
The certificate gives a 21-dimensional integer-coefficient kernel in
the 25-dimensional space of displayed sections modulo M_Q.

All 300 paired bases at this control have the same generic new height
block `24 I_2`, and all have rational points. Nevertheless, 168 pairs
give quotient rank two, 126 give rank one, and six give rank zero at
this specialization. In contrast, every one of the 19 retained pairs
on the three larger-jump fibres with at least two hits gives rank two.

This explains one source of overcounting: the number of distinct soluble
geometric covers can greatly exceed the number of independent
specialized directions. It does **not** explain why the whole fibre
has a large jump, since this atlas misses most directions on the
largest-jump controls.

## Mechanism ranking and next test

1. **Solubility: a common auxiliary curve with rational points remains
   a valid mechanism for blocks.** The retained genus-one pairs realize
   it for two directions, with generic independence already proved.
   These pairs explain only two of the ten known quotient directions
   at `2456/135`; they are not an explanation of that entire jump.
2. **Incidence: additional character spaces on a low-genus base remain
   a structural opening.** On a paired base, a section of the product
   twist `E^(q_i*q_j)` contributes a third orthogonal character without
   adding branch points or increasing the genus. Theorem F4 already
   proves this implication. What is missing is such a section and its
   survival as a new direction at the specialization.
3. **Weak explanation: accumulating soluble singleton covers.** The
   control with 25 lifts has the smallest quotient in this comparison,
   and 21 generic directions disappear modulo its generic subgroup.
   Neither split count nor generic character count predicts the jump.
4. **Visibility: multiplicities within a quotient line.** The eight
   covers exposing Q1 on the control are multiple ways to see a
   direction, not eight rank events. This information is useful for
   interpreting recovery yields, not for selecting high-rank fibres.

Do not extend the branch-overlap test to more atlas subsets: F3–F4
already settle it. Do not launch another pair-base point campaign;
the existing [paired-base census](../../elkies-k3/BISECTION_PAIR_COVER_GEOMETRY_2026-08-31.md)
already contains thousands of rationally soluble bases. The existing
[product-twist scores](../../elkies-k3/QUADRATIC_TWIST_RANK_CENSUS_2026-08-31.md)
also provide no certified extra section.

The next bounded retrospective test should check the **mixed-character
possibility** against the retained exceptional-point cover equations:
freeze the already fitted genus-one quartics for the historic exceptional
generators, factor their branch polynomials exactly, and ask whether any
is a product character of two atlas quadratics. An irreducible quartic
cannot be such a product; a factorization match still requires its rational
constant squareclass and section identity to agree. This test uses no new
parameters and makes no claim that a point's fitted cover is unique.

The fitted covers use oracle points, so this first comparison is strictly
retrospective. A positive result would identify a concrete section and
character block to reconstruct without oracle data. A negative result
would exclude only those retained presentations, not all covers through
the same points. The unresolved implication remains:

    a point-independent specialization condition
        => several globally soluble, independent classes
        => a large rational Mordell–Weil jump.
