# Native `074d9` record-twist Mordell--Weil audit

Date: 2026-09-04

<!-- status-consumer: EC-K3-R17-074D9-RECORD-TWIST-MW-OBSTRUCTION c794f827e9a8ac36 -->

## Status

The proposed full exceptional-quotient lift is impossible for both record
fibres.  Eight exact good reductions give

| twist | known rank | geometric rank upper bound | current `QQ(u)` rank |
|---|---:|---:|---|
| `q_tt04b07` | 1 | 4 | `UNKNOWN` in `[1,4]` |
| `q_tt11a44` | 1 | 4 | `UNKNOWN` in `[1,4]` |
| `q_tt11279` | 1 | 4 | `UNKNOWN` in `[1,4]` |
| `q_tt080fa` | 1 | 2 | `UNKNOWN` in `[1,2]` |

Consequently the two full twist groups can specialize with rank at most eight
at curve 356 and at most six at curve 385.  Both displayed exceptional
quotients have rank twelve.  Therefore

\[
Q_{356}\ne \operatorname{sp}MW(E^{q_{\tt04b07}})
            +\operatorname{sp}MW(E^{q_{\tt11a44}}),
\]

and likewise

\[
Q_{385}\ne \operatorname{sp}MW(E^{q_{\tt11279}})
            +\operatorname{sp}MW(E^{q_{\tt080fa}}).
\]

Thus no partition such as `5+7=12` can occur for these pairs.  At least four
directions at 356 and six directions at 385 are not specializations of these
two global quadratic-character blocks.  The stronger conclusion that all four
twists have rank one, leaving ten fibre-specific directions at each record,
remains `UNKNOWN`.

## Exact good-reduction bounds

For each twist use the regular surface

\[
q(u)y^2=x^3+A(u)x+B(u),
\]

birational on the generic fibre to the short twist

\[
Y^2=X^3+q(u)^2A(u)X+q(u)^3B(u).
\]

The Newton polytope has vertices

\[
(0,0,0),(12,0,0),(0,3,0),(0,0,2),(2,0,2).
\]

At each of `p=131,137,151,157,167,173,181,193`, the pinned open-source
ToricControlledReduction commit
`74cda9e8148cd8e9a3928fc15a558c9a70b67cc1` certifies nondegeneracy and
computes the degree-28 primitive Frobenius polynomial with Hodge vector
`[2,24,2]`.  Restoring the six-dimensional algebraic toric complement and
subtracting the geometric trivial lattice `U+D4+D4` gives

\[
\operatorname{rank}MW_{\overline{\mathbf F}_p(u)}
\le m_p+6-10,
\]

where `m_p` is the multiplicity of primitive eigenvalues equal to `p` times a
root of unity.  The resulting upper-bound sequences, in the prime order above,
are

| twist | bounds by prime | minimum |
|---|---|---:|
| `q_tt04b07` | `4,4,4,4,4,4,4,4` | 4 |
| `q_tt11a44` | `4,4,4,6,4,6,6,6` | 4 |
| `q_tt11279` | `4,6,6,4,6,6,4,4` | 4 |
| `q_tt080fa` | `6,6,2,4,2,4,2,8` | 2 |

Only the geometric bounds are promoted.  The six omitted toric classes have
not been assigned an arithmetic Frobenius action, so no smaller `+p`
multiplicity is used as an arithmetic rank bound.

## Exact sections and record-fibre images

For each rigid branch point `P` over `z^2=q(u)`, the anti-invariant point

\[
R=P-\sigma(P)=2P-\operatorname{Tr}(P)
\]

descends to an exact section of the short twist.  Direct substitution proves
`deg(X)=6`, `deg(Y)=9`, `P.O=0`, and height six.  Specialization satisfies

\[
\operatorname{sp}(R)=2P_+-\operatorname{sp}(\operatorname{Tr}(P)).
\]

The trace lies in the specialized generic `MW17`, so the image in the displayed
exceptional quotient is twice the rigid branch class.  In the basis
`P18,...,P29`, the four exact images are

| twist | record | quotient image |
|---|---:|---|
| `q_tt04b07` | 356 | `(0,0,0,0,0,0,0,2,0,0,0,0)` |
| `q_tt11a44` | 356 | `(0,-2,-2,0,-2,0,0,2,0,2,0,0)` |
| `q_tt11279` | 385 | `(0,0,0,2,0,0,0,0,0,0,0,0)` |
| `q_tt080fa` | 385 | `(0,0,0,2,0,0,2,-2,0,0,0,0)` |

Each record pair spans rank two.  These are the currently certified images;
they are not asserted to generate either full twist group.

## Function-field 2-descent boundary

Magma 2.29-9 implements `TwoSelmerGroup` over finite rational function fields,
but not over characteristic-zero function fields.  At the first globally good
prime `p=131`, all four exact jobs certify an irreducible 2-division cubic and
then exceed the public calculator's 60-second limit before returning a Selmer
group.  They therefore produce no good-reduction descent bound.

For diagnostics, the same jobs complete at coefficient-wise valid smaller
primes where surface fibres collide.  The 2-division cubic remains irreducible
and the Selmer dimensions are respectively `6,4,5,4` at primes `19,19,19,31`.
These four values are discovery data only and do not bound the characteristic-
zero ranks.

## Increasing-`P.O` section search

For every `k=0,1,2`, the exact finite-field search compiler writes the
cleared equation

\[
Y^2=X^3+q^2AXH^4+q^3BH^6,
\]

with

\[
\deg H=k,\qquad \deg X\le6+2k,\qquad \deg Y\le9+3k.
\]

Here `H` is monic on a chart where the section is affine.  The exporter uses
`k+1` distinct smooth chart fibres without rational 2-torsion, an exhaustive
cover because a degree-`k` denominator cannot vanish at all `k+1` fibres.  The
high coefficients recursively eliminate `Y`, leaving `6+3k` variables and
`9+3k` equations per leading-point block.

The small-prime `k=0,1,2` systems are discovery sieves, not rank certificates.
The first complete solve, for `q_tt04b07` at `p=19,k=0`, classifies all twelve
blocks: four have unit ideal and eight have nonunit leading ideal.  Rational
point extraction identifies the already known reduction among the retained
blocks; exact lifting of every retained block and complete `k=1,2` solves have
not closed.  No new characteristic-zero section is promoted.

## Certificates

- `artifacts/generated-results/elkies-k3-r17-074d9-twist-good-reduction-bounds-v1.json`
- `artifacts/generated-results/elkies-k3-r17-074d9-record-twist-sections-v1.json`
- `artifacts/generated-results/elkies-k3-r17-074d9-twist-2descent-audit-v1.json`
- `artifacts/generated-results/elkies-k3-r17-074d9-record-twist-mw-contribution-v1.json`

The first, second, and fourth certificates prove the negative quotient tests.
The descent and section-ladder records preserve the unresolved exact-rank
boundary rather than inferring ranks from bounded searches.
