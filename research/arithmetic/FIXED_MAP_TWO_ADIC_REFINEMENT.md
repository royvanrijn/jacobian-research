# Two-adic refinement for integral targets of the fixed map

This note completes the local test at two on two large parts of the
matching-resolvent locus and narrows its remaining odd/odd branch. It uses
[ILO1, ILO4 and ILO5](FIXED_MAP_INTEGRAL_LOCAL_OBSTRUCTIONS.md), including the
proved local identification of a totally ramified cubic extension of `Q_2`.
It makes no integral Hasse-existence or nonexistence claim.

Write an integral target as `(u,b,c)`, put `pi=-u!=0`, and set

\[
 F(S)=32E(S)=32\pi^5S^5-48\pi^4S^4+48\pi S^3
                 -(40\pi+9b)S^2+18S-9c.                       \tag{1}
\]

Let `R` denote the rational matching-resolvent locus: `E` is squarefree,
its rational irreducible factors have degrees two and three, and their
quadratic resolvents agree. Every Hasse-failing target lies in `R` by FMA3.
The statements conditional on `R` below concern local solubility at two;
they do not assert local solubility at the other places. Put
`e=v_2(pi)`, `f=v_2(b)`, with `f=infinity` if `b=0`.

## Theorem TAR1: exact criterion for even pi

On `R`, when `e>=1`, the entire fibre has a `Q_2`-point precisely in the
following cases:

| `e` | Allowed `b` |
|---|---|
| 1 | `v_2(b)=1` or `8|b` |
| 2 | `b` odd or `4|b` |
| 3 | `b` even |
| `>=4` | every integer `b` |

Every allowed case gives a simple local root even without the hypothesis
that the target lies in `R`. The excluded pairs `(e,f)` are exactly
`(1,0),(1,2),(2,1),(3,0)`. The conclusions allow arbitrary integer `c`.

**Proof.** The excluded cases have the irreducible totally ramified cubic
Newton block of ILO4/ILO5. On `R`, it must be the cubic factor, whose
quadratic resolvent is the nonsplit field `Q_2(sqrt(-3))`. Matching the
quadratic factor therefore leaves neither factor a local root. This gives
necessity, not just a failure of one sufficient criterion.

We verify every allowed case. If `f=0`, the allowed values are `e=2` or
`e>=4`, covered by the explicit simple-root rescalings in ILO4. If
`f=1,e=1`, use the rescaling in ILO5. For `f=1,e>=3` or `f=2,e>=2`,
we have `v_2(40pi+9b)=f`, and there is a maximal adjacent Newton edge

\[
 (2,f)\longrightarrow(3,4+e).                                 \tag{2}
\]

Indeed its supporting line at indices `0,1,4,5` has respective heights
`3f-2e-8`, `2f-e-4`, `2e+8-f`, `3e+12-2f`. These are strictly below the
actual coefficient valuations `v_2(c)>=0`, `1`, `4+4e`, `5+5e` in the
stated ranges. ILO1 supplies a simple local root.

Finally if `f>=3`, including `b=0`, the polynomial `F(Z/2)` is integral
and reduces to `Z+c` modulo two: the valuations of its quadratic and
higher coefficients are positive. Its unique residue root is simple.
This covers all remaining entries of the table. Every simple local root
survives derivative localization, proving sufficiency for the actual
entire fibre. \(\square\)

## Theorem TAR2: exact criterion for odd pi and even c

On `R`, if `pi` is odd and `c` is even, the entire fibre has a `Q_2`-point
if and only if `4|b`. The sufficient direction again holds without the
hypothesis `R`.

**Proof.** For `f=0`, ILO4 excludes a local point on `R`; for `f=1`, use
the `(e,f)=(0,1)` case of ILO5 and its same local resolvent argument.
If `4|b`, the integral polynomial `F(Z/2)` reduces to

\[
 Z^5+Z^4+\epsilon Z^2+Z,
 \qquad\epsilon=(b/4)\bmod2.
\]

The residue root zero is simple, as in ILO5. \(\square\)

## Theorem TAR3: the remaining odd pi, odd c branch

Suppose `pi,c` are odd. A target in `R` with a `Q_2`-point, and hence any
Hasse-failing target, must have `8|b`. Define the integer

\[
 A_0=F(1/2)
     =\pi^5-3\pi^4-4\pi+9-9c-\frac{9b}{4}.                    \tag{3}
\]

It is nonzero on `R`, and at least one of the following patterns is necessary:

\[
 \begin{array}{c|c|c}
 &\text{conditions on }\pi,b&\text{condition on }A_0\\\hline
 \mathrm{I}&\pi\equiv5\pmod8,\quad16\mid b&8\mid A_0\\
 \mathrm{II}&\pi\equiv1\pmod{16},\quad32\mid b&32\mid A_0.
 \end{array}                                                   \tag{4}
\]

In Case I the entire fibre always has a simple `Q_2`-point, even without
`R`. The present theorem gives only necessary conditions in Case II.

**Proof.** ILO4/ILO5 exclude `f=0,1`. When `f=2` and `c` is odd, ILO5
makes the quintic irreducible over `Q_2`, which is incompatible with `R`.
Thus `f>=3`. If `A_0=0`, `S=1/2` is a rational root, contradicting the
irreducible `2+3` factorization on `R`.

Write `G(T)=F(T+1/2)=sum g_i T^i`. Exact expansion gives

\[
 \begin{aligned}
 g_0&=A_0,\\
 g_1&=2(\pi-1)(5\pi^4-7\pi^3-7\pi^2-7\pi-9)-9b,\\
 g_2&=8\pi(\pi-1)(5\pi^3-4\pi^2-4\pi-4)-9b,\\
 g_3&=16\pi(5\pi^4-6\pi^3+3),\\
 g_4&=16\pi^4(5\pi-3),\qquad g_5=32\pi^5.                 \tag{5}
 \end{aligned}
\]

The two bracketed factors in `g_1,g_2` are odd. For odd `pi`,
`v_2(g_3)=v_2(g_5)=5`; also `v_2(g_4)=5` if `pi=1 mod4`, and at least
six if `pi=3 mod4`. Since `c` is odd and `8|b`, `A_0` is even. Put
`t=v_2(pi-1)`, allowing `t=infinity` at `pi=1`, and
`h=v_2(A_0)>=1`.

The following valuation table exhausts all cases with `f>=3` except
`t=2,f>=4` and `t>=4,f>=5`. Each listed edge is a maximal lower edge;
its intermediate and exterior coefficients are strictly above its
supporting line. Endpoints use the index and valuation of `g_i`.

| Conditions | Range of `h` | Consequence from the lower hull |
|---|---|---|
| `t=1`, any `f>=3` | `h=1` | edge `(0,1)` to `(5,5)`, irreducible quintic |
| `t=1`, any `f>=3` | `h>=2` | edge `(1,2)` to `(5,5)`, irreducible quartic block |
| `t>=2`, `f=3` | `h=1` | edge `(0,1)` to `(5,5)`, irreducible quintic |
| `t>=2`, `f=3` | `h>=2` | edge `(2,3)` to `(5,5)`, totally ramified cubic block |
| `t>=3`, `f=4` | `h<=3` | edge `(0,h)` to `(5,5)`, irreducible quintic |
| `t>=3`, `f=4` | `h>=4` | edge `(2,4)` to `(5,5)`, totally ramified cubic block |
| `t=3`, `f>=5` | `h<=3` | edge `(0,h)` to `(5,5)`, irreducible quintic |
| `t=3`, `f>=5` | `h>=4` | edge `(1,4)` to `(5,5)`, irreducible quartic block |

To check the table directly from (5), the two terms in `g_1` before
subtraction have valuations `t+1,f`, and in `g_2` they have valuations
`t+3,f`. In the finite equality case `f=t+1`, cancellation can raise
only the former sum; then
`v_2(g_2)=f`. If `f=t+3`, only the latter sum can be raised; then
`v_2(g_1)=t+1`. These observations give all stated strict supporting-line
inequalities, including the equality cases in those two coefficient sums.

An edge of length four or five in the table has slope with denominator
four or five and therefore gives an irreducible factor of that degree
over `Q_2`. Such a factor cannot occur in a rational quadratic-times-cubic
polynomial. The cubic edges have slope with denominator three and are
excluded for a local point on `R` by the same totally ramified cubic
resolvent argument as ILO4. This proves all exclusions in the table.

For `t=2,f>=4`, formula (5) gives `v_2(g_1)=3`, `v_2(g_2)>=4`, and
`v_2(g_3)=v_2(g_4)=v_2(g_5)=5`. If `h=1` or `2`, the single edge from
`(0,h)` to `(5,5)` gives an irreducible quintic. Thus `h>=3` is necessary,
which is Case I. In that range there is an adjacent lower edge
`(0,h)` to `(1,3)`, producing a simple root by ILO1. If `A_0=0` outside
the hypothesis `R`, the root zero of `G` is simple because `g_1!=0`.

For `t>=4,f>=5`, both `g_1,g_2` have valuation at least five, while
`g_3,g_4,g_5` have valuation five. If `h<=4`, the single edge from
`(0,h)` to `(5,5)` again gives an irreducible quintic. Hence `h>=5`,
which is Case II. These two surviving cases are exactly (4).
\(\square\)

## What the refinement does and does not decide

Theorems TAR1/TAR2 are exact local criteria on the stated rational
factorization and resolvent locus. Their excluded cases sometimes have
local roots without that hypothesis; it is the compatibility with the
rational factorization that rules them out. No unconditional classification
of all integral target fibres is being substituted for the stated result.

Case I of TAR3 has an actual local point, but can still fail the global
factorization, another local condition, or rational-point exclusion.
Case II remains unresolved by this valuation analysis. Given an individual
factorized target, the exact finite local algorithm of FMA3 still applies.
Neither case proves that an integral Hasse-failing target exists, and
excluding the other cases does not prove integral nonexistence.

The Newton and valuation inputs are those already identified precisely in
the ILO note: Milne's retained *Algebraic Number Theory*, v3.00,
Proposition 7.44, Remark 7.43 and Corollary 7.42. Simple-root Hensel lifting
is Proposition 7.31. The lower-hull convention here is by ascending
coefficient index, as in ILO4.

## Verification

Run from the repository root:

```sh
.venv/bin/python research/arithmetic/verify_fixed_map_two_adic_refinement.py --write-receipt
```

The [checker](verify_fixed_map_two_adic_refinement.py) verifies the exact
translated coefficients, residue controls for automatic roots, bounded
Newton hulls in all stated even-coordinate cases, and an independent
bounded exclusion check for the odd/odd valuation table. The
[receipt](fixed_map_two_adic_refinement_checks.json) records its finite
domains and hashes. These checks do not prove the unrestricted valuation
arguments, the imported Newton theorem, or integral Hasse nonexistence.
No formal verification, independent whole-proof replay, external review,
or graph-status change is claimed by this source.
