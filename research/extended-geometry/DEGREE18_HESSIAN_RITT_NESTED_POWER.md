# The degree-eighteen nested Hessian--Ritt component

Work over a field of characteristic zero on the monic twice-integrated
Hessian chart.  There are four ordered factor loci

\[
 \mathcal C_{2,9},\quad \mathcal C_{3,6},\quad
 \mathcal C_{6,3},\quad \mathcal C_{9,2}.
\]

This note records the first higher-degree counterexample search after the
complete degree-eight and degree-twelve diagrams.  It gives a negative answer
to the most literal refinement-tree formulation, but supports a slightly
richer formulation in terms of collision-labelled Ritt diagrams.

## 1. Pairwise degree-pattern census

The prime refinement words and the factor cuts which they carry are

| word | factor loci |
|---|---|
| `2 o 3 o 3` | `C_(2,9)`, `C_(6,3)` |
| `3 o 2 o 3` | `C_(3,6)`, `C_(6,3)` |
| `3 o 3 o 2` | `C_(3,6)`, `C_(9,2)` |

Thus three pair intersections have literal common refinements.  Two more are
the degree-six `2`-versus-`3` collision transported on the right or left by a
cubic.  The extreme pair `C_(2,9) cap C_(9,2)` is a direct coprime
`2`-versus-`9` Ritt collision.  No degree-eighteen pair requires more than one
Ritt core.

Ritt's second theorem puts the extreme pair in one source-translated power
family.  Unlike the degree-twelve `3`-versus-`4` pair, it has no separate
Chebyshev component: when the smaller degree is two, the Chebyshev case is
absorbed by the power normal form.  This is the overlap described in the
normal-form theorem of
[von zur Gathen](https://arxiv.org/abs/1308.1135).

## 2. Exact restriction of the middle factor orders

Put `y=W-s` and write the extreme power core as

\[
 f(W)=\bigl(yR(y^2)\bigr)^2,
 \qquad
 R(t)=t^4+r_3t^3+r_2t^2+r_1t+r_0.               \tag{1}
\]

Pulling the nine canonical residual equations for `C_(3,6)` and `C_(6,3)`
back to `QQ[s,r_0,r_1,r_2,r_3]` gives, in both cases, the same unique minimal
prime

\[
 r_2=\frac{r_3^2}{3},\qquad
 r_0=\frac{r_1r_3}{3}-\frac{r_3^4}{81}.          \tag{2}
\]

Exact minimal-prime computation gives one component of dimension three.
Consequently the reduced all-four intersection inside the extreme Ritt core
is an irreducible threefold, not merely a Chebyshev curve.  This calculation
does not assert that the unreduced intersection scheme has no nilpotents or
embedded structure.

The meaning of (2) is transparent.  Set

\[
 u=\frac{r_3}{3},\qquad
 v=r_1-\frac{r_3^3}{27},\qquad
 B(y)=y^3+uy,\qquad A(z)=z^3+vz.                 \tag{3}
\]

Then (2) is exactly the identity

\[
 yR(y^2)=A(B(y)),
 \qquad
 R(t)=(t+u)\bigl(t(t+u)^2+v\bigr).               \tag{4}
\]

Hence the component has the simple parametrization

\[
 \boxed{f_{s,u,v}(W)=\bigl(A(B(W-s))\bigr)^2.}   \tag{5}
\]

It carries all three prime refinement words.  Indeed, if

\[
 G_c(t)=t(t+c)^2,
\]

then

\[
 A(z)^2=G_v(z^2),\qquad B(y)^2=G_u(y^2).         \tag{6}
\]

Equations (5)--(6) give the chains

\[
 2\circ3\circ3,\qquad
 3\circ2\circ3,\qquad
 3\circ3\circ2,                                  \tag{7}
\]

and therefore all four two-factor cuts.  Conversely, the normal form for the
extreme pair and the exact pullback calculation show that every point of the
all-four intersection has (5).  On the affine monic chart its dimension is
three; the normalized endpoint incidence cuts its marked dimension to two.

The component meets the marked clean open.  One rational witness in (5) is

\[
 s=-3,\qquad u=5,\qquad v=-\frac{721476}{31}.     \tag{8}
\]

Its endpoint incidence vanishes, while the endpoint derivative gap, marked
Hessian value, both Hessian and primitive-root discriminants, and the
`H''(1)+2` condition are all nonzero.

## 3. Counterexample verdict

There is no pairwise counterexample in degrees `18`, `24`, or `30`: every pair
in the degree-pattern census is either a common refinement or one transported
coprime Ritt core.  The numbers are

| degree | factor loci | common-refinement pairs | one-core pairs | multi-core pairs |
|---:|---:|---:|---:|---:|
| 18 | 4 | 3 | 3 | 0 |
| 24 | 6 | 9 | 6 | 0 |
| 30 | 6 | 6 | 9 | 0 |

Degree eighteen nevertheless disproves the narrow higher-intersection
picture in which the exceptional power/Chebyshev pieces remain isolated
pairwise additions to a common-refinement lattice.  No single ordered degree
word carries all four cuts, but (5) is a positive-dimensional all-four
component obtained by two compatible power moves.

The reusable formulation should therefore be phrased in terms of **Ritt
diagrams**, not just common refinement trees:

> Irreducible components of multiple Hessian-composition intersections are
> compatibility loci of finite connected diagrams of normalized
> decompositions.  Vertices are common refinements and edges are labelled by
> power or Chebyshev Ritt moves.  Components are indexed only after diagrams
> with the same dense image are identified.

This is consistent with the global description of complete decompositions by
[Mueller and Zieve](https://arxiv.org/abs/0807.3578), but the irreducibility and
diagram-identification statement is additional algebraic geometry, not a
formal restatement of Ritt's theorems.

The later [general-intersection note](GENERAL_HESSIAN_RITT_INTERSECTIONS.md)
sharpens this verdict.  Ziegler's tame multi-collision theorem already gives
the relation-graph normal form for synchronized ordinary polynomial
compositions.  The new issue specific to this project is Hessian transfer:
one must prove that projection away from the linear term creates no
unsynchronized minimal primes.  This issue is now closed in degree eighteen:
canonical linear-lift ideal membership holds scheme-theoretically on all six
pairs of the four factor loci, and therefore on every higher intersection.
The family (5) is exactly the exponential block predicted by that relation
graph, with no remaining Hessian-synchronization caveat.

Degree `24` is the next compatibility test: its prime-refinement shadow is the
four-vertex path obtained by moving the cubic through three quadratic
factors.  Degree `30` is the first braid test: the six permutations of
`(2,3,5)` form the `S_3` Ritt graph, so distinct reduced move paths can have
the same Chebyshev image.  These are the correct searches after the
degree-eighteen nested-power calculation.

## 4. Reproduction

Run

```bash
.venv/bin/python scripts/explore_hessian_ritt_component_census.py
.venv/bin/python scripts/explore_degree18_hessian_ritt_power_core.py
```

The first command audits every factor-locus pair in degrees `18`, `24`, and
`30`.  The second constructs the canonical degree-eighteen residuals, checks
the factorization (4), and uses exact Singular minimal-prime computation to
certify (2).
