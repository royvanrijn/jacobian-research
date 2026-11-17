# Exact real chamber spectra of explicit Keller maps

This note strengthens the [all-degree rational-fiber theorem](ALL_DEGREE_RATIONAL_FIBERS.md).
In the terminology of
[Migus's generic-degree classification](https://arxiv.org/abs/2607.21572),
it refines the realization side inside the explicit weighted family:
Migus proves that non-dense real Keller images have even generic degree at
least four and that every such degree occurs, while the result below
determines the complete parity chamber spectrum, rational witnesses, and a
fold-adjacency chain for one explicit map in every degree.  The all-degree
weighted family and the first quartic empty-real-fiber example are due to
[Gallagher](https://jacobianfun.org/jacobian-explained); all-degree families
also occur in the earlier base manuscript.

The input is elementary real univariate geometry: once one member of a
degree-`N` pencil has `N` simple real roots, a generic parallel pencil line
loses those roots two at a time.

## Pencil theorem

Let `H in R[W]` have degree `N>=2` and put

\[
 E_{s,t}(W)=H(W)-sW+t,
 \qquad
 U_H=\{(s,t)\in\mathbb R^2:\operatorname{Disc}_W(E_{s,t})\ne0\}.
\]

On every connected component of `U_H`, the number of real roots is constant.
If one member `E_(s_0,t_0)` has `N` distinct real roots, then the set of root
counts occurring on components of `U_H` is exactly

\[
 \boxed{\{N,N-2,N-4,\ldots,N\bmod2\}}.                 \tag{1}
\]

More precisely, arbitrarily close to `s_0` there is a slope `s_1` for which
the vertical line `s=s_1` meets chambers with all the counts in (1), in that
order, and consecutive members are separated by one smooth fold of the
pencil discriminant.

If `H` and `(s_0,t_0)` are rational, `s_1` and one witness `t` in every one
of these chambers may all be chosen rational.

### Necessary parity

Off the discriminant all `N` complex roots are simple.  Nonreal roots occur
in conjugate pairs, so the real-root count `r` satisfies

\[
 0\le r\le N,
 \qquad r\equiv N\pmod2.                               \tag{2}
\]

This proves that no counts outside (1) can occur.

### A generic vertical line

Translate `t` so that `t_0=0` and write

\[
 P_s(W)=H(W)-sW.
\]

The property that `P_s` has `N` distinct real zeros is open, so it holds for
all `s` in a real interval about `s_0`.  Inside that interval avoid the
following finite set of slopes:

1. `s=H'(r)` with `H''(r)=0`;
2. slopes for which two distinct critical points of `P_s` have the same
   critical value.

The second set is finite because the discriminant normalization

\[
 r\longmapsto\bigl(H'(r),rH'(r)-H(r)\bigr)             \tag{3}
\]

is finite and birational onto its image.  Two equal critical values are
precisely two distinct normalization points over one discriminant point, and
a finite birational map of curves fails to be one-to-one over only finitely
many points.  Thus a suitable `s_1` exists.  If the data are rational, density
of `Q` supplies `s_1 in Q`.

Rolle's theorem now gives `N-1` real critical points of `P_(s_1)`.  There can
be no others, and the choice of `s_1` makes all of them nondegenerate with
distinct critical values.

### Losing pairs one at a time

Let `a` be the leading coefficient of `H`, put `d=sign(a)`, and vary

\[
 E_q(W)=P_{s_1}(W)+d q,\qquad q\ge0.                   \tag{4}
\]

At `q=0` the polynomial has `N` simple real roots.  After replacing the
whole polynomial by its negative when `a<0`, its local minima are negative,
its local maxima are positive, and (4) moves the zero level downward through
the minima.  Their absolute critical values are positive and pairwise
distinct.  Each crossing removes exactly two real roots.  There are
`floor(N/2)` such minima when `N` is even and `(N-1)/2` when `N` is odd.
For large `q`, an even-degree polynomial has no real root and an odd-degree
polynomial has one.  Hence the successive open intervals between critical
values realize every count in (1).

Every such interval contains a rational `q`.  With rational `s_1`, equation
(4) therefore gives rational pencil witnesses.  This also gives a direct
algorithm: compute the discriminant of `P_(s_1)+t` as a polynomial in `t`,
isolate its real roots, and choose rational numbers between the ordered
roots on the ray `d t>0`.

## Weighted Keller maps

For an admissible real weighted seed with constant `c!=0`, the inverse
pencil on `C!=0` is

\[
 E_{A,B,C}(W)=H(W)-BCW+cAC^2.                           \tag{5}
\]

The change of target coordinates

\[
 (A,B,C)\longmapsto(s,t,C)=(BC,cAC^2,C)
\]

is a real analytic isomorphism on each slice `C=C_0!=0`.  At a simple root,
reconstruction gives

\[
 \gamma=-E'(W)/c,
 \qquad x=-cC/E'(W),
\]

and all remaining source coordinates are finite and uniquely determined.
Thus, on

\[
 C\ne0,
 \qquad \operatorname{Disc}_W(E_{A,B,C})\ne0,          \tag{6}
\]

the real source-fiber cardinality is exactly the real-root count of (5), and
the complex fiber is complete of cardinality `N`.

Apply the pencil theorem to the explicit rational seed `H_N` of the
[all-degree construction](ALL_DEGREE_RATIONAL_FIBERS.md).  At

\[
 (s_0,t_0)=\left(-\frac{G_N'(0)}{\lambda_N},0\right)
\]

the pencil is `G_N/lambda_N` and has `N` distinct integral roots.  Taking
`C=1`, `B=s`, and `A=t` because this family has `c=1`, proves:

> **Exact real-sheet spectrum.** For every `N>=3`, the explicit Keller map
> `F_N` has nonempty open real target chambers carrying exactly
> `N,N-2,...,N mod 2` inverse sheets.  These are all possible counts on the
> complete regular locus (6).  In particular, its minimum there is zero in
> even degree and one in odd degree.  Every count has a rational target
> witness.

For example, for the explicit quartic seed the checker takes

\[
 (B,C)=\left(-\frac{2999}{1000},1\right)
\]

and certifies that `A=0,1,2` give respectively `4,2,0` real preimages.  For
the explicit quintic seed it takes `(B,C)=(-1499/1000,1)` and `A=0,1,2`,
giving respectively `5,3,1` real preimages.

The folds here are folds of the univariate pencil projection, not critical
points of the Keller map: the latter has everywhere nonzero Jacobian.  At the
wall the two colliding root reconstructions run into the reconstruction pole.

## Wall crossing and the full chamber graph

At a smooth real point of the discriminant, with double root `r`, the local
equation has the form

\[
 u+v(W-r)^2+\text{higher terms},\qquad uv\ne0.
\]

A transverse crossing therefore changes the real-sheet count by exactly two.
The construction above exhibits the adjacency chain

\[
 N\;--\;N-2\;--\;\cdots\;--\;(N\bmod2).               \tag{7}
\]

For a fixed rational seed, the complete chamber adjacency graph is an exact
real-algebraic computation.  Form

\[
 D_H(s,t)=\operatorname{Disc}_W(H(W)-sW+t).
\]

A cylindrical algebraic decomposition with projection polynomials including
the leading coefficients and `Disc_t(D_H)` cuts the `s`-axis at all cusp,
bitangent, and higher-collision events.  On every complementary `s`-interval,
isolate and order the real roots of `D_H(s,t)`; Sturm counting for
`E_(s,t)` labels every band by its sheet number.  Gluing the bands across the
event fibers and taking the dual graph gives the complete chamber adjacency
graph.  All coefficients are rational for `H_N`, so algebraic real-root
isolation makes this algorithm exact and rational sample points can be chosen
in every two-dimensional cell.

The closed all-`N` result proved here is the spectrum and the canonical chain
(7).  A closed formula for every additional node-induced adjacency of the
particular `H_N` discriminant is a separate, seed-specific enumeration; it is
not needed for (1).

## Verification

Run

```bash
.venv/bin/python scripts/verify_real_fiber_spectrum.py
```

Add `--show-witnesses` to print the exact rational target pairs found in each
audited degree.

The checker constructs the rational seeds, chooses an exact rational generic
slope, proves squarefreeness of its critical-value discriminant, isolates the
ordered critical values, and certifies by exact Sturm counts one rational
target for every permitted sheet count in the audited degree range.
