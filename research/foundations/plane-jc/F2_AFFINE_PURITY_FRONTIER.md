# F2 affine-purity frontier

> **Status.**  This note compiles everything presently forced about the
> affine purity row in F2.  A hypothetical counterexample must contain a
> boundary divisor with transverse index `e>1` dominating an affine
> nonproperness curve, together with positive affine-sheet degree over the
> same curve.  None of the components in the certified `27/48` source graphs
> can serve: they are already log-etale, map to target infinity, or are
> contracted to a point.  Purity therefore forces at least one new source
> boundary component, raising the rigorous component floors to `28/49`.
> The target curve is parametric of degree at most `124`, and its generic
> ledger is finite, with `2<=e<=d-1<=9374`.  Conversely, the coarse purity
> axioms admit a signature at every currently possible geometric degree, so
> they neither improve the floors `d>=6/12` nor determine `(e,f)`, the curve,
> or its Chern correction.  Those require new global map data.

The graph increment, Bézout range, exact row bounds, and all-degree coarse
witness family are replayed by
[`verify_f2_affine_purity_frontier.py`](../scripts/verify_f2_affine_purity_frontier.py).

## 1. The universal affine purity row

Let `B` be the normalization of `k[P,Q]` in `k(x,y)` and

\[
 \pi:\operatorname{Spec}B\longrightarrow\mathbb A^2_{P,Q}
\]

the canonical finite flat cover of geometric degree `d`.  Its restriction to
the original affine source is étale.  Zariski--Nagata purity and the absence
of nontrivial connected finite étale covers of `A^2` imply that, when `d>1`,
some missing-boundary prime `E` has transverse ramification index

\[
 e(E/C)>1                                      \tag{1.1}
\]

over an irreducible component `C` of the affine nonproperness set.

For all primes over the generic point of `C`, write the boundary rows as
`(e_i,f_i)` and the affine residue degrees as `a_j`.  Finite flatness gives

\[
 \boxed{d=\sum_i e_if_i+\sum_j a_j},
 \qquad \sum_j a_j\ge1,                         \tag{1.2}
\]

and at least one `e_i>=2`.  Consequently every ramified row satisfies

\[
 2\le e_i\le d-1,
 \qquad 1\le f_i\le\left\lfloor\frac{d-1}{e_i}\right\rfloor. \tag{1.3}
\]

The affine contribution in (1.2) is the companion forced by a height-one
factor of a prime equation of `C` after substitution `(P,Q)`.

## 2. F2 numerical range

The terminal `A_6` packet lies over target infinity and has row `(e,f)=(1,6)`.
In the squarefree case it gives

\[
 d=6+\rho_T,qquad \rho_T\ge0,                  \tag{2.1}
\]

while the two same-target packets in the double case give

\[
 d=12+\rho_T,qquad \rho_T\ge0.                 \tag{2.2}
\]

The purity ledger (1.2) is over a different, affine target curve.  It is a
second equality with the same `d`, not an additional summand in (2.1) or
(2.2).  Thus purity does not change these lower bounds.

Bézout gives the finite upper bound

\[
 d\le\deg(P)\deg(Q)=75\cdot125=9375.             \tag{2.3}
\]

Hence

\[
\begin{array}{c|c}
\text{F2 row}&\text{geometric-degree interval}\\ \hline
\text{squarefree}&6\le d\le9375,\\
\text{double}&12\le d\le9375.
\end{array}                                      \tag{2.4}
\]

The Jelonek--Lasoń parametrization theorem further supplies a polynomial
parametrization of every `C` of degree at most

\[
 \max(75,125)-1=124.                             \tag{2.5}
\]

Equations (1.2)--(2.5) make the generic purity search finite once `d` and the
target curve are specified.

## 3. The new component is unavoidable

The certified squarefree/double graphs currently have `27/48` boundary
components.  Their complete generic classifications are:

1. terminal, interior-attachment, carrier, aligned-arm, spectator, and
   outgoing components are log-etale;
2. the extraction-root components are contracted to one target point and
   carry the cyclic root term `27`; and
3. the terminal principal components dominate the target divisor at infinity
   with transverse index one.

No component in either graph dominates an affine target curve with
`e>1`.  The purity divisor from (1.1) is therefore a new irreducible source
boundary component.  This proves the strengthened lower bounds

\[
 \boxed{
 N_{\rm source}^{\rm squarefree}\ge28,
 \qquad N_{\rm source}^{\rm double}\ge49.}       \tag{3.1}
\]

Additional resolution of this row may raise the bounds further.  Without
knowing its attachment point, no universal increase in the number of leaves
is asserted.

## 4. Exact underdetermination at the coarse level

For every integer `d>=3`, the formal signature

\[
 \text{boundary rows }((2,1)),qquad
 \text{affine contribution }d-2                 \tag{4.1}
\]

satisfies all generic numerical conditions (1.1)--(1.3).  Its boundary
normalization may be assigned one puncture, for which the coarse
Riemann--Hurwitz residual cost is

\[
 f+s-2=1+1-2=0.                                 \tag{4.2}
\]

Thus such a coarse ledger exists for every integer in both intervals (2.4).
This is not construction of a finite cover or Keller map.  It proves a sharp
negative statement: generic purity, finite-flat degree accounting, and
puncture arithmetic alone cannot select an F2 degree, ramification index,
residue degree, or contradiction.

## 5. Minimal missing data and next theorem target

To turn (1.2) into an actual F2 row, one must supply:

1. an equation or polynomial parametrization of the target curve `C`;
2. the source boundary component and proximity chain dominating it;
3. every boundary pair `(e_i,f_i)` and every affine factor of `g(P,Q)`;
4. punctures and special residue ramification on each boundary normalization;
5. full logarithmic matrices at singular or colliding attachment points; and
6. self-intersections and kernel-line data needed for the localized-Chern
   contribution.

The first datum cannot be recovered from the terminal `(5,2)` valuation:
that row is centered at target infinity, whereas `C` lies in the affine
target.  The subsequent
[`target-curve atlas`](F2_AFFINE_TARGET_CURVE_ATLAS.md) nevertheless reduces
it to 24 normalization charts `(deg p,deg q)=(3k,5k)`, `1<=k<=24`, and forces
a nonunit divided-difference collision/critical ideal.  Selecting a chart and
factoring the resulting implicit equation still requires global or lower
Laurent data.

<!-- status-consumer: PF2ATC1 9ab722c45c586b73 -->

The later
[`puncture-attachment theorem`](F2_AFFINE_PURITY_PUNCTURE_ATTACHMENT.md)
locates every such curve on the extracted `(5,2)` target divisor with contact
`k`.  For `k=1`, the third Belyi value `lambda=125/729` selects the special
carrier point, but the terminal neighborhood is already resolved and cannot
extract the affine divisor.  Its index must be computed at another source
boundary locus.

<!-- status-consumer: PF2PPA1 b24c4d80c2f8230e -->

Accordingly, the purity obligation is now finite and precisely typed, and
the source-component bound improves, but the row itself is not constructed.
This does not exclude `(75,125)` or prove `JC(2)`.

The later complement-monodromy theorem proves that one purity row is not
enough on every immersed, distinct-image `k=1` collision partition and on
the generic `A_2+3A_1` and `2A_2+2A_1` cusp strata.  Their affine
complements are `Z`, while the positive affine-sheet remainder fixes a
sheet, so transitivity requires a second ramified affine component and gives
conditional source floors `29/50`.  The first genuine escape is the
noncyclic `E_6+A_1` stratum, which admits a transitive degree-six
fixed-sheet action.  The unconditional all-chart floors remain `28/49`
because that severe cusp locus and `k=2,...,24` remain.

<!-- status-consumer: PF2K1M1 fafcbb3c2e6ceb2b -->

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_purity_frontier.py
```
