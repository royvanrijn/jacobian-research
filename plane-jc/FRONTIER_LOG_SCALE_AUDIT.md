# Log-scale audit of the `(72,108)` frontier case

## What the primary proof actually supplies

The proof of “Case `(8,28)`” in the pinned source
[`external/arxiv-2204.14178/Increasing_Lower_Bound_20_04_2022.tex`](external/arxiv-2204.14178/Increasing_Lower_Bound_20_04_2022.tex)
contains more local boundary information than the final two Newton polygons,
but less than a global proximity graph.  After interchanging `x,y`, it proves

\[
 \operatorname{Pred}_P(1,0)\in\{(1,-2),(1,-3)\}.
\]

Its transformation tree is:

```text
Pred = (1,-2)
  -> y -> y + lambda_1*x^-2
  -> successor (2,-7), or successor (1,-3)
                              |
Pred = (1,-3)                 |
  -> one linear factor -------+
  |    y -> y + lambda*x^-3
  |
  -> two linear factors
       (x^3*y-alpha_1)(x^3*y-alpha_2)
       choose one factor and translate by alpha_1*x^-3

  -> intermediate polygon a, b, or c
  -> common edge y*(x^4*y-alpha)^7
  -> y -> y + alpha*x^-4
  -> final monomial map x -> x^-1, y -> x^4*y
  -> the two published Laurent-polygon cases
```

The three translation orders are therefore

\[
 q=2,\qquad q=3,\qquad q=4,                            \tag{1}
\]

where `q=2` and `q=3` occur on alternative/continued predecessor branches,
while `q=4` is common after the three intermediate cases.

## Local toroidal meaning

At the SNC crossing

\[
 x=0,\qquad w=1/y=0,
\]

the translation `y -> y+lambda*x^-q` selects a branch with

\[
 w\sim \lambda^{-1}x^q.
\]

In parameters `(u,v)=(w,x)`, this is the positive scale

\[
 [w:x^q],
\]

whose equality ray is `(q,1)`.  Thus (1) supplies the local rays

\[
 (2,1),\qquad(3,1),\qquad(4,1).                        \tag{2}
\]

The functions `frontier_72_108_translation_records()` and
`laurent_translation_branch_certificate()` in
[`cas/log_boundary_compiler.py`](cas/log_boundary_compiler.py) record this
case-labelled skeleton and compile the regular fan of each visible equality
branch.  For order `q`, that fan has centers

\[
 (Y_\infty,X_0),\quad
 (Y_\infty,E_1),\quad\ldots,\quad
 (Y_\infty,E_{q-1}),
\]

and `E_q` is the equality divisor.  For example, `(3,1)` produces

\[
 (Y_\infty,X_0),\quad(Y_\infty,E_1),\quad
 (Y_\infty,E_2)
\]

and the equality divisor `E3`.  The regression suite compiles all three
regular branch fans of lengths `2,3,4`.

The elementary arithmetic labels are:

| `q` | equality ray | common base order | branch/residue degree | projection differents | branch conductor |
| ---: | --- | ---: | ---: | --- | ---: |
| 2 | `(2,1)` | 2 | 1 | `(1,0)` | 0 |
| 3 | `(3,1)` | 3 | 1 | `(2,0)` | 0 |
| 4 | `(4,1)` | 4 | 1 | `(3,0)` | 0 |

In particular, each visible Laurent branch `w-cx^q=0` is smooth.  Its local
semigroup conductor and residue degree provide no obstruction; any useful
different/conductor contradiction must involve the globally glued map and
its target pullbacks, not these isolated binomial germs.

## The branch fan is not the map resolution

This distinction is essential.  In local coordinates, the rational
translation itself has target homogeneous `y`-coordinates

\[
 [x^q+\lambda w:\;wx^q].
\]

For `lambda != 0`, with the regular parameter

\[
 t=x^q+\lambda w,
\]

and after multiplying the second generator by the unit `lambda`, its base
ideal satisfies the exact equality

\[
 (x^q+\lambda w,wx^q)
   =(t,x^q(t-x^q))
   =(t,x^{2q}).                                        \tag{3}
\]

Thus the branch fan has length `q`, while the monomial base chain of the
rational map in adapted coordinates has length `2q`:

| `q` | branch-fan length | map-base-ideal length |
| ---: | ---: | ---: |
| 2 | 2 | 4 |
| 3 | 3 | 6 |
| 4 | 4 | 8 |

The curve `t=0` is not one of the original boundary components, even though
it passes through the same initial crossing.  After reaching the equality
divisor, the remaining base point is residue-dependent.  Therefore the
branch fan alone cannot be substituted for the graph resolution.

The adapted ideal also supplies the missing local chart.  After the first
`q` blowups, `t=0` has separated from `Yinf` and meets `E_q` at a smooth
point; the residual ideal there is

\[
 (t_q,x^q).
\]

A nested cluster rooted at `E_q` with scale `[x^q:t_q]` contributes exactly
the remaining `q` blowups.  The implementation
`laurent_translation_graph_certificate()` compiles these two clusters and
therefore the complete isolated source graph of lengths `4,6,8`.  What it
does not by itself do is construct compatible target modifications and
pullbacks.  The later additive-composition audit identifies the source graph
across the alternative Newton cases.  The exact algebra check is
`laurent_translation_base_ideal_audit()`.

The same audit does compute the generic pullback orders on the exceptional
sequence.  The unchanged target coordinate `x` has order one on every
exceptional.  After removal of the common base factor, the target
`Yinf`-coordinate has orders

\[
\begin{array}{c|l}
q=2&(2,2,1,0),\\
q=3&(3,3,3,2,1,0),\\
q=4&(4,4,4,4,3,2,1,0).
\end{array}                                             \tag{4}
\]

Because the first components have positive order for both target boundary
coordinates, several are contracted to the target corner.  The entries in
(4) are pullback valuations, not ramification indices of a finite
divisor-to-divisor map; applying the formula `different=e-1` to all of them
would be invalid.

There is one further global simplification.  Every Laurent step fixes `x`,
so the translations compose additively.  Because the order-four coefficient
is nonzero in all three intermediate cases, their total translation has the
form

\[
 y\longmapsto y+x^{-4}F(x),\qquad F(0)\ne0.
\]

Putting `t=x^4+wF(x)` gives the same base ideal `(t,x^8)` for cases `a,b,c`,
regardless of whether the lower orders are `(2)`, `(3)`, or `(2,3)`.
Consequently all cases share the same eight-blowup translation graph.  The
orders two and three survive as coefficient and equality-residue labels, not
as additional graph branches.  The exact identity is checked by
`laurent_translation_composition_audit()`.

The final monomial step is also absorbed by the correct completion.  The map

\[
 (x,y)\longmapsto(x^{-1},x^4y)
\]

is an involutive automorphism of `K[x,x^-1,y]`.  It is the transition between
the two base charts of the line bundle `O(-4)`, so on the Hirzebruch
completion `F_4` it adds no blowups.  Its valuation matrix is

\[
\begin{pmatrix}-1&0\\4&1\end{pmatrix},
\qquad \det=-1,\qquad M^2=I.
\]

The compiler now supports `GmA1_Fn`; the common order-four graph is emitted
on `GmA1_F4`.

The transition swaps the base divisors:
pre-transition `X0` becomes post-transition `Xinf`, while pre-transition
`Xinf` becomes post-transition `X0`.  The affine chart therefore fills
pre-transition `Xinf` (equivalently post-transition `X0`), not the component
through which the eight blowups run.  The ten post-transition boundary
primes are

\[
 X_\infty,\ Y_\infty,\ E_1,\ldots,E_8.
\]

Their `10 x 10` class matrix has determinant one and Smith diagonal
`(1,...,1)`.  Hence the filled open has unit rank zero and trivial Picard
group, exactly as an affine-plane source must.  Its intersection matrix has
determinant `-1`, inertia `(1,9,0)`, and the adjunction reconstruction gives
`K_X^2=0`, so `K_X^2+rho(X)=10`.  The report emits both matrices and the
intrinsic `A2` audit together with the divisor correspondence.  This
completes the source-boundary lattice translation
for the common Laurent graph; it does not yet prove that this source graph is
the complete boundary data of a hypothetical Keller map.

The symbolic residue tree itself is now encoded: `lambda_1` labels the
order-two step, `lambda` the one-factor order-three step,
`alpha_1 != alpha_2` the split factors with `alpha_1` selected, and nonzero
`alpha` the common order-four step.  The unselected factor can be continued
directly.  Put `beta=alpha_2-alpha_1`, which is nonzero.  After the selected
order-three translation it is `x^3y-beta`.  After the common translation
`y -> y+alpha*x^-4`, it becomes

\[
 x^{-1}(x^4y+\alpha-\beta x).                         \tag{5}
\]

At the order-four center `x=x^4y=0`, the numerator in (5) is the nonzero
value `alpha`; hence this factor is a unit and adds no infinitely-near base
point.  Under `X=x^-1`, `Y=x^4y`, its equation is
`X(Y+alpha)-beta`, which is `-beta != 0` on the filled divisor `X=0`.
The symbolic check is `frontier_72_108_unselected_factor_audit()`.

The primary leading edge also determines the target pole orders on every
exceptional in the common graph.  With `w=1/y`, put

\[
 R=x^{24}w^{-8}(x^4+\alpha w).
\]

The forced edge has `P`-leading term `R^2` and `Q`-leading term `R^3`.
The forced low vertices contribute `x^-1` to `P` and `x^2/w` together with
a constant to `Q`.  On `E1,...,E4`, the orders of `R` are
`(17,10,3,-4)`.  After writing `x^4+alpha*w=x^4*tau`, its orders on the
four residual exceptionals are `(-3,-2,-1,0)`.  Taking the exact competing
minima gives:

| divisor | `E1` | `E2` | `E3` | `E4` | `E5` | `E6` | `E7` | `E8` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| pole of `P` | 1 | 1 | 1 | 8 | 6 | 4 | 2 | 1 |
| pole of `Q` | 0 | 0 | 1 | 12 | 9 | 6 | 3 | 2 |
| pullback of target infinity | 1 | 1 | 1 | 12 | 9 | 6 | 3 | 2 |

Whenever a displayed minimum is negative, its competing orders are distinct,
so coefficient cancellation cannot change the table.  The executable check
is `frontier_72_108_exceptional_pole_audit()`.  This supplies the target-line
pullback on all eight exceptionals, but not on the original completion
components and not an exhaustiveness proof excluding other global clusters.

The two original-boundary entries are forced as well.  Generic
pre-transition `X0`, which is post-transition `Xinf`, has `(P,Q)` pole
orders `(1,0)`; generic `Yinf` has orders `(16,24)`.  In post-transition
boundary order `(Xinf,Yinf,E1,...,E8)`, the full common-graph target pole
vector is

\[
 p=(1,24,1,1,1,12,9,6,3,2).                         \tag{6}
\]

The intrinsic boundary gate gives

\[
 Qp=(0,12,0,0,11,10,0,0,2,1),\qquad p^TQp=427.
\]

The pullback is nef and its ordinary and logarithmic ramification vectors
are effective.  But every entry of `p` is positive, so no boundary prime
satisfies the dicritical criterion `p_i=0` and `(Qp)_i>0`.  Consequently the
common eight-blowup graph is **provably not** the complete resolution of a
hypothetical nonproper Keller map.  At least one additional dicritical
boundary cluster is forced.  The exact check is
`frontier_72_108_common_graph_pole_audit()`.

The first exact Laurent block now determines what actually happens near
`E3`.  At the crossing `E3 intersect E4`, use

\[
 s=Y,\qquad t={1\over XY},\qquad T={s\over t},\qquad z={1\over s}.
\]

Thus `s=0` and `t=0` are local equations for `E3` and `E4`.  Write

\[
 A(T)=T\,U(T),\qquad D(T)=T^2V(T),
\]

where `deg(U)=7` and `deg(V)=10`.  After multiplying the target sections
`[1:P:Q]` by `s*t^12`, their lowest pieces at the crossing are

\[
 [s t^{12},\ t^{11}U(s/t),\ t^{10}V(s/t)].            \tag{7}
\]

The crossing has common base order ten.  Its blowup creates an exceptional
`B0` of pole `1+12-10=3`.  The weighted-Wronskian equation

\[
 2AD'-3A'D=T^2
\]

reduces exactly to

\[
 UV+T(2UV'-3U'V)=1.                                   \tag{8}
\]

Consequently `V(0)=1`, `gcd(U,V)=1`, and every root of `V` is nonzero and
simple.  There are ten distinct roots.  Each is a free basepoint on `B0`
with local ideal `(epsilon,v)`; blowing it up creates a pole-two exceptional
of target-line degree one.  The extended pole vector is

\[
 p_{\rm src}=(p,3,\underbrace{2,\ldots,2}_{10}),
\qquad p_{\rm src}^{\,2}=317.                          \tag{9}
\]

It still has no dicritical.  This cluster is checked by
`frontier_72_108_forced_first_block_cluster_audit()`.

This also disposes of the previously apparent free-`E3` alternative.  On
the smooth torus of `E3`, the forced low monomial `X` is a unit after its
pole is removed, so it cannot be a source basepoint.  If one ignores
residues and tests all ten smooth boundary components and all nine crossings
with a new zero pole, a smooth `E3` point is still the unique numerical
candidate:

\[
\begin{aligned}
p'&=(1,24,1,1,1,12,9,6,3,2,0),\\
Q'p'&=(0,12,0,0,10,10,0,0,2,1,1),\\
(p')^TQ'p'&=426.
\end{aligned}
\]

That counterfactual is retained in
`frontier_72_108_minimal_dicritical_extension_audit()` because it shows that
the intrinsic lattice gate alone cannot recover coefficient residues.

The terminal nonvertical polygon, Case 2, supplies a dicritical cluster
directly.  At the smooth point `(X,z)=(0,0)` of `Yinf`, with `z=1/Y`, the
sections `[1:P:Q]z^24` have common order twelve.  After the first blowup
`X=e*u, z=e`, a monomial `X^iY^j` in Laurent band `k=2i-j` has
second-stage total order `12+k`.  Case 2 has only nonnegative bands, so the
second common order is again twelve.  The two new poles are

\[
 24-12=12,\qquad 12-12=0.
\]

On the zero-pole exceptional, with `u=e*r`, the residue map is

\[
 r\longmapsto [1:C(r):G(r)],
\qquad \deg C=8,\quad\deg G=12.                       \tag{10}
\]

It has target-line degree twelve.  Combining this target cluster with (9)
gives a 23-component intrinsic affine-plane boundary, remaining
self-intersection

\[
 427-10^2-10\cdot1^2-12^2-12^2=29,
\]

and one reconstructed dicritical, of degree twelve.  Every intrinsic
boundary, nefness, adjunction, and boundary divisor-class test passes.
Because the transformed pair has bracket `X^2`, the boundary-supported
`K+3H` vector used in this test is not yet the actual ramification divisor;
the correction is made below.
The exact support enumeration and matrix calculation are
`frontier_72_108_case2_boundary_package_audit()`.

The remaining corner is also forced by the final Newton edge.  Put

\[
 z=1/Y,\qquad s=X/z=XY.
\]

On the parallel upper edge write

\[
 P_{\rm edge}=Y^8A(s),\qquad Q_{\rm edge}=Y^{12}C(s).
\]

A direct bracket calculation gives

\[
 [P_{\rm edge},Q_{\rm edge}]_{\rm top}
 =4Y^{20}(3A'C-2AC').                                \tag{11}
\]

The right side cannot occur in `[P,Q]=X^2`, so `3A'C-2AC'=0`.  Unique
factorization in characteristic zero then gives

\[
 A=a\,r^2,\qquad C=c\,r^3,\qquad \deg r=4.            \tag{12}
\]

For terminal Case 2 the vertical vertices `(0,8),(0,12)` are absent.  The
support therefore forces `r=s^4`, a single root of multiplicity four at
`s=0`; this recovers the two-blowup package (10).  If one uses only the
Case-1 edge, both the constant and leading coefficients of `r` are nonzero.
Its roots are nonzero, and their multiplicities have one of the five
edge-only partitions

\[
 (4),\ (3,1),\ (2,2),\ (2,1,1),\ (1,1,1,1).          \tag{13}
\]

For the homogeneous edge model, the first blowup of the smooth `Yinf` point
has base order twelve and exceptional `D0` of pole twelve.  At a root of
multiplicity `m`, with local parameter `u`, the homogeneous ideal is

\[
 (e^{12},\,e^4u^{2m},\,u^{3m}).                       \tag{14}
\]

Let `g=gcd(m,4)`.  The three valuation orders in (14) agree on the primitive
ray

\[
 (m/g,4/g).
\]

Inserting that ray by regular star subdivisions principalizes the ideal.
The equality divisor is dicritical with residue

\[
 [1:t^{2g}:t^{3g}],                                  \tag{15}
\]

while every other new toric divisor is contracted.  The complete local data
are:

| `m` | equality ray | inserted rays in creation order | pole orders | residue exponents | hyperplane degree |
| ---: | --- | --- | --- | --- | ---: |
| 1 | `(1,4)` | `(1,1),(1,2),(1,3),(1,4)` | `9,6,3,0` | `(2,3)` | 3 |
| 2 | `(1,2)` | `(1,1),(1,2)` | `6,0` | `(4,6)` | 6 |
| 3 | `(3,4)` | `(1,1),(1,2),(2,3),(3,4)` | `3,0,0,0` | `(2,3)` | 3 |
| 4 | `(1,1)` | `(1,1)` | `0` | `(8,12)` | 12 |

The multiplicity-three row is worth stressing: two of its pole-zero
divisors are contracted.  Treating the row as two successive free blowups
misidentifies the dicritical and gives the wrong global degree.

Every root branch `u=0` is smooth, so its source conductor is zero.  The
image of (15) is the cusp with semigroup `<2,3>` and conductor two.  The map
from the equality divisor to the cusp normalization has degree `g` and tame
different exponent `g-1` at the marked point; its two displayed coordinate
projections have different exponents `2g-1` and `3g-1`.

Attaching these homogeneous fans to the common graph and the forced
ten-point first block gives:

| partition | boundary components | dicritical hyperplane degrees | dicritical `K+3H` coefficients | remaining self-intersection |
| --- | ---: | --- | --- | ---: |
| `(4)` | 23 | `12` | `0` | 29 |
| `(3,1)` | 30 | `3,3` | `3,3` | 29 |
| `(2,2)` | 26 | `6,6` | `1,1` | 29 |
| `(2,1,1)` | 32 | `6,3,3` | `1,3,3` | 29 |
| `(1,1,1,1)` | 38 | `3,3,3,3` | `3,3,3,3` | 29 |

For every edge-only row the full boundary-class matrix is square and
unimodular; the
intersection matrix passes the affine-plane adjunction and Noether gates;
the target pullback is nef; and the ordinary and logarithmic ramification
class representatives are effective.  These representatives agree with
actual ramification only before a nonconstant-Jacobian coordinate
transformation is made.  The machine-readable matrices, pole vectors,
differents, and conductors are emitted by
`frontier_72_108_plane_return_partition_audits()`.  These five rows classify
what the leading edge and the intrinsic gate alone permit; for a multiple
root they do not, by themselves, control transverse terms.

The primary split-factor formula supplies the missing stronger information.
After the selected order-three and common order-four shifts, it gives

\[
\begin{aligned}
P_{\rm edge}&=(Y+\alpha)^8
 (X(Y+\alpha)-\beta)^8,\\
Q_{\rm edge}&=(Y+\alpha)^{12}
 (X(Y+\alpha)-\beta)^{12},
\end{aligned}
\qquad \beta=\alpha_2-\alpha_1\ne0.
\]

Therefore `r=(s-beta)^4`: the source selects partition `(4)`, not an
arbitrary member of (13).  Transverse terms are fixed by the symmetry between
the two legal choices of order-three factor.  Repeating the construction with
`alpha_2` selected gives a second final plane chart related by

\[
 Y_2=Y_1-{\beta\over X}+\delta,\qquad
 {1\over Y_2}={Xz\over X-\beta z+\delta Xz}.           \tag{16}
\]

On the first blowup `X=e u,z=e`, the exact adapted parameter is

\[
 v=u-\beta+\delta e u.
\]

Blowing up `e=v=0` gives `Y_2=r/u` on the new exceptional.  Hence its
residue is the restriction of the alternate polynomial pair to `X=0`, of
degrees eight and twelve.  The second base order is twelve, its pole is zero,
and its target-line degree is twelve.  No assumption that transverse
coefficients vanish is used.  This source-selection check is
`frontier_72_108_case1_boundary_package_audit()`.

## The transformed-Poisson correction

There is an important distinction between the intrinsic divisor-class audit
and the actual ramification of the final Laurent pair.  Proposition Case
`(8,28)` ends with

\[
 [P,Q]_{X,Y}=X^2,
\]

not a nonzero constant.  If `c=K+3H` is the boundary-supported class vector
emitted by the generic Keller gate and

\[
 \operatorname{div}(X)=L_X+\sum_i\nu_iD_i,
\]

then the actual ramification divisor is

\[
 R_{\rm act}=K+3H+\operatorname{div}(X^2)
 =2L_X+\sum_i(c_i+2\nu_i)D_i.                        \tag{17}
\]

For the two source-selected 23-component packages, the common boundary
ordering is

```text
X0,Yinf,E1,...,E8,B0,B1,...,B10,(two return divisors).
```

The common part of the `X`-valuation vector is

```text
-1,0,-1,-1,-1,-1,-1,-1,-1,-1,-2,-2,...,-2,
```

with ten `-2` entries after `B0`.  The two return entries are `(1,1)` in
Case 1 and `(1,2)` in Case 2.  Applying (17) gives the common actual
boundary coefficients

```text
3,70,2,1,0,32,24,16,8,6,3,1,...,1,
```

followed by `(37,2)` in Case 1 and `(37,4)` in Case 2.  The principal-divisor
identity `L_X.D=-Q*nu` says that `L_X` meets only `D0` in Case 1 and only the
dicritical `T2` in Case 2, once in either case.

The local charts independently verify the same coefficients:

| case | final dicritical chart | pulled-back Jacobian | coefficient of the dicritical in `R_act` | generic normal index | boundary part of `R_act.D` | affine part | total |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |
| 1 | `X=e*u`, `u=(beta+e*r)/(1+delta*e)`, `Y=e^-1` | `e^2*u^2/(1+delta*e)` | 2 | 3 | 35 | 0 | 35 |
| 2 | `X=e^2*r`, `Y=e^-1` | `e^4*r^2` | 4 | 5 | 33 | 2 | 35 |

Thus the old intersection number `35` was correct, but interpreting the
zero entry of the boundary-supported `K+3H` vector as zero ramification was
not.  The corrected ledgers are checked by
`frontier_72_108_poisson_ramification_audits()`.

## Residue degree and the lower-jet reduction

Let `delta_cov` be the degree from the dicritical `P^1` to the normalization
of its image.  Since the exact residue has coordinate degrees `(8,12)`,
polynomial Lüroth theory gives

\[
 \delta_{\rm cov}\mid\gcd(8,12)=4.
\]

Consequently the only cover/image-degree splits and total Hurwitz differents
are

| `delta_cov` | image degree | total normalization different |
| ---: | ---: | ---: |
| 1 | 12 | 0 |
| 2 | 6 | 2 |
| 4 | 3 | 6 |

The homogeneous quartic initial form realizes `delta_cov=4`, but it does not
fix the lower coefficients of the exact alternate-chart residue in Case 1.
It would therefore be incorrect to declare degree four forced there.  With
the actual normal indices from (17), the contribution `e*delta_cov` to the
degree-29 valuation identity and the unaccounted remainder are:

| case | normal index `e` | possible contributions | remaining degrees |
| --- | ---: | --- | --- |
| 1 | 3 | `3,6,12` | `26,23,17` |
| 2 | 5 | `5` | `24` |

There is a concrete data boundary on the Case-1 decomposition test.  The
archived descent is a complete necessary *upper truncation*: it contains
`P` only through `z^-5` and `Q` only through `z^-4`, whereas the full
Newton polygons continue through `z^-8` and `z^-12`.  The alternate residue
`[1:P_2(0,r/u):Q_2(0,r/u)]` depends on those omitted lower bands.  Thus its
exact coefficient vector cannot be reconstructed from the archive, and a
Case-1 right-component sieve must not insert guessed zero coefficients.
That sieve becomes executable only after the missing bands are derived, or
after a lemma proves that its remainders depend solely on the archived
truncation.

The Case-2 row is now strictly smaller than its a priori degree split.
If its normalization cover had degree `delta=2` or `4`, then the polynomial
parametrization would factor through a monic polynomial right component

```text
delta=2: h=t^2+a*t
delta=4: h=t^4+a*t^3+b*t^2+c*t.
```

The top coefficients of `C` reconstruct `a`, and then `b,c`, uniquely.
After solving `(J3),(J2)` and only the determined part of `(J1)`, division
of `C,G` by powers of this general `h` gives respectively 9 and 12 distinct
positive-degree remainder equations.  Both exact ideals are the unit ideal.
No `(J1)` compatibility equation and no `(J0)` equation is used.  Thus

\[
 \boxed{\delta_{\rm cov}=1,\qquad \deg(\operatorname{image})=12}
\]

in Case 2.  This is a general polynomial-decomposition test, not an
evenness or parity assumption.  The pinned checker is
[`cas/audit_case2_residue_strata.py`](cas/audit_case2_residue_strata.py).

The remaining row is empty too.  Degree twelve forces `G_12 != 0`.
After adjoining `w*G_12-1`, the seven residual `(J1)` compatibility
equations--each an eight-term cubic in the three remaining parameters--
generate the unit ideal over the exact degree-35 field.  This uses no
coefficient of `(J0)`.  The pinned checker is
[`cas/case2_infinity_resolution.py`](cas/case2_infinity_resolution.py).
Before compatibility, its generic infinity branch would have characteristic
orders `(4,13)` and a seven-ray toric resolution; the endpoint certificate
shows that this target graph has no compatible Case-2 input.

For Case 2, the two bottom Laurent equations now admit a useful exact
reduction.  Put

\[
 H=\gcd(C',G'),\qquad C'=Hc,\qquad G'=Hg,\qquad \gcd(c,g)=1.
\]

Then `(J0)` is equivalent to

\[
 B=Kc,\qquad F=Kg.
\]

The support bounds imply `deg(K)<=deg(H)+1`; because `B(0)=F(0)=0`,
coprimality further gives `K=0` or `t` divides `K`.  Substitution in `(J1)`
leaves the single equation

\[
 2H(Ag-cE)+K^2(cg'-c'g)=0.                           \tag{18}
\]

The first coefficients sharpen this once more.  Triangularly, `(J4)--(J0)`
give

```text
D_2=1, B_1=2*E_2, F_1=F_2=0, C_1=0, G_1=G_2=0.
```

In particular `t` divides `H` in every hypothetical solution, so the
degree-zero gcd stratum is impossible.  Before the decomposition audit, a
cover of degree `1,2,4` forced respectively `deg(H)>=1,1,3`; the last two
cover strata are now empty.  The endpoint-compatibility audit also removes
the last degree-one cover.  The exact gcd stratification below is therefore
a geometric analysis of the pre-compatibility family, not a list of live
Case-2 rows.  The exact symbolic check is
`frontier_72_108_case2_lower_jet_audit()`.

The smallest surviving stratum is already rigid at the origin.  If
`deg(H)=1`, then `H` is proportional to `t`, coprimality forces
`ord(C')=1`, and a leading-order comparison in `(J2)` leaves only

```text
(ord(K),ord(g))=(1,2) or (2,1).
```

The second pair makes the `A*G'` term uniquely lowest in `(J1)` and is
impossible.  Therefore every degree-one-gcd solution would have the forced
pattern

```text
ord(K)=ord(B)=1, ord(E)=2,
ord(g)=2, ord(G')=ord(F)=3.
```

This removes another local branch before any number-field coefficient
calculation.

It also fixes part of the boundary linear series.  On this row
`[1:C:G]` is a basepoint-free `g^2_12`.  At the forced affine point its
vanishing sequence is `(0,2,4)`, of Pluecker weight `3`; at infinity the
coordinate degrees `(0,8,12)` give `(0,4,12)`, of weight `13`.  The total
weight is `3(12-2)=30`, hence exactly `14` remains:

```text
C'*G''-C''*G' = t^3*W_14(t),  deg(W_14)=14.
```

This turns the next birational-row test into a finite allocation of fourteen
flex/singularity weights on the resolved target graph.

The Wronskian is not independent of (18).  Exact differentiation gives

\[
 C'G''-C''G'=H^2(cg'-c'g),
\]

and multiplying (18) by `H^2` gives

\[
 K^2(C'G''-C''G')=-2H^3(Ag-cE).                    \tag{19}
\]

Since the affine Wronskian has degree `17`, an exact gcd degree
`h=1,...,7` leaves relative-Wronskian degree
`17-2h=15,13,11,9,7,5,3`.  Thus the sole birational cover is now seven
finite gcd/Wronskian rows; on the first row the forced extra factor at the
origin reduces `15` to the displayed `14` free weights.

The maximal row is empty.  When `h=7`, the degree-seven polynomial `C'`
divides `G'`.  Divide `G'` by `C'` after the exact `(J3),(J2)` solution and
the determined part of `(J1)`.  Of the seven remainder coefficients, those
of degrees `0,1,2`, together with the coefficient of `t^19` in
`(J0)=BG'-C'F`, generate the unit ideal over the exact degree-35 first-block
field.  Their term counts are `(155,155,155,9)` and their parameter degrees
are `(13,13,13,4)`.  No residual `(J1)` compatibility equation is used.
The Singular input is pinned by SHA-256
`1ac0b4db7ddd0b50fcbef6c93d49c28f7f80cb4133a73b5f2158af6c78f3b069`;
run

```bash
make verify-plane-case2-maximal-gcd
```

Hence only the six exact gcd degrees `1,...,6` remain.  Their
relative-Wronskian degrees are `15,13,11,9,7,5`.

For context, arbitrary numerical positive-pole preparatory blowups on the
unextended common graph give four connected witnesses:

| first center | first pole | resulting degree | dicritical degree | dicritical canonical coefficient |
| --- | ---: | ---: | ---: | ---: |
| `Yinf` | 12 | 139 | 12 | 0 |
| `E4` | 2 | 323 | 2 | 0 |
| `E7` | 1 | 422 | 1 | 3 |
| `E8` | 1 | 425 | 1 | 4 |

They are not asserted to be source basepoints; they delimit what intersection
theory alone can exclude.  They are checked by
`frontier_72_108_two_step_dicritical_witnesses()`.

The complete machine-readable local package is emitted by

```bash
.venv/bin/python plane-jc/cas/log_boundary_compiler.py --frontier-72-108
```

It contains both matrices for every branch fan and source graph, the
composition collapse, the `F_4` transition, the affine-plane fill, the
symbolic residue tree, the unselected-factor continuation, the common-graph
pole audit, the forced cluster (7)--(9), the Case-2 package (10), and all
five edge-only models (11)--(15), together with the source-selected Case-1
package (16), the transformed-Poisson correction (17), all three residue
degree splits, and the lower-jet reduction (18).

## Source-selected conclusion

The source proof works through successive Laurent automorphisms and calls its
last monomial substitution a morphism rather than a polynomial automorphism.
On the Laurent ring it is the involution above, and the source completion is
fixed by the `F_4` graph and affine-plane fill.  The first block now supplies
the missing `E3 intersect E4` source cluster and rules out the numerical
smooth-`E3` option.  The final Poisson-square edge supplies every remaining
dicritical chain.  Case 2 selects `(4)` at the zero root.  The split-factor
symmetry selects `(4)` at the nonzero root in Case 1 and supplies its exact
transverse parameter.  Both source-selected terminal packages therefore have
23 boundary components, remaining self-intersection 29, and one
degree-twelve dicritical.

None of these identifications is determined by the final corners
`(8,28)` and `(11/4,7)`.  The common order-four graph is now justified by the
composition identity, but treating every intermediate factor as the same
residue point—or as a separate global boundary component—would still be an
unsupported assertion.  The implementation emits the source-selected
architecture and the larger edge-only comparison family.  It rejects only
the legacy aggregate `NewtonBoundaryCertificate`, because that generic IR
cannot serialize the nonmonomial first-block root cluster; this is an IR
limitation, not missing boundary geometry.

## Consequence

The earlier “find either a free `E3` residue or some longer chain” gap is
closed: the first Wronskian forces the longer cluster (7)--(9), and no smooth
`E3` basepoint exists.  The nonvertical terminal polygon is compiled through
its dicritical, and split-factor symmetry gives the same numerical package
for the vertical terminal polygon.  The other four edge-only partitions also
pass the intrinsic gates.  The log-boundary lead has therefore become a
completed prefilter for this frontier: it does not exclude `(72,108)`;
instead it proves that the current intrinsic package is realizability-neutral
and that a stronger invariant than the present lattice/ramification gate is
needed before coefficient elimination can be replaced.

The next target is now finite and correctly normalized.  The dicritical has
hyperplane degree twelve and total ramification intersection `35`, but its
actual normal index is three in Case 1 and five in Case 2.  Case 1 retains
cover degrees `1,2,4`, and its decomposition sieve awaits the omitted lower
bands.  Case 2 has no compatible degree-twelve endpoint and is closed before
`(J0)`.  The immediate one-pair task is therefore solely Case 1: derive the
missing lower bands or prove the truncation lemma, then run the three
composition strata.  The Case-2 equation (18) and its gcd rows remain useful
only for compressing the already established exclusion geometrically.
