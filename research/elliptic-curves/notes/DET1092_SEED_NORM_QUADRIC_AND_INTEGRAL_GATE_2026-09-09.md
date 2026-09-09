# The seed norm equations: one rational quadric, two integral lattices

## Result and unresolved endpoint

**New explicit application of elementary norm/trace algebra.** The two
unresolved seed-relative principality equations have the same rational
geometric model: a smooth rational quadric constructed from the cubic field
alone. Explicit maps in both directions and two exact integrality matrices
are independently checked. All missing rational chart points are accounted
for; none of those exceptional points solves either integral equation.

**New explanatory deduction.** The relative ideal class is not encoded in
the rational isomorphism type of these norm surfaces. It survives in their
**integral lattices**. Dropping that lattice makes every case soluble, whether
or not its ideal is principal. A uniform equation-only formula also makes
the corresponding norm quadric rational on302 and all eight frozen controls.
Rational solubility of this smaller variety is therefore not a seed predictor.

Neither principality equation is solved. No new unit, ideal-class relation,
elliptic point, rank gain or prospective302 selector is claimed. The
[previous arithmetic fork](DET1092_SEED_VIRTUAL_UNIT_AND_IDEAL_PARITY_2026-09-09.md)
remains the endpoint; this is an explicit reduction preserving that fork,
not a replacement success criterion.

## 1. Audit of the two equations

**Previously verified input.** For the302 cubic field `K`, two integral ideals
`R_j` remain after the seven Artin characters and all compatible generic
corrections. Their exact matrices, elements and norms satisfy

\[
 R_j^2=(\beta_j),\qquad N(\beta_j)=N_j^2,\qquad N_j=N(R_j)>0.
\]

For an integral basis `e_{j,0},e_{j,1},e_{j,2}` of `R_j`, the unresolved test is

\[
 F_j(x)=N\left(\sum_i x_i e_{j,i}\right)/N_j=1,
 \qquad x\in\mathbf Z^3.                                \tag{1}
\]

Existence of an integral point is equivalent to principality: inclusion of
its principal ideal in `R_j` and equal norms force equality. Conversely,
change a principal generator's sign if necessary to give norm `N_j`.

The existing rational witness is `a_j=N_j/beta_j`, with `N(a_j)=N_j`.
Both equations are also soluble over every `Z_p`; that was proved in the
previous note, not inferred from a finite prime sample.

**New explicit isomorphism.** Define the norm-one torus
`T={u∈K: N(u)=1}`. Multiplication by `a_j` gives a rational linear isomorphism

\[
 T\longrightarrow\{F_j=1\},\qquad
 u\longmapsto \operatorname{coord}_{R_j}(a_j u).          \tag{2}
\]

This works for every square ideal satisfying the displayed norm identity,
not just these two examples. Even the projective norm cubic closures become
isomorphic over `Q` under the same linear change. Their integral models are
not identified. Thus their rational geometry alone cannot tell whether a
relative ideal is principal.

## 2. A common quadratic model, with explicit inverse

**Established algebra, explicitly applied.** Write the cubic in depressed
form

\[
 K=\mathbf Q(\eta),\qquad \eta^3+A\eta+B=0,
 \qquad \Delta=-4A^3-27B^2\ne0.
\]

Take the trace-zero basis `e1=eta`, `e2=eta²+2A/3`, and write
`v=v1 e1+v2 e2`. Since `Tr(v)=0`,

\[
 N(v+s)=s^3-\tfrac12\operatorname{Tr}(v^2)s+N(v).
\]

The smooth projective quadric

\[
 \boxed{Q_K:\quad
 A v_1^2+3Bv_1v_2-\frac{A^2}{3}v_2^2+a^2+ab+b^2=0}     \tag{3}
\]

therefore satisfies `N(v+a)=N(v+b)`, away from the factor `a-b`. The forward
map is

\[
 u=(v+a)/(v+b).                                         \tag{4}
\]

For a nonidentity norm-one element, the inverse is explicit:

\[
 z=\frac1{u-1},\quad b=\operatorname{Tr}(z)/3,\quad
 a=b+1,\quad v=z-b.                                    \tag{5}
\]

These identities prove birationality, not just a dominant parametrization.
For the302 field, `u≠1` makes `u-1` invertible. Conversely `N(v+b)` cannot
vanish at a rational point of (3): a zero norm in the field would give
`v+b=0`, then trace zero gives `b=v=0`, and (3) forces `a=0` as well.
The locus `a=b`, where defined, maps only to the identity.

**Verified discriminant identity.** The Gram determinant of (3) is

\[
 \det Q_K=\Delta/16.
\]

Its discriminant squareclass is the quadratic resolvent squareclass of the
cubic, not a newly discovered ideal or point class. The302 value is nonsquare;
the quadric has real signature `(2,2)`.

### A uniform rational point and parametrization

**New explicit application.** A rational point is

\[
 p=(0,-3,2A,-A).
\]

For `A B≠0`, (4) sends it to `u0=1+(A/B)eta`, whose norm is1.
Set

\[
 q=A r^2+s^2+sh+h^2,\qquad \ell=-9Br+3As.
\]

Projection from `p` gives the homogeneous quadratic map

\[
 \boxed{[r:s:h]\longmapsto
 [-\ell r:-3q:2Aq-\ell s:-Aq-\ell h].}                  \tag{6}
\]

Direct substitution verifies (3). The standard line-projection argument is
also checked algebraically below; no conic solver or rational-point search
is required. Rationality of a quadric with a rational point is established
geometry; see [Skorobogatov's algebraic-geometry notes](https://www.ma.ic.ac.uk/~anskor/AG.PDF).
The contribution here is the explicit formula and its faithful integral
transport, not a new general rationality theorem for tori.

For consistency with the old integral cubic basis, the two ideal certificates
use the equivalent trace-zero basis
`theta-Tr(theta)/3`, `theta²-Tr(theta²)/3`, and projection source
`theta³/N(theta)`. Those were frozen before the relative ideals were read.
The short-cubic formula (3)--(6) is a separately verified simplification;
it does not replace or retune the original certificate's chart.

## 3. The integral information that must not be discarded

**Verified construction.** Let `B_j` be the matrix whose columns are the
ideal basis in cubic power coordinates, and `M(a_j)` the multiplication
matrix of `a_j`. The exact transport matrix is

\[
 L_j=B_j^{-1}M(a_j),\qquad x=L_j\operatorname{coord}(u).
\]

The two `L_j` are recorded in full. Clear their denominators as
`L_j=C_j/d_j`, with an integral matrix `C_j` and positive integer `d_j`.
To make (4) explicit, form the multiplication matrix `M(v+b)` and put

\[
 D=N(v+b),\qquad
 (N_0,N_1,N_2)^t=\operatorname{adj}(M(v+b))
                   \operatorname{coord}(v+a).
\]

All four expressions are homogeneous cubics in the quadric coordinates.
Multiply all four by a common positive integer to clear their rational
coefficient denominators; their quotients are unchanged. For an integral
projective representative, (1) is then exactly the divisibility condition

\[
 \boxed{d_jD\mid(C_jN)_i\quad(i=0,1,2).}                \tag{7}
\]

After substituting (6), clear the fixed rational coefficient denominators
before applying (7). Homogeneous scaling does not change the three quotients.
This is an explicit integral gate, not a claim that checking some fixed
finite set of congruences suffices or that a small parameter box is complete.

**Independent boundary audit.** For the certificate's projection point `p`,
the nonsquare Gram determinant excludes a rational line on the quadric.
A second rational point in its tangent plane would span a totally isotropic
plane with `p`, forcing square determinant. Thus the only rational tangent
exception is `p` itself. The inverse (5), line projection and its inverse
therefore cover every nonidentity norm-one point except its explicitly
listed projection source. The identity and projection source are tested
separately on both lattices. All four resulting norm-form points are
nonintegral. No possible integral solution is lost in the stated chart
reduction; none has been found by it either.

If a principal generator exists, all norm-`N_j` generators form its orbit
under norm-one units. The quadratic parametrization does not calculate that
unit group or supply a bound on the integer parameters required in (7).
This is the remaining global arithmetic issue, not a further rational-point
existence problem on the quadric.

## 4. Equation-only comparison with the eight controls

**Verified application.** Specialize only the saved parent equation at the
unchanged nine original parameters. For its Weierstrass invariants take
`A=-c4/3`, `B=-2c6/27`; this is the depressed2-division cubic in
`eta=4x+b2/3`. All nine have `A B Delta≠0`. The point `(0,-3,2A,-A)` and the
norm-one matrix `1+(A/B)eta` verify exactly on every case.

| Cases | Cubic discriminant sign | Rational norm quadric |
| --- | --- | --- |
|302 and six controls|Positive|Yes, by the same formula|
|scale-0131232 and scale-0487239|Negative|Yes, by the same formula|

No exceptional points or relative seed ideals enter this panel computation.
The sign variation does not affect this rationality result. Neither the
quadric nor its rational point separates302 from these controls.

**Object distinction.** The norm surfaces in (1) are already trivial torsors
under `T`, because `a_j` is a rational point. This is not a statement that an
elliptic2-cover is trivial or has a rational point. In the present302
retrospective experiment, the marked seed's zero Sha image comes from its
independently certified elliptic rational point, not from this norm-quadric
calculation. Principality, unit squareclasses, point admissibility and
elliptic Selmer/Sha membership remain separate questions.

## 5. Certificates, replay and next step

The [immutable packet](../../artifacts/generated-results/elliptic-curves/det1092_seed_norm_quadric_v1/)
contains a generic-only quadric, quadratic parametrization, cubic torus maps,
the two integral matrices, all chart exceptions, and the nine-control replay.
Polynomial records store coefficient/exponent pairs; matrices use the printed
basis order. The two relative matrices remain openly seed-derived.

The independent checker uses only rational companion-matrix arithmetic,
not number-field or ideal backends or producer imports. It verifies norms,
inverse maps, both directions of the projection identity, determinant and
every integral transport. The separate universal check proves the identities
with symbolic `A,B` before replaying the nine control equations.

```sh
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/construct_det1092_seed_norm_quadric.sage
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_seed_norm_quadric.sage
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_norm_quadric_family.sage
```

All finish below one second in Sage10.9, within25-second caps. No class/unit
group, factorization, norm/point search, new parameter, later cascade input,
production mutation or detached process was used. No formal verification
or external review is claimed.

The next unresolved requirement is still an integral solution of one of
(1), or a proof excluding both; equivalently, solve the exact gate (7).
Even success there would classify the known seed's unit-versus-ideal origin,
not yet produce a seed on a new parameter without its retrospective ideal
input. The active seed-construction goal is not complete.
