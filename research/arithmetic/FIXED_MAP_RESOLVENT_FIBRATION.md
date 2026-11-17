# A genus-two fibration for the matching-resolvent target locus

The general quadratic-times-cubic branch of the fixed quintic map has a
concrete arithmetic description beyond the translated-pure-cubic surface.
Its matching-resolvent equation is a family of quintic hyperelliptic
curves, generically of genus two. The centered quadratic exception consists
exactly of the three previously classified affine lines.

The result organizes all possible Hasse-failing targets of this one map.
It does not give a rational parametrization of the entire locus, a rational
point classification of the curves, or an ambient asymptotic.

## 1. The general factorization chart

Retain `Phi`, `pi=-u` and `W=pi S` from
[the ambient criterion](FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md).
Its normalized inverse polynomial is

\[
 P(W)=W^5-\frac32W^4+A W^3+\beta W^2+C W+\delta,
 \quad A=\frac3{2\pi^2},\quad C=\frac9{16\pi},
\]
\[
 \beta=-\frac5{4\pi}-\frac{9b}{32\pi^2},\qquad
 \delta=-\frac{9c}{32}.                              \tag{1.1}
\]

For rational `pi,ell,v` with `pi ell!=0`, put

\[
 r=-\frac32-\ell,\qquad s=A-\ell r-v,
 \qquad t=\frac{C-vs}{\ell},\qquad
 h(W)=W^3+rW^2+sW+t.                                 \tag{1.2}
\]

The associated target is

\[
 \boxed{\left(-\pi,
 -\frac{32\pi^2}{9}\left(t+\ell s+vr+\frac5{4\pi}\right),
 -\frac{32vt}{9}\right).}                            \tag{1.3}
\]

Then `P(W)=(W^2+ell W+v)h(W)`. This is the exact coefficient chart of
Proposition 5.1 in [the rigidity note](FIXED_QUADRATIC_LINE_RIGIDITY.md),
not an assumption that every cubic field has a translated pure cubic
*polynomial* as this specified factor.

Define

\[
 D=\ell^2-4v,\qquad
 \Delta=r^2s^2-4s^3-4r^3t-27t^2+18rst,
 \qquad F_{\pi,\ell}(v)=D\Delta.                     \tag{1.4}
\]

Thus `D` and `Delta` are respectively the quadratic and cubic
discriminants. Restrict the factorization incidence to

\[
 \pi\ell D\Delta\operatorname{Res}_W(W^2+\ell W+v,h)\ne0.
                                                               \tag{1.5}
\]

This is an algebraic open subset. Irreducibility of each factor over `Q`
is a further arithmetic condition and is not being described as a
Zariski-open condition on rational parameters.

## 2. The resolvent double cover

**Theorem 2.1 (exact matching-resolvent cover).** On (1.5), the two
discriminants have the same rational squareclass if and only if

\[
 \boxed{Y^2=F_{\pi,\ell}(v)=
   (\ell^2-4v)\bigl(r^2s^2-4s^3-4r^3t-27t^2+18rst\bigr)}          \tag{2.1}
\]

has a nonzero rational `Y`. When both factors are irreducible, the
resulting fibre is Hasse-failing precisely when, at every finite prime
ramified in the cubic splitting field, either `D` is a square over `Q_p`
or `h` has a root over `Q_p`.

**Proof.** With both discriminants nonzero,
`Delta/D` is a rational square if and only if `D Delta` is a rational
square: multiply or divide by the square `D^2`. This proves the cover
equation. When the factors are irreducible, their quotient algebras are a
quadratic field and a cubic field; matching the nonsquare discriminant
also makes the cubic non-Galois. The finite local condition is exactly
Theorem FMA2/FMA3 of the ambient note, applied to the complete fibre
identified by (1.3). \(\square\)

The test is effective in the precise earlier sense. After choosing a
monic integral model `H` of `h`, only primes dividing `Disc(H)` need be
checked. At such a prime, a cubic root can be tested modulo
`p^{floor(v_p(Disc(H))/2)+1}` by the optimal cubic root test in Lemma
FMA3a. No new analytic or local-solubility oracle is assumed here.

## 3. The generic fibres have genus two

**Theorem 3.1 (generic genus-two fibres).** Regard (2.1) as a fibration
over the base with coordinates `(pi,ell)` and `pi ell!=0`.
The polynomial `F_{pi,ell}` has degree five in `v`, with leading
coefficient `108/ell^2`. Its discriminant in `v` is not the zero rational
function. Consequently, over the nonempty base open subset where that
discriminant does not vanish, the smooth projective fibre has genus two.

**Proof.** In (1.2), `s` is linear in `v` with leading coefficient `-1`,
and `t` is quadratic with leading coefficient `1/ell`. The only degree
four term of the cubic discriminant (1.4) is the leading part of
`-27t^2`. Hence `Delta` has degree four and leading coefficient
`-27/ell^2`, and `D Delta` has the claimed degree and leading coefficient.

At the rational specialization `pi=ell=1`, the quintic is

\[
 F_{1,1}(v)=\frac{(4v-1)
   (6912v^4-67840v^3+205216v^2-194256v+59043)}{256}.    \tag{3.1}
\]

Its exact discriminant is

\[
 \operatorname{Disc}_v(F_{1,1})
    =-2^{24}3^3 7^5 11^3 29^3\ne0.                  \tag{3.2}
\]

A smaller squarefreeness certificate is available modulo five. Since
the denominator in (3.1) is a unit modulo five, its reduction is

\[
 f(v)=3v^5+3v^4+4v^3+3v+2\in\mathbb F_5[v],
\]

and it satisfies

\[
 (v^2+2)f(v)+(v^4-2v-1)f'(v)=1\pmod5.                \tag{3.3}
\]

Thus this specialization is squarefree. Either (3.2) or (3.3) proves
that the general discriminant does not vanish identically. For any
specialization with squarefree degree-five right-hand side, the
projective double cover of the `v`-line branches at the five distinct
finite roots and at infinity, whose pole order is five. Characteristic
zero makes all six branches tame. Riemann--Hurwitz gives
`2g-2=2(-2)+6=2`, hence `g=2`. \(\square\)

The exceptional base divisor `Disc_v(F_{pi,ell})=0` is retained. Its
singular or lower-genus fibres can still contain Hasse targets and may
support useful special families. The generic-genus statement does not
exclude or classify those fibres. Nor does it prove that the total
three-dimensional double-cover variety has no rational parametrization.

**Corollary 3.2 (rational two-torsion on the fibre Jacobian).** Every
smooth genus-two fibre of (2.1) has a nonzero rational two-torsion class
on its Jacobian.

**Proof.** The curve has the rational branch point
`P_0=(ell^2/4,0)` and the unique rational point `P_infinity` above infinity.
The function `v-ell^2/4` has divisor `2P_0-2P_infinity`, so the rational
degree-zero class `P_0-P_infinity` is killed by two. It is nonzero: if
principal, a function with that divisor would define a degree-one map to
the projective line, making the curve rational, contrary to genus two.
\(\square\)

**Proposition 3.3 (height and all targets on a fixed factor curve).** Fix
rational `pi,ell` with `pi ell!=0`. Before imposing matching resolvent or
factor irreducibility, the target map (1.3) satisfies

\[
 H((-\pi,b(v),c(v)))\asymp_{\pi,\ell} H(v)^3.
                                                               \tag{3.4}
\]

Consequently the number of its rational target images of height at most
`B` is `Theta_{pi,ell}(B^{2/3})`. Removing the finite set where (1.5)
fails does not change this order. This is a count of all factorized
targets on that curve, not of its Hasse-failing subset.

**Proof.** The polynomial `b(v)` has degree two and leading coefficient
`-32pi^2/(9ell)`; `c(v)` has degree three and leading coefficient
`-32/(9ell)`. Homogenizing the coordinates of (1.3) to degree three gives

\[
 [N^3:-\pi N^3:N b_2(M,N):c_3(M,N)],                 \tag{3.5}
\]

where `b_2,c_3` are the respective degree-two and degree-three
homogenizations. These forms have no common projective zero: for
`N!=0` the first coordinate is nonzero, and for `N=0,M!=0` the last
coordinate equals `-32M^3/(9ell)`. Thus (3.5) is a morphism on the whole
projective line.

Here is an elementary height justification, including projective content.
Clear a fixed common denominator `D_0>0`, so the first homogeneous form
is `D_0 N^3` and the last is an integral cubic with nonzero leading
coefficient `a_3`. For coprime integers `m,n`, their common content is
bounded by `D_0|a_3|^3`. Indeed, at a prime not dividing `n`, the first
form bounds its valuation by `v_p(D_0)`. At a prime dividing `n`, put
`k=v_p(n)` and `L=v_p(a_3)`. If `k>L`, the leading term of the last
form has strictly smallest valuation `L`; otherwise the first form
bounds it by `v_p(D_0)+3L`. This proves the content bound prime by prime.
On the real compact set `max(|m|,|n|)=1`, the maximum absolute value of
the homogeneous forms has a positive minimum because they have no common
zero. Homogeneity and the content bound give the lower bound in (3.4).
The upper bound follows directly from their fixed coefficients.

The number of rational parameters of height at most `T` has order `T^2`,
by the elementary coprime integer-pair count. A target has at most two
parameters, because its `b` coordinate is a nonconstant quadratic.
Together with (3.4) this gives both bounds of the stated target order.
In fact `Q(b,c)=Q(v)`: the intermediate extension degree divides both
the polynomial degrees two and three. This confirms that the image is
a rational cubic curve. Finally `D`, `Delta`, and the quadratic--cubic
resultant are nonzero polynomials in `v`; the first two have degrees one
and four, and the resultant has leading term `v^4/ell^2`. Their zeros
remove only finitely many parameters. \(\square\)

**Corollary 3.4 (finitely many matching-resolvent targets on each smooth
fixed-base curve).** Fix rational `pi,ell` with `pi ell!=0` and
`Disc_v(F_{pi,ell})!=0`. There are only finitely many admissible rational
matching-resolvent targets (1.3), and hence only finitely many Hasse-failing
targets on that fixed-base curve. Their count is `O_{pi,ell}(1)` as `B`
varies. No effective bound or bound uniform in `(pi,ell)` is asserted.

**Proof.** Theorem 3.1 gives a smooth projective geometrically integral
genus-two curve over `Q`. The Mordell conjecture proved by Faltings says
that such a curve has finitely many rational points. The required
primary statement is G. Faltings, *Endlichkeitssatze fur abelsche
Varietaten uber Zahlkorpern*, Inventiones Mathematicae **73** (1983),
349--366, Satz 7 on p. 365 (PDF p. 17):
[primary scan](https://math.uchicago.edu/~drinfeld/Deligne%27s_conjecture_Manin_conf/Faltings_argument/Faltings.pdf),
[retained copy](sources/faltings-1983-endlichkeitssaetze.pdf).
The exact statement, whose genus hypothesis is `g>=2`, was visually
inspected in that primary source; its deep proof is adopted as an import,
not independently replayed here. Rational points of the admissible cover
are a subset of those projective rational points, and their target image
is therefore finite. \(\square\)

It follows that infinitely many Hasse targets with a fixed nonzero pair
`(pi,ell)` can occur only on the exceptional base divisor. The corollary
does not make their union over varying pairs finite, and it does not
provide an effective procedure for listing every rational point of each
smooth fibre. The earlier local test remains an effective decision for a
supplied rational target.

**Proposition 3.5 (exact exceptional base equations).** Write `p=pi` and
`l=ell` just in this proposition. Define

\[
 \mathcal A=(2l+3)^2\bigl((20l^2+12l-27)p^2+54\bigr)-243p,
\]

and `mathcal B=sum B_j(l)p^j`, `mathcal R=sum R_j(l)p^j` by the following
table; all omitted coefficients are zero.

| `j` | `B_j(l)` | `R_j(l)` |
| --- | --- | --- |
| 8 | `l^3(2l+3)^4(5l+8)` | `0` |
| 7 | `-3l(7l+12)(23l^2+48l+18)` | `0` |
| 6 | `376l^6+2088l^5+4302l^4+3888l^3+1296l^2+729` | `l^5(l+2)(25l^2+30l-36)` |
| 5 | `-162l(7l+12)` | `-6l(5l^3-72l^2-108l+108)` |
| 4 | `6l(2l+3)(76l^2+141l+36)` | `260l^6+696l^5+144l^4-432l^3+729` |
| 3 | `-1458` | `324l(l+6)` |
| 2 | `18l(71l+108)` | `12l^2(83l^2+132l-36)` |
| 0 | `729` | `1152l^2` |

On `pl!=0`, one has the exact identities

\[
 \operatorname{Disc}_v(\Delta)
   =-\frac{3\mathcal A^3\mathcal B}{256l^8p^{14}},\qquad
 \operatorname{Res}_v(D,\Delta)
   =-\frac{3\mathcal R}{l^2p^6},                    \tag{3.6}
\]
\[
 \boxed{\operatorname{Disc}_v(F)
   =-\frac{27\mathcal A^3\mathcal B\mathcal R^2}
                {256l^{12}p^{26}}.}               \tag{3.7}
\]

Consequently the exceptional base is exactly
`mathcal A mathcal B mathcal R=0` in the nonzero `(pi,ell)` base.
The equation `mathcal R=0` records collision of the rational branch
`v=ell^2/4` with a root of `Delta`; the other two equations record
repeated roots of the quartic `Delta`. These equations do not yet
classify the genera or rational points of every exceptional fibre.

**Proof.** Substitute (1.2) into (1.4), keeping the rational denominators.
The quartic numerator `N=256l^2p^6 Delta` has discriminant

\[
 \operatorname{Disc}_v(N)
  =-2^{40}3l^4p^{22}\mathcal A^3\mathcal B.
\]

Scaling a degree-four polynomial multiplies its discriminant by the sixth
power of the scale. This gives the first identity in (3.6). Evaluation
at the root of `D` gives
`Res_v(D,Delta)=(-4)^4 Delta(l^2/4)`, which expands to the second
identity. The product-discriminant identity
`Disc(D Delta)=Disc(D)Disc(Delta)Res(D,Delta)^2`, with `Disc(D)=1`,
then gives (3.7). These are exact polynomial identities, checked before
and after clearing denominators in the accompanying checker. At
`p=l=1`, the values are `mathcal A=1232`, `mathcal B=22736`,
`mathcal R=7424`, agreeing with (3.2). \(\square\)

**Proposition 3.6 (the first exceptional component has finitely many
rational base pairs).** The exceptional component `mathcal A=0` is
birational over `Q` to the genus-two curve

\[
 y^2=59049-216x^5(5x-24).                            \tag{3.8}
\]

It therefore has only finitely many rational pairs `(pi,ell)` with
`pi ell!=0`. This is a statement about the base pairs on this one
component, not finiteness of all target parameters above those pairs.

**Proof.** Put `x=2ell+3`. Then `mathcal A=0` becomes

\[
 x^3(5x-24)\pi^2-243\pi+54x^2=0.
\]

On `x(5x-24)!=0`, completing the square gives (3.8) with

\[
 y=2x^3(5x-24)\pi-243,\qquad
 \pi=\frac{y+243}{2x^3(5x-24)},\qquad
 \ell=\frac{x-3}{2}.                                \tag{3.9}
\]

These are inverse rational maps on the indicated opens; the original
base additionally excludes `ell=0`, or `x=3`, and `pi=0`.
Writing `G(x)=59049-216x^5(5x-24)`, one has
`G'(x)=-6480x^4(x-4)`. Its only possible common roots with its
derivative are `0,4`, and `G(0)=59049`, `G(4)=943785` are both
nonzero. Thus `G` is squarefree of degree six. The smooth projective
double cover has six finite branch points and no branch at infinity,
so Riemann--Hurwitz again gives genus two. Faltings's theorem, with
the exact primary import used in Corollary 3.4, makes its rational-point
set finite.

The two cases removed in completing the square are explicit. At `x=0`,
the equation forces `pi=0`, outside the base. At `x=24/5`, it is linear
and has the unique solution `pi=(2/9)x^2=128/25`; this gives the
admissible base pair `(pi,ell)=(128/25,9/10)`. Adding this one pair
preserves finiteness. No effective list or height bound for the rational
points of (3.8) is asserted. \(\square\)

**Proposition 3.7 (an explicit elliptic degeneration).** At the rational
exceptional base pair `(pi,ell)=(128/25,9/10)`, the resolvent cover has
the exact equation

\[
 Y^2=\frac{(400v-81)(819200v-243483)^2
              (20480000v^2-110592000v+143467857)}
             {41231686041600000000}.                \tag{3.10}
\]

Its smooth projective normalization is an elliptic curve over `Q`.

**Proof.** This pair satisfies `mathcal A=0` by Proposition 3.6.
Direct substitution into (1.4) gives (3.10). Put

\[
 J(v)=(400v-81)(20480000v^2-110592000v+143467857),\quad
 K=41231686041600000000.
\]

On `819200v-243483!=0`, the substitutions
`Z=Y/(819200v-243483)` and `Y=(819200v-243483)Z` are inverse maps
to the cubic cover `Z^2=J(v)/K`. The exact discriminant

\[
 \operatorname{Disc}_v(J)
   =2^{33}3^6 5^{14}7^3\cdot1889\cdot2709173^2\ne0
                                                               \tag{3.11}
\]

makes this cubic squarefree. Its projective double cover therefore has
four branch points and genus one. It has the rational point
`(v,Z)=(81/400,0)` and a unique rational point at infinity. With the
latter chosen as origin, the former is nonzero rational two-torsion.
The two covers have the same function field, so this smooth genus-one
curve is the normalization of (3.10).

The deleted value `819200v=243483` already has `Delta=0` and lies
outside the admissible factorization open (1.5). Thus the normalization
describes all admissible matching-resolvent targets at this base.
Factor irreducibility and the finite local tests remain necessary to
deduce Hasse failure.
\(\square\)

**Proposition 3.8 (infinitely many matching-resolvent targets on the
elliptic exceptional fibre).** At the fixed base pair
`(pi,ell)=(128/25,9/10)`, there are infinitely many distinct rational
targets (1.3) satisfying (1.5) and the matching-squareclass equation
(2.1). No assertion is made that their quadratic and cubic factors are
irreducible, or that infinitely many of these targets are Hasse-failing.

**Proof.** The normalization in Proposition 3.7 is isomorphic over `Q`
to the integral short Weierstrass model

\[
 E:\quad T^2=X^3-10758382X+6600759744                 \tag{3.12}
\]

by the inverse affine changes

\[
 X=\frac{6400}{3}v-3984,\quad
 T=\frac{20971520000}{3}Z,\qquad
 v=\frac{3(X+3984)}{6400},\quad
 Z=\frac{3T}{20971520000}.                          \tag{3.13}
\]

Direct substitution verifies the equation, and its discriminant is
`2^9 5^2 7^3 *1889*2709173^2`, which is nonzero. The rational point

\[
 P=\left(\frac{199200}{49},\frac{59459616}{343}\right)\in E(\mathbb Q)
                                                               \tag{3.14}
\]

also satisfies the equation exactly. Its `X` coordinate is not an
integer. The Lutz--Nagell theorem for a nonsingular short Weierstrass
equation with integral coefficients says that every rational torsion
point has integral coordinates. The primary source used here is
A. Dujella, *Number Theory* (2021), Theorem 15.5, p. 488, in
[the retained section 15.3](sources/dujella-number-theory-2021-section15-3.pdf),
with [author's book page](https://web.math.pmf.unizg.hr/~duje/numbertheorybook.html).
Thus `P` is nontorsion. This uses the general integral-model theorem,
not only its earlier specialization to `y^2=x^3+1`.

The distinct multiples of `P` supply infinitely many rational points
on `E`. The coordinate `v` in (3.13) has degree two as a map to the
projective line, so their `v` values form an infinite set. The polynomial
conditions excluded by (1.5) remove only finitely many values, as proved
in Proposition 3.3. The additional deleted linear factor in Proposition
3.7 already occurs in that excluded set. Hence infinitely many points
give admissible rational points on (2.1). Finally each rational target
has at most two such `v` parameters by its quadratic `b(v)` coordinate.
The resulting target images are therefore infinite. The irreducibility
and finite local conditions of Theorem 2.1 remain separate arithmetic
requirements. \(\square\)

The displayed point also illustrates why the local tests remain essential.
It gives `v=73953/19600`, with irreducible quadratic factor of
discriminant `-6(54/35)^2`. Its cubic factor satisfies

\[
 h(Z+1)=Z^3+\frac35Z^2-
       \frac{134708181}{40140800}Z+
       \frac{290236959237}{78675968000}.             \tag{3.15}
\]

Each nonleading coefficient has three-adic valuation exactly one, so
this cubic is Eisenstein at three. It is irreducible over `Q` and has
no `Q_3`-root. Its discriminant has the same squareclass as the
quadratic by the cover equation. The fibre is consequently a full
rank-five product of `Q(sqrt(-6))` and a non-Galois cubic field with that
quadratic resolvent. Such a cubic field cannot be a pure cubic field,
whose quadratic resolvent is `Q(sqrt(-3))`. Both components lack a
three-adic point, so this particular matching-resolvent fibre is not
Hasse-failing.

## 4. The complete centered exception and coverage of all Hasse targets

**Theorem 4.1 (centered factors and complete coverage).** Every rational
target of `Phi` with a monic irreducible quadratic inverse factor whose
linear coefficient in `W` is zero lies on one of the three fixed-quadratic
lines of the rigidity theorem. Only the original line `y_a` and second
line `z_r` among them can contain Hasse failures.

Every other Hasse-failing target has exactly two rational lifts to (2.1),
with opposite nonzero `Y`, satisfying the factor irreducibility and finite
local conditions of Theorem 2.1. This includes the exceptional base
divisor from section 3.

**Proof.** If the quadratic is `W^2+v`, comparison with (1.1) gives

\[
 v(3/(2\pi^2)-v)=9/(16\pi),\qquad
 16\pi^2v^2+9\pi-24v=0.                              \tag{4.1}
\]

This calculation applies to an individual target; it needs no affine
line hypothesis. The elliptic-curve reduction and complete rational
point calculation in Theorem 3.1 of the rigidity note show that
`(pi,v)` is `(1,3/4)`, `(-2,3/4)`, or `(-2,-3/8)`.
The remaining cubic constant is free and gives exactly its three lines.
Theorem 4.1 there excludes Hasse failures on the third line.

By FMA1 every Hasse target has `pi!=0` and squarefree quintic, and by
FMA2/FMA3 its rational irreducible factor degrees are exactly `2+3` with
matching resolvent. Its monic quadratic factor is unique over `Q`. If
noncentered, its coefficients recover a unique triple `(pi,ell,v)` in
the chart (1.2); condition (2.1) then has exactly two rational values
`Y` and `-Y`. Conversely every point satisfying Theorem 2.1 gives a
Hasse target. This proves the claimed coverage and multiplicity.
\(\square\)

The translated-pure-cubic rational surface in
[the moving-quadratic note](FIXED_MAP_MOVING_QUADRATIC_FAMILY.md) is one
explicit subfamily within this general resolvent locus. Its three-adic
integral obstruction concerns that subfamily and cannot be applied to
arbitrary cubic generators in (1.2).

## 5. Verification and the remaining arithmetic problem

Run `.venv/bin/python research/arithmetic/verify_fixed_map_resolvent_fibration.py`.
The [checker](verify_fixed_map_resolvent_fibration.py) verifies the
factorization chart and discriminants, both leading coefficients, the
displayed genus-two specialization, its exact discriminant and independent
modulo-five Bezout certificate, the rational branch, the height-map
degrees and its nonzero value at infinity, the exceptional-base
discriminant identities, the first-component birational equation and
squarefree sextic, the explicit elliptic degeneration and its nonzero
cubic discriminant, its integral elliptic model and nontorsion point,
and the centered coefficient constraint.
The [receipt](fixed_map_resolvent_fibration_checks.json) records exact
source hashes and results, including the retained primary Faltings scan.
The genus, two-torsion, and height statements use the written arguments;
the Faltings finiteness and general Lutz--Nagell theorems are explicitly
imported. Finite
specializations alone are not a proof about all rational fibres.

The remaining parametrization and counting work includes rational points
on this fibration and its exceptional base divisor, together with the
finite local conditions. The existing pointwise ambient decision
criterion is complete; a uniform rational parametrization or an ambient
counting asymptotic is a stronger, separate requirement. No formal
verification, independent whole-proof replay or external review is claimed.
