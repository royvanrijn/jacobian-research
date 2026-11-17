# A rational surface of moving quadratic and pure cubic fibres

The fixed quintic Keller map has Hasse-failing targets outside its two
previously classified affine lines. This note gives a rational surface of
factorized targets, an exact local criterion within that surface, and an
obstruction to integral Hasse targets in this particular branch. The cubic
factor here is a translated pure cubic polynomial; this condition is
stronger than requiring its field to be isomorphic to a pure cubic field.

Retain `Phi`, `pi=-u`, and the inverse coordinate `W=pi S` from
[the fixed-quadratic line note](FIXED_QUADRATIC_LINE_RIGIDITY.md).
The normalized inverse polynomial is

\[
 P(W)=W^5-\frac32W^4+\frac3{2\pi^2}W^3
 +\beta W^2+\frac9{16\pi}W+\delta,
 \quad\beta=-\frac5{4\pi}-\frac{9b}{32\pi^2},
 \quad\delta=-\frac{9c}{32}.                         \tag{1.1}
\]

## 1. The explicit quadric family

Let `Q` be the affine quadric

\[
 Z^2-5T^2+3\pi T+\frac34\pi^2=2.                    \tag{1.2}
\]

Work on its rational open subset where `pi Z(2T-pi)!=0`, and define

\[
 \alpha=T/\pi,\qquad d=Z/\pi,\qquad
 \ell=3(\alpha-1/2),\qquad
 v=\frac3{2\pi^2}+6\alpha^2-\frac92\alpha.            \tag{1.3}
\]

Set

\[
\begin{aligned}
 a&=\frac{10\alpha^4-8\alpha^3+3\alpha^2/\pi^2-3/(8\pi)}
                  {2\alpha-1},\\
 h_0&=-\alpha^3-a,\qquad
 \beta=h_0+3\ell\alpha^2-3\alpha v,\\
 b&=-\frac{32\pi^2}{9}\left(\beta+\frac5{4\pi}\right),
 \qquad c=-\frac{32v h_0}{9}.                        \tag{1.4}
\end{aligned}
\]

These formulas give the normalized target `(-pi,b,c)` of the same fixed
map `Phi`; they do not define a new Keller map.

**Theorem 1.1 (moving quadratic factorization).** At every rational point
of the specified open subset, (1.1) factors as

\[
 \boxed{P(W)=(W^2+\ell W+v)\bigl((W-\alpha)^3-a\bigr),
 \qquad \ell^2-4v=-3d^2.}                            \tag{1.5}
\]

The quadratic is irreducible and represents `K=Q(sqrt(-3))`. Conversely,
every factorization (1.5) of (1.1), with irreducible quadratic of field
`K` and `alpha!=1/2`, is obtained from these formulas, with the two signs
of `d` giving the same target.

**Proof.** The `W^4` coefficient forces `ell=3(alpha-1/2)`. The `W^3`
coefficient then forces the value of `v` in (1.3). The condition
`ell^2-4v=-3d^2` is exactly

\[
 d^2=5\alpha^2-3\alpha-3/4+2/\pi^2,
\]

which becomes (1.2). Matching the `W` coefficient gives

\[
 3v\alpha^2-\ell\alpha^3-\ell a=9/(16\pi).
\]

Solving for `a` gives (1.4); the remaining two coefficients give `b,c`.
This also proves the converse. Nonzero rational `d` makes the quadratic
discriminant nonzero of squareclass `-3`. \(\square\)

**Proposition 1.2 (rational parametrization and image dimension).** The
quadric (1.2) is a rational surface. One explicit dense parametrization is

\[
 D=n^2-5m^2+3m+3/4,\qquad
 \lambda=-\frac{2n-2m+3}{D},
\]
\[
 \pi=1+\lambda,\qquad T=1/2+m\lambda,\qquad Z=1+n\lambda,
                                                               \tag{1.6}
\]

on the open subset where its denominators and the conditions of Theorem
1.1 hold. The closure of the target image has dimension two.

**Proof.** The quadric contains `(pi,T,Z)=(1,1/2,1)`. Intersecting the
line from this point in direction `(1,m,n)` with the quadric gives the
second point (1.6). Conversely, off `pi=1`, recover
`m=(T-1/2)/(pi-1)`, `n=(Z-1)/(pi-1)`; substitution recovers (1.6).
The projective quadratic form of (1.2) is nonsingular, so this is a
birational parametrization on dense open subsets of a geometrically
integral surface.

Over a squarefree inverse quintic there are only ten choices of a monic
quadratic factor. Such a choice recovers `alpha` from `ell`, and then
recovers `T` and the two possible values of `Z`. Thus the target map has
finite fibres of size at most twenty on this open subset. That subset is
nonempty by the example in section 3 below, proving image dimension two.
\(\square\)

## 2. Exact arithmetic criterion and exceptional fibres

Let `H` be the already classified set

\[
 \mathcal H=\{a\in\mathbb Q^*\setminus\mathbb Q^{*3}:
       3\mid v_p(a)\ (p\equiv2\pmod3),\quad
       a\in\mathbb Q_3^{*3}\}.                       \tag{2.1}
\]

**Theorem 2.1 (local classification within the quadric family).** The
target (1.4) is Hasse-failing if and only if its computed parameter `a`
belongs to `H`. If `a` is a noncube, its full fibre is finite etale of
rank five with algebra

\[
 K\times\mathbb Q[U]/(U^3-a).                        \tag{2.2}
\]

For a nonzero cube `a` the fibre always has a rational point. At `a=0`
the entire fibre is `Spec K` and has no real point.

**Proof.** For a noncube `a`, the cubic is irreducible of degree three;
it is separable and coprime to the irreducible quadratic of degree two.
The derivative-unit reconstruction therefore gives (2.2). Apply the
[existing rational local criterion](../../unpublished/fixed-map-hasse-failures/sections/03-local-global.tex)
to this same algebra. It yields precisely (2.1).

If `a=r^3!=0`, the cubic has the simple rational root `W=alpha+r`, which
cannot lie in the irreducible quadratic. It survives derivative
localization and gives a rational point. More explicitly the fibre is
`Q x K x K` if the quadratic and cubic are coprime; if their quadratic
factors coincide, localization removes that repeated factor and leaves
`Q`. Thus no coprimality claim for cube parameters is being assumed.
At `a=0`, the triple root `W=alpha` is disjoint from the quadratic and
is removed, leaving exactly `K`. \(\square\)

This theorem describes an algebraic surface of factorized targets and its
Hasse criterion. It does not assert that the Hasse targets are Zariski
dense on that surface or have a positive counting density.

## 3. A Hasse-failing target outside the two affine lines

In (1.6), take `m=1,n=-2/3`. Then

\[
 \pi=17/29,\quad \alpha=5/34,\quad d=37/17,\quad
 \ell=-18/17,\quad v=4431/1156,\quad
 a=\frac{26263}{34^3}.                              \tag{3.1}
\]

Its target is

\[
 \boxed{\left(-\frac{17}{29},\frac{46972}{128673},
                           \frac{12991692}{1419857}\right).}    \tag{3.2}
\]

The integer `26263` is prime and is `1 (mod 9)`. The denominator is a
cube, so `a` satisfies every inert valuation condition and is a
three-adic cube; its valuation at `26263` is one, so it is not a rational
cube. Theorem 2.1 proves that (3.2) is Hasse-failing. Its first coordinate
is neither `-1` nor `2`, so it lies on neither previously classified line.
The quadratic factor in the coordinate `W` is
`W^2-(18/17)W+4431/1156`, and the cubic is
`(W-5/34)^3-26263/34^3`. This verifies the claim independently of any
unbounded search.

## 4. A three-adic obstruction to integral targets

**Theorem 4.1 (no integral Hasse targets in this branch).** No
Hasse-failing target in the quadric family has all three coordinates in
`Z_3`. In particular it has no integral Hasse-failing target.

**Proof.** Suppose first that `pi` is a three-adic unit and `b,c in Z_3`.
Then the monic polynomial (1.1) belongs to `Z_3[W]`. Every root of either
monic factor (1.5) is therefore integral over `Z_3`. Hasse failure implies
`a=r^3` for some nonzero `r in Q_3`. In `Q_3(zeta_3)`, the three cubic
roots are `alpha+r, alpha+zeta_3 r, alpha+zeta_3^2 r`.
The identity `(zeta_3-1)^2=-3zeta_3` gives valuation `1/2` for
`zeta_3-1`. Integral differences of roots imply `v_3(r)+1/2>=0`;
since `r` has integral valuation, `r in Z_3`, and then `alpha in Z_3`.

Consequently `ell in 3Z_3`, while (1.3) gives

\[
 v/3=1/(2\pi^2)+2\alpha^2-(3/2)\alpha
       \equiv2(1+\alpha^2)\not\equiv0\pmod3.
\]

Thus `v_3(v)=1`. The constant term `delta=v h_0=-9c/32` forces
`h_0 in 3Z_3`. Now `beta=h_0+3ell alpha^2-3alpha v` is divisible by
three, whereas `beta=-5/(4pi)-9b/(32pi^2)` is a unit. This contradiction
excludes the unit case.

Suppose next `e=v_3(pi)>=1`. In
`d^2=5alpha^2-3alpha-3/4+2/pi^2`, a valuation `v_3(alpha)<-e` would
make the leading unit squareclass five, which is nonsquare modulo three;
a valuation `v_3(alpha)>-e` would instead give the nonsquare unit class
two from `2/pi^2`. Hence `v_3(alpha)=-e`, so `T=pi alpha` is a unit.
Formula (1.4) becomes

\[
 a\pi^3=
 \frac{10T^4-8\pi T^3+3T^2-3\pi^3/8}{2T-\pi}
       \equiv2T^3\equiv\pm2\pmod9.                  \tag{4.1}
\]

For the congruence, subtract `2T^3(2T-pi)` from the numerator. The
difference is `3T^2(2T^2-2pi T+1)-3pi^3/8`, which is divisible by nine
when `pi` is divisible by three and `T` is a unit. Thus `a pi^3` is not
a three-adic cube, and neither is `a`. This already excludes Hasse
failure, without any assumption on `b,c` in this case.

An integral first target coordinate has `v_3(pi)>=0`, so the two cases
exhaust the assertion. \(\square\)

In the excluded centered case `alpha=1/2`, the coefficient conditions
give `pi=1` or `-2`, with `v=3/(4pi)`. The first is the original
`K`-factor line; the second has quadratic field `Q(sqrt(6))` and is the
Hasse-free third line of the rigidity classification. Thus only `pi=1`
belongs to the present `K` branch. It has no integral targets by Corollary
4.2 of [the rigidity note](FIXED_QUADRATIC_LINE_RIGIDITY.md), so adding
this centered exception preserves the integral obstruction. General
cubic generators, including the Cardano representatives on the second
`K`-factor line, require their own analysis.

## 5. Verification and remaining scope

Run `.venv/bin/python research/arithmetic/verify_fixed_map_moving_quadratic.py`.
The [checker](verify_fixed_map_moving_quadratic.py) verifies the exact
coefficient identities, quadric parametrization, the displayed target and
its admissible parameter, the finite primality test, the exceptional cube
and zero algebra examples, and the residue computations used in Theorem
4.1. The [receipt](fixed_map_moving_quadratic_checks.json) records exact
source hashes and controls. Finite examples do not prove rational-point
density or an asymptotic. The fixed-map reconstruction and original
local criterion retain their own proof and assurance scopes.

Moving quadratic factors with general cubic generators and other
quadratic fields remain outside this particular surface. No ambient
counting asymptotic, classification of integral Hasse targets in those
other branches, formal verification or external review is asserted.
