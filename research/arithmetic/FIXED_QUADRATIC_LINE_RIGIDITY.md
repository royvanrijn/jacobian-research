# Rigidity of affine target lines with a fixed quadratic factor

The two known Hasse-failing target lines of the fixed quintic map exhaust
the affine-line construction that retains one fixed irreducible quadratic
factor. A third line retains a different quadratic field and has no Hasse
failures. Moving quadratic factors remain outside this classification.

Use the fixed map `Phi` and the unnormalized coordinate `pi=-u` from
[the ambient arithmetic criterion](FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md).
Its scalar inverse polynomial is

\[
 E_{\pi,b,c}(S)=\pi^5S^5-\frac32\pi^4S^4+\frac32\pi S^3
 +\left(-\frac54\pi-\frac9{32}b\right)S^2
 +\frac9{16}S-\frac9{32}c.                            \tag{1.1}
\]

A **fixed-quadratic line** means a nonconstant rational affine target
line `lambda -> (u(lambda),b(lambda),c(lambda))`, with `pi` generically
nonzero, together with one monic irreducible `q(S) in Q[S]` of degree two
which divides (1.1) over `Q(lambda)`. The polynomial `q` itself is fixed,
not merely the isomorphism class of its quotient field.

## 1. The inverse coefficient constraint

**Lemma 1.1 (a fixed factor forces a constant first coordinate).** Every
fixed-quadratic line has constant nonzero `pi`. After setting `W=pi S`,
its monic quadratic factor is `W^2+v`, where `v!=0`, and

\[
 16\pi^2v^2+9\pi-24v=0.                              \tag{1.2}
\]

Conversely any rational `pi,v` satisfying (1.2), with `pi v!=0` and
`-v` nonsquare, defines precisely the affine line

\[
 \boxed{(u,b,c)=\left(-\pi,
 -\frac{32\pi^2}{9}\left(t-\frac{3v}{2}+\frac{5}{4\pi}\right),
 -\frac{32vt}{9}\right),\qquad t\in\mathbb Q.}        \tag{1.3}
\]

**Proof.** If `pi=pi_0+pi_1 lambda` with `pi_1!=0`, the coefficient of
`lambda^5` in (1.1) is `pi_1^5 S^5`. Reduction modulo the fixed irreducible
quadratic would force `q` to divide `S^5`, which is impossible. Thus `pi`
is constant, and is nonzero by the definition.

In the coordinate `W`, the monic polynomial becomes

\[
 P(W)=W^5-\frac32W^4+\frac{3}{2\pi^2}W^3
       +\beta W^2+\frac9{16\pi}W+\delta,
 \quad \beta=-\frac5{4\pi}-\frac{9b}{32\pi^2},
 \quad \delta=-\frac{9c}{32}.                        \tag{1.4}
\]

Subtract the polynomials at two targets on the line. Their difference is
`Delta beta W^2+Delta delta`. The coefficient `Delta beta` cannot vanish:
otherwise a nonzero constant would be divisible by the quadratic, or the
line would be constant. Hence a monic quadratic dividing this difference
has no linear term and equals `W^2+v`.

Polynomial division now gives

\[
 P(W)=(W^2+v)
       \left[W^3-\frac32W^2+
              \left(\frac3{2\pi^2}-v\right)W+t\right] \tag{1.5}
\]

if and only if
`v(3/(2pi^2)-v)=9/(16pi)`, `beta=t-3v/2`, and `delta=vt`.
These are (1.2)--(1.3). Irreducibility of the quadratic is exactly the
condition that `-v` be nonsquare. \(\square\)

## 2. The rational points of the auxiliary elliptic curve

**Lemma 2.1 (the complete auxiliary point set).**

\[
 E:\ Y^2=X^3+1,\qquad
 E(\mathbb Q)=\{O,(-1,0),(0,1),(0,-1),(2,3),(2,-3)\}. \tag{2.1}
\]

**Proof.** We spell out the rank computation and adopt the classical
two-isogeny descent theorem. For an integral equation
`y^2=x^3+a x^2+b x`, with `b(a^2-4b)!=0`, let `alpha` take the
`x`-squareclass, with values `1` at infinity and `b` at `(0,0)`.
For its two-isogenous curve the coefficients are `a'=-2a`,
`b'=a^2-4b`, with analogous map `alpha'`. The theorem states that

\[
 2^{r+2}=|\operatorname{im}\alpha|
                 |\operatorname{im}\alpha'|.          \tag{2.2}
\]

The images are subgroups of rational squareclasses represented by signed
squarefree divisors `d` of `b` or `b'`; the remaining classes occur exactly
when the corresponding quartic
`N^2=d M^4+a M^2 e^2+(b/d)e^4` has a primitive integral solution.
The exact imported statement and quartic construction are in Andrej
Dujella, *Number Theory*, Skolska knjiga, Zagreb, 2021, section 15.5,
printed pp. 507--509, equations (15.29)--(15.31) and the rank formula on
p. 508: [author PDF](https://web.math.pmf.unizg.hr/~duje/pdf/section15-5.pdf),
[retained copy](sources/dujella-number-theory-2021-section15-5.pdf).
This is a classical theorem adopted from that treatment; the underlying
general descent theorem is not newly proved here.

Translate `x=X+1` on (2.1). The resulting curve and its two-isogenous
curve are

\[
 E_0:y^2=x^3-3x^2+3x,\qquad
 E_1:y^2=x^3+6x^2-3x.                                \tag{2.3}
\]

For `E_0`, a nonzero rational `x` must be positive, since
`x^2-3x+3` is positive on the real line. Thus the descent image is
contained in `{1,3}`; both occur at infinity and `(0,0)`.

For `E_1`, the possibilities are `{1,-1,3,-3}`, with `1,-3` present.
The class `-1` would require

\[
 N^2=-M^4+6M^2e^2+3e^4,\qquad (M,e)=1.               \tag{2.4}
\]

Modulo three this forces `3|M,N`; hence `3` does not divide `e`. Modulo
nine the right-hand side is then three, whereas the left-hand side is
zero. This contradiction excludes `-1`. The image is a group containing
`-3`, so `3` is excluded too: its product with `-3` is the squareclass
`-1`. The two images in (2.2) therefore both have order two, proving rank
zero.

By the Mordell--Weil theorem all rational points are now torsion.
Lutz--Nagell and its discriminant bound give integral coordinates and
either `Y=0` or `Y^2|27`. The possibilities are `Y=0,+/-1,+/-3`, giving
exactly the five affine points in (2.1); all satisfy the curve equation.
These imported torsion statements are Dujella's Theorem 15.5 and
Proposition 15.6, printed p. 488, with Mordell--Weil at Theorem 15.3,
p. 486: [author PDF](https://web.math.pmf.unizg.hr/~duje/pdf/section15-3.pdf),
[retained copy](sources/dujella-number-theory-2021-section15-3.pdf).
The proposition's displayed proof and the quartic application above were
inspected; the general Mordell--Weil and Lutz--Nagell proofs remain
classical imported inputs. \(\square\)

The [source record](sources/dujella-number-theory-2021-context.json) gives
the exact edition, author-hosted files and hashes. The author identifies
the book as *Number Theory*, ISBN 978-953-0-30897-8, translated by Petra
Svob; “Elliptic curves” is chapter 15, not the title of a separate book.

## 3. Exactly three fixed-quadratic affine lines

**Theorem 3.1 (complete fixed-quadratic line classification).** The only
fixed-quadratic lines are the following three. The last two columns give
the fixed monic factor in the original coordinate `S` and its quadratic
field.

| `(pi,v)` | Normalized target line for `Phi` | Fixed factor | Field |
|---|---|---|---|
| `(1,3/4)` | `y_a=(-1,32a/9,(8a+1)/3)` | `S^2+3/4` | `Q(sqrt(-3))` |
| `(-2,3/4)` | `z_r=(2,56/3-1024r/9,-7/6-64r/3)` | `S^2+3/16` | `Q(sqrt(-3))` |
| `(-2,-3/8)` | `w_a=(2,8/3+128a/9,-1/6-4a/3)` | `S^2-3/32` | `Q(sqrt(6))` |

In particular the first two are exactly the lines retaining one fixed
quadratic factor with quotient field `Q(sqrt(-3))`.

**Proof.** Equation (1.2) has the transformation

\[
 X=\frac{8v}{3},\qquad Y=\frac{32v^2\pi+9}{9},\qquad
 Y^2-X^3-1=\frac{64v^2}{81}(16\pi^2v^2+9\pi-24v).    \tag{3.1}
\]

Thus it gives a finite rational point on (2.1). Since `v!=0`, `X!=0`;
the only possible points are `(2,3),(2,-3),(-1,0)`. Conversely at these
points the inverse formulas

\[
 v=3X/8,\qquad \pi=2(Y-1)/X^2                         \tag{3.2}
\]

give the three pairs displayed in the table. Each quadratic is irreducible.
Formula (1.3) gives the three target lines. For the table's parameters,
take respectively `t=-a-1/8`, `t=7/16+8r`, and `t=-a-1/8` in (1.5).
The factor in `S` is obtained by dividing `q(pi S)` by `pi^2`.
\(\square\)

The `pi=0` plane is not a hidden fourth case. Its inverse polynomial has
degree at most two, with fixed nonzero linear coefficient `9/16`. If a
fixed irreducible quadratic divided it along a nonconstant affine line,
that coefficient would fix the scalar multiple and then fix both target
coordinates, a contradiction. Independently, every target in this plane
has a rational source point by the ambient note's boundary reconstruction.

## 4. Exactly two of these lines can contain Hasse failures

**Theorem 4.1 (the third line is locally soluble only at rational cubes).**
The fibre at `w_a` is everywhere locally soluble if and only if `a` is a
nonzero rational cube. These are exactly its rational-point parameters.
Consequently the line `w_a` has no Hasse-failing fibres.

**Proof.** In the coordinate `W=-2S`, its exact inverse factorization is

\[
 E_{-2,b_a,c_a}(S)
   =\left(W^2-\frac38\right)
          \left[\left(W-\frac12\right)^3-a\right].    \tag{4.1}
\]

The quadratic has discriminant `3/2`, of squareclass six. The resultant
with the cubic is

\[
 \left(a+\frac{11}{16}\right)^2-\frac{243}{512}
 =\left(a+\frac{11}{16}\right)^2-6\left(\frac9{32}\right)^2,
                                                               \tag{4.2}
\]

which is nonzero for every rational `a`. For `a!=0` the cubic is
separable, so the entire fibre has algebra

\[
 \mathbb Q(\sqrt6)\times\mathbb Q[Z]/(Z^3-a).         \tag{4.3}
\]

If `a` is a nonzero cube this has a rational point. Otherwise the cubic
is irreducible with quadratic resolvent `Q(sqrt(-3))`, which differs from
`Q(sqrt6)`. The ambient note's rank-five theorem FMA2 excludes everywhere
local solubility for that wrong-resolvent product. At `a=0`, derivative
localization removes the triple cubic root and leaves the nonsplit
quadratic field `Q(sqrt6)`, which is not everywhere locally soluble.
This proves the assertion. \(\square\)

The first line has the existing exact parameter criterion `a in H`, and
the second has the criterion `r=a+(3/32)^3/a` for some `a in H`, as in
FMA4 of [the ambient note](FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md).
Both contain Hasse failures. Therefore these are exactly the affine lines
with a fixed irreducible quadratic factor that contain any Hasse failures.
This statement does not classify lines on which the quadratic factor moves.

**Corollary 4.2 (no integral targets on the fixed-quadratic lines).** None
of the three lines contains a target in `Z^3`.

**Proof.** Their last two target coordinates satisfy, respectively,

\[
 12c-9b=4,\qquad 48c-9b=-224,\qquad 96c+9b=8.          \tag{4.4}
\]

For integral `b,c` every left-hand side is divisible by three, whereas
none of the right-hand sides is. \(\square\)

Integer *parameters* in the earlier counting theorem therefore still
give rational, nonintegral target coordinates. Integral ambient target
counts require a different family.

## 5. The remaining quadratic-factor incidence chart

The fixed-factor rigidity has a precise limit: allowing the quadratic
itself to move gives a three-parameter incidence chart.

**Proposition 5.1 (rational chart for noncentered quadratic factors).** For
`pi,u,v in Q` with `pi u!=0`, put

\[
 A=\frac3{2\pi^2},\quad C=\frac9{16\pi},\quad
 r=-\frac32-u,\quad s=A-ur-v,\quad t=\frac{C-vs}{u}.
                                                               \tag{5.1}
\]

Then the target

\[
 \left(-\pi,-\frac{32\pi^2}{9}
           \left(t+us+vr+\frac5{4\pi}\right),
            -\frac{32vt}{9}\right)                   \tag{5.2}
\]

has normalized inverse factorization

\[
 P(W)=(W^2+uW+v)(W^3+rW^2+sW+t).                   \tag{5.3}
\]

Every rational factorization of this type with nonzero quadratic linear
coefficient occurs uniquely from the displayed coefficient data. In
particular the quadratic field is `Q(sqrt(-3))` precisely when
`v=(u^2+3d^2)/4` for some nonzero rational `d`.

**Proof.** Compare coefficients of `W^4,W^3,W` in (5.3) with (1.4), in
that order. They give exactly `r,s,t` in (5.1). The coefficients of
`W^2,1` then give (5.2). The field condition is the quadratic discriminant
condition `u^2-4v=-3d^2`. \(\square\)

To obtain Hasse failures from this chart one must still impose
irreducibility and coprimality of the two factors, the matching cubic
resolvent, and the finite local conditions of FMA3. In the specified
`Q(sqrt(-3))` branch the matching condition is
`disc(W^3+rW^2+sW+t)=-3z^2` with nonzero rational `z`.
The chart supplies those equations explicitly; it does not assert their
rational parametrization or an ambient counting asymptotic.

## 6. Verification scope

Run `.venv/bin/python research/arithmetic/verify_fixed_quadratic_line_rigidity.py`.
The [checker](verify_fixed_quadratic_line_rigidity.py) verifies the exact
polynomial coefficients, elliptic transformations, three line formulas,
third-line resultant and discriminants, incidence chart, the complete
primitive residue obstruction modulo nine in (2.4), and the finite
Lutz--Nagell candidate list. The [receipt](fixed_quadratic_line_rigidity_verification.json)
pins its source and results. These controls do not formalize the imported
general descent, Mordell--Weil, Lutz--Nagell, or arithmetic local-global
theorems. No formal verification, independent whole-proof replay,
external review or publication is claimed.
