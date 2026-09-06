# The marked minimal carrier of a successful two-direction block

For the certified two-direction subblock on the +7 R17 fibre `08234-003`,
the minimal carrier of the **specified native character directions** is
an explicit genus-one curve of degree four over the parameter line:

\[
\boxed{C:\quad u^2=f(t),\quad v^2=g(t)}
\]

where

\[
\begin{aligned}
f(t)&=16496921457+11037654810t+1654807609t^2,\\
g(t)&=4897771825+2856794060t+580965316t^2.
\end{aligned}
\]

These are `0b2d0` and `19e45`, with the rational square scalings removed
from their native cover equations. The carrier has a rational origin and
its Jacobian has exact rank three. This does not make every rational t
soluble: the question is whether the **fibre of C over t** contains a
rational point.

The exact necessary and sufficient condition away from the two branch
divisors is

\[
\boxed{f(t),g(t)\in\mathbf Q^{\times2}.}
\]

More substantively, an elimination using the generic three-cover relation
produces a degree-12 finite intersection scheme with factor degrees
**1+11**. Its sole rational point lies over

\[
\boxed{65t+288=0.}
\]

At that parameter the generic equations construct all three rational
lifts, and a separate exact witness certificate proves that their quotient
span has dimension two. The degree-11 factor is irreducible modulo 73,
which proves irreducibility over Q. Thus the extra rationality here is a
rational component of a finite intersection scheme, not a forced
degree-one intersection as in the preceding +8 pair.

This characterizes a **subblock**, not the entire observed +7 jump. Absolute
minimality over other auxiliary constructions or undisplayed sections,
and a carrier for the whole seven-dimensional quotient, remain UNKNOWN.

## Minimality and the isogeny lifting condition

The two quadratics are squarefree with disjoint geometric branch divisors.
Their classes in Q(t)*/Q(t)*² are independent. Thus

\[
\mathbf Q(C)=\mathbf Q(t)(\sqrt f,\sqrt g)
\]

has degree four, and Riemann–Hurwitz gives genus one. Equivalently, C is
the smooth intersection of two quadrics in P3 after homogenizing the
displayed equations.

Each native point has x-coordinate \(x_0(t)+\sqrt q\,x_1(t)\), with
\(x_1\ne0\). Its coordinates therefore recover its quadratic root. Making
both marked points rational over a function-field extension requires both
roots, hence an extension containing Q(C). The same conclusion holds for
the entire two-dimensional native character block: its two independent
nonzero characters must both become trivial, regardless of a change of
basis within that block. No smaller parameter extension carries this
marked block. A curve dominating C cannot have genus zero, so genus one
is minimal as well. This argument does not cover alternative blocks in
the full Mordell–Weil group.

There is a tempting smaller cover,

\[
C':\quad w^2=f(t)g(t),\qquad
\pi:C\longrightarrow C',\quad (t,u,v)\longmapsto(t,uv).
\]

It also has genus one, but its rational points need not lift to C. The
involution \((u,v)\mapsto(-u,-v)\) has no geometric fixed points because
the branch divisors are disjoint. After choosing the known rational origin
on C and its image on C', the map is a degree-two isogeny. An independent
Jacobian calculation verifies that the corresponding elliptic models are
indeed related by a degree-two isogeny.

For a non-branch rational point \((t,w)\in C'(\mathbf Q)\), the exact
obstruction to lifting is

\[
\delta_\pi(t,w)=[f(t)]\in\mathbf Q^*/\mathbf Q^{*2}.
\]

It vanishes precisely when u is rational, in which case \(v=w/u\) is
rational too. The origin used below has square f-value, so this formula
has the normalization required by the pointed isogeny descent. A product
square alone only makes the two squareclasses equal; it does not make
them trivial. This is a **global solubility** condition, not a visibility
feature or an original-curve rank prediction.

At \(t_0=-288/65\), the rational origin can be taken as

\[
(t_0,u_0,v_0)=
\left(-\frac{288}{65},\frac{44253}{5},\frac{3924473}{65}\right).
\]

The auxiliary Jacobian has minimal Weierstrass model

\[
Y^2+XY=X^3
-235950017004391353572674922721546566X
+33693852184918766381840507850394292045967182353721700.
\]

A bounded PARI 2-descent gives rank interval [3,3], rational 2-torsion
dimension one, and 2-Selmer dimension four. It returns two explicit
independent points; its rank lower bound three is the descent result,
not a claim that three generators were returned. This distinction follows
the documented semantics of
[PARI ellrank](https://pari.math.u-bordeaux.fr/dochtml/html-stable/Elliptic_curves.html#ellrank).
Consequently C has infinitely many rational points. Its projection to t
still has degree four, and imposing another native square condition remains
an additional lifting problem.

## Three lifts, one relation, and two directions

Write A=`01333`, B=`0b2d0`, and D=`19e45`. The retained relation is

\[
P_A-P_B+P_D=S(t_0),
\]

where S is a rational generic section. In the published generic basis its
word is

```text
[-1,-1,0,2,-2,2,-1,0,0,-1,-1,-1,0,1,1,0,0].
```

The third primitive square polynomial is

\[
h(t)=3255283501715844+1254304425186516t+125950947365881t^2.
\]

Requiring all three displayed native lifts over a variable parameter means
using

\[
C_3:\quad a^2=h(t),\quad u^2=f(t),\quad v^2=g(t).
\]

Its degree over t is eight and its genus is five. The extra lift
\(C_3\to C\) is a genuine quadratic condition; quotient dependence at
one specialization does not remove it generically. Thus there are three
different objects:

| Object | Degree over t | Genus | What rational points supply |
|---|---:|---:|---|
| C', product-square quotient | 2 | 1 | A candidate common squareclass, still subject to isogeny lifting |
| C, marked two-direction carrier | 4 | 1 | The two chosen rational native points |
| C3, all three displayed lifts | 8 | 5 | All three native points, without imposing their relation |

The original group has marked generic rank 17. The two native characters
give two independent generic directions over Q(C), so its displayed generic
subgroup has rank 19. This is not the rank of the auxiliary Jacobian and
does not assert the full generic pullback rank. On the specific fibre,
the [exact witness replay](PAIRED_SOLUBILITY_AND_SPECIALIZATION_COLLAPSE.md)
proves that B and D are independent modulo the generic subgroup. The three
points above have quotient rank exactly two. The fourth successful cover
`13109` contributes another independent direction, leaving the previously
certified quartet quotient rank three and four retained +7 directions
unexplained.

## Why the intersection has degree twelve

Let \(\tau_A,\tau_B,\tau_D\) be the native traces; each has height ten.
Map the genus-one pair carrier C to the original K3 by

\[
Q=S+P_B-P_D.
\]

The image D has degree four over t. To see that the map is birational to
its image, consider the three nonidentity elements of the biquadratic
Galois group. A single sign change changes a nonzero anti-invariant
section; a double sign change could fix Q only if the two different
nonzero characters cancelled. Their independence rules this out. Hence
the generic stabilizer is trivial.

Put

\[
z=2S+\tau_B-\tau_D,\qquad R=z-\tau_A.
\]

Over C the old generic heights multiply by four. Each native anti-invariant
section has height 24 there, by the previous quadratic-cover height-12
calculation and the further degree-two base change. Character orthogonality
therefore gives

\[
\widehat h_C(Q)=h(z)+12.
\]

The degree-four pullback has Euler characteristic eight. The branch
divisors are disjoint from the old nodal fibres, so no reducible-fibre
correction is introduced. Thus \(Q\cdot O=(h(z)-4)/2\). The image class is

\[
D=\left(\frac{h(z)}2+2,\,4,\,2z\right)
\quad\text{in }U\oplus M(-1).
\]

Consequently

\[
D^2=16,\quad p_a(D)=9,\qquad
B_{\tau_A}\cdot D=h(R)+2.
\]

The exact lattice calculation gives \(h(S)=10\), \(h(R)=10\), and
\(D\cdot O=6\). Hence the intersection degree is **12**. The image has
normalization genus one but arithmetic genus nine; its singularities have
total delta-invariant eight. Counting only its normalization genus would
miss this contribution. More generally, in this construction the
intersection number is even, so a degree-one explanation is impossible.

## The condition on t from exact elimination

The [frozen triple protocol](NATIVE_TRIPLE_INTERSECTION_PROTOCOL.json)
forms Q over \(\mathbf Q(t)(\sqrt f,\sqrt g)\) and imposes A's residual
chord. Its equation is linear in the second quadratic root. Substituting
that expression in its square equation gives an equation linear in the
first root. Eliminating this root and imposing A's residual x-quadratic
produces the saturated intersection polynomial

\[
F_S(t)=\left(t+\frac{288}{65}\right)H_{11}(t).
\]

The full coefficients and the rational functions giving all three roots
are in the immutable certificate. The initial norm numerator has degree
72; the residual-quadratic gcd has degree 12. The final polynomial is
squarefree of degree 12 and coprime to every recorded branch, singular-fibre,
and denominator exclusion.

Independent exact arithmetic verifies all three square equations and
\(P_A-P_B+P_D=S\) over \(\mathbf Q[t]/(F_S)\). Thus the certificate
exhibits twelve distinct geometric intersection points. Their count equals
the independent intersection number, so it exhausts the proper
intersection, including a check that none remain at excluded places or
infinity. Reduction modulo 73 proves \(H_{11}\) irreducible. Therefore the
finite Q-scheme is

\[
\operatorname{Spec}\mathbf Q\;\sqcup\;
\operatorname{Spec}\bigl(\mathbf Q[t]/(H_{11})\bigr).
\]

For this fixed generic relation, rational solubility is equivalent to
\(t=-288/65\). At this t, the primitive third root is
\(a=848451138/65\), alongside u0 and v0 above. No new parameters were
searched or evaluated; the only cohort hit among the 32 frozen values is
the already known `08234-003`.

The stronger condition \(F_S(t)=0\) is **sufficient** for C to have the
required rational fibre and for the third native point to satisfy this
particular relation. It is **not necessary** for C itself to have a
rational fibre: C has infinitely many rational points. Nor does it
classify other rational points of C3 with a different group relation.

## What this explains and what it does not

The successful chain now has the form

\[
\text{rational factor of the translated intersection scheme}
\Longrightarrow\text{three simultaneous rational lifts}
\overset{\text{independent certificate}}{\Longrightarrow}
\text{a two-dimensional quotient block}.
\]

The rational factor is a **solubility** event; the final dimension check
is **incidence**. No point-search visibility input enters this derivation.
The pair carrier and its isogeny give a precise smaller lifting problem
before the genus-five condition is imposed.

The generic translate S was selected using the retained exceptional
relation. Although elimination uses only generic equations thereafter,
this is not a prospective explanation of why that rational factor should
exist. The missing implication is a pre-point arithmetic criterion that
forces such a rational component while retaining at least two independent
directions. The norm-six pair criterion from the previous note cannot
supply it: this three-cover intersection has degree twelve and no
degree-one guarantee.

The next bounded comparison should keep this same triple and intersection
degree, choose alternative generic translates by the trace lattice alone,
and examine the rational-factor/Galois obstruction without searching new
specializations. That would test whether the linear component reflects
additional generic structure or the oracle choice of S. For Agent 1, the
current usable mathematical statements are an exact isogeny-lifting test
and a certified block at one retained parameter, not a new rank selector.

## Reproduction and evidence

Both workers had a 60-second bound. The auxiliary calculation used Sage
10.9 and PARI 2.17.3, effort zero and a 64 MiB initial PARI stack. From
the repository root:

```sh
sage -python elliptic-curves/rank-jump/native_triple_intersection.py check
sage -python elliptic-curves/rank-jump/verify_native_triple_carrier.py check
```

The verifier uses a separate explicit group law and modular irreducibility
check. The auxiliary rank is replayed with the same PARI 2-descent, not a
second independent descent implementation.

Immutable evidence: [triple input](../../artifacts/generated-results/elliptic-curves/rank_jump_native_triple_intersection_inputs_v1.json),
[intersection and root maps](../../artifacts/generated-results/elliptic-curves/rank_jump_native_triple_intersection_v1.json),
[minimal carrier input](../../artifacts/generated-results/elliptic-curves/rank_jump_minimal_native_block_carrier_inputs_v1.json),
[auxiliary geometry and descent](../../artifacts/generated-results/elliptic-curves/rank_jump_minimal_native_block_carrier_v1.json),
and [joint verification](../../artifacts/generated-results/elliptic-curves/rank_jump_native_triple_carrier_verification_v1.json).
All are hash-bound. Active search protocols, populations, scripts, and
mathematical status entries were not changed.
