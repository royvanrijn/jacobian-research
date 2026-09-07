# Explicit projection fibres isolate the remaining solubility obstruction

The common Jacobian construction now has explicit equations for its
projection fibres, Mumford divisor representatives, and a degree-four map
to the specialized elliptic curve. Six independently constructed quartics
at \(u=-1\) have exactly the retained Kummer labels. Three CT comparisons
reproduce \(0,1,0\), with 55 verified local Hilbert-symbol witnesses.

The main structural consequence is stronger than this small check:
**the common conic already splits for every inherited Selmer class**,
by the Hasse principle for conics. Its rational parametrization cannot
distinguish a rational Mordell–Weil class from Sha. The remaining
genus-one double cover carries that distinction.

This is a **solubility** analysis of the fixed-cubic control family,
whose arithmetic generic rank is zero. It is not a rank predictor for
MW17/MW16, and it uses public anchor points only retrospectively.
The [initial comparison panel and paired studies](ANALYSIS.md) remain the
place for high-gain versus low-gain evidence in those families.

## Equations for a rational divisor-class lift

Use the notation and labelled maps of
[the common-Jacobian construction](JACOBIAN_GLUING_AND_SHA_BLOCKS.md):
\[
E_0:b^2=f(a)=a^3+Aa+B,\quad D=1+Au^2+Bu^3,\quad uD\ne0,
\]
\[
C_u:\quad Y^2=g_u(z)
=uD-(3u+Au^3)z^2+3uz^4-uz^6,\qquad J_u=\operatorname{Jac}(C_u).
\]
Write \(q_u=\pi_{0*}:J_u\to E_0\), and fix a non-torsion anchor point
\(P=(p,q)\). Here \(q\) denotes its ordinate, not the map \(q_u\).
Put \(\kappa(p)=1+up+u^2(A+p^2)\).

In projective coordinates \((s:t:v:w)\), the fibre \(q_u^{-1}(P)\)
has the following smooth genus-one model:
\[
\begin{aligned}
Q_1&=v^2+us^2-2utw-u(up+2)w^2=0,\\
Q_2&=ut^2+(1-up)v^2+2u^2qvw-u\kappa(p)w^2=0.
\end{aligned}
\]
More precisely, the following open divisor chart gives a birational
identification; its smooth projective completion is the projection fibre.
Set \(w=1\), and
\[
h(z)=z^2-sz+t,\qquad n=v(up-1-t)-u^2q.
\]
The Mumford representative is
\[
h(z)=0,\qquad Y\equiv svz+n\pmod {h(z)}.
\]
Subtract the rational degree-two divisor
\(D_\infty=\pi_0^*(O)\) to obtain a class in \(J_u\).
The two quadric equations imply
\(g_u(z)-(svz+n)^2\equiv0\pmod {h(z)}\).

For the norm to \(E_0\), the line through the two image points has
slope \(-v/u\). The sum of those points has coordinates
\[
x_0=\frac{v^2}{u^2}-\frac{2-s^2+2t}{u},\qquad
y_0=-\frac n{u^2}-\frac{v(1+t)}{u^2}+\frac{vx_0}{u}.
\]
Reduction by \(Q_1,Q_2\) gives \(x_0=p,\ y_0=q\).
Conversely these chord identities recover \(Q_1,Q_2\) on the open chart.
A general degree-two divisor class on a genus-two curve has a unique
effective representative, so this is a divisor-class fibre model.

There is no assumption here that a rational point exists on that fibre.
The displayed Mumford representative is the rational divisor class one
would obtain **if** the two quadrics had a rational point.

## The other norm and the exact 2-cover label

On the open set \(st\ne0\), put
\[
e=v(up-1)-u^2q,\qquad
\mu=\frac{vt^2+e(s^2-t)}{ust}.
\]
The other norm is
\[
\begin{aligned}
x_u&=\mu^2+\frac2u-\frac{D(s^2-2t)}{ut^2},\\
y_u&=\frac{Des}{u^2t^2}
-\mu\left(\frac D{ut}+Au+\frac1u+x_u\right).
\end{aligned}
\]
It satisfies
\[
y_u^2=x_u^3+2Au x_u^2+(A+3Bu+A^2u^2)x_u
       +B+ABu^2-B^2u^3.
\]
These rational expressions extend to the smooth projective curves.
They arise by adding the two images under \(\pi_u\);
\(\pi_{u*}(D_\infty)\) has sum \(O\).

Let \(2R=P\) over an algebraic closure and use the origin
\(\Phi(R,0)\) on the fibre. Its other norm is \(O\), and the restriction
of the other norm to its embedded \(E_u\) is multiplication by two.
Its 2-cover cocycle is therefore exactly the transported Kummer class
\[
\beta_P=p-\theta,\qquad \theta^3+A\theta+B=0.
\]
Thus the equations identify the labelled cover, as well as its underlying
genus-one torsor. For the six retained quartics, the certificate separately
checks this label through an exact cubic square-root identity.

## Why the common conic cannot close solubility

The second quadric is a conic in \((t:v:w)\). Its Gram determinant is
\[
\det\begin{pmatrix}
u&0&0\\
0&1-up&u^2q\\
0&u^2q&-u\kappa(p)
\end{pmatrix}
=-u^2D.
\]
Indeed
\((1-up)\kappa(p)+u^3q^2=D\).
The determinant is independent of the chosen anchor point.

If \(\beta_P\in W_u\), its genus-one cover has a point in every completion.
Projection to \((t:v:w)\) gives a point on this conic in every completion:
the three coordinates cannot all vanish on the intersection.
The Hasse–Minkowski theorem then gives a rational point on the conic.

Consequently every inherited Selmer class admits this first
parametrization, including the CT-certified Sha classes. At \(u=-1\),
\(\dim W_{-1}=18\), while CT bounds the rational-lift subspace by dimension
two. All those locally soluble classes pass the common-conic step.
This implication uses the theorem, not an enumeration of their conics.

For \(u\ne0\), rank reduction of this conic requires \(D=0\), where
the elliptic specialization is singular. The excluded \(u=0\) limit is
the previously identified stable splitting event at the anchor.
This determinant therefore supplies no new smooth-fibre rank-jump
condition in the six controls.

The unresolved step remains a square on the parametrized conic:
a rational point on its remaining double cover, equivalently
\[
P\in q_u(J_u(\mathbb Q)).
\]
Simultaneous rational lifts of \(k\) independent anchor Kummer classes
would give \(k\) independent classes in \(E_u(\mathbb Q)/2E_u(\mathbb Q)\).
There is no rational 2-torsion here, so they would force rank at least
\(k\). A common rational conic parametrization alone gives none of them.

## Two quartic models can have different cover labels

At \(u=-1\), the equations become
\[
s^2=v^2+2tw+(2-p)w^2,\qquad
t^2=(p+1)v^2+2qvw+(p^2-p+A+1)w^2.
\]
Parametrizing the second conic gives the certified quartic with label
\(\beta_P\).

There is another tempting construction. With \(w=1\), put
\(a=s-v,\ b=s+v,\ t=(ab+p-2)/2\).
Eliminating \(b\) gives \(y^2=h_p(a)\), where
\[
\begin{aligned}
h_p(a)={}&(p+1)a^4-4qa^3+(6p^2-6p+4A)a^2\\
&+4q(2-p)a+p^3-3p^2+4B-4A.
\end{aligned}
\]
The raw discriminant is \(4h_p\). In the convention where
\(I=a_2^2-3a_4,\ J=-2a_2^3+9a_2a_4-27a_6\) for the cubic of \(E_{-1}\),
the normalized quartic \(h_p/4\) has these invariants. Its cubic invariant is
\[
\beta_P\kappa,\qquad
\kappa=1-\theta+A+\theta^2.
\]
This differs from \(\beta_P\) by the affine class
\(\eta=D(1+\theta)=\kappa(1+\theta)^2\).
The difference is the Kummer class of the already known rational point
\[
Q_\eta=(A+1,\ A-B+1)\in E_{-1}(\mathbb Q).
\]
The genus-one torsor is unchanged, but the degree-two model has changed
its 2-cover label by a rational Kummer class. Exact identities verify
this for all 20 retained anchor generators.

In particular, three alternate quartics for \(P,P',P+P'\) have product
class \(\eta\), not the trivial squareclass required by the three-cover
Fisher formula. Blindly reusing their labels would be an error even
though all three quartics describe the intended underlying torsors.

## Bounded retrospective certificate and replay

The [original protocol](PROJECTION_FIBRE_PROTOCOL.json) selects the first
three retained \(W_{-1}\) basis masks \(1,6,10\), together with their sums
\(7,11,12\). It constructs six conic parametrizations and quartics,
retaining every change of variables and cubic label witness in
[the fibre certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_projection_fibres_v1.json).

The original three CT calls each reached their 60-second limit.
Those attempts remain recorded as UNKNOWN. The
[separate retry protocol](PROJECTION_PAIR_PROTOCOL.json) changes the
arithmetic construction: the required square root of the product of the
three cubic invariants is multiplied from the six already certified
roots and the repeated basis factors. Support factorization first divides
by the known discriminant primes. Each retry has the same 60-second cap.

The [pairing certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_projection_pairings_v1.json)
reproduces the retained submatrix
\[
\begin{pmatrix}0&0&1\\0&0&0\\1&0&0\end{pmatrix}.
\]
Its three pairs retain 17, 19 and 19 local terms. All 55 Hilbert symbols
also agree with PARI. This independently constructs the covers; it
reuses the repository's Fisher formula and is not an independent CT
algorithm. The restricted one-dimensional radical is forced by parity
and does not certify a rational class.

From the repository root:
```sh
sage -python elliptic-curves/rank-jump/projection_pairings.py verify
sage -python -m unittest discover -s elliptic-curves/rank-jump -p test_projection_fibres.py
```
Replay performs exact polynomial, model, class and retained-witness
checks without conic solving, integer factorization or witness search.
It also verifies generic symbolic Mumford and determinant identities,
tests the norm maps at five other small parameters, and rejects corrupted
quartics, labels, square roots and conic maps through the separate tests.
Construction used Sage 10.9 and PARI 2.17.3. Checkpoints and failed draft
logs stay under the two ignored `artifacts/local/rank-jump-projection-*`
directories. No point or parameter search was performed.

## Consequence for the mechanism ranking

1. **Solubility: a large rational projection image in one auxiliary
   Jacobian remains the strongest precise target.** The extension and its
   torsion/local conditions exist before exceptional points are supplied.
   Splitting at the anchor is a proved simultaneous-solubility event.
   A comparably useful event in MW17/MW16 remains unidentified.
2. **Incidence: the relative full-Selmer and ramification calculations
   remain valid necessary structure.** These equations do not enlarge
   those groups or convert their CT radicals into rational points.
3. **Weak explanation: a shared conic, its splitting, or an easy conic
   parametrization.** Every inherited Selmer class already passes that
   step. Small coefficients can help subsequent point-search visibility;
   they are not evidence for a larger rational image.

The exact missing implication is a usable specialization condition
forcing a large subspace into \(q_u(J_u(\mathbb Q))\), beyond local
lifting and CT compatibility. For Agent 1, this currently supplies a
mathematical exclusion from future selectors: do not score conic
splitting as evidence that Selmer classes are rational. No active
selection policy or search output has been changed.

