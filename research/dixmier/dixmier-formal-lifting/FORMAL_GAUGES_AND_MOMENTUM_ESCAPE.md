# Formal fixed-R gauges: finite momentum and the remaining escape

14 November 2026. **QFG1.** The comparison gap after
[QPL1](PAINLEVE_GLOBAL_LONG_LIFT.md) can be narrowed substantially.

Every identity-mod-h formal Weyl automorphism fixing Long's R,
whose two position images have uniformly finite momentum order,
is a formal cotangent coordinate change followed by a positional
gradient shift. Its action on the R centralizer preserves every
finite momentum order, even when its position degrees are unbounded.

Consequently such a gauge cannot carry the explicit Painleve
lift into a spatially bounded full lift. If a bounded exact-R
Long lift exists, its unique target formal comparison with QPL1
must send at least one original position to a series of unbounded
momentum orders.

This is a necessary condition, not a nonlifting theorem. There
are explicit global formal gauges with precisely that unbounded
momentum behavior. In fact they can fix both R and QPL1's exact
Euler connection K, at arbitrarily high identity contact. Thus
preserving the scaling connection alone does not force the
comparison into the finite-momentum category.

## The completed algebra and the finite-momentum condition

Use the normalization of [QGR1](BOUNDED_GAUGE_RIGIDITY.md):
\[
 [p,x]=[z,y]=h,\quad
 r=R/2=x(1+xy),\quad u=1+xy,\quad
 b=x^2p-(1+2xy)z,
\]
and
\[
 C_h=K[h,x,y][b;h\theta],\qquad
 \theta=x^2\partial_x-(1+2xy)\partial_y,\quad\theta(r)=0.
                                                               \tag{1}
\]
Let Ahat be the h-adic completion of the full Weyl algebra,
and Chat its R centralizer. Every h coefficient is a finite
polynomial, but momentum orders may grow across coefficients.

Write O=K[x,y][[h]]. The finite-b-order part of Chat is
O[b;h theta]. Its b order is additive on nonzero products:
the associated order-graded ring is the domain O[b].
A uniform ambient momentum bound for an element of Chat
is equivalent to finite b order, since b is a nonzero linear
form in the original momenta and the coefficientwise
centralizer decomposition preserves momentum degree.

Let tau be an h-continuous formal endomorphism of Ahat,
or of Chat, with tau=id modulo h and tau(r)=r. Such an
endomorphism is automatically a formal automorphism by the
convergent inverse of id+O(h). Assume only that tau(x) and
tau(y) have finite b order. No spatial bound is assumed.

## Factorization forces the positions to remain positional

Put X=tau(x), Y=tau(y), U=tau(u)=1+XY. All have finite
b order and
\[
                              XU=r.                           \tag{2}
\]
Additivity forces X and U to have b order zero. Since
U-1 has nonzero classical reduction xy, additivity in
XY=U-1 forces Y to have order zero too. Thus X,Y lie in O.

Reduce (2) modulo x. Since U mod x=1+O(h) is a unit,
X is divisible by x coefficientwise. Write X=x a,
a=1+O(h). Then U=u/a. The relation U=1+XY gives
u-a=x a^2Y; reducing modulo x makes a=1+x f with f in hO.
We obtain the complete position formula
\[
 X=x(1+xf),\qquad
 Y={y-f\over(1+xf)^2},\qquad f\in hO.                          \tag{3}
\]
Every such formula preserves r and is an identity-mod-h
formal coordinate automorphism of O. Its inverse is obtained
coefficientwise by formal composition. Denominators in
(3) are formal units; their individual h coefficients are
polynomials.

Its ordinary position Jacobian is the formal unit
\[
 j_f=\det{\partial(X,Y)\over\partial(x,y)}
       ={1+2xf+\theta(f)\over(1+xf)^2}=1+O(h).                \tag{4}
\]
Equivalently theta(X)=j_f X^2 and
theta(Y)=-j_f(1+2XY). Formula (3) allows genuinely moving
positions. For example f=h already has unbounded position
degrees on expansion, unlike QGR1's finite-spatial group.

## Complete centralizer classification and admission

Let L_f be the formal cotangent lift of (3), with all
coefficients on the left in the momentum images. If M is
the two-by-two Jacobian matrix of (3), its i-th momentum
image is sum_j (M^-1)_(ji)*p_j. The dual vector fields
commute, so these are exact Weyl operators, and
\[
                              L_f(b)=j_f^{-1}b.
\]
For a centralizer automorphism tau as above, compose with
L_f^-1 to fix the positions. The difference between its
b image and b then commutes with both x,y, and is positional.
This proves the full classification
\[
       \tau(x)=X,\quad\tau(y)=Y,\quad
                    \tau(b)=j_f^{-1}b+F,\qquad F\in hO.      \tag{5}
\]
No momentum-bound assumption on tau(b) is needed once the
two position images have finite order. Every (f,F) gives
such an automorphism of Chat.

For a full Weyl automorphism, position-fixing maps are
exactly the gradient shifts gamma_G, with
p -> p+G_x, z -> z+G_y, G in hO modulo scalar series.
This follows directly from the joint position centralizer
and polynomial Poincare lemma coefficientwise. Compare
[Milicic, Chapter 1, Lemmas 1.1--1.2](https://www.math.utah.edu/~milicic/Eprints/book.pdf);
the formal assertions here also have this direct proof.
Thus every full tau in this category has the form
\[
                \tau=\gamma_G\circ L_f,\qquad
                \tau(b)=j_f^{-1}(b+\theta G).                \tag{6}
\]
In particular all its momentum-generator images, and those
of its inverse, have momentum order at most one.

A centralizer map (5) extends to a full formal Weyl
automorphism exactly when
\[
                             j_f F\in\theta(O).               \tag{7}
\]
Indeed it is gamma_(j_f F) composed with L_f at the
centralizer level; a full gradient realizes the translation
precisely when theta G=j_f F. By
[QRP1](RELATIVE_PERIOD_TORSION.md), this is the full global
relative-cohomology condition coefficient by coefficient.
Generic fiber exactness still discards its exceptional torsion.

The same criterion decides whether applying (5) to the
fiber images of any admitted exact-R Long tuple admits
a polynomial-coefficient fourth image. Sufficiency is the
full extension (6). For necessity, after the admitted L_f
change, take the first coefficient of the translation
j_f F with a nonzero global class. Remove the finitely
many earlier exact coefficients by gradient shifts.
QRP1's first-difference full-frame argument then excludes
a fourth image at that order. Higher coefficients cannot
repair it. This uses no common spatial bound.

For instance the formal cotangent base f=h is admitted.
The maps with the same base and F=h^N/j_f have translation
class h^N*[1], which is globally nonzero by QRP1. They give
nonadmissible triples agreeing with an admitted one to
any prescribed finite order.

## Finite momentum is preserved on every centralizer element

If P=sum_(i=0)^m a_i b^i in O[b;h theta] is nonzero,
equation (5) gives leading coefficient
tau(a_m)*j_f^(-m), which is nonzero. Therefore
\[
                            \operatorname{ord}_b\tau(P)=m.    \tag{8}
\]
This holds even if f,F and the a_i have unbounded position
degrees. It is not a statement about arbitrary formal
gauges whose position images already have unbounded
momentum orders.

Apply (8) to QPL1: its fiber orders remain exactly two
and three under every full gauge in this category.
If the resulting full tuple had a common spatial bound,
its generic scalar coefficient field would be algebraic,
and QSO2's quadratic-mate theorem would exclude that pair.
Hence none of these gauges gives a bounded full tuple.

Allowing a finite regular source Weyl automorphism with a
finite inverse does not change the conclusion: if the
whole source-changed tuple were bounded, composing back
with that finite inverse would bound the original tuple.

## What a hypothetical bounded comparison must do

For completeness, the identity-mod-h target comparison
between two full polynomial-coefficient formal lifts of
the same classical symplectic map is unique. At a first
difference, invert the classical Jacobian polynomially
to obtain its ambient vector field. The linearized full
relations make it symplectic, and polynomial H1=0 makes
it Hamiltonian. Successively removing these differences
gives the comparison, as in [QFL1](PROOF.md).

If a nonidentity target automorphism fixed all the first
tuple's images, its first nonzero vector field would
annihilate every classical image coordinate. The unit
Jacobian forces that vector field to vanish, a contradiction.
This proves uniqueness of the target automorphism, modulo
the irrelevant scalar ambiguity of its Hamiltonian generators.

If a bounded exact-R lift Phi_b existed, its comparison
tau with the PI lift would satisfy tau(R)=R. Equation (8)
and QSO2 show that tau(x),tau(y) cannot both have finite
momentum order:
\[
 \Phi_b=\tau\circ\Phi_{\rm PI}\quad\Longrightarrow\quad
 \sup_n\operatorname{ord}_{p,z}[h^n]\tau(x)=\infty
 \quad\hbox{or}\quad
 \sup_n\operatorname{ord}_{p,z}[h^n]\tau(y)=\infty.              \tag{9}
\]
Thus the missing comparison is more than an unbounded
position change. It requires unbounded position-momentum
mixing in the image of an original position coordinate.
No bounded exact-R lift is proved impossible by (9).

## Such mixing can preserve the exact scaling connection

The category excluded in (9) cannot simply be imposed
from R and the normalizer K. Here is a global all-orders
control, with every full relation intact.

Use QPL1's original uppercase positions X,Q, operators
R,B,w and its exact rational Euler connection K. Define
\[
 I={X^3w+wX^3\over2},\quad H=I^2,\qquad
 \tau_N=\exp\left(h^N{\operatorname{ad}H\over h}\right),
                         \quad N\ge1.                        \tag{10}
\]
H is a fixed global polynomial operator of momentum order
two, so each h coefficient of tau_N is polynomial and
the exponential is a full formal Weyl automorphism.
It is identity modulo h^N and fixes R.

QPL1 gives [K,X]=-hX/(5R) and [K,w]=3hw/(5R).
The opposite weights in X^3w prove [K,I]=[K,H]=0.
Consequently tau_N fixes K exactly as well. Both its
forward and inverse K images are the same finite
rational operator.

Nevertheless its position images have unbounded momentum.
The first-order operator I acts on X by
\[
                            [I,X]/h=-3RX^6.
\]
Let theta_I=-3RX^6*partial_X on the position subring K[R,X].
The leading I-order part of the n-th normalized H commutator
on X is 2^n*theta_I^n(X)*I^n, and
\[
 \theta_I^n(X)=(-3R)^n
       \left(\prod_{j=0}^{n-1}(5j+1)\right)X^{5n+1}.
\]
Since the leading B coefficient of I is -RX^4,
\[
 \operatorname{top}_B
       ((\operatorname{ad}H/h)^nX)
   =6^n R^{2n}\left(\prod_{j=0}^{n-1}(5j+1)\right)
                           X^{9n+1}B^n\ne0.                 \tag{11}
\]
In the h^(Nn) coefficient of tau_N(X), no later iterate
can occur and no earlier one reaches B order n.
Thus (11), divided by n!, cannot cancel. The inverse
flow has the same unbounded-order property.

These flows also produce full Long lifts with R and K
unchanged but unbounded momentum orders. For S_PI, the
classical leading I coefficient is -X^-6/(3R) times I^2.
Repeated theta_I differentiation has factors -6+5j,
none zero at integral j. The h^(Nn) coefficient after
transport therefore has nonzero B order n+2. Earlier
iterates and the positional positive PI coefficients
cannot cancel that order.

This is not a bounded candidate. It disproves the inference
that exact R, exact K, every full relation and arbitrary
identity contact force a finite-momentum comparison or
the original quadratic/cubic PI form. Additional hypotheses
on an actual bounded target tuple could still matter.

## Research consequence and verification boundary

Within uniformly finite momentum, the formal fixed-R gauge
and its global admission criterion are now explicit. A
bounded candidate, if it exists, must use the genuinely
unbounded-momentum comparison left outside this category.
Exact preservation of the Euler connection does not
eliminate that escape.

[Exact controls](check_formal_gauges.py) verify the arbitrary
base Jacobian identities, a full cotangent/gradient example,
the exact R,K invariance of H and selected normalized
commutator coefficients. The classification, coefficientwise
admission, comparison uniqueness and all-n growth use written
proofs. The finite samples are controls, not an all-orders
search or independent theorem replay.

No bounded full Long lift, all-representative nonlifting,
DC2 solution, formal verification, external review or
novelty priority is claimed.
