# The four rational cubic directions on the marked genus-eleven cover

The four norm functions in
[MO3](QUADRATIC_MIRROR_CLASS_RANK_FOUR_FROM_MARKED_OCTIC.md)
exhaust the rational \(3\)-torsion of the genus-eleven
Jacobian on which they live. The involution of its even
equation gives two smaller Jacobians on which finite-field
group computations detect enough points of order nine.
A separate Prym calculation works for every nonzero scaling
of the cubic base change and expresses the rational
\(3\)-torsion in terms of divisibility of ramification-point
differences on the original Jacobian.

Both the full-Jacobian and Prym assertions hold for
every nonzero rational scale. Chebotarev repeats the
two checked actions on three-torsion at primes away
from any fixed scale, removing restrictions at the
original certificate primes. In particular,
[BS2–BS4](MARKED_OCTIC_BETA_CUBIC_SCALE_CLASSIFICATION.md#the-allowable-cube-classes-and-the-cover-isomorphisms)
classify the scales allowed by virtuality of the fourth
mark. Changing scale supplies no fifth generic direction.
The four scales also allowing local cubes above three give two
inequivalent covers, each with rational \(3\)-rank four.

## A ramification divisor criterion for the cubic Prym

**Proposition CK1.** Let \(C\) be a smooth projective
geometrically connected curve over a number field
\(k\not\supset\mu_3\). Suppose that \(t\in k(C)\) has
nonempty simple zeros and poles, all individually rational.
List them as \(Q_1,\ldots,Q_r\) with orders
\(\varepsilon_i\in\{1,-1\}\). For any \(c\in k^\times\),
let \(\pi:C_c\to C\) be the normalization of
\(t=cu^3\), and let \(P_c=\ker(\pi_*:\operatorname{Jac}(C_c)\to J)\),
where \(J=\operatorname{Jac}(C)\).
Write \(\chi\) for the nontrivial character of \(k(\mu_3)/k\).
Assume \(J[3](\overline k)\) contains no line with character
\(\chi\). Define
\[
 H=\{(m_i)\in\mathbb F_3^r:\textstyle\sum_i m_i=0\}
             /\langle(\varepsilon_i)\rangle,\qquad
 \partial:H\longrightarrow J(k)/3J(k),
 \qquad
 (m_i)\longmapsto [\textstyle\sum_i\widetilde m_iQ_i],
\]
using integer lifts \(\widetilde m_i\) of sum zero.
Then \(P_c\) is connected and there is an exact sequence
\[
 0\longrightarrow J(k)[3]\xrightarrow{\pi^*}P_c(k)[3]
                \longrightarrow\ker\partial\longrightarrow0.
 \tag{1}
\]
In particular its rational \(3\)-rank is independent of \(c\).

Each branch point \(Q_i\) has one rational point \(Q_i'\)
above it with ramification index three. Hence pullback
on Jacobians is injective and its dual norm has connected
kernel by [Li–Mulase, Proposition 1.7](sources/li-mulase-prym-integrable-systems.pdf).
In this case the argument is also immediate from
\(\pi_*\pi^*=[3]\). If a nonzero \(3\)-torsion line
bundle pulled back to zero, its connected étale cubic
cover over \(\overline k\) would factor the map \(\pi\),
contradicting the total ramification.

Let \(\tau\) be the automorphism \(u\mapsto\zeta_3u\)
over \(k(\mu_3)\). Every \(\tau\)-invariant divisor class
over \(\overline k\) has a \(\tau\)-invariant divisor
representative. Indeed, scale a chosen isomorphism
\(\tau^*L\simeq L\) so that its threefold composite is the
identity. Hilbert's Theorem 90 on the generic fiber gives
a nonzero invariant rational section of this linearized
line bundle. The divisor of that section is invariant.
An invariant divisor is a sum of a pullback divisor and
integer multiples of the \(Q'_i\). For degree-zero classes
its residue coefficients at \(Q_i'\), modulo three, sum
to zero. Every vector of sum zero occurs using lifts
whose integer sum is zero. A principal invariant divisor
comes from an eigenfunction because \(\tau(f)/f\) is
constant. Every eigenfunction has the form \(u^j g\)
for \(g\in\overline k(C)^\times\) and \(j=0,1,2\).
The residues of its divisor are \(j(\varepsilon_i)\),
giving precisely the displayed relation. If the residues
of an invariant divisor vanish, it is a pullback. Therefore
\[
 \operatorname{Jac}(C_c)^\tau/\pi^*J\simeq H
\]
on geometric points, with trivial Galois action on \(H\).
Intersecting the numerator with \(P_c\) gives the extension
\[
 0\longrightarrow\pi^*J[3]\longrightarrow
       P_c[1-\tau]\longrightarrow H\longrightarrow0.
 \tag{2}
\]
For surjectivity one can subtract a pullback whose norm
matches the norm of a representative, since \([3]\) is
surjective on geometric points of \(J\).

To lift \(m\in H\) rationally through (2), start from the
rational representative \(\sum_i\widetilde m_iQ_i'\).
One can subtract \(\pi^*A\) to bring its norm to zero
exactly when \([3]A=[\sum_i\widetilde m_iQ_i]\).
For a rational lift, \(A\) is rational by injectivity of
\(\pi^*\). Changing the lifts adds a multiple of three
on \(J\); changing residues by \((\varepsilon_i)\) adds
the principal divisor of \(t\) modulo such a multiple.
This gives the rational image \(\ker\partial\).

On \(P_c\), \(1+\tau+\tau^2=0\), so
\(\kappa=\tau-\tau^2=\tau(1-\tau)\) satisfies
\(\kappa^2=-3\), and \(P_c[1-\tau]\subset P_c[3]\).
For \(v\in P_c(k)[3]\), \(\kappa v\) belongs to
\(P_c[1-\tau]\) and transforms by \(\chi\).
By (2) its image in \(H\) is zero, and the hypothesis
on \(J[3]\) forces \(\kappa v=0\).
Thus every rational \(3\)-torsion point of \(P_c\) is
already in the fixed kernel of (2), proving (1).

## Exact rational torsion on the marked cover

**Theorem CK2.** Put
\[
 F(t)=49t^8+750t^6-2403t^4+27000t^2+63504
\]
and let \(C_c:R^2=F(cu^3)\) for \(c\in\mathbb Q^\times\).
For every \(c\) its norm kernel \(P_c\) over
\(C:R^2=F(t)\) has
\[
 \operatorname{Jac}(C_c)(\mathbb Q)[3]=P_c(\mathbb Q)[3]
    =\pi^*J(\mathbb Q)[3]\oplus
       \langle[2B_{3,c}]\rangle
    \simeq(\mathbb Z/3)^4 ,
\]
using the boundary divisor \(B_{3,c}\) of
[MO5](QUADRATIC_MIRROR_CLASS_RANK_FOUR_FROM_MARKED_OCTIC.md).
There are no local restrictions on \(c\).

To apply CK1, the four ramification points are
\(P_0^\pm,P_\infty^\pm\) in MO5.
The degree-\(2^9\) isogeny of MO1 and
[HG5](GENUS_THREE_THREE_MARKED_ELLIPTIC_QUOTIENTS.md#exact-transfer-of-three-primary-descent)
identifies \(J[3]\) with the \(3\)-torsion of the three
elliptic factors. Their rational marked lines are trivial
Galois modules, and the irreducible residual cubic
\(3\)-division factors exclude any \(\chi\)-line, since a
\(\chi\)-eigenpoint has a rational \(X\)-coordinate.
Hence \(J(\mathbb Q)[3]\simeq(\mathbb Z/3)^3\).

In CK1 the group \(H\) has dimension two, with basis
given by the two signed differences at \(0\) and \(\infty\).
Their sum lifts to \(B_0\), whose class has order two by
MO5 and is therefore divisible by three. The difference
at \(0\) is not divisible by three. Its pushforward
to the third elliptic factor is \(-2P_0\) for
\(P_0=(75/121,1750/1331)\). At the good prime \(17\),
\(E_3(\mathbb F_{17})\) has order \(18\) and the
reduction of \(-2P_0\) has order nine. It cannot be
three times a point of this group. Thus \(\ker\partial\)
has dimension one. The class \([2B_{3,c}]\) maps to
its nonzero direction and gives the displayed basis.

For the full Jacobian, put
\[
   G(s)=49s^4+750s^3-2403s^2+27000s+63504,\qquad
   f_c(v)=G(c^2v^3).
\]
The involutions \((u,R)\mapsto(-u,R)\) and
\((u,R)\mapsto(-u,-R)\) have respective quotient maps
\((v,Z)=(u^2,R)\) and \((v,Z)=(u^2,uR)\), giving
\[
 A_c:Z^2=f_c(v)\quad(g=5),\qquad
 B_c:Z^2=vf_c(v)\quad(g=6).
\]
The pullbacks give the plus and minus eigenspaces
for the first involution, since the hyperelliptic
involution acts as \(-1\) on a Jacobian. Off-diagonal
push-pull maps are zero by these signs. Diagonal
push-pull is multiplication by two, so the sum of
pullbacks is an isogeny with kernel killed by two. It
identifies the rational \(3\)-torsion with the direct
sum of that of the two factors.

First fix \(c_0=3/7\). The following good reductions
of its quotient models give bounds on their rational
three-torsion. Prime-to-\(p\) torsion specializes
injectively, and reduction identifies the geometric
three-torsion with that of the special fiber, respecting
the unramified Frobenius action. For the reduction
injection see
[Müller–Reitsma, Section 1.1](sources/muller-reitsma-genus-three-torsion-2211.03372v4.pdf).

| factor | \(p\) | Jacobian order | affine points giving order-nine witnesses |
|---|---:|---:|---|
| \(A_{3/7}\) | 83 | \(27\cdot169907840\) | \((0,30)\) in the odd model below |
| \(B_{3/7}\) | 137 | \(81\cdot90130506752\) | \((1,7),(3,24)\) |

In the first row, \(v=62\) is a simple root of \(f_{3/7}\)
modulo \(83\). The change \(v=62+1/x\) gives the degree
eleven model with right side \(x^{12}f_{3/7}(62+1/x)\).
In the second row the degree thirteen equation is
\(Z^2=77v^{13}+69v^{10}+88v^7+16v^4+73v\bmod137\).
Take the divisor of each affine point minus the unique
point at infinity and multiply by the prime-to-three
cofactor shown. All resulting points have order nine.
In the second row the triples of the two points are
independent, since their reduced Mumford first
coordinates \(x^6+68x^3+8\) and \(x^6+135x^3+128\)
are different. A group of order with \(3\)-part \(27\)
and an element of order nine has \(3\)-rank at most
two. In the second row the independent triples imply
that the two order-nine points generate a subgroup
of order \(81\), exhausting the \(3\)-primary subgroup.
Thus the three-torsion fixed by Frobenius at \(83\)
on \(\operatorname{Jac}(A_{c_0})[3]\), and at \(137\)
on \(\operatorname{Jac}(B_{c_0})[3]\), has dimension
at most two in each case.

To transport each bound, let \((J_0,p_0)\) be either
\((\operatorname{Jac}(A_{c_0}),83)\) or
\((\operatorname{Jac}(B_{c_0}),137)\), and put
\(L=\mathbb Q(J_0[3],\mu_3)\). This is a finite Galois
extension unramified at \(p_0\): multiplication by three
on the abelian scheme of good reduction has finite
étale kernel, since \(p_0\ne3\). Its geometric fibers
identify the three-torsion and its Frobenius action.
Let \(\gamma\) be its
Frobenius conjugacy class at \(p_0\). By
[Lenstra, Theorem 3.1](sources/lenstra-chebotarev-density-notes.pdf),
infinitely many primes have Frobenius class \(\gamma\).
Its restriction to \(\mathbb Q(\mu_3)\) is nontrivial,
so all these primes satisfy \(p\equiv2\bmod3\).
Conjugate actions on \(J_0[3]\) have the same fixed-space
dimension. At every such good prime the reduced
Jacobian therefore has three-torsion rank at most two.
Only this mod-three action is transported; the full
group order need not agree with the certificate row.

For any fixed \(c\in\mathbb Q^\times\), choose such a
prime away from \(3,7\), the bad reduction primes and
the numerator and denominator of \(c\). Then
\(c/c_0\) is a \(p\)-adic unit. Cubing is bijective
on \(\mathbb F_p^\times\), and Hensel's lemma lifts
each cube root uniquely to a unit of \(\mathbb Q_p\).
Write \(d^3=c/c_0\). The map
\((u,R)\mapsto(du,R)\) identifies \(C_c\) with
\(C_{c_0}\) over \(\mathbb Q_p\), compatibly with
both involutions. On the quotients it gives
\((v,Z)\mapsto(d^2v,Z)\) for \(A_c\), and
\((v,Z)\mapsto(d^2v,dZ)\) for \(B_c\).
Consequently the corresponding Jacobian for \(c\)
has good reduction and rational three-torsion rank at
most two. The primes for the two factors may be
chosen separately. Their direct sum bounds the full
rational three-torsion rank by four for every \(c\).
The Prym already supplies four independent points,
giving equality.

The [exact checker](verify_marked_cubic_base_change_rational_torsion.py)
uses SageMath 10.9 to check squarefreeness, recompute
the two Frobenius polynomials, and verify the order-nine
divisors and the elliptic boundary reduction. The proof
uses the finite group structures as well as the orders.
An internal replay with SageMath 10.9 verified these
certificates unchanged for the unrestricted-scale
deduction. Chebotarev transport is the written argument
above, not a computation of torsion fields or an
independent mathematical review.

## Consequence for cubic norm functions

**Corollary CK3.** For every \(c\in\mathbb Q^\times\),
let \(a\in\mathbb Q(C_c)^\times\) have all divisor
multiplicities divisible by three and quadratic norm
\(a\bar a\in\mathbb Q(u)^{\times3}\).
Modulo cubes in \(\mathbb Q(C_c)^\times\), \(a\) lies
in the span of the three pullbacks of \(\alpha_i\) of MO1
and the pullback of \(\beta\) of MO3.

The divided divisor gives a rational point of
\(\operatorname{Jac}(C_c)[3]\). The kernel of this map
on such function classes is trivial. Indeed, a principal
divided divisor can be represented by a rational function
by descent of principal divisors, leaving a constant
\(d\in\mathbb Q^\times\). The norm of this constant is
\(d^2\), which is a cube in \(\mathbb Q(u)\) only when
\(d\) is a rational cube. The four given functions have
cube norms, and their divided divisors are the basis
in CK2.

Consequently every such function gives a class already
in the four marked directions at every specialization
where all functions in the cube identity are defined
and nonzero. The result concerns functions with
divisors divisible by three on this cover. The
arithmetic routes in
[MO7](QUADRATIC_MIRROR_CLASS_RANK_FOUR_FROM_MARKED_OCTIC.md#the-real-unit-line-is-the-remaining-arithmetic-loss)
remain open, including proving that the marked space
avoids the unit line on an infinite subfamily.

## Virtuality on a positive proportion of parameters

**Proposition CK4.** Let \(C/\mathbb Q\) be the smooth
projective model of \(Y^2=f(u)\), where \(f\) is a
nonconstant squarefree rational polynomial. Let
\(a\in\mathbb Q(C)^\times\) and \(m>1\). For a rational
parameter \(s\) giving a quadratic field
\(K_s=\mathbb Q(\sqrt{f(s)})\), write \(a_s\) for either
conjugate evaluation when it is defined and nonzero.
If \(\operatorname{div}(a)\) is not a multiple of
\(m\), then
\[
 \#\{s\in\mathbb Q:H(s)\le H,\ [K_s:\mathbb Q]=2,\
                 v_{\mathfrak p}(a_s)\equiv0\bmod m
                 \text{ for all finite }\mathfrak p\}
      =o(H^2),
 \tag{3}
\]
where \(H(A/B)=\max(|A|,B)\) for coprime \(A,B\) and
\(B>0\). The same conclusion holds if one ignores any
fixed finite set of primes. Conversely if the divisor
is a multiple of \(m\), all evaluations have the stated
divisibility outside a fixed finite set.

The specialization-inertia theorem of
[Conrad, Theorem 0.2](sources/conrad-inertia-groups-fibers-1999.pdf)
relates geometric ramification of a cover to its
intersection multiplicity in a specialization. For
this assertion local parameters give the required
valuation directly, and allow the excluded neighborhoods
to be counted in a rational height box.

Choose a geometric point \(Q\) where
\(n=\operatorname{ord}_{Q}(a)\) is nonzero modulo \(m\).
Let \(s_0\) be its image on the parameter line, and
choose a finite Galois number field \(L\) over which
\(Q,s_0\) are rational. At \(Q\) choose a uniformizer
\(z\) and write \(w=z^e U,\ a=z^n V\), where \(w=u-s_0\)
for finite \(s_0\), \(w=1/u\) otherwise, \(e=1\) or
\(2\) is the ramification index of the hyperelliptic
map, and \(U,V\) are units at \(Q\). These equations
extend to a neighborhood of the section \(Q\) in a
smooth model after excluding finitely many primes of
\(L\). In particular the reductions of \(U,V\) are
nonzero at that section, the reduced hyperelliptic map
has the same index there, and the other point over
\(s_0\), if any, has distinct reduction.

For every rational prime \(p\) that splits completely
in \(L\) outside the excluded set, choose an embedding
\(L\to\mathbb Q_p\). Consider the parameters for which
\(v_p(w)=1\). If \(e=1\), one point over the parameter
reduces to \(Q\), and is defined over \(\mathbb Q_p\)
by Hensel's lemma. If \(e=2\), the point over the
parameter is defined over a ramified quadratic extension
of \(\mathbb Q_p\). Indeed in a coordinate at the
simple branch point the quadratic equation, up to a
square factor, has right side \(w\) times a \(p\)-adic
unit. This works at infinity using \(1/u\). In both
cases at the corresponding prime \(\mathfrak p\) of
a quadratic \(K_s\) the normalized valuation of \(w\)
is \(e\). The identity \(w=z^eU\) shows that the
valuation of \(z\) is one, and \(v_{\mathfrak p}(a_s)=n\).
Conjugating the evaluation simply permutes the two
places in the split case. Thus the depth-one parameter
neighborhood is forbidden for virtuality.

Modulo \(p^2\) there are \(p-1\) such projective
parameter classes out of \(p^2+p\) classes. For
integer pairs \((A,B)\) in the height box, each class
has \(\varphi(p^2)\) distinct lifts with at least one
coordinate a unit, so the forbidden fraction among
all pairs is \((p-1)^2/p^3\). Every pair counted in
(3) avoids all such classes, at every selected prime.
For any fixed finite subset \(\mathcal P\) of them,
the Chinese remainder theorem and counting each
lattice residue class give
\[
 \limsup_{H\to\infty}\frac{1}{2H^2}
   \#\{s\text{ counted in (3)}\}
 \ \le\ \prod_{p\in\mathcal P}
              \left(1-\frac{(p-1)^2}{p^3}\right).
\]
Chebotarev's density theorem gives a positive density
of completely split primes in the fixed field \(L\),
so their reciprocal sum diverges. The product tends
to zero as \(\mathcal P\) grows. The finite exceptional
set can be enlarged at will.

For the converse, extend the smooth curve and the
equality \(\operatorname{div}(a)=mD\) to a smooth proper
model over \(\mathbb Z[1/S]\) after enlarging \(S\).
The closure of \(D\) is Cartier on this regular model,
and only finitely many vertical components need be
removed from the equality. An evaluated point extends
over each valuation ring away from \(S\) by properness.
Locally \(a\) is a unit times the \(m\)-th power of
an equation for the Cartier divisor. Pulling back gives
the required divisibility whenever \(a_s\) is defined
and nonzero.

**Corollary CK5.** For every fixed \(c\in\mathbb Q^\times\),
a fixed function \(a\in\mathbb Q(C_c)^\times\) can
supply a fifth independent virtual cubic class at only
\(o(H^2)\) rational parameters of height at most \(H\).
In particular, no finite list of fixed rational functions
supplies a fifth class throughout a nonempty admissible
open slice defined at finitely many places.

If the divisor of \(a\) is not divisible by three,
apply CK4. If it is divisible by three, taking the norm
of the divisor shows that \(a\bar a=d\,b(u)^3\) for a
constant \(d\in\mathbb Q^\times\) and
\(b(u)\in\mathbb Q(u)^\times\). At any virtual evaluation
the rational norm is a cube, hence \(d\) is a cube.
Now CK3 places the class in the marked span wherever
the identity can be evaluated. The excluded parameter
values form a finite set. A finite union of the
zero-density parameter sets still has density zero.

For the last assertion, a nonempty open slice at
finitely many places contains one given by fixed
congruence classes on primitive pairs and an open
real interval in the ratio \(A/B\). If needed, one can
also force the specialization to be quadratic by
imposing a valuation-one condition on \(f(s)\) near a
simple root at one new split prime. The Chinese remainder
theorem and primitive lattice point count give
\(\gg H^2\) parameters in the resulting slice.
The bound in (3) is in parameter height. It leaves open
which virtual classes can occur on thinner parameter
sets or after selecting new functions separately for
each field. The sharper
[VB4 bound](QUADRATIC_VIRTUALITY_BINARY_FORM_POWER_OBSTRUCTIONS.md#applications-to-the-existing-field-counts)
is \(O(H^{3/2}\log H)\) for any fixed finite list of
functions, giving \(O(X^{1/16}\log X)\) exceptional
fields through the height-bounded witnesses in MO4.
It gives no upper count for extra classes obtained
from witnesses of unrestricted height.
