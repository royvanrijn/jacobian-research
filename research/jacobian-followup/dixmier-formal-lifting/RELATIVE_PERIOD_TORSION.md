# Global Hamiltonian corrections and the invisible Broughton torsion

14 November 2026. **QRP1.** Long's undeformed first coordinate has an
explicit, all-degree test for a polynomial vertical derivation to be
Hamiltonian. Generic-fiber periods see only one part of that test.
There is also an infinite tower supported at the exceptional value
R=0, with elements requiring arbitrarily high powers of R to clear.
The same obstruction holds exactly in the polynomial quantum
centralizer. It also produces nonadmissible quantum triples with
Long's classical reduction that agree with an admitted triple to
any prescribed finite order. This does not decide whether Long
admits a bounded full quantum lift.

## The actual centralizer has Broughton's base

Use [QCC1's polynomial coordinates](CUBIC_CORE_NONCOMPLETION.md):
\[
 R=x(2-3xq),\quad B=3x^2p+2(1-3xq)z,\quad
 \{B,-\}|_{K[x,q]}=3x^2\partial_x+(2-6xq)\partial_q.
\]
Put
\[
 y=-3q/2,\quad r=R/2=x(1+xy),\quad b=B/3,\quad u=1+xy.
\]
Then the classical centralizer of R is
\[
 C=K[x,y,b],\qquad \{b,g\}=\theta(g),\qquad
 \theta=x^2\partial_x-(1+2xy)\partial_y,\quad \theta(r)=0. \tag{1}
\]
Here K is any characteristic-zero field. The map r is smooth:
\[
 (1-2xy)r_x+4y^2r_y=1.                              \tag{2}
\]
For r nonzero its fibers have coordinate ring K[r,r^{-1}][x,x^{-1}],
whereas its zero fiber is the disjoint union x=0 and u=0.
Smoothness does not make this nonproper family a product.

For the normal-order quantum algebra with [p,x]=[z,q]=h, the same
polynomial R and B satisfy
\[
 [B,x]=3hx^2,\quad [B,q]=2h(1-3xq),\quad [B,R]=0.
\]
Their polynomial subalgebra is exactly the Ore algebra
\[
 C_h=K[h,x,y][b;h\theta],\qquad [b,g]=h\theta(g).      \tag{3}
\]
The PBW powers of B are independent: their highest momentum symbols
are powers of the nonzero linear form B. After inverting h, this
is also the whole polynomial centralizer of R. Indeed ad(R)/h is
a derivation with a slice D_0, kills x,q,B, and has associated
momentum-graded derivation -partial_(D_0); successively subtracting
the leading invariant symbols proves exhaustion. Only the explicit
Ore subalgebra (3) is needed below.

## A complete global obstruction space

The relative one-form cohomology is
\[
 H=\Omega^1_{K[x,y]/K}/(dK[x,y]+K[x,y]\,dr).         \tag{4}
\]
Every relative one-form is closed, since (2) makes the relative
cotangent bundle a line bundle. Evaluation on theta identifies
relative one-forms with K[x,y]. In particular
\[
 \alpha=4y^2dx+(2xy-1)dy,\qquad \alpha(\theta)=1.
\]
It follows that
\[
                 H \simeq K[x,y]/\theta K[x,y],
                 \quad F\longmapsto[F\alpha].        \tag{5}
\]
This is a quotient of vector spaces and K[r]-modules, not a
quotient ring or a finitely generated K[x,y]-module.

There is a finite exact normal-form algorithm for every polynomial
form, without any degree cutoff. First remove dy terms using
\[
 [x^i y^jdy]=-\frac{i}{j+1}[x^{i-1}y^{j+1}dx],
 \qquad [x^a dx]=0.
\]
Write e_(a,b)=[x^a y^b dx], b>=1, and assign weight d=a-b.
The remaining relations, obtained from x^i y^j dr, are precisely
\[
 e_{a,b}+\frac{2b-a}{b+1}e_{a+1,b+1}=0,\quad b\ge1,
 \qquad e_{a,1}=0\quad(a\ge2).                       \tag{6}
\]
They never mix different d. At each integer weight the quotient
is one-dimensional:

| Weight d | Nonzero basis representative |
|---|---|
| d<=-1 | [y^(-d) dx] |
| d=0 | [xy dx] |
| d>=1 | [x^(2d+1)y^(d+1) dx] |

For d<=0 every coefficient in the recurrence (6) is nonzero.
For d>=1 it kills the initial chain through b=d, has a zero
coefficient at b=d, and leaves one free chain starting at b=d+1.
This proves both spanning and independence, rather than assuming
that a finite monomial calculation exhausts (4).

Consequently a polynomial F is theta of another polynomial if and
only if every coefficient in the normal form of F alpha is zero.
This is an all-degree linear membership test for the specified
derivation, not a search through bounded candidate primitives.

## The connection makes the missing tower explicit

A polynomial lift of partial_r is
\[
 V=(1-2xy)\partial_x+4y^2\partial_y,\qquad V(r)=1.
\]
Lie differentiation induces a well-defined connection nabla on H,
with [nabla,r]=1. A different lift differs by a vertical vector
field, whose action on relative cohomology is zero by Cartan's
formula.

Let eta=[y dx]. Direct differentiation and (6) give
\[
 \nabla\eta=3[y^2dx],\qquad r\nabla\eta=0,
 \qquad
 \nabla[y^m dx]=\frac{2m(2m+1)}{m+1}[y^{m+1}dx].     \tag{7}
\]
For k>=1 another application of (6) gives
\[
                  r^k\eta=-[x^{2k-1}y^k dx].        \tag{8}
\]
The vectors r^k eta (k>=0) and nabla^j eta (j>=1) have
distinct weights and are exactly the bases just listed. Hence,
for the first Weyl algebra D_r=K<r,partial_r> with
[partial_r,r]=1,
\[
                      H\simeq D_r/D_r(r\partial_r).  \tag{9}
\]
The map sends the cyclic vector to eta. The relation and the
displayed independent basis prove that there are no extra
relations.

Its r-power torsion is
\[
 T=K[\nabla]\nabla\eta\simeq D_r/D_r r,\qquad
                         H/T\simeq K[r].             \tag{10}
\]
The latter is the usual module with its derivative connection.
Thus T is not finitely generated over K[r], although it is cyclic
as a D_r-module. The extension (10) does not split as D_r-modules:
nabla is injective on H by the basis above, whereas K[r] has a
nonzero horizontal vector. Localizing at r kills exactly T and
leaves K[r,r^{-1}] eta with nabla eta=0.

## Explicit invisible classes of every order

For each integer n>=1 put
\[
 \omega_n=y^n du,\qquad
 E_n(u)=\int_0^u s^n(s-1)^n\,ds.
\]
These are polynomial forms and polynomial functions, and
\[
 r^n\omega_n=dE_n(u),\qquad
 E_n(1)-E_n(0)=(-1)^n\frac{(n!)^2}{(2n+1)!}\ne0.     \tag{11}
\]
The class of omega_n is nonzero and has exact r-torsion order n.

Here is an independent all-n proof of exact order. If
r^(n-1) omega_n=dG+Q dr, equation (11) would imply
d(E_n-rG) wedge dr=0. The polynomial constants along theta are
exactly K[r]: after inverting r, the ring is
K[r,r^{-1},x,x^{-1}], theta=x^2 partial_x, whose constants
are K[r,r^{-1}]; intersection with K[x,y] is K[r], since a
negative power of r has a pole on x=0. Therefore E_n-rG=P(r).
Restrict to the two components x=0 and u=0. The same value P(0)
would have to equal both E_n(1) and E_n(0), contradicting (11).

On every nonzero fiber omega_n=d(E_n/r^n). Its restriction to
each component of the zero fiber is zero as well, since u is
constant there. Thus every fiber period vanishes, including
those of the exceptional fiber, but the global polynomial
relative class is nonzero.

In contrast eta restricts on r=c!=0 to
\[
                  (c/x^2-1/x)\,dx
\]
and has residue -1. Generic periods detect the torsion-free
part of (10) and discard the entire tower T. Smoothness,
generic exactness, or checking each individual fiber cannot
supply a polynomial primitive in the total family.

## Exact classical and quantum correction obstructions

For any F in K[x,y], the vertical rule
\[
 \Delta_F(x)=\Delta_F(y)=0,\qquad \Delta_F(b)=F
\]
is a Poisson derivation of (1) and an ordinary algebra derivation
of (3), with Delta_F(h)=0. In the quantum algebra over K((h)),
\[
 \Delta_F=h^{-1}\operatorname{ad}(H)
 \quad\Longleftrightarrow\quad F=-\theta(H)
                \text{ for some }H\in K((h))[x,y].   \tag{12}
\]
To prove necessity, commutation of H with x rules out positive b
order: the highest commutator coefficient for b order m>0 is
h*m*x^2 times the nonzero leading coefficient. Thus H depends
only on x,y, and the remaining equation is (12).
The identical argument with the Poisson bracket proves the
classical assertion. Membership in theta K[x,y] is preserved
under characteristic-zero field extension, as is also evident
from the normal form (6).

In particular take
\[
                         F_n=r y^n.
\]
Then Delta_(F_n) is not globally Hamiltonian, classically or
quantum mechanically. On r nonzero it is exactly Hamiltonian:
\[
 H_n=E_n(u)/r^n,\qquad
 \{H_n,-\}=\Delta_{F_n},\qquad
 h^{-1}[H_n,-]=\Delta_{F_n}.                          \tag{13}
\]
Indeed theta(u)=-r. The least nonnegative k such that
r^k Delta_(F_n) is globally Hamiltonian is exactly n.
This is a genuine unbounded denominator requirement in the
actual R centralizer, not a bound suggested by finite samples.

The assertions concern vertical derivations fixing the base
x,y. They do not classify every derivation of C_h, nor prove
that a particular fourth-image defect of Long lies in this
vertical subspace. No fixed quantum S,T is completed or
excluded by (12) alone; the later argument (16) also uses
full-frame infinitesimal rigidity.

## Which fiber translations are actual canonical changes

The same normal form gives an exact extension test for finite
automorphisms, not just derivations. Put w=-2z/3, so
\[
 [p,x]=[w,y]=h,\qquad b=x^2p-(1+2xy)w.
\]
The fiber translation
\[
                   \sigma_F(x,y,b)=(x,y,b+F(x,y))     \tag{14}
\]
is always an automorphism of the polynomial Ore algebra (3).
It extends to a polynomial canonical transformation of the
full rank-two Weyl algebra fixing x,y if and only if
\[
                              [F\alpha]=0\text{ in }H. \tag{15}
\]
The same statement holds classically.

For necessity, any full Weyl assignment fixing x,y must have
p'=p+A(x,y) and w'=w+C(x,y), by the joint centralizer of
the positions. The last Weyl relation is C_x=A_y.
Polynomial integration gives a polynomial G with
A=G_x and C=G_y. Hence b'=b+theta(G), so (5) proves
(15). Conversely these gradient shifts give the required
polynomial canonical automorphism whenever F=theta(G).
The argument works over K((h)) and for any spatial degree;
the classical proof uses the identical partial-derivative
equations. No classification of arbitrary canonical changes
that move the positions is asserted.

For F_n=r y^n, (14) therefore has no such polynomial
canonical extension, although it has one after inverting r,
using G=-E_n/r^n. Replacing F_n by r^k F_n permits an
extension exactly when k>=n. Thus the exceptional torsion
measures an actual canonical-transport boundary in Long's
centralizer. A fiber automorphism, even one canonically
extendible on every nonzero fiber, cannot automatically be
used as a global canonical gauge.

## Nonadmissible quantum triples arbitrarily close to admitted ones

There is a direct consequence for Long's actual admission problem.
Let (R,S_h,T_h,D_h) be any admitted full formal Weyl tuple with
Long's classical reductions and with R undeformed. Such tuples
exist by [QFL1](PROOF.md) and
[QCR1's formal R-fixing](CENTRALIZER_NORMAL_FORM.md).
For every N>=1 apply the exact Ore automorphism
\[
 \tau_N(x)=x,\quad\tau_N(y)=y,\quad\tau_N(b)=b+h^N
\]
to the h-adic Ore completion containing S_h,T_h, and put
\[
                 S'_h=\tau_N(S_h),\qquad T'_h=\tau_N(T_h).
                                                               \tag{16}
\]
These images satisfy all three partial Weyl relations exactly,
have the same classical reductions, and agree with the original
triple modulo h^N. Nevertheless **no polynomial-coefficient
formal fourth image completes this fixed new triple**.

First note that Delta_1=partial_b is not inner: its class is
\[
                         [\alpha]=3[y^2dx]\ne0.       \tag{17}
\]
It is r-torsion of order one, with theta(-u)=r and
theta(-1/x)=1 after localization. Thus tau_N is a particularly
simple instance of (14) whose leading derivation is outer.

Suppose a completing D'_h existed. Its classical image differs
from D only by P(R), by QCR1; subtract that polynomial in R.
For every subsequent order k<N, if the first difference of
D'_h and D_h is h^k C, the Weyl relations and agreement of the
partial tuples imply
\[
                      \{C,R\}=\{C,S\}=\{C,T\}=0.
\]
Again QCR1 gives C in K[R]. Subtracting h^k C(R) preserves all
relations with the new triple. Repeat finitely often. We may
therefore assume the two full tuples agree modulo h^N.

Their difference at order N is now a polynomial solution of
the linearized full canonical-frame equations. The H^1=0
argument in QFL1 makes that solution Hamiltonian: there is
a polynomial H whose brackets with the four classical
images are these four differences. This uses the polynomial
inverse Jacobian of the full symplectic differential frame,
not an inverse polynomial map.

In particular {H,R}=0, and on S,T the derivation {H,-}
equals partial_b. Both derivations act on C=K[x,y,b].
They agree on r,S,T, which are algebraically independent:
contracting a relation among dr,dS,dT with the Hamiltonians
of S,T first kills its dT,dS coefficients, then its dr
coefficient. Since K(x,y,b) is a separable algebraic extension
of K(r,S,T), a derivation agreeing on those three functions
agrees everywhere. Hence {H,-}|_C=partial_b.

But H belongs to C by {H,R}=0. Commutation with x then
forces H_b=0, and the bracket with b would require
-theta(H)=1. Equation (17) rules this out. This contradiction
proves the assertion for every N and every admitted starting
tuple. The same argument works for any non-inner leading
vertical derivation in this centralizer.

More generally, for any F in K[x,y], replace b+h^N in
(16) by b+h^N F. The resulting fixed quantum triple admits
a formal polynomial-coefficient fourth image **if and only
if [F alpha]=0 in H**. Necessity is the same first-difference
argument, with Delta_F in place of partial_b. For sufficiency,
write F=theta(G) and use the exact full gradient automorphism
\[
 p\mapsto p+h^N G_x,\qquad
 w\mapsto w+h^N G_y,\qquad x,y\mapsto x,y.
\]
It restricts to the prescribed fiber translation and transports
D_h as well. This transformation and its inverse have finite
spatial degree, so they preserve the existence of a common
spatial bound for the transported full tuple. For nonzero
[F alpha], no completion exists even without such a bound.

This is a complete admission test for this specified family
of h-small fiber translations of an admitted Long triple.
It does not characterize all admissible triples or arbitrary
fiber automorphisms.

This does not contradict full formal unobstructedness: the
new partial triple is fixed, and its fourth-image equation
fails admission. Its finite jets below N still extend, using
the original full tuple. Nor does it obstruct Long's formal
class, which has the admitted starting tuples. It shows
precisely why exact pair lifting and arbitrarily long
matching full jets do not ensure admission of a prescribed
all-orders pair.

## Research consequence and the next admission gate

Keep all four images and their relations during a constructive
lift. Equations (16)--(17) show that a pair-only correction can
leave the admissible locus at any chosen order while preserving
every partial relation exactly.

A proposed completion argument through generic differential
operators must retain both parts of (10). After showing its
defect is of the form Delta_F, apply the global normal form
of F alpha. A nonzero generic period is an obstruction; zero
generic period leaves the exceptional torsion test. Both must
vanish to produce a polynomial Hamiltonian in this class.

The next useful test is to identify this defect for an actual
compatible quantum triple, or prove that a chosen construction
keeps it outside T at every order. Merely computing periods
on the generic algebraic curve is insufficient. Even if the
class vanishes coefficient by coefficient, a common spatial
bound for the resulting primitives remains separate.

This complements [QCR1's normalized conjugate](CENTRALIZER_NORMAL_FORM.md):
uniqueness after admission does not establish global admission.
It also separates the present global exactness issue from
[QSO1's finite-order string obstruction](ALGEBRAIC_CURVE_STRING_OBSTRUCTION.md)
and [QFL1's formal lifting](PROOF.md). The latter provides full
formal lifts through choices of all images and does not say
that every independently chosen fiber pair is admissible;
(16) now gives counterexamples with the actual Long reduction.

## Literature and verification boundary

The relative-de-Rham framework and its dependence on behavior
at infinity are classical; see Philippe Bonnet,
[Relative cohomology of polynomial mappings](https://arxiv.org/html/math/0602273v1),
especially the introduction and Theorem1.3. That theorem's
finiteness hypotheses are not automatic here: for any positive
weights the leading form of r is x^2 y, whose ideal is not
radical. We do not import its finite-free conclusion.
Bonnet--Dimca's earlier
[Relative differential forms and complex polynomials](https://doi.org/10.1016/S0007-4497(00)01055-1)
is a bibliography pointer, not an independently audited proof
input. The quotient, connection, exact torsion orders and
quantum vertical criterion above have self-contained proofs.
No novelty priority is asserted.

[Exact controls](check_relative_period.py) verify the coordinate
translation, quotient relations, connection identities, explicit
primitives, endpoint differences, normal forms and torsion orders
at the stated sample indices. The universal exhaustion and
nonmembership statements use the written proofs, not those
finite controls. No independent replay, formal verification,
external review, bounded Long lift or DC2 solution is claimed.
