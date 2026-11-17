# Every bounded fixed-R change is a translation

14 November 2026. **QGR1.** In rank two, every finite polynomial
endomorphism of Long's quantum R centralizer that fixes R is a
fiber translation. Every finite Weyl endomorphism fixing R is a
gradient automorphism fixing both positions. The same statements
hold classically.

Consequently QRP1's admission obstruction covers all bounded
changes of this centralizer, rather than one selected shear family.
When their correction has a global Hamiltonian primitive, the
normalized primitive has degree exactly one less than the correction.
Thus a bounded admissible centralizer change always has a bounded
full canonical extension. No bounded lift of Long is constructed.

## Algebra and meaning of bounded

Use the coordinates of [QRP1](RELATIVE_PERIOD_TORSION.md):
\[
 y=-3q/2,\quad w=-2z/3,\quad
 r=R/2=x(1+xy),\quad u=1+xy,\quad
 b=B/3=x^2p-(1+2xy)w.
\]
Let K be any characteristic-zero field. With [p,x]=[w,y]=h,
the centralizer is the Ore algebra
\[
 C_h=K[h,x,y][b;h\theta],\qquad
 \theta=x^2\partial_x-(1+2xy)\partial_y,\qquad
 [b,g]=h\theta(g).                                  \tag{1}
\]
We use either a field with a fixed nonzero h, or the finite PBW
algebra over K[[h]]. The latter permits scalar power series,
but only finitely many spatial monomials. An endomorphism is
linear over that scalar ring and fixes h.

The centralizer description works over K[[h]] without division
by h. Indeed put
\[
 d=D_0=(1-2xy)p/2+2y^2w.
\]
The polynomial inverse formulas are
\[
 p=2(1+2xy)d+4y^2b,\qquad
 w=2x^2d-(1-2xy)b.
\]
Normal ordering gives
\[
 A_{2,K[[h]]}(h)=C_{K[[h]]}[d;h\delta],\quad
 \delta x=(1-2xy)/2,\quad\delta y=2y^2,\quad
 \delta b=-3yb,\quad\delta r=1/2.                   \tag{2}
\]
The linear momentum change has determinant -1/2; its associated
graded is a polynomial coordinate change over K[[h]][x,y].
It supplies spanning and independence for the finite PBW
expansion in b,d. If F=sum c_j d^j, c_j in C, then
\[
                  [F,r]/h=\tfrac12\sum j c_jd^{j-1}.
\]
Hence its centralizer is exactly C. Finite conversion between
these PBW presentations uses polynomial coefficients and finitely
many reorderings. It preserves the property of having a common
finite spatial bound; it need not preserve its numerical value.

## All center-preserving endomorphisms of C are translations

Let tau be any endomorphism of C with tau(r)=r. No inverse is
assumed. First work over a field with h nonzero. Put
\[
                 X=\tau(x),\quad Y=\tau(y),\quad
                 U=\tau(u)=1+XY,\quad Z=\tau(b).
\]
Since XU=r and the b-order of a product is the sum of its
orders, both X and U have b-order zero. The same is true of
Y, by XY=U-1. If U-1 were zero then X=r, which is
excluded by the commutator relation below. Thus X,Y,U lie
in the commutative polynomial ring K[x,y].

Neither X nor U is a unit. If X were scalar, or if U were
scalar so that X were a scalar multiple of r, centrality would
contradict
\[
                              [Z,X]=hX^2\ne0.       \tag{3}
\]
The factors x and u of r are distinct irreducible polynomials.
Unique factorization therefore leaves only
\[
                 X=cx,\ U=c^{-1}u,
       \quad\hbox{or}\quad X=cu,\ U=c^{-1}x
                                                        \tag{4}
\]
for a nonzero scalar c. The second case would make u divide
c^{-1}x-1 in K[x,y], impossible after setting y=-1/x.
In the first case the equation U=1+c xY, evaluated at x=0,
gives c=1. Then Y=y.

It follows from (3) that Z-b commutes with x. An element of
positive b-order m has leading commutator with x equal to
h*m*x^2 times its leading coefficient, so it cannot commute.
Consequently
\[
             \boxed{\tau(x)=x,\quad\tau(y)=y,\quad
                         \tau(b)=b+F(x,y).}         \tag{5}
\]
Every F defines such an automorphism; its inverse replaces F
by -F. This proves exhaustion in every b-order and every
coefficient degree. It also proves that every such endomorphism
is automatically invertible.

For scalar coefficients in K[[h]], pass to K((h)) and apply
the same proof. The resulting identities descend, and F
retains its original regular scalar coefficients and finite
spatial support. The proof does not require tau to be the
identity modulo h. Replacing the commutator by the Poisson
bracket gives the identical classification in characteristic
zero for the classical algebra K[x,y,b].

## Every full finite endomorphism fixing r is an automorphism

Let Psi be a finite Weyl endomorphism of A_2 with Psi(r)=r.
It preserves the centralizer, so (5) forces
\[
                            \Psi(x)=x,\qquad\Psi(y)=y.
\]
The remaining position commutators force
\[
                  \Psi(p)=p+A(x,y),\qquad
                  \Psi(w)=w+C(x,y).
\]
The last relation is C_x=A_y. Polynomial integration in
characteristic zero gives a polynomial G with
\[
                   A=G_x,\quad C=G_y,\qquad
                   \Psi(b)=b+\theta(G).             \tag{6}
\]
This is an exact automorphism, with inverse given by -G.
Its scalar coefficients can also lie in K[[h]]; integration
introduces no spatial denominators or parameter poles.
The classical symplectic assertion follows in the same way.

The part of (6) classifying automorphisms once they fix the
positions is classical: see Milicic,
[Localization and Representation Theory of Reductive Lie Groups,
Chapter1, Lemmas1.1--1.2](https://www.math.utah.edu/~milicic/Eprints/book.pdf).
Those lemmas describe differential-operator automorphisms
fixing functions by closed one-forms. Here the factorization
argument proves that merely fixing this particular r already
forces the positions to be fixed. The elementary proof above
does not import an unverified general automorphism theorem.

Thus no finite quantum or classical canonical change preserving
R can secretly mix its positions with momenta. Searching for
higher-degree or higher-order exceptions within that stabilizer
cannot enlarge the group.

## A sharp bound for a global Hamiltonian primitive

The only obstruction to extending (5) is F in theta K[x,y],
which [QRP1](RELATIVE_PERIOD_TORSION.md) decides by the full
relative cohomology normal form. This membership test includes
the exceptional torsion discarded by generic periods.

For G in K[x,y], remove its r-polynomial ambiguity by
\[
 {\cal N}_r(G)=G-\sum_{k\ge0}
               [x^{2k}y^k]G\,r^k.                  \tag{7}
\]
Each r^k has exactly one pivot monomial x^(2l)y^l,
namely l=k, of coefficient one. Thus (7) is an idempotent,
has kernel K[r], never increases degree and leaves theta(G)
unchanged.

Write
\[
 \theta=\theta_+ -\partial_y,\qquad
 \theta_+=x^2\partial_x-2xy\partial_y.
\]
On a monomial x^a y^c,
\[
                  \theta_+(x^a y^c)
                         =(a-2c)x^{a+1}y^c.
\]
Different monomials have different images. Its kernel consists
exactly of the pivot monomials x^(2k)y^k. A nonzero normalized
G has no such monomial, including its constant term. Therefore
its highest homogeneous part is not killed by theta_+.
The lower operation -partial_y has degree two less, so
\[
 \boxed{\quad \theta(G)=F\ne0,\ {\cal N}_r(G)=G
                 \quad\Longrightarrow\quad
                         \deg G=\deg F-1.\quad}     \tag{8}
\]
There is a unique normalized primitive: two primitives differ
by K[r], as proved in QRP1. Equations (7)--(8) hold in all
degrees, not just for the displayed controls.

In particular, if F_h has a common spatial bound D and each
coefficient has a global polynomial primitive, its normalized
primitive G_h has a common bound D-1. It belongs to the finite
PBW coefficient ring K[[h]][x,y]. The resulting full gradient
automorphism has degree at most max(1,D-2), as does its
inverse. For F_h=0 choose G_h=0. This is a positive bounded
extension theorem within the complete fixed-R group.

The same argument works for Laurent scalar series: the
normalized primitive has the same parameter valuation as
F_h. In particular no hidden parameter pole is required
when F_h is regular at zero.

## Degree-controlled division computes both obstruction and primitive

The bound above gives a constructive decomposition for every
polynomial F, whether or not it is globally exact:
\[
 F=\theta(G)+Q,\qquad {\cal N}_r(G)=G,\qquad
 Q\in\operatorname{span}_K\{y^j,\ x^{2j+1}y^j:j\ge0\}. \tag{11}
\]
The decomposition is unique, with deg G<=deg F-1 and
deg Q<=deg F whenever the corresponding output is nonzero.

For a current monomial c*x^i*y^j with i>=1 and i-1-2j
nonzero, use
\[
 x^iy^j=
 \theta\!\left({x^{i-1}y^j\over i-1-2j}\right)
       +{j\over i-1-2j}x^{i-1}y^{j-1}.             \tag{12}
\]
Record the displayed primitive and continue with the residual
term, whose degree is two lower. A monomial with i=0 or
i=2j+1 is retained in Q. Processing highest degree first
terminates for every finite input and produces (11).
No recorded primitive has a pivot monomial.

A nonzero Q cannot belong to theta K[x,y]: its highest
homogeneous part consists entirely of monomials excluded
from the image of theta_+, whereas (8) forces the top of
every nonzero exact input to lie in that image. This proves
independence, uniqueness and the global obstruction test
Q=0. Thus (11) is a direct-sum decomposition, not a bounded
search for a primitive. It works with arbitrary scalar
series coefficients of finite spatial support as well.

This supplies both an exact rejection remainder and an actual
polynomial primitive when the test passes. It agrees with
QRP1's cohomology test while avoiding a separate linear
system at each degree.

## Complete admission test for every bounded fixed-R change

Fix an admitted full formal Long tuple (R,S_h,T_h,D_h) with
R exact. Such tuples exist by [QFL1](PROOF.md) and
[QCR1](CENTRALIZER_NORMAL_FORM.md); they are not known to have
bounded spatial degree.

Let tau be any finite-PBW endomorphism of C over K[[h]]
fixing R, with tau equal to the identity modulo h.
By (5), tau is b to b+F_h(x,y), where F_h is in hK[[h]][x,y]
and has a common finite degree bound. The following are
equivalent:

1. The fixed triple R,tau(S_h),tau(T_h) admits a
   polynomial-coefficient formal fourth image.
2. Every coefficient of F_h has zero global class in
   K[x,y]/theta K[x,y].
3. Tau extends to a finite-PBW full Weyl automorphism fixing R.

The implication 2 to3 follows from (8) and the exact gradient
extension (6), and 3 to1 transports D_h.

For 1 to2, suppose the first nonzero obstruction class occurs
at order N. The finitely many lower coefficients are global
primitives, so their exact gradient automorphisms remove them.
The first remaining change is h^N Delta_(F_N), with a nonzero
global class. QRP1's first-difference argument now forbids
completion: subtracting lower K[R] conjugate differences
would turn a hypothetical completion into a full-frame
Hamiltonian realizing that non-inner derivation. Higher
coefficients cannot change the contradiction at order N.

Every admitted bounded change of this centralizer therefore
has a bounded full extension and inverse. It preserves the
existence of a common spatial bound for the original full
tuple. The group of all these changes, rather than only a
preselected translation ansatz, is now covered.

Over a field the whole restriction and obstruction sequence is
\[
 0\longrightarrow K[r]/K
 \longrightarrow \operatorname{Stab}_{A_2}(r)\simeq K[x,y]/K
 \xrightarrow{\ G\mapsto\theta(G)\ }
 \operatorname{Aut}_{K[r]}(C_h)\simeq K[x,y]
 \longrightarrow H
 \longrightarrow0.                                \tag{13}
\]
These are additive groups; the stabilizer consists of the
gradient automorphisms (6). The kernel is the r-potential
freedom, and the cokernel is exactly QRP1's global obstruction
space. The same description holds classically.

## Why the finite spatial bound is essential

The formal Hamiltonian gauge
\[
                         \exp(h^2\operatorname{ad}(b)/h)
\]
fixes r, but on positions it acts by
\[
 x\mapsto {x\over1-h^2x},\qquad
 y\mapsto(1-h^2x)(y-h^2u).
                                                        \tag{14}
\]
Indeed theta(x)=x^2, theta(y)=-(1+2xy), and theta(u)=-r.
The new u is u(1-h^2x), so r is preserved exactly.
Every coefficient is polynomial, but the coefficient of
h^(2j) in the first image is x^(j+1). Its spatial degree
is unbounded.

This is an actual full formal canonical gauge, since it is
an exponential of an h-adically small normalized inner
derivation. Thus (5)--(6) cannot be extended to arbitrary
formal gauges. A formal factorization can acquire nonconstant
units such as 1-h^2x; the finite polynomial factorization
argument does not apply to them.

## A bounded first observable need not admit bounded target normalization

A second control distinguishes formal target normalization from a
bounded one. For any nonzero c(h) in hK[[h]], put
\[
                            r_h=r+c(h).
\]
It has finite spatial degree three and the bounded mate d:
[d,r_h]=h/2. Over K((h)), however, r_h is irreducible.
As a polynomial in y it is x^2 y+x+c; its two coefficients
are coprime in K((h))[x], so Gauss's lemma proves
irreducibility. A factorization in the full Weyl algebra
would have momentum order zero in both factors, hence would
already be a factorization in that commutative polynomial ring.
By contrast r=xu is a product of two nonunits.

It follows that no finite full Weyl automorphism can send
r_h to r, irrespective of its degree. The formal gauge
\[
                   \exp(-2c(h)\operatorname{ad}(d)/h)
\]
does send r_h to r, since [d,r]/h=1/2. Its coefficients are
polynomial at every order, but its spatial degree must be
unbounded by the irreducibility argument.

One may also shift the first image of any admitted full
formal Long tuple by the scalar 2c(h), leaving the other
images unchanged. It remains admitted with the same classical
reduction, and the same target-normalization obstruction
applies. A bounded full tuple is not supplied by this construction.

This control does not obstruct normalization using both source
and target changes: subtracting that scalar from an assigned
first image is itself a bounded source translation. It rules
out inferring a bounded target Hamiltonian normalization
merely from a bounded first observable, a bounded mate or
formal normalizability. The choice of allowed side matters.

## Consequence for the remaining Dixmier problem

The bounded fixed-R canonical group and its admission gate
are now exhausted. Within this group there is no additional
momentum-mixing freedom and no primitive-degree escape after
global exactness has been verified.

This does not prove that all bounded full lifts are connected
by bounded centralizer automorphisms. Formal gauge equivalence
does not establish that stronger assertion. Nor does it give
a bounded change sending a deformed first image back to R,
or cover canonical changes that move R.

The hypothesis Psi(R)=R is a stabilizer condition. Long's
candidate Weyl map instead requires its first generator
image to be R; it does not fix R as a polynomial observable.
The fixed-observable endomorphism theorem therefore does
not exclude Long's bounded quantization.

A useful next step must control a full tuple outside this
exhausted bounded gauge freedom, establish bounded reachability
between full lifts with the allowed source and target changes, or construct a bounded full scalar-anomaly
carrier. No larger fixed-R shear or operator-order census is
needed. [QFL2](BOUNDED_REALIZATION.md) remains the endpoint:
one admitted full lift with a common bound gives a genuine
Weyl counterexample for a noninvertible classical input.

[Exact controls](check_bounded_gauge.py) verify the Ore
presentation, normalized primitive identities and degree
equalities, gradient shifts and the formal moving-position
control. The universal endomorphism exhaustion, sharp bounds
and full admission equivalence use the written proofs.
No independent replay, formal verification, external review,
novelty priority, bounded Long lift or DC2 solution is claimed.
