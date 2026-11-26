# Effective rational-bias neighborhoods for two opposed factorial heights

Fix an ordinary two-monomial core, fixed positive integer gaps
\(d,e\), and a fixed signed tail gap \(h=dr_+-er_-\).
Append \(r_+\) copies of \((L,L+d)\) and \(r_-\) copies of
\((Y+e,Y)\), with ordinary heights and counts. Write
\(r_\pm=c_\pm r+s_\pm\) with fixed integer \(c_\pm,s_\pm\)
and \(u=1/r,\ p=r/L,\ q=r/Y\), as in the
[local Puiseux theorem](MULTIVARIATE_OPPOSED_DIFFUSE_PUISEUX_HEIGHT_BRANCHES.md).

**Theorem.** Given a positive Gaussian contact at adjacent
orders with specified rational bias, one can compute a
neighborhood of the contact and a height cutoff beyond
which no ordinary collision occurs in that neighborhood.
The contact is specified by its two positive algebraic
clocks. This holds for both interior and boundary contacts,
with no simplicity or transversality hypothesis.

The local argument uses the finite analytic algebra of the
two Weierstrass polynomials for the collision branches.
It constructs a single integer auxiliary polynomial in both
heights and the count. A finite-jet ideal calculation then
forces every sufficiently high ordinary collision onto a
computable list of algebraic height curves. This calculation
only searches for a finite Taylor approximation to a
meromorphic relation.

## A finite analytic algebra at the contact

Interchange the two columns and height names if necessary
so that \(p_*>0\) at the contact. The nonconstant count is
positive at all sufficiently high samples. An exact
collision is necessary for all degree-at-least-two Sylvester
minors to vanish, by strict positivity of the even moment
on the real axis and conjugation. Call this list of analytic
minors \(M_1,\ldots,M_\ell\) after the exact physically common
normalization. Also let \(A,D\) be the angular quotient
and physical resultants. These are the analytic equations
in the [central isolation argument](MULTIVARIATE_OPPOSED_DIFFUSE_PUISEUX_HEIGHT_BRANCHES.md#the-exact-complex-one-parameter-germ).

Use centered logarithmic clock coordinates
\(W=\log(w/w_*), V=\log(v/v_*)\) and put \(\theta=T_*\).
Both \(p,q\) are rational linear functions of \(W,V\) plus
rational multiples of \(\theta\), since the bias is specified
rational. Hermite–Lindemann gives the transcendence of
\(\theta\). After isolating all algebraic central
coefficients into a number field \(K\), Taylor calculation
can be done over \(K(\Theta)\) with a symbol \(\Theta\)
that is specialized to \(\theta\). As in the
[one-clock lift](MULTIVARIATE_OPPOSED_AUXILIARY_POLYNOMIAL_RATIONAL_BIAS.md#one-transcendental-coefficient-parameter),
the resultants and the minors lift to
\(K[[U_0,U_1,U_2,W,V]]\) with \(U_i=u\Theta^i\).
The lift applies to the minors by taking sums of products
of the exact coefficient expansions.

Prepare \(A\) in \(W\), obtaining a monic polynomial \(P\)
of degree \(a>0\). Take the norm of \(D\) in this quotient
and prepare it in \(V\), obtaining a monic polynomial \(Q\)
of degree \(b>0\). Reduce the coefficients of \(P\) in
\(V\) modulo \(Q\) and call the resulting polynomial
\(P^\sharp\). The algebra
\[
       \mathcal B=\mathbb C\{u\}[V,W]/(Q,P^\sharp)
\]
is free over the convergent series ring, with basis
\(V^i W^j\), \(0\le i<b,\ 0\le j<a\).
Every nearby ordinary collision is an evaluation on this
algebra. To represent an analytic function on it, first
take the Weierstrass remainder in \(W\) modulo \(P\),
and then take each coefficient's remainder in \(V\)
modulo \(Q\). All these operations preserve the stated
formal lift, and the resulting finite jets belong to
\(K(\Theta)\).

For clarity about effectivity, both central one-variable
preparations have finite orders, known to be finite by
Gaussian isolation. Their orders can be found by searching
successive exact Taylor coefficients. The first central
germ is nonzero because its angular eliminant is a nonzero
polynomial in the varying clock \(w\). For the second one,
the central norm is a power of the physical resultant at
\(W=0\), which has finite order in \(V\). Each finite jet of
the monic factors is computable by formal preparation.
Equivalently use Newton identities on the contour power
sums and take the residues of their finite parameter
expansions at the central roots. Weierstrass division
computes finite remainder jets by the same operations.
Uniqueness identifies the formal jets with the analytic
ones. None of these calculations tests an infinite series
for identity.

Analytic radii and supremum bounds in these calculations
are computable as well. For each preparation choose two
small concentric circles about its central root whose
closed outer disc contains no other central zero. A bound
from the first nonzero Taylor coefficient at the center
gives strict positive lower bounds on both circles.
Shrink the other variable discs by explicit holomorphic
coefficient estimates so that Rouché's theorem keeps
all prepared roots inside the inner circle. The norm
is bounded by the supremum of the original function on
the prepared roots, to the power of the monic degree.
For a remainder of \(f\) modulo a monic \(P_0\), the
formula
\[
 \frac{1}{2\pi i}\int_{|z|=R}
 \frac{f(z)}{P_0(z)}
       \frac{P_0(z)-P_0(w)}{z-w}\,dz
\]
bounds each remainder coefficient on smaller parameter
discs. Here \(|P_0(z)|\) is bounded below on the contour
by the corresponding power of the gap between the two
circles. Cauchy's estimate then turns any certified jet
vanishing into a computable uniform remainder bound.

## One effective auxiliary polynomial

Shrink the real neighborhood so that \(p\) is bounded
above and bounded away from zero. In the interior \(q\)
is also bounded away from zero. At the boundary, the
[exact mark-product height constraint](MULTIVARIATE_OPPOSED_MARK_PRODUCT_HEIGHT_CONSTRAINT.md#consequences-for-arbitrarily-separated-scales)
with the taller row marked gives \(L\ge cY^\beta\) at
any collision, with computable \(c,\beta>0\).
Since \(L=1/(up)\), there are a computable positive
constant \(c_0\) and integer \(k\ge0\) such that
\[
                              pq\ge c_0 u^k                     \tag{2}
\]
at all ordinary collisions under consideration. One may
take \(k=0\) at an interior contact.

For a polynomial \(R(L,Y,r)\) of degree at most \(D_0\),
use its cleared expression
\[
     R^\#=u^{D_0}p^{D_0}q^{D_0}R(1/(up),1/(uq),1/u).
\]
This is polynomial in \(u,p,q\). Take an unknown
\(H\in\mathbb Q[L,Y,r]\) of degree at most \(D_0\)
and impose vanishing of the \(\mathcal B\)-remainder
of \(H^\#\) through order \((k+2)D_0-1\) in \(u\).
The coefficient of \(u^j\) in any remainder coordinate
is polynomial in \(\Theta\) of degree at most \(2(D_0+j)\).
Indeed all dependence in the initial powers of \(p,q\)
has degree at most \(2D_0\), and the formal lift allows
at most \(2j\) more. The preparation, norm and division
are performed before specialization in the lifted
parameter ring. They are linear in the dividend
\(H^\#\), and introduce no other occurrences of the
symbol \(\Theta\). Thus the imposed constraints amount
to \(O(D_0^2)\) rational linear equations on
\(\binom{D_0+3}{3}\) unknowns, with fixed constants
depending on the preparation degrees and \(K,k\).
A finite exact search finds a nonzero \(H\).

At small positive \(u\), the certified remainder bound
and (2) give \(|H(L,Y,r)|\le C u^{D_0}\) at any ordinary
collision in the neighborhood, with computable \(C\).
Since \(H\) has fixed rational coefficients, its integer
evaluations lie in a fixed rational lattice. For all
sufficiently small \(u\), effectively, every such
collision therefore satisfies \(H=0\).

## An effective reduction to the finite catalogue

Take the finite rational height-curve list from the
[affine group catalogue](MULTIVARIATE_AFFINE_GROUP_HEIGHT_CURVE_CERTIFICATE.md)
for the two single row types. Formal divisor independence
holds by the
[opposed character argument](MULTIVARIATE_OPPOSED_AX_BRANCH_RANK.md#characters-of-the-factorial-ramps).
Discard geometrically irreducible curves on which the
restricted lattice of character relations has rank at
most one, and keep the geometrically reducible entries
for the moment. Let \(F(L,Y)\) be the product of the
remaining defining integer polynomials, with value one
for an empty list. Put \(f=pqF^\#\), using the clearing
degree for \(F\).

Work over the meromorphic series field \(\mathbb C((u))\).
The ideal generated by \(H^\#,M_1,\ldots,M_\ell\) in
the finite algebra \(\mathcal B\) after base change has
all its geometric points on \(f=0\). In fact these
points are convergent Puiseux germs by successive
preparation in \(V\) and \(W\). If \(pq=0\) identically,
the assertion is automatic. Otherwise their meromorphic
height and count functions satisfy \(H=0\), and the minors
show that the gcd has degree at least two. Their height
functions are nonconstant because \(p,q\) are bounded
near the contact and \(u\to0\). The Ax rank argument and
the analysis of a single character relation imply at least
two independent character relations. These give the
catalogue curve by the separated-function and Lüroth
arguments. For a geometrically irreducible catalogue
curve its function field injects into the branch field,
so its rank defect is the one tested in the filter.
This proves the assertion for every geometric point,
without making an integrality assumption on that point.

An element of a finite-dimensional algebra that vanishes
at every geometric point of a quotient is nilpotent in
the quotient. In this case the nilpotence index is at
most \(ab=\dim\mathcal B_{\mathbb C((u))}\): successive
nonzero powers of its multiplication operator must have
strictly decreasing image dimensions before reaching
zero. Thus \(f^{ab}\) is an ideal combination of the
\(H^\#,M_i\) over the meromorphic field. Solvability of the
finite linear system descends to \(K(\Theta)((u))\),
because all the multiplication data are over this
field and a field extension does not alter linear ranks.

An exact meromorphic identity is not needed to get a
computable height threshold. Choose
\[
                         E>ab\bigl(\deg F+(1+\deg F)k\bigr).
\]
Allow arbitrary finite Laurent polynomials over \(K(\Theta)\)
as coefficients of an ideal-combination approximation,
including all \(\mathcal B\) coordinates of those
coefficients. A finite search over pole and degree bounds,
using exact linear algebra on Taylor jets, finds a combination
congruent to \(f^{ab}\) modulo \(u^E\). Indeed any
meromorphic solution can be truncated so as to retain
each required jet. Polynomial dependence on \(\Theta\)
or rational-function denominators in it are checked before
specialization, using its transcendence. Multiplying
by a sufficient power of \(u\) to clear Laurent poles,
the analytic remainder estimate gives a computable
constant for the error \(O(u^E)\). At ordinary collisions
all ideal generators vanish. But if \(F(L,Y)\ne0\), then
it is an integer and
\[
                         |f|\ge u^{\deg F}(pq)^{1+\deg F}
                             \ge c_0^{1+\deg F}
                                 u^{\deg F+(1+\deg F)k}.
\]
The choice of \(E\) makes these bounds incompatible
beyond a computable cutoff. Thus every sufficiently high
ordinary collision near this contact lies on one of the
retained catalogue curves. This argument obtains a radius
from finite jets without deciding identities on a Puiseux
branch.

## Excluding the catalogue curves effectively

Consider the finitely many branches of each retained curve
where both heights diverge with positive real values.
Their Puiseux expansions and leading coefficients are
algebraic and computable by Newton–Puiseux. A convergence
radius and remainder bounds follow from a finite
truncation and implicit isolation of the resulting
simple root. Exclude branches whose limiting ratio
\(L/Y\) is incompatible with the contact, shrinking
the neighborhood and increasing the height cutoff if
needed. Geometrically reducible rational curves can be
handled by their finite geometric intersections, since
any rational point belongs to all conjugates of its
component. Singular points of an irreducible plane
curve also give finitely many directly computable pairs.

For an interior contact let \(\xi_*=q_*/p_*>0\) be its
rational reciprocal height ratio. On a remaining
algebraic branch let \(s=L^{-1/M}\). If \(L/Y-\xi_*\)
has a first nonzero term of order \(b/M<1\) in \(1/L\),
use the
[fractional clock deformation](MULTIVARIATE_AFFINE_GROUP_HEIGHT_CURVE_CERTIFICATE.md#fractional-changes-of-the-height-ratio-are-excluded).
Its prepared norm after putting
\(\zeta=p_*s^b\) is a nonzero analytic germ with
coefficients polynomial in the transcendental \(p_*\)
over a number field. Its first nonzero coefficient is
found by an exact finite search and bounded away from
zero by certified numerical evaluation of the positive
algebraic clocks and their logarithms. The same contour
and Cauchy estimates give a computable punctured
exclusion radius. If there is no such fractional term,
the expansion gives
\[
                         Y=\xi_*^{-1}L+\kappa+o(1)
\]
with algebraic \(\kappa\). The expansion bounds
effectively force the exact rational affine relation or
exclude lattice membership: \(Y-\xi_*^{-1}L\) lies in a
fixed discrete rational lattice. The
[effective opposed affine-height theorem](MULTIVARIATE_OPPOSED_INDEPENDENT_HEIGHT_RATIONAL_BIAS.md#first-marks-exclude-the-forced-affine-family)
handles the exact affine relation uniformly in the counts.

At a boundary contact the remaining compatible branches
have \(Y\asymp L^\beta\) with \(\beta>1\).
For \(1<\beta<2\) the same effective clock deformation
applies, now using the first term of \(L/Y\).
For \(\beta\ge2\), let \(C\) denote the normalization of
a geometrically irreducible catalogue curve. If the
geometric points at infinity form a single Galois orbit
over \(\mathbb Q\), the pole orders of each height
function are equal at all points of this orbit. Since
both heights have a pole at the selected branch,
\(\deg Y/\deg L=\beta\). The
[factorial divisor calculation](MULTIVARIATE_OPPOSED_BOUNDARY_POLE_DEGREES.md#divisors-of-the-factorial-relations)
depends only on that degree ratio and on at least two
character relations. It gives \(\beta<2\), so this
case cannot occur. This use of the divisor calculation
does not need an ineffective integral-point bound.

With at least two Galois orbits at infinity, choose a
positive integer \(a_0\) such that \(L+a_0Y\) has a pole
at every point at infinity. The finite choices that
cancel a leading pole coefficient can be computed and
avoided. At a positive integer height pair this
function is an integer. Runge's theorem on the
normalized projective curve over \(\mathbb Q\), with
just the archimedean place, effectively bounds all
such values, since there are at least two pole orbits.
Positivity bounds both integer heights by that sum.
[Levin, Theorem 1.2](https://jtnb.centre-mersenne.org/item/10.5802/jtnb.634.pdf)

Every stage uses finite inputs and a terminating search
or an explicit bound. This proves the effective
neighborhood theorem. The first few adjacent orders
may instead use the [ordinary low-terminal theorem](MULTIVARIATE_TERMINAL_FOURTEEN_ALL_PROFILES.md).

For the balanced case on a mirrored core, the
[effective radial and affine reduction](MULTIVARIATE_BALANCED_TWO_TYPE_HEIGHT_CUTOFF.md)
gives a global finite decision without assuming regularity of the
core. The present local theorem allows general cores and
nonzero signed twists at rational-bias contacts.

## Effective exclusion given a relation between irrational-bias clocks

The same calculation is effective at a specified positive
interior irrational-bias contact if an exact nonzero
polynomial relation over \(\mathbb Q\) between its two
logarithmic clocks is supplied. This is the effective form
of the [two-clock dependence theorem](MULTIVARIATE_OPPOSED_TWO_CLOCK_DEPENDENCE.md).
The relation is part of the input.

In this case use \(p_*,q_*\) as the two coefficient-field
parameters, which are invertible rational linear
combinations of the two logarithms. At least one logarithm
is transcendental, and adjoining a finite number field
\(K\) for the central algebraic coefficients preserves
transcendence degree one of the resulting field. Factor
the supplied relation over \(K\). Exactly one irreducible
factor vanishes at \((p_*,q_*)\): distinct factors are
relatively prime and their intersections have algebraic
coordinates, inconsistent with a transcendental coordinate.
Certified approximation of the logarithms detects
nonvanishing of all other factors in finite time. The
kernel of evaluation of \(K[P_*,Q_*]\) at the two
parameters is precisely the prime ideal generated by
that irreducible factor. Indeed the kernel has height one,
and the exhibited irreducible polynomial generates a
height-one prime. Thus arithmetic and exact equality in
the coefficient field can be carried out as rational
function arithmetic on that polynomial curve.

The formal lift now has six variables \(U_{ij}=uP_*^iQ_*^j\)
for \(i+j\le2\), as in the two-clock dependence theorem.
At the \(u^j\) coefficient of the remainder of \(H^\#\),
the total parameter degree is at most \(2(D_0+j)\).
Dividing by the irreducible relation in a total-degree
monomial order preserves this bound. The normal form uses
only monomials not divisible by the leading monomial of
that relation. Their number through total degree \(m\)
is \(O(m+1)\), because it is the difference of the numbers
of monomials through degrees \(m\) and \(m-\deg P\), where
\(P\) is the relation polynomial. Hence choosing \(H\)
again imposes only \(O(D_0^2)\) rational conditions.
The contour bounds and the integer lattice argument
force \(H=0\) at sufficiently high ordinary collisions.

Here \(p,q\) are bounded away from zero even on the
geometric branches used in the ideal argument. At
irrational bias the limiting height ratio is
transcendental by Gelfond–Schneider. Any nontrivial
character relation would separate into nonconstant
rational functions of \(L\) and \(Y\), whose leading
powers at infinity make that ratio algebraic. Thus the
characters have full rank on each geometric branch.
The identity \(H=0\) gives base transcendence degree
at most two, so the Ax bound excludes the degree-two
gcd minors on every branch. Consequently the finite
meromorphic ideal is the unit ideal. Approximate a
meromorphic representation of \(1\) by finite Laurent
jets as above. The generators are zero at an ordinary
collision, so the computable jet error excludes it
for all small positive \(u\). This proof makes no
assumption about which logarithmic relations actually
exist at irrational-bias contacts.
