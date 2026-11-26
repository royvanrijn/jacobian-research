# Finite height-curve certificates for two affine groups of factorial rows

Consider an ordinary binomial core and two independently varying
affine groups of repeated rows, in the setting of the
[two-group branch theorem](MULTIVARIATE_OPPOSED_DIFFUSE_PUISEUX_HEIGHT_BRANCHES.md#two-groups-of-affine-row-types).
The smaller exponents in each group have the form
\(a_i+\alpha_i L\) or \(a_i+\alpha_i Y\), with
\(a_i,\alpha_i\) rational and \(\alpha_i>0\). Fix the positive integer
gaps \(d_i\), orientations \(\epsilon_i\), and multiplicities
\(c_i r+s_i\), with \(c_i\ge0\) and \(c_i,s_i\) integers. Assume
\(\sum_i\epsilon_i c_i d_i=0\) and that the two Gaussian
curvature and bias directions are linearly independent.
The ordinary parameters \(L,Y,r\) are integers. Both heights
tend to infinity in the comparisons below.

At adjacent orders \(n,N=n+1\), normalize the two parents by
their constant coefficients and the common argument change in
the branch theorem. Denote their nonconstant coefficients by
\[
                   G_i(L,Y) F_i(L,Y)^r,\qquad 1\le i\le n+N,
                                                               \tag{1}
\]
where \(G_i\) is a nonzero rational function including the core
and the offset multiplicities, and the characters \(F_i\) are
products of ramp factors separated in \(L,Y\). Normalize the
factors to have value one at infinity. Assume the \(F_i\)
are multiplicatively independent modulo constants as rational
functions of the two independent formal height variables.
This is an exact rational divisor-rank condition on the row
groups and the two orders.

**Theorem.** A finite list of rational plane curves, computable
from the row groups and the adjacent orders, contains every
accumulating branch of ordinary common roots approaching a
positive rational-bias Gaussian contact. This includes boundary
contacts and allows arbitrary multiplicities or tangencies of the
reference roots. At an interior rational-bias contact the
accumulating branch must in fact lie on a rational affine line
\(Y=\rho L+\kappa\). On the boundary, any branch with
\(Y\sim c L^\beta\) and \(1<\beta<2\) is excluded.
The curve list is independent of the core and of the ordinary
heights and count. For each listed curve there is a
finite algebraic test that excludes accumulation if the
common-factor locus has codimension at least two in the
coefficient subtorus prescribed by that curve. In particular a
successful test for each curve is a terminal-uniform exclusion
certificate when applied through the absolute cutoff.
It does not require a numerical bound on the height at which
the exclusion begins.

The list is finite despite the unbounded exponents in possible
multiplicative relations. Two independent relations generate
a rational intermediate field with at least three marked fibers.
All preimages in those fibers lie in fixed factorial grids.
This bounds the degrees of the two maps to the intermediate
field and makes the list explicit.

## Two relations at every accumulating branch

At a positive rational-bias contact put \(\theta=T_*>0\).
In the analytic coordinates \(u=1/r,\ p=r/L,\ q=r/Y\), the
central \(p_*,q_*\) are rational multiples of \(\theta\):
\(T,\lambda\) are independent rational linear combinations of
\(p,q\), and the bias is \(\lambda_*/T_*\in\mathbb Q\).
The same fact holds when one share is zero. The exact
[normalization for affine heights](MULTIVARIATE_OPPOSED_DIFFUSE_PUISEUX_HEIGHT_BRANCHES.md#two-groups-of-affine-row-types)
expresses the coefficient logarithms before their constant
terms are removed as sums of
\[
                  (c/u+s)\log(1+\kappa u p)
            \quad\text{or}\quad (c/u+s)\log(1+\kappa u q),
\]
where the rational number \(\kappa\) includes the affine offset.
The central coefficients themselves are algebraic, since
\(\exp(-\theta/(2nN))\) is algebraic. After centering \(p,q\),
the \(u^a\) terms of the logarithms have degree at most \(2a\)
in the single coefficient parameter \(\theta\). This is the
hypothesis of the
[one-clock auxiliary construction](MULTIVARIATE_OPPOSED_AUXILIARY_POLYNOMIAL_RATIONAL_BIAS.md#one-transcendental-coefficient-parameter).
Its formal lift, two Weierstrass preparations, Puiseux coefficient
bounds and integer auxiliary-polynomial argument therefore
apply unchanged to these factors. The isolated central fiber is
the same Gaussian one, since the two group directions are
independent. Thus an accumulating branch satisfies a nonzero
polynomial identity \(J(L,Y,r)=0\) over \(\mathbb Q\).

The common-factor minors vanish identically on such a branch.
The common degree is at least two because the exact polynomials
are real and the even moment is strictly positive on the real
axis. The field \(\mathbb C(L,Y,r)\) has transcendence degree
at most two by \(J=0\). The
[Ax rank theorem](MULTIVARIATE_OPPOSED_AX_BRANCH_RANK.md#a-bound-from-differential-transcendence)
excludes full rank for \(F_i\) on the branch. The
[one-relation argument](MULTIVARIATE_OPPOSED_AX_BRANCH_RANK.md#boundary-branches-need-at-least-two-character-relations)
uses only \(J=0\), the limiting value one for the characters,
and the degree-two common-factor locus of a pair of full
parents. It therefore also excludes exactly one relation here.

A nonzero integral relation for the characters in (1) can be
written as \(A(L)=B(Y)\). Indeed every character separates in
the two heights, and its limit is one. Each of \(A,B\) is
nonconstant: otherwise, since both heights vary, it would
give a forbidden relation in the two formal variables.
The same formal independence makes the mapping from relation
vectors to \(A\) injective modulo constants. Choose two
independent integral vectors and write their equalities as
\(A_1(L)=B_1(Y)\) and \(A_2(L)=B_2(Y)\). As functions of the
single formal variable \(X\), each of the pairs
\((A_1,A_2)\) and \((B_1,B_2)\) is multiplicatively independent
modulo constants. All zeros and poles of the \(A_i\)'s lie
in the fixed finite rational grid \(D_L\) from the ramp
factors of the first group. All zeros and poles of the \(B_i\)'s
lie in the analogous grid \(D_Y\). The grids include zero
if needed and exclude infinity because all factors tend to one.

## Enumerating the height curves

Within the field of branch germs take the common intermediate
field
\[
               \mathbb Q(A_1(L),A_2(L))
                   =\mathbb Q(B_1(Y),B_2(Y)).
\]
It is contained in \(\mathbb Q(L)\) and in \(\mathbb Q(Y)\).
By [Milne, Theorem 9.19 on Lüroth's theorem over the original field](https://www.jmilne.org/math/CourseNotes/FT.pdf)
it has one nonconstant rational generator \(W\). Thus
\(W=U(L)=V(Y)\) for rational \(U,V\) over \(\mathbb Q\);
there are \(H_1,H_2\in\mathbb Q(X)\) with
\(A_i=H_i\circ U,\ B_i=H_i\circ V\).
The union of the zeros and poles of the \(H_i\)'s contains at
least three points: on a set of at most two points the space of
degree-zero divisors has rank at most one. Multiplicative
independence of \(H_1,H_2\) follows from that of \(A_1,A_2\).

Take any three marked values in that union. For each value all
points of its \(U\)-fiber lie in \(D_L\), since a rational map
pulls back a nonzero divisor order to a nonzero divisor order.
All points of its \(V\)-fiber similarly lie in \(D_Y\).
The three marked values are rational points of \(\mathbb P^1\)
because each is the image of a rational grid point.
Apply the same rational fractional linear change to \(U,V\)
sending those values to \(0,1,\infty\).
For a degree-\(D\) map of projective lines, the total
ramification is \(2D-2\). Over three distinct values, the
total of fiber sizes counted with multiplicity is \(3D\).
The number of distinct points in the three fibers is therefore
at least \(3D-(2D-2)=D+2\). Hence
\[
                      \deg U\le |D_L|-2,\qquad
                      \deg V\le |D_Y|-2.                         \tag{2}
\]
If the full relation space has rank \(\nu\ge2\), using all its
generators gives at least \(\nu+1\) marked values simultaneously.
The same argument then bounds each degree by
\((|D_L|-2)/(\nu-1)\) or \((|D_Y|-2)/(\nu-1)\), respectively.

For a grid \(D_0\) the degree bound (2) makes the normalized
maps explicitly enumerable. Enumerate disjoint zero and pole
divisors on \(D_0\) with equal positive degrees at most
\(|D_0|-2\). Their monic polynomials \(P,Q\) determine the
map up to a multiplier \(v\), as \(vP/Q\). Choose one grid
point \(t\) in neither divisor as a preimage of one. This
forces \(v=Q(t)/P(t)\). Test that \(vP-Q\) has the full
degree and all its roots, counted with multiplicity, in \(D_0\).
Keep exactly these rational maps. The equality
\(U(L)=V(Y)\) on the branch is thus taken from a finite list
of explicitly known fiber products. Factor their cleared
equations over \(\mathbb Q\) and keep irreducible components
with both heights nonconstant. This is the claimed list.
The point \(X=\infty\) of the source is not needed
for the enumeration: it is simply outside the three specified
fibers, with any ramification there allowed.

## Fractional changes of the height ratio are excluded

On any of the algebraic candidate curves through
\(L,Y=\infty\), put \(t=1/L\) and
\(\xi=L/Y\). Choose a convergent Puiseux expansion in
\(s=t^{1/M}\) on the ordinary branch. Its coefficients are
algebraic. When both shares are positive, the limiting
\(\xi_*=q_*/p_*\) is positive rational at rational bias,
and at the boundary with \(L\) active it is zero. Suppose
\(\xi-\xi_*\) has first nonzero term \(\kappa s^b\) with
\(0<b<M\).

The role of the normalized count in the exact coefficients
is \(\sigma=r/L\). After stripping the common scale and
normalizing constants, the vector of coefficient logarithms
from its repeated blocks is \(\sigma H(s)\). At \(s=0\),
\(H(0)\) is \(H_L+\xi_*H_Y\), where \(H_L,H_Y\) are the two
rational Gaussian log directions in the normalized coefficients.
All contributions to \(H\) from the nonlinear terms of the
factorial logarithms are \(O(t)\) when \(\xi\) is bounded.
Consequently \((H(s)-H(0))/s^b\) is analytic, has algebraic
coefficients, and is \(\kappa H_Y\) at \(s=0\). The offset
multiplicities contribute an additional analytic log vector
with algebraic coefficients that vanishes at zero.

The limit \(\sigma_*=p_*\) is a nonzero rational multiple of
\(\theta\), hence transcendental. The coefficient prefactors
containing \(\exp(\sigma_*H(0))\) are algebraic. Apply the
[boundary clock deformation](MULTIVARIATE_OPPOSED_BOUNDARY_CLOCK_EXCLUSION.md#the-two-independent-gaussian-directions)
with \(\eta=\sigma-\sigma_*\) and the independent parameter
\(\zeta\) replacing \(\sigma_*s^b\). Here the two central
directions are \(H_L+\xi_*H_Y\) and \(\kappa H_Y\), which are
independent by the two-group hypothesis. The central
physical and off-diagonal quotient resultants therefore
have an isolated simultaneous zero, just as in the boundary
proof. Preparation in \(\eta\) applies because the Gaussian
strength changes nontrivially with \(\eta\): its coefficient
is \(A_1+\xi_* A_2>0\) in the notation of the two-group branch
theorem. Its norm gives a nonzero analytic germ \(E(s,\zeta)\)
with algebraic coefficients. Transcendence of \(\sigma_*\)
prevents \(E(s,\sigma_*s^b)\) from vanishing identically.
Thus the branch cannot have accumulating collisions.

On an interior branch, only
\(\xi=\xi_*+O(t)\) remains. The convergent Puiseux expansion
then gives \(Y=\rho L+\kappa_0+o(1)\), with
\(\rho=1/\xi_*\in\mathbb Q_{>0}\). Since \(L,Y\) are integers,
\(Y-\rho L\) belongs to a fixed discrete rational lattice.
On an accumulating sequence it eventually equals the rational
constant \(\kappa_0\). The irreducible algebraic height curve
containing that sequence is therefore the line
\(Y=\rho L+\kappa_0\). On the boundary the excluded fractional
changes correspond precisely to \(1<\beta<2\) in
\(Y\sim cL^\beta\).

## The algebraic exclusion test

For each candidate curve one may discard components without
a branch where both heights diverge. For a remaining component
\(C\) use its geometric function field. An accumulating
sequence forces \(C\) to be geometrically irreducible:
otherwise rational points would belong to pairwise
intersections of conjugate components and would be finite.
On the smooth compactification of \(C\), compute the
divisors of the restricted characters \(F_i\). This involves
only valuations of explicit rational functions at the
points over the finite grids. One can compute these using
algebraic Puiseux series, taking the normalization of the
curve into account. The integral kernel is the lattice of
character relations modulo constants. The constants in those
relations are one, by evaluation at a branch with both
heights divergent. The lattice is saturated since the divisor
group has no torsion.

An accumulating branch on \(C\) therefore has coefficients
in the translated coefficient subtorus defined by
\[
                \prod_i a_i^{m_i}=\prod_i G_i^{m_i}
                         \quad\text{for every lattice row }m.
                                                               \tag{3}
\]
If the lattice has rank at most one, the preceding Ax argument
excludes that curve. Otherwise one can use a basis in (3).
Impose degree-at-least-two gcd by the Sylvester rank minors
on the two coefficient lists, with endpoints nonzero.
This is a finite Laurent polynomial calculation over the
algebraic function field of \(C\). Gröbner elimination, for
example, computes the dimension after adjoining inverse
variables for the coefficient coordinates, so the resulting
dimension is taken in the torus. If the closed locus has
codimension at least two in (3), its exponential coordinates
contradict the codimension form of the Ax bound, which allows
at most \(t-1\le1\) for \(t=\operatorname{trdeg}_{\mathbb C}
\mathbb C(L,Y,r)\le2\). This is the advertised certificate.
A codimension-one locus requires further argument.

The [absolute terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
restricts the list and the test to finitely many orders for
the same grouped family. For two single opposed types the
grids have sizes \(2nd+1\) and \(2ne+1\); the list bound is
thus \(2nd-1,2ne-1\) on the rational maps. The
[pole-degree theorem](MULTIVARIATE_OPPOSED_BOUNDARY_POLE_DEGREES.md)
can additionally restrict the ratio of the two map degrees
at boundary contacts, and the boundary clock exclusion then
closes that entire regime for two single opposed types. The finite
curve certificate is stated more generally for affine groups. It uses the formal divisor
rank condition and a positive rational-bias contact.
