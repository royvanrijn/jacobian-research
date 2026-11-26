# Quadratic diffuse contacts for wide ordinary factorial cores

Fix any finite ordinary core, possibly empty, and append \(r\)
groups of rows with prescribed nonzero integer signed gaps and
prescribed positive integer multiplicities. The \(r\) smaller
exponents \(A_i\) vary independently. Let \(J\ge4\) be an integer
at least the [absolute terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
and suppose \(r\ge2J\).

**Theorem.** There is an effectively described proper rational
algebraic subset of the base space such that every ordinary
core outside it has the following properties.

- Every positive Gaussian physical reference contact of order
  \(3\le N=n+1<J\) in the allowed curvature and integer-gap
  range has exactly one nonreal conjugate pair of common roots.
  The roots are simple in both parents and their angular
  coincidence is transverse as the strength varies.
- There is a computable relative-gap tolerance for appended
  ordinary shared tails of arbitrary dimension and curvature.
  Below this tolerance, any exact ordinary collision, in any
  degree, has an irreducible quadratic gcd over \(\mathbb Q\).
  It occurs at positive curvature bounded away from zero and
  within \(O(\eta)\) of the finite list of reference contacts,
  where \(\eta\) is the maximum relative gap of the tail.

The tolerance and the coefficient bounds for the nearby
quadratic depend on the fixed core. The exceptional fraction
among the base vectors in \(\{0,\ldots,H\}^r\) is
\(O((H+1)^{-1})\). The variable part of the core may consist
entirely of unit gaps of one orientation. The fixed part may
contain private rows.

There is also a local consequence for every finite ordinary
core, with no lower bound on its row count. A specified simple
transverse interior conjugate-pair contact is a limit of
contacts of wider finite ordinary cores having exactly one
common conjugate pair. The added rows can have unit gaps,
zero total signed gap and arbitrarily small curvature.

## Rank charts with two free Gaussian parameters

For a terminal \(N=n+1\), \(N\ge4\), put \(D=n+N\), so
\(r\ge D+3\). Strip the endpoint factorial factors in each row
and normalize both coefficient arrays to constant term one.
Their consecutive binomial-stripped ratios are rational
coordinates \(y_{j,k}\), for \(1\le k\le j\) and \(j=n,N\).
Work in the complex open base set where all the sampled
rising factors are nonzero. The fixed core multiplies the
ratios by nonzero constants and does not change their
logarithmic derivatives.

The [paired coefficient rank proof](MULTIVARIATE_EXCEPTIONAL_CODIMENSION_TWO.md#differential-of-the-paired-coefficient-map)
gives the increasing-row column
\[
                  g_{j,k,d}(A)
                   =\sum_{\ell=(k-1)d+1}^{kd}(A+\ell/j)^{-1}.
\]
A decreasing row reverses its indices and changes the sign.
A fixed multiplicity multiplies the column by that integer
and does not change its annihilators. For \(d\ge2\) the row
functions span their full \(D\)-dimensional dual space.
The increasing and decreasing unit-gap columns have one
identically annihilating line each. Away from these two
annihilators, a projective kernel vector forces each base
to lie in a finite set. The rank-defective configurations
arising this way have dimension at most \(D-1\).
At one of the two distinguished annihilators only the
corresponding unit-gap bases can vary freely.

Write \(s_+,s_-\) for the two unit-gap group counts.
If \(s_+,s_-\le r-3\), the rank-defective set therefore has
dimension at most \(r-3\). Outside it the coefficient map
is a submersion to the full pair of coefficient arrays.
If some orientation has \(s\ge r-2\) unit-gap groups, reverse
both monomials if necessary so that this orientation is
increasing. Keep the other bases as coordinates \(B\), of
number at most two. For the \(s\) selected columns the rank
is \(D-1\) whenever they take at least \(D-1\) distinct
base values. Indeed the \(D-1\) distinct reciprocal
linear functions in these columns have linear combinations
with numerator degree at most \(D-2\). The excluded
diagonals have dimension at most \(D-2+2\le r-3\).
The map with \(B\) retained is a submersion to the
chart \(y_{N,N}/y_{n,n}=C(B)\) in the
[dominant unit orientation proof](MULTIVARIATE_EXCEPTIONAL_CODIMENSION_TWO.md#the-case-of-a-single-dominant-unit-orientation).
For fixed \(B\) it leaves all lower coefficients free and
imposes only one fixed ratio between the upper coefficients
of \(z^N,z^{N-1}\).

Introduce \(q,v\in\mathbb C^\times\) and an integer \(h\).
The two Gaussian reference polynomials are \(H_n(q,z)\)
and \(H_N(q,(N/n)^h v z)\), with
\[
                 H_j(q,z)=\sum_{k=0}^j
                    [z^k]Q_j(z)\,q^{(nN/j)k(j-k)} z^k .
\]
At fixed \(q,v\) the deformation multiplies each coefficient
by a nonzero constant. Both submersion charts and their
dimension bounds persist. In the dominant chart the upper
top ratio can depend on \(q,v\) and the lower coefficients,
while the coefficients at \(z,\ldots,z^{N-2}\) are free.

## Excluding extra and multiple common roots

At fixed complex \(q,v\), consider four distinct common
nonzero roots. This is impossible for \(n<4\). Otherwise
fix the lower coefficients and select four of their roots,
which gives finitely many choices. Three upper root conditions
already give independent linear equations on the first three
upper coefficient variables. Their evaluation determinant is
the product of the root arguments and a Vandermonde
determinant. These variables are free in both charts when
\(N\ge5\). The resulting locus in base space has dimension
at most \(r-3\).

Next consider two distinct common nonzero roots that are
both multiple in one parent. If the lower parent has the
two multiple roots, this requires four independent Hermite
conditions in its free coefficients at fixed roots. The
two root variables therefore cost at least two dimensions
in the lower coefficient chart. Requiring the upper parent
to vanish at these roots costs two more dimensions, using
its free coefficients at \(z,z^2\). For \(n<4\) the two
lower multiple roots are impossible. If the upper parent
has the two multiple roots, fix the lower coefficients
and choose the roots from them. For \(N\ge5\), three of the
four upper Hermite conditions are independent on the free
coefficients at \(z,z^2,z^3\). To see this, a polynomial in
these three monomials vanishing twice at one nonzero root
and once at the other would give a quadratic with three
zeros after division by \(z\). The same proof works in the
full quartic chart. In the dominant quartic chart the upper
polynomial is uniquely fixed by the two prescribed double
roots and its constant coefficient. The upper chart has
three free coefficients, so the codimension is at least
three there too.

Allowing the two clock parameters increases dimension by
at most two. The rank-chart exclusions and the two loci
just considered therefore project to subsets of base space
whose closures have dimension at most \(r-1\).
At a real physical contact with positive strength the even
parent has no real roots, by [positive lognormal averaging](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#the-finite-set-of-physical-reference-contacts).
Extra common roots would give four distinct roots once
the multiple cases are removed. A multiple nonreal common
root in either parent makes its conjugate multiple in that
parent as well. Thus every physical contact outside the
exclusions consists of a single simple conjugate pair.

## Transversality from the root-pair incidence

Consider the complex incidence of two distinct simple
nonzero common roots \(z_1,z_2\) of the two parents over
the good rank charts. Locally solve the lower equations
for \(z_1,z_2\) by simplicity. The two upper equations
are independent because the upper variables at \(z,z^2\)
are free with the lower coefficients held fixed. The
incidence is smooth of dimension \(r\), starting from the
\(r+2\) base and clock variables after eliminating the lower
root equations.

Project the incidence to the \(r\) base variables. An
irreducible component with nondense image contributes a
proper exceptional set. On any dominating component the
map is generically étale. Indeed its function field is a
finite extension of the base function field, separable in
characteristic zero, so the differential has full rank at
the generic point. The nonsurjective differential locus is
therefore a proper algebraic subset of that component and
its image has dimension at most \(r-1\). This is the
equal-dimension case of generic smoothness. It also follows
from the [separable generic-point criterion](https://stacks.math.columbia.edu/tag/07ND).

For the physical interpretation take the upper root branches
\(w_i(q)\) before the multiplier \((N/n)^h v\), and lower
branches \(\zeta_i(q)\). Put \(r_i(q)=w_i(q)/\zeta_i(q)\).
In the solved incidence equations
\(r_i(q)-(N/n)^h v=0\), the determinant in \(q,v\) is
nonzero exactly when \(r_1'(q)\ne r_2'(q)\). For a real
physical conjugate pair the two derivatives are conjugate,
and \(r_1(q)\) is positive real at the contact. Therefore
this determinant is nonzero exactly when the derivative of
the angular difference is nonzero. The real Gaussian
strength is \(-2nN\log q\) with nonzero derivative, so the
two notions of transversality agree.

All the loci used above can be specified by polynomial
equations and inequations over \(\mathbb Q\), using implicit
differentiation and clearing nonzero denominators for the
last determinant. Their closures can be computed by
elimination. The proof gives a proper algebraic exceptional
set for each fixed \(n,N,h\), without specializing the
positive real clock parameters. A simultaneous reversal
only changes \(h,v\) into \(-h,v^{-1}\) and preserves the
conclusions.

## The curvature endpoints and exact ordinary gcd

The lower endpoint \(\tau=0\) has \(q=v=1\).
For any fixed \(h\), the two-distinct-root argument on the
same charts gives codimension at least two at this fixed
clock, using the upper coefficients at \(z,z^2\).
The bad rank-chart configurations have codimension at least
three. Thus adding these algebraic sets excludes physical
contacts at zero strength. For real ordinary cores, both
here and at positive strength, a common root would have
to come with its nonreal conjugate. Terminal two is excluded
by this real positivity. Terminal three has no Gaussian
physical contact, by the [endpoint ratio argument](MULTIVARIATE_DIFFUSE_RATIONAL_SLOPE_GENERIC_CLOSURE.md#a-single-terminal-integer-twist-and-rational-slope).

The integer range can be taken independently of the base
variables because the total unsigned core gap \(L\) is fixed.
The [center bound](MULTIVARIATE_DIFFUSE_RATIONAL_SLOPE_GENERIC_CLOSURE.md#the-finite-union-and-the-ordinary-tail)
allows \(|h|\le L+N(n+N+24)\) for a contact of strength
at most \(12\), or an actual collision with tail curvature
at most \(12\) and relative gaps at most \(1/2\).
Take the finite union of the proper exceptional sets at
these integers and terminals. A nonzero rational polynomial
vanishes on the union, which gives the integer-grid bound.

Now fix an ordinary core outside the union. The reference
resultants at zero strength are nonzero for every allowed
integer and terminal. Rational interval evaluation on a
small compact strength range therefore excludes ordinary
collisions there for all sufficiently diffuse tails, using
the [effective coefficient and center approximation](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#uniform-approximation-of-the-physical-pair).
If the tail curvature is at least \(12\), the
[absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
already excludes exact collisions. The eventual cutoff
excludes the high terminals. On the remaining compact
positive strength range, the [finite contact confinement theorem](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md)
places every sufficiently diffuse collision in one of
finitely many designated contact neighborhoods.

At such a contact choose small disjoint conjugate disks
around the common pair. Each disk contains one root of
each reference parent, with multiplicity one, and there
are no common roots outside them. By compact root
separation and Rouché these facts give an open coefficient
neighborhood in which any common roots must belong to
those disks and each disk has at most one root of either
parent. All these neighborhoods can be certified by root
isolation and coefficient enclosures at the algebraic
reference arrays. Hence an ordinary collision consists
of exactly the two nonreal simple roots. Its monic gcd is
rational and quadratic, and irreducible because both roots
are nonreal.

For an effective rational coefficient comparison, normalize
the common argument by the lower tail first-factor scale
of the [projective factor sieve](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md#consequences-at-finite-physical-diffuse-contacts).
The true reference pair in this normalization is
\(H_n(q_*,q_*^{-N(n-1)}z)\) and
\(H_N(q_*,(N/n)^h v_*q_*^{-N(n-1)}z)\).
The [simple transverse confinement bound](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#a-local-bound-when-the-contacting-roots-are-simple)
gives a parameter error \(O(\eta)\). The same Rouché
estimate for simple roots gives an \(O(\eta)\) error in
each common root and thus in the two coefficients of its
monic quadratic. The constants can be bounded from the
isolating disks, the nonzero angular derivative and the
coefficient and center errors. This completes the effective
all-degree reduction for the fixed core.

## Keeping an arbitrary transverse contact through the extension

Given a specified simple transverse interior conjugate-pair
contact of any finite ordinary core, add an even number
\(r\ge2N+2\) of independent unit-gap groups with multiplicity
one, half in each orientation. The total signed gap of
the added core rows is zero. At the specified orders and signed gap, the algebraic
incidence proof applies with the original core fixed and
with both clock parameters unrestricted. Thus outside a
proper algebraic subset every physical reference contact of
the widened core at these orders and signed gap has only a
single simple conjugate pair. The complement
contains ordinary integer choices with every added base
over any given height, by Zariski density of each such
integer grid.

The exact added core rows have curvature tending to zero
uniformly as the least base tends to infinity. Their
endpoint-normalized log multipliers and the correction to
their relative center are \(O(\sum_i A_i^{-1})\), by the
ordinary factorial chord and center bounds. With the
original finite core fixed, these changes make the
combined Gaussian reference coefficients, including their
physical scale, tend to the original ones. For the specified
simple pair, the two real equations of the physical contact
have nonsingular Jacobian in \(q,v\), as shown above.
The implicit function theorem for the coefficient
perturbation therefore gives a nearby real physical contact.
The original cone inequalities are strict, so they remain
strict. The combined curvature and total integer gap also
converge to the given values, with the total gap preserved
exactly. At every generic widened core this nearby contact
has the asserted quadratic-only form. The argument concerns
reference contacts and does not assert persistence of an
exact ordinary collision.

## Excluding prescribed quadratic coefficient loci

The algebraic exceptional set can be enlarged to handle any
fixed proper rational algebraic subset \(\mathcal Z\) of the
coefficient plane of monic quadratics. In the common rational
normalization above, every sufficiently diffuse ordinary collision
for a generic core then has its quadratic gcd outside \(\mathcal Z\).
This holds uniformly in the added rows for the fixed core.

Indeed at fixed \(q,v\) consider two distinct nonzero lower
roots \(z_1,z_2\). The map from these two roots to the two
coefficients of their monic quadratic in the
\(q^{-N(n-1)}\) normalization has finite fibers. Requiring
the quadratic point to be in \(\mathcal Z\) therefore leaves
at most one dimension of possible ordered root pairs.
Each prescribed pair costs two independent lower coefficient
conditions, and then two independent upper coefficient
conditions. The two coefficients at \(z,z^2\) are free
for each of those operations on the rank charts, taking the
lower polynomial first. At fixed \(q,v\) the resulting
codimension is at least three, so the projection after
allowing the two clocks is still a proper algebraic set.
A finite collection of prescribed loci can be imposed at
the same time. At a core outside the enlarged exceptional
set the finite limiting quadratics lie outside the closed
locus. Their effective coefficient neighborhoods therefore
exclude actual factors in that locus.

For example write \(z^2-sz+t\) for the common quadratic.
Its reciprocal trace is \(s^2/t-2\). A bound \(L_0\) on
the reduced denominator of this trace allows only finitely
many rational values in \([-2,2]\), since its roots are
nonreal conjugates and their ratio has modulus one.
Take for \(\mathcal Z\) the union of the corresponding
curves \(s^2=(2+\theta)t\).
The generic-core theorem then excludes every sufficiently
diffuse ordinary extension for which the exact quadratic
denominator budget is at most \(L_0\). The budget can be
supplied by the
[projective valuation identity](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md#exact-projective-denominators).
In particular \(L_0=1\) covers common quadratic slope
packages with equal valuations of the two roots at every
prime. This argument still permits the reduced denominator
to grow without bound along other ordinary profiles.
