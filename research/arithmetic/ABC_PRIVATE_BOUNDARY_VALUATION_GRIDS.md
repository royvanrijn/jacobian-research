# Primitivity at the boundary of high-power components

A low-degree component cannot have a private affine zero of
a high-powered root. On primitive integer points, an arithmetic
Nullstellensatz then confines the prime support of that root.
This gives a determinant count with overlapping roots within
each letter and treats translated grids of exact valuations.

## PV1. Bounded releases with roots on the same letter

Fix \(m\ge2\). Take three nonzero signed terms
\[
       A_T \prod_{j\in I_T} x_j^{e_j},\qquad 1\le T\le3 ,
\]
with fixed integer coefficients and fixed positive integer exponents on active coordinates. Their sum is zero, the \(x_j\) are positive integers
at most \(H\), and there are at most \(m\) active coordinates.
Require that the supports of the moving variables on different
roles be disjoint. On each role designate a set of coordinates
of exponents at least \(E\). Every other active coordinate on
that role has support in the product of its designated roots.
There is no restriction on overlap within one role.

For \(E>4m\), the number of tuples is
\[
 O_{m,(e_j),E,\varepsilon}
       \left(H^{\,m^2(4m/E)^{1/(2(m-1))}+\varepsilon}\right).
                                                        \tag{PV1}
\]
The constant is uniform in the coefficients. With one active
coordinate the bound is constant, and with all active
coordinates on just one letter the divisor bound is available.
The exponent in (PV1) is explicit in \(m,E\); the other constants
are allowed to depend on \(m\) at each fixed budget.

With no active coordinates there is at most one tuple.
Otherwise the initial polynomial is nonzero and its geometric components
have dimension at most \(m-1\). In the latter one-letter case
any actual tuple bounds the fixed product by
\(H^{\sum_j e_j}\), so the coefficient uniformity of the divisor
bound follows. Otherwise work with the positive torus points on
its geometric components, using the following low-degree count.

### Private affine zeros force high degree

Let \(Z\) be a positive-dimensional geometrically integral
component or subvariety of projective degree \(d\), meeting the
open torus. Denote its affine part by \(Z_{\rm aff}\). Fix a
designated \(x_i\) of role \(T\) and put
\(P_{\ne T}=\prod_{j\notin I_T}x_j\).
Suppose \(Z_{\rm aff}\) contains a point with
\[
                     x_i=0,\qquad P_{\ne T}\ne0 .
\]
On the normalization of \(Z_{\rm aff}\), there is a prime
divisor with \(x_i\) vanishing and all other-role coordinates
units. Indeed this is the zero set of one nonzero regular
function on that open subset, and a component through the
chosen fibre has codimension one. All same-role exponents are
nonnegative. Thus the role-\(T\) monomial ratio \(F\) to either
other term, with signs chosen for \(1-F\), vanishes to order
at least \(E\) there.

Take a general linear section to a curve, as in the bounded
high-root release proof. A selected divisor meets this
section transversely on its smooth open part on the
normalization. One can avoid its intersections with the
other coordinate zero and pole divisors. The image curve
has degree \(d\), its normalization has genus at most
\((d-1)(d-2)/2\), and the orders of the three terms on the
selected branch have the asserted inequalities. The zeros
and poles of \(F\) and \(1-F\) lie in the union of the
coordinate hyperplane sections and the hyperplane at
infinity, with at most \(2md\) points in total. Riemann–Hurwitz
for \(F\) over \(0,1,\infty\) gives
\[
                        E\le\deg F\le 2g-2+2md<4m d^2 .
\]
Consequently, if \(d<D_0=\sqrt{E/(4m)}\), every designated
\(x_i\) has empty private affine zero set of this type.
[General sections and the three-value Riemann–Hurwitz argument](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#releasing-a-bounded-number-of-high-degree-coordinates)

### A fixed integer supports every root on a low-degree component

Use the rational polynomial space of degree at most \(d\)
vanishing at all height-\(H\) integer tuples under consideration
on \(Z\). A maximal independent set of evaluation rows gives
an integral basis of this rational space, by minors of size
bounded solely in terms of \(m,d\). The coefficients are
\(H^{O_{m,d}(1)}\). Let \(Y\) be their common zero set.
It lies in \(Z_{\rm aff}\). For any projective point outside
\(Z\), choose a linear space through that point of complementary
dimension minus one disjoint from \(Z\), and project from a
hyperplane in that linear space avoiding the point. The image
of \(Z\) is a hypersurface of degree at most \(d\) and the
external point lies outside it. Pullback thus separates any
external point by a form of degree at most \(d\).
The complex polynomials of that degree vanishing on the
integer sample lie in the complexification of the rational
evaluation kernel. This also proves the containment for
components defined over a nonrational field.

For each designated \(x_i\), append \(x_i=0\) and
\(zP_{\ne T}-1=0\) to the equations of \(Y\), with new variable
\(z\). The system has no common complex zero. The arithmetic
Nullstellensatz supplies an integer combination of these
polynomials equal to a nonzero integer
\(B_i\) with \(|B_i|\le H^{O_{m,d}(1)}\). For an actual tuple
and a prime \(p\mid x_i\), primitivity across roles gives
\(p\nmid P_{\ne T}\). Evaluate the certificate modulo \(p\)
using \(z=P_{\ne T}^{-1}\) to get \(p\mid B_i\).
The given height bound now permits only \(H^{o(1)}\)
choices of each designated root by fixed-support enumeration.
The other roots use only the union of these same supports
on their role, so they also have \(H^{o(1)}\) choices.
No bound on the size of the original coefficients \(A_T\)
was used in producing the certificates.
[Krick–Pardo–Sombra, Theorem 1](https://www.maia.ub.edu/~sombra/publications/arith-nss/maino.pdf)

For \(d\ge D_0\) use the finite birational projection and
Salberger covering in the bounded-root release proof. At
dimension \(j\ge1\) the new components of lower dimension
are paid at \(H^{(j+1)/d^{1/j}+\varepsilon'}\), with their
degrees bounded at fixed parameters. The upper bound on
this exponent apart from \(\varepsilon'\) is
\(m D_0^{-1/(m-1)}\), and at most \(m\) steps occur.
At any stage, the rational points of a geometric component
that is not defined over \(\mathbb Q\) lie in its intersection
with a distinct conjugate. Descend to the components of
that intersection before applying the cover, with degrees
bounded in terms of the original ones.
At dimension zero the number of points is bounded by degree.
Choose the temporary epsilon smaller than the final allowance
divided by the bounded number of operations. This proves (PV1).
[Salberger, Theorem 0.10 and the projection argument](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#releasing-a-bounded-number-of-high-degree-coordinates)

## PV2. Adaptive translated grids and growing numbers of grids

At height \(X\), divide the retained primes into disjoint
sectors on their original roles. In each sector allow an
arithmetic progression of valuations
\[
                         e=a+qt,\qquad t\ge0,
                      \qquad a,q\ge K .
\]
Both the choice of sectors and their orders may vary with the
row. On an ordinary divisibility sector one can just use
\(e=qt\) for positive \(t\). The exceptional prime powers
have combined radical at most \(X^f\). Use at most \(m\)
variables in total, charging a shifted sector two slots
and an ordinary divisibility sector one. Then at fixed \(m\) and \(K>4m\)
the primitive count is
\[
           O_{m,K,\varepsilon}
               \left(X^{f+C_m K^{-1-1/(2(m-1))}+\varepsilon}\right).
                                                          \tag{PV2}
\]

For each shifted sector put its prime support in a squarefree
head \(R\) and its positive digits in an integer \(W\) supported
on \(R\). Its full product is \(R^a W^q\).
Choose the exceptional full coefficients via their radicals
at cost \(X^{f+o(1)}\). For a cutoff \(T_0\ge K\), a factor
raised to an order greater than \(T_0\) can be fixed at
\(X^{1/T_0+o(1)}\) cost, enumerating the exponent up to
\(\log_2 X\) when the root is nonunit. A fixed head gives a containing
support for its digit. Record the full powered digit
value up to \(X\) at height-subpower cost even when the
order is unbounded. Fixing a digit alone
preserves the high head. Drop unit variables. The remaining
exponents all belong to \([K,T_0]\) and the roots are bounded
by \(H=X^{1/K}\). In (PV1) take \(E=K\). Altogether the
exponent is at most
\[
              f+\frac{m}{T_0}
                 +\frac{m^2}{K}(4m/K)^{1/(2(m-1))}
                 +\varepsilon .
\]
The choice \(T_0=\lceil K^{1+1/(2(m-1))}\rceil\) gives (PV2).
A grid of step \(q\) and residue \(0\le r<q\) with quotient
at least one is represented by \(a=q+r\). Thus the count
applies also to shifts relatively prime to their steps.
It retains the exact original valuations.

For example one shifted grid per role uses \(m\le6\). Take
\(K=L^{-\theta}\) with fixed \(10/11<\theta<1\).
Then \(LK\to0\) and the loss in (PV2) is \(o(L)\).
In particular, exception cap \(X^{(1/3-\delta)L}\) for
fixed \(\delta>0\) gives count at most
\(X^{(1/3-\delta/2)L+\varepsilon}\) for sufficiently
small fixed \(L\), without a purchased complementary
divisor. The same result allows a growing number of grids.
For
\[
       K=\left\lceil\frac{1}{L\log(1/L)}\right\rceil,\qquad
       m\log\log(1/L)=o(\log(1/L)),
\]
the explicit bound above with \(T_0=\lceil K^2\rceil\)
has both losses \(o(L)\). Indeed \(m/T_0=o(L)\) and the
negative power of \(K\) in the determinant term dominates
\(m^2\log(1/L)\). At every fixed \(L\), the values of
\(m,K,T_0\) are fixed before sending \(X\) to infinity.
Constants depending on this entire finite list are allowed.

The argument improves the bounded-grid residual estimate
for at least four disjoint ordinary root slots as well.
The three-single-root estimate of Browning–Verzobio is
stronger in its scope. A diffuse profile supported on prime
valuation indices can still avoid even these growing menus:
a grid of step at least \(K\) has \(O(1/(LK)+1)\) positions
in a fixed scaled band, while the prime degree indices
there number \(\asymp (L\log(1/L))^{-1}\).
Obtaining enough grid mass from the additive relation
and all three radical caps remains a separate extraction.
[Existing primitive Fermat and adaptive grid bounds](ABC_FINE_VALUATION_DIVISIBILITY.md)

## PV3. Aggregate private valuations on a component

The argument of PV1 can account for all designated roots
on one curve at once. Let \(Z\) be a geometrically integral
component of projective degree \(d\) in the signed monomial
equation, with \(m\) active coordinates. Let \(J_{\rm priv}\)
consist of the coordinates \(x_j\) of exponent \(e_j\) that
have a private affine zero on \(Z\), in the sense of PV1.
For each such coordinate choose one prime divisor on the
normalization where it vanishes and the other roles are
units. Take one general curve section that meets all the
chosen divisors and avoids intersections of distinct
boundary components. Coordinates of one role can vanish
on the same divisor. Their orders then add. Different
roles cannot share a private divisor.

On the normalization of this curve the signed ratios can
be written \(F\) and \(1-F\) with denominator the third
role. The multiplicities of its three fibres over
\(0,1,\infty\) sum to \(3\deg F\). Every chosen private
coordinate contributes at least \(e_j\) to this sum, with
same-role contributions at a marked point added before
the bound. For any nonempty private set, \(F\) is
nonconstant. Riemann–Hurwitz and the genus and boundary
bounds of PV1 then give
\[
                 \sum_{j\in J_{\rm priv}} e_j
                    \ \le\ C\bigl(d^2+(m+1)d\bigr) .      \tag{PV3}
\]
A zero-dimensional component meeting the open torus has
empty private set. For every remaining coordinate, the
arithmetic certificate from PV1 confines the primes on
the primitive sample to a prescribed polynomial-height
integer.

In particular consider the exact valuation roots
\(r_{T,e}=\prod_{p^e\parallel T}p\) of the original three
letters. Fix the full factors at orders \(e>B/L\), of
combined raw cost at most \(3/B\) in base \(n=X^L\), and
keep at most \(m\le 3\lfloor B/L\rfloor\) root coordinates.
Suppose a selected band \(\mathcal J\) has
\(e\ge a/L\), total actual raw mass at least \(\rho\),
and every one of its roots is at most \(n^\omega\).
The private mass of the band on \(Z\) is at most
\[
                (C L/a)\,\bigl(d^2+(m+1)d\bigr)\omega.
\]
The rest of that mass is supported on the fixed product
of certificates. Under total radical cap \(n^{B_0}\),
EB4 prices the surviving rows by at most one third of the
remaining outside radical cap plus an arbitrarily small
allowance at sufficiently small fixed \(L\). This gives a
strict sector saving on any sufficiently cheap catalogue of
bounded-degree components, with the polynomial equation
holding identically on each component. The chosen degrees,
ambient dimensions and certificate exponents are fixed at
each height limit. Their constants may depend on \(L\).
[Effective price of a prescribed prime support EB4–EB5](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## PV4. Dense rational parameter boxes and simultaneous common factors

The following obstruction allows a rational parameter map
whose image has large degree. Use original primitive shell
rows with three original caps and a total cap \(n^{B_0}\).
Fix the factors above \(B/L\), retain the \(m\) exact
valuation roots, and select a band \(\mathcal J\) as in
PV3. On a chart assume each retained root is
\[
              r_j=P_j(\mathbf t)/Q_j(\mathbf t)>0,
                  \qquad \mathbf t\in\prod_i I_i ,
\]
where the two integer polynomials are reduced over
\(\mathbb Q\), \(Q_j\ne0\) on the actual inputs and the
\(I_i\) are integer intervals. Translate the parameter intervals. Let \(n^v\) bound the box volume and let
\(n^\kappa\), with \(\kappa>0\), be a lower bound on
each side. The dimensions and degrees are bounded at
each fixed \(L\) before the height limit. Put
\[
 \mathcal D=\sum_{j=1}^m(\deg P_j+\deg Q_j),\qquad
 \mathcal K=
 \max_{\substack{J\subseteq\mathcal J\\
                 \sum_{j\in J}e_j\le3(\mathcal D+1)}}
                 \ \sum_{j\in J} b_j ,
\]
where \(r_j\le n^{b_j}\) on the chart. In particular
\(\mathcal K\le 3(\mathcal D+1)\omega L/a\).

For any \(\chi>0\) and fixed normalized allowance
\(\eta>0\), the count per chart is at most
\[
 O_{L,\eta,\varepsilon}\!\left(
    X^\varepsilon n^{\,\max\{v-\min(\chi,\kappa),
                          (B_0-(\rho-\chi-\mathcal K)_+)/3\}
                            +\eta}\right)              \tag{PV4}
\]
at sufficiently small fixed budgets. Count distinct actual
primitive rows, with all of their eligibility conditions.
Charge the chart catalogue and the cost of fixing the
over-cutoff factors in addition to this exponent. If the
rational additive equation fails identically as a formal
equality, the single bound \(v-\kappa\) applies.

For this last assertion clear denominators. The resulting
nonzero polynomial has bounded degree at fixed \(L\) and
vanishes on at most its degree times the volume divided by
the shortest side. This grid zero bound follows by induction
on the number of variables, separating the exceptional
values for the leading coefficient in the last variable.
A parameter vector determines the row from the fixed tail
and the exact root tuple.

The coefficient heights can be reduced before the
specialization argument. If the parameter sample has
size at most \(O_L(n^{v-\kappa})\), count it directly.
Otherwise evaluation of all polynomials of degree at
most twice the maximum displayed degree is injective
on the sample, by the grid zero bound. Solve the linear
equations \(P_j(\mathbf t)=r_j Q_j(\mathbf t)\) for each
slot in rational coefficients at the displayed degrees.
There is a nonzero solution and integer minors give a
solution of polynomial height, since the translated
parameters and the exact roots have polynomial height.
The denominator polynomial cannot be zero, again by
injectivity. The cross product with the original
presentation vanishes on the sample, so is identically
zero by injectivity at twice the degree. Reducing the
new fractions gives the same rational functions with
polynomial-height coefficients. Reduction preserves
the bounds on the degrees of reduced presentations.

In the identity branch factor all \(P_j,Q_j\) into primitive
irreducibles over \(\mathbb Z\), separating each content.
Include contents in a prescribed prime pool. These factors
have polynomial-height coefficients, uniformly at fixed
degree and dimension. For example substitute increasing
fixed powers of one variable to keep all monomials distinct
in factors and products, and use the univariate integer
factor height bound. Call a nonconstant irreducible
*exclusive* if it occurs in exactly one numerator and in
no denominator. Choose a generic affine line in the
parameter space preserving degrees of all the factors and
avoiding every common zero set of two distinct factors.
These sets have codimension at least two. On this line
the identity gives rational functions \(F,1-F\) on
\(\mathbb P^1\). The set of their zeros and poles has size
at most \(\mathcal D+1\), including the place at infinity.
An exclusive factor in a slot \(j\) gives a private fibre
contribution at least \(e_j\), as in PV3. Therefore
\(\sum_j e_j\le 3(\mathcal D+1)\) over the slots having an
exclusive factor, by Riemann–Hurwitz. If the ratios are
constant there are no exclusive factors. This bounds the
total selected-band price of these slots by \(\mathcal K\).
[The three-value Riemann–Hurwitz argument in the generalized Fermat setting](https://arxiv.org/html/2608.24512v1#S2.SS1)

For specialization at integer parameters, two distinct
factors can share a prime even though they share no
polynomial factor. The following count treats the union
of such collision primes. Translate the parameter
intervals and shear along an integer direction with
coordinate 1 on a shortest side. Choose a bounded
direction at the fixed degrees and dimension on which
none of the top homogeneous forms of the factors
vanishes. The condition is possible because their
restrictions to that coordinate being 1 are nonzero.
Use \(t\) on the shortest side, of size \(H_s\), and the
other coordinates \(\mathbf y\) as base variables. The
sheared box fits inside a product with side \(H_s\)
in \(t\) and with each base side of size between its
original size and a bounded multiple. All primitive
irreducible factors have constant nonzero leading
coefficient in \(t\). Enlarge the fixed pool by the
primes dividing these coefficients.

The product \(\mathcal R(\mathbf y)\) of the pairwise
resultants in \(t\) of the distinct factors is a nonzero
polynomial. The leading coefficients are constant and
the factors are coprime over \(\mathbb Q(\mathbf y)\).
The base points at which \(\mathcal R=0\) cost
\(O_L(n^v/H_s)\) after reinserting the \(t\) coordinate,
by the same grid bound. In one parameter the product
is a nonzero constant. At each other base point all
collision primes outside the fixed pool divide the
nonzero polynomial-height integer \(\mathcal R(\mathbf y)\).
Their exact squarefree product \(Q_{\rm col}\) has
height-subpower choices. For every prime in that product
the coordinate \(t\bmod p\) lies among a bounded number
of roots of the union of the factor polynomials.
All leading coefficients are units there. CRT gives
\(X^\varepsilon(1+H_s/Q_{\rm col})\) choices of \(t\)
for a prescribed \(Q_{\rm col}\), with allowances resized
to absorb the divisor choice. The branch
\(Q_{\rm col}\ge n^\chi\) therefore costs at most
\(n^{v-\min(\chi,\kappa)}X^\varepsilon\).

On the remaining branch consider a prime outside the
pool and all collisions. It divides at most one factor
value. If that factor occurs in a denominator it cannot
divide any retained root: the reduced fraction having
that denominator would violate integrality. A factor
in two distinct numerators and no denominator would
divide two exact valuation roots, contrary to their
pairwise coprimality. Hence the prime must be accounted
for by an exclusive slot if it divides any root.
The selected raw mass on the prescribed pool is therefore
at least \(\rho-\chi-\mathcal K\), since the selected
roots are pairwise coprime. Apply EB4 to the remaining
outside radical cap. Each pool has a polynomial-height
prime product uniformly at fixed \(L\), so the exact
dictionary prices and their height-subpower errors are
admissible for EB2–EB4. This proves (PV4).
[Effective root dictionary prices EB2–EB4](ABC_EFFECTIVE_RADICAL_BUDGET.md)

In the balanced limiting problem GE1 supplies nonatomic
mass at least \(1/2\). Choose compact positive-degree
bands of mass \(\rho\) arbitrarily close to this bound,
with boundaries of zero limiting mass. Inside the bands
discard root indices of large individual raw price.
For any fixed price threshold, the discarded indices in
each band can only converge to a finite sum of atoms.
A slow diagonal decrease of the threshold \(\omega\)
thus retains the required nonatomic mass. Take \(B\)
large enough for a prescribed small fixed tail price
\(\tau\) with \(\tau\ge3/B\). If the carrier degrees are
bounded per slot and \(\omega\to0\), then \(\mathcal K=o(1)\).
The following extraction suffices on such strata. Cover
the rows by \(n^b\) charts with total catalogue and
parameter volume \(b+\tau+v\le 1/3+o(1)\), short-side
bound \(\kappa>0\), and \(b+\tau<\rho/3\) by a fixed
margin. Choose
\(\chi=\min\{\kappa,(\rho-3b-3\tau)/4\}\).
Both terms of (PV4), with catalogue price included,
are then strictly below \(1/3\) by a fixed margin.
The extraction has to control the common parameter box
and the short-side bound, in addition to the rational
functions individually.

## PV5. Exact rational carriers from separate deep moduli

One way to supply the rational functions in PV4 is by
simultaneous polynomial congruences. Translate the common
parameter box so that each coordinate has magnitude at
most \(n^v\). Suppose at every actual point the retained
root \(r_j\le n^{C}\) is congruent to a polynomial of
total degree at most \(D\) in the common parameters,
modulo a fixed numerical modulus \(M_j\ge X^\nu\).
The moduli and polynomial coefficients can vary with
the slot and the chart. The degree, parameter dimension
\(r\), and bounds \(C,v,\nu>0\) can all be fixed before
\(L\) decreases.

Use the columns of all degree-at-most-\(D\) monomials
on the actual parameter sample. Choose a rational
basis of \(s\) columns, and an invertible
\(s\)-square evaluation minor. For each root the
augmented minors of these columns with the root
column are divisible by \(M_j\), by replacing the
root column by its polynomial combination of all
the monomials. Each replacement determinant is
zero even if the reduction of the basis loses rank.
The augmented minor size is at most
\((s+1)!\,n^{C+s Dv}\), hence less than \(M_j\) once
\(L(C+sDv)<\nu\) with a fixed positive margin and
\(X\) is large. This also gives a sufficient
condition for varying degree or dimension at the
fixed-budget height limit.

The minors vanish, so each root column is in the
exact rational span of the same degree-\(D\) columns.
Cramer's rule using the chosen minor gives a common
integer denominator of polynomial height and
numerators with polynomial-height coefficients.
These polynomials are valid on the entire sample,
and PV4 applies after reducing the fractions. The
denominator comes from one prescribed matrix minor
for the chart, with no extra catalogue price.
In this formulation the moduli need neither have
common prime support nor supply independent rank
directions. The open step for general endpoint
families is a cheap, common parameter model for the
exact valuation roots.

## PV6. Short polynomial parameters of arbitrary box volume

Keep the selected exact valuation roots and band of
PV3. Suppose the entire retained tuple is given by
integer polynomials of degree at most \(D\) in \(r\)
common integer parameters, each lying in an interval
of length \(H=n^C\). The coefficients can vary without
height bounds. Neither a density hypothesis in the
box nor a bound on its volume is needed here. For a
chosen degree cutoff \(d_*\ge2\), put
\[
 \xi_{\rm cov}=C r(r+1)d_*^{-1/r},\qquad
 \xi_{\rm live}=C_0(Ld_*^2+(BD+L)d_*)\omega/a
\]
with a sufficiently large absolute \(C_0\).
The count per chart is at most
\[
 O_{L,\eta,\varepsilon}\!\left(
 X^\varepsilon n^{\,\xi_{\rm cov}
                 +(B_0-(\rho-\xi_{\rm live})_+)/3+\eta}
                         \right).                     \tag{PV6}
\]
All parameters and degree bounds are fixed before the
inner height limit. The price of the tails and the chart
catalogue is added afterwards. At a bounded \(r,C,B,D\)
one can take \(d_*\to\infty\) slowly as \(\omega,L\to0\),
making both losses \(o(1)\). This holds even if \(D\)
grows, when \(D\omega\to0\) and \(d_*\) grows slowly
enough. If \(r\) also grows, the displayed losses specify
a sufficient condition.

Translate the parameter intervals, with the induced
translation of the polynomials. Intersect parameter
space with the pulled-back signed trinomial, or start
with the whole affine space if it vanishes identically.
The initial algebraic set has degree bounded at fixed
\(L,r,D,B\) uniformly in the coefficients. On a
geometrically integral component \(Z\) of parameter
degree \(d<d_*\), call a slot private if its polynomial
vanishes somewhere affine on \(Z\) with every other
coordinate polynomial nonzero. Choose a general curve
through the selected private divisors in the sense of
PV1. Its normalization has genus \(O(d^2)\). Each
coordinate polynomial or its reciprocal has at most
\(O(Dd)\) zero and pole places on this curve by
intersection with a hypersurface, including infinity.
Apply the three-fibre argument of PV3 to the signed
letter ratios. The sum of the exponents of private
slots is \(O(d^2+(mD+1)d)\). Thus the total private
selected-band mass is at most \(\xi_{\rm live}\).

For any other slot use the closed affine graph
\[
           \{(\mathbf t,\mathbf x):
               \mathbf t\in Z_{\rm aff},\
               x_j=P_j(\mathbf t)\ \text{for every }j\}.
\]
Its projective degree is at most
\(d\max(1,D)^{\dim Z}\). Indeed a generic system
of linear slices on the graph pulls back to equations
of degree at most \(\max(1,D)\) on \(Z\), and
Bézout bounds the number of isolated intersections.
Use the integers in the actual graph sample to
construct the bounded-degree evaluation kernel of
PV1. The translated parameters and exact roots both
have polynomial height. The graph has no private
affine zero for this slot. Append the slot zero and
an inverse of the product of the other exact roots.
The arithmetic certificate fixes the root's prime
support because all exact roots are pairwise coprime,
including those on one original letter. A component
of dimension zero gives a bounded fibre directly.
The resulting fixed pool captures at least
\((\rho-\xi_{\rm live})_+\) of the selected radical
mass, and EB4 gives the outside-budget price in (PV6).
The argument uses no size bound on the polynomial
coefficients defining the graph.

On components of degree at least \(d_*\), use the
coefficient-uniform determinant covers of PV1
in parameter space. At dimension \(j\le r\) the
price is at most \(H^{(j+1)d_*^{-1/j}+\varepsilon'}\).
There are at most \(r\) descents, with degrees
bounded at the fixed parameters. A geometric
component outside \(\mathbb Q\) is descended by
intersecting its distinct conjugates before using
the count on rational points. The total covering
price in base \(n\) is at most \(\xi_{\rm cov}+\eta'\).
This yields (PV6) after the low-degree count and
absorbing the fixed number of allowances.

For two polynomial parameters one can alternatively
factor the pulled-back nonzero equation and use
the coefficient-uniform Bombieri–Pila bound on
components of degree at least \(d_*\), with
\(n^{C/d_*+\eta}\) for that branch. If the equation
is an identity the whole parameter plane is already
a degree-one component. The low-degree proof applies
in either case. Thus two bounded-degree polynomial
parameters in arbitrary boxes of side \(n^C\) cannot
support balanced saturation on a catalogue of price
below \(\rho/3\) by a fixed margin, including its fixed
tail cost. This holds across the balanced dictionary
range. The restriction is on the arithmetic model
of the exact roots.
[Browning–Verzobio, Lemma 2.2](https://arxiv.org/html/2608.24512v1#S2.SS1)
