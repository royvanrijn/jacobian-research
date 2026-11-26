# Hermite splitting at simultaneous large factorial heights

Fix an ordinary two-monomial core \(C\), which may be empty. Append
finitely many independently variable smaller exponents \(A_i\),
with fixed positive gaps \(d_i\), fixed positive multiplicities
\(r_i\), and orientations \(s_i=1\) for an increasing row and
\(s_i=-1\) for a decreasing row. In an ordinary profile
\(A_i\in\mathbb Z_{\ge0}\). Write
\[
 D=\sum_i s_i r_i d_i,\qquad
 C_j(z)=Q_j^C(j^D z),\qquad j=n,N=n+1,\quad n\ge2 .
\]
Here \(Q_j^C\) is the factorial moment polynomial of the fixed
core, with its binomial coefficients included. For the empty
core it is \((1+z)^j\). The
[weak-tail angular theorem](MULTIVARIATE_WEAK_TAIL_MULTIPLE_ROOT_RESPONSE.md)
gives the earlier Hermite separation of noncentral branches.
Here the common physical scaling retains the signed twist
\(j^D\), and the exact polynomial expansion also applies
to noncancelling complex escape directions.

**Theorem.** Let \(g\) be the number of distinct complex roots
of the two comparators \(C_n,C_N\) that have odd multiplicity
in both. There is an effectively computable \(L\) such that,
if every \(A_i\ge L\), the ordinary consecutive moment
polynomials have a common factor of degree at most \(g\).
If the core is empty, \(g=0\) for every signed total gap \(D\).
This includes balanced tails with \(D=0\), where the limiting
comparators themselves have a full degree-\(n\) common factor.

The same bound holds along any complex sequence with all
\(A_i\to\infty\) and
\[
 \liminf\frac{|\sum_i r_i d_i^2/A_i|}{\max_i|1/A_i|}>0 .
                                                               \tag{1}
\]
Here the parents over the complex height space are defined
by finite rising factorials. In particular a branch at
infinity of the degree-at-least-two locus with all heights
polar and \(g<2\) must have cancellation among the leading
reciprocal heights if the two comparators are not coprime.
If the comparators are coprime, ordinary resultant
continuity excludes such a branch without (1).

**Consequences.**

1. For any two independently varying row heights with empty
   core, fixed gaps, fixed multiplicities and either
   orientation, all-degree coprimality is effectively
   decidable by a finite list of ordinary height pairs.
   The same is true with a fixed core if \(g\le1\) at every
   terminal below the absolute cutoff.
2. For three independently varying individual unit-gap rows,
   with empty core and any orientation pattern, there is
   an effective finite list of ordinary height triples
   containing every collision in any degree.
3. For the mixed orientations in part 2, the relevant
   components of the saturated degree-at-least-two collision
   locus have dimension at most one. Every geometrically
   irreducible curve among them that is defined over
   \(\mathbb Q\) satisfies the Runge pole condition over
   \(\mathbb Q\).

## Exact rescaling at a comparator root

Put \(\epsilon_i=1/A_i\) for positive \(A_i\).
Strip the factorial \((jA_i)!\) from each copy of a tail row.
The remaining row factor is \((jA_i+1)_{d_i k}\) for an
increasing row and \((jA_i+1)_{d_i(j-k)}\) for a decreasing
row at coefficient index \(k\). Divide further by the
scalar \(\prod_{s_i=-1}(jA_i)^{j d_i r_i}\) and substitute
the common argument scale \((\prod_i A_i^{-s_i d_i r_i})z\).
The rescaled parents are polynomials in \((\boldsymbol\epsilon,z)\).
They have the exact expression
\[
 F_j(\boldsymbol\epsilon,z)
       =\mathbb E[U_j^j C_j(zV_j/U_j)],\qquad
 U_j=\prod_{\substack{i,h\\s_i=-1}}Y_{j,i,h}^{d_i},\quad
 V_j=\prod_{\substack{i,h\\s_i=1}}Y_{j,i,h}^{d_i}.                \tag{2}
\]
The index \(h\) ranges from \(1\) to \(r_i\). The independent
\(Y_{j,i,h}\) have Gamma shape \(jA_i+1\) and scale \(1/(jA_i)\).
Only their moments are used. At any formal or complex
\(\epsilon_i\), define the same expression by
\(\mathbb E[Y_{j,i,h}^k]=\prod_{a=1}^k(1+a\epsilon_i/j)\)
and multiplicativity on separate copies. The parent in
(2) is polynomial because each monomial in \(C_j\) has
degree at most \(j\).

Write
\[
 \lambda^2=\sum_i r_i d_i^2\epsilon_i,\qquad
 \omega_i=\epsilon_i/\lambda^2,\qquad
 \sum_i r_i d_i^2\omega_i=1 .
\]
For ordinary positive heights take \(\lambda>0\). Its direction
\(\boldsymbol\omega\) lies in a compact real simplex. For (1)
choose either square root. Then \(|\lambda|\to0\) and the
complex directions remain bounded.

Let \(Z_{j,i,h}=(Y_{j,i,h}-1)/\lambda\). These centered
variables have polynomial cumulants in \(\lambda,\omega_i\).
The exact formulas are
\[
 \kappa_1(Z_{j,i,h})=\lambda\omega_i/j,\qquad
 \kappa_k(Z_{j,i,h})=(k-1)!\lambda^{k-2}
   (\omega_i/j)^{k-1}(1+\lambda^2\omega_i/j)\quad(k\ge2).
\]
For positive heights these follow by differentiating the
Gamma logarithmic moment generating function. They are
also polynomial identities for the prescribed complex
moments by the moment-cumulant formula. In particular
every fixed mixed moment of the \(Z_{j,i,h}\) is polynomial
in \((\lambda,\boldsymbol\omega)\). At \(\lambda=0\) the
limiting formal cumulants are those of independent centered
Gaussians of variance \(\omega_i/j\).

Let \(\zeta\) be any root of \(C_j\) of multiplicity \(m\);
it is nonzero because the core has nonzero constant term.
Write \(C_j(w)=\sum_{a=m}^j b_{j,a}(w-\zeta)^a\). At
\(z=\zeta(1+\lambda x)\), (2) gives the exact polynomial
identity
\[
 \lambda^{-m} F_j(\boldsymbol\epsilon,\zeta(1+\lambda x))
 =\sum_{a=m}^j b_{j,a}\zeta^a\lambda^{a-m}
  \mathbb E\!\left[
       U_j^{j-a}\left(xV_j+(V_j-U_j)/\lambda\right)^a\right].
                                                               \tag{3}
\]
Indeed \(V_j-U_j\) is formally divisible by \(\lambda\)
after substituting \(Y_{j,i,h}=1+\lambda Z_{j,i,h}\).
The polynomial in parentheses at \(\lambda=0\) is
\(x+\sum_{i,h} s_i d_i Z_{j,i,h}\), whose limiting cumulants
are those of \(x\) plus a centered Gaussian of variance
\(\sum_i r_i d_i^2\omega_i/j=1/j\). Thus the coefficientwise
limit of (3) is
\[
 b_{j,m}\zeta^m G_{j,m}(x),\qquad
 G_{j,m}(x)=
 \sum_{\ell=0}^{\lfloor m/2\rfloor}
 \binom m{2\ell}(2\ell-1)!!\,j^{-\ell}x^{m-2\ell}.
                                                               \tag{4}
\]
Use \((-1)!!=1\). The difference from (4) is uniformly
\(O(|\lambda|)\) on bounded direction sets, with effective
coefficient bounds from (3).

## No nonzero shared scaled Hermite nodes

The roots of \(G_{j,m}\) are \(-i/\sqrt j\) times the
roots of the probabilists' Hermite polynomial
\(\operatorname{He}_m\). They are simple. Schur's
irreducibility theorem for the nonzero squared Hermite
nodes applies to each parity. It is unchanged by the
rational scaling and sign in \(x^2\). The resulting monic
polynomial for the nonzero squared roots of \(G_{j,m}\)
has degree \(\lfloor m/2\rfloor\) and root sum
\(-m(m-1)/(2j)\).
[Schur, §§1(D) and 4](https://doi.org/10.1515/crll.1931.165.52)

Suppose \(G_{n,m}\) and \(G_{N,b}\) share a nonzero root,
where \(m\le n\) and \(b\le N=n+1\). The two irreducible
square-variable polynomials must then coincide, so
\(\lfloor m/2\rfloor=\lfloor b/2\rfloor\ge1\). If \(m=b\),
their root sums differ because \(n\ne N\). If \(m=2a\) and
\(b=2a+1\), equality of root sums would force \(2a-1=2n\),
which is impossible. If \(m=2a+1\) and \(b=2a\), both the
larger numerator and the smaller denominator give a
strict inequality. Thus the local limits share only
the simple root zero when both multiplicities are odd.

Choose disjoint small disks around the distinct roots
of both comparators, including one disk per shared root.
Coefficient convergence of the parents at
\(\boldsymbol\epsilon=0\) preserves each multiplicity
count in these disks by Rouché's theorem. In each shared
disk use (3) to choose much smaller disks in the
common \(x\)-coordinate around the distinct roots of
both local limits. These disks map into the original
\(z\)-disk for small \(|\lambda|\), and Rouché's theorem
puts one parent root in each of its local disks.
Their counts exhaust the original \(z\)-disk. Only
the disk for \(x=0\) can contain a root of both parents,
and it contains at most one root of either parent.
Summing over the shared disks proves the degree bound.

For effectivity all comparator roots and local limit
roots are algebraic. Choose isolating disks with
rational radii and lower bounds for the comparator
or local-limit modulus on each relevant circle by
polynomial factorization and root separation. Enforce
the Rouché inequalities using the sums of absolute
coefficients of the exact perturbation polynomials in
(2) and (3), with \(|\omega_i|\le1\) in the ordinary
simplex. This gives an effective threshold for
\(\lambda\), and
\(\lambda^2\le(\sum_i r_i d_i^2)/\min_i A_i\) gives \(L\).
For the empty core with \(D\ne0\), the two comparator
roots \(-n^{-D},-N^{-D}\) differ. If \(D=0\), their sole
shared root has multiplicities \(n,N\) of opposite
parities. In either case \(g=0\).

## Finite ordinary decisions

At a consecutive pair of ordinary orders the even
parent is positive for all real arguments. Indeed it
is the integral of an even power of a real nonzero
binomial against the positive factorial density.
Thus any common root at an ordinary height has its
distinct conjugate as another common root. If \(g\le1\),
a collision therefore forces some height below \(L\).

For two independent heights, fix that bounded height.
The remaining height is a single affine ray, with
all other rows fixed. The
[all-direction affine-ray theorem](MULTIVARIATE_AFFINE_BASE_RAY_ALL_DIRECTIONS.md)
gives an effective finite parameter bound on each
slice for each terminal. Apply the
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
and take the finite union over smaller terminals. The
order-one window is excluded by positivity of order two.

For three individual unit-gap rows, first bound one
height using the empty-core result. Exchange the two
monomials in that slice if needed to make the retained
row increasing. The other two have signed total gap
in \(\{-2,0,2\}\), so the
[one-retained-row velocity proposition](MULTIVARIATE_THREE_HEIGHT_EFFECTIVE_DECISION.md#one-retained-unit-row-with-tails-of-either-orientation)
applies and bounds one of their heights. Fix that
second bounded height and finish by the affine-ray
theorem on the remaining height.

## Runge for the mixed three individual unit rows

Work at a fixed adjacent pair, on the open chart of
full degrees and nonzero constant terms. Impose common
factor degree at least two by the Sylvester rank
condition. By a relevant geometric component we mean
one containing an ordinary integer point in this
chart. In a mixed orientation pattern one row has
the lone sign. Its height cannot be constant on any
positive-dimensional relevant component. Otherwise
it is an ordinary fixed height and the other two rows
are cooriented, contradicting the
[same-orientation two-height specialization for a fixed core](MULTIVARIATE_TWO_HEIGHT_FACTOR_SPECIALIZATION.md).
If the component has dimension at least two, its
projection to the lone-sign height has
positive-dimensional fibers generically, including at
some ordinary integer values. The same two-height
theorem excludes those fibers. Thus only curves
and finitely many zero-dimensional components matter.

Let a relevant curve be geometrically irreducible
and defined over \(\mathbb Q\). At least two heights
are nonconstant. Otherwise the other two fixed
heights are ordinary integers, and the affine-ray
resultant is a nonzero polynomial on the remaining
coordinate line. If all three vary, then they cannot
all be polar at one point on the smooth projective
normalization: the total signed gap \(D\) is odd,
and the empty-core comparator roots are distinct.
Full-degree resultant continuity excludes the
common-pole point.

If exactly two vary, their orientations are opposite
because the lone-sign height varies. The fixed row
has an ordinary integer height. Exchange monomials
if needed to orient it increasingly. Assume both
remaining heights are polar at one point. Their
reciprocal coordinates \(\epsilon_+,\epsilon_-\)
tend to zero. The one-row comparators have only
simple roots, by their second-order differential
equation. Any common roots of the nearby parents
must converge to common comparator roots, and at
most one can converge to each such root. There
are at least two distinct limits \(\zeta_1,\zeta_2\)
after passing to a subsequence through the
full-degree chart. Divide the two reciprocal
coordinates by their maximum modulus and pass to
a nonzero limiting direction \((u_+,u_-)\).
By formula (3) of the
[one-retained-row velocity proposition](MULTIVARIATE_THREE_HEIGHT_EFFECTIVE_DECISION.md#one-retained-unit-row-with-tails-of-either-orientation),
the first-order root differences at both limits
must vanish:
\[
 \frac{n^{-1}-N^{-1}}2(u_++u_-)
       +(1/n-1/N)\zeta_i u_- =0,\qquad i=1,2 .
\]
The same analytic expansion holds for complex
reciprocal heights, since these are simple roots.
Subtracting the equations gives \(u_-=0\), and
then \(u_+=0\), a contradiction. Thus the two
nonconstant heights have no simultaneous pole.

In both cases, on the projective normalization
the sets of labeled height coordinates polar at a
point cannot be identical at all points at infinity:
each nonconstant height has a pole, and no point
carries every nonconstant height. These labeled
sets are Galois invariant. There are therefore at
least two Galois orbits at infinity. Take an
integer linear combination of the heights with a
pole at every point at infinity by avoiding the
finitely many proper linear conditions for leading
pole cancellation. It is a rational function with
at least two distinct \(\mathbb Q\)-irreducible poles
and integer value on ordinary height triples.
Levin's curve form of Runge then gives effectivity.
Singular affine points, intersections of conjugate
components and zero-dimensional components give a
finite algebraic remainder.
[Levin, Theorem 1.2](https://www.numdam.org/item/10.5802/jtnb.634.pdf)

The three-row decision here concerns fixed gaps and
multiplicities as stated. The Hermite limit controls
all fixed gaps and multiplicities at once near
simultaneous escape, while more independent heights
require control on the slices with fewer escaping
types. For balanced mixed complex directions,
condition (1) identifies a separate possibility for
the pole analysis.
