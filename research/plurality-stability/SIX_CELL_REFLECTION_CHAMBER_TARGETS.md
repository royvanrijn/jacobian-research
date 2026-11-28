# Reflection-chamber targets for the simplex standard kernel

The [five-cell chamber argument](FIVE_CELL_A3_BRIDGE.md) compares
averages of reflected Mehler kernels at the same noise parameter
under a radius-preserving change of chamber.

**Theorem (Six-cell reflection target).** In six cells, \(A_4\)
is the only possible finite real reflection chamber for the
unweighted pointwise comparison (CT2) under the stated chamber
homeomorphism. Angular volume and the diagonal image
multiplicities give the restriction.

## The comparison and its necessary conditions

Let \(k=n+2\), \(H=\mathbf1^\perp\cap\{x_1=x_2\}\), and
\(C=\{x\in H:x_2\ge x_3\ge\cdots\ge x_k\}\). Use the normalization
in (A34) of the five-cell proof, with \(k\) in place of five.
Against conditional Gaussian probability on \(C\), its native kernel
is a common factor
\[
 \sigma^{-n}\exp\!\left[-\frac{\rho^2}{2\sigma^2}
                              (|x|^2+|y|^2)\right],
 \qquad \sigma^2=1-\rho^2 ,
\]
times the average of \(\exp(\rho x\cdot gy/\sigma^2)\) over
\[
 \{g:g(1)=1\}\;\sqcup\;
 \{g:g(1)=1,\ g(2)=2\}.
 \tag{CT1}
\]
Indeed, in (S49) of the [interface reduction](SIMPLEX_INTERFACE_SYMMETRY.md),
the \(n!\) folded self terms occur twice and the \(n\,n!\) folded
adjacent terms occur once. The latter and one copy of the self terms
make up the subgroup on positions \(2,\ldots,k\). The total size in
(CT1) is \((n+2)n!\).

Suppose \(W\) is a finite real reflection group on \(\mathbb R^n\),
with closed chamber \(C_*\), and \(T:C\to C_*\) is a radius-preserving
homogeneous homeomorphism. Assume the native kernel dominates
\[
 \frac1{|W|}\sum_{g\in W} M_\rho^{(n)}(Tx,gTy)
 \tag{CT2}
\]
pointwise at some fixed \(0<\rho<1\), where \(M_\rho^{(n)}\) is the
standard \(n\)-dimensional Mehler kernel. No comparison of individual
summands is assumed. The common factor cancels in (CT2).
Use \(x=ru\) and \(y=rv\) for unit vectors in \(C\) and let \(r\to\infty\).
The leading native exponent is \(u\cdot v\), by rearrangement of the
ordered coordinates \(2,\ldots,k\). In the target it is
\(Tu\cdot Tv\). For \(Tu\) in the interior, an orbit point outside
\(C_*\) violates a chamber wall, and reflecting across that wall
strictly increases its dot product with \(Tu\). Thus the maximum
is at the orbit representative \(Tv\), with boundary cases following
by continuity. Taking the logarithm of the finite positive sums,
dividing by \(\rho r^2/\sigma^2\), and passing to the limit gives
\(u\cdot v\ge Tu\cdot Tv\). Hence the inverse on the unit section is
\(1\)-Lipschitz in spherical distance, and the target spherical area
is at least the native one.

Write \(p_n\) for the Gaussian mass of the full reference facet
before ordering the \(n\) untied coordinates. The native spherical
fraction is \(p_n/n!\), whereas the target fraction is \(1/|W|\).
Consequently
\[
                         |W|\le n!/p_n .                     \tag{CT3}
\]
For \(u=v\) unit, division by the common maximal exponential and
the same limit show a second necessary condition. If the consecutive
equal-score blocks of \(u\) on positions \(2,\ldots,k\) are
\(b_1,\ldots,b_s\), then
\[
 \frac{|\operatorname{Stab}_W(Tu)|}{|W|}
 \ \le\ 
 \frac{(b_1+1)(b_1-1)!\prod_{j=2}^s b_j!}{(n+2)n!}.           \tag{CT4}
\]
Equality in \(u\cdot gu\le1\) means that \(g\) fixes \(u\).
The first multiset in (CT1) then contributes \(b_1!\prod_{j=2}^s b_j!\)
and the second contributes \((b_1-1)!\prod_{j=2}^s b_j!\).
At a generic diagonal point the right side of (CT4) is
\(2/[(n+2)n!]\).

## The chamber restriction

The facet probability is \(p_n=\mathbb E Q^n\) for
\(Q=\Phi(Z/\sqrt2)\), with standard normal \(Z\). Indeed a normal
of variance \(1/2\) represents the tied coordinate, independently of
the remaining \(n\) standard normals, and removing the common trace
does not affect their ordering. Set \(U=Q-1/2\). Symmetry and the
two-dimensional Gaussian orthant formula give
\(\mathbb E U^2=\arcsin(1/3)/(2\pi)\): the two Gaussian differences
from the tied coordinate have correlation \(1/3\). The orthant formula
also follows by measuring the angle of the intersection of two
halfplanes for an isotropic planar Gaussian. Since \(|U|\le1/2\),
\[
 \frac{25}{176}
 <\frac1{16}+\frac1{4\pi}
 <p_4
 =\frac1{16}+\frac32\mathbb E U^2+\mathbb E U^4
 \le \frac1{16}+\frac{7}{8\pi}\arcsin(1/3)
 <\frac16 .                                                   \tag{CT5}
\]
For the lower bounds use \(\pi<22/7\) and
\(\arcsin(1/3)>1/3\). For the upper bound use
\(\arcsin(1/3)\le 1/(2\sqrt2)\) by integrating the increasing
derivative of arcsine, followed by \(\pi>3\) and \(7<5\sqrt2\).
Also \(p_{j+1}\ge p_j/2\), because the increasing functions
\(Q^j\) and \(Q\) have nonnegative covariance and \(\mathbb E Q=1/2\).
So \(p_n>2^{-(n-1)}\) for every \(n\ge4\).

Coxeter's classification and order table leave only \(A_n\) among
irreducible targets in rank \(n\ge4\). Every other irreducible
rank-\(n\) order is at least \(2^{n-1}n!\), which is too large
by (CT3). [Coxeter, *The Complete Enumeration of Finite Groups*,
Theorem 1 and table on p. 25](sources/coxeter-finite-reflection-groups-1935.pdf)

In six cells \(n=4\), (CT3) and (CT5) give \(|W|<180\).
The four native unit rays with exactly two score blocks have
right sides of (CT4), in order of increasing \(b_1\),
\[
                           \frac13,\quad\frac18,\quad
                           \frac19,\quad\frac5{24}.
\]
All other unit points have at least three score blocks and
right side at most \(1/12\), by the same factorial formula.
A target with a nonzero fixed subspace would have a unit point
of stabilizer fraction one. A reducible essential target with
an \(A_1\) factor has a unit point along this factor with stabilizer
fraction \(1/2\). Both violate (CT4). The only other reducible
essential cases are \(I_2(a)\times I_2(b)\), \(3\le a\le b\).
The product chamber has four unit rays with fractions
\(1/a,1/a,1/b,1/b\), since a boundary ray in one rank-two factor
is fixed by its wall reflection and by the whole other factor.
If \(a\le4\), two distinct rays would require two native points of
fraction at least \(1/4\), and only one exists. If \(a\ge5\), the
order \(4ab<180\) requires \(b\le8\). The four preimages would then
all require fraction at least \(1/8\), and only three native points
qualify. This excludes the reducible targets.

Thus \(A_4\) is the sole target surviving (CT2) for six cells.
Moreover (CT5) makes the six-cell worst diagonal row bound
\(2/(6p_4)>2\). An \(A_4\) comparison using just its degree-two
reference energy and that row bound cannot yield a positive gap,
since \((1-\rho^2)/(1-\rho)=1+\rho<2\). A transport and a sharper
row estimate would still be needed within this proof route.
