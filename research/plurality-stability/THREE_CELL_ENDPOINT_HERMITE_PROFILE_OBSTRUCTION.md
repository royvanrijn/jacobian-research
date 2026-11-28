# A Hermite profile with two simplex contacts

Fix \(0<a<b<1\). The scalar endpoint inequalities, Parseval
normalization and nonnegative Hermite weights permit one profile
that exceeds the balanced simplex exactly on \((a,b)\). The
profile fits in a real three-function Gaussian \(L^2\) tuple
with balanced means, pairwise orthogonal coordinates, and an
orientation-independent covariance response. The tuple lacks
the indicator constraint. This gives an obstruction to an
endpoint proof using these spectral conditions alone.

## Exact normalization and contact signs

Let \(w_n=W_n(T)\) be the balanced simplex weights, so
\(S_z(T)=\sum_{n\ge0}w_nz^n\), \(w_0=1/3\), and \(\sum_nw_n=1\).
All \(w_n\) are strictly positive. Indeed the ray formula
\(S_z'(T)=6I_1(z)+3I_{-1/2}(z)\) shows that each coefficient
of the derivative is a positive multiple of \(6+3(-1/2)^j\).
For even powers use the binomial expansion of
\((1-z^2)^{-1/2}\), and for odd powers use its product with
\(\arcsin z\).
[Exact simplex derivative and ray formula (D1)–(D2)](THREE_CELL_DEFICIT_AND_PERIMETER.md#exact-correlation-derivative-for-the-three-cell-simplex)

With the branch of the square root positive for real \(z<1\), put
\[
 F(z)=z(z-a)(z-b)\sqrt{1-z}=\sum_{n\ge1}d_n z^n,
 \qquad V(z)=S_z(T)-\varepsilon F(z).                   \tag{EP1}
\]
For \(n\ge1\), the coefficient \(g_n\) of \(z^n\) in
\(\sqrt{1-z}\) is
\[
 g_n=-\frac{1}{2n-1}\binom{2n}{n}4^{-n}
     \sim -\frac{1}{2\sqrt\pi}n^{-3/2}.
\]
For \(n\ge4\) we have \(d_n=ab g_{n-1}-(a+b)g_{n-2}+g_{n-3}\).
Thus \(d_n\sim-(1-a)(1-b)n^{-3/2}/(2\sqrt\pi)\), so all but
finitely many \(d_n\) are negative. Choose \(\varepsilon>0\)
small enough that \(v_n=w_n-\varepsilon d_n>0\) for every
\(d_n>0\), and put \(v_0=w_0\). All coefficients of \(V\) are
then nonnegative. Absolute convergence at one and \(F(1)=0\)
give \(\sum_n v_n=1\). The first-degree loss is
\(w_1-v_1=\varepsilon ab>0\).

The real deficit \(S_r(T)-V(r)\) is positive for
\(0<r<a\) and \(b<r<1\) and negative for \(a<r<b\).
Its two contact derivatives are
\(\varepsilon a(a-b)\sqrt{1-a}<0\) and
\(\varepsilon b(b-a)\sqrt{1-b}>0\), respectively. Its limit
divided by \(\sqrt{1-r}\) as \(r\uparrow1\) is the positive
number \(\varepsilon(1-a)(1-b)\). The same construction has
strict eventual cumulative domination by the simplex, since
\[
 \sum_{n>N}(v_n-w_n)\sim
       \frac{\varepsilon(1-a)(1-b)}{\sqrt\pi\sqrt N}.      \tag{EP2}
\]
This is the sign and order required by the perimeter tail
comparison for a nonsimplex partition. For \(-1\le r<0\) one
has \(V(r)>S_r(T)\ge0\) because \(F(r)<0\), and
\(|V(r)|\le\sum_n v_n=1\).
[Hermite tails and perimeter (HC2)](THREE_CELL_HERMITE_CONTACT_BOUNDS.md#eventual-cumulative-domination-for-a-fixed-partition)

Arbitrary prescribed odd contact orders are possible. Replace
\(z(z-a)(z-b)\) by \(P(z)=z(z-a)^j(z-b)^k\) with positive
odd integers \(j,k\). The first coefficient is \(ab\) in the
simple-root case and \(a^j b^k\) in general. The sign pattern
on the real intervals is unchanged, and the coefficient
asymptotic becomes \(-P(1)n^{-3/2}/(2\sqrt\pi)\) by the
same finite convolution. In particular both contact slopes
can agree exactly with the simplex slopes.

## A Gaussian tuple and its matrix covariance response

Identify each real planar Hermite chaos with its symmetric
tensor space, with the isometric tensor normalization. In
degree \(n\ge1\) choose a two-dimensional subspace invariant
under planar rotations and reflections as follows. For odd
\(n\), take the harmonic linear polynomials multiplied by
\(|x|^{n-1}\) in the homogeneous tensor representation. For
even \(n\), take the harmonic quadratic polynomials multiplied
by \(|x|^{n-2}\). Use an orthonormal basis \(h_{n,1},h_{n,2}\)
of the corresponding Hermite subspace. The tensor polynomials
here describe the leading homogeneous terms, with lower terms
chosen by the chaos isometry.

Let \(u_1,u_2\) be an orthonormal basis of the label plane
\(\mathbf 1^\perp\subset\mathbb R^3\) and take the vector tuple
\[
 f=\frac{\mathbf1}{3}+
       \sum_{n\ge1}\sqrt{v_n/2}
                   (u_1h_{n,1}+u_2h_{n,2}).                \tag{EP3}
\]
The series converges in Gaussian \(L^2\). Its coordinate sum
is one and each coordinate has mean \(1/3\). Orthogonality
and \(\sum_{n\ge1}v_n=2/3\) give
\(\langle f_i,f_j\rangle=\delta_{ij}/3\). For each real scalar
correlation in \([-1,1]\) the transition pairings have diagonal
entries \(V(r)/3\) and off-diagonal entries \((1-V(r))/6\).
They therefore have the nonnegativity and row-sum properties
of a balanced label transition matrix. Its nonconstant
eigenvalues are \((3V(r)-1)/2\), both nonnegative for \(r\ge0\).

The entire tuple has only angular frequencies zero, one and
two on almost every radius. This follows by taking the
Gaussian \(L^2\) limit in the closed sum of the angular
frequency subspaces; the Hermite correspondence commutes
with rotations. Thus its angular sections have continuous
representatives. An indicator tuple would be constant on
each such circle because its values lie almost everywhere
in the three-point label set. That would contradict the
positive first-degree weight of (EP3). Moreover a tuple
with all entries in \([0,1]\) would already be an indicator
tuple, since \(\mathbb E f_i^2=\mathbb E f_i=1/3\).

For matrix covariance, each degree energy depends on the
coefficient Gram operator \((v_n/2)P_n\), where \(P_n\) is the
projection onto the two-dimensional subspace just chosen.
It commutes with rotations, so \(S_R(f)\) is invariant under
orthogonal conjugation of \(R\) and its first covariance
derivative at \(rI\) is scalar. The positive Hermite-tensor
calculations in (AM5)–(AM7) apply unchanged to this \(L^2\)
tuple: the eigenbasis weights are nonnegative, and on a
compact positive covariance domain the derivatives of every
fixed order are bounded multipliers. The trace-free covariance
Laplacian is strictly positive. For example its degree-three
input contains \(|x|^2\) times a nonzero linear harmonic mode,
whose Euclidean Laplacian is a nonzero linear polynomial.
Smoothing cannot remove that degree-one contribution. The
fixed-trace maximum-principle comparisons (AM8)–(AM9) also
hold for the envelope over this tuple, the simplex, and their
rotations, since the proofs use the same positive Laplacian
and an orientation-invariant comparison class.
[Covariance convexity and trace-free subharmonicity (AM5)–(AM9)](GAUSSIAN_ANISOTROPIC_MAXIMIZER_DIMENSION.md#logarithmic-covariance-comparison)

The example imposes neither pointwise indicator values nor
stationary winner equations. A proof for actual partitions
must use additional constraints beyond the two scalar interval
comparisons and these covariance identities.

## Accumulation and width of actual failure intervals

Write \(M(r)=\sup_A S_r(A)\) for the ordinary balanced Gaussian
optimum and \(d(r)=M(r)-S_r(T)\). The dimension theorem gives
planar maximizers, and the uniform operator Taylor bound makes
\(M\) continuous. A strict failure component is a connected
component \((a,b)\) of \(\{r:d(r)>0\}\). By the initial and
near-one results, actual failure components lie in a compact
subinterval of \((0,1)\) and both their endpoints are
simplex-optimal.
[Exact near-one optimality and the endpoint alternative](THREE_CELL_NEAR_ONE_EXACT.md#two-endpoints-of-any-strict-failure-interval)

Suppose distinct failure components accumulate at a scalar
correlation \(q\) in an interior compact band. Their lengths
tend to zero because their interiors are disjoint. In each
component choose \(r_j\) maximizing \(d\) on its closure and
choose a planar maximizer \(A_j\) there. The point \(r_j\) is
interior since the endpoint values of \(d\) are zero. For the
fixed shape \(A_j\), the analytic advantage \(S_r(A_j)-S_r(T)\)
is bounded by \(d(r)\) and equals \(d(r_j)\) at \(r_j\).
It therefore has a local maximum at \(r_j\), with zero first
derivative and nonpositive second derivative. The hypotheses
require only genuine interior improvement and endpoint
optimality, with no open interval of optimality adjacent to
an individual component.

Along a subsequence \(r_j\to q\), strong planar maximizer
compactness gives an ordinary maximizing limit \(A\) at
\(q\), and the adjacent simplex-optimal endpoints imply
\(d(q)=0\). Uniform local growth about the simplex orbit
excludes a simplex limit, because \(A_j\) improves on \(T\).
The first two correlation derivatives converge by strong
\(L^2(\gamma_2)\) convergence and the uniform differentiated
Hermite multiplier bounds on the band. Consequently
\[
 S_q(A)=S_q(T),\qquad S'_q(A)=S'_q(T),\qquad
                  S''_q(A)\le S''_q(T).                  \tag{EP4}
\]
Thus excluding nonsimplex tangential ties would already imply
that the strict failure set has only finitely many components.
[Maximizer compactness and local growth](GAUSSIAN_SIMPLEX_DIMENSION_FREE_LOCAL_GROWTH.md)

There is also a uniform bound on the width of a component
with a separated endpoint slope. Let \(K\) bound
\(|S''_r(A)-S''_r(T)|\) on the compact band, uniformly over
ordinary balanced partitions. For example one may take
\(K=2\sup_{n\ge2} n(n-1)u^{n-2}\) if \(u<1\) is the upper
correlation of the band. If a fixed competitor ties at \(a\)
with \(S'_a(A)-S'_a(T)\ge\eta>0\), then by Taylor's formula
and exact simplex optimality at the far endpoint,
\[
 0\ge S_b(A)-S_b(T)\ge\eta(b-a)-\tfrac12 K(b-a)^2 .       \tag{EP5}
\]
Hence \(b-a\ge 2\eta/K\). The same width bound holds if a
competitor ties at \(b\) with \(S'_b(A)-S'_b(T)\le-\eta\),
by applying the expansion toward \(a\). Both estimates need
only simplex optimality at the opposite endpoint and the
stated slope advantage toward the component. Optimality
accumulating on the exterior side of the given endpoint
forces that slope's weak sign by differentiation, while the
positive gap \(\eta\) is an additional condition.

## The abstract profile obstruction persists under compactness

Choose pairwise disjoint intervals \((a_j,b_j)\) in a compact
interior band, with \(a_j,b_j\to c\in(0,1)\). Use (EP1) with
the same \(\varepsilon>0\) for every interval. This is possible:
the coefficient ratios \(g_{n-k}/g_n\) tend to one for each
\(k=1,2,3\), so \(d_n/g_n\) tends uniformly in \(j\) to
\((1-a_j)(1-b_j)\), bounded away from zero. Thus the indices
of positive \(d_n\) lie in one common finite set, where their
sizes are bounded and the simplex weights are positive.

The perturbed weight sequences converge in \(\ell^1\) to the
weights of the profile with deficit
\(\varepsilon z(z-c)^2\sqrt{1-z}\). In fact multiplication of
the absolutely summable square-root series by the degree-three
polynomials is continuous in their coefficients. Using the
same Hermite bases in (EP3) gives strong \(L^2\) convergence
of the tuples, since \((\sqrt x-\sqrt y)^2\le |x-y|\) for
nonnegative weights. These tuples stay uniformly separated
from the simplex orbit: the latter has a nonzero fourth
angular Fourier component whose norm is independent of
rotation and relabeling. Every tuple in (EP3) has zero
fourth angular component. For instance the integral of
\(e^{4i\theta}\) on a \(120^\circ\) sector centered at angle
zero is nonzero.

Take the pointwise supremum over these profiles, their limit,
and the simplex profile. Its strict failure components are
exactly the prescribed intervals. The limiting profile ties
tangentially at \(c\) and never exceeds the simplex. Thus even
strong compactness, uniform separation from the simplex,
and a uniform first-degree loss among the abstract competitors
allow infinitely many failure components with an accumulation
tie of the type in (EP4).
