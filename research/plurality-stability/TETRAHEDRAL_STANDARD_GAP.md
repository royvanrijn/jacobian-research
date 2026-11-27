# The balanced tetrahedral interface gap through correlation 23/25

For every \(0<\rho\le23/25\), the tetrahedral standard sector of the
Gaussian simplex interface form has gap at least \(\ell/1000\). Here
\(\ell\) is the distance between two of the four Gaussian cell-moment
vectors. The constant-coefficient comparison below reaches \(17/25\).
The [variable-coefficient comparison](#extension-using-the-variable-local-coefficient)
extends it to \(23/25\). The other tetrahedral sectors are covered by
[(S56), (S60) and (S61)](SIMPLEX_INTERFACE_SYMMETRY.md#four-cells-domination-settles-the-remaining-sign-block).
It follows that for every fixed correlation in this range and every
ambient dimension \(d\ge3\), the balanced full-interface form has a
strictly positive gap off the geometric rotations, with a gap constant
independent of \(d\) by [(S42)](SIMPLEX_INTERFACE_SYMMETRY.md#every-cell-count-the-active-gap-lifts-to-all-ambient-dimensions).
The remaining standard-sector correlations are \(23/25<\rho<1\).
The interface result gives the local comparison under
[(B10)–(B11)](POSITIVE_CORRELATION_CONTINUATION.md) in the proved range.
A distant nonsimplex global tie is a separate open obligation.

## Ordinary balance and the tetrahedral angle

Use \(E=\mathbf1^\perp\subset\mathbb R^4\) and the reference facet
\(F_{12}\). An orthonormal basis of \(E\) is
\[
 p=(1,1,-1,-1)/2,\quad q=(0,0,1,-1)/\sqrt2,\quad n=(1,-1,0,0)/\sqrt2.
\]
Write \(x=up+\sqrt2\,tq\) on \(F_{12}\), so \(u\ge|t|\).
In the orthonormal coordinates \((u,\sqrt2\,t)\) the wedge has opening
\(2\alpha=\arccos(-1/3)\), with \(\alpha=\arctan\sqrt2\).
Equivalently, its two unit rays are the normalizations of
\((1,1,1,-3)\) and \((1,1,-3,1)\), whose dot product is \(-1/3\).

Put \(\phi_0=(2\pi)^{-1/2}\), \(w=\nu_{12}(F_{12})=\phi_0\alpha/\pi\)
and \(\mu=\nu_{12}/w\) on \(F_{12}\). Gaussian divergence gives
\(\ell=4w\) as in [(S77)](SIMPLEX_INTERFACE_SYMMETRY.md#the-exact-initial-standard-spectral-minimum).
The standard multiplicity space consists of functions even in \(t\).
The permutation group has order \(24\). Its unoriented-facet stabilizer
has order \(4\): \((12)\) fixes the facet pointwise and reverses the
oriented speed, while \((34)\) reflects \(t\). The standard scalar
interaction on the full reference wedge has weights \(2,1,1\) from
[(S49)](SIMPLEX_INTERFACE_SYMMETRY.md#every-cell-count-the-standard-sector-is-a-scalar-positive-kernel).
Using the full wedge with probability measure \(\mu\) accounts for the
facet mass without a factor of \(24\).

Let \(K_\rho\) on \(L^2(\mu)\) denote integration against the kernel
\(2M_\rho(x,y)+M_\rho(x,(23)y)+M_\rho(x,(24)y)\) on the even subspace.
Here \(M_\rho\) is the three-dimensional Mehler kernel relative to
Gaussian measure. Let \(P\) be orthogonal projection off the constants
on that subspace and define
\[
  A_\rho=\frac{P K_\rho P}{4\rho},\qquad \sigma_\rho=\sqrt{1-\rho^2}.
\]
The balance constraint is \(\int h\,d\mu=0\). If
\(q_{\rm std,\rho}=\mathfrak q_\rho(f^{v,h})/(4|v|^2)\) in (S51), its
lower bound from [(B7)](POSITIVE_CORRELATION_CONTINUATION.md) is
\[
 q_{\rm std,\rho}(h)
 \ \ge\ \ell\left(\sigma_\rho^{-1}-\|A_\rho\|\right)
                   \|h\|_{L^2(\nu_{12})}^2.
 \tag{T1}
\]
The factor \(4\) in \(A_\rho\) is from \(\ell=4w\). Dividing both norm
and interaction by \(w\) converts the comparison to \(L^2(\mu)\).
This uses ordinary mass balance before making the constant lower bound
on the variable local coefficient \(a_\rho/\rho\).

## Finite Hermite coefficients and the tails

The Hermite expansion [(S11)](SIMPLEX_INTERFACE_SYMMETRY.md) of the full
facet interactions restricts to the standard sector and then to its
balanced subspace. In particular \(A_\rho=\sum_{m\ge1}\rho^{m-1}C_m\)
with \(C_m\ge0\). Each coefficient is the restriction of a sum of squares
of the degree-\(m\) Hermite moments of the signed cell boundary measures.
The constant coefficient vanishes after \(P\). This proves monotonicity
of \(A_\rho\) in quadratic form order. If \(A_{\rho,N}\) is the sum
through \(m=N\), then for \(0<\rho<b<1\)
\[
 0\le A_\rho-A_{\rho,N}\le(\rho/b)^N A_b,\qquad
 \|A_\rho-A_{\rho,N}\|\le
  (\rho/b)^N\,\frac{\pi}{2\alpha b\sigma_b}.
 \tag{T2}
\]
The first inequality bounds each omitted power separately. The positive
scalar kernel at \(b\) has row sum \((a_b/b)/w\) on \(\mu\). Symmetry
and the upper bound \(a_b/b\le2\phi_0/\sigma_b\) in (B7) give the Schur
bound after division by \(4b\).

For the full Hilbert–Schmidt norm a second bound is
\[
 0\le\|A_\rho\|_{\rm HS}-\|A_{\rho,N}\|_{\rm HS}
      \le(\rho/b)^N\,\frac{\pi}{\alpha b\sigma_b^3}.
 \tag{T3}
\]
In fact the squared Mehler kernel at \(b\), times the two full-plane
standard Gaussian densities, has a Gaussian exponent with precision
blocks \((1+b^2)\sigma_b^{-2}I\) on the diagonal and
\(-2b\sigma_b^{-2}D\) off diagonal, where \(D\) is the compression of
the corresponding permutation to the pair plane and \(\|D\|\le1\).
For each singular value \(s\) of \(D\), the determinant factor is
\(((1+b^2)^2-4b^2s^2)/\sigma_b^4\ge1\). The full-plane
Hilbert–Schmidt norm is therefore at most \(\sigma_b^{-3}\).
Conditioning both variables to the wedge increases this upper bound by
at most \(\pi/\alpha\). Apply the triangle inequality to the weighted
Mehler kernels and then compress to obtain
\(\|A_b\|_{\rm HS}\le\pi/(\alpha b\sigma_b^3)\).
For positive Hilbert–Schmidt operators \(0\le X\le Y\), one has
\(\|X\|_{\rm HS}\le\|Y\|_{\rm HS}\), since the difference of squared
norms is \(\|Y-X\|_{\rm HS}^2+2\operatorname{tr}(X(Y-X))\ge0\).
This proves (T3) using the first inequality of (T2).

To compute \(C_m\), use the polynomials
\(H_{j+1}^{(v)}(z)=zH_j^{(v)}(z)-jv H_{j-1}^{(v)}(z)\), with
\(H_0^{(v)}=1\) and \(H_1^{(v)}(z)=z\). In the coordinates
\((p\cdot x,\sqrt2\,q\cdot x,\sqrt2\,n\cdot x)\), their variances are
\(1,2,2\). The coordinates of \(x\) on the facet are \((u,2t,0)\).
The coordinates of \((23)y\), where \(y=u'p+\sqrt2\,t'q\), are
\((t',u'+t',u'-t')\). The \((24)\) contribution is the conjugate of
the \((23)\) contribution by \((34)\).

For each \(i+j+l=m\) with \(j,l\) even, form the polynomial product
\[
 \frac{H_i^{(1)}(u)H_j^{(2)}(2t)H_l^{(2)}(0)}
      {2^{1+j+l}i!\,j!\,l!}
 \left[
 H_i^{(1)}(u')H_j^{(2)}(2t')H_l^{(2)}(0)
 +H_i^{(1)}(t')H_j^{(2)}(u'+t')H_l^{(2)}(u'-t')
 \right].
 \tag{T4}
\]
The parity restriction on the first variable removes odd \(j\), and
\(H_l^{(2)}(0)=0\) for odd \(l\). Keep only the even part in \(t'\).
The factor \(2^{-1}\) includes the two identical contributions in
each of the self and adjacent terms after even compression and division
by \(4\). The remaining powers of \(2\) are from the variances in the
Mehler expansion. Centering each monomial, which kills all constant
terms in either variable, produces the exact rational coefficient
matrix of \(C_m\). As low-degree checks, the degree-one matrix has
coefficient \(1/2\) on the centered \(u\) pair. In degree two the
centered \((u^2,t^2)\) matrix is
\(\left(\begin{smallmatrix}1/4&1/4\\1/4&5/4\end{smallmatrix}\right)\).

## Ball Gram matrix and the Rayleigh bound

For even \(j\), the \(\mu\)-moment of \(u^i t^j\) is
\[
 \frac{2^{i/2}}{\alpha}\Gamma(1+(i+j)/2)
       \int_0^\alpha \cos^i\theta\sin^j\theta\,d\theta .
 \tag{T5}
\]
Here direction and radius are independent, the folded angle is uniform
on \([0,\alpha]\), and the radius has density \(r e^{-r^2/2}\).
For the angular integral \(I_{ij}\), integration of
\(\frac{d}{d\theta}(\cos^{i-1}\theta\,\sin^{j+1}\theta)\) gives
\((i+j)I_{ij}=(i-1)I_{i-2,j}+\cos^{i-1}\alpha\,\sin^{j+1}\alpha\)
for \(i\ge2\). For \(i=1\), \(I_{1j}=\sin^{j+1}\alpha/(j+1)\).
For \(i=0\), \(I_{00}=\alpha\) and integration of
\(\frac{d}{d\theta}(\cos\theta\,\sin^{j-1}\theta)\) gives
\(j I_{0j}=(j-1)I_{0,j-2}-\cos\alpha\,\sin^{j-1}\alpha\).
The base moments include
\(\mathbb E_\mu u=\sqrt\pi/(\sqrt3\alpha)\) and
\(\mathbb E_\mu u^2=1+\sqrt2/(3\alpha)\).

Let \(G\) be the Gram matrix of the centered monomials and \(Q\) the
rational coefficient matrix of \(A_{\rho,N}\) from (T4).
The action on polynomial coefficient vectors is \(J=QG\), so the
squared Hilbert–Schmidt norm is \(H^2=\operatorname{tr}(J^2)\).
The trial coefficient vector \(z\) is obtained by applying \(J\) fourteen
times to the centered \(u\) coordinate. Its Rayleigh quotient \(r\)
and squared residual \(d^2\) are
\[
 r=\frac{z^\top GJz}{z^\top Gz},\qquad
 d^2=\frac{(Jz)^\top G(Jz)}{z^\top Gz}-r^2.
\]
In particular the trial quotient for the full \(A_\rho\) is between
\(r\) and \(r\) plus the operator tail in (T2).

If \(s=\sqrt{H^2-r^2}<r\), the top eigenvalue \(\lambda_1\) of the
truncation is at least \(r\) and all other eigenvalues are at most \(s\),
since their squared sum is at most \(H^2-r^2\). Integrating
\((\lambda_1-x)(x-s)\le0\) against the spectral measure of the normalized
trial gives
\[
 \lambda_1\le r+\frac{d^2}{r-s}.
 \tag{T6}
\]
This is the one-sided Temple argument used in
[(DG42)](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md). The isolation check
\(s<r\) is part of the certificate. For classical context see
[Kuroda, Theorem A, p. 506](https://www.kurims.kyoto-u.ac.jp/~prims/pdf/43-2/43-2-24.pdf),
which requires spectral isolation and a residual condition. Heilman's
[Hyperstable Sets, Theorem 1.5](https://arxiv.org/pdf/2209.11216) gives
the three- and four-cell global conclusion under a hyperstability
hypothesis on maximizing partitions. The planar-facet calculation here
uses the prescribed tetrahedral geometry and supplies its local
second-variation estimate.

## Certified domain and replay

The [checker](check_tetrahedral_standard_gap.py) generates (T4) over
the rationals and evaluates (T5)–(T6) at 256 bits using Arb ball
arithmetic. The upper bound on each row uses the operator tail (T2),
with \(N=32,\ b=19/20\). The first row uses \(H\) for the truncated
operator and the later rows use (T6).

| Left endpoint \(a\) | Right endpoint \(c\) | Strict upper bound \(U\) on \(\|A_c\|\) |
|---:|---:|---:|
| \(0\) | \(1/2\) | \(19/20\) |
| \(1/2\) | \(3/5\) | \(113/100\) |
| \(3/5\) | \(16/25\) | \(31/25\) |
| \(16/25\) | \(33/50\) | \(13/10\) |
| \(33/50\) | \(67/100\) | \(133/100\) |
| \(67/100\) | \(169/250\) | \(269/200\) |
| \(169/250\) | \(17/25\) | \(339/250\) |

For each row the checker verifies \(U+1/1000<(1-a^2)^{-1/2}\).
For every \(\rho\in[a,c]\cap(0,1)\), monotonicity gives
\(\|A_\rho\|\le\|A_c\|\) and \(\sigma_\rho^{-1}\ge(1-a^2)^{-1/2}\).
Thus the seven endpoint calculations prove (T1) with the stated gap
throughout the interval. At \(17/25\), (T3) bounds the full
Hilbert–Schmidt norm between \(1.8275268\) and \(1.8300908\).
The full trial quotient is between \(1.3553863\) and \(1.3555114\).
The upper bound on the operator norm from (T6) and (T2) is less than
\(1.355512\).

Run from the repository root with Python 3.12.3 and python-flint 0.9.0:

~~~sh
.venv/bin/python research/plurality-stability/check_tetrahedral_standard_gap.py
~~~

The [verification record](tetrahedral_standard_gap_verification.json)
retains the ball enclosures and replay parameters. This is an internal
computer-assisted proof and replay. Independent mathematical review and
formal verification remain open.

## Extension using the variable local coefficient

The constant in (T1) is the local coefficient at the wedge vertex. Away
from the vertex, both gaps to the untied coordinates can increase the
coefficient. Retaining that increase gives a much better comparison.
Here and below the scale is \(D_\rho=\sigma_\rho a_\rho/(\rho\ell)\).
In \(L^2(\mu)\), the scaled standard form divided by \(w\ell\) is
\[
 \int D_\rho h^2\,d\mu-\langle h,\sigma_\rho A_\rho h\rangle_\mu
 \qquad \left(\int h\,d\mu=0\right).
 \tag{TW1}
\]

### A reciprocal bound from the two facet gaps

Put \(q_0=\ell/(2\phi_0)=2\alpha/\pi\) and
\(g_\pm=\rho(u\pm t)/\sigma_\rho\). Both gaps are nonnegative on the
facet. In the cone integral (B7), the \(n\) coordinate has tilt
\(z_N\,1_{\{z_N>0\}}\). After normalizing this integral, let \(N\) be
Rayleigh with density \(n e^{-n^2/2}\) for \(n>0\), and let \(U,T\)
be independent centered Gaussians of respective variances \(1,1/2\),
independent of \(N\). The shifted cone inequalities in the two untied
coordinates give
\[
 q_0D_\rho(x)=
 \Pr\{U-T+N/\sqrt2\ge-g_-,\ U+T+N/\sqrt2\ge-g_+\}.
 \tag{TW2}
\]
In particular \(D_\rho\ge1\), and \(D_\rho\) is pointwise
nondecreasing in \(\rho\). Both statements follow from monotonicity in
the two gaps and the value \(q_0\) at zero shifts, or directly from
(B7).

Write \(p(g)=\Pr(U-T+N/\sqrt2<-g)\). The joint success probability
in (TW2) is at least \((1-p(g_-))(1-p(g_+))\). To see this, represent
\(U-T,U+T\) in law as \(S/\sqrt2+E_-,S/\sqrt2+E_+\) with independent
standard normals \(S,E_-,E_+\). Conditional on \(S,N\), the success
probabilities are two increasing functions of \((S+N)/\sqrt2\).
Their covariance is nonnegative by the identity expressing it as
half the expected product of differences at two independent copies.
This supplies the required association without an additional
Gaussian inequality assumption.

For the standard Gaussian upper tail \(\overline\Phi\) and \(x,y\ge0\),
\(\overline\Phi(x+y)/\overline\Phi(x)\le 2\overline\Phi(y)\).
Indeed the ratio is the tail above \(y\) of the density proportional
to \(e^{-xv-v^2/2}\) on \(v\ge0\), whose extra decreasing factor
\(e^{-xv}\) only reduces the tail from its value at \(x=0\).
Condition on \(N=n\) in \(p(g)\) and take \(x=n/\sqrt3\),
\(y=\sqrt{2/3}\,g\). The unshifted marginal is
\[
 p(0)=\int_0^\infty n e^{-n^2/2}\overline\Phi(n/\sqrt3)\,dn
     =\frac12-\frac1{\sqrt{6\pi}}\int_0^\infty e^{-2n^2/3}\,dn
     =\frac14
\]
by integration by parts. It follows that
\(p(g)\le \tfrac14\,2\overline\Phi(\sqrt{2/3}\,g)\). The latter
Gaussian tail has the polar formula
\[
 2\overline\Phi(\sqrt{2/3}\,g)
   =\frac2\pi\int_0^{\pi/2}e^{-g^2/(3\sin^2\theta)}\,d\theta .
\]
This is obtained by imposing the coordinate cutoff on a standard
two-dimensional Gaussian, integrating the Rayleigh radial tail, and
using symmetry. The integrand increases with \(\theta\).
Using the right ends of the three equal angular intervals gives
\[
 Z(g)=\tfrac13\big(e^{-4g^2/3}+e^{-4g^2/9}+e^{-g^2/3}\big)
       \ \ge\ 2\overline\Phi(\sqrt{2/3}\,g).
\]
On \(0\le z\le1\), \((1-z/4)^{-1}\le 1+z/3\).
The reciprocal envelope used in the certificate is therefore
\[
 \frac1{D_\rho(x)}
 \ \le\ E_\rho(x):=q_0(1+ Z(g_-)/3)(1+ Z(g_+)/3).
 \tag{TW3}
\]
All of its terms are Gaussian tilts along the two facet gaps.

### Uniform spectral comparison on a correlation interval

There is a direct tail bound using the entire pair plane. For one
Hermite degree \(m\), let \(v_m\) be the vector of orthonormal
ambient Hermite moments of \(h\,d\mu\). The same-facet quadratic
coefficient in \(K_\rho\) is \(\|v_m\|^2\) with weight two.
Each adjacent coefficient has the form
\(\langle v_m,Rv_m\rangle\) with weight one, where \(R\) is the
orthogonal action of the corresponding transposition on degree \(m\).
Its absolute value is bounded by \(\|v_m\|^2\), so division by the
normalizing factor four of \(A_\rho\) leaves a tail bounded by the
same-facet tail.

To bound the latter, extend \(h\) by zero to the pair plane with its
standard Gaussian probability measure. For a planar Hermite degree \(k\),
the unnormalized ambient Mehler tail has multiplier
\[
 \rho^k\sum_{j\ge0:\ k+2j>N} a_j\rho^{2j},\qquad
 a_j=\frac{\operatorname{He}_{2j}(0)^2}{(2j)!}
     =\frac{1}{4^j}{2j\choose j}.
\]
Indeed the coordinate normal to the plane is zero, so odd normal
degrees vanish and the planar degree is unchanged. The sequence
\(a_j\) is decreasing. Reindexing from the first included \(j\)
bounds the multiplier by
\(\rho^{N+1}\sum_{j\ge0}a_j\rho^{2j}=\rho^{N+1}/\sigma_\rho\).
The wedge probability in the plane is \(p_F=\alpha/\pi\).
The two facet measures contribute \(p_F^{-2}\) in the interaction
and the squared norm of the zero extension contributes \(p_F\).
After the remaining factor \(1/\rho\) in \(A_\rho\) and the scale
\(\sigma_\rho\), this proves the operator bound
\[
 0\preceq\sigma_\rho (A_\rho-A_{\rho,N})
    \preceq \frac{\rho^N}{p_F} I .
 \tag{TW4}
\]
The lower bound also follows from the positivity of the coefficients
in (T2). Centering by \(P\) restricts all these bounds to the ordinary
mean-zero space.

The coefficient of the full-plane bound in a low planar degree is
strictly smaller than the scalar cap in (TW4). After scaling by
\(\sigma_\rho/\rho\), the exact multiplier of that bound is
\[
 \lambda_k(\rho)=\rho^{k-1}
    \left(1-\sigma_\rho\sum_{j<J}a_j\rho^{2j}\right),
 \qquad J=\max\{0,\lfloor(N-k)/2\rfloor+1\}.
\]
For \(k\ge1\) this is nondecreasing. For \(J=0\) this follows
from the factor \(\rho^{k-1}\). For \(J\ge1\), differentiating
the bracket and using \(2(j+1)a_{j+1}=(2j+1)a_j\) leaves
\((2J-1)a_{J-1}\rho^{2J-1}/\sigma_\rho\ge0\).
The degree-zero planar moment vanishes by mass balance.
Consequently on a slab \([a,c]\), one can subtract from the scalar
cap the terms
\[
 \sum_{k=1}^N (c^N-\lambda_k(c))
       \sum_{i+j=k}\left|\int h\,
       \frac{\operatorname{He}_i(u)\operatorname{He}_j(\sqrt2t)}
            {\sqrt{i!j!}}\,d\mu\right|^2.
\]
Each term is nonnegative by the bound used in (TW4).
Here the two factors of \(p_F\) in the planar moments of the zero
extension cancel the \(p_F^{-2}\) conversion of measures. Odd \(j\)
vanishes by reflection symmetry. The checker subtracts rational
lower bounds on these deficits, obtained by rounding the multipliers
\(\lambda_k(c)\) upward.

On a slab \([a,c]\), the coefficient of \(C_m\) in
\(\sigma_\rho A_\rho\) is \(F_m(\rho)=\sigma_\rho\rho^{m-1}\).
Its maximum is at an endpoint or at
\(\sqrt{(m-1)/m}\) if that point lies in the interval.
For \(m\le N\), let \(f_m\) be a rational upper bound on that maximum.
The checker rounds the square-root factors upward to precision
\(10^{-12}\). Thus
\(\sigma_\rho A_\rho\preceq K_{a,c}+\tau I\), where
\(K_{a,c}\) is \(\sum_{m=1}^N f_m C_m\) minus the certified planar
deficit terms and \(\tau=c^N/p_F\).

Fix the target gap \(\gamma=1/1000\). Since \(D_\rho\ge D_a\ge1\),
write \(e=E_a/(1-\tau-\gamma)\) whenever \(1-\tau-\gamma>0\).
Then \(D_\rho-\tau-\gamma\ge e^{-1}\).
It is enough to show that \(K_{a,c}\), evaluated on ordinary
mean-zero \(h\), has norm below \(1\) relative to
\(\|h\|_{L^2(e^{-1}\mu)}\).

For a monomial \(v_i\) from the even basis in (T4), the functional
\(h\mapsto\int h v_i\,d\mu\) in this weighted space is represented on
the mean-zero subspace by \(e(v_i-\int ev_i\,d\mu/\int e\,d\mu)\).
The same functional results from the ordinary-centered monomial
when \(\int h\,d\mu=0\). The appropriate Gram matrix is exactly
\[
 G_{ij}^{(e)}
 =\int e v_i v_j\,d\mu-
   \frac{\int e v_i\,d\mu\int e v_j\,d\mu}{\int e\,d\mu}.
 \tag{TW5}
\]
If \(Q\) is the symmetric rational coefficient matrix of \(K_{a,c}\),
the finite matrix \(J=QG^{(e)}\) has the nonzero eigenvalues of the
self-adjoint weighted operator. Thus for any positive even integer \(p\),
\(\|K_{a,c}\|_{\rm weighted}\le(\operatorname{tr}J^p)^{1/p}\).
Here \(p\) is chosen from powers of two so the trace is checked by
successive matrix squarings. This extends spectral validation to slabs
where the ordinary Hilbert–Schmidt norm cannot isolate one eigenvalue.

### Evaluation of the weighted moments

Every term in (TW3) is of the form
\(e^{-k_-(u-t)^2-k_+(u+t)^2}\), up to a constant. On the folded
wedge \(0\le t\le u\) the combined Gaussian has precision entries
\[
 P_{11}=1+2(k_-+k_+),\qquad P_{22}=2+2(k_-+k_+),\qquad
 P_{12}=2(k_+-k_-).
\]
With \(A=\sqrt{P_{11}}\) and \(B=\sqrt{P_{22}-P_{12}^2/P_{11}}\),
change coordinates to \(z_1=Au+(P_{12}/A)t,\ z_2=Bt\).
The resulting Gaussian is isotropic and the angular interval is
\([0,\beta]\) with \(\beta=\arctan(B/(A+P_{12}/A))\).
The denominator is positive because \(P_{11}+P_{12}=1+4k_+>0\).
The transformed measure is \( (\sqrt2/(\alpha AB))r e^{-r^2/2}dr\,d\theta\).
Expand
\(u=z_1/A-P_{12}z_2/(P_{11}B),\ t=z_2/B\) in the monomials.
Each angular integral is given by the recurrences in (T5), with
\(\alpha\) replaced by \(\beta\) and the additional base case
\(I_{01}=1-\cos\beta\). The radial moments satisfy \(R_0=1\),
\(R_1=\sqrt{\pi/2}\) and \(R_j=j R_{j-2}\).
Reflection exchanges \(k_-,k_+\). The weight in (TW3) is even after
adding every reflected term. The two terms in a reflected pair are
evaluated individually on the folded wedge before adding their moments.

The [weighted checker](check_tetrahedral_weighted_gap.py) builds the
rational Hermite matrices of (T4) with a common integer denominator
per degree, then uses ball arithmetic for (TW5) and the trace
comparisons. The low-degree coefficients, the untilted wedge mass and
the strict order of the two marginal tilt masses serve as controls.
The accompanying [replay record](tetrahedral_weighted_gap_verification.json)
retains the exact slabs, degrees, upper tail bounds and verified
trace bounds.

### Certified weighted domain and replay

| Left endpoint \(a\) | Right endpoint \(c\) | Degree \(N\) | Trace power \(p\) | Strict upper bound on \(\operatorname{tr}(J^p)^{1/p}\) |
|---:|---:|---:|---:|---:|
| \(17/25\) | \(77/100\) | 25 | 4 | 0.951298 |
| \(77/100\) | \(41/50\) | 28 | 8 | 0.933372 |
| \(41/50\) | \(17/20\) | 34 | 8 | 0.962949 |
| \(17/20\) | \(7/8\) | 39 | 16 | 0.973899 |
| \(7/8\) | \(89/100\) | 44 | 16 | 0.981231 |
| \(89/100\) | \(9/10\) | 48 | 16 | 0.986667 |
| \(9/10\) | \(91/100\) | 54 | 32 | 0.994334 |
| \(91/100\) | \(183/200\) | 56 | 16 | 0.997069 |
| \(183/200\) | \(23/25\) | 60 | 32 | 0.995045 |

These bounds already reserve \(\gamma=1/1000\) in the scaled form
(TW1). Thus the full, unscaled standard form on these slabs has gap at
least \(\ell/(1000\sigma_\rho)\). The constant-coefficient comparison
supplies \(\ell/1000\) below the first slab. The endpoint concerns the
tetrahedral interface form, with the same ambient dimension lift.

Run from the repository root with Python 3.12.3 and python-flint 0.9.0:

~~~sh
.venv/bin/python research/plurality-stability/check_tetrahedral_weighted_gap.py
~~~

The full replay uses 640-bit Arb arithmetic. A separate internal
coefficient comparison through degree eight and a check of the summed
reflected masses against an independent full-wedge formula also passed.
Independent mathematical review and formal verification remain open.
