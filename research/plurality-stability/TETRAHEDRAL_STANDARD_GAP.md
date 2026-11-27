# The balanced tetrahedral interface gap through correlation 17/25

For every \(0<\rho\le17/25\), the tetrahedral standard sector of the
Gaussian simplex interface form has gap at least \(\ell/1000\). Here
\(\ell\) is the distance between two of the four Gaussian cell-moment
vectors. The other tetrahedral sectors are covered by
[(S56), (S60) and (S61)](SIMPLEX_INTERFACE_SYMMETRY.md#four-cells-domination-settles-the-remaining-sign-block).
It follows that for every fixed correlation in this range and every
ambient dimension \(d\ge3\), the balanced full-interface form has a
strictly positive gap off the geometric rotations, with a gap constant
independent of \(d\) by [(S42)](SIMPLEX_INTERFACE_SYMMETRY.md#every-cell-count-the-active-gap-lifts-to-all-ambient-dimensions).
The remaining standard-sector correlations are \(17/25<\rho<1\).
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
