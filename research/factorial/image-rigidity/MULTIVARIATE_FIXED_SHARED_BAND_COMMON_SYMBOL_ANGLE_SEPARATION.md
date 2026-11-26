# Uniform angular separation of blended symbols on fixed shared bands

Let the two distinct ordinary exponent columns be fully shared,
omitting equal rows. Write \(d_i=b_i-a_i\), \(m_i=(a_i+b_i)/2\),
and \(W=\sum_i d_i^2/m_i\). Define
\[
 G(q)=\sum_i[(a_i+d_iq)\log(a_i+d_iq)]_{\rm c},\qquad
 \rho(q)=q(1-q)G''(q),
\]
where the subscript subtracts the chord through zero and one.
Let \(L_j\) be the endpoint-normalized log-factorial symbol
per degree at order \(j\).

**Theorem.** Fix positive rational \(W_-\le W_+\) and rational
\(0<\kappa<1\). There are a computable constant \(g>0\) and
terminal threshold such that whenever \(W\in[W_-,W_+]\)
and either \(\max\rho\le1-\kappa\) or \(\max\rho\ge1+\kappa\),
the root arguments of
\[
 P_m(z)=\sum_{k=0}^m\binom mk e^{m\widehat L_m(k/m)}z^k,
 \qquad
 \widehat L_m=(1-\lambda_m)L_n+\lambda_m L_{n+1},
 \qquad m=n,n+1
\]
have circular distance at least \(g/n^2\) for every larger
terminal and every independent choice
\(\lambda_n,\lambda_{n+1}\in[0,1]\).
The coefficients need not be mirrored.
The constants are unevaluated.

The localization argument below gives a compact root annulus
for the blended symbols. The existing fixed-band contour
charts can then be used with the same stationary values.
At positive-ray coexistence the exact packet ratios are
real on the real axis, which controls their imaginary
errors even when the blends differ between orders.

## A uniform annulus for blended factorial symbols

Suppose \(W\le W_+\). Put \(\ell_i=\min(a_i,b_i)\),
\(V=\sum_i d_i^2/\ell_i\) and
\(V_+=W_+(1+\max(2,W_+)/2)\). Then \(V\le V_+\).
Indeed if \(|d_i|\le2\ell_i\), the ratio
\(|d_i|/\ell_i\) is at most two. Otherwise \(m_i<|d_i|\)
implies \(|d_i|<W_+\), so \(|d_i|/\ell_i\le W_+\).
Multiplying each \(d_i^2/m_i\) by \(m_i/\ell_i\) proves
the stated bound. Write \(\mathcal W_+=1+V_+\) and
take \(r=(10^5\mathcal W_+)^{-1}\).

We claim that every sufficiently large-order \(P_m\)
has all root moduli in \(r<|z|<1/r\).
Let \(K=\lceil16rm\rceil+8\).
Write \(c_k\) for the coefficients of \(P_m\),
\(B=c_1\), and \(B_r=r B\).
Convexity and endpoint normalization give
\(c_k\le\binom mk\) and \(B_r\le rm\).
For \(f(x)=\log(c_x x!/B^x)\) on \(x\in[0,K]\)
use the analytic interpolation
\[
  f(x)=\log\Gamma(m+1)-\log\Gamma(m+1-x)-x\log m
                 +m(\widehat L_m(x/m)-x\widehat L_m(1/m)).
\]
The same function applies to the coefficients after
rescaling \(z\) by \(r\). In one orientation let
\(L=\max(1,\max_i |d_i|/a_i)\) and \(\alpha=2L/m\).
The positive Gamma kernel is bounded by \(1/u\),
as in the
[matched high-curvature interpolation](MULTIVARIATE_ABSOLUTE_CURVATURE_27_COMMON_SYMBOL.md#matched-logarithms-of-the-two-halves).
Since \(x_i(x/m)\ge a_i/2\) and \(m-x\ge m/2\),
for \(2\le s\le K\) it gives
\[
           |f^{(s)}(x)|\le H_s
              :=(s-2)!\alpha^{s-2}4(1+V)/m .
\]
This holds for either constituent symbol and its
convex blend. Bound the \(s\)-th forward difference
by the same \(H_s\) using the cube integral.
The exponential edge weight is then bounded by
\(2H_s\) because \(H_s\le H_2\le1/2\).

Use the
[connected-hypergraph log bound](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#sharper-connected-sum)
with \(\beta=8e(1+V)/m\). If \(m\ge1000\mathcal W_+\),
then \(L\le1+V\), \(\alpha K\le1/4\) and
\(\beta K\le1/4\). Also \(\alpha/\beta\le1/(4e)\).
The weighted tree formula and extra-edge bound
therefore give the truncated logarithm \(A\) through
degree \(K\) of \(P_m(rz)\) with
\[
   |[z^k]A|\le \frac{B_r}{k}(\Xi B_r)^{k-1},\qquad
   \Xi=\frac{3\beta}{1-\alpha K}
                  e^{\beta K/(2(1-\alpha K))}\le5\beta .
\]
The \(k=1\) term is \(B_r\) exactly.
Here \(\Xi B_r\le40e\mathcal W_+r<1/8\).
Consequently \(\|A\|_{|z|\le1}\le2B_r\) and
\(\|A\|_{|z|\le2}\le3B_r\).

The polynomial tail beyond \(K\) on \(|z|\le1\)
is bounded by \(4^{-K}\) using
\(c_k r^k\le(e rm/k)^k\) and \(k>K\ge16rm\).
Cauchy on radius two bounds the same tail of
\(e^A\) by \(e^{3B_r}2^{-K}\). Since \(P_m(rz)\)
and \(e^A\) agree through \(K\),
\[
 |P_m(rz)e^{-A(z)}-1|
     \le e^{2B_r}4^{-K}+e^{5B_r}2^{-K}
     \le e^{-7K/8}+e^{-3K/16}<1
                      \qquad (|z|\le1).
\]
We used \(B_r\le K/16\), \(\log2\ge1/2\), and
\(K\ge8\). This proves the inner zero exclusion.
Column reversal gives the reciprocal polynomial
with the same \(W,V\) and both symbols reflected,
so it gives the outer exclusion as well.

## The same stationary geometry for both blends

Use \(w=-\log z\) at the endpoint-normalized root.
Positive radial scaling changes only \(\Re w\).
The annulus puts all the fields in
\[
                  |\Re w|\le -\log r,\qquad -\pi\le\Im w\le0
\]
after conjugation. On fixed contour neighborhoods
avoiding the row cuts the Stirling expansion is
\[
                   L_j(q)=G(q)+A_0(q)/j
                                    +A_1(q)/j^2+O(j^{-3}).
\]
Uniformity follows from the same sector and
row-weight estimates as in the
[fixed-band contour cover](MULTIVARIATE_EFFECTIVE_SHARED_REGULAR_COVER.md#uniform-expansion-and-the-common-root-certificate).
The sum of inverse row midpoint powers is
bounded by \(W_+\). Subtracting the endpoint
chords preserves uniformity. Use the actual
finite-profile \(G,A_0,A_1\) in each chart.

At \(m=n,n+1\), \(m(\widehat L_m-G)=A_0+O(1/n)\).
For a relative complex order disk about \(n\), extend
\(\lambda(m)=\lambda_n+(m-n)(\lambda_{n+1}-\lambda_n)\)
and keep \(L_n,L_{n+1}\) fixed. On that disk one has
\[
 m(\widehat L_m-G)
   = b_m A_0(q)+O(1/n),\qquad
 b_m=\frac mn-\frac{m\lambda(m)}{n(n+1)}.
                                                               \tag{2}
\]
Indeed \(b_m\) is bounded even though \(\lambda(m)=O(n)\).
The common \(A_1\) terms have coefficient \(O(1/n)\),
using \(1/(n+1)^2-1/n^2=O(n^{-3})\), and the
remaining errors are \(O(1/n)\) as well.
In particular the row amplitude is analytic and
bounded away from zero on the complex-order saddle disks.

The normalized Rice integrand has the same rational
binomial kernel and power of \(z\) as the
[exact ordinary residue formula](MULTIVARIATE_RICE_CONTOUR_AND_REAL_BOUNDARY.md#exact-contour-and-its-phases),
with row factor \(e^{m\widehat L_m(q)}\).
This factor is holomorphic on the original cut-free
winding contours and throughout the region swept to
the coefficient poles. Indeed the Gamma factors at
the fixed orders \(n,n+1\) have holomorphic logarithms
with the chosen real normalization on the common
row-slit domain. The exact coefficient residues at
\(q=k/m\) give \(P_m(z)\).
Its leading phase is
\[
               G(q)-q\log q-(1-q)\log(1-q)-wq
\]
on the upper half-plane, with the same lower
continuation and real-gap reflection terms as before.
The uniform analytic bounds (2) preserve all strict
exterior losses and one-saddle exclusions.
For the positive-real boundary the corresponding
coefficient lattice-panel contours use the same
entropy phase, and the real coefficient tails
remain exponentially smaller on the same fixed disk.

Thus the paired-field witnesses of the
[fixed-band angular theorem](MULTIVARIATE_FIXED_SHARED_BAND_ANGULAR_SEPARATION.md#the-paired-field-geometry)
apply with the blended row factor. At a regular
two-saddle chart the two integer orders have the
same leading amplitude ratio, with relative
\(O(1/n)\) analytic errors. At a simple fold, (2)
and the usual binomial Gamma continuation give
the shared complex-order finite Airy integral
with nonzero analytic central amplitude. The
division by \(u^2-\lambda(w)\) on the Airy contour
uses the same fixed amplitude bounds. The integer
residue identities and strict exterior loss give
the same exponentially small integer-order
remainder as in the
[fixed-band fold argument](MULTIVARIATE_FIXED_SHARED_BAND_ANGULAR_SEPARATION.md#the-subcritical-fold-and-complex-order).
No common coefficient drift is required.

## Quantitative separation, including coexistence

Assume the root angles are within \(g/n^2\),
where \(0<g\le1\) is fixed small. By conjugating
individual roots one may take both angles in
\([0,\pi]\) without enlarging their difference.
Write \(w_j=x_j+i y_j\). All near-equal-height
field pairs at large orders are in one of the
same local alternatives as the equal-height
compact cover: a one-field exclusion or a common
two-saddle, fold, or real-coexistence chart.
Indeed the cover consists of finitely many
open witness neighborhoods of the compact
equal-height pair space, so it also covers
pairs uniformly sufficiently close to that space.

In a regular common chart the real log equations
first give \(|x_{n+1}-x_n|=O(n^{-1})+O(|y_{n+1}-y_n|)\).
Subtracting them improves this to
\(O(n^{-2})+O(|y_{n+1}-y_n|)\), by the
[real modulus comparison](MULTIVARIATE_FIXED_SHARED_BAND_ANGULAR_SEPARATION.md#regular-interior-and-lower-boundary-charts).
Here one adds the \(O(|y_{n+1}-y_n|)\) change
from moving one field vertically to the other
height inside the common chart.
The two root equations then give
\(e^{D(w_n)}=1+O(n^{-1})\), against the positive
regular phase margin.
At a fold the same vertical move in the Airy
parameter is \(O(g/n^2)\); the real displacement
bound from its imaginary derivative is \(O(1/n)\).
The Airy spacing test forces a common label as
in the fixed-band comparison. Its complex-order
branch has imaginary displacement at least a
positive multiple of \(n^{-5/3}\) by the
nonzero imaginary inverse derivative, exceeding
\(g/n^2\) for all large orders.

At the supercritical coexistence chart on the
positive ray, take actual retained left and right
lattice sums on the two saddle intervals from
the
[positive boundary construction](MULTIVARIATE_EFFECTIVE_SHARED_REAL_BOUNDARIES.md#upper-boundary).
Assign the exterior sum to either of the two
packets. The strict exterior loss makes this
addition exponentially small relative to both
packet scales on a fixed symmetric complex
disk. The two exact packets, their stationary
values, and their leading amplitudes are
analytic and real on the real axis. The saddle
expansion gives an exact root equation in
the chart
\[
             e^{mD(w_m)}=L(w_m)(1+R_m(w_m)),
                \qquad \|R_m\|=O(1/n),
\]
where \(D\) is real on the real axis, \(L\) is
real negative there, and \(R_m\) is analytic
and real there. This form uses the complete
packet error, without assuming its first
correction is common to the two orders.
Shrink the chart so that \(|D|<1/4\) and the
real derivative of \(D\) has fixed sign and
magnitude at least \(c>0\).

The real modulus comparison still gives
\(|x_{n+1}-x_n|=O(n^{-2})\). Write
\(\ell=\log(-L)\), \(r_m=\log(1+R_m)\), using
branches real on the real axis. The exponent
in the divided root equations is
\[
 D(w_n)+(n+1)(D(w_{n+1})-D(w_n))
     =\ell(w_{n+1})-\ell(w_n)
          +r_{n+1}(w_{n+1})-r_n(w_n).
\]
There is no phase-period ambiguity, since both
displayed sides have modulus less than \(1/2\)
by the small \(D\) and field displacement
bounds. For \(y=y_n\),
\(|\Im D(w_n)|\ge c|y|\) by integrating the
real derivative over the vertical segment.
The horizontal change in \(\Im D\) is bounded
by \(C|y|\,|x_{n+1}-x_n|\) because \(D'\) is
real on the axis; the added vertical change
costs \(C|y_{n+1}-y_n|\). The same bounds
hold for \(\ell\). Finally the imaginary
parts of \(r_m\) at the two fields are at
most \(C(|y|+|y_{n+1}-y_n|)/n\) by reality
on the axis and Cauchy's estimate. Thus
\[
                  c|y|\le C|y|/n+C n |y_{n+1}-y_n|.
\]
A polynomial of degree \(n\) with positive
coefficients has no root with
\(|\arg z|<\pi/n\): after rotation its
summands would all be in an open half-plane.
Consequently \(|y|\ge\pi/n\). Taking \(n\)
large and \(g\) sufficiently small contradicts
the displayed inequality. This gives the
claimed gap also at the coexistence chart.

All annulus and blend estimates are uniform in
profile dimension and use explicit majorants.
The retained-and-diffuse compact cover has the
same stationary tests as the fixed-band proof,
with the larger real cap from the annulus.
The new upper-ray tests use uniform Cauchy
bounds on its finite lattice-panel witnesses.
Thus choosing the common gap constant small
and the terminal cutoff large uses the same
finite witness search, and yields computable
constants without evaluating them.
