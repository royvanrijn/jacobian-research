# Circle roots along mirrored physical Gaussian mergers

Let \((X,Y)\) be a correlated gamma pair with common shape
\(\gamma>0\), unit marginal scale, and transform
\[
 \mathbb E e^{uX+vY}
   =\bigl(1-u-v+(1-\rho)uv\bigr)^{-\gamma},\qquad 0\le\rho<1 .
\]
At \(\gamma=1\) these are the squared moduli of the two complex
Gaussian coordinates in the physical row merger, with squared
correlation \(\rho\). For each integer \(r\ge0\), gap \(d\ge1\),
and degree \(m\ge1\), define
\[
          P_{m,r,d}(z)=\mathbb E[(XY)^r(X^d+zY^d)^m].
\]
Every \(P_{m,r,d}\) has simple unit-circle roots. In particular this
holds with \(r=am\) for ordinary mirrored rows
\((a,a+d),(a+d,a)\), uniformly in each physical correlation
strictly below one. The same conclusion holds after combining any
number of independent such pairs, with independent choices of the
correlations, shapes, smaller exponents, and gaps. Additional
independent equal rows contribute only a nonzero scalar.

For unit gap, every noncentral root moves strictly toward the
negative real axis as \(\rho\) increases. With zero smaller
exponent and unit gap, consecutive orders strictly interlace
throughout \(0\le\rho<1\) for a single correlated pair.
The root-location theorem applies at every degree without
requiring a common order comparison across the blocks.

## The exact angular polynomial for unit gap

Write \(\delta=1-\rho\), \(A=1-\delta u\), \(B=1-\delta v\),
and \(D=AB-\rho\). The transform is
\(\delta^\gamma D^{-\gamma}\). Its positivity as a real measure
also follows by taking \(J\) with probabilities
\(\delta^\gamma(\gamma)_j\rho^j/j!\), then using conditionally
independent gamma variables of shape \(\gamma+J\) and scale
\(\delta\).

Differentiate \(r\) times in each of \(u,v\) to tilt the law by
\((XY)^r\). The resulting transform is
\[
 (\gamma)_r^2\delta^{\gamma+2r} D^{-\gamma-2r}
               \sum_{\ell=0}^r\frac{(-r)_\ell^2}
                   {(\gamma)_\ell\ell!}(AB)^{r-\ell}\rho^\ell.
                                                               \tag{1}
\]
To see this, expand \((AB-\rho)^{-\gamma}\) in \(\rho/(AB)\)
near \(u=v=0\). The differentiated series is
\((\gamma)_r^2(AB)^{-\gamma-r}
 {}_2F_1(\gamma+r,\gamma+r;\gamma;\rho/(AB))\).
Euler's transformation extracts \((1-\rho/(AB))^{-\gamma-2r}\)
and leaves the terminating series in (1). The transformation
follows by substituting \(y(w)=(1-w)^{\gamma-a-b}v(w)\)
into \(w(1-w)y''+[\gamma-(a+b+1)w]y'-ab y=0\).
The resulting parameters for \(v\) are \(\gamma-a,\gamma-b,\gamma\).
The power series solution normalized at zero is unique by the
coefficient recurrence, since \(\gamma>0\).

For \(0\le j\le r\), the coefficient of \(D^{r-j}\) in the
numerator of (1) is
\[
    \binom rj \frac{(\gamma+r)_j}{(\gamma)_j}\rho^j.
\]
Indeed expand \(AB=D+\rho\), and after taking out
\(\binom rj\rho^j\) the remaining sum is
\(\sum_{\ell=0}^j(-r)_\ell(-j)_\ell/((\gamma)_\ell\ell!)\).
Its value is \((\gamma+r)_j/(\gamma)_j\), as follows by taking
the coefficient of \(t^j\) in
\((1-t)^{-\gamma}(1+t/(1-t))^r=(1-t)^{-\gamma-r}\)
and dividing by \((\gamma)_j/j!\).

Put \(\alpha=\gamma+r\) and \(q=\rho/\delta\). At \(u=t,v=zt\)
we have \(D=\delta(1-(1+z)t+\delta zt^2)\). For \(z=e^{2i\theta}\)
put \(x=\cos\theta/\sqrt{\delta}\).
Using the generating function of \(C_m^\alpha\), the Gegenbauer
polynomial, the tilted moment of order \(m\) becomes
\[
 P_{m,r,1}(e^{2i\theta})
   =m!(\gamma)_r^2\delta^r(\sqrt{\delta}\,e^{i\theta})^m H_m(x),
 \qquad
 H_m(x)=\sum_{j=0}^r
       \binom rj \frac{(\alpha)_j}{(\gamma)_j}q^j
                                      C_m^{\alpha+j}(x).       \tag{2}
\]
All parameter dependence due to the correlation is explicit here.

## Meixner factors and contraction of real roots

For \(E=x\,d/dx\) the Gegenbauer generating function gives
\[
      (\alpha+m/2+E/2) C_m^\alpha=\alpha C_m^{\alpha+1}.
\]
In fact apply \(\alpha+(t\partial_t+x\partial_x)/2\) to
\((1-2xt+t^2)^{-\alpha}\) and read the coefficient of \(t^m\).
Consequently \(H_m=T(E/2)C_m^\alpha\), where
\[
 T(v)=\sum_{j=0}^r\binom rj
             \frac{(\alpha+m/2+v)_j}{(\gamma)_j}q^j
     = M_r(-\alpha-m/2-v;\gamma,\delta).
                                                               \tag{3}
\]
Here \(M_r(s;\gamma,c)={}_2F_1(-r,-s;\gamma;1-1/c)\).
For \(0<c<1\), it is orthogonal on the nonnegative integers
with weight \((\gamma)_k c^k/k!\). One direct verification uses
\[
 \sum_{r\ge0}\frac{(\gamma)_r}{r!}M_r(k;\gamma,c)t^r
        =(1-t/c)^k(1-t)^{-k-\gamma}.
\]
Expand the terminating series and interchange the two sums near
zero to check the identity. Summing the product of two such
generating functions against the weight multiplied by
\((1-c)^\gamma\) gives \((1-ts/c)^{-\gamma}\).
Coefficient comparison gives orthogonality. The usual
sign-change argument now places all \(r\) zeros in \((0,\infty)\),
simple: if there are fewer than \(r\) sign-changing zeros
in that interval, multiplying by the linear factors at those
zeros gives a nonzero sum of one sign against a polynomial of
degree less than \(r\). Infinite support makes the sum strict.

For \(0<\rho<1\), the polynomial \(T\) therefore has factors
\(v+\alpha+m/2+s_i\), with all \(s_i>0\). At \(\rho=0\),
the multiplier is one. At \(r=0\) it is one for any \(\rho\).
Each corresponding Euler operator \(E+B\), \(B>0\), preserves
simple real symmetric roots in \((-1,1)\). To prove this,
write a polynomial of parity \(\epsilon\in\{0,1\}\) as
\(x^\epsilon f(x^2)\). If the positive squared roots of \(f\)
are \(0<u_1<\cdots<u_p<1\), the transformed squared roots
are the zeros between the poles of
\[
                 B+\epsilon+2\sum_{i=1}^p\frac{s}{s-u_i}.
                                                               \tag{4}
\]
The expression decreases strictly on each interval without
a pole because its derivative is
\(-2\sum_i u_i/(s-u_i)^2\).
It is positive at zero and tends to minus infinity at the
first pole from below. The same change of sign occurs between
each successive pair of poles. Thus it has exactly one simple
zero in every interval \((0,u_1),(u_1,u_2),\ldots\).
The central simple root in odd degree remains simple.
The base Gegenbauer polynomial has simple symmetric roots in
\((-1,1)\) by its positive orthogonality weight. Applying the
factors in (3) proves the same statement for \(H_m\).
Equation (2) then gives all \(m\) simple unit-circle roots.

The same computation gives monotonicity for unit gap. The
monic Meixner polynomial \(p_r(s;c)\) has
\[
 \partial_c p_r=-\frac{h_r}{c\,h_{r-1}}p_{r-1},\qquad
 h_j=\sum_{k\ge0}\frac{(\gamma)_k c^k}{k!}p_j(k;c)^2 .
\]
Differentiate its orthogonality to \(p_j\) for \(j<r\).
Only \(j=r-1\) contributes, since the derivative of the
weight gives \(k/c\), and \(kp_{r-1}=p_r+\) lower terms.
The opposite-degree monic orthogonal polynomials interlace
strictly by their three-term recurrence. Implicit differentiation
at a root of \(p_r\) now shows that each \(s_i\) increases
with \(c\). Thus every Euler parameter in (3) decreases as
\(\rho\) increases. In (4) each zero is increasing in \(B\)
and in each old \(u_i\), because the derivative in \(s\) is
negative and the corresponding parameter derivatives are
positive. Applying the commuting factors in increasing index
order, the squared positive roots of \(H_m\) decrease as
\(\rho\) increases. This comparison extends to zero correlation
by the coefficient formulas and continuity of the simple
limits. The positive factor \(\sqrt{1-\rho}\) in the map from
\(x\) to \(\cos\theta\) makes the physical angular motion
strict away from the central root.
For zero tilt, (2) reduces to one fixed-parameter Gegenbauer
polynomial at each degree. Strict consecutive interlacing in
that case follows directly from orthogonality and the same
angle map in both degrees.

## Passing to arbitrary gaps and multiple independent pairs

For \(M=dm\), take a real reciprocal degree-\(M\) polynomial
\(F(w)=\sum_{j=0}^M\binom Mj f_j w^j\) with unit-circle roots.
The linear map
\[
                    F\ \longmapsto\
        \sum_{k=0}^m\binom mk f_{dk}z^k                         \tag{5}
\]
preserves unit-circle location. Indeed, for fixed \(0<|z|<1\),
form \(B_z(w)=(1+z w^d)^m\). The right side of (5) is the
degree-\(M\) Schur–Szegő composition of \(F\) with \(B_z\),
evaluated at \(w=1\). Every root \(\zeta\) of \(B_z\) has
\(|\zeta|=|z|^{-1/d}>1\). By the
[polarization argument for circular composition](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#mirrored-profiles-at-every-fixed-order),
this evaluation is a nonzero scalar times the polarization
of \(F\) evaluated at the points \(-1/\zeta\), all in the
open unit disk. The Grace–Walsh–Szegő theorem gives a nonzero
value there. At zero the output is \(f_0\ne0\). Reciprocity
of the real output places all its roots on the circle.

The map (5) is surjective on the real reciprocal coefficient
spaces, by assigning the coefficients at multiples of \(d\)
and putting the others to zero. It is therefore open between
those spaces. Simple circle-rooted polynomials are precisely
the relative interior of the circle-rooted set in each
fixed-degree reciprocal space with nonzero leading coefficient.
This follows by removing the structural \(1+w\) factor in odd
degree and using the real polynomial in \(w+w^{-1}\).
The interior condition is that its roots are distinct and in
\((-2,2)\). Thus (5) preserves simple circle roots.
Apply it to \(P_{dm,r,1}\), retaining the same tilt \(r\).
Its output is exactly \(P_{m,r,d}\).

For independent blocks the combined coefficients, divided
by \(\binom mk\), are the products of the corresponding
block coefficients. Circular composition preserves this
location and simplicity, as in the
[ordinary mirrored factorization](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md).
At correlation one a block has \(X=Y\) in law and contributes
only a positive scalar times \((1+z)^m\), which is the
identity under this composition after scalar normalization.
Thus simple circle roots persist at physical merger endpoints
provided at least one other unequal block remains at a
correlation below one.

The circle-location theorem also preserves any radial
separation already available for an independent residual binary
profile. At each order, composition with the correlated mirror
blocks can only keep the residual roots inside their original
closed annulus. Apply the
[disc comparison and its reversed form](MULTIVARIATE_PAIRED_GAMMA_ROOT_ANNULUS.md#matched-factor-bounds-for-the-ordinary-profiles)
since each additional factor has both radii equal to one.
Thus disjoint physical radius intervals for the residual profile
remain disjoint throughout all these correlations, including
merger endpoints.

For general positive tilt, gap, and multiple pairs, strict
consecutive interlacing requires a further comparison across
the degrees. The Meixner factorization provides exact
degree-dependent Euler parameters for that comparison.

## Strict coefficient flattening for every balanced-gap pair

The physical correlation kernel gives a separate comparison
without requiring mirror symmetry of the smaller exponents.
Suppose the two rows have opposite integer gaps of the same
nonzero magnitude. At a fixed moment order their two exponent
arguments \(p,q\) then have fixed total \(T=p+q\) across the
coefficient array. For this block the logarithmic second
differences of the binomial-normalized coefficients decrease
strictly with \(\rho\) on \([0,1]\) and vanish at one.
Every interior coefficient after endpoint normalization strictly
increases to one. These conclusions persist if an independent
unchanged block is included, with the unchanged block's curvature
and endpoint-normalized coefficients as the respective limits.

In fact,
\[
 \mathbb E X^p Y^q=(\gamma)_p(\gamma)_q F_S(\rho),\qquad
 F_S(\rho)={}_2F_1(-p,-q;\gamma;\rho),\qquad S=pq .
\]
This follows by taking \(p,q\) derivatives of the transform
and using the same power-series expansion and Euler transformation
as in (1). For fixed \(T\), define
\[
        Y_S(\rho)=\rho\, F_S'(\rho)/F_S(\rho).
\]
At each \(0<\rho<1\), \(Y_S\) is strictly increasing and
strictly concave on the ordered set of distinct products \(S\)
available at that \(T\). The hypergeometric equation gives
\[
 (1-\rho)Y_S'=S-\left(\frac{\gamma-1}{\rho}+T\right)Y_S
                        -\frac{1-\rho}{\rho}Y_S^2 .          \tag{6}
\]
The expansion at zero is
\[
 Y_S(\rho)=\frac{S}{\gamma}\rho
  -\left(\frac{S^2}{\gamma^2(\gamma+1)}
          +\frac{(T-1)S}{\gamma(\gamma+1)}\right)\rho^2
   +O(\rho^3).
\]
This follows from the first two nonconstant coefficients of
\(F_S\), which are \(S/\gamma\) and
\(S(S-T+1)/(2\gamma(\gamma+1))\).
Strict increase and strict concavity hold at small positive
\(\rho\). At a first equality for an increasing pair, subtract
(6) for the two products. The derivative on the right has the
strictly positive forcing from their difference in \(S\),
and the other terms cancel. At a first equality in the
concavity comparison, subtract the weighted average of (6)
at the two endpoint products from (6) at their weighted
product. The linear terms cancel and the quadratic term gives
strictly positive weighted variance, since the endpoints
have distinct \(Y_S\). Both signs preclude a first equality.

For the two rows at hand \(p=p_0+kd\) and \(q=T-p\)
up to reversal, so \(S\) is a strictly concave quadratic
function of \(k\). Interpolate \(Y_S\) linearly between the
available products. Its interpolant is strictly increasing
and concave, hence \(Y_S/\rho\) is strictly concave along the
integer \(k\) array. Because this is
the \(\rho\) derivative of the log mixed moment, each log
second difference decreases strictly along the path. At
\(\rho=1\), \(X=Y\) and the mixed moment only depends on \(T\),
so the log second differences of this block are zero. Subtracting
the affine endpoint interpolation from \(Y_S/\rho\) gives
a strictly positive interior derivative for the log normalized
coefficient. At the merger every normalized coefficient is one.
The derivative comparison uses ordinary lattice products at
fixed total and holds for every order and every integer gap.

## Matching physical corrections across homogeneous orders

Fix positive \(u,v\) and let the order \(j\) run over positive
integers for which \(ju,jv\) are integers. Take
\(\tau=u+v,\ s=uv\), and put
\[
       Z_j(\rho)=\frac{\rho}{j}
             \frac{d}{d\rho}\log{}_2F_1(-ju,-jv;\gamma;\rho).
\]
For \(\gamma\ge1\) and \(0<\rho<1\), \(Z_j(\rho)\) strictly
increases with \(j\). It converges to
\[
 B(\rho)=\frac{2s}{\tau+
                      \sqrt{\tau^2+4s(1-\rho)/\rho}}.         \tag{7}
\]
This is a comparison at the same fractional arguments in
both orders.

To prove it, divide (6) by \(j\). The equation is
\[
 (1-\rho)Z_j'
    =j\left(s-\tau Z_j-\frac{1-\rho}{\rho}Z_j^2\right)
                               -\frac{\gamma-1}{\rho}Z_j .
\]
The positive zero of the parenthesis is \(B(\rho)\), which
strictly increases. At small positive \(\rho\), \(Z_j=O_j(\rho)\)
and \(B\) is of order \(\sqrt{\rho}\). At a first contact
\(Z_j=B\), the derivative of \(Z_j\) is nonpositive whereas
\(B'\) is positive, so \(Z_j<B\) throughout the interval.
At a first contact in the order comparison the parenthesis
is therefore strictly positive, and its multiplier \(j\)
gives the higher order a strictly larger derivative.
The comparison starts at zero since
\(Z_j=js\rho/\gamma+O_j(\rho^2)\).

For the limiting statement subtract the equation from the
derivative of \(B\). The error \(D_j=B-Z_j\) satisfies
\[
 D_j'=B'+\frac{\gamma-1}{\rho(1-\rho)}Z_j
        -j\left(\frac{\tau}{1-\rho}
                +\frac{B+Z_j}{\rho}\right)D_j.
\]
On every compact subinterval of \((0,1)\), use a larger
compact interval starting to its left. The forcing term is
bounded there, since \(Z_j<B\). The damping is at least
\(j\tau\), so the integrating factor gives \(D_j=O(j^{-1})\)
on the smaller interval. By dominated convergence, using
\(0\le Z_j/\rho\le B/\rho\) and the integrable square-root
bound at zero, the correction per degree has the exact
limit
\[
 \lim_{j\to\infty}\frac1j
            \log{}_2F_1(-ju,-jv;\gamma;\rho)
                    =\int_0^\rho B(t)\,\frac{dt}{t}
                    =-u\log(1-B(\rho)/u)-v\log(1-B(\rho)/v).
                                                               \tag{8}
\]
For the evaluation use
\(t=B(t)^2/[(u-B(t))(v-B(t))]\), which follows from
the positive root equation. Thus \(B(t)\,d\log t\) is
\([u/(u-B)+v/(v-B)]\,dB\), and \(B(0)=0\).
For opposite equal gaps and positive arguments \(u,v\), their
sum \(\tau\) is fixed over the fractional coefficient grid.
The right side of (7) is an increasing, strictly concave
function of their product \(s\) at \(0<\rho<1\), because
it is the positive inverse of
\(z\mapsto\tau z+(1-\rho)z^2/\rho\).
As \(s\) is a strictly concave quadratic in the fractional
index, the rate of change of the limiting per-degree log
curvature under correlation is strictly negative as well.

On compact fractional intervals with both arguments positive,
the derivatives of the integrand in (8) through order two are
bounded by a constant times \(t^{-1/2}\) near zero, as is seen
from the square root in (7). This justifies differentiation of
the limiting curvature. For a fully shared balanced block,
the maximum of its stationary curvature over the fractional
interval decreases strictly with correlation, because the
pointwise decrease is strict in the interior and the
stationary curvature vanishes at each endpoint.
