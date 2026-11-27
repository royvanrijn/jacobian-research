# Exponential localization in the negative voting interval

Let \(T_k\) denote the cylindrical regular Gaussian simplex and
let \(f:\mathbb R^d\to\Delta_k\) be a measurable fractional assignment,
with cell probabilities \(p_i\). For all sufficiently large \(k\),
every \(0<r\le 1/(k-1)\), and every finite \(d\ge k-1\), a competitor
satisfying \(S_{-r}(f)\le S_{-r}(T_k)\) obeys
\[
 \inf_T\sum_i\|f_i-1_{T_i}\|_{L^1(\gamma_d)}
 \le C\exp\!\left[-c/(r^2 k(\log k)^2)\right].
 \tag{VL1}
\]
The infimum ranges over orientations of the cylindrical simplex.
For all assignments on the same parameter range,
\[
 S_{-r}(f)\ge S_{-r}(T_k)
            -Cr\exp\!\left[-c/(r^2 k(\log k)^2)\right].
 \tag{VL2}
\]
Constants are absolute. The statements have no prescribed mass
constraint. The improvement uses the simplex score's reflection
identity in the linear higher-degree term and a relative
quadratic Hermite bound.

## Notation and inputs

Write \(L=\log k\), \(z=\Phi^{-1}(1-1/k)\asymp\sqrt L\), and use
\[
 \Delta=F_k-W_1(f)+\delta^2/r,\qquad
 \delta=\|p-\mathbf1/k\|_2.
\]
For a possible improver, the [covariance and centroid-spread estimates](NEGATIVE_VOTING_INTERVAL.md#dimension-uniform-localization-on-the-voting-scale)
give
\[
 \Delta\ge 0,\qquad
 \|Q-Q_*\|_F^2\le Ck\Delta,\qquad
 H\le CkL\Delta,\qquad p_i\le C/k .
 \tag{VL3}
\]
Here \(H\) is the conditional centroid-length variance from (V11),
\(Q\) is the normalized centroid Gram matrix, and \(Q_*\) is its
regular counterpart. These are (V3a), (V18) and the spread bound
following (V18) in that source. Choose the same polar alignment
of the simplex \(T\) and write \(b_i^*\) for its moment vectors.
Set \(e_i=f_i-1_{T_i}\) and let
\[
 \mathcal R = \int\left(\max_i b_i^*\cdot x
                         -\sum_i f_i(x)b_i^*\cdot x\right)
                          d\gamma_d(x).
\]
The polar comparison separates the score scale from covariance
shape. On the active label space \(Q_*^{1/2}\) is scalar, so
\(\|Q^{1/2}-Q_*^{1/2}\|_F^2\le(k-1)\|Q-Q_*\|_F^2\).
Also \((\sqrt{W_1(f)}-\sqrt{F_k})^2\le|W_1(f)-F_k|\le\Delta\).
With \(W_1(f)\asymp L/k\), the triangle inequality for the moment
matrices yields
\(\|B_f-B_*\|_F^2\le C L\|Q-Q_*\|_F^2+C\Delta\le CkL\Delta\).
Since the exact regret is
\((F_k-W_1(f)+\|B_f-B_*\|_F^2)/2\), it follows that
\[
 \mathcal R \le C k L\Delta .
 \tag{VL4}
\]
All the estimates below concern that aligned simplex. Bounds on
active variables therefore have the same meaning for all ambient
dimensions under consideration.

We use a basic normal density estimate. For fixed integers
\(a\ge0\), \(j\ge1\), \(l\ge0\) and fixed \(D\ge0\),
\[
 \int_0^\infty (1+y)^a e^{D y/z}\phi(y)^j
                     \Phi(y)^{k-l}\,dy
       \le C_{a,j,l,D}\frac{z^{a+j-1}}{k^j}.
 \tag{VL5}
\]
In the change of variable \(w=k(1-\Phi(y))\), one factor \(\phi\)
cancels the Jacobian. Mills' bound gives
\(\phi(y)\le C(1+y)w/k\). The Gaussian bound
\(y\le\sqrt{2\log(k/w)}\) yields
\((1+y)/z\le C(1+\sqrt{\log_+(1/w)})\) and
\(e^{D y/z}\le C_D\max(1,w^{-1/2})\), the latter by Young's
inequality in \(\sqrt{\log_+(1/w)}\).
For \(k\ge 2l\) the remaining power of \(\Phi\) is at most
\(e^{-w/2}\). The resulting majorant after division by
\(z^{a+j-1}/k^j\) is integrable at zero and infinity for \(j\ge1\),
which proves (VL5).

## The higher-degree terms

Represent the regular moment scores on the active label space
\(E=\mathbf1^\perp\) as \(a_* x_i\), where
\(a_*=\sqrt{F_k/(k-1)}\asymp z/k\).
The regular score strip bound in
[(V17)](NEGATIVE_VOTING_INTERVAL.md#dimension-uniform-localization-on-the-voting-scale)
and the definition of \(\mathcal R\) give
\[
 \left(\sum_i\|e_i\|_1\right)^2\le C k\mathcal R .
 \tag{VL6}
\]
Indeed a transferred label outside the top score strip of
width \(t\) costs at least \(t\). The strip bound and the
total regret bound the transferred mass by
\(C kt+\mathcal R/t\).

Let \(P_j\) be Gaussian Hermite projection and
\(R_r=\sum_{j\ge3}(-r)^jP_j\) the signed higher-degree
operator. Its contribution to the difference of stabilities is
\[
 2\sum_i\langle e_i,R_r1_{T_i}\rangle
       +\sum_i\langle e_i,R_re_i\rangle .
 \tag{VL7}
\]
The
[relative quadratic estimate (RT2)](NEGATIVE_VOTING_RELATIVE_TAIL.md)
applies since \(p_i\le C/k\) in (VL3).
Using (VL4) gives the bound
\[
 \left|\sum_i\langle e_i,R_re_i\rangle\right|
 \le C r^3 L^3\mathcal R\le C r^3 k L^4\Delta .
 \tag{VL8}
\]
This bounds the signed quadratic term at every transfer scale.

For the linear term use the actual simplex response
\[
 (T_\rho1_{T_i})(x)=\Pi_i(\kappa x),
 \quad \kappa=\frac{\rho}{\sqrt{1-\rho^2}},
 \qquad
 \Pi_i(w)=\int\phi(y-w_i)
                  \prod_{l\ne i}\Phi(y-w_l)\,dy .
\]
This is the conditional winning probability of independent normals
with mean shifts \(w_i\). Projecting those normals to \(E\) does
not change their winner. For \(1\le m\le3\),
\[
 |D^m\Pi_i(w)[v_1,\ldots,v_m]|
 \le C_m\frac{z^m}{k}\prod_{\ell=1}^m\|v_\ell\|_\infty
 \qquad(\|w\|_\infty\le 1/z).
 \tag{VL9}
\]
To prove this, distribute the \(m\) derivatives among the factors
under the integral. If the leading density receives \(a\) derivatives
and \(j\) distinct CDF factors receive \(b_1,\ldots,b_j\)
derivatives, then \(a+\sum b_l=m\).
The density derivatives produce polynomial multiples of the
density and each differentiated CDF factor produces a density
times a polynomial of degree \(b_l-1\).
Hence an integral of this type contains \(1+j\) density factors,
a polynomial bounded by \(C_m(1+|y|)^{m-j}\), and \(k-1-j\)
undifferentiated CDF factors. There are at most \(C_m k^j\)
such terms after bounding each directional coordinate by its
supremum norm. For \(y\ge0\), put \(t=y+1/z\).
Each CDF is bounded by \(\Phi(t)\), and each density is bounded
by \(\phi(t)e^{2t/z}\). This follows by translating by at most
\(2/z\) from \(t\). Applying (VL5) with \(a=m-j,\ j\) there
replaced by \(1+j\), bounds each such integral by
\(C_m z^m/k^{1+j}\). For \(y<0\), the CDF factors are at most
\(\Phi(1/z)\le3/4\) for large \(k\). A remaining density factor
gives integrability of the polynomial and an exponentially small
contribution even after the factor \(C_mk^j\).
Dominated differentiation on compact shift sets proves (VL9).

Fix \(M\ge1\) with \(2rM\le1/z\), and let \(\mathcal C_M\) be
the cap \(\max_i|x_i|\le M\) in the active coordinate representation.
For a unit normal \(n=(e_i-e_j)/\sqrt2\) to the pair plane,
\(\partial_n T_\rho1_{T_l}=\kappa D\Pi_l(\kappa x)[n]\).
Taylor expansion of \(D\Pi_l\) to first order and (VL9) give
\[
 \left|\partial_n R_r1_{T_l}(x)\right|
                  \le C r^3 z^3 M^2/k
                  \qquad(x\in\mathcal C_M).
\]
The first two Taylor coefficients in \(\rho\) are the derivatives
of \(P_1 1_{T_l}\) and \(P_2 1_{T_l}\), by Mehler diagonalization.
This may be checked by differentiating the conditional Gaussian
kernel at zero. On compact \(x\)-sets its first \(x\)-derivative
and each of the requisite correlation derivatives are dominated
by a Gaussian times a polynomial in the integration variables.
Replacing \(\kappa\) by \(\rho\) in the first two coefficients
costs only \(O(r^3 z/k+r^4 z^2M/k)\), covered by the displayed
error.

By reflection interchanging labels \(i,j\), \(R_r(1_{T_i}-1_{T_j})\)
is odd across their pair plane. The normal projection of \(x\) to
the plane averages coordinates \(i,j\), preserving \(\mathcal C_M\)
along the projection segment. Thus for \(x\in T_i\cap\mathcal C_M\)
the absolute difference of the two higher-degree scores is at most
\(C r^3 z^2 M^2\,a_*(x_i-x_j)\).
Summing against the ordered transfers, the cap contribution to
the linear term in (VL7) is at most \(C r^3 L M^2\mathcal R\).

Outside the cap the Gaussian probability is at most
\(2k e^{-M^2/2}\), since each \(x_i\) has variance at most one.
At comparison scale \(s_0=1/L\), Hermite diagonalization gives
\(\|R_r1_{T_i}\|_2\le (r/\sqrt{s_0})^3\|T_{\sqrt{s_0}}1_{T_i}\|_2
\le C r^3 L^{3/2}/k\) by the small-correlation smoothing estimate.
Cauchy–Schwarz and \(|e_i|\le1\), summed over labels, therefore give
\[
 \left|\sum_i\langle e_i,R_r1_{T_i}\rangle\right|
 \le Cr^3 L M^2\mathcal R
       + Cr^3 L^{3/2}\sqrt{k}\,e^{-M^2/4}.
 \tag{VL10}
\]

## Absorption and the resulting bounds

The second-degree comparison (V13)–(V15) in
[the voting-interval argument](NEGATIVE_VOTING_INTERVAL.md#centroid-length-spread-and-the-second-hermite-deficit)
can use the truncation height \(M\) instead of \(R\).
The ordered coupling after the common Gaussian shift gives
\(B_*-B\le C(M/\sqrt L)(k\Delta+2H)+Ck e^{-M^2/4}\) by (V14),
using uniform \(\alpha=(k-1)/k\ge1/2\).
Dropping the nonpositive matrix penalty from the
second-degree deficit and using (VL3) gives
\[
 W_2(T_k)-W_2(f)
 \le C M L^{3/2}\Delta + C L e^{-M^2/4}.
 \tag{VL11}
\]
The three terms in stability are the zeroth- and first-degree
difference \(r\Delta\), the second-degree difference
\(r^2(W_2(f)-W_2(T_k))\), and (VL7).
Combining (VL4), (VL8), (VL10) and (VL11), for a possible improver,
\[
 \Delta \le
 C(rM L^{3/2}+r^2 k L^2 M^2+r^2 k L^4)\Delta
 +C(rL+r^2 L^{3/2}\sqrt{k})e^{-M^2/4}.
 \tag{VL12}
\]

Choose \(M^2=b/(r^2 k L^2)\) for a fixed small absolute \(b>0\).
For sufficiently large \(k\), uniformly in the voting interval,
\(M\ge1\) and \(2rM\le1/z\).
The three coefficients before the common \(C\) in (VL12) are
\(\sqrt{bL/k}\), \(b\), and \(O(L^4/k)\).
Choose \(b\) small and then \(k\) large to absorb their sum.
This gives \(\Delta\le C e^{-M^2/4}\), since the polynomial
prefactor in the last term of (VL12) is bounded for large \(k\).
By (VL6) and (VL4) the label distance is at most
\(C k\sqrt L\,e^{-M^2/8}\). As
\(M^2\gtrsim k/L^2\gg\log k\), this proves (VL1).
For (VL2), use the same estimates in the actual stability
expansion. After multiplication by \(r\), every error in
(VL12) is at most \(C r e^{-M^2/4}\), using the bound on
\(\Delta\). A positive actual deficit already satisfies (VL2).

The [polynomial local radius](NEGATIVE_VOTING_POLYNOMIAL_NEIGHBORHOOD.md)
contains every such localized assignment for sufficiently large \(k\).
Together they prove the
[exact same-partition inequality](NEGATIVE_VOTING_INTERVAL.md#exact-gaussian-simplex-minimization).
