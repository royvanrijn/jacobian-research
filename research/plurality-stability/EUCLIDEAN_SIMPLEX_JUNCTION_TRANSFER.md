# Euclidean comparison for Gaussian simplex junctions

Let \(T\) be the regular simplex partition of its active space
\(E=\mathbf 1^\perp\), with dimension \(d=k-1\), \(k\ge3\).
For a symmetric nonnegative kernel \(K\) and ordinary partition
\(U=T\) outside \(B_R\subset E\), define \(h=1_U-1_T\) and
\[
 L_K(U)=\sum_{i\ne j}\int_{T_i\cap U_j}
     \int K(p,q)(1_{T_i}-1_{T_j})(q)\,dq\,dp,
 \qquad
 Q_K(U)=\iint K(p,q)h(p)\cdot h(q)\,dp\,dq .
\]
The finite relative cut is \(\Delta_K=2L_K-Q_K\). Write
\(\Delta_b,L_b,Q_b\) for the convolution kernel
\(\phi_b(p-q)\), where \(\phi_b\) is the centered Gaussian density
with covariance \(b^2 I_d\). Write \(\Delta_e^{\mathrm G},L_e^{\mathrm G}\)
for the kernel (BC7) in the
[Gaussian scaling formula](TWO_BLOCK_SIMPLEX_JUNCTION.md).
The physical stability deficit and score cost are \(e^d\) times
the corresponding scaled quantities.

## Comparing a weighted Gaussian cut to ordinary Euclidean cuts

For \(R\ge1\), \(e^2R\le1\), and \(\rho\) sufficiently close to one, set
\[
 a^2=(1+\rho)/2,\qquad b=a/\sqrt{\rho},\qquad
 w(p)=\exp[-e^2|p|^2/(4a^2)],\qquad C_e=(2\pi a^2)^{-d/2}b^d .
\]
For \(0<u<1\), let \(U_u=U\) on \(\{w^2>u\}\) and \(U_u=T\)
elsewhere. The partitions \(U_u\) are ordinary and their possible
transfers are all inherited from \(U\). Uniformly over \(U\),
\[
 \begin{split}
 \left|\frac{L_e^{\mathrm G}(U)}{C_e}
                   -\int_0^1 L_b(U_u)\,du\right|
   &\le C_k e^2(R+1)\int_0^1 L_b(U_u)\,du,\\
 \left|\frac{\Delta_e^{\mathrm G}(U)}{C_e}
                   -\int_0^1 \Delta_b(U_u)\,du\right|
   &\le C_k e^2(R+1)\int_0^1 L_b(U_u)\,du .
 \end{split}                                                     \tag{ET1}
\]
No balance or transfer-degree restriction enters (ET1).

Indeed, the midpoint identity
\(|(p+q)/2|^2=(|p|^2+|q|^2)/2-|p-q|^2/4\) gives the exact
factorization
\[
                 K_e(p,q)=C_e w(p)w(q)\phi_b(p-q).
                                                               \tag{ET2}
\]
For \(p\in T_i\), use \(s=p\cdot n_{ij}\ge0\), where \(n_{ij}\)
is the unit oriented normal to the central \(ij\)-plane.
Write \(t=q\cdot n_{ij}\). Reflection of \(q\in T_i\) through
that plane sends it to \(T_j\) and preserves its radial weight.
In the reflected difference the integrand is
\(\phi_b(p-q)(1-\exp[-2st/b^2])\) for \(q\in T_i\).
In particular the Euclidean pointwise score \(\ell_{b,ij}\)
obeys
\[
                c_k\min(s,1)\le\ell_{b,ij}(p)
                               \le C_k\min(s,1).                 \tag{ET3}
\]
For the lower bound restrict \(q-p\) to a fixed ball strictly
inside \(T_i\) on which \(n_{ij}\cdot(q-p)\ge1\). The sum of
two points of the cone belongs to the cone. For the upper
bound, \(1-\exp[-2st/b^2]\le C\min(s,1)(1+|q-p|)\), since
\(t\le s+|q-p|\) when \(s\le1\).

Replacing \(w(q)\) by \(w(p)\) in the reflected row incurs an
error bounded by \(C_k e^2(R+1)w(p)\min(s,1)\).
For \(|p|\le R\) and \(z=q-p\), apply \(|\exp(u)-1|\le
|u|\exp(u_+)\) to
\(u=-e^2(2p\cdot z+|z|^2)/(4a^2)\). This gives
\[
 \left|\frac{w(q)}{w(p)}-1\right|
   \le C e^2(R+|z|)|z|\exp(Ce^2 R |z|).
\]
Since \(e^2R\le1\), the reflected Gaussian integrand absorbs
these displacement factors. Multiplying by \(w(p)\) from
(ET2) and using (ET3) proves the first estimate of (ET1),
since \(\int_0^1 1_{\{w(p)^2>u\}}du=w(p)^2\).

For the quadratic term, the two threshold indicators integrate
to \(\min(w(p)^2,w(q)^2)\). On \(B_R^2\) its difference from
\(w(p)w(q)\) is at most
\[
 w(p)w(q)|\log w(p)-\log w(q)|
                 \le C e^2R |p-q|w(p)w(q).
\]
Write \(D=\{U\ne T\}\) and \(\bar L=\int_0^1 L_b(U_u)\,du\).
The required weighted interaction estimate is
\[
 \iint w(p)w(q)1_D(p)1_D(q)|p-q|\phi_b(p-q)\,dp\,dq
                   \le C_k\bar L .                              \tag{ET4}
\]
The estimate also holds without the displacement factor. Dominate the
kernel by \(C\phi_{\sqrt2 b}\). This Gaussian convolution is positive as
a quadratic form, so Cauchy--Schwarz reduces the sum over transfers to
their individual terms. For a transfer \(i\to j\) with \(s\ge1\), the
\(L^2\) contraction bounds that term by \(\int w^2 1_{F_{ij}}\),
paid by (ET3). For \(0\le s<1\), integrate \(w1_{F_{ij}}\)
in the normal direction to obtain a section integral \(m_w(z)\).
Bounding the normal kernel by its supremum leaves the
tangential Gaussian convolution, whose quadratic term is
at most \(C\int m_w(z)^2\,dz\). In these orthogonal coordinates,
\(w(z,s)=w(z,0)\exp(-e^2s^2/(4a^2))\).
The section inequality \(m(z)^2\le2\int_0^1 s 1_{F_{ij}}(z,s)\,ds\)
and the variation of the normal factor on \([0,1]\) give
\(m_w(z)^2\le C\int_0^1 s w(z,s)^2 1_{F_{ij}}(z,s)\,ds\).
This is paid by (ET3). One more Cauchy--Schwarz combines thick
and thin parts. The same argument with \(w=1\) bounds the
unweighted interaction by \(C_k L_b(U)\).
These bounds give the second estimate of (ET1).

## Unrestricted changes on a small noise-scale ball

For every fixed \(k\ge3\), there is \(r_k>0\) such that arbitrary
ordinary competitors supported in \(B_r\), \(0<r\le r_k\), satisfy
\[
         \Delta_1(U)\ge L_1(U),\qquad
         \Delta_e^{\mathrm G}(U)\ge L_e^{\mathrm G}(U)
                                                                 \tag{ET7}
\]
for all sufficiently small \(e\), with no mass constraint.
Indeed \(s\le r\le1\) on each transfer. The projections of the
transfer sets onto their pair planes have area at most \(C_kr^{d-1}\).
Cauchy--Schwarz on each projection and the section inequality from
(ET4) give \(|D|^2\le C_k r^{d-1}L_1(U)\).
The same estimate with \(L_e^{\mathrm G}\) follows from the first
part of (ET1), with the ambient support bound \(R=1\), and (ET3).
On these supports the two kernels are uniformly bounded, so
\(|Q|\le C_k|D|^2\). Choose \(r_k\) with \(C_k r_k^{d-1}\le1\)
and use the definition of the deficit. Here the physical Gaussian
radius is \(r e\).

## A sufficient Euclidean coercivity estimate

Suppose a class of ordinary changes is closed under spatial
restriction to balls and positive dilations. For every \(R\ge1\),
assume the Euclidean estimate
\[
       \Delta_1(U)\ge c R^{-\alpha}L_1(U)                       \tag{ET5}
\]
on that class when \(U=T\) off \(B_R\), for constants \(c>0\) and
\(0\le\alpha<5\). This hypothesis applies without mass restrictions.
Dilation gives the same inequality for \(b\) near one with a
uniform constant, by allowing radius \(\max(R/b,1)\) after
rescaling. The support of each \(U_u\) is in a ball of radius
at most \(R\), because \(w\) is radial and decreasing. Apply
(ET5) to the integral in (ET1) and absorb the error when
\(e^2(R+1)^{\alpha+1}\) is sufficiently small. The result is
\[
          \Delta_e^{\mathrm G}(U)\ge c' R^{-\alpha}
                                      L_e^{\mathrm G}(U).
                                                               \tag{ET6}
\]
For every fixed \(C>0\) the absorption is uniform through the
rescaled radius \(R=C[e^{-1}\log(1/e)]^{1/3}\) near one, since
\(e^2(R+1)^{\alpha+1}\to0\) there. This is physical radius
\(C e^{2/3}[\log(1/e)]^{1/3}\). On that scale \(R^{-\alpha}\)
dominates \(e^\beta\) for any \(\beta>\alpha/3\), so a balanced
truncation on this scale with deficit error \(o(e^\beta L)\)
would transfer (ET5) to the original maximizer.

The available [CT19–CT20 localization](THREE_CELL_EXTERIOR_TUBES.md)
has cutoff radius of order \([e^{-1}\log(1/e)]^{1/3}\) in
physical coordinates. Its resulting noise-rescaled radius is
of order \(e^{-4/3}[\log(1/e)]^{1/3}\). At that radius the
error bound in (ET1) is absorbable against \(R^{-\alpha}\) only
when \(\alpha<1/2\). The
[balanced ordered-ball examples](TWO_BLOCK_SIMPLEX_JUNCTION.md#sharpness-of-the-ordered-ball-exponent)
give ratios at most \(C R^{-2}\) even in the Euclidean problem,
so a uniform estimate (ET5) cannot have \(\alpha<1/2\).
A separate localization or a sharper comparison sensitive to
the geometry of the defects is needed at the larger radius.

Conversely, even nonnegative sign in the unrestricted Euclidean
comparison is necessary for exact balanced Gaussian optimality
along any sequence with \(e\downarrow0\). To see this, take a fixed
compact ordinary Euclidean defect \(U\) and scale it to \(ep\).
Its Gaussian mass differences are \(O(e^d)\). They can be repaired
in \(k-1\) fixed, mutually disjoint small patches in the interiors
of simplex facets whose label pairs form a tree. In each patch
transfer the signed required flow between the two cells, using a
one-sided normal strip of thickness \(O(e^d)\) and a fixed positive
tangential section. These strips are disjoint from the scaled
defect for small \(e\), and the required masses are attained by
continuity and strictly positive normal density on each patch.
The repair score cost is \(O(e^{2d-1})\), from the reflected
\(O(e^{-1})\) Lipschitz score bound and the strip width and mass.
The [near-one uncentered transfer estimate](THREE_CELL_EXTERIOR_TUBES.md)
bounds the pure repair interaction by the same order. Alternatively,
the diagonal section proof of (ET4) applied to the scaled repair
strips in (ET2) gives this estimate directly. In scaled coordinates
the strip widths are \(O(e^{d-1})\) and their volumes are \(O(1)\),
so the total normal moment is \(O(e^{d-1})\). The physical scaling
factor is \(e^d\).
The core interaction is \(O(e^d)\) by (ET4) and scaling.
Cauchy--Schwarz for the Gaussian interaction gives a mixed error
\(O(e^{(3d-1)/2})=o(e^d)\), since \(d\ge2\). Dividing the deficit
of the repaired partition by \(e^d\) and passing to the limit
in (ET1) on a fixed support gives
\(\phi_1(0)\Delta_1(U)\ge0\).
