# A distant bounded pocket requires exterior interaction near its contracted center

Fix \(0<\rho<1\), write \(\sigma=\sqrt{1-\rho^2}\), and
let \(A\) be an ordinary planar three-cell Gaussian
noise-stability maximizer at positive prescribed masses.
Let \(S\) be a bounded positive-area component surrounded by
one label \(j\), as in the bounded-component analysis.
The component may contain both non-\(j\) labels. For
every fixed \(D<\infty\), if \(S\subset B_D(p)\) and
\(|p|\) is sufficiently large, some remaining non-\(j\)
material must lie in a ball of radius depending only on
\(\rho,D\) about \(\rho p\). In fact the ball contains at
least half of a specified exterior interaction measure.
Both the large-radius threshold and the interaction radius
can be chosen uniformly on compact positive-correlation bands.

The estimate uses the first extra-coordinate stability
inequality with a direction-dependent speed weight. The
weight cancels the Gaussian drift in the radial projection
and makes the residual same-pocket interaction nonpositive.
This locates the exterior interaction at bounded distance
from the contracted center even when the distance of the
pocket from the origin tends to infinity.
[Unweighted exterior covariance bound, IC1–IC4](THREE_CELL_ISOLATED_SUBCLUSTERS.md#covariance-bound-from-the-first-extra-coordinate-block)

## Setup and statement

Fill \(S\) with label \(j\) to obtain the partition \(B\),
and set \(H=1_A-1_B\). The affected interfaces \(D_S\)
include the interior walls of \(S\) and its exterior
boundary. On them the jump of \(H\) is the original
label jump. This holds for bounded components of
\(\mathbb R^2\setminus\overline{O_j}\), where \(O_j\)
is the open strict-winner set. The argument also applies
to any bounded fill with this distributional jump property.
Let \(E\) be the non-\(j\) material of \(B\). Write
\(p_\rho(x,y)=\phi_2(x)\phi_2(y)M_\rho(x,y)\) for the
joint density in the two planar variables. Set
\[
 W_{\rm in}(x,y)=H(x)\cdot H(y),\qquad
 W_{\rm out}(x,y)=H(x)\cdot1_B(y)+1_S(x).
\]
The first weight is nonnegative and supported on \(S^2\).
The second is supported on \(S\times E\) and equals
one or two there, according to the two inside labels.
[Bounded fills and their interface measure](THREE_CELL_ISOLATED_SUBCLUSTERS.md)

Take a unit vector \(v\), its perpendicular unit vector
\(w\), and \(0<a\le b<\infty\) and \(L<\infty\) such that
\[
  a\le x\cdot v\le b,\qquad |x\cdot w|\le L
          \quad(x\in S\text{ a.e.}).
\]
Define
\[
 c_\rho=\frac{1+\rho^2}{2\rho},\quad
 \alpha_\rho=\frac{\sigma^2}{2\rho},\quad
 \zeta=a-\rho b,\quad
 C_{\rho,L}=\frac{\rho}{\sigma^2}
                  +\frac{(1+\rho)^2L^2}{\sigma^4}.
\]
Assume
\[
                        \zeta\ge\sigma,\qquad b\le c_\rho a.  \tag{DC1}
\]
Let \(g(s)=1/\phi_1(s)\) and use the
exterior interaction law proportional to
\(W_{\rm out}(x,y)g(x\cdot v)p_\rho(x,y)\,dx\,dy\).
Its total mass \(M_{\rm out}\) is positive under
(DC1). With \(\mathbb E_{\rm out}\) expectation under
this law and \(Z=Y-\rho X\), the bounds are
\[
 \mathbb E_{\rm out}(Z\cdot v)^2<\sigma^2,\qquad
 \sqrt{\mathbb E_{\rm out}(Z\cdot w)^2}
     \le \alpha_\rho L+
             \sqrt{\sigma^2+\alpha_\rho^2 L^2+
                       \frac{\sigma^5 C_{\rho,L}}{\rho\zeta}}
             =:K_{\rho,L,\zeta}.                              \tag{DC2}
\]
In particular for \(S\subset B_D(p)\), take
\(R=|p|\), \(v=p/R\), \(a=R-D\), \(b=R+D\), and \(L=D\).
Condition (DC1) holds if \(R>D\) and
\[
 R\ge\max\left\{\frac{\sigma+(1+\rho)D}{1-\rho},
                 \frac{(c_\rho+1)D}{c_\rho-1}\right\}.
\]
Minkowski's inequality and (DC2) give
\[
 \bigl(\mathbb E_{\rm out}|Y-\rho p|^2\bigr)^{1/2}
       \le \sqrt{\sigma^2+K_{\rho,D,\zeta}^2}+\rho D.
\]
The ball centered at \(\rho p\) with radius \(\sqrt{2}\)
times the right-hand side contains at least half of the
exterior interaction law in \(Y\). Under (DC1) the radius
bound is uniform in \(p\).

## The exact transport identity

Put \(s=x\cdot v\) and let \(\Phi\) denote the standard
Gaussian cumulative distribution. If \(S\) is positively
separated from \(E\), transport the whole pocket by
\(\Phi(s)\mapsto\Phi(s)+t\) in its first coordinate,
leaving its transverse coordinate unchanged. For small
positive or negative \(t\) the transported pocket has no
overlap with \(E\). The transport preserves Gaussian volume
of every inside label exactly, since \(\Phi(s)\) is a uniform
coordinate for the one-dimensional marginal. Therefore the
derivative of stability at \(t=0\) vanishes. The transport
velocity is \(g(s)v\). Symmetry of the self-interaction
reduces its derivative to twice the derivative in its first
endpoint. Adding \(1_S(x)\) to the exterior label pairing
integrates to a constant under this transport. If \(Z_v
=(y-\rho x)\cdot v\), the resulting equation is
\[
 \iint (W_{\rm in}+W_{\rm out})g(s) Z_v p_\rho(x,y)\,dx\,dy=0.
                                                               \tag{DC3}
\]

The same equation holds when the closure of the pocket
meets \(E\). Its infinitesimal proof uses stationarity and
does not require a disjoint finite transport in that case.
Indeed, with the stationary prices \(\beta\), the
distributional derivative \(D_vH\) gives
\(\int (T_\rho1_A-\beta)\cdot g(s)\,d(D_vH)\,\phi_2=0\)
because the price-adjusted pair scores agree on each
affected wall. In Cartesian coordinates,
\(D_v(g(s)\phi_2)=0\). Integration by parts therefore
gives \(\int H\cdot g(s)D_v T_\rho1_A\,d\gamma_2=0\).
The derivative of the conditional Gaussian kernel supplies a
factor \(\rho Z_v/\sigma^2\). Substitute \(1_A=1_B+H\)
and use \(T_\rho1=1\) to add \(1_S(x)\) in the exterior
weight, giving (DC3). Compact cutoffs equal to one near
the pocket justify testing the compactly supported
distribution \(D_vH\) with these coefficients.
[Stationary pair scores and normal translation identity, G2](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

On \(S^2\), (DC1) gives \(Z_v\ge\zeta>0\).
Writing \(M_{\rm in}=\iint W_{\rm in}g(s)p_\rho\),
we have \(M_{\rm in}>0\). In particular (DC3) forces
\(M_{\rm out}>0\) and
\[
                -M_{\rm out}\mathbb E_{\rm out} Z_v
                         \ge \zeta M_{\rm in}.                 \tag{DC4}
\]

## A weighted spectator variation

Let \(f(x)=\sqrt{g(x\cdot v)}\) and let \(u\) be a constant
unit direction. Use the wall speed \(f(x)(u\cdot N)\) on
\(D_S\) and zero elsewhere in the first extra-coordinate
form \(Q^{(1)}\). This speed is bounded on the bounded
affected wall. It can be approximated by compact smooth
speeds on regular interiors in \(L^2(a_e\,d\nu)\), using
local finiteness of the measure and the full
surface-measure coverage by regular pair patches.
There is no linear volume constraint on these speeds
in \(Q^{(1)}\), because the independent spectator
factor is an odd function. Continuity of the form in
the weighted space follows from positivity of the
interaction and the operator bound in (G1).
[First extra-coordinate form, G1, and weighted domain](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)
[Weighted approximation for a bounded fill, IC1 and PC4](THREE_CELL_ISOLATED_SUBCLUSTERS.md)

On the affected wall, the normal translation identity
applied to \(1_A=1_B+H\) reads
\[
 a_e(u\cdot N)=\rho K_{D_S}(u\cdot N)
                         -d_e\cdot D_u T_\rho1_B.
\]
Multiply by \(g(s)(u\cdot N)\) and integrate on \(D_S\).
Symmetrizing the same-pocket interaction, then subtracting
the interaction of the weighted speed, leaves the factor
\((f(x)-f(y))^2/2\). For the \(B\) term use
\(D_u H=-\sum_{e\in D_S}d_e(u\cdot N)\,ds_e\)
and integrate by parts. Thus
\[
\begin{split}
 0\le Q^{(1)}(f\,1_{D_S}(u\cdot N))
   ={}&-\int H\cdot
       \left[g(s)D_{uu}T_\rho1_B+
           (D_ug(s)-(x\cdot u)g(s))D_uT_\rho1_B\right]d\gamma_2\\
      &+\mathcal R_u,\\
 \mathcal R_u={}&\frac{\rho}{2}
       \iint W_{\rm in}(x,y)\,
       D_u^xD_u^y\left[(f(x)-f(y))^2p_\rho(x,y)\right]dx\,dy .
\end{split}                                                    \tag{DC5}
\]
The expression for \(\mathcal R_u\) follows by integrating
the two distributional jump measures of \(H\) by parts.
All test kernels can be cut off to a neighborhood of
\(\overline S^2\) before taking these derivatives.

### Radial projection

Take \(u=v\) in (DC5). Since \(D_vg(s)=s g(s)\), its
first integral equals
\[
 \frac{\rho^2}{\sigma^4} M_{\rm out}
                    \bigl(\sigma^2-\mathbb E_{\rm out}Z_v^2\bigr).
                                                               \tag{DC6}
\]
Here and below \(1_S(x)\) can be added to the exterior
weight because it corresponds to differentiating
\(T_\rho1=1\).

For the residual, write \(t=y\cdot v\) and expand
\((f(x)-f(y))^2p_\rho\) as the sum of \(g(s)p_\rho\)
and \(g(t)p_\rho\) minus \(2f(x)f(y)p_\rho\).
The ratios of their respective mixed derivatives
to these three positive kernels, with the factor
\(-2\) taken out, are
\[
 \frac{\rho}{\sigma^2}
           -\frac{\rho(t-\rho s)^2}{\sigma^4},\qquad
 \frac{\rho}{\sigma^2}
           -\frac{\rho(s-\rho t)^2}{\sigma^4},\qquad
 \frac{\rho}{\sigma^2}
           +\frac{\rho^2(t-c_\rho s)(s-c_\rho t)}{\sigma^4}.
\]
Indeed the radial log derivative of \(p_\rho\) in
its first variable is \((\rho t-s)/\sigma^2\);
derivatives of \(g\) and \(f\) add \(s\) and \(s/2\),
respectively, to that log derivative. By (DC1), the
first two ratios are nonpositive and the third is
strictly positive. Thus the mixed derivative of
the squared-difference kernel is strictly negative
throughout \(S^2\). As \(W_{\rm in}\) is positive
there, \(\mathcal R_v<0\). Equations (DC5)–(DC6) give
the first inequality in (DC2). Together with
(DC4), it also gives
\[
                       M_{\rm in}/M_{\rm out}\le\sigma/\zeta . \tag{DC7}
\]

### Transverse projection

For \(u=w\), \(D_wg=0\). The first integral in (DC5)
is
\[
 \frac{\rho^2}{\sigma^4}M_{\rm out}
   \left[\sigma^2-\mathbb E_{\rm out}Z_w^2+
              \frac{\sigma^2}{\rho}
                 \mathbb E_{\rm out}((X\cdot w)Z_w)\right],
\]
where \(Z_w=Z\cdot w\). In \(\mathcal R_w\) the tangential
mixed derivative only acts on \(p_\rho\). Its ratio to
\(p_\rho\) is at most \(C_{\rho,L}\) on \(S^2\):
the mixed log derivative is \(\rho/\sigma^2\), and the
two ordinary log derivatives have magnitudes at most
\((1+\rho)L/\sigma^2\). Since \(f\ge0\),
\((f(x)-f(y))^2\le g(s)+g(t)\). Symmetry of
\(W_{\rm in}p_\rho\) therefore gives
\(\mathcal R_w\le \rho C_{\rho,L}M_{\rm in}\).
Use (DC7) and Cauchy–Schwarz on the bracket to get,
for \(q_w=\sqrt{\mathbb E_{\rm out}Z_w^2}\),
\[
 q_w^2-2\alpha_\rho L q_w
    \le \sigma^2+\frac{\sigma^5 C_{\rho,L}}{\rho\zeta}.
\]
Solving this quadratic gives the second bound in (DC2).

The same estimate applies to pockets in balls \(B_D(p)\)
with \(D=o(|p|)\). In that case (DC1) holds eventually,
\(\zeta\) grows linearly in \(|p|\), and the half-mass
radius in (DC2) is \(O(1+D)=o(|p|)\), uniformly on a
compact positive-correlation band. A pocket deep inside
a fixed exterior sector thus cannot have all its other
nonhost material confined to the other sectors, because
the distance from \(ho p\) to those sectors grows
linearly in \(|p|\).

The condition requires inward exterior nonhost material
for a far-out bounded pocket. It allows radial chains of
such material and imposes no connectedness hypothesis on
the exterior. In particular it does not exclude an
arbitrary compact loop configuration at a correlation
endpoint.
