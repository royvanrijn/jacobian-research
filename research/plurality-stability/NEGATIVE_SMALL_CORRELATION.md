# Gaussian simplex optimality at small negative correlation

Write \(T=(1_{T_i})_{i=1}^k\) for a central regular-simplex assignment
in \(\mathbb R^d\), extended in the unused coordinates. Let
\(u=(1/k,\ldots,1/k)\), and let \(f:\mathbb R^d\to\Delta_k\) be a
measurable fractional assignment. Its Gaussian noise stability is
\[
S_\rho(f)=\sum_i\langle f_i,T_\rho f_i\rangle_{L^2(\gamma_d)},\qquad
T_\rho g(x)=\mathbb E\,g(\rho x+\sqrt{1-\rho^2}\,Z).
\]

**Fixed-cell theorem.** For every fixed \(k\ge2\) there is \(r_k>0\)
such that, for every \(d\ge k-1\), \(0<r\le r_k\), and measurable
\(f:\mathbb R^d\to\Delta_k\),
\[
S_{-r}(f)\ge S_{-r}(T).
\tag{N1}
\]
Equality holds exactly for cylindrical regular-simplex assignments
up to orthogonal transformations, labels and null sets. The masses
\(\int f_i\,d\gamma_d\) are free. This is an internal written deduction
from the [equal-cell first-moment input and moment alignment estimate
(M2)](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md#from-moments-to-regular-score-regret).

For three cells, Heilman's *Sharp Hardness for MAX-3-CUT and Quantum
MAX-CUT*, [arXiv:2608.00333](https://arxiv.org/abs/2608.00333),
Theorem 7.1, gives the full negative voting interval
\([-1/2,0)\) for ordinary measurable partitions with free masses in
every dimension \(d\ge2\). Its proof minimizes a bilinear functional
over pairs with corresponding masses equal and applies planar
dimension reduction and circular rearrangement. Section 1.4 identifies
the rearrangement obstruction for four or more cells. The argument
below uses first-moment geometry in the \(k-1\) active coordinates.
The threshold in (N1) may depend on \(k\).

## Localization with free masses

Let \(B_f\) have rows \(\int xf_i(x)\,d\gamma_d(x)\), and write
\(p=\int f\,d\gamma_d\). The total Hermite energies \(W_m(f)\) satisfy
\(W_0=\|p\|^2\), \(W_1=\|B_f\|_F^2\), and
\(\sum_m W_m=\int\sum_i f_i^2\,d\gamma_d\le1\).
The Mehler eigenvalues give
\[
\big|S_{-r}(f)-\|p\|^2+r\|B_f\|_F^2\big|\le r^2.
\tag{N2}
\]
Denote the simplex's first energy by
\(F_k=(\mathbb E\max_{i\le k}Z_i)^2/(k-1)\), with independent standard
normals \(Z_i\). A hypothetical sequence with \(r\downarrow0\) and
\(S_{-r}(f)\le S_{-r}(T)\) satisfies
\[
\|p-u\|^2+r(F_k-\|B_f\|_F^2)\le2r^2.
\tag{N3}
\]
In particular \(\|p-u\|^2\le r+2r^2\) and
\(\|B_f\|_F^2\ge F_k-2r\).

Rebalance \(f\) by removing the fraction \((p_i-1/k)/p_i\) of \(f_i(x)\)
for each \(p_i>1/k\). Distribute the pointwise sum of the removals
among labels with deficient mass in proportion to their deficits.
This gives an exactly balanced assignment \(g\) with
\(\sum_i\int|f_i-g_i|\,d\gamma_d=\sum_i|p_i-1/k|\).
If \(p=u\), take \(g=f\). Hermite projection gives
\(\|B_f-B_g\|_F^2\le\sum_i\|f_i-g_i\|_2^2\le\sum_i|p_i-1/k|\).
The feasible-label inequality in (M1) holds for balanced fractional
assignments by averaging the pointwise winning-score inequality
with their label weights. Thus the same deduction (M2) applies to
\(g\), giving \(\|B_g\|_F^2\le F_k\) and moment alignment.
By (N3), \(\|B_g\|_F^2\to F_k\).

Choose a regular moment matrix \(B_*\) closest to \(B_f\).
The triangle inequality and (M2) give \(\|B_f-B_*\|_F\to0\)
with estimates independent of \(d\). The associated simplex has score regret
\[
\int\left(\max_i b_i^*\cdot x-\sum_i f_i(x)b_i^*\cdot x\right)
\,d\gamma_d=F_k-\langle B_f,B_*\rangle\longrightarrow0 .
\]
Outside the strip of width \(a>0\) between the top two regular scores,
the regret controls \(a\) times the wrong-label weight. The strip
has measure at most \(C_k a\) by the union bound on one-dimensional
Gaussian score differences. Consequently
\(\sum_i\|f_i-1_{T_i}\|_1\to0\) for the aligned simplex.

Rotate its active space to \(E=\mathbf1^\perp\subset\mathbb R^k\)
and write points as \((x,z)\in E\oplus\mathbb R^{d-k+1}\).
The regular matrix is \(B_*=a_*[P,0]\), where \(P\) projects onto \(E\)
and \(a_*=\sqrt{F_k/(k-1)}\).
Differentiating the squared distance to the moment orbit shows that
the active block of \(B_f-B_*\) is symmetric and its transverse block
vanishes. Active rotations yield every skew active matrix as a
tangent, and mixed rotations yield every centered transverse matrix
as a tangent. Put \(\ell=\sqrt{2}a_*\).

## Exact negative-correlation scores

For \(i\ne j\), write \(n_{ij}=(b_i^*-b_j^*)/\ell\),
\(H_{ij}=n_{ij}^\perp\) in \(E\), \(F_{ij}=T_i\cap T_j\) in \(E\),
and \(t=n_{ij}\cdot x\). Reflection \(G_{ij}\) in \(H_{ij}\) exchanges
\(T_i,T_j\). Set
\[
V_{r,ij}(x)=-T_{-r}(1_{T_i}-1_{T_j})(x),\qquad m=-rx/\sqrt{1-r^2}.
\]
For \(x\in T_i\), put \(q=-m\cdot n_{ij}\ge0\).
The cones are invariant under positive dilations. Changing variables
by \(G_{ij}\) in the \(T_j\) integral gives
\[
V_{r,ij}(x)=\int_{T_i}\phi_{k-1}(y-G_{ij}m)
\big(1-e^{-2q(y\cdot n_{ij})}\big)\,dy .
\tag{N4}
\]
The integrand is nonnegative. On the ball
\(D_i=B(b_i^*/|b_i^*|,1/4)\subset T_i\),
\(y\cdot n_{ij}\ge1/\sqrt{2}-1/4>0\) and \(|y|\le5/4\).
There the density is bounded below by
\(\phi_{k-1}(y)\exp(-5|m|/4-|m|^2/2)\).
For \(r\le1/16\) the exponent is at least \(-|x|^2/32-C\)
by completing the square. Since \(q\ge rt\), the restriction to \(D_i\)
gives
\[
\frac{V_{r,ij}(x)}r\ge c_k e^{-|x|^2/32}\min(t,1).
\tag{N5}
\]
The score vanishes on \(H_{ij}\), and \(V_{r,ij}/r\to\ell t\)
locally in \(C^1\). This follows by differentiating the Mehler kernel
at zero correlation, where its derivative is \(x\cdot y\).
On compact \(x\)-sets the correlation derivative and its first
\(x\)-derivatives have integrable Gaussian majorants in \(y\)
for \(|\rho|\le1/16\).

Put \(e=f-T\), \(h_{ij}=1_{T_i}f_j\), and
\(L=\sum_{i\ne j}\int h_{ij}V_{r,ij}\,d\gamma_d\).
Here \(0\le h_{ij}\le1\), its support has \(t\ge0\), and
\(e_i=\sum_{j\ne i}(h_{ji}-h_{ij})\) almost everywhere.
Expanding the stability gives
\[
S_{-r}(f)-S_{-r}(T)
 =2L+\|p-u\|^2-r\|B_f-B_*\|_F^2+\sum_{m\ge2}(-r)^m W_m(e).
\tag{N6}
\]
By (N5), \(L>0\) whenever \(f\) differs from the aligned simplex
on a set of positive measure.

## Transfer estimates in the active coordinates

Use \(w(x)=(1+|x|^2)e^{|x|^2/16}\) on \(E\). For an active transfer
\(h:E\to[0,1]\) supported on \(t\ge0\), set
\[
J(h)=\int \min(t,1)e^{-|x|^2/32}h(x)\,d\gamma_{k-1}(x).
\]
There are \(C_k<\infty\) and \(\beta_k(R)\to0\) such that
\[
\left(\int w h\,d\gamma_{k-1}\right)^2\le C_kJ(h),\qquad
\left(\int_{|x|>R} w h\,d\gamma_{k-1}\right)^2\le\beta_k(R)J(h).
\tag{N7}
\]
For \(t\ge1\), Cauchy–Schwarz bounds the square by \(J(h)\) times
the integral of \(w^2e^{|x|^2/32}\) on the integration region.
For \(0\le t<1\), write \(x=y+tn_{ij}\) and
\(m_0(y)=\int_0^1 h(y+tn_{ij})\phi_1(t)\,dt\).
A density bounded by \(\phi_1(0)\) on the positive axis has first
moment at least its squared mass divided by \(2\phi_1(0)\), by filling
closest to zero. Hence \(m_0(y)^2\le 2\phi_1(0)\int_0^1
t h(y+tn_{ij})\phi_1(t)dt\).
Here \(e^{-|x|^2/32}\ge e^{-(|y|^2+1)/32}\) and
\(w(x)\le2e^{1/16}w(y)\). Apply Cauchy–Schwarz in \(y\) with the
normal cost weight, leaving an integral of \(w(y)^2e^{|y|^2/32}\)
up to a constant. The exterior integral uses only
\(|y|>\sqrt{R^2-1}\). The resulting Gaussian tails tend to zero.
In the zero-dimensional tangential space the last tail is zero for \(R>1\).

Set \(s^2=L/r\). Along the localized sequence, \(s\to0\).
Indeed, the normalized score bound \(|V_{r,ij}(x)|/r\le C_k w(x)\)
follows from the Mehler derivative bound used in
[the positive-correlation argument](SMALL_POSITIVE_CORRELATION_INTERVAL.md#weighted-tails-in-dimension-k-1).
It holds for both signs since the kernel bound uses \(|\rho|\).
The wrong-label weights tend to zero in \(L^1\), and their integrals
against \(w\) tend to zero by truncation in the active coordinates.

Apply (N7) to \(\bar h_{ij}(x)=\mathbb E_z h_{ij}(x,z)\).
The measures \(\bar h_{ij}\gamma_{k-1}/s\) have bounded weighted mass
and uniformly vanishing weighted tails. Applying the global (N7)
on each active slice at fixed \(z\) and integrating the squares
in \(z\) also gives \(\sum_i\|J_i\|_2^2\le C_k\), where
\(J_i(z)=s^{-1}\int w(x)|e_i(x,z)|\,d\gamma_{k-1}(x)\).
These estimates use \(L/r\ge c_k\sum_{i\ne j}\mathbb E_z J(h_{ij}(\cdot,z))\)
from (N5).

Let \(u_{i,p,q}\) be the part of \(e_i/s\) with active Hermite degree
\(p\) and transverse degree \(q\). At \(b=1/16\) the active Mehler
kernel is bounded by \(C_k w(x)w(y)\). Applying it on each fixed
transverse slice gives
\[
\sum_{i,p,q} b^p\|u_{i,p,q}\|_2^2\le C_k\sum_i\|J_i\|_2^2 .
\]
For \(p+q\ge2\) and \(0<r\le b\),
\(r^{p+q-1}\le(r/b^2)b^p\), so
\[
\left|\sum_{m\ge2}(-r)^m W_m(e)\right|\le C_k rL .
\tag{N8}
\]
The first-moment term in (N6) is \(O_k(L)\) by the weighted bound
on the averaged transfers and the vanishing transverse block.
If the difference in (N6) is nonpositive, then
\(\|p-u\|^2\le C_k L\), so \((p-u)/s\to0\).

## Limiting flux on the facets

Here \(\nu_{ij}=\phi_1(0)\gamma_{k-2}\) on \(H_{ij}\).
On each compact part of the plane the normalized normal integrals
\[
\theta_{r,ij}(y)=\frac1{s\phi_1(0)}\int_0^\infty
\bar h_{ij}(y+tn_{ij})\phi_1(t)\,dt
\]
are bounded in \(L^2(\nu_{ij})\). For \(t<1\) use the normal
mass-cost inequality in (N7), with the cost weight bounded below
on the compact cylinder. For \(t\ge1\) use Cauchy–Schwarz in \(t\)
with reciprocal weight \(e^{t^2/32}\).
Pass to a subsequence of locally weak limits \(\theta_{ij}\ge0\).
For each fixed \(b_0>0\), the part with \(t\ge b_0\) has mass
\(O(s^2)\) on each fixed ball by (N5).
The weighted tail bound from (N7) then shows that its normalized
mass tends to zero on all of \(E\).
Uniform continuity on compact sets identifies the normalized ordered
transfer limit with \(\theta_{ij}\nu_{ij}\), supported on \(F_{ij}\)
because the transfers originate in the closed cone \(T_i\).
Junctions have zero \(\nu_{ij}\)-measure when \(k\ge3\).
For \(k=2\) the plane is a point of mass \(\phi_1(0)\) and there is no junction.

On any compact tangential cylinder with \(0\le t\le b_0<1\),
the local score convergence gives \(V_{r,ij}/r\ge(\ell-o(1))t\).
Normal integrals truncated to \(t\le b_0\) have the same local
weak limits: the normalized difference tends to zero in \(L^1\)
and is locally bounded in \(L^2\), so bounded tests suffice by density.
The one-dimensional first-moment inequality for the truncated
integrals gives the squared norm bound after division by \(s^2\).
Sum the ordered transfer contributions using their nonnegative
costs, apply weak lower semicontinuity, and exhaust the planes to obtain
\[
1\ge\frac{\ell}{2}\sum_{i\ne j}\|\theta_{ij}\|_{L^2(\nu_{ij})}^2
\ge\frac{\ell}{2}\sum_{i<j}\|v_{ij}\|_{L^2(\nu_{ij})}^2,\qquad
v_{ij}=\theta_{ji}-\theta_{ij}.
\tag{N9}
\]
The last inequality uses \(\theta_{ij}\theta_{ji}\ge0\).
The normalized signed transfers have limits
\(\mu_i=\sum_{j\ne i}v_{ij}\nu_{ij}|_{F_{ij}}\), with \(v_{ji}=-v_{ij}\).
Weighted tails pass masses and first moments to the limit.
Thus \(\mu_i(E)=0\) by \((p-u)/s\to0\), and the moment matrix \(H\)
with rows \(\int x\,d\mu_i(x)\) is symmetric on \(E\) in the moment gauge.

The [convex facet covariance estimate (L2)](SMALL_POSITIVE_CORRELATION_INTERVAL.md#an-interface-estimate-in-the-moment-gauge)
applies to these balanced speeds in the symmetric gauge. It gives
\[
\|H\|_F^2\le(1-1/k)\ell\sum_{i<j}\|v_{ij}\|_{L^2(\nu_{ij})}^2
\le2(1-1/k).
\]
Divide (N6) by \(L=r s^2\), discard the nonnegative mass term,
and use (N8) and convergence of the normalized moment matrix.
Nonpositivity of the original difference would give \(\|H\|_F^2\ge2\).
The contradiction excludes every such sequence and proves (N1)
with equality rigidity and a threshold depending only on \(k\).

## The large-alphabet endpoint

The localization uses fixed-\(k\) moment alignment, and the exterior
and smoothing bounds also depend on \(k\). The remaining target is
\[
\exists M\ge3\quad \forall k\ge M,\ \forall d\ge k-1,\
\forall\,0<r\le 1/(k-1),\ \forall f:\mathbb R^d\to\Delta_k
\quad S_{-r}(f)\ge S_{-r}(T).
\]
A lower bound \(r_k\ge1/(k-1)\) for all sufficiently large \(k\)
or a separate comparison on the portion outside the fixed-cell
neighborhood would complete this interval.
