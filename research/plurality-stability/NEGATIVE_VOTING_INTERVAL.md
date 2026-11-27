# The negative voting interval for large alphabets

Let \(k\ge3\), \(d\ge k-1\), and \(0<r\le1/(k-1)\). For a measurable
assignment \(f:\mathbb R^d\to\Delta_k\), set \(p_i=\int f_i\,d\gamma_d\) and
write
\[
 S_{-r}(f)=\sum_i\mathbb E[f_i(X)f_i(Y)],
 \qquad \operatorname{Cov}(X,Y)=-rI_d.
\]
Let \(T_k\) denote the central regular simplex partition, extended
cylindrically to dimension \(d\). The exact target is the existence of
one \(M\) such that \(S_{-r}(f)\ge S_{-r}(T_k)\) for all \(k\ge M\),
all \(0<r\le1/(k-1)\), all \(d\ge k-1\), and all such \(f\).
The comparison below is uniform in every one of these parameters.
It leaves an additive error of order \(r/k\).

## V1. Uniform full-interval comparison and mass control

Write \(u=1/k\), let \(\Phi,\phi\) denote the standard normal CDF and
density, and let
\[
 h_r(p)=J_{-r}(p,p),\qquad
 J_{-r}(p,q)=\Pr[Z\le\Phi^{-1}(p),\,W\le\Phi^{-1}(q)]
\]
for a standard normal pair with correlation \(-r\). Use continuity at
\(p,q=0,1\). Then
\[
 S_{-r}(f)\ \ge\ k h_r(u)+\frac1{72}\sum_i(p_i-u)^2.
 \tag{V1}
\]
Moreover, uniformly for \(0<r\le1/(k-1)\) as \(k\to\infty\),
\[
 S_{-r}(T_k)-k h_r(u)
   =\bigl(2(1-\gamma_{\rm E})+o(1)\bigr)\frac r k,
 \tag{V2}
\]
where \(\gamma_{\rm E}\) is the mean of the standard Gumbel law.
In particular,
\[
 0\le S_{-r}(T_k)-\inf_{f:\mathbb R^d\to\Delta_k}S_{-r}(f)
       \le \bigl(2(1-\gamma_{\rm E})+o(1)\bigr)\frac r k .
 \tag{V3}
\]
Every \(f\) with \(S_{-r}(f)\le S_{-r}(T_k)\) satisfies
\(\sum_i(p_i-u)^2=O(r/k)\) with an absolute constant. The improvement
of the simplex over independent labels is
\[
 1/k-S_{-r}(T_k)=(2+o(1))r\log(k)/k
\]
on the same interval, so (V3) is smaller by a factor of order
\(1/\log k\).

For (V1), the negative-correlation functional Borell inequality gives
\(\mathbb E[J_{-r}(f_i(X),f_i(Y))]\ge J_{-r}(p_i,p_i)\).
Also \(J_{-r}(a,b)\le ab\). Indeed, condition on the first coordinate:
the conditional second-coordinate lower-tail probability increases in
that coordinate, while the first-coordinate lower-tail indicator decreases.
Their covariance is nonpositive, by the independent-copy formula for
covariance. This proves the cellwise lower bound for the actual bilinear
product. The imported statement is [Mossel–Neeman, *Robust optimality of
Gaussian noise stability*, Theorem 1.2 in arXiv:1210.4126v3, the reversal
of (1.2) at negative correlation](https://arxiv.org/pdf/1210.4126v3).

Set \(a=\sqrt{(1+r)/(1-r)}\). Differentiating both thresholds in the
bivariate normal distribution gives
\[
 h_r'(p)=2\Phi(a\Phi^{-1}(p)),\qquad
 h_r''(p)=2a\exp\bigl(-(a^2-1)\Phi^{-1}(p)^2/2\bigr).
\]
On \([u/2,1-u/2]\), the Chernoff bound gives \(\Phi^{-1}(p)^2\le
2\log(2k)\), and \(a^2-1=2r/(1-r)\le 2/(k-2)\). Hence \(h_r''(p)\ge
2(2k)^{-2/(k-2)}\ge1/18\), using \(2k\le6^{k-2}\). For each
\(p\in[0,1]\), the half of the segment from \(u\) to \(p\) nearest \(u\)
lies in this curvature interval. The Taylor integral on that half has
length \(|p-u|/2\) and weight at least \(|p-u|/2\). The rest of the
segment has nonnegative curvature. Therefore
\[
 h_r(p)\ge h_r(u)+h_r'(u)(p-u)+(p-u)^2/72 .
\]
Summing and using \(\sum_i p_i=1\) proves (V1).

## V2. Hermite remainder at the voting scale

For any list of \(k\) Gaussian indicators, each of measure \(u\),
let \(W_m\) be the sum of their squared degree-\(m\) Hermite norms.
The list may repeat the same halfspace. Take \(s=1/\log k\), so
\(r\le s<1\). By the positive-correlation case of Borell's inequality
(Mossel–Neeman, Theorem 1.1 of the same source),
\[
 \sum_{m\ge2}r^m W_m
 \le (r/s)^2\sum_{m\ge0}s^m W_m
 \le (r/s)^2 k J_s(u,u).
\]
Put \(z=\Phi^{-1}(1-u)>0\) and \(Q=1-\Phi\). The bivariate density
at correlation \(s\) is at most \((1+s)/\sqrt{1-s^2}\) times the
product of the two \(N(0,1+s)\) densities, by \(2xy\le x^2+y^2\).
Integrating over the positive quadrant starting at \(z\) gives
\[
 J_s(u,u)\le \frac{1+s}{\sqrt{1-s^2}}
                Q(z/\sqrt{1+s})^2 .
\]
Integration by parts gives the Mills bounds
\(\phi(x)x/(1+x^2)\le Q(x)\le\phi(x)/x\) for \(x>0\). Thus
\[
 \frac{Q(z/\sqrt{1+s})}{Q(z)}
 \le \sqrt{1+s}(1+z^{-2})\exp\bigl(sz^2/(2(1+s))\bigr)
 \le C_0
\]
with an absolute \(C_0\), since \(z\ge\Phi^{-1}(2/3)\),
\(z^2\le2\log k\), and \(s=1/\log k\).
As \(s\le1/\log 3<1\) we obtain, for an absolute \(C\),
\[
 \sum_{m\ge2}r^m W_m\le C r^2(\log k)^2/k .
 \tag{V4}
\]

Put \(z=\Phi^{-1}(1-u)\), \(\nu=k\phi(z)\),
\(\mu=\mathbb E\max_{i\le k}Z_i\) for independent standard normals,
and \(F_k=\mu^2/(k-1)\).
The degree-one mass of the simplex is \(F_k\). To see the normalization,
take scores \(Z_i-\bar Z\) in \(\mathbf 1^\perp\). The centroid of the
\(i\)-th cell is \(b(e_i-\mathbf 1/k)\) by symmetry and the sum of the
winning-score integrals is \(b(k-1)=\mu\). The halfspace list has
degree-one mass \(k\phi(z)^2=\nu^2/k\). Both lists have degree-zero
mass \(1/k\). Thus the exact OU expansion and (V4) give
\[
 \left|S_{-r}(T_k)-k h_r(u)
        -r\left(\frac{\nu^2}k-\frac{\mu^2}{k-1}\right)\right|
 \le 2 C r^2(\log k)^2/k .
 \tag{V5}
\]

## V3. The first-moment difference between the two benchmarks

Write \(Q=1-\Phi\). For each fixed real \(t\), as \(z\to\infty\),
\[
 \frac{Q(z+t/z)}{Q(z)}\longrightarrow e^{-t}.
\]
For example, substituting \(x=z+v/z\) in the two normal tail
integrals gives the ratio of the integrals of
\(\exp(-v-v^2/(2z^2))\) starting at \(t\) and \(0\), respectively.
For \(t\ge0\) the ratio is bounded by \(e^{-t}\), by shifting the
normal density in the tail. Because \(kQ(z)=1\), this limit gives
\(z(\max_i Z_i-z)\) converging in law to the standard Gumbel
distribution, with CDF \(\exp(-e^{-t})\). Its positive-tail
probability is bounded by \(e^{-t}\) for \(t\ge0\). For
\(0\le t\le z^2\), shifting the tail in the other direction gives
\(kQ(z-t/z)\ge \exp(t-t^2/(2z^2))\ge e^{t/2}\). The negative-tail
probability is then at most \(\exp(-e^{t/2})\). The remaining
negative-tail integral is bounded by
\[
 \int_{z^2}^\infty
 \Pr\bigl(z(\max_i Z_i-z)<-t\bigr)\,dt
 =z\int_0^\infty\Phi(-v)^k\,dv
 \le z\,2^{1-k}/\sqrt{2\pi}\longrightarrow0 .
\]
Here \(z^2\le2\log k\) by the Chernoff bound. The tail-integral formula
and dominated convergence now give \(z(\mu-z)\to\gamma_{\rm E}\).
Likewise \(\nu=\mathbb E[Z\mid Z>z]\), so that
\[
 z(\nu-z)=\int_0^\infty Q(z+t/z)/Q(z)\,dt\longrightarrow1 .
\]
It follows that
\[
 \nu^2-\frac{k}{k-1}\mu^2
 =(\nu-\mu)(\nu+\mu)-\frac{\mu^2}{k-1}
 \longrightarrow 2(1-\gamma_{\rm E}).
\]
Since \(r(\log k)^2\le(\log k)^2/(k-1)\to0\), (V5) proves (V2).
In addition \(z^2/(2\log k)\to1\), from the same tail integral
\(Q(z)\sim\phi(z)/z\). Substitution in (V4) for \(T_k\) proves
the claimed expansion of \(1/k-S_{-r}(T_k)\).

## Exact results and the remaining comparison

[Heilman, *Sharp Hardness for MAX-3-CUT and Quantum MAX-CUT*,
arXiv:2608.00333v1, Theorem 7.1](https://arxiv.org/pdf/2608.00333v1)
proves exact simplex minimization for \(k=3\) throughout
\([-1/2,0)\), in every \(d\ge2\), with no mass premise. Its
negative proof proceeds through the paired minimum with only
corresponding masses equal, after reflecting the second partition.
Section 1.4 identifies the spherical rearrangement used in the
planar case as an obstacle to this route for \(k\ge4\).

The value inequality in Theorem 7.1 extends to fractional
assignments in the same dimension. For a sequence of increasing
finite Gaussian measurable subdivisions generating the Borel
sigma-algebra, split each atom into \(k\) sets of measures
\(\int_{\rm atom} f_i\,d\gamma_d\), using atomlessness. The resulting
partition indicator tuples converge weakly to \(f\) in \(L^2\),
since integrals against every bounded test constant on a fixed
subdivision coincide from that index onward and those tests are dense.
For \(r<1\), \(T_{-r}\) is compact on Gaussian \(L^2\) in fixed finite
dimension, as its degree-\(m\) eigenvalues are \((-r)^m\) with finite
multiplicities. Thus the actual stability converges. The approximation
gives the value inequality for fractional assignments. It does not
transfer equality classification by itself.

The large-alphabet target requires eliminating the loss in (V3).
The scalar Borell comparison alone leaves the positive gap (V2) even
on the simplex. Control of the joint geometry at the \(r/k\) scale,
or a sharper comparison that vanishes at \(T_k\), is required. The [exact small-negative interval](NEGATIVE_SMALL_CORRELATION.md)
for each fixed \(k\) has a threshold that can depend on \(k\). A bound
at least \(1/(k-1)\) for all sufficiently large \(k\) remains open here. The estimates and source applications above
have been checked internally.
