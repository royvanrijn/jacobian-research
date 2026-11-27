# The negative voting interval for large alphabets

Let \(k\ge3\), \(d\ge k-1\), and \(0<r\le1/(k-1)\). For a measurable
assignment \(f:\mathbb R^d\to\Delta_k\), set \(p_i=\int f_i\,d\gamma_d\) and
write
\[
 S_{-r}(f)=\sum_i\mathbb E[f_i(X)f_i(Y)],
 \qquad \operatorname{Cov}(X,Y)=-rI_d.
\]
## Exact Gaussian simplex minimization

There is an integer \(M\ge3\) such that for every \(k\ge M\),
every \(d\ge k-1\), every \(0<r\le1/(k-1)\), and every measurable
\(f:\mathbb R^d\to\Delta_k\),
\[
 S_{-r}(f)\ge S_{-r}(1_{T_k}),
 \tag{VG1}
\]
where \(T_k\) is the regular central simplex extended cylindrically
to dimension \(d\). The cell masses are unrestricted. Equality
holds exactly for a regular central simplex assignment, up to
orthogonal maps, labels and null sets.

To prove this, consider an assignment with
\(S_{-r}(f)\le S_{-r}(1_{T_k})\).
[Exponential localization (VL1)](NEGATIVE_VOTING_EXPONENTIAL_LOCALIZATION.md)
bounds its summed label \(L^1\) distance to a simplex by
\(C\exp[-c/(r^2 k(\log k)^2)]\).
Since \(r\le1/(k-1)\), this is at most
\(C\exp[-c' k/(\log k)^2]\), uniformly in \(d,r\).
[The polynomial local comparison (PN1)](NEGATIVE_VOTING_POLYNOMIAL_NEIGHBORHOOD.md)
holds when the distance of a pair to a common simplex is at
most \(k^{-A}\), with one fixed \(A\). Apply it to \((f,f)\).
Its pair-distance bound is twice the single-assignment bound,
and \(2C\exp[-c' k/(\log k)^2]\le k^{-A}\) for all
sufficiently large \(k\). This gives the inequality and the
equality classification in (VG1).

The localization uses the prescribed-mass covariance estimate,
which derives its flux identities directly and uses the published
Gaussian multi-bubble theorem.
The second-degree estimate requires a Gaussian maximum comparison
after a common Gaussian shift, in the scope of
[Mulgund's Theorem 2.1](https://arxiv.org/pdf/2607.14087v2)
or [Friedland–Guédon–Souli, Corollary 1.5](https://arxiv.org/pdf/2609.38904v1).
The local step uses the signed weighted Hessian throughout the
voting interval. The proof has internal mathematical audit.
The value of \(M\) is qualitative.

## Discrete low-influence voting

For \(k\ge M\) and any \(\epsilon>0\) there are a low-degree cutoff
\(D\) and an influence threshold \(\tau>0\), depending on
\(k,\epsilon\), such that every \(f:[k]^N\to\Delta_k\)
with all coordinate influences through degree \(D\) at most
\(\tau\) satisfies
\[
 S_\rho(f)\ge
       \lim_{m\to\infty} S_\rho(\mathrm{PLUR}_{m,k})-\epsilon
 \qquad\left(-\frac1{k-1}\le\rho<0\right).
 \tag{VG2}
\]
This is the negative-correlation inequality of Plurality Is Stablest,
with the same low-degree thresholds for the entire interval.
There is no balance premise.

For a fixed \(\rho\), [Isaksson–Mossel, arXiv:0903.3362v3,
Theorem 7.1 and Lemma 7.3](https://arxiv.org/pdf/0903.3362v3)
give, respectively, Gaussian replacement in dimension \(N(k-1)\)
with arbitrarily small error in stability, and the simplex
value as the plurality limit. Their range includes the endpoint
\(-1/(k-1)\). Apply (VG1) in dimension \(N(k-1)\) to the
replacement. To take thresholds uniform in \(\rho\), use a
finite sufficiently fine net in the closed interval
\([-1/(k-1),0]\). The zero value is bounded below by \(1/k\)
by Cauchy–Schwarz. Both the discrete stability of any \(f\)
and the Gaussian simplex stability are Lipschitz uniformly
on this interval. Their power-series coefficients are
nonnegative with sum at most 1, and
\(\sup_{j\ge1}j(k-1)^{-(j-1)}<\infty\). Take the maximum
of the degree cutoffs and the minimum of the influence
thresholds over the net, allocating part of the error to
the Lipschitz differences.

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
\[
 \sum_i(p_i-u)^2\le \bigl(2(1-\gamma_{\rm E})+o(1)\bigr)r/k .
 \tag{V3a}
\]
In particular every such assignment has all cell probabilities
between \((1-\sqrt{2(1-\gamma_{\rm E})}-o(1))/k\) and
\((1+\sqrt{2(1-\gamma_{\rm E})}+o(1))/k\).
The lower endpoint is positive for all sufficiently large \(k\)
since \(\gamma_{\rm E}>1/2\).
The improvement
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

For the sharper asymptotic coefficient in (V3a), take
\(\varepsilon=k^{-2}\). On \([\varepsilon,1-\varepsilon]\) the
same derivative estimate gives
\(h_r''\ge2\exp(-4\log k/(k-2))\). The discarded part of the
Taylor integral from \(u\) to \(p\) has length and Taylor weight
at most \(\varepsilon\). If any part is discarded, then
\(|p-u|\ge u-\varepsilon\). The resulting tangent inequality is
uniform for \(p\in[0,1]\) with quadratic coefficient at least
\[
 e^{-4\log k/(k-2)}
 \left(1-\varepsilon^2/(u-\varepsilon)^2\right)=1-o(1).
\]
Summing and applying (V2) proves (V3a). The Euler constant satisfies
\(\gamma_{\rm E}>\sum_{j=1}^6 1/j-\log 7>1/2\).
The first bound follows from
\(\gamma_{\rm E}=\lim_m(\sum_{j=1}^m 1/j-\log m)\) and the
strict increase of \(\sum_{j=1}^m 1/j-\log(m+1)\).
For the second, \(\sum_{j=1}^6 1/j=49/20\) and
\(e^{39/20}>7\) by its Taylor polynomial of degree six.

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
The same tails give \(\operatorname{Var} M=O(1/z^2)\).
For the second moment of \(z(M-z)\), the bounded right and
left tails up to the same cutoffs are integrated against \(2t\,dt\).
The omitted left tail, after \(t=z^2+zv\), is at most
\(2z\,2^{1-k}\int_0^\infty (z^2+zv)\Phi(-v)\,dv=o(1)\).
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

## V4. First moments of possible improving assignments

Let \(B_i=\int x f_i(x)\,d\gamma_d(x)\), \(W_1(f)=\sum_i\|B_i\|^2\)
and \(\delta^2=\sum_i(p_i-u)^2\).
A competitor satisfying \(S_{-r}(f)\le S_{-r}(T_k)\) obeys
\[
 W_1(f)\ \ge\ F_k+\delta^2/r-C_1 r/k
 \tag{V6}
\]
for all sufficiently large \(k\), with an absolute \(C_1\).
Moreover, with \(I(p)=\phi(\Phi^{-1}(p))\) extended by zero at the
endpoints, it satisfies
\[
 0\le \sum_i\bigl(I(p_i)^2-\|B_i\|^2\bigr)
 \le \bigl(2(1-\gamma_{\rm E})+o(1)\bigr)/k
 \tag{V7}
\]
uniformly over the whole voting interval and ambient dimensions.

Indeed, (V1)–(V2) give \(\delta^2=O(r/k)\), so \(p_i\le b/k\) for
an absolute \(b>1\). The positive functional Borell inequality and
\(J_s(a,b)\ge ab\) for positive \(s\) give
\(\sum_m s^m W_m(f)\le \sum_i J_s(p_i,p_i)\le kJ_s(b/k,b/k)\).
The bivariate density comparison above now applies at
\(z_b=\Phi^{-1}(1-b/k)\), with \(z_b\ge\Phi^{-1}(2/3)\) for
\(k\ge3b\) and \(z_b^2\le2\log k\). It bounds the last term by
\(C_2/k\). The same Hermite argument as in (V4) bounds the absolute
higher-degree remainder of \(f\) by \(C_2 r^2(\log k)^2/k\).
Since \(W_0(f)=1/k+\delta^2\), the first comparison gives
\[
 \Delta:=F_k-W_1(f)+\delta^2/r\le C_3 r(\log k)^2/k .
\]
We next sharpen this estimate using the sign of the second degree.

For a cell with \(p_i>0\), let
\(A_i=\int (xx^T-I_d)f_i(x)\,d\gamma_d(x)\) and let \(v_i\) be
a unit vector in the direction of \(B_i\), arbitrary when \(B_i=0\).
By Cauchy–Schwarz in the measure \(f_i\,d\gamma_d\),
\[
 p_i v_i^TA_i v_i\ge\|B_i\|^2-p_i^2 .
\]
For zero-mass cells both sides are zero. If \(P_2=\sum_i p_i^2\),
then \(W_2(f)=\tfrac12\sum_i\|A_i\|_F^2\) gives
\[
 W_2(f)\ \ge\ L(W_1(f),P_2),\qquad
 L(w,p)=\frac{(w-p)_+^2}{2p}.
\]
Indeed \(\sum_i p_i v_i^TA_i v_i\le
\sqrt{2P_2W_2(f)}\). The preliminary bound on \(\Delta\) gives
\(W_1(f)\ge F_k-o(1/k)\), while \(P_2=1/k+O(r/k)\).

For the simplex, let \(M=\max_i Z_i\) and \(\mu=\mathbb E M\).
Its winner's axial second moment is
\(k(\mathbb E M^2-1/k)/(k-1)\), because the centered scores
are independent of the common sample mean of variance \(1/k\).
The other active second moments are all equal by the permutations
fixing that label. Conicity gives conditional trace \(k-1\) in
the active space. Consequently
\[
 W_2(T_k)=\frac{k(\mathbb E M^2-1)^2}{2(k-1)(k-2)} .
\]
This follows by summing squared centered moment matrices with
axial eigenvalue \(\lambda-1\), \(k-2\) other eigenvalues
\(-(\lambda-1)/(k-2)\), and mass \(1/k\) per label.
Using \(\mu^2=O(\log k)\) and \(\operatorname{Var} M=O(1/\log k)\)
shows
\[
 W_2(T_k)-L(F_k,1/k)=O(1/k).
\]
In fact \(L(F_k,1/k)=k(\mu^2-1+1/k)^2/(2(k-1)^2)\) for large \(k\).
The difference of the prefactors contributes \(O((\log k)^2/k^2)\),
and the variance shift in the squared terms contributes \(O(1/k)\).
When \(w\le F_k\) and \(p\ge1/k\), the first derivative of \(L\)
in \(w\) is at most \(O(\log k)\), and the absolute derivative
in \(p\) is at most \(O((\log k)^2)\). Using monotonicity in \(w\)
for \(W_1(f)>F_k\) gives
\[
 W_2(T_k)-W_2(f)\le
 C_4\left(\frac1k+(\log k)\Delta_+
                +(\log k)^2\delta^2\right).
\]
Here \((F_k-W_1(f))_+\le\Delta_+\).

The bound for Hermite degrees at least three is
\(O(r^3(\log k)^3/k)\) for both assignments, by the same
positive smoothing with \((r/s)^3\). Since the second Hermite
coefficient in \(S_{-r}\) is \(+r^2\),
\[
 \Delta\le r\bigl(W_2(T_k)-W_2(f)\bigr)
                 +O(r^2(\log k)^3/k).
\]
With \(\delta^2=O(r/k)\) and \(r\le1/(k-1)\), the right side
is at most \(C_5 r/k + C_4 r\log(k)\,\Delta_+\).
If \(\Delta>0\), absorb the last term for all sufficiently
large \(k\). If \(\Delta\le0\), the same upper bound
\(\Delta\le C_1 r/k\) holds. This proves (V6).

For every cell, \(\|B_i\|\le I(p_i)\). If \(B_i\ne0\), compare \(f_i\)
with the indicator of the upper halfspace of measure \(p_i\) in the
direction of \(B_i\). Against the coordinate in that direction minus
the threshold, their difference integrates to a nonpositive value,
and the threshold times their mass difference is zero. This proves
the first inequality in (V7).

On \((0,1/2]\), the second derivative of \(I(p)^2\) is
\(2(\Phi^{-1}(p)^2-1)\le4\log(1/p)\). For \(u\le p\le b/k\le1/2\),
the Taylor remainder at \(u\) is at most
\(2\log(k)(p-u)^2\). For \(0\le p<u\), write
\(\log(1/t)=\log k+\log(u/t)\) in the Taylor integral with weight
\(t-p\). The inequality \((t-p)\log(u/t)\le u-t\) gives an upper
bound of \(2(\log k+1)(p-u)^2\), also at \(p=0\) by continuity.
Summing cancels the first derivative term. Therefore
\[
 \sum_i I(p_i)^2\le k I(u)^2+2(\log k+1)\delta^2.
\]
Insert (V6), \(\delta^2=O(r/k)\), and the limit for
\(k(k I(u)^2-F_k)\) proved above. Both
\((\log k+1)\delta^2\) and \(r(\log k)^2/k\) are \(o(1/k)\)
on the voting interval. This proves (V7).

## Orthoplex competitors and moment localization

For even \(k=2m\ge4\), use \(m\) orthogonal Gaussian coordinates and
assign the label \(+i\) or \(-i\) according to the sign and index of
the largest absolute coordinate. Denote this balanced orthoplex
partition by \(O_k\) and extend it to every \(d\ge k-1\).
With \(M=\max_{i\le k}Z_i\) and \(U=\max_{i\le m}|Z_i|\) for independent
standard normal coordinates in the two definitions, one has
\[
 W_1(T_k)-W_1(O_k)
   =(2+o(1))\frac{\log k}{k^2},\qquad
 \inf_{T}\sum_i\|1_{O_{k,i}}-1_{T_i}\|_{L^1(\gamma_d)}
   \ge 1-o(1),
 \tag{V8}
\]
where the infimum is over all orientations and labels of the regular
simplex. Despite this small first-moment deficit, the orthoplex loses
on the actual voting interval, uniformly in \(0<r\le1/(k-1)\):
\[
 S_{-r}(O_k)-S_{-r}(T_k)
   =(2+o(1))\frac{r\log k}{k^2}.
 \tag{V9}
\]

For \(t\ge0\) and \(q=Q(t)\), the distribution functions satisfy
\[
 0\le \Pr(M\le t)-\Pr(U\le t)
 =(1-q)^k-(1-2q)^m
 \le m q^2(1-q)^{k-2}.
\]
This is the difference of the \(m\)-th powers of
\((1-q)^2\) and \(1-2q\). Put \(z=Q^{-1}(1/k)\).
For \(t\ge z/2\), the Mills inequality gives \(\phi(t)\ge tQ(t)\).
Changing variables to \(q\), the integral of the last bound
against \(t^j\,dt\) for \(j=0,1\) is bounded above by
\[
 m(z/2)^{j-1}\int_0^1 q(1-q)^{k-2}\,dq
 =\frac{(z/2)^{j-1}}{2(k-1)}.
\]
For \(0\le t\le z/2\), the difference is bounded by
\(\Pr(M\le z/2)\le\exp(-kQ(z/2))\).
Translating the Gaussian tail integral by \(z/2\) gives
\(kQ(z/2)\ge e^{3z^2/8}\).
For \(t\ge0\), \(\Pr(M<-t)=\Phi(-t)^k\le2^{1-k}\Phi(-t)\),
which controls the negative-tail contributions to the first two
moments. Consequently
\[
 0\le \mathbb E U-\mathbb E M=O(1/(kz)),\qquad
 |\mathbb E U^2-\mathbb E M^2|=O(1/k).
 \tag{V10}
\]
The conditional centroid of an orthoplex cell is along its winning
direction and has length \(\mathbb E U\), by symmetry. Thus
\(W_1(O_k)=(\mathbb E U)^2/k\).
Writing \(\mu=\mathbb E M\) and \(\omega=\mathbb E U\), the difference
in first energies is
\[
 \frac{\mu^2}{k-1}-\frac{\omega^2}{k}
 =\frac{\mu^2}{k(k-1)}-\frac{\omega^2-\mu^2}{k}.
\]
Here \(\mu/z\to1\), \(z^2\sim2\log k\), and
\(\omega^2-\mu^2=O(1/k)\) by (V10). This proves the first part of (V8).

To prove the distance bound, pair the labels \(+i,-i\) of \(O_k\).
For each pair, at least one of its two unit score directions has
inner product at most \(\sqrt{k/(2(k-1))}\) with the corresponding
unit simplex direction. Indeed, if both projections onto the
corresponding signed orthoplex axis exceed this number, the difference
of the two simplex directions has length exceeding
\(\sqrt{2k/(k-1)}\), its known length. For large \(k\) the bound
is at most \(3/4\). The maxima of the unit scores for \(O_k\) and for
a simplex have distributions \(U\) and \(\sqrt{k/(k-1)}(M-\bar Z)\).
Both exceed \(0.95z\) with probability tending to one, uniformly
over the orientations, by the bounds above and \(\bar Z/z\to0\)
in probability. Agreement on one of these \(k/2\) labels then requires
two standard Gaussian scores of correlation at most \(3/4\) to
exceed \(0.95z\). The Gaussian tail bound for their sum gives the
uniform union bound
\[
 k\exp\big(-(0.95z)^2/(1+3/4)\big)=o(1).
\]
The other \(k/2\) labels carry total orthoplex mass \(1/2\).
Agreement is at most \(1/2+o(1)\), and twice the disagreement is the
sum of the indicator \(L^1\) distances in (V8).

For the actual comparison, in an orthoplex cell the signed
permutations fixing its label force a scalar transverse second
moment and zero mixed entries. The axial second moment is
\(\mathbb E U^2\), and the trace is \(m\) by conicity. Applying
the same axial second-moment formula as above, this time with
\(m-1\) transverse coordinates, gives
\[
 W_2(O_k)=\frac{(\mathbb E U^2-1)^2}{2(k-2)} .
\]
Both second moments \(\mathbb E U^2,\mathbb E M^2\) are \(O(\log k)\)
by integrating the Gaussian union bound for the largest absolute
coordinate. By (V10), the displayed energies differ by
\(O((\log k)^2/k^2)\). The proof of (V4), with the factor
\((r/s)^3\) for degrees at least three, bounds the remaining
absolute energies of both partitions by \(O(r^3(\log k)^3/k)\).
Thus the first-degree difference in (V8) gives (V9). Indeed the
ratios of the degree-two and higher error bounds to
\(r\log k/k^2\) are \(O(r\log k)\) and
\(O(r^2 k(\log k)^2)\), respectively, and both tend to zero
uniformly on the voting interval.

The sharper competitor criterion (V6) excludes this family throughout
the voting interval for large enough even \(k\), since \(r/k=O(k^{-2})\).
In particular a localization criterion requiring only
\(W_1(T_k)-W_1(f)=o(1/k)\) cannot force convergence to the simplex
in label measure for growing \(k\). The exact sign in (V9)
excludes \(O_k\) as an improving competitor.

## Centroid-length spread and the second Hermite deficit

Use the full negative voting interval and write
\(\Delta=F_k-W_1(f)+\delta^2/r\).
For a possible improver, meaning \(S_{-r}(f)\le S_{-r}(T_k)\), all
\(p_i>0\) for sufficiently large \(k\) by (V3a). Set
\[
 m_i=\|B_i\|/p_i,\qquad m=\sum_i p_i m_i,\qquad
 H=\sum_i p_i(m_i-m)^2 .
\]
Let \(v_i\) be the unit direction of \(B_i\), taking \(v_i\) in the
span of the nonzero \(B_j\) if \(B_i=0\). Put
\[
 G=\left\|\frac1k\sum_i v_i v_i^T\right\|_F^2,\qquad
 L=\log k,\qquad \Lambda=\sqrt{\log k\,\log(k/r)} .
\]
The centroid span is nonzero by (V6). Since \(\sum_i B_i=0\), its
dimension is at most \(k-1\), so \(G\ge 1/(k-1)\).
There are absolute \(c,C>0\) such that
\[
 \Delta+c\,\frac{rL^2}{k}\left(G-\frac1{k-1}\right)
 \le C\left(\frac{r\Lambda}{k}H+\frac{r^2 L^3}{k}\right).
 \tag{V11}
\]
In particular, if the conditional centroid lengths \(m_i\) are
equal, with no prescribed cell masses, then the actual noise
stability is at least
\[
 S_{-r}(T_k)-C r^3(\log k)^3/k .
 \tag{V12}
\]
In that class, every possible improver has
\(\delta^2\le C r^3(\log k)^3/k\) and
\(G-1/(k-1)\le C r\log k\). Equation (V11) applies to arbitrary
centroid lengths and quantifies the extra term introduced by
their spread.

**The scalar comparison.** Let \(X\) be standard Gaussian and
conditionally draw the label \(I\) with probabilities \(f_i(X)\).
Write \(Y=v_I\cdot X\) and let \(Y_*\) be the maximum of
the unit regular-simplex scores. Take \(\alpha=(k-1)/k\).
Add independent common Gaussian noise \(N\) of variance
\(1/(k-1)\) to the scores \(v_i\cdot X\). If \(G_0\) is their
Gram matrix, multiplying the shifted vector by \(\sqrt\alpha\)
gives a unit-variance Gaussian vector with covariance
\(\alpha G_0+\mathbf1\mathbf1^T/k\).
Its maximum is stochastically dominated by the maximum of \(k\)
independent standard normals by
[Mulgund, *Stochastic Domination of Gaussian Maxima: A Resolution of the Weak Simplex Conjecture*, arXiv:2607.14087v2, Theorem 2.1](https://arxiv.org/pdf/2607.14087v2).
The theorem permits singular covariance and requires
\(R-\mathbf1\mathbf1^T/k\succeq0\) for the compared
correlation matrix \(R\), as holds here.
The same comparison follows from
[Friedland–Guédon–Souli, *On partition functions of Gaussian random variables*, arXiv:2609.38904v1, Corollary 1.5 and Remark 1.6](https://arxiv.org/pdf/2609.38904v1)
with the common log-concave profile a half-line indicator.

For the regular unit simplex scores this common-noise
construction gives independent coordinates after scaling.
Since \(Y\le\max_i v_i\cdot X\) pointwise, it follows that
\[
 \widetilde Y=Y+N\ \le_{\rm st}\
   \widetilde Y_*=Y_*+N',\qquad
 \widetilde Y_*\ \stackrel{\rm law}=\ \alpha^{-1/2} M_k ,
 \tag{V12a}
\]
where \(N'\) is also independent of the regular scores and has
the same variance. The common shift leaves both original means
and their second-moment difference unchanged. Denote
\(a=\mathbb E Y_*=\sqrt{kF_k}\), \(B_*=\mathbb E Y_*^2\),
and \(B=\mathbb E Y^2\). Then \(m=\mathbb EY\le a\).
Also \(m\ge c\sqrt{\log k}\), because \(m_i\le C\sqrt{\log k}\)
by the halfspace centroid bound and (V3a), and
\[
 k W_1(f)=k\sum_i p_i^2 m_i^2
       \le C\sum_i p_i m_i^2\le C\sqrt{\log k}\,m .
\]
The lower bound (V6) gives \(kW_1(f)\ge (2+o(1))\log k\)
as a lower estimate.

Put \(q_i=p_i-1/k\). Expanding \(\sum_i p_i m_i^2=m^2+H\)
and \(kW_1=k\sum_i p_i^2m_i^2\) gives
\[
 a^2-m^2=k(F_k-W_1)+H+k m^2\delta^2+
                  k\sum_i q_i p_i(m_i^2-m^2).
\]
As \(|m_i^2-m^2|\le C\sqrt L\,|m_i-m|\) and \(p_i\le C/k\),
Cauchy–Schwarz bounds the last term in absolute value by
\(C\delta\sqrt{k L H}\). Young's inequality and \(r\le1/(k-1)\)
therefore give, for sufficiently large \(k\),
\[
 0\le a^2-m^2\le k\Delta+2H-\frac{k}{2r}\delta^2 .
 \tag{V13}
\]
Couple the smoothed variables in (V12a) in nondecreasing
quantile order. For \(R>0\), put \(\psi_R(t)=\min(t,R)^2\).
Its slope is bounded above by \(2R\) on the entire real line,
so for \(x\le y\) the difference
\(\psi_R(y)-\psi_R(x)\) is at most \(2R(y-x)\). This one-sided
bound suffices even for negative scores.
The difference between \(\widetilde Y_*^2\) and
\(\psi_R(\widetilde Y_*)\) has expectation at most
\((2k/\alpha)e^{-\alpha R^2/2}\) by the Gaussian union bound,
while the corresponding difference for \(\widetilde Y\)
is nonnegative. Thus for any \(R>0\),
\[
 B_*-B \le 2R(a-m)+(2k/\alpha)e^{-\alpha R^2/2}.
 \tag{V14}
\]
Taking \(R=\sqrt{(2/\alpha)\log(k/r^2)}\) and using
\(a+m\ge a\) in (V13) bounds this by
\(C(R/a)(k\Delta+2H)+Cr^2\).

**The matrix constraint.** Let \(V_i=v_i v_i^T\) and
\(\bar V=k^{-1}\sum_i V_i\). The centered second moments \(A_i\)
defined above satisfy \(\sum_i A_i=0\), as the assignments sum to
one. Therefore
\[
 B-1=\sum_i\langle A_i,V_i-\bar V\rangle_F,\qquad
 \sum_i\|V_i-\bar V\|_F^2=k(1-G).
\]
Since \(B\ge m^2\ge cL\), one has \(G<1\) for sufficiently large
\(k\), and Cauchy–Schwarz gives
\[
 W_2(f)\ge \frac{(B-1)^2}{2k(1-G)} .
 \tag{V15}
\]
The same formula is exact at the simplex with \(G=1/(k-1)\).
Indeed its active second-moment matrices are proportional to
\(V_i-\bar V\) by the axial calculation in (V4).

Let \(g_*=1/(k-1)\). The bound
\((1-G)^{-1}\ge (1-g_*)^{-1}+G-g_*\) follows from \(G\ge g_*\).
If \(B\le B_*\), the difference of the squares
\((B_*-1)^2-(B-1)^2\) is at most \(C L (B_*-B)\).
If \(B>B_*\), the difference is nonpositive.
Use \(B_*=O(L)\), \(B\ge cL\), and the nonnegative upper bound
in (V14) to obtain
\[
 W_2(T_k)-W_2(f)\le
 C\left(\Lambda\left(\Delta+2H/k\right)+r^2 L/k\right)
       -c\,\frac{L^2}{k}(G-g_*).
\]
Here \(LR/a=O(\Lambda)\).
The Hermite remainder of degrees at least three from (V4),
together with nonpositive actual deficit, gives
\(\Delta\le r(W_2(T_k)-W_2(f))+C r^2 L^3/k\).
The coefficient \(r\Lambda\) tends to zero uniformly on the
voting interval. This follows by differentiating
\(r^2\log(k/r)\), which is increasing on that range.
Absorbing the \(\Delta\) term proves (V11).

For \(H=0\), (V13) gives \(\delta^2/(2r)\le \Delta\).
In particular \(\Delta\ge0\), so (V11) proves the asserted mass
and matrix bounds. To obtain (V12), consider only a possible improver.
The second-moment bound above and \(\Delta=O(r^2L^3/k)\)
show that \(r\Delta+r^2(W_2(f)-W_2(T_k))\) is bounded below
by \(-C r^3 L^3/k\), using \(\Delta\ge0\) and small \(r\Lambda\).
The higher-degree remainder has the same absolute bound.
An assignment with positive deficit already satisfies (V12).

## Dimension-uniform localization on the voting scale

The [prescribed-mass covariance estimate](NEGATIVE_VOTING_COVARIANCE.md)
and the [mass estimate for possible improvers](NEGATIVE_VOTING_MINIMIZERS.md)
strengthen the global conclusions of (V3) and (V6). For all sufficiently
large \(k\), for \(0<r\le1/(k-1)\), and uniformly in \(d\ge k-1\),
\[
 \inf_{f:\mathbb R^d\to\Delta_k} S_{-r}(f)
   \ge S_{-r}(T_k)-C\,r^3(\log k)^3/k .
 \tag{V16}
\]
For every possible improver \(f\), there is a cylindrical regular
simplex \(T\) such that
\[
 \sum_i\|f_i-1_{T_i}\|_1
     \le C\,r\sqrt{k}\,(\log k)^2.
 \tag{V17}
\]
In particular this distance tends to zero throughout the voting
interval as \(k\to\infty\). The constants here are absolute.

Let \(B_f\) be the matrix with rows \(B_i\), and write
\(Q=B_f B_f^T/W_1(f)\) and \(Q_*=(I-\mathbf1\mathbf1^T/k)/(k-1)\).
The mass estimate gives \(\max_i |kp_i-1|=o(1)\) uniformly on
the voting interval for possible improvers. The halfspace centroid
bound and (V6) give \(\max_i Q_{ii}\le C_0/k\) for an absolute
\(C_0\). The covariance estimate consequently gives
\[
 W_1(f)\le F_k+C_1(\log k)^2\delta^2
        -c_1\left(\|\operatorname{diag}(Q-Q_*)\|_2^2
                     +\|Q-Q_*\|_F^2/k\right).
\]
Since \(r(\log k)^2\) is small, this implies
\[
 \Delta\ge \frac{\delta^2}{2r}
        +c_1\left(\|\operatorname{diag}(Q-Q_*)\|_2^2
                     +\|Q-Q_*\|_F^2/k\right)\ge0 .
 \tag{V18}
\]
The centroid-length spread in (V11) satisfies
\(H\le C_2 k(\log k)\Delta\).
To see this, upper-bound the variance by the squared distance to
the constant \(\sqrt{kW_1(f)}\). Use
\(m_i=\sqrt{W_1(f)}\sqrt{Q_{ii}}/p_i\) and split the difference
into that of the square roots of the diagonal entries and that
of the reciprocal probabilities. Since \(p_i\asymp1/k\) and
\((\sqrt{Q_{ii}}-1/\sqrt{k})^2\le k(Q_{ii}-1/k)^2\), the two
terms give
\[
 H\le C_3 k(\log k)
        \left(\|\operatorname{diag}(Q-Q_*)\|_2^2+\delta^2\right).
\]
Applying (V18) proves the bound on \(H\).
Now \(r(\log k)\Lambda=O((\log k)^2/k)\) uniformly.
Absorb the spread term in (V11) and discard its nonnegative
matrix term to get \(\Delta\le C_4 r^2(\log k)^3/k\).
In particular (V18) yields
\(\delta^2=O(r^3(\log k)^3/k)\) and
\(\|Q-Q_*\|_F^2=O(r^2(\log k)^3)\).
The signed second-degree comparison used in (V11) bounds the
second-degree loss in the actual stability by
\(C_5 r^2(\log k)\Lambda\Delta+O(r^4\log(k)/k)\).
This is absorbed by \(r\Delta\) up to the
\(O(r^3(\log k)^3/k)\) Hermite remainder of degrees at least
three. Assignments of positive deficit already satisfy (V16).

For the label distance, choose a regular moment matrix \(B_*\) in
a matching polar gauge with \(B_f\). Such a choice exists even
for singular \(B_f\), since the active label space has dimension
\(k-1\) and \(d\ge k-1\). On that label space
\(Q_*^{1/2}\) is scalar, and
\[
 \|Q^{1/2}-Q_*^{1/2}\|_F^2\le(k-1)\|Q-Q_*\|_F^2.
\]
It follows that
\(\|B_f-B_*\|_F^2=O(r^2(\log k)^4)\),
using \(W_1(f)=O(\log k/k)\) and
\(|W_1(f)-F_k|\le\Delta\) by (V18).
The exact regret against the regular centroid scores is
\(F_k-\langle B_f,B_*\rangle_F\), which equals
\((F_k-W_1(f)+\|B_f-B_*\|_F^2)/2\).
It is therefore \(O(r^2(\log k)^4)\).

The Gaussian measure of the strip where the top two regular
scores differ by at most \(t\) is at most \(C_6 kt\).
Indeed, for \(k\) independent standard normals the probability
that the two largest differ by at most \(h\) is
\[
 k(k-1)\int_{\mathbb R}\phi(x)\Phi(x)^{k-2}
                       \int_x^{x+h}\phi(y)\,dy\,dx .
\]
On \(x\ge0\) this is bounded by \(h k(k-1)\int\phi^2\Phi^{k-2}
=h\,\mathbb E M\) using Gaussian integration by parts.
On \(x<0\) the bound is \(h\,k\phi(0)2^{1-k}\).
The mean subtraction in the regular score representation
does not change the gap, and its scale factor
\(\sqrt{F_k/(k-1)}\asymp \sqrt{\log k}/k\) proves the
strip bound. Off that strip a wrong label costs at least \(t\)
in regular-score regret. Minimizing the resulting bound
\(C_7 kt+O(r^2(\log k)^4/t)\) proves (V17).

The signed higher-degree comparison can use the error itself.
Let \(\mathcal D=\sum_i\|f_i-1_{T_i}\|_1\) for the polar
alignment just used, and write \(P_j\) for Gaussian Hermite
projection. Define the remainder operator
\(R_r=\sum_{j\ge3}(-r)^jP_j\). The difference of the two
higher-degree remainders in the objective is exactly
\(\sum_i\langle f_i-1_{T_i},R_r(f_i+1_{T_i})\rangle\).

Take \(s=1/\log k\), \(q=\log k\), and \(t=(q-1)^{-1/2}\)
for sufficiently large \(k\). Factoring out \(T_t\) from
\(R_r\), Gaussian hypercontractivity
from \(L^2\) to \(L^q\) gives, for \(v_i=f_i+1_{T_i}\),
\[
 \|R_r v_i\|_q
 \le (r/(t\sqrt{s}))^3\|T_{\sqrt{s}}v_i\|_2
 \le C_8\,r^3(\log k)^3/k.
\]
The first inequality follows degree by degree in \(L^2\), since
\(r/(t\sqrt s)\le1\). The second uses hypercontractivity
from \(L^{1+s}\) to \(L^2\) on the two summands and their
probabilities \(p_i,u\le C/k\).
Use
[Nelson's inequality as stated in Gross, Theorem 1.1](https://doi.org/10.1007/BF02392573).
Hölder on each label and concavity in their \(L^1\) errors give
\[
 \left|\sum_i\langle f_i-1_{T_i},R_r v_i\rangle\right|
 \le C_9 r^3 L^3 \mathcal D^{1-1/L}/k .
 \tag{V19}
\]
Here \(|f_i-1_{T_i}|\le1\) gives its \(L^{q/(q-1)}\) norm
at most \(\|f_i-1_{T_i}\|_1^{1-1/q}\), and summing those
powers costs at most \(k^{1/q}=e\).

Repeat the absorption argument using (V19) in place of the
absolute remainder. The discarded second-moment tail in (V14)
contributes \(O(r^3 L/k)\) to the bound on \(\Delta\).
Meanwhile the same polar-alignment and strip estimates show
\(\mathcal D^2\le C_{10} k^2 L\Delta\) from (V18).
It follows that
\[
 \Delta\le C_{11}\left(
       r^2 L^3 \mathcal D^{1-1/L}/k+r^3 L/k\right),\qquad
 \mathcal D^2\le C_{12}\left(
            k r^2 L^4\mathcal D^{1-1/L}+k r^3 L^2\right).
\]
Split according to which term dominates on the right of the
second inequality to obtain
\[
 \mathcal D \le C_{13}
      \left((k r^2 L^4)^{1/(1+1/L)}+\sqrt{k r^3}\,L\right)
   =O((\log k)^4/k)
 \tag{V20}
\]
uniformly for \(r\le1/(k-1)\). This uses
\(k^{1/(1+L)}\le e\), and then
\(\mathcal D^{1-1/L}\le C_{14} L^4/k\).

Substituting (V20) in the true deficit estimate gives the
uniform refinement of (V16)
\[
 \inf_f S_{-r}(f)\ge S_{-r}(T_k)-C_{15}\,r^3(\log k)^7/k^2 .
 \tag{V21}
\]
In this estimate the second-moment truncation contributes
\(O(r^4 L/k)\), which is absorbed for \(r\le1/(k-1)\).
Also (V18) now gives
\[
 \|Q-Q_*\|_F^2
     =O\left(r^2(\log k)^7/k+r^3\log k\right).
\]
In particular \((k-1)\|Q-Q_*\|_{\mathrm{op}}=o(1)\)
uniformly on the voting interval. Every possible improver
has full centroid rank \(k-1\) for sufficiently large \(k\),
and its normalized covariance has all active eigenvalues
\((1+o(1))/(k-1)\).

The polynomial-distance estimate here is strengthened in
[exponential localization](NEGATIVE_VOTING_EXPONENTIAL_LOCALIZATION.md).
The stronger version fits the local radius in (PN1).

## Earlier interval and large-alphabet comparisons

[Khot–Kindler–Mossel–O'Donnell, *Optimal Inapproximability Results
for MAX-CUT and Other 2-Variable CSPs?*, Proposition 11.4, final
version, pp. 31–32](https://www.stat.berkeley.edu/~mossel/publications/max_cut_final.pdf)
derive a large-alphabet lower bound for low-influence discrete
fractional functions at the endpoint \(-1/(k-1)\). It has leading
correction \(2\log(k)/k^2\) and error of order \(\log\log(k)/k^2\).
The analytic comparison (V1)–(V7) concerns the Gaussian functional
throughout the negative voting interval.

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

The three-cell result and (VG1) have different ranges of \(k\).
The exact result (VG1) proves the requested eventual Gaussian
inequality and its discrete consequence (VG2).
