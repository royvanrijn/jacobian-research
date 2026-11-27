# An explicit alphabet cutoff for negative voting

For every \(k\ge\lceil e^{3000}\rceil\), every finite \(d\ge k-1\),
every \(0<r\le 1/(k-1)\), and every measurable
\(f:\mathbb R^d\to\Delta_k\),
\[
 S_{-r}(f)=\sum_i\mathbb E[f_i(X)f_i(Y)]\ \ge\ S_{-r}(1_{T_k}),
 \qquad \operatorname{Cov}(X,Y)=-rI_d.
 \tag{EC1}
\]
Here \(T_k\) is the regular central simplex extended cylindrically.
Equality holds precisely for simplex assignments up to orthogonal
maps, labels and null sets. Cell masses are unrestricted. The same
alphabet cutoff applies to the discrete low-influence conclusion
[(VG2)](NEGATIVE_VOTING_INTERVAL.md#discrete-low-influence-voting).

## Preliminary bounds for a minimizer

Fix \(k,d,r\). The weak compactness argument leading to
[(M5)](NEGATIVE_VOTING_MINIMIZERS.md#existence-purity-and-pointwise-optimality)
gives a minimizer in the full class of fractional assignments.
All minimizers have objective at most the simplex value. Let
\(L=\log k\), \(u=1/k\), \(p_i=\int f_i\,d\gamma_d\),
\(\delta^2=\sum_i(p_i-u)^2\), \(z=\Phi^{-1}(1-1/k)\), and
\(W_j(f)=\sum_i\|P_jf_i\|_2^2\) for Gaussian chaos projection \(P_j\).
Write \(F_k=W_1(1_{T_k})\), and let \(M_k\) be the maximum of \(k\)
independent standard normals. By symmetry,
\(F_k=(\mathbb E M_k)^2/(k-1)\) and \(\ell=\sqrt{2F_k/(k-1)}\)
is the distance between distinct simplex cell moments.

The scalar comparison in
[(V1)](NEGATIVE_VOTING_INTERVAL.md#v1-uniform-full-interval-comparison-and-mass-control)
gives \(S_{-r}(f)\ge k h_r(u)+\delta^2/4\) for \(L\ge20\).
Indeed \(h_r''(v)\ge 2e^{-2\log(2k)/(k-2)}\ge1\) on
\([u/2,1-u/2]\), and the half of each Taylor segment nearest
\(u\) gives at least \((p_i-u)^2/4\).

Let \(Q_1(x)=1-\Phi(x)\), \(\mu=\mathbb EM_k\), and \(\nu=k\phi(z)\).
Mills and Chernoff bounds give \(\sqrt L\le z\le\sqrt{2L}\) for
\(L\ge20\), since \(kQ_1(\sqrt L)\ge e^{L/2}/(8\sqrt L)>1\).
Mills also gives \(\nu\le z+1/z\). Taking the expectation of
\(M_k\le z+\sum_i(Z_i-z)_+\) gives \(\mu\le\nu\).
For \(0\le t\le z\), the hazard bound
\(\phi(y)/Q_1(y)\ge y\) for \(y\ge0\) gives
\(kQ_1(z-t)\ge e^{zt/2}\). Thus
\[
 \int_0^z\Pr(M_k\le z-t)\,dt
 \le \int_0^\infty e^{-1-zt/2}\,dt=\frac{2}{ez}.
\]
The remaining part of \(\mathbb E(z-M_k)_+\) is at most
\(2^{1-k}/\sqrt{2\pi}\) by integrating
\(\Phi(-t)^k\le 2^{1-k}\Phi(-t)\) for \(t\ge0\).
It follows that \(\mu\ge z-2/z\). In particular
\(\mu\ge z/2\), \(L/(4k)\le F_k\le4L/k\), and
\[
 k\phi(z)^2-F_k\le(\nu^2-\mu^2)/k\le12/k.
\]

Set \(s=1/L\). Nelson's inequality and \(0\le f_i\le1\) give
\(\sum_{j\ge0}s^j\|P_j f_i\|_2^2\le p_i^{2/(1+s)}\).
For the simplex and for a list of \(k\) indicators of halfspaces
of mass \(u\), the sum of the degree-two-and-higher terms in
absolute value at correlation \(-r\) is bounded by
\(e^2 r^2L^2/k\) per list. The stability of the halfspace list is
\(k h_r(u)\), and its first-degree weight is \(k\phi(z)^2\).
Comparison with the simplex gives
\(S_{-r}(1_T)-k h_r(u)\le 13r/k\), using \(2e^2 rL^2\le1\)
for \(L\ge20\). Therefore \(\delta^2\le52r/k\) and \(p_i\le12/k\).
With this mass bound the total hypercontractivity estimate is
\(\sum_i p_i^{2/(1+s)}\le12e^2/k\). The exact stability expansion
for \(f\) and \(1_T\) gives
\[
 \Delta=F_k-W_1(f)+\delta^2/r\le110rL^2/k.
\]
This yields \(W_1(f)\ge L/(5k)\) for \(L\ge20\).

Write \(b_i=\int x f_i(x)\,d\gamma_d(x)\).
Gaussian exponential moments and Jensen under \(f_i/p_i\) give
\(\|b_i\|\le p_i\sqrt{2\log(1/p_i)}\), with the same bound by
continuity at zero. Monotonicity of \(p\sqrt{\log(1/p)}\) and
\(p\log(1/p)\) for \(p\le12/k\) then gives
\(\max_i\|b_i\|\le17\sqrt L/k\) and \(W_1(f)\le24L/k\).
In particular the normalized Gram matrix
\(Q_{ij}=b_i\cdot b_j/W_1(f)\) satisfies \(\max_i Q_{ii}\le2000/k\).

For the sharper mass estimate choose \(i\) with largest \(p_i\).
By (M5), \(\langle f_i,T_{-r}(f_i-f_j)\rangle\le0\) for every \(j\).
The absolute centered pairing with \(T_{-r}-P_0\) is at most
\[
 \frac r s\left(\|T_{\sqrt s}f_i\|_2^2+
             \|T_{\sqrt s}f_i\|_2\|T_{\sqrt s}f_j\|_2\right)
 \le \frac{4096rL}{k^2}.
\]
The second bound uses \(p_l^{1/(1+s)}\le12e/k\) for each
label \(l\). The constant pairing is \(p_i(p_i-p_j)\), and
\(p_i\ge1/k\). This proves
\[
 \max_i|kp_i-1|\le4096rL,\qquad
 \frac{L}{5k}\le W_1(f)\le\frac{24L}{k},\qquad
 \max_i Q_{ii}\le\frac{2000}{k}.
 \tag{EC2}
\]
In particular for \(L\ge3000\), the mass error in (EC2) is below
\(e^{-50}\).

## Quantitative localization

Apply [(VC9)](NEGATIVE_VOTING_COVARIANCE.md#numerical-constants-for-the-covariance-estimate)
with \(B_0=2000\), \(c=e^{-1201}\) and \(C=e^{1500}\).
The calculation below quantifies the localization estimates from
[(VL3)–(VL12)](NEGATIVE_VOTING_EXPONENTIAL_LOCALIZATION.md) with
these inputs. It also applies with other \(c,C>0\) satisfying (VC1)
on the given \(Q,p\), as long as \(L\ge1000\),
\(p_i\in[1/(2k),3/(2k)]\), and \(2CrL^2\le1\).
Set \(D_0=1+1/c\) and \(K=2^{76}D_0\).

### First and second moments

The preliminary maximum bounds give
\[
 a_*=\sqrt{F_k/(k-1)}\ge\frac{\sqrt L}{2k},\qquad
 \ell=\sqrt2 a_*\ge\frac{\sqrt L}{2k},\qquad z/\ell<3k .
\]
From (VC1),
\[
 \Delta\ge\frac{\delta^2}{2r}
  +c\bigl(\|\operatorname{diag}(Q-Q_*)\|_2^2
                   +\|Q-Q_*\|_F^2/k\bigr).
\]
For polar alignment, the scaled triangle inequality gives
\(\|B_f-B_T\|_F^2\le 48L\|Q-Q_*\|_F^2+2\Delta\).
Indeed \(Q_*^{1/2}\) is scalar on \(E\) and
\((\sqrt{W_1}-\sqrt{F_k})^2\le |W_1-F_k|\le\Delta\).
Here \(B_f\) and \(B_T\) denote the moment matrices in the polar
alignment of (VL4). Thus the regular score regret satisfies
\(\mathcal R\le 2^6D_0 kL\Delta\).

Write \(m_i=\|b_i\|/p_i\), \(m=\sum_i p_i m_i\) and
\(H=\sum_i p_i(m_i-m)^2\). For this centroid-length variance use \(\sqrt{kW_1}\) as the
comparison constant. Splitting the difference
\(\sqrt{W_1}(\sqrt{Q_{ii}}/p_i-\sqrt{k})\) into a diagonal and
a reciprocal-mass term gives
\[
 H\le4W_1 k^2\bigl(\|\operatorname{diag}(Q-Q_*)\|_2^2
                        +\|p-u\|_2^2\bigr)
  \le 2^8 D_0 kL\Delta .
\]
The preliminary centroid bound gives \(m_i\le\sqrt{2\log(1/p_i)}\le2\sqrt L\).
For example, Jensen under the normalized density \(f_i/p_i\)
gives \(e^{t m_i}\le e^{t^2/2}/p_i\) in the centroid direction.
Since \(W_1\le (\max_i p_i)(\max_i m_i)m\), one gets
\(m\ge\sqrt L/15\) and \(B\ge m^2>1\).
Write \(B\) and \(B_*\) for the second moments in (V12a)–(V15).
The unit simplex axial moment satisfies \(B_*\le6L\) by the upper
tail union bound for \(\max_i|Z_i|\) and the factor \(k/(k-1)\).

The mixed term in (V13) is bounded in modulus by
\(5\delta\sqrt{kLH}\le H+(25/4)kL\delta^2\).
The additional \(km^2\delta^2\) costs at most \(4kL\delta^2\).
Thus (V13) holds with its stated coefficients as \(22rL\le1\).
The smoothed coupling (V14) at cap \(R\) gives
\[
 B_*-B\le \frac{4R}{\sqrt L}(k\Delta+2H)+4k e^{-R^2/4}.
\]
From (V15), dropping the nonpositive \(G\) penalty and using
\(B>1\) and \(B_*\le6L\), the \(W_2\) difference is at most
\((12L/k)(B_*-B)_+\). Hence
\[
 W_2(T_k)-W_2(f)
 \le 2^{16}D_0 RL^{3/2}\Delta+64L e^{-R^2/4}.
\]

### Explicit normal integral controls

For \(a\ge0,\ j\ge1,\ a+j\le4,\ 0\le l\le4,\ 0\le D\le8\), (VL5)
holds in the form
\[
 \int_0^\infty (1+y)^a e^{Dy/z}\phi(y)^j
                 \Phi(y)^{k-l}\,dy
 \le 2^{40} L^{(a+j-1)/2}/k^j.
\]
To see this, with \(w=k(1-\Phi(y))\) and \(u=\log_+(1/w)\),
use \(1+y\le 3\sqrt L(1+\sqrt u)\),
\(\phi(y)\le3(1+y) w/k\),
and \(e^{D y/z}\le e^{13}\max(1,w^{-1/2})\) by
\(D\sqrt{2u/L}\le u/2+D^2/L\).
The product constants before integration are at most \(2^{30}\).
Use \(\Phi(y)^{k-l}\le e^{-w/2}\). For \(w\ge1\) the remaining
integral is at most \(2^j\Gamma(j)\le96\). For \(w<1\) use
\((1+\sqrt u)^3\le16(1+u^2)\) and the worst factor \(w^{-1/2}\),
giving at most \(16\int_0^\infty e^{-u/2}(1+u^2)du=288\).
The total is below \(2^9\).

For the simplex response (VL9) with \(1\le m\le3\), the ordered allocations of
derivatives to \(j\) distinct untied factors and the leading
density number at most \(2^8 k^j\). For \(y\ge0\), the shifted
Hermite factors cost at most \(2^m(1+t)^{m-j}\) after
\(t=y+1/z\) (each involved Gaussian argument differs from \(t\)
by at most \(2/z\), and \(|H_a^{\rm unnorm}(x)|\le(1+|x|)^a\)
for \(a\le3\)). Each density contributes at most
\(\phi(t)e^{2t/z}\). The preceding integral bound applies with
\(a=m-j,\ j\) there replaced by \(1+j,\ l=1+j,\ D=2(1+j)\).
For \(y<0\) all remaining CDF factors are \(\le3/4\).
Use the one leading density and bound all other Gaussian densities
by one, and the degree-at-most-three polynomial in the shifted
normal variable by \((3+|x|)^3\), whose Gaussian integral is below
\(128\). The resulting extra contribution is at most
\(2^{17} k^3(3/4)^{k-4}\le 1/k\).
Summing over \(j\), (VL9) follows with \(2^{60}L^{m/2}/k\)
in place of \(C_m z^m/k\).

On the cap \(\max|x_i|\le R\), assume \(R\ge1,\ 2rR\le1/z\).
Using \(|\kappa|\le2r,\ |\kappa-\rho|\le2r^3\),
\(|\kappa^2-\rho^2|\le2r^4\), Taylor's formula bounds each
normal derivative of the remainder by \(2^{65}r^3 L^{3/2}R^2/k\).
Apply this to each term of the score difference for the two cells. Use odd reflection and
divide by \(\ell\). The cap contribution is bounded by
\(2^{67}r^3 L R^2\mathcal R\). Outside the cap,
\(\|R_r1_{T_i}\|_2\le e r^3L^{3/2}/k\) by hypercontractivity
at the comparison scale \(1/L\). Summing Cauchy–Schwarz in \(i\)
with the cap-complement probability \(2ke^{-R^2/2}\) costs
at most \(4 r^3 L^{3/2}\sqrt{k}e^{-R^2/4}\).

### Explicit quadratic transfer bound

In (RT3), shift by \(1/(\sqrt2 z)\) and apply the explicit
integral bound with \(a=0,j=2,l=2,D=2\) on the positive part.
The negative part costs at most \(2^{-(k-2)}\).
The density factor \(2\sqrt\pi\) is at most four, giving
\(\gamma(\Omega_{ij})\le2^{43}\sqrt L/k^2\).
Consequently \(\gamma(\Omega_{ij})^{1-2/L}\le2^{49}\sqrt L/k^2\):
the displayed upper bound is at least \(k^{-2}\), and its
negative power costs at most \(e^4\).

Since \(q=L\ge1000\), the bound \(e^{3j/2}\) on normal Hermites
in (RT4) is at most \(q^{j/2}\). Cauchy–Schwarz across the
normal degrees costs \(\sqrt{m+1}\), so the squared tangential
norm control is at most \((2L)^m\). Thus the squared near
pairing (RT5) is at most \(2^{50}(2L)^m g_{ij}/k\).
For the far density grouped by label, use
\(\int H_{\rm far}\le2/k\), so \((\int H_{\rm far})^{1-2/L}\le2^{10}/k\).
With \(z/\ell<3k\), the squared far pairing costs at most
\(2^{12} L^m\) times the corresponding grouped regret.
The two near and two far contributions cost a factor four
after squaring their sum, and each transfer is counted twice
over labels. Therefore
\[
 \sum_i \|P_m e_i\|_2^2\le 2^{54}(2L)^m\mathcal R .
\]
Since \(2rL\le1/2\), the geometric tail gives (RT2) with
\(2^{60}r^3L^3\mathcal R\).

### Absorbing the remainder

Combining the estimates, including the factor two on the linear
tail, gives
\[
 \Delta\le
 K(r R L^{3/2}+r^2 k L^2 R^2+r^2 k L^4)\Delta
 +(64 rL+8r^2 L^{3/2}\sqrt{k})e^{-R^2/4}.
\]
Take \(R^2=1/(16K r^2 k L^2)\) and \(k\ge64K L^4\).
The three coefficients in parentheses after multiplication by
\(K\) are each at most \(1/16\). Also \(R\ge1\) and \(2rR\le1/z\).
The prefactor of the exponential is below one, proving
\(\Delta\le2 e^{-R^2/4}\).

The strip estimate in (V17) can use the constant \(4kt\).
For independent normals the integral on positive runner-up
height is at most \(h\,\mathbb E M_k\) by integrating
\(k(k-1)\phi^2\Phi^{k-2}\), and the negative height contributes
at most \(h k\phi(0)2^{-(k-1)}\). Rescale by
\(h=t/a_*\le 2kt/\sqrt L\).
Hence minimizing \(2(4kt+\mathcal R/t)\) gives summed label
distance at most \(8\sqrt{k\mathcal R}\). Substituting the regret bound
and \(r\le1/(k-1)\) gives
\[
 2^8\sqrt{D_0}\,k\sqrt L\,
     \exp\!\left[-\frac{k}{512K L^2}\right].
 \tag{EC3}
\]
If the local radius is \(k^{-A}\), this is less than \(k^{-A}\)
whenever \(k/L^4\ge2^{10}K(A+3+D_0)\). This stronger inequality
also implies the earlier \(k\ge64K L^4\).

## Comparing the two radii

For \(L\ge3000\), one has \(D_0\le e^{1202}\), \(K\le e^{1255}\),
\(2CrL^2\le4e^{1500}L^2/k<1\), and
\[
 k/L^4\ge e^{2960}
     >2^{10}K(2003+D_0).
\]
The first bound on \(k/L^4\) follows from monotonicity of
\(L-4\log L\). Take \(A=2000\) in (EC3). Every minimizer then
lies inside the explicit same-assignment local radius
[(LC6)](NEGATIVE_VOTING_LOCAL_CONSTANTS.md#an-explicit-single-assignment-local-radius).
That comparison forces the minimizer to be a regular simplex
assignment. Any assignment tying the simplex is itself a minimizer
and obeys the same argument. This proves (EC1).

The numerical cutoff is conservative. The argument uses first
variation only for minimizers, while the inequality and equality
classification apply over the full measurable fractional class.
All estimates use either the \(k-1\) active directions or
Gaussian bounds independent of the remaining dimension.
