# Explicit constants for the negative voting local comparison

Write \(L=\log k\) and take
\[
 k\ge \left\lceil e^{3000}\right\rceil,\qquad
 0<r\le\frac1{k-1},\qquad d\ge k-1 .
\]
The following estimates quantify the Hessian and polynomial local
neighborhood in
[VH4](NEGATIVE_VOTING_HESSIAN.md#the-full-negative-paired-gap) and
[PN1](NEGATIVE_VOTING_POLYNOMIAL_NEIGHBORHOOD.md).
The [explicit cutoff proof](NEGATIVE_VOTING_EXPLICIT_CUTOFF.md)
uses the local radius derived below.

## A weighted Hessian gap

Let \(\nu\) be the Gaussian surface measure of a simplex facet in
the active space, \(v_k=\nu(F_{12})\), \(\ell=kv_k\), and
\(\mu=\nu/v_k\). The identities, representation blocks and normalization
of the interaction are as in (VH1), (VH4)–(VH5) and (S82).
In particular \(D_r\) is the local form with multiplier \(a_r\),
and \(E_j\) is the nonnegative interaction form of Hermite degree
\(j\).

Here are finite bounds for the tied height \(Y\) with density
\(\phi(y)^2\Phi(y)^{k-2}/I_2(k)\) and the mean \(m_k=\mathbb E_\mu U\)
along the unit tie axis:
\[
 \frac{\sqrt L}{2k}\le\ell\le\frac{4\sqrt L}{k},\qquad
 1.74L\le\mathbb E Y^2\le2.26L,\qquad
 m_k^2\le4.6L,\qquad
 \gamma_k\le1-\frac{3.4L}{k}.
 \tag{LC1}
\]
For the first bound use \(\ell=\sqrt2\,\mathbb E M_k/(k-1)\),
with \(M_k\) the independent Gaussian maximum. Its upper expectation
bound \(\sqrt{2L}\) follows from the exponential moment at
\(\sqrt{2L}\). Its lower expectation is at least
\(\sqrt L[1-\exp(-\sqrt k/(8\sqrt L))]-1\): the Mills lower bound
\(1-\Phi(t)\ge \phi(t)t/(1+t^2)\ge e^{-t^2/2}/(8t)\) for \(t\ge1\)
gives the probability estimate at \(\sqrt L\), and
\((M_k)_-\le |Z_1|\).

Now \(I_2(k)=\ell/(\sqrt2 k)\ge\sqrt L/(4k^2)\).
With \(t_-=\sqrt{7L/4}\), the unnormalized integral of
\((1+y^2)\phi(y)^2\Phi(y)^{k-2}\) on \(y<t_-\) is at most
\(\exp[-k^{1/8}/(32\sqrt L)]\), by the same Mills lower bound.
On \(y>t_+=\sqrt{9L/4}\) the integral is at most
\((t_++1)e^{-t_+^2}\), by integrating the Gaussian square.
Each bound divided by \(I_2(k)\) is less than \(10^{-3}\) on the
stated range. On \([t_-,t_+]\) this puts the second moment between
the bounds in (LC1). For example the tail ratios are bounded by
\(\frac{4k^2}{\sqrt L}\exp[-k^{1/8}/(32\sqrt L)]\) and
\(12k^{-1/4}\), respectively. The mean identity in (VH6) gives
\(m_k^2=\frac{2k}{k-2}(\mathbb E Y)^2\le4.6L\).
Finally, \(\lambda(y)=\phi(y)/\Phi(y)\) satisfies
\((k-2)\mathbb E[Y\lambda(Y)]=2\mathbb E Y^2-1\) by integrating
the derivative of \(y\phi(y)^2\Phi(y)^{k-2}\).
The conditional contrast variance is
\(\gamma_k=1-\mathbb E[Y\lambda(Y)+\lambda(Y)^2]\), which yields
the final bound in (LC1).

Put \(\sigma=\sqrt{1-r^2}\), \(\beta=r^2/(2\sigma^2)\),
\(\gamma=r m_k/(2\sigma)\), \(X=\beta|x|^2+\gamma U\) and
\(A_r=a_r/(r\ell)\). By (VH5), \(A_r^{-1}\le e^X\).
The centered linear moment bound in that proof gives
\(\mathbb E e^{tZ}\le e^{t^2/2}\) for \(Z=U-m_k\) and for
any unit remaining-label contrast \(Z=V\). Its second derivative
at zero gives \(\|Z\|_2\le1\). Integrating the resulting tail
bound \(2e^{-t^2/2}\) gives \(\|Z\|_4\le2,\ \|Z\|_8\le3\). Also
\(\||x|^2\|_4\le2k\) from the \(\chi^2_{k-2}\) law on the facet.
On our interval,
\[
 \beta\le .52/k^2,\qquad
 0\le\gamma\le 1.1\sqrt L/k,\qquad
 \gamma m_k\le 2.4L/k.
\]
Consequently \(\|X\|_4\le3L/k\) and \(\|e^X\|_4\le2\).
For the latter, Cauchy–Schwarz bounds \(\mathbb E e^{4X}\)
using \((1-16\beta)^{-(k-2)/2}\le2\) and
\(\exp(8\gamma m_k+32\gamma^2)\le2\).
Since \(U\ge0\) on the facet, the bound
\(0\le e^X-1-X\le X^2e^X/2\), Hölder and (LC1) give
\[
\begin{aligned}
 \mathbb E_\mu[V^2/A_r]
 &\le \gamma_k+2.4L/k+4.2/k+4.4\sqrt L/k+81L^2/k^2
 \le 1-.9L/k,\\
 \mathbb E_\mu[(U-m_k)^2/A_r]
 &\le1+24L/k\le1.01.
\end{aligned}
\tag{LC2}
\]
The \(4.2/k\) term bounds
\(\beta\mathbb E[V^2|x|^2]\) by Cauchy–Schwarz.
The centered linear term uses
\(|\mathbb E[V^2(U-m_k)]|\le4\).
The second line uses \(e^X-1\le Xe^X\).

For the degree-two form, the probability of the facet within its
supporting plane is \(v_k/\phi(0)\ge k^{-2}/2\).
A unit degree-two Hermite polynomial restricts on that plane to
a degree-two polynomial of Gaussian \(L^q\) norm at most \(q\)
for \(q\ge2\), by the decomposition and hypercontractivity in
the proof of (VH8). For \(q=8L\), Hölder bounds its conditional
\(L^4(\mu)\) norm by \(16L\). Thus
\(\mathbb E_\mu[P^2/A_r]\le512L^2\), and summing the coefficient
pairings over incident facets as in (VH8) gives
\[
 r^2 E_2(f)\le 1024 r L^2\, D_r(f,f).
 \tag{LC3}
\]

For the higher-degree tail use \(b=1/2\) in the dilation argument
of (VH8). Young's inequality in (VH5) gives
\(a_r(x)^{-1/2}\le e(r\ell)^{-1/2}e^{\delta |x|^2/2}\)
with \(\delta=\beta+1/k\le2/k\).
The new correlation \(b'\) solves
\((b/b')(1+(b')^2)=1+b^2-2(1-b^2)\delta\).
On \([1/2,3/5]\) the derivative of the left side is at most
\(-8/9\), so \(b'\le1/2+4/k\).
The dilation factor is at most \((1+8/k)^{k/2}\le e^4\).
The absolute row bound at \(b'\) is less than 5, as in (S81).
The positive comparison kernel therefore has normalized operator
norm at most \(3000/(r\ell)\).
Multiplying by \((r/b)^3\) and using (LC1) gives
\[
 \left|\sum_{j\ge3}(-r)^j E_j(f)\right|
       \le \frac{50000}{k\sqrt L}\,D_r(f,f)
       \le .31\frac{L}{k}\,D_r(f,f).
 \tag{LC4}
\]

For \(D_r+K_r\), the two-row bound in (LC2) leaves a gap at least
\(.5 L/k\) after (LC4). In the standard block the centered
degree-one coefficient is bounded by \((5/8)1.01\).
The mean term costs at most \((5/2)m_k^2\le 12L\), which is
absorbed by \(1/r\) when the constant term is retained. The
whole standard block therefore has gap at least \(1/4\).
The rotation-hook rank argument following (VH8) gives gap at
least \(3/4\) off the rotations, since (LC3) and (LC4) bound
its remainder by \(1/4\). The rest of the symmetric blocks has
gap at least \(3/4\), by (LC4).
In \(D_r-K_r\) with matched masses the constant term vanishes
and degree one has favorable sign. Equations (LC3) and (LC4)
give gap at least \(1/2\). These inequalities and the lift
(NL7) give a weighted paired Hessian gap of at least
\[
 \frac{L}{4k}
 \tag{LC5}
\]
modulo common rotations in every \(d\ge k-1\).
With centered interaction, the same bound holds on the
active symmetric standard block after subtraction of
\(12L\|Cu\|^2/\ell\) in the norm weighted by \(A_r\), as in
(PN11). Removing the constant term from the active
antisymmetric form gives the \(1/2\) gap on all speeds there.

## An explicit single-assignment local radius

For the global application only the same-assignment specialization
of (PN1) is needed. On the range above, a sufficient radius is
\[
 \inf_T\sum_i\|f_i-1_{T_i}\|_1\le k^{-2000}.
 \tag{LC6}
\]
Here and below \(T\) ranges over cylindrical regular simplices.
Inside this radius \(S_{-r}(f)\ge S_{-r}(T)\), with equality only
on the simplex orbit. We quantify the estimates in
[the measurable projection argument](NEGATIVE_VOTING_POLYNOMIAL_NEIGHBORHOOD.md)
with ample slack.

Use the active coordinates \(x\in E=\mathbf1^\perp\), the unused
coordinates \(z\), and the notation \(A_r,U_{ij}\) of (PN2)–(PN9).
Set \(J=\sum_h\int c_h\,d\gamma(z)\) summing over the ordered
transfers of the single assignment. Its linear score cost
in the stability expansion is \(2r\ell J\).
For \(y\) anywhere in the tie plane, (LC1)
and the Jensen bound (PN3) give \(A_r(y)\ge e^{-1}e^{-2|y|^2/k}\).
Indeed \(r m_k/(2\sigma)\le1.1\sqrt L/k\), and
\(1.1\sqrt L\,|y|/k\le |y|^2/k+L/k\) by Young's inequality.
Also \(\beta\le .52/k^2\). Thus
\(U_{ij}(x)\ge e^{-1}e^{-2|x|^2/k}\min(t,1)\)
by (PN2). Differentiating the conditional Gaussian in a unit
normal gives \(|V_{ij}|\le2rt\) on \(t\ge0\) and
\(A_r\le2/\ell\le4k\).

Take \(R=16\sqrt{k}\) as the active cap and
\(W(x)=(1+|x|^2)e^{4|x|^2/k}\).
The elementary moment bounds we use, in any dimension \(m\le k\),
are
\[
 \int W(x)^2 e^{2|x|^2/k}\,d\gamma_m(x)\le e^{30}k^2,\qquad
 \int_{|x|>\sqrt{R^2-1}} W(x)^2 e^{2|x|^2/k}\,d\gamma_m(x)
        \le e^{30}k^2 e^{-k}.
\]
For both bounds use \((1+u)^2\le4k^2e^{u/k}\) for \(u\ge0\).
For the first integrate \(e^{11|x|^2/k}\) using the chi-square
moment generating function. For the second insert
\(e^{(|x|^2-(R^2-1))/8}\) on the indicated region and use the
same formula.

For a transfer \(h\), the active mass \(j_h(z)=\int W h\) and
the normalized active cost \(c_h(z)=\int U_{ij}h\) satisfy
\(j_h(z)^2\le e^{40}k^2 c_h(z)\) by normal filling for \(t<1\)
and weighted Cauchy–Schwarz for \(t\ge1\).
In the first case \(W(y+tn)\le 3 W(y)\) and
\(e^{2|y+tn|^2/k}\le 2 e^{2|y|^2/k}\).
Summing the squared bounds for each incident label costs at
most \(2k\), by Cauchy–Schwarz. Thus the \(\ell^2\) norm over
labels of the active mass \(L^2(d\gamma(z))\) norms is at most
\(e^{25}k^{3/2}\sqrt J\). For the exterior the same bound gains
\(e^{-k/2}\) by the weighted radial estimate above.

Discard transfers in those exteriors, in \(t>\xi\) on the cap,
and in \(0\le t\le\xi\) with the normal foot outside the genuine
facet. On the cap, \(W\le e^{1040} k\) and the reciprocal of
the weight in the lower bound for \(U_{ij}\) is at most
\(e^{513}\). Hence the sum of the \(L^1(d\gamma(z))\) masses
for long-normal transfers, including incidences, is at most
\(e^{1600} k J/\xi\). Interpolate this bound with the \(L^2\)
bound label by label and apply Hölder in the label sum.
For off-facet projections, the Gaussian measure of the plane
strips is at most \(2k\xi\). Using the cap bound on \(W\) in
weighted Cauchy–Schwarz gives squared slice mass at most
\(e^{2700} k^3\xi c_h(z)\) per transfer. In particular the
\(\ell^2\) norm over labels of the discarded \(L^{3/2}\)
slice norms divided by \(\sqrt J\) is bounded by
\[
 E_{\rm bad}=
  e^{3000} k^2
   \left(e^{-k/2}+\sqrt{\xi}+(\sqrt J/\xi)^{1/3}\right).
 \tag{LC7}
\]

For clarity, the constants in the kernel and gauge comparison
can be bounded as follows. For \(0\le r\le1/(k-1)\), the active
Mehler prefactor is at most 2 and its exponent is at most
\((|x|^2+|y|^2)/k\). Differentiating its formula in \(r\)
bounds \(|M^x_{-r}-1|/r\) by \(32W(x)W(y)\).
Differentiation in \(x,y\) on the cap bounds the spatial
derivatives by \(e^{600}r\sqrt{k}\).
Use (PN6), the contraction bounds there, and
\(\ell^{-1}\le 2k\). The normalized centered error incurred
by discarding and projecting is bounded by
\[
 e^{4000}k^4(E_{\rm bad}+E_{\rm bad}^2)
      + e^{4000}k^6\xi .
\]
In the projection term, the transverse centered kernel is
unchanged because each retained slice mass is preserved. The
other term uses the active derivative bound, the displacement
at most \(\xi\), and the two mass norms from (PN7).
The normalized moment-gauge error is at most
\(e^{30}(E_{\rm bad}+k^{3/2}\xi)\) by testing each unused
direction with its \(L^3\) Gaussian norm and each active
direction against \(W\). The normalized mass error is at most
\(E_{\rm bad}\).

The weighted facet moment map has norm at most \(e^{10}k\).
For each row this follows by testing against a unit direction,
using \(\int(1+|x|^2)A_r^{-1}\,d\nu\le e^{10} k v_k\),
then adding at most \(k-1\) interface contributions by
Cauchy–Schwarz. The transverse unit test has Gaussian second
moment 1. The estimate uses the radial \(\chi^2_{k-2}\) law
of the normalized facet measure. The inverse of this map on
geometric rotations, from the tangent moment matrix to
weighted facet speeds, has norm at most \(e^{10} k^2\).
Indeed the regular moment singular value \(\ell/\sqrt2\)
controls the active and mixed skew blocks. The weight upper
bound \(2/\ell\), the surface moments, and
\(\sum_{i<j}n_{ij}n_{ij}^T=(k/2)I_E\) bound the norm of the
corresponding rotation speed by \(2\sqrt{k}\) times the norm
of those blocks. A weaker \(e^{10}k^2\) bound on the inverse
therefore suffices. Applying these two map bounds to the
moment gauge puts
\(\operatorname{dist}_{D_A}(f_{\rm flux},\mathcal R)^2
 \ge (1+e^{40} k^6)^{-1}D_A(f_{\rm flux},f_{\rm flux})\)
for the rotation-corrected speed \(f_{\rm flux}\).

The centered facet form divided by \(r\ell\) also has operator
norm at most \(e^{40}k^3\). To obtain this directly use (PN6)
and the centered kernel bound \(32r W(x)W(y)\) for boundary
measures. For each facet \(\int W^2 A_r^{-1}d\nu
\le e^{40} k^2 v_k\) by the radial law. Weighted
Cauchy–Schwarz and the \(k-1\) incidences per label give
the bound, since \(kv_k=\ell\). Consequently rotating to the
exact moment gauge changes the squared norm and the
interaction by at most \(e^{4000}k^8(E_{\rm bad}+\xi)\)
when that quantity is small.

Use the weighted Hessian gap and the centered mass penalty
following (LC4). After enforcing the rotation gauge the
coercivity is at least \(2k^{-9}\) for \(k\ge e^{3000}\), by
the preceding norm comparison and the full gap \(L/(4k)\).
The centered penalty is at most \(12L\|Cf_{\rm flux}\|^2/\ell\).
The normalized original mass term is
\(\|p-u\|^2/(r\ell J)\), and
\(\|Cf_{\rm flux}\|^2\le 2\|p-u\|^2/J+2E_{\rm bad}^2\).
Thus the true mass term pays this penalty because \(24rL\le1\).
Using (PN9) and the preceding error bounds, the normalized
deficit for the same-partition objective, in units \(r\ell J\),
is at least
\[
 2k^{-9}-e^{5000}k^9(E_{\rm bad}+E_{\rm bad}^2+\xi).
 \tag{LC8}
\]
All bounds assume the error factors are at most one, as
verified next.

Write \(\eta=d_T(f,f)/2\). First Hermite projection puts
the moment of \(f\) within \(2\sqrt\eta\) of any initial
alignment at label distance at most \(2\eta\). The nearest
moment simplex changes the initial reference moment matrix
by at most \(4\sqrt\eta\). The union bound over Gaussian
wedges, with pair normal change at most \(8\sqrt\eta/\ell\),
gives aligned label distance
\(\eta_*\le 10 k^3\sqrt\eta\). The global score upper bound
\(2rt\) and Cauchy–Schwarz imply
\(J\le (2/\ell)\sqrt{k\eta_*}\).
Consequently \(\sqrt J\le 8 k^2\eta^{1/8}\le k^{-247}\)
under (LC6). Take \(\xi=k^{-80}\).
Formula (LC7) gives \(E_{\rm bad}\le k^{-36}\) for
\(k\ge e^{3000}\). The error in (LC8) is then less than
\(k^{-9}\). For \(J>0\) the deficit is strictly positive.
If \(J=0\) the positive score margin forces all transfers
to vanish. This proves (LC6).
