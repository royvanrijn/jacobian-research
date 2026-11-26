# All-degree interlacing from one large-gap shared mirrored pair

Let ordinary exponent columns be fully shared and mirrored after
equal rows are omitted. Suppose one mirrored pair is
\[
                     (\ell,\ell+d),\quad(\ell+d,\ell),
                          \qquad d\ge7\ell,\quad \ell\ge1 .
\]
Then the consecutive ordinary factorial moment polynomials strictly
interlace in symmetric angular order in every degree. Other mirrored
pairs may have arbitrary positive smaller exponents and integer
gaps. In particular the two monomials satisfy SFC(2) over every
characteristic-zero field.

The selected pair supplies a stronger endpoint scale than its
quadratic curvature lower bound. The proof uses that scale in an
exact finite geometric estimate and in a paired analytic estimate
of the two endpoint logarithms. It uses no individual-degree
root tests.

## A hard-pair bound on the full coefficient array

Use endpoint normalization, with coefficients \(c_{m,k}\) at order
\(m\), including the binomial coefficient. Thus the two endpoints
are one and the full polynomial is reciprocal. Let \(s_0\) be the
total absolute curvature of all rows other than the chosen pair,
where each row \((a,b)\) contributes \((b-a)^2/\max(a,b)\). Put
\[
 s_h=\frac{2d^2}{\ell+d},\qquad S=s_h+s_0,\qquad
 q=\left(\frac{\ell}{\ell+d}\right)^d e^{-s_0/2},\qquad
 L_*=\frac{d^2}{\ell}+\frac{s_0}{2}.
\]
At \(k\le m/2\), the exact factorial multiplier from the chosen
pair is
\[
 \prod_{h=1}^{dk}\frac{m\ell+h}{m(\ell+d)-dk+h}
       \ \le\ \left(\frac{\ell+d k/m}{\ell+d}\right)^{dk}.
\]
The denominator starts above the numerator, so each ratio is at
most its value at \(h=dk\). The ordinary curvature bound for the
other rows supplies \(\exp[-s_0 k(1-k/m)/2]\).
[Endpoint curvature bound](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)
Writing \(w=k/m\), we have
\[
 \begin{aligned}
 c_{m,k}^{1/k}
    &\le \frac e w\left(\frac{\ell+dw}{\ell+d}\right)^d
                            e^{-s_0(1-w)/2},\\
 c_{m,k}
    &\le \frac{(m q)^k}{k!}\exp(L_* k^2/m).
 \end{aligned}                                                   \tag{1}
\]
The right side of the first inequality has its maximum on any
closed subinterval of \(0<w\le1/2\) at an endpoint. Indeed, with
\(r=d/\ell\), the logarithmic derivative has the sign of
\(-1+(dr-r+s_0/2)w+(s_0/2)rw^2\), which is strictly increasing.
At \(w=1/2\) this bound is less than
\(2e(9/16)^7<1/10\).

We use the symmetric half polynomial
\[
 H_m(z)=\sum_{0\le k<m/2}c_{m,k}z^k
                  +\tfrac12\mathbf 1_{m\ {\rm even}}c_{m,m/2}z^{m/2}.
\]
Then \(P_m(z)=H_m(z)+z^m H_m(1/z)\).
On any disk centered at zero on which \(H_m\) is zero-free,
take its logarithm to vanish at zero. If
\(\phi_m(\theta)=\Im\log H_m(e^{i\theta})\), then the centered
circle value is \(2|H_m(e^{i\theta})|\cos(m\theta/2-\phi_m)\).
Also \(\phi_m(0)=\phi_m(\pi)=0\) since that logarithm is real
on the real diameter.

Let \(N=n+1\) and \(D=Nq\). The case \(n=1\) follows from
\(P_1=1+z\) and the simple conjugate unit roots of \(P_2\).
Otherwise it suffices to show the two phases
\(\psi_m=m\theta/2-\phi_m\) strictly increasing on \([0,\pi]\),
and
\[
                 |\phi_N(\theta)-\phi_n(\theta)|<\theta/2
                             \qquad(0<\theta\le\pi).             \tag{2}
\]
Then each positive root is labeled by
\(\psi_m=(i-1/2)\pi\), and at a lower root the upper phase lies
strictly between the same and the next labels. Since each
\(\psi_m(\pi)=m\pi/2\), reflection and parity finish interlacing.

## Direct bounds in the short endpoint ranges

If \(N\le30001\), the first bound in (1) is less than \(1/10\)
at \(w=1/m\) as well. Drop the extra exponential and use
\(d/\ell\ge7\), \(d\ge7\) to bound that endpoint by
\(e m[(1+7/m)/8]^7\). Its maximum on \(2\le m\le30001\)
is at an endpoint because its logarithmic derivative has the sign
of \(m+7-49\). At \(m=2\) it is \(2e(9/16)^7<1/10\).
At \(m=30001\) use \((1+7/30001)^7<2\) and \(e<3\).
Consequently all nonconstant terms in either half have
coefficients bounded by \(10^{-k}\).

For larger \(N\), if \(S\ge2\log N+5\), the high-curvature
coefficient estimate gives the half bound \(8^{-k}\).
[Large-curvature geometric cap](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)
Both geometric bounds imply \(\|H_m-1\|_{|z|\le1}<1\) and
\(\|zH_m'/H_m\|_{|z|\le1}<1/4\), by their geometric sums.
These derivative bounds give (2) and the phase monotonicity.

We now assume \(N>30001\) and \(S<2\log N+5\), and put
\(R=\log N\). Since \(d\le s_h\le S\),
\[
 \frac{L_*}{n}\le
 \frac{(2R+5)^2+(2R+5)/2}{e^R-1}<1/32.                           \tag{3}
\]
The bound follows by decrease of the displayed function for
\(R\ge10\), using \(R>10\) and \(e^{10}>20401\).
If \(D\le1/10\), the second bound in (1) gives
\(c_{m,k}\le(1/8)^k/k!\) for \(1\le k\le7\), since
\(e^{7/32}<5/4\).
On \(8\le k\le m/2\) use the endpoint maximum in (1).
At \(k=8\) the bound is at most \(e D e^{8/32}/8<1/10\),
and at \(k=m/2\) it is less than \(1/10\) as above.
The exponential sum to order seven and the geometric tail
therefore give \(\|H_m-1\|,\|zH_m'\|<3/20\) on \(|z|\le1\).
This gives the same \(1/4\) logarithmic derivative bound.

## A paired estimate for the two long endpoint sums

In the remaining range \(D>1/10\), take
\(K=\lceil32D\rceil+14\). Write \(V\) for the sum of
\(d_i^2/\ell_i\) over all unequal rows and \(W=1+V\).
The bounds \(S\le V\le S(S+2)\) hold since every
\(\ell_i\ge1\) and \(d_i/\ell_i\le1+d_i^2/(\ell_i+d_i)\).
In particular
\[
 W q\le (s_h+s_0+1)^2(d/\ell+1)^{-d}e^{-s_0/2}
         \le (2d+1)^2 8^{-d}\le 225/8^7 .
\]
The first upper expression decreases in \(s_0\ge0\), since
\(s_h\ge d\ge7\). The last expression decreases in \(d\ge7\).
As \(q\le8^{-7}\) and \(D>1/10\), \(n>200000\).
Therefore
\[
 \frac{D}{n}<\frac1{720000},\quad
 \frac{WD}{n}<\frac1{9000},\quad
 K\le182D,\quad \frac{WK}{n}<\frac{182}{9000},\quad K\ge18 .
 \tag{4}
\]
Both halves contain all coefficients through \(K\), with
\(K+1<n/2\). For \(K\le k\le m/2\), use the endpoint maximum
in (1). At \(k=K\) the bound is at most
\((e/32)e^{L_* K/n}<1/8\), since \(L_*\le W/2\) and (4)
holds. At the midpoint the bound was already below \(1/10\).
This gives the required tail \(c_{m,k}\le8^{-k}\).

For the endpoint interpolation at real order \(m\in[n,N]\),
the first coefficient is
\[
 B_m=m \prod_i\prod_{h=1}^{d_i}
               \frac{m\ell_i+h}{m(\ell_i+d_i)-d_i+h},
\]
with the index \(i\) now counting each mirrored pair once.
The log derivative of each ratio is
\(-d_i(\ell_i+h)/[(m\ell_i+h)(m(\ell_i+d_i)-d_i+h)]\).
It follows that
\[
      0\le 1/m-(\log B_m)'\le V/m^2 .
\]
Here the second denominator factor is at least
\(m(\ell_i+d_i)/2\).
By (4), \(V/n<1\), so \(|(\log B_m)'|\le2/n\).
Also \(B_m/m\) decreases and the integer bound (1) at \(n\)
gives \(B_m\le D\,e^{L_*/n}<101D/100\).
For the last bound, \(R>12\) in the present range, so the
decreasing function in (3) is less than \(1/150\) using
\(e^{12}>128326\).

We give the inputs to the paired endpoint logarithm estimate.
The endpoint interpolation uses
\(\alpha=2\max(1,\max_i d_i/\ell_i)/n\) and \(\beta=8eW/n\).
Here \(\alpha K\le2WK/n<1/2\), \(4W/n<1/2\), and
\(n/2>K+1\), so the finite-difference and derivative bounds on
all edge weights apply.
[Finite-difference endpoint interpolation](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md#coefficient-halves-and-endpoint-interpolation)
For the connected-sum estimate the combined exponent is
\(\Xi=3\beta(1-\alpha K)^{-1}
          \exp(\beta K/[2(1-\alpha K)])\).
Indeed \(\alpha K<1/20\), \(\beta K<1/2\), and
\(\beta B_m<1/300\) by (4) and \(e<11/4\). The exponential
in \(\Xi\) is at most \(e^{1/3}<3/2\), so \(\Xi B_m<1/60\).
Thus the formal logarithm coefficients through \(K\) satisfy
\[
 |\ell_k(m)|\le\frac{B_m}{k}(1/60)^{k-1},\qquad
 |\ell_k'(m)|\le\frac{4B_m}{n}(1/60)^{k-1}.
\]
These are the paired connected-sum inequalities, using the
first-coefficient derivative just proved for \(k=1\).
[Connected-sum estimate](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#sharper-connected-sum)

Put \(A_m(z)=\sum_{k=1}^K\ell_k(m)z^k\).
The geometric series gives \(\|A_m\|_{9/5},\|A_m\|_{6/5}\le2B_m\)
and \(\|A_N-A_n\|_{6/5}\le5D/n\), with norms on the indicated disks.
At integer \(m\), \(e^{A_m}\) and \(H_m\) agree through degree
\(K\). On \(|z|\le6/5\), the exponential tail relative to
\(|e^{A_m}|\) is at most \(2(2/3)^K e^{4B_m}\) by Cauchy's
estimate on the larger disk. The tail of the actual polynomial
is at most \((3/20)^K e^{2B_m}\) relative to the same value.
Using \(K\ge32D+14\), \(B_m<101D/100\) and
\(\log(3/2)>2/5\), the first is less than \(2e^{-6}\),
and the second is below \(e^{-K}\). Their sum is less than
\(1/130\). Therefore the half polynomial is zero-free on this
disk and its analytic logarithm differs from \(A_m\) by less
than \(1/120\). In particular
\[
             \|\log H_m\|_{6/5}\le3(D+1),\qquad
             \|\log H_N-\log H_n\|_{6/5}<1/55 .
\]

On the unit circle the angular derivative of the logarithm
difference has modulus at most \(5/55\) by Cauchy's formula.
Its imaginary part is zero at \(\theta=0\), giving (2).
Each single phase correction has angular derivative at most
\(15(D+1)\), which by \(D>1/10\) and (4) is less than \(m/2\).
This gives phase monotonicity and hence the stated interlacing.
The integer coefficient resultant is nonzero over \(\mathbb Q\),
so the ordinary coprimality conclusion extends to every
characteristic-zero field.

## A criterion using all the endpoint factors

There need not be a designated large-gap pair. For any fully shared
ordinary mirrored profile, with the index counting each pair once,
put
\[
 Q(w)=\prod_i\left(\frac{\ell_i+d_i w}{\ell_i+d_i}\right)^{d_i},
 \qquad q=Q(0),\qquad W=1+2\sum_i d_i^2/\ell_i .
\]
The consecutive polynomials strictly interlace in every degree if
\[
          Q(1/30001)\le10^{-6},\qquad Q(1/2)\le1/55,\qquad
                               Wq\le10^{-4}.                    \tag{6}
\]
These are rational comparisons on the exact row data.

Multiplying the factorial ratio estimate over all pairs gives
\(c_{m,k}^{1/k}\le eQ(w)/w\) at \(w=k/m\le1/2\) and
\(c_{m,k}\le (mq)^k\exp(L k^2/m)/k!\), with
\(L=V/2=\sum_i d_i^2/\ell_i\).
Writing \(r_i=d_i/\ell_i\), the derivative of
\(\log(eQ(w)/w)\) after multiplication by \(w\) is
\(-1+\sum_i d_i r_i w/(1+r_i w)\), strictly increasing.
Thus the interval maximum is again at an endpoint.
Likewise the derivative of \(\log(mQ(1/m))\) after multiplication
by \(m\) is \(1-\sum_i d_i r_i/(m+r_i)\), which is increasing.
For \(N\le30001\) the two extreme checks are therefore
\(2eQ(1/2)\le2e/55<1/10\) and
\(e\cdot30001\,Q(1/30001)\le e\cdot30001/10^6<1/10\).
This proves the same geometric cap at both orders.

For \(N>30001\) the high-curvature branch \(S\ge2\log N+5\)
again uses the existing slope bound. In its complement,
\(L\le S(S+2)/2\). Since \(2R+5>1\), the right side of (3)
still bounds \(L/n\), by
\((2R+5)(2R+7)/2\le(2R+5)^2+(2R+5)/2\).
For \(D=Nq\le1/10\), the estimate for \(k\le7\) and the
endpoint maximum for \(k\ge8\) now apply as written, with
\(eQ(1/2)/(1/2)<1/10\) at the midpoint.

For \(D>1/10\), the condition \(q\le Q(1/30001)\le10^{-6}\)
gives \(n\ge100000\), hence \(D/n<1/720000\) and
\(WD/n<1/9000\). Choose \(K=\lceil32D\rceil+14\) as above.
All bounds (4) hold, and \(L/n\le W/(2n)<1/1800\) gives
\(B_m\le D e^{L/n}<101D/100\) from the integer left endpoint
and the monotonicity of \(B_m/m\). The bound on
\(L K/n\) and the midpoint estimate again give
\(c_{m,k}\le8^{-k}\) for \(K\le k\le m/2\).
These inequalities verify the same inputs to the inherited
finite-difference endpoint interpolation and connected-sum
estimate, including \(\alpha K<1/2\), \(4W/n<1/2\), the
first-coefficient derivative, and \(K+1<n/2\). Thus the
paired logarithm bounds and the phase comparison apply with
no new estimates of individual degrees.

The large-gap hypothesis implies (6). The first two conditions
follow from its selected factor, bounded respectively by
\(2/8^7\) and \((9/16)^7\). For the last, write
\(W\le w_0+s_0(s_0+2)\), \(w_0=1+2d^2/\ell\ge99\), and
\(q\le(\ell/(\ell+d))^d e^{-s_0/2}\). The function
\((w_0+s_0(s_0+2))e^{-s_0/2}\) decreases for \(s_0\ge0\),
so \(Wq\le(1+2d^2)/8^d\le99/8^7<10^{-4}\).
As examples with no designated large-gap pair, (6) holds for
twenty copies of \((1,2),(2,1)\), or three copies of
\((2,7),(7,2)\).
