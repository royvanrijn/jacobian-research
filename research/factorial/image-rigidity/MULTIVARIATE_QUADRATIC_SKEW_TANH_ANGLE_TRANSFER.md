# Quadratic signed skew in the exact matched angular transfer

Let \(d_i=b_i-a_i\), \(m_i=(a_i+b_i)/2\), and
\(A_i=m_i/|d_i|\) for the unequal rows of an ordinary
two-monomial profile. Put
\[
 W=\sum_i d_i^2/m_i,\qquad B_\infty=-\sum_i d_i^3/m_i^2,\qquad
 D=\sum_i d_i^4/m_i^3 .
\]
For each \(j\ge2\), let \(t_j\) be half the exact
order-\(j\) midpoint curvature, and set
\[
 \gamma_j=\min(\pi,\sqrt{32t_j/\pi}),\qquad
 M_j=\pi/\gamma_j,\qquad r_j=1/(128M_j^2),\qquad K_j=M_j+2 .
\]
These are the reference quantities of the
[midpoint-inflation root chart](MULTIVARIATE_MIDPOINT_INFLATION_POLYNOMIAL_ANGLE_TRANSFER.md#a-root-chart-along-a-short-complex-deformation).
If all \(A_i>K_j\), write
\[
 b_j=|B_\infty|+D/j,\qquad
 S_j=4r_j^{-4}\sum_i\frac{|d_i|}{(A_i-K_j)^3},\qquad
 E_j=2r_j^{-3}b_j+S_j .
\]

**Theorem.** If \(E_j<r_j/j\), the actual order-\(j\)
factorial roots and their exact midpoint Gaussian
reference roots can be paired, up to positive radial
scaling, with circular angular displacement at most
\[
                Q_j=2j r_j^{-5}b_j^2+(1+r_j^2/2)S_j.         \tag{1}
\]
Consequently the exact adjacent pair, at \(n\ge2\) and
\(N=n+1\), has disjoint root arguments if both root
charts apply and
\[
                   Q_n+Q_N<\frac{t_n}{36nN(1+t_n)} .         \tag{2}
\]
The proof allows off-circle roots and retains the
signed cubic curvature contribution quadratically.

The criterion is exactly expressible in matched tanh
data. For each matched Gauss factor put
\[
                 u_i=\left(\frac n{T'_i}-\frac{n-1}{T_i}\right)^{-1}.
\]
The
[matched Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)
gives \(u_i=|d|/(2m)\) on a row of midpoint \(m\)
and signed gap \(d\). Thus with orientation sign
\(\epsilon_i\) and sums over Gauss factors,
\[
       W=2\sum_i u_i,\qquad
       B_\infty=-4\sum_i\epsilon_i u_i^2,\qquad
       D=8\sum_i u_i^3,\qquad
       \sum_{\rm rows}\frac{|d|}{(A-K)^3}
               =\sum_i\left(\frac{2u_i}{1-2K u_i}\right)^3 .
\]
This retains exact sign cancellation in the projected
midpoint skew before estimating angular movement.
The match removes the internal Gauss shift from the
row rate, while the exact finite-order clocks still
use their factorial midpoint curvatures.

## Projecting the root chart onto the circle

Fix \(j\) and suppress the subscripts of \(r,M,K,t\).
Use the normalized midpoint inflation at parameter
\(\lambda\ge1\) and the root chart \(Z(w)\) from the
[midpoint-inflation proof](MULTIVARIATE_MIDPOINT_INFLATION_POLYNOMIAL_ANGLE_TRANSFER.md#interpolation-by-inflating-all-midpoints).
The hypothesis on \(E_j\) gives its uniform bound
\(|w_\alpha|\le j E_j<r\) for the root coordinates
\(w_\alpha=j\log(z_\alpha/\zeta_\alpha)\)
throughout that deformation.
Here \(\zeta_\alpha\) are the circle reference roots.
Set \(\kappa=\max_\alpha|\Re w_\alpha|\).
In the exact product for \(Z(w)\) replace each
\(w_\alpha\) in its cotangents by \(i\Im w_\alpha\),
leaving the argument \(w\) unchanged.
Call the resulting chart \(Z^\circ(w)\).
All its cotangents \(q^\circ_\beta\) are real.

At any one root \(\alpha\), the perturbed cotangents
satisfy
\[
      \sum_{\beta\ne\alpha}|q_\beta-q^\circ_\beta|
                           \le16M^2\kappa j .
\]
In fact the segment from each projected cotangent
argument to its actual argument has length at most
\(\kappa/j\), and the mesh estimate of the root chart
bounds the derivative there by \(4M^2j^2/k^2\) at
cyclic distance \(k\). Sum using at most two roots
at each distance. Similarly
\(|q^\circ_\beta|\le2Mj/k\) and
\(\sum |q^\circ_\beta-q^0_\beta|\le16M^2rj\), where
\(q^0\) denotes the Gaussian reference cotangents.

For \(|w|\le r\) we have \(|Z^\circ(w)|\le2|w|\).
Indeed the prefactor in the product satisfies
\(|2j\sinh(w/(2j))|\le|w|e^{|w|/(2j)}\).
Use the logarithm expansion at zero for the
remaining factors. The absolute linear exponent
is at most \(r\pi/(4t)+8M^2r^2\),
by the reference signed force bound and the
cotangent differences. The quadratic term costs
at most \(17M^2r^2\), using
\(\sum_{\beta\ne\alpha}(1+|q^\circ_\beta|^2)
 \le17M^2j^2\).
The log expansion is valid since
\(|q^\circ_\beta w/(2j)|\le Mr\le1/128\).
With \(r\pi/(4t)\le1/(16\pi^2)\) from the
reference mesh and \(r=1/(128M^2)\), the sum
of all three exponents including the prefactor
is less than \(1/8\).

Each projected factor
\(\cosh(w/(2j))-i q^\circ_\beta\sinh(w/(2j))\)
has modulus at least \(1/2\) for \(|w|\le r\)
by the same smallness. The relative change
of that factor when replacing \(q^\circ_\beta\)
by \(q_\beta\) has modulus at most
\(2|w|\,|q^\circ_\beta-q_\beta|/j\).
Therefore
\[
 |Z(w)-Z^\circ(w)|
 \le2|w|(e^{32M^2\kappa |w|}-1)
 \le128M^2\kappa |w|^2\le\kappa |w|
                                     \qquad (|w|\le r).
\]
Here the middle inequality uses \(\kappa\le r\)
so the exponent is small. Reality of the
projected cotangents gives
\(Z^\circ(-\bar w)=-\overline{Z^\circ(w)}\),
and hence its odd Taylor coefficients are real.
Cauchy on the radius-\(r\) circle now gives
\[
                       |\Im [w^3]Z|\le\kappa/r^2 .          \tag{3}
\]

## Angular integration

In the exact inflation velocity, the term of
order three in the root chart is \([w^3]Z\)
times the signed first jump moment
\(\int h\,\nu_\lambda(dh)\).
Its integral in absolute value over \(\lambda\)
is at most \(b_j\) by the signed curvature
calculation in the inflation proof.
The integrated absolute remainder, using the
fourth-order chart bound, is at most \(S_j\).
Taking the imaginary part of the velocity and
applying (3) with \(\kappa\le jE_j\) therefore
bounds the final angle change by
\[
                     j E_j r_j^{-2}b_j+S_j .
\]
Since \(E_j<r_j/j\) implies \(j b_j/r_j^2<r_j^2/2\),
substituting \(E_j=2r_j^{-3}b_j+S_j\) proves (1).
The exact reference clock ratio is in the
linked physical interval, by the
[matched midpoint clock test](MULTIVARIATE_MATCHED_MIDPOINT_CLOCK_ANGLE_TEST.md#the-gaussian-clocks-and-a-bound-on-root-rays).
The
[polynomial Gaussian gap](MULTIVARIATE_GAUSSIAN_POLYNOMIAL_PHYSICAL_GAPS.md)
then proves (2).

## One sufficient condition at every degree

For a chosen positive rational \(w\) assume \(W\ge w\) and write
\[
 r=\min\{1/128,w/(16\pi^3)\},\qquad
 h=\frac{w}{72(w+4)},\qquad H=|B_\infty|+D .
\]
Let \(J\) be any integer at least four, at least
the
[absolute ordinary eventual angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md),
and large enough that \(J\ge8h/r^3\).
A sufficient condition for angular separation
at every degree is
\[
                      2Jr^{-5} H^2+33r^{-4}D<h/J^2 .        \tag{4}
\]
In detail \(r_j\ge r\), and \(K_j\le K_w\)
where \(K_w=(128r)^{-1/2}+2\).
The bound on \(D\) in (4) is smaller than
\(r^4<(2K_w)^{-3}\), using \(K_w\le3(128r)^{-1/2}\).
Since each summand of \(D\) is at least
\(A_i^{-3}\), all \(A_i>2K_j\).
It follows that \(S_j\le32r^{-4}D\).
The bound on \(H\) gives
\(2r^{-3}H<\sqrt{2h}r^{-1/2}J^{-3/2}\le r/(2J)\).
Also \(32r^{-4}D<h/J^2\le r/(2J)\).
Thus the root-chart hypothesis holds for
all \(j\le J\), and (1) is at most
\(2Jr^{-5}H^2+33r^{-4}D\) there.
The sum of the two angular bounds meets
(2) at every \(N\le J\), because \(t_n\ge w/4\).
At \(N=2\) strict positivity of the even
moment on the negative ray gives separation.
All larger terminals use the absolute
ordinary cutoff. The constants in (4)
are unevaluated and give a sufficient
condition in the matched row moments.

## Cancellation of the first finite-order skew correction

Keep \(A_i=m_i/|d_i|\), write \(\epsilon_i=\operatorname{sign}(d_i)\),
and denote the signed cubic factor in the inflation velocity by
\(I_j(\lambda)=\int h\,\nu_\lambda(dh)\), \(\lambda\ge1\). Define the
row sums
\[
 H_3=\sum_i\frac{\epsilon_i}{A_i^3},\qquad
 R_4=\sum_i\frac1{|d_i|A_i^4}.
\]
Uniformly in the rows and the degree,
\[
 \left|I_j(\lambda)-\frac{2B_\infty}{\lambda^3}
                    -\frac{3H}{j\lambda^4}\right|
       \le\frac{2R}{j^2\lambda^5}.
 \tag{5}
\]
The same signed correction can be retained in the adjacent difference.
For \(N=n+1\),
\[
 \left|I_N(\lambda)-I_n(\lambda)+\frac{3H}{nN\lambda^4}\right|
       \le\frac{2R}{\lambda^5}(n^{-2}-N^{-2}).
 \tag{6}
\]

To see this, put \(f(s)=s/(e^s-1)=1-s/2+e(s)\). For \(s>0\),
\[
 e(s)=\frac{s}{2}\coth\frac{s}{2}-1,\qquad
                 0\le e'(s)\le s/6.
\]
The lower derivative bound is equivalent to \(\sinh s\ge s\).
The upper is equivalent to \(3\sinh s\le s(\cosh s+2)\).
In the latter difference the coefficients at \(s\) and \(s^3\)
vanish, and the coefficient at \(s^{2k+1}\) for \(k\ge2\) is
\((2k-2)/(2k+1)!\). Since \(e(0)=0\), for \(s\ge t\ge0\)
one has \(0\le e(s)-e(t)\le(s^2-t^2)/12\).
The inflation measure gives
\[
 I_j(\lambda)=-\sum_i\epsilon_i |d_i|A_i
     \int_0^\infty q^2 e^{-\lambda A_iq}f(q/(j|d_i|))\,dq .
\]
Integrating \(1-s/2\) gives the two terms in (5).
The integral of \(q^4 e^{-\lambda A_i q}\) is
\(24/(\lambda A_i)^5\). Applying the bound on \(e\) with
\(t=0\) proves (5). Applying it with \(s=q/(n|d_i|)\) and
\(t=q/(N|d_i|)\) proves (6).

Both quantities can be formed at a single reference pair. In the exact
joint Gauss grid let \(u_i=1/(2A_i)\) and \(c_{i,h}=2h/|d_i|-1\).
Using \(\sum_h c_{i,h}=1\) and
\(\sum_h (3c_{i,h}^2-1)=2/|d_i|\) gives
\[
 H_3=8\sum_{i,h}\epsilon_i c_{i,h}u_i^3,\qquad
 R_4=8\sum_{i,h}(3c_{i,h}^2-1)u_i^4.
 \tag{7}
\]
The identities keep the complete integer Gauss grid of each row.
The same formulas apply to every reference pair.

The proof of (1) uses \(b_j\) as an upper bound on
\(\int_1^\infty|I_j(\lambda)|\,d\lambda\), both for continuation
and for angular integration. Substituting \(x=1/\lambda\) in (5)
therefore permits the replacement
\[
 b_j\ \longleftarrow\
 \min\left\{|B_\infty|+D/j,
       \int_0^1|2B_\infty x+3H x^2/j|\,dx+\frac{R_4}{2j^2}\right\}.
 \tag{8}
\]
The integral splits at most once on \((0,1)\). In the all-degree
condition (4), the quantity \(|B_\infty|+D\) can accordingly be
replaced by
\[
   |B_\infty|+\min\{D,\ |H_3|/2+R_4/8\}.
\]
The original fourth-order term \(33r^{-4}D\) in (4) still applies.

The two signed moments can vanish for genuinely asymmetric ordinary
rows. For any positive odd \(L\), take increasing unit-gap rows
with smaller exponents \((3L-1)/2\) and \((7L-1)/2\) with
multiplicities \(27\) and \(343\), and decreasing unit-gap rows
with smaller exponent \((5L-1)/2\) with multiplicity \(250\).
Then \(B_\infty=H_3=0\) since
\(27/9-250/25+343/49=27/27-250/125+343/343=0\).
Taking \(m\) copies with \(L=108m+1\) makes the curvature
\(W=216m/(108m+1)\) tend to two, with both cancellations exact.
For such cancellations the signed cubic bound is \(O(j^{-2})\)
at one order, and its matched difference has integrated size
\(O(n^{-3})\). The comparison of the fourth and higher phase terms
for general ordinary rows remains open.
