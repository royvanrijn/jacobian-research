# Exact angular separation from a uniform matched Gauss remainder

Let \(N=n+1\ge3\), and let the ordinary exponent columns be distinct
and fully shared on their unequal rows. At order \(n\) let
\(X=\sum_i x_{n,i}\) be the total of the matched Gauss reciprocals.
The index includes every Gauss factor within a row.

**Theorem.** If
\[
                         X\le \frac{10^{-16}}{16}N^{-29/9},
 \tag{1}
\]
the two factorial moments have simple roots and disjoint root
arguments under independent positive rescalings. For any ordinary
distinct exponent columns, the sufficient condition
\[
          S=\sum_i \frac{(b_i-a_i)^2}{\max(a_i,b_i)}
                    \le 10^{-18}N^{-20/9}
 \tag{2}
\]
implies (1) and hence coprimality over every characteristic-zero field.

This gives an exact angular range by controlling the analytic
remainder of the
[matched fourth-order model](MULTIVARIATE_MATCHED_SECOND_ORDER_ANGLES.md).
For sufficiently large degrees, (2) extends the earlier
\(10^{-16}N^{-8/3}\) range. The two explicit bounds can be used together.

The [sixth-order strengthening](#sixth-order-strengthening) enlarges
the sufficient bounds to \(X\le10^{-16}N^{-19/6}/16\) and
\(S\le10^{-18}N^{-13/6}\). The quartic proof below establishes the
same complex disk using elementary bounds and is independent of
the sixth-order calculation.

## A uniform complex disk for the matched deformation

Write \(\xi_{j,i}=x_{j,i}/X\) at \(j=n,N\) and
\(R_j=\sum_i\xi_{j,i}\). Thus \(R_n=1\) and
\(n/N\le R_N\le1\). At each order deform the reciprocals to
\(\lambda^2\xi_{j,i}\) with the same complex parameter \(\lambda\).
The actual common parameter is \(\sqrt X\). Put
\[
   c=10^{-8},\qquad L=c N^{-3/2},\qquad
   r=\frac1{80\sqrt N}.
\]
For every \(|\lambda|\le L\), the exponential-coordinate polynomial
\[
 C_j(z,\lambda)=\lambda^{-j}
    \sum_{k=0}^j(-1)^k\binom jk e^{-\lambda v z} F_k(\lambda^2),
                 \qquad v=k-j/2,
 \tag{3}
\]
has exactly one zero in each radius-\(r\) disk about a root of
\[
                  H_j(z)=j![u^j]e^{zu+R_j u^2/2}.
\]
Here \(F_k\) is the normalized Gauss product defined in the
[matched angular expansion](MULTIVARIATE_MATCHED_SECOND_ORDER_ANGLES.md#exact-angular-coefficient).
Each zero is simple and analytic in \(\lambda\) for \(|\lambda|<L\).
The singularity at zero in (3) is removable by alternating
cancellation. We prove the disk assertion by bounding the exact
Gauss deformation on the Hermite root circles.

Let \(A_{j,a}(v)\) be the centered single-factor log coefficient of
\(\lambda^{2a}\) for reciprocal \(\xi=1\) in the increasing orientation:
\[
 A_{j,a}(v)=\frac{(-1)^{a+1}}{a}
 \left(\sum_{h=0}^{k-1} h^a-\frac{k}{j}\sum_{h=0}^{j-1}h^a\right),
                        \qquad k=v+j/2.
\]
The sums in the first term are interpreted as the Faulhaber
polynomial in \(k\). Reflection replaces \(v\) by \(-v\).
The bounds in the proof apply to either sign.

For \(a\le2N\) and complex \(v\), one has
\[
                   |A_{j,a}(v)|\le9(|v|+4N)^{a+1}.
 \tag{4}
\]
Indeed the generating function \(u/(e^u-1)\) has modulus at most
\(4\) on \(|u|\le1\). This follows from
\(\left|(e^u-1)/u-1\right|\le e-2<3/4\).
The Bernoulli coefficients consequently satisfy
\(|B_m|\le4m!\). For \(m\le2N+1\), this is at most \(4(3N)^m\).
Expanding the Bernoulli polynomials gives
\[
 \left|\frac{B_{a+1}(k)-B_{a+1}(0)}{a+1}\right|
                \le\frac{8}{a+1}(|v|+4N)^{a+1}.
\]
The displayed quotient equals the first power sum by the
Bernoulli generating function, obtained by multiplying
\(\sum_{h=0}^{k-1}e^{hu}\) by \(e^u-1\) at integer \(k\).
The chord term is at most \((|v|+N/2)N^a\), using the integer
sum at \(j\). These bounds give (4).

Truncate the log at \(a=2N\), and substitute \(w=\lambda v\).
Write \(\mathcal L(w,\lambda)\) for the resulting polynomial for
the whole profile, summing the factors with their orientation
and powers \(\xi_{j,i}^a\). The substitution is regular at
\(\lambda=0\) since \(\deg A_{j,a}\le a+1\le2a\).
The first log term is \(R_jw^2/2-R_j j^2\lambda^2/8\).
Write \(\alpha=|\lambda|\) and suppose \(|w|\le13\sqrt N\).
Since \(\sum_i \xi_{j,i}^a\le1\), (4) and
\(|w|+4N\alpha\le14\sqrt N\) bound the remaining terms by
\[
 \sum_{a=2}^{2N}9\alpha^{a-1}(14\sqrt N)^{a+1}
                        \le49392\alpha N^{3/2}.
 \tag{5}
\]
The geometric ratio \(14\alpha\sqrt N\) is less than \(1/2\).
The first scalar term \(R_j j^2\lambda^2/8\) has modulus at most
\(\alpha N^{3/2}/8\).

## Centered finite differences give a relative circle bound

Let \(C_j^{\rm tr}\) be (3) with the truncated log multiplier.
For any entire function \(G\), a centered difference has the
exact integral formula
\[
 G(x+\lambda/2)-G(x-\lambda/2)
            =\lambda\int_{-1/2}^{1/2}G'(x+\lambda s)\,ds.
\]
This holds for complex \(\lambda\) by integrating along the straight
segment. Iterate \(j\) times and use
\(G(w)=e^{zw+\mathcal L(-w,\lambda)}\). With \(Y_j\) the sum
of \(j\) independent uniforms on \([-1/2,1/2]\), one obtains
\[
       C_j^{\rm tr}(z,\lambda)
          =\mathbb E\,G^{(j)}(\lambda Y_j)
          =\mathbb E\,j![u^j]G(u+\lambda Y_j).
 \tag{6}
\]
This is an exact identity, and remains valid at zero by continuity.

Use \(R_0=12\sqrt N\) as the coefficient contour radius.
For \(|u|\le R_0\), \(h=\lambda Y_j\) satisfies
\(|h|\le\alpha N/2\) and \(|u+h|\le13\sqrt N\).
On the Hermite root circles \(|z|\le3\sqrt N\).
The error exponent relative to the Gaussian generating function is
\[
 zh+R_j(uh+h^2/2)
       +\mathcal L(-(u+h),\lambda)-R_j(u+h)^2/2 .
\]
The first two terms are at most \(8\alpha N^{3/2}\) in modulus.
Combining (5) and the first scalar log term bounds the whole
exponent by \(50000\alpha N^{3/2}\).
This is at most \(1/2\), so the exponential minus one has modulus
at most \(10^5\alpha N^{3/2}\) on \(|u|\le R_0\).
Cauchy bounds its \(u^m\) coefficient by that quantity divided by
\(R_0^m\).

On the same circles the
[Hermite derivative bound](MULTIVARIATE_EFFECTIVE_WEAK_CURVATURE_ANGULAR_PARITY.md#individual-root-disks)
gives
\[
      \left|H_j^{(m)}(z)/H_j(z)\right|
                             \le30(6\sqrt N)^m,\qquad m\ge1.
\]
This bound only requires the Hermite variance in \([n/N,1]\)
and the stated same-order radius. Extracting \(u^j\) in (6)
against the Gaussian exponential gives
\[
  \left|C_j^{\rm tr}(z,\lambda)/H_j(z)-1\right|
                  \le 31\cdot10^5\alpha N^{3/2}\le31/1000 .
 \tag{7}
\]

It remains to bound the discarded log terms at the discrete
coefficient indices. Set \(\tau=N\alpha^2<1/2\).
For integer \(0\le k\le j\) the full and truncated logs have
modulus at most \(4N\tau<1\), by bounding the two power sums in
each chord by \(2N N^a\). The discarded log tail has modulus
at most \(4N\tau^{2N+1}\). Thus the two Gauss multipliers differ
by at most \(12N\tau^{2N+1}\).
On the root circles the exponential in (3) has modulus at
most two. The same-order Hermite separation is at least
\(1/(20\sqrt N)=4r\), so the monic \(H_j\) has modulus
at least \(r^j\) on the circle at each root.
Therefore for \(\alpha>0\),
\[
 \begin{aligned}
 \left|\frac{C_j-C_j^{\rm tr}}{H_j}\right|
  &\le24N(160\sqrt N)^j\alpha^{4N+2-j}N^{2N+1}\\
  &\le24c^2(160c^3)^N N^{-2N-1}<1/100.
 \end{aligned}
\]
The inequalities also hold at zero by continuity.
Together with (7), Rouché proves the asserted simple analytic
root branches on the complex disk.

## Remainder after the fourth-order interlacing model

For a root \(iy\) of \(H_j\), let \(z_{j,y}(\lambda)\)
be its exact branch. It remains within distance \(r\)
of \(iy\), so the Taylor coefficients of its displacement
have moduli at most \(r/L^m\) in degree \(m\).
The reflection identity
\(C_j(z,-\lambda)=(-1)^j C_j(-z,\lambda)\) and real coefficients
give \(z_{j,y}(-\lambda)=-\overline{z_{j,y}(\lambda)}\)
for real \(\lambda\). The odd-degree displacements are real.
The imaginary part, at the common real \(\lambda=\sqrt X\),
is the fourth-order angular model plus an error at most
\[
       r\frac{\eta^6}{1-\eta^2},
                     \qquad \eta=\frac{\sqrt X}{L}.
 \tag{8}
\]
In (1), \(\eta\le\tfrac14 N^{-1/9}\), so the sum of two
possible errors is at most \(N^{-7/6}/153600\).

For the fourth-order model the correction maps have slope
at least \(1/2\) on the whole Hermite height range.
At a common \(y\) with \(y^2\le4N\), the difference between
their corrections is at most
\[
 \left(\chi+25\chi^2+80000\chi^2(1+\chi)\right)\frac{|y|}{n^2},
                  \qquad\chi=nX.
\]
These are the explicit matched bounds in the
[fourth-order model proof](MULTIVARIATE_MATCHED_SECOND_ORDER_ANGLES.md#the-fourth-order-angular-model).
Here \(\chi\le1/1000\), so the coefficient in parentheses
is less than \(1/12\).
For either positive adjacent pair, the underlying Hermite gap
is at least \(y/(4n^2)\), with \(y\) the order-\(n\) height,
by the Jacobi comparison in that proof.
The model gap is consequently at least one sixth of the
underlying gap. Pairs involving the zero height have
at least half the underlying gap, since both correction
maps fix zero. In particular, the model gaps are at least
\(N^{-7/6}/19200\) by the
[scaled Hermite edge estimate](MULTIVARIATE_EFFECTIVE_WEAK_CURVATURE_HERMITE_EDGE.md).
The remainder (8) is smaller, so the positive angular
heights of the exact branches retain their disjoint order.
The odd-degree zero height stays real in the exponential
coordinate, and the nonzero exact heights keep their signs
by the same-order spacing and (8). Conjugation covers the
negative heights.

These branches account for all physical roots. The exponential
map \(t=-e^{-\lambda z}\) is injective on their disks since
\(6\lambda\sqrt N<2\pi\), and its derivative is nonzero.
Endpoint normalization rescales \(t\) by a positive factor
at each order, preserving the angular comparison.

Finally (2) forces every unequal row to be shared.
For \(S<1\), each row has \(d_i/(a_i+d_i)\le S\), with \(a_i\)
the smaller exponent, and hence
\[
      X\le\sum_i\frac{d_i^2}{n a_i}\le\frac{S}{n(1-S)}.
\]
Here \(n\ge2N/3\) and \(S\le1/2\), so (2) gives
\(X\le3\cdot10^{-18}N^{-29/9}\), which satisfies (1).
Rational Bézout transfers the resulting coprimality to
every characteristic-zero field.

## Sixth-order strengthening

The [sixth-order checker](verify_multivariate_matched_sixth_order_weight.py)
extends the exact formal coefficient recurrence of the
[matched angular calculation](MULTIVARIATE_MATCHED_SECOND_ORDER_ANGLES.md#the-fourth-order-angular-model)
through degree six, initially at variance \(1\).
In the centered alternating moments the additional coefficient
is \(j/322560+j(j-1)/46080+j(j-1)(j-2)/82944\) in degree six
after the leading factor. It comes from the sixth-order term of
\((\sinh(s/2)/(s/2))^j\). The log expansion needs only the
power sums through degree seven. Higher log terms have scaling
degree at least seven because the degree in \(v\) at log order
\(a\) is at most \(a+1\). Root reversion uses exact rational
polynomial arithmetic in \(j,z\) and the signed and unsigned
moments. It checks symbolically that the sixth root coefficient
divided by \(z\) is even in \(z\) and has degree at most three in
\(j,z^2\). Its coefficient \(\ell^1\) norm is
\(57500100121/195955200<300\), and the degree in the moment
variables is at most six. The same checker tests the resulting
polynomial by 44 independent finite Gauss expansions and exact
quotient-ring root reversions.

For general \(R\), replace \(\xi_i\) by \(\xi_i/R\) and use the
unit-variance exponential coordinate with parameter
\(\lambda\sqrt R\) and height \(z/\sqrt R\).
The sixth coefficient divided by \(z\) is then \(R^3\) times
the unit-variance polynomial evaluated at \(j,z^2/R\).
Writing \(J=jR\) and \(Y=z^2\), it is homogeneous of degree
three in \(J,R,Y\), with the same coefficient norm and moment
degree. Denote this polynomial by \(E_j(Y)\).

On the common range \(y^2\le4N\), the physical variables
\(J,R,Y\) have moduli at most \(6n\) and the normalized moments
have moduli at most one. The estimate
\[
    \left|\partial_y\bigl(y E_j(-y^2)\bigr)\right|
                               \le453600 n^3
\]
uses the coefficient bound \(300\), the physical degree three
and the factor \(1+2b\le7\) for a term of degree \(b\) in \(Y\).
At fixed \(Y\), the physical variables change by
\(|\Delta J|\le X\), \(|\Delta R|\le1/n\).
Differentiating each monomial bounds this contribution to
\(|\Delta E|\) by \(32400(n+n^2X)\).
The normalized weights change by a relative amount at most
\(X/n\). Each signed or unsigned moment up to power seven
therefore changes by at most \(80X/n\), because
\(m(1+X/n)^{m-1}<80\) for \(m\le7,\ n\ge2,\ X\le1\).
A moment monomial changes by at most \(480X/n\).
Using the homogeneous physical bound \((6n)^3\) and
coefficient norm \(300\) gives
\[
                 |\Delta E|\le32{,}000{,}000 n(1+nX).
\]
For \(\chi=nX\le1/2000\), the slope bound for the sixth-order
model \(y(1+XB_j(y)+X^2 C_j(-y^2)+X^3 E_j(-y^2))\) is at least
\[
               1-\tfrac{10}{3}\chi-3600\chi^2-453600\chi^3
                           \ge 1/2 .
\]
The three correction differences at a common height have sum
bounded by
\[
 \left(\chi+25\chi^2+80000\chi^2(1+\chi)
                  +32{,}000{,}000\chi^3(1+\chi)\right)\frac{|y|}{n^2}
                         <\frac{|y|}{12n^2}.
\]
Thus the sixth-order model has the same lower bound of one sixth
of the original Hermite gap.

The exact branch remainder now starts with the eighth degree
in \(\lambda\). Under
\[
                        X\le\frac{10^{-16}}{16}N^{-19/6},
\]
one has \(\sqrt X/L\le\tfrac14 N^{-1/12}\).
The sum of the two possible angular errors is at most
\[
 \frac{2}{80\sqrt N}
       \frac{(N^{-1/12}/4)^8}{1-(N^{-1/12}/4)^2}
                       \le\frac{N^{-7/6}}{2457600},
\]
which is smaller than the lower bound for the sixth-order model
gap. The root count and simplicity argument is unchanged.
For arbitrary ordinary distinct exponent columns, it therefore
suffices that
\[
                             S\le 10^{-18}N^{-13/6}.
\]
The same conversion \(X\le S/(n(1-S))\) proves this last assertion.

Run the symbolic sixth-order check from the repository root:

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_matched_sixth_order_weight.py
~~~
