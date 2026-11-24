# Ordinary Strong Factorial for the support \(\{4,14\}\)

For every characteristic-zero field \(K\), every nonzero
\(f=X^4(A+BX^{10})\) in \(K[X]\), and every integer \(N\ge2\), the
ordinary factorial moments of orders \(N-1\) and \(N\) cannot both
vanish. The proof uses an explicit complex-contour estimate for
\(N\ge4000\) and exact checks at the intervening orders. The estimates
use the weighted two-packet argument of the
[lower-two theorem](ORDINARY_LOWER_TWO_STRONG_FACTORIAL.md) with
new inequalities for this support.

## 1. The two actual moments

It suffices to prove the complex assertion. With both coefficients
nonzero, put \(n=N-1\), \(w=(A/B)N^{-10}\), and
\(g(s)=s^4(w+s^{10})\). Substitution \(X=Ns\) in the factorial integral
writes the two moments as nonzero multiples of
\[
 I=\int_0^\infty g(s)^N e^{-Ns}\,ds,\qquad
 J=\int_0^\infty g(s)^n e^{-Ns}\,ds.
\]
The common scale is \(N\) in both integrals. Put
\[
 H=14^{14}e^{-14},\qquad L=4^4e^{-4},\qquad W=H/L.
\]
Split each integral at \(9\), using subscripts \(U,L\) for the intervals
\([9,\infty)\), \([0,9]\). The lower-packet subscript is distinct from the
constant \(L\).

Write \(\Phi(s)=\log g(s)-s\) where a local logarithm is defined.
The stationary equation is \(w(s-4)=s^{10}(14-s)\).
For \(|w|\le6W/5\), Rouché's theorem on \(|s-14|=1/4\) selects a unique
stationary point \(s_U\) inside the circle. Indeed,
\(\lvert w/s^{10}\rvert<1/100\) there and
\((1/100)(41/4)<1/4\).
For \(|w|\ge4W/5\), the same argument on \(|s-4|=1/16\) selects \(s_L\).
There \(\lvert s^{10}/w\rvert<1/1200\) and
\((1/1200)(161/16)<1/16\).
The zeros are simple by the counts with multiplicity. The first branch
extends to \(w=0\) with \(s_U=14\), and the second extends to \(1/w=0\)
with \(s_L=4\).

Set \(S_X=\Re\Phi(s_X)\) and \(A_X=-\Phi''(s_X)\) for \(X=U,L\).
At a stationary point
\[
 -\Phi''(s)=\frac{s^2-17s+56}{s^2}.
\]
Subtracting \(1/b\), for \(b=14,4\), factors the difference as
\((s-b)((b-1)s-56)/(bs^2)\).
The bounds on the two saddle disks give
\[
 \Re A_U>1/20,\quad |A_U|<1/10,\qquad
 \Re A_L>1/5,\quad |A_L|<3/10.
\]
To bound heights, integrate the stationary values from \(w=0\) for the
upper packet and from \(1/w=0\) for the lower packet after factoring out
\(w\) from \(g\). The stationary-point terms vanish under differentiation.
Along the radial paths, the bounds on the derivatives multiplied by the
radial parameter are \(1/99\) and \(1/1199\) at the respective ends and
scale at most linearly with the radial parameter. Equivalently the
derivative integrals are bounded by \(q/(1-q)\) with the respective
bounds \(q<1/100,1/1200\). This proves
\[
 |S_U-\log H|<1/80,\qquad |S_L-\log(|w|L)|<1/800 .
\]

## 2. Full packet estimates

Deform the upper ray to height \(\Im s_U\), including the vertical segment
at \(9\). Keep \([0,2]\) real in the lower packet and deform \([2,9]\)
to height \(\Im s_L\), including the two vertical ends. These are
deformations of entire integrands. The upper vertical at infinity
vanishes by exponential decay for each fixed \(N,w\).

On the central horizontal ranges \(12\le\Re s\le16\) and
\(3\le\Re s\le5\), respectively, we have
\[
 \Re\Phi''(s)<-1/30,\qquad \Re\Phi''(s)<-1/8 .
\]
For verification, use the representations of \(\Phi\) with leading
terms \(14\log s-s\) and \(4\log s-s+\log w\), and perturbation
\(\log(1+q(s))\) with \(q=w/s^{10}\) or \(q=s^{10}/w\).
Throughout the two strips, \(|q|\) is bounded by \(1/25,1/150\),
respectively. The inequality
\(\Re(1/s^2)\ge (1-3(\Im s/\Re s)^2)/(\Re s)^2\) gives at least
\(99/(100(\Re s)^2)\) in each strip. If \(y=|q|/(1-|q|)\), the
second derivative correction has modulus at most
\((110y+100y^2)/|s|^2\). Inserting the strip bounds gives the displayed
negative curvatures. Integration twice from the saddle bounds the
central moduli by
\(\exp(m S_X-m b_{\rm strip}v^2)\) for
\(v=\Re(s-s_X)\), with \(b_{\rm strip}=1/60,1/16\).

On the local horizontal intervals \(|v|\le\rho_X\) with
\(\rho_U=1/2,\rho_L=1/4\), the Taylor remainder after the quadratic
term has modulus at most \(B_X |v|^3\), where
\(B_U=1/250,B_L=1/25\). The derivative estimate used here is
\[
 |\Phi'''(s)|\le
 \frac{2b+1320y+3300y^2+2000y^3}{|s|^3}
 \qquad (b=14,4).
\]
It follows from the third derivative of \(\log(1+q)\) and applies on the
whole disks \(|s-14|\le3/4\) and \(|s-4|\le5/16\).
Together with the saddle curvatures, the estimates give local decay
\(\exp(-m b_Xv^2)\) relative to the saddle height, with
\(b_U=1/45,b_L=1/12\), for both the original integrand and its
quadratic approximation.

For the remaining pieces set
\[
 A(x)=(x/4)^4e^{4-x},\qquad B(x)=(x/14)^{14}e^{14-x}.
\]
The real upper envelope is \(B(x)+(6/5)A(x)\), relative to \(H\).
The real lower envelope is \(A(x)+(5/4)B(x)\), relative to \(|w|L\).
On \([9,10],[10,11],[11,12]\) and at \(16\), the first is at most
\(189/200\), by bounding \(A\) and \(B\) individually at their maximal
endpoints. For \(x\ge16\), both terms decrease relative to their values
at \(16\) by at least the factor \(e^{-(x-16)/8}\).
For the other packet, on \([0,3],[5,6],[6,7],[7,8],[8,9]\) the second
envelope is at most \(93/100\), by the same endpoint comparison.
The necessary monotonicity is that each term increases up to its
center and decreases thereafter.
For a complex point of real part \(x\), each envelope acquires at
most the factor \(\exp(7(\Im s)^2/x^2)\). On the nonreal pieces of
the upper and lower contours, respectively, this is bounded by
\(\exp(1/180)\) and \(\exp(1/140)\). The inequalities
\[
 (189/200)e^{1/180+1/80+1/30}<1,\qquad
 (93/100)e^{1/140+1/800+1/20}<1
\]
follow from \(e^z\le(1-z)^{-1}\) for \(0\le z<1\).
Thus all noncentral pieces lose \(1/30\) or \(1/20\), respectively,
relative to the selected saddle height, with the retained far-right
decay. This also bounds the vertical pieces and the origin interval.

The two leading terms in each packet are
\[
 I_X^0=e^{N\Phi(s_X)}\sqrt{\frac{2\pi}{NA_X}},\qquad
 J_X^0=e^{n\Phi(s_X)-s_X}\sqrt{\frac{2\pi}{nA_X}},
\]
with the same branch of the square root in the two terms. After
removing the constant amplitude in \(J_X^0\), the weight on the
translated contour is \(e^{-(s-s_X)}\).
For the constant-amplitude \(I\) use the same upper error bounds.
The following estimates give relative errors \(<1/4\) in all four
packets for \(N\ge4000\).

For either \(m=N\) or \(n\), take out \(e^{mS_X}\) and the appropriate
constant amplitude. On the local interval, interpolation between the
two exponents bounds their integrand difference by
\(m B_X |v|^3 e^{-m b_X v^2}\). The extra weight is bounded by
\(e^{\rho_X}\) and differs from one by at most \(e^{\rho_X}|v|\).
Integration over the real line bounds the local error by
\(e^{\rho_X}(B_X/b_X^2+1/b_X)/m\).
For the omitted Gaussian tails one uses
\[
 2\int_\rho^\infty e^{-mbv^2}\,dv
       \le \frac{e^{-mb\rho^2}}{mb\rho}.
\]
The off-local central tails use the same inequality with
\(b=b_{\rm strip}\). The weight on those ranges is bounded by \(10\)
and \(3\), respectively, using \(|\Re(s-s_U)|\le9/4\) and
\(|\Re(s-s_L)|\le17/16\). On the noncentral contours its modulus
is bounded by \(200\) and \(60\).
After allowing for the two finite lengths and the far-right integral
\(8/m\), the noncentral contributions are at most
\(800e^{-m/30}\) and \(600e^{-m/20}\).

The normalized leading moduli exceed \(7/\sqrt m\) and \(4/\sqrt m\).
Consequently the two relative errors are bounded by
\[
 \begin{aligned}
 E_U(m)&=
 \frac{(5/3)(81/10+45)+90e^{-m/180}+1200e^{-m/240}}{7\sqrt m}
          +\frac{800}{7}\sqrt m\,e^{-m/30},\\
 E_L(m)&=
 \frac{(4/3)(144/25+12)+48e^{-m/192}+192e^{-m/256}}{4\sqrt m}
          +\frac{600}{4}\sqrt m\,e^{-m/20}.
 \end{aligned}
\]
For \(m\ge3999\), the \(\sqrt m\,e^{-m/c}\) terms decrease and all
terms are \(<1/4\) in total for each packet. The
[exact checker](verify_lower_four_gap_ten.py) checks the numerical
inequalities in this section by rational arithmetic. It uses
\(163/60<e<87/32\) from the series for \(e\) and explicit rational
upper bounds for the error budgets.

## 3. Consecutive zero exclusion

For \(|w|<4W/5\), the real lower-packet integrand envelope for the
upper-order moment is \((4/5)A(x)+B(x)\), relative to \(H\).
Its modulus is \(<82/100\) throughout \([0,9]\), by the endpoint
bounds on \([0,5],[5,6],[6,7],[7,8],[8,9]\).
Combining this with the lower bound on \(|I_U|\) from Section 2 gives
\[
 |I_L/I_U| < 2\sqrt N\,(9/10)^N<1 .
\]
For \(|w|>6W/5\), the real upper envelope relative to \(|w|L\) is
\(A(x)+(5/6)B(x)\), at most \(9/10\) on \([9,12],[12,16]\) and at
\(16\), with the further factor \(e^{-(x-16)/8}\) thereafter.
The lower bound on \(|I_L|\) gives
\[
 |I_U/I_L| < 3\sqrt N\,(92/100)^N<1 .
\]
The last two comparisons hold for \(N\ge4000\) by monotonicity and
their value at \(4000\). Thus \(I\) has no zeros off the overlap.

In the overlap, the stationary values satisfy
\[
 |g(s_U)|>14^{14}(55/56)^{14}(99/100),\qquad
 |g(s_L)|<14^{14}(6/5)(65/64)^4(1201/1200)e^{-10}.
\]
The first bound exceeds four times the second. Each leading ratio
\(J_X^0/I_X^0\) equals \(\sqrt{N/n}/g(s_X)\). Each actual ratio
\(J_X/I_X\) has modulus between \(3/5\) and \(5/3\) times its leading
modulus by the error bounds. The two actual ratios therefore differ
in modulus. If both moments vanished, \(I_U=-I_L\) and \(J_U=-J_L\)
would instead give \(J_U/I_U=J_L/I_L\). This proves the analytic tail.

## 4. Finite completion and descent

The [terminal-through-twelve theorem](STEPANOV_AUXILIARY_DETERMINANTS_SFC2.md)
covers \(2\le N\le12\). The exact checker covers every
\(13\le N\le3999\) using the actual integer polynomials
\[
 P_m(t)=\sum_{k=0}^m\binom mk (4m+1)_{10k}t^k .
\]
Every factorial argument is at most \(55986\).
The checked primes \(1000003,1000033\) preserve both degrees,
and the checker verifies a full Bézout identity at a
degree-preserving prime for each pair. The
[replay receipt](lower_four_gap_ten_verification.json)
records the completed range. The replay command from the repository root
is `python3 research/factorial/image-rigidity/verify_lower_four_gap_ten.py`
with Python 3.12.3 and python-flint 0.9.0.

For \(A\ne0\), the ordinary moment is \(A^m(4m)!P_m(B/A)\).
The complex proof implies rational coprimality in the analytic tail,
and the finite witnesses imply rational coprimality in the checked range.
Rational Bézout transfers all these pairs to every characteristic-zero
field. If either coefficient is zero, the nonzero monomial has every
positive factorial moment nonzero. This proves the stated support
theorem.
