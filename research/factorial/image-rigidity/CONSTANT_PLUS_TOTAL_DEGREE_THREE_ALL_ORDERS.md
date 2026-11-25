# Strong Factorial for a constant and a monomial of total degree three

For \(d=(1,1,1)\), \((2,1)\), or \((3)\), set
\[
 R_n(z)=\Lambda((z+X_1^{d_1}\cdots X_r^{d_r})^n)
       =\sum_{k=0}^n\binom nk\prod_{i=1}^r(d_i k)!\,z^{n-k}.
\]
Then \(R_n\) is squarefree for every \(n\ge1\).
The identity \(R_n'=nR_{n-1}\) gives consecutive-moment
separation at every terminal.
Coordinates with equal exponents in both monomials contribute
a nonzero scalar to each moment and can be removed.

The [exact replay](verify_constant_plus_total_degree_three.py)
checks both multivariate prefixes through terminal \(99{,}999\).
A quantitative contour argument covers \(n\ge100{,}000\) for
those profiles, and a \(3\)-adic argument covers \((3)\) at
every order.

## The profile \((3)\)

The following is the \(d=3\) instance of the
[univariate prime-power argument](../../../unpublished/strong-factorial-binomials/constant-support.tex).
Scale the reversed polynomial by \(t\mapsto t/3\) to obtain
\[
 \sum_{k=0}^n \frac{n!}{(n-k)!}\,U_k t^k,\qquad
 U_k=\frac{(3k)!}{3^k k!}
     =\prod_{\substack{1\le j\le3k\\3\nmid j}} j
     \equiv(-1)^k\pmod3.
\]
Write \(n=\sum_i a_i3^i\) in base three and
\(s_i=\sum_{j\le i}a_j3^j\), with \(s_{-1}=0\).
The coefficient ordinate is \(v_3(n!)-v_3((n-k)!)\).
In the interval \([s_{i-1},s_i]\) with \(a_i>0\), let
\(k=s_i-\Delta\). The difference between the ordinate and
the chord through the two endpoints is
\((\operatorname{digitsum}_3(\Delta)-\Delta/3^i)/2\),
by Legendre's formula and the disjoint digit positions of
\(n-s_i\) and \(\Delta\). This is nonnegative, with
equality exactly at multiples of \(3^i\) away from \(s_i\).
The chord has slope \(h_i/e_i\), where
\(e_i=3^i,\ h_i=(3^i-1)/2\). These slopes are strictly
increasing and in lowest terms, so the chords are the
edges of the lower Newton polygon.

The unit of \(l!\) after removing its powers of three
is congruent to \((-1)^{v_3(l!)}\) times the product of
the factorials of its base-three digits. This follows by
recursively removing the multiples of three and multiplying
the remaining residues in each block. Along the edge,
only the digit \(a_i\) in \((n-k)!\) changes. Hence after
dividing the edge coefficients by their powers of three
and normalizing the constant term, the residual polynomial
is
\[
 \sum_{j=0}^{a_i}\frac{a_i!}{(a_i-j)!}(\eta_iY)^j,
 \qquad \eta_i=(-1)^{h_i+e_i}.
\]
It is a scaled reversal of \(\sum_{j=0}^{a_i}Z^j/j!\).
The latter is squarefree in \(\mathbb F_3[Z]\), since
\(a_i<3\) and its difference with its derivative is
\(Z^{a_i}/a_i!\).

For clarity, at a fixed slope \(u/e\) take \(\pi^e=3\)
in a valued extension. Substitute \(t=\pi^{-u}Z\),
divide by a power of \(\pi\) to give minimum coefficient
valuation zero, and reduce. This reduction multiplies
under polynomial products. The surviving degrees are
congruent modulo \(e\), so removing the lowest monomial
gives the same identity of residual polynomials in
\(Y=Z^e\). A repeated factor would therefore give a
repeated residual factor at any one of its slopes.
The squarefree residuals just proved exclude this.

## Finite prefix from nearby primes

The coefficient of \(t^k\) in \(t^n R_n(1/t)\) is
\(n!(k!)^2/(n-k)!\) for \(d=(1,1,1)\), and
\(n!(2k)!/(n-k)!\) for \(d=(2,1)\).
For \(d=(1,1,1)\) a prime terminal \(n=p\) is immediate:
\(R_n\equiv z^n\pmod p\) and \(R_{n-1}(0)\not\equiv0\pmod p\).

Otherwise write \(n=p+a\) with \(p\) prime and \(a>0\).
For \(d=(1,1,1)\) assume \(2a<p-2\). The valuations
of the displayed coefficients are zero on \(0\le k\le a\),
one on \(a<k<p\), and three on \(p\le k\le n\),
since \(n<2p\). The positive Newton slopes of the reverse
polynomial are \(1/(p-1-a)\) and \(2/(a+1)\),
from the vertices \((a,0),(p-1,1),(n,3)\).
At order \(n-1\) replace \(a\) by \(a-1\).
The two positive-slope sets are disjoint, using \(2a<p-2\).

For \(d=(2,1)\) assume \(3a<p-3\).
Here \(n<4p/3\), \(2n<p^2\), and \(p\) is odd.
Writing \(m=(p-1)/2\), the valuation is zero up to \(k=a\),
then one through \(k=m\), two through \(k=p-1\), and
three through \(k=n\), by Legendre's formula.
The intermediate point \((m,1)\) lies on or above the
segment from \((a,0)\) to \((p-1,2)\).
The positive slopes are \(2/(p-1-a)\) and \(1/(a+1)\),
in this order by \(3a<p-3\). Again the sets at \(n\) and
\(n-1\) are disjoint.

A root of a monic integral common divisor of the
reciprocals would have valuation zero at \(p\), by these
Newton polygons. Such a divisor would retain its degree
and have nonzero constant term on reduction. At both
profiles the initial coefficient recurrences give
\[
 R_n(z)\equiv z^p R_a(z),\qquad
 R_{n-1}(z)\equiv z^p R_{a-1}(z)\pmod p .
\]
Indeed all terms at \(k>a\) have positive valuation at
order \(n\), and for \(k\le a\) the falling factorial
\(n!/(n-k)!\) may be replaced by \(a!/(a-k)!\) modulo \(p\).
The same argument applies at order \(n-1\) with \(a-1\).
Thus a constant modular gcd at the smaller orders proves
characteristic-zero coprimality.

For each terminal the replay tries the four nearest primes
at most \(n\). For \(d=(1,1,1)\) there are \(9{,}592\)
prime terminals and \(90{,}403\) near-prime certificates,
with three additional full modular checks. For \(d=(2,1)\)
there are \(99{,}987\) near-prime certificates and eleven
full modular checks. Each total accounts for all \(99{,}998\)
terminals from \(2\) to \(99{,}999\). Order one is immediate.

## Radial density with an explicit error

Take independent unit-rate exponential random variables \(X_i\).
Put \(B_0=\prod_i d_i^{d_i}\) and
\(T=(\prod_i X_i^{d_i}/B_0)^{1/3}\).
Then \(R_n(z)=\mathbb E(z+B_0T^3)^n\).
Write \(X_i=d_i t e^{u_i}\) with \(\sum_i d_i u_i=0\),
and use \(x_i=\sqrt{d_i}\,u_i\) on the Euclidean hyperplane
\(W=\{x:\sum_i\sqrt{d_i}x_i=0\}\).
The \(t\) column of the change to
\(\sqrt{d_i}\log(X_i/d_i)\) is orthogonal to \(W\) with
length \(\sqrt3/t\). The density of \(T\) is therefore
\[
 f(t)=C_0 t^m e^{-3t}\int_W e^{\ell(x)-tG(x)}\,dx,
 \quad m=r-1,\quad C_0=\sqrt3\prod_i\sqrt{d_i},
 \quad
 \ell=\sum_i u_i,\quad G=\sum_i d_i(e^{u_i}-1-u_i).
 \tag{1}
\]
It is holomorphic for \(\Re t>0\) since the unsubtracted
phase grows exponentially on \(W\). Indeed the largest
\(u_i\) is at least a positive constant times \(\|x\|\)
on the weighted-zero hyperplane. For \(d=(1,1,1)\)
one has \(\ell=0\), and for \(d=(2,1)\) its norm as a
functional on \(W\) is \(1/\sqrt6\), by orthogonal projection.
In both cases \(|\ell(x)|\le\|x\|\).

Let \(\tau=\Re t\ge50{,}000\) and \(|t|/\tau\le51/50\).
We prove
\[
 f(t)=C_2 t^{m/2} e^{-3t}(1+\epsilon(t)),\quad
 C_2=C_0(2\pi)^{m/2},\qquad |\epsilon(t)|<1/40 .
 \tag{2}
\]
For \(\rho=\|x\|\), one has
\(G(x)\ge\rho^2/[2(1+\rho)]\).
For negative \(u=-v\), the scalar bound follows by
integrating \(1-e^{-v}\ge v/(1+v)\). For nonnegative \(u\),
use the exponential series. Sum over coordinates using
\(|u_i|\le\rho\).
For \(\rho\le1/10\), Taylor's formula bounds
\(\left|G(x)-\rho^2/2\right|\le(3/16)\rho^3\),
since \(e^{1/10}<9/8\).

The absolute difference between the integrand in (1) and
\(e^{-t\rho^2/2}\) on this small ball is bounded by
\[
 (9/8)(\rho+\tau\rho^3/5)\,e^{-(12/25)\tau\rho^2}.
\]
This uses \(|e^z-1|\le |z|e^{|z|}\) and
\(|t(G-\rho^2/2)|\le\tau\rho^3/5\le\tau\rho^2/50\).
Relative to the modulus \((2\pi/|t|)^{m/2}\) of the
complex Gaussian integral, the bound integrates to
less than \(4/\sqrt\tau\). In detail, the ratio of the
Gaussian normalization constants, including \(9/8\),
is at most \((9/8)(51/48)^{m/2}<6/5\).
For the real Gaussian with coordinate variance
\(25/(24\tau)\) and \(m\le2\), its first radial moment
is less than \((3/2)/\sqrt\tau\), and its third is at most
the \(3/4\) power of its fourth radial moment, hence less
than \((11/2)/\tau^{3/2}\).

For \(1/10\le\rho\le1\), the true integrand is at most
\(e^{-\tau\rho^2/5}\), using \(G\ge\rho^2/4\) and
\(\rho\le\tau\rho^2/20\). Splitting the Gaussian bound
into two equal factors gives relative integral at most
\(6e^{-\tau/1000}\).
For \(\rho\ge1\), use \(G\ge\rho/4\) to get
\(e^{-\tau/10}e^{-\tau\rho/10}\).
In orthonormal coordinates
\(\rho\ge\sum_j|x_j|/\sqrt m\), so the integral here
is at most \(e^{-\tau/10}(20\sqrt m/\tau)^m\).
Dividing by \((2\pi/|t|)^{m/2}\) leaves at most
\(e^{-\tau/10}\).
The omitted part of the comparison Gaussian outside
\(\rho\le1/10\) similarly gives relative error at most
\(3 e^{-\tau/400}\).
Altogether \(4/\sqrt\tau+10e^{-\tau/1000}<1/40\).
Evaluating the full complex Gaussian then proves (2).

## Scaling, root radius and contour geometry

In the remainder let \(n\ge100{,}000\), and put
\[
 A=(3/e)^3,\qquad w=Av,\qquad z=B_0(n/e)^3v,\qquad L=3/2.
\]
The change \(s=3T/n\) shows that the scaled polynomial
is a nonzero scalar times
\[
 I_n(v)=K_n(v)+J_n(v)
 =\int_0^L(w+s^3)^n(n/3) f(ns/3)\,ds
 +\int_L^\infty(w+s^3)^n(n/3)f(ns/3)\,ds .
 \tag{3}
\]
Both parts are polynomials in \(v\).

The roots have \(|v|<7/6\). For the scaled monic polynomial
let \(b_k\) be the coefficient of \(v^{n-k}\).
When \(\theta=k/n\le1/20\), use \(\binom nk\le(en/k)^k\)
and \((d_i k)!\le(d_i k)^{d_i k}\) to obtain
\(b_k^{1/k}\le e^4\theta^2<3/20\).
For \(k\ge n/20\), Stirling and the binomial entropy
bound give
\[
 b_k^{1/k}\le q(\theta)(2ek)^{3/k},\qquad
 q(\theta)=\theta^2(1-\theta)^{-(1-\theta)/\theta}.
\]
Here \(l!\le el(l/e)^l\) and \(\prod_i ed_i k\le(2ek)^3\).
The unique maximum of \(q\) is at
\(3\theta+\log(1-\theta)=0\), by differentiation.
The solution exceeds \(9/10\), and at it
\(\log q=2\log(1-u)+3u<u<1/10\) with \(u=1-\theta\).
Thus \(q<10/9\). Also
\((2ek)^{3/k}\le e^{60\log(2en)/n}<101/100\),
since the exponent decreases and \(\log(100{,}000)<12\).
After division by \(7/6\) the \(k\)-th-root bound is
less than \(97/100\). The first range contributes a
geometric sum less than \(1/4\), and the second a
geometric tail starting at \(k\ge5{,}000\) less than
\(1/4\). Hence the leading term dominates on \(|v|\ge7/6\).

On \(|v|\le6/5\), the saddle \(s_0\) solves
\(w=s_0^2(3-s_0)\) and \(|s_0-3|<\delta=9/(2e^3)<9/40\).
For existence and uniqueness, on that circle divide the
modulus of \(s^2(3-s)\) by \(A\) and use
\((3/2)(1-\delta/3)^2>6/5\) in Rouché's theorem.
The saddle is simple and analytic in the disk.
Set \(\Psi=3s_0^2 e^{-s_0}\), \(H=|\Psi|\), and
\(\kappa=1-2/s_0\). The inequality
\(|3\kappa-1|=|2(s_0-3)/s_0|<6/37\) gives
\(\Re\kappa>4/15\) and \(|\kappa|<2/5\).
Also \(H/A>9/10\) by expanding
\(2\log(s_0/3)+3-s_0\): its linear term has modulus below
\(3/40\), and the remaining terms have modulus at most
\((3/40)^2/(37/40)\).

Put \(g(s)=e^{-\Re s}|w+s^3|\) and
\(\Phi(s)=\log(w+s^3)-s\) where the logarithm is used
locally. On \(s=x+i\Im s_0\), \(x\ge L\), factoring
the derivative gives
\[
 \Phi'(s)=\frac{\Re s_0-x}{s}
  \frac{1+(s_0-3)s^{-1}(1+s_0/s)}{1+w/s^3}.
\]
Both the numerator and denominator of the fraction
differ from one by less than \(1/2\), using
\(|s|\ge3/2,\ |s_0|<129/40,\ |w|<81/50\).
Their absolute arguments are below \(\pi/6\),
and \(|\arg s|<3/20<\pi/6\).
The real part of the factor multiplying \(\Re s_0-x\)
is positive. Thus the strict maximum on the shifted
ray is at \(s_0\).
For the connector \(s=L+iy,\ |y|\le\delta\), use
\[
 \frac{g(s)}{A}\le \frac65 e^{-3/2}
               +(e^{1/2}/2)^3 e^{2\delta^2/3}
 <\frac3{11}+\frac{(5/6)^3}{1-27/800}<\frac{22}{25}.
 \tag{4}
\]
The scalar bounds here follow from \(e^3>20\),
\(\sqrt{20}>22/5\), and \(e^{1/2}<5/3\).
In particular \(g/H<44/45\) on the connector.

## Explicit saddle error and the tail derivative

Deform \(J_n\) to the connector and the shifted ray.
The density is holomorphic there, and (2) controls the
closing segment as its real part tends to infinity.
On the path \(t=ns/3\) satisfies \(\Re t\ge n/2\)
and \(|t|/\Re t\le51/50\).

For \(|s-s_0|\le1/4\), \(|s|\ge101/40\).
Write \(\Phi(s)=3\log s+\log(1+G_1(s))-s\)
with \(G_1(s)=w/s^3\).
One has
\[
 |G_1|<51/500,\qquad
 |G_1'|\le(6/5)|G_1|,\qquad
 |G_1''|\le(19/10)|G_1|,\qquad
 |G_1'''|\le(15/4)|G_1|.
\]
Taking three derivatives of the logarithms then gives
\(|\Phi'''|<9/10\). The stated rational bounds imply
this directly by the formula
\((\log(1+G_1))'''=G_1'''/(1+G_1)
 -3G_1'G_1''/(1+G_1)^2+2(G_1')^3/(1+G_1)^3\).
Consequently for real \(|h|\le l=1/4\),
\[
 \Phi(s_0+h)-\Phi(s_0)=-\kappa h^2/2+E(h),\quad
 |E(h)|\le(3/20)|h|^3,\quad
 g(s_0+h)/H\le e^{-h^2/12}.
 \tag{5}
\]
Here the last inequality uses
\(\Re\kappa/2>2/15\) and \(|E(h)|\le(3/80)h^2\).
The amplitude ratio \((1+h/s_0)^{m/2}\)
differs from one by at most \(|h|\), by differentiating
along the segment, and its modulus is below two.

The complex Gaussian main term is
\[
 J_0(v)=C_2(n/3)^{(r+1)/2}s_0^{m/2}\Psi(v)^n
                \sqrt{\,2\pi/(n\kappa)\,}.
\]
It is analytic and nonzero on the disk.
Normalize by the amplitude and \(\Psi(v)^n\).
Relative to the full Gaussian modulus, which is at
least \(\sqrt{5\pi/n}\), the local phase and amplitude
errors are less than \(6/\sqrt n\) and \(4/\sqrt n\).
These bounds follow by integrating respectively
\((3/20)n|h|^3e^{-nh^2/12}\) and \(|h|e^{-nh^2/12}\)
on the real line. The relative density error is below
\((2/40)\sqrt{12/5}<2/25\), by (2).

Outside the local segment the shifted ray has
\(g/H\le e^{-1/192}\), by the strict monotonicity
and (5) at \(h=\pm l\). Equation (4) gives the same
bound on the connector. On the ray up to \(x=6\)
the normalized amplitude modulus is below three and
the density ratio below two. This gives contribution
at most \(27e^{-n/192}\). The connector contributes
at most \(e^{-n/192}\), since its length is below \(9/40\)
and its amplitude ratio is at most one.
For \(x\ge6\),
\[
 \Re\Phi'(s)\le\frac{3}{|s|(1-|w|/|s|^3)}-1<-2/5,
\]
and the normalized amplitude is bounded by \(x+1\).
Integration bounds that part by
\((35/n+25/(2n^2))e^{-n/192}\).
The omitted Gaussian tail is at most
\((30/n)e^{-n/120}\).
Thus the total relative estimate is
\[
 \left|\frac{J_n}{J_0}-1\right|
 <\frac{10}{\sqrt n}+\frac2{25}+10\sqrt n\,e^{-n/192}<\frac18 .
 \tag{6}
\]
All decreasing terms can be evaluated at \(100{,}000\),
using \(\sqrt{100{,}000}<317\) and \(e^{100{,}000/192}>2^{20}\).

Use Cauchy's estimate with radius \(1/40\) at
\(|v|\le7/6\) to bound the logarithmic derivative
contributed by the ratio in (6) by \(40/7\).
Also \(|s_0'|=|A/[3s_0(2-s_0)]|<21/100\).
The logarithmic derivative of \(s_0^{m/2}\kappa^{-1/2}\)
has modulus below one, since
\(\kappa'=2s_0'/s_0^2\) and \(|\kappa|>4/15\).
Together with
\(\Psi'/\Psi=A/(3s_0^2)=(3-s_0)/(3v)\), this gives
\[
 \left|\frac{J_n'}{J_n}-\frac{n(3-s_0)}{3v}\right|<7
 \quad (|v|\le7/6).
 \tag{7}
\]
The first expression for \(\Psi'/\Psi\) is regular at zero.

## Endpoint integral and simple roots

On every real subinterval of \([0,L]\), \(g\) attains
its maximum at an endpoint. In fact the derivative of
\(g(x)^2\) has the sign of \(Q(x)-|w|^2\), where
\[
 Q(x)=x^5(3-x)+a_w x^2(3-2x),\quad
 Q'(x)=3x\{x^3(5-2x)+a_w(2-2x)\},\quad a_w=\Re w.
\]
When \(a_w\ge0\) the bracket is positive for \(0<x\le1\),
and for \(1\le x\le3/2\) its first term is at least two
whereas the second has magnitude below \(81/50\).
If \(a_w<0\), the bracket is strictly increasing up to
one and positive beyond one. Thus \(Q\) first decreases
and then increases from \(Q(0)=0\), or is increasing.
In either case \(g^2\) has at most an interior minimum.

Let \(S=\sum_i X_i\), so \(3T\le S\) by weighted
AM–GM. The bound
\[
 W_n=\mathbb E[e^{3T}\mathbf1_{\{3T\le nL\}}]
   \le (3n)^3+64
 \tag{8}
\]
follows by dividing into \(S\le2nL\) and \(S>2nL\).
On the first event integration of \(e^S\) against the
joint exponential density gives the simplex volume
\((2nL)^r/r!\). On the second event use \(e^{3T}\le e^{nL}\)
and \(\mathbb E e^{3S/4}=4^r\).

If \(|w|\le(199/200)H\), the endpoint comparison and (4)
give \(|K_n|\le W_n((199/200)H)^n\). This is less than
\(|J_n|\) by (6) and \(|J_0|\ge H^n\).
Thus any root lies at \(|v|<7/6,\ |w|>(199/200)H\).

For \(|w|\ge(99/100)H\), one has \(|w|>4/5\), using
\(H/A>9/10\) and \(A>1\).
On \(s\le1/10\) the ratio \(s^3/|w|\) is at most \(s/4\).
At \(s=1/10\) this gives \(g(s)/|w|\le e^{-3/40}<99/100\).
By (4) the same bound holds at \(L\).
Therefore it holds on \([1/10,L]\) by the endpoint
comparison, and this region contributes at most
\(W_n e^{-n/100}\) to the modulus of \(K_n/w^n\).

For \(s=3T/n\le1/10\), the normalized multiplier differs
from one by at most
\[
 n(s^3/|w|)e^{ns^3/|w|}
 \le (5/4)S^3 n^{-2}e^{S/4}.
\]
Since \(\mathbb E[S^3e^{S/4}]
 =r(r+1)(r+2)(4/3)^{r+3}<340\), its integral is
less than \(425/n^2\). The missing probability for
\(s>1/10\) is at most \(8e^{-n/20}\) by
\(\mathbb E e^{S/2}\le8\). Combining the three bounds,
\[
 \left|\frac{K_n}{w^n}-1\right|
 \le W_n e^{-n/100}+\frac{425}{n^2}+8e^{-n/20}
 <\frac1{1000}.
 \tag{9}
\]
The polynomial factors in (8) are dominated by the
displayed exponentials throughout \(n\ge100{,}000\).

On \(|v|\le6/5\) one has \(|(w/\Psi)'|<2\):
use \(A/H<10/9,\ |w|/H<4/3\), and
\(|\Psi'/\Psi|=A/(3|s_0|^2)<3/50\).
Every disk of radius \(1/1000\) around a root therefore
stays inside the region of (9), since the root satisfies
\(|w|/H>199/200\). Cauchy's estimate gives
\(\left|K_n'/K_n-n/v\right|<2\) at the root.
There \(J_n=-K_n\ne0\), and (7) gives
\[
 \left|\frac{J_n'}{J_n}-\frac{K_n'}{K_n}
                +\frac{ns_0}{3v}\right|<9.
\]
The main term \(ns_0/(3v)\) has modulus exceeding
\(3n/4\) since \(|s_0|>111/40\) and \(|v|<7/6\).
Thus \(I_n'\) is nonzero at every zero of \(I_n\).
The scaled polynomial and \(R_n\) have the same
multiplicities. This closes the analytic range and hence
all terminal orders for both profiles. Their rational
Bézout identities transfer the result to every
characteristic-zero field.
