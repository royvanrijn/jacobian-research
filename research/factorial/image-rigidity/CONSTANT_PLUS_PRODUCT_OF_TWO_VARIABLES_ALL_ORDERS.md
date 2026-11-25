# Strong Factorial for a constant and the product of two variables

For the ordinary factorial functional \(\Lambda(X^iY^j)=i!j!\), put
\[
 R_n(z)=\Lambda((z+XY)^n)
       =\sum_{k=0}^n\binom nk(k!)^2z^{n-k}.
\]
Then \(R_n\) is squarefree for every \(n\ge1\).
In particular \(R_n'=nR_{n-1}\) gives coprimality at every
terminal. The characteristic-zero statement follows from the
rational Bézout identities. Equal-exponent rows may be removed
from the original support, since they only multiply each moment
by a nonzero factorial scalar.

The proof uses exact arithmetic up to terminal \(99{,}999\)
and a quantitative complex saddle estimate from terminal
\(100{,}000\). The
[verification program](verify_constant_plus_product_of_two_variables.py)
replays the finite prefix and the rational contour comparison
using only the Python standard library.

## Exact checks using a prime near the terminal

At a prime terminal \(n=p\), \(R_n\equiv z^n\pmod p\) and
\(R_{n-1}(0)=((p-1)!)^2\not\equiv0\pmod p\).
These two monic reductions are coprime.

Suppose next that \(n=p+a\) for a prime \(p\) and
\(0<2a<p-2\). The coefficient of \(t^k\) in the reverse
polynomial \(t^nR_n(1/t)\) is \(n!k!/(n-k)!\).
Since \(n<2p<p^2\), its \(p\)-adic valuation is zero for
\(k\le a\), one for \(a<k<p\), and two for \(p\le k\le n\).
The lower Newton hull of the reverse polynomial therefore
has its edges through
\[
 (0,0),\quad(a,0),\quad(p-1,1),\quad(n,2).
\]
The two positive slopes are \(1/(p-1-a)\) and \(1/(a+1)\),
in increasing order by \(2a<p-2\). At order \(n-1\),
replace \(a\) by \(a-1\); the slopes are \(1/(p-a)\) and
\(1/a\). No positive slope occurs at both orders. A common
root of the reciprocals must therefore have \(p\)-adic
valuation zero. Any nonconstant monic rational common divisor
would be integral by Gauss's lemma and would have unit constant
term at \(p\). Its reduction would retain its degree and have
nonzero constant term. But
\[
 R_n(z)\equiv z^p R_a(z),\qquad
 R_{n-1}(z)\equiv z^p R_{a-1}(z)\pmod p,
\]
by the same valuation computation and
\(n!/(n-k)!\equiv a!/(a-k)!\pmod p\) for \(k\le a\).
Consequently checking the degree-\(a\) modular pair
\(R_a,R_{a-1}\) suffices.

At \(2\le n\le99{,}999\) the verifier applies this argument
with the largest prime at most \(n\). It finds \(9{,}592\)
prime terminals and \(90{,}403\) successful tests using a
strictly smaller prime. At \(n=50{,}533\) it uses the
alternative prime \(50{,}513\), so \(a=20\). At \(n=4,10\)
it checks the full reciprocal pairs modulo \(5,11\),
respectively. This accounts for all \(99{,}998\) terminals
in that interval. Order one is immediate.

## Density, scaling and the outer root radius

Let \(X,Y\) be independent unit-rate exponential random
variables and \(T=\sqrt{XY}\). The change of variables
\(X=t e^u,\ Y=t e^{-u}\) has determinant of absolute value
\(2t\). The density of \(T\) is
\[
 f(t)=4t\int_0^\infty e^{-2t\cosh u}\,du
     =4t e^{-2t}\int_0^\infty
         e^{-t q^2}(1+q^2/4)^{-1/2}\,dq .
\]
In the last identity take \(q=2\sinh(u/2)\).
The formula extends holomorphically to \(\Re t>0\).
The Gaussian integral and
\(0\le1-(1+q^2/4)^{-1/2}\le q^2/8\) give
\[
 f(t)=2\sqrt\pi\, t^{1/2}e^{-2t}(1+\epsilon(t)),\qquad
 |\epsilon(t)|\le\frac{|t|^{1/2}}{16(\Re t)^{3/2}}.
 \tag{1}
\]
The power uses the right-half-plane branch. For \(t>0\)
the main term itself is an upper bound for \(f(t)\).

In the rest of the proof let \(n\ge100{,}000\). Write
\[
 A=4/e^2,\qquad w=Av,\qquad z=(n/e)^2v.
\]
Changing variables \(s=2t/n\) in \(R_n(z)=\mathbb E(z+T^2)^n\)
gives, up to a nonzero scalar,
\[
 I_n(v)=K_n(v)+J_n(v)
       =\int_0^1 (w+s^2)^n \frac n2 f(ns/2)\,ds
         +\int_1^\infty (w+s^2)^n \frac n2 f(ns/2)\,ds .
 \tag{2}
\]
Both integrals are polynomials in \(v\).

All roots have \(|v|<5/4\). In fact the monic polynomial
\((e/n)^{2n}R_n((n/e)^2v)\) has coefficient
\(b_k=\binom nk(k!)^2(e/n)^{2k}\) at \(v^{n-k}\).
For \(1\le k\le n/20\), \(k!\le k^k\) gives
\(b_k^{1/k}\le e^2 k/n\), so
\(b_k/(5/4)^k\le(1/3)^k\).
For \(k\ge n/20\), the binomial entropy bound and
\(k!\le ek(k/e)^k\) give
\[
 b_k^{1/k}\le
 q(k/n)(ek)^{2/k},\qquad
 q(\theta)=\theta(1-\theta)^{-(1-\theta)/\theta}.
\]
The maximum of \(q\) occurs at the interior solution to
\(2\theta+\log(1-\theta)=0\), by differentiating its logarithm.
That solution exceeds \(79/100\) since \(e^{79/50}>100/21\).
At the maximum write \(u=1-\theta\). Then
\(\log q=\log(1-u)+2u<u<21/100\), so \(q<1.234\).
Also \((ek)^{2/k}\le\exp(40(1+\log n)/n)<1.006\) in
this range, since the exponent decreases and
\(\log(100{,}000)<12\). Thus \(b_k/(5/4)^k<(199/200)^k\).
The total sum of the first bounds is at most \(1/2\).
The total sum of the second bounds is at most
\(200 e^{-25}<1/100\). Hence the leading term dominates
all other terms when \(|v|\ge5/4\).
The exponential comparisons here use the series lower
bound through degree six for \(e^{79/50}\), and
\(e^x\le1+x+x^2/[2(1-x/3)]\) for the upper bounds
at \(x=21/100\) and \(x=52/10000\).

## Saddle and contour geometry

Use the disk \(|v|\le13/10\), so \(|w|<c=71/100\) by
\(e>163/60\). Define
\[
 s_0=1+\sqrt{1-w}=1+u+ib,\qquad
 \kappa=1-1/s_0,\qquad \Psi(v)=2s_0e^{-s_0},\qquad
 H=|\Psi(v)|.
\]
The square root has positive real part. One has
\(53/100<u<131/100\) and \(|b|\le c/(2u)<67/100\),
using \(1-c\le u^2\le1+c\) and \(2ub=-\Im w\).
Also \(|\sqrt{1-w}|\le\sqrt{1+c}\), so
\(\Re\kappa>1/3\) and \(|\kappa|=|\sqrt{1-w}|/|s_0|<1\).
Set \(g(s)=e^{-\Re s}|w+s^2|\). Since \(w+s_0^2=2s_0\),
the saddle modulus is \(g(s_0)=H\).

For \(s=x+ib,\ x\ge1\), the logarithmic phase
\(\Phi(s)=\log(w+s^2)-s\) has derivative
\[
 \Phi'(x+ib)=(1+u-x)
 \frac{x+u-1+2ib}{x^2+1-u^2+2ib(x-u)}.
\]
The real part of the denominator is at least \(1-c\).
The real part of the quotient is positive, since multiplying
its numerator by the conjugate of its denominator gives real
part
\[
 (x+u-1)(x^2+1-u^2)+4b^2(x-u).
\]
For \(x\ge u\) this is positive. For \(x<u\), one has \(u>1\),
so the first term is at least \(1-c\), whereas the negative
term has magnitude below \(c^2(31/100)<1-c\).
Thus \(g(x+ib)\) increases strictly up to \(x=1+u\)
and decreases strictly thereafter. For \(x\ge5\) its
logarithmic derivative is at most \(-1/2\), by
\[
 \Re\Phi'(s)\le\frac{2}{|s|(1-c/|s|^2)}-1<-\tfrac12 .
\]

On \(s=1+iy\), the squared modulus of \(w+s^2\) is
\((2-u^2+b^2-y^2)^2+4(y-ub)^2\).
Its second derivative is \(12y^2+8-4(1+\Re w)>0\).
Thus the maximum on the vertical connector is at an
endpoint. At the real endpoint \(y=0\), writing \(t=b^2\),
\[
 \frac{g(1)^2}{H^2}
 =\frac{e^{2u}((2-u^2)^2+(4+2u^2)t+t^2)}
         {4((1+u)^2+t)},\qquad
 (1-u^2)^2+2(1+u^2)t+t^2\le c^2 .
\]
The verifier proves \(g(1)^2/H^2<24/25\) on
\(53/100\le u\le131/100\). For each rational subinterval
\([a,h]\) of length \(1/1000\), it uses
\(\Delta=\max(0,1-h^2,a^2-1)\le|1-u^2|\).
If \(\Delta>c\) the interval does not meet the constraint.
Otherwise \(t\le(c^2-\Delta^2)/(2(1+a^2))\).
The verifier bounds \((2-u^2)^2\) by its decreasing
endpoint value at \(a\) and \(4+2u^2\) by its increasing
endpoint value at \(h\). It uses the smaller denominator
\(4(1+a)^2\) and bounds \(e^{2u}\le e^{2h}\) by the
Taylor sum through degree sixteen plus a geometric tail. It checks \(770\) subintervals.

## Quantitative tail integral and its logarithmic derivative

In (2) deform \(J_n\) to the connector from \(1\) to \(1+ib\)
and the shifted ray from \(1+ib\) to infinity. Equation (1)
bounds the far closing vertical integral by a polynomial
in its real part times a decaying exponential, which
justifies the deformation. On the deformed path
\(\Re s\ge1,\ |\Im s|\le67/100\), so (1) gives
\(|\epsilon(ns/2)|\le1/n\).

For \(|h|\le1/4\), the phase has the local expansion
\[
 \Phi(s_0+h)-\Phi(s_0)
  =\log(1+h+h^2/(2s_0))-h
  =-\frac{\kappa h^2}{2}+E(h),\qquad |E(h)|\le |h|^3 .
 \tag{3}
\]
Indeed \(1/(2|s_0|)<1/3\). For \(q=h+h^2/(2s_0)\)
use \(|q|\le(1+|h|/3)|h|\) and the bound for the log-series
remainder beyond degree two,
\(|q|^3/[3(1-|q|)]\). At \(|h|=1/4\) the resulting
upper bound on \(|E(h)|/|h|^3\) is
\[
 \frac13+\frac1{72}
 +\frac{(13/12)^3}{3(1-13/48)}<1 .
\]
For real \(h\) with \(|h|\le l=1/24\), (3) gives
\(g(s_0+h)/H\le e^{-h^2/8}\) and
\(\left|(1+h/s_0)^{1/2}-1\right|\le|h|\).
The latter follows by differentiating along the line from
\(0\) to \(h\). Monotonicity on the ray gives
\(g(s)/H\le e^{-1/4608}\) outside that neighborhood
on the ray. This also bounds the connector, since it
bounds the endpoint \(1+ib\) and the real endpoint
satisfies \(g(1)/H<49/50<e^{-1/4608}\).

The complex Gaussian main term is
\[
 J_0(v)=2\sqrt\pi (n/2)^{3/2}s_0^{1/2}\Psi(v)^n
             \sqrt{\,2\pi/(n\kappa)\,}.
\]
It is analytic and nonzero throughout the disk.
Divide the local integral by the fixed saddle amplitude
and \(\Psi(v)^n\). To compare with the full Gaussian,
bound the three local error contributions from
\(e^{nE(h)}-1\), the amplitude ratio minus one, and the
density error. Their absolute integrals are at most
\[
 n\int_{\mathbb R}|h|^3e^{-nh^2/8}\,dh,\qquad
 \int_{\mathbb R}|h|e^{-nh^2/8}\,dh,\qquad
 \frac2n\int_{\mathbb R}e^{-nh^2/8}\,dh .
\]
These follow from \(\Re\kappa>1/3\),
\(|E(h)|\le |h|^3\le h^2/24\), and
\(|e^z-1|\le |z|e^{|z|}\).
The full complex Gaussian has modulus
\(\sqrt{2\pi/(n|\kappa|)}>\sqrt{2\pi/n}\).
Dividing the three bounds by this lower bound gives
less than \(26/\sqrt n\), \(4/\sqrt n\), and \(4/n\).

For the normalized integral outside \(|h|\le l\) up
through \(x=5\), bound the amplitude modulus
\(|s/s_0|^{1/2}\) by two, and the density ratio by two.
The contribution on the ray is at most
\(16 e^{-n/4608}\), and on the connector it is at most
\(2 e^{-n/4608}\) since \(|s|<|s_0|\) there.
On \(x\ge5\) use
\(|s/s_0|^{1/2}\le x+1\) and the exponential decrease
at rate \(n/2\) from \(x=5\). This contribution is
at most \(e^{-n/4608}(24/n+8/n^2)\).
The omitted Gaussian tail is at most
\(6/(nl)\,e^{-nl^2/6}\) by direct integration.
In total the relative error is bounded by
\[
 \left|\frac{J_n}{J_0}-1\right|
 \le \frac{30}{\sqrt n}+\frac4n+20\sqrt n\,e^{-n/4608}
 <\frac18\qquad(n\ge100{,}000).
 \tag{4}
\]
For the last comparison all terms decrease in the stated
range. At the cutoff use
\(\sqrt{100{,}000}<317\) and
\(e^{100{,}000/4608}>2^{20}\).

On \(|v|\le5/4\), apply Cauchy's estimate with radius
\(1/25\) to the analytic function \(J_n/J_0-1\).
Since (4) also ensures \(|J_n/J_0|>7/8\), its contribution
to the logarithmic derivative of \(J_n\) is less than four.
The logarithmic derivative of the prefactor of
\(\Psi(v)^n\) in \(J_0\) has modulus less than one.
Indeed \(|s_0'|=A/(2|\sqrt{1-w}|)<52/100\),
\(|s_0|>153/100\), and \(|\kappa|>1/3\).
As \(\Psi'/\Psi=A/(2s_0)\), the result is
\[
 \left|\frac{J_n'}{J_n}-\frac{n(2-s_0)}{2v}\right|<5.
 \tag{5}
\]
At zero interpret the expression \(n(2-s_0)/(2v)\) by
the analytic formula \(nA/(2s_0)\).

## Endpoint integral and exclusion of repeated roots

On any subinterval of the real segment \(0\le x\le1\),
the maximum of \(g(x)\) occurs at an endpoint. For a
direct check put \(a_w=\Re w\). The derivative of \(g^2\)
has the sign of \(Q(x)-|w|^2\), where
\[
 Q(x)=x^3(2-x)+a_w x(2-2x),\qquad
 Q'(x)=2\bigl[x^2(3-2x)+a_w(1-2x)\bigr].
\]
When \(a_w<0\), the bracket increases on \([0,1/2]\)
and is positive on \((1/2,1]\), so \(Q\) decreases and
then increases from \(Q(0)=0\). When \(a_w\ge0\),
the bracket is positive on \((0,1/2]\), and it stays
positive for \(x\ge1/2\) by \(a_w<1\) and
\(x^2(3-2x)-(2x-1)=(1-x)(2x^2-x+1)\ge0\).
In both cases \(g^2\) can have only an interior minimum.

If \(|w|\le (199/200)H\), the endpoint comparison and
\(g(1)/H<49/50\) give \(g(x)\le(199/200)H\)
throughout \([0,1]\). The real density bound from (1) gives
\[
 |K_n|\le 2n^{3/2}((199/200)H)^n,\qquad
 |J_0|=\pi n |s_0/\kappa|^{1/2}H^n>3nH^n.
\]
Thus \(|K_n|<|J_n|\) by (4) at \(n\ge100{,}000\).
Every root of \(I_n\) therefore lies at
\(|v|<5/4,\ |w|>(199/200)H\).

For \(|w|\ge (99/100)H\) we bound \(K_n/w^n\)
uniformly on the larger parameter disk. In this region
\(|w|>2/5\), because
\(H\ge 2(231/100)e^{-231/100}>21/50\).
Here \(2a e^{-a}\) decreases for \(a>1\), and
\(e^{231/100}<(11/4)^2/(1-31/100)<11\).
For \(x\le1/10\) use \(x^2/|w|\le x/4\).
At \(x=1/10\), this gives \(g(x)/|w|\le e^{-3/40}<99/100\).
At \(x=1\), \(g(1)/|w|<(49/50)/(99/100)<99/100\).
The endpoint comparison then gives the same bound
throughout \([1/10,1]\), so this part contributes
at most \(2n^{3/2}e^{-n/100}\) to \(|K_n/w^n|\).

In the remaining part put \(t=nx/2\). Then
\[
 |(1+x^2/w)^n-1|
 \le n|x^2/w|e^{n|x^2/w|}
 \le (10t^2/n)e^{t/2}.
\]
Its expectation against \(f(t)\) is less than \(50/n\),
using \(f(t)\le2\sqrt{\pi t}\,e^{-2t}\) and the
elementary Gamma integral. The missing probability at
\(t>n/20\) is at most \(3 e^{-3n/40}\), since
\(\sqrt t\le e^{t/2}\). Together,
\[
 \left|\frac{K_n}{w^n}-1\right|
 \le \frac{50}n+3e^{-3n/40}+2n^{3/2}e^{-n/100}
 <\frac1{1000}\qquad(n\ge100{,}000).
 \tag{6}
\]
Each bound decreases in the stated range, with the
exponential terms at the cutoff below \(1/10000\).

Every disk of radius \(1/1000\) about a root of \(I_n\)
is contained in the region \(|v|<13/10,\ |w|>(99/100)H\).
Indeed throughout \(|v|\le13/10\) one has
\(\left|(w/\Psi)'\right|<2\), from
\(A<55/100,\ H>21/50,\ |w|<71/100\), and
\(|\Psi'/\Psi|=|A/(2s_0)|<18/100\).
Thus \(|w|/H\) changes by less than \(2/1000\) on each
disk. Cauchy's estimate on (6), and
\(|K_n/w^n|>999/1000\), yield
\(\left|K_n'/K_n-n/v\right|<2\) at every root.
There \(J_n=-K_n\ne0\), so (5) gives
\[
 \left|\frac{J_n'}{J_n}-\frac{K_n'}{K_n}
       +\frac{n s_0}{2v}\right|<7.
\]
Since \(|s_0|>153/100\) and \(|v|<5/4\), the main
term has modulus greater than \(3n/5\).
Hence the difference of logarithmic derivatives is
nonzero, and the zero of \(I_n\) is simple.
This proves the large-terminal statement and completes
the squarefreeness proof for every order.
