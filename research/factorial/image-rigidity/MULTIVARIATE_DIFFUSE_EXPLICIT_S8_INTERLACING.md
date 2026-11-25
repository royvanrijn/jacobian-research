# A half-saddle phase bound for Gaussian coupling \(S\ge8\)

Write
\[
 B_j(z;S)=\sum_{k=0}^j\binom jk e^{-S k(j-k)/(2j)}z^k .
\]
For \(S\ge8\) and \(j\ge600\), the roots of \(B_j\) and
\(B_{j+1}\) are simple and strictly interlace on the unit circle.
Thus the Gaussian references at the slowest physical upper time
have disjoint arguments there. The
[inner comparison and outward motion](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md)
extend the separation to every physical upper time between
the slowest one and equal temperature.

It suffices to prove the new range \(8\le S\le13\), as the
[previous all-degree estimate](MULTIVARIATE_DIFFUSE_EXPLICIT_S13_INTERLACING.md)
covers \(S\ge13\). In the new range there is an explicit analytic
logarithm of the polynomial half on \(|z|\le6/5\).
We obtain the saddle estimate by bounding the far horizontal
base at its endpoints and cancelling the local cubic term
in the symmetric Gaussian integral.

## The analytic polynomial half

Put \(a=S/2\in[4,13/2]\), write
\[
 G(x,z)=e^{-x^2/(4a)}(1+ze^{x-a}),\qquad
 I_j(z)=\sqrt{\frac j{4\pi a}}\int_{-\infty}^a G(x,z)^j\,dx .
\]
Expansion and Gaussian integration show
\(B_j(z;S)=I_j(z)+z^jI_j(1/z)\).
Each half has positive coefficients. Let
\(F(x,z)=-x^2/(4a)+\log(1+ze^{x-a})\) where the
logarithmic series converges.

For \(|z|\le6/5\), use the saddle equation
\[
 q=z e^{-a}\exp\!\left(\frac{2a q}{1+q}\right),\qquad
 x_*=\frac{2aq}{1+q}.
\]
On \(|q|=2e^{-a}\) one has \(2a|q|/(1-|q|)\le1/3\):
use \(e^{-4}<1/50\) and the decrease of \(ae^{-a}\).
Thus the modulus of the right-hand side is at most
\((6/5)e^{-a}e^{1/3}<2e^{-a}\).
Rouché gives one simple analytic solution in the disk.
It has conjugation symmetry and \(|x_*|\le1/3\).
Also \(2e^{-a}\le2(7/19)^4<1/27\).
Writing \(F_*=F(x_*,z)\) and \(A_*=-F_{xx}(x_*,z)\) gives
\[
 |F_*|<1/20,\qquad
 |zF_*'(z)|=\left|\frac{q}{1+q}\right|\le1/26,\qquad
 |2aA_*-1|\le2/5 .
\tag{1}
\]
The first bound uses \(1/144+1/24<1/20\).
For the last use \(|A_*-1/(2a)|\le |q|/(1-|q|)^2\),
\(2e^{-4}<1/25\) and the decrease of \(ae^{-a}\).
A direct interval check improves the real lower bound to
\(\Re A_*\ge7/100\). On the successive half-width intervals
\([4,9/2], [9/2,5], [5,11/2], [11/2,6], [6,13/2]\),
use upper bounds \(1/27,2/81,1/72,1/108,1/192\) for \(|q|\).
They follow from \(2e^{-4}<1/27\), \(e^{-1/2}<2/3\),
and \(e^{-1}<3/8\). At each interval, the reciprocal of twice
its upper \(a\) endpoint minus \(q_{\max}/(1-q_{\max})^2\)
is at least \(7/100\). The amplitude \(C_*=(2aA_*)^{-1/2}\)
has an analytic logarithm with \(|\log C_*|\le1/3\) and
\(|C_*|>4/5\).

## Horizontal and connector estimates

Shift the integral horizontally through \(x_*(z)\) and then
vertically at \(\Re x=a\) back to its original endpoint.
This uses the entire integrand \(G^j\) and its Gaussian decay
as the left closing side tends to negative infinity.

On the horizontal line write \(x=x_*+t\) with \(t\) real.
For \(\Re x\le a/2\) put \(p=ze^{x-a}\).
Here \(|p|\le b=(6/5)e^{-a/2}<1/6\), and
\[
 \Re F_{xx}(x,z)
 =-\frac{1}{2a}+\Re\frac{p}{(1+p)^2}
 \le-\frac{1}{2a}+\frac{b}{(1+b)^2}
 \le-\frac{3}{200a}.                                        \tag{2}
\]
For the first inequality, with \(r=|p|\) and \(v=\cos\arg p\)
the real part is
\([r(1+r^2)v+2r^2]/(1+2rv+r^2)^2\).
It increases with \(v\in[-1,1]\) when \(r\le1/6\), since
the derivative numerator has sign bounded below by that of
\((1+r^2)^2-2r(1+r^2)-8r^2>0\).
Its maximum \(r/(1+r)^2\) increases with \(r<1\).
For the second inequality, \(2ab/(1+b)^2\) decreases with
\(a\ge4\). Its logarithmic derivative is
\(1/a-1/2+b/(1+b)<0\), and at \(a=4\) the bound
\(b\le(6/5)(7/19)^2\) gives \(2ab/(1+b)^2<97/100\).

On the horizontal interval \(|t|\le1/2\),
\(|p|\le 2e^{-a}e^{1/2}<1/15\). Direct differentiation gives
\[
 |F_{xxx}|=\left|\frac{p(1-p)}{(1+p)^3}\right|<9/100,\qquad
 |F_{xxxx}|=\left|\frac{p(1-4p+p^2)}{(1+p)^4}\right|<3/25 .
\]
Writing \(F(x_*+t,z)-F_*=-A_*t^2/2+E(t)\) and
\(c_3=F_{xxx}(x_*,z)/6\), we have
\[
 |E(t)|\le3|t|^3/200,\qquad |E(t)-c_3t^3|\le |t|^4/200.
\]
In particular \((\Re A_*)t^2/2-|E(t)|\ge t^2/37\),
since \(7/200-3/400>1/37\).
The integral of \(j c_3t^3 e^{-jA_*t^2/2}\) over
\([-1/2,1/2]\) vanishes. From the exponential series the
remaining difference is bounded by
\[
 \left(\frac j{200}|t|^4+
              \frac{9j^2}{80000}|t|^6\right)e^{-jt^2/37}.
\]
Its integral after multiplication by \(\sqrt{j/(4\pi a)}\)
is at most \(25/j\), by the Gaussian fourth and sixth
moments. For example the exact constant from these two
integrals is at most
\[
 \frac{61}{40}\left(\frac{3\cdot37^2}{800}
              +\frac{135\cdot37^3}{640000}\right)<25 .
\]

On both ends \(t=\pm1/2\) of the local interval, the
real phase is below \(\Re F_*\) by at least \(1/148\).
Its outward derivative is at most \(-1/50\):
Taylor expansion of \(F_x\) gives
\(-\Re A_*/2+(9/100)/8\le-7/200+9/800<-1/50\).
These bounds persist to the left and to \(\Re x=a/2\)
by (2), giving the loss
\(1/148+(|t|-1/2)/50\) throughout the nonlocal left
horizontal portions. The same bounds hold for the two
tails of the matching Gaussian. A pair of tails with
this loss costs at most
\((100/\sqrt{4\pi a})e^{-j/148}/\sqrt j
<15e^{-j/148}/\sqrt j\) after normalization.

For \(a/2\le\Re x\le a\) and \(|\Im x|\le1/3\), the
base modulus is at most \(e^{1/144}\) times the maximum of
\[
                 e^{-u^2/(4a)}(1+(6/5)e^{u-a})
                         \qquad (a/2\le u\le a).
\]
The maximum is attained at an endpoint. Indeed the derivative
has nonpositive sign precisely when
\((6/5)e^{-a}\le e^{-u}u/(2a-u)\), and the right-hand side
decreases because \(1/u+1/(2a-u)\le8/(3a)<1\).
The larger endpoint bound at \(a=4\) suffices for all \(a\ge4\).
Using \(e>19/7\), \(e^{-1/4}<32/41\) and \(e^{-\Re F_*}<20/19\),
we obtain
\[
 |G(x,z)e^{-F_*}|
 \le \frac{144}{143}\frac{20}{19}
           \max\left\{\frac{32}{41}\left(1+\frac65(7/19)^2\right),
                        \frac{11}{5}\frac7{19}\right\}
       <\frac{97}{100}.                                      \tag{3}
\]
Here \(e^{-\Re F_*}<20/19\) follows from (1) and the geometric
bound on \(e^{1/20}\). The bounded horizontal portion and
vertical connector have combined length less than \(4\).

## Phase comparison

The calculated Gaussian completion, local error, and all
three nonlocal contributions give
\[
 I_j(z)=e^{jF_*(z)} C_*(z)(1+\eta_j(z)),\qquad
 |\eta_j(z)|\le \frac{32}{j}
     +\frac{38}{\sqrt j} e^{-j/148}+5\sqrt j(97/100)^j
       \qquad(j\ge600).
\]
The three terms decrease on that range, and their sum is less than
\(1/12\) at \(j=600\), using \(\sqrt{600}>24\),
\(e^{-600/148}<1/54\), and \((97/100)^j\le e^{-3j/100}\).
For \(j\ge1700\) their sum is less than \(1/50\).

A separate adjacent estimate avoids adding the local errors at
both orders. Divide the local integral by \(e^{jF_*}C_*\)
and put \(y=\sqrt j\,t\). Its real-order continuation is
\[
 {\cal J}_j=\frac1{\sqrt{4\pi a}\,C_*}
       \int_{-\sqrt j/2}^{\sqrt j/2}
       \exp[-A_*y^2/2+j E(y/\sqrt j)]\,dy .
\]
The derivative of the exponent with respect to \(j\) is
\(E(t)-t E'(t)/2\), of modulus at most \(3|t|^3/80\)
by the local cubic bounds. The scaled exponent is bounded
above in real part by \(-y^2/37\).
Gaussian integration bounds the contribution of the
differentiated integrand by \(10/j^{3/2}\).
The two moving endpoints cost at most
\(e^{-j/148}/(10\sqrt j)\), using \(|C_*|>4/5\).
Integrate these estimates for \({\cal J}_{n+1}-{\cal J}_n\).
The remaining contributions to \(\eta_{n+1}-\eta_n\) come
from the actual nonlocal left horizontal portions and the
far bounded contour at the two orders. There is no Gaussian
completion term in this difference. The preceding tail bounds
give, for \(n\ge600\),
\[
 |\eta_{n+1}-\eta_n|
 \le\frac{10}{n^{3/2}}
   +\frac{381}{10\sqrt n}e^{-n/148}
   +10\sqrt n(97/100)^n<\frac{31}{1000}.
\]
For the last inequality use \(\sqrt{600}>24\),
\(e^{-600/148}<1/54\), and \(e^{18}>(19/7)^{18}\);
each term decreases for \(n\ge600\).

The half factorization is analytic and nonvanishing on \(|z|\le6/5\).
Let \(H_j=jF_*+\log C_*+\log(1+\eta_j)\), with the
branches real on the positive radius. The amplitude is common
to both orders. Since the line segment between \(\eta_n,\eta_N\)
stays in \(|w|<1/12\), their logarithm difference is bounded by
\(12|\eta_N-\eta_n|/11\). Cauchy at radius \(1/5\) gives on \(|z|\le1\)
\[
 |H_{n+1}'-H_n'| \le 5(1/20+(12/11)(31/1000))<1/2,\qquad
 |H_j'|\le5(j/20+1/3+1/11)<j/2 .
\]
The first bound and conjugation at the positive radius give
\[
       0<T_{j+1}(\theta)-T_j(\theta)<2\pi
               \quad(0<\theta<2\pi),\qquad
       T_j(\theta)=j\theta-2\Im H_j(e^{i\theta}).
\]
Indeed the error from \(\theta\) is less than
\(\min(\theta,2\pi-\theta)\).
The second bound makes \(T_j\) strictly increasing
from \(0\) to \(2j\pi\). The exact half decomposition
gives zeros at the odd multiples of \(\pi\) in \(T_j\).
These strict inequalities give all alternating roots in
argument and show they are simple.

