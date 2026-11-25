# Eventual Strong Factorial for a constant and a monomial

Let \(\Lambda(X^\gamma)=\prod_i\gamma_i!\) be the ordinary factorial
functional. Fix positive integers \(d_1,\ldots,d_r\), and write
\(D=\sum_i d_i\). The reciprocal moment polynomial
\[
 R_n(z)=\Lambda\bigl((z+X_1^{d_1}\cdots X_r^{d_r})^n\bigr)
       =\sum_{k=0}^n\binom nk\prod_i(d_i k)!\,z^{n-k}
\]
is squarefree for every \(n\) above an absolute threshold, uniformly
in the exponent profile and the number of variables. It is also
squarefree at every \(n\) with \(2^D\ge4n^2\).
Since \(R_n'=nR_{n-1}\), these statements give
consecutive-moment separation at terminal \(n\) over every
characteristic-zero field.

A coordinate with the same exponent on both sides contributes a
nonzero scalar to each moment and can be removed before applying
the theorem. Thus the theorem also covers a common monomial factor
in variables disjoint from those in the nonconstant monomial.
The theorem supplies an absolute existence cutoff. Below any
certified numerical cutoff there are only finitely many remaining
terminal orders and, up to permutation, finitely many exponent
profiles to check.

## Large total exponent at a fixed terminal

The inequality from the
[univariate constant-support proof](../../../unpublished/strong-factorial-binomials/constant-support.tex)
gives
\[
 \frac{(d_i k)!}{((d_i n)!)^{k/n}}
 \le (k/n)^{d_i k/2}.
\]
Multiplying over \(i\) bounds the normalized mixed coefficients
of \(1+t^n+\sum_{k=1}^{n-1}b_k t^k\) by
\(b_k\le\binom nk q^{\min(k,n-k)}\), where \(q=2^{-D/2}\).
Here one uses \(k\log(n/k)\ge(\log2)\min(k,n-k)\).
By symmetry of the majorant and the inequality
\(q^{\min(k,n-k)}\le q^k+q^{n-k}\), the two sums with respective
weights \(k/n\) and \(1-k/n\) are each at most
\((1+q)^n-1<1\) if \(q\le1/(2n)\).
At a repeated root with \(|t|\le1\), subtracting \(t/n\) times
the derivative from the polynomial contradicts the second sum
bound. At a repeated root with \(|t|\ge1\), dividing \(t/n\)
times the derivative by \(t^n\) contradicts the first.

## Uniform radial density

It suffices from now on to assume \(D<2\log_2(2n)\), since
the complementary range is settled above. This ensures
\(r\le D=O(\log n)\) with absolute implied constants.

Let \(X_i\) be independent unit-rate exponential variables and set
\(B_0=\prod_i d_i^{d_i}\). Then \(R_n(z)=\mathbb E(z+B_0 T^D)^n\),
where \(T=(\prod_i X_i^{d_i}/B_0)^{1/D}\).
For \(t>0\), write \(X_i=d_i t e^{u_i}\), with
\(\sum_i d_i u_i=0\). Put \(x_i=\sqrt{d_i}\,u_i\), and use
Euclidean measure on the hyperplane
\[
 W=\{x:\sum_i \sqrt{d_i}\,x_i=0\}.
\]
The change from \(t,x\) to the coordinates
\(\sqrt{d_i}\log(X_i/d_i)\) has determinant of absolute value
\(\sqrt D/t\), since its \(t\) column is orthogonal to \(W\).
Thus \(T\) has density
\[
 f(t)=C_0 t^m e^{-Dt}\int_W
       e^{\ell(x)-tG(x)}\,dx,\quad
 m=r-1,\quad C_0=\sqrt D\prod_i\sqrt{d_i},
\]
where \(\ell(x)=\sum_i u_i\) and
\(G(x)=\sum_i d_i(e^{u_i}-1-u_i)\).
The formula is holomorphic for \(\Re t>0\). Indeed the
unsubtracted phase \(\sum_i d_i e^{u_i}\) grows exponentially
in \(\|x\|\) on \(W\) for each profile, since
\(\max_i u_i>0\) on the unit sphere.

For an absolute \(c>0\), the following bounds hold:
\[
 |\ell(x)|\le\sqrt r\,\|x\|,\qquad
 G(x)\ge c\frac{\|x\|^2}{1+\|x\|},\qquad
 \left|G(x)-\tfrac12\|x\|^2\right|\le C\|x\|^3
 \quad(\|x\|\le1).
 \tag{1}
\]
The first bound is Cauchy–Schwarz. For the second,
\(e^u-1-u\ge c u^2/(1+|u|)\) for real \(u\).
On \(u\ge0\) this follows from its series. For \(-1\le u<0\)
use the integral remainder
\(u^2\int_0^1(1-y)e^{yu}\,dy\), and for \(u\le-1\) use
convexity of \(e^{-v}-1+v\) to obtain
\(e^{-v}-1+v\ge v/e\) at \(v=-u\ge1\).
Finally \(|u_i|\le\|x\|\), which gives the second bound
upon summation and gives the third by Taylor expansion.

With \(C_2=C_0(2\pi)^{m/2}\), one has
\[
 f(t)=C_2 t^{m/2}e^{-Dt}(1+o(1))
 \tag{2}
\]
uniformly for \(\Re t\ge n/2\) and \(|\arg t|\le3e^{-D}\).
For \(D=2\) it also holds on any fixed closed sector of the
right half-plane. To verify uniformity put \(\tau=|t|\).
On \(\|x\|\le\tau^{-2/5}\), the bounds (1) give relative
integrand error \(O(\sqrt r\,\tau^{-2/5}+\tau^{-1/5})\)
when replacing \(e^{\ell(x)-tG(x)}\) by
\(e^{-t\|x\|^2/2}\). The integral of the absolute value of
the latter Gaussian, divided by the modulus of its
integral, is \((\sec(\arg t))^{m/2}\). This stays bounded
since \(m\le D\) and \(|\arg t|\le3e^{-D}\). It also stays
bounded on a fixed sector when \(D=2\).

On \(\tau^{-2/5}\le\|x\|\le1\), (1) and completion of a
square bound the absolute integrand by
\(e^{-c_1\tau^{1/5}+O(r/\tau)}e^{-c_1\tau\|x\|^2}\)
for an absolute \(c_1>0\), or one depending on the fixed sector
for \(D=2\). This is negligible relative to the Gaussian
integral, since \(m=O(\log\tau)\). Beyond \(\|x\|\ge1\)
the integral is at most
\(e^{-c_2\tau}(C\sqrt{m}/\tau)^m\). This follows by
bounding the integrand by
\(e^{-c_2\tau}e^{-c_2\tau\|x\|}\), using
\(\|x\|\ge\sum_j|x_j|/\sqrt m\) in orthonormal
coordinates, and integrating separately. The Gaussian
outside \(\|x\|\le\tau^{-2/5}\) has a negligible integral
by the same quadratic estimate. The Gaussian on all of
\(W\) has integral \((2\pi/t)^{m/2}\), first for
positive real \(t\) and then by holomorphy for \(\Re t>0\).
This proves (2), with the zero-dimensional case immediate.

## Root radius and the integral split

If \(D=1\), the polynomial is \(n!\sum_{j=0}^n z^j/j!\).
It differs from its derivative by \(z^n\) and has nonzero
constant term, which finishes this case.

For \(D\ge2\), put
\[
 A=(D/e)^D,\qquad w=Av,\qquad z=(n/e)^D B_0v,\qquad L=D/2.
\]
All roots have \(|v|<R\) for sufficiently large \(n\), with
\(R=7/6\) for \(D\ge3\) and \(R=5/4\) for \(D=2\).
To see this, let \(b_{n,k}\) be the coefficient of
\(v^{n-k}\) after scaling and renormalizing to a monic polynomial.
For \(1\le k\le n\), with \(\theta=k/n\), Stirling's bound and the
binomial entropy bound give
\[
 b_{n,k}^{1/k}\le q_D(\theta)(2ek)^{D/k},
 \qquad
 q_D(\theta)=\theta^{D-1}(1-\theta)^{-(1-\theta)/\theta}.
\]
Indeed \(l!\le el(l/e)^l\), and the product of prefactors
\(ed_i k\) is at most \((2ek)^D\) because \(d_i\le2^{d_i-1}\).
The endpoint value at \(\theta=1\) is interpreted by continuity.
For \(j=2,3\) the maximum of \(q_j\) is at the unique interior
solution to \(j\theta+\log(1-\theta)=0\), by direct differentiation.
For \(j=3\) the solution exceeds \(9/10\), and for \(j=2\) it
exceeds \(79/100\). These follow from
\(e^{27/10}>10\) and \(e^{79/50}>100/21\).
Writing \(u=1-\theta\) at the maximum gives
\(\log q_j=(j-1)\log(1-u)+ju<u\).
Thus \(q_D\le q_3<e^{1/10}<7/6\) when \(D\ge3\), and
\(q_2<e^{21/100}<5/4\). The exponential lower bounds follow by
truncating the series at degree six. The two upper bounds follow
from \(e^x\le1+x+x^2/[2(1-x/3)]\) for the positive \(x<1\) here.

For \(\theta\) below a fixed sufficiently small \(\epsilon>0\),
the bound \(\binom nk\le(en/k)^k\) instead gives
\[
 b_{n,k}^{1/k}\le e\theta^{D-1}(2ek)^{D/k}
             \le 4e^3\theta .
\]
The last step uses \((2ek)^{1/k}\le2e\), \(D\ge2\) and
\(2e\theta\le1\). For \(\theta\ge\epsilon\), the extra factor
\((2ek)^{D/k}\) in the entropy bound tends uniformly to one,
since \(D=O(\log n)\). Together the two estimates give
\(b_{n,k}/R^k\le q^k\) with an absolute \(q<1\) for all \(k\)
once \(n\) is large. For each fixed positive \(k\), the same
ratio tends uniformly to zero by the small-\(\theta\) bound.
Hence \(\sum_{k=1}^n b_{n,k}/R^k\to0\) uniformly by geometric
domination. The leading monomial dominates on \(|v|\ge R\).

After changing variables \(s=Dt/n\), the scaled moment is a
nonzero constant times
\[
 I_n(v)=\int_0^\infty(w+s^D)^n \frac nD f(ns/D)\,ds
       =K_n(v)+J_n(v),
\]
where \(K_n\) integrates from \(0\) to \(L\) and \(J_n\) from \(L\)
to infinity. Both are polynomials in \(v\). Write
\(\Phi(s)=\log(w+s^D)-s\) where needed locally. The saddle equation
is \(w=s_0^{D-1}(D-s_0)\), so its phase satisfies
\[
 e^{\Phi(s_0)}=D s_0^{D-1}e^{-s_0},
 \qquad -\Phi_{ss}(s_0)=1-\frac{D-1}{s_0}.
 \tag{3}
\]
Put \(g_v(s)=e^{-\Re s}|w+s^D|\) and \(H_v=g_v(s_0)\).
We establish, on a larger disk \(|v|\le B\) with \(B>R\), a unique
strict maximum at \(s_0\) on the ray \(x+i\Im s_0\), \(x\ge L\).
The vertical segment connecting its start to \(L\) has
\(g_v<H_v\) throughout.

## Contour geometry for \(D\ge3\)

Take \(B=6/5\), \(\delta=(3/2)D e^{-D}\), and \(\xi=\delta/D\).
The saddle equation has exactly one solution \(|s_0-D|<\delta\).
Indeed on the circle the modulus of \(s^{D-1}(D-s)\) is at least
\(\delta(D-\delta)^{D-1}\). Dividing by \(A\) gives
\((3/2)(1-\xi)^{D-1}\ge(3/2)(1-(D-1)\xi)>6/5\), since
\((D-1)\xi\le3/e^3<3/20\). Rouché's theorem gives the saddle,
which is simple and analytic in \(v\). Its real part exceeds
\(D-1\), so the right side of (3) for the second derivative
has positive real part.

Put \(a_0=\Re s_0\) and \(s=x+i\Im s_0\) for \(x\ge L\).
Factoring \(w+s^D-Ds^{D-1}\) at \(s_0\) gives
\[
 \Phi'(s)=\frac{a_0-x}{s}
     \frac{1+(s_0-D)s^{-1}\sum_{j=0}^{D-2}(s_0/s)^j}
          {1+w/s^D}.
\]
The numerator differs from one by at most
\[
 2\xi\sum_{j=0}^{D-2}[2(1+\xi)]^j
 \le \tfrac32(2/e)^D e^{(D-1)\xi}<12/17<1/\sqrt2 .
\]
The denominator differs from one by at most
\((6/5)(2/e)^D<12/25\). The corresponding argument bounds
are \(\pi/4\) and \(\pi/6\), by the tangent to a disk centered
at one. Also \(|\arg(1/s)|\le\arctan(2\xi)<\pi/12\).
The factor multiplying \(a_0-x\) has positive real part, which
proves strict increase and then strict decrease along the ray.

For the connector, first use (3) and expand \(\log(s_0/D)\)
about one. The linear term in \(\log(H_v/A)\) has modulus at
most \(\delta/D\) and the remaining terms have modulus at most
\(\delta^2/[2D(1-\xi)]\). Their sum is less than \(1/10\), since
\(\xi<3/40\) and \(\delta^2/D\le27/1600\). Hence \(H_v>(9/10)A\).
On \(s=L+iy\) with \(|y|\le\delta\), one has instead
\[
 \frac{g_v(s)}A
 \le \frac65 e^{-D/2}
        +\left(\frac{e^{1/2}}2\right)^D e^{2\delta^2/D}
 < \frac{3}{11}+\frac{(5/6)^3}{1-27/800}
 < \frac{22}{25}.
\]
The scalar comparisons use \(e^3>20\), \(\sqrt{20}>22/5\),
\(e^{1/2}<5/3\), and the geometric-series upper bound for the
last exponential. In particular the connector has a strict
uniform gap below the saddle.

## Contour geometry for \(D=2\)

Take \(B=13/10\) so that \(|w|<c=71/100\). Indeed
\(e>163/60\) by summing its series through \(1/5!\) and
\(4B/(163/60)^2<71/100\). Write
\(s_0=1+\sqrt{1-w}=1+u+ib\), using the square root with positive
real part. Then \(53/100<u<131/100\) and \(|b|\le c/(2u)\).
Indeed \(1-c\le u^2\le1+c\), using the real part and modulus
of \(1-w\), and \(2ub=-\Im w\).
The saddle is analytic, and \(1-1/s_0\) has positive real part.

On \(s=x+ib\) with \(x\ge1\), factor the derivative as
\[
 \Phi'(s)=(1+u-x)\,
 \frac{x+u-1+2ib}{x^2+1-u^2+2ib(x-u)} .
\]
The denominator has real part at least \(1-c\).
After multiplication by its conjugate, the real part of the
numerator is
\[
 (x+u-1)(x^2+1-u^2)+4b^2(x-u).
\]
If \(x\ge u\) it is positive. For \(x<u\), one has \(u>1\),
so the first term is at least \(1-c\) and the magnitude of the
negative second term is less than \(c^2(31/100)<1-c\).
This again gives strict increase and then strict decrease
along the shifted ray.

The squared modulus of the polynomial factor on the connector is
\((2-u^2+b^2-y^2)^2+4(y-ub)^2\) at \(s=1+iy\).
It is convex as a function of \(y\), since its second derivative
is \(12y^2+8-4(1+\Re w)>0\). The endpoint \(y=b\) lies
strictly below the saddle by the ray comparison. For \(y=0\)
the squared ratio with the saddle is
\[
 \frac{g_v(1)^2}{H_v^2}
   =\frac{e^{2u}\bigl((2-u^2)^2+(4+2u^2)t+t^2\bigr)}
          {4((1+u)^2+t)},\qquad t=b^2.
\]
The condition \(|w|\le c\) gives
\((1-u^2)^2+2(1+u^2)t+t^2\le c^2\).
The following exact rational check bounds the squared ratio by
\(24/25\). For each subinterval \([a,h]\) it uses
\(\Delta\le|1-u^2|\) and \(T\ge t\), discards \(t\) in the
denominator, and majorizes the exponential by its series
through degree sixteen plus a geometric bound for the tail.

    from fractions import Fraction as Q
    from math import factorial

    c = Q(71, 100)
    for k in range(530, 1310):
        a, h = Q(k, 1000), Q(k + 1, 1000)
        Delta = max(Q(0), 1 - h*h, a*a - 1)
        if Delta > c:
            continue
        T = (c*c - Delta*Delta) / (2*(1 + a*a))
        x = 2*h
        E = sum(x**j / factorial(j) for j in range(17))
        E += x**17 / factorial(17) / (1 - x/18)
        bound = E*((2 - a*a)**2 + (4 + 2*h*h)*T + T*T)
        if bound >= Q(24, 25)*4*(1 + a)**2:
            raise ArithmeticError(k)

Together with convexity and compactness of the parameter disk,
this gives a strict uniform gap below the saddle on the connector.

## Uniform estimates for the two integrals

Deform the tail across the connector to the shifted ray.
The integrand is holomorphic in the right half-plane, and (2)
makes the far vertical side vanish for sufficiently large \(n\).
Formula (2) applies uniformly along the deformed path,
where \(\Re(ns/D)\ge n/2\). For \(D\ge3\) the angle is
at most \(2\delta/D=3e^{-D}\). For \(D=2\) a fixed closed
sector works by the bounds on the saddle.

For \(D\ge3\) the factor multiplying \(a_0-x\) in the
horizontal derivative has real part at least \(c_3/|s|\)
with an absolute \(c_3>0\). Indeed its modulus has this
lower bound, up to a constant, from the displayed bounds
on the numerator and denominator. Its argument is at most
\(\pi/4+\pi/6+3/20<\pi/2\) in absolute value.
Writing \(y=x-a_0\) and integrating the derivative shows,
for another absolute \(c_4>0\),
\[
 g_v(s_0+y)/H_v\le
 \begin{cases}
   e^{-c_4 y^2/D},& L\le a_0+y\le 2D,\\
   e^{-c_4(y+D)},& a_0+y>2D .
 \end{cases}
 \tag{3a}
\]
For the first bound use \(|s|\le3D\) on this interval.
Beyond \(2D\), one has \((x-a_0)/|s|\) bounded below,
and the first bound supplies a gap at \(2D\).
For \(D=2\), the same two types of bounds with adjusted
absolute constants follow from the explicit derivative.
Its real factor has a positive minimum on every bounded
part of the parameterized ray, by compactness, and
\(\Re\Phi'\to-1\) uniformly at infinity.

The amplitude on the ray satisfies
\(\left|s^{m/2}/s_0^{m/2}\right|\le e^{|y|}\),
since \(|(\log|s|)'|\le 2/D\) and \(m\le D\).
The ratio is bounded on the connector as well:
\(|s|\le |s_0|\) for \(D\ge3\), and \(D=2\) follows by
compactness. The connector has length bounded by an
absolute constant and a uniform strict gap below \(H_v\).
These facts give uniform control of the integrals outside
the saddle after rescaling \(y=\sqrt{D/n}\,Y\).
On the first part of (3a) the majorant is
\(e^{-c_4Y^2+|Y|}\) since \(D/n<1\). The second part and
the connector contribute exponentially small errors.

Let \(\kappa=1-(D-1)/s_0\). The numbers \(D\kappa\) lie
in a fixed compact subset of the right half-plane.
For \(D\ge3\) this follows from
\(D\kappa=(1+\zeta)/(1+\zeta/D)\) with \(\zeta=s_0-D\)
and \(|\zeta|<9/40\). For \(D=2\) use compactness.
Uniformly in a fixed small disk around the saddle, one
has \(|\Phi_{sss}|\le C/D^2\). For \(D\ge3\), write the
phase as \(D\log s+\log(1+w/s^D)-s\). In a fixed larger
disk, say \(|s-s_0|\le1/2\), the second logarithm
is \(O(e^{-D})\) as \(D\) grows because
\(|w/s^D|\le B e^{-D}|D/s|^D=O(e^{-D})\).
Use Cauchy's estimate there and differentiate \(D\log s\)
explicitly. For the finitely many smaller \(D\), including
\(D=2\), choose the smaller disk uniformly by compactness
and the nonvanishing of \(w+s_0^D\).
Taylor expansion at \(s_0\) therefore gives
\(n(\Phi(s)-\Phi(s_0))=-D\kappa Y^2/2+o(1)\) uniformly
on bounded \(Y\). Also \(s^{m/2}/s_0^{m/2}\to1\)
uniformly there. Formula (2), the integrable majorant,
and the complex Gaussian integral give
\[
 J_n(v)=C_2(n/D)^{(r+1)/2}s_0^{(r-1)/2}
       (w+s_0^D)^n e^{-ns_0}
       \left(\frac{2\pi}{n\kappa}\right)^{1/2}
       (1+o(1)).
 \tag{4}
\]
Uniform relative error follows because the Gaussian
integral in \(Y\) is bounded away from zero.
All powers in the prefactor use right-half-plane branches.
The quotient of \(J_n\) by its nonzero analytic main term
converges uniformly to one on the disk. The logarithmic
derivative of the prefactor is bounded on the smaller
disk \(|v|\le R\) by Cauchy's formula, since
\((r-1)\log(s_0/D)/2-\log(D\kappa)/2\) is uniformly
bounded on \(|v|\le B\). The saddle equation cancels the
contribution from the derivative of \(s_0\) in the phase.
Consequently, uniformly on \(|v|\le R\),
\[
 \frac{J_n'}{J_n}=\frac{nA}{w+s_0^D}+O(1)
                 =\frac{n(D-s_0)}{Dv}+O(1)
 \tag{5}
\]
with the first form used at \(v=0\).

On the real interval \(0\le x\le L\), the only candidates for a
maximum of \(g_v(x)\) are the endpoints, with a strict inequality
in the interior. To verify this, put \(a_w=\Re w\). The derivative
of \(g_v(x)^2\) has the sign of \(Q(x)-|w|^2\), where
\[
 Q(x)=x^{2D-1}(D-x)+a_w x^{D-1}(D-2x),\qquad
 Q'(x)=D x^{D-2}\bigl[x^D(2D-1-2x)+a_w(D-1-2x)\bigr].
\]
If \(a_w<0\), the bracket increases strictly up to
\((D-1)/2\) and is positive thereafter. The first term is
strictly increasing on \(0<x<L\) by differentiation.
Thus \(Q\) first decreases and then increases, starting at
zero. If \(a_w\ge0\), the bracket is positive before
\((D-1)/2\). For \(D\ge3\), beyond that point its first
term is at least \(((D-1)/2)^D(D-1)>(6/5)A\), whereas the
negative term has magnitude at most \(|w|\).
For \(D=3\) the strict bound uses \(2(e/3)^3>6/5\), and
for \(D\ge4\) it uses \((D-1)(3e/8)^D>6/5\).
For \(D=2\), use
\(x^2(3-2x)-(2x-1)=(1-x)(2x^2-x+1)\ge0\)
on \([1/2,1]\) and \(a_w<1\).
Thus \(Q\) is increasing when \(a_w\ge0\).
In all cases \(g_v^2\) can have at most an interior minimum,
which proves the strict endpoint comparison.

Choose absolute constants \(0<\rho_-<\rho_+<1\) with
the uniform connector bound \(g_v(L)/H_v<\rho_-\).
We use a uniform bound on the real integral from
\(0\le x\le L\). Put \(S=\sum_i X_i\).
Weighted AM–GM gives \(DT\le S\). If \(t=T\), the constraint
\(x=Dt/n\le L\) gives \(DT\le nL\). It follows that
\[
 \mathbb E[e^{DT}\mathbf 1_{\{DT\le nL\}}]
 \le (2nL)^r/r!+4^r e^{-nL/2}
 =e^{O((\log n)^2)} .
 \tag{6}
\]
For \(S\le2nL\), integrate \(e^S\) against the joint
exponential density to obtain the simplex volume
\((2nL)^r/r!\). For \(S>2nL\), use
\(\mathbb E e^{3S/4}=4^r\) and \(e^{DT}\le e^{nL}\).
If \(|w|\le\rho_+H_v\), the real endpoint comparison gives
\(g_v(x)\le\rho_+H_v\) on the integration interval.
Thus (6) shows that \(K_n\) is exponentially smaller than
the main term in (4). Its prefactor, after removing
\(H_v^n\), is at least a constant times \(n^{r/2}/\sqrt D\),
using \(C_2\ge1\), \(|s_0/D|^{m/2}\) bounded below, and
\(|D\kappa|\) bounded above.

Uniformly where \(|w|\ge\rho_-H_v\), one has instead
\(K_n(v)/w^n\to1\). To see this, choose a small absolute
\(x_0>0\). On \(0\le x\le x_0\) we have
\(\log(1+x^D/|w|)\le x/2\), choosing \(x_0\) sufficiently
small, because \(|w|\) is bounded below by a fixed
positive multiple of \(A=(D/e)^D\). On \(x_0\le x\le L\),
the ratio \(g_v(x)/|w|\) is uniformly below one.
For large \(D\) this follows from
\[
 \frac{g_v(x)}{|w|}
 \le e^{-x_0}+C_5(e^{1/2}/2)^D ,
\]
since \(x^D e^{-x}\) increases up to \(L\).
For the remaining finitely many \(D\), use the strict real
endpoint comparison, \(g_v(L)<|w|\), and compactness.
Bound (6) now makes the normalized integral over
\(x_0\le x\le L\) exponentially small.

On \(0\le x\le x_0\), the normalized integrand against
the probability measure is \((1+x^D/w)^n\), bounded in
modulus by \(e^{DT/2}\le e^{S/2}\).
The portion \(S>(\log n)^2\) has absolute integral at most
\(4^r e^{-(\log n)^2/4}\), by \(\mathbb E e^{3S/4}=4^r\).
On \(S\le(\log n)^2\), one has \(x\le(\log n)^2/n\).
Since \(D\ge2\) and \(|w|\) stays above a fixed positive
multiple of \(A\), the normalized integrand tends uniformly
to one there and the probability tends to one.
This proves \(K_n(v)/w^n\to1\).

For large \(n\) all roots satisfy \(|v|\le R\) and
\(|w|>\rho_+H_v\). These points lie a uniformly positive
distance inside \(|v|<B,\ |w|>\rho_-H_v\).
Indeed \(H_v/A\) is bounded above and below by absolute
positive constants. On each disk slightly smaller than \(B\)
it is uniformly Lipschitz, by Cauchy's formula on the
analytic function \((s_0/D)^{D-1}e^{D-s_0}\) whose modulus
is \(H_v/A\). For \(D=2\) this also follows by the
explicit saddle formula. Hence Cauchy's formula on the
uniform neighborhoods gives \(K_n'/K_n=n/v+o(1)\)
at every root, with an absolute error.
There \(J_n=-K_n\ne0\) and
\[
 \frac{J_n'}{J_n}-\frac{K_n'}{K_n}
   =-\frac{n s_0}{Dv}+O(1).
\]
Since \(s_0\) stays away from zero and \(|v|\le R\), this
expression is nonzero for large \(n\). Consequently \(I_n'\)
is nonzero at every zero of \(I_n\), and the original
polynomial \(R_n\) is squarefree. The coefficients of \(R_n\)
are integers, so its rational Bézout identity with
\(R_{n-1}\) also holds over every characteristic-zero field.
