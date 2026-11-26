# All degree physical separation through private unit curvature two thirds

Take an ordinary two-monomial profile with exactly one private
unit row in each orientation and arbitrary shared rows. Put
\[
 V=\sum_{\rm shared}\frac{(b_i-a_i)^2}{\min(a_i,b_i)},
                  \qquad \Delta=\sum_i(b_i-a_i).
\]

**Theorem.** If \(V\le2/3\), all ordinary consecutive physical
factorial moment pairs are coprime, at arbitrary integer
signed degree difference. The first positive angle box is
handled below for every \(n\ge13\) by one analytic comparison.
The existing higher-label theorem handles the other boxes,
and the established ordinary terminals through thirteen
supply the low orders. No eventual cutoff or new finite
moment computation is used.
[Higher angle labels](MULTIVARIATE_PRIVATE_UNIT_ALL_HIGHER_LABELS.md)
[Ordinary terminal thirteen theorem](MULTIVARIATE_TERMINAL_THIRTEEN_ALL_PROFILES.md)

The change to the weak-curvature comparison is to keep the
imaginary part of each endpoint packet in the Taylor error
instead of replacing it by the entire modulus. At a first-box
contact the two endpoint phases add to a specified number,
which bounds this error where the angle is smallest.
[Previous weak curvature criterion](MULTIVARIATE_PRIVATE_UNIT_WEAK_TWENTY_THREE_SIXTIETHS.md)

## Keep the endpoint phases in the error

Take \(n\ge13,\ N=n+1\), and suppose the two moments have a
common physical root in the first positive angle box. Let
\(\theta\) be its argument and \(x_n,x_N\) its log radii after
separate endpoint normalizations.
The private unit factorials cancel the binomial coefficient.
For each orientation, the normalized endpoint laws link as
\(\log T_n\stackrel{\rm law}=\log T_N+Z\) with independent noise \(Z\)
from the actual adjacent Beta factors.

Use \(F(s)=(1-e^{s+i\theta})^{-1}\),
\(\alpha=1/\sin\theta\), and \(b=n\sin\theta\).
At a fixed real shift \(x\) put \(E_N=\mathbb E F(x+\log T_N)\),
and define \(E_n\) at the same \(x\).
The Taylor argument for the linked resolvent gives the
following stronger form of its relative error, with
\(\phi_N=\arg E_N\):
\[
       |E_N-E_n|\le \epsilon \sin\phi_N\,|E_N|,\qquad
 \epsilon=\frac{\mathbb E Z}{\sin\theta}
                 +\frac{\mathbb E[Z^2e^{\alpha|Z|}]}{\sin^2\theta}.
                                                               \tag{1}
\]
Indeed the first and second variations are both bounded by
multiples of \(\operatorname{Im}E_N\), using \(F'\), \(F''\)
and the exponential bound on \(\operatorname{Im}F\) under
real shifts. The absolute bound of that argument follows
by replacing \(\sin\phi_N\) by one.
[Linked resolvent Taylor estimate](MULTIVARIATE_PRIVATE_UNIT_HIGHER_LABEL_LINKED_RESOLVENT.md#resolvent-variation-at-any-horizontal-position)

Write \(K_m(x,\theta)\) for the log modulus of the forward
packet at log radius \(x\) minus that of the reverse packet
at log radius \(-x\). Take the reverse packet with the same
positive angle before conjugating it in the geometric sum.
In particular both packet phases lie in \((0,\pi-\theta)\).
At the upper root of first label their sum is \(N\theta-\pi\).
This is the first-label phase equation of the
[geometric angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry).
Use (1) at \(x_N\) for the forward packet and \(-x_N\)
for the reverse. A common upper bound \(\epsilon<1\) in
both orientations gives
\[
            |K_N(x_N,\theta)-K_n(x_N,\theta)|
               \le \frac{\epsilon}{1-\epsilon}(N\theta-\pi).
                                                               \tag{2}
\]
Indeed each log-modulus change is at most
\(\epsilon\sin\phi_N/(1-\epsilon)\), and the sum of the
two sines is at most the sum of the phases.

The modulus equations are \(m x_m=K_m(x_m,\theta)\).
In the horizontal variable \(|\partial_x K_n|\le2/\sin\theta\)
by the exact first-derivative packet bound. Subtract the
two equations and use \(|x_n|\le V(n-1)/(2n)\)
from the log chord and Eneström–Kakeya.
The half-sized error in physical endpoint drift is
\(V/(2nN)\). Because the physical root is common,
\(x_N-x_n=\log(G_N/G_n)\), with \(G_m\) the normalization
scale. For \(\Delta\ne0\) this gives
\[
 |\Delta|<
 \frac{N\sin\theta}{N\sin\theta-2}
       \left[\frac V2+\frac{\epsilon}{1-\epsilon}(N\theta-\pi)\right].
                                                               \tag{3}
\]
Here the contraction factor is at least one, so it absorbs
the positive \(V/(2n)\) drift error and the negative
\(-V/(2n)\) in the bound on \(|x_n|\).
The strict inequality uses \(\log(N/n)>1/N\).
[Half-sized endpoint drift](MULTIVARIATE_PRIVATE_UNIT_WEAK_TWENTY_THREE_SIXTIETHS.md#integer-drift-with-the-sharper-packet-estimate)
[Eneström bound from the chord](MULTIVARIATE_PRIVATE_UNIT_WEAK_THREE_TENTHS.md#modulus-estimate-with-the-weak-row-cap)

## A common noise bound

Set \(v=2/3\). In either orientation the weak row cap
\(|d_i|/a_i\le V\) gives \(0\le\mu=\mathbb E Z\le V^2/(4nN)\).
For \(L_0=\Lambda''(0)\), where \(\Lambda(t)=\log\mathbb E e^{tZ}\),
we have \(L_0\le V/(nN)\). At \(\pm\alpha\), let
\(C_2=\mathbb E[Z^2\cosh(\alpha Z)]\).
The [evenized exponential moment estimate](MULTIVARIATE_PRIVATE_UNIT_WEAK_TWENTY_THREE_SIXTIETHS.md#evenized-exponential-moment)
gives \(C_2\le c_2\mathbb E Z^2\), with
\[
 c_2\le
 \frac{(1+\eta^2)/(1-\eta^2)^2+
       v^2/[2b_0 N_0(1-\eta)]+v/[b_0^2(1-\eta)^2]}
      {1-v^2/(4b_0N_0)-v/[2b_0^2(1-\eta)]},
 \quad \eta=v/b_0
\]
when \(b\ge b_0=31/10\), \(N\ge N_0=14\).
This follows from the finite Beta log-variance sum at
\(\pm\alpha\), the signed mean bound and
\(e^x\le(1-x)^{-1}\). The domain is open around
\([-\alpha,\alpha]\) since \(b>V\).
The displayed rational constant is less than \(4/3\).
Cauchy–Schwarz in the odd exponential term then yields
\[
 \mathbb E[Z^2e^{\alpha|Z|}]
 \le \big(c_2+\sqrt{2c_2(c_2-1)}\big)\mathbb E Z^2
 \le(23/10)V/(nN).
\]
Here \(\mathbb E Z^2\le[1+v/(16\cdot13\cdot14)]V/(nN)\)
and \(4/3+\sqrt{8/9}<16/7\).
Consequently (1) is bounded by
\[
                       \epsilon\le\frac{v^2}{56b}
                                       +\frac{(23/10)v}{b^2}.
                                                               \tag{4}
\]

## A single polynomial positivity check

Write \(u=N\theta/\pi\), so \(1<u<2\).
Define
\[
                  C=(24/25)(157/50),\qquad D=(13/14)C .
\]
Since \(\theta<2\pi/N\le\pi/7\), the cubic sine bound
gives \(N\sin\theta\ge Cu\) and \(b\ge Du\).
Also \(b\ge n\sin(\pi/n)>31/10\) for \(n\ge13\), by
\(\pi-\pi^3/(6n^2)\) and \(157/50<\pi<22/7\).
Thus the right side of (3) is bounded by
\[
      \frac{Cu}{Cu-2}
       \left[\frac v2+\frac{22}{7}(u-1)
                              \frac{e(u)}{1-e(u)}\right],
\quad
      e(u)=\frac{v^2}{56\max(b_0,Du)}
                       +\frac{(23/10)v}{\max(b_0,Du)^2}.        \tag{5}
\]
We show (5) is strictly less than one on the whole interval.

Let \(u_*=b_0/D\in(1,2)\). On \([1,u_*]\), \(e(u)=e(1)\),
so it suffices to check the positive gap
\[
              1-v/2-2/(Cu)
                    -(22/7)(u-1)e(1)/(1-e(1)).
\]
This is concave in \(u\); its endpoint values are positive.
For \(u\ge u_*\), write \(e(u)=\beta/u+\gamma/u^2\)
with \(\beta=v^2/(56D),\ \gamma=(23/10)v/D^2\).
Multiplying the same gap by the positive
\(u(u^2-\beta u-\gamma)\) gives the cubic
\[
 \big((1-v/2)u-2/C\big)(u^2-\beta u-\gamma)
                   -(22/7)u(u-1)(\beta u+\gamma).
\]
Every coefficient of this cubic in powers of \(u-u_*\)
is positive, as checked below.
Thus the integer \(\Delta\) in (3) must be zero.
The
[balanced Euclidean theorem](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md)
then excludes a common root. Positivity on the real axis
in even degree excludes real common roots, and conjugation
covers the lower half-plane.

## Exact scalar replay

~~~python
from fractions import Fraction as F

v, b0 = F(2,3), F(31,10)
eta = v/b0
K = v*v/(56*b0)+v/(2*b0*b0*(1-eta))
G = ((1+eta*eta)/(1-eta*eta)**2
                +v*v/(28*b0*(1-eta))+v/(b0*b0*(1-eta)**2))
assert K < 1 and G/(1-K) < F(4,3)
assert F(16,7)*(1+v/(16*13*14)) < F(23,10)
assert F(157,50)-F(22,7)**3/(6*13**2) > b0
assert 1-F(22,49)**2/6 > F(24,25)

C, W = F(24,25)*F(157,50), F(23,10)
D = F(13,14)*C
q = b0/D
e0 = v*v/(56*b0)+W*v/b0**2
assert 1<q<2 and e0<1 and C>2
for u in [F(1), q]:
    assert 1-v/2-2/(C*u)-F(22,7)*(u-1)*e0/(1-e0) > 0

beta, gamma = v*v/(56*D), W*v/D**2
A, B, J = 1-v/2, 2/C, F(22,7)
p3 = A-J*beta
p2 = -A*beta-B-J*(gamma-beta)
p1 = -A*gamma+B*beta+J*gamma
p0 = B*gamma
assert p3>0
assert p2+3*q*p3>0
assert p1+2*q*p2+3*q*q*p3>0
assert p0+q*p1+q*q*p2+q**3*p3>0
~~~
