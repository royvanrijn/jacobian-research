# Ordinary mirrored profiles can shrink the exact Gaussian angular gap

The exact physical Gaussian midpoint clocks need not give a lower
bound on the corresponding factorial cross-degree angular gap.
This fails already at orders two and three, for fully shared ordinary
mirrored rows. It also fails along arbitrarily diffuse mirrored
profiles with zero signed skew and fixed midpoint curvature sum.
Thus an angular transfer using the
[matched tanh hierarchy](MULTIVARIATE_MATCHED_TANH_ANGULAR_HIERARCHY.md)
must allow a loss in the reference gap. The result does not exclude
the smaller quantitative reserve in the
[polynomial Gaussian physical gap theorem](MULTIVARIATE_GAUSSIAN_POLYNOMIAL_PHYSICAL_GAPS.md).

## The two angular gaps

Take rows \((\ell,\ell+d),(\ell+d,\ell)\) and put \(m=\ell+d/2\).
The endpoint-normalized degree-two and degree-three polynomials
are \(1+2u z+z^2\) and \(1+3v z+3v z^2+z^3\), where
\[
 u=\frac{(2\ell+d)!^2}{(2\ell)!(2\ell+2d)!},\qquad
 v=\frac{(3\ell+d)!(3\ell+2d)!}{(3\ell)!(3\ell+3d)!}.
\]
Their positive arguments are
\(\theta_2=\arccos(-u)\) and \(\theta_3=\arccos((1-3v)/2)\),
with the structural degree-three root at \(-1\) as well.
For \(u<1/4\), \(v<1/3\) these satisfy
\(\pi/3<\theta_3<\pi/2<\theta_2<2\pi/3\).
The minimum circular cross-degree angle gap is
\(g(u,v)=\theta_2-\theta_3\). Indeed the structural root
is farther away, and reflection of one angle to the
opposite half-circle cannot reduce their distance.

The two
[exact Gaussian midpoint references](MULTIVARIATE_MATCHED_MIDPOINT_CLOCK_ANGLE_TEST.md)
replace these numbers by
\[
   u_G=e^{-C_2/4},\qquad v_G=e^{-C_3/3},\qquad
                C_j=2j d^2\psi'(1+jm).
\]
They use the actual linked physical clocks and satisfy the same
formula for their gap when \(u_G<1/4\), \(v_G<1/3\).

For the concrete ordinary pair \((1,7),(7,1)\) one has
\(g(u,v)\in(530/1000,531/1000)\) and
\(g(u_G,v_G)\in(532/1000,533/1000)\).
The shrinkage exceeds \(24/10000\) radians.
The following rigorous Arb checks use exact rational input and
\(\psi'(j)=\pi^2/6-\sum_{k=1}^{j-1}k^{-2}\) at positive integers.

~~~python
from math import factorial as f
from flint import arb, ctx, fmpq

ctx.prec = 128
u = arb(fmpq(f(8)**2, f(2)*f(14)))
v = arb(fmpq(f(9)*f(15), f(3)*f(21)))

def trigamma(j):
    return arb.pi()**2/6 - sum(arb(fmpq(1,k*k)) for k in range(1,j))

def gap(a, b):
    assert 0 < a < arb(fmpq(1,4)) and 0 < b < arb(fmpq(1,3))
    return (-a).acos() - ((1-3*b)/2).acos()

g = gap(u, v)
gG = gap((-36*trigamma(9)).exp(), (-72*trigamma(13)).exp())
assert arb(fmpq(530,1000)) < g < arb(fmpq(531,1000))
assert arb(fmpq(532,1000)) < gG < arb(fmpq(533,1000))
assert gG-g > arb(fmpq(24,10000))
~~~

## The same loss for arbitrarily diffuse inputs

Let \(d=6k\), \(m=d^2/12=3k^2\), and \(\ell=m-d/2\)
for integer \(k\ge2\). These are ordinary fully shared
profiles with \(W=2d^2/m=24\) and \(d/m\to0\).
The signed moments cancel by reflection. For every fixed
order, higher unsigned Gauss moments of degree greater than
one tend to zero, since there are \(2d\) factors with sizes
\(O(d/m)\).

In coordinate \(w=s-1/2\), the exact normalized
log-coefficient curvature is
\[
 K_j(w)=L_j''(1/2+w)
       =2d^2\int_0^\infty e^{-m x}h(x/j)\cosh(dwx)\,dx,
               \qquad h(y)=y/(e^y-1).
\]
This is the exact kernel in the matched hierarchy.
Rescaling \(x\) by \(m^{-1}\), using \(0<h(y)\le1\), gives
\(C_j=K_j(0)\to24\) and \(K_j''(0)\sim4d^4/m^3\).
For \(|w|\le1/2\) the remainder after the quadratic term in
the curvature is \(O(d^6/m^5)\).
Indeed the even exponential tail is bounded by a constant
times \((dx)^4e^{dx/2}\) inside the integral, and \(d/m\to0\).
Integrating twice with zero values at the endpoints yields
\[
 L_j(s)+C_j s(1-s)/2
   = \frac{K_j''(0)}{24}\bigl((s-1/2)^4-1/16\bigr)
                                      +O(d^6/m^5).
\]
Consequently with \(\varepsilon=d^4/m^3\),
\[
 \log(u/u_G)=-\varepsilon/48+o(\varepsilon),\qquad
 \log(v/v_G)=-5\varepsilon/162+o(\varepsilon).
\]
The limits of \(u_G,v_G\) are \(u_0=e^{-6}\) and \(v_0=e^{-8}\).
Differentiating the exact two-angle formula gives
\[
 \lim_{k\to\infty}\frac{g(u,v)-g(u_G,v_G)}{\varepsilon}
   =-\frac{u_0}{48\sqrt{1-u_0^2}}
      +\frac{5v_0}{108\sqrt{1-((1-3v_0)/2)^2}}<0 .
\]
For the sign, the second denominator's square root exceeds \(1/2\),
so its positive term is less than \(5v_0/54\).
This is less than \(u_0/48\) because \(e^2>40/9\).
The magnitude of the negative term exceeds \(u_0/48\).
Thus the exact ordinary gap is strictly smaller than the
Gaussian midpoint gap throughout a diffuse sequence.
