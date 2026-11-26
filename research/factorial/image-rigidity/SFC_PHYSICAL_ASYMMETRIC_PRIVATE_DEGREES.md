# Sharp private-degree condition for binary transfer over a positive base

Let \(U_0,V_0>0\) have finite mixed moments, and let
\[
   U=U_0\prod_i X_i^{a_i},\qquad
   V=V_0\prod_i Y_i^{b_i},
\]
where all private factors are independent rate-one exponentials
and independent of the base. The reserved exponents are positive
integers. Denote their two sums by \(d_U,d_V\).

**Theorem.** If \(d_U,d_V\ge1\) and \(d_U+d_V\ge4\), the
moments \(\mathbb E(U+tV)^n\) and \(\mathbb E(U+tV)^{n+1}\)
have disjoint root arguments at every \(n\ge1\). The coefficients
of the two degrees may be normalized by independent positive
changes of variable. In particular, a physical Gaussian
correlation or merger among the other coordinates preserves
binary coprimality when these private degrees stay reserved.
The condition is sharp for a theorem uniform over all positive
base laws and all ways to realize the private degrees. A
strictly positive finitely supported base with private degrees
\((1,2)\) gives an exact common root at orders four and five.

When each private degree is at least two, apply the
[existing private reserve theorem](SFC_PHYSICAL_BINARY_RESIDUAL_PRIVATE.md).
It remains to prove the case with one degree at least three and
the other at least one.

## The asymmetric coefficient envelope

Take \(d_U\ge1,d_V\ge3\), reversing the variable if necessary.
Write \(m\) for the order. Normalize the endpoints to one, so the
polynomial is
\[
       T_m(z)=1+z^m+\sum_{k=1}^{m-1} c_{m,k}z^k .
\]
Hölder bounds the endpoint-normalized base factors by one.
For a private exponent \(a\), the log convexity of
\((aj)!/(j!)^a\) in \(j\), as in the
[private reserve proof](SFC_PHYSICAL_BINARY_RESIDUAL_PRIVATE.md#a-distribution-independent-coefficient-bound),
bounds its endpoint-normalized factor at index \(k\) by
\((k!/(m!)^{k/m})^a\). The factor on the other side uses \(m-k\).
Since each such normalized bound is at most one, the exponents
one and three suffice. After multiplication by the binomial
coefficient, this gives
\[
        0<c_{m,k}\le A_{m,k}=\left(\frac{k!}{g_m^k}\right)^2,
        \qquad g_m=(m!)^{1/m}.                               \tag{1}
\]
The proof below uses this envelope independently at each degree.
It requires no additional coefficient relation across the base
moments.

## Separation for every lower degree at least forty

Split at \(k=\lfloor m/2\rfloor\), writing
\(T_m(z)=L_m(z)+z^m U_m(1/z)\), with both half polynomials having
constant term one. On \(|w|\le e^{1/m}\) the two half bounds are
\[
 \begin{array}{c|cc}
 & |H(w)-1| & |wH'(w)|\\ \hline
  L_m & 1/100 & 3/100\\
  U_m & 9/40 & 7/25 .
 \end{array}                                                  \tag{2}
\]

Here are exact derivations of the bounds. Integral comparison for
\(\log(m!)\) gives \(m/e\le g_m\le m e^{-1}(em)^{1/m}\).
The upper relative bound decreases in \(m\) and is below \(5/12\)
at \(m=40\). Use \(e^{1/m}\le 40/39\).
For the coefficient at index \(h\) in the left half, multiplied
by \(e^{h/m}\), use
\[
 (h!)^2\left[(87/1280)^2(40/39)\right]^h
 \quad(1\le h\le6).
\]
This follows from \(g_m\ge m/e\) and \(e<87/32\).
For \(7\le h\le m/2\), integral comparison for \(\log(h!)\)
bounds the \(h\)-th root by
\[
       (h/m)^2 \exp\bigl(2(1+\log h)/h+1/m\bigr)<2/5.
\]
For \(7\le h\le20\) use \(m\ge40\) and the increase of
\(h^2\exp(2(1+\log h)/h)\) in \(h\). For \(h\ge20\) use
\(h/m\le1/2\) and the decrease of \((1+\log h)/h\).
Both give at most \(\tfrac14 e^{1/8}20^{1/10}<2/5\).

For the right half the corresponding coefficient is bounded by
\[
      \left(\frac{g_m^h}{m(m-1)\cdots(m-h+1)}\right)^2 e^{h/m}.
\]
For \(1\le h\le6\) use
\[
        \frac{[(5/12)^2(40/39)]^h}
             {\prod_{j=0}^{h-1}(1-j/40)^2}.
\]
For \(7\le h\le m/2\), the average of
\(\log(1-j/m)\) for \(0\le j<h\) is at least the corresponding
integral average, which is at least \(\log 2-1\).
The \(h\)-th root is consequently bounded by
\((5/12)^2(40/39)(e/2)^2<1/3\).
Summing these bounds and their products with \(h\) gives (2).
All numerical inequalities in this calculation are rational
after replacing \(e\) by the bounds \(2718/1000<e<87/32\).
For instance, the left tail with first index seven uses
\(\sum_{h\ge7}q^h=q^7/(1-q)\) and
\(\sum_{h\ge7}h q^h=q^7(7-6q)/(1-q)^2\) at \(q=2/5\).
The checker below verifies the rational margins for both halves.

On \(|z|=e^{1/m}\), the term \(z^m U_m(1/z)\) dominates
\(L_m(z)\) by (2), since \(e(1-9/40)>1+1/100\).
All zeros of the dominating polynomial lie inside the contour,
as \(U_m\) has no zeros in \(|w|\le e^{1/m}\).
Apply the reciprocal argument using
\(e(1-1/100)>1+9/40\) to obtain
\(e^{-1/m}<|z|<e^{1/m}\) at every root.

The logarithms of the halves can be chosen to be real on the
positive radius. The sum of the absolute phase slopes at the
radii of any root is at most
\[
             \frac{3}{99}+\frac{7/25}{1-9/40}<\frac12 .
\]
This uses \(|d\arg H(re^{i\theta})/d\theta|\le |wH'(w)/H(w)|\).
At a root of argument \(0<\theta\le\pi\) the equality
\(z^m=-L_m(z)/U_m(1/z)\) therefore gives
\(m\theta=(2j+1)\pi+e_m(\theta)\) with
\(|e_m(\theta)|<\theta/2\).
A common argument in consecutive degrees would give
\(2(j_N-j_n)\pi=\theta+e_n(\theta)-e_N(\theta)\), strictly
between zero and \(2\theta\le2\pi\). Positive coefficients exclude
argument zero, and conjugation covers negative angles.

## The bounded degrees

For lower degree one, the first moment has a negative real root
and the second moment is strictly positive on the real line.
For lower degree two the normalized quadratic has coefficient
\(0<c_{2,1}\le1/2\) and both roots are in the open left half-plane.
The cubic has \(c_{3,1}c_{3,2}\le1/9<1\) from (1). It has a
negative real root \(-r\) and can be written as
\((z+r)(z^2+bz+1/r)\). If \(b\ge0\), then
\(c_{3,1}c_{3,2}=1+b(r^2+1/r)+rb^2\ge1\).
Thus \(b<0\). The quadratic factor has no positive real root,
since \(T_3\) has positive coefficients, so its two roots are
a nonreal conjugate pair in the open right half-plane. Neither
shares an argument with the order-two roots.

For the remaining finite band \(3\le n<40\), the
[interval checker](verify_sfc_asymmetric_private_phase.py)
excludes shared arguments on the entire product of the independent
coefficient envelopes (1). It uses the following compact domain
and separation test.

At either order \(m\in\{n,n+1\}\), the checker verifies
\(\epsilon_m=\sum_{k=1}^{m-1}A_{m,k}<7/10\).
For \(|z|\le1\), a root satisfies \(|1+z^m|\le\epsilon_m\).
For \(|z|\ge1\), apply the same estimate to the reciprocal
polynomial. Consequently, with \(q=m\log|z|\) and
\(s=\arg z/\pi\) after conjugation,
\[
 |q|\le3/2,\qquad
 s\in[0,1]\cap
       \bigcup_{j\ge0}
        \left[\frac{2j+1-1/4}{m},\frac{2j+1+1/4}{m}\right].
\]
The angle estimate uses \(\epsilon_m<\sin(\pi/4)\).
For a rational box of \((s,q)\), the checker encloses
the real and imaginary projections of
\[
 e^{-i\pi k s}
  \left(1+e^{q+i m\pi s}
     +\sum_{h=1}^{m-1} c_{m,h}e^{h q/m+i h\pi s}\right),
                    \qquad 0\le k\le m ,
\]
using the full intervals \(0\le c_{m,h}\le A_{m,h}\).
A strictly signed projection excludes a root in that box.
To certify an angle interval, it partitions the full radius
interval until every part is excluded at one of the two orders.
If a tested order cannot certify the angle interval within its
radial subdivision limit, it tries the other order and then
bisects the angle interval if needed. The starting angle intervals
are the exact intersections of the displayed grids. All interval
evaluations use Arb at 96 bits. Every angle interval is certified,
with 413 angle boxes and 6827 radial boxes visited. The proof
covers the continuum of coefficients and radii at those degrees.

Replay in the repository environment with python-flint installed.

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_sfc_asymmetric_private_phase.py
~~~

## Sharpness for arbitrary positive base distributions

Reserve one exponential for \(U\) and two for \(V\).
For a point \((a,b)\) of the base, the moment divided by \(m!\)
is
\[
        F_m(a,b;z)=\sum_{k=0}^m k!\,a^{m-k}b^k z^k .
\]
First prescribe limiting base laws separately at degrees four
and five. At degree four take the axis atoms \((1,0),(0,1)\)
with unnormalized masses \(1,1/4!\). Their polynomial is
\(1+z^4\), with simple root \(\zeta=e^{i\pi/4}\).
At degree five use those same axis directions and an interior
atom \((1,r)\) with \(r=3/8\).
Put \(H=F_5(1,r;\zeta)\). Directly,
\[
 \Re H=1-24r^4+(r-6r^3-120r^5)/\sqrt2<0,\qquad
 \Im H=2r^2+(r+6r^3-120r^5)/\sqrt2>0 .
\]
For the signs use \(7/10<1/\sqrt2<3/4\); both are checked
rationally in the replay. As \(\zeta^5=-\zeta\), the axis masses
\(\lambda=\Im H-\Re H\) and \(\mu/5!\), where
\(\mu=\sqrt2\,\Im H\), are positive and give
\(\lambda+\mu\zeta^5+H=0\).

To realize a single positive base law, let \(\varepsilon>0\).
Multiply the coordinates of the degree-four atoms by
\(\varepsilon\) and their masses by \(\varepsilon^{-4}\).
Multiply those of the degree-five atoms by \(\varepsilon^{-1}\)
and their masses by \(\varepsilon^5\). Before scaling the
coordinates, replace each zero coordinate by \(\varepsilon\).
Treat the two degree-five axis masses as variables
\(\lambda,\mu/5!\) and keep the interior mass one.
By homogeneity, the fourth moment is the nominal lower
polynomial plus \(\varepsilon\) times the nominal upper
contribution, with all arguments now real analytic in
\(\varepsilon,\lambda,\mu\) near their limiting values.
The fifth moment has the analogous description with the roles
reversed.

The simple fourth-moment root continues real analytically
near \(\zeta\) in the three real parameters. Evaluating the
fifth moment at this branch, its derivative in \((\lambda,\mu)\)
at \(\varepsilon=0\) has columns \(1\) and \(\zeta^5\). They
are linearly independent over the reals. The implicit function
theorem therefore gives positive axis masses for all sufficiently
small positive \(\varepsilon\) with an exact shared nonreal root.
Normalize the five masses by their sum. The resulting base
features are strictly positive and have finite mixed moments
of every order. Reversing sides covers \((2,1)\). Absorbing one
of the reserved exponential factors into the base gives the
same counterexample for \((1,1)\).

For a physical Gaussian path the base law is more restricted.
The counterexample concerns the unrestricted positive-base
class. The sufficient theorem applies to the entire physical
path, including its ordinary endpoints. At those endpoints a
rational Bézout identity transfers the complex coprimality to
every characteristic-zero field.
