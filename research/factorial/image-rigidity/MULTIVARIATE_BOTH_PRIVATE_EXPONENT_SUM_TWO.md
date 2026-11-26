# Angular separation under the reciprocal binomial coefficient bound

For every \(m\ge1\) let
\[
 F_m(z)=1+z^m+\sum_{k=1}^{m-1}c_{m,k}z^k,\qquad
             0\le c_{m,k}\le\binom mk^{-1}.
\]
Then \(F_n\) and \(F_{n+1}\) have disjoint root arguments for every
\(n\ge1\). The arrays at different degrees can be chosen independently
and the coefficients need not be symmetric.

In particular, ordinary factorial binomials satisfy SFC(2) in every
degree if the sum of the private exponents on each side is at least
two. A private row has one zero exponent and one positive exponent.
The remaining shared or private rows may have either orientation.
An optional common Gaussian coefficient factor preserves the conclusion.

Mao, Wang and Zheng prove strict angular interlacing for the
coefficientwise reciprocals of palindromic polynomials with negative
real zeros. Their Newton bound gives a reciprocal binomial coefficient
cap. The argument below applies directly to that cap and compares the
phases at the actual radii, so it also covers asymmetric coefficients.
[Interlacing on the Unit Circle via Coefficientwise Reciprocals](https://arxiv.org/abs/2609.34551)

## A common angular interval for every half

Set \(\alpha=19/50\), \(\beta=3/25\), and \(R_m=e^{7/(10m)}\).
Write
\[
 F_m(z)=L_m(z)+z^m U_m(1/z),
\]
where the two halves give weight \(1/2\) to the middle coefficient
if \(m\) is even. For \(1\le k\le\lfloor m/2\rfloor\) define
\[
 b_{m,k}=\frac{R_m^k}{\binom mk}\begin{cases}
  1/2 & (2k=m),\\
  1 & (2k<m).
 \end{cases}
\]
We will use the three bounds
\[
 A_m=\sum_k b_{m,k}<\tfrac12,\qquad
 \sum_k(k-\alpha)b_{m,k}<\alpha,\qquad
 b_{m,1}+\tfrac{27}{10}\sum_{k\ge2}(k+\beta)b_{m,k}<1.             \tag{1}
\]
For \(m=1\) the sums and \(b_{m,1}\) are zero.

At \(|z|=R_m\), the real coefficients of \(U_m\) let us compare \(L_m(z)\)
with \(U_m(\bar z^{-1})\) on the same ray. Both terms in the coefficient
difference at index \(k\) lie in \([0,b_{m,k}]\) after their respective
radial factors. Hence
\[
 \big||L_m(z)|-|U_m(1/z)|\big|\le A_m,\qquad
 |U_m(1/z)|\ge1-A_m.
\]
The ratio is at most \(1/(1-A_m)<2<R_m^m\).
Rouché gives every root inside \(R_m\), because \(z^m U_m(1/z)\) has all
\(m\) zeros there. Indeed \(U_m\) is nonzero on \(|w|\le R_m\) by
\(A_m<1/2\), and any nonzero root of \(z^m U_m(1/z)\) has reciprocal
outside \(R_m\). Reversal gives every root outside \(R_m^{-1}\).

Now let \(H=L_m\) or \(U_m\). For \(0\le r\le R_m\) and \(0<\theta\le\pi\)
we claim
\[
                  -\beta\theta<\arg H(r e^{i\theta})<\alpha\theta.  \tag{2}
\]
Here \(\Re H>1/2\) by the bound on \(A_m\). For \(s\ge\alpha\) and
\(0<\theta\le\pi\), the decreasing sinc function on \([0,\pi]\)
and \(|\sin x|/x\le1/\pi\) for \(x\ge\pi\) imply
\[
                 |\sin(s\theta)|\le (s/\alpha)\sin(\alpha\theta).
\]
In \(\Im(e^{-i\alpha\theta}H(r e^{i\theta}))\), the constant term is
\(-\sin(\alpha\theta)\) and the rest are bounded above by the same sine
times \(\sum_k(k-\alpha)b_{m,k}/\alpha\). This proves the upper bound.

For the lower bound use \(\Im(e^{i\beta\theta}H(r e^{i\theta}))\) and its
positive constant term \(\sin(\beta\theta)\). The negative part of
\(\sin((1+\beta)\theta)\) is at most \(\sin(\beta\theta)\), because
\((1+\beta)\theta-\pi\le\beta\theta\) when that sine is negative.
For \(k\ge2\), a negative sine requires \(\theta>\pi/(k+\beta)\).
Since \(\beta\theta<\pi/2\), the ratio of its negative part to
\(\sin(\beta\theta)\) is at most \(1/\sin(\pi\beta/(k+\beta))\).
The inequalities \(25/8<\pi<22/7\) and
\(\sin x\ge x(1-x^2/6)\) for \(0\le x\le9/50\) give
\[
             \frac1{\sin(\pi\beta/(k+\beta))}<\frac{27}{10}(k+\beta).
\]
The last inequality in (1) now proves the lower bound in (2).
The rotated arguments used for both bounds lie in \((-\pi,\pi)\).

At any root of argument \(0<\theta<\pi\), the phase equation is
\[
 m\theta-C_m(\theta)\in\pi+2\pi\mathbb Z,\qquad
 C_m(\theta)=\arg L_m(r e^{i\theta})+\arg U_m(r^{-1}e^{i\theta}),
\]
with \(r\) the radius of that root. At a hypothetical common argument,
the two values of \(C_m(\theta)\) both lie in
\((-2\beta\theta,2\alpha\theta)\), an open interval of length \(\theta\).
Subtracting the equations would give an integer multiple of \(2\pi\)
strictly between \(0\) and \(2\theta\). At \(\theta=\pi\) the halves
are real positive, so only odd degree can have a root there. Positive
coefficients at the endpoints exclude \(\theta=0\), and reflection
handles the lower half-plane.

## Bounding the scalar sums in every degree

Put \(B_m=\sum_k k b_{m,k}\). For \(m\ge9\) we have
\(B_m\le\max(B_9,B_{10})<7/25\). This proves (1), since the second
sum is at most \(B_m\), and the third is at most
\((27/10)(1+\beta/2)B_m<1\).

Here is the comparison \(B_{m+2}\le B_m\) in these ranges.
Every term before the old last index decreases, since the radius
decreases and the binomial denominator grows. For even \(m=2h\),
the sum of the last two new terms divided by the old last term
is at most
\[
 \frac{h+2}{2h+1}+\frac{(h+1)^2}{2h(2h+1)}<1\qquad(h\ge5).
\]
The new last term has the same exponent of \(e\) as the old last
term. For odd \(m=2h+1\), the corresponding ratio is at most
\[
 \frac{(h+3)(h+2)}{(2h+3)(2h+2)}
      +\frac{21}{20}\frac{h+1}{h}\frac{h+2}{2(2h+3)}<1
                 \qquad(h\ge2).
\]
Each summand decreases in \(h\) and their sum at \(h=2\) is less
than one. The factor \(21/20\) bounds
\(e^{7/[10(2h+1)(2h+3)]}\le e^{1/50}\).
Thus it suffices to check the two initial weighted sums and (1)
through \(m=8\). The following exact arithmetic verifies these
scalar bounds without constructing any polynomials.

~~~python
from fractions import Fraction as F
from math import comb, factorial

def exp_upper(x):
    return sum(x**j / factorial(j) for j in range(13)) + x**13 / factorial(13) / (1 - x/14)

def caps(m):
    return [
        exp_upper(F(7*k, 10*m)) / comb(m, k) / (2 if 2*k == m else 1)
        for k in range(1, m//2 + 1)
    ]

alpha, beta = F(19, 50), F(3, 25)
assert 1 + F(7, 10) + F(7, 10)**2/2 + F(7, 10)**3/6 > 2
assert exp_upper(F(1, 50)) < F(21, 20)
assert F(10, 21) + F(21, 20)*F(3, 2)*F(2, 7) < 1
assert F(22, 7)*beta/(2 + beta) < F(9, 50)
assert F(3, 8)*(1 - F(9, 50)**2/6) > F(10, 27)
for m in range(2, 9):
    upper = caps(m)
    assert sum(upper) < F(1, 2)
    assert sum((k-alpha)*b for k, b in enumerate(upper, 1)) < alpha
    assert upper[0] + F(27, 10)*sum(
        (k+beta)*b for k, b in enumerate(upper, 1) if k >= 2
    ) < 1
for m in (9, 10):
    assert sum(k*b for k, b in enumerate(caps(m), 1)) < F(7, 25)
assert F(27, 10)*(1 + beta/2)*F(7, 25) < 1
print("Reciprocal binomial half-phase scalar checks passed.")
~~~

## The private exponent condition

Normalize a factorial polynomial by a positive scalar and positive
rescaling of \(z\) so that its endpoint coefficients are both one.
A private row contributing \((d k)!\), for integer \(d\ge1\), has
normalized row factor at most
\[
             \frac{(d k)!}{((d m)!)^{k/m}}
                 \le\left(\frac{k!}{(m!)^{k/m}}\right)^d .
\]
In fact \(T(k)=(d k)!/(k!)^d\) is log-convex in the integer index \(k\),
as its consecutive ratio is
\(\prod_{h=1}^d(d-(d-h)/(k+1))\), a product of positive increasing
functions. Its endpoint-normalized factor is therefore at most one.
A private row in the other orientation has the analogous bound with
\(m-k\). All other normalized factorial row factors are at most one,
by log-convexity of the factorial along the row. The normalized
Gaussian factor is also at most one.

If the two private exponent sums are at least two, multiply the two
private bounds and include the binomial coefficient. This gives
\[
 c_{m,k}\le\binom mk
       \left[\frac{k!(m-k)!}{m!}\right]^2=\binom mk^{-1}.
\]
The angular conclusion applies to the original factorial moments
because the normalizations are positive. Coprimality over
\(\mathbb C\) implies coprimality over every characteristic-zero field,
since the unnormalized ordinary polynomials have integer coefficients.
