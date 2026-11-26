# Both-private separation at residual curvature seven halves

Consider an ordinary two-monomial support with a private
variable on each side. Allow arbitrary shared row orientations
and set
\[
                     S_*=\sum_{\max(a_i,b_i)>0}
                              \frac{(b_i-a_i)^2}{\max(a_i,b_i)}-2.
\]
If \(S_*\ge7/2\), the consecutive factorial moments have
disjoint root arguments in every degree. In particular the
profile satisfies SFC(2) over every characteristic-zero field.
For upper order \(m\ge14\), the arguments strictly interlace
when listed from the positive real axis. This is total
ordinary curvature at least \(11/2\).

The proof sharpens the geometric endpoint remainder of the
[residual-curvature-four argument](MULTIVARIATE_BOTH_PRIVATE_RESIDUAL_CURVATURE_FOUR.md).
Terminals through thirteen are already covered by the
[complete ordinary results](MULTIVARIATE_ALL_DEGREE_REDUCTION.md).
The higher-order argument uses no tests of individual
moment polynomials.

## A narrow annulus from the private reserve

Reserve one factorial on each side and normalize both
endpoints to one, as in the
[both-private angle-box theorem](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md).
For the normalized polynomial \(P_m(z)=\sum_{k=0}^m c_{m,k}z^k\),
the [residual endpoint bound](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_CURVATURE_INTERLACING.md#a-common-radial-bound)
gives
\[
                    c_{m,0}=c_{m,m}=1,\qquad
                 0<c_{m,k}\le e^{-7k(m-k)/(4m)}.              \tag{1}
\]
For \(m\ge14\), pair \(k,m-k\) on \(|z|=e^{2/(3m)}\).
The paired radial multiplier is at most \(1+e^{2/3}<3\),
with only half that bound needed for the even middle term.
For the first five coefficient envelopes use \(m\ge14\)
and for the rest use \(k\le m/2\). Their sum is at most
\[
       \sum_{k=1}^5 \exp[-(7k/4)(1-k/14)]
                                      +\frac{e^{-42/8}}{1-e^{-7/8}}
                                  <\frac3{10}.
\]
Thus the sum of the interior terms has modulus below
\(9/10<e^{2/3}-1\). Rouché with \(1+z^m\) and then with
the reversed polynomial gives
\[
                         e^{-2/(3m)}<|z|<e^{2/(3m)}
                                                               \tag{2}
\]
at every root.

## The signed geometric remainder

At a root \(z=r e^{i\theta}\) use the
[two endpoint resolvents](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry).
For either one let \(s\) be its positive endpoint ratio
including the radius \(r\) or \(1/r\), and write
\[
      c=\cos\theta,\quad q=\sin\theta,\quad D=1-2cs+s^2,\quad
      A=\mathbb E[1/D],\quad B=\mathbb E[s/D],\quad
      u_j=\mathbb E s^j .
\]
By (1) and (2), both endpoints obey
\[
                    u_j\le a_m(j)=
                        \exp\left[-\frac{7j}{4}
                                     +\frac{7j^2/4+2j/3}{m}\right].
                                                               \tag{3}
\]
For \(m\ge14\) set \(\theta_0=\pi/(m-1)\). We bound
the two phases in the range \(\theta_0\le\theta\le\pi/2\).

For an integer \(L\) between three and \(\lfloor m/2\rfloor\),
take the imaginary part of the finite expansion of
\((1-se^{i\theta})^{-1}\) with its exact geometric remainder.
The contribution of the remainder to \(B\) is
\[
       \mathbb E\!\left[
       s^L\,\frac{\sin(L\theta)-s\sin((L-1)\theta)}{qD}
       \right]
       \le \frac{1+\cos((L-1)\theta)}{2q^2}\,u_L .             \tag{4}
\]
For the inequality set \(\beta=(L-1)\theta\) and
\(x=(s-c)/q\). The scalar quotient inside the expectation,
excluding \(s^L\), becomes
\((\cos\beta-x\sin\beta)/(q^2(1+x^2))\).
Its claimed upper bound follows by completing the square
in
\[
   \frac{1+\cos\beta}{2}(1+x^2)-(\cos\beta-x\sin\beta).
\]

Using \(A\ge1+2cu_1-u_2\) and (4) gives
\[
 \begin{split}
 (1+2c)B-A+1
  \le{}&u_1+\bigl(1+2c(1+2c)\bigr)u_2\\
       &+(1+2c)\sum_{j=3}^{L-1}
                             \frac{\sin(j\theta)}{\sin\theta}u_j
       +\frac{(1+2c)(1+\cos((L-1)\theta))}{2q^2}\,u_L .
 \end{split}                                                    \tag{5}
\]
Each scalar coefficient on the right has its maximum in
the stated angle range at \(\theta_0\). For the sine ratio,
\(\sin(j\theta)/\sin\theta\) decreases on
\((0,\pi/j)\), by differentiating and using the increase
of \(\tan x/x\) on \((0,\pi/2)\) before \(j\theta=\pi/2\).
Beyond \(\pi/j\) the ratio is at most \(1/\sin(\pi/j)\),
which is no larger than its value at \(\pi/(2j)\).
We have \(\theta_0\le\pi/(2j)\). Multiplication by the
positive decreasing \(1+2c\) preserves the bound.
For the last coefficient the claim follows by
monotonicity through \(2\theta_0\), since
\(2(L-1)\theta_0\le\pi\). Beyond \(2\theta_0\) the
[endpoint trigonometric bound](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_RESOLVENT_INTERLACING.md#bounding-the-resolvent-with-endpoint-moments)
bounds the coefficient by \(3/(4\theta_0^2)\).
At \(\theta_0\) it is at least \(1/\theta_0^2\), using
\(c\ge1/2\), \(\sin\theta_0\le\theta_0\), and
\((L-1)\theta_0\le\pi/2\).

## Uniform scalar separation

For \(14\le m\le19\), the following scalar upper budgets
bound (5) after substituting \(\theta_0\) and (3).
They use only the curvature envelope (1).
\[
\begin{array}{c|rrrrrr}
 m&14&15&16&17&18&19\\ \hline
 L&6&6&7&7&8&8\\
 \text{upper budget}&49/50&93/100&88/100&84/100&4/5&77/100
\end{array}
\]
The rational replay below verifies every inequality.

For \(m\ge20\), let \(L=\lfloor m/2\rfloor\) and
use \(q_0=e^{-7/8}\).
The first five terms of (5) are bounded by
\[
              a_{20}(1)+7a_{20}(2)+9a_{20}(3)
                        +12a_{20}(4)+15a_{20}(5)<71/100 .
\]
By (3), the sum of the terms with \(6\le j<L\) is at
most \(3e^{1/3}\sum_{j\ge6}j q_0^j<13/50\).
In the last term, the cosine factor evaluated at
\(\theta_0\) is at most \(5/8\), since
\((L-1)/(m-1)\ge3/7\) and
\(\cos(3\pi/7)=\sin(\pi/14)<1/4\).
Therefore the last term is at most
\[
                    \frac56 L^2 e^{1/3-7L/8}
                          \le \frac{250}{3}e^{-101/12}<1/50.
\]
Here use \((m-1)^2/3\) for the trigonometric bound
before its cosine factor and \(m-1\le2L\).
The sequence in \(L\) decreases for \(L\ge10\).
The sum of the three budgets is \(99/100\).

Thus \((1+2c)B<A\) at both endpoints throughout the
required angles. By the
[exact half-phase comparison](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_CURVATURE_INTERLACING.md#strict-adjacent-angular-order),
the corresponding endpoint phases are both below
\(\theta/2\). At an upper-order root with label \(k\),
the angle is then strictly below
\((2k-1)\pi/(m-1)\) whenever it lies between
\(\theta_0\) and \(\pi/2\). Below \(\theta_0\)
this bound is automatic. Above \(\pi/2\) the upper
angle box already places the angle below that lower
left endpoint. The angle boxes yield the other side
of interlacing directly. This proves the stated
higher-order strict angular separation. The complete
earlier terminal results and rational Bézout give
the all-degree ordinary conclusion.

## Exact replay of the scalar comparisons

The code uses positive exponential Taylor bounds
and alternating trigonometric Taylor bounds on the
small ranges shown. The rational bounds
\(333/106<\pi<22/7\) suffice. The lower bound is
also verified by the Machin formula as in the earlier
residual-curvature-four proof.

~~~python
from fractions import Fraction as F
from math import factorial

def lo(x):
    return sum(F(x)**j/factorial(j) for j in range(25))

def hi(x):
    x = F(x)
    return lo(x) + x**25/factorial(25)/(1-x/26)

def neg(x):
    return 1/lo(x)

def a(m,j):
    return neg(F(7*j,4)-F(F(7*j*j,4)+F(2*j,3),m))

def cosup(x):
    return 1-x*x/2+x**4/24

def sinup(x):
    return x-x**3/6+x**5/120

q = neg(F(7,8))
annulus = sum(neg(F(7*j,4)*(1-F(j,14))) for j in range(1,6))
assert annulus+q**6/(1-q) < F(3,10)
assert hi(F(2,3)) < 2 and lo(F(2,3)) > F(19,10)
x = F(1,5)
assert 4*(4*(x-x**3/3+x**5/5-x**7/7)-F(1,239)) > F(333,106)
assert F(22,7)/14 < F(1,4)

for m,L,cap in [
    (14,6,F(49,50)), (15,6,F(93,100)), (16,7,F(88,100)),
    (17,7,F(84,100)), (18,8,F(4,5)), (19,8,F(77,100)),
]:
    x = F(333,106*(m-1))
    y = F(22,7*(m-1))
    c = cosup(x)
    s = x-x**3/6
    head = a(m,1)+(1+2*c*(1+2*c))*a(m,2)
    head += (1+2*c)*sum(sinup(j*y)/s*a(m,j) for j in range(3,L))
    error = (1+2*c)/s**2*(1+cosup((L-1)*x))/2*a(m,L)
    assert head+error < cap

assert a(20,1)+7*a(20,2)+9*a(20,3)+12*a(20,4)+15*a(20,5) < F(71,100)
assert 3*hi(F(1,3))*q**6*(6-5*q)/(1-q)**2 < F(13,50)
assert F(250,3)*neg(F(101,12)) < F(1,50)
assert F(11,10)**2*q < 1
~~~
