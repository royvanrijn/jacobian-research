# Angular separation with curvature on doubled shared exponents

Let \(S\) be the absolute curvature of an ordinary two-monomial
profile, with optional common Gaussian strength \(\tau\ge0\)
included as additional curvature. Let \(S_H\) be the part of \(S\)
from shared unequal rows whose larger exponent is at least twice
their smaller exponent. Put \(S_L=S-S_H\), including private rows
and \(\tau\) in \(S_L\). Define
\[
                         B=\frac{61}{100}S_H+\frac12 S_L .
\]
If
\[
       S\ge14,\qquad (1+2S)e^{-B}\le29 e^{-427/50},             \tag{1}
\]
then the consecutive moment polynomials have disjoint root
arguments in every degree. At \(\tau=0\) this gives ordinary
SFC(2) over every characteristic-zero field. There is no
orientation condition.

In particular it suffices that the support is fully shared,
every unequal row has its larger exponent at least twice its
smaller exponent, and \(S\ge14\). Condition (1) also permits
shared rows with smaller relative gaps and private rows.
It uses the actual curvature proportions. The proof changes
the entropy input to the
[absolute curvature sixteen comparison](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md)
while keeping its short and long phase bounds.

## Entropy budgets

Retain every shared row in the
[oriented endpoint criterion](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md).
Write \(q=\max Q_\pm(0)\), \(V\) for the Gamma curvature,
\(W=1+V\), and \(L_0\) for the maximum of one and of the
shared relative gaps. An ordinary shared row of smaller
exponent \(\ell\) and relative gap \(r=d/\ell\) has
\(S_i=\ell r^2/(1+r)\), \(V_i=\ell r^2\).
Its smaller oriented endpoint cost \(-\log Q_i(0)\)
is \(\ell(r-\log(1+r))\), in the decreasing orientation.

For a row counted in \(S_H\), put \(c=61/100\) and
\[
 E_c(r)=r-\log(1+r)-\frac{c r^2}{1+r}.
\]
Let \(C=1+2S\ge29\). For \(r\ge1\),
\[
 E_c(r)\ \ge\ \log\left(1+\frac{(r-1)r^2}{C(1+r)}\right).      \tag{2}
\]
At \(r=1\) this follows from \(\log2<139/200\).
For derivatives, cancel \(r/(1+r)^2\) on both sides.
The inequality reduces to
\[
  [(1-c)r+1-2c]\,[C+(r-1)r^2/(1+r)]\ge2(r^2+r-1).
\]
Multiply by \(1+r\), and set \(t=r-1\ge0\).
The left side minus the right side is at least its
value for \(C=29\), since \((1-c)r+1-2c>0\) on \(r\ge1\).
That value is
\[
       t^2\left(\frac{39}{100}t^2-\frac{21}{20}t+
                 \frac{51}{25}\right)+\frac{343}{25}t+
                 \frac{293}{50}.
\]
The quadratic in parentheses has negative discriminant
and positive leading coefficient. This proves (2).

Set \(U_i=V_i-2S_i\ge0\) on the rows counted in \(S_H\).
Multiplying (2) by \(\ell\) gives a cost beyond \(cS_i\)
at least \(\log(1+U_i/C)\), by
\((1+x)^\ell\ge1+\ell x\).
For the other rows the oriented cost is at least \(S_i/2\)
and \(V_i\le2S_i\). This includes private rows and \(\tau\).
Summing the costs and multiplying the elementary factors gives
\[
                q\le e^{-B},\qquad Wq\le C e^{-B}.             \tag{3}
\]
The largest relative gap also has a separate bound
\[
                              L_0q\le5e^{-B}.                 \tag{4}
\]
If \(L_0=1\), this follows from (3).
Otherwise take a row attaining \(r=L_0\). Its cost beyond
\(cS_i\) is at least \(E_c(r)\ge0\), since \(\ell\ge1\).
As \(E_c(r)\ge(1-c)r-\log(1+r)\), it suffices to bound
\(r(1+r)e^{-(39/100)r}\) by five. The maximum is attained
at the positive root of
\((39/100)r^2-(161/100)r-1=0\), between \(23/5\) and \(24/5\).
Substituting the upper bound in the polynomial factor and
the lower bound in the exponential gives a value below five.

## The midpoint gain

For a shared row with \(r\ge1\) the midpoint cost
\(-\log Q_i(1/2)\) is at least \(7S_i/25\).
Indeed it suffices to take \(\ell=1\) after division by \(\ell\),
since the subtracted concave correction is nonnegative
and independent of \(\ell\). The cost at \(\ell=1\) is
\[
 f(r)=(r+5/4)\log(1+r)-(r+5/2)\log(1+r/2).
\]
At \(r=1\) it exceeds \(7/50\).
Also
\[
 f'(r)=\log\frac{2(1+r)}{2+r}
                     -\frac{r}{4(1+r)(2+r)}
\]
increases on \(r\ge1\), by differentiating.
On \(1\le r\le3/2\), it exceeds \(147/625\) by its
value at one, while the derivative of \(7r^2/[25(1+r)]\)
is at most \(147/625\).
On \(r\ge3/2\), it exceeds \(7/25\) by its value at
\(3/2\), while the latter derivative is below \(7/25\).
This proves the claimed cost. Every other row contributes
at least \(S_i/4\) to the midpoint cost, by the
[full retention chord estimate](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md),
including the direct curvature factor for private rows
and for \(\tau\). Hence
\[
 -\log Q_\pm(1/2)\ge\frac7{25}S_H+\frac14 S_L
                        =\frac3{11}B+\frac5{44}S.
\]
Since \(S\ge14\), (1) implies \(B\ge427/50\).
The displayed midpoint cost is therefore at least \(98/25\).

## Reduction to the short and long phase comparisons

Assume \(14\le S<16\), since the existing curvature-sixteen
result gives the remaining stronger curvature range.
Then \(V\le S(S+2)<288\), and set \(L=3V/4\).
The original oriented increment bound is
\(\log Q_\pm(w)\le\log q+Lw\).
Put \(u=e^{-427/50}\).
At \(m=256\), for \(j=1,2\), use
\[
 Q_\pm(j/256)
 \le u\min\!\left(1,\frac{29}{1+V}\right)e^{3jV/1024}.
\]
The second line of (3) gives the second factor in the minimum.
For \(0\le V\le28\) the maximum of this envelope is at
\(V=28\). For \(28\le V\le288\) the logarithmic derivative
of the envelope is \(3j/1024-1/(1+V)\), increasing in \(V\),
so the maximum is at \(28\) or \(288\).
These two values give \(256Q_\pm(1/256)<1/16\) and
\(128eQ_\pm(1/128)<1/8\).

For \(2\le m\le256\) the \(k=1\) estimate follows by
the maximum of \(mQ_\pm(1/m)\) at \(m=2,256\),
using \(2Q_\pm(1/2)<1/16\) for the first endpoint.
For \(k\ge2\), the bound \(c_{m,k}^{1/k}\le eQ_\pm(w)/w\)
maximizes at \(w=2/m\) or \(w=1/2\).
At the former use the maximum of \((m/2)Q_\pm(2/m)\)
at \(m=4,256\), by the same logarithmic derivative
argument as in the oriented criterion.
At the midpoint use \(2eQ_\pm(1/2)<1/8\).
This gives the short annulus and phase bounds from
the curvature-sixteen proof.

For \(N=n+1>256\) and \(D=Nq\le4/25\), we have
\(L/n<89/100\) and
\[
 \frac{LD}{n}
     \le\frac34\frac{257}{256}\,29u < \frac{43}{10000}.
\]
These are the identical short-half bounds of the
[curvature sixteen proof](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md#the-short-ranges).
Its first three coefficient estimates and the cap
\(8^{-k}\) for the remaining terms follow, using
\(2eQ_\pm(1/2)<1/8\) for its midpoint test.
Thus the same short phase contradiction holds.

Finally if \(D>4/25\), then \(n\ge800\) since \(q<1/5000\).
Equations (3) and (4), with (1), imply
\[
  D/n<337/10^6,\qquad WD/n<572/10^5,\qquad
                                L_0D/n<101/10^5 .
\]
They are the three inputs of the
[long curvature-sixteen argument](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md#the-long-range-with-sharper-gamma-derivative-estimates).
Here \(\log N>6\) and \(S<16<2\log N+5\), as required
for its first-factor derivative estimate. The midpoint
bound just proved suffices for the long geometric tail.
All Gamma interpolation, connected logarithm and paired
phase inequalities of that argument apply with the
same constants. Terminal two separates by positivity
of the even quadratic on the real axis.

## Exact scalar checks

~~~python
from fractions import Fraction as F
from math import factorial

def lo(x):
    return sum(F(x)**j/factorial(j) for j in range(81))

def hi(x):
    x=F(x)
    return lo(x)+x**81/factorial(81)/(1-x/82)

assert lo(F(139,200))>2
assert (F(21,20)**2)<4*F(39,100)*F(51,25)
a=F(39,100)
assert a*F(23,5)**2+(a-2)*F(23,5)-1<0
assert a*F(24,5)**2+(a-2)*F(24,5)-1>0
assert F(24,5)*F(29,5)/lo(a*F(23,5))<5
assert hi(F(14,25))*F(3,2)**14<2**9
assert hi(F(147,625)+F(1,24))<F(4,3)
assert hi(F(7,25)+F(3,70))<F(10,7)
u=1/lo(F(427,50))
assert u<F(1,5000)
assert 2/lo(F(98,25))<F(1,16)
assert 2*hi(1)/lo(F(98,25))<F(1,8)
for v in [28,288]:
    assert 256*u*F(29,1+v)*hi(F(3*v,1024))<F(1,16)
    assert 128*hi(1)*u*F(29,1+v)*hi(F(6*v,1024))<F(1,8)
assert F(216,256)<F(89,100)
assert F(3,4)*F(257,256)*29*u<F(43,10000)
assert F(801,800)*u<F(337,10**6)
assert F(801,800)*29*u<F(572,10**5)
assert F(801,800)*5*u<F(101,10**5)
assert hi(6)<801
~~~
