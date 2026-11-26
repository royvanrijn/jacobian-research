# All-degree Strong Factorial for a constant and a monomial of total degree four

Let \(M=X_1^{d_1}\cdots X_s^{d_s}\) with positive integer
\(d_i\) and \(\sum_i d_i=4\). Every polynomial
\(R_n(z)=\Lambda((z+M)^n)\) is squarefree. The identity
\(R_n'=nR_{n-1}\) proves consecutive coprimality over every
characteristic-zero field.

Together with the bounds for [degree five](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_FIVE.md),
[degree six](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_SIX.md),
[degrees at least seven](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_SEVEN.md)
and the cases of degree at most three given below, this proves
all-degree SFC(2) for a constant and any nonconstant monomial.
A common monomial factor on disjoint variables changes only
the positive scalar multiplying each moment.

## Endpoint estimates

Use the [coefficient and half-polynomial bounds for degree at least seven](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_SEVEN.md)
with \(D=4\), keeping only their formulas valid for arbitrary \(D\).
Write
\[
 M_k=\prod_i(d_i k)!,\qquad
 T_n(z)=\sum_{k=0}^n c_k z^k,\qquad
 c_k=\binom nk M_k/M_n^{k/n},\qquad
 T_n(z)=L(z)+z^n U(1/z),
\]
including the middle coefficient in \(L\).
Put \(B=c_{n-1}\), \(R=9/8\), \(R'=9/5\).
For the remainder of the analytic proof assume \(n\ge256\).
For a power series \(H\) let
\(\|H\|_r^+=\sum_k|[z^k]H|r^k\) when this sum converges.
Set \(C(n)=\log(2\pi n)/(2n)+1/(12n^2)\).
The explicit right endpoint bound at \(h=1\), using
\(\binom n1=n\), gives
\[
 B\le n e^{-4+4C(n)}<\frac{39}{2000}n,\qquad
                       4C(n)<29/500.                            \tag{1}
\]
The same right endpoint comparison, with \(w=h/n\), gives
\[
 c_{n-h}^{1/h}\le (e/w)e^{4C(n)}(1-w)^{4(1-w)/w}
                                    \quad (1\le h\le n/2).
\]
The right side takes its maximum over \([1/6,1/2]\) at an
endpoint, by logarithmic differentiation as in the degree
at least seven proof. At both endpoints, its product with
\(R\) is less than \(13/25\). Consequently, with \(K=\lfloor n/6\rfloor\),
\[
              \sum_{h>K}c_{n-h}R^h
                 \le (13/25)^{K+1}/(1-13/25),                  \tag{2}
\]
where the sum is restricted to the terms of \(U\).

At the left endpoint, for \(1\le k\le9\) use
\[
                  c_k^{1/k}\le e^4(k!)^{3/k}/n^3<1/1000 .
\]
For \(10\le k\le n/2\), the
[Stirling calculation for degree five](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_FIVE.md#the-left-half)
applies with the entropy term
\(3\log w-(1-w)\log(1-w)/w\le-2\log2\), \(w=k/n\).
Indeed its derivative has the sign of
\(4+\log(1-w)/w>0\) on \(0<w\le1/2\).
The remaining prefactors give
\[
 c_k^{1/k}
 \le \frac14 \exp\left(\frac{3\log(2\pi k)+\log2}{2k}
                            +\frac9{24k^2}\right)<49/100 .
\]
For the last inequality the exponent decreases for \(k\ge10\).
The two geometric sums at radius \(R\) give
\(\|L-1\|_R^+<1/100\), and hence
\(\sup_{|z|\le R}|\log L(z)|<1/50\).

## The logarithm at the monomial endpoint

As in the degree at least seven proof, interpolate
\(F(h)=\log(c_{n-h}h!/B^h)\) by Gamma functions on \(0\le h\le K\).
It satisfies \(F(0)=F(1)=0\). For each \(2\le k\le K\)
put \(f=k/n\), and define
\[
 \alpha_k=\frac{1}{n(1-f)},\qquad
 H_k=\frac{3}{n(1-f)^2},\qquad \beta_k=e^{H_k}H_k,\qquad
 u_k=\alpha_k k=\frac{f}{1-f}.
\]
For \(2\le r\le k\), \(0\le x\le k\), the differentiated
Gamma series has three terms per integer shift after
cancelling the binomial term exactly:
\[
 \frac{F^{(r)}(x)}{(r-1)!}
   =\sum_{j\ge0}\left[
      \sum_{i=1}^s\sum_{a=1}^{d_i}(n-x+j+a/d_i)^{-r}
                           -(n-x+j+1)^{-r}\right].
\]
Each row has its \(a=d_i\) term at shift \(1\), so the
subtraction cancels one such term. The remaining shifts are
positive. Bound each remaining series by its initial term
plus its integral, replacing \(n-x+a/d_i\) by \(n-k\).
Using \((r-1)/(n-k)\le k/(n-k)\) gives
\[
        0\le F^{(r)}(x)\le (r-2)!\alpha_k^{\,r-2}H_k \le H_k .
\]
In the last step, \((r-2)\alpha_k\le k/(n-k)\le1/5\).
The exact integral of \(F^{(r)}\) over the unit cube gives
\(\Delta^r F(0)\). Thus the hyperedge weights satisfy
\[
    0\le e^{\Delta^r F(0)}-1
                       \le (r-2)!\alpha_k^{\,r-2}\beta_k.
\]

Let \(\ell_k=[z^k]\log U(z)\) in the formal expansion.
Apply the [sharper discovery-tree count](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#sharper-connected-sum)
separately at each size \(k\). Only its absolute bound is used.
For clarity the two constants in its proof change as follows.
For the sum of all possible additional edges on \(k\) labels, use
\[
 \sum_{r=2}^k
       \binom kr(r-2)!\alpha_k^{r-2}\beta_k
 \le \beta_k k(k-1)\left(\frac12+\frac{u_k}{6(1-u_k)}\right).
\]
The discovery tree has the base
\((1+t)^{(1+t)/t}\), \(t=\alpha_k/\beta_k\le1/3\).
It is at most \((4/3)^4=256/81\), since its logarithmic
derivative has the sign of \(t-\log(1+t)>0\).
The tree prefactor is at most \(e/k\le1\) when \(k\ge3\).
For \(k=2\) the single edge gives the same majorant directly.
Thus \(\ell_1=B\) and, for \(2\le k\le K\),
\[
 |\ell_k|\le B p_k^{k-1}/k,\qquad
 p_k=\frac{256 B\beta_k}{81(1-u_k)}
        \exp\left(\beta_k k
                  \left(\frac12+\frac{u_k}{6(1-u_k)}\right)\right).
                                                                  \tag{3}
\]
This bound uses the actual number \(k\) of vertices throughout.

Using (1) in (3) gives the following bounds. In each interval
use the indicated upper value of \(k\) to bound \(k/n\), clipped
at \(1/6\), and use \(n\ge256\).
\[
 \begin{array}{c|rrrrrrr}
 k & 2\text{--}3 &4\text{--}6 &7\text{--}10&11\text{--}15&
                 16\text{--}23 &24\text{--}32 &33\text{ and up}\\ \hline
 p_k\ \text{at most}
 &1/5&21/100&23/100&13/50&31/100&19/50&51/100
 \end{array}
\]
For \(A(z)=\sum_{k=1}^K\ell_k z^k\), sum these scalar series to get
\[
             \|A\|_R^+<\frac{32}{25}B,\qquad
             \|A\|_{R'}^+<\frac{23}{10}B.                         \tag{4}
\]
The polynomial \(U\) and \(e^A\) agree through degree \(K\).
For the tail of \(e^A\), its coefficient norm on \(R'\) is at
most \(e^{\|A\|_{R'}^+}\). The tail on \(R\) therefore costs
the factor \((R/R')^{K+1}\). With (2) and (4) this proves
\[
 \sup_{|z|\le R}|U(z)e^{-A(z)}-1|
 \le (5/8)^{K+1}e^{(32/25+23/10)B}
          +\frac{(13/25)^{K+1}}{1-13/25}e^{(32/25)B}
 < \frac3{25}+\frac1{200}=\frac18 .
\]
Here \(K+1\ge n/6\), \(\log(8/5)>47/100\) and
\(\log(25/13)>13/20\). The two differences in exponential
rates are positive and the rational replay checks \(n=256\).
Hence \(U\) is zero-free on the disk and
\[
    \sup_{|z|\le R} \left(|\log U(z)|+|\log L(z)|\right)
      <(32/25)B+7/50+1/50<13n/500 .
\]
The same bound holds for the sum of the separate sup norms.
At \(R_0=17/16\), one has
\(\log R_0>1/17>13/500\) and \(R_0/(R-R_0)=17\).
The [two-half criterion](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_SEVEN.md#a-criterion-using-the-two-halves)
now puts each root in \(R_0^{-1}<|z|<R_0\).
By Cauchy's estimate the sum of the two logarithmic derivatives
in that annulus has modulus less than \(17(13n/500)<n\).
The same criterion gives simplicity at every root.

## Finite range and the complete constant support

The [grouped analytic comparison for orders fourteen through 256](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_FOUR_MIDDLE.md)
controls the endpoint logarithm after removing its first
exponential term. It handles all five partitions simultaneously.
For \(1\le n\le13\), check those five partitions modulo
the prime \(1021\) using the short independent replay below.
The polynomials \(R_n\) are monic and \(n<1021\), so a
modular gcd of degree zero with the derivative proves
squarefreeness in characteristic zero.

For total degree one, \(R_n=n!\sum_{k=0}^n z^k/k!\)
satisfies \(R_n-R_n'=z^n\), so no nonzero root is multiple.
For total degree two, the \((1,1)\) profile is covered by
the [two-variable product theorem](CONSTANT_PLUS_PRODUCT_OF_TWO_VARIABLES_ALL_ORDERS.md).
For the \((2)\) profile, reverse \(R_n\) and set \(t=y/2\).
Its coefficient at \(y^k\) is the integer
\[
           \frac{n!}{(n-k)!}\frac{(2k)!}{2^k k!}.
\]
The second factor is odd. Write \(n=\sum_i a_i2^i\).
The \(2\)-adic Newton polygon has one edge of width \(2^i\)
and slope \(1-2^{-i}\) for each bit \(a_i=1\), in increasing
order of \(i\). Indeed on \(k=s_i-\Delta\),
\(s_i=\sum_{j\le i}a_j2^j\), \(0\le\Delta\le2^i\), the
ordinate minus the edge chord is
\(\operatorname{digitsum}_2(\Delta)-\Delta/2^i\ge0\), by
Legendre's formula and the disjoint bit positions of \(n-s_i\)
and \(\Delta\). A repeated irreducible factor over \(\mathbb Q_2\)
would occupy at least twice the denominator of one of the
slopes in horizontal width, by additivity of Newton polygons
and the integrality of the factor's slope endpoint ordinates.
Each edge has width exactly its slope denominator. This
excludes all repeated factors. For total degree three use
the [all-degree cubic result](CONSTANT_PLUS_TOTAL_DEGREE_THREE_ALL_ORDERS.md).

## Exact replays

The scalar inequalities are checked with rational Taylor bounds.

~~~python
from fractions import Fraction as Q
from math import factorial

def low_exp(x):
    x = Q(x)
    return sum(x**k / factorial(k) for k in range(25))

def high_exp(x):
    x = Q(x)
    assert 0 <= x < 26
    return low_exp(x) + x**25 / (factorial(25) * (1-x/26))

def geometric_tail(ratio, start):
    return ratio**start / (1-ratio)

n0 = 256
R, outer = Q(9, 8), Q(9, 5)
b = Q(39, 2000)
delta = Q(29, 500)
assert low_exp(Q(37,5)) > 2*Q(22,7)*n0
assert 2*Q(37,5)/n0 + Q(1,3*n0**2) < delta
assert high_exp(delta) < b*low_exp(4)
assert high_exp(4) < 55 and factorial(9) < 72**3
assert Q(55*72,n0**3) < Q(1,1000)
assert low_exp(Q(83,20)) > 20*Q(22,7) and low_exp(Q(7,10))>2
assert high_exp((3*Q(83,20)+Q(7,10))/20 + Q(9,2400))/4 < Q(49,100)
assert geometric_tail(R/1000,1)+geometric_tail(R*Q(49,100),10)<Q(1,100)
assert low_exp(Q(1,50)) > Q(100,99)
groups = [
    (3,Q(1,5)), (6,Q(21,100)), (10,Q(23,100)),
    (15,Q(13,50)), (23,Q(31,100)), (32,Q(19,50)),
    (48,Q(51,100)),
]
for k, bound in groups:
    f = min(Q(k,n0),Q(1,6))
    u = f/(1-f)
    h = 3/(n0*(1-f)**2)
    beta = 3*high_exp(h)/(1-f)**2
    psi = Q(1,2)+u/(6*(1-u))
    p = b*Q(256,81)*beta/(1-u)*high_exp(beta*f*psi)
    assert p < bound

def norm(radius):
    start, total = 1, Q(0)
    for stop, p in groups:
        total += sum(radius*(p*radius)**(k-1)/k for k in range(start,stop+1))
        start = stop+1
    return total+radius*geometric_tail(groups[-1][1]*radius,start-1)/start

assert norm(R) < Q(32,25)
assert norm(outer) < Q(23,10)
v=Q(13,25)
assert R*6*high_exp(1+delta)*Q(5,6)**20 < v
assert R*2*high_exp(1+delta)/16 < v
assert high_exp(Q(47,100))<outer/R
assert high_exp(Q(13,20))<1/v
loss=Q(47,600)-b*(Q(32,25)+Q(23,10))
assert loss>0 and low_exp(loss*n0)>Q(25,3)
loss=Q(13,120)-b*Q(32,25)
assert loss>0 and (1-v)*low_exp(loss*n0)>200
assert Q(3,25)+Q(1,200)<Q(1,8)+Q(1,1000)
assert Q(3,25)+Q(1,200)==Q(1,8)
assert low_exp(Q(7,50))>Q(8,7)
assert b*Q(32,25)+(Q(7,50)+Q(1,50))/n0 < Q(13,500)
assert Q(13,500)<Q(1,17)
~~~

The finite modular check uses ascending coefficient lists.

~~~python
from math import comb, isqrt, prod

PARTITIONS = ((4,), (3,1), (2,2), (2,1,1), (1,1,1,1))

def remainder(a, b, p):
    a = list(a)
    last_inverse = pow(b[-1], -1, p)
    while len(a) >= len(b):
        shift = len(a) - len(b)
        lead = a.pop() * last_inverse % p
        if lead:
            for i in range(len(b)-1):
                a[shift+i] = (a[shift+i] - lead*b[i]) % p
        while a and a[-1] == 0:
            a.pop()
    return a

def squarefree(f, p):
    derivative = [k*f[k] % p for k in range(1,len(f))]
    while derivative:
        f, derivative = derivative, remainder(f,derivative,p)
    return len(f)==1

p=1021
assert all(p%i for i in range(2,isqrt(p)+1)) and p>4*13
factorials = [1]
for m in range(1,4*13+1):
    factorials.append(m*factorials[-1]%p)
for d in PARTITIONS:
    for n in range(1,14):
        coefficients = [
            comb(n,k)*prod(factorials[h*k] for h in d)%p
            for k in range(n,-1,-1)
        ]
        assert squarefree(coefficients,p), (d,n)
~~~
