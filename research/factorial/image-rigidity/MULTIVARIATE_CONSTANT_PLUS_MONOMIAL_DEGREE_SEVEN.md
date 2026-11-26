# All-degree Strong Factorial for a constant and a monomial of total degree at least seven

Let \(d_1,\ldots,d_s\) be positive integers, \(D=\sum_i d_i\ge7\),
and let \(\Lambda\) denote the ordinary multivariate factorial
functional. Then
\[
                   R_n(z)=\Lambda((z+X_1^{d_1}\cdots X_s^{d_s})^n)
\]
is squarefree for every \(n\ge1\). In particular \(R_n\) and
\(R_{n-1}\) are coprime because \(R_n'=nR_{n-1}\). This gives
SFC(2) for the support over every characteristic-zero field.
Shared equal-exponent coordinates in independent variables can
be removed as scalar factors. The theorem is uniform in the
number of factors and in the partition of \(D\).

The argument uses a mostly constant half at the constant-monomial
endpoint. At the opposite endpoint, a stronger factorial estimate
gives an effective first coefficient of size \(ne^{-D}\).
A direct comparison or a short logarithm approximation then makes
the logarithmic derivative of the dominant monomial larger than
both endpoint derivatives. The proof has no new tests of individual
factorial moment polynomials.

## An exact coefficient comparison

Put \(M_k=\prod_i(d_i k)!\) and write the endpoint-normalized
moment as
\[
             T_n(z)=\sum_{k=0}^n c_k z^k,\qquad
             c_k=\binom nk M_k/M_n^{k/n}.
\]
Reversal and positive scaling relate this polynomial to \(R_n\).
Write \(f_{n,k}=k!/(n!)^{k/n}\). For every positive integer \(d\),
the sequence \((dk)!/(k!)^d\) is log convex in \(k\) and starts at
one. Its successive ratios are products of the increasing factors
\(d-(d-h)/(k+1)\), \(1\le h\le d\). Thus
\[
                         c_k\le \binom nk f_{n,k}^{\,D}.         \tag{1}
\]
In particular, lowering \(D\) in this upper bound increases it.

Set \(C(x)=\log(2\pi x)/(2x)+1/(12x^2)\) for \(x>0\).
[Robbins's Stirling bounds](https://heuklyd.github.io/papers/pdf/Robbins-1955.pdf)
give
\(\log n!^{1/n}\le\log n-1+C(n)\) and \(n!^{1/n}\ge n/e\).
For \(1\le k\le n/2\), \(k\ge6\), use \(C(k)\le C(6)<31/100\)
to obtain
\[
       c_k^{1/k}\le e(k/n)^{D-1}e^{DC(k)}
                       \le e^{1+7(31/100)}/64<3/8.
\]
The bound decreases in \(D\ge7\), because
\((k/n)e^{C(k)}<1\). For \(1\le k\le5\), \(n\ge14\),
the arithmetic–geometric mean gives \(k!^{1/k}\le(k+1)/2\).
Thus
\[
       c_k^{1/k}\le(en/k)(e(k+1)/(2n))^7
                        \le (14e/5)(3e/14)^7<3/8,
\]
since this upper bound increases in \(k\ge1\) at fixed \(n\)
and decreases in \(n\) at fixed \(k\).

At the other end put \(T=ne^{-D}\) and \(\delta=D C(n)\).
For \(1\le h\le n/2\), let \(w=h/n\). We have
\[
 \begin{split}
   c_{n-h}^{1/h}&\le \frac e w\,e^\delta (1-w)^{D(1-w)/w},\\
   c_{n-h}&\le \frac{(Te^\delta)^h}{h!}
                               e^{D h(h-1)/n}.                 \tag{2}
 \end{split}
\]
Indeed \(f_{n,n-h}=(n!)^{h/n}/\prod_{i=0}^{h-1}(n-i)\).
Bound the sum of the \(\log(1-i/n)\) below by its integral on
\([0,w]\) for the first inequality, and use
\(\log(1-i/n)\ge-2i/n\) for the second.
Use \(\binom n h\le(en/h)^h\) or \(\binom n h\le n^h/h!\),
respectively. The first bound in (2) takes its maximum on any
closed interval of \(w\) at an endpoint. Its logarithmic
derivative multiplied by \(w\) is
\(-1+D[-\log(1-w)-w]/w\), which is strictly increasing.
At \(h=h_0\) the bound in (2) is at most
\((eT/h_0)\exp(\delta+D h_0/n)\), since
\(\log(1-w)\le-w\).
At \(w=1/2\) the exact bound is \(2e\,e^\delta 2^{-D}\).

## A criterion using the two halves

Split \(T_n(z)=L(z)+z^n U(1/z)\) at \(\lfloor n/2\rfloor\),
including the middle term in \(L\). Both halves have constant
coefficient one. If they are zero-free on \(|w|\le R_0>1\) and
\[
           R_0^n\min|U|>\max|L|,\qquad
           R_0^n\min|L|>\max|U|                                \tag{3}
\]
on that disk, then Rouché on \(|z|=R_0\) and on the
corresponding circle for the reciprocal puts all roots in
\(R_0^{-1}<|z|<R_0\).
At any root inside this annulus,
\[
                  \frac{zT_n'(z)}{L(z)}
                   = -n+\frac{zL'(z)}{L(z)}
                        +\frac{z^{-1}U'(z^{-1})}{U(z^{-1})}.
                                                               \tag{4}
\]
Thus the root is simple if the sum of the last two moduli
is less than \(n\).

For \(n\ge14\), the left coefficient cap gives at radius \(6/5\)
a geometric ratio at most \(9/20\). In particular
\(|L-1|\le9/11\) and \(\|\log L\|_{6/5}<2\), with the
logarithm chosen to be zero at zero. Whenever
\(\|\log L\|_{6/5}+\|\log U\|_{6/5}<n/12\) and the same sum
times \(11\) is less than \(n\), (3) and (4) hold with \(R_0=11/10\).
Indeed \(\log(11/10)>1/12\), and Cauchy's estimate on
the two nested disks bounds each logarithmic derivative
in (4) by eleven times its log norm.

## The direct ranges

Terminals through thirteen follow from the
[ordinary terminal results through thirteen](MULTIVARIATE_ALL_DEGREE_REDUCTION.md#the-fixed-terminal-obligation)
and the derivative identity for \(R_n\).

For \(14\le n\le64\) we may use \(D=7\) for the right coefficient
bounds, by (1). In the first line of (2) take \(w=1/n,1/2\).
The midpoint bound is below \(1/4\) since \(C(14)<81/500\).
For \(w=1/n\), the upper bound
\(en\exp[-7+7(C(n)+1/n)]\) has its maximum on this interval
at an endpoint. The derivative of its logarithm has the
sign of
\[
          n-7\left(\tfrac12\log(2\pi n)+\tfrac12+\tfrac1{6n}\right).
\]
The ratio of \(n\) to the quantity in parentheses increases
with \(n\). The endpoint bounds, using
\(C(14)<81/500\) and \(C(64)<47/1000\), are both below \(1/4\).
Thus \(U\) has geometric ratio at most \(3/10\) on \(6/5\).
The lower and upper norm bounds from the geometric sums give
(3) with \(R_0=6/5\), already using \(n\ge14\).
For a half with geometric ratio \(x<1/2\), the modulus of
its normalized logarithmic derivative is at most
\(x/[(1-x)(1-2x)]\).
At \(x=9/20,3/10\) the sum is below fourteen, proving (4).

For \(64\le n\le256\) again use \(D=7\) in (1) and in (2)
as bounds on the normalized coefficients of the profile.
Then \(\delta<1/3\), \(D/n<1/9\), and
\((6/5)Te^\delta<2/5\).
For \(1\le h\le7\) use the second bound in (2).
For \(h\ge8\) in the right half, the first bound in (2)
has endpoint bounds below \(3/10\), at \(h=8\) and at the
midpoint. On the disk of radius \(6/5\) the two sums give
\[
 |U-1|\le
       \sum_{h=1}^7\frac{(2/5)^h(9/8)^{h(h-1)}}{h!}
                +\frac{(9/25)^8}{1-9/25}<3/5.                  \tag{5}
\]
We used \(e^{1/9}<9/8\). Thus \(\|\log U\|_{6/5}<1\).
The log criteria following (4) apply with \(n\ge64\).

For \(n\ge256\) and \(D\ge2\log n+5\) squarefreeness follows
from the [large-curvature theorem](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md)
at terminal \(n\) and the derivative identity.
Otherwise \(\delta<1/4\) and \(D/n<1/15\).
Indeed \(C(n)(2\log n+5)\) and \((2\log n+5)/n\) decrease for
\(n\ge256\), as is seen by differentiating, and the stated
bounds hold at 256. If \(T\le2/5\), the same use of (2)
gives \((6/5)Te^\delta<5/8\) and endpoint tail bounds
below \(5/16\) for \(h\ge8\). Thus
\[
     |U-1|\le\sum_{h=1}^7\frac{(5/8)^h(15/14)^{h(h-1)}}{h!}
                     +\frac{(3/8)^8}{1-3/8}<19/20 .            \tag{6}
\]
The log norm of \(U\) is below three, and the log criteria
following (4) again apply.

## The remaining endpoint logarithm

Here \(n\ge256\) and \(T>2/5\). Put
\[
 K=\lceil13T\rceil+12,\qquad B=c_{n,n-1},\qquad \epsilon=K/n .
\]
Since \(n>(2/5)e^D\), we have
\[
 \delta <1/15,\qquad B<(107/100)T,\qquad
 \epsilon <21/500,\qquad D\epsilon<3/10 .
                                                               \tag{7}
\]
For the first inequality use
\[
 DC(n)\le \frac{5D(D+1)}{4e^D}+\frac{25D}{48e^{2D}}
\]
and evaluate at \(D=7\), since the terms decrease for \(D\ge7\).
The bound used \(\log(4\pi/5)<1\).
The other inequalities follow from \(B\le Te^\delta\) and
\(K\le(91/2)T\). The first line of (2) at \(K\) and at
the midpoint gives \(c_{n-h}\le(1/3)^h\) for \(K\le h\le n/2\).

We use the connected-sum estimate for the logarithm coefficients of
the factorial endpoint interpolation. Its inputs are sharper here
because the \(D\) Gamma shifts all point towards a constant row.
Write \(F(x)\) for the interpolated logarithm of
\(c_{n,n-x} x!/B^x\), with \(F(0)=F(1)=0\).
For \(r\ge2\), differentiating the Gamma interpolation
as in the
[endpoint proof](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md#coefficient-halves-and-endpoint-interpolation),
and grouping each Gamma series into its \(d_i\) residue classes,
gives
\[
 F^{(r)}(x)=(r-1)!\sum_{k\ge0}
     \left[\sum_i\sum_{a=1}^{d_i}(n-x+k+a/d_i)^{-r}
                              -(n-x+k+1)^{-r}\right].
\]
The last term is the binomial contribution.
The derivative is nonnegative, and for \(0\le x\le K\) and
\(2\le r\le K\) it is at most
\[
 H_r=(r-2)!\alpha^{r-2}H_2,\qquad
 \alpha=\frac1{n-K},\qquad H_2=\frac{D}{n(1-\epsilon)^2}.
\]
Indeed bound the series by its initial term and the integral
from zero to infinity. The relative correction to the integral
bound is at most \((r-1)/(n-K)\le\epsilon/(1-\epsilon)\).
Now \(H_r\le H_2<9/500\). Integrating the finite differences
over unit cubes gives \(0\le\Delta^r F(0)\le H_r\), so the
hyperedge weights satisfy
\[
 |w_r|=|e^{\Delta^r F(0)}-1|\le e^{H_2}H_r
         =(r-2)!\alpha^{r-2}\beta/e,\qquad \beta=e\,e^{H_2}H_2.
\]
This derivation uses no interpolation in the degree.

For this choice, \(\alpha/\beta\le1/(eD)\le1/(4e)\),
\(\alpha K<1/22\), \(\beta K<1\), and \(\beta B<11/500\)
by (7). The absolute part of the
[sharpened connected-sum bound](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#sharper-connected-sum)
gives log coefficients \(\ell_k\), \(1\le k\le K\), with
\[
 |\ell_k|\le Bp^{k-1}/k,\quad
 p=\frac{3\beta B}{1-\alpha K}
        \exp\left(\frac{\beta K}{2(1-\alpha K)}\right)<1/8 .
\]
The Taylor polynomial \(A(z)=\sum_{k=1}^K\ell_k z^k\)
therefore has
\(\|A\|_{9/5}+\|A\|_{6/5}<4B\) and \(\|A\|_{6/5}<2B\)
by the geometric sums.
Cauchy's estimate bounds the tail of \(e^A\) relative to
its own modulus on \(|w|\le6/5\) by
\(2(2/3)^K e^{4B}<2e^{-24/5}<1/50\),
using \(K\ge13T+12\) and (7).
The actual polynomial \(U\) agrees with \(e^A\) through
degree \(K\). Its remaining tail relative to the same modulus
is bounded by \((2/5)^K e^{2B}<e^{-54/5}<1/1000\).
It follows that \(U\) is zero-free on the disk, and
\(\|\log U\|_{6/5}\le2B+1/20\).
Together the log norms are bounded by
\[
             2+2B+1/20\le (107/50)ne^{-D}+41/20<n/100
\]
because \(n>(2/5)e^D\) and \(D\ge7\).
The two log criteria following (4) now give the conclusion.

## Exact scalar checks

The grouped estimates above use the following rational Taylor
bounds. They check the constants and the two short sums,
independently of any search over terminal orders.

```python
from fractions import Fraction as F
from math import factorial
def lower(x):
    x=F(x)
    return sum(x**k/factorial(k) for k in range(41))
def upper(x):
    x=F(x)
    return lower(x)+x**41/factorial(41)/(1-x/42)
def head(b,r):
    return sum(b**h*r**(h*(h-1))/factorial(h) for h in range(1,8))
R=F(6,5); eb=F(11,4); eb7=F(1096); eps=F(21,500)
assert upper(1)<eb and lower(7)>eb7
assert lower(F(92,25))>12*F(22,7)
assert F(92,25)/12+F(1,432)<F(31,100)
assert upper(1+7*F(31,100))/64 < F(3,8)
assert (14*eb/5)*(3*eb/14)**7 < F(3,8)
assert lower(F(9,2)) > 28*F(22,7)
assert F(9,2)/28+F(1,12*14*14) < F(81,500)
assert lower(6)>128*F(22,7)
assert F(6)/128+F(1,12*64*64)<F(47,1000)
assert 14/lower(6-7*(F(81,500)+F(1,14))) < F(1,4)
assert 64/lower(6-7*(F(47,1000)+F(1,64))) < F(1,4)
assert upper(1+7*F(81,500))/64 < F(1,4)
a=F(9,20); b=F(3,10)
assert (a/((1-a)*(1-2*a))+b/((1-b)*(1-2*b))) <14
assert R**14 > (1/(1-b))/((1-2*a)/(1-a))
assert R**14 > (1/(1-a))/((1-2*b)/(1-b))
assert R*256/eb7*upper(F(1,3)) < F(2,5)
assert eb*256/eb7/8*upper(F(11,9))<F(3,10)
assert 2*eb*upper(F(1,3))/128<F(3,10)
assert head(F(2,5),F(9,8)) + F(9,25)**8/(1-F(9,25)) < F(3,5)
assert lower(F(28,5))>256
assert lower(F(37,5))>512*F(22,7)
assert (2*F(28,5)+5)*(F(37,5)/512+F(1,12*256**2))<F(1,4)
assert (2*F(28,5)+5)/256<F(1,15)
assert F(12,25)*upper(F(1,4))<F(5,8)
assert eb*F(2,5)/8*upper(F(47,60))<F(5,16)
assert 2*eb*upper(F(1,4))/128<F(5,16)
assert head(F(5,8),F(15,14)) + F(3,8)**8/(1-F(3,8)) < F(19,20)
assert lower(1)>4*F(22,7)/5
assert F(70)/eb7+F(25*7,48)/eb7**2 < F(1,15)
assert upper(F(1,15))<F(107,100)
assert F(91,2)/eb7<eps and F(91,2)*7/eb7<F(3,10)
assert eb/13*upper(F(11,30))<F(1,3)
assert 2*eb*F(107,100)/128<F(1,3)
assert 7/(F(2,5)*eb7*(1-eps)**2)<F(9,500)
assert upper(F(9,500))<F(51,50)
assert eps/(1-eps)<F(1,22)
assert eb*F(51,50)*F(3,10)/(1-eps)**2<1
assert eb*F(51,50)*F(107,100)*7/eb7/(1-eps)**2<F(11,500)
assert 3*F(22,21)*upper(F(11,21))*F(11,500)<F(1,8)
assert F(9,5)/(1-F(9,5)/8) + F(6,5)/(1-F(6,5)/8) < 4
assert upper(F(2,5)) < F(3,2)
assert upper(F(9,10)) < F(5,2)
assert 2/lower(F(24,5)) + 1/lower(F(54,5)) < F(1,40)
assert 2*F(107,100)/eb7+F(41,20)/(F(2,5)*eb7)<F(1,100)
assert lower(1)>F(5,2) and lower(2)>F(11,2) and lower(3)>20
assert upper(F(1,15))<F(107,100)
assert upper(F(11,30))<F(3,2)
assert lower(F(1,10))>F(11,10) and upper(F(1,12))<F(11,10)
```
