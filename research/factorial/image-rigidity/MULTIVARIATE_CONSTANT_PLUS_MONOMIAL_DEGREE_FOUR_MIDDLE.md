# Grouped analytic proof for degree-four monomials at orders 14 through 256

Let \(M=\prod_i X_i^{d_i}\) have ordinary positive integer exponents
with \(\sum_i d_i=4\). Then the factorial moment \(R_n(z)=\Lambda((z+M)^n)\)
is squarefree for every \(14\le n\le256\). The derivative identity
\(R_n'=nR_{n-1}\) gives coprimality at those consecutive terminals.
The proof uses grouped bounds with rational arithmetic for a universal
Gamma comparison. No individual factorial moment polynomial is tested.

## A central-product bound for the left half

Use the normalization and split \(T_n(z)=L(z)+z^n U(1/z)\)
from the [degree at least seven proof](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_SEVEN.md).
In particular the coefficient at \(k\) is at most
\(\binom nk(k!/(n!)^{k/n})^4\), and the even midpoint belongs to \(L\).
For \(k\le n/2\),
\[
                  \left(\frac{k!}{(n!)^{k/n}}\right)^2
                                      \le\binom nk^{-1}.
\]
Indeed the geometric mean of the central factors
\(k+1,\ldots,n-k\) is at least the geometric mean of
\(1,\ldots,n\). Pair \(i\) with \(n+1-i\): their products
increase toward the center, and an odd middle factor squared
is at least each pair product. The midpoint case with
empty central product is equality.
Thus the left coefficients satisfy \(c_k\le1/\binom nk\).
For \(k\le4\) use \(1/\binom{14}{k}\), and for \(k\ge5\) use
\(2^{-k}\), since \(\binom nk\ge(n/k)^k\ge2^k\).
At \(R=28/25\), the sums give \(L\) zero-free and
\[
             \|\log L\|_{R}<3/10,\qquad
                            \|zL'/L\|_{R}<2 .
\]

## Positive coefficients after removing the first exponential

Write \(B=c_{n,n-1}\) and
\(c_{n,n-h}=B^h e^{F_n(h)}/h!\) at right-end indices.
Splitting each factorial into its positive Gamma factors
with shapes \(a/d_i\), \(1\le a\le d_i\), cancels one
factor of the binomial ratio with shape one. Thus, up to the
affine normalization \(F_n(0)=F_n(1)=0\), the log coefficient
is a sum of three terms \(\log\Gamma(n+\alpha-h)\), with
\(0<\alpha\le1\). For example the full expression for
\(F_n(x)\) is the sum over these remaining three factors of
\[
 \log\Gamma(n+\alpha-x)-\log\Gamma(n+\alpha)
                               +x\log(n+\alpha-1).
\]
For a group \(a\le n\le b\), compare it to
\[
 J_a(x)=3\bigl(\log\Gamma(a-x)-\log\Gamma(a)
                                    +x\log(a-1)\bigr).
\]
For every integer \(r\ge2\) and \(0\le x<a\),
\(0\le F_n^{(r)}(x)\le J_a^{(r)}(x)\) by the polygamma
series. Integrating over the unit cubes gives the same
inequalities on \(\Delta^r F_n(0)\) for \(r<a\).
Also \(e^{F_n(h)}\le G_a(h):=e^{J_a(h)}\) at every integer
\(0\le h<a\), by comparison of the factors
\((1-i/(n+\alpha-1))^{-1}\) and \((1-i/(a-1))^{-1}\).
The sequence \(G_a\) is rational:
\[
                       G_a(h)
              =\left((a-1)^h\frac{(a-1-h)!}{(a-1)!}\right)^3 .
\]

For any sequence \(f\) with \(f(0)=f(1)=0\), write
\(w_r=e^{\Delta^r f(0)}-1\). Newton's difference formula
expresses \(e^{f(h)}\) as the product of \(1+w_{|E|}\) over
subsets \(E\subset\{1,\ldots,h\}\) with \(|E|\ge2\).
Inclusion-exclusion then expresses \(\Delta^k e^f(0)\)
as the weight sum of hypergraphs on \(k\) labels with no
isolated vertices. In fact a graph using a given subset of
the \(k\) vertices acquires coefficient zero unless that
subset is the whole vertex set.
Consequently \(0\le\Delta^k e^{F_n}(0)\le\Delta^k G_a(0)\)
for \(2\le k<a\), by the monotonicity of the nonnegative
edge weights. This compares the coefficient after canceling
the exponential rather than bounding the absolute exponential error.

## Grouped relative error and simplicity

For a group \([a,b]\) take \(K\le\lfloor(a-1)/2\rfloor\)
and \(H=\lfloor(b-1)/2\rfloor<a\).
For \(x\ge RB\), set \(d_k=x^k\Delta^kG_a(0)/k!\) for
\(0\le k\le K\). Thus \(d_0=1\), \(d_1=0\).
The Taylor polynomial \(Q(z)\) through degree \(K\)
of \(e^{-Bz}U(z)\) has coefficients bounded in modulus on
the \(R\)-disk by these \(d_k\), from the finite difference
comparison. The polynomials \(U\) and the series \(e^{Bz}Q\)
agree through degree \(K\), so
\[
 \begin{split}
 \|e^{-Bz}U-1\|_R \ \le\ E:={}&
                 \sum_{k=2}^{K}d_k\\
      &+e^x\left[
           \sum_{h=K+1}^{H}\frac{x^h}{h!}G_a(h)
            +\sum_{k=0}^{K}d_k
                    \sum_{h>K-k}\frac{x^h}{h!}\right].          \tag{1}
 \end{split}
\]
The first sum is the centered jet. The other sums bound the tails
of \(U\) and \(e^{Bz}Q\). No value of an actual factorial moment
coefficient is required.

Take \(R_0=11/10\). The checks below choose \(x\) using
\(B\le n\exp[-4+4C(n)]\) from the Stirling bound in the
degree at least seven proof, and verify \(E<1\) and
\(A>-\log(1-E)\) in every group.
They give
\[
       x+A+3/10 < (9/100)a,\qquad
                          x+55A+2<a .
\]
The log norms of the two halves are at most \(x+A+3/10\).
Since \(\log R_0>2/21>9/100\), the two Rouché inequalities
put all roots into the annulus \(R_0^{-1}<|z|<R_0\).
On the \(R_0\)-disk the normalized logarithmic derivative
of \(U\) is bounded by \(x+55A\): use the \(Bz\) term
directly and use Cauchy's formula only on
\(\log(e^{-Bz}U)\), with \(R_0/(R-R_0)=55\).
Together with the left derivative bound this is less than
\(n\), proving simplicity by the half-polynomial criterion.

## Exact rational replay

Each block uses a rational upper bound on \(\log(2\pi a)\).
Every sum below bounds a scalar majorant throughout a block.

~~~python
from fractions import Fraction as F
from math import comb, factorial as f
def lo(x):
 x=F(x)
 return sum(x**k/f(k) for k in range(41))
def hi(x):
 x=F(x)
 return lo(x)+x**41/f(41)/(1-x/42)
def etail(x,k):
 return sum(x**h/f(h) for h in range(k+1,51))+x**51/f(51)/(1-x/52)
def roundup(x,d):
 return F(-(-x.numerator*d//x.denominator),d)
def g(a,k):
 return F((a-1)**k*f(a-k-1),f(a-1))**3
R=F(28,25)
r=R/2
left=sum(R**k/comb(14,k) for k in range(1,5))+r**5/(1-r)
left_deriv=sum(k*R**k/comb(14,k) for k in range(1,5))+r**5*(5-4*r)/(1-r)**2
assert left<1 and lo(F(3,10))>1/(1-left)
assert left_deriv/(1-left)<2
assert F(2,21)>F(9,100)
blocks=[(14,16,450),(16,20,462),(20,26,484),(26,34,510),(34,46,537),
        (46,64,567),(64,88,600),(88,116,632),(116,148,660),(148,184,684),
        (184,220,706),(220,256,724)]
for a,b,L in blocks:
 L=F(L,100)
 assert lo(L)>2*F(22,7)*a
 x=roundup(R*b/lo(4)*hi(4*(L/(2*a)+F(1,12*a*a))),1000)
 K=min((a-1)//2,40)
 J=(b-1)//2
 assert K<=J<a
 G=[g(a,k) for k in range(J+1)]
 d=[F(1),F(0)]
 d += [sum((-1)**(k-i)*comb(k,i)*G[i] for i in range(k+1))*x**k/f(k) for k in range(2,K+1)]
 assert min(d)>=0
 head=sum(d[2:])
 utail=sum(x**h/f(h)*G[h] for h in range(K+1,J+1))
 et=sum(d[k]*etail(x,K-k) for k in range(K+1))
 err=head+hi(x)*(utail+et)
 assert 0<=err<1
 A=next(F(j,100) for j in range(1,150) if lo(F(j,100))>1/(1-err))
 assert x+A+F(3,10)<F(9,100)*a
 assert x+55*A+2<a
 print(a,b,K,float(x),float(head),float(err),A)
~~~
