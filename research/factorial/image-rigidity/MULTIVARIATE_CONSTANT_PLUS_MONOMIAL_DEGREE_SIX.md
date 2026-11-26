# All-degree Strong Factorial for a constant and a monomial of total degree six

Let \(M=X_1^{d_1}\cdots X_s^{d_s}\) with \(d_i\ge1\) and
\(\sum_i d_i=6\). The ordinary factorial moment
\(R_n(z)=\Lambda((z+M)^n)\) is squarefree for every \(n\ge1\).
As \(R_n'=nR_{n-1}\), this proves SFC(2) over every
characteristic-zero field for these supports.
Combined with the
[degree at least seven theorem](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_SEVEN.md),
this covers total degrees at least six.

The proof keeps the coefficient comparisons and half-polynomial
criterion from that theorem, with refined bounds on the low
endpoint coefficients and a comparison with \(e^{Bz}\) in the
middle range. No individual factorial polynomials are tested.

## Normalization and the left half

Use \(T_n(z)=L(z)+z^n U(1/z)\) with
\(c_k=\binom nk\prod_i(d_i k)!/\prod_i((d_i n)!)^{k/n}\),
as in the degree at least seven theorem. Put
\[
 C(n)=\log(2\pi n)/(2n)+1/(12n^2),\quad T=ne^{-6},\quad
 \delta=6C(n),\quad B=c_{n,n-1},\quad R=6/5 .
\]
The coefficient comparison there gives
\(c_k\le \binom nk(k!/(n!)^{k/n})^6\), and the right endpoint
estimates apply with \(D=6\):
\[
 c_{n-h}^{1/h}\le (e/w)e^\delta(1-w)^{6(1-w)/w},\qquad
 c_{n-h}\le (Te^\delta)^h e^{6h(h-1)/n}/h!,
                 \qquad w=h/n\le1/2.                           \tag{1}
\]
At \(h=h_0\) the first bound is at most
\((eT/h_0)e^{\delta+6h_0/n}\); at \(h=n/2\) it is
\(2e e^\delta/64\). Both use the same endpoint-maximum
principle as in the cited proof.

The cases \(n\le13\) follow from the
[ordinary theorem through thirteen](MULTIVARIATE_ALL_DEGREE_REDUCTION.md)
and the derivative identity, so suppose \(n\ge14\).
For \(1\le k\le5\) the elementary estimate used in the degree
at least seven proof gives
\[
 c_k^{1/k}\le (14e/5)(3e/14)^6<3/10 .
\]
For \(6\le k\le n/2\), put \(w=k/n\). The single term
\(\binom nk w^k(1-w)^{n-k}\) in the binomial theorem is at most
one. Applying Stirling's bounds then gives
\[
   \log c_k^{1/k}
       \le 5\log w-\frac{1-w}{w}\log(1-w)+6C(k).
\]
The derivative of the two \(w\)-terms, multiplied by \(w\),
is \(6+\log(1-w)/w>0\) on \(0<w\le1/2\).
Since \(C(k)\le C(6)<31/100\), this gives
\(c_k^{1/k}\le e^{6(31/100)}/16<41/100\).
Consequently on the disk of radius \(R\),
\[
       \|L-1\|_R<3/5,\qquad \|zL'\|_R<6/5,\qquad
                            \|\log L\|_R<1 .                   \tag{2}
\]
These follow by summing the geometric ratios \(9/25\) in
the first five terms and \(123/250\) thereafter, including
their termwise derivative sums.

## Direct estimates through one hundred twenty-eight

For \(14\le n\le128\), the second inequality of (1) gives
\(R^h c_{n,n-h}\le (1/2)^h/h!\) for \(1\le h\le3\).
For \(h\ge4\) the first inequality in (1) gives the
geometric bound \((21/50)^h\) on the \(R\)-disk.
For the first bound use the values at \(n=14,128\) of
\(Rn\exp(-6+6C(n)+12/n)\).
For the right tail, use the values at the same two endpoints
of \((n/4)\exp(-5+6C(n)+24/n)\), and the midpoint bound
\(2e e^{6C(14)}/64\).
For any constant \(b\) the logarithmic derivative of
\(n\exp(6C(n)+b/n)\) has the sign of
\(n+3(1-\log(2\pi n))-1/n-b\), which is increasing for \(n\ge14\).
Thus there is no interior maximum in either endpoint estimate.
The scalar bounds give
\(\|U-1\|_R<7/10\) and \(\|zU'\|_R<11/10\).
With (2), the two Rouché inequalities for the halves hold at
radius \(R\) because \(R^{14}>17/4\) and \(R^{14}>16/3\).
The sum of logarithmic derivative moduli is at most
\((6/5)/(2/5)+(11/10)/(3/10)<14\), proving simplicity.

## Comparison with an exponential from one hundred twenty-eight to five hundred

For \(128\le n\le500\), write \(b=RB\) and
\(F(h)=\log(c_{n,n-h}h!/B^h)\).
The fixed-order Gamma differentiation in the degree at least
seven proof gives \(F(0)=F(1)=0\), and
\[
          0\le F''(x)\le\frac{6}{n(1-8/n)^2}
                      \le\frac{512}{75n}\qquad(0\le x\le8).
\]
Indeed the same derivative series and sum-integral comparison
apply with \(K=8,D=6\).
For integer \(h\le7\) summing the second differences gives
\(0\le F(h)\le \xi h(h-1)\) with \(\xi=256/(75 n)\).
The contribution of the first seven terms of
\(U-e^{Bz}\), divided by \(\min_{|z|\le R}|e^{Bz}|=e^{-b}\),
is consequently at most \(\xi b^2 e^{2b+42\xi}\).
Use \(e^v-1\le ve^v\) and sum the second factorial
moment of the exponential series.

For the remaining terms use the following uniform bounds:
\[
 \begin{array}{c|cc}
 n & b/n\le\alpha & R^h c_{n,n-h}\le v^h\quad(h\ge8)\\ \hline
 [128,256] & 7/2000 & 7/20\\
 [256,500] & 13/4000 & 59/100
 \end{array}
\]
The first column follows from \(B\le Te^\delta\) at the
left endpoint of each range for the factor \(e^\delta\).
For the other column use (1) at \(h=8\) and the midpoint.
As in the direct range it suffices to check each range's two
endpoint values for \(n\), and use its left endpoint for the
midpoint.

On each listed range the head bound is at most
\[
             (256/75)\alpha^2 n
                      \exp(2\alpha n+42(256/75)/n).
\]
This is increasing there, since
\(n+2\alpha n^2>42(256/75)\) already at the left endpoint.
Let \(n_1\) be the right endpoint and \(b_1=\alpha n_1\).
The two relative tail errors are at most
\[
      e^{b_1}v^8/(1-v),\qquad
      e^{b_1} b_1^8/(8!(1-b_1/9)).
\]
The head and tail sums are less than \(17/20\) on both
ranges by the scalar estimates below. Hence \(U\) is
zero-free on the disk and
\(\|\log U\|_R<b_1+2\). The sum of the log norms of \(L,U\)
is less than five by (2). The half-polynomial criteria at
radius \(11/10\) now give simplicity for \(n\ge128\) in this
range, using \(11\cdot5<128\) and \(5<128/12\).

## Endpoint logarithm from five hundred onward

For \(n\ge500\), put \(K=\lceil16T\rceil+6\) and
\(\epsilon=K/n\). Then
\[
        T>6/5,\quad \delta<1/20,\quad B<(263/250)T,\quad
                                      \epsilon<27/500 .
\]
Using (1) at \(h=K\) and at the midpoint gives
\(c_{n,n-h}\le(1/4)^h\) for \(K\le h\le n/2\).
Define the constants of the
[connected-sum argument](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_SEVEN.md#the-remaining-endpoint-logarithm)
with
\[
 \alpha=\frac1{n-K},\quad H_2=\frac{6}{n(1-\epsilon)^2},
                          \qquad \beta=e e^{H_2}H_2 .
\]
The Gamma derivative bound and finite-difference proof there
give exactly these inputs for the log expansion.
We have \(H_2<27/2000\), \(\alpha K<29/500\),
\(\alpha/\beta\le1/(6e)\), \(\beta K<1\) and
\(\beta B<49/1000\). Hence the connected-sum bound gives
\[
   |\ell_k|\le Bp^{k-1}/k,\qquad
    p=\frac{3\beta B}{1-\alpha K}
                     e^{\beta K/(2(1-\alpha K))}<27/100
                         \qquad (1\le k\le K).
\]
For \(A(z)=\sum_{k=1}^K\ell_k z^k\), sum the logarithmic
series at \(9/5\) and \(R\) to obtain
\(\|A\|_{9/5}+\|A\|_R<4B\) and \(\|A\|_R<(8/5)B\).
The relative tails when comparing \(U\) to \(e^A\) are at
most \(2(2/3)^K e^{4B}\) and \((3/10)^K e^{(8/5)B}\).
Use \(\log(3/2)>2/5\) and \(\log(10/3)>1\) together
with the displayed bounds on \(K,B,T\) to see that their sum
is less than \(1/50\).
Thus \(U\) is zero-free on the \(R\)-disk and
\(\|\log U\|_R<(8/5)B+1/40\).
Together with (2) this is a sum of log norms less than
\(n/100\). The half-polynomial criteria at radius \(11/10\)
prove simplicity for the remaining range.

## Rational scalar replay

The following checks use only grouped constants.

```python
from fractions import Fraction as F
from math import factorial as f
def lo(x):
 x=F(x)
 return sum(x**k/f(k) for k in range(41))
def hi(x):
 x=F(x)
 return lo(x)+x**41/f(41)/(1-x/42)
def tails(v,k):
 return v**k/(1-v),v**k*(k-(k-1)*v)/(1-v)**2
def c(n, log):
 assert lo(log)>2*F(22,7)*n
 return log/(2*n)+F(1,12*n*n)
e=F(68,25); R=F(6,5); R0=F(11,10)
assert hi(1)<e and 403<lo(6) and hi(6)<404
assert c(6,F(92,25))<F(31,100)
Cs={14:c(14,F(9,2)),128:c(128,F(67,10)),256:c(256,F(37,5)),500:c(500,F(81,10))}
assert 6*Cs[14]<F(972,1000)
assert 6*Cs[128]<F(158,1000)
assert 6*Cs[256]<F(87,1000)
assert 6*Cs[500]<F(49,1000)
assert hi(6*F(31,100))/16<F(41,100)
assert (14*e/5)*(3*e/14)**6<F(3,10)
l,ld=tails(F(123,250),6)
l+=sum(F(9,25)**k for k in range(1,6))
ld+=sum(k*F(9,25)**k for k in range(1,6))
assert l<F(3,5) and ld<F(6,5) and lo(1)>F(5,2)
print('left',float(l),float(ld))
for n in [14,128]:
 assert R*n/403*hi(6*Cs[n]+F(12,n))<F(1,2)
 assert e*n/(4*403)*hi(6*Cs[n]+F(24,n))<F(7,20)
assert 2*e/64*hi(6*Cs[14]) < F(7,20)
u,ud=tails(F(21,50),4)
u+=sum(F(1,2)**k/f(k) for k in range(1,4))
ud+=sum(k*F(1,2)**k/f(k) for k in range(1,4))
assert u<F(7,10) and ud<F(11,10)
assert R**14 > max((1+u)/(1-l),(1+l)/(1-u))
assert ld/(1-l)+ud/(1-u)<14
print('direct U',float(u),float(ud))
A=F(256,75)
for a,b,s,v in [(128,256,F(7,2000),F(7,20)),(256,500,F(13,4000),F(59,100))]:
 assert R*hi(6*Cs[a])/403<s
 for n in [a,b]:
  assert R*e*n/(8*403)*hi(6*Cs[n]+F(48,n))<v
 assert R*2*e/64*hi(6*Cs[a])<v
 assert a+2*s*a*a>42*A
 head=A*s*s*b*hi(2*s*b+42*A/b)
 utail=hi(s*b)*v**8/(1-v)
 etail=hi(s*b)*(s*b)**8/f(8)/(1-s*b/9)
 print('poisson',float(head),float(utail),float(etail))
 assert head+utail+etail<F(17,20)
assert lo(2)>F(20,3)
eps=F(27,500); eb=F(263,250)
assert 16/F(403)+F(7,500)<eps
assert hi(F(1,20))<eb
h2=6/(500*(1-eps)**2)
assert h2<F(27,2000) and hi(F(27,2000))<F(507,500)
assert eps/(1-eps)<F(29,500)
assert e*F(507,500)*6*eps/(1-eps)**2<1
assert e*F(507,500)*eb*6/(403*(1-eps)**2)<F(49,1000)
assert 3*F(49,1000)/(1-F(29,500))*hi(F(1,2)/(1-F(29,500)))<F(27,100)
P=F(27,100)
assert lo(4*P)*(1-F(9,5)*P)*(1-R*P)>1
assert lo(F(8,5)*P)*(1-R*P)>1
assert e/16*hi(F(1,20)+6*eps)<F(1,4)
assert 2*e/64*hi(F(1,20))<F(1,4)
assert hi(F(2,5))<F(3,2)
assert hi(1)<F(10,3)
T=F(6,5)
assert F(500,404)>T
assert 2/lo(F(12,5)+(F(32,5)-4*eb)*T) + 1/lo(6+(16-F(8,5)*eb)*T)<F(1,50)
assert lo(F(1,40))>F(50,49)
assert F(41,40)/500+F(8,5)*eb/403<F(1,100)
print('all degree-six grouped constants pass')
```
