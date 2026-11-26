# All-degree Strong Factorial for a constant and a monomial of total degree five

Let \(M=X_1^{d_1}\cdots X_s^{d_s}\) with positive integer
\(d_i\) and \(\sum_i d_i=5\). Then
\(R_n(z)=\Lambda((z+M)^n)\) is squarefree for every \(n\ge1\),
where \(\Lambda\) is the ordinary factorial functional.
Since \(R_n'=nR_{n-1}\), this proves SFC(2) over every
characteristic-zero field. Together with the theorems for
[degree six](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_SIX.md)
and [degree at least seven](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_SEVEN.md),
it covers every total degree at least five.

The proof uses the first-order logarithm estimate for fixed
degrees. The extra weight needed to differentiate that estimate
in the order parameter can be omitted here. This improves the
long-range bound at degree five. The remaining ranges use
grouped scalar inequalities.

## The left half

Use the endpoint normalization \(T_n(z)=L(z)+z^n U(1/z)\) and
the exact coefficient comparison
\(c_k\le\binom nk (k!/(n!)^{k/n})^5\)
from the degree at least seven theorem. Set
\[
 C(n)=\frac{\log(2\pi n)}{2n}+\frac1{12n^2},\quad
 T=ne^{-5},\quad \delta=5C(n),\quad B=c_{n,n-1},\quad R=6/5 .
\]
For \(n\le13\), the derivative identity and the
[ordinary theorem through thirteen](MULTIVARIATE_ALL_DEGREE_REDUCTION.md)
supply the claim. Assume \(n\ge14\).
For \(1\le k\le5\),
\[
             c_k^{1/k}\le e^5(k!)^{4/k}/n^4<19/100 ,
\]
using \(\binom nk\le n^k/k!\), Stirling's lower bound for
\(n!\), and \(120^{4/5}<47\).

For \(6\le k\le n/2\), put \(w=k/n\).
Robbins's upper and lower Stirling bounds give
\[
 \binom nk\le
 \frac{\exp(n[-w\log w-(1-w)\log(1-w)]+1/(12n))}
      {\sqrt{2\pi k(1-w)}}.
\]
The \(w\)-dependent expression \(4\log w-(1-w)\log(1-w)/w\)
has derivative of the sign of \(5+\log(1-w)/w>0\)
on \(0<w\le1/2\). Thus the entropy term combined with
\(5\log w\) is at most \(-3\log2\).
With \(n\ge2k\), the remaining Stirling factors give
\[
 c_k^{1/k}\le \frac18
       \exp\left(\frac{4\log(2\pi k)+\log2}{2k}
                           +\frac{11}{24k^2}\right)<23/50.
\]
The exponent decreases at \(k\ge6\), and the last bound uses
\(k=6\). Summing the geometric majorants at radius \(R\) gives
\[
          \|L-1\|_R<2/5,\qquad \|zL'\|_R<1,\qquad
                                      \|\log L\|_R<1.           \tag{1}
\]

## Direct estimates and the comparison with an exponential

The two right coefficient bounds from the degree at least seven
proof apply with five in place of seven:
\[
 c_{n-h}^{1/h}\le(e/w)e^\delta(1-w)^{5(1-w)/w},\qquad
 c_{n-h}\le (Te^\delta)^h e^{5h(h-1)/n}/h!,
                \qquad w=h/n\le1/2 .                          \tag{2}
\]
The first bound has its maximum at an endpoint of any closed
\(w\)-interval. At the left endpoint \(h=h_0\),
\(\log(1-w)\le-w-w^2/2\) bounds it by
\[
         \frac{eT}{h_0}\exp(\delta+5w(1+w)/2).
                                                               \tag{3}
\]
At the midpoint use \(2e e^\delta/32\).
To make these estimates uniform on intervals in \(n\), the
function \(n\exp(5C(n)+c/n+b/n^2)\), for constants \(c,b\ge0\),
has no interior maximum on \(n\ge14\).
Indeed the derivative of its logarithm has the sign of
\(n+(5/2)(1-\log(2\pi n))-c-(5/6+2b)/n\), which increases there.

For \(14\le n\le36\), the second inequality of (2) gives
\(R^h c_{n,n-h}\le(57/100)^h/h!\) for \(1\le h\le3\).
Equations (2)–(3) at \(h_0=4\) and the midpoint
give the bound \((47/100)^h\) for \(h\ge4\) on the \(R\)-disk.
The endpoint checks at \(n=14,36\) are recorded below.
The sums give \(\|U-1\|_R<43/50\) and
\(\|zU'\|_R<29/20\).
With (1) the Rouché bounds hold at radius \(R\), since
\(R^{14}\) exceeds both ratios of maximum to minimum half
moduli. The sum of the logarithmic derivative moduli is below
\(14\), so the derivative criterion of the degree at least
seven proof proves simplicity.

For \(36\le n\le128\) write \(b=RB\) and, as in the degree
six proof, \(F(h)=\log(c_{n,n-h}h!/B^h)\).
The fixed-order Gamma derivative bound there gives
\(F(0)=F(1)=0\) and
\(0\le F''(x)\le5/[n(1-7/n)^2]\) on \(0\le x\le7\).
In particular \(0\le F(h)\le\xi h(h-1)\) for integer
\(h\le6\) if \(\xi\ge5/[2n(1-7/n)^2]\).
The following rational bounds apply:
\[
 \begin{array}{c|ccc|c}
 n & b\le b_1&\xi&R^h c_{n,n-h}\le v^h\;(h\ge7)
       &\text{relative error}\\ \hline
 [36,64]&19/25&27/250&v=7/20&9/20\\
 [64,96]&99/100&1/20&v=11/25&3/4\\
 [96,128]&31/25&31/1000&v=27/50&19/25
 \end{array}
\]
For \(b\) use \(B\le Te^\delta\). The geometric tails use
(3) at \(h_0=7\) and the midpoint of the half.
The last column bounds \(|U-e^{Bz}|/|e^{Bz}|\) on
\(|z|\le R\), obtained by summing
\[
  e^{b_1}\sum_{h=2}^{6}\frac{b_1^h}{h!}
                           (e^{\xi h(h-1)}-1)
          +e^{b_1}\frac{v^7}{1-v}
          +e^{b_1}\frac{b_1^7}{7!(1-b_1/8)}.
\]
The first term uses the exact zeros \(F(0)=F(1)=0\), and
the other terms are the two tails.
Thus \(U\) is zero-free on the disk, and
\(\|\log U\|_R\) is less than \(19/25+3/5\)
on the first range and \(31/25+3/2\) on the other ranges.
With (1) the total log norm is below \(n/12\);
the half-polynomial test on radius \(11/10\) applies.

## The connected logarithm for \(n\ge128\)

Put \(K=\lceil14T\rceil+6\) and \(\epsilon=K/n\).
Then \(\delta<33/250\), \(B<(571/500)T\),
and \(\epsilon<3/20\).
Using (2) at \(h=K\) and the midpoint gives
\(R^h c_{n,n-h}\le(29/50)^h\) for \(h\ge K\).
Use the Gamma derivative estimate from the
[degree at least seven proof](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_SEVEN.md#the-remaining-endpoint-logarithm)
with
\[
         \alpha=1/(n-K),\quad H_2=5/[n(1-\epsilon)^2],
                              \qquad \beta=e^{H_2}H_2 .
\]
It gives \(|w_j|\le(j-2)!\alpha^{j-2}\beta\) for
the finite-difference hyperedge weights, \(2\le j\le K\).
Here the connected sum is only bounded in absolute value.
Thus the extra factor \(e\) used to control derivatives of
the weights in the
[sharpened connected-sum proof](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#sharper-connected-sum)
is unnecessary. The argument for its absolute estimate
still gives the same tree factor \(3\).
In fact with \(t=\alpha/\beta\le1/5\), its exact base
\((1+t)^{(1+t)/t}\) is increasing, by logarithmic
differentiation, and at \(t=1/5\) equals
\((6/5)^6<3\). Its remaining prefactor is still at most \(e/k\)
for \(k\) labels and hence at most one for \(k\ge3\).
The cases \(k=1,2\) have the same bound directly.

In this range \(H_2<11/200\), \(\alpha K<3/17\),
\(\beta K<11/10\), and \(\beta B<113/2000\).
Therefore the log coefficients satisfy
\[
 |\ell_k|\le Bp^{k-1}/k,\qquad
         p=\frac{3\beta B}{1-\alpha K}
                       e^{\beta K/(2(1-\alpha K))}<81/200 .
\]
For \(A(z)=\sum_{k=1}^K\ell_k z^k\), summing the logarithmic
series gives
\(\|A\|_{9/5}+\|A\|_R<(49/10)B\) and
\(\|A\|_R<(33/20)B\).
The relative errors comparing \(U\) to \(e^A\) on the
\(R\)-disk are bounded by
\[
      2(2/3)^K e^{(49/10)B},
            \qquad (29/50)^K e^{(33/20)B}.
\]
Using \(\log(3/2)>2/5\), \(\log(50/29)>27/50\),
\(K\ge14T+6\), and \(T\ge128/149\), their sum is below
\(19/100\). The first relative bound is at most
\(2e^{-12/5}\) because \((49/10)(571/500)<28/5\).
Consequently \(U\) is zero-free there and
\(\|\log U\|_R<(33/20)B+1/4\). The sum of the log norms
with (1) is below \(n/40\). The half-polynomial test at
radius \(11/10\) proves simplicity in this last range.

## Rational scalar replay

The following checks use only the grouped inequalities above.

```python
from fractions import Fraction as Q
from math import factorial as f
def lo(x):
 x=Q(x)
 return sum(x**k/f(k) for k in range(41))
def hi(x):
 x=Q(x)
 return lo(x)+x**41/f(41)/(1-x/42)
def tail(v,k):
 return v**k/(1-v),v**k*(k-(k-1)*v)/(1-v)**2
def cb(n,log):
 assert lo(log)>2*Q(22,7)*n
 return log/(2*n)+Q(1,12*n*n)
R=Q(6,5); e=Q(68,25)
assert hi(1)<e and 148<lo(5) and hi(5)<149 and lo(Q(7,10))>2
assert lo(Q(92,25))>12*Q(22,7)
assert 120**4<47**5 and Q(149*47,14**4)<Q(19,100)
assert hi((4*Q(92,25)+Q(7,10))/12+Q(11,24*36))/8<Q(23,50)
Cs={14:cb(14,Q(9,2)),36:cb(36,Q(109,20)),64:cb(64,Q(6)),
    96:cb(96,Q(13,2)),128:cb(128,Q(67,10))}
ds={14:Q(81,100),36:Q(19,50),64:Q(47,200),96:Q(17,100),128:Q(33,250)}
for n in ds: assert 5*Cs[n]<ds[n]
l,ld=tail(R*Q(23,50),6)
l+=sum((R*Q(19,100))**k for k in range(1,6))
ld+=sum(k*(R*Q(19,100))**k for k in range(1,6))
assert l< Q(2,5) and ld<1 and lo(1)>Q(5,3)
print('left',float(l),float(ld))
def tailtest(n,h,ratio):
 w=Q(h,n)
 return R*e*n/(148*h)*hi(ds[n]+Q(5,2)*w*(1+w))<ratio
for n in [14,36]:
 assert R*n/148*hi(ds[n]+Q(10,n))<Q(57,100)
 assert tailtest(n,4,Q(47,100))
assert R*2*e/32*hi(ds[14]) < Q(47,100)
u,ud=tail(Q(47,100),4)
u+=sum(Q(57,100)**k/f(k) for k in range(1,4))
ud+=sum(k*Q(57,100)**k/f(k) for k in range(1,4))
assert u<Q(43,50) and ud<Q(29,20)
assert R**14>max((1+u)/(1-l),(1+l)/(1-u))
assert ld/(1-l)+ud/(1-u)<14
print('direct U',float(u),float(ud))
for a,b,B,xi,v,limit in [(36,64,Q(19,25),Q(27,250),Q(7,20),Q(9,20)),
    (64,96,Q(99,100),Q(1,20),Q(11,25),Q(3,4)),
    (96,128,Q(31,25),Q(31,1000),Q(27,50),Q(19,25))]:
 assert R*b/148*hi(ds[a])<B
 assert 5/(2*a*(1-Q(7,a))**2)<xi
 assert tailtest(a,7,v) and tailtest(b,7,v)
 assert R*2*e/32*hi(ds[a])<v
 head=hi(B)*sum(B**h/f(h)*(hi(xi*h*(h-1))-1) for h in range(2,7))
 utail=hi(B)*v**7/(1-v)
 etail=hi(B)*B**7/f(7)/(1-B/8)
 print('poisson',float(head),float(utail),float(etail),float(head+utail+etail))
 assert head+utail+etail<limit
assert lo(Q(3,5))>Q(20,11) and lo(Q(3,2))>Q(25,6)
assert (1+Q(19,25)+Q(3,5))*12<36 and (1+Q(31,25)+Q(3,2))*12<64
assert Q(14,148)+Q(7,128)<Q(3,20)
assert hi(ds[128])<Q(571,500)
eps=Q(3,20); Hb=Q(11,200); eb=Q(571,500)
assert 5/(128*(1-eps)**2)<Hb and hi(Hb)<Q(1057,1000)
assert Q(1057,1000)*5*eps/(1-eps)**2<Q(11,10)
assert Q(1057,1000)*5*eb/(148*(1-eps)**2)<Q(113,2000)
assert Q(6,5)**6<3
assert 3*Q(113,2000)*Q(17,14)*hi(Q(11,10)*Q(17,28))<Q(81,200)
p=Q(81,200)
assert lo(Q(49,10)*p)*(1-Q(9,5)*p)*(1-R*p)>1
assert lo(Q(33,20)*p)*(1-R*p)>1
assert e/14*hi(ds[128]+5*eps)<Q(29,60)
assert 2*e/32*hi(ds[128])<Q(29,60)
assert hi(Q(27,50))<Q(50,29) and hi(Q(2,5))<Q(3,2)
assert Q(49,10)*eb<Q(28,5)
assert 2/lo(Q(12,5)) + Q(29,21)/lo(Q(81,25)+(Q(189,25)-Q(33,20)*eb)*Q(128,149)) < Q(19,100)
assert lo(Q(1,4))>Q(100,81)
assert Q(5,4)/128+Q(33,20)*eb/148<Q(1,40)
print('all degree-five grouped constants pass')
```
