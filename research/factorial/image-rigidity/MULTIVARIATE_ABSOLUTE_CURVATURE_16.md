# All-degree angular separation at absolute curvature sixteen

Let \(a,b\) be distinct ordinary nonnegative exponent columns.
Omit equal rows and put
\[
 S=\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)}+\tau ,
\]
where \(\tau\ge0\) is an optional common Gaussian strength,
with coefficient multiplier \(e^{-\tau k(m-k)/(2m)}\) at
index \(k\) of order \(m\).

**Theorem.** For \(S\ge16\), the consecutive moment polynomials
have disjoint root arguments in every degree. This includes
arbitrary private rows and orientations and independent positive
normalizations of the two orders. At \(\tau=0\) it proves
SFC(2) over every characteristic-zero field.

The proof strengthens the
[entropy and connected logarithm comparison at curvature twenty](MULTIVARIATE_ABSOLUTE_CURVATURE_20.md).
The endpoint cost controls both the sum of the row curvatures and
the largest relative gap. Keeping these controls separate permits
sharper Gamma derivative bounds. No individual moment polynomial
is tested.

## The largest relative gap costs endpoint mass

Keep \(V,W=1+V\), and the oriented envelopes \(Q,q\) of the
[previous proof](MULTIVARIATE_ABSOLUTE_CURVATURE_20.md#the-endpoint-cost-absorbs-the-excess-curvature).
Both bounds in its (1) hold at \(S\ge16\).
In its derivative comparison take \(C=1+S\ge17\).
The difference polynomial at \(C=17\) is
\(r^3-4r^2+7r+11>0\), by its strictly positive derivative.
All four envelopes in its (2) also hold with \(L=3V/4\).

Let \(L_0\) be the maximum of one and of \(d/\ell\) on shared
rows with gap \(d>0\) and smaller exponent \(\ell\).
In an orientation with all starting exponents positive,
\(\max(1,\max_i |b_i-a_i|/a_i)\le L_0\). We claim
\[
                       L_0q\le 3e^{-S/2}.                     \tag{1}
\]
If \(L_0=1\) this follows from \(q\le e^{-S/2}\).
Otherwise take a row attaining \(r=L_0\). Its exact logarithmic
cost is at least \(S_i/2+\ell E(r)\) in both orientations,
with \(E(r)\ge0\) as in the previous proof. Since \(\ell\ge1\),
\[
 r e^{-\ell E(r)}\le r e^{-E(r)}
                   \le r(1+r)e^{-r/2}<3 .
\]
The second inequality uses
\(E(r)=r/2+1/2-1/[2(1+r)]-\log(1+r)\).
The last function has its unique maximum at the positive root of
\(r^2-3r-2=0\), lying between \(7/2\) and \(18/5\).
Its value is bounded by \((414/25)e^{-7/4}<3\).
The other rows contribute at least their original costs \(S_i/2\),
giving (1).

## The short ranges

Write \(N=n+1\ge3\) and \(m=n,N\).
For \(N\le256\) use the discrete bound
\(c_{m,k}\le\binom mk e^{-8 k(1-k/m)}\) in either
orientation. At \(k=1\), the bound \(m e^{-8(1-1/m)}\)
takes its maximum at an endpoint of \([2,256]\),
and is below \(1/8\). For \(2\le k\le m/2\), the
bound on \(c_{m,k}^{1/k}\) has its maximum at \(w=2/m\)
or \(w=1/2\), using \(e e^{-8(1-w)}/w\).
The first expression at \(w=2/m\) takes its maximum for
\(4\le m\le256\) at an endpoint, and both this maximum
and the value at \(1/2\) are below \(1/8\).
For \(m\ge4\), \(e^{1/(2m)}/8<1/7\). For \(m=2,3\)
the only possible index in a half besides the constant is
\(k=1\), and the same geometric ratio bound follows from
\(m e^{-8(1-1/m)}<1/16\).
For \(N\ge256\) and \(S\ge2\log N+5\), the
[large-curvature coefficient estimate](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)
again gives ratio below \(1/7\) on \(|z|\le e^{1/(2m)}\).
Both cases are covered by the short phase argument of the
curvature-twenty theorem.

Otherwise \(N>256\), \(S<2\log N+5\), and \(D=Nq\).
Then \(L/n\le89/100\) by the decrease of
\(3S(S+2)/[4(e^{(S-5)/2}-1)]\) on \(S\ge16\).
Also \(L/n\le(43/10000)/D\),
since \(LD/n\le(3/4)WD/n\) and
\(WD/n\le (257/256)17e^{-8}\).
Suppose \(D\le4/25\).
For any positive integer \(k\) let
\[
       T_k=\max\left\{
           \frac{43}{8900}e^{89k/100},
           \frac4{25}e^{(43/10000)k/(4/25)}\right\}.
\]
This bounds \(D\exp(kL/n)\): for \(D\) below the switch
\(43/8900\) the upper envelope is linear in \(D\);
above the switch the logarithm of that envelope has derivative
\((D-k(43/10000))/D^2\), so its maximum is at one end.
The first three coefficients satisfy \(c_{m,k}\le T_k^k/k!\).
The endpoint comparison of the previous proof bounds
\(c_{m,k}^{1/k}\) for \(k\ge4\) by
\(\max(eT_4/4,\,2e^{-3})<1/8\).
On radius \(R_m=e^{1/(2m)}<501/500\),
the replay gives \(|H_m-1|<183/1000\) and
\(|zH_m'|<201/1000\) from these first three terms
and the geometric tail. Thus \(|zH_m'/H_m|<1/4\).
The ratio \((1183/1000)/(817/1000)<\sqrt e\) suffices
for Rouché on this radius and its reciprocal. The same
angular phase contradiction applies.

## The long range with sharper Gamma derivative estimates

Now \(D>4/25\) and \(K=\lceil32D\rceil+14\).
Then \(N>(4/25)e^8>476\), hence \(n\ge476\).
The endpoint bounds give
\[
 \frac{D}{n}<\frac{337}{10^6},\quad
 \frac{WD}{n}<\frac{572}{10^5},\quad
 \frac{L_0D}{n}<\frac{101}{10^5},\quad K\le126D.
\]
In particular \(WK/n<721/1000,\ L_0K/n<16/125\),
and \(K+1<n/2\).
Testing \(k=K\) using
\((e/32)e^{(3/4)WK/n}<1/6\), and at the midpoint,
gives \(c_{m,k}\le6^{-k}\) for \(k\ge K\) in the half.
The same first coefficient estimate used at curvature twenty
has numerator \(5V+S(3+\log N)\le50R^2\),
\(R=\log N>6\).
Since \(50R^2/(e^R-1)<9/2\), we have
\(|(\log B_m)'|<11/(2n)\) and
\[
                B_m\le De^{L/n+11/(2n)}<(26/25)D.
\]

For a positive starting endpoint use the same Gamma interpolation.
With \(\kappa=16/125\), the Gamma arguments in its derivatives
are at least \((1-\kappa)ma_i\) and \((1-\kappa)m\) for
the binomial term. Put \(\alpha=L_0/[(1-\kappa)n]\) and
\[
 H_2=\frac{W}{(1-\kappa)n}
                 \left(1+\frac{K}{(1-\kappa)n}\right),
 \qquad H_r=(r-2)!\alpha^{r-2}H_2.
\]
The polygamma estimate
\[
 |\psi_{r-1}(y)|\le(r-2)!y^{1-r}(1+(r-1)/y)
\]
applied to each positive argument proves the bound by \(H_r\).
The order derivative uses \(\psi_r(y)\) and replaces
\(r-1\) by at most \(K\) in the parenthesis. Its additional
factor relative to \(H_r\) is at most
\((r-1)/[(1-\kappa)n]\le2r/n\).
The Gaussian term is absorbed by its contribution to \(W\)
in \(H_2\), for the second difference and its order derivative.
Thus the same bounds apply to the finite differences and
their derivatives as at curvature twenty. Here
\(H_2<(121/100)W/n<11/250\) and \(\alpha K<147/1000\).

Apply the
[differentiated connected sum](MULTIVARIATE_ABSOLUTE_CURVATURE_20.md#a-smaller-differentiated-connected-logarithm-bound)
with \(\beta=e^{H_2+1/2}H_2\) and \(u=\alpha K\).
Since \(t=\alpha/\beta\le e^{-1/2}<2/3\), the tree base
is at most \((5/3)^{5/2}<18/5\). Use this bound in
place of \(256/81\) in \(\Xi\). The scalar bounds give
\[
 \beta K<151/100,\quad \beta B_m<1/80,\quad
 p_m=\Xi B_m<3/25,\quad
                   |\ell_k'(m)|\le7B_m p_m^{k-1}/n.
\]
The derivative constant uses \(11/2+4/e<7\).
For the truncated log \(A_m\) the absolute coefficient
norms on \(9/5\) and \(6/5\) are at most \((21/10)B_m\)
and \((7/5)B_m\), respectively. Its difference across
orders on \(6/5\) is below \(11D/n\).
The exponential tail and the actual polynomial tail on
\(6/5\), relative to \(e^{A_m}\), are bounded by
\((2/3)^{K+1} e^{(7/2)B_m}\) and
\((1/5)^K e^{(7/5)B_m}\). They are less than \(e^{-7}\)
and \(e^{-K}\), respectively, and \(K\ge20\).
Their sum is below \(1/1000\), giving
\(\|\log H_m-A_m\|_{6/5}<1/900\).

For a private starting row, the same factorial cancellation
gives a geometric ratio \((6/5)e^{3/2+LK/n}q<1/300\)
through \(K\), and \(1/5\) beyond \(K\), on radius \(6/5\).
The two tails bound \(\|\log H_m\|_{6/5}<1/240\).
In both cases the individual and matched log bounds
are \(M=2D\) and \(\delta=1/100\), respectively.
Use Rouché at \(e^{5D/m}\) and its reciprocal to obtain
\(|u_m|\le5D/m\). The half arguments lie in \(21/20\).
The same real radius comparison gives
\[
       n|u_N-u_n|\le \frac{5D/n+2\delta}{1-14M/n}<1/40 .
\]
The angular correction difference at a proposed common
argument is at most
\((14\delta+210(M/n)/40)\theta<\theta/3\).
As before, the phase difference rules out the
common argument. Terminal two follows from strict
positivity of the quadratic on the real line.

## Exact rational scalar checks

~~~python
from fractions import Fraction as F
from math import factorial

def lo(x):
    return sum(F(x)**k/factorial(k) for k in range(61))

def hi(x):
    x = F(x)
    return lo(x)+x**61/(factorial(61)*(1-x/62))

assert hi(F(1,8))<F(8,7) and hi(F(1,4))<F(4,3)
assert 2/lo(4)<F(1,16) and 256/lo(8-F(1,32))<F(1,8)
assert 2/lo(3)<F(1,8) and 128/lo(7-F(1,16))<F(1,8)
assert lo(F(1,2))>F(6,4)
assert F(3,4)*16*18/(lo(F(11,2))-1)<F(89,100)
assert F(3,4)*17*F(257,256)/lo(8)<F(43,10000)
assert hi(F(1,512))<F(501,500)
small=F(43,10000)/F(89,100)
d=F(4,25)
def short(k):
    return max(small*hi(k*F(89,100)),d*hi(k*F(43,10000)/d))
assert hi(1)*short(4)/4<F(1,8)
r=F(501,4000)
head=[(F(501,500)*short(k))**k/factorial(k) for k in range(1,4)]
assert sum(head)+r**4/(1-r)<F(183,1000)
assert sum(k*a for k,a in enumerate(head,1))+r**4*(4-3*r)/(1-r)**2<F(201,1000)
assert F(201,817)<F(1,4) and lo(F(1,2))>F(1183,817)
assert F(7,2)**2-3*F(7,2)-2<0 and F(18,5)**2-3*F(18,5)-2>0
assert F(414,25)/lo(F(7,4))<3
assert d*lo(8)>476 and hi(6)<477
base=F(477,476)/lo(8)
assert base<F(337,10**6) and 17*base<F(572,10**5) and 3*base<F(101,10**5)
v=F(721,1000)
assert 126*F(572,10**5)<v and 126*F(101,10**5)<F(16,125)
kap=F(16,125)
u=F(147,1000)
assert kap/(1-kap)<u and 126*F(337,10**6)+F(1,476)<F(1,2)
factor=(1+126*F(337,10**6)/(1-kap))/(1-kap)
assert factor<F(121,100)
wn=F(572,10**5)/d
assert F(121,100)*wn<F(11,250)
assert hi(1+F(3,4)*v)/32<F(1,6)
assert 22+131/F(6)+190/F(36)<50
assert 50*36/(lo(6)-1)<F(9,2)
assert hi(F(3,4)*wn+F(11,2*476))<F(26,25)
assert F(11,2)+4/lo(1)<7
assert F(5,3)**5<F(18,5)**2
inflate=F(121,100)*hi(F(1,2)+F(11,250))
assert inflate*v<F(151,100)
assert inflate*F(26,25)*F(572,10**5)<F(1,80)
assert F(18,5)*F(1,80)/(1-u)*hi(F(151,100)*(F(1,2)+u/(6*(1-u))))<F(3,25)
assert lo(F(21,10)*F(3,25))>F(125,98)
assert lo(F(7,5)*F(3,25))>F(125,107)
assert 7*F(6,5)*F(26,25)/F(107,125)<11
assert hi(F(2,5))<F(3,2) and F(2,5)*(32*d+15)-F(7,2)*F(26,25)*d>7
assert hi(F(11,10))<5 and F(7,5)*F(26,25)/32<F(1,10)
assert 1/lo(7)+1/lo(20)<F(1,1000) and lo(F(1,900))>F(1000,999)
assert F(6,5)*hi(F(3,2)+F(3,4)*v)/lo(8)<F(1,300)
assert F(1,299)+F(1,5)**20/(1-F(1,5))<F(1,250)
assert lo(F(1,240))>F(250,249)
assert F(7,5)*F(26,25)+F(1,900)/d<2
assert 11*F(337,10**6)+F(2,900)<F(1,100) and F(2,240)<F(1,100)
M=2*F(337,10**6)
R=5*F(337,10**6)
assert hi(R)<F(21,20)
assert (R+F(1,50))/(1-14*M)<F(1,40)
assert F(14,100)+210*M/40<F(1,3)
assert 3/lo(F(16,3))<F(1,16)
assert F(2,5)*(32*d+15)-F(7,2)*F(26,25)*d>7
assert lo(1)>F(8,3)
~~~
