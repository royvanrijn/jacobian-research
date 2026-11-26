# All-degree angular separation from two opposed shared gaps

An ordinary two-monomial profile has disjoint consecutive root
arguments in every degree if it contains two shared unequal
rows with opposite orientations, each with gap at least five
times its smaller exponent. The magnitudes and smaller exponents
of the two rows may differ. All remaining rows are arbitrary,
including private rows. The conclusion also holds with a
common physical Gaussian addition of strength \(\tau\ge0\).
At \(\tau=0\) this proves SFC(2) for the indicated ordinary
profiles over every characteristic-zero field.

The two retained rows control both oriented endpoints. Using
binomial entropy at the midpoint keeps the far ends of the
coefficient halves small even for the first case of relative
gap five. No terminal polynomials are tested.

## Two oriented endpoint estimates

Use the notation \(S,V,W=1+V\) of the
[oriented entropy criterion](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md),
including \(\tau\) in \(S,V\). The
[absolute curvature sixteen theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md)
permits us to assume \(S<16\).
In either orientation let \(d\) be the gap of the retained
increasing row and \(e_0\) the gap of the retained decreasing
row. Both are at least five. Let \(s_0\) be the curvature
of all other rows plus \(\tau\). Since the retained curvature
is at least \(5(d+e_0)/6\), we have \(d+e_0\le19\) and
\[
 V\le(d+e_0-5)^2+25+s_0(s_0+2)<222.
\]
For the last bound put \(s=d+e_0\). The resulting quadratic
upper bound after \(s_0\le16-5s/6\) is convex and is below
\(222\) at \(s=10,19\).

Retain exactly these two shared rows in the two functions
\(Q_\pm(w)\) of the
[oriented entropy criterion](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md).
Write \(q=\max(Q_+(0),Q_-(0))\).
For one retained row with ratio \(r=d/\ell\ge5\), its
increasing endpoint cost per unit gap is
\((1+1/r)\log(1+r)-1\ge (6/5)\log6-1\).
Its decreasing endpoint cost per unit gap is
\(1-\log(1+r)/r\ge1-\log6/5\).
Both lower bounds exceed \(2/5\), and their sum is \(\log6\).
It follows that in both orientations
\[
 q\le6^{-5}
       \exp\big[-(2/5)(d+e_0-10)-s_0/2\big].                 \tag{1}
\]

On an increasing row the derivative of \(A(w)/w\)
in the oriented criterion is at most \(d^2/(2\ell)\).
Its concave correction \(G(w)/2\) is at most \(dw/(2\ell)\).
For the decreasing row, on \(w\le1/2\) the arguments of its
shared affine function are at least \(e_0/2\).
The derivative of \(A(w)/w\) is then at most \(e_0\).
Its concave correction is at most
\(w\log(1+e_0/\ell)/2\le w e_0/2\).
Thus one may use the common increment bound
\[
 \log Q_\pm(w)\le \log q+Lw,\qquad
 L=\min\!\left\{\frac{3V}{4},
        \frac{\max(d,e_0)^2+\max(d,e_0)+3\min(d,e_0)+s_0}{2}\right\}.
                                                                    \tag{2}
\]
The value \(3V/4\) is the original oriented increment bound.

Let \(L_0\) be the maximum of one and all shared relative gaps
\(d_i/\min(a_i,b_i)\). Equations (1) and (2) give the useful
four budgets
\[
        q\le\frac1{7776},\quad Wq\le\frac{51}{7776},\quad
        L_0q\le\frac5{7776},\quad Lq\le\frac{25}{7776}.        \tag{3}
\]
Indeed \(W\le1+d^2+e_0^2+s_0(s_0+2)\) and
\(L_0\le\max(d,e_0,1+s_0)\).
For the \(Wq\) and \(Lq\) bounds, at fixed \(s_0\) each
of the two logarithmic derivatives of the polynomial
prefactors in \(d,e_0\) is at most \(2/5\).
This is immediate for \(W\). For the second bound in the
minimum in (2), in either ordering of the two gaps, use
\((2d+1)/(d^2+d+3e_0)\le2/d\) and
\(3/(d^2+d+3e_0)\le1/15\).
The maxima are therefore at \(d=e_0=5\).
Then the remaining expressions
\((51+s_0(s_0+2))e^{-s_0/2}\) and
\((45/2+s_0/2)e^{-s_0/2}\) decrease with \(s_0\).
For \(L_0q\) apply the same argument to \(d\) and \(e_0\)
separately, and use \((1+s_0)e^{-s_0/2}\le2\).
These comparisons hold for the larger of the two endpoint
values, regardless of which orientation attains it.

## A binomial entropy endpoint bound

Let
\[
 H(w)=-w\log w-(1-w)\log(1-w),\qquad B(w)=e^{H(w)/w}.
\]
The binomial identity gives \(\binom mk^{1/k}\le B(k/m)\).
For the two chosen rows the maximum of \(B(w)Q_\pm(w)\)
on a closed interval of \(0<w\le1/2\) is attained at an
endpoint. In fact the derivative of \(w(\log B)'(w)\) is
\[
       -\frac{w/(1-w)+\log(1-w)}{w^2}
                                  \ge-4(1-\log 2)>-5/4.
\]
The first inequality follows by expanding the numerator
as \(\sum_{j\ge2}(1-1/j)w^j\).
The derivative of \(w(\log Q_\pm)'(w)\) from any retained
row is nonnegative by the oriented entropy criterion.
The increasing row with relative gap \(r\ge5\) supplies
at least \(E(5)/2>5/4\), where \(E\) is the endpoint
derivative bound of that criterion. The comparison follows
directly from
\[
 E(5)=8\log(7/2)-40/7+20/49-\log6,\quad
                     \log(7/2)>5/4,\quad \log6<9/5 .
\]
Hence \(w(\log(BQ_\pm))'\) increases. This proves the
endpoint maximum in \(w\), and also proves the endpoint
maximum of \(B(2/m)Q_\pm(2/m)\) as a function of \(m\ge4\).
For the latter, multiply its logarithmic derivative by
\(m\), giving \(-w(\log(BQ_\pm))'\) at \(w=2/m\).

In either orientation,
\[
                         2Q_\pm(1/2)<1/16.                   \tag{4}
\]
For a row of gap \(d\) and ratio \(r\ge5\), define
\[
 h(r)=\frac{(r+2)\log(1+r/2)-(r+1)\log(1+r)}r,\qquad
 g(r)=\tfrac12\log(1+r/2)-\tfrac14\log(1+r).
\]
Its midpoint logarithm is \(d h(r)+g(r)\).
The numerator defining \(h(r)\) vanishes at zero
and has decreasing derivative
\(\log((1+r/2)/(1+r))\); hence \(h\) decreases.
The function \(g\) increases with \(g'(r)<1/[4(1+r)]\).
Since \(r\le d\), the logarithm is at most \(d h(5)+g(d)\).
It decreases in \(d\ge5\), using \(h(5)<-1/4\).
The largest value is attained at \(d=r=5\).
Thus both midpoint factors are bounded by those of the
single shared mirrored pair \((1,6),(6,1)\) together.
For that pair \(Q(1/2)^2=(7/2)^{30}/6^{25}<1/32^2\).
The unused curvature only decreases \(Q\). This proves (4).

## Separation at every order

Write \(N=n+1\), \(m=n,N\), \(D=Nq\), and
\(F_m(z)=L_m(z)+z^m U_m(1/z)\), assigning the middle
coefficient, if any, to just one half. Terminal two
separates by positivity of the quadratic on the real line,
also preserved by Gaussian averaging. Thus assume \(n\ge2\).

For \(N\le256\), first bound \(mQ_\pm(1/m)\) at the
endpoints \(m=2,256\), by the original oriented criterion.
For the coefficients of index \(k\ge2\), apply the entropy
endpoint maximum first on \([2/m,1/2]\), then at
\(m=4,256\) for \(B(2/m)Q_\pm(2/m)\).
The midpoint checks are (4). At \(m=256\) use \(B(w)\le e/w\)
for the second check and use the exact binomial at \(k=1\).
For \(j=1,2\) the required bound at this endpoint uses
\[
 q e^{jL/256}\le \frac1{7776}e^{j(45/2)/256}.
\]
This follows by combining (1) and the second bound
in (2). The logarithm on the right of this combined
majorant decreases in \(d,e_0,s_0\), since \(d,e_0\le14\)
and \(j\le2\). These endpoint bounds give \(c_{m,k}\le8^{-k}\)
throughout each half and \(c_{m,1}<1/16\).
Thus on radius \(e^{1/(2m)}\) their geometric ratio is below
\(1/7\), using the stronger first coefficient for \(m=2,3\).
The short phase comparison of the
[curvature sixteen proof](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md#the-short-ranges)
applies.

For \(N>256\) and \(D\le4/25\), we have
\(L/n\le2/3\) by \(V<222\) and \(n\ge256\), and
\(LD/n<(323/100000)\) by (3).
For each \(1\le k\le4\), the expression
\(D\exp(k\min(2/3,(323/100000)/D))\)
has its maximum at the switch of the two bounds or
at \(D=4/25\). In fact the second endpoint dominates
for these four indices. The normalized first three
coefficients are therefore bounded by
\[
           \frac{\big((4/25)\exp((323/100000)k/(4/25))\big)^k}{k!}.
\]
The endpoint bound at \(k=4\) is below \(1/8\) by
\(eD\exp(4L/n)/4\). Together with (4) the entropy
endpoint maximum gives \(8^{-k}\) for \(k\ge4\).
On radius \(501/500\), the three explicit powers
and geometric tail give \(|H_m-1|<183/1000\) and
\(|zH_m'|<201/1000\) for either half.
Thus the same short phase proof applies at \(e^{1/(2m)}\).

In the remaining case \(D>4/25\), we have \(n\ge1244\).
Put \(K=\lceil32D\rceil+14\le126D\).
The budgets in (3) give
\[
 \frac{D}{n}<\frac{129}{10^6},\quad
 \frac{WD}{n}<\frac{657}{10^5},\quad
 \frac{L_0D}{n}<\frac{644}{10^6},\quad
 \frac{LD}{n}<\frac{323}{10^5}.
\]
In particular \(L_0K/n<82/1000,\ K/n<163/10000\).
At \(k=K\) use \(e\exp(407/1000)/32<1/6\) and
at the midpoint use (4). The entropy maximum gives
\(c_{m,k}\le6^{-k}\) from \(K\) to the end of each half.

For a positive starting endpoint the Gamma derivative
bound for the real-order first coefficient \(B_m\) has
numerator \(5V+S(3+\log N)\le1110+16(3+\log N)\).
Since \(N\ge1245\), this divided by \(n\) is less than \(3/2\).
For example, its ratio to \(e^R-1\) decreases for \(R\ge7\)
and is below \(3/2\) at \(R=7\).
This gives \(|(\log B_m)'|<5/(2n)\) and
\(B_m\le De^{L/n+5/(2n)}<(103/100)D\).
The Gaussian term contributes \(-\tau/(2m^2)\)
to the derivative and is included in the same bound.

The
[sharper Gamma and differentiated connected estimates](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md#the-long-range-with-sharper-gamma-derivative-estimates)
now apply with \(\kappa=82/1000\),
\(\alpha=L_0/[(1-\kappa)n]\), and
\(H_2=W(1+K/[(1-\kappa)n])/[(1-\kappa)n]\).
The same argument includes the Gaussian quadratic
in the endpoint logarithm. In this case
\[
 H_2< (111/100)W/n<46/1000,\quad \alpha K<9/100,\quad
 \beta K<159/100,\quad \beta B_m<13/1000 ,
 \qquad \beta=e^{H_2+1/2}H_2.
\]
With the same tree base bound \(18/5\), the connected
parameter \(p_m=\Xi B_m\) remains below \(3/25\).
The first-coefficient derivative bound is stronger
than the one used at absolute curvature sixteen.
Thus the same logarithm and phase estimates apply
to all positive starting halves, with individual
bound \(2D\) and difference across orders below \(1/100\).

For a private starting half the factorial cancellation
gives coefficient-root bound \(e^{3/2}Q_\pm(k/m)\).
Up to \(K\) the ratio on radius \(6/5\) is below
\((6/5)e^{3/2+407/1000}/7776<1/500\).
Beyond \(K\) the ratio is \(1/5\).
The geometric sums give a logarithm of modulus less
than \(1/240\) on this disk, and so the same individual
and paired log bounds as above.
The
[final paired-radius phase comparison](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md#the-long-range-with-sharper-gamma-derivative-estimates)
uses only these two bounds and \(D/n<337/10^6\).
It excludes a common argument. This completes the proof.

## Exact rational scalar checks

~~~python
from fractions import Fraction as F
from math import factorial

def lo(x):
    x = F(x)
    return sum(x**j/factorial(j) for j in range(81))

def hi(x):
    x = F(x)
    return lo(x)+x**81/factorial(81)/(1-x/82)

assert hi(F(11,16))<2
assert lo(F(9,5))>6 and hi(F(5,4))<F(7,2)
assert (8*F(5,4)-F(40,7)+F(20,49)-F(9,5))/2>F(5,4)
assert F(7,2)**7*hi(F(5,4))<6**6
assert F(7,2)**30/6**25<F(1,32)**2
assert F(256,7776)*hi(F(45,512))<F(1,16)
assert F(128,7776)*hi(1+F(45,256))<F(1,8)
assert hi(F(1,8))<F(8,7) and hi(F(1,4))<F(4,3)

a,b,c,g=F(129,10**6),F(657,100000),F(644,10**6),F(323,100000)
d=F(4,25)
small=g/F(2,3)
assert F(3,4)*222/256<F(2,3)
assert F(257,256)*F(25,7776)<g
assert small*hi(4*(F(2,3)-g/d))<d
assert hi(F(1,512))<F(501,500)
assert d*hi(1+4*g/d)/4<F(1,8)
v=F(501,4000)
head=[(F(501,500)*d)**k*hi(g*k*k/d)/factorial(k)
      for k in range(1,4)]
assert sum(head)+v**4/(1-v)<F(183,1000)
assert sum(k*t for k,t in enumerate(head,1))+v**4*(4-3*v)/(1-v)**2<F(201,1000)
assert F(201,817)<F(1,4) and lo(F(1,2))>F(1183,817)

n=1244
assert d*7776>n
assert F(n+1,n)/7776<a
assert F(n+1,n)*F(51,7776)<b
assert F(n+1,n)*F(5,7776)<c
assert F(n+1,n)*F(25,7776)<g
assert 126*c<F(41,500) and 126*a<F(163,10000)
assert 126*g<F(407,1000)
assert F(1)/(1-F(41,500))*(1+126*a/(1-F(41,500)))<F(111,100)
assert F(41,500)/(1-F(41,500))<F(9,100)
assert hi(7)<1245 and F(1110+16*10,1)/(lo(7)-1)<F(3,2)
assert hi(g/d+F(5,2*n))<F(103,100)
assert F(111,100)*b/d<F(23,500)
beta=F(111,100)*hi(F(273,500))
assert beta*126*b<F(159,100)
assert beta*b*F(103,100)<F(13,1000)
alpha=F(9,100)
assert F(18,5)*F(13,1000)/(1-alpha)*hi(F(159,100)*(F(1,2)+alpha/(6*(1-alpha))))<F(3,25)
assert hi(1+F(407,1000))/32<F(1,6)
assert F(6,5)*hi(F(3,2)+F(407,1000))/7776<F(1,500)
assert F(1,499)+F(1,5)**20/(1-F(1,5))<F(1,250)
assert F(5,2)+4/lo(1)<7
assert F(129,10**6)<F(337,10**6)
~~~
