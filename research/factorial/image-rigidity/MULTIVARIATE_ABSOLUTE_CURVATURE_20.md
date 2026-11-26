# All-degree angular separation at absolute curvature twenty

For arbitrary distinct ordinary exponent columns \(a,b\), put
\[
 S=\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)}+\tau ,
\]
omitting equal rows, where \(\tau\ge0\) is an optional common
Gaussian strength. Its multiplier on the order-\(m\) normalized
coefficient at \(k\) is \(e^{-\tau k(m-k)/(2m)}\).

**Theorem.** If \(S\ge20\), the consecutive moment polynomials
have disjoint root arguments at every order, including after
independent positive rescaling. At \(\tau=0\) this proves
SFC(2) over every characteristic-zero field. The support can
have private variables and arbitrary orientations.

The proof uses the exact oriented endpoint cost of every
shared row. It absorbs the excess of the Gamma curvature
parameter over \(S\). A sharper connected logarithm bound then
gives the adjacent comparison using no tests of individual
moment polynomials.

## The endpoint cost absorbs the excess curvature

Use the two endpoint envelopes \(Q\) from the
[oriented endpoint comparison](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md),
retaining every shared row. Let \(q\) be the larger of
their two values at zero. Let \(V\) be the sum of the
squared gaps divided by the smaller exponent on shared
rows, of the private exponents, and of \(\tau\).
Put \(W=1+V\) and \(C=1+S\). We prove
\[
        q\le e^{-S/2},\qquad
                  Wq\le(1+S)e^{-S/2}.                         \tag{1}
\]

On a shared row its weaker oriented endpoint cost
\(-\log Q_i(0)\) is \(\ell(r-\log(1+r))\), where
\(\ell=\min(a,b)\) and \(r=|b-a|/\ell\).
Indeed this is the decreasing orientation by
\(\log(1+r)\ge2r/(2+r)\).
The absolute curvature is \(S_i=\ell r^2/(1+r)\),
and the excess \(V_i-S_i=\ell r^3/(1+r)\).
Let \(E(r)=r-\log(1+r)-r^2/[2(1+r)]\).
Then \(E(0)=0\) and \(E'(r)=r^2/[2(1+r)^2]\).
For \(C\ge21\) this derivative is at least the derivative
of \(\log(1+r^3/[C(1+r)])\). After multiplying by the
positive denominators, the required bound is
\[
 C(1+r)+r^3\ \ge\ 2(3+2r)(1+r).
\]
At \(C=21\) the difference is \(r^3-4r^2+11r+15>0\),
whose derivative is positive on \(r\ge0\).
Integrating and using the elementary product inequality gives
\[
 \ell E(r)\ge \ell\log(1+r^3/[C(1+r)])
                 \ge \log(1+(V_i-S_i)/C).
\]
Private rows and the Gaussian have cost \(S_i/2\) and
zero excess. Summing over rows yields
\(-\log q\ge S/2+\log(1+(V-S)/C)\), which proves (1).

We use also the coefficient envelopes from the same oriented comparison.
With \(w=k/m\), \(L=3V/4\), they are
\[
 c_{m,k}^{1/k}\le(e/w)Q(w),\quad
 c_{m,k}\le(mq)^k e^{L k^2/m}/k!,\quad
 Q(w)\le q e^{Lw},\quad Q(1/2)\le e^{-S/4}.                   \tag{2}
\]
The first coefficient envelope takes its maximum at an endpoint
of any closed \(w\)-interval in \((0,1/2]\).

## Short endpoint halves

Write \(N=n+1\), \(n\ge2\), and let \(m=n,N\).
Use the two oriented half-polynomials \(F_m(z)=L_m(z)+z^m U_m(1/z)\).
If \(N\le256\), the [discrete curvature bound](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)
gives \(c_{m,k}\le20^{-k}\) for \(k\le m/2\).
Indeed the scalar endpoint bounds are
\(em e^{-10(1-1/m)}\) and \(2e e^{-5}\); the first has
no interior maximum for \(2\le m\le256\), and both are
less than \(1/20\).
For \(N\ge256\) and \(S\ge2\log N+5\), use the
[large-curvature geometric estimate](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate)
\(c_{m,k}\le8^{-k}\), which uses the same coefficient
bound with \(S\) including the Gaussian.
In these cases at radius \(R_m=e^{1/(2m)}\) the geometric
ratio is strictly below \(1/7\), giving
\(|H_m-1|<1/6\) and \(|zH_m'/H_m|<7/30<1/4\)
on the disk for both halves \(H_m=L_m,U_m\).

Otherwise \(S<2\log N+5\) and \(n>1807\).
Indeed \(N>e^{(S-5)/2}\ge e^{15/2}>1808\).
Consequently \(L/n<183/1000\), by \(V\le S(S+2)\) and the
decrease of \(3S(S+2)/[4(e^{(S-5)/2}-1)]\) for \(S\ge20\).
Put \(D=Nq\). When \(D\le2/25\), on \(|z|\le R_m\)
the first seven half coefficients are majorized by
\[
              (101/1250)^k e^{183 k^2/1000}/k!.
\]
Beyond this use \(8^{-k}\), by the endpoint test in (2)
at \(k=8\) and at the midpoint. The respective bounds there
are at most \(e^{1+8(183/1000)}/100<1/8\) and \(2e^{-4}<1/8\).
The replay checks both geometric tails on this disk,
and shows that the total for \(|H_m-1|\) and for \(|zH_m'|\)
is below \(3/25\) in each case.
Thus \(|zH_m'/H_m|<1/4\) also in this region.

In the short cases the
[phase comparison on the short annulus](MULTIVARIATE_HARD_MIRRORED_PAIR_ARBITRARY_TAIL.md#short-halves-and-large-absolute-curvature)
now applies. In fact the half norms just proved give
Rouché at \(R_m\) and at its reciprocal, since \(R_m^m=\sqrt e\)
exceeds \(7/5\). The phase error from both halves has modulus
less than \(\theta/2\) at any root of argument \(0<\theta\le\pi\).
Subtracting the two phases at a common argument would put
an integer multiple of \(2\pi\) strictly between \(0\) and
\(2\theta\le2\pi\).

## A smaller differentiated connected logarithm bound

It remains to treat \(D>2/25\) in the same band \(S<2\log N+5\).
Put \(K=\lceil32D\rceil+14\). By (1), using \(n>1807\)
and the decrease of \((S+1)e^{-S/2}\) on \(S\ge20\), we have
\[
 \frac{D}{n}<\frac{46}{10^6},\qquad
 \frac{WD}{n}<\frac{955}{10^6},\qquad
 K\le220D,\qquad
 \frac{WK}{n}<\frac{421}{2000}.
                                                               \tag{3}
\]
In particular \(K+1<n/2\).
At \(k=K\) the first coefficient bound in (2) is at most
\((e/32)\exp((3/4)WK/n)<1/8\). At the midpoint its
bound is at most \(2e^{-4}<1/8\), giving the whole geometric
tail \(c_{m,k}\le8^{-k}\) for \(k\ge K\).

For a positive starting endpoint let \(B_m\) be the first
factor interpolated at real orders \(m\in[n,N]\).
The [Gamma first-factor derivative bound](MULTIVARIATE_HARD_MIRRORED_PAIR_ARBITRARY_TAIL.md#long-halves-and-the-change-of-order)
gives
\[
 |(\log B_m)'-1/m|\le (5V+S(3+\log N))/m^2 .
\]
The same augmentation by \(\tau\) applies to this inequality.
Here \(R=\log N>15/2\). The numerator is at most
\(22R^2+131R+190\le43R^2\), and
\(43R^2/(e^R-1)<3/2\) on that range.
Thus \(|(\log B_m)'|<5/(2n)\).
Equation (2) at \(m=n,k=1\), (3), and \(D>2/25\)
then give \(B_m\le De^{L/n+5/(2n)}<(51/50)D\).

The [exact endpoint interpolation](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md#coefficient-halves-and-endpoint-interpolation)
uses the logarithm of \(c_{m,k}k!/B_m^k\).
Its finite differences \(\nu_r(m)\) and their order
derivatives are bounded by \(H_r\) and \((2r/n)H_r\), where
\[
 H_r=(r-2)!\alpha^{r-2}H_2,\qquad
 H_2=4W/n,\qquad
              \alpha=2\max(1,\max_i|d_i|/a_i)/n\le2W/n .
\]
The Gaussian adds \(\tau k(k-1)/(2m)\) to this interpolation;
its second difference and derivative are absorbed by the
increase in \(W\). The bounds (3) give \(\alpha K<421/1000\)
and \(K+1<n/2\), ensuring the Gamma argument restrictions.
Also \(H_r\le H_2<6/125\) for \(2\le r\le K\).
The edge weights \(e^{\nu_r}-1\) and their derivatives
therefore have absolute bounds \(e^{H_2}H_r\) and
\((2r/n)e^{H_2}H_r\).

For a connected hypergraph on \(k\) labels, differentiating
its product of weights costs at most \(2X/n\) times the
product of the absolute weight bounds, with \(X\) the sum
of its edge sizes. Use
\[
                   X\le (2/e)k e^{X/(2k)} .
\]
Each edge has size at most \(k\). The factor \(e^{X/(2k)}\)
can thus be absorbed by inflating each edge weight by
\(e^{1/2}\). This gives a common absolute majorant
\((r-2)!\alpha^{r-2}\beta\) with
\(\beta=e^{H_2+1/2}H_2\), and the derivative costs only
an additional factor \(4k/(en)\).

Use the [discovery tree sum](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#sharper-connected-sum)
with \(u=\alpha K\). The tree base
\((1+t)^{(1+t)/t}\) at \(t=\alpha/\beta\le1/3\) is
at most \(256/81\).
For the additional edges use
\[
  \sum_{r=2}^k \binom kr(r-2)!\alpha^{r-2}\beta
          \le \beta k(k-1)\left(\frac12+
                                       \frac{u}{6(1-u)}\right).
\]
Combining these estimates shows that the absolute connected
sum is at most \((k-1)!\Xi^{k-1}\), with
\[
            \Xi=\frac{256\beta}{81(1-u)}
                         \exp\left(\beta K
                                   \left(\frac12+\frac{u}{6(1-u)}\right)\right).
\]
Its derivative is at most \(4k/(en)\) times the same
majorant, by the inflation already incorporated in \(\beta\).
The tree prefactor is at most \(e/k\le1\) for \(k\ge3\);
the cases \(k=1,2\) follow directly.
Let \(\ell_k(m)\) be the formal logarithm coefficients.
They equal \(B_m^k/k!\) times the connected sum. Hence for
\(p_m=\Xi B_m\),
\[
 |\ell_k(m)|\le B_m p_m^{k-1}/k,\qquad
 |\ell_k'(m)|\le4 B_m p_m^{k-1}/n ,
\]
where the derivative constant uses \(5/2+4/e<4\).
From (3) one has \(\beta K<73/50,\ \beta B_m<17/2500\),
and thus \(p_m<1/10\).

## Phase comparison in the long region

Write \(A_m(z)=\sum_{k=1}^K\ell_k(m)z^k\).
The absolute coefficient norm on \(9/5\) is at most \(2B_m\),
by \(10\log(50/41)<2\), and the same bound holds on \(6/5\).
The coefficient derivative bound gives
\(\|A_N-A_n\|_{6/5}<6D/n\).
The actual half and \(e^{A_m}\) agree through degree \(K\).
On radius \(6/5\) their two tails, relative to \(e^{A_m}\),
are bounded by
\[
          (2/3)^{K+1} e^{4B_m}
                       \quad\hbox{and}\quad
          (3/20)^K e^{2B_m}.
\]
The first uses the absolute coefficient norm of \(e^{A_m}\)
on radius \(9/5\); the second uses the factor \(8^{-k}\).
Since \(K\ge32D+14\) and \(D>2/25\), the first is less than
\(e^{-6}\) and the second less than \(e^{-K}\).
The sum is below \(1/130\); the analytic logarithm of the
half differs from \(A_m\) by less than \(1/120\) there.

For a private starting row, reserve one factorial as in
the [oriented private-start estimate](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md#uniform-comparison-across-degrees).
This gives \(c_{m,k}^{1/k}\le e^{3/2}Q(w)\).
For \(k\le K\) the ratio on \(6/5\) is at most
\((6/5)e^{3/2+LK/n}q<1/3000\), and beyond \(K\)
it is at most \(3/20\). Thus in that case
\(\|\log H_m\|_{6/5}<1/190\).
In both cases we obtain
\[
       \|\log H_m\|_{6/5}\le M:=3(D+1),\qquad
             \|\log H_N-\log H_n\|_{6/5}<\delta:=1/55 .
\]
Use the [paired radius and phase identities](MULTIVARIATE_HARD_MIRRORED_PAIR_ARBITRARY_TAIL.md#the-two-radii-do-not-affect-the-phase-contradiction).
Rouché confines root log radii \(u_m\) to \(|u_m|\le7(D+1)/m\),
so the half arguments lie in the disk \(21/20\).
Put \(J_m(u,\theta)=\log L_m(e^{u+i\theta})-\log U_m(e^{-u-i\theta})\).
Cauchy's estimate gives the paired angular derivative bound
\(14\delta\) and the mixed derivative bound \(210M\).
At a common argument the real equations give
\[
 n|u_N-u_n|
       \le \frac{7(D+1)/n+2\delta}{1-14M/n}<1/22 ,
\]
by (3) and \((D+1)/n<(27/2)D/n\).
Both phase functions are real along the positive radius.
Integration of the derivative bounds therefore gives
\[
 |\Im(J_N(u_N,\theta)-J_n(u_n,\theta))|
     \le (14\delta+210(M/n)/22)\theta<\theta/3 .
\]
Here \(M/n<(81/2)(46/10^6)\). Subtracting the two phase
equations at the common argument gives the same contradiction
as in the short case.
Terminal two separates by positivity of the quadratic moment
on the real line, including under the Gaussian averaging.

## Exact rational scalar checks

~~~python
from fractions import Fraction as F
from math import factorial

def lo(x):
    return sum(F(x)**k/factorial(k) for k in range(61))

def hi(x):
    x = F(x)
    return lo(x)+x**61/(factorial(61)*(1-x/62))

assert 2/lo(4) < F(1,20)
assert 256/lo(9-F(10,256)) < F(1,20)
assert lo(F(15,2)) > 1808
assert F(3,4)*20*22/(lo(F(15,2))-1) < F(183,1000)
assert hi(F(1,2)/1807) < F(101,100)
r=F(101,800)
head=sum(k*(F(101,1250)**k)*hi(F(183*k*k,1000))/factorial(k) for k in range(1,8))
tail=r**8*(8-7*r)/(1-r)**2
assert head+tail < F(3,25)
assert hi(1+8*F(183,1000))/100 < F(1,8)
assert F(1808,1807)/lo(10) < F(46,10**6)
assert 21*F(1808,1807)/lo(10) < F(955,10**6)
assert 220*F(955,10**6) < F(421,2000)
assert 220*F(46,10**6)+F(1,1807)<F(1,2)
assert hi(1+F(3,4)*F(421,2000))/32 < F(1,8)
assert 22+131/F(15,2)+190/F(15,2)**2 <43
assert 43*F(15,2)**2/(lo(F(15,2))-1)<F(3,2)
wn=F(955,10**6)/F(2,25)
assert 4*wn<F(6,125)
assert hi(F(3,4)*wn+F(5,2*1807)) < F(51,50)
assert lo(F(1,2))>F(3,2) and F(5,2)+4/lo(1)<4
u=F(421,1000)
inflate=4*hi(F(1,2)+F(6,125))
assert 2*F(421,2000) == u
assert inflate*F(421,2000)<F(73,50)
assert inflate*F(51,50)*F(955,10**6)<F(17,2500)
assert F(256,81)*F(17,2500)/(1-u)*hi(F(73,50)*(F(1,2)+u/(6*(1-u))))<F(1,10)
assert lo(F(1,5))>F(50,41)
assert (4*F(6,5)/(1-F(3,25)))*F(51,50)<6
assert hi(F(2,5))<F(3,2)
assert F(2,5)*(32*F(2,25)+14)-4*F(51,50)*F(2,25)>6
assert hi(F(11,10))<F(20,3) and 2*F(51,50)/32<F(1,10)
assert 1/lo(6)+1/lo(17)<F(1,130) and lo(F(1,120))>F(130,129)
assert F(6,5)*hi(F(3,2)+F(3,4)*F(421,2000))/lo(10)<F(1,3000)
assert F(1,2999)+F(3,20)**17/(1-F(3,20))<F(1,200)
assert lo(F(1,190))>F(200,199)
assert 6*F(46,10**6)+F(2,120)<F(1,55)
M=F(81,2)*F(46,10**6)
R=7*F(27,2)*F(46,10**6)
assert hi(R)<F(21,20)
assert (R+F(2,55))/(1-14*M)<F(1,22)
assert F(14,55)+210*M/22 <F(1,3)
~~~
