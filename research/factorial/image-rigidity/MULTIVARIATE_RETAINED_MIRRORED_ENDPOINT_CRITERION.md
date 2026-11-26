# A rational endpoint criterion from a retained mirrored subset

Choose any mirrored subset of the rows of an ordinary binomial.
Write its pairs as \((\ell_i,\ell_i+d_i),(\ell_i+d_i,\ell_i)\),
where \(\ell_i,d_i\ge1\). All remaining rows are arbitrary,
including private rows and asymmetric orientations. Let \(s_0\)
be the absolute curvature of those remaining rows and define
\[
 Q(w)=e^{-s_0(1-w)/2}
       \prod_i \left(\frac{\ell_i+d_iw}{\ell_i+d_i}\right)^{d_i},
 \quad q=Q(0),\quad
 W=1+\sum_{\rm shared}\frac{(a_i-b_i)^2}{\min(a_i,b_i)}
                 +\sum_{\rm private}\max(a_i,b_i).
\]
If
\[
          Q(1/30001)\le10^{-6},\qquad Q(1/2)\le1/55,\qquad
                                      Wq\le10^{-4},           \tag{1}
\]
the full ordinary profile has disjoint consecutive root arguments
in every degree, hence satisfies SFC(2) over every
characteristic-zero field. The same conclusion holds with a
physical Gaussian addition of strength \(\tau\ge0\) upon adding
\(\tau\) both to \(s_0\) and to \(W\).

The exponential factors in (1) can be bounded by directed rational
Taylor sums. A purely rational sufficient test is obtained by
replacing \(e^{-s_0(1-w)/2}\) at each tested \(w\) by any certified
rational upper bound. In particular one can use
\(\big(\sum_{j=0}^M (s_0(1-w)/2)^j/j!\big)^{-1}\)
at \(\tau=0\). The criterion always includes a designated
mirrored pair with \(d\ge7\ell\), and it can also use multiple
smaller retained pairs.

## The two coefficient estimates

Put \(V=W-1\), \(S=s_0+2\sum_i d_i^2/(\ell_i+d_i)\), and
\(L_*=\sum_i d_i^2/\ell_i+s_0/2\). Then
\(S\le V\le S(S+2)\) and \(L_*\le V/2\).
In either orientation let \(c_{m,k}\) be the coefficient of the
full endpoint-normalized order-\(m\) polynomial, including
the Gaussian factor if present. For \(w=k/m\le1/2\), the exact
mirrored row ratio on each retained pair and the ordinary
curvature bound on the rest give
\[
 c_{m,k}^{1/k}\le eQ(w)/w,\qquad
 c_{m,k}\le (m q)^k e^{L_* k^2/m}/k! .                         \tag{2}
\]
Indeed the retained pair contributes at most
\(((\ell_i+d_i w)/(\ell_i+d_i))^{d_i k}\). The other rows
contribute at most \(e^{-s_0 k(1-w)/2}\), by the
[discrete endpoint estimate](MULTIVARIATE_HIGH_CURVATURE_ANGULAR_SLOPE.md#an-endpoint-coefficient-estimate).
Combining with either binomial bound and
\(\log(1+(d_i/\ell_i)w)\le (d_i/\ell_i)w\) proves (2).

Writing \(r_i=d_i/\ell_i\), differentiation of
\(\log(eQ(w)/w)\) and multiplication by \(w\) gives
\(-1+\sum_i d_i r_iw/(1+r_iw)+s_0w/2\), an increasing
function of \(w\). Thus interval maxima lie at endpoints.
Similarly \(m\) times the derivative of
\(\log(mQ(1/m))\) is \(1-\sum_i d_i r_i/(m+r_i)-s_0/(2m)\),
which is increasing in \(m\).

## Short halves

For \(N=n+1\le30001\), test \(w=1/m\) and \(w=1/2\) for
each \(m=n,N\), \(m\ge2\). The bound at \(w=1/m\) maximizes
at \(m=2\) or \(m=30001\). The two values from (1),
\(2eQ(1/2)\) and \(30001eQ(1/30001)\), are below \(1/10\).
Consequently \(c_{m,k}\le10^{-k}\) on the interior of each
oriented half. If \(N>30001\) and \(S\ge2\log N+5\), the
large-curvature bound gives \(8^{-k}\).

In the remaining regime \(S<2\log N+5\), let \(D=Nq\).
With \(R=\log N>10\), the bound
\[
        L_*/n\le
        ((2R+5)^2+(2R+5)/2)/(e^R-1)<1/32
\]
follows from \(V\le S(S+2)\) as in the
[short endpoint estimate](MULTIVARIATE_MIRRORED_HARD_PAIR_INTERLACING.md#direct-bounds-in-the-short-endpoint-ranges).
When \(D\le1/10\), (2) gives \(c_{m,k}\le (1/8)^k/k!\)
through \(k=7\). At \(k=8\) the first bound of (2) is at
most \((e D/8)e^{8L_*/n}<1/10\) and at \(w=1/2\) it
is below \(1/10\) by (1), giving this geometric bound
through the rest of each half.

These are exactly the three coefficient bounds of the
[oriented short-half phase proof](MULTIVARIATE_HARD_MIRRORED_PAIR_ARBITRARY_TAIL.md#short-halves-and-large-absolute-curvature).
It uses \(R_m=e^{1/(2m)}\) and bounds each half's logarithmic
angular derivative by \(1/4\). Its Rouché and phase argument
therefore excludes a common root argument in these regimes,
independently of reflection symmetry.

## Long halves

Assume \(D>1/10,\ S<2\log N+5\). By \(q\le Q(1/30001)\le10^{-6}\)
and \(Wq\le10^{-4}\) we have \(n\ge100000\) and, with
\(K=\lceil32D\rceil+14\),
\[
 D/n<1/720000,\quad WD/n<1/9000,\quad
 K\le182D,\quad WK/n<182/9000,\quad K+1<n/2 .
\]
The two endpoints in (2) at \(k=K,m/2\) give
\(c_{m,k}\le8^{-k}\) for \(K\le k\le m/2\), since
\((e/32)e^{L_*K/n}<1/8\) at the first endpoint.

For a positive starting endpoint, the
[general first-coefficient derivative bound](MULTIVARIATE_ABSOLUTE_CURVATURE_27.md#the-endpoint-cluster-region)
at real \(m\in[n,N]\) gives
\[
          |(\log B_m)'-1/m|
                \le (5V+S(3+\log N))/m^2<1/n .
\]
The last bound follows from \(S<2\log N+5\), \(\log N>11\),
and \((5V+S(3+\log N))\le40(\log N)^2<n\).
The Gaussian contribution to this derivative is
\(-\tau/(2m^2)\), included in the augmented bound.
At \(m=n\), (2) gives \(B_n\le D e^{L_*/n}\). Thus
\(B_m\le D e^{(L_*+2)/n}<101D/100\).
The logarithmic endpoint interpolation after dividing by
\(B_m^k/k!\) gains only a term \(\tau k(k-1)/(2m)\)
from the Gaussian, and its second finite difference and
order derivative satisfy the augmented edge bounds.
Consequently the preceding long \(D\) inequalities are the
inputs of the
[paired endpoint cluster estimate](MULTIVARIATE_MIRRORED_HARD_PAIR_INTERLACING.md#a-paired-estimate-for-the-two-long-endpoint-sums).
As in the
[asymmetric hard-pair extension](MULTIVARIATE_HARD_MIRRORED_PAIR_ARBITRARY_TAIL.md#long-halves-and-the-change-of-order),
\(B_m/m\) need not be monotone for this use of that estimate.
The resulting half logarithm bounds on radius \(6/5\) are
\[
          \|\log H_m\|\le3(D+1),\qquad
                              \|\log H_N-\log H_n\|<1/55 .    \tag{3}
\]

If the half has a private starting row, it lies outside the
retained subset. The private bound
\(\binom mk k!/(m!)^{k/m}\le e^k\) consumes one unit of residual
curvature and eliminates the binomial coefficient. Therefore
\[
                c_{m,k}^{1/k}\le e^{1+(1-k/m)/2}Q(k/m).
\]
For \(k\le K\), \(Q(k/m)\le q e^{L_*K/n}\) by bounding
the logarithmic derivative of \(Q\) by \(L_*\).
At radius \(6/5\) the resulting geometric ratio is less
than \(10^{-5}\), using \(q\le10^{-6}\) and
\(L_*K/n\le WK/(2n)<1/90\).
The remaining terms on the same disk start beyond \(K\)
with geometric ratio \(3/20\).
Since \(K\ge18\), the sum gives \(|H_m-1|<1/200\).
This implies \(|\log H_m|<1/190\) and establishes (3) also
for the private orientation, which is the same at both orders.

The
[paired radius and phase argument](MULTIVARIATE_HARD_MIRRORED_PAIR_ARBITRARY_TAIL.md#the-two-radii-do-not-affect-the-phase-contradiction)
uses only (3) and \(D/n<1/720000\). It excludes a common
root argument in the long regime. The first-to-second pair
separates by positivity of the quadratic moment on the real
line, also preserved by Gaussian lognormal averaging.
This proves the all-degree assertion.

For the single designated pair \(d\ge7\ell\), the first
two conditions of (1) follow from \(2/8^7<10^{-6}\)
and \((9/16)^7<1/55\). For the last, put
\(w_0=1+2d^2/\ell\ge99\). Since
\(W\le w_0+s_0(s_0+2)\), the bound
\((w_0+s_0(s_0+2))e^{-s_0/2}\) decreases in \(s_0\ge0\).
Thus \(Wq\le(1+2d^2/\ell)8^{-d}\le99/8^7<10^{-4}\).
