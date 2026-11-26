# All-degree separation at large signed degree difference

Let \(a,b\) be distinct nonnegative integer exponent columns of
an ordinary factorial binomial in any dimension. Define
\(\Delta=\sum_i(b_i-a_i)\), and write \(Q_j\) for its degree-\(j\)
moment polynomial. Take \(J\ge4\) at least the
[computable absolute eventual angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
If
\[
                             |\Delta|\ge146J-73,
\]
then \(Q_n(z)\) and \(Q_{n+1}(\lambda z)\) have no common
complex root for every \(n\ge1\) when \(\lambda\ge1\) for
\(\Delta>0\), or \(0<\lambda\le1\) for \(\Delta<0\).
In particular the original profile satisfies SFC(2) at every
order over every characteristic-zero field. The bound on
\(\Delta\) allows shared rows of either orientation and any
allocation of the private rows. The numerical value of \(J\)
is unevaluated.

## An explicit linked radial condition

Delete equal rows. Let \(S=\sum_i(b_i-a_i)^2/\max(a_i,b_i)\)
be the total curvature, \(P\) the sum of private exponents,
and \(V\) the sum of \((b_i-a_i)^2/\min(a_i,b_i)\) over
shared rows. Put \(R_P=\max(100,P\max(4,P))\).
The [bounded-private annulus](MULTIVARIATE_BOUNDED_PRIVATE_ROOT_ANNULUS.md)
bounds the absolute log modulus of each endpoint-normalized
root by \(P\log R_P+V/2\). For \(P=0\) use the shared bound
from that result.

Write \(G_j=([z^j]Q_j/Q_j(0))^{1/j}>0\) for the positive
endpoint scale. The normalized root associated to a physical
root \(z\) is \(G_jz\). If two roots match at \(z,\lambda z\)
with \(N=n+1\), the annulus gives
\[
                   \big|\log(\lambda G_N/G_n)\big|
                                   \le 2P\log R_P+V .
\]
The two endpoint scales satisfy
\[
 \left|\log(G_N/G_n)-\Delta\log(N/n)\right|
                     \le \frac{V+P\log(1+N)}{nN}.               \tag{1}
\]
Together these inequalities give the necessary condition
\[
 |\log\lambda+\Delta\log(N/n)|
         \le 2P\log R_P+V+\frac{V+P\log(1+N)}{nN}.              \tag{2}
\]
This condition applies at any order, including when
\(\log\lambda\) and \(\Delta\) have opposite signs.

To prove (1), extend the order continuously to \(j>0\).
A row with entries \(a,b\) contributes to \(\log G_j\) by
\[
    \frac{\log\Gamma(1+jb)-\log\Gamma(1+ja)}{j}
       =\int_a^b\psi(1+jx)\,dx .
\]
The derivative is the signed integral of
\(x\psi'(1+jx)\). For \(u>0\), comparison of the trigamma
series \(\sum_{k\ge1}(u+k)^{-2}\) with its upper and lower
integrals gives
\[
                  0\le1-u\psi'(1+u)\le\frac1{1+u}.
\]
The same inequalities hold at \(u=0\) by continuity.
A shared row of gap \(d=|b-a|\) and smaller exponent
\(\ell\ge1\) therefore contributes \( (b-a)/j\) to the
order derivative with error at most
\(d/(j^2\ell)\le d^2/(j^2\ell)\).

For an increasing private row of exponent \(d\), the
derivative is exactly
\[
         \frac d j-\frac1{j^2}\int_0^{jd}(1-u\psi'(1+u))\,du .
\]
The decreasing orientation gives its negative. The
error relative to the signed \(d/j\) is bounded by
\(\log(1+jd)/j^2\). Since
\(\log(1+Nd)\le d\log(1+N)\) for integral \(d\ge1\),
summing the errors and integrating \(j\) from \(n\) to
\(N\) proves (1).

## The physical grid at fixed reflected mean

Assume all unequal rows are shared. At lower order \(n\), let
\(U_\pm\) be the sums of the row outer losses
\[
 U_i=\frac12\sum_{h=1}^{d_i}
     \log\left(1+\frac{(n-1)d_i}{n\ell_i+h}\right)
\]
in the increasing and decreasing orientations, respectively.
Thus \(U_++U_-=-\log u\) in the endpoint-normalized
reflected-mean convention. Put \(L_n=\log(2n/(n+1))\).

**Proposition.** The physical comparison scale obeys
\[
 \Delta\log(N/n)-\frac{U_+}{nN L_n}
       \ \le\ \log(G_N/G_n)\
       \le\ \Delta\log(N/n)+\frac{U_-}{nN L_n}.                  \tag{3}
\]
In particular, a physical common root requires the ratio of
the two normalized root moduli to lie in this interval
after taking logarithms. Only the integer grid
\(\Delta\log(N/n)\) contributes to the central term.

For the proof, the trapezoid inequality for the convex function
\((u+x)^{-2}\) gives, for \(u>0\),
\[
  \int_0^\infty\frac{dx}{(u+x)^2}
       \le\frac1{2u^2}+\sum_{k\ge1}\frac1{(u+k)^2}.
\]
The decreasing-function integral comparison gives the other
inequality, so \(0\le1-u\psi'(1+u)\le1/(2u)\).
The correction to the derivative of \(\log G_j\) for an
increasing row is therefore nonpositive and has magnitude
at most \(d_i/(2j^2\ell_i)\), by the shared-row integral above.
For a decreasing row the sign reverses.
On the other hand,
\[
 U_i\ \ge\ \frac{d_i}{2}
     \log\left(1+\frac{n-1}{n\ell_i+1}\right)
 \ \ge\ \frac{d_i}{2\ell_i} L_n .
\]
The first inequality uses \(h\le d_i\) and
\(d_i/(n\ell_i+d_i)\ge1/(n\ell_i+1)\).
For the second, \(\ell\log((n\ell+n)/(n\ell+1))
=\int_1^n \ell/(n\ell+x)\,dx\) increases with \(\ell\ge1\).
Integrating \(j^{-2}\) from \(n\) to \(N\) proves (3).

The signed slot information can give a one-sided refinement.
With \(T_{n,h}=(n-1)d_i/(2(n\ell_i+h)+(n-1)d_i)\),
write \(P=\tfrac12\sum_i\epsilon_i\sum_h T_{n,h}^2\), with
\(\epsilon_i\) positive on increasing rows. If
\(p_{\max}>0\) satisfies
\(T^2/(2\operatorname{artanh} T)\le p_{\max}\) for every slot,
then \(P\le p_{\max}U_+\).
Hence when \(s=P/(p_{\max}(-\log u))\ge0\), one has
\(U_-\le (1-s)(-\log u)\), giving the smaller upper error
in (3). Reversing the columns gives the symmetric conclusion.

## The absolute scalar bound

A candidate with \(N\ge J\) is already excluded by the
eventual angular theorem. The first and second moments have
disjoint root arguments by strict positivity of the quadratic.
We prove a parametric bound for an integer \(K\ge11\) whenever
all profiles of curvature \(S\ge K\) have angular separation.
The [absolute curvature twelve theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
gives \(K=12\). It remains to consider \(n\ge2\), \(N<J\),
and \(S<K\).

For a shared row put \(r=d/\ell\). Its contributions to
\(S\) and \(V\) are \(\ell r^2/(1+r)\) and \(\ell r^2\),
so the latter is \(1+r\) times the former.
If \(r>K-1\) and \(\ell\ge2\), the row alone has curvature
above \(2(K-1)^2/K>K\).
If \(\ell=1\), integrality and \(S<K\) leave only \(d=K\)
among gaps above \(K-1\), since the next gap gives
\((K+1)^2/(K+2)>K\).
Such a row contributes \(K^2/(K+1)>K-1\), so it
occurs at most once and is incompatible with \(P\ge1\).

If this row occurs, the remaining curvature is below
\(K/(K+1)<1\). Every other shared row then has \(r\le1\):
\(r>1\) requires \(d\ge2\) and would give
\(S_i=dr/(1+r)>d/2\ge1\).
Thus \(V<K^2+2\). If the row is absent, the bound
\(r\le K-1\) gives \(V\le K(K-P)\). With \(B=K^2+2\),
both cases give
\[
                       2P\log R_P+V\le B-P,\qquad V\le B.
\]
For \(P\ge1\), use \(P\le K-1\) and
\(2\log R_P\le4\log(K-1)\le K-1\).
The last inequality holds at \(K=11\) since \(\log10<5/2\),
and persists by differentiation for larger \(K\).

Under the stated sign condition on \(\lambda\),
\(|\log\lambda+\Delta\log(N/n)|\ge|\Delta|\log(N/n)\).
Now \(\log(N/n)>1/N\), \(\log(n+2)/n<1\) for \(n\ge2\),
and (2) imply
\[
  |\Delta|<N(B-P)+B/n+P
       \le BN+B/2\le B(J-1/2).
\]
Thus for any such integer \(K\), signed degree difference at least
\(\lceil (K^2+2)(J-1/2)\rceil\) suffices. At \(K=12\) this is
\(146J-73\).
