# Compactification of exponent profiles at a fixed terminal

At a fixed terminal, arbitrary dimensions and arbitrarily large
exponents can be treated together. The limiting data consist of
multiplicities of fixed unequal exponent rows, a Gaussian variance for
rows escaping to infinity, and a signed correction to the ratio of the
two moment centers. The last correction must be retained even though
the relative gap in every escaping row tends to zero.

We give both directions of this description and explicit truncation
errors. They give a finite algebraic certificate at each truncation
level. A sequence of these certificates terminates whenever the
compactified family is separated. The argument leaves open whether
the compactified family is separated at every terminal.

## Normalized polynomials and row data

Fix \(m\ge3\) and \(N=m+1\). Delete equal rows from the two exponent
columns, since they contribute a nonzero scalar in each order. An
ordered unequal row is denoted by \(w=(a,b)\). Write
\[
 r(w)=\max(a,b),\qquad d(w)=b-a,\qquad s(w)=d(w)^2/r(w).
\]
The full curvature is \(S=\sum_w s(w)\), with multiplicities, and the
total degree difference is \(\Delta=\sum_w d(w)\).
For \(n=m,N\), use
\[
 g(x)=\log\Gamma(x+1),\qquad
 J_{n,k}(w)
   =g(na+kd(w))-(1-k/n)g(na)-(k/n)g(nb).
\]
If \(G_n=(C_{n,n}/C_{n,0})^{1/n}\) is the endpoint center, the
normalized polynomial is
\[
 P_n(z)=\sum_{k=0}^n \binom nk
         \exp\!\left(\sum_w J_{n,k}(w)\right)z^k,\qquad z=G_n t.
 \tag{1}
\]
Here \(t\) is the ratio of the two monomial coefficients. Define
\[
 H_w(x)=\frac{g(bx)-g(ax)}{x},\qquad
 E(w)=H_w(N)-H_w(m)-d(w)\log(N/m).
\]
The two normalizations satisfy
\[
 C=\log(G_N/G_m)=\Delta\log(N/m)+\sum_w E(w).
 \tag{2}
\]
Consequently the original pair is coprime exactly when \(P_m(z)\)
and \(P_N(e^Cz)\) are coprime.

We use the [fixed-terminal curvature bounds](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md)
to restrict the possible values of \(S,\Delta\) at a common root.
Small and large curvature give \(0<s_0\le S\le K\) for suitable
\(s_0>0\) and any \(K>4\log N+10\). The low terminals are already
covered by the [low-terminal proof](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md).
All intermediate coefficients in (1) lie between zero and
\(\binom nk\) by log convexity. Every root of \(P_n\) therefore lies
in \(2^{-n}<|z|<2^n\). For the upper bound, when \(|z|\ge1\) the
root equation gives \(|z|^n\le(2^n-1)|z|^{n-1}\), and the
reciprocal gives the lower bound.

The trigamma estimates used in the endpoint-radius proof also give
a per-row bound
\[
 |E(w)|\le \frac{2\log(N+1)}{mN}s(w).
 \tag{3}
\]
Indeed \(xH'_w(x)=x^{-1}\int_{ax}^{bx} t g''(t)\,dt\), and
\(t/(1+t)\le t g''(t)\le1\). For \(x\in[m,N]\),
the resulting error in \(xH'_w(x)-d(w)\) is at most
\(\bigl|\log((1+bx)/(1+ax))\bigr|/x\).
For \(0\le a<b\) and \(d(w)\ge b/2\), the numerator is at most
\(d(w)\log(1+x)\le2s(w)\log(1+x)\). For \(0<d(w)<b/2\) it is at
most \(\log(b/a)\le d(w)/a\le2s(w)\). Reversing the row preserves
the bound. Integrating \(2s(w)\log(N+1)/x^2\) gives (3).
At a common root, the root annuli bound \(|C|\) by
\((m+N)\log2\). Equations (2) and (3) then give
\[
 |\Delta|<
 \frac{(m+N)\log2+2K\log(N+1)/(mN)}{\log(N/m)}.
 \tag{4}
\]
Thus the integer \(\Delta\) ranges over a finite set independent
of the largest row and the dimension.

## Uniform estimates for rows of large size

Assume \(S\le K\). At a cutoff \(B>K\), retain all rows with
\(r(w)\le B\). For the other rows use
\[
 \eta=\sqrt{K/(B+1)}<1,\qquad
 \tau=\sum_{r(w)>B}s(w),\qquad
 L=\sum_{r(w)>B}\log(b/a).
\]
Each large row satisfies \(|d(w)|/r(w)\le\eta\), so both of its
exponents are positive and \(\min(a,b)\ge(1-\eta)r(w)\).
In particular \(|L|\le\tau/(1-\eta)\), since
\[
 |\log(b/a)|\le |d(w)|/\min(a,b)\le s(w)/(1-\eta).
 \tag{5}
\]

For the coefficient contributions, a second-derivative comparison
on the segment between \(na\) and \(nb\) gives
\[
 \left|J_{n,k}(w)+\frac{k(n-k)}{2n}s(w)\right|
 \le
 \frac{k(n-k)}{2n}s(w)
 \left(\frac{\eta}{1-\eta}+\frac1{n(B+1)}\right).
 \tag{6}
\]
In fact
\[
 \frac1{nr(w)+1}\ \le\ g''(x)\ \le\
                 \frac1{nr(w)(1-\eta)}
\]
on that segment. Since the endpoint interpolation error for a
quadratic of second derivative \(A\) is exactly
\(-A k(n-k)d(w)^2/2\), subtracting the two bounding quadratics
gives (6). This proves the Gaussian correction to every coefficient
with an error that is uniform in the number of large rows.

For the center ratio we need one further term. For \(v>0\), the
trapezoidal sum for \(y\mapsto(v+y)^{-2}\) gives
\[
 \psi'(v+1)=\frac1v-\frac1{2v^2}+R(v),\qquad
                 0\le R(v)\le\frac1{4v^3}.
\]
On each unit interval the trapezoidal error is the integral of
the second derivative against \(t(1-t)/2\), with \(0\le t\le1\).
Summing gives the bound \(\tfrac18\int_0^\infty 6(v+y)^{-4}dy\).
Integrating the formula between \(mx\) and \(Nx\) gives
\[
 0\le \psi(Nx+1)-\psi(mx+1)-\log(N/m)+\frac1{2mNx}
      \le\frac1{4m^3 x^2}.
\]
Now integrate between \(a\) and \(b\), with the actual orientation.
Since \(|d(w)|\le d(w)^2\), summation yields
\[
 \left|\sum_{r(w)>B}E(w)+\frac{L}{2mN}\right|
 \le\frac{\tau}{4m^3(1-\eta)^2(B+1)} .
 \tag{7}
\]
This signed term in (7) can survive a limit of arbitrarily large
rows, independently of the variance term in (6).

## Exact compactified data

**Theorem.** Fix \(0<s_0\le K\) and an integer \(\Delta_0\).
The limit triples \((P_m,P_N,C)\) along ordinary profiles with
\(\Delta=\Delta_0\) and
\(s_0\le\liminf S\le\limsup S\le K\) are exactly the following.
Choose a countable multiset of ordered unequal rows, of
multiplicities \(c_w\in\mathbb Z_{\ge0}\), and choose real
numbers \(\tau\ge0,\lambda\) with
\[
 S_{\rm ret}=\sum_w c_w s(w),\qquad
 s_0\le S_{\rm ret}+\tau\le K,\qquad |\lambda|\le\tau .
 \tag{8}
\]
Then the triple is
\[
 \begin{aligned}
 \mathcal P_n(z)&=\sum_{k=0}^n\binom nk
       \exp\!\left(\sum_w c_w J_{n,k}(w)
                      -\frac{\tau k(n-k)}{2n}\right)z^k ,
                      && n=m,N,\\
 \mathcal C&=\Delta_0\log(N/m)+\sum_w c_wE(w)-\frac{\lambda}{2mN}.
 \end{aligned}                                                     \tag{9}
\]
Both sums converge absolutely. There is no further constraint
on \(\Delta_0\) from the retained multiplicities or \(\lambda\).

To prove necessity, take \(S\to S_*\) along a subsequence. At each
fixed cutoff \(B\) there are finitely many row types, each having a
bounded integer multiplicity since \(s(w)\ge1/B\) there. Diagonal
extraction gives eventual multiplicities \(c_w\) at every finite
cutoff. Thus \(S_{\rm ret}\le S_*\), and the curvature of the
large rows tends first to \(S_*-\sum_{r(w)\le B}c_ws(w)\)
and then to \(\tau=S_*-S_{\rm ret}\).

For the signed correction use the total log ratio over rows with
both entries positive. This quantity is bounded along the
sequence. For large rows use (5), and below a fixed large cutoff
use bounded multiplicities and finitely many types. Extract a
limit for the total. The series
\(\sum_{ab>0}c_w\log(b/a)\) converges absolutely by (5) on the
large rows. Subtract it from the extracted limit to get
\(\lambda\). Inequality (5) at finite cutoffs and then \(B\to\infty\)
give \(|\lambda|\le\tau\).
Equation (3) gives absolute convergence of the \(E\)-sum and (6)
of the \(J\)-sums. Applying (6) and (7) at successively larger
cutoffs proves (9).

For sufficiency take successively larger finite truncations of
the retained multiset. At a large integer \(A\), add
\[
 \left\lfloor\frac{(A+1)(\tau+\lambda)}2\right\rfloor
        \text{ rows }(A,A+1),\qquad
 \left\lfloor\frac{(A+1)(\tau-\lambda)}2\right\rfloor
        \text{ rows }(A+1,A).
\]
Their curvature tends to \(\tau\) and their signed log ratio
tends to \(\lambda\), since
\((A+1)\log(1+1/A)\to1\). The total integer degree difference
of this finite profile can be adjusted to \(\Delta_0\).
If the required correction is \(D\), add \(|D|\) unit-gap rows
of the appropriate orientation with size at least
\(A(1+|D|)^2\). These rows contribute \(O(1/A)\) to both curvature
and log ratio. Because \(s_0>0\), the supports are distinct for
large \(A\). Equations (6) and (7) now show that the constructed
profiles tend to (9). In this direction one needs the displayed
limiting bounds on \(S\). For example, when \(S_{\rm ret}=K\)
a degree adjustment can raise \(S\) above \(K\) by a term tending
to zero. \(\square\)

The [pure diffuse model](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#the-limiting-pair-when-each-relative-gap-tends-to-zero)
is the case of zero retained multiplicities. In that model the
limiting roots lie on the unit circle and a common root forces
\(\mathcal C=0\). Formula (9) permits nonsymmetric finite or
countable retained profiles, whose limiting roots need not have
modulus one.

## A finite certificate for compactified separation

**Theorem.** Fix positive rational numbers \(s_0\le K\), and
allow a finite set of integers \(\Delta\). There is a sequence
of finite tests over the real algebraic numbers with the following
properties.

- A successful test proves ordinary consecutive coprimality on
  the whole range \(s_0\le S\le K\) with the allowed degree
  differences. Neither the dimension nor the largest exponent
  needs to be bounded.
- Some test succeeds if and only if no triple (9) with this
  curvature range and degree set satisfies
  \(\mathcal P_m(z)=\mathcal P_N(e^{\mathcal C}z)=0\)
  for any complex \(z\).

For a fixed \(B>K\), enumerate the finitely many multisets of
rows with \(r(w)\le B\) and \(S_{\rm core}\le K\). In a given
core, put
\[
 \max(0,s_0-S_{\rm core})\le\tau\le K-S_{\rm core}.
\]
Approximate the remaining rows by (9) with this finite core.
For a true profile, replace the signed log ratio \(L\) of the
large rows by \(\lambda=(1-\eta)L\), so \(|\lambda|\le\tau\)
by (5). Use rational relative error bounds tending to zero with
\(B\). In logarithmic coordinates it suffices to use the
following respective bounds for each coefficient and for the
center:
\[
 \frac{NK}{8}\left(\frac{\eta}{1-\eta}+\frac1{m(B+1)}\right),
 \qquad
 \frac{K}{4m^3(1-\eta)^2(B+1)}
     +\frac{\eta K}{2mN(1-\eta)}.
 \tag{10}
\]
The extra term in the center bound accounts for replacing \(L\)
by \(\lambda\). These bounds follow from (6), (7), and
\(k(n-k)/(2n)\le n/8\). If \(0\le\delta<1/2\) is an upper
enclosure for a logarithmic error, admit the multiplicative
error interval \([1-\delta,1+2\delta]\). Take \(B\) large enough
that the two errors admit enclosures of this size.

Put \(u=e^{-\tau/(2mN)}\) and \(v=e^{-\lambda/(2mN)}\).
The variance factors in (9) are now \(u^{N k(m-k)}\) and
\(u^{m k(N-k)}\), and the center factor is
\[
 \exp(\mathcal C)
   =(N/m)^\Delta \exp\!\left(\sum_{\rm core}E(w)\right)v.
\]
All coefficients from the finite core and the center factor
apart from \(v\) are positive real algebraic numbers, given by
rational powers of ratios of factorials and of \(N/m\).
Use rational outer enclosures of the two exponential endpoints
for \(u\), with error tending uniformly to zero, and impose
\(u\le v\le u^{-1}\). The endpoints for \(u\) are bounded away
from zero in terms of \(K,m,N\). They may be enclosed by the
exponential series with its tail bounded by a geometric series.
Error variables for the nonleading coefficients and the center
range over the closed intervals above. The normalized endpoint
coefficients at each order are fixed at one.

For each core and each allowed \(\Delta\), test for the
existence of a common complex root after these perturbations.
Writing real and imaginary parts of the two polynomial equations
gives a finite test of polynomial equalities and inequalities
over the real algebraic numbers. Such a test is decidable by
quantifier elimination over real closed fields.
[Perrucci and Roy, Theorem 1](https://arxiv.org/abs/1609.02879)
If all the tests at a cutoff are infeasible, (6)–(10) imply
coprimality. Over every characteristic-zero field the conclusion
then follows by rational Bézout.

If no cutoff succeeds, choose feasible data with \(B\to\infty\).
Each inexact exponential interval can be clipped to the exact
interval with parameter changes tending to zero. More explicitly,
the logarithm is uniformly Lipschitz on the \(u\)-intervals,
which stay bounded away from zero. Clip \(\tau\) to the
curvature range for its finite core, then clip \(\lambda\) to
\([-\tau,\tau]\), with errors tending to zero in both.
The unperturbed data are now triples (9) with finite cores.
The convergence statements in the compactification theorem show
that their limit triples are again in (9): one can first
approximate each triple by an ordinary finite support with
errors tending to zero and with its curvature in the limiting
range, and then apply the necessity direction. The actual
perturbations also tend to zero by (10). Roots stay bounded
since the monic coefficients are uniformly bounded, and the
center factors stay bounded above and away from zero by (3),
(8), and the finite degree set. A subsequence therefore gives
a common root of a compactified pair. Conversely every
compactified triple is a limit of ordinary profiles. For a
fixed \(B\), take the corresponding limit of the core
multiplicities, tail curvature and normalized data of those
profiles. If the curvatures approach \(K\) from above, apply the
estimates with \(K+o(1)\) in place of \(K\) and use the strict
inequality \(B>K\). Passing to the limit gives the original
closed estimates (6)–(10) and the exact allowed tail interval.
Thus a compactified common root makes every test feasible.
This proves the criterion. \(\square\)

The finite test is a sufficient certificate for ordinary
coprimality. A compactified common root would prevent this
certificate from succeeding. It would require a separate argument
to obtain an actual ordinary support with a common root.
