# Uniform coefficient phases for supports that vary with the order

Consider ordinary factorial moments with integer exponent rows
\(a_i,b_i\ge1\), after removing equal rows. Define
\[
 m_i=(a_i+b_i)/2,\quad d_i=b_i-a_i,\quad r_i=d_i/(2m_i),\quad
 w_i=d_i^2/m_i,\quad W=\sum_i w_i,\quad \Delta=\sum_i d_i .
\]
The profile and the dimension may depend on the terminal order.
Within each consecutive pair use the same profile at both orders.
For \(x_i(u)=a_i+d_i u\), use \([f]_{\rm c}(u)=f(u)-(1-u)f(0)-uf(1)\).
The fractional coefficient log, after removing the binomial
coefficient and normalizing both endpoints, is
\[
 L_n(u)=\sum_i[\log\Gamma(1+n x_i(u))]_{\rm c}.
\]
At the coefficient of index \(k\), take \(u=k/n\).
Define two analytic functions with zero endpoint values by
\[
 G(u)=\sum_i[x_i(u)\log x_i(u)]_{\rm c},\qquad
 H(u)=\tfrac12\sum_i[\log x_i(u)]_{\rm c}.
\]

**Theorem.** Fix \(W_0>0\) and \(Y<\infty\), and allow any
such profile with \(W\le W_0\). Uniformly on
\(0\le\Re u\le1,\ |\Im u|\le Y\), with constants depending
only on \(W_0,Y\),
\[
 L_n=nG+H+O(W/n),\qquad L_{n+1}-L_n=G+O(W/n^2).
 \tag{1}
\]
For the positive endpoint center
\[
 \mu_n=\frac1n\sum_i\log\frac{(n b_i)!}{(n a_i)!}
\]
the exact cross-order shift satisfies
\[
 \mu_{n+1}-\mu_n=\Delta\log(1+1/n)+O(W/n^2).       \tag{2}
\]
The possible limiting pairs \(G,H\) on such bounded profile
bands are determined by retained integer rows and a single
nonnegative diffuse mass. The mass adds a quadratic term to
\(G\) and adds zero to \(H\).

## Uniform estimates with an order derivative

Since \(w_i=2|d_i||r_i|\), a row with \(|r_i|\ge1/2\) has
\(|d_i|\le W_0\), so \(|r_i|\le W_0/(W_0+2)\) by
\(\min(a_i,b_i)\ge1\). Thus \(|r_i|\le R<1\) uniformly,
where \(R=\max(1/2,W_0/(W_0+2))\).
For \(0\le\Re u\le1\), this gives
\(\Re x_i(u)\ge(1-R)m_i\). We have \(m_i\ge1\).

For \(\Re x>0\), the trigamma series gives the exact formula
\[
 n^2\psi'(1+nx)=n\int_0^\infty e^{-xu}K(u/n)\,du,\qquad
 K(v)=\frac{v}{e^v-1}.
\]
Indeed, for \(k\ge1\) integrate \(t e^{-(nx+k)t}\) and sum,
then substitute \(u=nt\).
The real function \(K\) is analytic at zero with
\(K(0)=1,\ K'(0)=-1/2\). For \(v\ge0\),
\[
 |K(v)-1+v/2|\le C v^2,\qquad
 |K(v)-vK'(v)-1|\le C v^2 .
\]
On \([0,1]\) these follow by Taylor's theorem, and on
\([1,\infty)\) by the displayed formula and boundedness
of \(K,K'\). Consequently for all positive real \(n\ge1\),
\[
 \begin{aligned}
 n^2\psi'(1+nx)&=\frac n x-\frac1{2x^2}
             +O\!\left(\frac{1}{n(\Re x)^3}\right),\\
 \partial_n[n^2\psi'(1+nx)]&=\frac1x
             +O\!\left(\frac{1}{n^2(\Re x)^3}\right).
 \end{aligned}
\]
The derivative follows by dominated integration on any compact
interval of positive orders. These formulas bound the second
\(u\)-derivatives of the errors in (1), namely by \(C W/n\)
and \(C W/n^2\) after integrating the order derivative across
\([n,n+1]\), because \(d_i^2/m_i^3\le w_i\). Both errors have
zero endpoint values. A bound on the second derivative in the
rectangle bounds the function there: first use
\(f'(0)=-\int_0^1(1-v)f''(v)\,dv\), then integrate twice along
a straight segment from zero. This proves (1), with no error
from independently subtracting the two estimates of order
\(1/n\). [Trigamma series](https://dlmf.nist.gov/5.15.E1)

For (2), allow positive real \(n\) and write
\(\mu_n=\sum_i\int_{a_i}^{b_i}\psi(1+nx)\,dx\).
For real \(x>0\) integral comparison of the trigamma sum
gives \(1/(1+nx)\le\psi'(1+nx)\le1/(nx)\).
It follows that
\[
 |\partial_n\mu_n-\Delta/n|
 \le n^{-2}\sum_i|\log(b_i/a_i)|
 \le \frac{W}{(1-R)n^2}.
\]
The last inequality uses
\(|\log((1+r_i)/(1-r_i))|\le2|r_i|/(1-R)\) and
\(2|r_i|\le w_i\). Integrate over \([n,n+1]\).

## Retained measures and the diffuse term

The phase derivatives are the transforms
\[
 G''(u)=\int \frac{d\nu(r)}{1+(2u-1)r},\qquad
 H''(u)=-\frac12\int\frac{d\xi(r)}{(1+(2u-1)r)^2},
 \quad
 \nu=\sum_i w_i\delta_{r_i},\quad
 \xi=\sum_i\frac{w_i}{m_i}\delta_{r_i}.
\]
At each cutoff \(|r|\ge\epsilon>0\) only finitely many row
types can occur on the band. Indeed
\(1\le |d_i|\le W_0/(2\epsilon)\) and
\(m_i=|d_i|/(2|r_i|)\) is bounded. Each nonzero multiplicity
is bounded by the mass budget. A diagonal subsequence makes
every individual row-type multiplicity eventually constant and
makes the total masses converge to some \(W_*\le W_0\). The limiting
retained list has total mass at most \(W_*\). The difference
\(\tau\ge0\) becomes the diffuse atom \(\tau\delta_0\)
in the limit of \(\nu\), in addition to the retained measure.
No diffuse mass enters \(\xi\) since \(1/m_i\le2|r_i|\).
The estimates on the discarded terms are uniform on every
rectangle above, by the positive lower bound on the real
parts of the denominators. The transform for \(G''\) sees
the diffuse mass as the constant \(\tau\), and the corresponding
term in \(G\) is \(-\tau u(1-u)/2\).

Conversely, take any summable list of retained integer rows and
a diffuse mass \(\tau\) consistent with the band \(W\le W_0\).
Truncate the list and approximate \(\tau\) from below by many
copies of the increasing unit-gap row \((A,A+1)\), of individual
mass \(1/(A+1/2)\). As \(A\to\infty\), the masses sum to
\(\tau+o(1)\), their \(r\) tend to zero, and their contribution
to \(\xi\) tends to zero. Passing to longer truncations proves
attainability. This identifies the coefficient phase and its
degree increment on bounded-curvature sequences. A uniform
root comparison still requires control of the oscillatory
binomial sum with these phases.
