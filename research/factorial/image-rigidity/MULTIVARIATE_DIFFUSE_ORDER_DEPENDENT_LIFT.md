# An order-dependent lift from diffuse to factorial moments

Let \(a,b\) be distinct nonnegative integer exponent columns.
After removing equal rows, put
\[
 S=\sum_i\frac{(b_i-a_i)^2}{max(a_i,b_i)},\qquad
 \eta=\max_i\frac{|b_i-a_i|}{max(a_i,b_i)} .
\]

**Theorem.** There is an integer \(N_0\) such that for every
terminal \(N\ge N_0\), the ordinary factorial moment
polynomials of orders \(N-1,N\) are coprime whenever
\[
                         \eta\le e^{-20N^6}.
\]
This holds in arbitrary dimension, for all \(S\) and all
total-degree differences, over every characteristic-zero field.
The relative-gap bound is explicit. The order threshold is
the unevaluated threshold of the
[uniform eventual diffuse theorem](MULTIVARIATE_DIFFUSE_VANISHING_COUPLING.md),
increased if necessary so that \(N\ge4\).

The proof bounds the minimum of the integer diffuse resultant
on a rational interval using its critical values. It then
compares the two Sylvester determinants at the actual ratio
of the positive endpoint centers.

## A quantitative minimum for the diffuse resultant

Write \(n=N-1\) and set
\[
 P(z;u)=\sum_{k=0}^n\binom nk u^{N k(n-k)}z^k,\qquad
 Q(z;u)=\sum_{k=0}^N\binom Nk u^{n k(N-k)}z^k,qquad
 R(u)=\operatorname{Res}_z(P(z;u),Q(z;u)).
\]
For \(0<u<1\), these are the diffuse polynomials with
\(s=-2nN\log u>0\). For all terminals under consideration
they are coprime there. At \(u=0\) the two polynomials
are \(1+z^n\), \(1+z^N\); their resultant is \(2\), by
evaluating \(1+z^N=1-z\) at a root of \(1+z^n\).
At \(u=1\) the resultant vanishes. Thus \(R\) is an integer
polynomial of degree \(D\ge1\), with no zero in \([0,1)\).

Each entry of the Sylvester matrix has degree at most
\(N^3/4\) in \(u\), hence \(D\le N^4\). On \(|u|=1\), the
coefficient row of \(P\) has Euclidean norm at most \(2^n\)
and that of \(Q\) at most \(2^N\), by the sum of the
binomial coefficients. Hadamard's inequality and Cauchy's
coefficient bound give
\[
 |R(u)|\le 2^{2nN}\quad (|u|=1),\qquad
 \max_j |[u^j]R|\le H=2^{2N^2}.
\]
Let \(C(y)=\operatorname{Res}_u(R(u)-y,R'(u))\).
This is a nonzero integer polynomial, since for a generic
\(y\) it has no common vanishing with the fixed derivative.
Its degree in \(y\) is at most \(D-1\). For \(|y|=1\), every
Sylvester row has norm at most \(2(D+1)^2H\) by summing
its coefficient magnitudes. The matrix size is \(2D-1\).
A second application of Hadamard and Cauchy's bound gives
height at most \(e^{8N^6}\), using
\(\log(D+1)\le N^2\) and
\(\log 2+2N^2\log 2\le2N^2\).
After removing the maximal power of \(y\), its constant
term has modulus at least one. Any nonzero root \(y\)
of modulus below one therefore satisfies
\[
 1\le D e^{8N^6}|y|,qquad |y|\ge e^{-9N^6}.
\]
This bounds the modulus of every nonzero critical value
of \(R\).

Put \(T=10^{17}N^8\). The endpoint \(1-1/T\) is rational,
and the nonzero value there has modulus at least \(T^{-D}\),
by clearing denominators. Since \(\log T\le4N^2\) for
\(N\ge4\), this is at least \(e^{-4N^6}\). The value at
zero is \(2\). An interior minimum of \(|R|\) on this
interval is a nonzero critical value in modulus, because
\(R\) is real with constant sign there. It follows that
\[
              |R(u)|\ge e^{-9N^6}\qquad
              (0\le u\le1-1/T).                            \tag{1}
\]

## The actual normalized pair

The existing
[weak-curvature bound](MULTIVARIATE_WEAK_CURVATURE_ALL_TERMINALS.md)
and
[large-curvature bound](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#large-curvature-and-the-endpoint-radii)
cover \(S\le10^{-16}N^{-6}\) and \(S\ge4\log N+10\),
respectively. For the remaining band, take
\(u=e^{-S/(2nN)}\). The inequality \(1-e^{-x}\ge x/(1+x)\)
gives \(1-u\ge1/T\) since \(S/(2nN)>1/(2\cdot10^{16}N^8)\).
Thus (1) applies throughout the band.

At order \(j=n,N\), let \(\widehat P_j(z)\) be the actual
polynomial with its endpoints normalized to one and its
variable rescaled by the positive center
\[
 G_j=\left(\frac{\prod_i(jb_i)!}{\prod_i(ja_i)!}\right)^{1/j},
 \qquad z=G_jt .
\]
The order-\(j\) limiting polynomial is \(P\) or \(Q\),
evaluated at the fixed \(u\) above. For each unequal row let
\(M_i=\max(a_i,b_i)\), \(\ell_i=\min(a_i,b_i)\) and
\(d_i=b_i-a_i\). Since \(\eta\le1/2\), all rows are shared,
\(M_i\ge1/\eta\), and \(\ell_i\ge(1-\eta)M_i\).
For \(g(x)=\log\Gamma(1+x)\), integral comparison of the
trigamma series gives \(1/(x+1)\le g''(x)\le1/x\).
On the row's interpolation segment
\(j\ell_i\le x=ja_i+kd_i\le jM_i\), this gives
\[
 1-\eta/j\le jM_i g''(x)\le1+2\eta .
\]
The second derivative in \(k\) of the coefficient logarithm
after removal of the binomial, summed over rows, differs
from \(S/j\) by at most \(2\eta S/j\).
Subtract the endpoint chord, as in the normalization.
The error relative to the diffuse logarithm has zero
endpoints and absolute value at most
\(\eta S k(j-k)/j\le\eta S j/4\), by comparison with
the quadratic having that second-derivative bound.
This last quantity is below \(1/2\) under the theorem's
hypothesis. The diffuse coefficients are at most the
binomial coefficients, so the coefficient errors satisfy
\[
 \bigl|[z^k](\widehat P_j(z)-B_j(z;S))\bigr|
            \le E\binom jk,qquad E=7N^2\eta .
 \tag{2}
\]
Here \(|e^v-1|\le2|v|\) for \(|v|\le1/2\), and
\(S<4\log N+10\le7N\).

Each actual normalized root has \(1/2<|z|<2\).
Indeed all diffuse roots are on the unit circle by the
[circle theorem](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#the-limiting-pair-when-each-relative-gap-tends-to-zero).
For \(|z|\ge2\), its degree-\(j\) monic polynomial
has modulus at least \((|z|-1)^j\), whereas (2) bounds
the error by \(E(1+|z|)^j\le E3^j(|z|-1)^j\) with
\(E3^N<1\). Apply the same argument after coefficient
reversal, since the diffuse polynomial is reciprocal
and actual endpoints are one. A common original root
therefore gives the relation \(z_N=c z_n\) between its
normalized variables, with \(c=G_N/G_n\in[1/4,4]\).

## Comparison at the actual center ratio

Write \(z_j,w_k\) for the diffuse roots at orders \(n,N\).
Their unit moduli give
\[
 |c z_j-w_k|^2=(c-1)^2+c|z_j-w_k|^2.
\]
Since \(P\) is monic, its resultant with \(Q(cz;u)\)
is the product of these pairwise factors in modulus.
Consequently
\[
 \left|\operatorname{Res}_z(P(z;u),Q(cz;u))\right|
    \ge c^{nN/2}|R(u)|\ge e^{-10N^6},\qquad c\in[1/4,4].
 \tag{3}
\]
No closeness of the two endpoint centers is needed.

The corresponding actual resultant is
\(\operatorname{Res}_z(\widehat P_n(z),\widehat P_N(cz))\).
Each diffuse and actual Sylvester row has norm at most
\(2\cdot5^N\), and the norm of its perturbation is
at most \(E5^N\), by (2) and
\(\sum_k\binom Nk c^k\le5^N\). Replace the rows one
at a time and use multilinearity and Hadamard's inequality.
The absolute difference of the two determinants is at most
\[
        2N E (2\cdot5^N)^{2N}\le E e^{6N^2}.
\]
For \(\eta\le e^{-20N^6}\) this is strictly below the lower
bound (3), so the actual resultant is nonzero at the
center ratio of any hypothetical common root. This
excludes common roots over \(\mathbb C\). The original
polynomials have integer coefficients, giving rational
Bézout identities which extend to every characteristic-zero
field. \(\square\)
