# Weak shared tails and integer drift for two private unit rows

Consider an ordinary binomial with one private unit row on each
side and arbitrary shared rows. Set
\[
          V=\sum_{\rm shared}\frac{(b_i-a_i)^2}{\min(a_i,b_i)},
          \qquad
          \Delta=\sum_{\rm shared}(b_i-a_i).
\]
At any reference adjacent pair, the
[exact matched tanh grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)
gives \(V=\sum_i2v_i/(1-v_i)\) over shared Gauss factors, with
\(v_i=(n/T_i'-(n-1)/T_i)^{-1}\). The factor orientations give
the integer \(\Delta\).

**Theorem.** At a consecutive pair with terminal \(N\), the
ordinary physical moments are coprime if \(NV\le1\).
In particular, for \(J\ge4\) at least the
[computable absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md),
the condition \(V\le1/J\) gives SFC(2) in every degree.
There is no restriction on the number, orientations, or signed
total gap of the shared rows. The result concerns physical
coprimality. Its radial estimate uses the fact that the
curvature of the coefficient skew is quadratic in \(V\).

## Normalization and a bound for the coefficient skew

After the two private factorials cancel the binomial factor,
normalize the degree-\(m\) polynomial to have endpoint
coefficients one. Write it as
\[
                     P_m(z)=\sum_{k=0}^m c_k z^k .
\]
For real index \(x\in[0,m]\), write \(c(x)=e^{f(x)}\)
for the continuous Gamma interpolation, with its affine
endpoint log subtracted. Then \(f(0)=f(m)=0\) and
\(0\le f''(x)\le V/m\). In particular \(f\le0\),
\(f\ge-mV/8\), and \(|f'|\le V/2\).

Put \(g(x)=f(m-x)\). For a row with signed gap \(d\)
and smaller exponent \(\ell\), the corresponding
trigamma arguments at \(x,m-x\) differ by at most
\(m|d|\), and both exceed \(m\ell\). The estimate
\[
           |\psi''(1+y)|
                =2\sum_{h=1}^{\infty}(y+h)^{-3}
                       \le y^{-2}\qquad(y>0)
\]
follows by integral comparison. Thus the contribution
to \(|f''-g''|\) is at most \(|d|^3/(m\ell^2)\).
Since \(|d|\ge1\), we have \(|d|/\ell\le d^2/\ell\),
and summing gives
\[
                       |f''-g''|\le V^2/m .
\]
The zero endpoints of \(f-g\) then imply
\(|f'-g'|\le V^2/2\) and \(|f-g|\le mV^2/8\).
For the coefficient skew \(\delta(x)=e^{f(x)}-e^{g(x)}\),
use the Lipschitz bound for the exponential on the
nonpositive half-line to get
\[
        |\delta'(x)|\le V^2/2+mV^3/16 .
\]
When \(mV\le1\), its sampled total variation satisfies
\[
                \sum_{k=1}^m |\delta(k)-\delta(k-1)|
                                  \le 9mV^2/16.             \tag{1}
\]

## A quadratic bound on the radial displacement

Average the coefficients of \(P_m\) and its reciprocal
to obtain \(\widehat P_m\). These averaged coefficients
are the moments of a positive random variable.
Indeed, the original normalized coefficients have the
form \(c_k=\mathbb E X^k\), \(c_m=1\). Weighting the law
of \(X\) by \(X^m\) and using \(1/X\) gives moments
\(c_{m-k}\). Mix these two laws in equal proportions.
The averaged polynomial is reciprocal, and its roots
lie on the unit circle: the
[private-unit angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
give one simple root per argument box for such moment
coefficients, and reciprocal conjugation preserves
the argument. In the geometric case the assertion
follows directly from the roots of unity.
Distinct reference arguments have circular distance
at least \(\pi/m\), and each has distance at least
\(\pi/m\) from the positive ray. This follows by comparing
successive box endpoints, with conjugation and the
negative root in odd degree.

Assume \(m\ge4\), \(0<V\), and \(mV\le1\).
Write
\[
                  (1-z)\widehat P_m(z)
                         =1-z^{m+1}+E(z).
\]
The averaged coefficients are log-convex, by their
moment representation, and lie between \(e^{-mV/8}\)
and \(1\). Thus their total variation is at most
\(2(1-e^{-mV/8})\le mV/4\). The coefficients of \(E\)
are these first differences, so
\[
                 |E'(z)|\le e^{1/4} m^2 V/4 <m/3
                       \qquad (|z|\le1+1/(4m)).
\]
Let \(\zeta\) be any reference root and put \(\eta=V^2\).
Then \(\eta\le1/(4m)\). On \(|z-\zeta|=\eta\), binomial
expansion gives
\[
              |z^{m+1}-\zeta^{m+1}|
                              \ge \tfrac34(m+1)\eta .
\]
Indeed \(\alpha=(m+1)\eta\le3/8\), and the sum of binomial
terms of order at least two is at most
\(\alpha^2/[2(1-\alpha/3)]\le\alpha/4\).
Subtracting the value at \(\zeta\) in the expression
for \((1-z)\widehat P_m(z)\) now gives a modulus of at
least \(5m\eta/12\).

On the other hand, \((1-z)(P_m-\widehat P_m)\) has
coefficients equal to half the first differences in
(1), including zero skew at the endpoints. Its modulus
on this circle is at most
\[
                      \frac{e^{1/4}}{2}\,\frac{9mV^2}{16}
                                   <\frac38 mV^2 .
\]
Rouché's theorem thus gives one actual root in each
reference disk. These disks are disjoint and avoid
\(1\) by the reference angle spacing. Since every
reference center has modulus one, every actual root
satisfies
\[
                                 |\log|z||\le 2V^2.           \tag{2}
\]
When \(V=0\) the polynomial is already geometric
and (2) is immediate.

## Integer drift at a physical common root

Let \(g_m\) be the positive scale carrying a physical
root in degree \(m\) to its endpoint-normalized root.
The
[linked center estimate](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md#an-explicit-linked-radial-condition)
gives
\[
       \left|\log(g_N/g_n)-\Delta\log(N/n)\right|
                                    \le V/(nN).
\]
At a common physical root with \(N\ge5\) and \(NV\le1\),
equation (2) applies at both orders. Hence
\[
             |\Delta|\log(N/n)\le4V^2+V/(nN).
\]
If \(\Delta\ne0\), using \(\log(N/n)>1/N\) would give
\[
               1\le|\Delta|<4NV^2+V/n
                                  \le 4/N+1/(Nn)<1.
\]
Thus any common root forces \(\Delta=0\). The
[positive-radius gate](MULTIVARIATE_PRIVATE_UNIT_LINKED_POSITIVE_RADIUS_GATE.md)
excludes a common root at \(\Delta=0\) when
\(V\le6/(n+2)\), which follows from \(NV\le1\).
The orders with \(N\le4\) are covered by the
[ordinary low-order theorem](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md).
This proves the pairwise result. When \(V\le1/J\),
all \(N<J\) satisfy it, and \(N\ge J\) follow from the
absolute angular theorem.
