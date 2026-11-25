# A variable-prime Eisenstein criterion for a marked exponent row

Fix a terminal \(N\ge2\). For an ordinary exponent profile
\(a,b\) let
\[
 M_N(t)=\sum_{k=0}^N\binom Nk
          \prod_i((N-k)a_i+kb_i)!\,t^k .
\]
Select one row with exponents \(0\le A<B\), and orient the columns
so this row increases. Let \(B_0\) be the maximum exponent among
all other rows, with \(B_0=0\) if there are none.

**Theorem.** If there exists a prime \(p\) with
\[
 \max\{N,NA,NB/2,NB_0\}<p\le (N-1)A+B,
 \tag{1}
\]
then \(M_N\) is irreducible over \(\mathbb Q\). In particular the
order-\((N-1)\) and order-\(N\) moments are coprime over every
characteristic-zero field.

At index zero the selected row contributes \((NA)!\), a unit
modulo \(p\). At every positive index its factorial argument
lies in \([p,2p)\), by (1), so it contributes exactly one power
of \(p\). Every other row has factorial arguments below \(p\),
and all binomial coefficients are units because \(p>N\).
The reversed polynomial \(t^N M_N(1/t)\) is Eisenstein at \(p\).
The degree argument gives coprimality over \(\mathbb Q\), and
rational Bézout transfers it to other characteristic-zero fields.

The same proof works if the condition \(p>NB_0\) is replaced by
the condition that the sum of the other rows' coefficient
valuations is affine in \(k\). If this sum is \(u+vk\), first
multiply the variable by \(p^{-v}\) and divide the polynomial
by \(p^u\). The resulting coefficients are integral over the
local ring at \(p\) and satisfy the displayed Eisenstein test.
Irreducibility over \(\mathbb Q_p\) implies irreducibility over
\(\mathbb Q\).

## Arbitrary subordinate rows and arbitrarily small marked curvature

For any prime \(p>N^2\), set
\[
 A=\lfloor p/N\rfloor,\qquad r=p-NA,\qquad B=A+r.
\]
Here \(1\le r\le N-1\), because \(p>N\).
Then \((N-1)A+B=p\) and \(NB<2p\), since
\(A\ge N>r\). Any subordinate profile with exponents at most \(A\)
satisfies (1). The new marked curvature is \(r^2/(A+r)\), which
can be arbitrarily small.

In particular, every fixed ordinary profile can be extended by
two fully shared rows \((A,A+r),(A,A-r)\) so that its terminal
pair is coprime, with its original private degrees and signed
total-degree difference preserved. The second row satisfies the
subordinate condition. The added total curvature is at most
\(2(N-1)^2/A\) and tends to zero over the primes. For instance,
with \(N=6\) the marked gap is always at most five, while the
remaining shared and private rows are arbitrary up to the bound
\(A\).

The pair of centered coefficient arrays and the ratio of the
two positive normalizing centers also converge to the original
ones under this balanced extension. To see this, write \(n=N-1,N\)
and let \(\mu_n\) be the logarithm of the order-\(n\) endpoint
center. The two new rows add
\[
 \frac{\log((n(A+r))!)+\log((n(A-r))!)-2\log((nA)!)}{n}
\]
to \(\mu_n\). This is \(O_N(1/A)\), by bounding the second
derivative of \(\log\Gamma(1+x)\) by \(1/x\) on the indicated
interval. The same bound gives \(O_N(1/A)\) for every added
normalized coefficient log, by comparison with the quadratic
having second derivative \(1/[n(A-r)]\) along its factorial
arguments. Both estimates are uniform in the bounded \(r\).

This criterion selects its prime from a single row and permits
arbitrary subordinate profiles. The simultaneous finite-prime
obstruction in the [digit analysis](MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md#a-simultaneous-obstruction-at-finitely-many-primes)
concerns certificates with primes fixed in advance.
