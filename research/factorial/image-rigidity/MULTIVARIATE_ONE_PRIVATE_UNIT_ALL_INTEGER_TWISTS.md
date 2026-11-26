# All-degree weak shared tails around a single private unit

Consider ordinary exponent columns with private variables on only
one side, and with total private degree one. Set
\[
                S_{\rm sh}=\sum_{\rm shared}
                   \frac{(a_i-b_i)^2}{\max(a_i,b_i)} .
\]

**Theorem.** There is a computable absolute positive bound such
that \(S_{\rm sh}\) below that bound implies physical SFC(2)
in every degree, for arbitrary signed total degree difference
and arbitrary mixed shared orientations. Specifically, if \(J\ge4\)
is any terminal cutoff for the
[absolute ordinary angular theorem](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md),
it suffices to have
\[
                             S_{\rm sh}< J^{-20J^4}.
\]
The value of \(J\) is unevaluated. No degree-by-degree moment
test is needed to establish the stated sufficient bound.

The exact core claim is that \(E_n(w)\) and
\(E_{n+1}(((n+1)/n)^h w)\) are coprime for every \(n\ge1\)
and \(h\in\mathbb Z\), where
\[
                           E_m(w)=\sum_{k=0}^m w^k/k!.
\]
To prove this, orient the private row as \((1,0)\); by itself
it gives the moment polynomial \(m!E_m(w)\). At \(h=0\) the
difference \(E_{n+1}-E_n=w^{n+1}/(n+1)!\) excludes any common
root. At \(n=1\) the lower root is real, while the upper
polynomial at any positive scaling is strictly positive there.

For \(n\ge2,\ h\ne0\), choose an odd prime \(p\) dividing
the odd member of \(n,N=n+1\). At this prime all roots of
every \(E_m\) have valuations in
\[
                              [0,1/(p-1)] .
\]
Indeed the coefficient valuations \(-v_p(k!)\) are
nonincreasing. All slopes of their lower convex hull are
therefore nonpositive. Also every coefficient point lies
on or above the line of slope \(-1/(p-1)\) through the
constant coefficient, since \(v_p(k!)=\sum_{j\ge1}
\lfloor k/p^j\rfloor\le k/(p-1)\).
The first Newton slope is at least this slope, and all
later ones are larger. The root valuations are the negatives
of the slopes. The dilation \((N/n)^h\) shifts all upper
root valuations by the same nonzero integer. Its shifted
interval is disjoint from the lower one since \(p\) is odd.
A complex common factor of the rational polynomials would
also give a common root over an algebraic closure of
\(\mathbb Q_p\), so this proves the core claim.

## An explicit tolerance at the bounded terminals

For a shared tail put \(\tau=S_{\rm sh}\), let \(\Delta\)
be its signed degree difference, and assume \(\tau\le1/10\).
After stripping its constant and normalizing its leading
factor by a positive change of variable, write the tail's
multipliers as \(t_{m,k}\), with \(t_{m,0}=t_{m,m}=1\).
Let \(D_m\) denote the tail's positive endpoint scale.
The
[endpoint log chord and radial drift estimates](MULTIVARIATE_ARBITRARY_CORE_TWISTED_TAIL_STABILITY.md#proof-of-stability)
give
\[
  0\le 1-t_{m,k}\le m\tau/4,\qquad
  D_N/D_n=(N/n)^\Delta e^\epsilon,\qquad |\epsilon|\le\tau.
                                                               \tag{1}
\]
The logs of \(t_{m,k}\) are convex in \(k\). Their continuous
interpolation comes from log Gamma with second derivative at
most \(\tau/[m(1-\tau)]\), using the same chord estimate.
Hence their consecutive differences have absolute value at most
\(\tau/[2(1-\tau)]\): the extremal differences occur at the ends,
where the quadratic chord bound gives that estimate.
Eneström–Kakeya applied to
\(\sum_k t_{m,k}w^k/k!\) therefore puts all its roots
in \(1/2\le |w|\le2m\).

Suppose a physical common root gives normalized variables
\(w=D_n z\) and \(w'=D_N z\).
For \(2\le n<N<J\), put \(\lambda=(N/n)^\Delta\).
The annuli and (1) imply
\[
 (8J)^{-1}\le\lambda\le8J,\qquad |\Delta|\le J^2 .
\]
For the second inequality use \(\log(N/n)>1/N\)
and \(\log(8J)<J\) for \(J\ge4\).
The weighted coefficient bound at \(w,w'\), followed by integration
of \(wE_N'(w)\) along the radial segment from \(w'\) to
\(\lambda w\), gives
\[
                |E_n(w)|,\ |E_N(\lambda w)|\le 5J e^{4J}\tau .
\]
Indeed the two coefficient errors are at most
\(J e^{2J}\tau/4\), and on the segment \(|w|\le4J\)
the derivative norm is at most \(4J e^{4J}\).

Write \(\lambda=u/v\) with positive coprime integers.
Both \(u\) and \(v\) are at most \(J^{J^2}\).
The integer polynomials \(A=n!E_n\) and
\(B=v^N N!E_N(\lambda\,\cdot)\) have coefficients bounded
by \(M=J^{2J^3}\), and their two values at \(w\) have absolute
value at most \(M\cdot5Je^{4J}\tau\).
Their Sylvester determinant is a nonzero integer by the
core claim. The adjugate identity gives integer polynomials
\(U,V_0\) of degrees at most \(2J-1\) satisfying
\(UA+V_0B=\pm\operatorname{Res}(A,B)\), with coefficients
bounded by \((2J)!\,M^{2J}\). In fact every coefficient
is a signed determinant of size less than \(2J\).
Since \(|w|\le2J\), evaluating gives
\[
    1\le 20J^2 e^{4J}(2J)^{2J}(2J)!\,M^{2J+1}\tau
       \le J^{20J^4}\tau .
\]
For the last inequality, \(e^4\le J^3\) for \(J\ge4\),
\((2J)!\le(2J)^{2J}\), \(2J\le J^2\), and \(20J^2\le J^5\).
Then the power of \(J\) is at most
\(4J^4+2J^3+11J+5\le20J^4\).
This contradicts the claimed bound on \(\tau\).
At terminal \(N\ge J\) use the absolute angular theorem,
and \(N=2\) follows from even-moment positivity.
Rational Bézout descent gives the result over every
characteristic-zero field.
