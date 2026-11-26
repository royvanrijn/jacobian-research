# Quadratic trace certificates from adjacent Newton polygons

Consider an ordinary two-monomial support with a private variable on
each side. The private exponents are arbitrary. The following
criterion resolves a possible common factor of degree two by
comparing its arithmetic trace with the actual angular boxes.
The Newton polygons may have any number of edges.

At an adjacent pair \(n,N=n+1\), write \(F_n,F_N\) for the moment
polynomials after dividing out any common coefficient scalars so
that both have integer coefficients. Let \(S_p\) be the set of
shared root valuations of the two Newton polygons at the prime
\(p\), with \(v_p(p)=1\). Thus the numbers in \(S_p\) are
negatives of shared slopes, with no change of variable between
the two polynomials. For a shared valuation \(s\), let
\(w_{m,p}(s)\) denote the horizontal width of the corresponding
edge in degree \(m\). Assume that the rational common factor degree
has been bounded by two. For example, it suffices that at some
prime
\[
           \sum_{s\in S_p}\min(w_{n,p}(s),w_{N,p}(s))\le2.      \tag{1}
\]

**Theorem.** Suppose that \(S_p\cap\mathbb Z\) contains at most one
value at every prime. The adjacent pair is then coprime unless
at least one of the following two conditions holds:
\[
 \begin{array}{ll}
 n\equiv5\pmod 8,& S_2\text{ contains a half-integer outside }\mathbb Z,\\
 n\equiv 7,8,\text{ or }9\pmod {12},&
                  S_3\text{ contains a half-integer outside }\mathbb Z.
 \end{array}                                                   \tag{2}
\]
Each condition requires both entries on its line. The hypothesis
on integral shared valuations needs checking only at primes dividing
the endpoint products of both polynomials. A wider trace sieve,
given below, can also use primes that have more than one integral
shared valuation.

If these certificates or other exact separations cover
\(3\le N<J\) for an original profile, the
[computable absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
handles \(N\ge J\). The first and second moments are always coprime
by strict positivity of the quadratic on the real line. The
argument in this note establishes physical coprimality.

## The possible quadratic and its angle

A common rational factor has even degree because the even
moment polynomial is strictly positive over the real line.
Under the degree bound a common factor, if present, is therefore
an irreducible quadratic with nonreal roots \(\alpha,\bar\alpha\).
Choose the one in the upper half-plane. Put
\[
   u=\alpha+\bar\alpha,\qquad c=\alpha\bar\alpha,\qquad
              t=\frac{\alpha}{\bar\alpha}
                      +\frac{\bar\alpha}{\alpha}
                 =\frac{u^2}{c}-2 .
\]
The
[both-private angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
give \(0<\theta=\arg\alpha<\pi/2\), and both roots at the two
orders have the same positive label \(k\) after choosing the
upper root. In particular \(u>0,c>0\), and
\[
       \frac{(2k-1)\pi}{n}<\theta<\frac{2k\pi}{n+2},
       \qquad 4k<n+2,\qquad t=2\cos(2\theta)\in(-2,2).          \tag{3}
\]

At each prime both quadratic roots have valuations in \(S_p\).
If their valuations differ, they are both integers. Indeed,
a quadratic irreducible over \(\mathbb Q_p\) would have roots
conjugate over \(\mathbb Q_p\) with the same valuation, by
uniqueness of the extended valuation. Distinct values therefore
force splitting over \(\mathbb Q_p\). Consequently, under the
integral-valuation hypothesis of the theorem, the roots have
equal valuation at every prime. The quotient
\(\alpha/\bar\alpha\) and its inverse are then unit everywhere,
and \(t\in\mathbb Z\).

Only \(t=-1,0,1\) lie strictly between \(-2\) and \(2\).
These give \(\theta=\pi/3,\pi/4,\pi/6\), respectively.
For \(\theta=\pi/3\), (3) would require
\(6k-3<n<6k-2\), impossible for integer \(n\).
The two other possibilities give precisely the congruences
in (2). If \(\theta=\pi/4\), the equality
\(c=u^2/(t+2)=u^2/2\) gives \(v_2(c)\) odd.
Since the root valuations agree at \(2\), they must be
half-integers outside \(\mathbb Z\).
For \(\theta=\pi/6\) the same proof uses \(c=u^2/3\)
and \(p=3\). This proves (2).
If one of the integer polynomials has unit endpoints at a
prime, all its roots have valuation zero there. This explains
the finite restriction on primes.

## An exact sieve with more than one integral shared slope

For the finitely many primes dividing both endpoint products, put
\[
 \delta_p=\max(S_p\cap\mathbb Z)-\min(S_p\cap\mathbb Z)
\]
if there are at least two such values, and \(\delta_p=0\)
otherwise. Let \(D=\prod_p p^{\delta_p}\). At any prime where
the root valuations differ, \(t\) has valuation exactly the
negative of their difference in absolute value. This is
because the two quotient summands in \(t\) have distinct
valuations there. At equal root valuations, \(t\) is integral.
It follows that \(D t\) is an integer.

Thus \(t\) can be restricted to the finite rational grid
\(D^{-1}\mathbb Z\) in the intervals
\[
     2\cos\!\left(\frac{4k\pi}{n+2}\right)
              <t<
     2\cos\!\left(\frac{2(2k-1)\pi}{n}\right),
                       \qquad 4k<n+2.                         \tag{4}
\]
The endpoints are real algebraic numbers, so these are
exact interval tests. At each remaining candidate, let
\(v_p(0)=+\infty\) if \(t=0\).
If \(v_p(t)<0\), any possible root slopes must be distinct
integers in \(S_p\) with difference in absolute value
\(-v_p(t)\). Then \(v_p(u)\) is the smaller root valuation,
since the two summands in \(u\) have different valuations.
If \(v_p(t)\ge0\), the two root valuations must agree at some
\(s\in S_p\), and
\[
                          v_p(u)=s+\tfrac12 v_p(t+2)           \tag{5}
\]
must be integral. These tests include every prime dividing a
numerator or denominator of \(t+2\) as needed in (5).
Outside the union of the joint endpoint support and the
prime divisors of the numerator and denominator of \(t+2\),
the trace valuation must be zero. In particular the choices
at the relevant primes give finitely many positive rational
values of \(u\). Each determines the only possible monic
quadratic,
\[
                          z^2-uz+\frac{u^2}{t+2}.
\]
The candidates can be tested by exact division or discarded
by additional necessary conditions. Thus the sieve confines
a possible degree-two gcd using the Newton polygons and angles
without computing either moment's complex roots.

## Native sources of the degree bound

The bound (1) is immediate from the valuation of a common
root and the root counts given by the edge widths. In
particular it does not require single-edge polygons.

A simple source of bound two is a native delayed mark.
At \(p>N\), strip the \((N\ell)!\) row scalars from the
upper polynomial. Suppose one increasing row has a unique
divisible factor \(N\ell+h\), first sampled at index
\(j=\lceil h/d\rceil\), and every other sampled factor
is a unit. At indices below \(j\) the upper valuations
are zero. At later indices they equal the positive
valuation \(E\) of the marked factor. The two upper
edges have widths \(j-1\) and \(W=N-j+1\), with slopes
zero and \(E/W\). If the lower roots are all units and
\(j=3\), (1) holds. A decreasing mark is handled by
reversing both polynomials.

For \(j=3,4\) one can instead proceed without a condition
on the lower polynomial at this prime. If \(\gcd(E,W)=1\),
one rational irreducible factor of the upper polynomial
must consume the entire edge of width \(W\), since
any nonzero width of a factor at that slope is a multiple
of its reduced denominator \(W\). The angular boxes bound
the common factor degree by \(2\lfloor N/4\rfloor\) on the
entire both-private class, since each common conjugate
pair has one of the labels \(k\) with \(4k\le N\).
Consequently \(W>2\lfloor N/4\rfloor\) rules out that
irreducible factor in the gcd. Every remaining common
root is a \(p\)-adic unit, and the common factor degree is
at most \(j-1\). For \(j=3,4\) its degree is therefore
at most two by parity.

The
[opposite marked-row test](MULTIVARIATE_OPPOSITE_MARKED_ROW_PRIME.md)
provides another source: a narrow interval of upper unit
coefficients and a lower polynomial with unit endpoints.
Allowing the upper interval to have width two, rather than
width at most one, gives (1), and the trace tests above can
then exclude the remaining quadratic. These sources use
more than one upper Newton edge.

## Arithmetic of the two-unit-root factor

For a delayed mark at \(j=3\), write the upper polynomial as
\(F_N(z)=A(z)+p^E H(z)\), where \(A(z)=a_0+a_1z+a_2z^2\)
and its endpoints are \(p\)-adic units. The reciprocal
polynomial reduces to \(z^{N-2}(a_2+a_1z+a_0z^2)\), with
coprime factors. Hensel factorization, followed by reversing
the degree-two factor and making it monic, gives the unique
monic polynomial \(g\) of the two unit roots.
Let \(g_0=A/a_2\). The same factorization modulo \(p^E\)
gives \(g\equiv g_0\bmod p^E\). In fact
\[
              g\equiv g_0+
                \frac{p^E}{a_2}\operatorname{rem}_{g_0}H
                                             \pmod {p^{2E}}. \tag{6}
\]
To check (6), write \(g=g_0+p^E g_1\), where
\(\deg g_1<2\). Division by \(g\) gives remainder
\(p^E(\operatorname{rem}_g H-a_2g_1)=0\).
Remainders by the two monic quadratics agree modulo \(p^E\),
so the formula follows. No assumption about a repeated
root in the quadratic reduction is needed.

If all possible common roots are \(p\)-adic units,
their common quadratic must be \(g\). In particular its
trace parameter satisfies
\[
                            t\equiv
                        \frac{a_1^2}{a_0a_2}-2
                                                           \pmod p.
\]
Formula (6) gives the next precision. These can be combined
with the integer trace restriction or the rational grid
without factoring either moment over \(\mathbb Q\).
