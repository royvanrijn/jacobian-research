# Two-prime coprimality for same-orientation unit gaps

Suppose the distinct ordinary monomials differ by one in each
unequal row and every unequal row points in the same direction.
There can be any finite number \(R\ge2\) of these rows, with
arbitrary fixed smaller exponents. Then for a set of even orders
of relative lower density at least \(\log 2\), the actual factorial
moment polynomial is irreducible over \(\mathbb Q\). The same
arithmetic certificate gives consecutive coprimality on a set
of terminal orders of lower density at least \(\log 2\).

This extends the two-prime mechanism of
[Filaseta and Rocks, Theorem 1.1 and Lemmas 2.3–2.4](https://people.math.sc.edu/filaseta/papers/ColloquiumMathPaper2016.pdf).
Their single-row constant-plus-power polygons are replaced
by the two polygons of the actual multivariate coefficients.
The argument does not use the general eventual angular cutoff.

## An explicit two-prime criterion

Remove the equal rows and orient all the unit gaps to increase.
Write \(A_i\ge0\) for their smaller exponents and \(A=\max_i A_i\).
At order \(m\) put
\[
                M_m(z)=\sum_{k=0}^m
                 \binom{m}{k}\prod_{i=1}^R(m A_i+k)!\,z^k .
\]
Suppose \(n=Kq\), where \(q\) is prime and
\[
                              q>\max\{R,(A+1)K\}.
 \tag{1}
\]
Take a second prime \(p\le n\), and let
\[
                  t=n-p,\qquad w=p-(A+1)t-1 .
\]
It suffices that
\[
      (A+1)t<p,\qquad p^2>(A+1)n,\qquad 2w>p,\qquad n-w<q-1.
 \tag{2}
\]
Under (1)–(2), \(M_n\) either is irreducible or is a product
of a linear polynomial and an irreducible polynomial of degree
\(n-1\). It is irreducible when \(n\) is even.
The conditions also certify coprimality at terminal \(n\) when
\(n\) is odd.

Use the constant-to-leading orientation for Newton polygons,
subtracting any common scalar valuation. The polygon product
rule of [Filaseta and Rocks, Theorem 2.1](https://people.math.sc.edu/filaseta/papers/ColloquiumMathPaper2016.pdf)
says that a rational factor uses widths from the same slopes
as the full polynomial. At a slope with reduced denominator
\(b\), the width contributed by a rational factor is a multiple
of \(b\), since its edge endpoints have integral coordinates.

### The prime dividing the order

By (1) every factorial argument is below \(q^2\).
At \(k=hq+s\), \(0\le s<q\), the change in \(q\)-valuation of
each row from its constant coefficient is \(h\).
The binomial valuation is \(\mathbf1_{s>0}\), by the factorial
floor formula with \(n=Kq<q^2\).
The ordinates are therefore \(Rh+\mathbf1_{s>0}\).

The first edge runs from \(0\) to \(q-1\) with slope
\(1/(q-1)\). Over the next complete blocks the last point of
each block lies on the line of slope \(R/q\) through \(q-1\)
and \(n-1\); the interior points and the starts of each block
lie above that line. For the starts, its value at \(hq\) is
\(Rh+1-R+R/q<Rh\). The final slope at width one is \(R-1\).
The middle edge is absent for \(K=1\).
The displayed slopes are strictly ordered, and \(R/q\) is
reduced. Thus any positive-degree rational factor of degree
strictly less than \(q-1\) is linear.

### The nearby prime

Again there are no contributions from levels above \(p\) by
(2). At index \(k\), the row \(i\) changes valuation once,
at \(k=p-A_i t\): its quotient in the factorial floor is
\(A_i+\lfloor(A_i t+k)/p\rfloor\), and \(A_i t+k<2p\).
The binomial valuation is \(\mathbf1_{t<k<p}\).
The full ordinate is therefore zero on \(0\le k\le t\),
one on \(t<k\le p-A t-1\), and at least two to its right.
At indices from \(p\) onwards this last assertion uses
\(R\ge2\). As \(2w>p\), every slope from the point at
\(t\) to an index beyond \(p-A t-1\) is strictly larger
than \(1/w\). Hence the first positive edge has width \(w\)
and slope \(1/w\).

Some rational irreducible factor uses this edge, which has
width equal to its reduced denominator. That factor has
degree at least \(w\). Its complement has total degree at
most \(n-w<q-1\), so the conclusion from \(q\) leaves only
degrees zero and one for the complement. At even \(n\)
the independent exponential representation gives
\(M_n(z)>0\) for every real \(z\): the integrand is a positive
even power except on a null set. Thus there is no linear
factor and \(M_n\) is irreducible.

Irreducibility certifies coprimality with the preceding
moment by degree, and with the following moment by
[nondivisibility](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#an-inequality-using-the-endpoint-coefficients).
For odd \(n\), the even preceding moment cannot share the
possible linear factor. Sharing the large factor of degree
\(n-1\) would make the preceding moment divide \(M_n\),
also excluded by nondivisibility. Rational Bézout transfers
the coprimality assertions to every characteristic-zero field.

## Density

Fix a real \(C>A+4\).
For all sufficiently large \(n\), conditions (1)–(2) hold
if \(n\) has a prime factor \(q>C\sqrt n\) and there is a
prime \(p\) in \((n-\sqrt n,n]\). Indeed \(q^2>C^2n\)
gives \(q>(A+1)K\), and
\(n-w=(A+2)t+1<q-1\) for large \(n\).
The remaining inequalities in (1)–(2) follow with \(A,R\)
fixed.

The argument for
[Filaseta and Rocks, equations (3.2)–(3.3)](https://people.math.sc.edu/filaseta/papers/ColloquiumMathPaper2016.pdf)
shows that the first prime condition has asymptotic density
\(\log 2\), for any fixed \(C>1\), and the second has density
one. The former proof uses the prime reciprocal sum on
\((\sqrt n,n]\), with a constant change at the lower
endpoint, and the latter uses the estimate for the total
length of prime gaps exceeding the square-root threshold.
The first condition has the same relative density among
even \(n=2h\), because \(q>C\sqrt{2h}\) is odd for large
\(h\), and then \(q\mid2h\) exactly when \(q\mid h\).
The second condition still has relative density one on
even \(n\) because its complement has density zero.
This proves the stated irreducibility density on even
orders. Every such even order certifies two consecutive
terminals, its own and its successor, and these terminal
pairs are disjoint as the even order varies.

