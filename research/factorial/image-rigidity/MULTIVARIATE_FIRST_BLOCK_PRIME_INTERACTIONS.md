# Large prime interactions forced by a consecutive factorial collision

In the original ordinary two-monomial problem, discard equal rows.
For an unequal row write \(L\) for its smaller exponent,
\(d\ge1\) for its gap, and retain its orientation. Multiplicities
are positive integers. A row of a given type can occur any number
of times. Fix a computable
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
\(J\ge4\). There are no restrictions on the number of heights,
the relative heights of different rows, or the multiplicities.

For \(m\ge2\) and a row \(w\) with data \((L,d)\), define its full
factor block and its selected first block by
\[
 {\cal R}_m(w)=\{mL+1,\ldots,mL+md\},\qquad
 {\cal A}_m(w)=\{mL+h:1\le h\le d,\ h\equiv1\pmod m\}.
\]
For an adjacent pair \(m,q\), say that another row blocks a prime
from the first block of \(w\) in order \(m\) if one of the following
holds.

- It has opposite orientation from \(w\), and its \(m\)-th full
  block contains a multiple of \(p\).
- It has the same orientation as \(w\), and its \(q\)-th full
  block contains a multiple of \(p\).

The role of order \(m\) can be taken by either of the two orders.
Identical copies of \(w\) are harmless in the criterion below.

**Theorem.** Suppose some adjacent orders \(n,N=n+1\) have a
common root. Necessarily \(3\le N<J\). At these orders the
following hold for every unequal row \(w\).

1. For each \(j=n,N\), every prime \(p>J^2d\) dividing a
   member of \({\cal A}_j(w)\) is blocked.
2. If \(d\ge 2^{4000J^3}\), put
   \[
       B(d)=\left\lfloor
         \frac{d\lfloor\log_2 d\rfloor}{2000J}\right\rfloor .
   \]
   At each of \(n,N\), strictly more than \(d/(8J)\) distinct
   blocked primes exceed \(B(d)\) and divide the selected first
   block. The two sets of primes are disjoint.
3. If all the blocked primes from both selected first blocks
   are at most \(B\), where \(B\ge J^2d\), then \(L\le M(B)\)
   for an effectively computable \(M(B)\) independent of the
   profile and the selected orders.

In particular, the original moments are coprime in every degree
if at each residual pair there is a row of sufficiently large
gap for which at one order at most \(d/(8J)\) selected first-block
primes above \(B(d)\) are blocked. The row may depend on the
adjacent pair. A bounded-gap row of sufficiently large height
likewise suffices when its blocked primes at both orders are
bounded uniformly over the residual pairs. These conclusions
are over every characteristic-zero field. They apply to the
physical common variable.

## The Newton polygon comparison

Take \(j\) as the marked order and \(q\) as the other order.
Reverse the two exponent columns if necessary so that \(w\)
increases. Strip every scalar factorial \((mL)!\) in order \(m\).
An increasing row then contributes the product from the first
\(dk\) entries of \({\cal R}_m\) at index \(k\). A decreasing row
uses the first \(d(m-k)\) entries. The corresponding valuations
are nondecreasing and nonincreasing in \(k\), respectively.

If \(p>J^2d\) divides \(jL+h\) with \(h\equiv1\pmod j\), it cannot
divide any \(qL+t\) with \(1\le t\le qd\) of the selected row.
Indeed, it would divide \(qh-jt\), whose absolute value is less
than \(J^2d\). This integer is nonzero because
\(\gcd(j,q)=\gcd(j,h)=1\). In particular identical copies of the
row cannot block \(p\). The first-block occurrence in order \(j\)
increases the marked valuation by at least one at every positive
index relative to index zero. No increasing row at that order
can decrease an ordinate relative to zero. If no opposite row
blocks the prime there, the order-\(j\) lower Newton polygon has
strictly positive slopes. On the other hand, if no same-orientation
row blocks at order \(q\), the valuation array of that order is
nonincreasing, so all its lower Newton slopes are nonpositive.
The binomial coefficients are units since \(p>J\). The Newton
polygon product rule excludes a common factor, proving part 1.
The nearest existing criteria use a marked first block in a fixed
or variable prime regime.
[Native upper chord theorem](MULTIVARIATE_NATIVE_PRIME_CHORD_ANGULAR_CERTIFICATE.md#checking-the-polygon-from-the-factorial-blocks) ·
[Variable-prime Eisenstein criterion](MULTIVARIATE_MARKED_ROW_VARIABLE_PRIME.md)

The comparison can use multiplicities at a candidate prime.
In the reflected orientation define \(A_m^+\) to be the sum
of the valuations of the entire first blocks of increasing
rows, with multiplicity, and \(E_m^-\) the sum of full-block
valuations of decreasing rows. If \(A_j^+>E_j^-\), every
positive-index ordinate in order \(j\) is strictly above the
ordinate at zero, whether or not an opposite row meets the
prime. In order \(q\), define \(A_q^-\) and \(E_q^+\) with the
roles of increasing and decreasing exchanged. If
\(A_q^-\ge E_q^+\), every ordinate before the leading one is
at least the leading ordinate, hence the lower polygon has
nonpositive slopes. Either this last inequality or the absence
of a same-orientation companion block suffices for the
comparison at \(q\).

## A uniform number of primes at a large gap

The
[absolute curvature bound](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
gives \(\sum_w d_w^2/(L_w+d_w)<12\) at any collision, with all
multiplicities included. Thus for the selected row of large
gap, \(L\ge d^2/24\). Write \(Y=jL+1\), \(X=jL+d\), and let
\(K=\lfloor(d-1)/j\rfloor+1\ge d/j\) be the number of terms of
\({\cal A}_j(w)\). Then \(X\le2Y\), all terms are at least \(Y\),
and any prime \(p\nmid j\) has at most \(K/p^a+1\) multiples of
\(p^a\) in that progression. No prime dividing \(j\) occurs.
Consequently, for \(B=B(d)\), the logarithm of the part of its
product supported on primes at most \(B\) is at most
\[
            K\sum_{p\le B}\frac{\log p}{p-1}
                      +\pi(B)\log X .                        \tag{1}
\]
For completeness, the elementary bounds
\[
 \sum_{p\le B}\frac{\log p}{p-1}\le\log B+12,\qquad
 \pi(B)\le 9B/\log B
\]
follow as follows. At integer \(b\ge1\) the prime powers in
\((b,2b]\) each contribute their \(\log p\) to
\(\log\binom{2b}{b}\le 2b\log 2\), with valuations counted
using the factorial floor formula. Summing in dyadic intervals
gives \(\psi(x)=\sum_{p^a\le x}\log p\le4x\). From the floor
formula for \(B!\),
\[
   B\sum_{p^a\le B}\frac{\log p}{p^a}
             \le \log(B!)+\psi(B)\le B(\log B+4).
\]
The geometric sums over \(p^a>B\), for \(p\le B\), contribute
at most \(2\sum_{p\le B}\log p/B\le8\). Finally split the prime
count at \(\sqrt B\) and use \(\psi(B)\le4B\) and
\(\log B\le\sqrt B\). This proves both bounds.

Our threshold and definition of \(B\) give
\(B\ge J^2d\) and \(B\le d\log d/(1000J)\). In particular,
\(\log B\ge\log d\) and
\(\pi(B)/K\le9/1000<1/100\).
We also have \(\log Y\ge2\log d-\log12\) and
\(\log B+12\le(2/3)\log Y\). For the latter use \(x=\log d\ge100\)
and \(\log x+12+(2/3)\log12\le x/3\), whose slack increases
for \(x\ge100\). Since \(\log X\le2\log Y\), (1) is strictly
less than \((3/4)K\log Y\). Therefore the sum of
\(v_p(\prod{\cal A}_j(w))\log p\) over \(p>B\) exceeds
\((1/4)K\log Y\).

Such a prime can occur in only one term of the progression:
the positive differences of its terms are at most \(d-1<p\).
Its logarithmic contribution is at most \(\log X\le2\log Y\).
Thus there are more than \(K/8\ge d/(8J)\) distinct primes
above \(B\). Every one is blocked by part 1. The cross-order
determinant from that proof also shows that none occurs in the
selected first block in the other order. This proves part 2.
The quantitative claim uses the entire selected first-block
product and remains valid when every individual marker is
composite.

## Bounded gaps and a varying number of heights

The first factors \(nL+1\) and \(NL+1\) are selected at their
respective orders. By part 1, if all their blocked prime
divisors are at most \(B\ge J^2d\), both factors are
\(B\)-smooth. The two \(B\)-smooth integers
\[
              x=n(NL+1),\qquad x+1=N(nL+1)
\]
are consecutive. Let \(M(B)\) be the largest first entry of a
consecutive \(B\)-smooth pair. It is computable by the finite
Pell procedure of
[Lehmer, Theorem 1](https://doi.org/10.1215/ijm/1256067456).
As \(x\ge L\), this proves part 3.
More generally, for bounded \(d\) the permitted bound \(B\)
can grow with \(L\) as
\(o(\log_2^* L\,\log_3^* L/\log_4^* L)\), where here
\(\log_k^* L\) denotes the \(k\)-fold iterated natural logarithm.
This notation is unrelated to the base-two logarithm used in
\(B(d)\). To see this, take \(y(t)=\max(3,B(t))\) nondecreasing
and apply
[Stewart, Theorem 1 and (4)](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/on_sequences-2.pdf)
to the consecutive integers above. The hypothesis implies
that the gaps between successive integers with greatest prime
factor at most \(y(t)\) tend to infinity. The factors \(n,N\)
are fixed within the finite menu of orders and
\(x,x+1\) are comparable to \(L\), so the small-prime cutoff
\(B(L)\) is within this hypothesis.

For example, fix any bounded-height retained profile and
a tagged bounded-gap row with variable smaller exponent \(L\).
Choose \(B\) greater than \(J\) times the retained maximum
height and greater than \(J^2d\) and \(JD\). Here \(D\) is a
fixed gap cap for the remaining rows. Let \(P\) be the product
of primes above \(B\) dividing any selected first block of
the tagged row at \(2\le j<J\). Every other varying row may
independently have any smaller exponent divisible by \(P\),
in either orientation, and any gap at most \(D\). Indeed its
whole factor block is a unit modulo each \(p\mid P\), since
its entries are congruent to integers in \([1,JD]\).
Therefore the bounded-gap conclusion holds with any number
of heights and any multiplicities for \(L>M(B)\), even if
their total curvature tends to a positive limit. A similar
congruence gives the large-gap conclusion by taking a tagged
row with \(B(d)\) exceeding the fixed core and other gaps in
the same sense. These conditions isolate one arithmetic
interaction in the exact arrays. Cases where sufficiently
many of the first-block primes are blocked remain possible
under the conclusion of this theorem.
