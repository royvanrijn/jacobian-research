# Radical counts on moving affine slices of primitive triples

A fixed linear relation between two letters of a primitive abc
triple gives a stronger exponent than the unrestricted count.
The coefficients and the right-hand side may move with height.
The argument retains the third radical cap to deal with small
primitive quotients. When the affine intercept is numerically
short compared with the two weighted values, a real-strip
argument gives a further saving.

## A fixed weighted shift inside the three-letter count

Fix \(K,C_0>0\). Let \(L>0\) be small and fixed before
\(X\to\infty\), and write \(n=X^L\).
Use positive primitive triples \(a+b=c\) of height
\(X/2<c\le X\) with individual radical caps
\[
                  \operatorname{rad}(T)\le n^{l_T},
            \qquad T=a,b,c,\qquad 0\le l_T\le C_0 .
\]
Choose two different letters \(U,V\). Prescribe nonzero
integers \(p,q,t\), each of size at most \(X^K\), and
require
\[
                              pU+qV=t.
\]
For every \(\zeta,\varepsilon>0\), at all sufficiently small
fixed \(L\) the count is
\[
       O_{K,C_0,L,\zeta,\varepsilon}
             \bigl(X^{((l_U+l_V)/3+\zeta)L+\varepsilon}\bigr).
                                                               \tag{AS1}
\]
The constants and budget threshold are uniform in \(p,q,t\)
and in the indicated radical caps. If, for a fixed
\(\tau>0\), the sample also satisfies
\[
          |t|\le X^{-\tau}\max(|p|U,|q|V),
\]
then the count is
\[
       O_{K,C_0,L,\zeta,\tau,\varepsilon}
             \bigl(X^{(\min(l_U,l_V)/2+\zeta)L+\varepsilon}\bigr).
                                                               \tag{AS2}
\]
The third letter can have any height and radical within its
stated cap. Its exact value is determined by the other two.

### The small-quotient case

Suppose (AS1) fails with a fixed margin. Refine anatomic
coordinates on all three original letters at each fixed \(L\),
then take a sequence of violating fixed budgets \(L\to0\).
Use normalized limsup count exponent \(d\) in base \(n\),
with \(d>(l_U+l_V)/3\) by a fixed positive amount, writing
\(l_T\) also for the limiting caps. The root-dictionary
exponents are bounded above by the original radical caps
up to \(o(1)\) in these units. All coefficient-height and
root lists are handled with the height limit first.
[Exact valuation measures in anatomic boxes](ABC_SECANT_ENDPOINT_TRANSFER.md#the-height-first-setting)

Apply (NS1) to the two positive weighted terms
\(|p|U,|q|V\), with \(|p|,|q|\) as fixed coefficients
and signs as prescribed. Common divisors are permitted in
(NS1). For any fixed \(a_0>0\), let \(b_U(a_0),b_V(a_0)\)
be limiting upper caps for the original root mass in
degrees \(eL<a_0\). The middle bound in (NS1) gives
\[
              d\ \le\
       \max\{(l_U+l_V)/3,\ \min(b_U(a_0),b_V(a_0))\}.
\]
Indeed its small-quotient limb counts each complete value
by the low positions after stripping the prescribed
shift supports. Stripping can only lower the original
dictionary caps. In the positive-quotient limb the sum
cost is used instead.
[Arbitrary-shift bound (NS1)](ABC_NONPRIMITIVE_FIXED_SHIFT_ENERGY.md#two-dictionary-bounds-for-any-nonzero-shift)

Take \(a_0\) through continuity cutoffs decreasing to zero in
the two original limiting valuation measures. Since the
left side exceeds the first limb, each original measure
has an atom at zero of mass at least \(d\). Their combined
positive-scale mass is therefore at most \(l_U+l_V-2d\).
Apply the prebalance positive-scale inequality (IP1) to
the actual primitive abc triples, with their whole third
letter tracked. This gives
\[
                       d< l_U+l_V-2d,
\]
contradicting the strict violation. The third cap is used
here even though the weighted shift involves only two
letters. This proves (AS1).
[Positive-scale pair inequality (IP1)](ABC_IMBALANCED_PRODUCT_RADICAL_GROWTH.md#a-degree-measure-bound-before-balance)

For (AS2), strip and fix the root factors supported on
the primes of \(t\), using subpower enumeration
uniformly at polynomial height. Fix
\(g=\gcd(|p|U,|q|V)\) at divisor cost, since \(g\mid t\).
After division by \(g\), the two primitive complete
terms retain the given root lists on the remaining
supports and arbitrary fixed integer coefficients. The
relative inequality on the shift is unchanged by the
division. Apply (FS1) in the real-gap case. This charges
one half of the smaller of the two remaining dictionary
costs, bounded by \(\min(l_U,l_V)/2\). This proves (AS2)
with the original radical caps.
[Fixed-shift real gap (FS1)](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#a-fixed-shift-count-without-a-radical-cap-on-the-shift)

### A restriction at equality in the balanced affine bound

Suppose \(l_U=l_V=1/3\) in a limiting sequence with
\(d=2/9\). For each fixed fibre, normalize the terms as
above by fixing their \(t\)-supported parts and dividing
out \(g\). On a subsequence retaining this exponent, the
maximum primitive quotient height has positive limiting
logarithm relative to \(X\). If it tends to zero, the
small-height branch of the proof of (NS1) makes each
original atom at zero at least \(d\), and (IP1) then
gives the strict contradiction \(d<2/3-2d=d\).
With a positive quotient height, the positive-height
branch of (NS1) gives \(3d\le w_U+w_V\) for the stripped
dictionary costs, so they both equal \(1/3\). In
particular no fixed positive radical mass on either free
letter can be lost in the support stripping or in its
dictionary caps.

The residue-concentration argument (FT3)–(FT4) applies
since \(d>\min(w_U,w_V)/2\). Neither residual variable
part can have vanishing numerical height. To see this,
if one does and has a positive-height fixed coefficient,
freeze the complement of a positive-weight root on the
other side at cost at most \(1/3-d<d\).
The first fixed coefficient is a modulus larger than
the remaining short root. The unit power congruence has
only subpower many residue branches, and therefore
subpower solutions. If both variable parts vanish,
a positive fixed coefficient on one of the two sides
directly reconstructs the opposite small variable part
as an integer in a specified residue class. If one
whole complete term tends to numerical height zero,
the other is tall in a fixed real strip around the
shift; the bounded-root concentration and (RUW1)
give subpower cost there. These are the respective
terminal arguments of (NS1), and all have strict
margins at \(d=2/9\).

Consequently on each side there is exactly one original
position with positive limiting numerical weight. Each
such position costs at least \(d\) by (FT4), and two
would exceed the cap \(1/3\). Fixing the complement and
releasing the two high roots by the irreducible affine
curve count in the proof of (FT1) gives
\[
                         d\le w_U+w_V-w_{U,\mathrm{high}}
                                         -w_{V,\mathrm{high}} .
\]
Since \(w_U+w_V=3d\), both high costs equal \(d=2/9\).
The raw caps agree with the dictionary costs on these
positions because the raw and charged sums agree.
Their scaled orders have finite positive limits along
a subsequence, since their individual numerical weights
are positive and bounded. The other positions have
vanishing numerical weight after stripping. Returning
to the original weak measures does not change their
mass, because stripping has vanishing normalized radical
cost. Thus each of the two original measures has one
positive atom of mass \(2/9\) and an atom of mass \(1/9\)
at zero.

If the third cap also is \(1/3\), its nonatomic part
must vanish. Otherwise a cheap full modulus there would
release the two atomic free letters via (GE1) at vanishing
normalized cost. Once all three measures are atomic, an
escaped positive moment on any one letter likewise gives
an arbitrarily cheap full modulus and a contradiction.
Every letter has limiting numerical height one: if one
were short, (RUW1) would release the two atomic tall terms.
Locking any positive-weight original position on the
third letter and using (GE1) now shows that its cost
is at least \(d\). There is only one such position.
If a noncoordinate line supports this sequence, use the
first two-letter equality argument with each pair. Thus
on each letter the measure is exactly
\[
                          \frac19\delta_0+\frac29\delta_{9/2}.
\]
The positive atom has order \(9/2\) because its first
moment is one and no moment escapes.

There must also be one prime of mass \(2/9-o(1)\)
inside that positive atom on each letter, on a subfamily
of the full limiting exponent. In fact a positive
fraction of that root split among primes, with largest
prime of mass at most \(2/9-\xi\) for fixed \(\xi>0\),
would supply a full divisor with normalized radical cost
bounded away from both zero and \(2/9\). To select it,
take a prime of intermediate mass if one exists, or add
small primes greedily to a fixed small threshold. The
root occupies one original order with scaled order
tending to \(9/2\). Such a full divisor therefore has
positive numerical weight and fullness of order \(1/L\).
Enumerate its exact value by its radical, and apply
(GE1) to the two atomic opposite letters. This costs less
than \(d\) by a fixed amount, uniformly over the
selected divisors. Discard these smaller families at
each fixed tolerance and then decrease the tolerance
along the budgets. This proves the single-prime
assertion.

These restrictions exclude equality on the affine line
in sectors separated from the displayed weak measures
or from their single-prime realization by a fixed margin.
[Two-root release (FT1)](ABC_FIXED_SHIFT_TWO_ROOT_BOUND.md#releasing-two-original-roots)
· [Atomic grouping and full locks (GE1)](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)
· [Real-strip release (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)
[Fixed-shift residue concentration (FT3)–(FT4)](ABC_FIXED_SHIFT_TWO_ROOT_BOUND.md#one-fixed-residue-forces-concentration)
· [Terminals with arbitrary coefficients](ABC_NONPRIMITIVE_FIXED_SHIFT_ENERGY.md#positive-height-quotients-with-arbitrary-coefficients)

## Noncoordinate lines and projection sizes

Let a single affine line in the plane \(a+b=c\) vary
with \(X\), and suppose that none of the three letters
is constant on it. If it is a line through the origin,
it contains at most one primitive positive triple. Otherwise
each pair of letters admits a fixed nonzero relation with both
coefficients nonzero. If the line contains two integral triples
in the shell, take a primitive equation for its two-letter
projection through those points. Its coefficient height is
\(O(X^2)\), irrespective of a supplied representation of the line.
Apply (AS1) to the two least expensive letters. With total
radical at most \(X^L\), finite radical bins then give
\(X^{(2/9+\zeta)L+\varepsilon}\) points on the line.
More generally, with unequal individual caps the exponent
is the sum of the two smallest caps divided by three.

An explicit inverse consequence concerns one projection rather
than one fibre. For a given pair \(U,V\) fix the nonzero
coefficients \(p,q\), which may vary with \(X\). A family
of \(n^{1/3-o(1)}\) triples in balanced input boxes has
at least \(n^{1/9-o(1)}\) distinct values of \(pU+qV\).
The zero level has at most one primitive triple. If the
whole projected family satisfies the quantitative shortness
condition of (AS2), then the lower projection exponent is
\(1/6\). Equivalently, a union of at most \(n^\gamma\)
noncoordinate affine lines gives an exponent saving below
\(1/3\) for \(\gamma<1/9\) by a fixed margin; in the
short-intercept case for the requisite cheapest pair,
\(\gamma<1/6\) suffices. This counts the lines by their
actual number and keeps the height and radical constraints
on all three letters.

## Application to the balanced Pell slice and moving offsets

On the family
\[
                     (a,b,c)=(m,m+\kappa,2m+\kappa)
\]
assume \(\kappa\ne0\), all three terms are positive, and the
triple is primitive. Every two letters have a weighted
difference \(\pm\kappa\), with weights in \(\{1,2\}\).
Consequently a single prescribed \(\kappa\) of any
polynomial height has, in a height shell, count
\(O(X^{(2/9+\zeta)L+\varepsilon})\) at total radical budget
\(X^L\). If \(|\kappa|\le X^{1-\tau}\) for fixed \(\tau>0\),
all three letters are comparable to the height there and
(AS2) gives \(O(X^{(1/6+\zeta)L+\varepsilon})\).
For a prescribed fixed nonzero \(\kappa\), primitivity
can be dropped by splitting the possible common divisor
\(g\mid |\kappa|\), dividing the letters by it, and using the
same estimate with \(\kappa/g\). The fixed height factor
and the adjustment from \(c^\lambda\) to a dyadic cap
are absorbed by an arbitrarily small error. In particular,
for fixed nonzero \(\kappa\) and every \(\zeta>0\), for all
sufficiently small fixed \(\lambda>0\) the positive family
with \(\operatorname{rad}(m(m+\kappa)(2m+\kappa))<
(2m+\kappa)^\lambda\) has count
\[
                   O_{\kappa,\lambda,\zeta,\varepsilon}
                           (X^{(1/6+\zeta)\lambda+\varepsilon})
\]
up to maximum height \(X\). This strengthens the coefficient
in Zeng's fixed unit-offset small-budget result and covers
the fixed nonunit offsets in his Remark 6.2. The argument
uses the simultaneous radical caps and the relative
height of the weighted shift. The square-Pell solution
may have a large fundamental unit, as in the example
in that remark.
[Zeng, Theorem 1.1 and Remark 6.2](https://hexagonmath.org/pdf/2611.00611v1)

For a varying primitive offset \(|a-b|\le X^{\theta L}\)
in a height shell, sum the short-intercept estimate over
the possible offsets. This gives coefficient at most
\(1/6+\theta+\zeta\), a strict saving when \(\theta<1/6\).
A large unrestricted family need not have such a short
offset or a small projection dictionary. The unrestricted
saving requires an arithmetic constraint controlling the
number or nature of its affine slices.
