# Independent real and finite-prime approximation of adjacent factorial coefficients

Fix coprime orders \(n,N\ge2\), a finite nonempty set of primes \(\mathcal P\),
and an integer precision \(L_p\ge1\) at each of those primes.
For an ordinary two-monomial exponent profile write
\[
 R_{j,k}=\frac{\prod_i(ja_i+k(b_i-a_i))!}{\prod_i(ja_i)!},
       \qquad j=n,N,\quad 1\le k\le j .
\]
Thus \(\binom jk R_{j,k}\) is the \(k\)-th coefficient after
dividing the \(j\)-th moment by its constant coefficient. For a
nonzero rational \(x\), use
\(\operatorname{unit}_p(x)=p^{-v_p(x)}x\) in the unit group
modulo \(p^{L_p}\).

**Theorem.** Starting from any ordinary profile, one can prescribe
arbitrary integers \(v_p(R_{j,k})\) and arbitrary unit residues
\(\operatorname{unit}_p(R_{j,k})\bmod p^{L_p}\), simultaneously
for all the indicated coordinates. To do so, append a finite
multiset of ordinary shared rows of signed total gap zero.
All added rows may use the same positive absolute gap, depending
only on \(n,N\), and may have arbitrarily large smaller exponents.
Their total curvature can be arbitrarily small. In addition,
all the new \(R_{j,k}\) can be made arbitrarily close to their
starting values over the reals.

The result applies to the original factorial coefficients in
the same variable. For example, the prescribed data can all be
zero valuations and unit residues one. At every selected prime
the two resulting constant-normalized moments then agree with
\((1+z)^n,(1+z)^N\bmod p^{L_p}\), while the real coefficient pair
approaches the starting pair and its signed degree difference
is preserved. The
[simultaneous prime-screen obstruction](MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md#a-simultaneous-obstruction-at-finitely-many-primes)
uses reflected rows to obtain binomial residues in a prescribed
curvature band. The construction here uses an integer basis of
coefficient valuations and units and allows an arbitrary starting
profile and arbitrary local targets.

## A basis of factorial factor indicators

Choose a common absolute gap \(d\) large enough that each interval
between consecutive elements of the ordered set
\[
       \{0,1\}\cup\{k/n:0<k<n\}\cup\{k/N:0<k<N\}
\]
contains a point of the form \(u/d\) with \(u\) integral.
Also take \(d>\max(n,N)\).
A row of smaller exponent \(a\) has the full order-\(j\) ramp
\(ja+h\), \(1\le h\le jd\). The increasing row uses the product
up to \(kd\) in \(R_{j,k}\). The decreasing row uses the inverse
of the product from \((j-k)d+1\) through \(jd\).

A marked factor at \(h=ju\) has exponent in the increasing
ratio equal to \({\bf1}_{u/d<k/j}\), since the chosen fraction
avoids every threshold. Order the internal coordinates by their
distinct fractions \(k/j\). The vectors obtained by varying \(u\)
in the chosen intervals are successive suffix indicators on
these coordinates and are one at both endpoint coordinates.
Consecutive differences and the last vector give the individual
internal coordinates and the sum of the endpoint vectors.
These form an integer basis of the vectors with equal endpoint
coordinates. The reflected row with its mark at \(h=j(d-u)\)
gives the negative of the same indicator.

To supply the remaining endpoint direction at a target prime \(p\),
choose an order \(j_*\in\{n,N\}\) with \(p\nmid j_*\).
Choose an integer \(h_*\) with
\((j_*-1)d<h_*<j_*d\) and \(j_*\nmid h_*\), and put
\(\alpha=h_*/j_*\).
This pole belongs to the ramp at order \(j_*\) only.
Indeed integrality at the other order would force
\(j_*\mid h_*\) since the orders are coprime.
Its indicator in the increasing ratios is one at that
endpoint and zero elsewhere. At pole \(d-\alpha\) the
decreasing ratios give its negative. Adding this private
endpoint indicator completes an integer basis at \(p\).

## Exact changes in valuations and units

A pole \(\alpha\) in this construction has denominator prime to \(p\).
For a large integer \(H\) and a unit \(c\bmod p^{L_p}\),
choose a positive integer base \(a\) such that
\[
                   a+\alpha\equiv p^H c\pmod{p^{H+L_p}} .
\]
At every other specified prime \(q\) take \(a\equiv0\bmod q^{K_q}\),
where \(K_q\) is large enough that the valuation and unit residue
of each \(ja+h\) are those of \(h\). This holds whenever
\(K_q\ge L_q+\max_{j,h}v_q(h)\), and the Chinese remainder
theorem gives all the required bases.

At \(p\), write
\[
                      ja+h=j(a+\alpha)+(h-j\alpha).
\]
For every nonzero offset on the right, fixed \(n,N,d,\alpha\)
allow \(H\) large enough that its valuation and unit residue
are unchanged by the first summand. For each zero offset the
valuation is \(H+v_p(j)\) and the unit residue is
\(\operatorname{unit}_p(j)c\).
Consequently the valuation vector of an increasing row
at depth \(H+1\) minus the vector at depth \(H\), both with
\(c=1\), is exactly the chosen basis indicator in the target
prime's coordinates. The unit data agree everywhere.
At fixed depth, changing \(c=1\) to \(c\) multiplies the
unit data by \(c\) on that indicator, with no valuation
change. The same rule covers each basis indicator, including
the private endpoint at every target prime.

Allow formal integer multiplicities for the moment.
The preceding differences can realize any integer changes
of all the valuations, since the indicators are an integer
basis. They can also realize arbitrary unit changes. The inverse integer
basis matrix acts on unit tuples by integer powers.
Every formal difference uses two increasing rows with the
same gap and hence has zero signed total gap.

## Passing from formal differences to actual rows

Record a row by the vector of all its valuations and its
signed gap divided by \(d\). Use the increasing and decreasing
marker rows for every basis direction at every prime, with
\(c=1\). Also take two rows of opposite orientation whose base
is large and divisible by the prescribed powers of all primes,
which keeps their valuations bounded independently of \(H\).

Divide the valuation coordinates by \(H\) and leave the
last coordinate unchanged. As \(H\) tends to infinity, the
increasing and decreasing vectors for a basis direction \(w\)
tend to \((w,1)\) and \((-w,-1)\), respectively. The two
additional rows tend to \((0,1),(0,-1)\). The increasing
limits form a basis, and the decreasing limits are their
negatives. For sufficiently large integral \(H\), the
increasing vectors still form a basis. The unique combination
of those vectors canceling the sum of the decreasing vectors
has all coefficients strictly positive, since those coefficients
tend to one. These coefficients are rational.

One can include every auxiliary row used for the preceding
formal differences in a positive rational combination with
zero valuation vector and zero signed gap. Give the auxiliary
rows weight one and cancel their sum by a rational combination
of the increasing basis vectors. Adding enough of the
strictly positive null combination just constructed makes
all the weights positive. Clear denominators to obtain a
positive integral multiset. Repeat it by a common multiple
of the orders of the unit groups modulo \(p^{L_p}\).
It then contributes zero to all the valuations and one to
all the unit residues, and has signed gap zero.
Adding sufficiently many copies of this multiset to the
chosen formal differences makes all multiplicities nonnegative,
without altering the prescribed changes.

Finally add the same large multiple \(T\) of a common modulus
to every smaller exponent just constructed. Choose the modulus
divisible by \(p^{L_p+\max v_p(ja+h)}\) at each selected prime,
where the maximum runs over all added rows and their sampled
factors before shifting. This keeps all valuations and unit
residues unchanged. The signed gap is also unchanged. Every
added row has curvature tending to zero. At order \(j\) an
added increasing row contributes \((jT)^{kd}(1+O(1/T))\)
to the raw ratio \(R_{j,k}\), and a decreasing row contributes
\((jT)^{-kd}(1+O(1/T))\). The finite number of powers cancel
because the added signed gap is zero, and the product tends to
one. This proves the real approximation and completes the proof.

## Simultaneous all-degree angular separation

The approximating profiles in the theorem can also be chosen to
have disjoint consecutive root arguments in every degree. Use the
same gap \(d\) in the interpolation and in an additional mirrored
pair, increasing \(d\) beforehand if needed to be at least \(J-1\).
Here \(J\) is the absolute angular cutoff in the
[universal opposed-row extension](MULTIVARIATE_UNIVERSAL_OPPOSED_PADIC_ANGLE_MARKING.md).

Choose a fixed positive reference base \(A_0\) for that mirrored
pair. Compute the pair's local factorial-ratio multipliers at each
of the original places. Prescribe preliminary targets for the
interpolation by dividing the desired local multipliers by these
reference multipliers, in the group of valuation and unit data.
The theorem gives an intermediate ordinary profile which is as
close as desired to the starting one over the reals. Its added
signed total gap is zero.

Now choose a prime \(r\) outside the prescribed places. Apply the
universal opposed-row extension to this intermediate profile, with
the chosen equal gaps and a common base \(A\equiv-1\bmod r^E\),
where \(E\) is the threshold for that profile. At each original
place \(p\), also require \(A\equiv A_0\bmod p^{K_p}\), with
\(K_p\ge L_p+\max_{j,h}v_p(jA_0+h)\) over the tested ramps of
the mirrored pair. These conditions preserve the precompensated
local data of the pair after it is inserted. The Chinese remainder
theorem gives arbitrarily large such \(A\). The pair adds no signed
gap, and its raw real multipliers tend to one at each fixed
coefficient as \(A\) increases. The opposed-row theorem supplies
the all-degree conclusion simultaneously.

## Consequence for a finite-cover approach

For a fixed adjacent pair, the joint closure of the real
normalized coefficient data and the local raw coefficient data
contains the product of the [real compactified closure](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md)
and all nonzero local coefficient tuples. The assertion holds at fixed
signed total gap. To see it, take ordinary profiles approximating
the given real compactified point, and apply the theorem to one
profile at a time with any finite collection of local boxes and
arbitrarily small real perturbation. This gives the defining
neighborhood condition in the product topology.

In particular, if a real compactified point has a common root,
adding finitely many coefficient-valuation or finite-unit tests
cannot exclude its product neighborhood by a finite precision
local coprimality criterion. Every neighborhood also admits
ordinary profiles with the local data of the binomial common-root
pair at the tested places. This result addresses coefficient-box
certificates. The possibility of exact constraints using the
integer row data, or of prime and precision choices growing
with the profile, remains open.

The coprimality hypothesis on the two orders is used to separate
their internal fractions and to produce a private endpoint.
For comparison, at multiple orders there are actual constraints.
For example, for a prime \(p\) and order \(j\),
\[
 v_p(R_{pj,pk})-v_p(R_{j,k})=k\sum_i(b_i-a_i)
\]
by \(v_p((p x)!)=x+v_p(x!)\). Thus independent interpolation
cannot be asserted for arbitrary collections of orders.

## Exact replay

The following test uses orders three and four and primes two and
three, so the private endpoint occurs at a different order at the
two primes. It checks the two integer bases, exact depth and unit
changes, the positive null combination, and a nontrivial prescribed
valuation and unit array after removing all signed multiplicities.
Factorials are evaluated through their short factor ramps.

~~~python
from collections import Counter
from fractions import Fraction as Q
from math import lcm, prod

import sympy as sp
from sympy.ntheory.modular import crt

orders, primes = (3, 4), (2, 3)
precision, gap, depth = {2: 3, 3: 2}, 24, 256
coordinates = [(j, k) for j in orders for k in range(1, j + 1)]
s = len(coordinates)
dimension = 2*s+1
target_val = [(-1)**i*(1+i % 4) for i in range(2*s)]
target_unit = [1, 3, 5, 7, 5, 3, 1, 1, 2, 4, 5, 7, 8, 4]


def marked_row(p, alpha, orientation, H=depth, c=1):
    moduli, residues = [], []
    for q in primes:
        if q == p:
            modulus = p**(H+precision[p])
            residue = (p**H*c-alpha.numerator *
                       pow(alpha.denominator, -1, modulus)) % modulus
        else:
            modulus, residue = q**(precision[q]+8), 0
        moduli.append(modulus)
        residues.append(residue)
    base = int(crt(moduli, residues)[0])+prod(moduli)
    return base, orientation


def row_data(row):
    a, orientation = row
    valuations, units = [], []
    for p in primes:
        modulus = p**precision[p]
        for j, k in coordinates:
            factors = (range(1, k*gap+1) if orientation == 1
                       else range((j-k)*gap+1, j*gap+1))
            total_val, total_unit = 0, 1
            for h in factors:
                y, v = j*a+h, 0
                while y % p == 0:
                    y //= p
                    v += 1
                total_val += orientation*v
                total_unit = total_unit*pow(y, orientation, modulus) % modulus
            valuations.append(total_val)
            units.append(total_unit)
    return sp.Matrix([*valuations, orientation]), units


fractions = sorted({Q(0), Q(1)}
                   | {Q(k, j) for j, k in coordinates if k < j})
common_poles = [Q(next(u for u in range(1, gap) if a < Q(u, gap) < b))
                for a, b in zip(fractions, fractions[1:])]
positive_rows, negative_rows, auxiliary = [], [], []
formal = Counter()
zero = sp.zeros(dimension, 1)

for place, p in enumerate(primes):
    j = next(j for j in orders if j % p)
    h = next(h for h in range((j-1)*gap+1, j*gap) if h % j)
    poles = common_poles+[Q(h, j)]
    W = sp.Matrix([[int(j*a % 1 == 0 and j*a < k*gap) for a in poles]
                   for j, k in coordinates])
    assert abs(W.det()) == 1
    inverse = W.inv()
    i0, modulus = place*s, p**precision[p]
    powers = inverse*sp.Matrix(target_val[i0:i0+s])
    for i, pole in enumerate(poles):
        c = prod(pow(target_unit[i0+k], int(inverse[i, k]), modulus)
                 for k in range(s)) % modulus
        plus, minus = marked_row(p, pole, 1), marked_row(p, gap-pole, -1)
        deeper = marked_row(p, pole, 1, H=depth+1)
        unit = marked_row(p, pole, 1, c=c)
        positive_rows.append(plus)
        negative_rows.append(minus)
        auxiliary.extend((deeper, unit))
        base_val, base_unit = row_data(plus)
        deep_val, deep_unit = row_data(deeper)
        unit_val, unit_unit = row_data(unit)
        indicator = sp.zeros(dimension, 1)
        indicator[i0:i0+s, 0] = W[:, i]
        assert deep_val-base_val == indicator and deep_unit == base_unit
        assert unit_val == base_val
        for t in range(2*s):
            q = primes[t//s]
            M = q**precision[q]
            assert unit_unit[t]*pow(base_unit[t], -1, M) % M == (
                pow(c, int(indicator[t]), M) if q == p else 1)
        formal[deeper] += int(powers[i])
        formal[unit] += 1
        formal[plus] -= int(powers[i])+1

generic = 2**20*3**20
positive_rows.append((generic, 1))
negative_rows.append((generic, -1))
P = sp.Matrix.hstack(*(row_data(row)[0] for row in positive_rows))
negative = sum((row_data(row)[0] for row in negative_rows), zero)
positive_null = -P.inv()*negative
assert all(x > 0 for x in positive_null)
extra = sum((row_data(row)[0] for row in auxiliary), zero)
correction = -P.inv()*extra
repeat = 1+max(0, int(max(-a/b for a, b in zip(correction, positive_null))))
leading = correction+repeat*positive_null
assert all(x > 0 for x in leading)
assert P*leading+repeat*negative+extra == zero
multiplier = lcm(*(int(x.q) for x in leading), 1)*lcm(4, 6)
null = Counter()
for row, count in zip(positive_rows+negative_rows+auxiliary,
                      list(leading)+[repeat]*dimension+[1]*len(auxiliary)):
    null[row] += int(multiplier*count)
add = 1+max(0, max(-count//null[row] for row, count in formal.items()))
actual = Counter({row: add*count for row, count in null.items()})
actual.update(formal)
assert min(actual.values()) >= 0
assert sum((count*row_data(row)[0] for row, count in actual.items()), zero) == (
    sp.Matrix([*target_val, 0]))
for t in range(2*s):
    p = primes[t//s]
    modulus = p**precision[p]
    assert prod(pow(row_data(row)[1][t], count, modulus)
                for row, count in actual.items()) % modulus == target_unit[t]
print("Exact two-prime prescribed coefficients and zero total signed gap passed")
~~~
