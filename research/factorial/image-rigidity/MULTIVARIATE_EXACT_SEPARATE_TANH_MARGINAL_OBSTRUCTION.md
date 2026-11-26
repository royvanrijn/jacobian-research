# Exact separate tanh marginals at finitely many degrees do not determine target moments

Let \(I\) and \(J\) be disjoint nonempty finite sets of integers
at least two. For a normalized ordinary binomial moment of order
\(j\), its separate tanh marginal means the multiset of all Gauss
values \(T_{j,h}\), keeping the increasing and decreasing
orientations separate.

**Theorem.** There are two distinct fully shared ordinary
two-monomial profiles whose separate tanh marginals agree exactly at
every order in \(I\). Their row counts, multisets of absolute gaps,
and two respective monomial total degrees agree. At each target
order in \(J\), their separate tanh marginals and every interior
endpoint-normalized factorial coefficient differ. Their full raw
moment polynomials can be made equal at every observed order by
adjoining equal rows. The profiles can have arbitrarily small
positive total curvature.

The two profiles can also be taken exactly mirrored, with degree
difference zero. In that case their normalized target polynomials
have different multisets of root arguments at every order in \(J\).
In particular the entire collection of separate marginals at a
prescribed finite set of orders, with the physical degree difference
given, need not determine the normalized polynomial or its argument
multiset at a previously unobserved order.

## Row differences and the cofactor construction

Write a row by its lower exponent \(\ell\) and positive gap \(d\).
Use increasing rows first. At order \(j\) its Gauss parameters
\[
                 x_{j,h}=(j\ell+h)/d,\qquad 1\le h\le d,
\]
give \(T_{j,h}=(j-1)/(j-1+2x_{j,h})\). Thus a separate tanh
marginal is equivalent to the multiset of the \(x_{j,h}\).
This change of coordinates is the one-order specialization of the
[joint Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid).

Record the marginal for the row as
\[
 S_{\ell,d}^{(j)}(X)=\sum_{h=1}^d X^{(j\ell+h)/d}.
\]
These are polynomials in a common fractional power of \(X\) for any
finite collection. The geometric sum identity gives
\[
 S_{\ell,d}^{(j)}-S_{\ell+1,d}^{(j)}
       =(1-X)\sum_{h=1}^j X^{(j\ell+h)/d}.                       \tag{1}
\]
This follows by factoring \(1-X^{j/d}\) from the difference and
using \(S_{\ell,d}^{(j)}=X^{(j\ell+1)/d}(1-X)/(1-X^{1/d})\).

Put \(m=|I|\). In the commutative ring of finite integer sums of
symbols \([q]\), \(q\in\mathbb Q_{>0}\), with
\([q][q']=[qq']\), define
\[
 B_{j,\ell}=\sum_{h=1}^{j}[j\ell+h].
\]
Take \(m+1\) nonnegative integer baselines
\(\ell_0,\ldots,\ell_m\). For each \(i\) let \(U_i\) be the signed
cofactor \((-1)^i\) times the determinant of the \(m\)-by-\(m\)
matrix obtained by deleting column \(i\) from
\((B_{j,\ell_i})_{j\in I,\,0\le i\le m}\). The row order is
fixed once. Expansion along a last appended row shows that
\[
 \sum_{i=0}^m U_i B_{j,\ell_i}=0\qquad(j\in I).                   \tag{2}
\]
For any target \(j\) the same expression, up to the global sign
\((-1)^m\), is the determinant after appending that order.

All \(q\) in the \(U_i\) are positive integers. Choose an integer
\(D>0\) divisible by every such \(q\). For each term \(c[q]\) of
\(U_i\), with \(d=D/q\), take \(c\) times the signed row difference
\[
                           (\ell_i,d)-(\ell_i+1,d),
\]
where the second coordinate here is the gap. The positive and
negative multisets of the resulting signed row array give the two
profiles, omitting any common cancellations. Apply (1) to see that
their marginal difference at order \(j\) is \(1-X\) times the image
of \(\sum_i U_i B_{j,\ell_i}\) under \([q]\mapsto X^{q/D}\).
This map is injective on finite sums, and multiplication by
\(1-X\) is injective on these polynomials. Therefore (2) gives
the observed equalities, and each nonzero enlarged determinant
gives a target difference.

Every signed row difference has zero net count at its gap.
The profiles therefore have the same row count and gap multiset.
Their difference in total lower exponent is
\(-\sum_i U_i(1)\), where \(U_i(1)\) is the augmentation
\([q]\mapsto1\). For \(m\ge2\) each augmented cofactor vanishes
because each row of the augmented matrix is constant. For \(m=1\)
the two augmented cofactors cancel. Adding the common total gap
gives equality of the upper exponent totals as well.

## Nonvanishing at the specified targets

Take \(s=\max(1,2m-2)\) and evaluate \([q]\mapsto q^s\).
After dividing each row \(j\) by \(j^{s+1}\) its entry in a column
with indeterminate baseline \(L\) becomes
\[
 \frac1j\sum_{h=1}^j(L+h/j)^s
        =\sum_{r=0}^s\binom sr L^{s-r}\mu_r(1/j),
 \qquad \mu_r(1/j)=\frac1j\sum_{h=1}^j(h/j)^r.
\]
Faulhaber's formula gives
\[
     \mu_r(x)=\frac{1}{r+1}
                   \sum_{a=0}^r\binom{r+1}{a} B_a(1)x^a,
\]
where \(B_a(1)\) is the Bernoulli polynomial evaluated at one.
The constant, linear, and successive even-top degrees have
nonzero leading coefficients. The odd coefficients above degree
one vanish. These facts follow from the Bernoulli polynomial
power sum, the parity of Bernoulli numbers, and the Fourier
formula at zero in even degree. [DLMF power sum](https://dlmf.nist.gov/24.4.E9) · [DLMF Bernoulli values](https://dlmf.nist.gov/24.2) · [DLMF even Fourier series](https://dlmf.nist.gov/24.8.E1)

For \(m=1\) use indices \(r=0,1\). For \(m\ge2\) use
\(r=0,1,2,4,\ldots,2m-2\). There are \(m+1\) indices in either
case. Index the chosen powers by \(r_0,\ldots,r_m\). By column
multilinearity, the coefficient of
\(\prod_{i=0}^m\ell_i^{s-r_i}\) is
\[
   \left(\prod_{i=0}^m\binom{s}{r_i}\right)
      \det_{j\in I\cup\{j_*\},\,0\le i\le m}\mu_{r_i}(1/j),
\]
for target \(j_*\), in the row order with \(j_*\) appended.
Each polynomial \(\mu_{r_i}\) has nonzero leading coefficient
\(B_{r_i}(1)\). All its lower powers lie among the preceding
indices because the odd Bernoulli values above index one vanish.
The determinant is therefore a nonzero scalar times the
evaluation determinant of the monomials \(x^{r_i}\) at the distinct
positive nodes \(1/j\). In fact a nonzero
linear combination of \(h\) distinct monomials has at most
\(h-1\) distinct positive zeros. Divide by its smallest power
and differentiate, then use Rolle's theorem inductively.

Thus every target gives a nonzero determinant polynomial in the
baselines. The product over the finite target set is nonzero,
so some tuple of nonnegative integers avoids its zero set.
Each coordinate degree of the product is at most \(s|J|\).
Successive univariate interpolation guarantees a nonzero value on
the grid \(\{0,1,\ldots,s|J|\}^{m+1}\).
At these baselines each target cofactor
sum is nonzero in the original ring as well.

## Target coefficients and arguments

Replace each row \((\ell,d)\) by \((\ell+Cd,d)\), with a common
positive integer \(C\). This translates all Gauss parameters by
\(jC\) in order \(j\), and preserves every observed equality,
every target difference, the gap counts, and the equal degree
totals. The resulting rows are fully shared and their total
curvature tends to zero as \(C\) tends to infinity.

The binomial-stripped normalized coefficient at index \(k\) is
the product over the shifted Gauss parameters \(jC+x\) of
\[
                        \frac{(jC+x)_k}{(jC+x)_j^{\,k/j}}.
\]
Indeed \((j\ell+dk)!/(j\ell)!\) factors over the \(d\)
residue classes into \(d^{dk}\prod_h(x_{j,h})_k\), and
the endpoint normalization cancels the powers of \(d\).
Since the gap multiset already agrees, the observed polynomials
even agree after dividing by their constant coefficients, with no
variable rescaling. To make the raw polynomials equal, append to
each profile the full list of lower exponents from the other
profile as equal rows. An added exponent \(b\) gives the row
\((b,b)\), multiplying the raw order-\(j\) polynomial by \((jb)!\).
The padded constant coefficients at every order are then the
product of the two original constant coefficients. The padding
adds the same total degree to each column of each profile,
because the original lower exponent totals were equal. It adds no
curvature and has no Gauss slots.

At a target order \(j\), let \(\nu_j\) be the difference of the
two finite parameter counting measures before translation.
It is nonzero and has mass zero. Let \(p\ge1\) be its first
nonzero power moment, denoted by \(M_p=\int x^p\,d\nu_j(x)\).
Such a \(p\) exists by the ordinary finite Vandermonde matrix.
The log of the ratio of the normalized coefficients is
\[
 \int \left(\sum_{r=0}^{k-1}\log(jC+x+r)
      -\frac{k}{j}\sum_{r=0}^{j-1}\log(jC+x+r)\right)d\nu_j(x).
\]
Expand the logarithms in \(1/(jC)\) for sufficiently large \(C\).
The pure power of \(x\) at each order cancels between the two
sums. All remaining terms through order \(p\) integrate to zero.
At order \(p+1\) only the term in \(r x^p\) contributes, so the
log ratio equals
\[
     \frac{(-1)^p M_p}{(jC)^{p+1}}\,\frac{k(k-j)}2
                              +O(C^{-p-2}).                     \tag{3}
\]
The expansion is convergent for large \(C\) on the entire finite
support. The leading term is nonzero at every interior index and
has the same sign for \(k\) and \(j-k\). Choose \(C\) large enough
for the finitely many targets and indices simultaneously.
The physical coefficients divided by their constant coefficient
have the analogous expansion without the second sum in the log
ratio. Its leading term is
\[
                  \frac{(-1)^{p-1} k M_p}{p(jC)^p}
                               +O(C^{-p-1})\qquad(1\le k\le j).
\]
Thus after the equal-row padding the raw target polynomials
differ at every nonconstant coefficient, with the same
constant coefficient. The physical endpoint scales for these
increasing profiles differ at each target as well.

To make both profiles mirrored, append the reflection of every
row to its own profile. The observed equality then holds in both
orientations. The target log ratio at \(k\) is the sum of (3) at
\(k\) and \(j-k\), with twice the same leading term.
The two monomials now have equal total degrees, with the common
degree also the same in both profiles. The normalized mirrored
moments are monic with all roots on the unit circle, by the
[mirrored circle theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md).
Distinct normalized target polynomials therefore have different
multisets of root arguments.

## An exact two-observation replay

Use \(I=\{2,3\}\), \(J=\{4,5\}\) and baselines \(0,1,2\).
Taking \(D\) as the least common multiple of the cofactor
supports gives \(D=151200\). After common-row cancellation the
profiles have 19 rows each and totals \(23,267826\) before
translation. The following replay computes the signed target
marginals and tests both normalized and physical interior
coefficient ratios by integer powers. The normalized ratios
raised to the order \(j\) are rational, and the physical
ratios are rational themselves.

~~~python
from collections import Counter
from fractions import Fraction as Q
from math import lcm

baselines = (0, 1, 2)
observed = (2, 3)

def block(order, base):
    return Counter(range(order * base + 1, order * (base + 1) + 1))

def add(left, right, sign=1):
    result = Counter(left)
    for key, value in right.items():
        result[key] += sign * value
    return Counter({key: value for key, value in result.items() if value})

def multiply(left, right):
    result = Counter()
    for x, a in left.items():
        for y, b in right.items():
            result[x * y] += a * b
    return Counter({key: value for key, value in result.items() if value})

cofactors = []
for i in range(3):
    a, b = [base for base in baselines if base != i]
    minor = add(multiply(block(2, a), block(3, b)),
                multiply(block(2, b), block(3, a)), -1)
    cofactors.append(Counter({q: (-1)**i * c for q, c in minor.items()}))
D = lcm(*(q for cofactor in cofactors for q in cofactor))
assert D == 151200
generators = [(base, D // q, count)
              for base, cofactor in zip(baselines, cofactors, strict=True)
              for q, count in cofactor.items()]
rows = Counter()
for base, gap, count in generators:
    rows[base, gap] += count
    rows[base + 1, gap] -= count
assert sum((+rows).values()) == sum((-rows).values()) == 19
for side in (0, 1):
    totals = [sum((base + side * gap) * count
                  for (base, gap), count in part.items())
              for part in (+rows, -rows)]
    assert totals == [[23, 267826][side]] * 2

expected = {4: Q(-1, 476280000), 5: Q(-1, 127008000)}
checks = 0
for order in (2, 3, 4, 5):
    marginal = Counter()
    for base, gap, count in generators:
        for h in range(1, order + 1):
            x = Q(order * base + h, gap)
            marginal[x] += count
            marginal[x + 1] -= count
    if order in observed:
        assert all(value == 0 for value in marginal.values())
    else:
        moments = [sum(value * x**r for x, value in marginal.items())
                   for r in range(4)]
        assert moments == [0, 0, 0, expected[order]]
    for translation in (1, 5):
        for k in range(1, order):
            normalized_ratio = Q(1)
            physical_ratio = Q(1)
            for base, gap, count in generators:
                for h in range(1, order + 1):
                    a = order * (base + translation * gap) + h
                    normalized_ratio *= Q(
                        a**(order - k) * (a + gap * order)**k,
                        (a + gap * k)**order,
                    )**count
                    physical_ratio *= Q(a, a + gap * k)**count
            assert (normalized_ratio == 1) == (order in observed)
            assert (physical_ratio == 1) == (order in observed)
            checks += 2
assert checks == 40
~~~

The agreement here is on the separate multiset at each observed
order. The matched joint data keep the pairing of the Gauss slots
between two orders and recover the rows by the
[exact integer multiplicity theorem](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md).
The construction gives a limitation of separate marginals even
when every power sum is available, and supplies no common root
of consecutive moments of either ordinary profile.
