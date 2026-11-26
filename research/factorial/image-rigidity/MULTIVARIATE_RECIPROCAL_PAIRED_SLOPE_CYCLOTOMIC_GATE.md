# Reciprocal Newton slope and cyclotomic certificates for mirrored moments

Let \(P_n,P_N\) be the monic ordinary moment polynomials at adjacent
orders \(n,N=n+1\) of a mirrored exponent profile. Each unequal row
\((\ell,\ell+d)\) occurs with the same multiplicity as
\((\ell+d,\ell)\), with \(\ell\ge0,d\ge1\).
The [mirrored-root theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#mirrored-profiles-at-every-fixed-order)
places all roots on the unit circle. The even moment is strictly
positive on the real line, so a common root would be nonreal.
The following certificate uses only the rational coefficients, their
Newton polygons, and polynomial arithmetic over finite fields.
It has no curvature restriction.

## Paired accounting for a possible factor

Fix a prime \(p\) and take \(v_p(p)=1\). For a root valuation \(s\)
write \(s=a/b\) in lowest terms, \(b>0\). The corresponding Newton
edge has slope \(-s\) and width divisible by \(b\).
Form its stepped residual polynomial in \(\mathbb F_p[X]\)
by reading the leading units of coefficients on that edge, starting
with its leftmost point and using index steps \(b\). Set a missing
edge's residual polynomial equal to one. Let \(G_{p,s}\) be the
monic gcd of the two stepped residual polynomials.

Here is a necessary degree test for any possible monic irreducible
common factor \(H\) of degree \(D\). At every prime \(p\) there must
be a monic divisor \(U_p\) of \(G_{p,0}\) and a monic divisor \(V_{p,s}\)
of \(G_{p,s}\) for each positive shared valuation \(s\), satisfying
\[
 D=\deg U_p+\sum_{s>0}2b\,\deg V_{p,s}.
 \tag{1}
\]
The divisor \(U_p\) is reciprocal with constant coefficient one,
and its multiplicities at residues \(1,-1\) are even, with the two
residues counted once if they are equal. Divisor one is allowed
throughout. Consequently a positive valuation of denominator \(b\)
whose stepped residual divisor has degree \(f\) consumes \(2bf\)
degrees of the common factor. The feasible degree lists in (1) can
be intersected across primes and with any separate degree bound.

Indeed the roots of \(H\) are on the circle, and complex conjugation
acts on them as inversion. No root is \(1\) or \(-1\).
In any common algebraic closure of \(\mathbb Q_p\), inversion
therefore bijects the roots at \(s\) with those at \(-s\),
with multiplicities. Roots of valuation \(s\) form a Galois-stable
set over \(\mathbb Q_p\), so their monic local factor is defined
over \(\mathbb Q_p\).

For the stepped residue rule, put \(\pi^b=p\), substitute
\(z=\pi^a x\), and divide the polynomial by the smallest power of
\(\pi\) occurring in its coefficients. The residue is a monomial
in \(x\) times the edge's stepped residual polynomial in \(x^b\).
The spacing follows from equality of \(b v_p(c_k)+ak\) on the
edge, where \(c_k\) is the coefficient at index \(k\).
The normalized reductions multiply, so the residue for the local
factor at that valuation divides the two edge residues after
removing their monomials. The same normalization with coefficients
in \(\mathbb Q_p\) gives stepped residues over \(\mathbb F_p\).
Since substitution \(X=x^b\) reflects divisibility between
polynomials in \(X\), the stepped divisor must divide \(G_{p,s}\).
Its width is \(b\) times its degree. At \(s=0\) the local roots
all have unit valuation, occur in distinct inverse pairs, and have
product one. Their monic factor reduces to \(U_p\).
Each inverse pair reducing to one of \(1,-1\) contributes twice
to its multiplicity. This gives (1), without an assumption on
tameness or simple residual roots.
[The general stepped face comparison](MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md#exact-recognition-of-the-newton-and-residue-data)

## Cyclotomic reduction

Suppose that for a feasible degree \(D\), every allocation (1)
of that degree at every prime has \(\deg V_{p,s}=0\) for every
\(s>0\). Then any common irreducible factor of degree \(D\) is
cyclotomic. Indeed all its roots have valuation zero at every
prime. Thus its monic coefficients, which are rational, are
integral everywhere and belong to \(\mathbb Z\).
All conjugates have absolute value one. By Kronecker's theorem
such a root is a root of unity. One elementary proof takes the
monic polynomials whose roots are the \(k\)-th powers of all
conjugates. Their coefficients are bounded integers, so the
polynomials and the union of their root sets are finite.
The sequence of powers of the original root must repeat.

Only finitely many primes are needed for the nonzero-slope
hypothesis. Outside the denominator support of the coefficients
of either monic reciprocal polynomial, its roots are integral at
that prime. Inversion then gives valuation zero for every root.
For an ordinary profile, every denominator prime is bounded by
the largest factorial argument used in the two arrays.

There are finitely many candidate cyclotomic orders \(h\), since
\(\varphi(h)=D\) implies \(h\le2D^2\). Indeed \(p-1\ge\sqrt p\)
for odd primes, so the prime-power formula for \(\varphi\) gives
\(\varphi(h)\ge\sqrt{h/2}\), with only \(p=2\) contributing the
factor \(1/\sqrt2\). The orders \(1,2\) are excluded by real
positivity of the even moment. For each remaining candidate and
each prime write \(h=p^e h_0\), \(p\nmid h_0\).
It can survive only if
\[
                 \overline{\Phi_{h_0}}^{\,\varphi(p^e)}
                         \mid G_{p,0}.
 \tag{2}
\]
For \(e=0\) this is direct, and otherwise the formula follows
from \(\Phi_{p^e h_0}(z)
=\Phi_{h_0}(z^{p^e})/\Phi_{h_0}(z^{p^{e-1}})\) and Frobenius.
A candidate can also be checked directly by exact division by
\(\Phi_h\).

This gives an exact certificate for coprimality and hence angular
separation of the two mirrored moments. If applied at all
\(3\le N<J\), the
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
supplies the remaining high orders, and the quadratic has no real
root at terminal two. Sparsity of common nonzero slopes is a
sufficient condition for the cyclotomic reduction, not an asserted
property of all mirrored profiles.

The accounting applies to any two monic rational polynomials
with all roots on the unit circle and with disjoint real roots.
For reciprocal moments it can also be performed on the
half-degree trace polynomials. Remove \(z+1\) from an odd moment
and write the remaining reciprocal polynomial as \(z^r A(z+z^{-1})\).
Each root pair gives one real root \(x=\alpha+\alpha^{-1}\) of
\(A\), with \(-2<x<2\). If \(v_p(\alpha)\ne0\), then
\(v_p(x)=-|v_p(\alpha)|\), and otherwise \(x\) is integral at \(p\).
Thus common negative trace valuations account for the paired
nonzero slopes. A totally real algebraic integral trace in this
range again gives a cyclotomic factor, by applying the preceding
argument to the roots of \(z^2-xz+1\).

## Exact example with two different gaps

Take the four shared rows \((7,9),(9,7),(11,12),(12,11)\) and
the pair of orders \(17,18\). The monic coefficient is
\[
 [z^k]P_m=\binom mk
           \prod_{(\ell,d)\in\{(7,2),(11,1)\}}
          \frac{(m\ell+kd)!\,(m\ell+(m-k)d)!}
               {(m\ell)!\,(m(\ell+d))!}.
\]
At \(p=2\) both polygons have only slope zero. The zero-slope
residues are \(1+z+z^{16}+z^{17}\) and \(1+z^{18}\).
Their gcd is \(1+z^2\), allowing degree at most two.
Among the 47 primes at most 216, the only common nonzero
valuations occur at \(p=7\), and are \(\pm1/7\).
They would consume at least fourteen degrees in (1).
At \(7\), the stripped zero-slope residues are
\(4+2z+2z^2+4z^3\) and \(6+4z+6z^2\),
giving \(G_{7,0}=1+3z+z^2\).

Any common factor must therefore be one of the degree-two
cyclotomic polynomials \(\Phi_3,\Phi_4,\Phi_6\).
The gcd at \(2\) permits only \(\Phi_4\), while the gcd at \(7\)
excludes it. Each of these two primes leaves degree two in its residual
degree test. The reciprocal
valuation constraint and the cyclotomic comparison close the
pair without factoring the factorial polynomials over \(\mathbb Q\).

The valuation and residue claims can be reproduced with the
following exact code. It calculates only this fixed example.
Coefficients and lists of polynomial coefficients are stored in
increasing order.

```python
from fractions import Fraction
from math import comb, factorial, isqrt

ROWS = ((7, 2), (11, 1))

def coefficients(m):
    out = [Fraction(comb(m, k)) for k in range(m + 1)]
    for ell, d in ROWS:
        denominator = factorial(m * ell) * factorial(m * (ell + d))
        for k in range(m + 1):
            numerator = factorial(m * ell + k * d)
            numerator *= factorial(m * ell + (m - k) * d)
            out[k] *= Fraction(numerator, denominator)
    return out

def valuation_and_unit(value, p):
    def split(integer):
        valuation = 0
        while integer % p == 0:
            integer //= p
            valuation += 1
        return valuation, integer
    a, u = split(value.numerator)
    b, v = split(value.denominator)
    return a - b, u * pow(v, -1, p) % p

def polygon_and_zero_residue(sequence, p):
    data = [valuation_and_unit(value, p) for value in sequence]
    hull = []
    for k, (valuation, _) in enumerate(data):
        while len(hull) >= 2:
            i, vi = hull[-2]
            j, vj = hull[-1]
            if (vj - vi) * (k - j) < (valuation - vj) * (j - i):
                break
            hull.pop()
        hull.append((k, valuation))
    slopes = {}
    for (i, vi), (j, vj) in zip(hull, hull[1:]):
        slopes[Fraction(vi - vj, j - i)] = (i, j)
    if 0 not in slopes:
        return set(slopes), [1]
    i, j = slopes[0]
    height = data[i][0]
    return set(slopes), [u if v == height else 0 for v, u in data[i:j+1]]

def monic_gcd(f, g, p):
    f, g = list(f), list(g)
    while g:
        inverse = pow(g[-1], -1, p)
        while len(f) >= len(g):
            index = len(f) - len(g)
            ratio = f[-1] * inverse % p
            for j, coefficient in enumerate(g):
                f[index + j] = (f[index + j] - ratio * coefficient) % p
            while f and f[-1] == 0:
                f.pop()
        f, g = g, f
    inverse = pow(f[-1], -1, p)
    return [c * inverse % p for c in f]

lower, upper = coefficients(17), coefficients(18)
primes = [p for p in range(2, 217)
          if all(p % j for j in range(2, isqrt(p) + 1))]
assert len(primes) == 47
nonzero = {}
for p in primes:
    slopes_lower, g = polygon_and_zero_residue(lower, p)
    slopes_upper, h = polygon_and_zero_residue(upper, p)
    shared = (slopes_lower & slopes_upper) - {0}
    if shared:
        nonzero[p] = shared
    if p == 2:
        assert slopes_lower == slopes_upper == {0}
        assert monic_gcd(g, h, p) == [1, 0, 1]
    if p == 7:
        assert g == [4, 2, 2, 4] and h == [6, 4, 6]
        assert monic_gcd(g, h, p) == [1, 3, 1]
assert nonzero == {7: {Fraction(-1, 7), Fraction(1, 7)}}
```
