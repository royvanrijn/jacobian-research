# Later Euclidean remainders for the linked beta pair

For every nonnegative integer exponent pair \(r,u\), the third
Euclidean remainder starting from \(J_{n+1},J_n\) has degree
\(n-3\) and positive leading coefficient for \(n\ge3\). The
fourth has degree \(n-4\) and positive leading coefficient for
\(n\ge4\). Each division uses the monic versions of the
previous two polynomials and the convention
\(\mathrm{dividend}=\mathrm{quotient}\cdot\mathrm{divisor}
+\mathrm{remainder}\).
These two additional steps advance
`GAP-FACTORIAL-LINKED-BETA-REGULAR-CHAIN`.

All Euclidean remainders down to the constant have strictly
positive leading coefficients for starts \(n=1,\ldots,9\).
Separate positive resultant certificates at starts \(9,10\)
complete the linked coprimality prefix through ten, beyond
`FACTORIAL-BINARY-SHIFTED-RESULTANT-POSITIVE-WINDOWS-THROUGH-FIVE`.

## 1. Exact recurrence for a short coefficient list

Recall the [physical monic linked moments](canonical-density/moment_pair_obstructions.md)
\[
 J_m(z;r,u)=\sum_{k=0}^m(-1)^k\binom{m}{k}
       \frac{(mr+1)_k}{(m(r+u)+2)_k}z^{m-k}.
\]
Let \(A_i\) and \(B_i\) be descending coefficient lists of two
polynomials of consecutive degrees, with the larger one listed
first. Each may have its own positive scale; in particular
\(A_0,B_0>0\). The actual monic polynomials use the lists
\(A_i/A_0\) and \(B_i/B_0\). The quotient of these polynomials is
\(z+A_1/A_0-B_1/B_0\). Thus the coefficient in slot \(j\) of the
expected monic-divisor remainder is
\[
  \frac{C_j}{A_0 B_0^2},\qquad
  C_j=A_{j+2}B_0^2-A_0B_{j+2}B_0-A_1 B_{j+1} B_0
                                      +A_0 B_1 B_{j+1}.             \tag{1}
\]
Slot zero denotes the expected leading coefficient, whether or
not that coefficient vanishes. Its strict positivity proves
the expected degree and supplies the next positive scale.

Every step uses two fewer descending coefficients. For \(h\)
steps, it suffices to begin with the first \(2h+1\) coefficients
of each polynomial, inserting zeros beyond its degree.
Write \(s=r+u\). A positive common denominator for these lists is
\[
 D=(2h)!\prod_{i=2}^{2h+1}\bigl(ns+i\bigr)\bigl((n+1)s+i\bigr).
\]
The exact scaled lists follow directly from the displayed formula
for \(J_m\). In particular, the falling factorial in \(m\) inserts
the required zero coefficients even when \(2h>m\).

We use (1) over the ring of integer polynomials in the parameter
charts. At each step we may divide all the \(C_j\) by their common
gcd with \(A_0B_0^2\) and keep the shortened coefficient list.
The common divisor is normalized to have positive constant
coefficient in the chart. The denominator \(A_0B_0^2\) is
positive throughout the chart by the earlier leading-coefficient
bounds. Since the common divisor divides it, the divisor has
no zero in the nonnegative real orthant and has positive sign
there by connectedness and normalization at the origin. In
particular the signs after cancellation still certify actual
leading signs. Equivalently the new list and the shortened
denominator retain the exact unnormalized remainder
coefficients. The argument does not assume in advance that a
degree persists.

By reflection it suffices to use \(r\le u\). In fact
\(J_m(z;u,r)=(-1)^mJ_m(1-z;r,u)\), and for any consecutive
monic dividend and divisor this operation multiplies its
remainder by the same sign as its dividend. When the
remainder has the expected degree, its leading coefficient
therefore keeps its sign. This permits induction over the
steps. We cover every sorted nonnegative integer profile
with the charts
\[
                       (r,u)=(0,v),\qquad
                       (r,u)=(1+y,1+y+v),
                  \qquad y,v\ge0.
\]
The positivity statements also hold for the real exponents
in these two charts.

## 2. All-order third and fourth pivots

The [first two leading-coefficient theorems](LINKED_BETA_UNEQUAL_PROFILE_PIVOTS.md)
supply the induction through step two. Applying (1) with
\(n=h+x,\ x\ge0\), the exact reduced leading numerators at
steps \(h=3,4\) have the following support and sign counts.
All listed nonzero coefficients are positive, including
the constant, so the polynomials are strictly positive
throughout their respective charts.

| Step \(h\) | Profile chart | Multidegrees in \((x,y,v)\) | Positive coefficients | Negative coefficients |
|---:|---|---|---:|---:|
| 3 | Axis | \((24,0,26)\) | 465 | 0 |
| 3 | Positive | \((34,33,26)\) | 19,845 | 0 |
| 4 | Axis | \((46,0,48)\) | 1,562 | 0 |
| 4 | Positive | \((62,61,48)\) | 117,306 | 0 |

The third-step result supplies the prior positive leading
coefficient when applying the fourth-step calculation.
For start \(n=4\), step four itself gives a nonzero constant,
so the whole Euclidean chain holds at that start independently
of the previously known resultant.

## 3. Finite extension of the actual coprimality window

At fixed starts \(n=5,6,7,8\), the same recurrence gives
positive leading results through the constant at every step.
The table records the number of positive coefficients for
each pivot in each chart. All counts follow coefficient
expansion, including at the last step of degree zero in the
polynomial variable. The remaining entries still depend on
the exponent parameters.

| Start | Axis counts in step order | Positive chart counts in step order |
|---:|---|---|
| 5 | \(10,14,24,36,45\) | \(65,147,468,1134,1980\) |
| 6 | \(12,16,27,42,57,67\) | \(90,184,567,1449,2850,4355\) |
| 7 | \(14,18,29,46,66,83,93\) | \(119,225,638,1679,3597,6059,8370\) |
| 8 | \(16,20,31,49,72,95,113,123\) | \(152,270,713,1862,4140,7505,11300,14637\) |

There is a separate coprimality certificate at start \(9\).
For the sorted positive chart define
\[
 G_m(z)=\sum_{k=0}^m\binom mk(mr+1)_{m-k}(mu+1)_kz^k.
\]
The exact change of chart is
\[
 J_m(z)=\frac{(z-1)^m}{(m(r+u)+2)_m}
                    G_m\!\left(\frac{z}{z-1}\right).
\]
It follows by the binomial expansion of the beta integral.
The projective endpoints give pure monomials and
have nonzero moments. The exact Sylvester resultant
\(\operatorname{Res}_z(G_9,G_{10})\) has 8,281 strictly
positive nonzero coefficients, including the constant,
in \(y=r-1,v=u-r\), and multidegrees \((135,90)\).
At start \(10\) the corresponding counts are 12,321
positive coefficients and multidegrees \((165,110)\).
The [all-order axis theorem](GAMMA_QUOTIENT_TRANSFER_AUDIT.md)
supplies the other chart at these orders. This proves
shifted coprimality at starts \(6,7,8,9,10\) for
every nonnegative integer profile. The all-order linked
coprimality and later Euclidean steps remain open.

For start \(9\), the first eight Euclidean remainders have
positive leading coefficients by the same expansion. Their
coefficient counts in step order are
\[
 (18,22,33,51,76,104,129,147)\quad\text{and}\quad
 (189,319,792,1989,4522,8684,13932,19257)
\]
for the axis and positive charts. The eighth remainder is
linear, and the coprimality at start \(9\) makes the final
remainder a nonzero constant. It is positive: the resultant
of the monic \(J_{10},J_9\) is positive, since \(J_{10}\)
is positive at every real root of \(J_9\) and products at
nonreal conjugate pairs are positive squares. The [regular
chain resultant product](canonical-density/moment_pair_obstructions.md),
equation (2), therefore makes the final coefficient positive
given the first eight signs.


## 4. Exact replay

The script below uses Python and python-flint 0.9.0. It
checks every coefficient of the listed integer polynomials.
The largest all-order coefficient check is the positive
chart for the fourth step. Fixed-order checks are smaller
because \(x\) is absent. No numerical sign estimate is used.

```python
from math import comb, factorial, prod
from flint import fmpz_mpoly_ctx

ring = fmpz_mpoly_ctx.get(("x", "y", "v"), ordering="lex")
x, y, v = ring.gens()

def leading_numerators(n, r, u, steps, take=None):
    top = 2 * steps
    s = r + u
    scale = factorial(top)
    lower = [n*s + j for j in range(2, top+2)]
    upper = [(n+1)*s + j for j in range(2, top+2)]

    def scaled(m, first, second):
        return [
            (-1)**k * (scale//factorial(k)) * prod(first[k:]) * prod(second)
            * prod((m-j)*(m*r+1+j) for j in range(k))
            for k in range(top+1)
        ]

    A = scaled(n+1, upper, lower)
    B = scaled(n, lower, upper)
    answer = []
    for _ in range(steps if take is None else take):
        a, b = A[0], B[0]
        den = a*b**2
        C = [
            A[j+2]*b**2 - a*B[j+2]*b
            - A[1]*B[j+1]*b + a*B[1]*B[j+1]
            for j in range(len(B)-2)
        ]
        common = den
        for coefficient in C:
            common = common.gcd(coefficient)
        if common[0,0,0] < 0:
            common = -common
        assert common[0,0,0] > 0
        C = [coefficient//common for coefficient in C]
        answer.append(C[0])
        A, B = B, C
    return answer

def positive(poly, count):
    assert len(poly.coeffs()) == count
    assert all(c > 0 for c in poly.coeffs())
    assert poly[0,0,0] > 0

for h, r, u, count, degrees in (
    (3, 0*y, v, 465, (24,0,26)),
    (3, y+1, y+1+v, 19845, (34,33,26)),
    (4, 0*y, v, 1562, (46,0,48)),
    (4, y+1, y+1+v, 117306, (62,61,48)),
):
    lead = leading_numerators(x+h, r, u, h)[-1]
    positive(lead, count)
    assert lead.degrees() == degrees

for n, axis_counts, positive_counts in (
    (5, (10,14,24,36,45), (65,147,468,1134,1980)),
    (6, (12,16,27,42,57,67), (90,184,567,1449,2850,4355)),
    (7, (14,18,29,46,66,83,93), (119,225,638,1679,3597,6059,8370)),
    (8, (16,20,31,49,72,95,113,123),
        (152,270,713,1862,4140,7505,11300,14637)),
):
    for r, u, counts in ((0*y, v, axis_counts),
                          (1+y, 1+y+v, positive_counts)):
        leads = leading_numerators(n+0*x, r, u, n)
        for lead, count in zip(leads, counts, strict=True):
            positive(lead, count)

for r, u, counts in (
    (0*y, v, (18,22,33,51,76,104,129,147)),
    (1+y, 1+y+v, (189,319,792,1989,4522,8684,13932,19257)),
):
    leads = leading_numerators(9+0*x, r, u, 9, take=8)
    for lead, count in zip(leads, counts, strict=True):
        positive(lead, count)


# Eliminate the polynomial variable x, here used as z.
r, u = 1+y, 1+y+v
def G(m):
    return sum(
        comb(m,k)*prod(m*r+j for j in range(1,m-k+1))
        *prod(m*u+j for j in range(1,k+1))*x**k
        for k in range(m+1)
    )

for n, count, degrees in ((9, 8281, (0,135,90)),
                          (10, 12321, (0,165,110))):
    res = G(n).resultant(G(n+1), 0)
    positive(res, count)
    assert res.degrees() == degrees
print("Exact pivot and resultant certificates passed.")
```

Independent rational Euclidean division of the full polynomials
checked the numerical third-pivot identity at 60 exact chart
points and the signs of all pivots on 245 fixed-profile samples
counted by step. An independently formed physical resultant at
start \(8\) in the sorted positive chart also had 5,329 positive
coefficients and no negatives.
