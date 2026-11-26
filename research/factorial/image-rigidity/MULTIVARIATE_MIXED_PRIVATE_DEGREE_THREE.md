# All-order closure with private exponent totals one and at least three

Suppose an ordinary two-monomial profile contains a private
variable on each side, and one side has total private exponent
at least three. There is no restriction on the number, orientation
or size of the other rows. In the [matched tanh grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md),
a slot is private exactly when
\(\bigl(n/T'-(n-1)/T\bigr)^{-1}=1\).
The count of such slots in each orientation equals its total
private exponent. Thus the criterion requires at least one
private slot in each orientation and at least three in one.

**Theorem.** Consecutive moment root arguments are disjoint whenever
the lower order is at least twelve. This conclusion also allows a
common nonnegative Gaussian coefficient strength. For ordinary
profiles the [known smaller terminals](MULTIVARIATE_TERMINAL_TWELVE_ALL_PROFILES.md)
then give all-degree SFC(2) over every characteristic-zero field.

## Two extra factorial factors

Orient the side with private exponent total at least three toward
the coefficient index \(k\). Normalize the constant and leading
coefficients in each degree \(m\) to one with a positive scalar
and positive change of variable. A private exponent \(b\) on the
chosen side supplies the factors \((h/b)_k\), \(1\le h\le b\),
up to a geometric power in \(k\). This follows by grouping the
factors in \((bk)!\) by residue modulo \(b\). A private exponent
on the opposite side supplies the analogous factors in \(m-k\).

Reserve one shape-one factor from a private row on each side.
These are \(k!\) and \((m-k)!\) and cancel the binomial factor
after endpoint normalization. At least two private factors with
shapes \(0<a\le1\) remain on the chosen side. The normalized
factor \((a)_k/(a)_m^{k/m}\) is at most \(k!/(m!)^{k/m}\).
Indeed the function \(\log\Gamma(a+x)-\log\Gamma(1+x)\) is
convex for \(x\ge0\) by the trigamma series. Subtracting its chord
on \([0,m]\) gives the comparison. Every other factorial factor
has normalized value at most one by log convexity. Thus the
actual normalized coefficients satisfy
\[
           0<c_{m,k}\le f_{m,k}^2,\qquad
           f_{m,k}=\frac{k!}{(m!)^{k/m}} .
\]
A common Gaussian multiplier \(e^{-\tau k(m-k)/(2m)}\)
with \(\tau\ge0\) respects the same bound.

## Uniform bounds for the two halves

For any fixed integer \(k\ge1\), both \(f_{m,k}\) and \(f_{m,m-k}\)
decrease with \(m\ge k\). For the first, the geometric mean
\((m!)^{1/m}\) increases. For the second use real \(m>k-1\) and
write
\[
 f_{m,m-k}=\frac{\Gamma(m+1)^{k/m}}{\prod_{h=0}^{k-1}(m-h)}.
\]
Its logarithmic derivative is
\[
 \frac{k}{m^2}\int_0^m s\psi'(1+s)\,ds
          -\sum_{h=0}^{k-1}\frac1{m-h}<0 .
\]
Here \(\psi'(1+s)=\sum_{j\ge1}(s+j)^{-2}<1/s\) for \(s>0\),
by comparison with an integral. This is the exact normalized
Gamma-factor monotonicity in the two orientations.

Split \(P_m(z)=L_m(z)+z^m U_m(z^{-1})\), where each half includes
the endpoint coefficient one and half of the middle coefficient
when \(m\) is even. For \(m\ge12\) the bounds at \(m=12\) suffice
for the first five coefficients, in either half. For \(k\ge6\)
the largest bound in either half is at \(m=2k\) and is
\[
                       f_{2k,k}^2=\binom{2k}{k}^{-1}<3^{-k}.
\]
The sequence \(\binom{2k}{k}^{1/k}\) is increasing, since the
ratios of successive central binomial coefficients are
\(4-2/k\) for steps from \(k-1\) to \(k\), and
\(\binom{12}{6}=924>3^6\).

Put \(R=25/24\). For a half \(H=1+\sum_{k\ge1}h_kz^k\)
let \(A=\sum h_kR^k\), \(B=\sum kh_kR^k\), and
\(C=\sum_{k\ge2}kh_kR^k\). The displayed coefficient bounds give
\[
\begin{array}{c|ccc}
 & A & B & C\\ \hline
 L_m & 1/20 & 83/1000 & 45/1000\\
 U_m & 69/250 & 393/1000 & 191/1000
\end{array}
\]
as strict upper bounds. They imply
\(R^{12}>(1+A_L)/(1-A_U)\) and
\(R^{12}>(1+A_U)/(1-A_L)\). Rouché on \(|z|=R\) and on
the reversed polynomial shows that all roots lie in
\(R^{-1}<|z|<R\). For example, on \(|z|=R\),
the term \(z^m U_m(z^{-1})\) has modulus at least
\(R^m(1-A_U)>1+A_L\) and has all its \(m\) zeros inside,
since \(U_m\) is nonzero on \(|z|\le R\).

## Subtracting the adjacent phase equations

The logarithm of each half with value zero at the origin is
analytic on \(|z|\le R\). For \(0<\theta\le\pi\) and
\(0<s\le R\), the first coefficient has nonnegative contribution
to \(\Im H(s e^{i\theta})\). All other terms have absolute
imaginary part at most \(kh_ks^k\theta\). Since
\(\Re H(s e^{i\theta})\ge1-A>0\) and
\(|\arctan x|\le |x|\), we have
\[
        -\frac{C}{1-A}\theta
        \ \le\ \Im\log H(s e^{i\theta})
        \ \le\ \frac{B}{1-A}\theta .
\]
At a root \(r e^{i\theta}\), the split identity gives
\[
    m\theta\equiv\pi+\phi_m(r,\theta)\pmod{2\pi},\qquad
    \phi_m(r,\theta)=\Im\log L_m(r e^{i\theta})
                  +\Im\log U_m(r^{-1}e^{i\theta}).
\]
The table gives \(-8\theta/25<\phi_m(r,\theta)<16\theta/25\).
Two adjacent roots on one ray would therefore satisfy
\(\theta\equiv\phi_{m+1}-\phi_m\pmod{2\pi}\) with
\(|\phi_{m+1}-\phi_m|<24\theta/25\). This is impossible
for \(0<\theta\le\pi\). Real coefficients give the lower
half-plane, and positive coefficients exclude \(\theta=0\).

## Exact scalar replay

~~~python
from fractions import Fraction as F
from math import comb, factorial

radius = F(25, 24)
ratio = radius / 3
tail_mass = ratio**6 / (1 - ratio)
tail_slope = tail_mass * (6 + ratio / (1 - ratio))
assert comb(12, 6) > 3**6

caps = (
    (358, 52, 17, 10, 9),
    (1943, 450, 126, 44, 19),
)
bounds = (
    (F(1, 20), F(83, 1000), F(45, 1000)),
    (F(69, 250), F(393, 1000), F(191, 1000)),
)
for side, (row, (mass_bound, slope_bound, higher_slope_bound)) in enumerate(
    zip(caps, bounds)
):
    scaled = []
    for k, numerator in enumerate(row, 1):
        index = k if side == 0 else 12 - k
        cap = F(numerator, 10_000)
        assert cap**6 * factorial(12)**index > factorial(index)**12
        scaled.append(cap * radius**k)
    assert sum(scaled) + tail_mass < mass_bound
    assert sum(k * x for k, x in enumerate(scaled, 1)) + tail_slope < slope_bound
    assert sum(k * x for k, x in enumerate(scaled, 1) if k > 1) + tail_slope < higher_slope_bound

left, right = bounds
assert radius**12 > (1 + left[0]) / (1 - right[0])
assert radius**12 > (1 + right[0]) / (1 - left[0])
assert sum(b / (1 - a) for a, b, c in bounds) < F(16, 25)
assert sum(c / (1 - a) for a, b, c in bounds) < F(8, 25)
~~~
