# First-digit exclusion of the three prime-front zero-rank residue roots

Let \(f=ax^2+bxy+cy^2+k\) with \(abck\ne0\), and let
\(M_j=\Lambda_2(f^j)\) with \(\Lambda_2(x^i y^j)=i!j!\).
In the
[prime-front zero-rank gate](BINARY_EVEN_QUADRATIC_PRIME_FRONT_ZERO_RANK_GATE.md),
let \(\ell=1\) denote the window \(p-1,p,p+1,p+2\)
and \(\ell=2\) the window \(p,p+1,p+2,p+3\).
After normalizing \(a=1\), call \(X\) the
residue of \((c-(b/2)^2)/k\).
The three residue roots
\[
             (p,\ell,X)=(11,2,-6),\ (19,1,-2),\ (37,2,-9)
\]
cannot lift to a four-zero window in characteristic zero.
The first window, \(11,\ldots,14\), is independently covered
for all coefficient residues by the
[prime-penultimate theorem at \(13\)](BINARY_EVEN_QUADRATIC_PRIME_PENULTIMATE_FOUR_WINDOW.md).
The other two are \(18,\ldots,21\) and \(37,\ldots,40\).

Each exclusion uses three original moments at the selected
orders and one additional valuation digit beyond the residue
test. The [exact integer replay](verify_binary_even_quadratic_prime_front_zero_lifts.py)
computes the needed coefficients directly from the native
multinomial expansion.

## The local scaling and its precision

Work over the algebraic \(p\)-adic closure with \(v(p)=1\).
Put
\[
 h=(p-1)/2,\quad \kappa=k^p/p^2,\quad
 \alpha=p/k^h,\quad x=(c-(b/2)^2)/k,\quad z=(b/2)/\alpha .
\]
The zero-rank gate gives \(v(\alpha)=1/p\),
\(v(\kappa+2)\ge(2\ell+3)/p\), and
\(\bar z=-2(-1)^h g_{h+1}(X)\) with the
[factorial polynomial](BINARY_EVEN_QUADRATIC_PRIME_FRONT_ZERO_RANK_GATE.md)
\(g_m(X)=\sum_{i=0}^m\binom mi(2i)!X^i\).
We have
\[
 p=\kappa^h\alpha^p,\qquad
 k=\kappa\alpha^2,\qquad b=2\alpha z,\qquad
                 c=\alpha^2(z^2+\kappa x).
\]
Both \(x,z\) are units.

For nonnegative \(I,B,J,K\) with \(I+B+J+K=n\),
let
\[
 N=\frac{n!}{I!B!J!K!}(2I+B)!(2J+B)!,
       \qquad w=v_p(N),\qquad m_0=B+2(J+K)+pw .
\]
Its contribution to \(M_n\) is
\[
 \alpha^{m_0}(N/p^w)\,2^B\kappa^{K+hw}
                      z^B(z^2+\kappa x)^J .                 \tag{1}
\]
Fix a target precision \(m+1\) in powers of \(\alpha\), and
let \(m_{\min}\) be the least \(m_0\) of the raw terms.
If
\[
                  m+1<m_{\min}+\min(p,2\ell+3),             \tag{2}
\]
coefficient integers can be reduced modulo \(p\) and
\(\kappa\) can be replaced by \(-2\) through that precision.
The two replacements gain valuation at least \(1\)
and \((2\ell+3)/p\), respectively, from each originating
term, so (2) accounts also for any cancellation between
terms of lower degree. Terms with \(m_0>m+1\) are higher
without any replacement.

Expanding (1) under these substitutions yields two first
nonzero layers \(E(x,z),F(x,z)\) at degrees \(m,m+1\)
over \(\mathbb F_p\) in each of the orders listed below.
The checker verifies that the intervening lower
coefficients vanish as polynomials over \(\mathbb F_p\)
and that (2) holds. The normalized equation for integer
coefficient lifts of those layers is
\[
       M_n/\alpha^m = E(x,z)+\alpha F(x,z)+\varepsilon,
                     \qquad v(\varepsilon)>1/p .           \tag{3}
\]

## The inconsistent first digit

For each row below the entries are the order, its
normalizing degree \(m\), and
\((E_x,E_z,F)\) evaluated at the specified \((\bar x,\bar z)\)
in the residue field. The value of \(E\) is zero in
every row.

| \(p\) | \(\ell\) | \((\bar x,\bar z)\) | \(n\) | \(m\) | \((E_x,E_z,F)\) |
|---|---|---|---|---|---|
| 11 | 2 | \((5,9)\) | 11 | 23 | \((1,2,2)\) |
| 11 | 2 | \((5,9)\) | 12 | 25 | \((9,1,4)\) |
| 11 | 2 | \((5,9)\) | 13 | 27 | \((3,4,1)\) |
| 19 | 1 | \((17,17)\) | 18 | 18 | \((14,0,16)\) |
| 19 | 1 | \((17,17)\) | 19 | 39 | \((18,2,5)\) |
| 19 | 1 | \((17,17)\) | 20 | 41 | \((11,16,3)\) |
| 37 | 2 | \((28,6)\) | 37 | 75 | \((18,2,20)\) |
| 37 | 2 | \((28,6)\) | 38 | 77 | \((5,15,2)\) |
| 37 | 2 | \((28,6)\) | 39 | 79 | \((34,28,20)\) |

To see why the first digit is required even in a ramified
extension, choose integer representatives \(x_0,z_0\) of the
displayed residues. The first two leading gradients have
nonzero determinant modulo \(p\). If the smaller valuation of
\(x-x_0,z-z_0\) were strictly below \(1/p\), the invertible
linear part of their first two equations (3) would give
that same valuation to at least one equation. The constants
\(E(x_0,z_0)\) are divisible by \(p\), the quadratic
Taylor remainders have strictly higher valuation than the
linear part, and the terms involving \(\alpha\) have
valuation at least \(1/p\). This is impossible.
Hence both corrections have valuation at least \(1/p\).

Put \((\xi,\zeta)=(\overline{(x-x_0)/\alpha},
\overline{(z-z_0)/\alpha})\). Dividing (3) by \(\alpha\)
and reducing gives \(E_x\xi+E_z\zeta+F=0\) in each row.
For \(p=11,19,37\), respectively, the first two rows
force \((\xi,\zeta)=(10,5),(7,1),(0,27)\).
The remaining row then gives the nonzero residual
\(7,1,36\), respectively. This proves the exclusions.

The uniform first corrections give a separate check of the two
windows with no prime-penultimate overlap. At \(p=19\), the
[second-position correction](BINARY_EVEN_QUADRATIC_PRIME_SECOND_ZERO_RANK_LIFT.md#the-first-digit-fixes-two-rational-residues)
requires \(X=359/300\equiv10\pmod {19}\), contradicting \(X=17\).
At \(p=37\), the
[first-position correction](BINARY_EVEN_QUADRATIC_PRIME_FIRST_ZERO_RANK_LIFT.md#boundary-data-and-the-first-correction)
has \(C=1,\ G=-\bar z/2=-3\) and hence \(W=X-2G^2=10\).
Its two quadratics give \(Q_1(W)=22\) and \(Q_2(W)=18\) modulo
\(37\), contradicting the implication \(Q_1(W)\ne0\Rightarrow
Q_2(W)=0\).

Combined with the prime-front residue results, the native
calculations remove every surviving zero-rank root in the stated
finite gcd ranges. For other primes the general finite-field
gate and the first-digit procedure remain available.
