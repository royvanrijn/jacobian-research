# Ordinary binomial Strong Factorial through terminal twenty-four

For every \(2\le N\le24\), every \(a\ge0,d\ge1\), and every
characteristic-zero field, a nonzero ordinary binomial \(X^a(A+BX^d)\)
has a nonzero factorial moment in orders \(N-1,N\).
The [terminal-through-twelve theorem](STEPANOV_AUXILIARY_DETERMINANTS_SFC2.md)
covers the earlier terminals. The new range \(13\le N\le24\) uses
the same analytic complements and a complete exact finite replay.

## Complement for unbounded supports

Use Theorem thm:rational-annulus-cutoff in
[the initial and boundary argument](lower_exponent_initial_and_boundary.tex).
For a fixed terminal \(N\), let \(C_m\) be positive integers satisfying
\((C_m+1)^m<2C_m^m\) for \(m=N-1,N\). Set
\(B_m=((d_0+1)m)!/m!\). The strict inequality
\[
 B_N^{N-1}>(C_{N-1}C_N)^{N(N-1)}B_{N-1}^{N}
\]
gives coprimality for every \(a\ge1,d\ge d_0\), by the monotone
center ratio in that proof. The checker chooses the least integer
\(C_m\) satisfying its inequality and verifies the center comparison
for the explicit \(d_0\) in the table below.

The [sharp shape separation](GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md#3-shape-separation-and-a-sharper-soft-complement)
gives coprimality for \(a\ge d(N-2)\), including equality.
For \(a=0\) the [large-gap squarefreeness theorem](ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex)
gives coprimality once \(2^d\ge4N^2\) using the exact fixed-shape Appell
identity. These complements apply before any finite modular check.

## Complete finite replay

The exact enumeration uses
\(1\le d<d_0,\ 1\le a<d(N-2)\), together with
\(a=0,\ 1\le d<d_{\mathrm{const}}\), where \(d_{\mathrm{const}}\)
is the least integer with \(2^{d_{\mathrm{const}}}\ge4N^2\).
The count at each \(N\) is
\((N-2)d_0(d_0-1)/2-(d_0-1)+d_{\mathrm{const}}-1\).

| \(N\) | \(d_0\) | \(d_{\mathrm{const}}\) | Pairs checked |
|---:|---:|---:|---:|
| 13 | 73 | 10 | 28,845 |
| 14 | 81 | 10 | 38,809 |
| 15 | 89 | 10 | 50,829 |
| 16 | 97 | 10 | 65,097 |
| 17 | 105 | 11 | 81,806 |
| 18 | 114 | 11 | 102,953 |
| 19 | 122 | 11 | 125,366 |
| 20 | 131 | 11 | 153,150 |
| 21 | 139 | 11 | 182,101 |
| 22 | 148 | 11 | 217,423 |
| 23 | 157 | 12 | 257,021 |
| 24 | 166 | 12 | 301,136 |

For every one of the \(1{,}604{,}536\) triples, the
[checker](verify_ordinary_binomial_terminals_13_24.py) constructs
both actual polynomials
\[
 P_{m,a,d}(t)=\sum_{k=0}^m
       \binom mk \frac{(am+dk)!}{(am)!}t^k,\qquad m=N-1,N.
\]
It checks exact degrees and a full polynomial Bézout identity over
\(\mathbb F_p\), using \(p=1000003\) first. Four pairs at that prime
have a linear gcd, so the corresponding identities use \(p=1000033\).
Both primes are verified by trial division and exceed the largest
factorial argument, \(91{,}056\), throughout the domain.
Thus good reduction makes each successful identity a certificate
of rational coprimality.

The [receipt](ordinary_binomial_terminals_13_24_verification.json)
records each completed terminal and all four fallbacks. From the
repository root, run the checker with Python 3.12.3 and
python-flint 0.9.0.
Rational Bézout transfers the result to every characteristic-zero field.
A nonzero monomial has every positive ordinary factorial moment nonzero.
