# Effective factorial moment windows with prescribed denominator primes

Let \(f\) be a polynomial in any finite number of variables with
algebraic coefficients and nonzero constant coefficient. Scaling
by its constant coefficient does not affect the zero indices, so
normalize the constant to one. Let \(S\) be any finite set of
rational primes containing the primes below all denominator
places of the normalized coefficients. For instance a common
positive integer \(D\) such that \(Df\) has algebraic-integer
coefficients is enough, using all prime factors of \(D\).

**Theorem.** Every positive index of a zero factorial moment
is \(S\)-smooth. For every prescribed \(S\) there is an effectively
computable finite list of adjacent pairs that contains an initial
pair of every possible zero window of length at least two, uniformly
in the dimension, degree, support, and coefficients with that
denominator-prime constraint. Thus for any fixed algebraic \(f\)
with nonzero constant and any requested window length at least two,
the assertion that every window of that length has a nonzero
factorial moment is decidable by exact algebraic arithmetic.

## A prime dividing the moment index

Write \(\mathcal L(x_1^{a_1}\cdots x_d^{a_d})=\prod_i a_i!\).
If \(p\notin S\) and \(p\mid j\), use any place over \(p\) of
the coefficient number field. All coefficients are integral
at this place. In its residue field,
\[
                \overline{f^j}=(\overline{f^{j/p}})^p .
\]
Every nonconstant monomial in this power has some exponent at
least \(p\), and its factorial multiplier is divisible by \(p\).
The constant of the power is one. Thus
\(\mathcal L(f^j)\) reduces to one and cannot vanish. In fact
at every place over a prime dividing a zero index, at least
one normalized coefficient is nonintegral.

## The explicit adjacent-index list

If \(2\notin S\), no adjacent pair can be \(S\)-smooth.
Otherwise apply the Størmer–Lehmer Pell reduction described below.
For clarity the recipe here only uses neighboring integers.
If \(n,n+1\) are both \(S\)-smooth, write \(n(n+1)=d y^2\)
with \(d>1\) squarefree and supported on \(S\). The value
\(d=1\) is impossible, because the two relatively prime
positive integers \(n,n+1\) would both be squares. Put
\(X=2n+1,\ Y=2y\), giving \(X^2-dY^2=1\), with \(Y\)
also \(S\)-smooth. For each of the finitely many \(d\)
compute the fundamental positive Pell unit \(X_1+Y_1\sqrt d\)
by continued fractions, and enumerate
\[
       X+Y\sqrt d=(X_1+Y_1\sqrt d)^h,\qquad
       1\le h\le \max(12,1+\max S).
\]
This suffices by the real Lucas primitive divisor theorem
in the form used by
[Luca and Najman, Section 2, pp. 430–431](https://web.math.pmf.unizg.hr/~fnajman/FLFNMC.pdf).
In detail \(Y_h/Y_1\) is the Lucas sequence for the two
conjugate Pell units. At \(h>12\) a primitive divisor
has rank of apparition \(h\) and is congruent to \(1\)
or \(-1\) modulo \(h\), so an \(S\)-smooth \(Y_h\)
forces \(h\le1+\max S\). Retain \(X\) odd and
\(n=(X-1)/2\ge1\) only if \(n,n+1\) both pass the
\(S\)-smooth test. This list is complete by the Pell
transformation. For a longer window discard starts
whose later indices have any prime outside \(S\), and
evaluate the exact moments at the remaining starts.

The finite list depends on \(S\). It gives no common
finite cutoff when the normalized coefficient denominator
primes range over all finite sets. A phase with zero constant
has no constant residue to survive Frobenius in this proof.
