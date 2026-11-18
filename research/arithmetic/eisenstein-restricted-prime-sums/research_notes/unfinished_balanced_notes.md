# Balanced Eisenstein sums: a preserved route, not a completed theorem

17 November 2026. These notes deliberately separate proved elementary
calculations from proposed uses of the quaternionic spectral machinery.

## Established finite statement

The companion exact sieve and witness replay checked every target
`N = 6k+2 < 100,000,000`. Under the condition

\[
N=a+b,\qquad a,b\in\mathcal A_3,
\qquad N/3\le a,b\le2N/3,
\]

the exceptions in this finite range are exactly **8 and 44**. Here
\(\mathcal A_3\) includes 1 and otherwise consists of integers whose prime
factors are all 1 modulo 3. This does not prove that there are no later
exceptions.

## An exact polynomial certificate

Put \(t=(a-b)/N\), \(z=t^2\), and define

\[
P(z)=18000(1/9-z)(z+1/36)
(z-3/10)^2(z-9/16)^2(z-4/5)^2(z-24/25)^2.
\]

The script `literature_balanced_polynomial.py` verifies, using rational
Bernstein coefficients on subdivided intervals, that

\[
-1/5\le P(t^2)\le\mathbf1_{|t|\le1/3}\qquad(-1\le t\le1).
\]

The sign outside the central interval also follows directly from the
factorization. Its Legendre expansion

\[
P(t^2)=c_0+\sum_{j=2,4,\ldots,20}c_jL_j(t)
\]

has exactly certified bounds

\[
c_0=\frac{3824117519863}{17459442000000}
=0.219028621869\ldots>7/32,
\qquad
\sum_{j\ge2}|c_j|<7/2.
\]

All rational coefficients and certificate receipts are preserved in
`literature_balanced_polynomial_receipt.json`. This is an elementary
certificate for the weight only; it does not establish the balanced theorem.

## Proposed archimedean spectral input

Write a quaternion as \(q=z+w j\), of norm \(n\), and let

\[
t_q=\frac{N_E(z)-N_E(w)}{n}.
\]

For even \(j>0\), the Legendre polynomial \(L_j(t_q)\) is the normalized
matrix coefficient of the weight-zero vector in the \((2j+1)\)-dimensional
representation of \(\mathrm{SU}(2)\). That vector is fixed by the twelve
units of the Eisenstein quaternion order: the complex units act through the
torus, and the other six units act by the Weyl involution, whose action on
this vector is \((-1)^j=1\).

This suggests the uniform moment estimate

\[
\left|\sum_{q\in\mathcal Q_n}L_j(t_q)\right|
\le12\tau(n)\sqrt n.\tag{M}
\]

The proposed justification is the same class-number-one quaternionic Hecke
operator, now in nontrivial archimedean weight. Jacquet--Langlands transfers
it to holomorphic weight \(2j+2\); the normalization divides its usual
Hecke eigenvalue by \(n^j\). Deligne then gives \(\tau(n)\sqrt n\).
There are no one-dimensional automorphic constituents in a nontrivial
archimedean type. **This application and its simultaneous finite-level
version still need to be written out and independently checked.**

For odd component norms, one subtracts the representations divisible by 2:

\[
\Theta_j(n)-\Theta_j(n/4).
\]

At \(n=2^e m\), the normalized local eigenvalue is
\(b_e=2^{e/2}U_e(\cos\theta)\). For \(e\ge2\),

\[
b_e-b_{e-2}
=2^{e/2}\left(u^e+u^{-e}
 +\tfrac12\sum_{r=1}^{e-1}u^{e-2r}\right),
\qquad |u|=1.
\]

Its modulus is at most \((e+3)2^{e/2}/2\), no larger than the usual
\((e+1)2^{e/2}\) bound. Thus, if (M) is established at the operator level,
the odd-norm subtraction need not increase its constant.

## Combining a finite-prime condition with the weight

Let \(Q(t^2)=P(t^2)+1/5\). Then \(0\le Q\le6/5\). The useful sieve
inequality is

\[
\#\{\text{balanced good representations}\}
\ge\sum_{\text{odd-norm representations}}P(t_q^2)
 -\sum_{\substack{p\equiv5\ (6)}}
   \sum_{\text{representations bad at }p}Q(t_q^2).
\]

The nonnegativity of Q is essential: P itself can be negative, so an
unqualified weighted union bound using P would be invalid.

The desired extension of (M) to the p-bad set comes from the tensor product
of the finite torus-normalizer orbit indicator with the same archimedean
vector. Its nonconstant Legendre moments should again have a bound
\(12\tau(n)\sqrt n\), independent of p. The constant moment contributes
the density \(2/[p(p-1)]\). This needs its own careful operator argument;
it is not a consequence of the scalar finite-quotient bound alone.

An exact elementary bound already available is

\[
\sum_{p\equiv5\ (6)}\frac2{p(p-1)}<\frac17.
\]

One proof replaces primes by all \(6k-1\), sums k=1 through 10 exactly,
and bounds the remainder by \(1/[3(60-2)]\). The script verifies the
resulting rational inequality. The unweighted odd-norm fraction is at
least 3/4. Thus the proposed weighted argument has a positive asymptotic
margin. Using \(c_0>7/32\), \(c_0<11/50\), and the bound above, the
constant part is at least

\[
\frac{21}{128}-\frac{(11/50+1/5)}7
=\frac{333}{3200}>0.
\]

This makes an eventual balanced theorem plausible once the moment lemmas
are justified. It does **not** close the gap above 100 million: directly
inserting the crude \(n^{1/5}\) divisor bounds can give a very large
effective cutoff. A sharper bound or additional argument would be needed
before claiming the complete exception set {8,44}.

## Literature position

Ramaré--Viswanadham explicitly discuss the same primitive Eisenstein
(Löschian) sifted set, but their additive theorem is ternary. Existing
sphere-sector methods might already imply an eventual balanced result.
Consequently neither this route nor eventual balance should be advertised
as novel until the relevant literature is checked more deeply. The exact
exception set, if proved, would be the more sharply defined target.
