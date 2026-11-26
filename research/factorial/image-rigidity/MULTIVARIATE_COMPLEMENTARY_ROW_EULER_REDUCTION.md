# Complementary factorial rows and an Euler comparison

A pair of complementary exponent rows has an exact arithmetic
comparison with an Euler operator on the retained factorial
moment polynomials. This gives an all-order ordinary family with
two private linear variables and an arbitrary fixed positive
overlap:
\[
 U=X Z^a W^{Q-a},\qquad
 V=Y Z^b W^{Q-b},\qquad a,b\ge1 .
\]
For each \(a,b\) and each prime \(p\) there is a computable
\(L_0\) such that these ordinary monomials satisfy SFC(2) for
every integer \(Q\) divisible by \(p^{L_0}\). The modulus can
be chosen so that \(Q>\max(a,b)\).
The proof applies to the actual integer exponents at every
order. It uses the absolute eventual binary cutoff and a direct
order-independent coprimality proof for the Euler comparison.

## A complementary pair over a retained core

Fix an ordinary support with coefficient variables \(c_i\).
Let \(M_j(c)\) be the moment polynomial of the retained rows,
and append the rows \(A_i\) and \(Q-A_i\), with \(A_i\ge1\).
For a coefficient index \(\nu\) of total degree \(j\) write
\(x=\sum_i\nu_i A_i\). If \(Q>\max_i A_i\), then \(1\le x<jQ\).
The appended factorial factor is \(x!(jQ-x)!\). Divide the
augmented polynomial by \((jQ)!/(jQ)\) and change every
\(c_i\) to \((-1)^{A_i}c_i\). Multiplying the polynomial by
one further factor \(-1\), its appended coefficient multiplier
is exactly
\[
          \frac{(-1)^{x-1}\,jQ}{\binom{jQ}{x}}
               =\frac{x}{\prod_{h=1}^{x-1}(1-jQ/h)}.
 \tag{1}
\]
The sign change of the coefficient variables is independent
of the order.

Let \(J\ge1\) be a fixed order bound, \(H=J\max_i A_i\),
and \(p^L\mid Q\). If \(L>\lfloor\log_p H\rfloor\), each
factor in the denominator in (1) is a \(p\)-unit congruent
to one modulo \(p^{L-\lfloor\log_p H\rfloor}\).
Consequently the normalized augmented polynomials, uniformly
at all orders \(j\le J\), are coefficientwise congruent to
\[
                            \sum_i A_i c_i\partial_{c_i}M_j
 \tag{2}
\]
at that precision. One may divide both sides by any common
\(p\)-power in the retained polynomial before making the
comparison. The formula follows coefficientwise from (1)
and thus allows any retained exponent rows.

For a distinct binary core, suppose the Euler polynomials
(2) are coprime at every consecutive terminal \(2\le N<J\),
where \(J\) is the
[absolute eventual ordinary cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md).
For a fixed prime \(p\), choose the precision larger than the
\(p\)-adic valuation of each of their finitely many nonzero
rational resultants, after clearing any scalar denominators.
Then the augmented resultants cannot vanish below \(J\), by
(1) and the Sylvester determinant. Every other terminal
is supplied by the eventual theorem. In particular no prime
needs to give a nonzero reduction at precision one. Sufficiently
large valuation of \(Q\) at any fixed prime is enough.

## Coprimality of the Euler-geometric comparison

For the retained core \(X+zY\), divide its order-\(j\)
moment by \(j!\) to get \(H_j(z)=1+z+\cdots+z^j\).
Take the first appended row entries \(a,b\), with \(a,b\ge1\).
After reversing both polynomials if necessary set \(d=b-a\ge0\).
The Euler comparison is
\[
                     G_j(z)=\sum_{k=0}^j(ja+dk)z^k .
\]
If \(d=0\), these are scalar multiples of \(H_j\), and
\(H_{j+1}-zH_j=1\) proves consecutive coprimality.
For \(d>0\), consider positive \(j=n\). All coefficients
are positive, so a common root is neither zero nor a
nonnegative real number. The identity
\[
 (1-z)^2G_j(z)=P_j(z)-z^{j+1}Q_j(z),\qquad
 P_j(z)=ja+(d-ja)z,\qquad
 Q_j(z)=ja+d(j+1)-j b z
\]
follows by differentiating the finite geometric series.
Eliminating \(z^{n+1}\) at a common root gives
\[
 \begin{aligned}
 0&=P_{n+1}(z)Q_n(z)-zP_n(z)Q_{n+1}(z)\\
  &=-(n+1)(z-1)^2\bigl[b(na-d)z-a(nb+d)\bigr].
 \end{aligned}
\]
This uses no division by either \(Q_j\).
There is no candidate if \(na=d\), and if \(na>d\) it is
positive. In the remaining case any common root is the
negative rational number with absolute value
\[
                         \frac{a(nb+d)}{b(d-na)}.
\]

Divide \(a,d\) by their gcd, which rescales each \(G_j\)
without changing its roots. Now \(\gcd(a,d)=\gcd(a,b)=1\).
By the rational root theorem applied to both \(G_n,G_{n+1}\),
the reduced numerator of the displayed candidate divides
\(\gcd(na,(n+1)a)=a\) and its reduced denominator divides \(b\).
Since \(\gcd(a,b(d-na))=1\), its reduced numerator is
in fact \(a\). If its denominator is \(v\), one therefore
has \(v(nb+d)=b(d-na)\).
But \(\gcd(nb+d,b)=1\), so \(nb+d\mid d-na\).
This contradicts \(0<d-na<nb+d\). The polynomials
\(G_n,G_{n+1}\) are therefore coprime for every \(n\ge1\).

For the original augmented moments, also divide by the
retained factor \(j!\). The congruence (1) becomes the
comparison to these exact \(G_j\)'s, with no loss of
precision from \(j!\). Choose \(L_0\) large enough that
\(p^{L_0}>\max(a,b)\) and
\[
       L_0-\lfloor\log_p(J\max(a,b))\rfloor
           >\max_{1\le j<J-1}v_p\operatorname{Res}(G_j,G_{j+1}) .
\]
The right side is bounded effectively even without computing
the resultants individually: all columns in their Sylvester
matrices have bounded coefficient magnitude, so the ordinary
determinant height bound gives an upper bound on the
valuation of each nonzero integer resultant.
The complementary row comparison handles every terminal
below \(J\). The eventual theorem handles the rest.
Rational Bézout transfers the exclusions to every
characteristic-zero field.

