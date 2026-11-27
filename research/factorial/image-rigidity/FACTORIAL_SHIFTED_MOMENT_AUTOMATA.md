# Automata for factorial moments after a common monomial shift

Let \(g=c+\sum_{\nu\ne0}u_\nu X^\nu/\nu!\) have arbitrary finite exponent
support, let \(\beta\in\mathbb N^d\), and put \(F=X^\beta g\).
Work over a finite extension of \(\mathbb Q_p\), with valuation ring
\(\mathcal O\) and \(v(p)=1\). Assume the \(u_\nu\) are integral and that
\(c=0\) or \(v(c)>1/(p-1)\). For every \(a\ge1\) the sequence
\[
             S_m=\frac{\Lambda(F^m)}{m!\,(m\beta)!}\pmod{p^a\mathcal O}
             \qquad(m\ge0)
 \tag{1}
\]
is computed by a finite automaton reading the digits of \(m\).
In particular, for any positive length \(s\), the starts of \(s\)
consecutive zeros at that precision are an effectively regular language.
This tests all starts at the chosen precision for a fixed valued phase.
The test concerns normalized residues.

## 1. Represent the shifted sequence in a rational coefficient array

Use \(W_N(T)=N![X^N]\exp(Tg)\) for \(N\ge0\).
The [unshifted factorial series construction](FACTORIAL_MOMENT_SERIES_PRIME_POWER_RATIONALITY.md#the-finite-matrix-at-precision-pa)
gives the recurrence for \(c=0\) by removing a block of exponent \(\nu\).
Its transition weights depend only on the recorded residue and
positive-coordinate flags for the target exponent \(N\). Keep a
variable \(Z_i\) for each coordinate and multiply the transition
for \(\nu\) by \(Z^\nu\). Then the same finite matrix, with its original
initial state, gives the rational series
\[
                   A(T,Z)=\sum_{N\ge0}W_N(T)Z^N
\]
over \(R=\mathcal O/p^a\mathcal O\). More explicitly the transition
matrix \(C(Z)\) has polynomial entries, and the inverse
\((I-TC(Z))^{-1}\) has a denominator with constant coefficient one.
A nonzero constant contributes the factor \(\exp(cT)\).
Its reduction is a polynomial by
\(v(c^j/j!)\ge j(v(c)-1/(p-1))\).

The shift acts on the divided-power basis by
\(\Lambda(X^{m\beta}X^N/N!)/(m\beta)!=
\prod_i\binom{m\beta_i+N_i}{N_i}\).
Its generating function in \(m,N\) is rational:
\[
 B(U,Z)=\sum_{m,N\ge0}U^m Z^N\prod_i\binom{m\beta_i+N_i}{N_i}
       =\frac{\prod_i (1-Z_i)^{-1}}
               {1-U\prod_i(1-Z_i)^{-\beta_i}} .
\]
Choose an integer \(M\ge1\) at least as large as every exponent entry
of \(g\), and let \(\mathbf M=(M,\ldots,M)\). Then \(S_m\) is the
coefficient at indices \((m,m,\mathbf M m)\) of
\[
              A(TZ^{\mathbf M},Z^{-1})B(U,Z).
 \tag{2}
\]
In fact the \(T\)-coefficient of order \(m\) from \(A\) carries the
exponent \(\mathbf M m-N\), and the \(U\)-coefficient of order \(m\)
from \(B\) supplies exponent \(N\). Every transition in the substituted
matrix introduces a monomial \(TZ^{\mathbf M-\nu}\) of nonnegative
exponents. Thus (2) is a rational power series in \(T,U,Z\) with
a denominator of unit constant coefficient. The construction poses
no bound on the size of \(\beta\) or on the overlap between its
coordinates and those of the exponent columns of \(g\).

## 2. Extract every start at the fixed precision

We give the finite-state extraction to handle ramified coefficient
rings and the direction \((1,1,\mathbf M)\) explicitly.
For the fixed-denominator Cartier construction over \(\mathbb Z_p\),
see [Rowland and Yassawi, Proposition 1.9 and Theorem 2.1](https://www.numdam.org/item/10.5802/jtnb.901.pdf).

Write (2) as \(P/Q\), let \(b\) be the cardinality of the residue field,
and denote all its variables by \(\mathbf z\). There is a computable
integer \(L\ge0\) such that \(D=Q^{p^L}\) satisfies
\(D^b=D(\mathbf z^b)\) over \(R\). Indeed \(Q^b=Q(\mathbf z^b)\)
modulo the maximal ideal \(\mathfrak m_R\). In \(R\), a congruence
modulo \(\mathfrak m_R^j\) for \(j\ge1\) lifts after taking \(p\)th
powers to a congruence modulo \(\mathfrak m_R^{j+1}\) by the
binomial theorem. The ideal \(\mathfrak m_R\) is nilpotent.

Use the denominator \(D\) and starting numerator \(P Q^{p^L-1}\).
For a residue vector \(\mathbf r\) with entries between \(0\) and
\(b-1\), let \(C_{\mathbf r}\) take the coefficients at indices
\(b\mathbf j+\mathbf r\) and place them at \(\mathbf j\).
Then
\[
             C_{\mathbf r}(V/D)=C_{\mathbf r}(V D^{b-1})/D.
\]
The numerator degrees in every variable remain in the box bounded
by the larger of its degrees in the starting numerator and in \(D\).
The ring \(R\) is finite, so the resulting numerator states are finite.

Set \(\kappa=(1,1,\mathbf M)\). A state also records a vector
\(\delta\) with \(0\le\delta_i\le\kappa_i\), initially zero, for
extracting coefficients at \(\kappa m+\delta\).
On reading a least-significant digit \(d\), take the remainder
\(\mathbf r\) and quotient \(\delta'\) of \(\kappa d+\delta\)
componentwise modulo \(b\) and apply \(C_{\mathbf r}\).
The new \(\delta'\) stays in the same box. At the end read the
coefficient at index \(\delta\) of the current rational state.
This proves (1). To handle offsets \(m+j\), attach the usual digit
addition carry, initially \(j\), with remaining carry digits read
after the input ends. A finite product for \(j=0,\ldots,s-1\)
recognizes the starts whose \(s\) outputs are all zero. Emptiness
of this language for positive starts is a finite reachability test.

## 3. The exact obstruction to a rational or bounded-gap transfer

For \(p=2,\ \beta=2,\ g=x\), the sequence (1) is
\(\binom{3m}{m}\bmod2\).
By the base-two binomial rule it equals one exactly when the
addition of \(m\) and \(2m\) has no carry, that is, when \(m\)
has no adjacent ones in binary. It is one at every order \(2^j\)
and zero throughout \([3\cdot2^j,\,4\cdot2^j-1]\).
An eventually periodic sequence with infinitely many nonzero
entries has bounded zero runs. A rational univariate series over a
finite field with unit constant denominator is eventually periodic.
Thus the normalized shifted sequence can be nonrational and can
contain arbitrarily long residue zero runs. The original moment is
\((3m)!\), nonzero at every order.

The exact transfer after a common monomial shift is therefore an
automaton for the indicated normalization. A sharp SFC window still
needs a normalization or a proof that detects valuations when this
residue is zero. No formal verification, independent replay or
external review is claimed.
