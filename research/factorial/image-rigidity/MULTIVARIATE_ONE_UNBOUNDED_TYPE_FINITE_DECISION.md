# A finite height bound for original profiles with one unbounded row type

For an ordinary two-monomial profile with distinct columns let
\(Q_j\) be its factorial moment in degree \(j\). Omit equal
rows and call \(\max(a_i,b_i)\) the height of an unequal row.
Fix an integer \(B\ge1\). Consider the profiles for which all rows
of height greater than \(B\) have the same ordered exponent pair.
There is no bound on the number of rows of that type or on the
other row multiplicities.

**Theorem.** There is a computable integer \(H_B\) such that every
physical common root of consecutive moments of one of these
profiles occurs on a profile whose unequal rows all have height
at most \(H_B\). Consequently all-degree ordinary SFC(2) is
decidable for the indicated family by finitely many exact
resultants. The bound applies to the original moments with the
same argument. No additional twist is required.

Let \(J\ge4\) be an integer at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
One can take
\[
 \begin{split}
    T_B&=12B+\min\{J(2J+24),146J-74\},\\
    C_B&=2^J(JB)^{12JB}[J(T_B+1)]^{JT_B},\\
    H_B&=B+T_B+1+(2J)!\,C_B^{\,2J}.
 \end{split}
\]
The numerical value of \(J\) is unevaluated.

**Corollary.** The all-degree original physical property for
profiles with at most two unequal rows has a computable finite
height reduction, with arbitrary shared or private exponents and
gaps. Equal rows and the original dimension are unrestricted.

## Bounding the gap of the varied row type

By the absolute cutoff it suffices to consider terminals \(N<J\).
Terminal two follows from positivity of the even moment on the
real line, so write \(2\le n<N=n+1<J\).
By the [curvature twelve theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
only
\[
 S=\sum_i \frac{(b_i-a_i)^2}{\max(a_i,b_i)}<12
\]
needs consideration.

Write \(\Delta=\sum_i(b_i-a_i)\), \(\Lambda=N/n\), and \(G_j\)
for the endpoint-normalizing argument scale in order \(j\).
The normalized coefficients lie between zero and their binomial
coefficients, with endpoints both one, by factorial log convexity.
Thus a root of the normalized order-\(j\) moment has modulus
between \(2^{-j}\) and \(2^j\), by the monic root bound and
reciprocity. A physical common root gives
\[
        |\log(G_N/G_n)|\le(n+N)\log 2 .
\]
The [fixed-terminal center comparison](MULTIVARIATE_BOUNDED_GAP_FINITE_DECISION.md#exact-algebraic-boxes-at-a-height-cutoff)
gives \(\log(G_N/G_n)=\Delta\log\Lambda+\sum_i E_i\)
with \(|E_i|\le2(b_i-a_i)^2/\max(a_i,b_i)\).
Since \(\log\Lambda\ge1/N\), we obtain
\[
                    |\Delta|\le N(n+N+24)\le J(2J+24).
\]
The [all-degree signed imbalance theorem](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
also gives \(|\Delta|<146J-73\) at any physical contact.

Take the rows of height at most \(B\) as the fixed core.
Each has gap at most \(B\) times its curvature, since its gap
is a positive integer. Their total absolute gap is at most
\(12B\). If there are any other rows, they have one common
orientation, gap \(d>0\), and smaller exponent \(A\). Let \(r\)
be their multiplicity. Their signed total gap is either
\(rd\) or \(-rd\). The two preceding bounds give \(rd\le T_B\).
In particular, there are finitely many core and varied gap
or multiplicity choices at any possible failure.

## A one-variable resultant bound

For each such choice, strip the smaller-end factorial of every
row in each order. The \(k\)-th coefficient in order \(j\) is
then an integer polynomial in \(A\). It is the product of the
binomial coefficient, fixed core rising factors, and \(r\)
copies of \((jA+1)_{dk}\) for an increasing varied row or of
\((jA+1)_{d(j-k)}\) for a decreasing one.

Its coefficient sum, equivalently its value at \(A=1\), is at
most \(C_B\). Indeed the core contains at most \(12JB\)
sampled factors, each at most \(JB\). The varied block contains
at most \(JT_B\) factors evaluated at one, each at most
\(J(T_B+1)\). Finally \(\binom jk\le2^J\).
The two-order Sylvester determinant has size \(n+N\le2J\).
The sum of absolute coefficients in its resultant polynomial
in \(A\) is therefore at most \((2J)!\,C_B^{2J}\).

This resultant is nonzero by the prime specialization underlying
the [single-row progression theorem](MULTIVARIATE_SINGLE_ROW_ALL_ORDER_ARITHMETIC_PROGRESSION.md#a-unit-resultant-at-one-selected-prime),
which also applies to repeated copies of the varied row.
Choose a prime \(p>\max(JB,N^2d,N)\), and choose
\(A\equiv-1/N\bmod p\). At order \(n\) every varied factor is
a unit, since
\(N(nA+h)\equiv N(h-1)+1\pmod p\) for \(1\le h\le nd\).
The core factors and binomial coefficients are units as well.
At order \(N\) only the constant term survives for an
increasing varied row, and only the leading term survives for
a decreasing one. Every other coefficient samples \(NA+1\).
Repeating the row leaves these assertions unchanged. The
reduced lower polynomial has unit endpoints. Its formal
resultant with either upper monomial is therefore a unit, by
the root product using the full lower degree. This proves
nonvanishing of the resultant as a polynomial in \(A\).

Its nonzero integer leading coefficient has absolute value
at least one. The Cauchy root bound now gives
\(A\le1+(2J)!\,C_B^{2J}\) at every nonnegative real zero,
so the height \(A+d\) is at most \(H_B\). This proves the
theorem. To implement the finite decision one can isolate
only the nonnegative integer zeros of the exact one-variable
resultants, checking the resulting finite profiles. A
rational Bézout identity transfers physical coprimality to
every characteristic-zero field.

For the corollary, take the computable positive parameter
\(\eta_0\) in the
[absolute diffuse angular theorem](MULTIVARIATE_DIFFUSE_ALL_ORDER_RESULTANT_TRANSFER.md)
and put \(B=\lceil12/\eta_0^2\rceil\). A failing profile
has some row of relative gap \(d_i/r_i>\eta_0\), since otherwise
that theorem excludes every angular contact. Because
\(S_i=r_i(d_i/r_i)^2<12\), that row has height at most \(B\).
A profile with at most two unequal rows then meets the
one-unbounded-type hypothesis. This supplies the stated
computable height reduction. It also shows where this
particular reduction stops in higher dimension: a failure
outside the tested family must have at least two distinct
ordered row types above \(B\).
