# All-order coprimality on single-row arithmetic progressions

Fix any ordinary two-monomial exponent profile except for the
smaller exponent \(B\) in one selected unequal row. Its positive
gap \(d\) and its orientation are fixed. The other rows may be
shared, private, or equal. Let \(J\ge3\) be an absolute terminal
cutoff as in the
[effective eventual theorem](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md).

**Theorem.** There is an explicitly constructible residue class
of \(B\) in which every nonnegative integer base gives SFC(2) at
all orders. The construction can refine any prescribed residue
class of \(B\). Irrespective of the size or number of the fixed
rows, at most
\[
                         \frac{d(J-1)J(J-2)}6
\]
nonnegative integer bases can fail at any order. The bounds and
construction apply to both orientations.

In particular, one can append just one shared unit-gap row of
arbitrarily small curvature in order to obtain an all-order
regular profile. One can also use a balanced pair of appended
rows by fixing the second with the reverse orientation and a
sufficiently large base before selecting the first.

## A unit resultant at one selected prime

For a terminal \(N=n+1<J\), take a prime
\[
                       p>\max\{N^2d,NM,N\},
\]
where \(M\) is the largest exponent in the other rows, or zero if
there are none. Choose
\(B\equiv-1/N\pmod p\). Normalize by the row-minimum factorials
as in the [translated resultant proof](MULTIVARIATE_TRANSLATED_RESULTANT_DEGREE.md).
Thus each factor in the selected row is \(jB+h\),
where \(1\le h\le d k\) for an increasing row and
\(1\le h\le d(j-k)\) for a decreasing row.

At the lower order the selected row contributes only units.
In fact
\[
                 N(nB+h)\equiv N(h-1)+1\pmod p
                    \qquad(1\le h\le nd).
\]
The right side is strictly between zero and \(p\).
All sampled factors in the other rows, and all binomial
coefficients, are also units. Thus the lower polynomial has
unit endpoints.

At order \(N\) the factor \(NB+1\) appears at every positive
index \(k\) in the increasing orientation, so only the constant
coefficient survives, and it is a unit. In the decreasing
orientation only the coefficient at \(k=N\) survives, again a
unit. Since the lower polynomial keeps its degree and constant
coefficient modulo \(p\), the formal degree-\((n,N)\) Sylvester
resultant is a unit. Indeed, its product formula multiplies the
values of the reduced upper polynomial at all lower roots by
the \(N\)-th power of the lower leading coefficient, regardless
of the upper degree drop. No lower root is zero. Hence the
ordinary polynomials have no common characteristic-zero root.

Take distinct primes with these properties for all \(2\le N<J\).
The Chinese remainder theorem gives one residue class of \(B\)
with the asserted property at every such terminal. Higher
terminals are covered by the effective eventual theorem. To
refine a prescribed residue class, choose the new primes
coprime to its modulus before applying CRT.

There is also a Cartesian residue class when all the signed
gaps are fixed and all the smaller exponents vary. Choose the
primes larger than \(N^2\) times the largest gap. At each prime
impose the same selected-base condition and set every other
base to zero modulo that prime. All the unselected factors
then reduce to nonzero integers \(h\) in their sampled ranges,
so the same proof applies at every representative of those
classes. This version does not require an upper bound on any
of the unselected bases.

## The uniform bound on exceptional bases

At a fixed terminal, the scalar-free resultant \(R_N(B)\)
has degree at most \(nNd\) in \(B\).
For an increasing selected row, put \(B=u^{-1}\) and replace
\(z\) by \(u^dz\). Its coefficient at \(z^k\) originally had
\(B\)-degree \(dk\), so the transformed coefficients are
polynomial in \(u\). The resultant scaling law multiplies the
original resultant by \(u^{nNd}\), proving the bound. In the
decreasing orientation, replace \(z\) by \(u^{-d}z\) and also
multiply order \(j\) by \(u^{jd}\). This gives the same net
power \(u^{nNd}\). The unit-prime construction proves that these
one-variable resultants are all nonzero.

For real \(B\) in an open neighborhood of \([0,\infty)\),
each even moment in this finite order range is a strictly
positive real even-power integral. For example \(B>-1/(2J)\)
keeps all Gamma arguments positive in the range. Each
resultant is nonnegative there: evaluate the odd polynomial
at the conjugate root pairs of the even polynomial and use
the positive leading coefficient. Swapping the two orders
does not affect the sign since \(nN\) is even.
Consequently every nonnegative real zero of the product of
the \(R_N\) has even multiplicity. Its degree is at most
\[
              d\sum_{N=2}^{J-1}N(N-1)
                =\frac{d(J-1)J(J-2)}3 .
\]
The exceptional count follows. Unlike a simultaneous bound
in all row bases, this one depends on the gap of the
single varied row.
