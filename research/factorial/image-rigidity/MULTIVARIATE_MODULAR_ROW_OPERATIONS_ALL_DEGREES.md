# All-degree coprimality under neutral modular row operations

Fix any ordinary two-monomial core that satisfies SFC(2) in every
degree. There is a prime \(p\) with the following simultaneous
closure properties.

- Add any nonnegative multiple of \(p\) to both entries of each
  core row. The multiples may be chosen independently.
- Append any collection of blocks of \(p-1\) shared rows with the
  same signed gap \(d\) in a block and with smaller exponents
  congruent modulo \(p\). If \(r\in[0,p-1]\) is their common
  residue, require \(J(r+|d|)<p\), where \(J\) is an absolute
  terminal cutoff. Their smaller exponents may be distinct
  positive integers with that residue.
- Multiply any resulting row by any nonnegative power of \(p\),
  independently. One may also pack rows by separating their base
  \(p\) digits to preserve all certificates below \(J\).

Every resulting ordinary profile satisfies SFC(2) in every degree.
The operations may be used together, with the indicated shift and
append operations carried out before the inflation or digit packing.
The added curvature and the number of added blocks are unrestricted.
This includes arbitrarily signed shifts of the total exponent gap
by multiples of \(p-1\), using blocks of unit-gap rows of either
orientation.

## One common prime

Take \(J\ge4\) at least the
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
The value of \(J\) is computable and unevaluated here.
For a row \((a,b)\), divide its factorial contribution at order
\(m\) by \((m\ell)!\), with \(\ell=\min(a,b)\).
These are positive scalar divisions independent of the argument.
The remaining factors are
\[
                 \prod_{h=1}^{|b-a|j}(m\ell+h)
\]
with \(j=k\) or \(j=m-k\), depending on orientation.
The resulting core polynomials \(F_m\) have integral coefficients
and the same roots as the ordinary moments.

Let \(M\) be the largest core exponent. Choose a prime
\(p>J(M+1)\) avoiding the nonzero integer resultants of
\(F_n,F_{n+1}\) for \(1\le n<n+1<J\), and their leading
coefficients. The hypothesis ensures the resultants are nonzero,
so all sufficiently large primes work.
At the chosen prime each of these pairs has the intended
degrees and coprime reductions.

A change \(\ell\mapsto\ell+pu\) with the orientation and gap
held fixed preserves each factor \(m\ell+h\) modulo \(p\).
Thus all independently shifted rows give the same reductions.
In each appended block, every factor in a row of gap magnitude
\(e=|d|\) lies, modulo \(p\), between \(1\) and
\(m(r+e)<p\) for \(m<J\).
Its value depends only on the common residue \(r\), the
orientation, and the sampled index \(k\).
Hence the \(p-1\) rows multiply each normalized coefficient
by a unit raised to the power \(p-1\), equal to one modulo
\(p\) by Fermat's theorem. The same reductions are obtained
after arbitrarily many such blocks.

## Inflation and digit packing at the same prime

For a factorial argument \(x\ge0\), the exact
[p-adic inflation identity](SFC_PADIC_ROW_INFLATION_HIGH_PRECISION.md#factorial-coefficients-at-all-orders)
gives
\((p^e x)!=x!\lambda_e^x U_e(x)\), with \(\lambda_e\) independent
of \(x\) and \(U_e(x)\) a \(p\)-adic unit congruent to one modulo
\(p\). The factor \(\lambda_e^x\) is removed by one diagonal
rescaling of the two coefficient variables, simultaneously
for all orders. Relative to the scalar normalization of the
original profile, the rescaled polynomials thus have the
same integral reductions even when many rows are inflated
independently. A p-adic scalar division here amounts to
multiplying the rational polynomial by a nonzero number.

The [digit packing](SFC_PADIC_ROW_MERGER_INITIAL_FORMS.md)
compares two rows \(A,B\) with the single row \(p^e A+B\).
Choose \(p^e>J\max_i B_i\) at each packing.
After the same diagonal rescaling the additional ratio of
factorials is
\[
        \frac{(p^e x+y)!}{(p^e x)!\,y!}
                           =\prod_{h=1}^y(1+p^e x/h)
                           \equiv1\pmod p ,
\]
since \(0\le y<p^e\).
This again preserves the normalized reductions in all sampled
orders. Digit separation also preserves distinctness of the
two columns. Iterating gives the stated packing operation
at the same prime.

A pair of rational polynomials of the same degrees as their
coprime reductions has nonzero resultant. The scalar and
diagonal changes are invertible over \(\mathbb Q\), so this
proves the assertion at \(N<J\). The absolute theorem
handles \(N\ge J\) for every resulting ordinary profile.
Rational Bézout identities transfer the conclusion to every
characteristic-zero field.

## Relation to integer twists

The same prime rules out every twist
\(((n+1)/n)^{\Delta}\) with \(\Delta\) divisible by \(p-1\):
for \(n+1<J\) it is a \(p\)-adic unit reducing to one.
Thus these integers avoid the
[finite exceptional set of the core](MULTIVARIATE_ARBITRARY_CORE_TWISTED_TAIL_STABILITY.md#the-finite-set-and-the-stability-constant).
The neutral row theorem gives a separate exact implementation
for the indicated residue blocks, without imposing a
small-curvature assumption.
