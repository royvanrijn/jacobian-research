# All-degree coprimality for borrowed Jacobi paths

Let \(J\) be the
[absolute eventual cutoff for ordinary binomials](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md).
There is a prime \(p\), chosen from \(J\) alone, with the
following property. Take a path \(A_0,\ldots,A_m\), \(m\ge1\),
with endpoints \(A_0=(1,2)\) and \(A_m=(0,1)\) and arbitrary
positive integer rows \(A_i\) at the internal vertices.
Choose each integer \(e_i\) so that
\[
                p^{e_i}>J\max_\nu A_{i+1,\nu},\qquad 0\le i<m .
\]
Use the \(m\) rows \(p^{e_i}A_i-A_{i+1}\) as the ordinary
exponent rows of a binomial. Its consecutive factorial moment
polynomials are coprime at every positive order over every
characteristic-zero field.

This gives fully shared ordinary examples in arbitrary
dimension without a retained regular core. The moment
comparisons below \(J\) are with a Jacobi family for which
the entire consecutive common factor is explicitly known.
There is no bound on the internal vertex exponents.

## The reduced arrays, including the endpoint

All entries of the actual exponent rows are positive. Their
columns are distinct, since the last row has signed difference
\(p^{e_{m-1}}(A_{m-1,2}-A_{m-1,1})-1\ne0\).
For order \(j\) and coefficient index \(k\), put
\(x_i=(j-k)A_{i,1}+kA_{i,2}\). Then \(x_0=j+k\), \(x_m=k\),
and all other \(x_i\) are positive. Assume \(1\le j<J\).
Use equation (3) of the
[borrowed factorial comparison](MULTIVARIATE_BORROW_CYCLE_ROW_TRANSFER.md#the-borrowed-factorial-identity)
with \(s=0\) and precision \(p\). Divide each row contribution
by \(-p^{-e_i}\Lambda_{e_i}^{x_i}(-1)^{x_{i+1}}\).
The product of these normalizations is a rational scalar
and one common rational diagonal change of coefficient
variables, independent of the order. For \(k\ge1\) each
normalized contribution is
\((x_i-1)!/(x_{i+1}-1)!\) times a unit congruent to one.
These quotients telescope.

At \(k=0\) only the final low argument is zero. Its actual
factor is then \((p^{e_{m-1}}x_{m-1})!\), so the inflation
identity in the same comparison, after the prescribed
normalization, gives
\[
          -p^{e_{m-1}}x_{m-1}(x_{m-1}-1)!
                      \cdot (1+p\mathbb Z_p).
\]
The prior path factors telescope with \((x_{m-1}-1)!\),
leaving a multiple of \(p\) times \((j-1)!\).
Thus the entire normalized ordinary array is \(p\)-integral
and reduces to the integer polynomial
\[
              Q_j(z)=\sum_{k=1}^j
                         \binom jk\frac{(j+k-1)!}{(k-1)!} z^k .
                                                               \tag{1}
\]
This calculation retains the ordinary constant coefficient:
it is a positive or negative nonzero rational number after
normalization, and vanishes only in the reduction.

## Orthogonality and the common width

For \(r\ge0\) put
\[
 R_r(t)=\frac1{r!\,t}\frac{d^r}{dt^r}
                  \bigl(t^{r+1}(1-t)^r\bigr).
\]
Expansion gives \(Q_j(z)=j!zR_{j-1}(-z)\), including \(j=1\).
These are shifted Jacobi polynomials. The coefficient
formula gives \(R_r(0)=r+1\) and \(\deg R_r=r\).
The \(R_r\) are orthogonal for the weight \(t\) on \([0,1]\):
multiply by \(t\) and a polynomial of degree less than \(r\),
integrate, and integrate by parts \(r\) times. The derivatives
through order \(r-1\) of the expression being differentiated
vanish at both endpoints.

The monic normalizations satisfy a three-term recurrence
with nonzero preceding coefficient. In the expansion
of \(tR_r\) in the monic basis, the coefficients below
degree \(r-1\) vanish by orthogonality, and the coefficient
of the monic \(R_{r-1}\) is the positive ratio of the
squared norms of the monic \(R_r\) and \(R_{r-1}\).
A common root of consecutive polynomials would therefore
propagate to the nonzero constant \(R_0\).
It follows, since \(R_r(0)\ne0\), that the gcd of each
\(Q_j,Q_{j+1}\) over \(\mathbb Q\) is exactly \(z\).

Choose \(p>2J\) avoiding the nonzero integer resultants
of \(Q_j/z,Q_{j+1}/z\) for \(1\le j<j+1<J\).
All these reduced polynomials then have their stated
degrees, their constant coefficients are units, and the
gcd of the two reduced arrays is exactly \(z\).
The normalized ordinary parents therefore have unit
leading coefficients at \(p\). By Gauss's lemma over
\(\mathbb Z_p\), a common factor can be made primitive
with unit leading coefficient, so retains its degree
on reduction. Its reduction divides both parents.
The common degree in characteristic zero is consequently
at most one.

One ordinary parent at any adjacent pair is an even moment.
For a real coefficient ratio its positive-orthant factorial
integral is strictly positive. Indeed, an even power is
nonnegative and the distinct monomials make the phase a
nonzero polynomial, which is nonzero on an open subset
of the positive orthant. The exponential density there
is positive. Thus the even parent has no real zero and
disallows a rational common linear factor.
This proves coprimality below \(J\); the absolute theorem
supplies every other terminal for the actual ordinary
profile. Rational Bézout identities give the result
over every characteristic-zero field.

An internal vertex may have equal coordinates. The
argument uses positivity of each coordinate and not
strict inequality of the coordinates inside any vertex.
Rows in the resulting ordinary path whose two entries
agree contribute only positive scalar factors.
