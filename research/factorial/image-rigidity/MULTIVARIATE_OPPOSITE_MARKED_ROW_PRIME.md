# A variable-prime test with opposite marked rows

For arbitrary terminal \(N\ge2\), two opposite exponent rows can
force coprimality when no single row can be isolated at the same
prime. The mechanism is a narrow interval of unit coefficients
in the higher order.

Let the two selected rows be \((A,A+d)\) and \((B+e,B)\) with
\(A,B\ge0\), \(d,e\ge1\). Set \(M=\max(A+d,B+e)\) and let \(M_0\)
be the largest exponent of every other row, or zero if there is
no other row. For \(j=N-1,N\) let \(Q_j\) be the ordinary factorial
moment polynomial for the two monomials.

**Theorem.** Suppose a prime \(p\) satisfies
\[
                    p>\max\{N,(N-1)M,NM_0\}.
                   \tag{1}
\]
Define the set of integers
\[
 I=\{0\le k\le N:\ NA+d k<p,\ N(B+e)-e k<p\}.
 \tag{2}
\]
If \(I\) consists of one or two indices, then \(Q_{N-1},Q_N\)
are coprime over every characteristic-zero field.

At the lower order every factorial and binomial coefficient
is a \(p\)-adic unit. At the upper order every unmarked factorial
is still a unit. The two marked factorial arguments are
less than \(2p\) since \(N M\le2(N-1)M<2p\). Each marked row
therefore contributes valuation zero below \(p\) and valuation
one otherwise. Thus the reduction of \(Q_N\) is supported exactly
on \(I\), with every coefficient there nonzero. The inequalities
in (2) make \(I\) an interval of integers.

For a possible rational common factor, use a primitive
representative integral at \(p\). Both its endpoint coefficients
are units because it divides \(Q_{N-1}\), whose endpoints are
units, and primitive products preserve minimum coefficient
valuation. Its reduction retains its degree, has nonzero
constant term, and divides the reduction of \(Q_N\). Removing
the initial monomial from the latter leaves a polynomial of
degree at most one. Hence the common factor has degree at most
one. The polynomial at the even order has no real zero, since
it is the exponential integral of an even power of a nonzero
real binomial. It follows that every rational common factor
has even degree. This proves coprimality over \(\mathbb Q\),
and rational Bézout transfers it to characteristic zero.

The condition on the unmarked rows can also be replaced by
an affine valuation condition. Allow their total valuation at
order \(j\) to be \(u_j+v k\), where \(v\) is the same integer
at both adjacent orders. Use \(t\mapsto p^{-v}t\) and divide
out the separate scalars \(p^{u_j}\). The normalized pair has
the same reductions used above. More generally the argument
uses only that the lower normalized pair has both endpoints
units and that the higher normalized polynomial has
valuation-minimizing support of width at most one.

## Balanced markings with arbitrarily small curvature

Take mirror rows \((A,A+d)\) and \((A+d,A)\), and set
\(c=\lceil N/2\rceil\), \(d=N\). Given a prime \(p>N\), write
\(p=N(A+c)+r\), with \(1\le r<N\).
For all sufficiently large \(p\) we have \(A>0\) and
\((N-1)(A+N)<p\). We may adjoin any unmarked profile with
exponents at most \(A+c\). Since
\[
                  NA+c d<p<NA+(c+1)d,
\]
the unit indices in (2) are exactly the integers
\(N-c\le k\le c\).
For even \(N\) this consists of the central index alone.
For odd \(N\) it consists of the two central indices.
The theorem proves coprimality in either case. The added
curvature of the two marked rows is \(2N^2/(A+N)\to0\)
as \(p\to\infty\). The signed total-degree difference is
unchanged, and the rows are both shared. Neither marked
row satisfies the [single-row criterion](MULTIVARIATE_MARKED_ROW_VARIABLE_PRIME.md)
at this prime because the other marked row has the same maximum
exponent \(A+N\).

At terminal six, the unequal opposing rows \((5,7)\) and
\((7,4)\) work at \(p=37\), with arbitrary subordinate rows
of maximum exponent at most six. The lower maximum marked
argument is \(35\), and the upper unhit indices are exactly
two and three. Thus the higher reduction is a monomial times
a degree-one polynomial, while the actual common factor would
have even degree. The selected rows by themselves have
curvature \(13/7\).

This test uses a prime chosen jointly from two exponent rows.
The simultaneous obstruction to a prescribed finite collection
of prime screens remains in force.
