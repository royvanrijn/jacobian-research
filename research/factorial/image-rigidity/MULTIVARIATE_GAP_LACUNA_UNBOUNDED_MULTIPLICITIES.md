# Gap lacunae with unrestricted row multiplicities

Fix an integer \(L\ge1\), and let \(J\ge4\) be at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Assume that for every consecutive pair \(n,N=n+1\) with
\(4\le N<J\), at least one of \(n,N\) does not divide \(L\).
In particular \(L\) may be one or any prime power.

In an original ordinary two-monomial profile, choose a nonempty
group of unequal rows and assign each a positive integer multiplier
\(c_i\mid L\). Write \(d_i=|b_i-a_i|\),
\(M_i=\max(a_i,b_i)\), \(g_i=c_i d_i\), and
\(X=\max_{\rm group}g_i\). Choose an integer \(M\) at least every
\(c_iM_i\) in the group.

**Theorem.** Fix an integer \(C\ge1\).
Suppose there is a gap value \(D\ge X/C\) in the group such that
no group gap \(g_i\) lies in \([D/12,D)\), and suppose
\[
                      M-D/12\le c_iM_i\le M
                            \qquad (i\text{ in the group}).
\]
There are explicit computable \(\delta,X_0>0\), depending
only on \(C,L,J\), such that the following condition for \(X\ge X_0\)
gives ordinary SFC(2) in every degree:
\[
                   \sum_{\rm outside}d_i\log(JM_i)
                                 \le \delta\,\frac XC\,\log M .
\]
One may omit from the outside sum any row with \(M_i\le JX\).
One may also omit outside rows whose two exponents, after
multiplication by a divisor of \(L\), lie in \([M-X/(6C),M]\).
The number and total gap of all omitted rows are unrestricted.
The result uses the physical argument of the original polynomials.

To prove this, select only the group rows with \(g_i\ge D\)
in the [rational-band theorem](MULTIVARIATE_RATIONAL_UPPER_BANDS_ALL_DEGREES.md),
with its parameters \(A=C\) and \(B=1\).
These rows have scaled gap in \([D,CD]\) and scaled upper
exponent in \([M-D/12,M]\).
Every remaining group row has \(g_i<D/12\), so both its
scaled exponents are in \([M-D/6,M]\). It is therefore
passive throughout the progression of marks in that theorem,
and costs nothing in the retained prime budget.
An outside row with \(M_i\le JX\) obeys
\(JM_i\le J^2X\le J^2CD\le HD\), where \(H\) is the
rational-band offset constant. The other permitted outside
rows satisfy the same scaled passive-strip condition.
Thus the retained budget is at most
\(\delta (X/C)\log M\le\delta D\log M\).
One can take \(\delta=1/(64J^3L)\) and \(X_0=CD_0\), with \(D_0\)
from the rational-band theorem at these parameters.

## A bounded number of distinct gaps

Suppose instead that at most \(K_0\) distinct scaled gaps occur
in the group, with no bound on their multiplicities.
Take \(C=12^{K_0-1}\) and suppose
\[
                 M-X/(12C)\le c_i M_i\le M
                       \qquad (i\text{ in the group}).
\]
Write the distinct scaled gaps as
\(d_1^*=X>d_2^*>\cdots>d_r^*>0\) and append \(d_{r+1}^*=0\).
At the first \(i\) with \(d_{i+1}^*<d_i^*/12\), set \(D=d_i^*\).
All previous successive ratios are at most twelve, so
\(D\ge X/12^{i-1}\ge X/C\). The gap interval \([D/12,D)\)
is empty, and the upper-spread bound above is at most \(D/12\).
The theorem therefore applies with the same outside budget.

For an ordinary example take \(L=1\), an integer \(X\to\infty\),
and \(M=2X^3\). Use \(X\) copies of each orientation of the row
\((M-X,M)\) and \(\lfloor M/4\rfloor\) copies of each orientation
of the row \((M-1,M)\). Keep the two reflected rows \((1,2),(2,1)\)
outside the group. There are only two group gap values and no
upper spread, but their multiplicities grow without bound.
The signed total gap is zero. The curvature is
\(1+2\lfloor M/4\rfloor/M+1\), which tends to \(5/2\), and
the outside rows have a fixed positive relative gap.
The arithmetic criterion proves all-degree physical coprimality
for all sufficiently large \(X\). Any bounded-size selected
subgroup leaves unboundedly many gap-\(X\) rows of upper exponent
\(M\) outside it, beyond both passive cutoffs and the small
retained gap budget in the
[bounded-row theorem](MULTIVARIATE_BOUNDED_UPPER_CLUSTER_GAP_SPREAD.md).
