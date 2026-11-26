# Projective nonvanishing from independent translations of existing rows

Let the \(s\ge2\) columns \(b_1,\ldots,b_s\) be distinct ordinary
exponent vectors. Suppose there is a column ordering and there are
pairwise disjoint, nonempty blocks \(R_1,\ldots,R_{s-1}\) of existing
rows with the following property. In block \(R_i\), column \(i\) is
coordinatewise minimal among columns \(i,\ldots,s\), and for each
\(j>i\) some row of \(R_i\) satisfies \(b_{\rho j}>b_{\rho i}\).
There is no ordering condition on the entries of columns \(j<i\).
A single row uniquely minimizing at \(i\) among the remaining
columns suffices for its block. Write \(R=\bigcup_i R_i\).

Independently increase every entry in row \(\rho\in R\) by the
nonnegative integer \(L_\rho\), and leave the other rows fixed.
This operation preserves every row difference and every difference
of total column degrees. In particular it preserves homogeneity
and can make the rows being varied fully shared, with arbitrarily
small relative differences.

**Theorem.** For each set \(I\) of \(s\) distinct positive orders,
the homogeneous resultant of the moment forms, in the original
coefficient variables and with the nonzero common factorials
removed, is a nonzero integer polynomial in the selected row
translations \(L_\rho\). For any finite collection of such sets, an
explicit product of integral residue classes gives regularity at
all its order sets simultaneously. Each coordinate in the product
can be taken arbitrarily large.

The same conclusion holds on any active subset of columns of size
at least two, with its own support size in place of \(s\). For a
fixed finite collection of order sets and active subsets, the
fraction of exceptional tuples in the grid \([0,H]^R\) is
\(O(1/H)\). In the binary case the absolute order cutoff gives
all-order conclusions for one such finite collection. For higher
support size the theorem requires the prescribed collection of
orders to be finite.

## Exact factorial forms and the arithmetic certificate

For \(m>0\) and \(|\alpha|=m\), write
\(\beta_\rho(\alpha)=\sum_j\alpha_j b_{\rho j}\). After removing
\(\prod_{\rho\in R}(mL_\rho)!\) at ordinary integer translations, the
order-\(m\) factorial moment form is
\[
 G_m(L,c)=\sum_{|\alpha|=m}\binom m\alpha c^\alpha
  \prod_{\rho\notin R}\beta_\rho(\alpha)!
  \prod_{\rho\in R}\prod_{h=1}^{\beta_\rho(\alpha)}(mL_\rho+h).
 \tag{1}
\]
Thus the family is defined polynomially outside the domain of
ordinary row translations as well. Its coefficients are integral.

Order \(I=\{m_1,\ldots,m_s\}\) so \(m_1>\cdots>m_s>0\), put
\(J=m_1,\ K=|R|\), and choose an integer \(B>\max(2J,2K)\).
Set \(e_i=B^{s-i}\) for \(1\le i<s\). Choose a prime
\[
                    p>J^2\max_{\rho,j} b_{\rho j}+J
\]
and, for \(\rho\in R_i\), impose
\[
 L_\rho\equiv-b_{\rho i}-m_i^{-1}+p^{e_i}
                                  \pmod {p^{e_i+1}}.          \tag{2}
\]
The inverse is computed modulo \(p^{e_i+1}\). For a coefficient
of \(G_m\), all multinomial and unchanged-row factors in (1) are
\(p\)-adic units. In a selected row \(\rho\in R_i\), the factor
\(mL_\rho+h\) is a \(p\)-adic unit unless
\[
                    m_i h-m_i m b_{\rho i}-m=0.
 \tag{3}
\]
Indeed its residue is \(m_i^{-1}\) times the displayed integer,
whose absolute value is at most \(J^2\max_{\rho,j}b_{\rho j}+J\).
In case (3) its valuation, normalized by \(v_p(p)=1\), is exactly
\(e_i\) by (2), since \(p>m\). At each order there is at most
one such factor per row, and it can occur only if \(m_i\mid m\).
In particular, a row block with a larger assigned order cannot
contribute a positive valuation at a smaller order.

Suppose there is a common nonzero coefficient vector over
\(\overline{\mathbb Q}_p\). Scale it so that
\(\min_j v_p(c_j)=0\), and let \(v_p(0)=+\infty\). Inductively,
at step \(i<s\), the earlier coordinates satisfy
\[
                   v_p(c_j)\ge e_j/(2J) \qquad (j<i).
 \tag{4}
\]
Every term at order \(m_i\) involving an earlier coordinate
therefore has valuation at least \(e_i\). All other coordinates
have nonnegative valuation and all coefficients in (1) are
integral, while \(e_j/(2J)>e_i\) for \(j<i\).

A term supported on the remaining columns which uses any column
\(j>i\) has, in some row \(\rho\in R_i\), one factor at
\(h=m_i b_{\rho i}+1\). All its contributing columns have row
entry at least \(b_{\rho i}\), and some contributing column has
strictly larger entry. That factor has valuation \(e_i\), so
the term has valuation at least \(e_i\). The pure term \(c_i^{m_i}\)
has no such factor in its own block. Its coefficient valuation
is at most
\(\sum_{j>i}|R_j|e_j\le K e_i/B<e_i/2\), because earlier blocks
have larger assigned orders. Vanishing at \(m_i\) forces the
pure term to have valuation at least \(e_i\), and hence
\(v_p(c_i)\ge e_i/(2J)\). This proves (4) at the next step.

At order \(m_s\) every coefficient in (1) is a unit, since
all assigned blocks have larger orders. The normalization
forces \(v_p(c_s)=0\). Its pure term has valuation zero and
every other term has positive valuation. This contradicts
vanishing at \(m_s\).

There are positive representatives of each congruence class (2),
so in particular the homogeneous resultant is nonzero at some
ordinary translation grid point. Its defining elimination
property then gives the conclusion over every characteristic-zero
field. The same argument yields nonvanishing at every grid point
satisfying (2). Distinct primes for a finite list of order sets
give compatible conditions by the Chinese remainder theorem.
A restriction to an active subset keeps the blocks belonging to
its retained columns except the last, in their induced order.
Each of those blocks still distinguishes its designated column
from every later retained column. The unused rows can be treated
as fixed for the proof of nonvanishing of the resultant polynomial
on the full selected parameter space, by assigning them fixed
ordinary values. For simultaneous residue certificates with
active subsets, prescribe each unused selected translation
modulo the corresponding prime as well. It suffices to prescribe
it to be zero modulo that prime and retain the same lower bound
on the prime, so every factor in its normalized product is a
unit.

Put \(A=\sum_{\rho\in R}\max_j b_{\rho j}\). Each order-\(m\)
coefficient in (1) has total translation degree at most \(mA\).
By multihomogeneity of the homogeneous resultant, its total
degree is at most \(s A\prod_{m\in I}m\). A nonzero polynomial
of degree \(d\) in \(|R|\) variables has at most
\(d(H+1)^{|R|-1}\) zeros in the Cartesian integer grid. This follows
by induction on the number of variables, using the leading
nonzero coefficient as a polynomial in the other variables.
Applying this count to the finite product of relevant
resultants proves the density conclusion.

## A clustered homogeneous test case

Take the two selected rows \((0,1,2)\), \((4,0,1)\), and the
unchanged balancing row \((0,3,1)\), so all three total degrees
are four before translation. Use orders \(4,2,1\).
The first two rows individually distinguish the first and second
columns in their stages. Directly specializing the two translation
parameters to \(-1/4,-1/2\) makes the order-four form identically
zero: the first row leaves only the first pure term, and the
second row kills that term. The distinct valuations in (2)
keep that pure term shallower than all rivals at order four.
In the proof the first coefficient then has sufficient
valuation at order two to be eliminated.

The nearby construction in
[triangular ordinary row additions](SFC_TRIANGULAR_INTEGER_ROW_TAGGING.md)
appends specially arranged rows for any core. It specializes
each new row without a valuation hierarchy because every earlier
column is minimal on every later row. Here all translated rows
belong to the original support, and earlier columns may exceed
the later marked minima. The homogeneous resultant property and
its multidegree formula are in
[Sturmfels, *Solving Systems of Polynomial Equations*, Theorem 4.4](https://www.math.ucdavis.edu/~deloera/MISC/LA-BIBLIO/trunk/Sturmfels/cbms.pdf).
The prescribed finite-window residue classes give no simultaneous
integer translation for all positive starting orders in the
higher-support case.
