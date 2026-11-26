# Prepared cores for endpoint-aligned factorial tails

Fix a finite ordinary two-monomial core \(C\), possibly empty.
Choose an integer \(J\ge4\) at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md),
and a positive integer \(D\). Let \(\mathcal P\) be any
nonempty finite set of primes greater than \(J^2D\).

At a prime \(p\in\mathcal P\), call a tail row admissible
if, in either orientation, it is obtained by multiplying
both entries of \((a,a+d)\) by a common power \(p^b\),
where \(1\le d\le D,\ a\ge0,\ b\ge0\), and one of the
following holds.

- The larger reduced exponent satisfies \(a+d\equiv0\bmod p\).
- All reduced ramp factors \(ja+t\), for \(2\le j<J\)
  and \(1\le t\le jd\), are \(p\)-units.

The second condition gives many residue classes for the
smaller exponent, since the union of excluded residues
has size at most \(d\sum_{j=2}^{J-1}j<p-1\).
The first condition singles out the larger endpoint.
The integer inflation powers can vary independently
and the actual gaps need no upper bound.

**Theorem.** Append a fixed balanced list of unit-gap
rows to \(C\), at arbitrarily small total curvature,
to form a prepared core \(C'\). One can choose the list
so that for every \(p\in\mathcal P\) and every \(2\le n<N=n+1<J\),
any finite collection of rows admissible at \(p\) leaves
the two Newton slope sets separated by more than an
arbitrary prescribed \(\delta>0\).
Both polygons have only width-one edges.
Consequently every ordinary appendage admissible at
one witness prime has coprime moments with \(C'\) in
every degree. The witness may also depend on the
adjacent pair. The construction requires no minimum
tail count or signed degree balance. It permits the empty
tail and uses at most \(2|\mathcal R|\) preparation rows,
where
\[
       \mathcal R=\{k/j:2\le j<J,\ 1\le k\le j\}.
\]
The same labels serve at all primes.

Two further assertions permit rows outside the admissible
classes. Arbitrary extra rows with sufficiently small
residual valuation oscillation at a witness prime preserve
the conclusion. A single arbitrary \(p\)-inflated unit-gap
row preserves it without a depth bound. Using two witness
primes of opposite slope ordering also permits a single
arbitrary reduced gap \(d\le D\), with an arbitrary
inflation at the prime used for its orientation. Any number
of synchronous remainder rows of the same orientation can
also occur when their cutoffs agree at the tested pair.

## Endpoint valuations and common slope removal

Strip the smaller-end factorial from each appended row.
Write \(j\) for its order and \(dk\) for the number of
reduced factors used at index \(k\) of an increasing
row, or \(d(j-k)\) in the decreasing case.
In an uninflated row of the first class, \(ja+t\) is
divisible by \(p\) exactly at \(t=jd\). Indeed \(jd<p\)
and \(a+d\equiv0\bmod p\). The depth there is
\(v_p(a+d)\), independent of \(j\), since \(p>j\).
It appears only at the endpoint using the larger exponent.
Every valuation from the second class is zero.

For inflation, if \(L\) is the number of reduced factors
used, the factorial valuation identity is
\[
 v_p\!\left(\frac{(p^b(ja+L))!}{(p^bja)!}\right)
   = L\,\frac{p^b-1}{p-1}
      + v_p\!\left(\frac{(ja+L)!}{(ja)!}\right).               \tag{1}
\]
This follows by separating the first \(b\) layers in the
factorial floor formula, and is the valuation part of the
[row inflation identity](SFC_PADIC_ROW_INFLATION_HIGH_PRECISION.md).
The extra term is an affine function of the coefficient
index with the same slope in all orders. The slope is
\(d(p^b-1)/(p-1)\) for an increasing row and its negative
for a decreasing row. Removing these common slopes is a
common rational variable change. Scalar intercepts may
differ between orders without affecting their polygons.

After the common change, the total tail valuation vector
at every order is
\[
                         (T_-,0,\ldots,0,T_+),               \tag{2}
\]
where \(T_\pm\) is the sum of the depths \(v_p(a+d)\)
over each orientation in the first class. Thus a
contributing inflated row adds depth
\(v_p(p^b(a+d))-b\) at its endpoint. In particular,
different allowed residue classes and any number of
neutral rows can occur together.

## An explicit preparation using unit gaps

Work first at one prime \(p\in\mathcal P\). For each
\(r\in\mathcal R\) reserve an increasing unit-gap row and
a decreasing unit-gap row, with respective smaller bases
\(a_r^+,a_r^-\). Prescribe
\[
                v_p(a_r^++r)=X_r,\qquad v_p(a_r^-+r)=Y_r,
\]
where \(X_r,Y_r\) are positive integers to be chosen.
All the fractions in \(\mathcal R\) have unit denominators
and distinct residues at \(p\): the numerator of their
nonzero pairwise difference has magnitude less than \(J^2\).
Consequently each label contributes a mark precisely
where its fraction appears in the factorial ratio.
If \(\rho_{j,k}\) is the valuation of the \(k\)-th
adjacent coefficient ratio of the original core, the
prepared ratio valuation is exactly
\[
       g_{j,k}=\rho_{j,k}+X_{k/j}-Y_{(j-k+1)/j}.
                                                                    \tag{3}
\]
The two signs in (3) follow from the increasing ramp
factor \(ja+k\) and the decreasing ramp factor
\(ja+j-k+1\) in the denominator.

Choose an orientation sign \(\epsilon\in\{1,-1\}\).
For \(r<1\), use
\[
 X_r=L+(\epsilon A+B)r,\qquad
 Y_r=L+\epsilon A(1-r),                                      \tag{4}
\]
where \(A,B\) are positive multiples of the least common
multiple of \(2,\ldots,J-1\), and \(L>A+B\) is an integer.
For interior ratio indices \(2\le k<j\), (3) becomes
\[
                 g_{j,k}=\rho_{j,k}+\epsilon A/j+Bk/j.        \tag{5}
\]
First take \(B\) sufficiently large so the interior
ratio valuations increase strictly in \(k\) for every
order. Then take \(A\) sufficiently large compared with
\(B\), the fixed \(\rho_{j,k}\)'s and \(\delta\).
At any adjacent pair, all central ratio valuations
of the lower order exceed all those of the upper
order by more than \(\delta\) if \(\epsilon=1\).
For \(\epsilon=-1\) the upper order exceeds the lower
one by more than \(\delta\). This follows from
\(1/n-1/N=1/(nN)\) in (5). The two first ratio
valuations have respective difference
\(\rho_{n,1}-\rho_{N,1}+(\epsilon A+B)/(nN)\), and the
two last ratio valuations have difference
\(\rho_{n,n}-\rho_{N,N}+\epsilon A/(nN)\).
Thus these differences have absolute values greater
than \(\delta\) and sign \(\epsilon\) as well.

Choose \(K\ge1\) bounding \(\delta\), all central ratio
valuations in absolute value, and the corresponding
first and last ratio differences at adjacent orders.
These differences do not depend on \(X_1,Y_1\).
Take \(X_1,Y_1\) sufficiently large so that every
first ratio valuation lies below \(-(2J+4)K\) and
every last one above \((2J+4)K\).
The prepared polygons now have strictly increasing
adjacent slopes. Adding (2) decreases the first
adjacent slope by \(T_-\) and increases the last
one by \(T_+\), preserving strict convexity.
The separation between all central and exterior
types follows from the buffers, and same-side
exterior differences are unchanged across the two
orders because \(T_\pm\) are common. This proves the
claimed margin after the common affine change.

Apply these prescriptions independently at every chosen
prime, using the same \(2|\mathcal R|\) labels. All
fractions have representatives with unit denominators
at each prime, so the congruence
\(a_r^++r\equiv p^{X_r}\bmod p^{X_r+1}\) and its
decreasing counterpart exist as ordinary residue
classes. The Chinese remainder theorem combines the
classes for each label. There is no varying label
count or padding requirement. All prepared bases can
be made arbitrarily large, so the added curvature
tends to zero while the chosen valuations remain exact.
Preparation has signed gap zero.

The Newton polygon product rule implies coprimality
at a pair with disjoint slope sets.
[Filaseta and Rocks, Theorem 2.1](https://doi.org/10.4064/cm6761-3-2016)
The terminal \(2\) is covered by even positivity, and
the absolute cutoff takes care of \(N\ge J\).

## A depth budget for arbitrary remainder rows

At a tested prime \(p\) and pair \(n,N\), allow any additional
ordinary rows outside the admissible classes. For a remainder
row \(i\) let \(L_{i,p}\) be the maximum \(p\)-valuation of
any factor in its smaller-end-stripped ramps at the two
orders. The [common slope removal theorem](MULTIVARIATE_RETAINED_COMMON_SLOPE_BUDGET.md)
subtracts a common slope from these contributions and
bounds the resulting ordinate oscillation at each order
by \(\sum_i L_{i,p}\). If their actual residual oscillations
are \(K_{n,p},K_{N,p}\), the two slope lists move by at most
these amounts respectively, after the common shift.
Thus \(K_{n,p}+K_{N,p}\le\delta\) preserves disjointness.
In particular the certificate permits
\(2\sum_i L_{i,p}\le\delta\). There is no bound on the
heights, gaps or row count of the remainder in the
oscillation formulation. Cancellation can be measured
using the [cyclic discrepancy bound](MULTIVARIATE_RETAINED_CYCLIC_DISCREPANCY.md).

For example, if the main tail is admissible at all
the witnesses, the same arbitrary remainder is permitted
at a pair whenever
\[
       2\sum_{p\in\mathcal P}\sum_i L_{i,p}\log p
                         \le \delta\sum_{p\in\mathcal P}\log p .
\]
A weighted average then gives a witnessing prime
for that pair. One can use a sufficient real size bound
without factoring the ramps, since
\(\sum_p L_{i,p}\log p\le (n+N)d_i\log(JM_i)\)
for a remainder row of absolute gap \(d_i\) and larger
exponent \(M_i\). The prepared margin \(\delta\) can be
chosen to accommodate any prescribed finite budget.

## A single arbitrary remainder at unbounded depth

At a given witness prime consider one arbitrary remainder
with reduced gap \(d\le D\), in either orientation, after
a possible common \(p\)-power inflation. Remove its
common affine inflation term by (1). Since each reduced
ramp has length at most \(JD<p\), there is at most one
divisible factor in each order. If the remainder is marked
in only one of the orders, it adds just one step to its
valuation vector. An edge that does not cross the
jump stays an original width-one edge: points within
either constant plateau form a strictly convex
chain, so any longer chord would lie above an
intermediate point. The only possible new shared
slope is the crossing edge, and the other order has
width-one slopes. The common factor degree is therefore
at most one. Even positivity excludes it. The same
reasoning works for any number of remainder rows whose
only marks at the pair aggregate to a single step in
one order.

If the remainder has marks in both orders, their
offsets \(t_n,t_N\) satisfy
\(Nt_n=nt_N\). Indeed the difference is divisible
by \(p\) and its absolute value is less than
\(nNd<p\). Thus \(t_j=ju\) with a common integer
\(1\le u\le d\), and both marks have the same depth
\(e=v_p(a+u)\). The case \(u=d\) is the admissible
endpoint case. In particular every unit-gap remainder
is already covered at a single witness with either
ordering sign.

For arbitrary reduced gap, take two witnesses with
opposite signs \(\epsilon\) in (4). For an increasing
remainder use the witness with \(\epsilon=1\).
For a decreasing remainder first reverse both columns
and use the witness with \(\epsilon=-1\).
In the reversed variable the two orders' slopes
are negated and reversed. Their central bands and
same-side exterior slopes therefore have the same
ordering as for the increasing case, and their large
buffers still apply. The base tail is assumed admissible
at both witnesses.

It remains to treat the increasing synchronous
mark with \(u<d\). Write
\(\ell_j=\lceil ju/d\rceil\), so the positive step of
depth \(e\) starts at index \(\ell_j\).
Write \(g'_{j,k}\) for the unperturbed slopes after
adding the admissible tail and removing its common
affine part. The first and last \(g'_{j,k}\) include
\(-T_-\) and \(T_+\), respectively. For each order
the crossing edge starts at index \(\ell_j-1\):
any earlier start would give a chord strictly above
the preceding convex plateau at \(\ell_j-1\).
Its slope \(s_j\) is the minimum of
\[
   S_j(w)=\frac{e+\sum_{k=\ell_j}^{\ell_j+w-1}g'_{j,k}}{w},
                  \qquad 1\le w\le j-\ell_j+1 .                \tag{6}
\]
All remaining edges retain their old width-one slopes.
Those after the crossing have slope strictly above \(s_j\),
taking the full maximal edge at slope \(s_j\).

The two cutoffs differ by zero or one. If both are at
least two, then \(N-\ell_N+1\ge n-\ell_n+1\).
For every width available at \(n\), the upper candidate
in (6) is strictly smaller than the lower one. This
uses the ordered central bands, and the fact that the
upper rightmost slope is included only when the lower
rightmost one is also included. Hence \(s_N<s_n\).
The same comparison holds if \(\ell_n=\ell_N=1\),
since the ordered first exterior slopes are paired.
In this latter case \(s_N\) cannot equal any old lower
slope, because they all lie after the crossing.
Thus at most one shared width-one branch remains.

Suppose both cutoffs are at least two.
The lower crossing slope exceeds all upper central
slopes and the upper left exterior slope.
It can meet only the old upper right exterior slope.
If it does, it is strictly below the lower right
exterior slope. Its edge cannot reach the lower
right endpoint: if the plateau contains the preceding
index, that point lies below such a chord, and if the
crossing width is one at the last index the slope
exceeds the old right slope by \(e\).
So the sum in its minimizing candidate (6) uses
central slopes alone. Therefore \(e\ge s_n-K\).
The upper crossing then satisfies
\(s_N\ge-K+e/J>K\), by the right buffer.
It exceeds all old lower central and left exterior
slopes. It lies below \(s_n\), which is below the
old lower right exterior slope. Hence it matches
no old lower slope, and again at most one shared
width-one branch remains.

Finally, if \(\ell_n=1,\ \ell_N=2\), both candidate
ranges in (6) have the same maximum width \(n\).
At each equal width the lower candidate is strictly
smaller: its first slope contains the negative left
buffer, the central differences total at most \(2JK\),
and the right exterior contribution, when included,
has difference bounded by \(K\). Hence \(s_n<s_N\).
The old upper slopes after its crossing are above
\(s_N\). Its only old slope before the crossing is
the left exterior. The lower crossing is at least its
own strictly larger left exterior slope, since it
averages an increasing slope list and adds \(e/w\).
Thus it cannot match the old upper left slope either.
This last case still allows at most one shared
width-one branch. Even positivity completes the
synchronous proof.

At a tested pair, any number of increasing synchronous
rows with the same cutoff pair
\[
                      \bigl(\lceil nu/d\rceil,\lceil Nu/d\rceil\bigr)
\]
fit into the same proof. After inflation removal the total
perturbation is a single step at each order, with the same
total depth \(\sum v_p(a+u)\). The arguments with (6) apply
to this total depth regardless of the number of rows.
Reversal gives the decreasing version. A simultaneous
congruence family across all the low orders is obtained
by taking all the reduced fractions \(u/d\) in one
left-open, right-closed interval of the grid
\[
                 \{0,1\}\ \cup\
                   \{k/j:2\le j<J,\ 1\le k<j\}.
\]
For each row use \(1\le u\le d\le D\) and
\(a+u\equiv0\bmod p\) in its reduced variables.
The cutoff \(\lceil ju/d\rceil\) is constant on the
chosen interval at every such order. Thus the added
gap sizes, heights, inflation powers and multiplicities
may vary independently within the family.

## Limitation for several synchronous cells

Fix any ordinary core, an adjacent pair \(n,N=n+1\) with
\(n\ge4\), a prime \(p>N^2\), and a precision \(L\ge1\).
Allow reduced gaps up to \(N\). With several synchronous
cells available, one can append increasing rows and
balanced endpoint buffers whose two central Newton
faces have a shared quadratic modulo \(p^L\).
All their heights can be arbitrarily large. As above, divide
each added row out by its smaller-end factorial at each order,
a scalar common to all coefficients of that order.

Use one increasing row for every distinct internal
fraction \(r=k/j\) of the two orders. Give it gap \(j\)
and choose its smaller base \(a_r\) so \(a_r+k\)
has positive valuation \(e_r\). At either parent order
\(i\), its sole marked factor is \(i(a_r+k)\), with
that same valuation. It enters at the coefficients
with index \(l/i\ge r\). If \(H(r)\) is the cumulative
sum of the depths in increasing \(r\), the resulting
interior ordinates at order \(j\) are \(b_{j,k}+H(k/j)\),
with \(b_{j,k}\) the original core valuations.

Choose \(H(k/j)=sk+c_j-b_{j,k}\) so both parents
have central slope \(s\). The two internal meshes
interlace in the order \(k/N<k/n<(k+1)/N\).
First choose \(c_n-c_N\) large enough for an increase
of at least \(L\) at every step \(k/N<k/n\);
then choose \(s\) large enough for that increase
at each remaining step. A common addition to both
\(c_j\)'s makes the first \(H\)-value at least \(L\).
Use the resulting successive differences as \(e_r\).
Append one endpoint-aligned increasing unit row and
one endpoint-aligned decreasing unit row.
Their two depths can be chosen large enough that the
endpoint valuations exceed the corresponding central
line by at least \(L\) at both orders. This balanced
pair changes no signed gap. If a specific signed gap
is wanted, neutral unit rows of the required sign
supply the remaining difference, with smaller bases
zero modulo \(p^L\).

In the physical row factors impose
\[
               a_r+k\equiv p^{e_r}w_r\bmod p^{e_r+L},
                          \qquad w_r\text{ a unit}.
\]
Every unmarked ramp factor at order \(i\) is
\(i(a_r+k)+(t-ik)\), a unit congruent to the
nonzero offset \(t-ik\) modulo \(p^L\), since
\(p>N^2\) and \(e_r\ge L\). The endpoint buffers
can also have depths at least \(L\), so their
interior units are fixed. After coefficientwise
removal of the powers of \(p\), the unit at \(k/j\)
is therefore a fixed unit times
\(\prod_{r\le k/j}w_r\). The mesh points are
distinct. Choose the \(w_r\)'s successively to
prescribe every interior unit in both orders.
For example after substituting \(z=p^{-s}x\) and
dividing the parent by \(p^{c_j}\), prescribe
the interior polynomial
\[
                   x(1+x+x^2)\sum_{i=0}^{j-4}x^i .
\]
All its designated coefficients are \(p\)-units.
The two full scaled polynomials are congruent to
these targets modulo \(p^L\).
Thus their central faces permit a common quadratic
at the specified precision.

The same labels work at any finite set of such
primes for this pair: choose the depths and units
separately and combine each row's base congruences
by CRT. The prepared slopes and depths from one
prime impose no conditions on the others.
The added rows can have arbitrarily small total
curvature. The same profiles can satisfy SFC(2) in every
degree. First prepare the core at an additional prime
\(q>J^2N\) outside the tested primes, using the theorem
above with \(D=N\). Then carry out the staircase
construction over this core and give all staircase,
buffer and gap-repair bases the class zero modulo \(q\).
Their entire ramps are \(q\)-units. Thus the prepared
certificate at \(q\) gives all-degree coprimality, while
the other primes have the prescribed local quadratic
compatibility. Both additions can be arbitrarily weak
in curvature.

## Positive diffuse limits

Among the admissible rows are ordinary unit gaps
whose larger exponents are multiples of
\(\prod_{p\in\mathcal P}p\). Their two orientation
counts can tend to infinity together while their
heights are chosen independently.
For fixed positive limiting curvature \(\tau\) and
signed bias \(|\lambda|\le\tau\), use \(t\) rows of
each orientation with respective upper exponent of
size \(t/((\tau\pm\lambda)/2)\), rounded up to a
multiple of that product. If a denominator is zero,
use larger exponents of order \(t^2\) in that
orientation. Any fixed signed gap is supplied by
finitely many additional unit rows at still larger
aligned exponents. Thus these arithmetic classes
can approach positive physical diffuse reference
parameters, including exceptional ones, while their
prepared ordinary profiles satisfy SFC(2) in every
degree.
