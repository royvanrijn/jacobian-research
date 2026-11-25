# Joint root dictionaries across letters

A full divisor on one letter identifies the other two letters' signed
ratio modulo that divisor. Taking the character average of the ratio
allows a group dictionary to use coordinates from both free letters.
Each group enters the bound at half its dictionary exponent. With two
full divisors, the two simultaneous ratios also allow groups that use
all three letters.

## One divisor and two free letters

Use the signed anatomic equation and the raw root caps of
[the grouped root dictionary count](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count).
Write \(L_T=\sum_e \alpha_{T,e}\), \(L=\sum_T L_T\le C_0\lambda\),
and \(X^{D+o_X(1)}\) for the number of triples in the box.
Choose one anatomic representation per triple. Fix the truncation
and all degrees before taking the height limit for each fixed
\(\lambda\). The subsequent budget limit is \(\lambda\to0\).
Coefficients are fixed at this height limit and roots in different
positions in each row are coprime.

Choose a prescribed root block on one letter \(T\) and fix its
exact values. Its full product \(M\) must have numerical weight
at least \(w_0>0\) and be \(\lfloor\kappa/\lambda\rfloor\)-full,
where \(w_0,\kappa\) are fixed independently of the small budget.
The choice of positions can depend on the box. A record may include
other root coordinates as well. Suppose there are at most
\(X^{f+o_X(1)}\) records in total. This applies with \(f\) equal
to the sum of the raw root caps of the frozen positions, or to the
sum of the actual dictionary costs of all groups recorded in full.

In each fixed record partition the **remaining** root positions
of the other two letters \(U,V\) into disjoint groups \(J\) of at
most \(m\) coordinates for a fixed \(m\). A group may meet both
\(U\) and \(V\). Suppose the dictionary of root vectors in group \(J\)
has at most \(X^{W_J+o_X(1)}\) elements, uniformly in the record.
For example, from an original partition on all three letters,
one may project each conditional dictionary onto the remaining
positions on \(U,V\). Discard empty projections. In particular a
group meeting the locked letter is charged at the cost of its
actual conditional projection, and the cost of any frozen values
is included in \(f\). If \(W=\sum_J W_J\) is a uniform upper bound
on the projected cost over records, then in the iterated limit,
\[
                           D\le f+\tfrac12 W+o(\lambda).              \tag{CG1}
\]
The \(o(\lambda)\) means that any fixed positive exponent allowance
times \(\lambda\) is valid at sufficiently small budgets with the
fixed margins. There is no lower degree requirement on the free
coordinates.

For the pairwise coprimality needed in the moments, first restrict
to a single prime colouring by full original coordinate position.
An independent uniform colouring retains each actual row with
probability at least \(c^{-\omega(abc)}=X^{-o_X(1)}\), with \(c\)
the number of positions at the fixed truncation. Some colouring
therefore retains the exponent. The dictionary upper bounds remain
valid. After this restriction, products of independently chosen
values in distinct positions are always coprime.

The free complete values are units modulo \(M\) and their ratio is
a fixed unit modulo \(M\). Absorb fixed factors and signs into the
target of the congruence. For the remaining factors, let \(s_i\)
be \(e_i\) on \(U\) and \(-e_i\) on \(V\). The group sum on the unit
characters of \(M\) is
\[
 A_J(\chi)=\sum_{\boldsymbol x\in\mathcal D_J}
                  \chi\!\left(\prod_{i\in J}x_i^{s_i}\right).
\]
The equation determines \(T\) from the two complete free values.
Consequently its count is at most the direct ratio congruence
count on the product of these dictionaries. By orthogonality this
is at most \(\mathbb E_\chi\prod_J |A_J(\chi)|\).

Write \(\alpha_J=\sum_{i\in J}\alpha_i\) and
\(A=\sum_J\alpha_J=O(\lambda)\). Omit zero-cap groups at inner
subpower cost. Set \(h_J=\lceil A/\alpha_J\rceil\) for the
remaining groups. When none remain, the fixed free values give at
most one triple for the record. At fixed height-exponent limits
these moment orders are fixed integers. One has \(h_J\alpha_J\le 2A\) and
\(\sum_J1/(2h_J)\le1\). In the \(2h_J\)-th group moment choose the
first \(h_J\) vectors at cost \(h_JW_J\) and collapse the other
\(h_J\) vectors to one product per position. All collapsed products
have height \(X^{O(\lambda)}\). They are coprime by the colouring,
and the signed monomial in them has a fixed unit residue modulo
\(M\). Apply (FC2), using the fixed positive numerical weight of
\(M\) and its fullness of order \(1/\lambda\). This counts these
product tuples at cost \(o(\lambda)\), uniformly in the first
choices. A given product tuple has divisor-many ordered
factorizations at each fixed moment order and truncation.
Thus
\[
    \left(\mathbb E_\chi |A_J(\chi)|^{2h_J}\right)^{1/(2h_J)}
            \le X^{\,W_J/2+o(\lambda)/(2h_J)+o_X(1)} .
\]
The cost \(o(\lambda)\) before division by \(2h_J\) is uniform
because the number of free positions per group and their collapsed
height constant are bounded. Hölder applies to the direct character
average and the sum of reciprocal moment orders is at most one.
Sum over the records to prove (CG1).
[Signed monomial congruence count (FC2)](ABC_OMITTED_COORDINATE_GRAPH.md#the-determinant-count-and-amplification)

## Two simultaneous divisors

One may instead fix full divisors \(F\mid A\) and \(G\mid B\) from
two distinct letters. They must each have a fixed positive
numerical weight and fullness at least \(\lfloor\kappa/\lambda\rfloor\).
The divisors are coprime by primitivity. Let the locks (and any
other fixed root positions) have at most \(X^{f+o_X(1)}\)
records. In a fixed record, partition **all remaining root
positions** into disjoint groups of size at most \(m\), allowing
each group to meet any number of letters. Let \(W\) bound the
sum of conditional group dictionary exponents over all records.
Then
\[
                           D\le f+\tfrac12 W+o(\lambda).              \tag{CG2}
\]

The two congruences are the signed ratios \(A/C\bmod G\) and
\(B/C\bmod F\). Take independent unit characters on the two
moduli and put each group's two signed ratios into its character
sum, so an \(A\)-coordinate appears only modulo \(G\), a
\(B\)-coordinate only modulo \(F\), and a \(C\)-coordinate in
both. The roots used in each congruence are units at its modulus
by primitivity. Enlarge the count to tuples satisfying these
two congruences with the fixed factors absorbed into their
targets. The direct character-average and Hölder step is the
same as for (CG1). In each \(2h_J\)-th group moment, choose
the first \(h_J\) vectors and collapse the other side. Solve
first for the collapsed \(A,C\) products modulo \(G\) by (FC2).
There are \(X^{o(\lambda)}\) possibilities. With those fixed,
apply (FC2) to the remaining \(B\) products modulo \(F\).
If one of these coordinate sets is empty the corresponding
congruence imposes no positive counting cost. The colouring
and exact-product factorization argument apply in both steps.
The moment costs are therefore the same as before, and prove
(CG2).

## Groups of roots in boundedly many shrinking degree windows

There is a second permitted type of group in (CG1) and (CG2).
The number of positions in this type of group can grow in the
budget limit, but must still be fixed before the height limit
for each fixed budget. In each letter the scaled degrees
\(\lambda e_i\) of its positions must be covered by at most
\(m_0\) windows of diameter \(\rho_\lambda\to0\), uniformly
in the groups. The centres may depend on the group and budget.
Equivalently, partition each letter part of the group into
at most \(m_0\) subblocks with integer representatives \(a_j\)
such that
\[
                          \lambda |e_i-a_j|\le \rho_\lambda
                  \qquad (i\text{ in subblock }j).
\]
Representatives may be taken positive by choosing a degree
in the relevant subblock.

In each group moment, after collapsing the second \(h_J\)
vectors, further group their products by subblock, each
a narrow root of height \(X^{O(\lambda)}\).
Write its original monomial contribution as the
\(a_j\)-th power of the new root times the product of the
original collapsed positions raised to \(e_i-a_j\).
After combining the correction terms into a signed
rational number the correction height is at most
\[
             X^{\,(\rho_\lambda/\lambda)h_J\alpha_J+o_X(1)}
                        =X^{O(\rho_\lambda)+o_X(1)}.
\]
Apply the wide-cofactor full-modulus count to the bounded
number of narrow roots and the wide correction. The
amplification for that count permits overlap of the wide
correction's prime support with the narrow roots. In this
application the fixed original-position colouring makes
the exact signed monomial determine each original
collapsed product: read its prime valuations and divide
by the signed original degree of its position. This also
holds for coordinatewise products of amplified tuples.
The exact fibre has size one on the compressed product tuples.
Sign cancellation is absent because all primes of a given colour
stay in the same signed original position. A product of a fixed
number of exact outputs has divisor-many ordered decompositions,
by factoring its numerator and denominator. Thus the usual
amplification applies to the compressed product tuples, with
the exponent constants fixed before the height limit.
Thus the wide-cofactor bound costs \(o(\lambda)\) as in
the bounded-coordinate moment. For (CG2), first apply
this argument to \(A,C\bmod G\), which determines their
original collapsed products, and then to \(B\bmod F\).
[Wide-cofactor count and overlap-compatible amplification](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)

## Joint entropy at saturation

Use the balanced limiting box with
\[
                 D/\lambda\to\tfrac13,\qquad
                 L_U/\lambda,L_V/\lambda,L_T/\lambda\to\tfrac13 .
\]
Suppose \(T\) supplies arbitrarily cheap positive-weight full
moduli. Precisely, for every \(a>0\), after taking the budget
sufficiently small one can choose a prescribed root block of raw
cost at most \(a\lambda\) with numerical weight bounded below and
fullness of order \(1/\lambda\), with the positive margins
allowed to depend on \(a\).
For any allowed partition of the positions on \(U,V\),
\[
             \sum_J H(\boldsymbol x_J)
                              =(L_U+L_V+o(\lambda))\log X             \tag{CG3}
\]
for the uniform row measure in the sequential limits. In
particular any disjoint matching of coordinates across these two
letters has total pair mutual information \(o(\lambda\log X)\).

Indeed, for a given row and group \(J\), take the surprise
\(-\log_X\Pr(\boldsymbol x_J)\) capped at its raw cost
\(\alpha_J\). Bin these scores into finitely many possibilities
at fixed truncation, with any desired small total upper
rounding error of order \(\lambda\). Each bin supplies a dictionary
of exponent at most its upper score cap. For the raw-cap last
bin this uses the raw box bound. Otherwise it follows by
counting the group vectors of probability at least the inverse
of the dictionary bound. If a score profile has total \(W\)
below \(L_U+L_V\) by a fixed positive multiple of \(\lambda\),
lock \(T\) with a smaller raw cost \(a\lambda\) and apply (CG1).
The profile's row count is then power-smaller. This holds also
for their union since the number of profiles is fixed at each
budget. Taking expectations, first of the capped scores and
then of the actual surprises, gives the entropy lower bound.
The raw group caps give the upper bound. Subtracting from
the sum of the singleton entropy bounds proves the matching
statement. The conclusion also measures aggregate correlations
of matched across-letter windows when the bounded-window
condition above holds.

The established balanced endpoint theory supplies a cheap
opposite modulus for every positive-height nonfinite-scale
letter. Thus (CG3) applies to at least two of the three
pairs at a positive-height balanced endpoint, and all three
when none of the letters is finite-scale. If the short height
vanishes, either tall letter provides the cheap opposite modulus.
The free letters themselves need no lower numerical height.
[Cheap moduli and finite-scale letters](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

If two letters supply arbitrarily cheap full moduli, use (CG2)
and the same information binning before choosing the two locks.
Freeze their raw root values at total cost at most \(2a\lambda\),
and for each group use its projection onto the unfrozen
positions. Projecting cannot increase dictionary cardinality.
For any allowed partition of the positions on all three letters,
\[
                       \sum_J H(\boldsymbol x_J)
                                    \ge (2/3-o(1))\lambda\log X.      \tag{CG4}
\]
This lower bound follows also from (CG3) by projecting every
group onto the two letters opposite one cheap lock. The two-lock
form (CG2) has a separate use when the charged conditional
dictionaries or lock records differ in the two applications.
Neither form gives the coefficient \(1/3\) for arbitrary
three-letter joint dictionary costs. In particular equality in
all pair entropy bounds can coexist with a three-way dependency
in abstract dictionaries. For instance, if independent
uniform labels \(u,v\) in a finite cyclic group determine
\(t=u+v\), every pair is uniform independent while the
triple's joint entropy is only twice the single-label
entropy. This example only describes the entropy obstruction;
it does not supply an arithmetic endpoint family.


## Conditional entropy under a prescribed lock

Let \(\mathcal R\) denote the full record of the frozen
coordinates in (CG1) or (CG2), with a prescribed
choice of its positions and uniform positive weight
and fullness margins. Fix the partition of the
remaining free positions in (CG1), or of all the
unfrozen positions in (CG2), into allowed groups
\(J\). On a uniform family of actual triples,
\[
                H(\mathcal R)+\frac12\sum_J
                H(\boldsymbol x_J\mid\mathcal R)
                         \ \ge\ (D-o(\lambda))\log X .        \tag{CG5}
\]
The error is in the iterated convention of (CG1).

For a given row, take the surprise of its record,
and the conditional surprise of each group in its
actual record, in base \(X\). Each score lies
between zero and \((D+o_X(1))\), since
probabilities refer to a finite uniform family
and its nonempty fibres. Bin all scores with a
fixed total upper rounding loss at most
\(\xi\lambda\). For fixed truncation the number
of score bins is bounded independently of \(X\).
A record surprise bin with upper cap \(s_0\) has
at most \(X^{s_0}\) different records. In each
of them the \(J\)-surprise bound \(s_J\) gives
a conditional dictionary with at most \(X^{s_J}\)
vectors, by normalization of conditional probability.
The corresponding subset of rows therefore has
count exponent at most \(s_0+\frac12\sum_J s_J+o(\lambda)\).
Whenever the realized score sum before binning
is below \(D-\zeta\lambda\), this gives a
power-smaller subset for small enough \(\xi\)
and sufficiently late budgets. Taking expectations
proves (CG5). At saturation in the one-lock case,
if \(H(\mathcal R)=o(\lambda\log X)\), this
recovers full aggregate conditional group entropy
on the two free letters for every allowed
partition that ignores the frozen positions.
