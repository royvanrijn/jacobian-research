# Sharp inverse degree for a parallel extension of the literal 58-variable map

**Theorem LiteralIndexFour58ParallelCeiling.** Let \(H\) be the
[literal 58-variable cubic](CUBIC_INDEX_FOUR_SYZYGY58.md) over any
characteristic-zero field. Adjoin any finite number of parallel
coordinates with homogeneous cubic components in all 58 old variables.
If the extended Jacobian has fourth power zero, the extended inverse
has degree exactly 207.

The [pre-\(S\) parallel restart theorem](CUBIC_INDEX_FOUR_COMPLETE_PARALLEL_RESTART.md)
applies to a family of cube-zero bases with a prescribed set of
permitted old variables. Here the specific literal cube-zero base and
its nine pre-\(S\) components and final \(S\) are fixed. The added
cubics may use every original base coordinate and the entire appended
block. The proof tests all high-weight cubics in that larger literal
space.

Let \(J=JH\) and write \(Q\) for the vector of new components. The
extended Jacobian and its fourth power are
\[
 M=\begin{pmatrix}J&0\\D Q&0\end{pmatrix},
 \qquad M^4=\begin{pmatrix}J^4&0\\D Q\,J^3&0\end{pmatrix}.
\]
Thus a scalar component \(q\) is admissible exactly when
\(dq\,J^3=0\). The old positive weights \(w_i\) give each
\((J^3)_{ij}\) weight \(w_i-w_j\) when that entry is nonzero.
For a weight-\(s\) component \(q_s\), the \(j\)-th entry of
\(dq_sJ^3\) has weight \(s-w_j\). Therefore different weight
components of a cubic cannot cancel one another in this condition.

For each nonempty weight space \(s>207\), take all ordinary cubic
monomials in the 58 coordinates of weight \(s\). Form one column
per monomial, indexed by a source index \(j\) and an ordinary
degree-eight monomial. Its entries are the coefficients of
\((dq_sJ^3)_j\). The
[exact checker](verify_literal_index_four_58_full_parallel_extension.py)
reconstructs \(J^3\) over \(\mathbb Q\) from the literal data.
It reduces coefficients modulo \(1\,000\,000\,009\), checking
the denominators are units, and performs exact column elimination
for each weight. All 128 nonempty spaces above weight 207 have full
column rank. Together they have 3,788 columns, including 129 at weight 209
and one at weight 621. These ranks
prove injectivity over \(\mathbb Q\) and hence over every
characteristic-zero field, since nonzero maximal minors survive
scalar extension.

Every admissible component consequently has weight at most 207.
With \(G=(I-H)^{-1}\), the extended inverse in old target
coordinates \(X\) and new ones \(Y\) is
\((G(X),Y+Q(G(X)))\). The inverse \(G\) preserves the old weights,
all at least one. Every monomial of \(Q(G(X))\) therefore has
ordinary degree at most 207. The old degree-207 component
remains present, and its nonzero third Jacobian row also prevents
a drop in nilpotency index. This proves the theorem for one
parallel layer over the literal map. It leaves iterated extensions
and the unrestricted weighted class \(U_{104}\) open.

## Size of a complete recursive grade-104 scan

**Proposition CubicIndexFourGrade104ScanSize.** The standard
recursive root-shuffle enumeration for grade 104 has at least
\(2\,829\,774\) candidate columns, and its complete loop enumerates
at least \(5\,211\,060\,703\,987\,899\,795\,082\,169\) basis
substitutions. These bounds already allow the usual quotienting
in all smaller grades.

Let \(d_j=\dim(A_{3,4})_j\). Exact rational reduction gives
\[
 d_0,\ldots,d_{12}
 =1,1,1,2,3,6,12,25,52,114,253,571,1304.
\]
The [four-dimensional cyclic Engel algebra](CUBIC_ENGEL_INDEX_FOUR_ONE_GENERATOR_OBSTRUCTION.md)
gives \(d_j\ge1\) at every larger grade. For a lower bound write
\(b_j=d_j\) through \(12\) and \(b_j=1\) for \(13\le j\le103\).
The multiset count for choosing lower-grade basis elements is
\[
 C(z,t)=\prod_{j=0}^{103}(1-tz^j)^{-b_j}.
\]
The usual recursive presentation first grafts triples of lower
basis elements of total grade 103. Their number is at least
\([z^{103}t^3]C\). Each root relation then takes a bottom basis
element, eight unordered off-path basis elements, and four new
internal vertices. Consequently the full loop visits at least
\[
 [z^{100}t^8]C(z,t)\sum_{j=0}^{100}b_jz^j
\]
substitutions. The same checker replays the grade-12 rational ranks
and computes the two coefficients by truncated integer
convolution. The substitution bound counts basis inputs to the
loop. It does not give their distinct nonzero row count or rank.

The [grade-104 problem](INDEX_FOUR_FINITE_WEIGHTED_TEST.md#the-grade-104-membership-problem)
concerns one weighted vector. To prove membership, a certificate
may use only a selected set of valid polarized root relations
and their grafted contexts. It suffices to give the exact rational
combination and verify its equality to the target. To prove
nonmembership by a dual, every relation must be killed, or one
must supply a finite Engel algebra with nonzero evaluation of the
target. The count above explains why deduplicating rows only
after a complete traversal of substitutions is insufficient to
make that traversal practical.

Replay from the repository root.

~~~sh
python3 research/jacobian/subtree-shuffling/verify_literal_index_four_58_full_parallel_extension.py
~~~
