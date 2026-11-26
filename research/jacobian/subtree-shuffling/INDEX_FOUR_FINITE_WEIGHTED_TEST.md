# A finite test for the sharp candidate at cubic index four

Let \(A_{3,4}\) be the symmetric ternary Engel algebra on one generator
with \(L_a^4(b)=0\), and let \(U_k\) be its weighted inverse class at
internal grade \(k\). The [58-variable example](CUBIC_INDEX_FOUR_SYZYGY58.md)
shows that any dimension-uniform inverse-degree ceiling is at least 207.
The following conditions are equivalent.

1. Every homogeneous cubic \(H\) with \((JH)^4=0\) over every
   characteristic-zero field has \(I-H\) invertible with inverse degree
   at most 207.
2. Every \(U_k\) with \(104\le k\le206\) is zero.
3. Every finite admissible restart tower over the retained
   exact-index-four fixed-coordinate scaffold has inverse degree at most
   207. The equivalence also holds with complete-basis saturation in
   every round.

By the [inverse-degree dictionary](SUBTREE_SHUFFLE_AND_INVERSE_DEGREE.md#inverse-degree-interpretation),
condition 1 is equivalent to \(U_k=0\) for every \(k\ge104\).
The [quadratic inverse gap](CUBIC_CUBE_ZERO_SHUFFLE_EXTENSION.md#2-a-quadratic-propagation-interval)
says the block \(U_K,\ldots,U_{2K-2}=0\) propagates to every later grade,
at any fixed arity and Engel length. Its premise at \(K=104\) is exactly
condition 2. The [all-round equivalence](SUCCESSIVE_RESTART_UNIVERSALITY.md#equivalence-of-the-all-round-and-unrestricted-uniform-bounds)
gives condition 3, including the complete-basis assertion.

If any test in the finite block fails, the graded quotient through
that grade produces an exact-index-four triangular map with inverse
degree \(2k+1\), between 209 and 413. Exact index follows because the
grade-three comb survives below the first Engel relations and is
retained by the truncation. If every test succeeds, the existing
degree-207 example makes the ceiling sharp. This finite block tests
the specific candidate ceiling 207. Finding any larger finite uniform
ceiling, or proving the inverse degrees unbounded, remains a separate
question.


## The grade-104 membership problem

The first obligation in this block is \(U_{104}=0\) in the rational
one-generator quotient. Its failure produces a degree-209 triangular
example. Its success alone leaves grades 105 through 206 to test.

The closest small-index tree theorem is [Wright, *The Jacobian
Conjecture: ideal membership and a combinatorial vanishing
conjecture*, §3.8, Theorem 3.19](https://arxiv.org/pdf/math/0511214).
It gives inverse degree at most 15 for homogeneous cubic \(I-H\) with
\((JH)^4=0\) and a **symmetric Jacobian matrix**. Its tree module uses
the scalar gradient potential. The symmetric ternary multiplication
in \(A_{3,4}\) comes from polarization of arbitrary cubic coordinates
and supplies no symmetric-Jacobian hypothesis. Wright's vanishing
interval thus does not give any grade of the present block.

An unreduced unordered-tree enumeration already has
\(27\,600\,147\,135\,309\,968\,703\,380\,177\,022\,909\,343\,380\,105\,889\)
columns in grade 104. The count follows from the cycle index for the
three unordered children at the root,
\[
T(z)=1+\frac{z}{6}\big(T(z)^3+3T(z)T(z^2)+2T(z^3)\big).
\]
The constant term represents the sole leaf. The identity, the three
transpositions, and the two 3-cycles give the three terms. The
unreduced count bounds the size of literal enumeration, while the
recursive quotient replaces the factor spaces by their lower-grade
quotients before grafting. Exact small replay with
`verify_subtree_shuffle.py --scan 3 4 12` gives quotient dimension
1304 out of 3057 unordered shapes at grade 12 and a dual evaluation
\(188811\) for the weighted inverse class. A 35-second capped attempt
through grade 15 did not finish. These figures provide no estimate
for the dimensions of the recursive quotient at grade 104.

A rational row combination for the recursively reduced grade-104 vector,
or a rational dual vanishing on every grade-104 relation and nonzero on
the vector, would decide that obligation. A finite-field calculation can
certify vanishing by proving full column rank in an integral presentation,
as in the cubic cube-zero computation. For a general rank-deficient
matrix, membership at a single prime does not certify rational membership.
A modular nonmembership can certify rational nonmembership when the
relation matrix is independently shown to have the same rank over
\(\mathbb Q\) and the residue field. The augmented matrix then has
rational rank at least its larger modular rank. The grade-104 class is
currently undecided.
