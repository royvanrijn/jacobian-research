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
