# Smaller cube-zero witnesses from the explicit 48-variable map

Use the [explicit rational cubic](CUBIC_CUBE_ZERO_SHARP48.md) and its
[exact inverse-line receipt](cubic_cube_zero_sharp48_verification.json).
The receipt records a nonzero vector at every grade from 0 through 14.
For example, the first nonzero coordinate in each recorded vector has
zero-based index
\[
             0,1,2,3,4,6,9,13,19,28,38,42,45,46,47 .
\]
These are the actual coefficients in the inverse at \(t e_0\), grouped by
weight \(2k+1\).

## Smaller triangular witnesses

For \(0\le k\le14\), keep the coordinates of grades at most \(k\).
The retained components use only coordinates of strictly smaller grade,
so projection intertwines the full map and the truncation. It also
intertwines their inverses by triangular uniqueness. The Jacobian on
the prefix is the leading principal block of a block lower triangular
Jacobian, so its cube is zero. For \(k\ge2\) the existing
\((JH(e_0))^2e_0=9e_2\) witness survives and gives exact index three.

The upper bound on inverse degree of the prefix is \(2k+1\) by the
positive weights. Its nonzero grade-\(k\) coefficient at \(t e_0\)
gives equality. In increasing order of \(k\), the dimensions of these
prefixes are
\[
              1,2,3,4,6,9,13,19,28,38,42,45,46,47,48.
\]
In particular, for every odd inverse degree from 5 through 29 there is
an explicit rational exact-index-three witness of the indicated
dimension. The same coefficients work over every characteristic-zero
field. Identity padding gives witnesses in every larger dimension.

## Back to the universal length-three calculation

By the [inverse-coefficient dictionary](SUBTREE_SHUFFLE_AND_INVERSE_DEGREE.md#inverse-degree-interpretation),
a nonzero degree-\(2k+1\) inverse term in a map with \((JH)^3=0\)
witnesses \(U_k\ne0\). Thus this one certified map supplies the full
nonvanishing side \(1\le k\le14\) of the length-three threshold.
The separately proved universal upper bound \(29\) supplies the
vanishing side, by the reverse graded-realization implication.
The same lower witness can be padded by identity coordinates to
dimension 195 and gives the earlier sharp-degree statement by an
alternative route. The exact whole-grade dimensions and upper-bound
certificates retain their own dependencies.
