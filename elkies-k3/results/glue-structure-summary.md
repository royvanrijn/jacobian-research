# Glue decomposition of the Elkies rank-17 MW lattice

Let M be the recovered rank-17 positive-definite Mordell-Weil lattice
of determinant 948.

A saturated rank-9 sublattice K was found from a maximal +2 clique.
It is an index-2 overlattice of 2A9 with:

- rank 9
- determinant 1280
- 160 signed norm-4 vectors
- automorphism group order 161280

Its orthogonal complement C in M has:

- rank 8
- determinant 303360
- Smith invariants [1,1,2,2,2,2,4,4740]
- 4 signed norm-4 vectors
- automorphism group order 8

The orthogonal direct sum K + C has index

    [M : K + C] = 640

with quotient

    M / (K + C) ~= (Z/2)^5 x Z/20.

Of the 640 glue cosets, 315 contain norm-4 vectors.

The 2622 signed norm-4 vectors decompose by glue order and
projection norms as follows:

| order | K norm | C norm | signed |
|------:|-------:|-------:|-------:|
| 1 | 0 | 4 | 4 |
| 1 | 4 | 0 | 160 |
| 2 | 1 | 3 | 56 |
| 2 | 2 | 2 | 288 |
| 2 | 3 | 1 | 160 |
| 4 | 5/4 | 11/4 | 100 |
| 4 | 9/4 | 7/4 | 252 |
| 10 | 4/5 | 16/5 | 36 |
| 10 | 6/5 | 14/5 | 162 |
| 10 | 9/5 | 11/5 | 64 |
| 10 | 11/5 | 9/5 | 504 |
| 10 | 14/5 | 6/5 | 60 |
| 20 | 21/20 | 59/20 | 66 |
| 20 | 29/20 | 51/20 | 210 |
| 20 | 41/20 | 39/20 | 180 |
| 20 | 49/20 | 31/20 | 320 |

Total:

    2622 signed norm-4 vectors
    = 1311 +/- pairs.

The denominator of the projection norms tracks the glue order:
integral for order 2, quarters for order 4, fifths for order 10,
and twentieths for order 20.

Thus the short-vector population is largely controlled by the glue
between K and C rather than by the component lattices separately.
