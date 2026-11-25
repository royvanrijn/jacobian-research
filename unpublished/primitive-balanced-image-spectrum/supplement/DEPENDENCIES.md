# Proof and verification dependencies

The proof of the main safety theorem uses three cases for the leading binary
form. The ordinary-reduction theorem is a separate result proved for the
stated geometric field action and signature. It is applied to the cyclic
factors appearing in the third case.

| Manuscript argument | Main inputs and where they are used |
| --- | --- |
| Transfer from contractions to the sphere | The kernel identity agrees with van den Essen, Wright and Zhao, Theorem 3.3. The proof is included using the action of \(\xi_i-\partial_{z_i}\) on monomials. Angular averaging and beta integration give the factorial transfer. Homogenization and differentiation handle the degree restriction and arbitrary fixed markers. [On the Image Conjecture](https://arxiv.org/abs/1008.3962v2) |
| Angular core and quadratic Laurent height | Müger–Tuset, Theorem 2.1, excludes a nonzero one-variable polynomial if all of its interval moments vanish. Apply it after a linear change of interval and an embedding of the finitely generated characteristic-zero coefficient field. The quadratic-height endpoint argument and polynomial shear are proved in the manuscript. [On the moments of a polynomial in one variable](https://arxiv.org/abs/1910.12148v1) |
| Joint endpoint construction | Long, Theorem 3.1, gives the quartic Hopf seed and its marked coefficient identity. The manuscript derives the joint endpoint identity, the factorial value for the common marker, the radial evaluation and the simple-root specialization. The finite replay checks samples of those identities in the original variables and in a separate Laurent expansion. [The Mathieu Property for Compact Connected Lie Groups](https://arxiv.org/abs/2609.16178v1) |
| Small leading root multiplicities | Projective adjunction, Weierstrass preparation and resolution of the pencil supply a proper phase map with regular top form. The proper-volume lemma in the manuscript combines positivity for the graph generator, proper direct image of the filtration, strictness, and Fourier transform. The graph input is Chen–Dirks–Mustaţă, Example 5.19, and the proper direct-image input is their Theorem 8.1. Their Proposition 3.19 and Corollary 3.20 give the strictness used in the quotient argument. [An introduction to \(V\)-filtrations](https://arxiv.org/abs/2402.07816v1) |
| Leading multiplicity equal to degree | The angular core of the top height is a nonzero power of a coefficient. The prime congruence isolates it in the interval denominator. The manuscript explains specialization to a number field and preservation of the required leading identity and nonzero coefficients. |
| Leading multiplicity above the degree | The manuscript derives the finite lattice bounds and treats one-sided, core, two-node, and higher-node faces. For the remaining faces, the normalized Kummer eigensheaf calculation identifies the isolated coefficient with a Frobenius scalar. The finite supporting-line and branch arithmetic is reproduced by the face checker. The written proof supplies the reduction, the Frobenius interpretation, and the collision shears. |
| Small-dimensional ordinary factors | Sawin, Theorem 1, Corollary 2 and Theorem 3, give density one after finite extension for surfaces and hence for elliptic curves by squaring. Fité, Theorem 1, gives positive density for threefolds with a geometric imaginary quadratic action. Where a threefold occurs with a smaller factor, apply Fité over the extension where the smaller factor is ordinary on a density-one set. Intersect with that set and the density-one absolute-degree-one set. [Sawin](https://arxiv.org/abs/1506.04784v2) and [Fité](https://arxiv.org/abs/2005.10185v3) |
| Rank-two ordinary-reduction theorem | Weil bounds for weight-one Frobenius, slope and duality properties of the split crystalline factors, and the degree-polynomial comparison for actual endomorphisms. The relative trace norm follows the approach in Cantoral Farfán–Li–Mantovan–Pries–Tang, Proposition 3.12. Adjoining torsion and roots of unity adapts the argument recorded by Ogus in Proposition 2.7.2, printed pp. 371–372. The manuscript proves the two-place divisibility for a nonordinary balanced pair and the contradiction at equality. [Ogus](https://publications.ias.edu/sites/default/files/Number45.pdf) and [Cantoral Farfán–Li–Mantovan–Pries–Tang](https://arxiv.org/abs/2508.11174v1) |

## Branch data for the largest cover orders

An entry \(d:g\) denotes the factor of conductor \(d\) and dimension \(g\).
Only positive-dimensional factors are listed. Branch multisets contain the
two endpoint multiplicities and the multiplicities at the roots of the face
polynomial, including their ramification multiplier.

| Order | Two distinct roots | Factor dimensions | Double root | Factor dimensions |
| --- | --- | --- | --- | --- |
| 7 | \(4,2,2,6\) | \(7:6\) | \(4,4,6\) | \(7:3\) |
| 9 | \(8,2,2,6\) | \(3:1,\ 9:6\) | \(8,4,6\) | \(9:3\) |
| 10 | \(5,3,3,9\) | \(2:1,\ 5:2,\ 10:4\) | \(5,6,9\) | \(10:2\) |
| 11 | \(7,3,3,9\) | \(11:10\) | \(7,6,9\) | \(11:5\) |
| 15 | \(10,4,4,12\) | \(3:1,\ 5:2,\ 15:8\) | \(10,8,12\) | \(15:4\) |

For two distinct roots, the primitive factors in this table have rank two
over their cyclotomic fields and satisfy the signature bound. Each occurs
as an isogeny factor cut out by a rational deck-action idempotent, and the
deck action preserves the Jacobian polarization. The double-root primitive
factors have rank one over their cyclotomic fields. At split good absolute
degree-one primes their crystalline factors have height one. The product
of the ordinary factors is ordinary, since the slopes on an isogeny product
are the multiset union of the slopes on its factors.

## Scope of the replays

The finite branch computation tests the proposed list of supporting lines,
all admissible node occupations on those crossing lines, and all root
partitions for each resulting face degree. Normalized eigenspace dimensions
and conductor factors are tested in the cases where their branch hypothesis
holds. The remaining branch partitions are routed to the binomial or
polynomial-shear arguments proved in the paper. The computation assumes the
mathematical reduction to the model lattice and the Kummer dimension
formula. The portable witness computation uses finitely many orders, and
its symbolic input is the integer family itself. The paper supplies the
all-order proof. Neither replay uses the optional stripped-chart Gröbner
calculations from the separate research notes.
