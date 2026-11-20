# Sharp inverse degree for cubic maps with cube-zero Jacobian

For every characteristic-zero field K and every cubic homogeneous polynomial
map H:K^n->K^n satisfying (JH)^3=0, the map I-H is a polynomial
automorphism and its inverse has degree at most 29. This bound is attained
over Q in dimension 195. The proof combines exact finite shuffle
certificates with the classical quadratic recurrence for inverse coefficients.

The finite certificates are checked by [verify_cubic_cube_zero_shuffle.py](verify_cubic_cube_zero_shuffle.py),
using only the Python standard library. The [input](cubic_cube_zero_shuffle_certificates.json.gz) contains the exact finite certificates. The retained orbit reduction
and finite graded inverse realization are in
[SUBTREE_SHUFFLE_AND_INVERSE_DEGREE.md](SUBTREE_SHUFFLE_AND_INVERSE_DEGREE.md).

## 1. Statement and classical scope

**Theorem CubicCubeZeroSharpInverseDegree.** Let K be any field of
characteristic zero, let n>=1, and let H in K[X_1,...,X_n]^n have homogeneous
cubic components, allowing zero components. If (JH)^3=0 as a polynomial
matrix, then I-H has a polynomial inverse of degree at most 29.
There is a homogeneous cubic H over Q in dimension 195 for which JH has
nilpotency index exactly 3 and the inverse degree is exactly 29.

**Theorem CubicLengthThreeEventualShuffle.** In the symmetric one-generator
ternary quotient A_(3,3) defined by polarized L_a^3(b)=0, where
L_a(b)=mu(b,a,a), the orbit-size-weighted total U_k is nonzero for
1<=k<=14 and is zero for every k>=15. Equivalently, the planar constant
belongs to the length-three subtree-shuffle span exactly for k>=15 among
positive internal grades.

These are fixed-nilpotency-index statements. Bisi, Dyszewski, Gantert,
Johnston, Prochno and Schmid, arXiv:2301.08221v3, Conjecture 1.3
(printed page 5), asks for an eventual constant-vector threshold for each
fixed arity d and shuffle length p. Their Theorem 1.4 would require all
p at one d>=3 to conclude the full Jacobian conjecture. Only d=p=3 is
settled here.

Several stronger bounds under additional hypotheses must be distinguished.

* Wright, arXiv:math/0511214v2, Theorem 3.2 (printed page 8), assumes a
  symmetric Jacobian and gives inverse degree <=2d-1 when (JH)^3=0.
* Yan--de Bondt, arXiv:1302.5864v3, Theorem 5.3(ii) (printed pages 11--12),
  additionally assumes dim(ker JH)=dim_K(ker_K JH). For homogeneous
  power-linear maps this applies and their inverse formula gives degree
  at most d^2.
* Adamus--Bogdan--Crespo--Hajto, arXiv:1508.02012v1, Theorem 3
  (printed pages 3--4), assumes that every component of H is a cube of a
  linear form and gives inverse degree at most 9.
* The bound d^(p-1) discussed by Bisi et al. in their related-work section
  assumes strong nilpotency: products JH at different points vanish.
  The single-point identity (JH)^p=0 used here is weaker.

The general recurrence used below is classical: Zhao,
arXiv:math/0403020v2, Theorem 2.4 and Proposition 2.5, equations (2.6)
and (2.10), printed pages 7--9. Symmetry is introduced only in his
Section 3 and is not a hypothesis of these statements. The paper explicitly
notes that the formulas and proofs hold over every Q-algebra.
The application to the present quotient and its finite certificates is
separate from this classical recurrence.

## 2. A quadratic propagation interval

Let H be any formal map of order >=2 over a characteristic-zero field.
Set F_t(x)=x-tH(x) and write its unique formal inverse as

    G_t(x)=x+tN_t(x)=x+sum_(k>=1) t^k N_k(x).

The inverse identities give N_t(F_t(x))=H(x). Differentiate with respect
to t and compose on the right with G_t. The formal chain rule gives

    partial_t N_t = (J N_t) N_t.

Comparing coefficients, for every k>=2,

    (k-1) N_k = sum_(i+j=k, i,j>=1) (J N_i) N_j.          (Q1)

For d-homogeneous H, N_k is homogeneous of degree (d-1)k+1.
This is a formal calculation; it uses neither analytic convergence nor a
gradient or nilpotency hypothesis.

**Lemma QuadraticInverseGap.** If N_K,...,N_(2K-2) are zero identically for
K>=2, then all N_k with k>=K are zero.

At k=2K-1, each i+j=k has i>=K or j>=K, and both i,j<k.
If N_i=0 then JN_i=0. If N_j=0 the product is zero. Equation (Q1),
division by k-1, and induction prove the claim.

**Corollary ShuffleQuadraticGap.** In A_(d,p), if
U_K,...,U_(2K-2)=0 for K>=2, then U_k=0 for all k>=K.

Evaluation of the universal tree coefficients in any finite-dimensional
symmetric d-ary algebra satisfying L_a^p(b)=0 gives the homogeneous inverse
coefficients N_k of its map I-H, H(a)=mu(a,...,a). The hypothesized U
identities therefore give the N identities as polynomial maps, so the
lemma applies in every such algebra. If some later U_k were nonzero, the
finite graded truncation A_(d,p)/(grades>k) would be such an algebra and
would have nonzero degree-(d-1)k+1 inverse coefficient at its generator.
This contradicts the lemma.

Thus for d=p=3 and K=15 it suffices to check grades 15 through 28.
The elementary ternary grafting interval 15 through 43 is unnecessary for
this inverse-coefficient conclusion.

## 3. Exact finite quotient calculation

Let A_k be the grade-k piece of A_(3,3), where grade counts internal
vertices. The certificates establish

    dim A_0,...,dim A_18 =
    1,1,1,1,2,3,4,6,9,13,19,26,34,40,35,16,13,9,2.

They give U_k!=0 for k=1,...,14, and U_k=0 for k=15,...,18.
They also establish A_k=0 for each k=19,...,28, using the same proof prime
ell=1000000009. The final zero-piece certificate has 37406 columns and
37406 independent valid root relations.

In each grade, first graft three previously established quotient basis
elements whose grades sum to k-1. Their unordered multisets form the
prequotient basis. This is the symmetric cube of the already reduced
lower pieces: every relation occurring below the root has already been
imposed. The remaining root relations are the polarized length-three
shuffle sums, with six off-path basis elements and one bottom basis
element. Their multilinearity makes substitutions of basis elements
sufficient. This is the retained quotient-completeness argument; no
unlisted relation type or additional identity is used.

For grades 1 through 18 the checker regenerates EVERY root relation over Q.
The input gives a reduced row basis R and a list S of selected original
rows. Let J be the pivot columns, B=S[:,J], and let C be the selected
free-column block. The following exact conditions are checked.

1. Each pivot row has coefficient 1 in its pivot, and every other nonzero
   term is in a free column strictly below that pivot index.
2. Every regenerated original row reduces to zero by R.
3. The square selected pivot-column block B has full rank modulo ell.
4. All reduced rational coefficients have denominator prime to ell.

Since S=B R and B is invertible over Q, R lies in the original relation
space, and condition 2 gives the converse inclusion. Thus this is an
exact rational quotient certificate. The same conditions also identify
the ACTUAL modular relation space in each initial grade: reduced
coefficients specialize, all original rows lie in the specialized R-span,
and the selected B remains invertible at this SAME ell. This last
invertibility check is essential; denominator invertibility alone would
not justify recursive base change.

The weighted U target is regenerated independently by the ternary inverse
recursion, including its multinomial factors. Its exact rational remainder
is nonzero through grade 14 and zero at grades15 through18. The coordinate
dual at grade 14 has the nonzero value

    -10529123457938566680915 / 34700041625447744.

For each grade19 through28, the certificate lists only selected root
substitutions. The checker validates the six off-path entries and bottom
entry against the actual earlier modular bases, checks the total grade,
regenerates the root sum, and verifies that each selected row increases
rank. Their number equals the ENTIRE prequotient dimension, proving
full rank. It is unnecessary to enumerate the unused root relations.

A previously proved zero modular piece may be replaced by zero while
forming later root sums. The checker requires its grade to be strictly
lower than the current grade and its established basis to be empty.
This shortcut is an induction step, not an assumption of eventual
vanishing.

## 4. Integral model and the passage from modular zeros to Q

Here are the base-change details behind the full-rank certificates.
Let R_0=Z[1/6!] and let F_(R_0) be the free symmetric ternary algebra on x,
with unordered trees as its basis and grafting as its product. Give each
tree its internal grade. Impose the ideal generated by the full
six-variable polarization of L_a^3(b).

At any fixed grade the free tree module is finite free and the relation
submodule has a finite matrix presentation: it is generated by root
polarized substitutions and their tree contexts, all within that grade.
The unlabelled root sum in the checker differs from a labelled
polarization value only by a product of factorials of repeated off-path
subtree multiplicities. This product divides 6!, hence is a unit in R_0.
The ternary symmetrization factors divide powers of 3!, also units.
Consequently the normalized root presentation is equivalent to the
polarized Engel presentation over R_0, over Q, and over F_ell.

Formation of this quotient commutes with base change. Indeed the
multilinear substitution and context generators can be expanded in the
tree basis, and the quotient of a finite free module by their span
commutes with tensor product. The recursive prequotient computation
commutes as well: its lower relations have been verified in both
fields, and its chosen pivot blocks are units at ell. Thus the initial
certificates really construct A_(R_0) tensor F_ell in grades<=18, and
the subsequent full-rank steps construct it in grades19,...,28.

For any finite relation matrix over R_0, its rank over Q is at least its
rank modulo ell: a nonzero minor modulo ell is a nonzero rational minor.
It follows that

    dim_Q A_k <= dim_(F_ell) (A_(R_0) tensor F_ell)_k.

In particular, the modular whole-piece zeros at grades19,...,28 prove
A_k=0 over Q. This argument never infers rational membership from a
modular zero of U alone.

## 5. Completion, equivalence and sharpness

**Theorem CubicCubeZeroUniformBoundEquivalence.** The following are equivalent:

1. In every dimension, every cubic homogeneous I-H over every characteristic-zero
   field with (JH)^3=0 has a polynomial inverse of degree at most 29.
2. U_k=0 in A_(3,3) for every k>=15.
3. U_k=0 in A_(3,3) for 15<=k<=28.

The quadratic gap proves 3=>2, and 2=>3 is immediate. Evaluation of universal
tree coefficients proves 2=>1. For 1=>2, if any U_k with k>=15 were nonzero,
the finite graded rational truncation at k would be a cube-zero cubic map
whose inverse degree is 2k+1>=31, contradicting1. This is a genuine reverse
implication; it does not identify a fixed nilpotency-index statement with
the unrestricted Jacobian conjecture.

Sections 3--4 give U_15=...=U_28=0 over Q. The quadratic gap lemma gives
U_k=0 for every k>=15. Rational tree identities remain true after extension
to any characteristic-zero field. Polarizing an arbitrary cubic H over
that field gives a symmetric ternary product with JH(a)=3L_a.
The hypothesis (JH)^3=0 therefore imposes exactly the Engel identity used
in the quotient. The formal inverse coefficients evaluate from the U_k,
so only grades<=14 remain. Their ordinary degrees are 2k+1, hence at most
29. Formal inverse identities then give two polynomial inverse identities.

For sharpness take

    B=A_(3,3)/(all grades>14).

The exact basis dimensions give dim_Q B=195. This is a finite-dimensional
graded symmetric ternary algebra, and H(y)=mu(y,y,y) is a cubic polynomial
map on its coordinate space. Grafting increases grade by one plus the
sum of child grades, so each output coordinate of H depends only on
strictly lower-grade input coordinates. Thus I-H is triangular with
identity diagonal.

The polarized identities give (JH)^3=27L_y^3=0. Its index is exactly 3:
the grade-two comb L_x^2(x) is nonzero since there is no relation below
grade 3, so (JH(x))^2(x)=9L_x^2(x)!=0. The inverse coefficient at grade 14
evaluated at the generator is U_14!=0. All products with more than 14
internal vertices vanish in B for grading reasons. Hence the inverse
degree is exactly 29.

This proves both stated theorems. No assertion of whole-quotient vanishing
after grade28 is required or made.

## 6. Reproduction and assurance

Run from the repository root:

    python3 research/jacobian/subtree-shuffling/verify_cubic_cube_zero_shuffle.py \
      --output /tmp/cubic_cube_zero_shuffle_verification.json

The standard-library replay takes about 30 seconds in the current environment.
The compressed input is 684,646 bytes; the uncompressed JSON is 3,987,765 bytes.
It contains no opaque CAS result and no generated polynomial in 195 named
coordinates. It stores rational pivot rows and selected genuine shuffle
substitutions; the checker reconstructs their source relations and all U
coefficients from the retained algebra definition.

The exploratory producer used FLINT's exact Dixon solve to obtain rational
pivot certificates. That CAS is not needed for checking them. Independent
replay and code review supply internal mathematical evidence, not formal
verification or external review. Historical novelty is not inferred merely
from the literature search.


Updated 20 November 2026. The [source ledger](cubic_cube_zero_shuffle_source_ledger.json) pins the retained proof inputs and five primary papers; the [receipt](cubic_cube_zero_shuffle_verification.json) records the exact replay.
