# A 48-variable cubic map with cube-zero Jacobian and inverse degree 29

The rational map specified in [cubic_cube_zero_sharp48.json](cubic_cube_zero_sharp48.json)
has 48 variables and 1,213 nonzero cubic monomials. It gives a smaller
explicit sharp example for the cubic cube-zero inverse-degree theorem.
The proof below uses only this literal polynomial data and its exact
derivatives; it does not import a shuffle quotient calculation.

**Theorem CubicCubeZeroSharp48.** There is a homogeneous cubic polynomial map
H over Q in dimension 48 such that (JH)^3=0, (JH)^2 is not zero, and I-H has
a polynomial inverse of degree exactly 29. The same rational formulas have
these properties over every characteristic-zero field.

## Definition of the map

Use coordinates y_0,...,y_47. The input lists their nonnegative grades g_i
in increasing order. The dimensions of the successive grades 0 through 14 are

    1,1,1,1,2,3,4,6,9,10,4,3,1,1,1.

A record with output o, inputs [i,j,k], and coefficient c defines the term

    c*y_i*y_j*y_k

in H_o. Repeated indices denote ordinary powers. The listed coefficient
already includes every multiplicity; no additional multinomial factor is
inserted when interpreting this file. Distinct records are distinct
commutative monomials, and all coefficients are nonzero exact rational
numbers. The data are the complete definition of H.

Every recorded term satisfies

    g_o=1+g_i+g_j+g_k.                                  (1)

In particular each input grade is strictly below its output grade.
Set w_i=2g_i+1. These weights are positive integers at most 29, and (1) gives

    w_o=w_i+w_j+w_k.                                    (2)

## Polynomial inverse and its degree

Order the coordinates by grade. Equation (1) makes I-H strictly triangular
with identity diagonal, so its inverse G exists by recursive polynomial
substitution. Equation (2) shows inductively that G_i is weighted
homogeneous of weight w_i. Since all variable weights are at least 1,
every ordinary monomial of G_i has ordinary degree at most w_i<=29.
Thus deg G<=29.

Let e_i denote the standard coordinate vectors. Specialize the inverse
input to s*e_0. Weighted homogeneity gives

    G_i(s*e_0)=v_i*s^(2g_i+1).

The constants v_i are obtained directly from G=x+H(G):
v_0=1, all other grade-zero inputs vanish, and at each positive grade the
coefficient v_o is the sum of the literal listed c*v_i*v_j*v_k.
The exact recurrence gives

    v_47=8803690708443294969/95329784685296 !=0.           (3)

The last coordinate has grade 14 and weight 29. Hence the actual inverse
G_47 contains the monomial v_47*x_0^29. Equation (3) and the upper bound
prove deg G=29.

**Lemma Sharp48LineCriterion.** Under the displayed strict grade and weight
profile, deg G=29 if and only if the x_0^29 coefficient of G_47 is nonzero,
equivalently the vector of weight-29 coefficients on the input line s*e_0
is nonzero.

Indeed x_0 is the unique weight-one input, coordinate 47 is the unique
weight-29 output, and every other output has weight below 29. A monomial
of ordinary degree 29 and weight at most29 must therefore be x_0^29,
in output47. Consequently absence of that coefficient rules out every
ordinary degree 29 term. The positive-weight upper bound already rules
out all higher ordinary degrees. This proves both directions, rather
than only a sufficient line test.

**Corollary Sharp48UniformLowerBound.** Every dimension-independent
inverse-degree upper bound C for all cubic homogeneous cube-zero Jacobians
satisfies C>=29, by this particular 48-dimensional map. Combined with the
[separate universal upper bound 29](CUBIC_CUBE_ZERO_SHUFFLE_EXTENSION.md), the least uniform upper bound is 29.
This implication does not assert that the abstract lower bound C>=29
forces the existence of this particular 48-dimensional realization.

## Exact Jacobian identity

The checker differentiates the literal 1,213 monomials, using the actual
multiplicity of a repeated variable. It obtains a polynomial Jacobian
with 749 nonzero entries and 3,198 monomial terms. Ordinary sparse polynomial
matrix multiplication over Q gives:

    JH squared: 551 nonzero entries, 3,219 monomial terms;
    JH cubed: the zero polynomial matrix.

The square and cube multiplication stages perform respectively 27,061 and
12,098 exact term products, with full coefficient collection and rational
cancellation. These are all-variable polynomial identities, not evaluations
at a finite set of points.

The same literal Jacobian also satisfies

    (JH(e_0))^2 e_0=9e_2 !=0.                           (4)

Thus its polynomial nilpotency index is exactly 3. The identity and
nonvanishing in (3)--(4) persist over every characteristic-zero field.
The same direct data give
JH(e_0)*JH(e_0)*JH(e_1)*e_0=27e_6, so strong cube-nilpotency fails.
Also (JH)_(1,0)=3y_0^2 while (JH)_(0,1)=0, and H_2=3y_0^2y_1
is not a cube of a linear form. These controls locate the example outside
the symmetric, strongly cube-nilpotent and power-linear special classes.

This proves the theorem.

## Reproduction

Run the standard-library checker with its input and source ledger beside it:

    python3 research/jacobian/subtree-shuffling/verify_cubic_cube_zero_sharp48.py \
      --output /tmp/cubic_cube_zero_sharp48_verification.json

The replay takes less than one second in the current environment. It
validates the 48 coordinates, grades, every index, all canonical nonzero
cubic monomials, strict triangularity, and the positive weight bound.
It then recomputes the full Jacobian cube, the exact point square in (4),
and all inverse-line coefficients up to (3). A valid single-coefficient
perturbation is required to produce a nonzero literal Jacobian cube.
Malformed grades, indices, zero coefficients and duplicate monomials are
also rejected.

The [source ledger](cubic_cube_zero_sharp48_source_ledger.json) pins this
proof, the literal data and the checker. The
[receipt](cubic_cube_zero_sharp48_verification.json) records the exact
replay. Optimized Python is refused. No CAS, quotient producer, source
algebra or external theorem is needed for this witness check. No claim
of minimum possible dimension is made.

Updated 20 November 2026.

