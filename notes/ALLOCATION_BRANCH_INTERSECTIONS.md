# Hensel products and allocation-branch intersections

Fix a collision polynomial

\[
M_0=\prod_{\rho=1}^{\ell}(W-r_\rho)^{m_\rho}
\]

with distinct roots, and two allocations in the same component,

\[
A=((i_\rho,j_\rho))_\rho,
\qquad
A'=((i'_\rho,j'_\rho))_\rho,
\]

where

\[
2i_\rho+3j_\rho=2i'_\rho+3j'_\rho=m_\rho,qquad
\sum i_\rho=sum i'_\rho,qquad
\sum j_\rho=sum j'_\rho.
\]

There are unique integers `k_rho` such that

\[
i_\rho-i'_\rho=3k_\rho,qquad
j_\rho-j'_\rho=-2k_\rho.                         \tag{1}
\]

The global degree constraints give

\[
\sum_\rho k_\rho=0.                              \tag{2}
\]

Thus every transfer in one direction is compensated at other collision
roots.

## Local transfer blocks

For `k>=1`, let `Z_k` be the completion along `U=S^3,V=S^2` of

\[
U^2=V^3,qquad \deg U=3k,quad\deg V=2k,quad\deg S=k.
\]

The cases currently known are

\[
\mathfrak Z_1=k[[s,\epsilon]]/(\epsilon^2)
\]

and

\[
\mathfrak Z_2=
k[[p,q,X,Y]]/(X^3,2XY-pX^2,Y^2-qX^2).
\]

They are finite flat of ranks `2` and `4` over their reduced factor spaces.

## Hensel-product theorem for strong equality

Let `D^=` denote the completed correspondence of the two allocations with the
strong equation

\[
M_A=M_{A'}.
\]

Then

\[
\boxed{
D^=\simeq R_{\mathrm{com}}
\widehat\otimes
\bigwidehat\otimes_{\rho:k_\rho\ne0}
\mathfrak Z_{|k_\rho|}.}                         \tag{3}
\]

Here `R_com` is a regular complete local ring containing the root-position
parameters and the deformations of factor blocks on which the two allocations
agree.

### Proof

Choose pairwise coprime discs around the distinct roots `r_rho`.  Formal
Hensel factorization uniquely decomposes every monic deformation of `M_0` as
a product of monic cluster polynomials of degrees `m_rho`.  Equality of the
two global monic products is therefore equivalent to equality of the
corresponding cluster polynomials root by root.

If `k_rho=0`, the two allocations of that cluster agree and its correspondence
is the diagonal, hence smooth.  If `k_rho!=0`, canceling the common double and
triple factors leaves

\[
U_\rho^2=V_\rho^3
\]

with degrees `3|k_rho|,2|k_rho|`; this is `Z_|k_rho|`.  Independence of the
Hensel clusters gives the completed tensor product (3).

The incidence equation `Phi=0` has nonzero derivative in a root-position
direction on a dense open subset of every exact collision stratum.  Formal
implicit-function elimination removes one smooth parameter from `R_com` and
does not alter the Artin multiplication of the transfer blocks.

## The affine equalizer

The normalized coefficient map forgets the constant and linear coefficients.
The actual branch correspondence `D^aff` is consequently defined by

\[
M_A-M_{A'}=\lambda W+\mu.                         \tag{4}
\]

There is a closed immersion

\[
D^=\hookrightarrow D^{\mathrm{aff}}.             \tag{5}
\]

Equation (4) is global: unlike strong equality, it does not formally split
root by root.  This is the sole coupling between the local transfer blocks.

### Equalizer criterion

Assume the local blocks in (3) are finite flat over their reduced factors and
put

\[
L(A,A')=\prod_{\rho:k_\rho\ne0}
\operatorname{rank}\mathfrak Z_{|k_\rho|}.       \tag{6}
\]

If a transverse slice of `D^aff` has length at most `L(A,A')`, then (5) is an
isomorphism.  Indeed, `D^=` is already a closed subscheme of that length, so
the induced surjection of local Artin rings has equal length.

This reduces the global affine problem to a jet or initial-ideal upper bound;
no further lower-bound construction is required.

When the finite-flat conjecture

\[
\operatorname{rank}\mathfrak Z_k=2^k             \tag{7}
\]

holds, the predicted transverse length is

\[
\boxed{L(A,A')=2^{\sum_\rho|k_\rho|}.}           \tag{8}
\]

## Proven cases

1. For transfers `(1,-1)` over the degree-twelve `(6,6)` collision, the two
   blocks are dual numbers.  The exact second-jet upper bound is four, so the
   equalizer criterion proves `D^aff=D^=`.
2. For an isolated `k=2` block, direct elimination proves that allowing
   `U^2-V^3` to be affine gives exactly the same ideal as strong equality.
3. For the first global higher-transfer vector `(2,-2)`, the quadratic-cubic
   initial ideal has colength sixteen, equal to the strong product.  Hence the
   affine and strong correspondences coincide in this case as well.

## Conductor once affine rigidity holds

Suppose the affine equalizer criterion is satisfied.  Let `B_A,B_A'` be the
two regular completed normalization branches and let

\[
B_A\twoheadrightarrow D,
\qquad
B_{A'}\twoheadrightarrow D
\]

be their quotient maps to the product algebra (3), after incidence
elimination.  Then the completed component ring is

\[
\widehat{\mathcal O}_{C,\bar\xi}
=B_A\mathop\times_D B_{A'},                      \tag{9}
\]

and its conductor in `B_A\oplus B_A'` is the direct sum of the two
quotient kernels.  Thus classifying the local blocks and proving the one
global upper bound determine both the completed intersection and conductor.

## The global `(2,-2)` theorem

Take transfer vector

\[
(k_1,k_2)=(2,-2),
\]

over a collision with two multiplicity-twelve roots.  Formula (3) gives the
strong transverse algebra

\[
\mathfrak Z_2\widehat\otimes\mathfrak Z_2
\]

of length sixteen.

Use the twenty nonleading coefficients of the four monic factors as ambient
coordinates.  Equality of normalized coefficients supplies the coefficients
of degrees `2,...,23`.  Their linearization has rank twelve, hence tangent
dimension eight.  Four tangent directions are the reduced `S_0,S_1`
factor-space directions; slice them away and use transverse coordinates
`e_0,...,e_3`.

Twelve equations and twelve normal coefficient variables have an invertible
linear block.  Formal implicit elimination through degree three leaves the
homogeneous local standard basis

\[
\begin{aligned}
e_0^2-9e_1^2,&\qquad e_0e_1+3e_1^2,\qquad e_1^3,\\
e_2^2,&\qquad e_2e_3,\qquad e_3^3.               \tag{10}
\end{aligned}
\]

The first line becomes `(X^3,XY,Y^2)` after taking
`X=e_0`, `Y=e_0+3e_1`; the second does so with `X=e_3`, `Y=e_2`.
Consequently (10) is the tensor product of two coincident `Z_2` fibers.  Its
standard basis has sixteen elements and Hilbert function `(1,4,6,4,1)`.

The initial ideal of the full affine correspondence contains (10), so its
transverse length is at most sixteen.  The strong Hensel product is already a
closed length-sixteen subscheme.  The equalizer criterion gives

\[
\boxed{
D^{\mathrm{aff}}_{(2,-2)}=D^=_{(2,-2)}
=\mathfrak Z_2\widehat\otimes\mathfrak Z_2.}     \tag{11}
\]

This is the first genuinely global higher-transfer branch-intersection
theorem.  The next new cases are mixed transfer vectors such as `(2,-1,-1)`
and the first unknown local block `Z_3`.
