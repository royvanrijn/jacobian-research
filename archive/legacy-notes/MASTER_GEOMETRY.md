# Master geometry: two universal constructions

The weighted theory in this repository is organized by two universal
constructions.  On the target side one projects the locus of polynomials with
no simple root modulo affine polynomials.  On the source side one pulls the
universal root cover back to an affine pencil, normalizes it, and retains the
open on which the marked root reconstructs a source point.

This note is a synthesis.  The target-side quotient statement is proved
precisely in [MASTER_QUOTIENT_THEOREM.md](../../experimental/transfer-and-cancellation/MASTER_QUOTIENT_THEOREM.md); only
its closed-immersion upgrade across all collision strata remains conjectural.
The ambient multiplication maps, their torus-normalized hyperplane slices,
the uniqueness of the cubic affine complement, and the obstruction to the
next `(2,3)` candidate are developed in
[UNIVERSAL_FACTORIZATION_GEOMETRY.md](../../experimental/transfer-and-cancellation/UNIVERSAL_FACTORIZATION_GEOMETRY.md).
Claim C01 remains independently certified by Dean Cureton's Lean 4
formalization; see [LEAN_C01.md](../../verified/LEAN_C01.md).

## 1. The target-side construction

Fix a characteristic-zero field `k` and put

\[
 V_n=k[W]_{\le n},\qquad L=\langle 1,W\rangle,
 \qquad \pi:V_n\longrightarrow V_n/L.
\]

A seed `H` determines the affine plane

\[
 H+L=\{H-sW+t:(s,t)\in\mathbb A^2\}.
\]

Thus the unnormalized inverse pencil depends only on `[H]`.  The conditions

\[
 H(0)=H'(0)=H(1)=0,\qquad H'(1)=-1
\]

choose a representative and a scale in the relevant open subset of
`V_n/L`.

Let `B_n^circ` be the exact-degree locus of polynomials all of whose roots
have multiplicity at least two.  Its projective closure is the union of the
classical coincident-root loci with partitions having no part equal to one.
Then:

- the ordinary discriminant meets `H+L` in the discriminant curve of C05;
- $B_n^\circ\cap(H+L)$ is the set of omitted inverse polynomials of C07; and
- Mason--Stothers rigidity says that a normalized admissible pencil contains
  at most one such polynomial.

The exact-degree qualification matters.  The projective closure contains
boundary forms and can meet the center of the linear projection, so the
strong projection statement below is only proposed on the relevant
quasi-affine open.

## 2. The normalization hyperplane

Write `M(W)=sum_(j=0)^n m_j W^j` and define the linear forms

\[
 \Phi(M)=M(1)-M(0)-M'(0)=\sum_{j=2}^n m_j,
 \qquad
 D(M)=M'(1)-M'(0)=\sum_{j=2}^n j m_j.
\]

On `Phi=0` and `D ne 0`, the scale-invariant normalized projection is

\[
 \rho(M)=\frac{-M+M'(0)W+M(0)}{D(M)}.             \tag{1}
\]

It satisfies

\[
 \rho(M)(0)=\rho(M)'(0)=\rho(M)(1)=0,
 \qquad \rho(M)'(1)=-1.
\]

Moreover the weighted forbidden divisor is exactly

\[
 M''(1)-2D(M)=0,                                  \tag{2}
\]

because `rho(M)''(1)=-M''(1)/D(M)`.  Consequently the exceptional seed locus
has the precise presentation

\[
 \mathcal N_n
 =\rho\!\left(B_n^\circ\cap V(\Phi)\cap
 D\bigl(D(M)(M''(1)-2D(M))\bigr)\right),          \tag{3}
\]

with the usual collision exclusions when one is naming an exact stratum.
Equivalently, quotient the set in parentheses by scalar multiplication, or
take its unique section `D=1`.  Formula (3) is the quotient/hyperplane form of
the full-contact classification already proved in C07--C10.

This explains the dimension calculation.  A coincident-root locus of type
`lambda` has affine monic dimension `ell(lambda)`.  The nonzero linear
condition `Phi=0` cuts its root normalization by one, giving

\[
 \dim \mathcal E_\lambda=\ell(\lambda)-1.
\]

The weighted-Vandermonde argument in
[UNIFORM_EXCEPTIONAL_SEEDS.md](../../experimental/geometry/UNIFORM_EXCEPTIONAL_SEEDS.md) additionally
proves that the map to seed coefficients does not lower this dimension.

## 3. Classical coincident-root loci

For a partition `lambda=(lambda_1,...,lambda_ell)` of `n`, let

\[
 \Delta_\lambda
 =\overline{\left\{\prod_i(W-r_i)^{\lambda_i}:r_i\ne r_j\right\}}.
\]

The standard root-parameter map, with equal parts divided by their symmetric
groups, is finite and birational onto `Delta_lambda`.  In quotient
coordinates it is expressed by monic polynomials `Q_mu` recording all roots
of multiplicity `mu`:

\[
 (Q_2,Q_3,\ldots)\longmapsto\prod_{\mu\ge2}Q_\mu^\mu. \tag{4}
\]

This is exactly the normalization map used in C08--C12.  In particular, for
a maximal partition it reduces to

\[
 (Q,R)\longmapsto Q^2R^3.                         \tag{5}
\]

The classical language accounts for the repository's combinatorics:

| Repository object | Coincident-root interpretation |
|---|---|
| contact partition | multiplicity partition `lambda` |
| collision order | coarsening order on partitions |
| common coarsenings | set-theoretic component intersections |
| equal-root quotient | standard symmetric normalization coordinates |
| allocation count | fiber of the normalization correspondence |

The existing proofs retain information not supplied merely by this
dictionary: the hyperplane section is nonempty and irreducible for maximal
`2/3` types, the admissible open is met, its quotient hypersurface is smooth,
and the local intersection multiplicities are computed.  Classical results
therefore give the ambient normalization framework; they do not by themselves
replace those slice and local-algebra arguments.

The complete finite-normalization, slice-nonemptiness, exact-closure, and
common-coarsening proof for C08--C11 is collected in
[COINCIDENT_ROOT_REBUILD.md](../../experimental/geometry/COINCIDENT_ROOT_REBUILD.md).

## 4. Why 2 and 3 are the cusp

The allowed multiplicity semigroup is

\[
 \{2,3,4,\ldots\}=\langle2,3\rangle.
\]

It is the value semigroup of the cusp

\[
 C=\operatorname{Spec}k[u,v]/(u^2-v^3),
 \qquad s\longmapsto(u,v)=(s^3,s^2).              \tag{6}
\]

Accordingly, `2` and `3` are both the maximal contact atoms and the two
generators of the cusp semigroup.  Formula (5) is the global monomial
factorization associated with these generators.  Comparing two allocations
and cancelling their common factors produces the polynomial cusp equation

\[
 \mathfrak Z_k=\{(U,V):U^2=V^3\},
 \qquad (U,V)=(S^3,S^2),\quad\deg S=k.             \tag{7}
\]

Thus `Z_k` is a polynomial mapping/factorization scheme through the cusp,
completed along maps which lift to its normalization.  Over a squarefree
`S`, its `k` separated roots carry `2^k` allocation choices.  The all-`k`
theorem identifies the complete block with the symmetric quotient of the
ordered-root square-zero Boolean thickening.  Collision changes the
multiplication without changing the length:

\[
 \operatorname{rank}\mathfrak Z_k=2^k
 \quad\text{for every }k\ge1.
\]

The special fibers have Hilbert series `(1+t)^k`, but need not
be Gorenstein.  They are flat degenerations of the squarefree subset algebra,
not literal tensor products of dual numbers.  The equality
`Z_k^aff=Z_k` says that permitting an affine difference between two
factorizations introduces no additional local structure for any block.  See
[ALL_K_TRANSFER_BLOCK_THEOREM.md](../../experimental/transfer-and-cancellation/ALL_K_TRANSFER_BLOCK_THEOREM.md).

## 5. The source-side construction

Fix an admissible seed `H`.  The universal inverse polynomial over target
space is

\[
 E_{A,B,C}(W)=H(W)-BCW+cAC^2.
\]

Its finite degree-`n` root incidence is

\[
 \mathcal I_H=V(E)\subset\mathbb A^3_{A,B,C}\times\mathbb A^1_W.
\]

This incidence is the pullback of the universal factorization map
`mu_(1,n-1)` along the three-parameter map
`(A,B,C) mapsto H-BCW+cAC^2`.  The coordinate `C` is the reconstruction scale,
so this is a three-dimensional parametrized pullback even though its image in
polynomial coefficient space lies in the affine pencil `H+<1,W>`.

Let `I_H^nu` be its normalization and let `R_H` be the locus where the
marked root reconstructs regular source coordinates.  C04 and C16 prove the
structural diagram

\[
 \mathbb A^3\ \cong\ \mathcal R_H
 \ \subset\ \mathcal I_H^\nu
 \ \longrightarrow\ \mathbb A^3_{A,B,C}.         \tag{8}
\]

The last arrow forgets the root and is the weighted Keller map under the
displayed isomorphism.  This gives a uniform reading of the main phenomena:

| Phenomenon | Universal-root explanation |
|---|---|
| noninjectivity | several surviving simple-root sheets |
| omission | the pencil meets `B_n^circ`, so no root is simple |
| discriminant | two or more root sheets meet |
| nonproperness | reconstruction has a polar divisor on `I_H^nu` |
| dicritical boundary | a removed polar divisor has nonconstant target image |

The cubic counterexample C01 is the first instance.  Its determinant and
three-point collision are also independently formalized in Lean 4 by Dean
Cureton; the later universal geometry is not asserted to be part of that
formalization.

## 6. The master quotient and global affine rigidity

On the monic exact-degree no-simple-root locus, the master quotient theorem
says that `pi` is a closed immersion onto its image. For a monic polynomial,
its constant and linear coefficients are therefore regular
scheme-theoretic functions of the coefficients of degrees `2,...,n-1`, even
where contact strata collide. After imposing `Phi=0`, `D!=0`, and weighted
admissibility, this identifies the resulting slice with the closed exceptional
seed locus. The exact-degree qualification is essential: the
degree-at-most-`n` locus contains all `(W-a)^2`, which have the same image
modulo affine polynomials.

The collision proof has two parts. C22 identifies every completed local
square/cube relation with its full Boolean thickening; Hensel factorization
tensors these relations over the distinct collision roots. The only remaining
coupling is one equation

\[
 Q^2R^3-S^2T^3=\lambda W+\mu.
\]

The Wronskian `S^2T^3(Q^2R^3)' - Q^2R^3(S^2T^3)'` is divisible by a monic
polynomial of degree greater than `n`, while the displayed affine equation
makes it have degree at most `n`. It vanishes, and its two leading
coefficients force `lambda=mu=0`. Thus affine and strong completed equalizers
coincide for every collision tree. Faithful flatness of completion upgrades
the local statement to the global closed immersion.

## 7. Reproducible coincident-root comparison

The official Macaulay2 `CoincidentRootLoci` package constructs
`Delta_lambda`, its ideal, parameterization, tangent space, and singular
locus.  The repository script
[`verify_coincident_root_slices.m2`](../../scripts/verify_coincident_root_slices.m2)
uses those package ideals through degree seven, substitutes the normalized
hyperplane section, saturates by the exact-degree and weighted-admissibility
factors, and eliminates the constant and linear coefficients.  In degree
eight it uses the package normalization maps and certifies maximal
differential rank on the admissible `Phi=0` slice; forcing a full ambient
implicitization there is needlessly memory-intensive.

Run it with a local Macaulay2 installation or the pinned Docker fallback:

```text
make verify-coincident-root-loci
```

It checks all maximal types in degrees five through eight:

\[
 (3,2),\quad(2,2,2),\quad(3,3),\quad(3,2,2),
 \quad(2,2,2,2),\quad(3,3,2).
\]

The asserted dimensions agree with C10. These are independent bounded-degree
ideal comparisons through degree seven and normalization-rank audits in degree
eight; they are regressions for, not proofs of, the all-degree
closed-immersion theorem.

## References and attribution

- Dean Cureton, [`deancureton/jacobian`](https://github.com/deancureton/jacobian),
  separately authored Lean 4 formalization of C01.
- Jaydeep Chipalkatti, [*On equations defining coincident root
  loci*](https://arxiv.org/abs/math/0110224).
- L. M. Feher, A. Nemethi, and R. Rimanyi,
  [*Coincident root loci of binary forms*](https://arxiv.org/abs/math/0311312).
- Simon Kurmann, [*Some remarks on equations defining coincident root loci*](https://doi.org/10.1016/j.jalgebra.2011.10.045).
- Macaulay2,
  [`CoincidentRootLoci` package documentation](https://macaulay2.com/doc/Macaulay2/share/doc/Macaulay2/CoincidentRootLoci/html/toc.html).
