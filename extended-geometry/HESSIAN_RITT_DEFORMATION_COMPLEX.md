# The Hessian--Ritt deformation complex

Work over a field \(k\) of characteristic zero.  This note replaces the
degree-by-degree Gröbner-growth formulation of `OP-RITT` by a deformation
problem attached to the coefficient-decorated Ritt 2-complex.

The replacement has three layers.

1. Ritt--Engström--Ziegler theory determines the reduced relation-graph
   components.
2. A finite tree/cellular cotangent complex determines tangent, excess, and
   obstruction groups along each such component.
3. The completed derived deformation algebra, not its tangent cohomology
   alone, determines the nilpotent thickening.

This separation is essential.  The degree-thirty braid has only one excess
first-order direction along the monomial divisor, but its three labelled
sectors have nilpotence indices \(4,3,4\) and different transverse Artin
algebras.  No vector-space tangent complex can recover those numbers without
its differential, filtration, and higher operations.

## 1. The composition differential

Let

\[
 \mathcal P_d=\{x^d+a_{d-1}x^{d-1}+\cdots+a_1x\}
\]

be the affine space of monic original degree-\(d\) polynomials.  For a
complete decomposition word
\(\mathbf d=(d_1,\ldots,d_r)\), write

\[
 \mu_{\mathbf d}:\prod_{i=1}^r\mathcal P_{d_i}\longrightarrow
 \mathcal P_N,\qquad
 (f_1,\ldots,f_r)\longmapsto f_1\circ\cdots\circ f_r.
                                                               \tag{1.1}
\]

Put

\[
 A_i=f_1\circ\cdots\circ f_{i-1},\qquad
 B_i=f_{i+1}\circ\cdots\circ f_r.
\]

The derivative of (1.1) is the sparse substitution operator

\[
 d\mu_{\mathbf d}(\dot f_1,\ldots,\dot f_r)
 =
 \sum_{i=1}^r
 \bigl(A_i'\circ f_i\circ B_i\bigr)
 \bigl(\dot f_i\circ B_i\bigr).                              \tag{1.2}
\]

Thus every linear calculation needed by the deformation theory is assembled
from differentiation, substitution, and addition along the decomposition
tree.  It does not require elimination in the \(N-1\) ambient coefficients.

If factors are not normalized, an internal edge carries an infinitesimal
affine coordinate change \(\xi=(ax+b)\partial_x\).  The intrinsic tangent
complex of a factorization is

\[
 C^{-1}_{\mathbf d}
   =\bigoplus_{\text{internal edges}}\mathfrak{aff}_1
 \longrightarrow
 C^0_{\mathbf d}
   =\bigoplus_i T_{f_i}\mathcal P_{d_i}
 \xrightarrow{\,d\mu_{\mathbf d}\,}
 T_f\mathcal P_N.                                            \tag{1.3}
\]

The first arrow differentiates postcomposition on the inner factor and
inverse precomposition on the outer factor.  The monic-original charts used
in this repository are slices for this gauge action, so their working
complex is obtained from (1.3) by deleting the contractible gauge summand.

For Hessian composition, replace the last term by
\(T_{\pi_N(f)}\mathcal K_N\) and the last arrow by
\(d\pi_N\circ d\mu_{\mathbf d}\), where \(\pi_N\) forgets the linear
coefficient.

## 2. Several cuts as a derived intersection

Let \(D=\{\mathbf d_0,\ldots,\mathbf d_s\}\) be the requested decomposition
words and put \(X_i=\prod_j\mathcal P_{d_{ij}}\).  Their ordinary polynomial
intersection is the classical truncation of

\[
 X_D^{\mathrm{der}}
 =
 X_0\times^{\mathbf R}_{\mathcal P_N}
 X_1\times^{\mathbf R}_{\mathcal P_N}\cdots
 \times^{\mathbf R}_{\mathcal P_N}X_s.                       \tag{2.1}
\]

Here the normalized triangular reconstruction identifies each factor chart
with its composition incidence subscheme.  This closed-immersion property is
what identifies the classical fiber product in (2.1) with the
scheme-theoretic intersection; injectivity on field-valued points alone
would not suffice.

At a synchronized factorization its normalized tangent complex is the
two-term complex

\[
 \mathbb T_D:\quad
 \bigoplus_{i=0}^s T X_i
 \longrightarrow
 \bigoplus_{i=1}^s T\mathcal P_N,\qquad
 (u_i)\longmapsto
 \bigl(d\mu_i(u_i)-d\mu_0(u_0)\bigr)_{i=1}^s.                \tag{2.2}
\]

Here \(H^0(\mathbb T_D)\) is the common tangent space and
\(H^1(\mathbb T_D)\) is the first excess/obstruction space.  Formula (2.2)
already explains why adding one composite cut can kill the unique excess
direction in the degree-thirty braid.

For the Hessian intersection use \(\mathcal K_N\) in (2.1)--(2.2).  The
canonical reconstructed linear coefficients
\(\lambda_i\) define the defect morphism

\[
 \Delta\lambda:X_{D,H}^{\rm der}\longrightarrow
 V_D:=\mathbb A^{s+1}/\mathbb A^1_{\rm diagonal},\qquad
 \Delta\lambda=(\lambda_i-\lambda_0)_{i=1}^s.                \tag{2.3}
\]

Reduced synchronization says that (2.3) vanishes after passing to the
reduction.  Scheme-theoretic synchronization says that its coordinate
functions vanish in \(H^0(\mathcal O_{X_{D,H}^{\rm der}})\), equivalently
that the classical defect map factors through the origin.  A compatible
null-homotopy in the completed derived algebra is the stronger datum needed
to transport this vanishing coherently through moves and 2-cells.  This is
the derived replacement for checking
\(\lambda_i-\lambda_j\in H_D\) by a new Gröbner basis in every degree.

## 3. The cellular complex of Ritt moves

A set of words is not merely a collection over \(\mathcal P_N\).  Adjacent
Ritt moves give coefficient correspondences, commuting moves give squares,
and the Coxeter relation gives braid hexagons.  Let \(K_D\) be the resulting
coefficient-decorated Ritt 2-complex.

Fix a reduced power or Dickson component \(B\) supplied by the
relation-graph theorem and complete every chart and correspondence along its
map from \(B\).  Applying the relative tangent complex cellwise gives a
double complex

\[
 C^{p,q}_{K_D/B}
 =
 \bigoplus_{\sigma\in K_D^{(p)}}
 \mathbb T^q_{\widehat X_\sigma/B},
 \qquad 0\le p\le2,                                          \tag{3.1}
\]

whose horizontal differential is the signed restriction map and whose
vertical differential is assembled from (1.2).  Its totalization

\[
 \mathfrak g^{\mathrm{lin}}_{K_D/B}
 =\operatorname{Tot} C^{\bullet,\bullet}_{K_D/B}             \tag{3.2}
\]

is the small linear deformation complex attached to the decomposition
diagram.

Its interpretation is:

* \(H^0\): infinitesimal deformations common to all requested
  decompositions;
* \(H^1\): path mismatch, excess conormal directions, and infinitesimal
  synchronization defects;
* \(H^2\): failure to fill commuting or braid cells coherently.

The complex is small in combinatorial width: it has one summand per factor,
move, and 2-cell, and its maps are the substitution operators (1.2).
Its polynomial modules still remember the factor degrees; no claim of a
degree-independent finite-dimensional bound is needed.

The 2-cell labels cannot be discarded.  In the degree-thirty braid, sectors
omitting cuts \(10,15,6\) have different filtered local algebras although
they lie on the same unlabelled hexagon.  Accordingly (3.1) is a complex of
coefficient correspondences, not the ordinary cellular cochain complex of a
permutohedron.

## 4. Why the cotangent complex must be completed

The linear complex (3.2) controls square-zero extensions only.  The full
formal neighborhood is controlled by the relative cotangent complex

\[
 i^*L_{\widehat X_D/B}
 \quad\text{or, dually, by}\quad
 \mathfrak g_{K_D/B}
 =
 \mathbf R\!\operatorname{Hom}
 (i^*L_{\widehat X_D/B},\mathcal O_B)[-1],                   \tag{4.1}
\]

where \(i:B\hookrightarrow\widehat X_D\).  Under the usual perfectness and
formal-moduli hypotheses this complex carries dg-Lie or \(L_\infty\)
brackets.  Its completed Chevalley--Eilenberg algebra,

\[
 \widehat{\mathrm{CE}}(\mathfrak g_{K_D/B})
 \simeq
 \widehat{\operatorname{Sym}}_{\mathcal O_B}
 \bigl(\mathfrak g^\vee_{K_D/B}[-1]\bigr),                   \tag{4.2}
\]

with the full differential, is the object whose classical \(H^0\) recovers
the nilpotent structure.  By perfect duality the generators in (4.2)
identify with \(i^*L_{\widehat X_D/B}\).  The displayed symmetric algebra is
derived and completed; its differential contains the brackets, so this is
not a claim that the thickening is split.

Filter (4.2) by its augmentation ideal \(F\).  Its associated graded algebra
is

\[
 \operatorname{gr}^m_F\widehat{\mathrm{CE}}(\mathfrak g)
 \cong
 \operatorname{Sym}^m_{\mathcal O_B}
 \bigl(\mathfrak g^\vee[-1]\bigr).                           \tag{4.3}
\]

The linear part of the differential is dual to (3.2), while higher brackets
give the higher-filtration differentials.  The obstruction to extending a
comparison from order \(m\) to \(m+1\) is a cohomology class in the
corresponding twisted associated-graded deformation complex; its precise
coefficient module is determined by (4.3), not assumed in advance to be a
split symmetric power of an ordinary normal bundle.  The degree-thirty
nilpotence and annihilator tables are finite samples of this tower.  The
degree-forty-two spectator calculation shows why the tower should be
relative to the full boundary: raw path ideals contain a common
spectator-dependent layer, while the path-to-boundary tangent difference
remains one-dimensional.

## 5. The all-degree theorem target

The proposed replacement for the old elimination programme is the following.

> **Deformation-complex Hessian--Ritt theorem (target).**  For every
> characteristic-zero multiple Hessian-composition intersection:
>
> 1. its reduced irreducible components are the power/Dickson
>    relation-graph components given by tame Ritt--Engström theory;
> 2. after completion along any component \(B\), its derived formal
>    neighborhood is the homotopy limit of the coefficient-decorated factor,
>    move, and 2-cell diagram, and its cotangent complex is computed by the
>    totalization (3.1);
> 3. the Hessian linear-lift defect (2.3) is null-homotopic in that completed
>    complex;
> 4. the obstruction classes for commuting squares and braid hexagons vanish
>    in \(H^2\) of every twisted filtered layer arising from (4.3);
> 5. \(H^0\) of (4.2) gives the scheme-theoretic intersection, including its
>    nilpotent synchronization defects and their gluing along power--Dickson
>    boundary divisors.

Items 2--4 are the new content.  Item 1 is reduced collision theory, and
item 5 is the formal consequence once the completed derived model is proved.
This formulation permits nontrivial nilpotent path schemes: coherence means
a specified derived comparison around a 2-cell, not literal equality of its
two underived half-braids.

## 6. Proof architecture

An all-degree proof can now be divided into structural lemmas.

1. **Tree model.**  Prove that the normalized factor chart presents the
   cotangent complex by (1.3), functorially under contraction and refinement
   of decomposition trees.
2. **Move model.**  Compute the relative cotangent complex of the universal
   power and Dickson Ritt correspondences.  Show that arbitrary degrees are
   obtained by base change and outer/inner composition.
3. **Cell descent.**  Prove that the derived intersection is the homotopy
   limit of the coefficient-decorated Ritt 2-complex, so that its cotangent
   complex is the totalization (3.1).
4. **Synchronization null-homotopy.**  Upgrade the missing-line theorem and
   common-right-factor top-jet theorem from radical vanishing to a
   null-homotopy of (2.3) on every filtered layer.
5. **Universal 2-cells.**  Establish obstruction vanishing for the universal
   commuting square and the labelled three-factor braid.  Composition and
   base change should then supply all spectators, provided the comparison is
   made relative to the full cell boundary.
6. **Boundary gluing.**  Compute the completed power--Dickson overlap once,
   including the \(z\)-adic filtration, and descend the sector modules along
   adjacent cells.

The reusable local calculations are concentrated in steps 2, 4, and 5.
They are universal factor-degree calculations, not a census over total
degrees \(N\); steps 1 and 3 are the structural comparison theorems that
make this reduction valid.

## 7. Immediate tests

Before claiming the theorem, the following tests are decisive.

* Reconstruct the three degree-thirty transverse Artin algebras from
  \(H^0\widehat{\mathrm{CE}}(\mathfrak g)\), including generator orders
  \((5)\), \((2,2)\), and \((4,2)\).
* The completed degree-forty-two ideal flag and first conormal layers are now
  computed: the prime-omitting path equals the boundary, the spectator
  quotient has minimal \(z\)-annihilation exponent one, and the relative
  sector quotient has exponent eight.  Compute their module presentations
  and the extension class in the cotangent transitivity triangle.
* Verify a commuting square with both power and Dickson labels; a bare
  topological square has no scheme-theoretic content.
* Show that the lift cocycle is null-homotopic, rather than merely zero on
  tangent spaces, on the existing degree-thirty and degree-forty-two charts.

Passing these tests would justify replacing Gröbner growth by a finite list
of universal cotangent-complex calculations.  Until then the deformation
complex is a precise theorem programme, not an established all-degree
result.

The first linear regression is now implemented:

```bash
.venv/bin/python scripts/verify_hessian_ritt_deformation_complex.py
```

It proves the tree formula (1.2) by literal dual-number differentiation,
computes the unlabelled filled-braid cellular cohomology as a baseline, and
reconstructs the three point-cotangent homology pairs from the exact
transverse complete-intersection presentations.  It also verifies that the
cut-\(15\) and cut-\(6\) sectors have identical linear cotangent ranks but
different completed Hilbert data.  This is a regression for the proposed
linearization and for the necessity of completion, not yet the
coefficient-decorated homotopy-limit theorem.

## Relation to the existing notes

The [general Hessian--Ritt note](GENERAL_HESSIAN_RITT_INTERSECTIONS.md)
supplies the reduced relation-graph theorem and the canonical lift cocycle.
The [restricted synchronization theorem](RESTRICTED_HESSIAN_SYNCHRONIZATION_THEOREM.md)
proves radical synchronization and records the primary frontier.  The
[Ritt move 2-complex calculation](RITT_MOVE_2_COMPLEX.md) supplies the
degree-thirty cotangent/Tor data and the degree-forty-two warning that the
correct comparison is relative to the full cell boundary.
