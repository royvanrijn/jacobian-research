# The Hessian--Ritt deformation complex

Let the ground field \(k\) carry the characteristic label
\(\chi(k)\in\{0,p\}\).  The all-degree theorem target in this note remains a
characteristic-zero statement, but the deformation complex itself is
defined fiberwise with this label.  In characteristic \(p\), its ordinary
composition part must be augmented by the Frobenius cell module introduced
below.  This note replaces the degree-by-degree Gröbner-growth formulation
of `OP-RITT` by a deformation problem attached to the
coefficient-decorated Ritt 2-complex.

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

### 1.1 Characteristic labels on composition edges

Every oriented composition edge \(H\circ R\), with \(R\) monic original of
degree \(r>1\), now carries the label

\[
 \bigl(\chi(k),\,\deg H',\,\tau_p(H)\bigr),\qquad
 \tau_p(H)\in\{\mathrm{Frob},\mathrm{aff\mbox{-}Frob},
                         \mathrm{ordinary}\}.               \tag{1.4}
\]

In characteristic zero the last two entries are just
\((\deg H',\mathrm{ordinary})\).  In characteristic \(p\), put

\[
\begin{array}{c|c|c}
\tau_p(H)&\text{outer polynomial}&H'\\ \hline
\mathrm{Frob}&H=G(x^p)&0\\
\mathrm{aff\mbox{-}Frob}&H=ax+G(x^p),\ a\ne0&a\\
\mathrm{ordinary}&\text{all remaining cases}&\deg H'\ge1.
\end{array}                                                  \tag{1.5}
\]

Thus a separable outer polynomial can still lie in the second row.  It is
not a purely inseparable degeneration, but Hessian projection cannot see
one normalized inner tangent there.

For an integral outer polynomial \(H=\sum h_jx^j\), reduction at \(p\) is
flagged by computing \(\overline{H}'=\sum\overline{j h_j}x^{j-1}\).
The Frobenius flag is raised exactly when every \(j h_j\) vanishes modulo
\(p\); the affine--Frobenius flag is raised exactly when all terms with
\(j\ge2\) vanish and \(\overline{h_1}\ne0\).  For monic \(H\) of degree
\(m>1\), the first flag can occur only when \(p\mid m\), but
\(p\mid m\) is not a pointwise classification.

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

### 3.1 The split Keller control complex

The
[separated-product Keller--Ritt theorem](../verified/KELLER_RITT_PRODUCT_THEOREM.md)
supplies a strict zero-defect control for this cellular architecture.
For \(\prod_iF_{n_i}\), the coordinate factors commute as polynomial Keller
maps, the intermediate-field lattice is Boolean, and the canonical split
composition series realize the Coxeter squares and braid hexagons
literally.  Thus their comparison modules and all higher path obstructions
vanish before reduction or completion.

This does not solve the coefficient-coupled cells studied here.  It fixes
the normalization of the obstruction theory: a universal commuting or braid
class should restrict to zero on the separated-product locus, and any
nonzero class must measure coupling of factor coefficients, affine
reconstruction opens, or boundary data rather than the bare Coxeter
relation.

### 3.2 The Frobenius cell module

Let \(\pi_{>q}\) retain the coefficients of degrees greater than \(q\).  The
one-sided invisible module on an oriented edge \(H\circ R\) is

\[
 \mathfrak F^{(q)}_p(H,R)
 =
 \ker\left(
 xk[x]_{<r}\xrightarrow{\ U\mapsto\pi_{>q}(H'(R)U)\ }
 k[x]_{>q}\right).                                         \tag{3.3}
\]

The complete positive-characteristic classification proved in
[the counterexample-search note](../NEW_COUNTEREXAMPLE_SEARCHES.md#theorem-a1-complete-one-sided-tangent-classification)
gives, for the Hessian cutoff \(q=1\),

\[
\mathfrak F^{(1)}_p(H,R)=
\begin{cases}
xk[x]_{<r},&\tau_p(H)=\mathrm{Frob},\\
kx,&\tau_p(H)=\mathrm{aff\mbox{-}Frob},\\
0,&\tau_p(H)=\mathrm{ordinary}.
\end{cases}                                                  \tag{3.4}
\]

More generally, if \(H'\ne0\) and \(d=\deg H'\), it is spanned by
\(x,\ldots,x^s\), where
\[
 s=\max(0,\min(r-1,q-rd)).
\]
This is an exact theorem, not a bounded-search heuristic.

For each factor, move, and labelled 2-cell in \(K_D\), take the direct sum
of (3.3) over its oriented composition edges and use the same signed
restriction maps as in (3.1).  The result is the **Frobenius cell module**
\(\mathfrak F^{\bullet,(q)}_{K_D/B,p}\).  Naturality of the composition
differential makes it a cellular subcomplex and gives

\[
 0\longrightarrow
 \mathfrak F^{\bullet,(q)}_{K_D/B,p}
 \longrightarrow
 \mathfrak g^{\mathrm{lin},(q)}_{K_D/B,p}
 \longrightarrow
 \mathfrak g^{\mathrm{vis},(q)}_{K_D/B,p}
 \longrightarrow0.                                        \tag{3.5}
\]

The quotient is the visible complex.  A vector-space splitting
of (3.5) may be chosen cellwise, but it is not canonical and does not erase
the attaching maps.  In particular, Frobenius tangents are not to be folded
into an unexplained increase of \(H^0\) or \(H^1\): they are a separate
labelled cell module whose restrictions around commuting and braid cells
must be checked.

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

### 5.1 Positive-characteristic extension target

The corresponding characteristic-\(p\) programme is deliberately separate
from the theorem target above:

1. stratify every factor, move, and relation-graph cell by the labels
   (1.4)--(1.5);
2. use the ordinary cellular complex on the ordinary stratum and attach
   \(\mathfrak F_{K_D/B,p}^{\bullet,(1)}\) on the Frobenius and
   affine--Frobenius strata;
3. prove synchronization only after computing the cohomology of the
   augmented complex, rather than by reducing a characteristic-zero tangent
   rank;
4. reprove the relation-graph input on the required tame locus.  When
   \(p\) divides a participating factor degree, the characteristic-zero
   Ritt--Engström component classification is not being asserted.

Thus (3.5) modifies synchronization, Hessian projection, and the
relation-graph deformation complex at the same place.  The theorem (3.4)
classifies the new one-sided tangent summand; it does **not** by itself prove
a positive-characteristic all-degree Ritt theorem or classify wild
relation-graph components.

## 6. Proof architecture

An all-degree proof can now be divided into structural lemmas.

1. **Tree model.**  Prove that the normalized factor chart presents the
   cotangent complex by (1.3), functorially under contraction and refinement
   of decomposition trees.
2. **Move model.**  Compute the relative cotangent complex of the universal
   power and Dickson Ritt correspondences.  Show that arbitrary degrees are
   obtained by base change and outer/inner composition.
3. **Cell descent.**  Cotangent descent for the full bar presentation of the
   actual derived intersection is now formal.  Prove coefficient
   effectivity of its compression to the finite decorated Ritt complex, so
   that its perfect dual is the totalization (3.1).
4. **Synchronization null-homotopy.**  Upgrade the missing-line theorem and
   common-right-factor top-jet theorem from radical vanishing to a
   null-homotopy of (2.3) on every filtered layer.
5. **Universal 2-cells.**  The split Keller control cells have zero
   obstruction by KRP1.  Establish the relative obstruction class for the
   coefficient-coupled commuting square and labelled three-factor braid,
   prove its required vanishing, and show that it restricts to the split
   zero section.  Composition and base change should then supply all
   spectators, provided the comparison is made relative to the full cell
   boundary.
6. **Boundary gluing.**  Compute the completed power--Dickson overlap once,
   including the \(z\)-adic filtration, and descend the sector modules along
   adjacent cells.
7. **Characteristic audit.**  For each integral model, record the
   characteristic label of every outer edge, attach (3.3), and separate
   good-reduction certificates from calculations made only over
   \(\mathbb Q\).

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
  sector quotient has exponent eight.  The completed ideal-module
  extension and the first-Postnikov conormal projection are non-split; the
  latter proves that the cotangent transitivity connecting morphism is
  nonzero.  The sector-to-total conormal map is injective after completion,
  and its apparent two-dimensional kernel modulo the base square is
  base-change Tor.  Compute the individual higher cotangent homology
  modules.
* Insert the exact sixth-jet embedded-support class \(c_6\) into the
  extension-retaining cellular model.  First determine whether it maps to
  the cut-\(6\) non-splitting class.  Then compute its restrictions through
  the completed cut-\(14\) and cut-\(21\) splittings and evaluate braid
  coherence.  The known modular order-seven lift is evidence only; the
  characteristic-zero theorem currently ends at order six.
* Compare a commuting square with both power and Dickson labels against the
  strict KRP1 split square; a bare topological square has no
  scheme-theoretic content.
* Show that the lift cocycle is null-homotopic, rather than merely zero on
  tangent spaces, on the existing degree-thirty and degree-forty-two charts.
* Reduce integral test charts at primes meeting all three rows of (1.5) and
  verify that the extra kernel is exactly the Frobenius cell module (3.4).

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

The next associated-graded regression is the
[cellular cotangent prototype](HESSIAN_RITT_CELLULAR_COTANGENT_PROTOTYPE.md).
It attaches explicit modules and signed matrices to vertices, moves,
commuting squares, and braid cells.  Its degree-thirty sectors have
\((\dim H^0,\dim H^1,\dim H^2)=(2,1,0)\).  In degree forty-two it separates
the first conormal \(H^1\) into the exact \(z^8\) sector and \(z\) spectator
layers and replays their finite-jet dimensions.  This is an exact
associated-graded model.  Its next verifier proves that the
sector--spectator module extension splits at order two and is non-split at
orders three and four, with explicit change-of-splitting obstruction
functionals.  A presentation-first finite tensor quotient further proves
that the completed ideal-module extension is non-split.  Killing the full
normal ideal then gives a non-split \(4\to6\to2\) conormal quotient, proving
that the actual cotangent transitivity connecting morphism is nonzero.
The quadratic overlap vanishes after completion, so its first homology
sequence is short exact; the two-dimensional kernel introduced by
base-square reduction is ordinary Tor.  What remains is to compute the
individual higher homology modules and prove the homotopy-limit comparison.

The general local algebra behind this conclusion is now separated in
[cellular Postnikov transitivity](CELLULAR_POSTNIKOV_TRANSITIVITY.md).
For an arbitrary ideal flag it proves the overlap formula and constructs the
successive conormal extensions; its executable finite-module tower accepts
any number of layers.  This closes the formal flag-transitivity part of the
programme.  It does not close the Ritt-specific cellular-descent comparison.

The two rotated degree-forty-two first-conormal extensions are now computed
separately from the non-split cut-\(6\) tower above.  Compatible order-five,
order-six, and order-seven sections come from one nested tensor
presentation.  Exact completed presentations give polynomial sections

\[
 e_4+(-3(1+\tau)^2+2\zeta)e_6,\qquad
 e_4+(-4(1+\tau)^3+8(1+\tau)\zeta)e_7
\]

for cuts \(14\) and \(21\).  Hence both rotated completed extension and
inverse-limit torsor obstruction classes vanish.  The one-dimensional
fourth-jet quotient-ring discrepancy is non-flat base-change Tor, while the
different correction polynomials retain the labelled sector data.

The categorical part of that comparison is now proved in
[cotangent descent for the Hessian--Ritt diagram](HESSIAN_RITT_COTANGENT_DESCENT_COMPARISON.md):
the full bar diagram computes the actual cotangent complex, and a complete
two-skeleton determines \(H^0,H^1\).  Local coefficient effectivity is now
proved on the rotated conormal towers; the remaining Ritt-specific statement
is coherence under the actual factor-chart braid and commuting-cell
restriction maps.  Genuine \(H^2\) also requires the next Coxeter cells; the
permutohedron three-cell kills the topological \(H^2\) of the four-factor
two-skeleton exactly.

## Good-reduction ledger for existing synchronization certificates

A characteristic-zero synchronization result survives reduction only when
both its algebraic certificate and its tangent model survive.  The current
ledger is:

| Input | Reduction statement |
|---|---|
| Composition differential (1.2) and literal integral polynomial identities | Survive in every characteristic after coefficient reduction. |
| Common-right-factor top-jet theorem `HRCF` | Survives over every reduction in which the outer degree \(m\) is a unit, hence over fields with \(p\nmid m\).  It makes no claim when \(p\mid m\). |
| A denominator-cleared membership \(c\Delta=\sum A_jg_j\) over \(\mathbb Z\) | Survives at primes with \(p\nmid c\), provided every localized chart denominator, unit pivot, and exact-degree leading coefficient remains invertible. |
| Published `HRSYNC` and transported degree-\(30\)/\(42\) Gröbner certificates over \(\mathbb Q\) | No blanket positive-characteristic claim is currently registered: their complete denominator/pivot exceptional sets have not been published.  A matching basis size modulo one prime is a regression, not a good-reduction theorem. |
| Characteristic-zero tangent exactness or null-homotopy | Transfers only on fibers where the cleared certificate is valid **and** every relevant edge has \(\mathfrak F^{(1)}_p(H,R)=0\).  The affine--Frobenius row must be excluded even though \(H\) is separable. |
| Tame power/Dickson relation-graph classification | May be reused only after its prime-to-degree and separability hypotheses are checked on the chosen fiber; it is not an automatic reduction of the characteristic-zero component theorem. |

Accordingly, “proved over \(\mathbb Q\)” never means “proved for all but
unspecified primes” in this programme.  A good-reduction claim must list the
integer certificate multiplier, localized denominators and pivots, relevant
outer-edge labels, and the resulting explicit exceptional-prime set.

## Relation to the existing notes

The [general Hessian--Ritt note](GENERAL_HESSIAN_RITT_INTERSECTIONS.md)
supplies the reduced relation-graph theorem and the canonical lift cocycle.
The [restricted synchronization theorem](RESTRICTED_HESSIAN_SYNCHRONIZATION_THEOREM.md)
proves radical synchronization and records the primary frontier.  The
[positive-characteristic tangent theorem](../NEW_COUNTEREXAMPLE_SEARCHES.md#theorem-a1-complete-one-sided-tangent-classification)
supplies the Frobenius cell module (3.3)--(3.4).  The
[Ritt move 2-complex calculation](RITT_MOVE_2_COMPLEX.md) supplies the
degree-thirty cotangent/Tor data and the degree-forty-two warning that the
correct comparison is relative to the full cell boundary.
