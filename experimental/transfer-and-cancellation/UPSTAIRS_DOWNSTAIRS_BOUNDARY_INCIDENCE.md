# Upstairs--downstairs boundary incidence

The target boundary-incidence diagram is only a quotient of the canonical
finite normalization.  This note defines a bipartite enhancement which keeps
the boundary primes upstairs, their finite maps to target divisors, all
multiple intersections on both sides, and the maps between those strata.

Use the admissible boundary setup

\[
 X\xhookrightarrow{j_F}\overline X_F
 \xrightarrow{\pi_F}Y                                      \tag{1}
\]

from
[the canonical boundary theorem](CANONICAL_BOUNDARY_INTERSECTION_INVARIANT.md).
Thus `X` and `Y` are normal integral affine varieties, `F:X -> Y` is dominant
and quasi-finite, `bar(X)_F` is the finite normalization of `Y` in `k(X)`,
and `j_F` is the canonical open immersion.  Write `partial_F` for its reduced
closed complement.

The maximal invariant is already the isomorphism class of the entire triple
(1).  The construction below is a finite boundary-component skeleton which
is much more computable while retaining substantially more information than
the target diagram alone.

There is no hidden contraction assumption in this terminology.  By
Proposition 1.2 of the canonical boundary theorem, `partial_F` is pure of
codimension one and every one of its irreducible components maps onto a prime
divisor of `Y`.  The first assertion uses that the distinguished open `X` and
the normal ambient scheme `bar(X)_F` are affine.  The second uses finiteness:
if a height-one prime `q` of the integral closure contracts to `p` in the
target ring, then `p` is nonzero and going down gives
`ht(p)<=ht(q)=1`, so `ht(p)=1`.  Thus a boundary divisor cannot be contracted
to a higher-codimension target stratum.  This would fail for a nonfinite
compactification map, such as a blow-up.

## 1. The bipartite vertices

Let

\[
 \mathcal E_F=\{E:E\text{ is a prime divisor of }\overline X_F,
                         \ E\subset\partial_F\}                \tag{2}
\]

and let `mathcal D_F` be the set of their images in `Y`.  The preceding
purity--noncontraction lemma says that these images are automatically prime
divisors.  The finite morphism gives a canonical surjection

\[
 q_F:\mathcal E_F\longrightarrow\mathcal D_F,
 \qquad E\longmapsto Z_E=\pi_F(E).                            \tag{3}
\]

The edge from `E` to `Z_E` records not merely `(e,f)`, but the complete finite
dominant morphism

\[
 \pi_E:E\longrightarrow Z_E,                                 \tag{4}
\]

including the induced extension
`kappa(Z_E) subset kappa(E)`.  The numerical labels

\[
 e(E/Z_E),\qquad f(E/Z_E)=[\kappa(E):\kappa(Z_E)]             \tag{5}
\]

are recoverable from (4), but (4) also retains branching and specialization
inside the residue-field cover.

Several vertices of `mathcal E_F` may map to the same target vertex.  They
remain distinct even when their numerical `(e,f)` labels agree.

The vertex set is therefore exhaustive at the level of irreducible boundary
components, rather than a choice to discard components with small images.
Higher-codimension boundary loci still occur inside these components.  The
multiple-intersection system below lists those forced by distinct boundary
components, while the complete edge morphisms `E -> Z_E` retain arbitrary
special loci inside a single component implicitly.  A further finite
stratification could be added only after specifying a canonical functorial
stratification rule; no arbitrary choice of strata is part of this invariant.

## 2. Upstairs and downstairs intersections

For every nonempty `S subset mathcal E_F`, define the upstairs
scheme-theoretic and reduced intersections

\[
 K_S(F)=\bigcap_{E\in S}^{\mathrm{sch}}E\subset\overline X_F,
 \qquad
 K_S^{\mathrm{red}}(F)=K_S(F)_{\mathrm{red}}.                 \tag{6}
\]

Empty intersections are allowed and retained as empty schemes.  If
`S subset T`, there is a canonical closed immersion `K_T -> K_S`.

For every nonempty `U subset mathcal D_F`, retain the downstairs intersections

\[
 J_U(F)=\bigcap_{Z\in U}^{\mathrm{sch}}Z\subset Y,
 \qquad
 I_U(F)=J_U(F)_{\mathrm{red}}.                               \tag{7}
\]

Put `q(S)={q(E):E in S}`.  Since `K_S` lies scheme-theoretically over every
`Z in q(S)`, restriction of `pi_F` gives a finite morphism

\[
 K_S(F)\longrightarrow J_{q(S)}(F).                          \tag{8}
\]

Let `W_S(F)` be its scheme-theoretic image.  Then (8) factors canonically as

\[
 K_S(F)\longrightarrow W_S(F)lhookrightarrow J_{q(S)}(F).   \tag{9}
\]

The first arrow is finite and schematically dominant.  If `S subset T`, then
`K_T -> K_S` induces a closed immersion `W_T -> W_S`.  Thus (6)--(9) form a
finite system of compatible upstairs strata, scheme images, downstairs
strata, and incidence arrows.

Recording `W_S` is essential.  The set-theoretic image of an upstairs
intersection forgets nilpotents and contact multiplicities, while the ambient
downstairs intersection `J_(q(S))` may contain components which are not met by
`K_S` at all.

## 3. Infinitesimal ramification data

Assume now that `k(X)/k(Y)` is generically separable, as it is for
characteristic-zero Keller maps.  The finite morphism has the canonical
relative different ideal

\[
 \mathfrak D_F=operatorname{Fitt}^0_{\mathcal O_{\overline X_F}}
                 \Omega_{\overline X_F/Y}.                   \tag{10}
\]

Record (10), its order along every `E in mathcal E_F`, and its restrictions to
the strata `K_S`.  For a Keller map, `pi_F` restricts to the étale map `F` on
`X`, so the different is supported in the canonical boundary.

One may strengthen (10) further by retaining the formal completion of the
finite map along each pair `K_S -> W_S`.  These formal maps encode transverse
contact and higher infinitesimal specialization which neither `(e,f)` nor the
scheme `K_S` alone determines.

There is generally no comparable canonical conductor between
`O_(bar(X)_F)` and `O_X`: the inclusion is an open immersion, not a finite
birational extension, and the affine ring inclusion may be a localization
with zero conductor.  Conductors attached to a chosen resolvent order inside
its normalization depend on that presentation.  The canonical replacements
are the open pair (1), the different (10), and the formal boundary maps.

## 4. Local monodromy as derived data

Let `V subset Y` be the maximal open over which `pi_F` is étale.  The finite
étale cover `pi_F^(-1)(V) -> V` is part of (1).  At a geometric point of a
stratum `W_S`, its restriction to a punctured strict-henselian neighborhood
gives a permutation representation of the corresponding local fundamental
group on the geometric generic fiber.  Its conjugacy class is therefore
derived functorially from the finite map and the chosen stratum; it need not
be added as coordinate-dependent branch cycles.

This local representation can distinguish covers with identical divisorial
`e,f` data.  A primitive resolvent supplies a valid computation only after its
local branches have been identified with the canonical normalization (1).

## 5. The enhanced invariant

### Definition 5.1

The **upstairs--downstairs boundary-incidence object**
`mathfrak I_(X->Y)(F)` consists of:

1. the bipartite finite sets `mathcal E_F -> mathcal D_F`;
2. every finite edge morphism `E -> Z_E`, including its function-field map;
3. all upstairs intersections `K_S` and their reductions;
4. all downstairs intersections `J_U` and their reductions;
5. every scheme image and finite factorization `K_S -> W_S -> J_(q(S))`;
6. all incidence arrows induced by inclusions of indexing subsets; and
7. in the separable case, the different ideal and, when desired, the formal
   completions along `K_S -> W_S`.

There is a canonical forgetful functor

\[
 \mathfrak I_{X\to Y}(F)\longrightarrow\mathfrak I_Y(F)       \tag{11}
\]

to the former target-only diagram.  It forgets the distinct upstairs
vertices, their intersections and scheme images, the actual residue covers,
and all infinitesimal ramification data.

## 6. Functoriality and stabilization

### Theorem 6.1 (stable upstairs--downstairs invariance)

Let `F:X -> Y` and `G:X' -> Y'` be admissible boundary setups.

1. A left--right isomorphism `G=beta F alpha` induces an isomorphism
   \[
   \mathfrak I_{X\to Y}(F)\simeq
   \mathfrak I_{X'\to Y'}(G).                                \tag{12}
   \]
2. For every `a>=0`, every scheme and morphism in
   `mathfrak I_(X->Y)(F times id_(A^a))` is the product with `A^a` of the
   corresponding algebraic scheme or morphism for `F`.  The optional formal
   completions are the completions of these products along the extended
   ideals.
3. Consequently the affine-cylinder class of the enhanced object is a stable
   polynomial left--right invariant.

**Proof.**  Source and target isomorphisms transport the function-field
inclusion, its integral closure, the distinguished affine-source open, and
the finite map.  Hence they transport the boundary primes, edge morphisms,
and all prime ideals defining (6)--(7).  They commute with sums of ideals and
radicals, so all full and reduced intersections are preserved.

The scheme-theoretic image of a finite morphism is defined by the kernel of
the map on coordinate rings.  Isomorphisms preserve that kernel, proving
functoriality of `W_S` and (9).  Kähler differentials and Fitting ideals are
also functorial under isomorphism, and formal completion is functorial for the
transported ideals.  This proves (12).

Under stabilization, integral closure commutes with polynomial extension and
the canonical triple becomes

\[
 X\times\mathbb A^a\hookrightarrow
 \overline X_F\times\mathbb A^a\longrightarrow
 Y\times\mathbb A^a.                                       \tag{13}
\]

Every boundary prime and every intersection is its old counterpart times
`A^a`.  Scheme-theoretic image commutes with this flat base change because
polynomial extension preserves the kernel defining the image.  Moreover,

\[
 \Omega_{(\overline X_F\times\mathbb A^a)/(Y\times\mathbb A^a)}
 \simeq\operatorname{pr}_1^*\Omega_{\overline X_F/Y},         \tag{14}
\]

so the different ideal and its restrictions pull back.  Formal completions
along extended ideals are functorial as well.  Affinely, if `I subset R`
defines a stratum, the stabilized completed ring is

\[
 \widehat{R[t_1,\ldots,t_a]}_{,I R[t_1,\ldots,t_a]},         \tag{15}
\]

not in general the ordinary polynomial ring over `hat(R)_I`; this completed
cylinder is the datum intended here.  This proves the stabilization and
stable-invariance assertions.  QED

## 7. First layer for C04 and C24

Theorem 5.1 of the
[boundary-exhaustion certificate](BOUNDARY_EXHAUSTION_CERTIFICATE.md)
determines the complete bipartite vertex set and the generic edge covers in
the canonical finite normalization.

For a generic degree-`n` weighted C04 seed over an algebraic closure:

- one boundary prime `E_Delta` maps birationally to the discriminant with
  ramification index two; and
- `n-3` boundary primes `E_rho` map birationally and unramifiedly to `C=0`.

For `C24_(m,r)` over an algebraic closure:

- one boundary prime `E_Delta` maps birationally to the discriminant with
  ramification index `r+1`; and
- `mr-1` boundary primes map birationally and unramifiedly to `P=0`.

Over the arithmetic coefficient field, the latter vertices are grouped by
the irreducible factors of the projective boundary polynomial after removal
of the distinguished affine root; the finite edge maps retain the resulting
residue extensions rather than only their degrees.

Thus the vertex and edge layer is already explicit.  The genuinely new C24
problem begins with the schemes

\[
 E_\Delta\cap E_i,
 \qquad E_i\cap E_j,
\]

their scheme images inside the thick target intersection, and the restrictions
of the different to those strata.

## 8. Why the enhancement is genuinely stronger

Two target-only diagrams may agree while the enhanced objects differ in any
of the following ways:

- the edge covers `E -> Z` may be nonisomorphic despite equal degree;
- equal-profile boundary primes over one `Z` may intersect differently over
  lower target strata;
- an upstairs intersection may map onto only one component or one thickening
  inside the larger downstairs intersection;
- the different or completed local maps may have different transverse
  structure; or
- the induced local monodromy representations may be nonconjugate.

For C24, the immediate next computation is the upstairs decomposition of the
thick target intersection

\[
 Q^{mr(m+1)}\{(r+1)RQ^m-C\}=0.                               \tag{16}
\]

The `mr-1` boundary primes over `P=0`, the discriminant boundary prime, and
their intersections upstairs should explain how the nilpotent affine-line
component in (15) is assembled.  This enhanced incidence is also a natural
candidate for distinguishing C24 maps attached to different parameter roots,
where the target-only numerical signatures agree.
