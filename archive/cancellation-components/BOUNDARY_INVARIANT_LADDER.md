# The boundary-invariant ladder

The canonical finite normalization supports a hierarchy of stable invariants,
not a single boundary diagram.  This note makes the three layers explicit:

\[
 \boxed{qquad
 \mathfrak I^{\mathrm{formal}}(F)
 \longrightarrow
 \mathfrak I^{\mathrm{sch}}(F)
 \longrightarrow
 \mathfrak I^{\mathrm{red}}(F).
 \qquad}                                                   \tag{1}
\]

Each arrow is a canonical forgetful functor.  The reduced layer is intended
for inexpensive first obstructions, the scheme layer for multiplicity and
contact, and the formal layer for completed transverse geometry.

## 1. Common canonical setup

Use an admissible boundary setup

\[
 X\hookrightarrow\overline X_F\xrightarrow{\pi_F}Y,        \tag{2}
\]

with boundary primes `mathcal E_F`, target divisors `mathcal D_F`, and map
`q_F:mathcal E_F -> mathcal D_F`.  For nonempty indexing sets retain

\[
 K_S=\bigcap_{E\in S}^{\rm sch}E,qquad
 J_U=\bigcap_{Z\in U}^{\rm sch}Z,                           \tag{3}
\]

and the finite factorization through the scheme image

\[
 K_S\longrightarrow W_S\lhookrightarrow J_{q(S)}.         \tag{4}
\]

Purity and finiteness ensure that every irreducible boundary component occurs
in `mathcal E_F` and maps onto a member of `mathcal D_F`.

## 2. The reduced layer

### Definition 2.1

The object `mathfrak I^red_(X->Y)(F)` consists of:

1. the bipartite vertices `mathcal E_F -> mathcal D_F` and the generic
   labels `(e,f)` and boundary sheet contribution;
2. the reduced finite edge morphisms `E -> Z_E`;
3. every reduced upstairs intersection `K_S^red`;
4. every reduced downstairs intersection `I_U=(J_U)_red`;
5. every reduced scheme image `W_S^red` and the induced incidence maps.

This layer remembers which strata and components meet, but forgets tangency,
nilpotents, embedded structure, and transverse infinitesimal neighborhoods.
It is usually the cheapest layer to compute.

The former target-only reduced diagram is the quotient

\[
 \mathfrak I_Y^{\rm red}(F)
 =\mathfrak I_{X\to Y}^{\rm red}(F)
   /\{\text{upstairs vertices and maps}\}.                 \tag{5}
\]

For weighted marked-root theorem versus every noncubic cancellation member, this first target quotient already
separates the maps: the marked reduced target intersection is `A^1` for the
generic weighted seed and `A^1 disjoint-union G_m` for the cancellation construction.

## 3. The scheme-theoretic layer

### Definition 3.1

The object `mathfrak I^sch_(X->Y)(F)` retains the reduced layer and replaces
its intersections and images by the full schemes

\[
 K_S,qquad J_U,qquad W_S,                                 \tag{6}
\]

together with the finite schematically dominant maps (4), all closed
immersions induced by inclusions of indexing sets, and the ideals defining
these arrows.

It therefore retains intersection multiplicities, nilpotent thickenings,
embedded components, and scheme-theoretic contact which reduction erases.
Termwise reduction gives the second arrow in (1):

\[
 \mathfrak I^{\rm sch}_{X\to Y}(F)longrightarrow
 \mathfrak I^{\rm red}_{X\to Y}(F).                        \tag{7}
\]

The former target-only scheme diagram `mathfrak J(F)` is precisely the target
quotient `mathfrak I_Y^sch(F)`.  For cancellation construction its marked target intersection has
nilradical index `mr(m+1)`, while the generic weighted intersection is
reduced.  Thus the scheme layer strengthens, but is not needed for, the first
weighted--cancellation separation.

## 4. The formal layer

Assume the function-field extension is generically separable.

### Definition 4.1

The object `mathfrak I^formal_(X->Y)(F)` augments the scheme layer by the
following canonical data.

1. **Completed stratum maps.** For every `S`, retain the adic formal morphism
   
   \[
    \widehat{\overline X_F}_{K_S}longrightarrow
    \widehat Y_{W_S},                                      \tag{8}
   \]
   
   together with the compatible completed maps for all inclusions of
   strata.  One may equivalently record the corresponding systems of
   infinitesimal neighborhoods.

2. **Different ideals.** Retain
   
   \[
    \mathfrak D_F=operatorname{Fitt}^0
       \Omega_{\overline X_F/Y},                            \tag{9}
   \]
   
   its orders along every boundary prime, and its images in the completed
   local rings of (8).

3. **Finite-stratum conductors.**  Write `g_S:K_S -> W_S`.  Since `W_S` is
   the scheme-theoretic image, the map
   
   \[
    \mathcal O_{W_S}\hookrightarrow(g_S)_*\mathcal O_{K_S}
   \]
   
   is injective.  Retain the canonical conductor sheaf
   
   \[
    \mathfrak c_S=operatorname{Ann}_{\mathcal O_{W_S}}
    \left((g_S)_*\mathcal O_{K_S}/\mathcal O_{W_S}\right), \tag{10}
   \]
   
   together with its extension to `(g_S)_*O_(K_S)`.  When useful, also
   retain the conductor of `W_S^red` into its finite normalization.

4. **Boundary valuation filtrations.**  Every `E in mathcal E_F` gives the
   canonical divisorial valuation `v_E` of `L=k(X)`.  Retain
   
   \[
    \mathcal F_E^nL=\{f\in L:v_E(f)\ge n\},
    \qquad
    \mathcal F^{\mathbf n}L=
       \bigcap_E\mathcal F_E^{n_E}L,                       \tag{11}
   \]
   
   and their induced filtrations on the relevant local and completed
   coordinate rings.  Computations may retain a finite range of indices, but
   the invariant is the filtered object.

5. **Derived local monodromy.**  On punctured strict-henselian formal
   neighborhoods of the strata, retain the conjugacy class of the local
   permutation representation of the finite étale locus.  This is derived
   from (2) and (8), rather than from coordinate-dependent branch cuts.

Forgetting completions, ideals, and filtrations gives

\[
 \mathfrak I^{\rm formal}_{X\to Y}(F)longrightarrow
 \mathfrak I^{\rm sch}_{X\to Y}(F).                        \tag{12}
\]

### Conductor warning

Equation (10) is a conductor for a **finite schematically dominant stratum
map**.  There is generally no useful conductor for the open immersion
`X -> bar(X)_F`: it is not a finite ring extension, and a localization can
have zero conductor.  A conductor of a chosen nonnormal resolvent order is
presentation-dependent and is not part of the canonical ladder unless that
order itself has first been specified functorially.

## 5. Stable functoriality of the ladder

### Theorem 5.1

Each layer in (1) is a polynomial left--right invariant, and its affine-
cylinder class is a stable polynomial left--right invariant.  The two
forgetful arrows commute with left--right isomorphisms and stabilization.

**Proof.**  A left--right isomorphism transports the function-field
extension, its integral closure, the distinguished affine open, all boundary
prime ideals, and the finite map.  It therefore transports every object in
the reduced and scheme layers, commuting with sums of ideals, radicals, and
scheme-theoretic images.

Formal completion is functorial for transported ideals.  Kähler
differentials and Fitting ideals are functorial under isomorphism.  The
conductor (10) is the annihilator of a canonically transported finite module,
and divisorial valuations are transported with the boundary primes.  Local
monodromy is functorial for the induced finite étale covers of punctured
formal neighborhoods.

Under stabilization, the canonical triple is its product with `A^a`.
Intersections, scheme images, and reductions take products.  The different
pulls back because relative differentials do.  If `A subset B` is one of the
finite stratum-ring inclusions, then

\[
 \operatorname{Ann}_{A[t]}((B/A)[t])
 =\operatorname{Ann}_A(B/A)[t],                            \tag{13}
\]

so conductor ideals extend polynomially.  Every boundary valuation acquires
the Gauss extension with `v_E(t_i)=0`, preserving (11).  Finally, completion
along an extended ideal gives the completed cylinder

\[
 \widehat{R[t_1,\ldots,t_a]}_{,IR[t_1,\ldots,t_a]},       \tag{14}
\]

which is the stabilized formal datum; it need not be the ordinary polynomial
ring over `hat(R)_I`.  This proves the theorem.  QED

## 6. Computation strategy

The ladder suggests a fixed workflow for new examples:

1. compute `I^red` and stop if topology or components separate the maps;
2. compute `I^sch` only when reduced incidence agrees;
3. compute selected formal neighborhoods, conductors, differents, or
   valuation-graded pieces only when the scheme diagrams still agree.

This prevents the invariant from changing ad hoc between examples.  Stronger
layers refine a fixed canonical object; they do not replace earlier layers.

For weighted--cancellation, the reduced target layer already separates the families.  The
next formal-layer computation is useful for a different purpose: determining
the cover-lifting congruence kernel in the target-boundary automorphism group
and its possible action on the cancellation parameter roots.
