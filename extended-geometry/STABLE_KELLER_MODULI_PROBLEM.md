# Degreewise stable Keller moduli

The stable-class results are now strong enough that producing another
explicit family is not the main classification problem.  The next object is
the degreewise moduli groupoid itself.

Work over an algebraically closed field `k` of characteristic zero and let
`N>=3` denote **geometric degree**, not maximum coordinate degree.  The
guiding question is:

> **Stable Keller-moduli problem (`OP-KMOD`).**  Construct a useful moduli
> stack (or an ind-algebraic stack with finite-type degreewise pieces)
> \(\mathfrak K^{\mathrm{st}}_{3,N}\) of geometric-degree-`N` Keller maps
> \(\mathbb A^3\to\mathbb A^3\) modulo stable polynomial left--right
> equivalence.  Determine its dimension, its irreducible components, and the
> geometric monodromy group of the generic map on every component.  Do the
> same for its stably atomic locus
> \(\mathfrak A^{\mathrm{st}}_{3,N}\).

“Generic monodromy” is necessarily componentwise: if `C` is an irreducible
component and a universal map exists after passage to an atlas, take the
geometric monodromy of the induced degree-`N` function-field extension over
the generic point of `C`.  This is distinct from monodromy of loops in the
parameter space.

## 1. The candidate groupoid

For a `k`-scheme `S`, start with families

\[
 F:\mathbb A^3_S\longrightarrow\mathbb A^3_S
\]

whose Jacobian determinant is a unit and whose geometric degree is
fiberwise `N` on the chosen stratum.  Two such families are stably
left--right equivalent if, after adjoining a common number of identity
coordinates, they differ by polynomial source and target automorphisms over
`S`.  Stackify this groupoid in a stated topology.

This definition is only a prestack-level target, not a representability
claim.  Three sources of unboundedness have to be handled explicitly:

1. the coordinate degrees of the maps;
2. the degrees of the source and target automorphisms; and
3. the number of identity variables used in a stable equivalence.

Thus a classical finite-type quotient by a fixed reductive group should not
be assumed.  One possible construction is an ind-object assembled from
finite-type coefficient schemes and finite-type bounded-equivalence
correspondences.  Another is a finite-type boundary-decorated rigidification
whose forgetful morphism to the stable groupoid is representable on a clean
locus.  Either approach must prove that the resulting dimension and
component notions are independent of the chosen filtration.

The word “stack” also matters.  Stable self-equivalences produce inertia;
on the clean quadratic-gauge receiver the physical inertia is already known
to be trivial modulo automorphisms of the added identity factors, but that
vertical gauge has not been classified globally.

## 2. What is already known inside the candidate moduli object

The established results give loci and lower bounds, not a construction of
the whole stack.

### 2.1 The weighted locus

For every `N>=4`, the Hessian-clean weighted seed space has dimension
`N-3`.  The marked Hessian-divisor invariant

\[
 (\mathbb P^1;\operatorname{div}(H''),0,\infty)
\]

gives a generically etale stable-separation map of exact rerooting degree
`N-2` onto an image of dimension `N-3`.  Adding the intrinsic affine root
sheet removes that generic rerooting ambiguity.  Every map in this locus has
geometric monodromy `S_N` and is absolutely and stably atomic.

Consequently any adequate realization of
\(\mathfrak A^{\mathrm{st}}_{3,N}\) has dimension at least `N-3`.  What is
not known is whether this locus is open in a component, dense in a
component, or contained in a component of larger dimension.

### 2.2 The clean quadratic-gauge receiver

For `N>=5`, the clean marked receiver is the represented scheme

\[
 \mathscr D_N^{\mathrm{quad},\circ}
 \simeq Q_N^\circ,
 \qquad \dim \mathscr D_N^{\mathrm{quad},\circ}=N-1,
\]

with the honest span

\[
 BS_N\longleftarrow
 \mathscr D_N^{\mathrm{quad},\circ}
 \longrightarrow
 \mathcal M_N^{\mathrm{quad}}.
\]

The target-forgetting fibers have dimension three and
`dim M_N^(quad)=N-4`.  This is a represented chart for one decorated
compiler stratum, not a global moduli stack and not the source of the
`N-3` lower bound above.  Its role in `OP-KMOD` is to test an eventual
stable-map stack: construct its morphism to that stack, decide whether the
image is locally closed, and compute its normal deformation theory.

### 2.3 Cubics and decomposition loci

Every degree-three Keller nonunit is atomic.  Fiber-invisible cubic gauge
lifts give infinitely many stable atomic classes with any prescribed
complete cubic fiber.  This proves neither positive-dimensional cubic
moduli nor infinitely many irreducible components: a discrete stable
boundary invariant can jump in a degeneration.

For general `N`, degree multiplicativity and the stable Keller--Ritt theorem
suggest composition incidence correspondences

\[
 \mathfrak K^{\mathrm{st}}_{3,a}\times
 \mathfrak K^{\mathrm{st}}_{3,b}
 \longleftarrow \mathfrak{Comp}^{\mathrm{st}}_{a,b}
 \longrightarrow
 \mathfrak K^{\mathrm{st}}_{3,ab}.
\]

Here the middle object must remember the identification of the intermediate
affine space; composition is not a well-defined operation on two arbitrary
left--right classes without that datum.  Its image belongs to an
imprimitive-monodromy stratum.  A decomposable
nonunit exists exactly when `N=ab` with `a,b>=3`; the first such degree is
nine.  The Ritt machinery controls factorization and intersections for the
known separated-product loci, but it does not yet show that these are whole
components of the stable moduli object.

## 3. The degreewise invariants to compute

Once a satisfactory stack or coarse receiver is constructed, the primary
questions are the following.

1. **Dimension.**  Is \(\dim\mathfrak K^{\mathrm{st}}_{3,N}\) finite?  Is
   the atomic locus equidimensional?  Is the known `N-3` lower bound sharp
   on an `S_N` component?
2. **Components.**  Is the degree-`N` stack locally of finite type, and does
   it have finitely many irreducible components?  Which boundary
   normalization types, affine-root markings, and Ritt types label generic
   points rather than only special strata?
3. **Generic cover monodromy.**  For every component `C`, determine the
   transitive subgroup \(G_C\leq S_N\) of its generic degree-`N` cover.
   Which components have primitive `G_C`, hence lie in the known
   stable-atomic region, and which are generically imprimitive composition
   components?
4. **Inertia and separatedness.**  Determine the generic stable
   self-equivalence group, including the vertical stabilization kernel, and
   decide whether boundary-decorated rigidification gives a separated or
   adequately separated quotient.
5. **GIT model.**  Find a boundary linearization or invariant ring on a
   finite-type rigidification that separates stable classes on a dense
   locus.  The full polynomial-automorphism ind-group is nonreductive, so
   ordinary reductive GIT applies only after such a reduction has been
   justified.

The monodromy and Ritt strata should be studied together.  On a component
with generic monodromy `S_N`, primitive-monodromy atomicity is automatic.
On an imprimitive component, the block systems give candidate intermediate
fields; the stable Keller--Ritt theorem then asks which of them algebraize
through affine-space Keller factors.

## 4. The first finite-type approximation

The generic weighted locus already gives a first finite-type approximation.
This is a consequence of the existing decorated-normalization and
affine-mark theorems, but it is useful to state it in moduli language.

### 4.1 A bounded coefficient realization

Let

\[
 \mathcal A_N=
 \left\{
 H(W)=\sum_{j=3}^N h_j(W^j-W^2):
 \sum_{j=3}^N(j-2)h_j=-1
 \right\}.                                           \tag{4.1}
\]

Thus `A_N` is an affine space of dimension `N-3`.  Remove the degree-drop,
`H''(1)=-2`, non-squarefree, boundary-collision, and nonordinary closed
failure loci, and then shrink to a rerooting-stable open
`A_N^(mod)`.  The weighted formula is regular on this open and has
determinant one.  Degree counting in the intrinsic formula gives the
uniform coordinate bound

\[
 d_N=5N-8,
 \qquad
 (\deg F_{H,1},\deg F_{H,2},\deg F_{H,3})
 \leq(5N-8,5N-9,4),                                  \tag{4.2}
\]

with equality on a nonempty open.  Hence the universal weighted family is a
finite-type morphism

\[
 \iota_N:\mathcal A_N^{\mathrm{mod}}
 \longrightarrow X(3,d_N).                           \tag{4.3}
\]

After a further nonempty shrinking, `iota_N` is a locally closed immersion.
Indeed, equality of two coefficient points implies equality of the two maps,
hence stable equivalence, and affine-mark faithfulness recovers the seed.
The same argument works after algebraically closed field extension.  A
tangent vector in the kernel of `d iota_N` maps to zero under the marked
Hessian-divisor invariant; generic unramifiedness of that invariant forces
the seed tangent to vanish.  Thus `iota_N` is generically radicial and
unramified.  Applying Zariski Main to its scheme-theoretic image and
shrinking to the normal locus of the reduced image gives the asserted
locally closed immersion.

This embeds the normalized family in a genuine bounded coefficient scheme.
It does not quotient the other coefficient directions by stable
left--right equivalence.

### 4.2 The finite rerooting groupoid

For a nonzero simple root `a` of `H`, put

\[
 \kappa_a=-\frac1{aH'(a)},
 \qquad R_aH(w)=\kappa_aH(aw).                        \tag{4.4}
\]

On the rerooting-stable open these arrows form a finite etale equivalence
relation

\[
 \mathcal R_N\rightrightarrows\mathcal A_N^{\mathrm{mod}}.
                                                               \tag{4.5}
\]

The identity is `a=1`.  If `b` is a root of `R_aH`, then `ab` is a root of
`H` and

\[
 R_b(R_aH)=R_{ab}H;                                   \tag{4.6}
\]

in particular `b=1/a` gives the inverse arrow.  Each source fiber has the
`N-2` nonzero simple roots of `H`.

The marked Hessian-divisor map forgets exactly this relation.  After
shrinking its image `V_N` in

\[
 \mathcal Q_{N-2}=[U_{N-2}/\mathbb G_m],
\]

the exact rerooting theorem identifies

\[
 \mathcal A_N^{\mathrm{mod}}/\mathcal R_N
 \simeq V_N,                                         \tag{4.7}
\]

and the map `A_N^(mod) -> V_N` is finite etale of degree `N-2`.  Its generic
monodromy is `S_(N-2)`, as seen from the selected-root quotient

\[
 [M_{0,N}/S_{N-3}]\longrightarrow[M_{0,N}/S_{N-2}]. \tag{4.8}
\]

This is **parameter-cover monodromy**.  It is different from the `S_N`
geometric monodromy of the generic inverse Keller cover.

The quotient `V_N` is only a coarse boundary-decoration receiver.  Its
points identify the rerooted seeds even though those seeds have different
intrinsic affine sheets and therefore represent different stable Keller
classes.

### 4.3 The affine-marked stable-class core

Retaining the unique unramified root sheet in the reconstruction open kills
the relation (4.5): preservation of that sheet forces `a=1`.  Consequently

\[
 \mathcal A_N^{\mathrm{mod}}(K)
 \longrightarrow
 \{\text{degree-}N\text{ Keller maps over }K\}/
       \text{stable LR}                               \tag{4.9}
\]

is injective for every algebraically closed characteristic-zero extension
`K/k`.  More generally it is injective on reduced families: if two maps
from a reduced base to `A_N^(mod)` become stably left--right equivalent,
their restrictions to every geometric point agree by affine-mark
faithfulness, hence their coefficient morphisms agree on the reduced base.

Thus `A_N^(mod)` is a represented, finite-type, `N-3` dimensional
**reduced coarse core** of the desired stable-moduli object.  On this core:

- the generic inverse-cover monodromy is `S_N`;
- every point is absolutely and stably atomic;
- the generic parameter-cover monodromy before affine rigidification is
  `S_(N-2)`; and
- the physical stable-class map is faithful after affine rigidification.

This is the precise gain from the present moduli reframing.  It upgrades an
inequivalent family to a bounded represented test chart, while keeping clear
that no ambient stable quotient has yet been constructed.

## 5. The next gate: normal deformation after boundary rigidification

The next step is not another root family and not the unrestricted Artin
quotient.  It is to construct a simultaneous-boundary deformation functor
near (4.3).

Let `X_N=X(3,d_N)`.  On reduced test schemes, define
`X_N^(bdy)` to parametrize a coefficient family together with:

1. a simultaneous finite normalization of the inverse graph;
2. the distinguished affine open and the complete boundary divisor;
3. the ordered ramified and zero-boundary target vertices;
4. the relative Fitting divisor and conductor; and
5. the distinguished unramified affine root sheet.

### 5.1 Reduction to simultaneous normalization

At fixed `d_N`, finite typeness of the decorations after normalization is
not the main uncertainty.  Take the biprojective closure of the universal
graph in `P^3_source x P^3_target x X_N` and apply flattening stratification
to its projective morphism over `X_N`.  On a stratum carrying a simultaneous
normalization of that graph which commutes with reduced base change, take
the Stein factor of its projection to `P^3_target`.  Over the affine target
chart its algebra is finite of generic rank `N` and recovers the
Zariski--Main normalization package.

Once that algebra is available, the remaining data are finite-type:

- the source affine chart is the complement of the pulled-back source
  hyperplane;
- the boundary components and their target images are closed supports of
  finitely presented modules after a further flattening stratification;
- relative differentials, Fitting ideals, and conductor quotients are
  finitely presented and can be stratified for base change; and
- the affine root sheet is a length-one point of a finite cover, hence is
  parameterized by its relative length-one Hilbert scheme.

Thus the first representability problem reduces to one precise gate:

> **Equinormalization gate.**  Construct the maximal locally closed reduced
> stratum through `iota_N(A_N^(mod))` on which the normalized projective
> graph and its Stein algebra commute with base change, or prove that no
> such finite-type neighborhood contains every reduced transverse arc.

If the gate succeeds, ordinary flattening and Hilbert/Quot constructions
produce `X_N^(bdy)` locally.  If it fails, the correct moduli object must be
stratified by normalized-graph and boundary Hilbert polynomials; one cannot
use a single boundary chart across those arcs.

The immediate theorem to prove is representability, after shrinking around
`iota_N(A_N^(mod))`, by a finite-type algebraic space or stack, with base
change for the normalization, Fitting, conductor, and affine mark.  Only
then is there a legitimate normal cone

\[
 C_{\mathcal A_N^{\mathrm{mod}}/X_N^{\mathrm{bdy}}}   \tag{5.1}
\]

whose components can decide whether the weighted core is a component, a
divisor in a larger component, or a higher-codimension stratum.

The existing coefficient calculations calibrate but do not replace this
normal cone.  In degree four, the minimal box has raw tangent dimension
`58`; the affine orbit has rank `22`; adding four target shears and the one
weighted parameter gives a reduced family of rank `27`; and its generic raw
tangent space still has dimension `49`.  The resulting transverse excess
may be other reduced branches, nilpotent thickening, or directions that
disappear after boundary rigidification.  It is not a stable-moduli tangent
dimension.

### 5.2 The quartic biprojective graph calibration

The first graph calculation can now be made completely explicit.  For

\[
 H_4(w)=\frac14w^2(w-1)(w-5),
\]

the weighted map has coordinate profile `(12,11,4)`.  Let
`G_0,G_1,G_2,G_3` be the common degree-twelve homogenization of
`[1:F_1:F_2:F_3]`; after harmless diagonal rescaling of the target, these
have respectively `1,16,14,3` terms.  The closure of the affine graph is
presented by

\[
 I_{\Gamma_4}=
 (Y_1G_0-Y_0G_1,\;Y_2G_0-Y_0G_2,\;Y_3G_0-Y_0G_3)
 :(X_0Y_0)^\infty .                                 \tag{5.2}
\]

The source base scheme has a simple reduced support hidden by its
multiplicities.  Direct restriction gives

\[
 G_1|_{X_0=0}=\frac34X_1^6X_2^4X_3^2,
 \qquad
 \nu_{X_0}(G_0,G_1,G_2,G_3)=(12,0,1,8).
\]

Consequently its reduced support is the coordinate triangle

\[
 V(X_0,X_1X_2X_3)=L_1\cup L_2\cup L_3,
 \qquad L_i=V(X_0,X_i),                              \tag{5.3}
\]

and the generic `(X_0,X_i)`-adic orders of the four coordinates are

\[
\begin{array}{c|c|c}
 & (\nu_{L_i}(G_0),\nu_{L_i}(G_1),
     \nu_{L_i}(G_2),\nu_{L_i}(G_3)) & \nu_{L_i}(I_B)\\ \hline
 L_1&(12,6,7,11)&6\\
 L_2&(12,4,4,8)&4\\
 L_3&(12,2,3,9)&2.
\end{array}                                         \tag{5.4}
\]

Thus the normalization problem first stratifies into three generic-line and
three vertex boundary problems.  This six-stratum decomposition is the
useful replacement for an ambient calculation with all graph generators at
once; each stratum may still require several affine Rees charts.

There are two complementary exact calculations.

- Over `Q`, two independent rational source-plane/target-codimension-two
  cards, after dehomogenizing and saturating by the restricted base ideal,
  are zero-dimensional radical schemes of length `28`.
- Over `F_32003`, a full saturation of (5.2) has affine-cone dimension `5`,
  hence biprojective dimension `3`.  General source/target hyperplane cards
  give the modular projective multidegrees

  \[
    (\delta_0,\delta_1,\delta_2,\delta_3)=(1,12,28,4). \tag{5.5}
  \]

Here `delta_i` uses `i` target hyperplanes and `3-i` source hyperplanes.
The endpoint values agree with the birational source projection and the
degree-four generic inverse cover.  If (5.5) lifts to characteristic zero,
the Segre degree is

\[
  \delta_0+3\delta_1+3\delta_2+\delta_3=125.          \tag{5.6}
\]

This is deliberately a **calibration**, not yet a characteristic-zero
multidegree theorem.  The two exact cards certify reduced length `28` in the
rational source chart, while the full graph computation is at one good
prime; a flat spread-out graph or an exact characteristic-zero saturation is
still required to identify the middle projective multidegree over `Q`.
The modular standard basis has `181` elements and its saturated homogeneous
ideal has `38` minimal generators.  Thus a brute-force ambient
Jacobian-minor calculation is the wrong next interface.

The useful next object is the normalized Rees algebra of the base ideal
`(G_0,G_1,G_2,G_3)`, computed chartwise with its bigrading retained.  One
must glue the normalized boundary charts, compute their conductor and
bigraded Hilbert polynomial, and produce an integer model on which this
normalization is flat.  That simultaneously tests the equinormalization gate
and decides whether (5.5) lifts.  The command

```bash
.venv/bin/python scripts/verify_quartic_biprojective_graph.py
```

replays the exact section and modular graph calculations; it requires
Singular.

### 5.3 Making normalization finite: surface clusters and colored facets

The six-stratum decomposition does make the next calculation finite.  Over
the function field `K=Q(t)` of an open side `L_i`, the completed source is a
regular surface `K[[u,v]]`, and the restriction `I_i` of the four-generator
base ideal is primary to `(u,v)`.  Repeated quadratic transforms give the
following complete point bases.  A pair `(e;m)` records residue degree `e`
over `K` and base-point multiplicity `m`:

\[
\begin{array}{c|l|c|c|c}
 &\text{point basis }(e;m)&\ell(K[[u,v]]/I_i)
 &\ell(K[[u,v]]/\overline{I_i})
 &\ell(\overline{I_i}/I_i)\\ \hline
L_1&(1;6),(1;4),(1;1),(3;1),(1;2),(1;2),(2;1)
 &52&43&9\\
L_2&(1;4),(1;3),(1;1),(1;1),(1;1)
 &25&19&6\\
L_3&5\times(1;2),(2;1)&18&17&1.
\end{array}                                           \tag{5.7}
\]

The first colengths come from exact Groebner bases over `Q(t)`.  The point
bases terminate because the final leading forms are coprime, including
after the irreducible cubic and quadratic residue extensions in `L_1` and
the quadratic extension in `L_3`.  The middle column is then the
Hoskin--Deligne formula

\[
 \ell(K[[u,v]]/\overline I)
 =\sum_P[k(P):K]\binom{m_P+1}{2}.                    \tag{5.8}
\]

Thus the normalized blowup is completely determined over every open side;
the defect lengths `(9,6,1)` replace three large Rees eliminations.

At a vertex, take the monomial support of the four local generators and its
compact positive Newton facets.  There are six, not an unbounded search.
For `V_12`, use local coordinates `(a,b,c)=(X_0,X_1,X_2)` with `X_3=1` and
put `q_12=a^2+bc`.  For `V_23`, use
`(a,b,c)=(X_0,X_2,X_3)` with `X_1=1` and put

\[
 q_{23}=5ab-3c,\qquad h_{23}=3a^3-5ab+3c.
\]

The compact-facet and generic-color ledger is

\[
\begin{array}{c|c|c|c|c}
\text{vertex}&\text{weight}&\nu_w(I_B)&
 \text{torus color in the face ideal}&\text{generic point basis}\\ \hline
V_{12}&(1,1,1)&10&q_{12}^3&(1,1,1)\\
V_{12}&(2,3,1)&22&q_{12},\ P_{18}&(1),\ (1)\\
V_{23}&(1,1,2)&8&q_{23}^2&(2,2)\\
V_{23}&(1,2,1)&9&\varnothing&\varnothing\\
V_{23}&(1,2,3)&11&h_{23}&(1)\\
V_{23}&(2,3,5)&21&q_{23}&(1,1,1).
\end{array}                                           \tag{5.9}
\]

The last column of (5.9) lists multiplicities over the color's own function
field, rather than the `(e;m)` notation of (5.7).  Here `P_18` is the second
weighted-degree-`18` factor of the `(2,3,1)`
initial form.  On the torus chart `c=1,a=t`, it is irreducible of degree
five in `b`.  Its generic color chart has multiplicity one and coprime final
leading forms over that degree-five function-field extension.  Moreover
`P_18(t,-t^2,1)=t^8`, so the two `V_12` colors do not meet in the torus.
The third vertex `V_13` has no new compact positive facet beyond its adjacent
side data.

Therefore the remaining normalization work is no longer the `181`-element
ambient standard basis.  It begins with the finite collection of
toric-boundary and facet-overlap corner charts left after (5.7)--(5.9).  The
exact replay for the generic strata is

```bash
.venv/bin/python scripts/verify_quartic_rees_stratification.py
```

It verifies all quadratic transforms, algebraic terminal branches, Newton
facets, face colors, and generic-color principalizations symbolically.

### 5.4 The finite corner atlas and its cusp conductors

The corner calculation can also be completed without normalizing the ambient
graph.  Let `e_a,e_b,e_c` be the coordinate rays and use `w_ijk=(i,j,k)` for
a compact facet ray.  The complete lower Newton fan has the following eight
maximal cones.  The support column is their common lower-support vertex; the
last column gives a transformed generator which is a unit at the closed
toric orbit after the common divisorial monomial has been removed.

\[
\begin{array}{c|c|c|c|c}
 &\text{support}&\text{rays}&\text{lattice index}&\text{unit coefficient}\\ \hline
V_{12}&(0,6,4)&e_a,w_{111},w_{231}&2&G_1:3/4\\
 &(8,2,0)&e_c,w_{111},w_{231}&1&G_1:3/4\\
 &(11,0,0)&e_b,e_c,w_{231}&2&G_1:-3/2\\ \hline
V_{23}&(1,3,2)&w_{112},w_{121},w_{235}&2&G_2:1\\
 &(3,5,0)&e_c,w_{112},w_{235}&1&G_2:25/9\\
 &(8,0,1)&e_b,w_{121},w_{123},w_{235}&\text{nonsimplicial}&G_3:1\\
 &(9,1,0)&e_c,w_{123},w_{235}&1&G_3:-5/3\\
 &(11,0,0)&e_b,e_c,w_{123}&1&G_3:1.
\end{array}                                                   \tag{5.10}
\]

The unique nonsimplicial cone has primitive ray relation

\[
  2e_b+2w_{235}=w_{121}+3w_{123}.                 \tag{5.11}
\]

Thus it is retained as one normal toric chart rather than hidden by a choice
of small triangulation.  More importantly, the unit column in (5.10) proves
that no toric closed orbit is a residual base point.  The only remaining
vertex phenomena occur on the strict closures of the colors.

The difficult color is `P_18`.  At the `b=1` boundary point put
`d=c+a^2`.  Its first new face and residual equation are

\[
 \operatorname{in}_{(1,4)}P_{18}=a^{12}+\frac34d^3,
 \qquad d=z a^4:\quad 3z^3+4=0.                  \tag{5.12}
\]

Hence there are three geometric branches, one degree-three closed packet
over `Q`, with mutual contact four.  Its conductor in
`Q[a,c]/(P_18(a,1,c))` is

\[
 (c^4,\ c^2(a^2+c),\ (a^2+c)^2),
 \qquad \delta=\ell(A/\mathfrak c)=12.            \tag{5.13}
\]

At the `c=1` point there are two tangent packets.  The three branches tangent
to `b=0` split at the separable irreducible cubic

\[
 25z^3+45z^2+54z+46,
\]

and have mutual contact two.  The two branches tangent to `5a-3b=0` split
at

\[
 -15625z^2+15750z+756,
\]

and have mutual contact three; branches in different packets meet once.
Consequently

\[
 \delta=\binom32 2+(3\cdot2)1+3=15.              \tag{5.14}
\]

An exact conductor ideal with that colength is

\[
\begin{split}
(&25a^2b^2-30ab^3+9b^4,\ 81b^5-625ab^3+375b^4,\\
 &27ab^4-125ab^3+75b^4,\ 5a^4b-5ab^3+3b^4,\ a^6).
\end{split}                                        \tag{5.15}
\]

These branch ledgers also recover the two original color-contact lengths

\[
 P_{18}(a,1,-a^2)=a^{12},\qquad
 P_{18}(a,-a^2,1)=a^8.                             \tag{5.16}
\]

The full Newton fan removes most of those singularities.  On the unimodular
`e_c,w_111,w_231` overlap chart the strict `q_12` equation is `s=0`; writing
the strict `P_18` equation as `F_3=A_3y^2+B_3y+C_3`, one has

\[
 \operatorname{in}F_3=\frac14(3s^3+4y^2),\qquad
 B_3^2-4A_3C_3=\frac34s^3(7s-4).                  \tag{5.17}
\]

The normalization is therefore obtained locally by

\[
 Z=(2A_3y+B_3)/s,
 \qquad Z^2=\frac34s(7s-4).                       \tag{5.18}
\]

Its conductor is `(y,s)`, of colength one.  On the index-two
`e_a,w_111,w_231` cover the strict color is `x^2+1=0`, and for
`F_1=A_1y^2+B_1y+C_1` the corresponding formulas are

\[
\begin{split}
 B_1^2-4A_1C_1
   &=\frac34x^{10}(x^2+1)^3(3x^2+7),\\
 Z&=(2A_1y+B_1)/(x^5(x^2+1)),\\
 Z^2&=\frac34(x^2+1)(3x^2+7).
\end{split}                                        \tag{5.19}
\]

Here the conductor is `(y,x^2+1)`, of `Q`-colength two.  In conductor-pushout
form the two local gluing maps are simply

\[
 \mathbf Q\longrightarrow\mathbf Q[\epsilon]/(\epsilon^2),
 \qquad
 \mathbf Q(i)\longrightarrow\mathbf Q(i)[\epsilon]/(\epsilon^2),          \tag{5.20}
\]

with constants embedded on the right.  Thus the residual `P_18` conductor
contribution is `1+2=3`, not a new large normalization problem.

The `V_23` colored overlaps are normal already.  On the
`w_112-w_235` overlap the two `q_23` descriptions both reduce to `5-3s`.
On `w_235-w_123` one obtains

\[
 q_{23}^{\rm str}=5-3s,\qquad
 h_{23}^{\rm str}=3s+3u-5,                         \tag{5.21}
\]

so their restrictions at `u=0` agree up to sign and glue as one smooth
color.  The remaining `w_235-w_121` overlap gives `5v^2-3`, which is a unit
at `v=0`; that color misses the corner.

The command

```bash
.venv/bin/python scripts/verify_quartic_rees_corner_atlas.py
```

reconstructs (5.10)--(5.21) exactly.  It uses Singular only for the four
conductor ideals and colengths.

### 5.5 Compactifying the side packets and matching the vertices

The nontrivial side-center covers are also small enough to compactify
directly.  If `t=T_0/T_1` is the side parameter, a degree-`d` direction
packet which becomes constant in

\[
 z=t^k v/u
\]

has normalized branches

\[
 [U:V]=[T_0^k:zT_1^k].                            \tag{5.22}
\]

For `f(z)=\sum_{i=0}^d c_i z^i`, its image in
`P^1_[T_0:T_1] x P^1_[U:V]` is

\[
 F_f=\sum_{i=0}^d c_iT_0^{ki}T_1^{k(d-i)}V^iU^{d-i}=0,             \tag{5.23}
\]

of bidegree `(kd,d)`.  All geometric branches meet `[U:V]=[0:1]`
over `T_0=0` and `[1:0]` over `T_1=0`, with pairwise contact `k`.
The three packets extracted from (5.7) are

\[
\begin{array}{c|c|c|c|c}
\text{packet}&k&f(z)&\text{bidegree}&\delta\text{ at each endpoint}\\ \hline
L_1^{(1)}&1&(z+1)(25z^3+45z^2+54z+46)&(4,4)&6\\
L_1^{(2)}&3&-15625z^2+15750z+756&(6,2)&3\\
L_3^{(2)}&3&3z^2-12z+5&(6,2)&3.
\end{array}                                                       \tag{5.24}
\]

The first endpoint singularity is an ordinary four-branch point.  In local
base/direction coordinates `(r,s)` its conductor is `(r,s)^3`.  Each
quadratic packet has two smooth branches of contact three and conductor
`(s,r^3)`.  Exact normalization gives the displayed colengths at both
endpoints.  Thus the packet conductor at either endpoint of `L_1` has total
colength `6+3=9`, numerically equal to its generic closure defect in (5.7).
This equality is a consistency check, not yet an identification of the two
quotient sheaves.

For a divisor of bidegree `(kd,d)`, with `m` the base degree and `n` the
direction degree, the image and normalization Hilbert polynomials are

\[
\begin{split}
 P_{\rm image}(m,n)&=dm+kdn+(kd+d-kd^2),\\
 P_{\rm norm}(m,n)&=dm+kdn+d.
\end{split}                                                       \tag{5.25}
\]

Their difference is `kd(d-1)=2\delta_endpoint`.  Hence the three rows of
(5.24) give respectively

\[
\begin{array}{c|c|c|c}
 &P_{\rm image}&P_{\rm norm}&P_{\rm norm}-P_{\rm image}\\ \hline
L_1^{(1)}&4m+4n-8&4m+4n+4&12\\
L_1^{(2)}&2m+6n-4&2m+6n+2&6\\
L_3^{(2)}&2m+6n-4&2m+6n+2&6.
\end{array}                                                       \tag{5.26}
\]

All twelve rational infinitely-near centers are graphs of explicit
monomial maps `P^1 -> P^1`, of degrees zero, one, or two, with no common
zero in their two homogeneous coordinates.  Their compactifications are
therefore already normal and add no packet conductor.

The facet Kummer degrees are obtained without another normalization.  If
the side parameter is the `j`-th vertex coordinate, the facet chart of
weight `w` uses `t=lambda^{w_j}`.  The complete endpoint table is

\[
\begin{array}{c|c|c|c}
\text{side endpoint}&\text{vertex}&\text{parameter}&
 \text{degrees in facet order}\\ \hline
L_1(0)&V_{12}&c&(1,1)\\
L_2(0)&V_{12}&b&(1,3)\\
L_2(\infty)&V_{23}&c&(2,1,3,5)\\
L_3(\infty)&V_{23}&b&(1,2,2,3).
\end{array}                                                       \tag{5.27}
\]

Here the facet orders are those of (5.9).  At `V_13` there is no compact
positive facet and hence no additional Kummer row.

The algebraic endpoint identifications are exact.  At `V_12`, the cubic
factor of `L_1^(1)` and the quadratic `L_1^(2)` are precisely the two
`c=1` tangent packets of `P_18` in (5.14).  The `L_2` endpoint reaches the
`b=1` chart through the unique cubic base change

\[
 a=\lambda^2A,\qquad b=\lambda^3,\qquad c=\lambda C,              \tag{5.28}
\]

under which `q_12/lambda^4=A^2+C` and
`P_18/lambda^18=P_18(A,1,C)`.

At `V_13`, the two quadratic packets have the same residue field
`Q(sqrt(21))`.  Their roots are identified by

\[
 z_{L_1}=\frac9{25}z_{L_3}-\frac{27}{125},\qquad
 f_{L_1}(z_{L_1})=-675f_{L_3}(z_{L_3}).            \tag{5.29}
\]

At `V_23`, in the `L_3` coordinates
`(a,b,c)=(u/t,1/t,v/t)`, one has

\[
 t^2q_{23}=5u-3tv,qquad
 t^3h_{23}=3u^3-5tu+3t^2v.                        \tag{5.30}
\]

Thus the first `L_3` direction is exactly `q_23=0`; on it the second
expression is `3u^3`, recovering the cubic `q_23-h_23` contact resolved by
the subsequent centers.  The `L_2` endpoint is the `c=1` fixed direction,
where `q_23=-3` and `h_23=3`, so it misses both colors.

The command

```bash
.venv/bin/python scripts/verify_quartic_rees_side_corner_matching.py
```

reconstructs (5.22)--(5.30), including six exact `normalConductor`
calculations.  Before the global pushout, one must still determine how these
normalized center curves sit inside the rational exceptional **surfaces**.
In particular, the divisor classes and normal-bundle degrees of the
successive centers are not determined by (5.22)--(5.30).  Only after that
surface-intersection ledger is known can one impose the incidence maps in the
bigraded Rees algebra and take the single conductor pushout.  No
characteristic-zero graph multidegree or simultaneous-normalization statement
is inferred yet.

### 5.6 Pause checkpoint: the two-coefficient Hilbert-polynomial gap

The following records a computational checkpoint from August 10, 2026.  It
is deliberately not promoted to a theorem or to `MATH_STATUS.json`.

A Singular 4.4.1 calculation over `F_32003`, using the saturated modular
graph ideal from `WQG1` with

\[
 \deg X_i=(1,0),\qquad \deg Y_j=(0,1),
\]

and Hilbert-series denominator `(1-s)^4(1-t)^4`, returned

\[
 P_\Gamma(m,n)=\frac{1}{6}\bigl(
 m^3+36m^2n+84mn^2+4n^3+6m^2-54mn-30n^2
 -19m-268n+582\bigr).                                      \tag{5.31}
\]

Its cubic part recovers the modular multidegrees `(1,12,28,4)`.  In the Rees
grading `d=m+12n`, (5.31) is equivalent to

\[
 P_{R/B^n}(d)=
 (58n^2+33n+5)d-
 \frac{1226}{3}n^3-247n^2+\frac{20}{3}n-96.                 \tag{5.32}
\]

The exact normalized open-side point bases give the corresponding
degree-in-`d` coefficient

\[
 (33n^2+10n)+(14n^2+5n)+(11n^2+6n)=58n^2+21n.              \tag{5.33}
\]

Consequently, within the finite side-and-corner conductor model, every term
of the normalization correction involving the source degree is fixed.  The
candidate normalized polynomial is reduced to two target-direction
coefficients:

\[
 P_{\widetilde\Gamma}-P_\Gamma
   =12mn+5m+\alpha n^2+\beta n-96,                          \tag{5.34}
\]

or equivalently

\[
 P_{\widetilde\Gamma}(m,n)=\frac{1}{6}\bigl(
 m^3+36m^2n+84mn^2+4n^3+6m^2+18mn
 +(6\alpha-30)n^2+11m+(6\beta-268)n+6\bigr).              \tag{5.35}
\]

The specialization `n=0` is correctly
`binomial(m+3,3)`.  Equations (5.34)--(5.35) are a reduction of the pending
calculation, not a proof that the finite model is already the global
normalization.  The integers `alpha` and `beta` remain unknown.

A direct ambient normalization is not a useful way to determine them.  Two
scratch runs of Singular `normalI` on the modular Rees algebra were stopped:
`normalI(I,4)` had run for more than 36 minutes, while the localized route
`normalI(I,8,0,1)` spent about 29 minutes on its main closure and then more
than 70 minutes in two local workers, each using roughly 2--2.5 GiB.  Neither
run produced a result used here.  These timings are failed exploratory
computations, not certificates, and this ambient route should not be resumed
without a new reason.

The correct restart point is a **normal-bundle ledger** on the three ruled
exceptional surfaces.  For each of the twelve rational centers and three
algebraic packets, record its residue or covering degree, divisor class,
Rees multiplicity, restriction degree of the current polarization, and
`deg N_(C/X)`.  Then apply, center by center,

\[
 \chi(\pi^*L-rE)=\chi(L)
 -\frac{r(r+1)}2\bigl(\deg(L|_C)+1\bigr)
 +\frac{r(r-1)(r+1)}6\deg N_{C/X},                          \tag{5.36}
\]

sum normalized algebraic branches, and impose the endpoint conductor
corrections `(12,6,6)` (twice the endpoint deltas `(6,3,3)`).  The incidence
maps (5.27)--(5.30) can then be used to perform one global pushout.  That
finite Riemann--Roch calculation, rather than ambient normalization, is the
next route to `alpha` and `beta`.

This construction directly joins invariant theory to `OP-CCDM`.  The first
deliverables are therefore:

1. prove local representability of `X_N^(bdy)` along the weighted core;
2. compute the normal-bundle and divisor-class ledger for the center curves
   (5.22)--(5.30), then assemble them inside their rational exceptional
   surfaces, perform the resulting global conductor pushout, identify the
   equinormalization stratum, and compute the normalized-graph Hilbert
   polynomial in the quartic box;
3. construct its stable-equivalence correspondence and quotient the
   vertical stabilization gauge;
4. compute (5.1), beginning with `N=4` and `d_4=12`;
5. determine the generic inverse-cover monodromy on every resulting normal
   component; and
6. compare any imprimitive component with the Keller--Ritt composition
   incidence.

At a fixed automorphism-degree cutoff there is no group action to quotient:
polynomial automorphisms of degree at most `b>1` are not closed under
composition.  Bounded left--right loci are finite-type correspondences, not
groupoids.  The boundary rigidification is what reduces the generic weighted
relation to the finite groupoid (4.5) and makes ordinary invariant theory
available on this first chart.

## 6. Claims not made

- No global algebraic or finite-type stable Keller-moduli stack is currently
  constructed.
- The `N-3` weighted image and the `N-4` clean quadratic-gauge image count
  different loci; neither identifies the dimension of the whole moduli
  object.
- Infinitely many cubic stable classes do not determine the cubic moduli
  dimension or its number of components.
- `S_N` is the generic monodromy on the established weighted and
  quadratic-gauge loci, not on every hypothetical component.
- The existing Ritt classification does not enumerate the irreducible
  components of degree-`N` stable moduli.
- The locally closed coefficient slice (4.3) is not asserted to be an
  irreducible component of `X(3,d_N)` or of stable moduli.
- The reduced coarse core does not control nonreduced bases or classify
  stable inertia; those are part of the boundary-deformation gate.

The endpoint is therefore the geometry of
\(\mathfrak K^{\mathrm{st}}_{3,N}\) and
\(\mathfrak A^{\mathrm{st}}_{3,N}\), not a larger catalogue of pairwise
inequivalent examples.
