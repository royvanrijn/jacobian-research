# Recursive resonance screens and the log-étale H1-STACK comparison

## Result

This document proves `H1-STACK` by identifying the recursive screen
construction with the fs-saturated logarithmic base change of the
admissible-cover stack.  The bounded chart families in the reproduction
section test the formulas; exhaustive coverage and effective descent follow
instead from the selected-factor normalization lemma and the log-étale
comparison theorem.

Work in characteristic zero on a labelled splitting chart of the polynomial
root space.  Let \(\mathcal N\) be the laminar family of critical-value
collisions indexing a stratum of the normalized wonderful branch graph.
After an étale localization that separates the relevant source points and
the Kummer extensions forced by the node indices, the selected polynomial
admissible cover has the following atlas.

1. Every target cluster \(S\in\mathcal N\) has a framed screen coordinate
   \(T_S\), a placement \(a_S\), and a scale \(q_S\), with
   \[
      T_{\operatorname{par}(S)}=a_S+q_ST_S.       \tag{0.1}
   \]
2. If a source node over that target node has index \(e\), choose its source
   center \(b\), a source scale \(s\), and a unit \(u\) with
   \[
      q_S=u s^e.                                  \tag{0.2}
   \]
3. If \(R(U)\) is the map on the parent source screen, the child map is the
   regular normalized initial form
   \[
      R_{\rm child}(V)
      =
      \frac{R(b+sV)-a_S}{q_S}.                    \tag{0.3}
   \]
   Regularity means coefficientwise divisibility after the scale
   substitution; it does not require \(R-a_S\) to contain
   \((U-b)^e\) over the whole smoothing base.
4. For every finite target flag with screen residue \(\alpha_h\), its
   weighted inverse-image divisor is
   \[
      D_h=V\!\left(R_{\rm child}(V)-\alpha_h\right). \tag{0.4}
   \]
   The parent-node fiber is the pole divisor.  Thus all third-fiber
   equations are identities rather than additional choices.
5. Deleting a screen composes two affine substitutions.  Consequently all
   source maps and every divisor (0.4) contract functorially.
6. The stack character at a nested stratum is obtained from the full deck
   centralizer on each source vertex, including permutations of isomorphic
   disconnected components.  Match the induced permutations of source-node
   cycles across every target edge, add the two half-node rotations, and
   impose the diagonal Kummer condition at every target node.

These charts cover every nested tame degeneration of the generic polynomial
cover.  Radial root collisions are the anchored instances already described
by the general radial atlas; nonradial Maxwell and caustic refinements are
the transposition-tree instances; mixed strata use the same recursion.

> **Recursive polynomial atlas theorem (`H1-STACK`).**
> The normalized closure of the labelled generic polynomial cover over the
> normalized wonderful target graph is the tame Deligne--Mumford stack
> obtained from (0.1)--(0.4), admissible node saturation, and the simultaneous
> centralizer character rule.  The construction is compatible with every
> nested-set contraction and descends from collision-adapted splitting charts.
> Its coarse space is the wonderful target graph.

The selected-factor lemma and log-étale comparison in Section 7 prove the
global identification, exhaustive coverage, and effective descent at once.
The previous finite-normalization theorem independently proves
`H1-COARSE`.  The original identification with the unmodified root-stable
quotient remains false, and the degree-five chart still requires the
normalized blowup of \((x^3,y^2)\).

The target-screen construction is the local “screen” form of the
Fulton--MacPherson/wonderful compactification.  Li's arrangement theorem
constructs such models by smooth blowups and includes the
Fulton--MacPherson spaces as examples:
[Wonderful compactification of an arrangement of
subvarieties](https://arxiv.org/abs/math/0611412).
The source Kummer equations are the standard normalization charts of
admissible covers; see Abramovich--Corti--Vistoli,
[Twisted bundles and admissible
covers](https://arxiv.org/abs/math/0106211).
Compatibility with the tropical edge-length and expansion-factor picture is
the classical/tropical branch compatibility of
Cavalieri--Markwig--Ranganathan,
[Tropicalizing the space of admissible
covers](https://arxiv.org/abs/1401.4626).

## 1. Framed nested screens

Let \(I\) be a finite set of labelled target flags and let
\(\mathcal N\) be a laminar family of subsets:

\[
 S,S'\in\mathcal N
 \quad\Longrightarrow\quad
 S\subset S',\ S'\subset S,\ \text{or }S\cap S'=\varnothing. \tag{1.1}
\]

The children of \(S\) are its maximal proper members in \(\mathcal N\)
together with the singleton labels not contained in one of those members.
Introduce a placement \(a_{S,C}\) for every child \(C\) and a scale \(q_C\)
for every cluster child.  Recursively put

\[
 z_{S,i}=
 \begin{cases}
   a_{S,i},&i\text{ is a singleton child},\\[2mm]
   a_{S,C}+q_Cz_{C,i},&i\in C.
 \end{cases}                                      \tag{1.2}
\]

At the root, \(z_{\varnothing,i}\) is the blowdown position of the original
target flag.  On the boundary \(q_C=0\), all labels in \(C\) have the same
position on the parent screen while their relative positions remain visible
on the \(C\)-screen.

The coordinates in (1.2) retain an affine frame on every screen.  Because a
non-root target component has its parent flag and at least two child flags,
one may quotient the affine frame, or set two child placements to \(0\) and
\(1\).  Keeping the frame is more efficient: all contraction maps are
polynomial and no case distinction is caused by changing normalizing flags.

This redundancy has an explicit gauge action.  Change the coordinate on
each screen \(S\) by

\[
 u'_S=c_S+d_Su_S,\qquad d_S\in\mathbf G_m.
\]

If \(C\) is a cluster child of \(S\), then

\[
 q'_C=\frac{d_S}{d_C}q_C,\qquad
 a'_{S,C}=c_S+d_Sa_{S,C}-q'_Cc_C.                 \tag{1.2a}
\]

For a singleton child \(i\), the rule is simply

\[
 a'_{S,i}=c_S+d_Sa_{S,i}.                          \tag{1.2b}
\]

Substitution in (1.2) cancels every internal \(c_C,d_C\) and changes a root
leaf position only by \(z'_i=c_{\rm root}+d_{\rm root}z_i\).  Thus the open
locus where child placements are pairwise distinct, modulo
\(\prod_S\operatorname{Aff}_1\), is the ordinary stable-screen chart.
Equations (1.2a)--(1.2b) make the framed presentation a smooth atlas rather
than extra moduli.

### Proposition 1.1: exact screen contraction

Suppose \(C\) is a child of \(S\), and delete the \(C\)-screen.  Promote each
child \(D\) of \(C\) to a child of \(S\) and set

\[
 \begin{aligned}
   a'_{S,D}&=a_{S,C}+q_Ca_{C,D},\\
   q'_D&=q_Cq_D
       &&\text{if \(D\) is itself a cluster}.       \tag{1.3}
 \end{aligned}
\]

Then every blowdown position \(z_{\varnothing,i}\) is unchanged.  Successive
contractions are independent of parenthesization.

### Proof

Equation (1.3) is the composition

\[
 u_S=a_{S,C}+q_Cu_C,\qquad
 u_C=a_{C,D}+q_Du_D.
\]

It gives

\[
 u_S=(a_{S,C}+q_Ca_{C,D})+(q_Cq_D)u_D.
\]

The equality of all leaf positions follows by recursion.  Associativity of
affine composition proves path independence. \(\square\)

The pulled-back target placements and scales are the requested recursive
root-residue coordinates: the branch values and their diagonal ideals are
polynomials in the labelled polynomial roots, and the normalized wonderful
pullback makes their successive ratios regular.

## 2. Source screens from normalized initial forms

Let \(A\) be the normal local ring of a target-screen chart and let

\[
 F(W)\in A[W]
\]

be the polynomial reconstructed from its labelled zero fiber.  The same
argument applies recursively with \(F\) replaced by a source-vertex map
\(R\).

Let \(a\) be the center of a target screen and let \(b\) be a source point
above it with expansion index \(e\).  The normalized admissible-cover chart
adjoins a source scale \(s\) satisfying

\[
 q=us^e,\qquad u\in A^\times.                       \tag{2.1}
\]

After the corresponding étale and Kummer localization, Taylor expansion has
the form

\[
 F(b+sU)-a=\sum_{k\ge0}d_kU^k.                      \tag{2.2}
\]

On the pullback of the selected normalized admissible-cover chart, regularity
of the universal map is equivalent to the coefficient relations

\[
 d_k=qh_k\qquad(k\ge0).                             \tag{2.3}
\]

These relations are the Rees-chart form of initial-form divisibility.  They
hold because the universal map to the expanded target is regular.  They are
not being used here to define the normalization: the comparison theorem in
Section 7 identifies the chart first from the standard admissible-cover
deformation ring, after which (2.3) extracts regular functions on that
chart.  This avoids assuming that coefficientwise divisibility by itself
generates the full integral closure.  Tame local admissibility says on the
boundary that

\[
 \bar h_0=\cdots=\bar h_{e-1}=0,
 \qquad \bar h_e\in A^\times.                       \tag{2.4}
\]

Absorb \(u\) into either framed coordinate.  The source-screen map is
therefore

\[
 R(U)=\frac{F(b+sU)-a}{q}
     =\sum_{k\ge0}h_kU^k.                           \tag{2.5}
\]

Equations (2.3)--(2.4), rather than an exact factorization of \(F-a\) over
the smoothing base, prove regularity and the required expansion index.  A
lower Taylor term may be nonzero but acquire an extra factor of \(s\) after
division by \(q\); it then vanishes in the initial form.  Applied at every
cycle of the boundary monodromy, (2.5) constructs all source vertices over
the new target screen.

This description includes:

- a degree-\(\mu\) radial connector, where the initial map is \(U^\mu\)
  up to a unit;
- a local polynomial tail, where the remaining initial terms retain the
  cluster-root residues;
- an identity strand, where \(e=1\);
- a Maxwell component, where disjoint transpositions give several quadratic
  initial maps; and
- a caustic component, where adjacent transpositions belong to one
  higher-degree initial map.

The monodromy-forest theorem determines which source point \(b\), degree,
and index occurs.  Formula (2.5) supplies its algebraic map and positions.

## 3. Flag divisors and the third-fiber equations

Let \(\alpha_h\) be the residue of a finite target flag \(h\) on its screen.
On a source component with map \(R(U)\), define

\[
 D_h=\operatorname{div}_0(R(U)-\alpha_h).           \tag{3.1}
\]

For the parent flag at infinity, use the pole divisor of \(R\).  If a
different target coordinate writes the map as \(A(U)/B(U)\), the equation is

\[
 A(U)-\alpha_hB(U)=0.                               \tag{3.2}
\]

These are explicit equations over the root-residue chart.  Factorization
after the same étale splitting cover records the individual source points;
without splitting, (3.1) is already the weighted effective divisor.

Two consequences close the former flag-extension problem.

1. Every flag divisor is pulled back from one and the same map \(R\), so the
   third-fiber compatibility equation is automatic.
2. Any two flag fibers and a third normalization recover \(R\) by
   source-vertex rigidity.  Thus (3.1) introduces no additional Hurwitz
   choice.

The node fibers from the two adjacent source screens are the same
Weierstrass factors of \(F-a\).  Their multiplicities are the common cycle
lengths, so the source points glue and satisfy admissibility.

## 4. Contraction compatibility

Consider two successive source substitutions

\[
 W=c+rU,\qquad U=b+sV
\]

and target normalizations

\[
 R(U)=\frac{F(c+rU)-a}{t},\qquad
 R_{\rm child}(V)=\frac{R(b+sV)-\beta}{u},
\quad \beta=R(b).
\]

Then

\[
\begin{aligned}
 R_{\rm child}(V)
 &=\frac{F(c+r(b+sV))-(a+t\beta)}{tu}.              \tag{4.1}
\end{aligned}
\]

But the composite source coordinate is

\[
 W=(c+rb)+(rs)V.                                    \tag{4.2}
\]

Equations (4.1)--(4.2) say that recursively taking normalized initial forms
is exactly the same as first contracting the intermediate screen and then
taking the direct normalized form.  Applying (4.1) to
\(R_{\rm child}-\alpha_h\) proves the same assertion for every flag divisor.

> **Initial-form transitivity theorem.**  
> The source maps, all weighted target-flag inverse images, node
> identifications, and expansion indices commute with every wonderful
> edge contraction.

This is an identity before specialization.  It is stronger than agreement
of the limiting dual graphs.

## 5. Automatic simultaneous characters

Choose branch paths adapted to the laminar family, so each collision cluster
is represented by an interval of the branch-cycle tuple

\[
 \gamma_1,\ldots,\gamma_m\in S_d.                  \tag{5.1}
\]

For a target vertex \(S\), replace every child cluster \(C\) by its boundary
product

\[
 \gamma_C=\prod_{i\in C}^{\rm path\ order}\gamma_i, \tag{5.2}
\]

and retain the direct branch cycles in \(S\) not lying in a child.  These
permutations are the local monodromy generators on that target screen.
Their full centralizer in \(S_d\) is the vertex deck group.

It is important to take the **full** centralizer, rather than a product of
centralizers on connected components.  At a nested boundary, a deck
transformation may exchange two isomorphic unlabelled degree-one components.
The exchange may be forced to continue across the adjacent vertex and can
contribute to inertia.

For every incident target node:

1. decompose its boundary product into oriented cycles;
2. record how a vertex deck transformation permutes those cycles;
3. record its rotation exponent on each cycle;
4. match the cycle permutation with the action from the adjacent target
   vertex; and
5. add the two half-node rotations modulo the cycle length.

The resulting rows are the actual global automorphism character

\[
 \chi:\operatorname{Aut}_{\rm lab}(C/P)
 \longrightarrow
 \prod_{\text{source nodes }j}\mu_{e_j}.            \tag{5.3}
\]

At each target node, intersect (5.3) with the diagonal
\(\mu_{\operatorname{lcm}(e_j)}\), simultaneously for all target nodes.
This is precisely the labelled-node saturation theorem.

The compiler retains more than the order of this intersection.  A global
automorphism is stored as the concrete tuple

\[
 (g_{\rm root},(g_S)_{S\in\mathcal N})\in
 \prod_{\text{target vertices}}S_d,
\]

and the inertia group is the subset whose phase row is diagonal at every
node.  Multiplication is vertexwise permutation composition.  Thus kernels,
component exchanges, and noncyclic inertia are not collapsed by recording
only the character image.

Changing a collision-adapted path system applies Hurwitz moves and
simultaneous conjugations to (5.1).  Centralizers, cycle matchings, and
rotation characters are transported isomorphically.  Hence the result is
independent of the paths and descends under braid and root-label actions.

### Pair--triple correction

For the nested degree-six matching

\[
 S_2=\{1,2\}\subset S_3=\{1,2,3\},
\]

the two node profiles are

\[
 (1,1,2,2),\qquad(2,2,2).                           \tag{5.4}
\]

The full compiler finds \(32\) compatible global cover automorphisms on the
special fiber, of which exactly \(4\) preserve one normalized Kummer branch.
Thus the inertia is

\[
 (\mu_2)^2.                                         \tag{5.5}
\]

A componentwise centralizer calculation incorrectly returns order \(2\):
it loses the automorphism that exchanges the two unlabelled identity
components on the inner bubble.  This is why the cycle-matching step is
essential.

## 6. Radial, resonance, and mixed charts

On a target screen containing the labelled zero flag, every sheet of its
zero fiber is anchored.  Its vertex deck action is therefore the identity.
This does **not** by itself kill a rotation on a different two-point power
connector: such a rotation can be an automorphism of the nodal special
fiber.  The correct test matches it across the whole chain and imposes the
diagonal Kummer condition at every node.  For the polynomial radial atlas,
the other endpoint is rigidified by the stationary monodromy.  This does not
always make normalized inertia trivial when cluster multiplicities are
unequal.  The full-chain answer has a closed form.  Write the ordered scale
partition, with \(B_0\) the first active block, as

\[
 B_0\mid B_1\mid\cdots\mid B_k,
\]

and put

\[
 M_j=\operatorname{lcm}_{i\in B_j}\mu_i,\qquad
 L_j=\operatorname{lcm}_{i\in B_j\cup\cdots\cup B_k}\mu_i. \tag{6.1}
\]

If \(t_j\in\mathbf Z/L_j\) is the diagonal phase at the node before level
\(j\), rigidification of a cluster \(i\in B_\ell\) gives the telescoping
condition

\[
 t_0+\cdots+t_\ell=0\pmod{\mu_i}.                  \tag{6.2}
\]

After the previous phases have been fixed, (6.2) leaves exactly
\(L_j/M_j\) choices at level \(j\).  Conversely, integrating the successive
phases along each power-connector chain recovers a unique tuple of connector
rotations because both endpoint actions are fixed.  Thus every solution is
realized and no additional kernel remains.  Hence

\[
 \boxed{\;
 |I_{\rm radial}|
   =\prod_{j=0}^{k}\frac{L_j}{M_j}.
 \;}                                                \tag{6.3}
\]

This is the optimized radial character formula.  Equal multiplicities make
every factor one, so all thirteen sextic \((2,2,2)\) scale types have
trivial normalized inertia even though their special fibers can have
\(2,4,\) or \(8\) compatible cover automorphisms.  Unequal multiplicities
can retain stack structure: the strict orders
\((2)\mid(3)\) and \((3)\mid(2)\) have inertia orders \(3\) and \(2\),
respectively, while \((2)\mid(4)\) has order \(2\).

For comparison, an arbitrary tame chain without the rigid polynomial root
vertex can retain ghost inertia; the checked profile
\((4,2)\rightsquigarrow(6)\) has order \(6\).  Thus “a downstream label kills
every connector rotation” is not a valid standalone principle.  Full
two-ended cycle matching is essential.

Away from zero, the generic branch cycles are transpositions.  Their
subforests give the source components and the character compiler gives all
nested Maxwell/caustic inertia.

At a mixed stratum, include the zero-fiber monodromy—with its arbitrary tame
cycle profile—among the local branch cycles and mark the vertex carrying the
zero flag as anchored.  The same full-centralizer fiber product then joins
the radial and resonance charts.  No separate mixed rule is necessary.

## 7. Global comparison and descent

This section proves the comparison with the normalized selected
admissible-cover closure.  In particular, it does not infer equality of two
atlases merely from agreement of their displayed formulas.

### 7.1 The globally defined stack

Let

\[
 B=B_N^{\mathrm{tgt}}
\]

be the normal wonderful target graph, and let

\[
 \operatorname{br}:
 \overline{\mathcal H}^{\mathrm{adm,lab}}_N
 \longrightarrow \overline M_{0,m}
\]

be the fully branch- and fiber-marked tame admissible-cover stack with the
polynomial discrete data.  Put

\[
 \mathcal Z
 =
 B\times_{\overline M_{0,m}}
 \overline{\mathcal H}^{\mathrm{adm,lab}}_N.
                                                               \tag{7.1}
\]

The branch morphism is proper with finite geometric fibers; after passing
to coarse spaces its branch morphism is finite.  (It need not be
representable at a boundary cover with relative inertia.)  Let
\(\mathcal W\subset\mathcal Z\) be the reduced closure of the generic
polynomial section and define

\[
 \mathcal G=\operatorname{Norm}(\mathcal W).        \tag{7.2}
\]

All stacks here are of finite type in characteristic zero.  Normalization is
therefore finite.  It follows that \(\mathcal G\) is a proper, separated,
tame Deligne--Mumford stack over \(B\), with finite coarse morphism.  Its
generic stabilizer is trivial.  The finite-normalization theorem proves
that its coarse space is \(B\).
Thus algebraicity, the DM property, separatedness, and the underlying coarse
comparison are consequences of the global construction (7.2), not of a
formal gluing assertion.

### 7.2 Log-étale comparison theorem

Give \(B\) the divisorial fs log structure of its wonderful boundary and
give the admissible-cover stack the divisorial fs log structure of its
nodal boundary.  At a geometric point, let \(P\) be the characteristic
monoid of \(B\), let \(E\) be the set of target nodes, and write
\(\mathbf N^E\to P\) for the target smoothing map.  If the source nodes
over \(\eta\in E\) have indices \(e_{\eta1},\ldots,e_{\eta r_\eta}\), set

\[
 Q_\eta
 =
 \left\langle
   \delta_\eta,\sigma_{\eta1},\ldots,\sigma_{\eta r_\eta}
   \ \middle|\
   \delta_\eta=e_{\eta j}\sigma_{\eta j}\ (1\le j\le r_\eta)
 \right\rangle
\]

and \(Q=\bigoplus_{\eta\in E}Q_\eta\), where the \(\eta\)-th generator of
\(\mathbf N^E\) maps to \(\delta_\eta\).  The normalized pullback monoid is

\[
 P_{\mathrm{cov}}
 =
 \left(P\oplus_{\mathbf N^E}Q\right)^{\mathrm{sat}}. \tag{7.2a}
\]

The following toroidal lemma is the point that turns the screen formulas
into a comparison proof.

> **Lemma 7.1 (selected-factor normalization).**
> Let \(X\) be a normal excellent fs logarithmic stack, let \(Y\to X\)
> be Kummer log-étale and finite on coarse spaces, and let
> \(U\subset X\) be the dense trivial-log locus.  If a section
> \(U\to Y_U\) selects a generic factor, then the normalization of its
> reduced closure in \(Y\) is Kummer log-étale over \(X\).  Strict étale
> locally, it is exactly the normal toroidal factor of the saturated
> pushout chart selected by the generic factor.

### Proof

Pass to a strict henselian toroidal chart
\(\operatorname{Spec}A\to X\).  Strict étale locally, a Kummer log-étale
map is strict étale over a Kummer toric chart.  Its unsaturated coordinate
ring has the form
\[
 A\otimes_{\mathbf Z[P]}\mathbf Z[P']
\]
for a Kummer map \(P\to P'\).  The normalization of an affine monoid
algebra is the algebra of the saturated monoid, so normalizing this chart
replaces the monoid pushout by its saturation.  Normalization commutes with
strict étale base change.  Finally, the normalization of a reduced ring is
the product of the normalizations in the function fields of its minimal
primes.  The generic section selects one of those fields; its reduced
closure has exactly that function field, and its normalization is therefore
the corresponding saturated toroidal factor.  Selection introduces neither
extra integral equations nor extra normalization branches.  The
construction is invariant under strict étale change and hence descends.
\(\square\)

For the two normalization steps used here, see the Stacks Project on
[normalization by generic factors](https://stacks.math.columbia.edu/tag/035P)
and [normalization under smooth, hence strict étale, base
change](https://stacks.math.columbia.edu/tag/082D).  The monoid statement is
the standard affine-toric normalization theorem; for example, it is
Proposition 1.3.8 in Cox--Little--Schenck,
[Toric Varieties](https://web.ma.utexas.edu/users/ikmartin/pages/blog/toric-page/documents/Toric-Varieties_Cox-Little-Schenck.pdf).

> **Theorem 7.2 (log-étale comparison).**
> Strict étale locally at every geometric point of \(\mathcal G\), each
> selected polynomial factor is the toroidal chart associated with
> \(P\to P_{\mathrm{cov}}\), together with one discrete Hurwitz class.
> These charts are Kummer log-étale over \(B\), are jointly surjective, and
> have no additional integral equations or normalization branches.

### Proof

Fix a geometric point \(x\to\mathcal G\), with target nodes indexed by
\(\eta\).  Write \(q_\eta\) for the corresponding target smoothing
functions on an étale neighborhood of its image in \(B\).  Suppose the
source nodes over \(\eta\) have expansion indices

\[
 e_{\eta1},\ldots,e_{\eta r_\eta}.
\]

After strict étale localization, pull the Harris--Mumford deformation
algebra back to the completed local ring \(A\) of \(B\).  It is

\[
 A[[s_{\eta j}]]\Big/
 \left(s_{\eta j}^{e_{\eta j}}-u_{\eta j}q_\eta\right)_{\eta,j},
 \qquad u_{\eta j}\in A^\times .                   \tag{7.3}
\]

There are no further relative deformation parameters: with all branch
points marked, a local Hurwitz class over a fixed target is discrete.  This
is the standard admissible-cover deformation theorem.  Its normalization
is the ACV twisted-cover chart.  Equivalently, the branch morphism is
étale-locally toroidal, and (7.3) is its characteristic-monoid chart.
See Abramovich--Corti--Vistoli,
[Twisted bundles and admissible
covers](https://arxiv.org/abs/math/0106211), and
Cavalieri--Markwig--Ranganathan,
[Tropicalizing the space of admissible
covers](https://arxiv.org/abs/1401.4626), Sections 4.2.1--4.2.2.
There the Harris--Mumford chart is explicitly possibly nonnormal toric, its
normalization has
\(\prod_j e_{\eta j}/\operatorname{lcm}_j(e_{\eta j})\) branches at
\(\eta\), and the branch morphism is toroidal.  Its cone map is an
isomorphism after forgetting the integral structures; hence the
characteristic-monoid map is Kummer, while the remaining strict part is
étale.  This is the Kummer log-étale assertion used below.

Normalization of a toroidal base change is saturation of its pushout
monoid.  Thus the completed local ring of the normalized pullback is the
completed monoid algebra of (7.2a), with the same strict-étale parameters
as \(A\).  In other words, it is the fs logarithmic fiber product.  The
admissible-cover chart is Kummer log-étale over the target chart, and this
property is preserved by fs base change.  After strict henselization, the
generic polynomial section selects a union of normal local factors.  That
selection changes neither the characteristic monoid nor the deformation
equations, by Lemma 7.1.

For one target node, put

\[
 L_\eta=\operatorname{lcm}_j(e_{\eta j}).
\]

Before this base change, over the universal target deformation ring, the
normalized factors of the one-node algebra are parametrized by

\[
 q_\eta=\rho_\eta^{L_\eta},\qquad
 s_{\eta j}
 =\omega_{\eta j}\rho_\eta^{L_\eta/e_{\eta j}},
 \qquad
 \omega_{\eta j}\in\mu_{e_{\eta j}},               \tag{7.4}
\]

modulo simultaneous reparametrization of \(\rho_\eta\).  Roots of the units
have been absorbed after strict henselization.  After base change to \(A\),
one must normalize
\(A[\rho_\eta]/(\rho_\eta^{L_\eta}-q_\eta)\); this incorporates any further
monoid saturation forced by the wonderful pullback.  For several target
nodes, take the tensor product of these node factors over \(A\) and
normalize.  A geometric point \(x\) selects one factor at every node and
one of the finitely many local Hurwitz classes.

Choose an actual étale toroidal chart of the admissible-cover stack at the
image of \(x\), pull it back to \(B\), take the reduced polynomial closure,
and normalize each component meeting the generic polynomial section.
Normalization commutes with étale localization, so these normalized schemes
map étale to \(\mathcal G\); allowing \(x\) to vary makes them jointly
surjective.  Excellence identifies their completed local rings with the
selected normalized factors of (7.3).  Conversely every selected factor
that meets the generic polynomial section occurs in \(\mathcal G\), by the
definition of reduced closure.  This proves chart exhaustion and rules out
additional integral equations or unrecorded normalization branches:
Lemma 7.1 identifies each normalized closure with its selected saturated
toroidal factor.  This proves Theorem 7.2. \(\square\)

### 7.3 Identification with the recursive formulas

On one of the étale charts above, split the relevant source points and
choose affine coordinates on every source and target component, with the
parent flag at infinity.  The stable-target chart then has

\[
 T_{\operatorname{par}(S)}=a_S+q_ST_S.
                                                               \tag{7.5}
\]

At a source node of index \(e\), (7.3)--(7.4) give
\(q_S=us^e\).  Regularity of the universal admissible-cover map to the
expanded target says coefficientwise that

\[
 F(b+sU)-a_S\in q_S A[U].
\]

Therefore its expression on the child screen is exactly

\[
 R_{\mathrm{child}}(U)
 =\frac{F(b+sU)-a_S}{q_S}.                          \tag{7.6}
\]

This proves (0.3) on an actual étale chart of \(\mathcal G\), rather than on
a separately postulated formal chart.  Pulling back a finite target flag
\(\alpha_h\) gives

\[
 D_h=V(R_{\mathrm{child}}-\alpha_h),                \tag{7.7}
\]

and the parent flag gives the pole divisor.  Thus the flag equations
(0.4), the third-fiber compatibility, and the vertex maps are restrictions
of the universal cover on \(\mathcal G\).  The monodromy-forest theorem
identifies the source components and their degrees; source-vertex rigidity
then shows that (7.7) leaves no second component map with the same
flag-complete data.

Every point of \(\mathcal G\) admits the construction above, so
(7.5)--(7.7), together with the saturated node equations (7.3), are an
étale logarithmic atlas of \(\mathcal G\).

### 7.4 Overlaps and effective descent

The charts of Section 7.3 are restrictions of étale charts of the single
global stack \(\mathcal G\).  Their overlap groupoids and cocycles are
therefore effective automatically.  The explicit transition functions are
as follows.

- A change of target or source frame is affine because it preserves the
  parent flag.  Under \(U'=c+dU\), \(T'=a+bT\), the vertex map becomes
  \[
   R'(U')=a+bR\!\left(\frac{U'-c}{d}\right).         \tag{7.8}
  \]
  Equations (1.2a)--(1.2b) give the corresponding screen transitions.
- Changing splitting labels permutes the indexed placements, source
  factors, and flag divisors.
- A change of collision-adapted paths changes only the topological
  trivialization of the same finite étale cover.  Hurwitz moves and
  simultaneous conjugation transport the monodromy centralizers and node
  cycles, so no new algebraic branch is introduced.
- Replacing a Kummer generator by a root-of-unity multiple is the phase
  action in (7.4).

Affine composition proves compatibility under contraction.  Since \(B\) is
defined as a normalized graph closure, and the building-set blowup and
normalization are functorial, this descent is independent of the ordering of
the wonderful centers and of every auxiliary splitting order.

### 7.5 Actual inertia

The fully marked stable target has no automorphisms.  On the normalization
of a source component over the punctured target vertex, an automorphism of
the cover is exactly a permutation of sheets centralizing the local
monodromy.  Taking the centralizer in the full symmetric group retains
permutations of isomorphic disconnected components.  A global
automorphism is exactly a tuple of such vertex automorphisms whose
permutations of incident node cycles agree on both sides of every edge,
together with the balanced half-node rotations.  This is the
full-centralizer fiber product of Section 5.
Restriction to the punctured source components is injective because those
opens are dense in the nodal source.  Conversely, compatible vertex
automorphisms with matched half-node actions glue uniquely across the
completed node rings.  Thus the fiber product is neither a subgroup estimate
nor merely an equality of orders: it is the full automorphism group of the
admissible cover.

Its action on the source smoothing parameters gives

\[
 \chi:
 \operatorname{Aut}_{\mathrm{lab}}(C/P)
 \longrightarrow
 \prod_{\eta,j}\mu_{e_{\eta j}}.                   \tag{7.9}
\]

An automorphism fixes the selected normalized factor (7.4) precisely when,
at every target node \(\eta\), its phase vector lies in

\[
 \Delta_{\eta}(\mu_{L_\eta})
 \subset\prod_j\mu_{e_{\eta j}}.
\]

Consequently the stabilizer of \(x\) in the actual normalized stack is

\[
 I_x
 =
 \chi^{-1}\!\left(
   \prod_\eta\Delta_\eta(\mu_{L_\eta})
 \right).                                          \tag{7.10}
\]

This is exactly the simultaneous character rule implemented in Section 5.
It is not merely a stabilizer calculation for a candidate chart: (7.10)
is the stabilizer of the selected factor in the local quotient presentation
of \(\mathcal G\).

Combining Theorem 7.2 and Sections 7.3--7.5 with the source-graph
finite-normalization theorem gives:

> **H1-STACK theorem.**
> The polynomial Hurwitz--LL graph is the normalized selected
> admissible-cover stack over the normalized wonderful pullback of the
> labelled branch-diagonal building set.  Its coarse space is that wonderful
> graph.  Its explicit logarithmic charts are (0.1)--(0.4), and its complete
> inertia is the simultaneous full-centralizer character of Section 5.

The computations below verify coordinate identities and sample the
character compiler; they are not used for exhaustive coverage or descent.
Those follow from Theorem 7.2 at arbitrary geometric points.  The coarse
comparison is also the independent `H1-COARSE` theorem.  Neither theorem
revives the original map from the unmodified root-stable space.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_recursive_resonance_atlas.py
```

The checker verifies:

- all \(534\) branch-only nested screen families with two through five
  moving labels and all \(534\) independent affine gauge changes;
- \(1{,}453\) one-screen and \(2{,}926\) two-screen contractions;
- \(84\) exact normalized initial forms, flag equations, and transitivity
  identities through degree seven;
- \(84\) source/target affine-frame transitions for the maps and flag
  divisors;
- \(63\) nonfactorized smoothing families satisfying coefficientwise
  initial-form divisibility;
- \(801\) automatically extracted single-node characters on every reduced
  polynomial factorization through degree five;
- every ordered pair of tame cycle profiles through degree six;
- the closed radial inertia formula on \(76\) ordered charts from twelve
  equal- and unequal-multiplicity profiles of total degree at most seven;
- all \(89\) interval-nested charts of a degree-six matching tree;
- all thirteen degree-six polynomial radial types, each with trivial
  normalized inertia despite up to eight special-fiber automorphisms;
- an unrigidified tame chain with genuine order-six ghost inertia;
- the order-four pair--triple Maxwell inertia;
- concrete global and inertia subgroup closure on \(180\) compiled atlases;
- invariance under all \(720\) permutations of the six sheet labels.

The machine-readable certificate is
[`recursive_resonance_atlas.json`](../artifacts/generated-results/recursive_resonance_atlas.json).
