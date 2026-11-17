# Completed local rings at nested collisions

## Verdict

The ambient collision-tree local algebra and the local algebra of the
repository's corrected selected polynomial graph are now determined.
The stronger assertion for an arbitrary, unmodified root-parameter
degeneration is false: there one must also retain the **integral contact
monoid** of the map from the root-cluster chart to the target-screen chart.

The collision forest and the special-fiber vertex maps determine:

- the target and source dual graphs;
- the algebraic cover on every vertex;
- the expansion-index profile over every target edge; and
- the phase matching and inertia of a normalized lift.

The forest data do not, without a primitivity hypothesis, determine the
orders with which the target smoothing parameters and the successive Rees
generators vanish in the root parameters.  Those orders determine the
normalized blowups.  They therefore cannot be discarded from a theorem
about the completed local ring of a graph constructed directly over that
root chart.

For such a direct presentation, the corrected invariant is

\[
 \boxed{\text{collision forest}+
        \text{contact monoid}+
        \text{vertex covers}+
        \text{node-cycle character}.}
\]

The repository avoids this ambiguity in its corrected graph.  It first
forms the normalized wonderful target graph \(B\), then takes the reduced
generic polynomial closure inside the finite admissible-cover pullback over
\(B\), and finally normalizes.  At a point of this selected closure, the
standard admissible-cover deformation ring pulled back to
\(\widehat{\mathcal O}_{B}\) determines the completed local algebra.
Consequently the contact monoid is the characteristic monoid already
present on \(B\); it is not reconstructed from collision tables.

Here "expansion indices" means the full multiset
\((e_{\epsilon,j})_j\) of source-edge indices over every target edge
\(\epsilon\), not one integer per target vertex.

The equivalence relation also needs to be fixed.  In this note, two complete
local \(k\)-algebras are **stably strictly-etale equivalent** if, after finite
separable residue-field extensions, adjoining finitely many formal
power-series variables makes them isomorphic.  Regularity, the toric lattice
indices, and the characteristic monoid are invariant under this equivalence.

## 1. The counterexample to the uncorrected statement

The degree-five calculation gives the local rational map

\[
 [\lambda_1:\lambda_2]=[u_1x^3:u_2y^2].
\]

Its graph is the normalized blowup of \((x^3,y^2)\).  The associated
two-dimensional fan inserts the primitive ray

\[
 \rho=(2,3).
\]

The two affine toric charts have lattice indices \(3\) and \(2\).

Now make the tame ramified base change \(x=X^2\).  The limiting cover,
collision forest, source-node expansion indices, normalized vertex maps,
and node-cycle permutations are unchanged.  The graph is instead the
normalized blowup of

\[
 (X^6,y^2),
\]

whose primitive ray is \((1,3)\).  Its chart indices are \(3\) and \(1\).
In particular, the chart on the \(y\)-side is a cyclic quotient singularity
of index \(2\) before base change and is regular afterwards.  A smooth
factor and a finite etale extension cannot turn a singular completed local
ring into a regular one.

Thus the listed decorations do not determine the completed local ring for
arbitrary degenerations.  What changed is exactly the contact homomorphism

\[
 q_1\longmapsto x^3
 \quad\rightsquigarrow\quad
 q_1\longmapsto X^6.
\]

If the conjecture is restricted to the universal labelled root chart with
primitive root-separation parameters, this counterexample is excluded.
Then the contact monoid may be reconstructible from the forest and the
coefficientwise initial-form filtration, but that reconstruction is a
lemma that must be stated and proved; it is not a consequence of the
special-fiber initial maps alone.  This is Route B of Section 6, not an
obstruction to the corrected normalized-graph construction.

## 2. Collision-tree local-algebra theorem

The finite collision tables are instances of one edgewise
complete-intersection calculation.

Let \(\Gamma_{\rm tgt}\) be a stable target tree and, for every target edge
\(\epsilon\), let

\[
 \mathbf e_\epsilon=(e_{\epsilon,1},\ldots,e_{\epsilon,r_\epsilon})
\]

be the expansion indices of the source edges above \(\epsilon\).  Assume
that the residue characteristic is prime to every expansion index.  Fix the
vertex covers and their marked branch positions, and suppose that they are
etale points of their fully marked Hurwitz spaces.  Write \(z_1,\ldots,z_m\)
for the resulting unobstructed vertex-deformation parameters.

The completed ambient admissible-cover algebra is

\[
 A_\Gamma =
 k[[z_1,\ldots,z_m,\,
        (q_\epsilon)_\epsilon,\,
        (s_{\epsilon,j})_{\epsilon,j}]]
 \big/
 \left(s_{\epsilon,j}^{e_{\epsilon,j}}-q_\epsilon\right)_{\epsilon,j}.
                                                               \tag{2.1}
\]

For each edge put

\[
 L_\epsilon=\operatorname{lcm}_j(e_{\epsilon,j}),\qquad
 n_{\epsilon,j}=L_\epsilon/e_{\epsilon,j},\qquad
 \Phi_\epsilon=
 \left(\prod_j\mu_{e_{\epsilon,j}}\right)/\mu_{L_\epsilon}.
                                                               \tag{2.2}
\]

> **Collision-tree local-algebra theorem.**
> After a finite separable residue-field extension containing the required
> roots of unity:
>
> 1. \(A_\Gamma\) is a reduced finite flat complete intersection over
>    \(B=k[[z_1,\ldots,z_m,(q_\epsilon)_\epsilon]]\).
> 2. Its normalization is
>    \[
>      \widetilde A_\Gamma
>      \cong
>      \prod_{\beta\in\prod_\epsilon\Phi_\epsilon}
>      k[[z_1,\ldots,z_m,(t_\epsilon)_\epsilon]],               \tag{2.3}
>    \]
>    where on the factor labelled by
>    \(\beta=(\beta_\epsilon)_\epsilon\),
>    \[
>      q_\epsilon=t_\epsilon^{L_\epsilon},\qquad
>      s_{\epsilon,j}
>      =\omega_{\beta,\epsilon,j}
>       t_\epsilon^{n_{\epsilon,j}}.                            \tag{2.4}
>    \]
> 3. The number of normalization branches is
>    \[
>      M_\Gamma
>      =\prod_\epsilon
>       \frac{\prod_j e_{\epsilon,j}}{L_\epsilon}.              \tag{2.5}
>    \]
> 4. On every normalization branch,
>    \[
>    \begin{aligned}
>      \operatorname{Fitt}_0\Omega_{A_\Gamma/B}\,
>        \widetilde A_{\Gamma,\beta}
>        &=
>        \left(
>          \prod_\epsilon
>          t_\epsilon^{
>            \sum_j(e_{\epsilon,j}-1)n_{\epsilon,j}}
>        \right),\\
>      \operatorname{Fitt}_0
>        \Omega_{\widetilde A_{\Gamma,\beta}/B}
>        &=
>        \left(
>          \prod_\epsilon t_\epsilon^{L_\epsilon-1}
>        \right),\\
>      \mathfrak c_\Gamma\widetilde A_{\Gamma,\beta}
>        &=
>        \left(
>          \prod_\epsilon t_\epsilon^{c_\epsilon}
>        \right),
>    \end{aligned}                                               \tag{2.6}
>    \]
>    with
>    \[
>      c_\epsilon
>      =\sum_j(e_{\epsilon,j}-1)n_{\epsilon,j}
>       -(L_\epsilon-1)
>      =(r_\epsilon-1)L_\epsilon
>       -\sum_jn_{\epsilon,j}+1.                                \tag{2.7}
>    \]
> 5. The dual tree and vertex covers enter the algebra through the edge
>    profiles and the smooth variables.  The node-cycle character does not
>    add equations to (2.1): it selects the subgroup of the vertex
>    centralizer stabilizing a branch \(\beta\), hence the stack quotient
>    and its coarse invariant ring.

### Proof

For one edge, eliminate any variables of index one.  The remaining ring is
finite free over \(k[[q_\epsilon]]\), with basis given by the monomials
\(\prod_j s_{\epsilon,j}^{a_j}\), \(0\leq a_j<e_{\epsilon,j}\).
It is a complete intersection.  Tameness makes the polynomials
\(X^{e_{\epsilon,j}}-q_\epsilon\) separable over \(k((q_\epsilon))\);
hence the ring is generically reduced and, being finite flat over the
reduced base, is reduced.

Every common Puiseux parametrization has the form (2.4).  Its choices of
roots of unity are \(\prod_j\mu_{e_{\epsilon,j}}\), while changing
\(t_\epsilon\) by an \(L_\epsilon\)-th root of unity acts diagonally.
Thus the branches are indexed by \(\Phi_\epsilon\), of cardinality
\(\prod_j e_{\epsilon,j}/L_\epsilon\).  The integral closure of each branch
is \(k[[t_\epsilon]]\): the gcd of the integers
\(L_\epsilon,n_{\epsilon,1},\ldots,n_{\epsilon,r_\epsilon}\) is one.

Distinct target nodes have disjoint smoothing variables.  The global ring
(2.1) is their completed tensor product, followed by the formally smooth
factor \(k[[z_1,\ldots,z_m]]\).  Normalization branches and their
parametrizations therefore take Cartesian products.  This proves
(2.3)--(2.5).

The relative Jacobian matrix for \(A_\Gamma/B\) is diagonal by edge and
source node, with diagonal entries
\(e_{\epsilon,j}s_{\epsilon,j}^{e_{\epsilon,j}-1}\).
The integers \(e_{\epsilon,j}\) are units, so its determinant gives the
first line of (2.6).  On a normalized branch, the map to the smoothing base
is the product of the maps \(q_\epsilon=t_\epsilon^{L_\epsilon}\); its
Jacobian gives the second line.

Finally \(A_\Gamma/B\) is a finite flat Gorenstein complete intersection.
Complete-intersection duality compares its Kähler different with the
different of its normalization:

\[
 \mathfrak D_{A_\Gamma/B}\widetilde A_\Gamma
 =
 \mathfrak c_\Gamma\,
 \mathfrak D_{\widetilde A_\Gamma/B}.             \tag{2.8}
\]

Subtracting the normalized different exponent \(L_\epsilon-1\) from the
raw exponent on every independent edge gives (2.7), and multiplying the
edge contributions gives the conductor in (2.6).  The last assertion is
the definition of the stabilizer of a chosen normalization factor under
the vertex centralizer. \(\square\)

Thus every Maxwell, caustic, radial, or mixed table is a substitution into
(2.2)--(2.7), not a separate case proof.  Tree topology matters for matching
the vertex automorphisms and node cycles, but the ambient local algebra
factorizes over its edges.

This theorem does **not** prove the corresponding assertion for a selected
graph cut out inside a root-parameter space.  The latter also remembers how
its primitive parameters map to the \(q_\epsilon\)'s.  That additional
integral contact datum is the subject of the next theorem.

## 3. Saturated contact-monoid theorem

The following gives a non-tautological version of the viable theorem.

Let \(x\) be a tame nested collision point on a labelled root-cluster chart.
Assume:

1. the target screen model is toroidal at the image of \(x\);
2. in a regular system of primitive root-cluster parameters, every generator
   of every transformed branch-diagonal and coefficient-divisibility ideal
   is a monomial times a unit;
3. the exponent matrix of those monomials is constant on the chosen
   stratum; and
4. the chosen vertex cover is an etale point of the relevant fully marked
   Hurwitz space after its branch positions are fixed.

Let \(P_x\) be the unsaturated contact/Rees monoid.  It records the
columns of the exponent matrix of every target scale \(q_\epsilon\), every
dominant transform of a branch-diagonal ideal, and every
coefficient-divisibility generator in primitive root-cluster parameters.

For each target edge \(\epsilon\), introduce the Harris--Mumford node monoid

\[
 H_\epsilon=
 \langle\delta_\epsilon,\sigma_{\epsilon,j}\mid
   \delta_\epsilon=e_{\epsilon,j}\sigma_{\epsilon,j}\rangle.
\]

Let \(T=\mathbb N^{E(\Gamma_{\rm tgt})}\), with its generator
\(\delta_\epsilon\) mapping both to the contact element of \(P_x\) and to
the element \(\delta_\epsilon\) of \(H_\epsilon\).  Set

\[
 Q_x=
 P_x\mathop{\oplus}_{T}
 \bigoplus_\epsilon H_\epsilon.                    \tag{3.1}
\]

The torsion characters of \(Q_x^{\rm gp}\) index the Kummer phase branches.
For a selected phase branch \(\beta\), remove the corresponding torsion and
saturate, obtaining \(Q_{x,\beta}^{\rm sat}\).

> **Monomial-contact completed-local-ring theorem.**
> After strict henselization and up to a formal power-series factor, the
> normalized completed atlas ring of the selected admissible-cover graph at
> \((x,\beta)\) is
> \[
>   k'[[Q_{x,\beta}^{\rm sat}]].                    \tag{3.2}
> \]
> Its stack chart is the quotient of this formal toric chart by the subgroup
> of the full vertex centralizer whose node-cycle character preserves
> \(\beta\).  Its coarse completed ring is the corresponding invariant
> ring.

### Proof

After the finite etale extraction of roots of the units, the completed
multi-Rees chart in hypothesis 2 is the completion of the semigroup algebra
\(k[P_x]\).  Kummer base change at all source nodes is tensor product over
the target-scale algebra \(k[T]\); for monoid algebras this tensor product is
the monoid pushout (3.1).  The torsion of its groupification gives the
connected phase factors.  On a selected factor \(\beta\), normalization of
the affine semigroup algebra replaces the torsion-free monoid by its
saturation.  Completion commutes with this finite normalization in the
excellent local rings at hand.  Hypothesis 4 contributes only formal
power-series parameters.  Finally, a cover automorphism acts on a chosen
factor exactly when its node-cycle character fixes \(\beta\), giving the
stated stack quotient. \(\square\)

Formula (3.2) specializes correctly in both directions:

- for the ambient ACV space, \(P_x=T\), and each selected normalized branch
  is smooth;
- for a graph presented directly over a root chart, \(P_x\) is the monoid
  supplied by the normalized multi-Rees construction, so weighted blowups
  and resonance subdivisions remain visible.

> **Corrected-graph corollary.**
> For the repository's selected graph constructed by normalized closure
> over the normal wonderful target graph \(B\), the conclusion (3.2) holds
> with \(P_x\) equal to the characteristic monoid of \(B\).  No
> coefficient-divisibility monomialization hypothesis is required.

Indeed, pull the Harris--Mumford algebra back to
\(\widehat{\mathcal O}_{B,x}\).  Selecting the factors met by the generic
polynomial section and normalizing gives exactly (3.1)--(3.2).
Normalization commutes with etale localization and with completion for
these excellent finite algebras.  The
[finite-normalization theorem](SOURCE_GRAPH_FINITE_NORMALIZATION.md)
constructs the selected closure, while the
[recursive resonance atlas](RECURSIVE_RESONANCE_ATLAS.md) records the same
local node pushout.  Hence simultaneous monomialization of
coefficient-divisibility ideals is not needed to determine these completed
local rings. \(\square\)

The recursive screen formulas still give valuable explicit coordinates.
On the already constructed normalized chart, regularity of the universal
admissible cover implies coefficientwise divisibility; normalized initial
forms and their contraction identities are consequences of the chart,
rather than equations used to define it.

## 4. The four smallest tests

### 4.1 Degree five

For \((x^3,y^2)\), \(P_x\) inserts the primitive ray \((2,3)\).  The
normalized blowup has chart indices \((3,2)\).  This test detects the contact
monoid and already disproves any version retaining only the limiting cover.

The one-node monomial curve

\[
 k[[s_1,s_2]]/(s_1^3-s_2^2)
\]

has normalization \(k[[t]]\), with \(s_1=t^2,s_2=t^3\), and conductor
\((t^2)\).  This conductor is determined by the unsaturated semigroup
\(\langle2,3\rangle\), not by its saturation \(\mathbb N\).

### 4.2 Pairwise Maxwell collision

Disjoint transpositions give the node profile

\[
 (1,1,2,2).
\]

Here \(L=2\) and \(M=2\).  After eliminating the index-one variables, the
non-normal factor is

\[
 k[[a,b]]/(a^2-b^2).
\]

It has two normalization branches, conductor \((t)\) on each branch, and
normalized branch map \(q=t^2\).  The relative different of a normalized
branch over \(k[[q]]\) is \((t)\).  The full cycle character leaves the
known diagonal \(\mu_2\) inertia on the polynomial Maxwell lift.

### 4.3 Simple caustic collision

Adjacent transpositions give

\[
 (1,1,1,3).
\]

After eliminating the index-one variables, (2.1) is already \(k[[t]]\),
with \(q=t^3\).  There is one normalization branch, unit conductor, relative
different \((t^2)\), and no Maxwell sign choice.  Thus Maxwell and caustic
have visibly different unsaturated node algebras even before vertex
automorphisms are considered.

### 4.4 Mixed nested collision

For two target edges carrying the preceding Maxwell and caustic profiles,
the ambient completed factor is their completed tensor product.  It has two
normalization branches, with

\[
 q_M=t_M^2,\qquad q_C=t_C^3.
\]

Its conductor on either normalized branch is generated by \(t_M\), while
the normalized relative different over the two target smoothing parameters
is generated by \(t_Mt_C^2\).  In a genuinely nested cover, the inertia is
not obtained by multiplying two local answers blindly: the shared
connector permutations must satisfy the simultaneous character condition.
That is exactly the role of the repository's cycle-matching compiler.

## 5. Conductor, Fitting ideals, and nilpotency

These three consequences require different inputs.

### Conductor

The conductor belongs to a finite map \(A\to\widetilde A\), not to the
normal ring \(\widetilde A\) alone.  In the monoidal model it is determined
by the **unsaturated** monoid \(Q_x\), its phase gluing, and its saturation.
On a domain branch its monomial exponents are

\[
 \mathfrak c(Q_x)=
 \{m\in Q_x^{\rm sat}:m+Q_x^{\rm sat}\subseteq Q_x\}.
                                                               \tag{5.1}
\]

For several branches one must additionally impose the phase-gluing
conditions.  This is finite affine-semigroup computation, but it cannot be
recovered from \(Q_x^{\rm sat}\) alone.

For the raw admissible-cover node ring (2.1), there is a closed formula.
Put

\[
 n_j=\frac{L}{e_j}.
\]

Then on **every** normalization branch the conductor is

\[
 \boxed{\;
 \mathfrak c_{\mathbf e}
   =(t^{c(\mathbf e)}),\qquad
 c(\mathbf e)
   =\sum_j(e_j-1)n_j-(L-1)
   =(r-1)L-\sum_j n_j+1.
 \;}                                                \tag{5.2}
\]

Indeed, \(A_{\mathbf e}\) is a finite flat complete intersection over
\(k[[q]]\), hence Gorenstein.  Its relative Kähler different pulls back to
order

\[
 d_{\rm raw}=\sum_j(e_j-1)n_j,
\]

while the different of the normalized branch
\(k[[q]]\subset k[[t]]\), \(q=t^L\), has order \(L-1\).  Complete-intersection
duality gives

\[
 \mathfrak D_{A_{\mathbf e}/k[[q]]}\widetilde A
 =
 \mathfrak c_{\mathbf e}\,
 \mathfrak D_{\widetilde A/k[[q]]},
\]

which proves (5.2).  Identity strands \(e_j=1\) cancel automatically from
the second expression.  For several independent target nodes, the conductor
on a normalization branch is the single monomial

\[
 \prod_\epsilon t_\epsilon^{c(\mathbf e_\epsilon)}. \tag{5.3}
\]

### Fitting ideals

One must name the differential module.  For the raw node ring (2.1),
tameness gives

\[
 \operatorname{Fitt}_0\Omega_{A_{\mathbf e}/k[[q]]}
 =
 \left(\prod_j s_j^{e_j-1}\right).                 \tag{5.4}
\]

On a normalized branch its pullback has \(t\)-order

\[
 \sum_j(e_j-1)\frac{L}{e_j}.                       \tag{5.5}
\]

For the normalized branch itself,
\(\operatorname{Fitt}_0\Omega_{k[[t]]/k[[q]]}=(t^{L-1})\).
Consequently the conductor identity can be read numerically as

\[
 \boxed{c(\mathbf e)=
   \operatorname{ord}_t\operatorname{Fitt}_0
       \Omega_{A_{\mathbf e}/k[[q]]}
   -
   \operatorname{ord}_t\operatorname{Fitt}_0
       \Omega_{k[[t]]/k[[q]]}.}                    \tag{5.6}
\]

For a selected graph, the remaining logarithmic Jacobian is determined by
the contact-monoid map; ordinary vertex-map Fitting factors must be retained
separately if they are part of the intended invariant.

### Collision nilpotency

A normal completed local ring has no nilpotents.  "Collision nilpotency"
therefore refers to a specified contraction or special fiber.  For the
coarse selected-root algebra of a block of \(\mu\) colliding roots, the
special fiber is

\[
 k[T]/(T^\mu),
\]

so its nilpotency index is \(\mu\).  More general nilpotency indices require
the quotient ideal defining the contraction.  They do not follow from the
normalized local ring without that extra morphism.

## 6. Two construction routes and present status

There are two logically different routes.

### Route A: corrected normalized graph

This route is complete for local algebras.

1. Construct the normal wonderful target graph \(B\).
2. Pull the finite admissible-cover stack back to \(B\).
3. Take the reduced closure of the generic polynomial section and normalize.
4. Etale locally, use the standard Harris--Mumford deformation algebra.
5. Form the node pushout (3.1), select its phase factors, and saturate.
6. Apply the simultaneous vertex-centralizer character to recover inertia.
7. Compute conductors and differents from (2.6)--(2.7), and treat contraction
   nilpotency only after specifying its quotient algebra.

This proves the collision-tree local-algebra theorem at every point of the
corrected selected graph.  It does not by itself provide one global explicit
coordinate atlas or prove effective descent of independently written screen
charts; those are separate `H1-STACK` tasks.

### Route B: direct equations in primitive root parameters

To recover the same normalization without first using the finite
admissible-cover pullback, one would:

1. define \(P_x\) as the value monoid of the multi-Rees algebra of all
   transformed branch-diagonal and coefficient-divisibility ideals;
2. prove simultaneous monomial contact after removing smooth parameters;
3. apply the saturated contact-monoid theorem.

This remains an optional presentation problem.  Its failure would mean that
two families with the same proposed contact monoid and vertex covers retain
different non-toric formal relations.  Such a failure would obstruct that
direct presentation, but not the completed-local-ring classification of
Route A.

## 7. Relation to the literature

Abramovich--Corti--Vistoli
([Twisted bundles and admissible covers](https://arxiv.org/abs/math/0106211))
identify the twisted-cover stack with the normalization of the
Harris--Mumford admissible-cover space.  The explicit local calculation in
Cavalieri--Markwig--Ranganathan
([Tropicalizing the space of admissible covers](https://arxiv.org/abs/1401.4626))
gives (2.1)--(2.3),
shows that the Harris--Mumford space is locally a possibly nonnormal toric
space, and that its normalization is smooth with toroidal boundary.  Their
warning is also essential: a tropical combinatorial type may represent
several algebraic strata, so the tropical cone does not select the
polynomial component.

Li's
[wonderful compactification theorem](https://arxiv.org/abs/math/0611412)
explains the target screens and their nested sets.  Affine-semigroup
normalization explains why the integral contact
monoid, rather than the real cone or forest alone, is the correct invariant:
normalization replaces a monoid by its saturation, while the conductor
remembers the inclusion of the original monoid.

The conceptual conclusion is therefore:

\[
\boxed{\text{collision tree}+\text{wonderful contact monoid}
       +\text{vertex cover}+\text{node-cycle character}
       \Longrightarrow\text{completed local algebra}.}
\]

Simultaneous monomial contact for the universal primitive root chart would
give a shorter direct equation-level presentation.  It is no longer needed
for the completed-local-ring classification itself.  If it fails, its
failure ideal is the extra decoration required only by that direct
presentation.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_combinatorial_completed_local_rings.py
```

The checker compiles the degree-five monomial blowups, the pairwise and
triple Maxwell node rings, the caustic node ring, and their mixed completed
tensor product.  It also checks the phase-quotient branch count and the
closed conductor/Fitting/different formula for all \(1{,}554\) ordered
profiles of width at most four with indices at most six.  An independent
phase-character calculation recovers the conductor for all \(155\) profiles
of width at most three with indices at most five.  The reusable compiler is
[`collision_local_rings.py`](../jcsearch/collision_local_rings.py).
