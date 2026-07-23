# Completed local rings at nested collisions

## Verdict

The proposed statement has a correct toroidal core, but it is false as
written for arbitrary tame degenerations.  The missing datum is not an
arbitrary analytic modulus.  It is the **integral contact monoid** of the
map from the root-cluster chart to the target-screen chart.

The collision forest and the special-fiber vertex maps determine:

- the target and source dual graphs;
- the algebraic cover on every vertex;
- the expansion-index profile over every target edge; and
- the phase matching and inertia of a normalized lift.

They do not, without a primitivity hypothesis, determine the orders with
which the target smoothing parameters and the successive Rees generators
vanish in the root parameters.  Those orders determine the normalized
blowups.  They therefore cannot be discarded from a theorem about the
completed local ring of the **selected polynomial graph**.

The corrected invariant is

\[
 \boxed{\text{collision forest}+
        \text{contact monoid}+
        \text{vertex covers}+
        \text{node-cycle character}.}
\]

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
special-fiber initial maps alone.

## 2. The exact ambient admissible-cover ring

At one target node with source-node indices \(e_1,\ldots,e_r\), the
Harris--Mumford completed factor is

\[
 A_{\mathbf e}
 =
 k[[q,s_1,\ldots,s_r]]/
 (s_1^{e_1}-q,\ldots,s_r^{e_r}-q).                 \tag{2.1}
\]

Put \(L=\operatorname{lcm}(e_1,\ldots,e_r)\).  After adjoining the required
roots of unity, a normalization branch has

\[
 q=t^L,\qquad s_j=\omega_jt^{L/e_j}.               \tag{2.2}
\]

The branches form

\[
 \left(\prod_j\mu_{e_j}\right)/\mu_L
\]

and hence have number

\[
 M_{\mathbf e}=\frac{\prod_j e_j}{L}.              \tag{2.3}
\]

For several target nodes these factors take completed tensor products.
After choosing one phase branch, the normalized ambient admissible-cover
ring is a power-series ring in one \(t_\epsilon\) per target edge, together
with the unobstructed vertex-deformation parameters.  This is the precise
sense in which the ambient ACV normalization is combinatorial.  The
matching-permutation character determines the stabilizer of the stack
chart, not additional equations in its smooth atlas.

This ambient statement does not prove the corresponding assertion for a
selected graph cut out inside a root-parameter space.  The latter also
remembers how its parameters map to the \(q_\epsilon\)'s.

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
- for the polynomial graph, \(P_x\) is the monoid supplied by the normalized
  multi-Rees construction, so weighted blowups and resonance subdivisions
  remain visible.

The repository's recursive screen formulas give the generators expected in
\(P_x\), and initial-form transitivity gives functoriality under contraction.
The remaining repository-specific problem is now exactly hypothesis 2:
prove simultaneous monomialization of all coefficientwise Rees-divisibility
ideals in primitive root parameters.  Equivalently, show that after removing
smooth vertex parameters the divisibility algebra has no non-toric formal
relations.

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

## 6. Proof program

The induction on nested-set depth can now be made precise.

1. Use tame Weierstrass preparation to split the source points and identify
   the smooth vertex-deformation factors.
2. Define \(P_x\) as the value monoid of the multi-Rees algebra of all
   dominant transforms of branch-diagonal and coefficient-divisibility
   ideals.
3. Prove a **simultaneous monomial-contact lemma**: after removing smooth
   parameters, all transformed coefficientwise-divisibility ideals are
   generated by monomials times units.  The monomial-contact theorem then
   supplies normalization by saturation.
4. Form the node pushout (3.1), separate its torsion characters, and
   saturate.  This recovers the Harris--Mumford/ACV branch count.
5. Use initial-form transitivity to show that the monoids and their phase
   branches commute with contraction of an intermediate screen.
6. Apply the full vertex-centralizer character to recover inertia and coarse
   invariant rings.
7. Compute node conductors by (5.2), general Rees conductors by (5.1),
   Fitting ideals from the monoid/Jacobian presentation, and nilpotency only
   after adjoining the specified contraction algebra.

The first likely genuine failure of this corrected theorem is no longer a
cross-ratio omitted from radial data: the full vertex initial maps already
retain that cross-ratio.  It would be a failure of simultaneous
monomialization, namely two selected families with the same contact monoid
and the same vertex covers but different non-toric formal relations.

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

The resulting breakthrough is therefore a sharpened target:

\[
\boxed{\text{prove simultaneous monomial contact for the universal
primitive root chart}.}
\]

That lemma would turn the repository's recursive atlas from a set of
compatible formulas into the desired completed-local-ring classification.
If it fails, its failure ideal is the precise "resonance decoration" that
must be added.

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
