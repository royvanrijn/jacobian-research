# Programme 4: the universal decomposition category

> **Status.** This is a research programme, not a classification theorem.
> It records two formal factorization consequences of the existing
> Keller results, gives a precise category-level formulation, and separates
> five open recognition and uniqueness problems.  In particular, it does
> **not** assert that every atomic Keller map is a marked-root suspension,
> that monodromy alone detects polynomial factorization, or that a
> Jordan--Hölder theorem currently holds.

Work over a characteristic-zero field \(k\).  Write

\[
 \mathcal K_d(k)=
 \{F:\mathbb A_k^d\longrightarrow\mathbb A_k^d:
   \det DF\in k^\times\}.
\]

The existing
[atomic-spectrum theorem](../verified/ETALE_MONOID_ATOMIC_SPECTRUM.md)
treats \(\mathcal K_d(k)\) as a monoid.  That is enough to define atoms and
to construct words of all prescribed degree patterns, but it forgets:

1. the intermediate affine spaces in a factorization;
2. stable polynomial left--right identifications between maps;
3. identifications of their finite marked fibers;
4. the block flags carried by a marked sheet of a composite; and
5. boundary, conductor, and Fitting data which can distinguish stable
   left--right classes.

Programme 4 keeps these data.  The correct first object is not an ordinary
quotient monoid.  Stable left--right equivalence is two-sided and therefore
does not descend naively through composition.  It belongs as a square, or
equivalently as an invertible \(2\)-morphism, in a decorated factorization
category.

## 1. The fixed-map decomposition category

Let

\[
 A_F=k[F_1,\ldots,F_d]\subseteq B=k[x_1,\ldots,x_d].
 \tag{1.1}
\]

Define the **polynomial-sandwich poset**

\[
 \operatorname{Sand}(F)=
 \left\{
 R:
 A_F\subseteq R\subseteq B,\quad
 R\simeq k[t_1,\ldots,t_d]
 \right\},
 \tag{1.2}
\]

ordered by inclusion.  The two endpoints are \(A_F\) and \(B\).

The
[polynomial-sandwich criterion](../verified/IMPRIMITIVE_KELLER_FACTORIZATION.md)
immediately gives the following dictionary.

> **Factorization--sandwich dictionary.**  Two factorizations
> \[
> F=G\circ H
> \]
> which differ only by a polynomial automorphism of the intermediate
> \(\mathbb A^d\) determine the same object
> \(R=k[H_1,\ldots,H_d]\) of \(\operatorname{Sand}(F)\), and every object of
> \(\operatorname{Sand}(F)\) arises this way.

Thus \(\operatorname{Sand}(F)\), rather than the full intermediate-field
lattice, is the intrinsic \(1\)-categorical factorization object.  A
factorization

\[
 F=F_r\circ\cdots\circ F_1
 \tag{1.3}
\]

gives a chain

\[
 A_F=R_0\subsetneq R_1\subsetneq\cdots
 \subsetneq R_r=B.                                  \tag{1.4}
\]

The chain is maximal exactly when every \(F_i\) is atomic.  Notice that
"maximal" is taken inside the polynomial-sandwich poset, not inside the
larger lattice of all intermediate fields.

This distinction is essential.  If \(M/k(F)\) is the Galois closure,
\(G=\operatorname{Gal}(M/k(F))\), and \(P=\operatorname{Gal}(M/k(x))\),
then generic intermediate fields form the subgroup interval

\[
 [P,G]=\{J:P\subseteq J\subseteq G\}.                \tag{1.5}
\]

Only the subgroups whose fields admit a compatible polynomial affine-space
model occur in (1.2).

## 2. A finite atomic factorization theorem

The following is a formal consequence of results already proved in the
repository.

> **Proposition 2.1.** Every noninvertible characteristic-zero Keller map
> over \(k\) admits a finite factorization into Keller maps which are atomic
> over \(k\).

**Proof.**  If \(F\) is not atomic, write \(F=G\circ H\) with both factors
noninvertible.  The chain-rule argument in
[primitive-monodromy atomicity](../verified/PRIMITIVE_MONODROMY_ATOMICITY.md)
shows that \(G\) and \(H\) are Keller maps.  Geometric degree is
multiplicative:

\[
 \operatorname{gdeg}F
 =\operatorname{gdeg}G\,\operatorname{gdeg}H.
 \tag{2.1}
\]

Both factor degrees are strictly smaller than \(\operatorname{gdeg}F\).
Induction on this positive integer terminates. \(\square\)

After base change to an algebraic closure, the same induction gives a
factorization into geometrically atomic factors there.  Descent of those
geometric factors to \(k\) is a separate question.  The proposition also
does not say that the \(k\)-atomic factors are absolutely or stably atomic.

The genuinely new generation question is therefore not mere existence of
atoms.  It is:

> **P4-A (marked-suspension generation).**  Is every geometrically atomic
> Keller map stably polynomially left--right equivalent to a
> geometrically atomic marked-root suspension?  If not, what additional
> atomic generators are required?

Here "marked-root suspension" must mean an intrinsic incidence
construction, not the existence of a primitive inverse equation after
choosing coordinates.  The universal maps
\(\mathcal U_N\), root-engineered quadratic gauges, weighted maps, and
fiber-invisible or power-shifted gauges provide all-degree test generators.
They do not constitute a completeness theorem.

## 3. The decorated Keller double category

The minimal global structure should be a double category
\(\mathbf{Kell}^{\mathrm{dec}}_{\mathrm{st}}\).

- Objects are affine spaces, with identity stabilization formally inverted.
- Horizontal arrows are Keller maps equipped with the intrinsic
  normalization package
  \[
  \mathcal B(F)=
  (\overline X_F\to Y,\,
    \mathbb A^d\hookrightarrow\overline X_F,\,
    \partial_F)
  \tag{3.1}
  \]
  and only those boundary or fiber marks which have been reconstructed
  intrinsically.
- Vertical arrows are stable polynomial automorphisms of source or target.
- A square is a stable left--right commutative square together with the
  induced isomorphism of intrinsic finite-cover decorations.
- A marked-fiber \(2\)-cell identifies finite étale fibers, or their
  incidence covers, compatibly with the selected affine sheet and all
  retained labels.

Horizontal composition is ordinary composition of maps.  A marked sheet of
a composite

\[
 X_0\xrightarrow{F_1}X_1\xrightarrow{F_2}\cdots
 \xrightarrow{F_r}X_r
 \tag{3.2}
\]

records the whole flag

\[
 x_0\longmapsto x_1\longmapsto\cdots\longmapsto x_r.
 \tag{3.3}
\]

On a geometric generic fiber this is the nested block flag associated with
the field tower (1.4).  A single unflagged identification of the composite
fiber loses this structure.

The use of a double category avoids a false operation.  If
\[
 F'\simeq_{\mathrm{LR}}F,\qquad G'\simeq_{\mathrm{LR}}G,
 \]
there is no reason for the middle automorphisms in the two equivalences to
match.  Hence \(G'\circ F'\) need not represent a well-defined product of
the two stable LR classes.  Compatible squares do compose.

The
[stable-normalization theorem](../verified/STABLE_NORMALIZATION_FUNCTORIALITY.md)
proves that (3.1), labeled valuations, intersections, Fitting ideals,
conductors, and intrinsic finite-stratum maps are transported by these
squares and simply acquire affine factors under stabilization.  The
[decorated-normalization invariant](DECORATED_NORMALIZATION_INVARIANT.md)
supplies the strongest current decoration on the weighted clean locus.

### 3.1 Marking discipline

A presentation root is not automatically a categorical mark.  A mark may be
retained only when it is:

1. the tautological sheet on a genuine marked incidence cover;
2. characterized uniquely from \(\mathcal B(F)\); or
3. supplied as explicit extra structure and required to be preserved by
   morphisms.

This is the same distinction already enforced by stable normalization
functoriality.  It prevents a chosen primitive element, target coordinate, or
root label from being promoted silently to a stable invariant.

## 4. What monodromy and normalization can recognize

Monodromy gives a finite list of candidate generic factorizations: the block
systems in (1.5).  It is decisive in the primitive case, where there are no
proper intermediate fields and the Keller map is atomic.  In the
imprimitive case it is only the first gate.

For \(J\in[P,G]\), let \(E=M^J\) and let

\[
 Z_E=\operatorname{Norm}_{\mathbb A^d_{\mathrm{target}}}(E).
 \tag{4.1}
\]

The polynomial-sandwich theorem gives an exact recognition criterion:

\[
\boxed{
\begin{array}{c}
E\text{ comes from a polynomial factorization}\\[2pt]
\Longleftrightarrow\\[2pt]
Z_E\text{ contains a compatible open }U_E\simeq\mathbb A^d
\text{ into which the source lift is regular.}
\end{array}}
\tag{4.2}
\]

For a chain of block systems, the affine reconstruction opens must be nested
compatibly.  Consequently:

- monodromy alone recognizes candidate block towers, not polynomial
  factorizations;
- monodromy plus coarse ramification or boundary graphs still need not
  recognize the affine reconstruction open;
- monodromy plus the full intermediate normalizations, their source lifts,
  and compatible affine opens recognizes the factorization by (4.2).

This suggests the executable problem.

> **P4-B (decorated algebraization).**  Construct an intrinsic finite test
> for the existence and uniqueness of the affine reconstruction open in
> (4.2), using valuation signs, conductor/Fitting data, and the marked
> affine sheet.

A positive answer would turn the group-theoretic block lattice into the
actual category \(\operatorname{Sand}(F)\).  A counterexample should be an
imprimitive Keller map for which an intermediate normalization is rational
but has no compatible polynomial sandwich.

### 4.1 Full-wreath order rigidity

The first composite laboratory already gives a nontrivial conclusion.

> **Proposition 4.1.**  Suppose the geometric monodromy of a degree-\(ab\)
> cover is the full wreath product \(S_a\wr S_b\) in its natural
> imprimitive action, with \(a,b\ge2\).  Its only proper nontrivial block
> system is the canonical partition into \(b\) blocks of size \(a\).
> Consequently every nontrivial polynomial factorization of the associated
> Keller map has inner degree \(a\) and outer degree \(b\).

**Proof.**  Let \(D\) be a block containing a point \(\omega\) in the
canonical block \(C_1\).  A base-group permutation supported on another
canonical block fixes \(\omega\), hence preserves \(D\).  The stabilizer of
\(C_1\) permutes the other canonical blocks transitively.  Therefore \(D\)
either meets none of them or contains all of them.  In the latter case, a
base permutation in \(C_1\) moving \(\omega\) produces a translate of \(D\)
which intersects \(D\) outside \(C_1\); the block axiom then forces equality
and hence \(D\) contains all of \(C_1\).  Thus \(D\) is the whole fiber.
In the former case, \(D\subseteq C_1\), and primitivity of the natural
\(S_a\)-action makes \(D\) either \(\{\omega\}\) or \(C_1\).

The factor-degree conclusion follows from the block-system/intermediate-field
dictionary and the polynomial-sandwich criterion. \(\square\)

The
[degree-twelve wreath calculation](../verified/IMPRIMITIVE_KELLER_FACTORIZATION.md)
proves

\[
 \operatorname{Mon}_{\mathrm{geom}}(C_{3,4})=S_3\wr S_4.
 \tag{4.3}
\]

Hence every nontrivial factorization of \(C_{3,4}\) has inner degree \(3\)
and outer degree \(4\).  In particular, \(C_{3,4}\) has no reverse
\(4\)-then-\(3\) factorization.  This proves ordered degree rigidity, not
uniqueness of the polynomial subalgebra with the unique intermediate
fraction field.

## 5. Boundary-decorated uniqueness

For an atomic factor \(A\), write

\[
 [A]_{\partial}
 =
 \bigl[
 A;\mathcal B(A),\operatorname{Fitt}_0\Omega,
 \text{conductor},\text{inertia},\text{intrinsic marks}
 \bigr]_{\mathrm{stable\ LR}}.                       \tag{5.1}
\]

A **decorated composition series** of \(F\) is a maximal chain (1.4)
together with:

1. the stable decorated class (5.1) of every adjacent atomic quotient;
2. the marked-sheet block flag of the composite generic fiber; and
3. specialization maps between the factor boundaries and the boundary of
   the composite.

The first uniqueness target should be local rather than global.

> **P4-C (diamond problem).**  Suppose two polynomial sandwiches
> \(R,S\in\operatorname{Sand}(F)\) are incomparable and both cover their
> intersection in the sandwich poset.  Determine when their join is again a
> polynomial sandwich and classify the resulting decorated diamond.

If every finite interval of \(\operatorname{Sand}(F)\) were lower
semimodular and every diamond were generated by a finite list of decorated
Ritt moves, maximal chains would have a common length and would be connected
by local replacements.  Neither assertion is currently proved for Keller
maps.

Degree data alone cannot supply this theorem.  The first total degree with
two different unordered nontrivial atom-degree patterns is \(24\):

\[
 24=3\cdot8=4\cdot6.
 \tag{5.2}
\]

The first total degree allowing different nontrivial atom-word lengths is
\(27\):

\[
 27=3\cdot9=3\cdot3\cdot3.
 \tag{5.3}
\]

These are only arithmetic possibilities among explicit words of known
atoms; they do not exhibit one Keller map with two such decompositions.
They show exactly where a degree-only Jordan--Hölder argument becomes
impossible.

## 6. Ritt moves as stable LR \(2\)-geometry

There are two Ritt structures in play and they must not be conflated.

1. A **Keller Ritt move** replaces two adjacent atomic factors by two others:
   \[
   A_2\circ A_1=B_2\circ B_1,
   \tag{6.1}
   \]
   together with a compatible stable LR square on the common composite and
   a marked-fiber block identification.
2. A **vertical seed Ritt move** is a second decomposition of
   \(H-sW\) inside one marked-root family.  The existing Hessian--Ritt
   programme studies these seed incidences and their stable boundary
   geometry.  It does not yet prove (6.1) for arbitrary Keller-map factors.

The
[Ritt move \(2\)-complex](RITT_MOVE_2_COMPLEX.md)
already gives the correct local model for coherence:

- vertices are complete normalized decompositions;
- edges are elementary adjacent Ritt correspondences;
- commuting squares and braid hexagons compare paths; and
- a \(2\)-cell identifies reduced images or normalizations, while retaining
  any nilpotent or derived defect of the full path schemes.

Programme 4 asks for a functor from the Keller factorization groupoid to this
coefficient-decorated \(2\)-complex on the marked-root subcategory.  On an
edge, the functor should transport:

1. the generic block flag;
2. the intermediate affine reconstruction open;
3. the decorated normalizations of all four adjacent atoms; and
4. the conductor, Fitting, and nilpotent comparison data on the common
   boundary.

The degree-thirty braid calculation warns against declaring two move paths
equal merely because they have the same reduction and normalization.  Full
stable LR geometry should therefore be correspondence-valued or derived,
not a bare graph quotient.

> **P4-D (Ritt realization).**  Determine which vertical seed Ritt edges
> lift to Keller Ritt \(2\)-cells, and compute the obstruction when the
> required intermediate affine reconstruction opens fail to glue.

The first useful tests are degrees \(12\), \(24\), \(27\), and \(30\).
Degree \(12\) has the two orders \(3\cdot4\) and \(4\cdot3\); degree \(24\)
is the first arithmetic competition between distinct two-factor multisets;
degree \(27\) is the first possible length competition; and degree \(30\)
already has a nontrivial braid coherence defect in the seed Ritt complex.

## 7. A Jordan--Hölder target

The direct analogue of Ritt's first theorem would be:

> **Decorated Jordan--Hölder conjecture.**  On a specified boundary-clean
> marked-root subcategory, any two decorated composition series of a Keller
> map:
>
> 1. have the same length;
> 2. have the same multiset of stable decorated atomic factors; and
> 3. are connected by finitely many decorated Keller Ritt moves whose
>    square and braid relations hold at the normalization level, with the
>    residual scheme structure recorded by coherent \(2\)-cells.

This is intentionally a subcategory conjecture.  For a general transitive
permutation group, maximal chains in the interval \([P,G]\) need not satisfy
a Jordan--Hölder theorem.  Classical polynomial decomposition is special:
the inertia at infinity supplies a long cycle, and the resulting
imprimitivity lattice is sufficiently semimodular for Ritt-type uniqueness.
The multivariable Keller setting has no established replacement for that
single boundary cycle.

The plausible substitute is the full decorated boundary inertia system:
the collection of tame inertia generators, their incidence graph,
normalization conductors, and the affine-sheet selector.  This leads to the
main structural question.

> **P4-E (boundary semimodularity).**  Find intrinsic hypotheses on the
> decorated boundary inertia system which force the polynomial-sandwich
> poset to be graded or lower semimodular.

A successful theorem would yield a Jordan--Dedekind length statement first.
Only after classifying decorated diamonds should one expect uniqueness of
the factor multiset or connectivity by Ritt moves.

## 8. Concrete first calculations

The programme should begin with finite, falsifiable calculations.

### 8.1 Composite laboratories

For
\[
 C_{a,b}=F_b\circ F_a,
 \]
compute:

1. the full monodromy group, not only its embedding in
   \(S_a\wr S_b\);
2. every subgroup between the point stabilizer and the monodromy group;
3. every compatible polynomial sandwich;
4. whether the displayed block system is unique; and
5. the decorated boundary maps of the two factors inside the composite.

The case \(C_{3,4}\) already has full wreath monodromy, and Proposition 4.1
settles its ordered factor degrees.  The remaining degree-twelve problem is
uniqueness of the polynomial affine model over that intermediate field.
The next cases should be \(C_{3,5}\), \(C_{4,3}\), and the competing
degree-\(24\) pairs.

### 8.2 Stable LR versus factor order

Compare
\[
 C_{a,b}=F_b\circ F_a,\qquad C_{b,a}=F_a\circ F_b.
 \tag{8.1}
\]

Their natural block sizes are different.  Test whether the full decorated
normalization of the composite reconstructs that ordered block flag and
therefore separates the two maps stably.  This is the smallest direct test
of whether factor order is visible in stable LR geometry.

### 8.3 The first diamond

Degree \(12\) is the smallest possible adjacent interchange.  Search for a
single Keller map having both a \(3\)-then-\(4\) and a \(4\)-then-\(3\)
polynomial sandwich.  Proposition 4.1 excludes the existing
\(C_{3,4}\) full-wreath composite, so a positive example must lie on a
proper-monodromy locus.  Such an example would produce the first Keller
Ritt square.

### 8.4 Fiber recognition

For a target with a complete squarefree fiber, retain the finite étale
algebra together with its block partitions under every candidate
intermediate field.  Compare this finite datum with the normalization-open
test (4.2).  The goal is to find the smallest example where:

- the fiber partition exists but does not globalize to a polynomial
  sandwich; or
- the marked fiber plus boundary decoration uniquely determines the
  sandwich.

## 9. Decision tree

For a supplied Keller map \(F\), the current rigorous workflow is:

1. compute geometric degree and generic monodromy;
2. if monodromy is primitive, certify atomicity;
3. otherwise enumerate block systems and intermediate fields;
4. normalize the target in each intermediate field;
5. test rationality and the compatible affine reconstruction open;
6. assemble all successful opens into \(\operatorname{Sand}(F)\);
7. decorate its maximal chains by intrinsic normalization and marked-fiber
   data; and
8. test whether different chains are connected by decorated diamonds.

Steps 1--4 are generic-cover problems.  Step 5 is the polynomial
algebraization gate.  Steps 6--8 are Programme 4.

## 10. External structural anchors

The classical comparison uses only the following general inputs.

- Finite étale covers form a Galois category, so connected generic covers
  are controlled by finite transitive monodromy sets and their block
  systems:
  [Stacks Project, Sections 58.5 and 58.7](https://stacks.math.columbia.edu/tag/03SF).
- Muzychuk--Pakovich prove a Jordan--Hölder theorem for imprimitivity
  systems when the monodromy contains a cyclic subgroup with at most two
  orbits; this includes the Ritt-type equality of maximal-decomposition
  lengths and degree multisets:
  [arXiv:0712.3869](https://arxiv.org/abs/0712.3869).
- Ritt transformations connect complete decompositions in the classical
  polynomial setting, while analogous statements fail in more general
  composition monoids without additional hypotheses:
  [Bavula, arXiv:0711.0913](https://arxiv.org/abs/0711.0913).

These results motivate P4-E.  They do not imply semimodularity of
\(\operatorname{Sand}(F)\) for multivariable Keller maps.

## 11. Current answers to the five questions

1. **Is every Keller map a composition of geometrically atomic marked-root
   suspensions?**  Every geometric Keller map has a finite atomic
   factorization, but completeness of marked-root suspension atoms and
   descent of the factors are open.
2. **Is factorization unique after boundary decoration?**  Unknown.  The
   first necessary problem is classification of decorated diamonds in the
   polynomial-sandwich poset.  Full-wreath monodromy already makes the
   ordered factor degrees of \(C_{3,4}\) rigid, but not its affine
   polynomial model.
3. **Does monodromy plus normalization recognize factorization?**  Monodromy
   gives candidate fields.  Full intermediate normalization plus a
   compatible affine reconstruction open recognizes polynomial factors.
   Coarse normalization data alone is not known to suffice.
4. **How are Ritt moves reflected in stable LR geometry?**  They should be
   invertible, correspondence-valued \(2\)-cells between adjacent factor
   words, preserving marked block flags and boundary decorations.  The
   existing Hessian--Ritt \(2\)-complex is the local seed model, not yet the
   global Keller factorization category.
5. **What is the Jordan--Hölder analogue?**  First prove that the decorated
   sandwich poset is graded under explicit boundary-inertia hypotheses;
   then classify diamonds and prove maximal-chain connectivity.  A global
   theorem for all Keller maps is presently unsupported.

## 12. Exact finite regression

Run

```bash
python3 scripts/verify_universal_decomposition_category.py
```

The dependency-free checker enumerates ordered degree words with factors at
least three through degree \(200\).  It verifies the four arithmetic
thresholds \(9,12,24,27\) used above.  It deliberately makes no claim that
two words occur as factorizations of the same Keller map.
