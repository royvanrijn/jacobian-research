# Multiple protected rank-one obstructions and the coloring budget

The local four-strand obstruction can be packed on disjoint sets of
top strands. The stronger [two-sided theorem](TWO_SIDED_PROTECTED_OBSTRUCTIONS.md)
now proves that `k` protected blocks force BOTH a root defect and a
covector defect of at least `k`, so `ar_C<=b-k`. The weaker
`b-ceil(k/2)` argument below is retained as a consequence of the
original local disjunction. Explicit protected knot closures are supplied;
a matching meridional lower bound and prime gap-two realization remain open.

There is a separate obstruction to double-transposition colorings.
Repeating the displayed alternating eight-symbol seed loses one distinct
transposition edge per block. A complete local tree-component census below
proves the same support budget for EVERY coloring of disjoint protected
blocks by double transpositions: `v<=2b-k+1`. Together with the weaker
bound this certified at most a gap of one; the stronger two-sided bound
leaves room for a gap of `ceil(k/2)` if compatible colorings and knot
closures can be constructed.

## 1. The local hypothesis really needs only two protected meridians

Use the left-to-right Hurwitz convention
`T(a,c)=(c,c^-1*a*c)` and put

```
U=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3].
```

The signed trace in [the released proof](../../../results/prime-knot-augmentation-gaps/representations.tex),
Lemma `lem:spans`, has input arcs `(1,0,8,13)` and output arcs
`(1,6,16,13)`. The repeated first and fourth labels come from their
closure relations. That argument proves the following LOCAL statement.

**Local lemma (weaker consequence of the two-sided theorem).** Suppose the matrices propagated through one copy of
`U` have nonzero displacement rank one and common determinant
`mu!=0`, over any field and in any finite ambient dimension. If its
first and fourth output matrices equal the corresponding input
matrices, then its four input roots and four input covectors cannot
both be independent.

For clarity, this imports the seven-relation proof in the released
source, rather than an inference about a tail's total product. If both
four-element families were independent, replace the first seed by its
conjugate under the second, a triangular change with diagonal one on
both families. Write `D` for that second seed and `H,F` for the other
three root/covector spans. Its root and covector are outside `H,F`.
The first protected closure and the first fold force `M_10=M_2`.
The three-crossing second fold and the other protected closure force
`M_5=M_13`. All seven remaining internal relations lie in `H,F`.
Lemma `lem:seven` of the cited proof excludes their simultaneous
three-dimensional root and covector spans. No identification of the
two middle output matrices, no tail product cancellation and no
relations involving other blocks enter this argument.

The assertion includes `mu=1`; it uses inverses of `mu`, not of
`1-mu`. Reciprocal rescaling of root/covector pairs preserves every
dependence statement below.

## 2. Independent defects on disjoint strand sets

**Theorem 1 (disjoint-block defect; weaker bound).** Let `b>=4k`, `k>=1`. Embed
`k` copies `U_i` of `U` on the pairwise disjoint sets

```
I_i={4i+1,4i+2,4i+3,4i+4},       0<=i<k.
```

Put `P_k=U_0...U_(k-1)`. Suppose a braid `tau` fixes, as FREE
meridian words, the positions `4i+1,4i+4` for every `i`, and that
`K=closure(P_k*tau)` is a knot. For every finite-dimensional
representation of its group over any field in which meridians have
nonzero displacement rank one, let `r` be the root-span dimension
and `q` the covector-span dimension. Then

```
r<=b, q<=b,          r+q<=2b-k.                         (1)
```

Every irreducible such representation therefore has

```
dim(V)<=b-ceil(k/2).                                   (2)
```

In the conventional complex augmentation setting `U=1,mu!=1`,

```
ar_C(K)<=b-ceil(k/2).                                  (3)
```

**Proof.** Choose nonzero pairs `(u_j,f_j)` for the `b` input
meridians, with `M_j=I-u_j*f_j`. Let

```
A:k^b->V,        A(e_j)=u_j,
C:k^b->V*,       C(e_j)=f_j.
```

Their image dimensions are exactly `r,q`: the input meridians
generate the group, and conjugation by these matrices and their
inverses preserves the input root and covector spans. Consequently

```
dim ker(A)=b-r,       dim ker(C)=b-q.                    (4)
```

The other prefix blocks do not touch block `i`. Because `tau` fixes
both of its outer free meridians, closure identifies precisely the two
outer output matrices of `U_i` with its inputs. The local lemma
therefore supplies a nonzero vector supported in `I_i` in `ker(A)`
or in `ker(C)`. Assign each block to one available side, choosing
just one side if both are available. If there are `a` root-side
blocks and `c` covector-side blocks, then `a+c=k`.

Vectors assigned to one side have disjoint nonempty supports, so are
linearly independent on that side. Thus `b-r>=a` and `b-q>=c`.
Adding proves (1). In an irreducible representation the nonzero root
span is all of `V`. The common fixed subspace is the common kernel
of the input covectors and is zero: otherwise irreducibility would
make the representation trivial. Hence `r=q=dim(V)`, giving (2).
Cornwell's augmentation--KCH correspondence gives (3). This proves
all assertions. The proof uses neither a preserved form nor local
complete reducibility. QED.

Four blocks give `ar_C<=b-2`; two give only `ar_C<=b-1` by this
weaker argument. From a local DISJUNCTION alone it is invalid to count
every block as a root dependence: the supplied alternatives could lie
on opposite sides. The actual signed `U` relations satisfy the stronger
[two-sided local theorem](TWO_SIDED_PROTECTED_OBSTRUCTIONS.md), which
excludes that possibility for these blocks. The checker includes a rank-one irreducible linear-algebra
control with two disjoint four-element blocks and `r=q=7` at `b=8`.
It meets the local dependence ALTERNATIVES, not the signed `U`
relations, so it diagnoses that particular counting inference only.

**Corollary (a preserved form doubles the bound).** If the same
representation preserves a nondegenerate bilinear form, or a
nondegenerate sesquilinear form for a field involution, then

```
r=q<=b-k.
```

Indeed a rank-one isometry `M=I-u*f` fixes `ker(f)` pointwise.
For `v` in that hyperplane and `w` with `f(w)!=0`, invariance of
the form gives `B(u,v)=0` by subtracting `B(Mw,Mv)=B(w,v)`.
Nondegeneracy makes `B(u,-)` nonzero, so it is a nonzero scalar
multiple of `f`. Thus on EVERY block the root and covector ranks
are equal, including for unipotent meridians. The local lemma
then supplies a root dependence on each block, and the same
disjoint-support proof gives `r<=b-k`. The form identifies the
two global spans, giving `r=q`. For a sesquilinear form use the
semilinear identification with the dual; it preserves dimensions
and dependence. General KCH representations have no such
preserved-form premise. The corollary does not improve (3) for
unrestricted complex augmentations.

## 3. Explicit knot closures satisfying all the protection hypotheses

For positions `i<j`, define the band braid

```
g_(i,j)=sigma_i ... sigma_(j-2) sigma_(j-1)
                  sigma_(j-2)^-1 ... sigma_i^-1.         (5)
```

On an interval of entries `(a,c_1,...,c_m,e)`, put
`W=c_1...c_m`. Successive Hurwitz moves in (5) give exactly

```
(W*e*W^-1, c_1,...,c_m, e^-1*W^-1*a*W*e).              (6)
```

Indeed the positive moves first shift `a` rightward by conjugation
through `c_1,...,c_m,e`; the negative moves move `e` leftward while
restoring each `c_l`. Formula (6) proves that ALL interior meridian
words are fixed, and its strand permutation is the transposition
`(i j)`. Exterior entries are unchanged. This includes adjacent
positions, where `W` is the empty product.

Let `m_1<...<m_s`, `s=b-2k`, list all positions except the `2k`
protected outer positions. Define

```
R_(k,b)=g_(m_1,m_2) ... g_(m_(s-1),m_s),
K_(k,b,l)=closure(P_k * R_(k,b)^(1+s*l)),    l>=0.        (7)
```

Each band fixes every protected free meridian by (6). The permutation
of `R_(k,b)` is an `s`-cycle on the mobile positions. Each prefix
block has permutation `(4i+1 4i+2)(4i+3 4i+4)`; hence each of its
two protected outer positions is paired with a different mobile
position. Composing this matching with the mobile cycle inserts
those `2k` outer positions into that cycle. It is a `b`-cycle.
The exponent in (7) has the same mobile permutation as `R_(k,b)`.
Thus every (7) is a knot and Theorem 1 applies.

In particular, (7) with `k=4,b=16,l=0` is an explicit knot with
`ar_C<=14`. No equality of braid index, bridge number or meridional
rank with 16 follows, nor is primeness, hyperbolicity, or parameterwise
distinctness asserted.

## 4. A rigorous obstruction to repeating the old alternating colors

The displayed old four-meridian seed is

```
(01)(23), (14)(23), (05)(16), (12)(67).                   (8)
```

Its eight transposition occurrences have just SEVEN distinct edges,
because `(23)` is repeated. Its support graph is connected on eight
symbols. Every simultaneous conjugate of this seed has the same
one-edge deficit.

**Theorem 2 (edge budget for literal seed stacking).** Suppose the
`b` initial meridian colors are double transpositions, with each of
the `k` disjoint four-strand blocks colored by a simultaneous
conjugate of (8). Allow arbitrary overlaps between their symbol
sets, and arbitrary other double-transposition colors. If their
generated action is transitive on `v` symbols, then

```
v<=2b-k+1.                                             (9)
```

**Proof.** There are `b-4k` other colors, contributing at most
`2(b-4k)` distinct edges. The union from the specified blocks has
at most `7k` edges, even when symbol sets overlap. The whole graph
therefore has at most `2b-k` edges. Transitivity of the generated
permutation action implies connectedness of its transposition-edge
graph: each generator preserves each graph component. A connected
graph on `v` vertices has at least `v-1` edges, proving (9). QED.

In particular, `k>=2` prevents a full `A_(2b)` coloring by this
literal recipe, irrespective of the tail. The usual support-graph
lower bound for a transitive image on `v` symbols is

```
mr(K)>=ceil((v-1)/2),
```

because each meridian supplies at most two edges. Under (9) this
particular lower bound is at most `b-floor(k/2)`, whereas (3) is
`ar_C<=b-ceil(k/2)`. Combining THESE TWO BOUNDS can certify a
gap of at most one (zero at even `k`, one at odd `k`). This does
not bound the actual gap of a knot, exclude a stronger invariant,
or classify different local colorings of the same braid blocks.

Matching the full displayed braid size by the previous literal
alternating seed is impossible. The next result shows that different
double-transposition colors obey the same budget on this disjoint-block
architecture. The stronger [two-sided local defect](TWO_SIDED_PROTECTED_OBSTRUCTIONS.md)
now improves the upper bound to `b-k`, making a gap of two compatible
with this budget for three blocks. Constructing a compatible coloring
and prime/hyperbolic realization remains a further obligation.

### 4.1. The support budget holds for every local coloring

All permutation products in this section mean composition from right
to left. No particular alternating seed, whole-tuple return, or tail
is assumed in the local computation.

**Theorem 3 (universal protected-block coloring budget).** Let
`(g_1,g_2,g_3,g_4)` be four double transpositions on a finite set.
Suppose the literal Hurwitz action of `U` returns its first and fourth
entries to `g_1,g_4`. The multigraph with one edge for EACH of its eight
transposition occurrences is not a forest. Parallel edges from different
colors remain distinct edges.

Consequently, if `k` disjoint protected `U` blocks occur on a `b`-strand
braid as in Theorem 1, EVERY transitive quotient coloring of its input
meridians by double transpositions on `v` symbols satisfies

```
v<=2b-k+1.                                             (10)
```

In particular no such quotient can be the natural `A_(2b)` action when
`k>=2`, regardless of the local colors or their overlapping symbol sets.
The standard support lower bound and Theorem 1 can therefore certify
an augmentation gap at most one by these two inequalities. This does
not bound a knot's actual gap or rule out a different lower invariant.

**Complete finite reduction for the local assertion.** Suppose the
occurrence multigraph were a forest. Each of its nontrivial connected
components is a tree with `n` vertices and `n-1` edges, where `2<=n<=9`:
there are only eight edges in the entire graph. Each color preserves
every such component, and restricts to an involution with zero, one or
two transpositions. Write its edge-count vector as `(a_1,a_2,a_3,a_4)`;
the counts sum to `n-1`. The two outer-return identities restrict to
each component because all propagated colors remain in the generated
permutation group.

The ordered product of these four restricted involutions is an `n`-cycle.
Indeed, expand them into their transposition edges. In ANY order of a
forest's distinct edges, the next edge joins two components of the
preceding forest. Inductively the partial permutation has one cycle on
each preceding component. Multiplication by this next transposition
merges precisely those two cycles. After every edge is inserted the
product is a single cycle. Conversely an `n`-cycle product with exactly
`n-1` transposition occurrences forces a connected graph and hence a
tree, without parallel edges.

Normalize that product to `c_n=(0 1 ... n-1)`. The centralizer of `c_n`
is exactly its powers: a commuting permutation is determined by the image
of zero and sends `c_n^i(0)` to `c_n^i(g(0))`. Thus we may choose the
first involution to be the lexicographically least simultaneous rotation
of its entries. We do NOT claim that this also removes the stabilizer
of the first involution; duplicate full conjugacy orbits are harmless.

For every `n=2,...,9` the checker explicitly generates the entire set
of involutions with at most two transpositions: one identity,
`binom(n,2)` single transpositions, and the three perfect matchings on
every four-element subset. It then enumerates the normalized first
color and all second and third colors, forcing the fourth to be

```
g_4=g_3*g_2*g_1*c_n.                                   (11)
```

It tests fourth-color membership, the total edge count and the literal
fourteen-letter signed braid action. Its exact census is:

| `n` | First-color representatives passing the capacity bound | Minimal factorizations tested | Both outer entries return |
|---:|---:|---:|---:|
| 2 | 2 | 4 | 0 |
| 3 | 2 | 12 | 3 |
| 4 | 5 | 49 | 0 |
| 5 | 6 | 163 | 0 |
| 6 | 13 | 647 | 0 |
| 7 | 19 | 1736 | 0 |
| 8 | 33 | 4395 | 0 |
| 9 | 42 | 6561 | 0 |

The ONLY returning edge-count vector is `(0,0,1,1)`. Its three
normalized representatives all have identity first and second colors,
with two distinct transpositions on three symbols as the last colors.
One is `(1,1,(01),(12))`. Since every nontrivial tree component would
have identity first and second restricted colors, the whole first and
second permutations would be identities. This contradicts their being
double transpositions, and proves the local assertion.


The finite packing formulation is equivalent to the local forest question,
in both directions. A forest tuple restricts to admissible connected tree
components whose edge-count vectors sum to `(2,2,2,2)`. Conversely, given
any such finite collection of admissible components, put their symbol
sets disjointly and take the disjoint union of each of the four colors.
Each resulting color is a double transposition because its total count is
two. The Hurwitz action and both outer-return equalities act componentwise,
so the assembled tuple satisfies them and its occurrence graph is a
forest. Thus infeasibility of that vector packing is exactly the absence
of a local forest tuple, rather than only a sufficient test for absence.
Isolated fixed symbols may be omitted in either direction.

For an efficient exact replay, let `ell(h)=n-number_of_cycles(h)`.
This is the least number of transpositions in a factorization of `h`:
one transposition changes the number of cycles by one, and each cycle
has the usual length-minus-one factorization. The triangle inequality
for this length shows that in a factorization of `c_n` of total length
`n-1`, EVERY successive residual after removing `g_i` has length reduced
by exactly `a_i`. The portable checker caches these complete residual
successor lists before forcing the last factor. This loses no candidate.
The capacity tests merely check that the remaining colors can supply the
remaining length (at most two per color). No sampled tuple, modular
calculation or heuristic group-generation test occurs.

A separate reference search through `n=6` uses a different normalization:
fix the FIRST color to an identity, `(01)` or `(01)(23)` under the full
symmetric group, enumerate all three other colors with the exact edge
budget, and check graph connectivity directly. It uses no product-cycle
normalization or residual-length pruning. It finds the same sole component
vector. This smaller cross-check does not replace the complete main census
through nine.

**Proof of the global consequence.** Form the full occurrence multigraph
of the `b` input meridians, so it has exactly `2b` edges. Its vertices
are the `v` symbols of the transitive action. It is connected because
each generator preserves every connected component of its edge graph.
Over `F_2`, the kernel of the edge-to-vertex boundary map has dimension
`2b-v+1`: the boundary matrix of a connected graph has rank `v-1`, as
seen on a spanning tree, and every column has even coordinate sum.

By the local assertion each protected block has a nonzero cycle vector.
Choose one for every block. These vectors are linearly independent in the
full edge space because their nonempty supports lie in disjoint sets of
OCCURRENCE edges. This remains true if different blocks reuse the same
geometric transposition or the same symbols. Thus `2b-v+1>=k`, proving
(10). For one block the old eight-symbol seed has cycle rank exactly one
and reaches the bound, so this local loss cannot be increased without
additional hypotheses. QED.

The result in particular excludes any transitive nine-symbol local
coloring of `U` by four double transpositions with these two outer
closures. It does not assert that disconnected local actions have total
support at most eight.

### 4.2. The Hurwitz budget and coalescence at genus zero

Consider a connected degree-`n` cover of the projective line of genus
`g`, with `b` finite branch values whose permutations are double
transpositions. Let the permutation over infinity have `s` cycles,
counting fixed points. Suppose the finite branch tuple contains `k`
copies of the four-factor local block of Theorem 3, with both outer
entries returned by `U` and with pairwise disjoint factor positions.
Then

\[
                         k\le 2g+s-1.
\]

The permutation index is two for each finite branch and `n-s` over
infinity. Riemann–Hurwitz gives `2g-2=-2n+2b+(n-s)`, so the cycle rank
of the connected occurrence graph is `2b-n+1=2g+s-1`. The local
cycles from Theorem 3 have disjoint edge supports in the occurrence
graph, even when different blocks reuse sheets. They give `k`
independent cycle vectors. This proof uses the finite branch tuple
and the two outer returns. A knot closure and a returning tail play
no role. For the genus-zero setting, see [Pakovich, equation (1.1)
and Theorem 1.1](https://pakovich.github.io/lau.pdf).

A genus-zero polynomial tuple has `s=1` and admits no such block.
A genus-zero Laurent tuple has `s=2` and admits at most one. If four
consecutive entries form that block, its occurrence graph has rank
one, so it must be one of the [two rank-one local classes, equations
(2) and (3)](PROTECTED_U_DOUBLE_TRANSPOSITION_CLASSIFICATION.md#2-complete-finite-proof-of-theorem-1).
In the first class the first two factors share `(23)` and their
secondary transpositions are `(01),(14)`. In the second they share
`(01)` with secondaries `(23),(24)`. Thus in either case the adjacent
product cancels the shared transposition and is a three-cycle. The
common edge is the unique doubled edge of the full occurrence graph.

After replacing the adjacent pair by its product, expand that
three-cycle into the two secondary edges. The resulting graph is
the original unicyclic graph with both occurrences of its doubled
edge removed. It has two tree components. The tree-product argument
in Section 4.1 shows that the unchanged total permutation has exactly
one cycle on each component. For pole type `(n-1,1)`, one component
is therefore a singleton. Deleting it gives the genus-zero polynomial
tuple on `n-1` sheets immediately. For other two-cycle pole types the
same replacement gives the two tree components with those cycle
lengths. The connected rank-one class in equation (2) has eight
sheets and pole type `(7,1)`, so `k=1` is attained. This local test
applies to a displayed protected block. The general `(n-1,1)`
Laurent reduction is established by the [cycle-folding theorem](../../arithmetic/ALTERNATING_HURWITZ_CYCLE_FOLDING.md).

## 5. Primary input, checks and assurance

The indispensable local signed-matrix input is Lemmas `lem:spans`
and `lem:seven` of the linked current released proof. It is an
internal written result with the assurance recorded by that paper.
The correspondence input is C. R. Cornwell, *KCH representations,
augmentations, and A-polynomials*, J. Symplectic Geom. 15 (2017),
983--1017, Theorem 1.2 and Definition 1.8, inspected in
[arXiv:1310.7526v3](https://arxiv.org/html/1310.7526v3#S1).
Its Theorem 1.9 is the existing connected-sum additivity result;
the present note establishes protection and disjoint-support
counting directly and does not invoke meridional-rank additivity.
No priority assertion is made.

The [standard-library checker](verify_multiple_protected_obstructions.py)
tests exact free-group identities for (6), the literal prefix
permutation and the knot permutations in (7), the disjoint-support
dimension calculation, an exact irreducible linear-algebra negative
control for the unjustified two-block/two-defect inference, and the
edge deficit in (8). The written arguments above establish the
unrestricted assertions; finite controls do not replace them.

```
.venv/bin/python research/meridional-rank/representations/verify_multiple_protected_obstructions.py \
  --write-receipt research/meridional-rank/representations/multiple_protected_obstructions_verification.json
```

The receipt records source pins and exact scopes. No local signed
obstruction replay, independent mathematical review, formal
verification, prime realization, or matching meridional lower bound
is inferred from this checker.

The separate [complete forest-coloring checker](verify_protected_block_forest_colorings.py)
and [receipt](protected_block_forest_colorings_verification.json) supply
the exhaustive finite part of Theorem 3. It reconstructs the census
through nine, the independent first-color normalization through six,
all admitted component vectors and their bounded packing. Its controls
retain the valid old eight-symbol coloring, reject an explicit
nine-symbol path coloring, distinguish one outer closure from both,
and verify inversion of the signed Hurwitz action. The finite proof
does not establish a new prime knot or an augmentation gap two.

```
.venv/bin/python research/meridional-rank/representations/verify_protected_block_forest_colorings.py \
  --write-receipt research/meridional-rank/representations/protected_block_forest_colorings_verification.json
```

The elementary tree factorization and cycle-space arguments are proved
above; no external classification of Hurwitz orbits is imported. Primary
comparison with Cornwell's results motivates the rank question but does
not supply this local census. The exact finite replay and written
reduction are internal evidence; formal verification, independent
external review and historical priority are not claimed.
