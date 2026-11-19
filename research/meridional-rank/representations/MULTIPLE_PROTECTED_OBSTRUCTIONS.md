# Multiple protected rank-one obstructions and the coloring budget

The local four-strand obstruction can be packed on disjoint sets of
top strands. The resulting loss is a sum of root and covector defects:
`k` protected blocks force an irreducible dimension loss of at least
`ceil(k/2)`. This gives explicit knots with such augmentation upper
bounds. It does not yet give a matching meridional-rank lower bound or
a prime knot with augmentation gap two.

There is a separate obstruction to the most immediate coloring attempt.
Repeating the displayed alternating eight-symbol seed loses one distinct
transposition edge per block. Its support-graph lower bound cannot certify
an augmentation gap greater than one from the new dimension bound.

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

**Local lemma.** Suppose the matrices propagated through one copy of
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

**Theorem 1 (disjoint-block defect).** Let `b>=4k`, `k>=1`. Embed
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
argument. It is invalid to count every local obstruction as a root
dependence: different blocks can supply their defects on opposite
sides. The checker includes a rank-one irreducible linear-algebra
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

Thus four protected blocks genuinely supply an augmentation upper
bound two below the displayed braid size, but matching that size by
the previous literal alternating seed is impossible. A gap-two
construction needs new colorings, a stronger local two-sided defect,
a more efficient packing, or a different meridional lower invariant.
Its prime/hyperbolic realization is a further obligation.

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
