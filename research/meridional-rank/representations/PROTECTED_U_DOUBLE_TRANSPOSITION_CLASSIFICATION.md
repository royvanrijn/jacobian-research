# Protected U double-transposition colors and the tail obstruction

Every double-transposition coloring of disjoint protected copies of `U`
obeys a stronger support budget after a knot-closing tail:

```
v<=2b-2k+2.
```

Here `b` is the strand count, `k` the number of blocks, and `v` the
support size of the quotient colors. The tail fixes the outer meridians
as exact free words. Combined with `ar_C<=b-k`, the standard support
lower bound can therefore certify a gap of at most one. This holds for
all local colorings, all numbers of blocks and arbitrary extra strands.
It does not bound a knot's actual meridional rank or actual gap.

The proof has two parts: a complete classification of 31 local color
tuples, and an unrestricted permutation argument controlling how a
knot-closing tail can match their mobile colors. The classification does
not assume connected local support.

## 1. Local classification and the component-slot certificate

Use the left-to-right Hurwitz convention `T(a,c)=(c,c^-1*a*c)` and

```
U=[1,2,-3,-3,2,-3,-3,2,1,1,2,2,2,-3].
```

A local tuple consists of four double transpositions `(A,X,Y,C)` such
that the first and fourth output entries of `U` return to `A,C`.
Its occurrence graph has one edge for EACH of the eight transposition
occurrences. Parallel edges remain distinct. Delete symbols fixed by
every color. Simultaneous conjugation may relabel symbols but does not
permute the four color positions.

**Theorem 1 (complete local classification).** There are exactly 31
simultaneous-conjugacy classes of local tuples. Their occurrence-cycle
ranks have the following distribution:

| Cycle rank | 1 | 2 | 3 | 4 | 5 | 6 |
|---:|---:|---:|---:|---:|---:|---:|
| Local classes | 2 | 4 | 11 | 10 | 3 | 1 |

For every class let `G_0` be the transposition support graph of its two
INPUT mobile colors `X,Y`, and let `t` count full-support symbols absent
from `G_0`. Each transposition in either OUTPUT mobile color has both
endpoints in a single component of `G_0`. Order the two transpositions
inside each input and output color. There exist an output ordering and
a permutation `tau` of the four occurrence slots such that

```
component(output slot i)=component(input slot tau(i)),
ell(tau)+t<=2,                                         (1)
```

where `ell(tau)=4-number_of_cycles(tau)` is its minimum transposition
length. Exact certificate counts are:

| `ell(tau)` | Omitted symbols `t` | Classes |
|---:|---:|---:|
|0|0|17|
|0|1|9|
|0|2|2|
|1|0|1|
|1|1|2|

The [receipt](protected_u_double_transposition_classification_verification.json)
contains a representative of every class, its literal prefix output,
component partition, chosen output-edge ordering and slot permutation.
All permutations in it are exact finite arrays.

Geometric edge equality is stronger than (1) and is unnecessary. Seven
classes change actual mobile transposition edges. In each of those
classes the local input-component signatures stay unchanged, so its
component-slot cost is zero. The proof below uses component membership,
not equality of geometric edges.

## 2. Complete finite proof of Theorem 1

Restrict the tuple to each nontrivial connected component of its full
occurrence graph. Every restricted color is an involution with zero,
one or two transpositions. The protected outer identities restrict
componentwise, since the Hurwitz operations are group operations.

The complete tree-component theorem in
[Multiple protected obstructions](MULTIPLE_PROTECTED_OBSTRUCTIONS.md),
Theorem 3, shows that the ONLY admissible nontrivial tree component has
three vertices and edge-count vector `(0,0,1,1)`: its last two colors
are distinct transpositions sharing a vertex. Its proof and complete
forest checker are indispensable inputs here.

Every cyclic component with `n` vertices has at least `n` occurrence
edges. There are at most eight edges in total, so `2<=n<=8`.
For every such `n`, enumerate its restricted first color up to `S_n`:
identity, `(01)`, or `(01)(23)`, omitting unavailable classes.
Enumerate the other three involutions without quotienting: one
identity, every transposition, and all three pairings of each
four-element subset. Retain tuples whose total edge count is from
`n` through eight. Test the literal fourteen signed Hurwitz moves,
both outer returns and graph connectivity. Finally quotient successful
tuples by the full centralizer of their fixed first color.

This is complete: every cyclic component normalizes its first color
in this way, and any conjugacy between two normalized tuples lies in
that centralizer. No candidate is dropped based on a guessed group
or sampled representative. The complete census is:

| `n` | First-color edges | Tuples tested | Outer returns | Connected | Orbits |
|---:|---:|---:|---:|---:|---:|
|2|0|4|1|1|1|
|2|1|7|2|2|2|
|3|0|27|0|0|0|
|3|1|54|3|2|1|
|4|0|540|3|0|0|
|4|1|864|9|2|1|
|4|2|981|18|14|3|
|5|0|10125|0|0|0|
|5|1|15300|45|24|2|
|5|2|17200|49|24|3|
|6|0|91125|0|0|0|
|6|1|182250|237|24|1|
|6|2|218700|196|88|6|
|7|1|1157625|0|0|0|
|7|2|1852200|573|48|1|
|8|2|9261000|1906|192|1|

The absent first-color rows cannot meet the edge budget or have no
such involution. This gives 22 cyclic component classes. Together
with the sole tree class they give every possible nontrivial component.
Each has at least two edges, so a full tuple has at most four components.
Enumerate multisets of at most four component classes whose edge-count
vectors sum to `(2,2,2,2)`; take their disjoint unions color by color.
This produces exactly 31 tuples and the rank distribution in Theorem 1.

Both directions of this assembly matter. Every original tuple has this
component multiset. Conversely a disjoint union of admitted components
with that total vector is a tuple of four DTs, and both protected
relations hold componentwise. Distinct multisets give distinct full
classes, because a simultaneous conjugacy must carry connected
components to isomorphic ordered-color components.

For each full representative, replay `U` and compute the components
of its two input mobile colors. Check every output transposition
stays within one such component. Try both orders for the two edges
of each output color and all 24 slot permutations. Recording a
minimum-length admissible slot matching gives exactly the table in
Theorem 1, including (1). This is a finite exact certificate check,
not an extrapolation from small support. It establishes Theorem 1.

For comparison, the two rank-one classes are

```
((01)(23),(14)(23),(05)(16),(12)(67)),                   (2)
((01)(23),(01)(24),(35)(67),(02)(78)).                   (3)
```

The first is connected on eight symbols. The second has a connected
six-symbol part and a three-symbol tree, hence disconnected support
of size nine. Thus connected-support assumptions would omit a real
class. The classification contains support up to ten symbols.

## 3. An unrestricted tail theorem

Let `b=4k+s`, with `k>=1,s>=0`. Put `k` copies of `U` on the first
`4k` strands, in disjoint blocks, and call their product `P_k`.
Protect the first and fourth free meridians of every block. All other
positions, including the `s` extras, are mobile; their count is
`m=2k+s`.

For consecutive mobile positions `i<j`, use the canonical band

```
g_(i,j)=sigma_i ... sigma_(j-1) sigma_(j-2)^-1 ... sigma_i^-1.
```

It fixes every intervening free meridian. On its endpoints with
intervening word `W` its exact action is
`(a,e)->(W*e*W^-1,e^-1*W^-1*a*W*e)`.

**Theorem 2 (universal mobile-tail support budget).** Suppose a tail
`beta` belongs to the subgroup generated by these consecutive-mobile
bands. Suppose a DT coloring of the `b` input meridians returns under
`P_k*beta`, and its braid closure is a knot. If `v` is its full support
size, then

```
v<=6k+2s+2=2b-2k+2.                                   (4)
```

Arbitrary overlaps between block supports and arbitrary extra DT
colors are allowed.

**Proof.** The tail fixes the protected free meridians, so each local
block has both outer returns. Theorem 1 applies separately to every
block.

Frame the mobile entries to turn consecutive-mobile bands into
ordinary Hurwitz generators on `m` entries. For block `j` with outer
colors `A_j,C_j`, its two mobile entries share the conjugating frame

```
w_0=1,  w_j=C_0*A_1*C_1*A_2*...*C_(j-1)*A_j.
```

Across the extra positions use the last frame multiplied by the last
protected color. The band formula proves this equality of actions
exactly. The frames are the same before and after the prefix, since
the protected outer colors return.

Let `G` be the support occurrence graph of the FRAMED INPUT mobile
colors, with `2m=4k+2s` edge occurrences. Its components are invariant
under the group `H` generated by those colors: each generator preserves
every component. A local input-mobile component, after framing, is a
connected subgraph of `G`. Hence membership in the same local component
implies membership in the same global component, under all overlaps.
In particular the local certificates (1) remain valid for component
labels in `G`, including for classes that change geometric edges.

For block `j` choose its certificate permutation `tau_j`, of length
`c_j`, and omitted-symbol count `t_j`, so `c_j+t_j<=2`.
Take their disjoint product `tau` on the `2m` mobile edge slots, fixing
all extra-color slots. Put `C=sum c_j,T=sum t_j`. Thus

```
ell(tau)<=C,     C+T<=2k.                              (5)
```

Write `L(q)` for the component of `G` containing input slot `q`.
Order output slots as in the local certificates. Every output slot
lies within `G`, and its label equals `L(tau(q))`.

Each prefix-output color preserves every component of `G`, because its
two transposition edges stay inside those components. Every subsequent
Hurwitz conjugating word therefore preserves the same components. This
is all the invariance needed below; membership in the subgroup generated
by the original local mobile colors is unnecessary.

Let `rho` be the tail permutation on the `m` mobile positions.
It is one `m`-cycle. Indeed the prefix pairs each protected outer
position with a distinct mobile position and fixes every extra
position. Composing that matching with the tail inserts the outer
positions into the cycles of `rho`, preserving the number of cycles.
A knot closure therefore has precisely one mobile cycle.

Match the two transposition occurrences of each prefix-output color
to those of the final corresponding input color using its component-preserving conjugation. This gives a slot permutation `pi` lifting the one-cycle color
permutation `rho`. Such a lift has at most two cycles: after one trip
around all `m` pairs it either preserves each pair's edge index or
interchanges the two indices. Component invariance gives

```
L(tau(q))=L(pi(q)).
```

Thus `L` is constant on cycles of `pi*tau^-1`. Every nonempty component
of `G` contains a slot, so its component count satisfies

```
c(G)<=cycles(pi*tau^-1)<=cycles(pi)+C<=2+C.              (6)
```

The middle inequality follows because multiplying by one transposition
changes the number of cycles by exactly one, and `tau` has a
factorization of length at most `C`. The mobile support size obeys

```
v_mobile<=2m+c(G)<=4k+2s+2+C.                           (7)
```

It remains to compare framed full supports with original full supports.
The cumulative union is unchanged. For each new block strip the final
`A_j` from `w_j`, since it preserves its own full support; the remaining
word is supported on earlier block sets. It permutes their cumulative
union and fixes its complement. Inductively adjoining the framed block
therefore gives the same cumulative full union. The extra frames are
supported on that block union, so adjoining extra-color supports also
preserves the total union.

Each block's full support omits at most `t_j` symbols from its mobile
pair, and extras omit none. Consequently

```
v<=v_mobile+T<=4k+2s+2+C+T<=6k+2s+2,
```

by (5)--(7). This proves the theorem. QED.

## 4. Exact free-meridian stabilizers and the limitation of support bounds

**Corollary 3 (all exact protected tails).** Theorem 2 holds whenever
`beta` fixes each protected free meridian as an exact word. This uses
[Theorem FreeMeridianStabilizer](FREE_MERIDIAN_STABILIZERS.md): its full
pointwise stabilizer is exactly the canonical mobile-band subgroup.
This imported geometric theorem is a separate indispensable premise;
its proof and primary-source ledger accompany that note.

For every knot in this protected architecture, the
[two-sided rank-one defect theorem](TWO_SIDED_PROTECTED_OBSTRUCTIONS.md)
gives `ar_C<=b-k`. A transitive DT quotient on `v` symbols supplies the
usual support lower bound `mr>=ceil((v-1)/2)`. By (4),

```
ceil((v-1)/2) - (b-k) <= 1.                             (8)
```

Therefore THESE TWO BOUNDS cannot certify a gap of two, for any `k`,
any number of extra strands, or any local DT coloring. In particular
adding a local cycle-rank-two or cycle-rank-three block cannot repair
this route. The statement leaves open different finite meridian
classes, different local braids or protection hypotheses, overlapping
obstructions, and stronger meridional lower invariants. It does not
exclude an actual gap-two knot within the architecture.

## 5. Exact replay and assurance

The [checker](verify_protected_u_double_transposition_classification.py)
uses exact finite permutations and integer graph components. NumPy
only batches the exhaustive tuple arrays; there is no floating-point
calculation or numerical tolerance. It replays the existing complete
forest census, enumerates every cyclic component, constructs all 31
full classes and their slot certificates, rejects crossing-sign
mutations, and checks 65,128 finite cases of the slot-cycle bound,
including extra mobile colors. A non-knot control violates that bound
when the color permutation has several cycles, confirming that the
closure hypothesis is used. The proof in Section 3 establishes the
all-`k` claim; the small slot controls do not replace it.

```
.venv/bin/python research/meridional-rank/representations/verify_protected_u_double_transposition_classification.py \
  --write-receipt research/meridional-rank/representations/protected_u_double_transposition_classification_verification.json
```

A separately implemented C++ exhaustive component search and an
independent full-component assembly agreed with the NumPy census
and all 31 local properties during internal review. The receipt pins
current proof sources and records its exact computational scope.
No formal verification, external review, historical priority or new
prime/hyperbolic gap-two example is asserted.
