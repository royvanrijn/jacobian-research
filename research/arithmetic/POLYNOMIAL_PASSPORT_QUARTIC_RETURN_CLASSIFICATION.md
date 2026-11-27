# Return braids and a complete four-sheet gluing criterion

Fix the [polynomial-passport conventions](POLYNOMIAL_PASSPORT_COMPONENTS.md):
\(g_1\cdots g_r=P_N\), total finite index \(N-1\), full braids, and a
preserved congruence partition into \(d\) blocks of size \(b\).
The [literal return construction](POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md)
assigns \(M=(r-1)d+1\) common-fibre returns with product \(P_b\).
Identity returns and identity quotient factors are retained.

The induced action is always an inner braid action, although only a
particular subgroup of inner braids is available. This distinction
gives a complete four-sheet gluing criterion. When \(b=4\), the only
ambiguous return pattern consists of three transpositions. Its extra
components are precisely the orbits of an explicit subgroup of \(A_4\)
on four labels. The quotient degree and number of original factors are
unrestricted.

This resolves the return-action problem for fibre size four by a
four-state criterion over the decorated quotient graph. It does not
classify arbitrary larger-fibre actions or primitive polynomial
passports.

## Primary source and the remaining distinction

Artin's [*Theory of Braids*, Theorems 15--16 and their proof,
printed pp113--115](sources/artin-theory-of-braids-1947.pdf), characterizes
braid substitutions on a free group: each basis generator is sent to a
conjugate of a positive basis generator, and their ordered product is
fixed. The proof on pp114--115 constructs a braid by shortening the
conjugating words. Formula (14) uses the same generator convention
\((x_i,x_{i+1})\mapsto(x_{i+1},x_{i+1}^{-1}x_ix_{i+1})\) as here.
Theorems 11--14 provide faithfulness. These actual statements and the
sufficiency proof were inspected in the retained primary source.

Artin's theorem determines which free-group substitutions are braids.
It does not say that a quotient stabilizer supplies every inner braid.
That subgroup restriction remains essential: the degree-eight
obstruction in L8 has one inner component and two glued components.
No historical novelty is asserted.

## Theorem R1. Universal realization of return transport by inner braids

For a fixed literal quotient tuple \(q\), replace its \(M\) returns by
the free generators \(x_1,\ldots,x_M\) of \(F_M\), in the boundary order
of L2. Write \(Z=x_1\cdots x_M\).
The lift, original signed Hurwitz moves, and extraction of L1--L3
can all be performed with free-group words as fibre transports.
For an original braid word \(\beta\), from \(q\) to \(q'\), let
\[
       \Phi_{q,\beta}(x_j)=W_j(x_1,\ldots,x_M)                 \tag{1}
\]
be its extracted words. Then
\[
 \begin{split}
  W_1\cdots W_M&=Z,\\
  W_j&=U_j^{-1}x_{\sigma(j)}U_j
               \quad(\sigma\in S_M).                        \tag{2}
 \end{split}
\]
Consequently (1) is an Artin braid substitution. It is independent
of the fibre degree, and its evaluation at any admissible return
tuple is exactly the original return transport.

In particular the quotient stabilizer acts through a specific
subgroup of \(B_M\). For successive words, the precise substitution
rule is
\[
 \Phi_{q,\beta\gamma}
       =\Phi_{q,\beta}\circ\Phi_{q\beta,\gamma}.               \tag{3}
\]
This convention agrees with executing tuple moves from left to right.
It gives a homomorphism on a fixed quotient stabilizer, with the
corresponding Artin composition convention.

**Proof.** The tree normalization and product straightening use only
group multiplication and inverse. They therefore work over every
group, without an index hypothesis on its elements. With returns
\(x_j\), straightening makes the full product have identity transport
between consecutive blocks, except for transport \(Z\) at the last
block. Original Hurwitz moves preserve this product. Its \(d\)-th
power on root block zero is \(Z\), so the boundary calculation of
L2 proves the first identity in (2) after extraction.

A Hurwitz move carries each quotient cycle either unchanged to the
other position or to its conjugate by the adjacent quotient
permutation. The corresponding whole-cycle return is conjugated by
a fibre transport. Changing the normalized tree charts conjugates
it again by words in the \(x_j\). Every old cycle slot occurs exactly
once among the new slots. Initially its return is the corresponding
\(x_j\). This proves the second identity in (2), including positivity
and the fact that \(\sigma\) is a permutation.

Artin Theorem 16 now supplies a braid. Equivalently, one can first
see invertibility by applying the opposite original word and using
the literal inverse construction L2 over the same free group.
All computations commute with group homomorphisms from \(F_M\), so
evaluation gives the exact permutation transport for every \(b\).
Inserting the normalized literal inverse between successive original
words proves (3). Original braid relations hold in wreath coordinates;
Artin faithfulness makes these equalities equalities of inner braids.
\(\square\)

Deleting identity returns has the usual braid meaning. Set the
corresponding \(x_j\) to one and delete their image slots. The remaining
images still satisfy (2) on the remaining, naturally ordered
generators. Artin's theorem therefore gives a braid on just the
nonidentity return slots. This also follows from the string-deletion
operation described by Artin on printed p116.

## Corollary R2. The inner component and its monodromy

The full inner braid orbit of the return tuple is an invariant of
every original component preserving the specified blocks.
The subgroup \(J=\langle R_q(g)\rangle\le S_b\) is literally invariant,
and equals the group induced on root block zero by its setwise
stabilizer in \(G=\langle g_1,\ldots,g_r\rangle\).

**Proof.** R1 realizes each original move by an inner braid. Inner
Hurwitz moves preserve the literal generated subgroup, proving the
first two assertions.

For the last assertion, consider the directed graph on quotient
blocks with one edge \(a\to q_i(a)\) for every pair \((i,a)\).
The retained edges of L1 are a spanning tree. Their paths from zero
give words in the original factors that carry root block zero to
the indicated blocks. Traversing a deleted edge and returning along
tree paths gives a word in its setwise stabilizer; its induced
permutation is precisely that normalized return.

Conversely any word stabilizing root block zero traces a closed
walk in this graph. Insert the tree path to zero at each successive
edge. The walk becomes a product of these deleted-edge loops and
their inverses. Its root-block action is therefore generated by the
returns. Since the root chart is the original chart, this is equality
as literal subgroups of \(S_b\). \(\square\)

The invariant is necessary, not generally sufficient. L8 already
has the same inner full braid orbit on both glued components.

## Lemma R3. The four-sheet return fibres

For \(b=4\), a decorated quotient has one of the following nonidentity
return patterns. Here \(T\) is a transposition, \(C\) a 3-cycle, and
\(D\) a double transposition.

| Return pattern | Literal fibre size | Orbits under conjugation by \(P_4\) |
|---|---:|---:|
| one 4-cycle | 1 | 1 |
| \(C,T\), in either specified order | 4 | 1 |
| \(D,T\), in either specified order | 2 | 1 |
| \(T,T,T\) | 16 | 4, each of size 4 |

Every other return is the identity. Thus, except for the last row,
the decorated quotient component is sufficient for an original
component. In the last row there are at most four components above
one decorated component.

**Proof.** The total return index is three, giving exactly the
displayed partitions of the index budget.
A sole return is forced to be \(P_4\). For two returns, choose the
transposition. The four adjacent chords of the \(P_4\)-cycle leave
a 3-cycle on completion; the two diagonals leave a DT. Conjugation
by \(P_4\) is transitive on each such set.

For three transpositions, choosing the first adjacent transposition
leaves a 3-cycle, which has three ordered two-transposition
factorizations. Choosing a diagonal leaves a DT, which has two.
There are therefore \(4\cdot3+2\cdot2=16\) triples.
No triple is fixed by a nontrivial power of \(P_4\): a quarter-turn
fixes no transposition, while the only transpositions fixed by a
half-turn are the two diagonals, which cannot generate \(P_4\).
So there are four rotation classes.

The \(d\)-th original full twist supplies simultaneous \(P_4\)
conjugation by L3. Transport a connecting decorated quotient braid
and then correct this phase. This proves sufficiency in the first
three rows and the four-component bound in the last. \(\square\)

## Theorem R4. A complete \(A_4\) criterion for four-sheet gluing

Let \(Y\) be one full decorated quotient braid component whose returns
have pattern \(T,T,T\). Order the three active slots by the boundary
order of L2. Label their four simultaneous-\(P_4\)-conjugacy classes
by the representatives
\[
 \begin{array}{c|ccc}
 0&(2\,3)&(1\,3)&(0\,3)\\
 1&(2\,3)&(0\,1)&(1\,3)\\
 2&(2\,3)&(0\,3)&(0\,1)\\
 3&(1\,3)&(1\,2)&(0\,3).
 \end{array}                                                 \tag{4}
\]
Each row has product \(P_4=(0\,1\,2\,3)\).
Every signed decorated quotient edge acts on these four labels by
an element of \(A_4\), computable from L3 on the four representatives.

Choose a base vertex \(y_0\in Y\) and paths to every vertex. Compose
the edge permutations along each path to obtain \(S_y\).
For each signed edge \(e:y\to y'\), with label permutation \(E_e\),
form
\[
           S_{y'}^{-1} E_e S_y,\qquad
 H_Y=\left\langle S_{y'}^{-1}E_eS_y:e\right\rangle\le A_4.     \tag{5}
\]
Original components above \(Y\) are exactly the orbits of \(H_Y\)
on the four labels. In particular:

| \(H_Y\), up to conjugacy | Label-orbit sizes | Original component sizes |
|---|---|---|
| \(1\) | \(1,1,1,1\) | \(4\lvert Y\rvert,4\lvert Y\rvert,4\lvert Y\rvert,4\lvert Y\rvert\) |
| \(C_2\) | \(2,2\) | \(8\lvert Y\rvert,8\lvert Y\rvert\) |
| \(C_3\) | \(1,3\) | \(4\lvert Y\rvert,12\lvert Y\rvert\) |
| \(V_4\) or \(A_4\) | \(4\) | \(16\lvert Y\rvert\) |

Together with R3, this is a complete component criterion for all
block-size-four sectors of polynomial passports. It requires the
decorated quotient graph and four-state edge actions, rather than
an original-degree braid search. It does not assert that the
decoration alone, or the abstract group \(H_Y\) without its orbits,
labels an individual component.

**Proof.** On the four classes (4), direct substitution in the inner
Hurwitz formula gives
\[
                  Q_1=(0\,3\,2),\qquad Q_2=(0\,2\,1).         \tag{6}
\]
They generate \(A_4\): the subgroup is transitive on four letters,
has an element of order three, and is contained in \(A_4\); its order
is therefore twelve.

For an arbitrary original edge, R1 followed by deletion of identity
returns gives a braid on the three active slots. Its action on (4)
is consequently in this \(A_4\). Transport commutes with simultaneous
conjugation, so these four-state edge permutations are well defined.
Formula (5) is the full Schreier action on these labels, by the same
path insertion argument as L3. The root phase connects all four
literal triples within each label. Thus no additional component is
lost on passing to four states, and no distinct \(H_Y\)-orbits can
be connected.

An \(H_Y\)-orbit of size \(h\) gives \(h\) labels and \(4h\) literal
returns over each vertex of \(Y\), proving the count \(4h|Y|\).
The subgroups of \(A_4\) are \(1,C_2,C_3,V_4,A_4\): a subgroup with
no 3-cycle lies in \(V_4\); a proper subgroup containing a 3-cycle
can contain no nonidentity double transposition or different
order-three subgroup, since those generate \(A_4\).
Their point orbits give exactly the table. \(\square\)

Deck groups are not discarded by this criterion. After choosing any
label in an \(H_Y\)-orbit, reconstruct its literal tuple with L2 and
compute its commuting powers of \(P_N\). This exact cyclic stabilizer
gives the absolute orbit size as the literal size times the stabilizer
order divided by \(N\), as in H4 and L3. Different block systems are
handled by the canonical-divisor filter L9; no uniqueness of a
polynomial decomposition chain is assumed.

## Theorem R5. Every inner component obstruction survives a cyclic outer lift

Let \(h=(h_1,\ldots,h_s)\) be a degree-\(b\) polynomial tuple, with
literal group \(J\), and let \(d\ge2\). There is a degree-\(bd\)
polynomial tuple with passport
\[
      (d^b),\qquad
      \operatorname{type}(h_i)\cup 1^{(d-1)b}
                            \quad(1\le i\le s),               \tag{7}
\]
whose literal group is the standard embedded \(J\wr C_d\).
Its quotient is \((P_d,1,\ldots,1)\), and deleting the identity
entries from its return tuple gives exactly \(h\).

If two such inner tuples have the same ordered cycle partitions
but distinct full inner braid components, their outer lifts belong
to distinct original components of the same passport (7).
Thus every inner component split gives a split for infinitely many
higher-degree passports with one additional finite branch value.
This is an obstruction-inheritance theorem, not a claim that these
are all the components of (7).

**Proof.** Use the normal coordinates of L1 with
\(q=(P_d,1,\ldots,1)\). Give every transport of \(P_d\) the identity
value. At the fixed root-block slot of identity factor \(i\), put
\(h_i\); all its other slots have identity return. In the boundary
order of L2, deleting identities leaves \(h_1,\ldots,h_s\).
Its product is \(P_b\), so L2 reconstructs the literal full product
\(P_{bd}\) and gives exactly (7).

Before product straightening, the generators are a pure block
rotation and the local \(h_i\) on root block zero. The latter generate
the full local \(J\); conjugating by block rotations gives independent
copies of \(J\) on every block. The group is therefore \(J\wr C_d\).
The straightening charts are powers of \(P_b\), since the only
nonidentity full-product transport is \(h_1\cdots h_s=P_b\).
As \(P_b\in J\), these charts belong to the base group and preserve
this standard embedded wreath group.

If an original braid connected the two lifted tuples, it would
preserve their chosen congruence partition. R1 would connect their
return tuples by an inner braid. Deleting identity strands would
then connect the original inner tuples, a contradiction. \(\square\)

For example H6 gives six distinct degree-eleven inner components
with the same literal \(A_{11}\). Thus, for every \(d\ge2\), the
passport
\[
             (d^{11}),\quad(3^2 1^{11d-6}),\quad
                              (3^3 1^{11d-9})                \tag{8}
\]
has at least six components generating the same standard embedded
\(A_{11}\wr C_d\). All six have the same decorated quotient and the
same inner cycle partitions. For odd \(d\), every factor has odd
order, and H5 gives the common standard spin \(+1\): its exponent is
\(11(d^2-1)/8+5+((11d)^2-1)/8=(33d^2+7)/2\), which is even.
The assertions concern these six exhibited components; other
components or block descriptions of (8) are not excluded.

## Proposition R6. All four component patterns occur

In the table below, a slot \((i,O)\) uses a one-based branch index
and zero-based quotient sheets. Put a transposition return at each
of the three indicated slots and identity returns elsewhere.
The return values can be any row of (4). Let \(Y\) be the resulting
decorated quotient component.

| Literal quotient \(q\) | Three active slots | \(\lvert Y\rvert\) | \(H_Y\) | Literal component sizes above \(Y\) |
|---|---|---:|---|---|
| \(((0\,1),1)\) | \((1,(0\,1)),(2,(1)),(2,(0))\) | 2 | \(C_2\) | \(16,16\) |
| \(((0\,1),(1\,2))\) on three sheets | \((1,(0\,1)),(2,(1\,2)),(1,(2))\) | 6 | \(1\) | \(24,24,24,24\) |
| \(((0\,1),1,1)\) | \((1,(0\,1)),(2,(1)),(3,(1))\) | 12 | \(A_4\) | \(192\) |
| \(((0\,1),(1\,2),(2\,3))\) | \((1,(3)),(2,(3)),(3,(1))\) | 64 | \(C_3\) | \(256,768\) |

In particular the four-component bound is attained. The \(1+3\)
case occurs with three nonidentity quotient factors, rather than
only in a rigid two-factor problem.

More fully, for the degree-sixteen passport
\[
                         (2^5 1^6)^3,                       \tag{9}
\]
the entire block-size-four sector has four components of literal
sizes \(256,256,768,768\), with trivial deck groups and absolute
sizes \(16,16,48,48\). Each generated group has block-count
divisors exactly \(1,4,16\); the selected fibre in L9 is therefore
really four. This statement classifies this sector and does not
exclude components of (9) without that block system.

**Proof.** For the four specified inputs, enumerate the entire
decorated quotient graph by both signed generators at every
position and compute (5). This gives the displayed vertex counts
and groups. Independently reconstruct the original tuples over
all four labels and all four root phases at every quotient vertex.
Their direct original Hurwitz closures agree exactly with the
four-state components and have the displayed sizes.

For (9), each quotient factor has index at most
\(\lfloor5/4\rfloor=1\), and the quotient index total is three.
Thus all three quotient factors are transpositions on four sheets.
Every such quotient tuple is in the sixteen-element component of
minimal three-transposition factorizations of \(P_4\), already
enumerated in R3--R4. Each original factor has return index one.
It cannot have its transposition return on its nontrivial quotient
cycle, since that would produce a 4-cycle upstairs. It therefore
has one transposition return at one of its two fixed quotient
points. There are exactly \(16\cdot2^3=128\) decorated tuples.
Their full graph has two components of 64 vertices, each with
\(H_Y=C_3\). Formula R4 gives the four asserted literal counts,
and direct original closure checks all \(128\cdot16=2048\) tuples.
The congruence-divisor and commuting-rotation tests verify the
stated block lattices and deck groups for every component.
\(\square\)

## Exact controls and remaining scope

The [checker](verify_polynomial_return_braids.py),
[receipt](polynomial_return_braids_verification.json), and
[source ledger](polynomial_return_braids_source_ledger.json) retain the
free-group substitutions, boundary and peripheral checks, inverses,
braid relations, identity deletion, finite \(A_4\) action, and direct
original component controls. The uniform assertion R1 uses the
written argument and Artin's actual sufficiency theorem; finite word
tests do not replace that proof.

The component criterion resolves fibre size four. General primitive
\(S_N/A_N\) components, return actions for fibres at least five outside
the previously classified cyclic family, and reduced lifting
invariants remain open. All present proof audit and computation are
internal; independent review, independent replay and formal
verification are not asserted.
