# Alternating subgroup covers and a uniform Hasse family

The regular rational degree-ten pencil produces infinitely many finite étale
Q-algebras of rank 301 which have a point over every completion of Q and have
no Q-point. Their splitting fields can be chosen linearly disjoint with group
A10. The construction controls ramified decomposition groups, as well as
individual Frobenius elements. Each algebra is a complete Keller fibre by the
retained finite-étale compiler theorem. A [single fixed four-dimensional
map](../jacobian/monodromy-theory/FIXED_KELLER_MAP_FINITE_ETALE_FAMILIES.md)
now realizes the entire product family; three-dimensional realization
remains a separate question.

The first input is the [degree-ten regular A10 pencil](RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md),
Theorem 1:

\[
\begin{split}
H(z)&=9z^{10}-864z^8+28161z^6-350150z^4+1127500z^2,\\
P_s(z)&=H(z)-sz+250000,\qquad f(z)=\frac{H(z)+250000}{z}.
\end{split}
\tag{1}
\]

Its splitting field M/Q(s) is regular with group A10. Its five finite branch
inertia generators are double transpositions; the inertia at infinity has
cycle type (9,1). The input monodromy proof is indispensable: a square
discriminant by itself would not establish this group.

An alternative construction uses exactly three irreducible factors of total
rank 1075. Three is the minimum number of subgroup classes even when all
ramified local groups are included; the rank-301 construction uses a fourth
factor to lower the total degree.

The independent finite controls are in the
[checker](verify_alternating_subgroup_hasse.py) and its
[actual replay receipt](alternating_subgroup_hasse_verification.json).
The arithmetic argument below imports two precise primary statements,
recorded in the [source ledger](alternating_subgroup_hasse_source_ledger.json).
The finite checker does not certify these imported arithmetic theorems.

## Proposition 1. Subset stabilizers alone require every subset size

For G=A10 define

\[
H_k=(S_k\times S_{10-k})\cap G,\qquad 1\leq k\leq5.
\]

The conjugates of these five subgroups cover G. No proper subfamily does.
Their coset actions have total degree

\[
                 \sum_{k=1}^5\binom{10}{k}=637.
\tag{2}
\]

More generally, for n=2b the b subset sizes 1 through b are all necessary
when the allowed subgroups are the stabilizers of subsets of these sizes.

Proof. An even permutation of an even number n of letters cannot be an
n-cycle. It has at least two cycles, and hence preserves a nonempty proper
subset; pass to its complement to obtain a size at most n/2. Conversely the
even permutation of cycle type (k,n-k) has only the proper invariant subset
sizes k and n-k. It is covered by precisely the size-k family in the range
1 through n/2. This includes k=n/2. Each H_k has index binomial(n,k).

These are stabilizers of individual subsets. Replacing a middle-size subset
by an *unordered partition into two equal parts* gives a different subgroup.

## Theorem 2. An optimal three-class element cover of A10

Let W5 be the stabilizer in A10 of an unordered partition of the ten letters
into two blocks of five. Thus

\[
W_5=(S_5\wr S_2)\cap A_{10},\qquad [A_{10}:W_5]=126.
\]

The three subgroup classes

\[
                       H_1,\quad H_3,\quad W_5
\tag{3}
\]

cover every element of A10. No two conjugacy classes of proper subgroups
cover A10. The sum of the three coset degrees is

\[
                         10+120+126=256.
\tag{4}
\]

Proof of the cover. A permutation with a fixed point lies in a conjugate of
H1. A permutation with a 3-cycle preserves its support and lies in a
conjugate of H3. Suppose neither occurs. If any cycle is odd, the only
possible odd-cycle case on ten letters is (5,5): lengths 7 and 9 would leave
a 1- or 3-cycle, and a 5-cycle must be accompanied by a second 5-cycle.
Such a permutation preserves the corresponding unordered equipartition.
If all cycles have even length, alternate two colours around every cycle.
There are five letters of each colour and the permutation interchanges the
two colour classes. It again fixes an unordered equipartition.

A10 is transitive on k-subsets and on unordered 5+5 partitions: the
corresponding S10 stabilizer contains an odd permutation, so the S10 orbit
does not split on restriction to A10. This also explains why splitting of
some S10 cycle types into two A10 classes causes no exception to the proof.

The lower bound uses the inspected primary classification:
Bubboloni, *Coverings of the symmetric and alternating groups*,
[Theorem 3.9, printed page 13](sources/bubboloni-coverings-1009.3866v1.pdf).
It states that An is covered by conjugates of two proper subgroups exactly
for 4<=n<=8. We import this theorem; our finite checker proves the displayed
upper bound and does not reprove the full classification.

The elementary proof here concerns the number of subgroup classes.
Proposition 4.2 additionally proves minimum total degrees with an explicitly
identified finite subgroup-classification dependency.

## Proposition 3. Element covers can fail at tame ramified primes

Use letters 0 through 9 and set

\[
d=(0\,1)(2\,3),\qquad
h=(0\,1)(4\,5\,6\,7\,8\,9),\qquad
D=\langle d,h\rangle.
\tag{5}
\]

Then D is an even permutation group isomorphic to C2 x C6. It has normal
double-transposition inertia I=<d> and cyclic quotient D/I, but is contained
in no conjugate of any subgroup in (3).

Proof. The generators commute, have orders two and six, and generate a
group of order twelve. Its orbits have sizes 2,2,6, so it has no invariant
subset of size one or three. For an invariant unordered 5+5 partition,
both d and h would have to preserve each block: each generator has fixed
letters, while a permutation exchanging the two blocks has none. Each
block would therefore be a union of D-orbits, impossible for size five.
This proves the claim.

The local group pattern itself is possible at an odd prime: take the
compositum of a ramified quadratic extension and an unramified extension
of degree six. Its decomposition group is C2 x C6, with the indicated
quadratic inertia. This does not claim that a particular specialization of
(1) realizes (5). It establishes that a proof based only on the three-class
element cover leaves a genuine ramified-group gap.

## Theorem 4. Local covers with three or four subgroup classes

The four classes

\[
                         H_1,\ H_2,\ H_3,\ W_5
\tag{6}
\]

cover every cyclic subgroup of A10 and every subgroup D<=A10 having a
normal subgroup I generated by a double transposition with D/I cyclic.
Their coset degrees sum to

\[
                         10+45+120+126=301.
\tag{7}
\]

Every one of the four classes is necessary within this displayed family
and this domain of local groups.

There is also a cover of the same local-group domain using exactly three
classes:

\[
                         H_1,\ H_3,\ W_2,
\tag{8}
\]

where W2=(S2 wr S5) cap A10 stabilizes a partition into five pairs, or
equivalently a perfect matching. Its index is 945, so the total degree is
10+120+945=1075. Three classes is optimal among arbitrary proper subgroups,
by the lower bound in Theorem 2. The degree 301 of (6) is smaller. Proposition 4.2 establishes both minimum
total degrees using a separate finite subgroup-classification input.

Proof of the four-class cover. Cyclic groups are covered by Theorem 2. For the other case, conjugate
I so that it is generated by d=(0 1)(2 3). Since Aut(C2) is trivial, I is
central in D. Lift a generator of D/I to h in D. Then D=<d,h> and h
commutes with d.

The moved set of d has four letters and is D-invariant. Its D-orbits are
either two pairs or a single orbit of size four. In the first case one
pair is D-invariant, giving a conjugate of H2. In the second case h
interchanges the two pairs of d. On the remaining six letters, D-orbits
are simply cycles of h. If one such cycle has length one or three, use H1
or H3. Otherwise all those cycles have even length: an odd length five
would force a fixed letter.

Colour the two d-pairs oppositely on the four-letter support. On each of
the even h-cycles on the remaining six letters, alternate colours.
This gives five letters of each colour. The element d preserves both
blocks and h exchanges them, so D fixes the unordered equipartition.
This gives W5 and proves the cover.

Deleting H1 loses the cyclic subgroup of type (9,1); deleting H3 loses
type (7,3); deleting W5 loses type (5,5); deleting H2 loses (5).
Consequently the cover is inclusion-minimal. The three-class alternative
below uses a larger total degree.

Proof of the three-class cover. First consider an individual element. If
it has a 1-cycle or 3-cycle use H1 or H3. Otherwise its cycles are all even,
or it has type (5,5), as in Theorem 2. On each even cycle pair opposite
letters. For type (5,5), pair corresponding letters in the two cycles.
This gives an invariant perfect matching.

For a local group D=<d,h> with normal double-transposition inertia, inspect
h on the six letters fixed by d. A cycle of length one or three gives H1
or H3. Otherwise all those cycles are even. Pair opposite letters in each
such cycle, and use the two transposition pairs of d on its four-letter
support. Since h commutes with d, it preserves that pair partition of the
support. The resulting perfect matching is invariant under both
generators, so D lies in a conjugate of W2. A10 is transitive on perfect
matchings for the same stabilizer-parity reason used in Theorem 2.
This proves the three-class local cover.

The checker exhausts the centralizer of d: C_S4(d) has eight elements,
and after combining with S6 and imposing even sign there are 2,880 possible
h. It tests (6) for every pair <d,h>. Exactly 450 of these generator pairs
fail (3), and each is repaired by H2. The count is of lifts h, not of
distinct subgroups D. It also generates all 945 perfect matchings and checks
the alternative three-class cover for all 22 even cycle types and all 2,880
centralizer lifts. In particular W2 covers the local counterexample (5).

## Proposition 4.1. A uniform local cover in every even degree

Let n=2b>=8. In An take the odd subset-size classes H_k with 1<=k<b,
together with W2=(S2 wr S_b) cap An, the stabilizer of a perfect matching.
These floor(b/2)+1 classes cover all cyclic subgroups and every subgroup
having normal double-transposition inertia and cyclic quotient.

Proof for cyclic groups. If a cycle has odd length less than b, its support
gives one of the listed H_k. If an odd cycle has length greater than b,
its complement has odd size less than b and does the same. The only case
with odd cycles and no such invariant subset is (b,b), with b odd: if
one b-cycle's complement had more than one cycle, one of its odd cycles
would have size less than b. For (b,b), pair corresponding points of the
two cycles. If all cycles have even length, pair opposite points in each
cycle. Both constructions give a fixed perfect matching.

For a double-transposition local group write D=<d,h> with h centralizing
d as before. If an odd h-cycle on the n-4 fixed letters of d has length
less than b, use its support. If it has length at least b, the complement
inside those n-4 letters has positive odd size less than b and is
D-invariant. If there is no such odd cycle, pair opposite points on each
even cycle, and use d's own two transposition pairs on its support.
The matching is D-invariant. This proves the uniform claim. The finite
checker also replays the partition consequences for every even degree
from 8 through 30; these finite controls supplement the written all-degree
proof.

Consequently, any rational function f in Q(z) of degree n with regular
An monodromy and only double-transposition finite inertia admits the same
conditional Hasse-specialization argument as Theorem 5. Choose a rational
source point away from its poles and critical points. One uses that source point to obtain a simple rational root
at a rational parameter, chooses a rational arithmetic progression within
the finitely many required Hensel neighbourhoods, and keeps its parameters
integral outside an enlarged exceptional prime set. The actual degree-ten
cover supplies these hypotheses explicitly. Existence of rational covers
in every arithmetically admissible degree is a different moonshot.

## Proposition 4.2. Minimum total degrees for the entire local-group domain

Using the classified maximal-subgroup list of A10, the following minima
hold among families of conjugacy classes of arbitrary proper subgroups:

* To cover all cyclic subgroups, the minimum sum of coset degrees is 256.
* To cover all cyclic groups and all abstract double-transposition-inertia,
  cyclic-quotient groups, the minimum sum is 301.
* With exactly three subgroup classes covering the latter domain, the
  minimum sum is 1075.

These minima are attained by (3), (6), and (8), respectively. This is a
classification-backed finite computation, not a classification-free proof.

The complete maximal-subgroup class list has indices

\[
                   10,\ 45,\ 120,\ 126,\ 210,\ 945,\ 2520.
\]

They are H1,H2,H3,W5,H4,W2 and a primitive subgroup of order 720.
The [Sage/GAP checker](verify_alternating_hasse_minimal_degree.py) and
[actual receipt](alternating_hasse_minimal_degree_verification.json) compute
all seven classes and their exact coset fixed-point counts for four
witnesses. The relevant incidence table is:

| Maximal class | Index | (9,1) | (7,3) | (5,5) | D from (5) |
| --- | ---: | ---: | ---: | ---: | ---: |
| H1 | 10 | 1 | 0 | 0 | 0 |
| H2 | 45 | 0 | 0 | 0 | 2 |
| H3 | 120 | 0 | 1 | 0 | 0 |
| W5 | 126 | 0 | 0 | 1 | 0 |
| H4 | 210 | 0 | 0 | 0 | 1 |
| W2 | 945 | 0 | 0 | 5 | 1 |
| primitive order 720 | 2520 | 0 | 0 | 5 | 0 |

A positive entry means that the whole witness group fixes a coset.
Every proper subgroup lies in a maximal subgroup. Replacing it by such a
maximal subgroup preserves fixed-point coverage and decreases its index;
repeated maximal classes can then be merged. It therefore suffices to
optimize over the seven listed classes.

The first two witness types force H1 and H3, of total degree 130. The
least expensive class covering (5,5) is W5, giving 256 for cyclic groups.
The extra witness D requires H2, H4 or W2. The cheapest completion is H2,
giving 301, while any choice using W2 already has degree at least 1075.
For exactly three classes, the third class must cover both (5,5) and D;
only W2 does, giving 1075. The elementary covers prove matching upper
bounds for the entire respective domains. The checker exhausts all
2^7-1 possible nonempty maximal-class subfamilies as a second control.

The replay uses Sage 10.9, GAP 4.15.1 and PrimGrp 4.0.1. We inspected the
[retained GAP source](sources/gap-4.15.1-gpprmsya.gi),
MaximalSubgroupsSymmAlt beginning at line 2196. Its ANMAXPRIMS[10]=[6]
entry supplies the primitive maximal class. We also inspected the
[PrimGrp 4.0.1 manual](sources/primgrp-4.0.1-manual.pdf), Section 1.1,
pages 4-5, and the nine primitive groups of degree ten in the library.
The proper groups in that list contained in A10 have orders 60,360,720.
Completeness of the classified library is imported; it is not established
by checking its group orders. The receipt records this dependency, source
hash, generators, exact action counts and all software versions.

The minimum 301 concerns a cover of **all abstract groups in the stated
local domain**. An actual infinite family of number fields could have
fewer kinds of ramified decomposition groups. This theorem does not rule
out a rank-256 Hasse family obtained by ensuring that all its actual
ramified groups fit (3).

## Theorem 5. Infinitely many rank-301 Hasse counterexamples

There are infinitely many pairwise linearly disjoint A10 extensions
M_j/Q arising as splitting fields of specializations of (1), such that

\[
A_j=M_j^{H_1}\times M_j^{H_2}\times M_j^{H_3}\times M_j^{W_5}
\tag{9}
\]

satisfies

\[
 \dim_{\mathbb Q}A_j=301,\qquad
 \operatorname{Spec}(A_j)(\mathbb Q)=\varnothing,\qquad
 \operatorname{Spec}(A_j)(\mathbb Q_v)\ne\varnothing
 \quad\text{for every place }v.
\tag{10}
\]

The factors have degrees 10,45,120,126, respectively. Their normal closures
are M_j, and their product is a finite étale algebra. This is a uniform
specialization construction, although it does not print degree-120 and
degree-126 resolvent coefficients.

Proof.

**Local fixed-point criterion.** For a finite Galois extension M/Q with
group G and a subgroup H, the embeddings of M^H into an algebraic closure
form the G-set G/H. A Q_v-point on Spec(M^H) is equivalent to a fixed
coset for a decomposition group D_v, and hence to D_v being contained in
a conjugate of H. A point on the spectrum of a product chooses one
factor. Thus (10) is exactly the decomposition-group covering condition
for (6). Since each factor in (9) has degree greater than one, there is
no Q-point.

**Good-prime ramification.** Apply Conrad,
*Inertia groups and fibers*,
[Theorem 0.2, page 2](sources/conrad-inertia-groups-fibers-1999.pdf),
to the regular Galois cover of P1 associated to M/Q(s).
Outside a finite set S of rational primes, specialization inertia is
trivial or conjugate to a subgroup of a power of a geometric branch
inertia group. Enlarge S to contain the primes at which finite branch
coordinates fail to be integral, any collisions with the infinity section,
and primes dividing |A10|. For integral s, specialization cannot meet
the infinity section outside S. Consequently every nontrivial inertia
group outside S is generated by a double transposition; the order-nine
inertia at infinity never enters this argument. The residue-field quotient
of a nonarchimedean decomposition group by inertia is cyclic. Theorem 4
therefore handles all primes outside S, ramified or unramified.

**All exceptional primes at once.** At z=1 one has the exact values

\[
 s_* = f(1)=1054656,\qquad
 P_{s_*}(1)=0,\qquad P_{s_*}'(1)=f'(1)=-38112.
\tag{11}
\]

For each p in S put e_p=v_p(38112), and choose

\[
 M_0=\prod_{p\in S}p^{\,2e_p+1},\qquad s=s_*+M_0n,\quad n\in\mathbb Z.
\tag{12}
\]

The product is one if S is empty. For every p in S,

\[
P_s(1)=-M_0n,\qquad P_s'(1)=-38112-M_0n.
\]

If n is nonzero, then v_p(P_s(1))>2e_p and
v_p(P_s'(1))=e_p. The strict Hensel inequality gives a root of P_s in Q_p.
For n=0 the displayed exact root already suffices. The natural degree-ten
action of D_p therefore fixes a sheet: D_p lies in a conjugate of H1.
This works for every integer n, and so does not depend on an unsupported
density assertion for an arbitrary p-adic open set.

**Full global group in this progression.** Substitute (12) into the
regular A10 cover. The parameter change is affine and invertible over Q,
so its generic group is still A10. Integral Hilbert irreducibility gives
full A10 for infinitely many integers n. The precise imported formulation
is Zywina, *Hilbert's irreducibility theorem and the larger sieve*,
[Theorem 1.2, page 2](sources/zywina-quantitative-hit-2010.pdf), which records
the Cohen--Serre estimate O(B^(1/2) log B) for exceptional integer parameters
in one variable. Monic normalization P_s/9 has rational-function
coefficients, as allowed by that statement. Discard the finite branch
and presentation exceptions.

At a real place the decomposition group has order at most two and is
cyclic, so Theorem 2 handles it. Together with the preceding finite-prime
arguments and the local fixed-point criterion, this proves (10).

**Infinitely many disjoint splitting fields.** At an inductive step let
L/Q be the compositum of previously selected splitting fields. Regularity
gives M cap L(s)=Q(s), so M L/Q(s) has group
A10 x Gal(L/Q). Apply the same integral Hilbert theorem to a primitive
element of this compositum, keeping the progression (12), and exclude the
finite parameters where its chosen presentation or specialization fails.
A full product-group specialization has an A10 splitting field linearly
disjoint from L. This proves the stronger inductive assertion. In
particular the resulting algebras are pairwise nonisomorphic: the
degree-ten factor is distinguished by its degree, and its normal closure
determines M_j. Each H in (6) has trivial core because A10 is simple and H
is proper, so all the factors have that same normal closure.

For the same fields M_j, the product
M_j^{H1} x M_j^{H3} x M_j^{W2} has the same everywhere-local,
globally-pointless property, now with exactly three field factors and
rank 1075. The local proof uses the second cover in Theorem 4.

The finite set S, and hence M0, are obtained from the cover's finite
exceptional reduction data. No explicit numerical value for S or M0 has
been computed here. Their existence and effective determination are part
of the imported inertia framework. The proof establishes the infinite
family without an explicit list of admissible integer parameters.

## Corollary 6. Whole Keller fibres with a Hasse failure

For every A_j in Theorem 5, the retained
[finite-étale polynomial realization theorem](../../unpublished/common-arithmetic-fibers/sections/01-keller-fibers.tex),
label cor:finite-etale-realization, produces a polynomial Keller map

\[
                   F_j:\mathbb A^3_{\mathbb Q}\longrightarrow
                         \mathbb A^3_{\mathbb Q}
\]

and y_j in A3(Q) whose entire scheme-theoretic fibre is Spec(A_j).
Its generic degree is 301, its coordinate degrees are at most

\[
                              6\cdot301+2=1808,
\]

and its determinant is one. The retained
[arithmetic transfer theorem](../../unpublished/common-arithmetic-fibers/sections/03-realization.tex),
label cor:transfer, transfers exactly (10). By the separately retained
[quadratic-tilt theorem](../jacobian/universal-keller-structure/PROOF.md),
UAQTILT1, these compiler maps have generic arithmetic and geometric S301.
The A10 group belongs to the splitting field of the distinguished finite
algebra.

This is an application of the existing compiler theorem and its stated
assurance; it is not a fresh independent verification of that theorem.
The maps in this corollary may depend on j. The [fixed-family
theorem](../jacobian/monodromy-theory/FIXED_KELLER_MAP_FINITE_ETALE_FAMILIES.md)
instead uses one fixed four-dimensional polynomial Keller map for infinitely
many of these complete product fibres. A fixed three-dimensional realization
remains open.

## Priority, scope, and remaining smallest questions

The general existence of intersective polynomials with a given nonsolvable
realizable group is established: Sonn,
*Two remarks on the inverse Galois problem for intersective polynomials*,
[Proposition, printed page 436](sources/sonn-intersective-two-remarks-2009.pdf).
Thus broad Hasse existence for A10 is not a new inverse-Galois theorem.
The concrete additions here are the chosen regular pencil, its
double-transposition local groups, the obstruction (5), the economical
local covers, their precise finite degree optimization, and the uniform
arithmetic-progression construction.
The real-place condition is included explicitly; finite-prime
intersectivity alone would not state it.

The [zero-branch ramification theorem](ALTERNATING_PENCIL_RAMIFICATION_DENSITY.md)
shows that locally soluble rank-256 integer parameters have density zero
in every arithmetic progression; a thinner infinite family remains possible.
One fixed A3 Keller map for the complete product fibres also remains open,
while the fixed A4 realization is proved. Three local classes are sufficient and optimal in
number. The total-degree minima 301 and 1075 are established for the
whole abstract local domain, with the explicit classification dependency
in Proposition 4.2.

The elementary cover claims are replayed exactly with the Python
standard library. The receipt checks 22 even cycle types, all 945 perfect matchings, and all
2,880 even centralizer lifts, as well as the essential-component negative
controls and numerical ranks. A separate Sage/GAP receipt establishes
the classification-backed minimum-degree statement. It is an internal replay, not external
human review or formal verification. The infinite existence statement
additionally depends on the written proof and the explicitly imported
inertia and integral Hilbert theorems.
