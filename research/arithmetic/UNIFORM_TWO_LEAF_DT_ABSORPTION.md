# Every two-fresh-sheet DT extension of the canonical component is absorbed

For every \(r\geq3\), append a double transposition with two fresh
leaf sheets to the canonical degree-\(2r+1\) polynomial tuple.
Its two old anchors may be any distinct sheets. The extension
always reaches the canonical tuple in degree \(2r+3\).

The proof retains the old cyclic distance between the anchors,
chooses a useful representative for each distance, and isolates
one of six fixed thirteen-sheet cores. Every outside attachment
pair is checked against the actual lower antipodal pairing of the
core. This supplies the extension control needed after a removable
DT extraction; it does not prove that arbitrary tuples admit such
an extraction.

The threshold is sharp. At \(r=2\), all ten extensions still have
full \(A_7\), but the five adjacent-anchor choices belong to the
other full-\(A_7\) orbit and have no wreath-leaf representative.

## The statement and metric

Let
\[
 \begin{split}
 W_r&=(A,D_{r-1},\ldots,D_1),\\
 A&=(x\,a_1)(01),\quad
 D_j=(0\,a_j)(1\,b_j),\\
 x&=2r,\quad a_j=2j,\quad b_j=2j+1 .
 \end{split}                                                \tag{1}
\]
This is the canonical tuple of the
[wreath normal-form theorem](IMPRIMITIVE_LEAF_DELETION_NORMAL_FORM.md).
Its product is
\[
             P=(0\,x\,a_1\,a_2\,\cdots\,a_{r-1}\,1\,
                           b_1\,b_2\,\cdots\,b_{r-1}).       \tag{2}
\]
Add two fresh sheets \(u,v\), and put
\[
                             \eta=(a\,u)(b\,v),
                    \qquad a\ne b\text{ old sheets}.        \tag{3}
\]
We count ordinary signed generators
\(Q_i(c,d)=(d,d^{-1}cd)\). Common sheet relabeling is free.

**Theorem 1.** For every \(r\geq3\) and every pair (3),
\((W_r,\eta)\) reaches \(W_{r+1}\), up to common sheet relabeling.
A uniform ordinary-word bound is
\[
                     M+r^3+r^2+5r+1,\qquad M=446265984.      \tag{4}
\]
Thus the old-pair choices of this two-fresh-sheet extension form
one absolute orbit over the canonical component.

The fixed finite core-controller words used below exist by retained
complete finite orbit theorems. They have not been recovered or
replayed. The checker replays the explicit cyclic transport,
core extraction and final wreath normalization, and substitutes
only those certified controller endpoints.

## 1. Cyclic distance and six representatives

A full twist on the old \(r\) factors conjugates them by \(P^{-1}\)
and has length \(r(r-1)\). If \(P^k(a_0)=a\), use the inverse
\(k\)-th power of that partial full twist and then relabel the old
sheets by \(P^{-k}\), fixing \(u,v\). The old tuple returns to
\(W_r\), and its two anchors are cyclically translated by \(P^{-k}\).
The new factor was untouched by the ordinary word.

After interchanging the two anchor/fresh-sheet roles if necessary,
their directed old-cycle distance is \(d\in\{1,\ldots,r\}\).
Any pair at that distance is cyclically translated to any other.
Choose \(-r\leq k\leq r\), so the transport costs at most
\[
                               r^2(r-1).                   \tag{5}
\]
This translation preserves \(d\); transitivity on individual old
sheets alone would not identify the pair choices.

For \(r\geq6\), use these representatives:

| Case | Distance | Ordered anchors \((a,b)\) | Third selected satellite index |
| --- | --- | --- | --- |
| G | \(1\leq d\leq r-5\) | \((a_j,a_{r-2})\), \(j=r-d-2\) | \(j\) |
| E1 | \(r-4\) | \((a_1,a_{r-3})\) | \(r-3\) |
| E2 | \(r-3\) | \((x,a_{r-3})\) | \(r-3\) |
| E3 | \(r-2\) | \((0,a_{r-3})\) | \(r-3\) |
| E4 | \(r-1\) | \((0,a_{r-2})\) | \(r-3\) |
| E5 | \(r\) | \((0,a_{r-1})\) | \(r-3\) |

The displayed distances follow directly from (2).
In G, \(3\leq j\leq r-3\). Retain the six factors
\[
                     A,D_{r-1},D_{r-2},D_j,D_1,\eta.        \tag{6}
\]
For the exceptional cases use \(j=r-3\).
Move these factors to the left in their original relative order,
using positive Hurwitz moves. Each selected factor remains unchanged;
an outside factor that it crosses is conjugated.

There are \(r-5\) outside factors. Each crosses at most five selected
factors, so the cost is at most \(5(r-5)\).
The selected factors have thirteen distinct sheets.
Their old first five factors are a relabeling of \(W_5\);
the new DT adds two leaves. Hence the core product is a 13-cycle.
Its group is full \(A_{13}\): the old \(A_{11}\) supplies a 3-cycle,
the new product is transitive, and the odd-degree DT block argument
proves primitivity. Jordan's consequence in
[Jones, Theorem1.1](sources/jones-primitive-cycle-1209.5169v1.pdf)
then applies.

## 2. The actual outside anchoring pairs

Relabel only the thirteen core sheets by
\[
 \begin{array}{c|rrrrrrrrrrrrr}
 \text{core label}&0&1&2&3&4&5&6&7&8&9&10&11&12\\ \hline
 \text{sheet}&0&1&a_1&b_1&a_{r-1}&b_{r-1}&a_{r-2}&b_{r-2}&
 a_j&b_j&x&u&v .
 \end{array}                                                \tag{7}
\]
The common first five core factors are
\[
              (01)(2\,10),\ (04)(15),\ (06)(17),\
                           (08)(19),\ (02)(13).            \tag{8}
\]
The sixth factor and required outside pairs are:

| Case | Sixth core factor | Outside core-anchoring pairs | Chosen leaf |
| --- | --- | --- | --- |
| G | \((8\,11)(6\,12)\) | \(\{11,9\}\), \(\{2,3\}\) | 8 |
| E1 | \((2\,11)(8\,12)\) | \(\{11,3\}\) | 0 |
| E2 | \((10\,11)(8\,12)\) | \(\{2,3\}\) | 0 |
| E3 | \((0\,11)(8\,12)\) | \(\{2,3\}\) | 0 |
| E4 | \((0\,11)(6\,12)\) | \(\{2,3\}\) | 0 |
| E5 | \((0\,11)(4\,12)\) | \(\{2,3\}\) | 0 |

Here is the derivation, including the effect of the moves.
In G, an outside factor before \(D_j\) is first crossed by \(D_j\),
so its old root pair \(\{0,1\}\) becomes \(\{a_j,b_j\}\).
The later \(D_1\) fixes that pair, and \(\eta\) sends it to
\(\{u,b_j\}\). An outside factor after \(D_j\) is first crossed
by \(D_1\), giving \(\{a_1,b_1\}\), which \(\eta\) fixes.
The two outside ranges have sizes \(d-1\) and \(j-2\).
An empty range causes no problem.

In E1-E5 the three selected satellites are consecutive at the
start of the old star. Every outside factor is therefore first
crossed by \(D_1\). The final crossing by \(\eta\) turns
\(\{a_1,b_1\}\) into \(\{u,b_1\}\) in E1 and leaves it fixed
in E2-E5. All the outside factor's other two sheets belong to
its own untouched old satellite block. Thus it remains a block
swap between the listed core pair and an outside pair.
This is genuine conjugation of old factors, not an inverse
bringing move that leaves the old odd subgroup unchanged.

For direct verification, the six core products are:
\[
 \begin{array}{ll}
 G:&(0\,10\,2\,8\,11\,6\,12\,4\,1\,3\,9\,7\,5),\\
 E1:&(0\,10\,2\,11\,8\,12\,6\,4\,1\,3\,9\,7\,5),\\
 E2:&(0\,10\,11\,2\,8\,12\,6\,4\,1\,3\,9\,7\,5),\\
 E3:&(0\,11\,10\,2\,8\,12\,6\,4\,1\,3\,9\,7\,5),\\
 E4:&(0\,11\,10\,2\,8\,6\,12\,4\,1\,3\,9\,7\,5),\\
 E5:&(0\,11\,10\,2\,8\,6\,4\,12\,1\,3\,9\,7\,5).
 \end{array}                                                \tag{9}
\]
Delete the chosen leaf in the table. Every required outside pair
is then antipodal in the resulting 12-cycle. In G, sheet8 is
the unique choice satisfying both pair requirements.
Thus a wrong core wreath leaf would not suffice.

## 3. Finite core controllers and the whole wreath deletion

The retained
[complete polynomial census](DEGREE_FOURTEEN_HURWITZ_ORBIT.md)
has one absolute DT orbit in each of the following degrees:

| Degree | Independent literal total | Complete absolute total | Literal controller bound |
| --- | ---: | ---: | ---: |
| 9 | 59049 | 6561 | 6656 |
| 11 | 14992384 | 1362944 | 1363143 |
| 13 | 5801453125 | 446265625 | 446265984 |

For degrees9 and11, these are the explicit retained lower-degree
character totals and complete BFS rows of the
[exact receipt](degree_fourteen_hurwitz_verification.json),
rather than an assertion inferred from the degree13 theorem's title.
The cyclic action is free: a nonidentity power of an odd full
cycle has all its orbits of the same odd length greater than1,
whereas a DT's four-sheet support cannot be a union of such orbits.
Division of the literal total by the degree gives the absolute total.
The BFS from a full-alternating cyclic seed reaches that whole bound.

For degree \(n=2m+1\), an absolute path of length at most \(V-1\)
can be lifted literally with at most \(n-1\) old full twists,
each of length \(m(m-1)\). Indeed all cyclic normalizations
along a path commute with Hurwitz moves and accumulate into
one final cyclic conjugation. This gives the last column.

For each of the six cores, take a conjugate of \(W_6\) with
the same product and chosen leaf from the table. Its lower
pairing is precisely the antipodal pairing of (9) with that
leaf omitted. The degree13 finite theorem supplies a literal
core controller of length at most \(M\).
It acts only in the first six factor positions, leaving every
outside factor unchanged.

The leaf is outside every outside anchoring pair, so it is a
leaf of the whole tuple. Use the six lower core pairs together
with each untouched outside two-sheet pair as a partition of
all remaining sheets. The core preserves its six pairs. Each
outside DT swaps one core pair with its own pair. These new
blocks are leaves on the connected lower core block tree.
Consequently the lower group is the full wreath group, either
by the
[deletion dichotomy](IMPRIMITIVE_LEAF_DELETION_NORMAL_FORM.md)
or directly by the connected swaps and a block flip.

Apply the constructive wreath normal form to the whole tuple.
For \(m=r+1\) factors, it uses at most
\[
                         2(m-1)(m-2)+2m-1=2r^2+1          \tag{10}
\]
ordinary moves. To check this bound, the recursive transposition
star word on \(m\) blocks costs at most \((m-1)(m-2)\):
at each deletion it brings one input and one target edge to the
end, costing at most \(2(m-2)\). The normal-form proof uses
this word twice, at most \(m-1\) moves to bring the leaf owner
first, at most \(m-2\) moves for the root adjustment, and two
moves for its square. The final relabeling is free.

Adding (5), the extraction cost, the controller and (10) is
at most the bound (4). This proves Theorem1 for \(r\geq6\).
For \(r=3,4,5\), the whole extension already has degree9,11,13.
Its expanded graph is a tree, so its product is a full cycle;
the old alternating subgroup supplies a 3-cycle and the
odd-degree argument gives full alternating group.
The corresponding finite row above handles it directly.
All three bounds are at most \(M\). \(\square\)

## The sharp degree-seven boundary

**Theorem 2.** Let \(r=2\), so \(W_2\) has group D5 on five
old sheets. All ten unordered old-pair extensions (3) have
full \(A_7\). Exactly the five pairs of old cyclic distance2
reach \(W_3\); the five pairs of distance1 have no wreath-leaf
representative.

**Proof.** The [checker](verify_uniform_two_leaf_dt_absorption.py)
enumerates the complete signed Hurwitz orbit of each of the
ten inputs. Each has147 literal tuples and588 directed signed
edges. Exact generated-group closure has order2520 in every case.

For distance2 it retains and replays a literal word to a
fixed-product conjugate of \(W_3\). For example the anchors
\(\{0,1\}\) admit \((1,-2,-2,1,-2)\) with the stated conventions.
For distance1, every leaf deletion at every one of the147
tuples is checked and is symmetric, never wreath.
Normalizing all products to the same 7-cycle gives exactly
two disjoint147-element orbits, one for each distance class.

This is a complete finite certificate of nonreachability to
any wreath leaf, with no fixed chosen-leaf restriction.
It also certifies full \(A_7\) in the negative class, so the
failure is an orbit distinction rather than a smaller generated
group. Old cyclic translations act transitively on each of
the two distance classes. \(\square\)

## Corollary 3. Recursive extensions and the buffered families

The canonical component in every degree \(2r+1\), \(r\geq3\),
is closed under adjoining a DT with any two distinct old anchors
and two fresh sheets. This also holds when the old tuple is any
representative of that component: apply its old normalizing
word while leaving the new factor untouched, extend the final
old relabeling to fix the two fresh sheets, and use Theorem1.

In particular all previously displayed odd-factor fresh-buffered
D5 bouquets lie in the canonical wreath component. Normalize
their old unbuffered prefix using the
[root-preserving bouquet theorem](D5_LOBE_BOUNDED_CORE_ABSORPTION.md).
Its final relabeling carries the buffer's old anchors to some
distinct pair; Theorem1 handles that pair. The same argument
applies after normalizing an
[ordered one-vertex D5 lobe tree](NATURAL_D5_LOBE_TREE_HURWITZ_REACHABILITY.md).

The extension theorem itself does not extract a factor from an arbitrary
tuple. That input is supplied by
[Theorem 20 of the terminal-mass note](TERMINAL_MASS_HURWITZ_DESCENT.md#theorem-20-every-residual-pure-dt-orbit-admits-extraction):
every pure-DT orbit on odd \(N\ge15\) sheets contains a removable factor.
Induct from the complete degree-thirteen pure-DT base: move the removable
factor to the end by positive moves, delete its two exclusive leaves,
normalize the lower tuple, then apply Theorem 1 at its two old anchors.
Thus all residual pure-DT tuples lie in the canonical component. The separate
mixed polynomial connectivity problem, and hence the original degree-sixteen
and all-degree even-passport conclusions, remain open.

## Exact evidence and limits

The [checker](verify_uniform_two_leaf_dt_absorption.py) verifies
all six core products, chosen-leaf pairings and their conjugated
\(W_6\) targets. It replays every normalized distance case for
\(6\leq r\leq12\), selected cases at \(r=25\), and all490 old
unordered anchor pairs for \(6\leq r\leq9\), including the
signed partial full twists and fresh-sheet interchange.
It checks all112 small-degree extension endpoints for \(r=3,4,5\).

Every positive whole-tuple wreath endpoint is actually normalized
by the frozen normal-form implementation, with bound (10) checked.
A deliberately wrong generic core leaf produces a symmetric
global deletion and is rejected. The complete degree-seven
boundary has both positive word replays and exhaustive negative
leaf checks.

The [receipt](uniform_two_leaf_dt_absorption_verification.json)
and [source ledger](uniform_two_leaf_dt_absorption_source_ledger.json)
distinguish those actual operations from substitutions of the
unrecovered finite-controller endpoints. Uniformity follows from
the distance representatives, exact conjugation formulas and
finite fixed-core proof, not the bounded control range.

There is no new large census, recovered finite core-controller
word, optimal distance claim, arbitrary-DT orbit classification,
formal verification, independent external replay, external human
review or novelty claim.
