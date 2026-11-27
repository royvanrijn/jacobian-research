# The complete two-heavy first deficit

Let the two non-transposition finite factors have nontrivial cycle
profiles
\[
 \mathbf a=(a_1,\ldots,a_u),\qquad
 \mathbf b=(b_1,\ldots,b_v),\qquad a_i,b_j\ge2,
\]
each of permutation index at least two. Put
\[
 p=u+v,\quad s=p-3\ge1,\quad
 A_0=\sum_i a_i,\quad B_0=\sum_j b_j,\quad N=A_0+B_0-2.
                                                               \tag{1}
\]
Consider the minimal full-cycle passport
\[
 (\,\mathbf a,1^{B_0-2}\,),\quad
 (\,\mathbf b,1^{A_0-2}\,),\quad T_1,\ldots,T_s.             \tag{2}
\]
For every \(p\ge5\), it has exactly one component, with actual group
\(S_N\). This holds under full braids and under pure braids on each
fixed ordered passport. The \(p=4\) cases are exactly the previously
classified two-heavy, one-simple sector: one primitive component,
and a second quadratic component only when both heavy factors
repeat their cycle length.
All separate sizes, block systems and deck groups are explicit in Z5.

The proof first creates a sheet fixed by both heavy factors. It uses
this sheet to reduce the heavy forest to three edges, then connects
all three-edge forests by leaf transfers. A pure loop independently
rotates each heavy component. This removes the phase left when a
one-heavy normalization is lifted through the two heavy factors.

## Sources and conventions

The closest general comparison is Wajnryb,
[Theorem 1, retained p2](sources/wajnryb-two-special-fibers-1995.pdf):
at most two nonsimple critical values are allowed. Infinity and
both heavy finite values in (2) give three, so that conclusion
cannot be imported. Jones--Zvonkin,
[Theorem 1.1 and Sections 1.3, 2.2, 5.4](sources/jones-zvonkin-braid-cacti-2002.pdf),
give the polynomial/cactus/braid correspondence, count conventions
and quotient-copy examples. These actual statements and their
limits were inspected. No novelty claim is made.

The tree facts are H1 of the
[polynomial passport note](POLYNOMIAL_PASSPORT_COMPONENTS.md).
The positive one-transposition base is
[Y1--Y5](POLYNOMIAL_TWO_CYCLE_PAIR_COMPONENTS.md).
Z4 uses the labeled-root version of the logarithmic-derivative
bridge F1 and its primary component input,
Chen--Gendron--Prado--Tahar,
[Theorem 1.7, printed p895](sources/chen-gendron-prado-tahar-isoresidual-curves-2026.pdf).
The bridge, nonresonance and removal of critical-value collisions
are audited in the
[single-cycle-tail note](POLYNOMIAL_SINGLE_CYCLE_TAIL_COMPONENTS.md).
Here the zeros of that differential all have order one, so neither
exceptional disconnected family applies. The primary theorem
labels both poles and zeros; this is essential in Z4.

Composition is right to left, with
\(Q_i(C,D)=(D,D^{-1}CD)\). Work first modulo simultaneous
sheet relabeling, letting the product range over its full-cycle
conjugacy class. Keep the two heavy strands in order and put
the transpositions last:
\[
                         (A,B,T_1,\ldots,T_s).
 \]
Intermediate braids may permute the transposition strands.
The final pure assertion is proved separately.

## Lemma Z1. A transposition hub with the heavy factors fixed

The hub assertion requires only minimality and a full-cycle product;
it applies with any positive number of transpositions.
Delete the transposition edges from the incidence tree.
Its \(s+1\) components are the orbits of \(\langle A,B\rangle\),
and \(R=AB\) is a full cycle on each. Choose any sheet \(r\)
in one component \(H\). There is an absolute braid representative
with \(A,B\) literal and
\[
                    T_j=(r,x_j),                          \tag{3}
\]
one \(x_j\) in each component other than \(H\).

**Proof.** The transposition edges form a tree after contracting
the heavy components. To increase the degree of \(H\), choose
an edge \(T\) not incident to \(H\). Gather it last by positive
moves through the other transpositions. Every crossed factor
is unchanged; only the selected transposition is conjugated.
Each crossed edge lies in one side of the original cut at \(T\),
so its two endpoints stay on their respective cut sides.

Let \(C\) be the product of the remaining prefix. Its restrictions
to the two cut trees are full cycles. A prefix full twist followed
by simultaneous relabeling by \(C\) restores every prefix factor
and conjugates \(T\) by \(C\). A power puts its endpoint on the
side containing \(H\) at \(r\); its other endpoint stays outside
that side. All other transpositions and both heavy permutations
are literal. The degree of \(H\) has increased by one.
Repeat until every transposition is incident to \(H\).

The contracted tree is now a star. Repeat the same gathering
and prefix rotation for each edge to put its \(H\)-endpoint
at \(r\). Its other cut side is just one heavy component, so
it stays in that component. Previously normalized edges are
unchanged. This proves (3). \(\square\)

## Lemma Z2. An outside sheet and edge deletion

Let \(e\) be the number of edges of the bipartite heavy-cycle
forest, and let \(z\) count sheets fixed by both heavy factors.
Then
\[
                            e=z+2.                        \tag{4}
\]
For \(p\ge5\), every orbit contains a representative with \(z\ge1\).
If \(z\ge1\), any chosen heavy-core edge can be deleted by a
braid, changing \((e,z)\) to \((e-1,z-1)\).

**Proof.** The heavy-cycle forest has \(p\) cycle vertices and
\(e\) edges. Including the \(z\) isolated sheets, its component
count is \(p-e+z=s+1=p-2\), proving (4).

If \(z=0\), then \(e=2\). Two edges use at most four cycle
vertices, so some heavy cycle is an isolated component.
Suppose it is an \(A\)-cycle; the other color is symmetric.
Use Z1 with a root outside this component to make it a leaf.
Gather its transposition last. Its endpoint in the isolated
cycle is always \(A\)-only. Rotate its other endpoint to a
\(B\)-only sheet. Such a sheet exists: \(B\)'s support has
size at least three, while only \(e=2\) sheets are shared.
Gather the selected transposition immediately after \(B\)
and use the literal identity
\[
 Q^2(B,(u,v))=((u,v)B(u,v),(u,B^{-1}u)),
       \quad B(u)\ne u,\quad B(v)=v.                       \tag{5}
\]
Here \(u\) was \(B\)-only and \(v\) was \(A\)-only.
The sheet \(u\) becomes outside both supports, and \(v\)
becomes shared. Thus \(z=1\).

For deletion, let \(q\) be the shared sheet of the chosen edge
and choose an outside sheet \(w\). In its heavy component
choose a root \(r\ne q\), and use Z1. Its \(w\)-edge is the
only factor moving \(w\), so \(w\) is an exclusive leaf.
Gather that factor last and use the complementary full cycle
to put its other endpoint at \(q\). Every other transposition
fixes both \(w,q\): its root is \(r\), and its other endpoint
is in another heavy component. It therefore commutes past
the selected factor when gathering it next to \(B\).
Apply (5). It replaces \(q\) by \(w\) in its \(B\)-cycle.
The old \(q\) remains in \(A\), the old outside sheet becomes
\(B\)-only, and the chosen core edge is removed.
No other heavy-core edge changes. \(\square\)

## Lemma Z3. Buffered leaf transfers normalize the heavy forest

Suppose \(z\ge1\). A leaf \(A\)-cycle can be detached from its
neighbor and reattached to any \(B\)-cycle having a free sheet
after detachment. Any such sheet can be chosen as the new
attachment. The whole factor \(A\) stays literal.
The symmetric statement holds for a leaf \(B\)-cycle.

Consequently, when \(u,v\ge2\), every orbit reaches a heavy forest
consisting of one four-vertex path and all other heavy cycles
isolated, with exactly one outside sheet. The path's two cycles
of each color, their internal/leaf roles and both oriented
attachment distances can be prescribed.

**Proof of the transfer.** Let \(L\) be the leaf \(A\)-cycle and
\(q\) its shared sheet. Choose the Z1 root in its heavy component
outside \(L\); the neighboring \(B\)-cycle has a sheet other
than \(q\), so this is possible. Every other transposition
then fixes \(L\). Take the exclusive outside leaf \(w\)
and move its transposition to \((w,q)\), as in Z2.
Capture (5) detaches \(L\), replacing \(q\) by \(w\) in its
old neighbor. The selected transposition now joins the pure
\(A\)-cycle \(L\) to the rest of the tree.

Gather it last, keeping every other factor literal. The two
complementary restrictions are full cycles. Rotate the endpoint
in the rest to any chosen \(B\)-only sheet \(y\) of the target
cycle. The other endpoint is some \(x\in L\), hence \(A\)-only.
Gather next to \(B\) and capture again. It replaces \(y\) by
\(x\), reattaching \(L\) and making \(y\) an outside sheet.
Every other heavy incidence is retained. The transposition
factors may change during gathering, which is permitted.
Thus the transfer preserves \(e,z\). Choosing the original
neighbor allows any freed slot there. Exchange the heavy
factors and undo the exchange to obtain the other-color version.

**Reduction of the graph.** Use Z2 to delete edges until \(e=3\).
Throughout, \(z=e-2\ge1\) at each completed step. A three-edge
bipartite forest has, besides isolated vertices, just one of:
three disjoint edges, a two-edge path and a disjoint edge,
a three-edge path, or a three-edge star.

A two-edge path plus an edge becomes a three-edge path by
transferring the appropriately colored leaf from the separate
edge to an endpoint of the two-edge path.
Three disjoint edges first become that configuration.
For a three-edge star, transfer one leaf to an isolated vertex
of the center's color; one exists since that color has at
least two cycles. This gives a two-edge path plus an edge.
These moves require only new degrees at most two, so every
cycle length at least two suffices.

Write the resulting path as \(A_1-B_1-A_2-B_2\).
Transferring \(A_1\) to \(B_2\) exchanges the two \(B\) roles;
transferring \(B_2\) to \(A_1\) exchanges the two \(A\) roles.
If \(A_3\) is an isolated cycle, the three transfers
\[
 B_2:A_2\to A_3,\qquad A_1:B_1\to B_2,\qquad
 B_1:A_2\to A_3                                    \tag{6}
\]
replace internal \(A_2\) by \(A_3\), leaving \(A_2\) isolated.
All target degrees are at most two. Together with role exchange,
this selects any desired two \(A\)-cycles. The other color is
symmetric. Finally detach and reattach to the same internal
neighbor at any free slot to prescribe each of the two oriented
distances, exactly as in Y3. The two internal cycles each have
only two attachment sheets. There is no other heavy-forest
cyclic-order parameter. \(\square\)

## Lemma Z4. Independent phases and labeled attachment normalization

Assume \(s\ge2\), and fix a heavy forest having exactly one outside sheet \(w\).
Modulo absolute braids, all transposition attachments and their
orders are equivalent for that same heavy-forest isomorphism type.

**Independent rotations.** First use Z1 with root \(w\), so that
each transposition is \((w,x_j)\), one for each other \(R\)-cycle
\(C_j\). Let \(L_j=|C_j|\). There is a pure braid loop fixing
every transposition and conjugating \(A,B\) by the single
component rotation \(R_j\), with identity on all other components.

Indeed gather the selected transposition first among them.
A full twist on the prefix \(A,B,T_j\) conjugates it by
\((RT_j)^{-1}\). Its power \(L_j+1\) is identity on
\(C_j\cup\{w\}\), and is \(R_k^{-(L_j+1)}\) on each other
component \(C_k\). In particular \(R,T_j\) are restored.
Undo gathering, which restores every transposition literally.
Now apply the inverse \((L_j+1)\)-st power of the two-heavy
full twist. It conjugates \(A,B\) by \(R^{L_j+1}\).
The combined conjugator is \(R_j\) on \(C_j\) and identity
elsewhere. The loop is pure, including the gathering conjugation.

**Attachment orders.** Treat \(R=AB\) as a single marked branch
factor, and label each of its cycles. A polynomial with that
labeled zero fiber has
\[
 f(z)=C\prod_j(z-x_j)^{|C_j|},\qquad
 \frac{df}{f}=\sum_j\frac{|C_j|\,dz}{z-x_j}.                 \tag{7}
\]
There are \(s\ge2\) labeled zeros of order one of \(df/f\),
corresponding to the transposition branches. Its labeled
simple poles have residues \((|C_j|)_j,-N\).
No proper subset sums to zero. The primary Theorem 1.7,
through the exact F1 bridge, says that this labeled
isoresidual fiber is connected: its zero orders are not all
even, and there are at least four simple poles.
F1 also justifies deleting critical-value collisions.
Thus all attachments and orders are connected even with
the \(R\)-cycles and the simple critical points marked.

Lift the resulting braid by replacing the \(R\)-strand by the
consecutive block \(A,B\). Crossing that block with a
transposition is realized by the two usual crossings:
\[
 (A,B,T)\longmapsto(T,TAT,TBT).
 \]
Throughout the lift, \(A,B\) are conjugated by the same
permutation. At its endpoint, after the one-heavy simultaneous
relabeling, this permutation centralizes \(R\) and fixes each
labeled \(R\)-cycle setwise. It is therefore a product of
independent component rotations.
Choose a target with all transpositions at the hub \(w\).
The pure loops just proved cancel those component rotations
while fixing every target transposition. This restores the
desired heavy forest as well as its attachments.
This last correction is necessary: connectivity after bundling
\(A,B\) alone would not establish connectivity of its refinements.
\(\square\)

## Theorem Z5. Complete two-heavy first-deficit classification

For (1)--(2), every \(p\ge5\) has exactly one full component and
one pure component on every fixed ordered passport, in both
literal and absolute conventions. Its actual group is \(S_N\),
its only block-count divisors are \(1,N\), and its deck group
is trivial.

For \(p=4\), Y4--Y5 give one primitive \(S_N\) component,
with a second quadratic component precisely when
\(u=v=2\), \(a_1=a_2\), \(b_1=b_2\). There are no other
components in the complete positive two-heavy first deficit.

For exact sizes, let \(m_{\mathbf a}(d)\) and
\(m_{\mathbf b}(d)\) be the nontrivial-length multiplicities.
The fixed-order total literal count is
\[
 H=N^{s+1}
 \frac{(u+B_0-3)!}{(B_0-2)!\prod_d m_{\mathbf a}(d)!}
 \frac{(v+A_0-3)!}{(A_0-2)!\prod_d m_{\mathbf b}(d)!}.        \tag{8}
\]
For \(p\ge5\), this is the single pure literal size, and the
pure absolute size is \(H/N\). Full sizes multiply these by
\[
 W=\frac{(s+2)(s+1)}
          {1+\boldsymbol1_{\mathbf a=\mathbf b}},           \tag{9}
\]
where profile equality means multiset equality.
For \(p=4\), all separate sizes and quadratic group data are Y4--Y5.

**Proof.** For \(p\ge5\) and \(u,v\ge2\), Z2 gives an outside sheet.
Z3 makes the same prescribed three-edge heavy forest in every
orbit, with one outside sheet. Z4 removes all attachment
phases and orders. This proves connectivity while keeping
the heavy branch strands in their original order.
If one heavy factor has a single nontrivial cycle, the
other factor together with that single cycle and the \(s\)
transpositions is F2's single-cycle-tail family. A
transposition tail gives one component. This also preserves
the two distinct heavy branch types.

Z1 with an outside root \(w\) also proves the actual group.
The heavy group is transitive on each of its components
and fixes \(w\). Conjugating the star transposition for that
component yields every \((w,x)\) in it. These transpositions
generate \(S_N\). For the single-cycle-factor case, Z2
still creates an outside sheet, so the same argument applies.

To pass to pure braids, first normalize to a hub
at an outside sheet. Any adjacent pair of distinct hub
transpositions shares exactly that sheet and satisfies
\(Q_i^3(T_i,T_{i+1})=(T_i,T_{i+1})\).
This literal stabilizer exchanges their two strand markings.
Its transpositions generate the entire symmetric group
on the \(s\) simple strands. Since the orbit was connected
with only those strand permutations allowed, the kernel,
the pure braid group on this fixed heavy order, is transitive.
Conjugating treats any fixed ordered passport.
The full twist is pure and realizes conjugation by the
full product, so the conclusion is literal as well as absolute.

Formula (8) is the corrected Goulden--Jackson minimal
cactus count: the total number of finite factors is \(s+2\),
each transposition contributes one, and the two heavy
factorials are displayed explicitly. The number of type
placements is (9). The symmetric group has the stated
blocks and deck, and its absolute/literal conversion is \(N\).
The \(p=4\) assertion is Y's already proved complete base.
\(\square\)

## Exact evidence and remaining scope

The [checker](verify_polynomial_two_heavy_first_deficit.py)
checks the hub construction, creation and deletion of outside
sheets, both colors of leaf transfer, the cabled crossing
identity, the independent component-rotation word and the
odd stabilizers on transposition markings. It separately
enumerates complete finite three-edge forest graphs and
complete pure-braid passport orbits against independent
character totals. Negative controls reject a wrong rotation
exponent and inadmissible forest moves.

The [receipt](polynomial_two_heavy_first_deficit_verification.json)
and [source ledger](polynomial_two_heavy_first_deficit_source_ledger.json)
give exact bounds, seeds, hashes and assurance limits.
The finite controls support the written all-profile argument;
they are not an extrapolation.

Run from the canonical checkout:

    python3 research/arithmetic/verify_polynomial_two_heavy_first_deficit.py

Every positive first-deficit passport with \(k=2\) is now
classified. The first-deficit problem for \(k\ge3\), deeper
deficits, zero-simple passports, larger noncyclic return
fibers and the full classification of genuine lifting
obstructions remain open outside already resolved families.
No geometric irreducibility of every auxiliary determinant
factor before its finite label quotient is asserted.
Independent mathematical review and replay remain pending.
