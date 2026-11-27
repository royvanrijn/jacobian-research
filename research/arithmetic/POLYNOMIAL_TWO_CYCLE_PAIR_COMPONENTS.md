# Two cycle pairs: complete components and pure braid normalization

Let \(a_1,a_2,b_1,b_2\ge2\), put
\[
 A_0=a_1+a_2,\quad B_0=b_1+b_2,\quad N=A_0+B_0-2,
\]
and fix an \(N\)-cycle \(P\). Consider the finite passport
\[
 (a_1,a_2,1^{B_0-2}),\qquad
 (b_1,b_2,1^{A_0-2}),\qquad T.                             \tag{1}
\]
Its index sum is \(N-1\).
There is exactly one primitive component, with group \(S_N\).
There is one additional component precisely when
\(a_1=a_2\) and \(b_1=b_2\); it is the known quadratic lift.
These statements hold for full braids and, on each fixed ordered
passport, for the pure braid group. All blocks, deck groups and
separate orbit sizes are given in Y4.

Two moves extend the
[duplicated-cycle normalization](POLYNOMIAL_DUPLICATED_CYCLE_PATH_NORMALIZATION.md).
A matching with different support types at its two transposition
endpoints can escape through either a path or a star. At a path,
reattaching a detached leaf cycle to the other cycle of the opposite
color changes which cycle is internal. This removes the role
obstruction created by unequal lengths.

## Sources, scope and conventions

Jones--Zvonkin,
[Theorem 1.1 and Sections 1.3, 2.2, 5.4](sources/jones-zvonkin-braid-cacti-2002.pdf),
supply the polynomial/cactus/braid correspondence and count conventions;
their actual constructions of quotient-copy examples do not supply
the uniform normalization here.
Wajnryb,
[Theorem 1, retained p2](sources/wajnryb-two-special-fibers-1995.pdf),
requires at most two nonsimple values. Infinity and both heavy finite
values in (1) are nonsimple, so that theorem does not apply.
These primary statements were inspected; no novelty claim is made.

The tree input is H1 of the
[polynomial passport note](POLYNOMIAL_PASSPORT_COMPONENTS.md).
The literal rotation and capture identities are X2 of the
[duplicated-cycle note](POLYNOMIAL_DUPLICATED_CYCLE_PATH_NORMALIZATION.md);
they require no equality of cycle lengths. The quadratic group's
classification and the corrected Goulden--Jackson count are U2 of the
[first-deficit note](POLYNOMIAL_FIRST_DEFICIT_COMPONENT_REDUCTION.md).
Y5 separately invokes F2 of the
[single-cycle-tail theorem](POLYNOMIAL_SINGLE_CYCLE_TAIL_COMPONENTS.md),
whose precise primary isoresidual theorem and label hypotheses
are recorded there.

Write the ordered tuple as \((A,B,T)\), with \(T\) last.
The two nontrivial cycles within each heavy factor are unordered.
Composition is right to left and
\(Q_i(C,D)=(D,D^{-1}CD)\).
Absolute equivalence allows simultaneous sheet relabeling and lets
the product range over the conjugacy class of \(P\).
For a fixed literal product it is the quotient by \(\langle P\rangle\).
Pure braids keep all three branch strands individually marked.
Other ordered placements are obtained by braid conjugation;
normality of \(PB_3\) transfers the statements to each one.

## Theorem Y1. General four-core description

Deleting the transposition cycle vertex from the cycle-incidence
tree leaves exactly two components. The permutation \(R=AB\)
is a full cycle on each component. Suppressing unshared sheet
vertices and replacing a shared sheet by an edge gives precisely:

1. two disjoint \(A\)-\(B\) edges;
2. an \(A\)-centered two-edge star and an isolated \(A\)-cycle;
3. a \(B\)-centered two-edge star and an isolated \(B\)-cycle;
4. a path \(A_\ell-B_i-A_i-B_\ell\) and one isolated sheet.

In a path, put
\[
 p=A_\ell\cap B_i,\quad q=A_i\cap B_i,\quad r=A_i\cap B_\ell,
 \qquad
 \alpha=\operatorname{dist}_{A_i}(q,r),\quad
 \beta=\operatorname{dist}_{B_i}(p,q).                       \tag{2}
\]
Once the internal/leaf assignment of the two lengths of each
color is chosen, the absolute heavy forest is determined by
\[
 1\le\alpha<|A_i|,\qquad 1\le\beta<|B_i|.                   \tag{3}
\]
Every choice occurs, and its transposition joins the isolated
sheet to any sheet of the other component.

**Proof.** There are four heavy cycle vertices. If their
bipartite core has \(e\) edges and there are \(z\) sheets
outside both heavy supports, its component count is
\(4-e+z=2\). The core is a forest in \(K_{2,2}\), so
\(e\le3\) and \(z\le1\). Its two-edge and three-edge forests
give the list. This repeats the incidence proof of X1
without its unnecessary equal-length specialization:
the incidence graph now has \(N+5\) vertices and
\(A_0+B_0+2=N+4\) edges. It is a tree by the minimal
index sum and full-cycle product.
The tree-product argument gives the two full restrictions of \(R\).
Each internal cycle has two distinguished sheets; their oriented
distance determines its cyclic order. A leaf cycle has one.
Gluing such cycles along the displayed path and adding one
fresh sheet with \(T\) proves existence and uniqueness.
\(\square\)

## Lemma Y2. Every nonaligned matching escapes

Every star reaches a path by pure braids in the absolute quotient.
A matching reaches a path unless its two components have the same
\(A\)-length, the same \(B\)-length, and aligned transposition
endpoints under their unique color-preserving cyclic identification.

**Proof.** We use X2's absolute synchronized rotation
\[
             (A,B,T)\longmapsto(A,B,R^kTR^{-k})              \tag{4}
\]
and its literal capture identity
\[
 Q_2^2(A,B,(u,v))=(A,(u,v)B(u,v),(u,B^{-1}u)),
       \quad B(u)\ne u,\ B(v)=v.                           \tag{5}
\]
The first is \(Q_1^{2k}\) followed by simultaneous relabeling
by \(R^k\). Both braid operations are pure.

In an \(A\)-centered star, rotate the main endpoint of \(T\)
to a \(B\)-only sheet. Its other endpoint is on the isolated
\(A\)-cycle and is \(A\)-only. Capture (5) gives a path.
For a \(B\)-centered star, exchange the heavy factors,
do this construction, and undo the exchange.
The combined braid is pure.

In each matching component, read the \(R\)-cyclic support word
\[
                   O\,B^{\,b-1}A^{\,a-1},                 \tag{6}
\]
where \(O\) denotes the shared sheet. Under (4), compare the
letters at the two transposition endpoints. Suppose a pair
of different letters occurs.

- For \(A,B\), (5) directly creates a path.
- For \(O,A\), (5) moves the shared \(B\)-attachment to the
  other \(A\)-cycle. This gives an \(A\)-centered star,
  which already reaches a path.
- For \(O,B\), first use \(Q_1^{-1}\), giving
  \((ABA^{-1},A,T)\). At the old shared sheet \(u\), the
  first factor now fixes \(u\): its predecessor in the
  \(A\)-cycle was \(A\)-only. The other endpoint was
  \(B\)-only, so it is still moved by the first factor
  and fixed by the second. The two endpoints are now
  exclusive and of opposite colors. Apply (5) to the
  second factor and \(T\), then undo the heavy exchange.
  This gives a path with the original branch order.
  The word \(Q_1^{-1}Q_2^2Q_1\) is pure.

If no different-letter pair ever occurs, the two doubly infinite
periodic words observed from the endpoints are identical.
Each component word has exactly one \(O\) per period, so its
least period is its component size. Equality forces equal sizes
and aligned \(O\) positions. The subsequent lengths of the
\(B\) and \(A\) runs must then agree. Thus both within-color
lengths are repeated and the two endpoints are aligned.
These matchings are the quadratic lifts already identified in U2.
No other matching is left. \(\square\)

## Theorem Y3. Pure path normalization with independent cycle roles

All path configurations for a fixed pair of unordered heavy profiles
are connected under \(PB_3\) in the absolute quotient.
They have actual monodromy \(S_N\).

**Proof.** Rotate \(T\)'s main endpoint to \(p\) and write
\(T=(p,v)\), where \(v\) is the isolated sheet. Applying
(5) replaces \(p\) by \(v\) in \(B_i\); it detaches the
leaf cycle \(A_\ell\), leaving an \(A\)-centered star.
The factor \(A\), including \(q,r\), is literal.

Choose either \(B\)-cycle of this star. It has one shared
sheet with \(A_i\), denoted \(s\), and all its other sheets
are \(B\)-only. Rotate the main endpoint of the new
transposition to any such sheet \(y\). The other endpoint
is some sheet \(x\) of the isolated \(A_\ell\)-cycle.
Capture replaces \(y\) by \(x\) in the chosen \(B\)-cycle,
reattaching \(A_\ell\) and making \(y\) isolated.

Choosing the old \(B_i\)-cycle preserves both roles and
\(\alpha\), while prescribing any new \(\beta\) in
\(1,\ldots,|B_i|-1\). The other \(B\)-cycle stays literal.
Choosing the old \(B_\ell\)-cycle exchanges the two \(B\) roles.
The \(A\) roles remain fixed and
\(\alpha\) becomes \(|A_i|-\alpha\).
The new \(\beta\) can again be prescribed arbitrarily.
In either construction the entire factor \(A\) stays literal.

Exchange the heavy factors by \(Q_1\).
The core adjacencies remain the same, and (2) becomes
\[
        (\alpha',\beta')=(|B_i|-\beta,\ |A_i|-\alpha).       \tag{7}
\]
Indeed the second factor is \(B^{-1}AB\): it carries the
old \(A\) cyclic orders with it, while on each old
\(B\)-cycle the shared sheets move back by one.
The new path orientation reverses the two attachment orders.
Apply the preceding construction and undo the exchange.
This permits independent control of \(\alpha\), or exchange
of the two \(A\) roles while preserving the \(B\) roles.
The full operation is pure and keeps the original \(B\)
factor literal.

Choose the smaller length of each profile as its internal
cycle length; role exchanges make this possible. Prescribe
\(\alpha=\beta=1\), then rotate the main \(T\)-endpoint to \(q\).
Y1 gives one absolute configuration with these data.
Equal lengths introduce no extra role distinction because
the cycles within a heavy factor are unordered.
Every construction is pure: its atomic braids are squares,
or a pure operation between a heavy exchange and its inverse.

Finally \(A,B\) fix the isolated sheet \(v\) and generate
a transitive group on the other \(N-1\) sheets.
Conjugates of \(T=(v,q)\) therefore include every
\((v,w)\), \(w\ne v\). They generate \(S_N\).
\(\square\)

## Theorem Y4. All components, markings, groups and sizes

For (1), put
\[
 e_A=1+\boldsymbol1_{a_1=a_2},\quad
 e_B=1+\boldsymbol1_{b_1=b_2},\quad
 H=\frac{N^2(A_0-1)(B_0-1)}{e_Ae_B},\quad
 W=3\bigl(2-\boldsymbol1_{\{a_1,a_2\}=\{b_1,b_2\}}\bigr),   \tag{8}
\]
where equality of the braces means equality of multisets.

If either within-color pair has unequal lengths, there is exactly
one full literal and one full absolute orbit. On each fixed ordered
passport there is exactly one pure literal and one pure absolute orbit.
Their group is \(S_N\), their block-count divisors are \(1,N\), and
their deck group is trivial. Their fixed-order pure sizes are
\[
                    H\quad\hbox{and}\quad H/N,             \tag{9}
\]
respectively. Full sizes are \(W\) times these.

If both pairs repeat, write \(a_1=a_2=a\), \(b_1=b_2=b\),
\(m=a+b-1\). There are exactly two components under either convention:

| Component | Actual group | All block-count divisors | Deck order | Fixed-order pure literal size | Fixed-order pure absolute size |
|---|---|---|---|---|---|
| Primitive | \(S_{2m}\) | \(1,2m\) | \(1\) | \(H-m^2\) | \((H-m^2)/(2m)\) |
| Quadratic lift | \(S_2\wr K\) | \(1,m,2m\) | \(2\) | \(m^2\) | \(m\) |

Here \(K=A_m\) for \(a,b\) both odd and \(K=S_m\) otherwise.
Again multiply by \(W\) for full literal and absolute sizes.

**Proof.** Y1--Y3 put every nonaligned configuration in one
primitive pure absolute orbit. The exceptional aligned matching
can occur only when both pairs repeat. It is imprimitive,
so cannot join that orbit. It has one pure absolute orbit:
rotation (4) carries its common \(T\)-anchor to every one
of the \(m\) identified component sheets. Thus it exhausts
all aligned configurations. U2 supplies its group, blocks,
deck and fixed-order literal size \(m^2\).

The full twist in \(PB_3\) realizes simultaneous conjugation
by the full product. Every relabeling preserving the fixed
literal product is a power of that product. Consequently
each pure absolute orbit lifts to just one pure literal orbit.
Full braids cannot merge the primitive and imprimitive orbits,
and Y2--Y3 already show no other full orbit can occur.

The corrected Goulden--Jackson ordered count gives (8):
the heavy-factor contributions are
\((B_0-1)/e_A\) and \((A_0-1)/e_B\);
the transposition contribution is one. Its leading factor
is \(N^2\). If no quadratic component exists, this is the
primitive size. Otherwise subtract \(m^2\).
The absolute/literal conversion is \(N/|\mathrm{Deck}|\).
There are \(W\) distinct type placements, proving all the
displayed counts. \(\square\)

## Corollary Y5. The complete two-heavy, one-simple first deficit

Let a minimal full-cycle polynomial passport have two
non-transposition finite factors, each of index at least two,
and one literal transposition. Suppose it is at the first
deficit \(s=p-k-1\), so \(p=4\).
It has exactly one primitive \(S_N\) component.
It has a second component exactly for the repeated pairs
in Y4. This is a complete classification under full braids
and under pure braids on each fixed ordered passport.
Every group, block system, deck group and separate size
is explicit.

**Proof.** The nontrivial-cycle counts of the two heavy
factors split as \(2+2\) or \(3+1\).
The former is Y4.
In the latter let the factor with three nontrivial cycles have profile
\[
        (\ell_1,\ell_2,\ell_3,1^{q-2}),\qquad \ell_j\ge2,
\]
and let the other factor have one cycle of length \(q\ge3\).
Then \(N=\ell_1+\ell_2+\ell_3+q-2\).
This is F2's one-arbitrary-factor family with two
single-cycle tails of lengths \(q,2\).
The even tail gives one component. U1 identifies its group
as \(S_N\); equivalently F3 gives this because the heavy
partition contains a fixed point and a tail is a transposition.
Its blocks are \(1,N\), and its deck is trivial.

All three factor types are distinct in this \(3+1\) case.
A full braid taking one fixed-order tuple to another must
therefore have trivial strand permutation. Hence the
single full component gives one pure component on each
ordered passport; the full twist supplies the literal conclusion.
If \(M=\prod_d(\#\{j:\ell_j=d\})!\), its ordered literal count is
\[
                      H_{3+1}=N^2\,q(q-1)/M.               \tag{10}
\]
This is the corrected Goulden--Jackson count, since the
factor with three nontrivial cycles contributes \(q!/((q-2)!M)\), and the
single-cycle factors each contribute one.
Its ordered absolute size is \(H_{3+1}/N\); full sizes
are six times the ordered sizes.
\(\square\)

## Exact evidence and remaining scope

The [checker](verify_polynomial_two_cycle_pair.py) tests general
four-core recognition, synchronized rotations, every mismatch escape,
both distance changes, both role exchanges and exact canonical
normalization. It includes an unequal matching that cannot escape
by opposite-exclusive support alone but does escape through a
shared sheet. Its component checks use only pure braid generators,
including when the two heavy profiles are equal.

The [receipt](polynomial_two_cycle_pair_verification.json)
and [source ledger](polynomial_two_cycle_pair_source_ledger.json)
record finite domains, literal seeds, orbit hashes, actual groups,
all blocks and deck groups, and independent character totals.
They control the written all-length proof; no infinite assertion
is inferred from a finite range.

Run from the canonical checkout:

    python3 research/arithmetic/verify_polynomial_two_cycle_pair.py

This resolves the full \(k=2,s=1,p=4\) first-deficit sector.
The remaining first-deficit problem includes \(k=2,p\ge5\)
with \(s=p-3\ge2\), and \(k\ge3\), outside already resolved
families. Deeper deficits, zero-simple passports, larger-fibre
lifting transport and the complete obstruction classification
remain in the overall goal.
No assertion about every auxiliary residue curve before its
finite label quotient, or about its genus or rationality, follows.
Independent mathematical review and replay remain pending.
