# Two components for all duplicated-cycle lengths

Let \(a,b\ge2\), \(m=a+b-1\), \(N=2m\), and fix an \(N\)-cycle \(P\).
The minimal finite passport
\[
        (a^2 1^{2b-2}),\qquad (b^2 1^{2a-2}),\qquad T       \tag{1}
\]
has exactly two full Hurwitz components, for every pair of lengths.
One has monodromy \(S_N\), and the other is the quadratic lift
with monodromy \(S_2\wr K\), where \(K=A_m\) if both lengths
are odd and \(K=S_m\) otherwise. The complete group, block,
deck and orbit-size statements are in X4.

The argument uses the four possible forests obtained by removing \(T\).
All configurations except one aligned matching reach a path with four
cycle vertices and one isolated sheet. Two excursions through a star
change the path's two cyclic distances independently. Thus the proof
closes the component problem without requiring geometric irreducibility
of every auxiliary residue equation before its finite label quotient.

## Primary comparison and conventions

Jones--Zvonkin,
[Theorem 1.1, Sections 1.3, 2.2 and 5.4](sources/jones-zvonkin-braid-cacti-2002.pdf),
give the polynomial/cactus/braid correspondence, the relevant count
convention, and examples of constructing imprimitive cacti from copies
of a quotient. Those actual statements and constructions were inspected.
The four-core normalization below is proved directly; their individual
exceptional examples are not a component theorem for all of (1).
Wajnryb,
[Theorem 1, retained p2](sources/wajnryb-two-special-fibers-1995.pdf),
allows at most two nonsimple critical values. Infinity and the two heavy
finite values here give three, so its connectivity conclusion is not
imported. No novelty claim is made.

The group and count input for the quadratic component is
[U2](POLYNOMIAL_FIRST_DEFICIT_COMPONENT_REDUCTION.md#theorem-u2-the-complete-imprimitive-component).
Its count uses the corrected Goulden--Jackson formula, with the printed
hypothesis discrepancy explicitly audited in that source.
The tree input is
[H1](POLYNOMIAL_PASSPORT_COMPONENTS.md).
All new connectivity steps are proved below independently of finite
orbit computations.

We compose permutations right to left and use
\(Q_i(C,D)=(D,D^{-1}CD)\). First work modulo simultaneous
sheet relabeling, allowing the product to range over the conjugacy
class of \(P\). This is the absolute fixed-product quotient.
Move the transposition last and write the tuple as \((A,B,T)\).
The two cycles within either heavy factor are unordered. Branch-type
placements are all allowed; exchanges used below need not preserve
their intermediate order.

## Theorem X1. Four complementary cores

The cycle-incidence graph of \(A,B,T\) is a tree. Removing the
cycle vertex of \(T\) leaves two components, and \(R=AB\) is
a full cycle on each one. Suppress the sheet vertices belonging
to exactly one nontrivial cycle, and suppress each shared sheet
to an edge between its \(A\)- and \(B\)-cycle vertices.
Keep sheets outside both supports as isolated vertices.
Exactly four possibilities remain:

1. a matching of two \(A\)-\(B\) pairs;
2. an \(A\)-centered two-edge star and one isolated \(A\)-cycle vertex;
3. a \(B\)-centered two-edge star and one isolated \(B\)-cycle vertex;
4. a path \(A_{\ell}-B_i-A_i-B_{\ell}\), and one isolated sheet.

The subscripts \(\ell,i\) mean leaf and internal cycle vertices.
In the path case, let
\[
 p=A_{\ell}\cap B_i,\qquad q=A_i\cap B_i,\qquad
 r=A_i\cap B_{\ell}.
\]
The intersections denote their unique sheet. Apart from the
endpoint of \(T\) on the nonisolated component, its absolute
configuration is determined by
\[
 \alpha=\operatorname{dist}_{A_i}(q,r)\in\{1,\ldots,a-1\},
 \qquad
 \beta=\operatorname{dist}_{B_i}(p,q)\in\{1,\ldots,b-1\}.       \tag{2}
\]
Every such pair of distances occurs.

**Proof.** The incidence tree has the \(N\) sheet vertices and five
nontrivial-cycle vertices, with \(2a+2b+2=N+4\) incidences.
Its connectedness follows from the full-cycle product; the minimal
index sum makes it a tree. Deleting the degree-two \(T\) vertex
therefore gives exactly two components.
Expanding the remaining cycles into consecutive transpositions
shows that their restrictions multiply to a full cycle on each
tree component.

The four heavy cycle vertices form a bipartite subgraph of
\(K_{2,2}\). Two cycles of different colors cannot share two
sheets, since that would give an incidence cycle.
If \(e\) is its edge count and \(z\) the number of sheets outside
both heavy supports, the number of components is \(4-e+z=2\).
The core is a forest, so \(e\le3\); hence \(z\le1\).
When \(z=0\), the two edges form a matching or a two-edge star.
When \(z=1\), the three edges form the stated path.

Each internal path cycle has two distinguished sheets; its
oriented distance between them determines its cyclic order up to
relabeling the intervening sheets. A leaf cycle has just one
distinguished sheet. This proves the uniqueness assertion in (2).
Conversely build these four cyclically ordered vertices along a
path and join any sheet of that component to a fresh isolated
sheet by \(T\). The incidence tree gives a full-cycle product,
so every pair in (2) occurs. \(\square\)

## Lemma X2. Synchronized rotation and opposite-color capture

The following transformations are available in the absolute braid orbit.

**Rotation.** Keep \(A,B\) literal and replace \(T=(x,y)\) by
\[
                           (R^k x,R^k y),\qquad R=AB.       \tag{3}
\]
In the path case this puts its nonisolated endpoint at any desired
sheet. In a star case it puts that endpoint at any desired sheet
of the star component.

**Capture.** If \(T=(u,v)\), \(B(u)\ne u\) and \(B(v)=v\), then
\[
 Q_2^2(A,B,T)=
        \bigl(A,\;TBT,\;(u,B^{-1}u)\bigr).                 \tag{4}
\]
In particular, if \(u\) is \(B\)-only and \(v\) is \(A\)-only,
then \(u\) becomes an isolated sheet and \(v\) becomes a shared
sheet in the new heavy forest.

**Proof.** The pair move \(Q_1^2\) simultaneously conjugates \(A,B\)
by \(R^{-1}\), and keeps \(T\) literal. Follow its \(k\)-th power
by simultaneous relabeling by \(R^k\); this restores both heavy
permutations and gives (3).
For (4), the first move gives \((T,TBT)\). The second gives
\[
      \bigl(TBT,\;TB^{-1}TBT\bigr).
 \]
The second entry is the transposition
\((TB^{-1}u,TB^{-1}v)=(B^{-1}u,u)\), since
\(B^{-1}u\notin\{u,v\}\) and \(B^{-1}v=v\).
This proves the identity, including its unchanged first heavy
factor. Conjugation by \(T\) replaces \(u\) by \(v\) in its
chosen \(B\)-cycle, leaving the other \(B\)-cycle literal.
\(\square\)

## Theorem X3. Uniform normalization of every primitive tuple

Every star configuration reaches a path configuration.
Every matching configuration reaches a path unless its two
transposition endpoints agree under the unique color-preserving
cyclic identification of its two components. All path configurations
belong to one absolute braid orbit. They have group \(S_N\).
The excluded aligned matchings form the quadratic imprimitive
component.

**Proof: the stars.** In an \(A\)-centered star, the other component
is a single \(A\)-cycle, so its \(T\) endpoint is \(A\)-only.
The star has \(B\)-only sheets because \(b\ge2\).
Use (3) to put the other endpoint at any such sheet and use
(4). One new shared sheet attaches the isolated \(A\)-cycle
to a \(B\)-cycle, giving the path; the removed \(B\)-only
sheet is now isolated. A \(B\)-centered star is symmetric:
exchange the heavy factors, perform this construction, and
exchange back.

**The matchings.** Each component is an \(a\)-cycle and a
\(b\)-cycle meeting at one sheet, and has size \(m\).
There is a unique color-preserving cyclic identification between
the components, carrying their shared sheets to each other.
In its \(R\)-cyclic order, each component has the word
\[
          O\,B^{\,b-1}A^{\,a-1},                          \tag{5}
\]
where \(O\) is the shared sheet and the other letters are
the two exclusive support types. Write the relative phase of
the two \(T\) endpoints as \(\delta\in\mathbb Z/m\).
Rotation (3) traverses all equal-phase pairs.

Every nonzero phase has a pair with one \(A\)-only and one
\(B\)-only endpoint. Indeed put \(B\)-positions at
\(1,\ldots,b-1\) and \(A\)-positions at \(b,\ldots,m-1\).
Their positive differences fill \(1,\ldots,m-2\);
the negative differences fill \(2,\ldots,m-1\) modulo \(m\).
Their union is every nonzero residue.
Equation (4) therefore gives a path whenever \(\delta\ne0\).
For \(\delta=0\), exchanging the two identified copies
commutes with \(A,B,T\). Its two-point orbits give an invariant
block system. These are precisely free quadratic lifts of the
two-cycle quotient, and U2 gives one imprimitive component.

**Changing \(\beta\) while keeping \(\alpha\).**
Start with a path as in X1. Let \(v\) be the isolated sheet,
and use (3) to set \(T=(p,v)\).
Apply (4) on the internal \(B\)-cycle. It replaces its sheet
\(p\) by \(v\), removing the edge \(A_{\ell}-B_i\) of the core.
The result is an \(A\)-centered star plus the isolated
\(A_{\ell}\)-cycle, with new transposition \((p,B^{-1}p)\).
Both attachment sheets \(q,r\) of \(A_i\) are untouched,
so \(\alpha\) is unchanged.

In this star, rotate the main endpoint of \(T\) to a chosen
\(B\)-only sheet \(y\) of the same internal \(B\)-cycle.
The other endpoint is some sheet \(x\) of \(A_{\ell}\).
Apply (4) again. It replaces \(y\) by \(x\) in that \(B\)-cycle
and restores the edge from \(A_{\ell}\), with \(y\) now
isolated. As \(y\) ranges over the \(b-1\) sheets other than
\(q\) in the cut \(B_i\)-cycle, the new distance from its
attachment \(x\) to \(q\) ranges over all
\(1,\ldots,b-1\). The permutation \(A\) and the entire
other \(B\)-cycle stay literal throughout these absolute
transformations. This proves independent freedom in \(\beta\).

**Changing \(\alpha\).** Exchanging the heavy factors by \(Q_1\)
gives \((B,B^{-1}AB,T)\). It preserves the core adjacencies:
on each \(B\)-cycle every shared sheet moves back by one
place; the \(A\)-cycles are conjugated with their cyclic orders.
In the conventions (2), the new distances are
\[
                         (\alpha',\beta')=(b-\beta,a-\alpha).
                                                               \tag{6}
\]
The isolated sheet remains outside both supports.
Apply the preceding \(\beta'\) construction to obtain any
desired value \(a-\alpha_{\rm new}\), then undo the exchange.
Formula (6) shows that \(\beta\) is retained and \(\alpha\)
has the prescribed value.

Thus all paths reach \((\alpha,\beta)=(1,1)\).
Rotation (3) fixes their remaining \(T\)-anchor freedom, and
X1 identifies their absolute configurations. This proves a
single path orbit.

Finally, in a path the group \(\langle A,B\rangle\) is transitive
on its nonisolated component and fixes \(v\).
Conjugating \(T=(v,x)\) by this group gives every
\((v,y)\) with \(y\ne v\). These transpositions generate
\(S_N\). Every star and nonaligned matching, being braid
equivalent to a path, has the same group. Aligned matchings
are imprimitive, so no such tuple is in the path orbit.
These four core types exhaust all tuples by X1.
\(\square\)

## Corollary X4. Complete groups, blocks, decks and orbit sizes

Put \(W=3(2-\boldsymbol1_{a=b})\) and \(D=(2a-1)(2b-1)\).
For all \(a,b\ge2\), passport (1) has exactly two full
absolute and two full literal fixed-product braid orbits:

| Component | Actual group | All block-count divisors | Deck order | Full literal size | Full absolute size |
|---|---|---|---|---|---|
| Primitive | \(S_{2m}\) | \(1,2m\) | \(1\) | \(Wm^2(D-1)\) | \(Wm(D-1)/2\) |
| Quadratic lift | \(S_2\wr K\) | \(1,m,2m\) | \(2\) | \(Wm^2\) | \(Wm\) |

Here \(K=A_m\) for \(a,b\) both odd and \(K=S_m\) otherwise.
The statement includes equal heavy profiles and every factor-type placement.

**Proof.** X3 proves two absolute orbits and determines the primitive
group directly. For any full-cycle product, the only relabelings
preserving that literal product are its cyclic centralizer.
The full braid twist realizes conjugation by the product, so
an absolute full orbit is also a single literal full orbit.
The quadratic group, blocks, deck and size are U2.
The primitive group gives its blocks and deck immediately.
The corrected Goulden--Jackson ordered total is \(m^2D\);
there are \(W\) type placements. Subtract the quadratic
size from that total to obtain the primitive literal size.
Divide literal sizes by \(N/\lvert\mathrm{Deck}\rvert\)
to get the absolute sizes. \(\square\)

## Exact controls and remaining scope

The [checker](verify_polynomial_duplicated_cycle_paths.py) verifies the
literal braid identity, the twist followed by relabeling, the complete
core classification and the two independent distance changes. Its
domains include every ordered \(2\le a,b\le7\), plus four
asymmetric larger pairs. It tests all path-distance pairs in the square
domain, both star types, all matching phases, and the unchanged
complementary cycles. Separate complete colored-braid closures and
independent hook totals check ten passports, including lengths with
both entries at least four. Their finite scopes and exact hashes are
recorded in the [receipt](polynomial_duplicated_cycle_paths_verification.json)
and [source ledger](polynomial_duplicated_cycle_paths_source_ledger.json).

Run from the canonical checkout:

    python3 research/arithmetic/verify_polynomial_duplicated_cycle_paths.py

This closes the entire duplicated-cycle family, including all lengths
left by the [rational](POLYNOMIAL_DT_DUPLICATED_CYCLE_COMPONENTS.md)
and [quadratic](POLYNOMIAL_DUPLICATED_THREE_CYCLE_COMPONENTS.md)
residue proofs. Those remain separate geometric descriptions.
No claim is made that every pair-ordered auxiliary residue curve is
irreducible, rational or of specified genus.
Other first-deficit profiles, deeper deficits and zero-simple passports
outside the already classified families remain open.
The [general two-cycle-pair theorem](POLYNOMIAL_TWO_CYCLE_PAIR_COMPONENTS.md)
now handles unequal within-color lengths by exchanging the internal
and leaf cycle roles. It also proves pure braid connectivity on
each fixed ordered passport, including in the duplicated sector,
and closes the entire \(k=2,s=1,p=4\) first deficit.
Independent mathematical review and replay remain pending.
