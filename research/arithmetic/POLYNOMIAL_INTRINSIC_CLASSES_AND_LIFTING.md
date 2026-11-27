# Intrinsic branch classes and genuine polynomial lifting obstructions

A cycle partition in \(S_N\), a conjugacy class in the actual generated
group \(G\), and a lifting class in a central extension of \(G\) are
different data. Component classification must retain that distinction.
The two-bundled-involution family needs an intrinsic conjugacy-class
refinement. The seven-sheet \(A_7\) pair and the eight-sheet coupling
pair admit genuine cubic lifting separation even after that refinement.
A fifteen-sheet primitive \(A_{15}\) passport admits quadratic lifting
separation with three finite branch values.

We use the [fixed-product polynomial conventions](POLYNOMIAL_PASSPORT_COMPONENTS.md):
\(\prod_i g_i=P_N\), total finite index \(N-1\), full braids
\(Q_i(X,Y)=(Y,Y^{-1}XY)\), and absolute quotient by powers of \(P_N\).
Every assertion below specifies its actual group and lift choices.
No general component completeness theorem is inferred from these cases.

## Primary sources and the applicable scope

Roberts--Venkatesh, [Section 4B, equations (4-1)--(4-2),
printed p518](sources/roberts-venkatesh-hurwitz-monodromy-2015.pdf),
define the reduced Schur multiplier by killing commutators of lifts of
\(x\in C\) and \(y\in C_G(x)\), for the specified **\(G\)-classes**.
Their Sections 3D and 4C--4D distinguish these classes, lift torsors,
and braid invariants. Their asymptotic connectivity statement requires
sufficiently many branch points of every prescribed class. It is not
applied to the minimal polynomial passports here.

Akita--Kawasaki--Satoh,
[Section 2.1 and Theorem 2.2](sources/akita-kawasaki-satoh-schur-covers-2024.pdf),
give \(H_2(S_N,\mathbb Z)=C_2\) for \(N\ge4\) and Schur-cover
presentations in which lifts of disjoint transpositions anticommute.
The retained primary author manuscript contains the actual statement
and the explicit cocycle calculation. We use its unit-root realization:
\(t_i=(e_i-e_{i+1})/\sqrt2\) in the real Clifford algebra
\(e_i^2=1,\ e_ie_j=-e_je_i\) for \(i\ne j\).
These roots satisfy that Schur presentation with \(t_i^2=1\)
and central \(z=-1\); their alternating subgroup is the spin double
cover used in D4 and D7.

The [ATLAS primary group page](sources/atlas-a7-group-2026.html)
records the multiplier of \(A_7\) as \(C_6\). Its raw
[45-point generator A](sources/atlas-3a7-45a-generator-1.mtx) and
[generator B](sources/atlas-3a7-45a-generator-2.mtx), together with the
[seven-point A](sources/atlas-a7-7-generator-1.gap) and
[seven-point B](sources/atlas-a7-7-generator-2.gap), supply exact
permutation data. The checker verifies the entire central extension;
it does not rely on the database's historical verification flags.

## Lemma D1. Canonical coprime lifts and class-defined lifts

Let \(1\to C\to\widetilde G\to G\to1\) be a finite central extension.
If \(x\in G\) has order prime to \(|C|\), it has a unique lift
\(\widehat x\) of the same order. These lifts are conjugacy-equivariant.
In particular their commutators with arbitrary lifts of \(C_G(x)\)
are trivial.

If the finite factors and their product \(P\) have such orders, then
\[
                 \lambda(g)=
              \widehat g_1\cdots\widehat g_r\widehat P^{-1}
                              \in C                         \tag{1}
\]
is invariant under the full braid action. More generally, the same
conclusion holds for any specified conjugacy-equivariant lift sections
on the branch classes and at infinity. When the extension is stem and
these classes are lift-split, (1) is a reduced Schur lifting invariant.

**Proof.** If \(x\) has order \(k\), any lift \(X\) has \(X^k\in C\).
The power map \(c\mapsto c^k\) is an automorphism of \(C\), so exactly
one \(cX\) has \(k\)-th power one. Its order is \(k\) by projection.
Conjugation preserves that order, proving equivariance and the
centralizer assertion.

Under one Hurwitz move the lifted pair is
\((\widehat Y,\widehat Y^{-1}\widehat X\widehat Y)\), whose product
is \(\widehat X\widehat Y\). This proves (1), and the same calculation
works for the stated general sections. The centralizer assertion
means that none of this kernel is killed in the corresponding
class-splitting quotient. A stem extension with that property is a
quotient of the reduced Schur cover of the specified branch data.
\(\square\)

One may also pull back an extension along \(\phi:G\to H\) and use
canonical lifts of the images \(\phi(g_i)\). It is the orders of
these images that must be coprime to the kernel. Stemness of a
pullback requires its own argument; it is not automatic.

## Proposition D2. The involution-family partition is intrinsic class data

For the [two-bundled-involution family](POLYNOMIAL_INVOLUTION_COMPONENT_CLASSIFICATION.md),
the actual embedded group together with the unordered list of its
intrinsic branch conjugacy classes is a complete component invariant.
In particular the two components of
\[
                       (2^p),(2^{p-3}1^6),T,T                \tag{2}
\]
for an odd prime \(p\ge5\) have different \(G\)-classes of the second
bundled factor, despite having the same \(G=S_p\wr C_2\).
They are different Nielsen data before a lifting invariant is chosen.

**Proof.** The normal transposition-support blocks determine \(n,m\)
from the actual group. When \(m\) is odd, B1 already gives one component
for that group. When \(m\) is even, B3 reads the extra partition from
the two fixed quotient cycles of the smaller bundled factor.
Conjugation in \(G\) permutes those cycles and conjugates their return
permutations. Their unordered fixed-point counts are therefore
invariant under intrinsic \(G\)-conjugacy. Different degree partitions
give different \(G\)-classes. B1 proves sufficiency of the resulting
labels.

For \(p=5\), the second factor has return types \(({\rm DT},1)\) or
\((T,T)\) on the two blocks. Their \(S_5\wr C_2\)-class sizes are
\(30\) and \(100\), respectively, by their centralizer orders
\(8\cdot120\) and \(2\cdot12^2\). The exact group conjugacy computation
agrees. \(\square\)

This does not apply to the examples D4, D5 and D7: each of those
pairs has identical intrinsic branch classes.

## Theorem D3. A transposition kills the reduced multiplier of \(S_N\)

For \(N\ge4\), if the prescribed \(S_N\)-classes include the
transposition class, then
\[
                         H_2(S_N,C)=0.                      \tag{3}
\]
Thus reduced Schur-multiplier data cannot separate components of a
primitive \(S_N\) polynomial passport containing a transposition.
This assertion does not imply connectivity.

**Proof.** The ordinary multiplier is \(C_2\), by the primary theorem
quoted above. The transpositions \(x=(0\,1)\), \(y=(2\,3)\) commute,
and \(y\in C_{S_N}(x)\). In either Schur cover their lifts have
commutator the nontrivial central element, by the disjoint-generator
relation in Theorem 2.2. The defining reduced quotient therefore kills
the whole multiplier. Adding other branch classes cannot restore it.
\(\square\)

## Theorem D4. The seven-sheet split is a cubic reduced lifting split

For \(G=A_7\), let \(C\) consist of the DT class and the class of
\(P_7^{-1}\) at infinity. Then
\[
                         H_2(A_7,C)=C_3.                    \tag{4}
\]
In its triple cover, use the unique order-two lifts of the DTs and
the unique order-seven lift of \(P_7\).
Among all \(392\) literal triples of DTs with product \(P_7\),
the values of (1) have counts
\[
 \begin{array}{c|ccc}
 \lambda&1&z&z^{-1}\\ \hline
 \text{literal tuples}&98&147&147 .
 \end{array}                                                 \tag{5}
\]
Here \(z\) is a chosen generator of the central \(C_3\).
The \(z\) and \(z^{-1}\) sets are exactly the two \(A_7\) braid
components. The neutral set consists of the two \(49\)-tuple
components with different embedded \(\operatorname{PSL}_3(2)\)
groups. The DTs form one intrinsic \(A_7\)-class, and infinity is
the same fixed \(P_7\) throughout.

**Proof.** The ATLAS multiplier bound is \(C_6\). Its 2-primary part
is killed by a commuting Klein four pair:
\[
                 x=(0\,1)(2\,3),\qquad y=(0\,2)(1\,3).
\]
Their standard spin lifts have commutator \(-1\). This can be seen
in the quaternion preimage of the Klein four group, or directly
from the eight unit-root factors retained by the checker.

The 3-primary part is not killed: DTs have order two, and the
infinity permutation has order seven, both prime to three.
D1 shows that every relevant centralizer commutator vanishes in
the triple cover. Hence the quotient is exactly \(C_3\).

For an explicit model, close the paired ATLAS generators in degrees
45 and seven. There are \(7560\) upper permutations, \(2520\)
lower permutations, and three upper elements over every lower
one. The three-element kernel is central. The upper generators
have orders \(3,5\), and their product has order seven.
These orders force trivial abelianization: in an abelian quotient,
the product of the 3- and 5-primary generators cannot have order
dividing seven unless both vanish. Thus this is a stem triple cover.
The imported cyclic multiplier identifies it with the 3-primary
Schur cover.

Independently enumerate the \(105\) DTs and complete a triple from
its first two factors using the fixed product. This gives exactly
\(392\) tuples. Their full braid closures have sizes
\(147,147,49,49\); direct generated-group closure gives orders
\(2520,2520,168,168\). Exact unique-order lifts on every tuple
give (5), with nonzero values on the two \(A_7\) components.
The checker also verifies the invariant on every signed braid edge.
Together with these complete closures, the values prove the stated
separation. \(\square\)

The choice of \(z\) labels the two nonzero values; replacing \(z\)
by its inverse exchanges their names. The two components themselves
are independent of that convention. No quotient by outer
automorphisms is silently added.

## Theorem D5. A cubic area invariant explains the eight-sheet coupling

Let
\[
 H=\{(u,v;\sigma)\in S_4\wr S_2:
                 \operatorname{sgn}(u)\operatorname{sgn}(v)
                         =\operatorname{sgn}(\sigma)\}.       \tag{6}
\]
The two \(16\)-tuple components of L8 have the same literal \(H\),
the same intrinsic class of their \((4,2,2)\) factor, and the same
intrinsic class of their DT factor. Nevertheless a class-splitting
stem central \(C_3\)-extension of \(H\) gives them the two nonzero
values of (1).

**Construction and proof.** Let \(S_4\) act on its three perfect
matchings. This gives \(S_4\to S_3\) with kernel \(V_4\), preserving
sign. Identify \(S_3\) with the affine permutations
\(t\mapsto\epsilon t+c\) of \(\mathbb F_3\).
The induced affine action of (6) on \(\mathbb F_3^2\) has determinant
one. Its quotient is
\[
          Q=\mathbb F_3^2\rtimes\langle J\rangle,\qquad
              J(v_1,v_2)=(-v_2,v_1),\quad |J|=4,             \tag{7}
\]
of order 36, with kernel \(V_4^2\).

For \(\omega(v,w)=v_1w_2-v_2w_1\), define a central extension of
\(Q\) by
\[
 (v,j,z)(w,k,t)=
 \bigl(v+J^jw,\ j+k,\ z+t+\tfrac12\omega(v,J^jw)\bigr).       \tag{8}
\]
The last coordinate lies in \(\mathbb F_3\), where \(1/2=2\);
the middle coordinate lies in \(\mathbb Z/4\).
Bilinearity and \(J\)-invariance of \(\omega\) prove associativity.
This is the order-108 Heisenberg semidirect product. Pulling it
back to \(H\) gives a central extension of order 1728.

Two commuting 3-cycles, one in each four-sheet block, have
independent translation images in (7). The commutator of any
lifts is the nonzero central area
\(\omega((1,0),(0,1))\). Hence the central kernel is in the
derived subgroup, so the pullback is stem.
The finite branch orders are four and two, and infinity has
order eight. D1 gives canonical lifts and shows that all these
classes are lift-split.

Exact conjugation in \(H\) verifies that the two components use
the same factor classes, of sizes \(72\) and \(36\).
Applying (8) to their canonical lifts gives \(z\) and \(z^{-1}\).
The complete \(16\)-tuple closures from L8, and all their signed
edges, have those respective values. Thus this is a genuine cubic
lifting distinction after fixing the intrinsic branch classes.
\(\square\)

Only this central \(C_3\) quotient is asserted for \(H\); no
calculation of its entire reduced Schur multiplier is claimed.

## Theorem D6. An infinite family with identical intrinsic classes

For every \(b\ge2\), consider the block-size-\(b\) sector of the
degree-\(7b\) passport
\[
                  (2^{2b}1^{3b})^3,\qquad(b\,1^{6b}).         \tag{9}
\]
The two components with literal group \(C_b\wr A_7\) have identical
intrinsic branch classes. Their two nonzero cubic quotient-lifting
values distinguish them. Each has literal size \(4116\), absolute
size \(588\), and deck-group order \(b\), as in L7.

**Proof.** L5--L7 give one lifted component over each seven-sheet
quotient component and prove the exact embedded wreath group.
Every copied DT factor has identity cycle returns, so it is
conjugate by cyclic base elements to its pure quotient permutation.
All quotient DTs belong to the same \(A_7\)-class. The additional
factor is a local \(b\)-cycle; transitivity of \(A_7\) on the
blocks gives it one intrinsic class. Infinity is the fixed
literal \(P_{7b}\). Thus the two components have identical
intrinsic branch data.

Pull back the triple cover of \(A_7\) through the block projection
\(C_b\wr A_7\to A_7\). It is stem, because its restriction to the
diagonal \(A_7\) is the perfect triple cover, whose central kernel
is in its derived subgroup. All finite quotient images have order
two or one; infinity has quotient order seven. The corresponding
canonical quotient lifts define conjugacy-equivariant split lift
classes in the pullback. Their invariant is exactly D4's value
on the quotient, so the two components have values \(z,z^{-1}\).
\(\square\)

This proof includes \(3\mid b\). In that case the order of the local
factor does not choose a unique lift. The specified lift is instead
chosen through its identity quotient image. The conclusion is
restricted to the stated block sector; other components of the
same ambient passport are not excluded.

## Theorem D7. A primitive three-value quadratic lifting split

On fifteen sheets the passport
\[
             (5,4,3,2,1),\qquad(3\,1^{12}),\qquad
                                           (3\,1^{12})       \tag{10}
\]
has exactly two literal full braid components, each of size \(8100\).
Both generate the identical \(A_{15}\), have the same intrinsic
branch classes, trivial deck group and absolute size \(540\).
A class-defined double-cover invariant takes the two values
\(+1,-1\) on these components.

The lift section for the even-order factor has a general
justification: if an even permutation has distinct cycle lengths
and a positive even number of even-length cycles, its class is
lift-split in the spin double cover of the alternating group.

**Proof of the class assertion.** Write its cycles as
\(\gamma_j\), and their unit-root lifts as \(X_j\).
Disjoint cycle lifts anticommute exactly when both have odd
parity. There are a positive even number of such factors.
Consequently the commutator of \(X=\prod_jX_j\) with a lift of
\(\prod_j\gamma_j^{a_j}\) is
\[
                         (-1)^{\sum_{\ell_j\ {\rm even}}a_j},
                                                                  \tag{11}
\]
the sign of the centralizer element. Because the lengths are
distinct, these powers form the entire symmetric centralizer.
Its alternating subgroup therefore centralizes \(X\).
Thus the preimage class splits in the double cover of \(A_N\).
Choosing either of its lift classes gives an equivariant section.
An odd centralizer element does exist, so the ordinary
\(S_N\)-class itself is one \(A_N\)-class.

**Finite component proof.** Put
\[
 \begin{split}
 G_0&=(0\,1\,2\,3\,4)(5\,6\,7\,8)(9\,10\,11)(12\,13),\\
 U&=(4\,8\,11),\qquad V=(11\,13\,14).
 \end{split}                                                  \tag{12}
\]
The first tuple is \((G_0,U,V)\), whose product is \(P_{15}\).
For the second, use \((G_0,U,V^{-1})\) and relabel its full
product to \(P_{15}\) by sending its orbit of zero to
\(0,1,\ldots,14\) in order. The checker retains both literal tuples.
Their full closures are disjoint, each has \(8100\) tuples, and
their total equals the independent exact hook-character count
of \(16200\).

All generators are even and include a 3-cycle. The group is
primitive: any proper block size in degree 15 is at least three,
so each 3-cycle fixes the quotient. The remaining quotient factor
would be a full cycle, yet it has a fixed block because \(G_0\)
has a fixed sheet. This is impossible. A primitive group
containing a 3-cycle contains the alternating group: take all
conjugate 3-cycle supports; their connected components are blocks,
and connected overlapping 3-cycles generate the alternating group
on their union. To see the last step inductively, two 3-cycles
sharing two or one letters generate \(A_4\) or \(A_5\),
respectively; adjoining further overlapping supports gives the
alternating group on the enlarged union. Thus both groups are
\(A_{15}\).

Choose the lift class of \(G_0\) represented by its adjacent-root
cycle expansions. For another factor of type \((5,4,3,2,1)\),
conjugate this lift by an even permutation. Equation (11) makes
the choice of conjugator immaterial. The 3-cycles and infinity
have their unique odd-order lifts. Their product (1) is therefore
a braid invariant. Exact Clifford multiplication gives opposite
signs on (12)'s two normalized tuples. The extension is stem:
the commuting Klein four pair used in D4 has commutator \(-1\)
in the spin cover. All prescribed classes are lift-split, so
this is a reduced Schur lifting obstruction.

Finally the strict partition class is one intrinsic \(A_{15}\)-class
by (11)'s odd-centralizer observation, and the 3-cycle class is
unsplit because its symmetric centralizer contains an odd
transposition on fixed sheets. Infinity is the same fixed
\(P_{15}\). The natural centralizer of \(A_{15}\) is trivial,
giving the deck and absolute counts. \(\square\)

No all-odd-order spin-rigidity theorem applies to (10): its first
factor has order 60, and its lift was chosen by a split class,
not by an arbitrary sign. That factor has index ten, outside
the scope of the index-two mixture connectivity theorem.

## Evidence and the remaining classification problem

The [checker](verify_polynomial_lifting_obstructions.py),
[receipt](polynomial_lifting_obstructions_verification.json), and
[source ledger](polynomial_lifting_obstructions_source_ledger.json)
retain the exact cover data, full paired closure, centralizer tests,
Clifford witnesses, finite component domains, cocycle checks and
the bounded cyclic-lift controls. The all-degree claims D1--D3,
D6 and D7's strict-class criterion use their written proofs.
Finite executions do not prove other passports connected or
classify their lifting data.

The full arbitrary-passport goal remains open. These results
separate intrinsic conjugacy-class refinement from genuine
lifting separation, and show that primitive three-value
connectivity itself fails with a higher-index factor.
General primitive component classification, larger-fibre return
actions and lifting data for other actual groups still require
a complete classification. Present proof audit and replay are
internal; independent review, independent replay and formal
verification are not asserted.
