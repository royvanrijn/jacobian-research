# Complete components for a five-cycle partition and two three-cycles

Let
\[
  1\le \ell_0\le\ell_1\le\ell_2\le\ell_3\le\ell_4,\qquad
  N=\sum_{j=0}^4\ell_j,\qquad D=\gcd(\ell_0,\ldots,\ell_4).
                                                               \tag{1}
\]
Consider the unordered finite passport
\[
     (\ell_0,\ell_1,\ell_2,\ell_3,\ell_4),\
     (3,1^{N-3}),\ (3,1^{N-3}),                                \tag{2}
\]
with fixed product \(P_N=(0\,1\,\ldots\,N-1)\), excluding the
identity first factor \((1,1,1,1,1)\).
Its index sum is \((N-5)+2+2=N-1\).
We use full braids \(Q_i(X,Y)=(Y,Y^{-1}XY)\), and quotient by
simultaneous powers of \(P_N\) for absolute components, as in the
[polynomial passport conventions](POLYNOMIAL_PASSPORT_COMPONENTS.md).
The two three-cycle branch values are interchangeable. Initially
label the five cycles; only for distinct lengths are these labels
intrinsic. Put \(M=\prod_{a\ge1}m_a!\), where \(m_a\) is the
multiplicity of the length \(a\).

There are exactly two full braid components if the five lengths
are distinct, and exactly one otherwise. Each distinct-length
component has \(36N\) absolute classes and \(36N^2\) literal tuples.
For a repeated length the sole component has \(72N/M\) absolute
classes and \(72N^2/M\) literal tuples, except for
\((1,1,1,1,3)\), where all three factors are three-cycles and
these counts are \(7\) and \(49\).
Every deck group is trivial. The actual embedded monodromy
group, explicitly determined below, is the same throughout.
For \(D=1\) this gives a
complete primitive \(S_N/A_N\) family in all degrees covered by (1).
This is a classification of (2), not of arbitrary passports.

## Literature and proof dependencies

Jones--Zvonkin,
[Section 1.1, Theorem 1.1 and Section 1.3](sources/jones-zvonkin-braid-cacti-2002.pdf),
give the cactus braid conventions, the connection to polynomial
topological classes, and the Goulden--Jackson count. Their Section 1.3
gives \(N^2(4!)/M=24N^2/M\) literal tuples with the five-cycle factor
first; the two three-cycle class factors each contribute one.
There are three distinct placements unless all three factors
have the same cycle type. This supplies a count comparison,
not the two-orbit proof. Their detailed classifications concern the
specified exceptional primitive groups; their historical discussion
does not establish a current literature completeness claim for (2).
No novelty or priority is asserted.

The argument below uses the
[cycle-incidence tree and full block-divisor lemma](POLYNOMIAL_PASSPORT_COMPONENTS.md#lemma-h1-the-cycle-incidence-tree).
Its alternating normal-closure step is proved here. The
[degree-fifteen quadratic example](POLYNOMIAL_INTRINSIC_CLASSES_AND_LIFTING.md#theorem-d7-a-primitive-three-value-quadratic-lifting-split)
is the smallest strict specialization \((1,2,3,4,5)\); its lifting computation
is separate from the present all-degree component argument.

## Lemma E1. Exact tripod coordinates

Put the five-cycle factor first and call the tuple \((g,u,v)\).
Contract each of the five cycles of \(g\), including a fixed point
as a one-sheet component. The contracted incidence tree has five
cycle vertices and two degree-three vertices. Thus \(u\) and \(v\)
each meet three different \(g\)-cycles, and their cycle-label sets
have exactly one common label \(c\).

Read \(u=(x_c,x_a,x_b)\) starting at its point on cycle \(c\), and
\(v=(y_c,x_d,x_e)\) starting at its point on the same cycle.
Then \((c,a,b,d,e)\) is a permutation of \(0,\ldots,4\).
Set
\[
          y_c=g^\delta x_c,\qquad 0\le\delta<\ell_c.
                                                               \tag{3}
\]
With the five cycles individually labeled, modulo simultaneous
\(S_N\)-relabeling the tuple is determined uniquely by
\[
                     (c,a,b,d,e;\delta).                      \tag{4}
\]
Every such coordinate occurs. Their number is \(24N\).
Forgetting labels quotients by
\(\prod_a S_{m_a}\), acting freely by renaming cycle labels in (4).
Thus the unlabelled count is \(24N/M\).
Every tuple, even after forgetting labels, has trivial
simultaneous centralizer.

**Proof.** Contracting the five disjoint \(g\)-components preserves
the tree property. Its six incidences force all five cycle
vertices to occur, with one of degree two and four of degree one.
Rotate each leaf cycle so that its unique attachment is at its
chosen origin, and rotate the central cycle so that \(x_c\) is
at its origin. The remaining attachment is exactly (3).
The labels leave no further cycle permutations. Conversely,
attach two ordered tripods at the prescribed origins and offset;
the incidence graph is a tree, so the product is a full cycle.

A permutation centralizing all three factors must preserve the
unique central \(g\)-cycle, identified by meeting both supports.
It acts there by a rotation and fixes the unique point of
\(\operatorname{supp}(u)\) on that cycle. Commuting with \(u\)
then fixes its other two support points. The same argument with
\(v\) fixes the other two leaf attachments. Thus each leaf cycle
is preserved individually, even if its length repeats, and all
five cycle rotations are trivial.
For each central label \(c\) there are \(4!\ell_c\) coordinates,
giving \(24N\). An equal-length relabeling fixes (4) only if it
fixes all five entries, so the label-forgetting action is free.
Normalizing the product to \(P_N\) identifies this
simultaneous-\(S_N\) quotient with the fixed-product absolute
quotient. \(\square\)

## Theorem E2. Two components exactly for distinct lengths

On (4) define
\[
               \varepsilon=\operatorname{sgn}(c,a,b,d,e).
                                                               \tag{5}
\]
On the cycle-labeled cover the two signs are exactly the two
components. For distinct lengths this is the actual component
classification of (2), with \(36N\) absolute states and
\(36N^2\) literal tuples per component. If a length repeats,
there is one component, with the counts stated above.

**Proof.** First restrict to braids returning the distinguished
factor to the first position. This subgroup is generated by
\(T=Q_1^2\) and \(B=Q_2\).
For completeness, use coset representatives
\(1,Q_1,Q_1Q_2\), according to that factor's position.
The nontrivial Schreier generators are \(B,T\) and
\(Q_1Q_2^2Q_1^{-1}\); the last is \(B^{-1}TB\), by the braid
relation. With repeated cycles their labels are transported by
the conjugations in the braid formula. If all three factors have
the same type, temporarily distinguish the factor \(g\); the
same subgroup argument applies to that marked branch value.

Here is the exact action on coordinates. If \(\delta>0\), then
\[
 B(c,a,b,d,e;\delta)
       =(c,d,e,a,b;\ell_c-\delta).
                                                               \tag{6}
\]
The supports of \(u,v\) are disjoint in this case. If \(\delta=0\),
their supports meet at their central point, and
\[
 B(c,a,b,d,e;0)=(e,c,d,a,b;0).                                 \tag{7}
\]
Indeed the new pair is
\(((x_c,x_d,x_e),(x_e,x_a,x_b))\).

To compute \(T\), put \(H=gu\). The pair \((g,u)\) is conjugated
by \(H^{-1}\) under \(Q_1^2\). Relabeling the whole tuple by \(H\)
therefore restores \((g,u)\) and replaces \(v\) by \(HvH^{-1}\).
The three \(g\)-cycles meeting \(u\) form one \(H\)-cycle.
Thus \(T\) moves the attachment of \(v\) one step along that
merged cycle. In coordinates,
\[
 \begin{split}
 T(c,a,b,d,e;\delta)
     &=(c,a,b,d,e;\delta+1\bmod\ell_c),&&\delta>0,\\
 T(c,a,b,d,e;0)
     &=(a,b,c,d,e;1\bmod\ell_a).&
 \end{split}                                                  \tag{8}
\]
For example, when the attachment is at \(x_c\), the next
position is \(g x_a\), rather than \(g x_c\).
The \(T\)-cycle has length \(\ell_c+\ell_a+\ell_b\);
its three zero-offset states have central labels \(c,a,b\).

Every permutation of labels in (6)--(8) is even.
This proves invariance of (5) under the entire stabilizer of the
first type placement. Any full braid path between two tuples with
that placement lies in this stabilizer, so opposite signs cannot
join through another placement.

It remains to prove sufficiency. Equation (8) moves every state
to offset zero. Starting at offset zero, \(T^{\ell_a}\) gives
\[
       (c,a,b,d,e;0)\longmapsto(a,b,c,d,e;0).                  \tag{9}
\]
Equation (7) gives the second operation
\[
       (c,a,b,d,e;0)\longmapsto(e,c,d,a,b;0).                 \tag{10}
\]
On the five positions these are a three-cycle and a five-cycle.
Conjugating the former by the latter gives a three-cycle with
exactly one common letter. These two three-cycles generate
\(A_5\), by the elementary one-letter overlap calculation used
in E3 below and retained explicitly by the checker.
Thus (9)--(10) act transitively on
each of the two parity classes of the \(5!\) zero-offset states.
Every state of the same sign is connected.

For distinct lengths and each \(c,\delta\), twelve of the \(4!\)
orders of the other labels have either sign. This gives \(12N\).
Each component meets all three type placements, equally: any
fixed braid moving the distinguished factor induces a bijection.
E1 gives free \(P_N\)-conjugation orbits of size \(N\); the full
twist acts by this conjugation and belongs to every literal
component. Hence absolute components lift to single literal
components with the displayed counts.

Relabeling the five cycle marks commutes with braid transport.
If a length repeats, a transposition of two equal-length labels
interchanges the two cycle-labeled components. Forgetting these
labels therefore leaves one component. E1 gives \(24N/M\)
absolute classes with the distinguished factor first. Usually
there are three different type placements, giving \(72N/M\);
free product conjugation multiplies this by \(N\).
The only case where the distinguished type equals the three-cycle
type is \((1,1,1,1,3)\), in degree seven. Forgetting which of
the three factors was temporarily distinguished divides the
three-placement count by three. This gives \(7\) absolute
classes and \(49\) literal tuples. \(\square\)

## Theorem E3. Actual group, every block system, and deck symmetry

Let \(n=N/D\), let
\[
 \mathcal B_D=\{B_j:j\in\mathbb Z/D\},\qquad
 B_j=\{j,j+D,\ldots,j+(n-1)D\},
\]
and put \(K=\prod_{j=0}^{D-1}A(B_j)\).
For every fixed-product tuple in (2), its actual generated group is
\[
                     G=K\langle P_N\rangle.                  \tag{11}
\]
In particular this is the same embedded group throughout the passport.
Its order is
\[
       |G|=(n!/2)^D
             \begin{cases}D,&n\text{ odd},\\2D,&n\text{ even}.
             \end{cases}                                    \tag{12}
\]
Its block-count divisors are exactly
\[
                         \operatorname{div}(D)\ \cup\ \{N\}.
                                                               \tag{13}
\]
Its centralizer in \(S_N\), and hence its deck group, is trivial.
For \(D=1\), this says \(G=A_N\) when \(N\) is odd and \(G=S_N\)
when \(N\) is even.

**Proof.** A three-cycle preserving a proper block system must
fix every block: moving any block of size \(b>1\) in a nontrivial
three-cycle would move at least \(3b\) sheets. Its support is
therefore contained in one block. In particular \(b=2\) is
impossible. On the block quotient only \(g\) acts, and it must
be a full cycle on the \(m\) blocks. Every \(g\)-cycle length
is consequently divisible by \(m\), so \(m\mid D\).

Conversely, if \(m\mid D\), give each \(g\)-cycle a phase in
\(\mathbb Z/m\) increasing by one under \(g\). On the central
cycle choose phase zero at \(x_c\), and hence phase \(\delta\)
at \(y_c\). Set phase zero at both \(u\)-leaf attachments and
phase \(\delta\) at both \(v\)-leaf attachments. These choices
are consistent since all lengths are multiples of \(m\).
Every phase occurs \(N/m\) times, \(g\) advances it, and \(u,v\)
preserve it. This constructs a block system for every \(m\mid D\).
After product normalization it is necessarily the residue
partition \(\mathcal B_m\). This proves (13).

Let \(L\) be the normal closure of \(u,v\) in \(G\). The connected
components of the hypergraph of their conjugate three-cycle
supports are its orbits and form a \(G\)-block system.
A connected family of three-cycle supports generates the
alternating group on its union. Here is an elementary argument
for the step used repeatedly: starting with a three-cycle gives
\(A_3\). If an added support meets the old union in two points,
the two relevant three-cycles generate \(A_4\). If it meets it in
one point, they generate \(A_5\). In either case, conjugating by
the alternating group already obtained supplies all three-cycles
on the enlarged union. Iterating along a spanning sequence of
intersecting supports proves the assertion. Disjoint connected
components give independent alternating factors.

If there are \(k\) such components, they are proper blocks,
so \(k\mid D\); the one-component case has \(k=1\).
Since \(L\) fixes every block of \(\mathcal B_D\), also \(k\ge D\).
Therefore \(k=D\) and \(L=K\). As \(u,v\in K\) and \(guv=P_N\),
the group is exactly (11).

Now \(P_N^j\in K\) first requires \(D\mid j\).
The restriction of \(P_N^D\) on each block is an \(n\)-cycle,
even exactly when \(n\) is odd. This gives (12).
The five positive lengths divided by \(D\) imply \(n\ge5\).
A permutation centralizing \(K\) must preserve the support of
each of its independent alternating factors. On each block it
centralizes \(A_n\), whose natural centralizer is trivial:
commuting with the point stabilizer \(A_{n-1}\) forces each
point to be fixed. This proves the centralizer assertion.
\(\square\)

## Scope of the lifting information

The orientation (5) is already a complete component label for
(2) with distinct lengths, including imprimitive cases.
Repeated lengths instead give one component.
For \((1,2,3,4,5)\), the
existing D7 computation identifies the two labels with opposite
quadratic reduced lifts and identical intrinsic \(A_{15}\)-classes.
This specializes E2 to two literal components of \(8100\).

More generally, for distinct lengths, when \(D=1\), \(N\) is odd,
and the number of even lengths is positive, that number is two
or four. D7's
written strict-class criterion supplies a class-splitting stem
double cover of \(A_N\) for the prescribed finite and infinity
classes. It does not by itself determine whether its two
values separate E2's components. No uniform value formula is
asserted here. The all-odd case also requires retaining the
intrinsic alternating-group class of the first factor.
For \(D>1\), (11) identifies the actual group but a full reduced
multiplier computation is not supplied.

## Verification and remaining scope

The [checker](verify_polynomial_two_tripods.py) and
[receipt](polynomial_two_tripods_verification.json) compare
the constant-size coordinate formulas with raw permutation
braids, exhaust the indicated compact state spaces, replay
selected complete literal components and independent character
totals, and test block and three-cycle normal-support certificates.
The fixed five-position \(A_5\) calculation, the \(A_4/A_5\)
overlap controls and deliberately invalid coordinates are
checked separately. These bounded checks support the written
unrestricted proofs E1--E3; no degree cutoff proves them.

The [source ledger](polynomial_two_tripods_source_ledger.json)
pins exact dependencies and primary literature. Review and
replay are internal. Independent replay, external review and
formal verification are not asserted. More three-cycle factors,
arbitrary higher-index factors and
the full polynomial component classification remain outside
this theorem.
