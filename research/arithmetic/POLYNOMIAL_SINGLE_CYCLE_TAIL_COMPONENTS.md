# One arbitrary factor and a tail of single cycles

Fix \(s\ge2\), integers \(q_i\ge2\), and a partition
\[
 \lambda=(\ell_1,\ldots,\ell_k),\qquad
 k=1+\sum_{i=1}^s(q_i-1),\qquad
 N=\sum_{j=1}^k\ell_j,\qquad D=\gcd(\ell_1,\ldots,\ell_k).
                                                               \tag{1}
\]
The factor of type \(\lambda\) is nonidentity. Consider all
tuples with unordered finite passport
\[
          \lambda,\quad(q_1,1^{N-q_1}),\ldots,
                            (q_s,1^{N-q_s})                   \tag{2}
\]
and fixed product \(P_N=(0\,1\,\ldots\,N-1)\).
The index sum is \(N-1\). A tail factor means exactly one
nontrivial cycle, with its other sheets fixed.
Full braids allow all distinct type placements. Absolute
classes quotient by simultaneous powers of \(P_N\).

There are exactly **two components** if every \(q_i\) is odd
and the parts of \(\lambda\) are distinct. In every other case
there is exactly **one component**. The two-component case
has a geometric parity label. The actual embedded monodromy,
all block systems and deck group are also determined below.
This classifies (2); arbitrary passports with two or more
factors having disconnected nontrivial support remain outside
the result. A uniform formula for the sizes of the two separate
orbits is not claimed.

## Primary theorem and the exact bridge to it

Chen--Gendron--Prado--Tahar,
[Isoresidual curves](sources/chen-gendron-prado-tahar-isoresidual-curves-2026.pdf),
*Journal de l'École polytechnique — Mathématiques* **13** (2026),
891--963, DOI [10.5802/jep.339](https://doi.org/10.5802/jep.339),
prove Theorem 1.7 for fixed residues outside the resonance
arrangement. With at least two zeros and only simple poles,
their fiber has two components exactly when all zero orders
are even; otherwise it is connected.
Their poles and zeros are labeled. The pertinent primary
statements and arguments are:

- Definition 1.1, printed p892: resonance means a proper
  nonempty subset of residues sums to zero.
- Theorem 1.7, p895, and its proof in Section 8.4, pp955--957.
- Lemma 8.1, pp948--949: merging any chosen zeros reaches
  every component; each boundary component is adjacent to a
  unique component before merging.
- Lemma 8.6, pp952--953: with two even-order zeros, one
  positive residue and all other residues negative, the
  components are the two parities of cyclic orders of the
  poles around the one-zero boundary.
- Section 8.1.2, p947: after transport to an appropriate
  nonresonant residue configuration the separating invariant
  is a surgical Arf invariant, independent of zero labels.

The downloaded published edition is retained and pinned.
The component theorem is an explicit imported mathematical
input. The passage to polynomial passports, the label
quotients, and the monodromy proof are supplied here.
The [two-tripod argument](POLYNOMIAL_TWO_TRIPOD_COMPONENTS.md)
gives an independent elementary proof when \(s=2\), \(q_1=q_2=3\).
No novelty or historical priority is asserted.

## Lemma F1. Logarithmic derivative correspondence

Put \(a_i=q_i-1\) and take the coefficient-residue vector
\[
                 \rho=(\ell_1,\ldots,\ell_k,-N).              \tag{3}
\]
Polynomials with a labeled zero fiber of type \(\lambda\)
and labeled remaining critical points of orders \(a_i\)
correspond, up to multiplication of the polynomial by a
nonzero constant, to the isoresidual fiber
\[
             \mathcal H(a_1,\ldots,a_s,-1^{k+1})_\rho.         \tag{4}
\]
Source changes of coordinates on the polynomial side match
biholomorphisms on the differential side. The locus where
the remaining critical values are distinct is a nonempty
open subset of each component, with the same components.

**Proof.** Move the pole of residue \(-N\) to infinity.
For roots \(x_j\), the polynomial
\[
                 f(z)=C\prod_{j=1}^k(z-x_j)^{\ell_j}
\]
satisfies
\[
                    \frac{df}{f}
                        =\sum_{j=1}^k\frac{\ell_j\,dz}{z-x_j}.
                                                               \tag{5}
\]
It has simple poles with residues (3). At a nonzero finite
critical value a zero of order \(a_i\) of \(df/f\) is a
ramification point of local degree \(a_i+1=q_i\).
The degree of its zero divisor is \(k-1=\sum_i a_i\).

Conversely a differential on the sphere with only those
simple poles is exactly the right side of (5): subtraction
leaves a holomorphic one-form on the sphere and hence zero.
Integrating gives the displayed polynomial, uniquely up
to \(C\in\mathbb C^*\), since its residues are the stated
positive integers and their negative sum. The distinguished
negative-residue pole is unique, so it introduces no choice
of infinity. This proves the correspondence.

No proper subset of (3) sums to zero. A subset omitting
\(-N\) has positive sum. A proper subset containing it has
negative sum. The primary theorem's nonresonance hypothesis
therefore holds even when positive residues repeat or have
a common divisor. Its convention uses periods rather than
coefficient residues: use \(2\pi i\rho\), or equivalently
the differential \((2\pi i)^{-1}df/f\) with vector \(\rho\).
Uniform nonzero rescaling does not change the component
or nonresonance assertion.

Outside resonance, relative periods between the \(s\)
zeros are local coordinates on the fixed-residue fiber.
Locally the ratio of two critical values is the exponential
of the corresponding relative period. Thus equality of
two such values is a proper analytic condition and cannot
hold identically on a component. It is also an algebraic
condition, using the product formula for \(f\).
The fiber is smooth in these period coordinates; its
connected components are irreducible. Deleting these
finitely many proper complex algebraic loci leaves each
component nonempty and connected. Multiplication by
\(\mathbb C^*\) and source affine-coordinate changes have
connected fibers. This proves the final assertion.
\(\square\)

## Theorem F2. Complete component criterion for the single-cycle tail

For (1)--(2), the number of full braid components is
\[
 \#\pi_0=
 \begin{cases}
   2,&\text{all }q_i\text{ odd and all }\ell_j\text{ distinct},\\
   1,&\text{otherwise}.
 \end{cases}                                                 \tag{6}
\]
When the first line applies, the parity of a one-zero
boundary cyclic order is a complete component label.
Transported surgical Arf parity is an equivalent geometric
label. It is not identified here with a reduced Schur lift
in the actual finite monodromy group.

**Proof.** Apply the primary Theorem 1.7 to (4).
There are \(s\ge2\) zeros and \(k+1\ge4\) simple poles.
The two-simple-pole exceptional family is inapplicable.
The other exceptional family occurs exactly when all
\(a_i\) are even, equivalently all \(q_i\) are odd.
Lemma F1 identifies these components with polynomial
components having labeled roots and labeled tail critical
points, after deleting critical-value collisions.

Suppose all \(a_i\) are even. Multiply the differential
by \(-1\), so that infinity has the unique positive
residue \(N\) and the roots have negative residues.
Merge the zeros into two nonempty groups. Their total
orders are even. By the primary Lemma 8.1 and proof of
Theorem 1.7, the two components of this two-zero boundary
are adjacent to the two distinct components of (4).
Lemma 8.6 reads their labels by the cyclic order of the
\(k\) finite poles in the further one-zero boundary.
Here \(k=1+\sum a_i\) is odd. The parity of a list of
these \(k\) labels is therefore unchanged by cyclic
rotation of the list. Every cyclic order occurs at the
one-zero boundary. Two orders have the same adjacent
component exactly when their list parities agree.

Permuting equal-order zero labels does not change the
surgical Arf invariant: the surgery and its quadratic
form do not use these labels. The component-preserving
action can equally be transported back from the residue
configuration in Section 8.1.2, because residue transport
commutes with zero relabeling. Thus forgetting the labels
on equal-length tail cycles preserves the two components.

If two parts of \(\lambda\) agree, their pole labels may
be interchanged while the residue vector stays fixed.
On the one-zero boundary this transposition reverses
list parity. Boundary adjacency is unique and equivariant
under pole relabeling, so the transposition interchanges
the two components of (4). Forgetting equal-residue root
labels therefore merges them. If the parts are distinct,
their lengths identify the roots, and no such quotient
is present. These arguments prove (6).

If \(\lambda\) is itself a single-cycle type occurring
among the tail factors, temporarily distinguish that
branch value before applying F1. Such a partition has
repeated parts: in the only possible strict one-cycle
case it would have at most two parts, whereas \(k\ge3\).
The distinguished space is connected, and forgetting
which equal-type branch was distinguished leaves it
connected. Finally the polynomial/cactus correspondence
identifies these connected components with full braid
orbits for the stipulated type placements.
\(\square\)

## Theorem F3. Exact monodromy, every block system, and trivial deck

Set \(n=N/D\), let \(B_j=\{j,j+D,\ldots,j+(n-1)D\}\), and put
\(K=\prod_{j\in\mathbb Z/D}A(B_j)\). Every tuple in (2) has
block-count divisors exactly
\[
                         \operatorname{div}(D)\cup\{N\}.
                                                               \tag{7}
\]
Its actual embedded group is
\[
 G=
 \begin{cases}
     K\langle P_N\rangle,&\text{all }q_i\text{ odd},\\
     \bigl(\prod_j S(B_j)\bigr)\rtimes C_D,
                    &\text{at least one }q_i\text{ even}.
 \end{cases}                                                 \tag{8}
\]
In the second line \(C_D\) is the induced cyclic block action;
the expression denotes the full cyclic-top wreath group.
Its orders are respectively
\[
 (n!/2)^D D
       \begin{cases}1,&n\text{ odd},\\2,&n\text{ even},\end{cases}
       \qquad (n!)^D D.                                      \tag{9}
\]
The deck group is trivial. In particular, both components
in the split case have the same literal subgroup of \(S_N\).

**Proof.** Use the
[quotient and return index budget](POLYNOMIAL_PASSPORT_COMPONENTS.md#theorem-h3-the-quotient-and-fibre-index-budget).
For a proper block system with block size \(b>1\), total
return index is \(b-1\). A single-cycle factor acting
nontrivially on blocks has exactly one nontrivial quotient
cycle, covers all sheets in those moving blocks, and has
a full \(b\)-cycle as its return. It alone consumes \(b-1\).
There is at least one other tail factor. That other
nonidentity factor contributes positive return index
whether its quotient is trivial or nontrivial. This is
impossible. Hence every tail factor fixes all blocks
and has its support inside one block.

Only \(g\), the factor of type \(\lambda\), acts on the
quotient. Its quotient is a full cycle, so the block count
\(m\) divides every \(\ell_j\), and hence \(D\).
Conversely contract the \(g\)-cycles in the incidence
tree. The resulting bipartite tree has \(k\) cycle
vertices and \(s\) tail vertices of degrees \(q_i\).
For any \(m\mid D\), give each \(g\)-cycle a phase
in \(\mathbb Z/m\) increasing by one under \(g\).
Start on one cycle and propagate the phase origins
through the tree so that every tail support has a
common phase. A tree gives no compatibility loop.
This produces equally sized phase blocks preserved
by every factor, proving (7). In fixed-product
coordinates they are the residue partitions.

There are two tail vertices sharing a cycle vertex
along a path in the contracted tree. Their sets of
\(g\)-cycle neighbors meet in exactly one cycle.
Conjugate one tail by a power of \(g\) so that their
attachments on that cycle agree. The two single-cycle
supports then meet in exactly one sheet. Their
commutator is a literal three-cycle (on the shared
sheet and its two images). In particular this
construction works even when no prescribed factor
is a three-cycle.

Let \(L\) be the normal closure of this three-cycle
in \(G\). Its support-hypergraph components are
equal-size \(G\)-blocks. Connected three-cycle
supports generate the alternating group on their
union, by the elementary \(A_4/A_5\) overlap proof
in [E3](POLYNOMIAL_TWO_TRIPOD_COMPONENTS.md#theorem-e3-actual-group-every-block-system-and-deck-symmetry).
If \(L\) has \(h\) orbits, (7) gives \(h\mid D\);
if \(h=1\) this is immediate. Since \(L\) fixes
every \(D\)-block, also \(h\ge D\). Thus \(h=D\)
and \(L=K\).

Each tail lies in one \(D\)-block. If all its cycles
have odd length, every tail is even there and
belongs to \(K\). The full product yields
\(G=K\langle P_N\rangle\). The intersection with
\(\langle P_N\rangle\) gives the first formula (9)
as in E3. If a tail has even length, it supplies
an odd permutation in one block. Together with
the alternating factor it gives the full symmetric
group there. Conjugation by \(P_N\) supplies the
independent full symmetric group on every block,
giving the second line of (8).

In the all-odd case \(n\ge k\ge5\), so the natural
centralizer of each alternating block factor is
trivial. In the other case \(n\ge k\ge3\), and
the same is true for each symmetric block factor.
Independence of the factors prevents a centralizer
from interchanging blocks. Thus the full centralizer,
and hence the deck group, is trivial. \(\square\)

## Exact counts and the remaining lifting question

Let \(m_a\) be the multiplicity of part \(a\) in
\(\lambda\), set \(M=\prod_a m_a!\), and let \(L_0\)
be the number of distinct placements of all finite
cycle types, including any equality with \(\lambda\).
The total number of literal tuples in (2) is
\[
             L_0\,N^s\,\frac{(k-1)!}{M},                    \tag{10}
\]
by the Goulden--Jackson formula in
[Jones--Zvonkin Section 1.3](sources/jones-zvonkin-braid-cacti-2002.pdf).
All deck groups are trivial, so the total number of
absolute classes is (10) divided by \(N\).
In every connected case these are the exact orbit sizes.
When (6) gives two components, (10) gives their sum.
No uniform equality-of-sizes formula is inferred from
the two-tripod specialization.

The surgical parity supplies a sufficient component
label in the split case. It is not automatically
the standard odd-order spin product, which can be
forced by the genus-zero passport, or a class-defined
lift in the actual group (8). The degree-fifteen
case retains its separately proved genuine reduced
quadratic invariant. Uniform identification of finite
lifting data for every member of (2) remains open.

## Evidence and assurance

The [checker](verify_polynomial_single_cycle_tail.py),
[receipt](polynomial_single_cycle_tail_verification.json)
and [source ledger](polynomial_single_cycle_tail_source_ledger.json)
retain exact small passport closures, independent
ordinary character totals, all residue block-divisor
checks, normal-three-cycle support certificates,
commutator witnesses, logarithmic derivative examples
and deliberately invalid inputs.
The unrestricted component count depends on the
explicit primary theorem, together with F1--F2's
bridge and label quotient arguments. Finite examples
are supporting controls, not substitutes for that
theorem. The unrestricted group statement has the
written proof F3. Review and replay are internal;
independent replay, external review and formal
verification are not asserted.
