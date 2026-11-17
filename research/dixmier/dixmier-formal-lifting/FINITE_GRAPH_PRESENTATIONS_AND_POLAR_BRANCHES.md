# Finite graph presentations and parameter-polar branches

For an actual polynomial symplectic graph, the h-complete graph module is
already free over the completed source Weyl algebra when its classical
bundle is trivial. The remaining effectivity question is an uncompleted
source lattice stable under all target actions. A finite monic module
reduction certificate can prove that property for a supplied presentation.
It does not construct such a presentation or bound the degree needed to
find one. Coherence and h-flatness alone are insufficient, even in a
commutative example with identity-graph reduction.

## Which classical graph is available

The canonical full polynomial Long frame is retained in
[the exact transcription](../dixmier-audit/verify_long_pbw.py), function
`long_images()`. Its source coordinates are \((X,Q,P,Z)\), and its
ordered target coordinates are \((R,T,D,S)\), with
\(\{D,R\}=\{S,T\}=1\). That script constructs all four polynomials,
including the corrected fourth image of degree 23 and 139 terms, and
checks every one of the six Poisson brackets. A fresh internal replay on
17 November 2026 passed all six checks, the three-point collision and
the complete support-size controls. This is a replay of an existing
classical input, not a new classical completion theorem.

The [centralizer normalization](CENTRALIZER_NORMAL_FORM.md) uses that
admitted fourth image. The [quintic phase chart](QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md),
equation (3), identifies its generic expression. Thus the canonical full
Long map is an eligible classical input to a graph theorem.

By contrast, [the residue-transport inputs](../dixmier-residue-transport/inputs.json)
contain a pair and selected pair corrections; they do not contain a
four-image frame, and their coefficients differ from the untransformed
canonical transcription. Their identification with a full classical
frame would require a checked transport or an explicit fourth image.
[QFL1](PROOF.md) does not infer such a completion for an arbitrary
retained fiber pair or a prescribed pair-only quantum jet. Every graph
argument below assumes an actual polynomial symplectic map, rather than
silently promoting a partial tuple to one.

## Formal source freeness follows from the actual graph

Let \(k\) have characteristic zero, let \(V=\mathbb A_k^{2n}\), and
let \(\phi:V\to V\) be polynomial symplectic. In target times opposite
source, its graph \(j:Y\hookrightarrow\mathcal X\) is affine and its
source projection is an isomorphism. Let \(E_h\) be an h-flat,
h-complete module over the product formal Weyl quantization, with
\[
                         E_h/hE_h\simeq j_*\mathcal O_Y^r.
                                                               \tag{1}
\]
For example, the hypotheses in
[the audited formal graph theorem](NEW_FORMAL_GRAPH_CARRIERS_AND_FINITE_ADMISSION_GAP.md)
give such a module. Write \(A_s\) for the finite-PBW source algebra over
\(k[[h]]\), \(\widehat A_s\) for its h-completion, and
\(M=\Gamma(\mathcal X,E_h)\).

Then
\[
                         M\simeq\widehat A_s^{\,r}
                   \quad\text{as complete right modules}.       \tag{2}
\]
Here is a direct proof, including the global-section step. Flatness
identifies every layer \(h^jE_h/h^{j+1}E_h\) with the coherent affine
sheaf in (1). Its positive cohomology vanishes. The short exact
sequences between the finite h-quotients therefore give surjective
maps of global sections, and reduction at the first level is
\[
                         M/hM\simeq k[V]^r.                    \tag{3}
\]
More explicitly, each global section at a finite level lifts to the
next because the intervening layer has zero first cohomology. Choose
compatible lifts; completeness identifies their inverse limit with
global sections of \(E_h\). Flatness identifies the kernel at the
first level with \(hM\), justifying (3).

Lift the standard basis of (3) to \(m_1,\ldots,m_r\in M\). Source
right multiplication defines a continuous map
\(\widehat A_s^r\to M\). Successive reduction and subtraction, followed
by completeness, prove surjectivity. If a nonzero tuple were in the
kernel, divide by its least common h-power. Flatness permits division
of the zero equality, and its reduction would be a nontrivial relation
between the basis vectors in (3). This is impossible. Separation
ensures that a nonzero tuple has a finite first nonzero h-order.
This proves (2).

The resulting target-action matrices belong to
\(M_r(\widehat A_s)\). Their coefficients may have arbitrarily growing
PBW degree. Formula (2) supplies neither a finite free \(A_s\)-lattice
nor matrices in \(M_r(A_s)\).

## A sufficient finite monic presentation certificate

For a concrete algebraic candidate take
\(R_0=k[h]_{(h)}\). Let \(A_s,A_t\) be the source and target finite
Weyl algebras over \(R_0\), and put
\[
                       B=A_t\otimes_{R_0}A_s^{\rm op}.           \tag{4}
\]
A left \(B\)-module is a target-left, source-right bimodule. The
opposite source Weyl relations have the opposite sign; the two factors
commute. Suppose a module is supplied by finitely many generators and
relations over \(B\).

The following finite certificate is sufficient.

1. Supply monic oriented rules for its relations, together with the
   ambient PBW rules, and a compatible well-founded monomial/module
   order for which every right side is smaller than its left side.
2. Supply finite reductions resolving every overlap and inclusion
   ambiguity. No nonunit of \(R_0\), especially \(h\), is divided out.
3. Show that the irreducible module words are exactly
\[
 \{\text{source-opposite PBW word}\}\,\mathcal T,
                   \qquad \mathcal T=\{\tau_1,\ldots,\tau_r\},   \tag{5}
\]
   where the \(\tau_i\) are a fixed finite set of target/module tails.

Condition (5) is not an infinite computational sample. For a PBW module
Groebner certificate it can be established from the finite initial
monomial module: its module generators involve only target monomials,
and their target complement is the listed finite set. Pure-power
bounds for each target variable give a finite box in which that
complement can be checked. An alternative finite rewriting certificate
must prove the same entire normal-word description.

For completeness, the normal-form implication needs no parameter
specialization theorem. Termination follows from the supplied order.
If two first reductions act on disjoint places, they commute. If they
overlap or include one another, the supplied ambiguity reduction
resolves their difference using smaller words. Induction on the order
therefore gives a unique normal form. Monicity makes this argument
valid over \(R_0\), without a field-division step. Normal words form
an \(R_0\)-basis of the quotient. In view of (5), its underlying right
source module is exactly
\[
                              M_0\simeq A_s^r.                 \tag{6}
\]
Independence here is as important as generation: no source polynomial
relation between the tails remains.

Reduce each target generator times each \(\tau_i\) to its finite
normal form. Grouping source words gives actual matrices
\(H_j\in M_r(A_s)\). Since they are induced by the \(A_t\)-action,
they satisfy every target Weyl relation. Their finite list of entries
has a computable common Bernstein PBW bound, namely the maximum
source degree in these normal forms. All their scalar coefficients
are regular at \(h=0\).

For the desired graph, additionally check
\[
                         \overline H_j=\phi_j I_r              \tag{7}
\]
in this source basis, for every target generator. This verifies the
classical graph, not a union of distinct classical branches. If the
claim is algebraization of one particular formal graph module, also
give a completed isomorphism or verify its full chosen formal
projective and scalar connection data. Graph reduction alone does not
identify that formal isomorphism class.

This is a sufficient certificate over an effective coefficient ring.
It is not asserted that rational functions of \(h\) cover every
bounded lift with arbitrary \(k[[h]]\) coefficients. Explicitly
presented algebraic-series coefficient rings can also be used when
their equality and unit certificates are supplied. No completeness of
a rational-parameter coefficient search follows.

An equivalent abstract filtration formulation is often less useful:
an uncompleted, target-stable source module with an exhaustive
nonnegative good source filtration and finite homogeneous free source
associated graded admits finite free source generators by leading-term
lifting and subtraction.
Finite filtration degrees for the finitely many target generators
then give finite PBW matrices. Conversely a finite free source module
with such matrices has its standard source filtration. Merely asking
for the existence of this structure restates the finite-free carrier
problem. A supplied finite presentation and the certificate above
make it a verifiable assertion about one candidate.

## Coherence and h-flatness do not remove polar branches

Already in commutative algebra, consider
\[
 M_0=R_0[x,t]/(ht^2+t-x)\simeq R_0[t],
                  \qquad x\text{ acts as }t+ht^2.               \tag{8}
\]
It is cyclic and coherently presented over \(R_0[x,t]\), and it is
h-flat. Its reduction is the identity graph \(t=x\). The map
\(x\mapsto t+ht^2\) has a unique inverse in the h-complete polynomial
ring, obtained order by order starting with \(t=x\). Thus its
completion is free of rank one over the completed source ring.

Nevertheless \(M_0\) is not finite over the uncompleted \(R_0[x]\).
Its quotient by \(x\) maps onto \(k(h)\) by
\[
                            t\longmapsto-1/h.
\]
If \(M_0\) were source finite, this quotient would be a finite
\(R_0\)-module. But multiplication by \(h\) is surjective on
\(k(h)\), which a nonzero finite module over the local ring \(R_0\)
cannot satisfy: the finite Nakayama argument would force it to be
zero. This is the extra parameter-polar branch invisible to the
selected completion. It is a commutative countercontrol to an
algebraization inference, not a constructed Weyl bimodule.

## Curvature and common triangular flags

There is a separate failure before the finite-source test for a naive
graph presentation. Let \(f_j\) be an actual classical symplectic frame,
let \(F_j\) be chosen uncorrected finite Weyl representatives, and in
(4) impose on a vector \(e\) the constraints
\[
 g_j e=0,\qquad
 g_j=y_j\otimes1-1\otimes F_j^{\rm op}.                         \tag{9}
\]
With target Poisson matrix \(J_{ij}\), their commutator is the
opposite-source curvature
\[
 [g_i,g_j]=-C_{ij}^{\rm op},\qquad
                         C_{ij}=[F_i,F_j]-hJ_{ij}.              \tag{10}
\]
Hence \(eC_{ij}=0\). If
\(C_{ij}=h^m(E_{ij}+O(h))\), \(\overline E_{ij}\ne0\), h-flatness
would give \(\overline e\,\overline E_{ij}=0\). This is impossible
when \(\overline e\ne0\) lies in the free classical graph module
over the polynomial source domain. In Weyl ordering the first defect
of the canonical Long frame is already nonzero. The
[independent Fraction control](check_long_graph_curvature.py) imports
only `long_images()` as its input and recomputes the polynomial algebra.
Its [retained coefficient receipt](long_graph_curvature_controls.json)
gives 30 terms in each of
\(E_{ST}=\Pi^3(S,T)/24\) and \(E_{DR}=\Pi^3(D,R)/24\), with
\[
 E_{ST}(1,0,0,0)=199/2,\qquad E_{DR}(1,0,0,0)=-27/2.
\]
Here \(C=h^3E+O(h^5)\); higher terms are not asserted zero.
Thus these unmodified graph constraints cannot be an h-flat
quantization of its graph.

Adding higher Koszul or BRST terms while retaining these exact induced
degree-zero constraints cannot remove (10). A successful construction
must change the induced target action or its constraints. Exterior
nilpotence by itself does not cancel a nonzero scalar curvature on a
free classical graph module. This argument does not address arbitrary
BRST constructions with different degree-zero actions.

Likewise, suppose actual finite target matrices are upper triangular
in a common admitted source-free frame. Each diagonal tuple separately
satisfies all scalar Weyl relations, because the diagonal of a product
of upper triangular matrices is the product of their diagonal entries.
Keeping the uncorrected \(F_j\) on a diagonal retains its curvature;
correcting that diagonal already supplies a finite scalar lift. With
scalar graph reduction (7), its reduction is the desired classical
map. Strict triangular or nilpotent corrections preserving this flag
therefore do not bypass the scalar finite gate. This is not a claim
that arbitrary creation/annihilation matrices, projective carriers or
matrix actions have such a common admitted flag.

## Exact remaining task

For canonical Long the classical full frame is an established input,
and formal source freeness follows from an actual formal graph module.
Neither fact constructs a finite presentation satisfying (5)--(7).
If one starts with unknown finite action matrices and writes only the
linear rules for their actions on a prescribed source basis, resolving
their ambiguities is exactly the bounded matrix Weyl ansatz again.
The criterion provides a finite verifier for another supplied module
presentation, not a new existence result or an a priori degree bound.

This note constructs no bounded Long lift, algebraization of the
nonconstant graph carrier, scalar realization or DC1/DC2 result.
Formal verification, independent whole-proof replay and external review
are not asserted. The finite certificate proves properties of a supplied
presentation; it is not an existence theorem for that presentation.
