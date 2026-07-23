# Restricted Hessian synchronization: the proved normal-form theorem and the square-free gap

Work over an algebraically closed field \(k\) of characteristic zero.  This
note audits the proposed restricted theorem:

> If the normalized relation graph has only singleton blocks and at most one
> tame power/Dickson strongly connected block with square-free degree
> multiset, then all canonical Hessian linear lifts agree modulo the summed
> Hessian residual ideal.

The conclusion must be separated into two levels:

\[
\begin{array}{rcl}
\text{generic synchronization}
&\Longleftrightarrow&
\lambda_{\mathbf d}-\lambda_{\mathbf e}\in\sqrt{H_D},\\[2mm]
\text{scheme-theoretic synchronization}
&\Longleftrightarrow&
\lambda_{\mathbf d}-\lambda_{\mathbf e}\in H_D.
\end{array}                                                   \tag{0.1}
\]

The relation-graph normal forms prove the first statement once component
completeness is known.  They do not, by themselves, prove either component
completeness or the second statement.  In particular, the proposed
restriction contains all five presently open degree-thirty pair
memberships.  Thus it is a good theorem target, but it is not yet a theorem
of the repository.

## 1. Canonical setup

Put

\[
 R_N=k[c_2,\ldots,c_{N-1}],\qquad
 F=W^N+\sum_{j=2}^{N-1}c_jW^j.
\]

For a proper outer cut \(a\mid N\), let \(H_a\subset R_N\) be the canonical
Hessian residual ideal for a normalized \(a\circ(N/a)\) decomposition.
Triangular reconstruction uses only the displayed coefficients and produces
a regular function

\[
 \lambda_a\in R_N
\]

whose value is the forced coefficient of \(W\) in the corresponding monic
original polynomial lift.  For a collection \(D\) of requested cuts or
words, write

\[
 H_D=\sum_{\mathbf d\in D}H_{\mathbf d}.
\]

The full synchronized polynomial ideal is

\[
 H_D+(\lambda_{\mathbf d}-c_1:\mathbf d\in D).                \tag{1.1}
\]

Consequently equality of all lifts modulo \(H_D\) is exactly the assertion
that (1.1) is the graph of one regular function over the entire, possibly
nonreduced, Hessian intersection.

## 2. What the power/Dickson charts prove formally

The useful theorem available without any further elimination is the
following.

> **Normal-form vanishing theorem.**  Let \(X\) be a reduced irreducible
> \(k\)-scheme and let
> \[
>  f_X=W^N+\sum_{j=1}^{N-1}u_jW^j\in\Gamma(X,\mathcal O_X)[W]
> \]
> be one monic original polynomial family.  Suppose that, for every
> requested word \(\mathbf d\in D\), the same \(f_X\) has a normalized
> decomposition of degree word \(\mathbf d\).  Let
> \(\phi_X:X\to\operatorname{Spec}R_N/H_D\) retain the coefficients of
> degrees at least two.  Then
> \[
>  \phi_X^*(\lambda_{\mathbf d}-\lambda_{\mathbf e})=0
>  \qquad(\mathbf d,\mathbf e\in D).                           \tag{2.1}
> \]

**Proof.**  Apply canonical triangular reconstruction to any one of the
given normalized decompositions of \(f_X\).  Uniqueness of the reconstruction
forces its reconstructed linear coefficient to be \(u_1\).  This is true for
every requested word, so both terms in (2.1) pull back to \(u_1\).  \(\square\)

This theorem applies integrally to every synchronized exponential or
Dickson family supplied by the tame multi-collision normal form.  It also
survives arbitrary common outer and inner composition decorations, because
the chart starts with one polynomial family carrying all requested
decompositions.

The precise dense-component consequence is:

> **Reduced transfer criterion.**  Suppose every minimal component of
> \(V(H_D)\) is dominated by a synchronized power/Dickson relation-graph
> chart.  Then
> \[
>  \lambda_{\mathbf d}-\lambda_{\mathbf e}\in\sqrt{H_D}
> \]
> for every requested pair.  If \(H_D\) is radical, then the stronger
> membership
> \[
>  \lambda_{\mathbf d}-\lambda_{\mathbf e}\in H_D             \tag{2.2}
> \]
> follows.

**Proof.**  Equation (2.1) gives vanishing at the generic point of every
minimal component, hence membership in every minimal prime and therefore in
\(\sqrt{H_D}\).  If \(H_D=\sqrt{H_D}\), this is (2.2).  \(\square\)

This is the strongest restricted synchronization theorem that follows
formally from the current relation-graph theory.  It makes the two missing
inputs explicit:

1. **component completeness:** no extra Hessian-projection component misses
   all synchronized normal-form charts;
2. **nilpotent control:** the lift difference vanishes in each primary
   thickening, not only on its reduced support.

## 3. Why Ziegler's theorem does not supply the missing inputs

Ziegler's tame multi-collision theorem starts with one polynomial

\[
 f\in\bigcap_{\mathbf d\in D}\mathcal D_{\mathbf d}.
\]

Its relation graph classifies the decompositions of that already
synchronized \(f\).  The compatibility of the original-shift parameters in
the proof is likewise deduced from decompositions of one \(f\).  A point of
the Hessian intersection gives instead a collection

\[
 F+\lambda_{\mathbf d}W\in\mathcal D_{\mathbf d},
\]

with no initial equality among the \(\lambda_{\mathbf d}\).  Applying the
ordinary theorem before proving their equality is circular.

The square-free assumption simplifies the gcd-refinement step in the
ordinary collision theorem.  It does not change this logical boundary:
forgetting the linear coefficient still need not commute with intersection.
Nor does a valuation at infinity distinguish the missing linear terms
directly, since all of them occur \(N-1\) orders below the leading term.
An infinity argument must therefore be a new lemma about two distinct
decomposable lifts whose difference is linear, not a quotation of the
ordinary multi-collision theorem.

There is also a strict scheme-theoretic obstruction to concluding ideal
membership from reduced normal forms.  The map

\[
 k[\epsilon]/(\epsilon^2)\longrightarrow k,\qquad\epsilon\longmapsto0,
\]

kills \(\epsilon\) on the reduced chart although
\(\epsilon\notin(\epsilon^2)\).  The degree-thirty braid calculations show
that analogous nonreduced path thickenings genuinely occur in the Ritt
geometry.  Relation-graph combinatorics alone cannot rule them out.

## 4. The square-free degree-thirty audit

The proposed restriction does not leave degree thirty as a distant braid
test.  It contains every currently open degree-thirty pair.

| open cuts | normalized nontrivial block | remaining singleton decoration |
|---|---|---|
| \(\{2,3\}\) | \(\{2,3\}\) | common inner degree \(5\) |
| \(\{2,5\}\) | \(\{2,5\}\) | common inner degree \(3\) |
| \(\{3,5\}\) | \(\{3,5\}\) | common inner degree \(2\) |
| \(\{2,15\}\) | \(\{2,3,5\}\), with words \((2,3,5)\) and \((3,5,2)\) | none |
| \(\{3,10\}\) | \(\{2,3,5\}\), with words \((3,2,5)\) and \((2,5,3)\) | none |

Every displayed multiset is square-free and pairwise coprime.  In the last
two rows the union relation graph is strongly connected: the moved vertex
has two-way relations with both other vertices, while the other two retain
their order.

Thus the proposed theorem would prove exactly the five memberships

\[
\begin{aligned}
\lambda_2-\lambda_3&\in H_2+H_3,&
\lambda_2-\lambda_5&\in H_2+H_5,\\
\lambda_3-\lambda_5&\in H_3+H_5,&
\lambda_2-\lambda_{15}&\in H_2+H_{15},\\
\lambda_3-\lambda_{10}&\in H_3+H_{10}.             \tag{4.1}
\end{aligned}
\]

These are precisely the five pairs left open by the permanent exact
degree-thirty audit.  Degrees eighteen and twenty-four are not open tests:
all their pair memberships have already been proved over characteristic
zero.  The all-six degree-thirty intersection is also synchronized by the
certified spanning tree

\[
 2-6-3-15-5-10,
\]

but that argument does not imply any of the five smaller-ideal memberships
in (4.1).

## 5. A noncircular proof target

A usable local theorem would have to begin before synchronization:

> **Primitive Hessian edge lemma (needed).**  Let \(r,s\ge2\) be coprime and
> let
> \[
> P=A\circ B,\qquad Q=C\circ D
> \]
> be monic original polynomials of degree \(rs\), with outer degrees \(r\)
> and \(s\), respectively.  If \(P''=Q''\), then \(P=Q\), at least at the
> generic points of every tame power/Dickson component; for
> scheme-theoretic synchronization the same assertion must hold in the
> corresponding infinitesimal families.

Since \(P-Q=\delta W\), the differentiated equation is

\[
 A'(B)B'-C'(D)D'=\delta.                              \tag{5.1}
\]

Ritt's theorem classifies (5.1) only after \(\delta=0\) has been proved.
The proposed ramification/infinity strategy must show that a nonzero
constant in (5.1) is incompatible with the two decomposition structures.
For transported rows of the table, it must also prove that stripping a
common decoration is legitimate for the Hessian equation; this is not
formal because second derivatives do not commute with composition.

The minimal successful proof package is therefore:

1. prove the primitive edge lemma on the exponential and Dickson generic
   points;
2. prove a transport lemma through arbitrary singleton blocks;
3. prove that the resulting charts exhaust every minimal prime of the
   Hessian intersection;
4. lift generic vanishing through the primary structures, or construct an
   explicit syzygy for each difference in (4.1).

Steps 1--3 prove generic synchronization.  Step 4 is additionally necessary
for the ideal-membership statement.

## 6. Exact falsification and certification protocol

For each pair in (4.1), pull the target residuals to the canonical factor
chart of the source cut and let \(I\) be the resulting exact ideal and
\(\Delta\) the pulled-back lift difference.

1. Saturate by the exact-degree and reconstruction denominators.
2. Compute modular minimal primes at several good primes.
3. Match their dimensions and degrees with the transported exponential and
   Dickson graph ideals.
4. Lift the candidate primes to characteristic zero and prove both ideal
   containments.
5. Compute a characteristic-zero primary decomposition
   \(I=\bigcap Q_i\), or an equivalent certified primary filtration.
6. Reduce \(\Delta\) modulo every \(Q_i\).

A single \(Q_i\) with \(\Delta\notin\sqrt{Q_i}\) is a genuinely
unsynchronized component.  If \(\Delta\in\sqrt{Q_i}\setminus Q_i\), generic
synchronization holds there but the claimed ideal membership fails because
of a nilpotent defect.  If \(\Delta\in Q_i\) for every \(i\), then
\(\Delta\in I\).

For discovery, a bounded-degree modular syzygy

\[
 \Delta=\sum_j q_jR_j
\]

is likely cheaper than a full primary decomposition.  Rational
reconstruction followed by direct expansion over \(\mathbb Q\) gives the
smallest permanent certificate.

## 7. Status

The rigorous conclusion is:

> Power/Dickson normal forms force equality of the canonical linear lifts on
> every synchronized normal-form component.  If those components exhaust the
> minimal primes, the restricted Hessian intersection is generically
> synchronized.  Ideal membership additionally requires radicality or
> primary-thickening control.  The square-free graph hypothesis alone does
> not currently provide either input, and in degree thirty it specializes to
> all five open pair certificates.

The ordinary normal-form source is Konstantin Ziegler,
[Tame Decompositions and Collisions](https://arxiv.org/abs/1402.5945),
especially the relation-graph splitting and strongly connected block
classification.  Its complete-decomposition background is Peter Müller and
Michael Zieve,
[On Ritt's Polynomial Decomposition Theorems](https://arxiv.org/abs/0807.3578).
