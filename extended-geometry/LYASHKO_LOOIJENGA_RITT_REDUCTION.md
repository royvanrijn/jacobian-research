# Lyashko--Looijenga profiles and the Ritt component reduction

Work over an algebraically closed field of characteristic zero.  This note
connects the polynomial Lyashko--Looijenga branch data to the existing
Hessian--Ritt relation graph.  Its main concrete consequence is that the
degree-thirty all-cut intersection has no additional reduced component away
from the six complete-decomposition charts.

The result concerns reduced component completeness.  It does not identify
the different nilpotent path structures around a Ritt two-cell; those remain
the subject of the Hessian--Ritt deformation complex.

## 1. The generic LL profile of a factor word

Let

\[
 f=f_1\circ\cdots\circ f_r,\qquad \deg f_i=d_i,
 \qquad N=\prod_i d_i,
\]

and put

\[
 s_i=\prod_{j>i}d_j.
\]

On the open where every factor is simply ramified, all displayed critical
values are distinct, and no critical value from one factor level coincides
with one from another, a critical point of \(f_i\) has \(s_i\) inverse images
under the inner suffix.  These are simple critical points of \(f\), and they
have one common critical value.  The corresponding branch cycle is a product
of \(s_i\) disjoint transpositions.

Therefore the generic LL branch-cycle signature is

\[
 \boxed{\prod_{i=1}^r (2^{s_i})^{\,d_i-1}.}          \tag{1.1}
\]

Here \((2^s)^m\) means \(m\) distinct finite branch values, each with local
monodromy a product of \(s\) disjoint transpositions.  In particular, for a
two-factor composition of degrees \((a,b)\),

\[
 \boxed{(2^b)^{a-1}(2)^{b-1}.}                      \tag{1.2}
\]

The Riemann--Hurwitz count is automatic:

\[
 \sum_i(d_i-1)s_i=N-1.                              \tag{1.3}
\]

The number of finite branch values is

\[
 \sum_i(d_i-1).
\]

After the one-dimensional target scaling used in the marked seed chart, the
LL stratum has dimension

\[
 \sum_i(d_i-1)-1,                                   \tag{1.4}
\]

exactly the normalized factor-chart dimension.  For \((a,b)\), subtracting
(1.4) from the ambient seed dimension \(N-3\) gives

\[
 \boxed{\operatorname{codim}\mathcal C_{a,b}=N-a-b.} \tag{1.5}
\]

Thus the codimension part of the Hessian atlas is already visible in the
generic LL collision profile.

## 2. Block systems and prime refinements

A decomposition \(f=A\circ B\) is equivalently an intermediate field

\[
 k(f)\subset k(B)\subset k(W),
\]

or a block system of the polynomial monodromy action.  The block size is
\(\deg B\).  Taking gcd refinements of several block systems produces prime
factor words; equal prime occurrences retain their relative order, while
coprime adjacent occurrences are related by Ritt moves.

For a requested outer cut \(a\mid N\), call a prime word a **carrier** if a
prefix has product \(a\).  A polynomial carrying several decompositions
chooses at least one carrier for every requested cut.  The union of their
order relations is precisely the normalized relation graph used by tame
multi-collision theory.

This observation supplies the missing coverage statement: complete
decomposition charts are not merely test charts.  After synchronization,
every point of the multiple intersection admits a prime refinement and hence
lies in one of the carrier systems enumerated below.

## 3. Degree twelve

The prime multiset is \((2_1,2_2,3)\), with \(2_1\) preceding \(2_2\).
The three prime words are

\[
 2_1\,2_2\,3,\qquad 2_1\,3\,2_2,\qquad 3\,2_1\,2_2.
\]

They carry the four proper cuts \(2,3,4,6\).  There are four ways to select
one carrier for each cut.  Every selection gives the same relation graph:

\[
 2_1\longrightarrow2_2,\qquad
 2_1\leftrightarrow3\leftrightarrow2_2.             \tag{3.1}
\]

Its two-way graph is connected, so it is one collision block.  The tame
normal form has one dense Dickson family; its power specialization lies at
the monomial boundary \(z=0\).  Consequently the reduced all-four
degree-twelve intersection is the Dickson component, agreeing with the
existing exact degree-twelve atlas.

## 4. Degree thirty

The prime multiset is \((2,3,5)\).  For the six requested cuts

\[
 2,3,5,6,10,15,
\]

each cut has two carrier words.  Hence there are \(2^6=64\) possible choices
of one carrier per cut.  Direct finite enumeration shows that every one of
the 64 choices gives the same relation graph:

\[
 2\leftrightarrow3,\qquad
 2\leftrightarrow5,\qquad
 3\leftrightarrow5.                                  \tag{4.1}
\]

Thus the normalized relation graph is the fully bidirected triangle,
independently of which complete decomposition refines each initial cut.
Tame multi-collision theory gives two normal-form labels on this strongly
connected block:

- the Dickson family
  \[
  D_{30}(W+t,z)-D_{30}(t,z);
  \]
- its power degeneration at \(z=0\).

The latter is contained in the closure of the former, so there is one reduced
irreducible component.

The all-six Hessian intersection is already known to be globally
scheme-theoretically synchronized by the spanning-tree lift certificates.
Synchronization identifies it with the graph of one linear coefficient over
the ordinary composition intersection, without changing its nilpotent
structure.  Combining that theorem with (4.1) proves:

> **Degree-thirty reduced component-completeness theorem.**
> The reduced all-six degree-thirty Hessian--Ritt intersection is the single
> Dickson component.  Every point has a prime-refinement carrier system, so
> there is no additional synchronized reduced component missing all six
> complete-decomposition charts.

This removes the ambient minimal-prime calculation from the degree-thirty
reduced problem.  The unreduced schemes attached to different half-braids
remain distinct and are not collapsed by this theorem.

## 5. What the reduction does and does not optimize

The LL/block-system front end now performs three tasks before coefficient
algebra:

1. it recovers the expected codimension from branch packets;
2. it enumerates the finite carrier systems and their relation graphs; and
3. after synchronization, it proves reduced component coverage.

Large Gröbner calculations remain necessary only for:

- synchronization when it is not already structural;
- nilpotent thickenings and path-to-boundary comparison;
- affine-sheet boundary intersections; and
- components whose relation graph contains a genuinely new power/Dickson
  overlap.

In particular, no ambient degree-thirty minimal-prime computation is needed.
The next all-degree problem is not reduced component completeness for the
degree-thirty braid.  It is lifting terminal refinements and Ritt
two-cell coherence through nonreduced primary thickenings.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_ll_ritt_reduction.py
```

The checker verifies (1.1)--(1.5), enumerates all degree-twelve and
degree-thirty carrier choices, constructs their relation graphs, checks the
single strongly connected block, and verifies directly that \(D_{12}\) and
\(D_{30}\) have two finite critical values while the power family is their
\(z=0\) specialization.
