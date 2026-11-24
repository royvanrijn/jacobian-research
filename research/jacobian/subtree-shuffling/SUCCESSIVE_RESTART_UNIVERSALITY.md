# Arbitrarily many cubic restart rounds recover the full fixed-index problem

The complete one-layer inverse ceiling does not turn an arbitrary number of
subsequent admissible layers into a modest extension problem. A fixed coordinate
already permits every finite graded inverse witness to grow as a separate branch.
This note makes that obstruction to a shortcut precise.

Work over a characteristic-zero field. A **scaffold** is a fixed polynomial
automorphism \(F_*=I-H_*\), where \(H_*\) is homogeneous cubic,
\((JH_*)^p=0\), and one coordinate \(u\) has \((H_*)_u=0\). Its inverse
degree is \(D_*\). A restart round adjoins finitely many independent coordinates
whose homogeneous cubic components use previously available coordinates and
whose addition preserves \(p\)-th Jacobian power zero. Components in one round
are parallel. No later component is required to involve a particular existing
high-degree output.

The argument applies in particular to the retained index-four scaffold,
including the fixed nine pre-\(S\) components and the complete first helper
layer. It leaves those components unchanged. The appended branch uses only
\(u\) and its own new coordinates, so a prohibition on using other coordinates
of the original base does not exclude this branch.

## Gluing along a fixed coordinate

**Lemma FixedCoordinateAmalgamation.** Suppose two homogeneous cubic maps on
\((u,v)\) and \((u,z)\), respectively, have zero \(u\)-component and satisfy
\(J^p=0\). Their union on \((u,v,z)\), sharing the coordinate \(u\), also
satisfies \(J^p=0\). A nonzero \((p-1)\)-st power in either original map
remains nonzero. If both identity-minus-map maps are automorphisms, the union
is an automorphism and its inverse degree is the maximum of their two inverse
degrees.

Write the combined Jacobian in the form
\[
 M=\begin{pmatrix}0&0&0\\a&B&0\\c&0&D\end{pmatrix}.
\]
For every \(k\ge1\), direct induction gives
\[
 M^k=\begin{pmatrix}
 0&0&0\\ B^{k-1}a&B^k&0\\D^{k-1}c&0&D^k
 \end{pmatrix}.
\]
The asserted nilpotence and its exactness follow by inspecting the two old
principal blocks. The map itself has the form
\[
 (u,v,z)\longmapsto (u,F_v(u,v),F_z(u,z)),
\]
and its inverse is \((u,G_v(u,v),G_z(u,z))\), proving the degree assertion.

The fixed-coordinate hypothesis is essential. Joining the cubic chain
\(a\mapsto u\), with component \(H_u=a^3\), to the chain
\(u\mapsto z_1\mapsto z_2\mapsto z_3\), with components
\(u^3,z_1^3,z_2^3\), produces a chain of length four and Jacobian index five.
The latter chain had index four when \(u\) was fixed.

## Every universal inverse witness occurs in a restart tower

Let \(A_{3,p}=\bigoplus_{j\ge0}A_j\) be the free symmetric ternary algebra
on one generator \(x\), modulo the polarized identity \(L_a^p(b)=0\), with
\(L_a(b)=\mu(b,a,a)\). Grade by internal vertices, so
\[
 \mu(A_i,A_j,A_k)\subseteq A_{i+j+k+1}.
\]
Each graded piece is finite-dimensional, and \(A_0=Kx\). Define
\[
 U_0=x,\qquad
 U_N=\sum_{i+j+k=N-1}\mu(U_i,U_j,U_k)\quad(N\ge1).
\]
These are the weighted inverse coefficients, with ordered triples in the
sum. The precise dictionary and graded realization are proved in
[ShuffleInverseDegreeEquivalence](SUBTREE_SHUFFLE_AND_INVERSE_DEGREE.md).

**Theorem RestartRealizesEveryGrade.** If \(U_N\ne0\), the scaffold admits
a tower of at most \(N\) restart rounds whose final inverse degree is exactly
\[
 \max\{D_*,2N+1\}.
\]
The old scaffold is retained unchanged. If its Jacobian already has exact
index \(p\), so does every enlarged map.

Take the finite-dimensional quotient
\[
 B_N=A_{3,p}/\bigoplus_{j>N}A_j.
\]
The discarded subspace is an ideal, since grafting adds one to the sum of
nonnegative grades. The quotient retains the Engel identity. In a basis
ordered by grade, put \(T(v)=\mu(v,v,v)\). The coordinate of its generator
has zero component, and every output of grade \(j\ge1\) uses only input
coordinates of strictly smaller grade. Therefore \(I-T\) is triangular.
Polarization gives \(JT(a)=3L_a\), hence \((JT)^p=0\).

The full set of coordinates of grades at most \(j\) is itself the quotient
\(B_j\), so its map still satisfies \(J^p=0\). Adjoin the basis coordinates
of grade \(j\) in the \(j\)-th parallel round, identifying the generator
coordinate with the fixed scaffold coordinate \(u\). The preceding lemma
proves admissibility of every round and preservation of an old exact index.

Give a coordinate of grade \(j\) weight \(2j+1\). Every monomial of \(T\)
has the weight of its output. Triangular inversion preserves these positive
weights, so the inverse degree is at most \(2N+1\). Its degree-\(2N+1\)
coefficient at the generator is \(U_N\ne0\), proving equality. The fixed
coordinate gluing then gives the displayed maximum.

This proof constructs an admissible family from an actual nonzero universal
coefficient. It does not assert that such coefficients exist at any new
index-four grades.

## Equivalence of the all-round and unrestricted uniform bounds

**Corollary AllRoundUniformityEquivalence.** For every integer \(D\ge D_*\),
the following statements are equivalent:

1. Every finite admissible restart tower over the fixed scaffold has inverse
   degree at most \(D\).
2. Every weighted class \(U_N\) with \(2N+1>D\) is zero in \(A_{3,p}\).
3. In every finite dimension over every characteristic-zero extension, every
   homogeneous cubic map \(H\) with \((JH)^p=0\) has polynomial inverse
   \((I-H)^{-1}\) of degree at most \(D\).

The theorem proves \(1\Rightarrow2\). Evaluation of the universal inverse
coefficients in the polarization algebra of any \(H\) proves
\(2\Rightarrow3\): the unique formal inverse terminates, so both inverse
identities are polynomial identities. Statement 3 applies to every enlarged
scaffold, proving \(3\Rightarrow1\).

It is enough to test the towers over the rationals when the scaffold is
rational. The free algebra and its relations are defined over \(\mathbb Q\),
and extension of a characteristic-zero coefficient field is faithful. In
particular nonzero rational universal coefficients cannot disappear under
such scalar extension.

For the complete index-four scaffold, \(D_*=207\). Consequently an
all-round ceiling of 207 is equivalent to
\[
 U_N=0\qquad\text{for every }N\ge104
\]
for **all** cubic index-four maps. The existing complete parallel-restart
theorem concerns a fixed finite layer and does not prove these identities.
A dimension-uniform bound is also distinct from the assertion that each
individual map has some polynomial inverse with a bound depending on that
map or its dimension.

## Complete-basis saturation does not avoid the embedding

The same implication holds if each round is required to adjoin a basis of
*all* admissible cubics on the currently allowed coordinates. At each grade,
the finitely many desired free-algebra components are admissible cubics
using only already selected branch coordinates. They are linearly
independent: by polarization, a linear relation among their diagonal cubic
components annihilates all products in that positive grade, hence the entire
grade, because the free algebra is generated by \(x\). Extend those components to a basis of the admissible cubic space.
If a different basis was chosen, an invertible constant linear change among
that round's new coordinates realizes the same selection.

Unselected helper coordinates may depend on the selected branch, but the
selected branch never depends on them. Projection onto the branch therefore
intertwines the full map and its inverse. The branch's inverse degree is a
lower bound for the full inverse degree. Thus allowing every cubic at every
round does not restrict away the universal witnesses.

## Consequence for further research

A stable algebra or module proving a uniform inverse-support bound under
*every* subsequent admissible cubic extension would settle the fixed-index
weighted inverse problem itself. Additional restrictions that rule out the
independent branch must be explicit and proved stable; removing duplicate
outputs alone does not do so. Conversely, a repeatable degree-increasing
block would provide an unbounded inverse family at the fixed index and
would refute that uniform bound.

The argument uses the existing universal inverse realization, rather than
new finite computations. Its classical context is Umirbaev, *Polarization
algebras and their relations*, MPIM 2014-60, Lemmas 7--8 and Corollary 2,
printed pages 7--8, together with the fixed-arity/fixed-length quantifiers in
Bisi et al., *Random planar trees and the Jacobian conjecture*,
arXiv:2301.08221v3, Conjecture 1.3. Those sources do not resolve the
index-four weighted inverse identities. No new bound, nonvanishing grade,
external review, or formal verification is claimed by this equivalence.
