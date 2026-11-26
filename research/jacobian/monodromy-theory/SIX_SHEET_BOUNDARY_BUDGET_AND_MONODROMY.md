# The boundary budget for five- and six-sheet plane Keller maps

Let \(F:\mathbb A^2_{\mathbb C}\to\mathbb A^2_{\mathbb C}\) be polynomial
with nonzero constant Jacobian, and let \(N=[\mathbb C(x,y):\mathbb C(F)]\)
be its generic sheet degree. Write \(W\) for the entire affine
normalization of the target in \(\mathbb C(x,y)\). Monodromy below
means the action on the \(N\) sheets of this finite map off its reduced
branch locus. The Zariski main theorem identifies the polynomial
source with an open subset of \(W\).

**Theorem (five- and six-sheet branch profiles).** If \(N=5\), every
irreducible branch curve has a transposition meridian and the
monodromy is \(S_5\). There are at most two ramified target curves.
If \(N=6\), the only possible profiles and monodromy groups are

| Distinct ramified target curves | Monodromy | Residual budget |
|---|---|---:|
| One transposition curve | \(S_6\) | 3 |
| Two transposition curves | \(S_6\) | 1 |
| One transposition curve and one 3-cycle curve | \(S_6\) | 0 |
| One 3-cycle curve | \(A_6\) | 2 |

In particular the natural five-sheet action of \(A_5\) is excluded.
Each ramified curve in the table is dominated by one boundary prime
with residue degree one. The table gives necessary conditions for
the original polynomial map.

The six-sheet census agrees with Sections 2–6 of
[the refined six-sheet budget](https://github.com/alok/jacobian-two/blob/ded8e67fb47c155f83ec9bd68af6014499fc2d61/docs/refined-six-sheet-budget.md).
Here the size of the support of each meridian is bounded before
classifying its generating group. This removes the need for a separate
imprimitive group exclusion. The pinned upstream Markdown has
SHA-256 \(9868eb4ce42287f990caf7a9034a55eb4acd1d6a834a3086da7fb4d1a3a1857d\).

## The weighted boundary identity

Use the resolution in Orevkov's Lemma 2.1. A connected component of
the part of the boundary with affine image is a chain with one
dicritical end \(E\). The dicritical curve minus its point over infinity
is an affine line \(E^\circ\). After passing to the Stein factor, a
constant chain meets \(E^\circ\) at just one point. Thus the map from
each \(E^\circ\) to its image in the contracted space is injective.
On the affine target the Stein factor is \(W\). Indeed its finite
normal algebra is the integral closure in the function field of the
resolution, and its fibres under the birational first map are exactly
the connected constant chains and the singleton fibres. This identifies
the affine local multiplicities \(\mu_{\pi(c)}\) of Orevkov with the
degrees of the local finite branches of \(W\) at \(\pi(c)\).

Denote the generic normal index of \(E\) by \(e_E\), and write
\(r_E:E^\circ\to B_E^\nu\) for its map to the normalization of its
affine target image. This is finite because any nonconstant target
coordinate is a polynomial on \(E^\circ=\mathbb A^1\) over which its
coordinate ring is finite. Passing to the integral closure in the
image field preserves this finiteness. After projective completion the
target is rational and its points at infinity have preimages among the
unique point at infinity of \(E^\circ\). Thus \(B_E^\nu\cong\mathbb A^1\).
Write \(d_E=\deg r_E\) and \(a_c\) for its local degree at \(c\).

For \(c\in E^\circ\), choose a nearby regular point in \(B_E^\nu\) mapping
to the smooth locus of the branch, with \(a_c\) preimages converging to
\(c\). Each of these preimages contributes \(e_E\) sheets after a small
normal displacement by Orevkov's generic local form in Lemma 3.1.
The normal displacement can be chosen within the neighbourhood defining
the finite local degree at \(\pi(c)\). It follows that
\(\mu_{\pi(c)}\ge e_Ea_c\), including at the attachment point of a
contracted chain. Distinct local branches are counted by disjoint
neighbourhoods before taking the normal displacement.

**Proposition (refined Orevkov identity).** In this resolution,
\[
 \sum_E(e_Ed_E+\delta_E)=N-1,\qquad
 \delta_E=\sum_{c\in E^\circ}(\mu_{\pi(c)}-e_Ea_c)\ge0.
\]
Indeed Orevkov's Lemma 4.2 sums \(e_E+\sum_{x\in\pi(E^\circ)}
(\mu_x-e_E)\) over \(E\). The injectivity just proved permits summation
over \(c\) instead of \(x\). Finite Riemann–Hurwitz for the polynomial
\(r_E\) gives \(\sum_c(a_c-1)=d_E-1\) because the missing point has
local degree \(d_E\). The formula follows term by term. In particular
an unramified dicritical counts in this identity with \(e_E=1\).

See [Orevkov's author retypesetting](sources/orevkov-three-sheet-author-1987.pdf),
Lemmas 2.1(a)–(c), 3.1, 4.2 and Corollary 4.3, PDF pp. 2, 4, 7.
The argument uses these statements in arbitrary degree, before the
separate specialization to three sheets in Section 6.
[Author copy and bibliographic record](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf).

## The sole branch curve and its normalization

For a reduced branch component \(B_i\) write \(c_i\) for the number of
sheets moved by a generic positive meridian. Above a generic smooth
point the boundary primes dominating \(B_i\) contribute \(d_E\)
disjoint cycles of length \(e_E\). Lemma 3.1 gives this generic local form on each branch. Thus
\[
 c_i=\sum_{E\to B_i,\ e_E>1} e_Ed_E,\qquad \sum_i c_i\le N-1.
\]
At a point \(y\in B_i\), each point of \(B_i^\nu\) over \(y\) contributes
at least \(c_i\) to the local degree sum. For each \(E\), sum
\(e_Ea_c\) over the finite preimage of that normalization point to get
\(e_Ed_E\), and use the disjointness of the neighbourhood counts.
Distinct preimages on a given \(E^\circ\) remain distinct under the
contractions. Consequently if \(2c_i>N\), the map \(B_i^\nu\to B_i\)
is injective.

**Proposition (sole branch support bound).** If the reduced branch
locus is one irreducible curve \(B\), then \(2c_B\le N\).

Suppose \(B^\nu\to B\) were injective. Since \(B^\nu\cong\mathbb A^1\)
and the map is finite, \(B\) is topologically contractible. By
Lin–Zaidenberg, in the formulation of Gurjar and Miyanishi, Theorem 2,
a polynomial target change gives \(B=\{X^m=Y^n\}\) with coprime positive
\(m,n\). Weighted dilations \((X,Y)\mapsto(r^nX,r^mY)\) give the
complement the product description by a sphere complement and a radial
interval. The norm on each positive dilation orbit is strictly
increasing from zero to infinity, so inclusion of the complement in
a ball about the origin induces an isomorphism on fundamental groups.

Over the origin, choose a point on a ramified dicritical. Orevkov's
unrefined bracket for that component contains that point's local
degree, with all other terms nonnegative. Hence its degree is at
most \(N-1\), and a second point of the finite fibre exists. By properness,
the inverse image of a sufficiently small ball separates into disjoint
neighbourhoods of the fibre points, each contributing sheets off the
branch. Its local monodromy is intransitive. The finite étale cover
over the global complement is connected by integrality of \(W\).
Purity for the finite map from the normal surface to the smooth plane
allows the entire nonétale image to be read from its ramified divisors.
The local-to-global fundamental-group isomorphism gives a contradiction.

See [Gurjar–Miyanishi, Theorem 2](sources/gurjar-miyanishi-contractible-curves-1996.pdf)
and [the normal-to-regular branch-purity criterion](https://stacks.math.columbia.edu/tag/0BMB).
In particular the argument concerns the reduced **ramified** target curve
of the finite normalization. It also applies when unramified omitted
primes have images elsewhere.

## Monodromy and the six-sheet cases with an unramified prime

A nontrivial meridian moves at least two sheets. For \(N=5\), any
meridian moving at least three sheets would be the only branch within
budget four, contradicting the sole branch support bound. For \(N=6\),
the same argument excludes all meridians moving at least four sheets.
The branch locus is nonempty for \(N>1\): otherwise purity would give a
connected étale cover of the simply connected affine plane of degree
\(N\). The possible positive supports in degree six are therefore
\(2\), \(2+2\), \(2+3\), and \(3\). A support of two is a single
transposition and a support of three is a single 3-cycle.

The meridians normally generate monodromy, by spanning a loop with a
disk in the plane and perturbing it to meet smooth branch points
transversely. If all generators are transpositions, the graph of their
supports, including conjugates, is connected by transitivity. Its
edges generate the full symmetric group. If a transposition and a
3-cycle occur in degree six, any nontrivial block size would be two or
three. A nontrivial permutation of blocks moves at least four sheets,
so every generating conjugate fixes every block setwise. This
contradicts transitivity. In the resulting primitive group the graph
of conjugate transpositions has components forming a block system and
has an edge. It is connected, so the group is \(S_6\).

If only 3-cycles occur, their supports form a connected hypergraph.
Beginning with one triple on a set \(S\), the generators already
contain \(A_S\). A triple with two old points and one new point grows
this to the alternating group on the enlarged set: conjugation by
\(A_S\), which is transitive on unoriented pairs, and inversion give
all triples made of the new point and pairs from \(S\). If a triple
has one old point \(v\) and two new ones \(x,y\), conjugate
\((v\,x\,y)\) to \((w\,x\,y)\) with \(w\ne v\) old. The product
\((v\,x\,y)(w\,x\,y)^{-1}\) has two old points and one new point,
and then the original triple adds the other new point. Connectedness
finishes the proof. All generators are even, so in this last case the
group is \(A_6\).

**Corollary (absence of a proper intermediate field).** For \(N=5,6\)
the monodromy action is primitive. The Galois closure correspondence
identifies an intermediate field with a subgroup containing a sheet
stabilizer, and such a subgroup defines a block containing that sheet.
Thus \(\mathbb C(F)\subset\mathbb C(x,y)\) has no proper intermediate
field. Every hypothetical map in these degrees belongs to the
[branch without an intermediate](../affine-complements/NORMAL_FORM_COVERAGE_AND_COMMON_LEADING_DIVISORS.md#1-three-global-coverage-obligations).

**Corollary (six sheets with an unramified boundary).** If an actual
omitted unramified boundary prime is present, a transposition curve
and a distinct 3-cycle curve are impossible. Two transposition curves
then force exactly three boundary primes, with
\((e_E,d_E)=(2,1),(2,1),(1,1)\) and \(\delta_E=0\) for each.
The boundary prime at the exceptional polar fold of type
\((p,r,n)=(2,1,5)\) uses two units of the budget because its own
normal index is two and its residue degree is one. The fold exponent
\(n\) is a local parameter of that model.

For the first statement, the ramified curves and any unramified
prime spend at least \(2+3+1=6\). In the two-transposition case
\(2+2+1=5\) forces equality in each summand. The
[unramified residue-degree theorem](../../hessian/hc4-triangular-graphs/UNRAMIFIED_BOUNDARY_RESIDUE_DEGREE_AND_TARGET_NORMALIZATION.md)
identifies its tangential degree as one even away from saturation.
The [exceptional fold chart](../../hessian/hc4-triangular-graphs/EXCEPTIONAL_POLAR_FOLD_WHOLE_CHARTS_AND_BIRATIONAL_IMAGES.md)
gives the local data used here. The two saturated transposition
curves still need an actual polynomial-admission or marked-trace
argument. Empty unramified boundary leaves the distinct transposition
and 3-cycle profile available to this count.

In the \(A_6\) case the sole branch normalization has a collision by
the local-to-global argument. Any such value has exactly two
normalization points above it, each contributing three sheets.
This exhausts the finite fibre degree. Every unramified omitted
prime has image disjoint from those collision values. The remaining
\(A_6\) admission problem includes the geometry away from them.

For either ramification index, an actual polar source velocity has
a [detectable marked trace at the original prime](../../hessian/hc4-triangular-graphs/EXCEPTIONAL_FOLD_POLYNOMIAL_ADMISSION_AND_MARKED_TRACE.md).
For pole order \(q>0\) along a prime of normal index \(e\), a marker
in the full finite order has trace pole order \(\lfloor q/e\rfloor+1\)
under a transverse target derivation. This locates the potentially
obstructing trace in both the index-two and the index-three cases.
Automatic ordinary codifferent membership for a first velocity
in a source-finite frame remains a separate condition.

The [finite cusp classification](KELLER_BOUNDARY_CUSP_PROFILES_AND_OREVKOV_INVARIANT.md)
determines the analytic singular branches from the six-sheet budget.
It also connects the saturated cases to the Orevkov invariant of
the reduced ramified curve.

## Finite check and source boundary

The [GAP census](verify_six_sheet_boundary_budget.py) enumerates
conjugacy classes of subgroups of \(S_5,S_6\), selects transitive
ones and enumerates all multisets of nontrivial meridian types
permitted by the budget and sole-branch bound. It checks normal
generation inside the candidate subgroup, including identical
cycle types in different conjugacy classes. Run it with
\(\texttt{sage -python}\) using SageMath 10.9 and GAP 4.15.1.
The geometric arguments above give the result directly. The census
controls their finite group consequences.

The [publisher record for Żołądek](https://www.sciencedirect.com/science/article/pii/S0040938308000207)
was checked for the bibliographic identity. The accessible abstract
and author bibliography did not supply the text of Theorem 6.12,
so the reported all-group degree-five exclusion still has a
primary-source verification task. The natural \(A_5\) mathematical
exclusion follows from the proof here using the separate sources
identified above.
