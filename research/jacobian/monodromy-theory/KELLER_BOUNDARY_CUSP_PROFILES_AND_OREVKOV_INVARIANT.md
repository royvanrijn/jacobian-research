# Local cusp types and the Orevkov invariant of Keller branch curves

Let \(F:X=\mathbb A^2_{\mathbb C}\to Y=\mathbb A^2_{\mathbb C}\) have
nonzero constant Jacobian. Let \(\pi:W\to Y\) be the full finite
normalization in the source function field and let \(B\subset Y\) be
its reduced branch curve. Write \(N=\deg\pi\). In the relatively minimal
finite-boundary resolution of [Orevkov, Lemma 2.1 and §5](sources/orevkov-three-sheet-author-1987.pdf),
the omitted primes in \(W\) are disjoint, each with affine-line
normalization. Each is the image of a smooth dicritical attached at
the endpoint of at most one constant rational chain. Contraction is
injective on the affine dicritical, so its image is analytically
unibranch. A point obtained from a nonempty constant chain is a
singular point of \(W\). See the
[boundary identity](SIX_SHEET_BOUNDARY_BUDGET_AND_MONODROMY.md#the-weighted-boundary-identity)
for the identification with the full affine normalization.

Write \(e\) for the generic normal ramification index of an omitted
prime \(E\), and \(\mu\) for the local degree of \(\pi\) at a point of
\(E\). These refer to the full finite map, including when the target
has other branches through the same value.

## Local degree and the finite branch quotient

**Proposition KellerSmoothBranchLocalModel.** For \(e>1\) the image
branch at the point of \(E\) is smooth if and only if \(\mu=e\).
In this case the local finite map is analytically equivalent to
\((u,v)\mapsto(u^e,v)\). Every omitted prime with \(e=1\) is a smooth
affine line in the étale locus of \(W\); its map to the target is
immersive and has residue degree one.

For the forward implication, separate the preimage of a small target
ball into the finite germs at its distinct preimages. The branch locus
in the target of the selected germ comes only from divisors ramified
through its source point, by purity. For \(e>1\) there is just one such
germ of divisor, by disjointness and analytic unibranchness of \(E\).
Its image is the selected branch, regardless of the other germs of
\(B\) through the target point. If that branch is smooth, the
complement of a small coordinate bidisk cut by this germ has
fundamental group \(\mathbb Z\). The cover on this complement is
connected because the source is an irreducible normal germ. It is
the unique connected cyclic cover of degree \(\mu\); uniqueness of
normal extension identifies the whole finite germ with the power
map. Its generic index is then \(\mu\). The reverse implication is
Orevkov's Lemma 5.2, proved in §5 before the specialization to three
sheets.

At a point of an unramified omitted prime, purity gives an unbranched
germ, hence local degree one. Apply the same Lemma 5.2 with \(e=1\):
the parametrization of the image branch is locally smooth. Since
\(\pi\) is a local isomorphism, \(E\) is smooth in \(W\). The induced
finite map between the two affine-line normalizations is unramified
at every point, hence is an isomorphism. This also recovers the
[residue degree calculation](../../hessian/hc4-triangular-graphs/UNRAMIFIED_BOUNDARY_RESIDUE_DEGREE_AND_TARGET_NORMALIZATION.md).
The target curve can have self-identifications.

**Theorem KellerOrevkovInvariantObstruction.** If every analytic
branch of \(B\) at every affine point is smooth, then the monodromy
of \(\pi\) factors through the Orevkov invariant \(\mathcal O(B)\).
Consequently for \(N>1\) at least one of the following is necessary:
\(\mathcal O(B)\) is nonabelian, or \(B\) has a singular analytic branch.

The invariant is the quotient of \(\pi_1(Y\setminus B)\) by the normal
closure of the images of the commutators of the local complement
groups at singular points of \(B\). At a given target value the finite
germs act on disjoint sets of sheets. Each local action is cyclic:
it either comes from a power map at a ramified prime, or comes from a
degree-one germ. This kills all the local commutator subgroups. The
resulting quotient action is still transitive. If it is abelian, it
is regular, since an element fixing a sheet then fixes every sheet.
But a generic branch meridian is nonidentity and fixes a sheet: its
moved support is bounded by \(N-1\) in the
[refined budget](SIX_SHEET_BOUNDARY_BUDGET_AND_MONODROMY.md).
For \(N>1\) there is at least one branch by purity and simple
connectivity of \(Y\).

The abelianity input can be checked from the link at infinity in
specific cases. In the conventions of
[Neumann–Norbury, Theorem 2](https://researchers.ms.unimelb.edu.au/~norbury%40unimelb/papers/orevkov.pdf)
([retained PDF](sources/neumann-norbury-orevkov-invariant.pdf)),
the splice condition of Theorem 5.2 on each component, together with
\(d_n>1\) as used in the proof of Theorem 2, suffices for a
polynomially parametrized curve. Here their \(d_n\) is the terminal
splice invariant \(-K\cdot\overline{B_i}-1\), different from residue
degree. Theorem 5.2 uses \(d_n\ge p_1q_2\cdots q_n-1\), with weights
taken on the path of their Figure 2. The curve tested in this
application is the reduced ramified curve \(B\). That theorem's
infinity condition must be checked on the union, with the mutual
link data of its components. For a single polynomially parametrized
component with one Puiseux pair at infinity, their Proposition 5.5
and the argument of §5.3 supply the positive braid case.

## Small ramification indices at a contracted endpoint

**Proposition KellerIndexTwoThreeEndpointTypes.** An omitted prime
of index two is everywhere on the smooth locus of \(W\). At a
singular endpoint on an index-three prime, the germ of \(W\) is the
quotient by \(-I\), and its local degree \(\mu\) is even. If \(\mu=4\),
the image branch has Puiseux pair \((2,3)\).

This uses the Hirzebruch–Jung characterization of the singularity
of a contracted rational linear chain as a cyclic quotient.
[Popescu-Pampu, Propositions 6.2–6.3 and §6.2](sources/popescu-pampu-math-0506432v2.pdf)
([arXiv v2](https://arxiv.org/abs/math/0506432v2))
gives the two directions of the resolution statement. If the
constant chain has a \((-1)\)-curve, contract it on the smooth
resolution first. The curve meets at most two other curves of
the chain with its dicritical attached, at distinct transverse
points. The contraction preserves the chain and the transverse
dicritical at its end. Repeating leaves the minimal chain with
all self-intersections at most \(-2\). The order and
parity argument follows the pinned
[one-dicritical analysis](https://github.com/alok/jacobian-two/blob/ded8e67fb47c155f83ec9bd68af6014499fc2d61/docs/one-dicritical-source-smoothness.md#1-the-cyclic-endpoint-obstruction).
The hypotheses of its local argument apply at any endpoint of the
affine boundary, because no other omitted prime passes through
that point.

For clarity, the lift of the dicritical to the cyclic cover is a
single smooth branch. In a quotient with weights \(1,a\) modulo
\(n\), write \(b a\equiv1\bmod n\), \(0<b<n\).
The chart adjacent to the first coordinate ray in the minimal
toric resolution has coordinates \(s=x/y^b\) and \(t=y^n\). Its rays
are \((1,0)\) and \((b/n,1/n)\) in the quotient lattice. The affine
part of the end exceptional curve in this chart is \(t=0\).
A transverse curve meeting it at a finite value of \(s\) has equation
\(s=\phi(t)\). Its lift is \(x=y^b\phi(y^n)\), a single smooth curve,
which can be straightened by the equivariant subtraction of this
right side. At the opposite end interchange the axes. In a
one-vertex chain the two endpoint charts cover the entire exceptional
curve, so the same argument applies at every curvette.

Let \(h=(u,v)\) be the composite of the cyclic cover and \(\pi\) at
such a point. The cover is étale away from the origin and there is
only one ramification divisor of the composite. By the usual
coordinate Jacobian order at its generic point,
\(\det dh=x^{e-1}\epsilon\) with unit \(\epsilon\). Both coordinates
are invariant, and have zero linear term because a small cyclic
quotient has no trivial tangent summand. Thus \(e=2\) is impossible
at a singular endpoint.

For \(e=3\), the initial forms of \(u,v\) are independent quadratics
whose Jacobian is a nonzero multiple of \(x^2\). If the two forms
had no common root, they would give a degree-two morphism of the
projective line whose two ramification contributions are at
distinct points. Its quadratic Jacobian could not be a doubled
line. If the forms share a root, the common line is \(x=0\):
multiplying two independent linear forms by a common line gives
a Jacobian proportional to that line squared. After linear target
changes and coordinates preserving \(x=0\), the initial forms are
\(x^2,xy\). The invariant products of \(x\) with two independent
linear forms show that the cyclic action is scalar and has order
two. In fact if \(x\) has character \(\chi\), the two residual
linear factors have character \(\chi^{-1}\), and one is \(x\)
itself. Thus \(u,v\) are even series for simultaneous negation.

The germ \(v=0\) consists of two smooth transverse branches. On
the branch tangent to \(y=0\), take \(x=t\). Differentiating \(u\)
on this branch gives \(J_h/v_y\), of order one, so the intersection
contribution to \(I_0(u,v)\) is two. On the other branch take
\(y=t,\ x=\beta(t)\), and put \(r=\operatorname{ord}\beta\). Here
the derivative is \(-J_h/v_x\), of order \(2r-1\), giving contribution
\(2r\). Finiteness rules out \(\beta=0\) identically: then \(J_h=0\)
and \(v=0\) along \(x=0\), forcing \(u\)
constant there. Negation preserves each branch, so \(\beta\) is
odd. Since \(2\mu=I_0(u,v)=2+2r\), \(\mu=r+1\) is even.

When \(\mu=4\), \(r=3\). Write
\(v(0,y)=c y^4+O(y^6)\), with \(c\ne0\) by the equation for
\(\beta\), and write \(u_x(0,y)=A y^3+O(y^5)\).
The coefficient at \(y^3\) in \(\partial_x J_h(0,y)=0\) is
\(8c-2A\). This follows from \(u_{xx}(0,0)=2,\ v_{xy}(0,0)=1\),
\(v_x(0,y)=y+O(y^3)\), and \(u_y(0,y)=O(y^3)\). Thus \(A=4c\).
Now \(J_h(0,y)=0\) gives
\(u_y(0,y)=16c^2 y^5+O(y^7)\), so
\(u(0,y)=(8c^2/3)y^6+O(y^8)\).
The quotient parameter along the omitted prime is \(y^2\).
Its two image coordinates have orders three and two in that
parameter, giving the claimed ordinary cusp.

The singular local case occurs among analytic finite maps. An
explicit model on the double cover is
\[
(u,v)=\bigl(x^2+4xy^3+\tfrac83y^6,\;xy+y^4\bigr),
\qquad J_h=2x^2.
\]
Both functions are invariant under simultaneous negation. The
branches of \(v=0\) are \(y=0\) and \(x=-y^3\), where the orders
of \(u\) are two and six. Thus the local degree after taking the
quotient is four, and the restriction to \(x=0\) has orders
three and two in \(y^2\). This local example leaves the global
polynomial source question open.

## Consequences in generic sheet degree six

**Theorem SixSheetFiniteCuspProfiles.** A six-sheet plane Keller
map with group \(S_6\) has smooth full finite normalization. An
index-two dicritical is immersed except at points where the
local degree \(\mu>2\); the corresponding singular image branch
has Puiseux pair \((\mu-1,\mu)\).

For \(A_6\) let \(k\) be the number of unramified omitted primes.
The possibilities are as follows.

| \(k\) | Full finite normalization | Singular analytic branches of \(B\) |
|---:|---|---|
| 0 | Smooth | One cusp \((2,5)\) from local degree five |
| 1 | One singularity of type \(A_1=\mathbb C^2/\{\pm I\}\) | One cusp \((2,3)\) from local degree four |
| 2 | Smooth | None |

The two saturated \(S_6\) profiles (mixed indices with no unramified
prime, or two transposition curves with one unramified prime) and
the \(A_6\) profile with \(k=2\) require a nonabelian Orevkov invariant
of \(B\). This necessity also holds in any unsaturated configuration
whose branch normalizations are immersive.

To check the local smooth-source cases, a smooth point of \(W\) on
a dicritical has smooth \(E\) in Orevkov's minimal contraction.
For \(e=2\), \(J_\pi\) has order one along \(E\), so the map has rank
one. For \(e=3\) and local degree four or five, rank zero is
impossible. In degree four the quadratic leading forms would
be coprime and their determinant a square of a linear form, giving
the impossible projective-line map above. In degree five the
quadratics must be independent, because cancellation to order
three implies local degree at least six by the local intersection
bound. Some linear combination of the quadratics has distinct
linear factors: the conic of doubled lines contains no projective
line. Use that combination for the first target coordinate.
Its two local smooth transversal zero branches have contributions
at least two each to local degree five, hence two and three.
On the branch giving three the Jacobian has order three:
differentiate the second coordinate along that branch and add
one for the transverse factor of the first coordinate. This
contradicts the even order of a unit times the square of the
local equation of \(E\).

In rank one, use coordinates \(\pi=(a,f(a,v))\) and make \(E\)
the graph \(a=v^m\) by reparametrizing \(v\). Here \(m\) is the
finite order of \(a|_E\); if it vanished identically the
determinant would force \(f|_E\) constant too, contrary to
finiteness. Then \(f_v\) is a unit times \((a-v^m)^{e-1}\) and
\(\mu=m(e-1)+1\). Subtract \(f(a,0)\) in the target. Along \(E\)
the first nonzero term of \(f\) at degree \(m(e-1)+1\) has
nonzero coefficient, a unit at the origin times
\(\int_0^1(1-s^m)^{e-1}\,ds\). Thus the branch orders
are \(m,m(e-1)+1\), coprime. This proves the \(S_6\) cusp assertion
and the smooth degree-five \(A_6\) case, in agreement with the
[one-dicritical local calculation](https://github.com/alok/jacobian-two/blob/ded8e67fb47c155f83ec9bd68af6014499fc2d61/docs/a6-one-dicritical-local.md).

For \(S_6\), all index-two primes and all unramified primes are
smooth by the preceding results. In the mixed profile the
index-three prime has zero excess, hence is smooth by the
power-map model. For \(A_6\), the single index-three prime \(E\)
has residue degree one. Each unramified prime spends one unit,
with zero excess. Thus \(\sum_{x\in E}(\mu_x-3)=2-k\).
Local degree four occurs only at a singular source point, and
there is at most one such point on \(E\) by the endpoint-chain
description. This excludes the \(1+1\) partition when \(k=0\).
Its unique local degree five point is smooth by the evenness
result for singular endpoints. When \(k=1\), the unique excess
is one and the singular endpoint calculation applies. When
\(k=2\), there is no excess. Cusps in the \(A_6\) table cannot
share a target point with another branch of the normalization:
\(4+3>6\) and \(5+3>6\). All remaining singularities of the
ramified curve arise by identifying distinct smooth branches.
The assertions about \(\mathcal O(B)\) now follow from the
finite branch quotient theorem.

The cusp classification is a necessary condition on a
hypothetical map. The remaining global obstruction must
exclude the indicated branch embeddings or their compatible
finite normalizations with the original polynomial source.
