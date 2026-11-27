# The two poles of a width-five alternating genus-two slice

The [coprime-width obstruction](KELLER_WIDTH_FIVE_ALTERNATING_ORDINARY_CUSP_EXCLUSION.md)
excludes this width from the singular alternating profile. The following
arguments give the separate pole and critical-fibre restrictions at width five.

Consider a hypothetical plane Keller map of generic sheet degree
six in the singular \(A_6\) profile with one unramified omitted
prime. The affine nonproper curve has normalization \(\mathbb A^1\)
and a unique ordinary cusp. The full affine source normalization
has one \(A_1\) point over that cusp. The ramified omitted prime
has normal index three and residue degree one. Write \(p,q\) for
a polynomial target frame. Suppose \(p\) has width five and the
generic compactified \(p\)-slice has genus two. Its two poles of
\(q\) have orders \(s_1,s_2\) with \(s_1+s_2=6\).

**Theorem SixSheetWidthFiveAlternatingGenusTwoPoleMarks.**
The two pole points on the generic compactified \(p\)-slice are
not exchanged by its hyperelliptic involution. The conclusion holds
over an algebraic closure of the generic base field.

**Proposition SixSheetWidthFiveAlternatingGenusTwoCriticalFibres.**
At an immersed point of the nonproper curve the projection on its
normalization has local degree at most three. If it has local
degree three, its projection value is distinct from those of the
affine cusp and of all other critical points. Its compactified
fibre has two rational normalizations, one with each pole. On the
fibre itself one component is smooth and the other has a node.
The two components meet at that node, forming an ordinary triple
point of the fibre. The pole orders are \(1,5\) or \(3,3\).

At the projection value of the affine cusp, if the strict fibre
has two components, their two branches on the resolution of the
\(A_1\) point lie on different components. For transverse
projection each pole order is odd. For tangent projection the
pole orders are \(2,4\). Consequently the tangent profile with
an immersed local-degree-three point has integral strict cusp
fibre. If an immersed local-degree-two point shares the tangent
cusp projection value, the strict cusp fibre splits, and the
component with the order-four pole contains the immersed cusp.

## Compactification and local degrees

Use the [singular relative compactification and local cusp
model](KELLER_SIX_SHEET_CUSP_SLICE_GENUS.md#singular-normalization-with-one-ordinary-affine-cusp).
Let \(S\) be the smooth relative surface obtained by normalizing
\(\mathbb A^1_p\times\mathbb P^1_q\) in the source field and
resolving the \(A_1\) point. Write \(D\) for the exceptional
\((-2)\)-curve, \(H_1,H_2\) for the disjoint pole sections, and
\(E,U\) for the strict transforms of the ramified and unramified
omitted primes. Denote the projection of the affine cusp by
\(p_0\). The finite-base fibres of \(S\to\mathbb A^1_p\) are
connected, reduced and of arithmetic genus two. Every component
other than \(D\) maps finitely and surjectively to
\(\mathbb P^1_q\) and contains a pole point. There are at most
two such components, and each meets exactly one pole section if
there are two. The pole sections meet the fibres at smooth
points, with the indicated orders of \(q\) on the corresponding
components' normalizations.

Off the \(A_1\) resolution, the immersed ramification germ has
normal index three. If \(m\) is the local degree of the projection
on the normalized image branch, the slice has local equation
\(t^m+w^3=0\). Indeed take \(t\) along the image branch and its
normal coordinate independently, choosing \(t\) from \(q\)
when \(m>1\). After an analytic change of source normal coordinate
the pulled-back equation is \(w^3+p|_E(t)-p(E(0))=0\).
An invertible change of \(t\) gives the displayed equation and
keeps \(q-q(E(0))\) a local coordinate in \(t\). Thus for
\(m=1,2,3,4\) the respective slice germs are smooth, ordinary
cuspidal, ordinary triple, and unibranch of delta invariant three.
On the normalization of these slice germs the local degrees of
\(q\) for \(m=1,2,3\) are \(3,3,1\), the last on each of its
three branches. For \(m=1\) this follows directly from the
index-three normal coordinate and transversality of the
\(p\)-slice with the ramification curve. Elsewhere the
unramified germs introduce no singularity or ramification on
finite slice points.

Near the intersection of the ramified prime with \(D\), the
resolved cusp model is
\[
 u=z(s^2+4sz+\tfrac83 z^2),\qquad v=z(s+z).
\]
Here \(D\) is \(z=0\), \(E\) is \(s=0\), and the linear part of
\(p-p_0\) is \(A u+B v\). When \(B\ne0\), the two strict
branches meet \(D\) at distinct points and the fibre has an
ordinary node at each. On the branch through \(s=0\) one has
\(s=-z+O(z^2)\), so \(u\) has order three with leading term
\(-z^3/3\). A complementary \(q\), modulo adding a multiple
of \(p\), has linear part a nonzero multiple of \(u\).
Consequently its local degree there is three. On the other
branch its linear term in the normal coordinate to \(D\) is
nonzero. In the two charts of the resolution this follows from
the quadratic tangent forms before resolution, \(x^2,xy\):
the second branch has direction \(Ax+By=0\) with \(x\ne0\),
and \(x^2\) is nonzero in this direction. Its local degree is
one.

When \(B=0\), the strict branches have expansions
\(s=\alpha z+O(z^2)\) for the two distinct roots of
\(\alpha^2+4\alpha+\tfrac83=0\). Indeed after dividing \(p-p_0\)
by \(z\), the terms from degree at least two in \(u,v\) have
order at least three in \(s,z\). The two branches and \(D\)
form an ordinary triple point. A complementary \(q\) has
nonzero coefficient of \(v\), and \(v=(\alpha+1)z^2+O(z^3)\)
on each branch, with \(\alpha+1\ne0\). The local degree of \(q\)
is two on both strict branches. On \(E\), the corresponding
local degrees of \(p\) at the affine cusp are two when \(B\ne0\)
and three when \(B=0\).

## Restrictions on the critical fibres

The projection on the normalization of the nonproper curve
is a degree-five polynomial, so its finite critical
multiplicities sum to four. The affine cusp contributes at
least one. An immersed degree-five point is therefore impossible.
An immersed degree-four point uses three and leaves only the
transverse affine cusp. Its value is different from \(p_0\),
since \(4+2>5\). Its slice has an isolated unibranch
singularity of delta invariant three and must be integral by
connectivity. This contradicts its arithmetic genus of two.

For a reduced connected fibre \(Z\), let \(k\) be the number
of its irreducible components, \(g_j\) the genera of their
normalizations and \(\delta_x\) the delta invariants of its
singularities. Then
\[
                  2=\sum_j g_j+\sum_x\delta_x-k+1 .
\]
If an immersed degree-three point shares its value with an
immersed simple critical point away from \(p_0\), the sum of
delta invariants is at least four and \(k\le2\). This is
impossible. It cannot share a tangent cusp value since
\(3+3>5\). At a transverse cusp value the two nodes on \(D\)
and the ordinary triple point give total delta at least five,
also impossible with \(k\le3\). It follows that the
degree-three point is the only directional critical point at
its value, and that its resolved fibre has exactly two
components, both with rational normalization. Connectivity
forces a \(1+2\) distribution of the three branches. The
single-branch component is smooth, the other has an ordinary
node, and their intersection multiplicity is two. Each has
one pole. By the local degrees above, all finite ramification
indices of \(q\) on the normalizations are odd. Applying
Riemann–Hurwitz on each component shows that its one pole
also has odd order, because every finite term in the
ramification sum is even.

At \(p_0\), all singularities off \(D\) are ordinary unibranch
cusps. If two strict components occur, connectivity therefore
puts one \(D\)-attachment branch on each. In the transverse
case the components form a chain with \(D\) in the middle
and the other two components each meeting it transversely
once. If the two strict components have arithmetic genera
\(a,b\), then \(a+b=2\). Both \(q\)-pole orders are odd by
Riemann–Hurwitz on each normalization, since the finite local
degrees include one or three at \(D\), and only odd indices
elsewhere. In the tangent case the two strict components and
\(D\) have pairwise intersection multiplicity one, all at
the ordinary triple point. Their arithmetic genera satisfy
\(a+b=1\). Each strict normalization has exactly one finite
point of even local degree, namely two at the \(D\)-attachment.
Riemann–Hurwitz forces both pole orders even, giving \(2,4\).
This excludes a split tangent cusp fibre in the presence of an
immersed degree-three point.

If an immersed simple critical point shares the tangent
cusp value, the total delta at that value is four (three at
the triple point and one at the immersed cusp). The genus
formula forces two strict components in addition to \(D\),
and both strict normalizations are rational. The immersed
cusp has \(q\)-local degree three on its normalization, so it
lies on the component of degree four. This proves the
proposition.

## The hyperelliptic fibre test

The Picard group of \(S\) is free on \(E,U,D,H_1,H_2\).
The complement of their union is the original affine plane,
whose Picard group vanishes and whose regular units are
constant. Divisor localization gives the generators and
also rules out relations: a rational function with divisor
supported on these curves restricts to a unit of the plane,
so is constant. The pullback of the target two-form gives
the canonical class
\[
                   K_S=2[E]+[D]-(s_1+1)[H_1]-(s_2+1)[H_2].
\]
The coefficient of \(E\) is its ramification defect, and the
pole coefficients follow by differentiating the pole of \(q\)
in each normal parameter with \(p\) a coordinate along the
section. At \(D\) the local equation above gives
\(du\wedge dv=z s^2\,ds\wedge dz\). A change to the target
frame \(p,q\) has nonvanishing Jacobian. This gives coefficient
one at \(D\), and the two-form has no other divisorial zeros
or poles.

The restriction of \(K_S\) to the smooth generic fibre \(C\)
is its canonical class. On a smooth curve of genus two the
canonical linear system is its degree-two hyperelliptic map:
it has dimension one and no base point by Riemann–Roch.
The two pole points could be conjugate only if
\[
                     L=K_S-[H_1]-[H_2]
\]
were generically trivial. Triviality over the algebraic
closure suffices for triviality over \(\mathbb C(p)\):
by extension of scalars the space of sections of the
degree-zero line bundle on \(C\) then has dimension one
over \(\mathbb C(p)\), so a nonzero section trivializes it.
A trivialization over the generic fibre presents \(L\) as
a vertical divisor class on \(S\). Whole fibres are principal
on the affine base. In particular the intersection at any
fibre component is computed from an integral combination of
components of that fibre alone, since components of other
fibres are disjoint from it.

If an immersed degree-three critical point exists, let \(A\)
be the smooth rational component at its projection value.
Its intersection with the other component is two, so
\(A^2=-2\) by intersection with the full reduced fibre.
Adjunction gives \(K_S\cdot A=0\) and hence \(L\cdot A=-1\).
The vertical intersection at \(A\) is always even, a
contradiction.

In the absence of an immersed degree-three point only the
fibre at \(p_0\) can be reducible, since every other fibre
singularity is unibranch. If its strict fibre is integral,
the only possible vertical class modulo whole fibres is an
integral multiple of \([D]\). The coefficient two of \([E]\)
in \(L\) rules this out.

Suppose the tangent cusp fibre splits, and write \(A,B\)
for its strict components, with \(A\) of arithmetic genus
zero (their genera sum to one). All three self-intersections
are \(-2\), by their pairwise intersections of one and the
full reduced fibre. Thus \(L\cdot A=-1\) and \(L\cdot D=0\).
A vertical divisor \(iA+jB+kD\) with zero intersection on
\(D\) satisfies \(i+j=2k\), and so its intersection on \(A\)
is \(3(k-i)\). This cannot be \(-1\).

Finally suppose the transverse cusp fibre splits into the
chain \(A+D+B\) and let \(a=p_a(A)\). We have \(A^2=B^2=-1\)
and \(D^2=-2\). By adjunction and the pole section on \(A\),
\(L\cdot A=2a-2\), and \(L\cdot D=0\). Subtracting an
integer multiple of the whole fibre, any vertical
representative is \(iA+kD\). The two intersection equations
give \(i=2k\) and \(k=2-2a\). Its class is
\((2-2a)(2[A]+[D])\), whose coefficient of \([E]\) in the
free basis is divisible by four. The coefficient in \(L\)
is two, a contradiction. This completes the hyperelliptic
fibre test in all the width-five genus-two profiles.
