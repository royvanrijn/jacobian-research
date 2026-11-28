# Affine symmetries at a three-cell Gaussian junction

Fix \(0<\rho<1\). Let \(A\) be a balanced planar Gaussian
noise-stability maximizer, and suppose an invertible affine map
\(F(x)=p+L(x-p)\) preserves it up to a permutation of the
labels and null sets. If \(L\) is contracting in spectral
radius, then \(A\) is the centered regular simplex. The same
holds if \(L^{-1}\) is contracting. More generally it suffices
that \(p\) is a genuine triple junction and at least one
eigenvalue of \(L\) has modulus different from one.
In these cases \(p=0\), and \(L\) is a positive scalar times an
orthogonal symmetry of the simplex.

The proof uses the ordinary analytic jump geometry of
actual maximizers. It covers a discrete affine symmetry and
does not require the linear part to be symmetric, normal, or
diagonalizable.
[Planar regular arcs and the two triple-junction cases, G20](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)
[Cones centered at an arbitrary point](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#cones-centered-at-an-arbitrary-point)

## The induced linear action on the junction

After passing to a positive power of \(F\), its label
permutation is trivial. If \(L\) is contracting and any cell
has positive mass, it has a positive-volume subset in some
finite ball. Its images under large powers of \(F\) are still
in the same cell up to null sets, have positive volume by
invertibility, and lie in every sufficiently small ball
about \(p\). Thus all three open strict winners approach
\(p\). For expanding \(L\) apply this argument to \(F^{-1}\).
It remains to work at a genuine triple junction.

Because the symmetry acts by an invertible affine
transformation, it preserves the essential boundaries
of the labeled cells. In particular, it permutes the finite
set of analytic half-arcs of the actual jump near the fixed
point, as local germs at that point. Pass to another positive
power so that it fixes every half-arc germ. Such a power
preserves each oriented tangent ray. It may also be taken
even to remove negative eigenvalues in the real
eigenvalue case. These powers preserve the assumption that
some eigenvalue is off the unit circle.

If the triple junction is noncollinear, its three tangent
rays have distinct unoriented lines. A planar linear map
fixing three distinct lines is scalar, and since it fixes
the oriented rays the scalar is positive. The scalar is
different from one by the eigenvalue hypothesis.

In the collinear case the map preserves the common tangent
line. If the map is scalar we again obtain a nontrivial
positive scalar symmetry. Otherwise write coordinates
\((x,y)\) about \(p\), with the tangent equal to the
\(x\)-axis. Every incident regular wall continues locally
as an analytic graph \(y=h(x)\) with \(h(0)=h'(0)=0\),
by nonvanishing of its pair score gradient. The genuine
collinear cusp has at least one nonzero such graph.

If both eigenvalues are the same positive number
\(\lambda\), a nonscalar map has the form
\((x,y)\mapsto(\lambda x+a y,\lambda y)\) with \(a\ne0\).
Invariance of a wall germ gives
\(h(\lambda x+a h(x))=\lambda h(x)\) on its active
half-interval. If the first nonzero term of \(h\) has
order \(q\ge2\), comparing that term gives
\(\lambda^q=\lambda\). This contradicts the assumption
that an eigenvalue is off the unit circle.

If the eigenvalues are distinct, take the possibly
oblique eigenbasis, with the first direction tangent.
The map is \((x,y)\mapsto(\lambda x,\mu y)\) with
\(\lambda,\mu>0\). Each wall equation is
\[
                         h(\lambda x)=\mu h(x).               \tag{AS1}
\]
The identity holds on an active half-interval and hence
as an analytic germ. Comparison of Taylor coefficients
shows that every wall is either zero or \(h(x)=a x^q\)
for the same integer \(q\ge2\), with \(\mu=\lambda^q\).
A nonzero wall exists, so \(\lambda\ne1\) under the
eigenvalue hypothesis. Both eigenvalues are simultaneously
contracting or expanding. Invert the map if necessary,
so \(0<\lambda<1\).

## A nontrivial power-law wall is impossible

On each sign half-plane of \(x\), put \(z=y/|x|^q\).
Iteration into the junction and invariance of \(z\)
show that the ordinary labels on this half-plane
depend only on \(z\), with the finitely many thresholds
of the local monomial walls. For clarity, for almost
every point there is equality of labels under every
iterate, and some iterate reaches the local chart.
Within a fixed sign half-plane in that chart, the
thresholds are constants in \(z\) by (AS1), so the same
thresholds hold at every initial point, modulo null
sets.

Choose one half-plane with a nonzero threshold \(c\).
If it has positive thresholds, take the largest
positive one. Otherwise take the smallest negative
one. Its incident pair has a constant unadjusted
score contrast on the entire wall
\((x,c|x|^q)\) with the chosen sign of \(x\). This
follows either from the reconstructed pair of
winning regions or from real analyticity of the
score contrast restricted to the wall.

At an observation \(p+xv_t+c|x|^q v_n\) on this
wall, where \(v_t,v_n\) are the two eigenvectors,
the conditional source is
\(\rho(p+xv_t+c|x|^qv_n)+\sqrt{1-\rho^2}\,Z\).
In the same oblique coordinates about \(p\), its
horizontal sign eventually agrees with that of \(x\)
and its invariant coordinate tends almost surely to
\[
                                \rho^{1-q}c .
\]
Indeed the coordinates of the bounded affine drift
and a fixed Gaussian sample contribute \(O(1)\),
while the leading horizontal and transverse scales
are \(|x|\) and \(|x|^q\). This computation does not
require independent Gaussian coordinates in the
oblique basis.

The limit lies beyond the chosen extreme threshold.
The outside label for that wall consequently has
conditional probability tending to one, by dominated
convergence, and the other incident label tends to
zero. Their contrast tends to \(1\) or \(-1\). Since
the contrast is constant on the wall and both cells
have positive Gaussian mass, its value at each finite
observation is strictly between \(-1\) and \(1\).
This is impossible. Thus the collinear case cannot
give a nonscalar affine symmetry under these hypotheses.

## Scalar invariance and uniform separation

An actual stationary partition invariant under a
nontrivial positive scalar dilation about \(p\) is
conical about \(p\). Indeed, on almost every ray from
\(p\) all pair-score differences have some first
nonzero Taylor coefficient at \(p\) in the ray direction.
For each pair the exceptional directions lie in the
zeros of a nonzero homogeneous polynomial. Thus the
ordering of the scores on that ray is constant at
sufficiently small positive distances from \(p\).
Dilation invariance, after inversion if necessary,
propagates the small-distance winner to almost every
point of the ray. Fubini justifies using the invariant
indicator identities simultaneously along almost
every ray. The balanced translated-cone classification
now gives the centered simplex. Its unique triple
junction is the origin, so the fixed junction \(p\)
equals zero. The original \(L\) maps the three oriented
simplex rays to positive multiples of one permutation
of them. Their only linear relation has three equal
coefficients, so those positive multiples agree and
\(L\) is a scalar times a simplex symmetry.

There is a uniform version for contracting maps.
Fix a compact positive-correlation band, a bounded
set of centers \(p\), and a compact set of invertible
linear maps with spectral radius bounded above by
a number less than one. Then a nonsimplex balanced
planar maximizer has a positive uniform lower bound
on its Gaussian label distance from its affine
transform by any of these maps, even after relabeling.
Otherwise take convergent affine parameters and
correlations, fix the permutation on a subsequence,
and use strong compactness of actual maximizers.
The limit remains nonsimplex by uniform local
simplex isolation on the band. Affine pullback is
continuous in this limit: the transformed Gaussian
laws have uniformly small tails and uniformly
bounded densities on compact sets, and bounded
continuous approximations to the limiting indicators
give strong convergence under the convergent maps.
Thus the limit would have the exact affine symmetry
already excluded.
[Strong compactness and uniform local comparison](POSITIVE_CORRELATION_CONTINUATION.md)

The result excludes exact contracting affine
self-similarity of a nonsimplex tied maximizer.
A partition can have non-self-similar interfaces,
which are not classified by this argument.
