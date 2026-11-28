# Ordinary compact partitions with all three tilted slice identities

For every fixed \(0<r<1\), there is an ordinary balanced planar
three-cell partition that agrees with a centered simplex off a
compact annulus, obeys all three identities (CD1) at \(r\), and has
nonzero spectator imbalance for every reference pair. Its three
angular masses are pairwise distinct on a positive-measure set of
radii. This gives an ordinary realization of the fractional
construction (CD7).
[The fractional construction and its nonzero witnesses](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#projection-constraints-allow-fractional-asymmetry)

## Which constraints can be preserved

Let \(D\) be a closed disk centered at the origin, and let
\(f_1,\ldots,f_k:D\to[0,1]\) satisfy \(\sum_i f_i=1\) almost
everywhere. There is an ordinary assignment on \(D\) preserving all
of the following data of the \(f_i\):

- The pushforward of \(e^{-\alpha|x|^2}f_i(x)\,dx\) by \(x\mapsto
  e\cdot x\), for every pair \((e,\alpha)\) in any prescribed finite
  list of unit vectors and real exponents.
- The pushforward by \(x\mapsto |x|\) of \(f_i(x)1_C(x)\,dx\), for
  every sector \(C\) of any fixed finite angular division of the
  disk.
- Any prescribed finite list of linear integrals of the tuple
  against \(L^1(D)\) test functions.

The data are preserved label by label in the first two items.
Assignments outside \(D\) can be left as given. The sectorwise
radial pushforwards already preserve all angular masses by circle
and the Gaussian masses in each such sector at any variance.

## A bounded perturbation with zero directional pushforwards

Let \(E\) be a bounded measurable subset of the plane with positive
area. Index the prescribed directional pushforwards by
\(j=1,\ldots,m\) and put \(w_j(x)=e^{-\alpha_j|x|^2}\). Choose small
nonzero vectors \(d_j\perp e_j\) and define the weighted translations
\[
 U_j h(x)=\frac{w_j(x-d_j)}{w_j(x)}h(x-d_j).
\]
For a compactly supported bounded \(h\), the \(j\)-th projected
measure of \((U_j-I)h\) vanishes. Multiplying by \(w_j\) makes it the
difference of a translate of \(w_jh\) by \(d_j\) and \(w_jh\) itself,
and the longitudinal coordinate is fixed by this translation.

The operators commute up to constant factors:
\[
 U_iU_j=\exp(2(\alpha_i-\alpha_j)d_i\cdot d_j)\,U_jU_i.
\]
Use the \((N+1)^m\) ordered monomials
\(U_1^{n_1}\cdots U_m^{n_m}\), \(0\le n_i\le N\), with \(N=m\).
For each \(j\), commute \(U_j\) to the left. In each group with
all other exponents fixed, require that the resulting polynomial
in \(U_j\) have coefficient sum zero. Each such group imposes one
linear condition and then has the left factor \(U_j-I\).
There are at most \(m(N+1)^{m-1}\) conditions, so there is a
nonzero real coefficient array satisfying all of them. The
resulting operator \(P\) has a left factor \(U_j-I\) for every
\(j\), separately.

The step lengths can be chosen arbitrarily small so that the
vectors \(D_n=\sum_i n_i d_i\) for these monomials are distinct.
For each two different exponent tuples the forbidden coincidence
is a proper linear subspace of the step-length parameters, and
there are finitely many such coincidences. Continuity of
translations in \(L^1\) ensures that the intersection
\(\bigcap_n(E-D_n)\) has positive measure for sufficiently small
steps. Choose a positive-measure subset \(F\) of that intersection
contained in a ball of diameter less than the smallest distance
between distinct \(D_n\). The supports of the monomials applied to
\(1_F\) are then disjoint subsets \(F+D_n\) of \(E\). Their
multipliers are positive and bounded on these supports. Therefore
\(h=P1_F\) is a nonzero bounded function supported in \(E\) and
all its directional pushforwards vanish. For an empty list one
can use \(h=1_E\).

## Radial constraints and ordinary realization

To preserve the radial pushforwards as well, shrink \(E\) to a
positive-measure set \(E_0\) in an open ball \(B\) strictly within one
angular sector and away from the origin. Choose \(B\) sufficiently
small that it contains no whole circle centered at the origin.
Take a rotation \(R\) through a sufficiently small irrational
multiple of \(2\pi\). Continuity of rotations in \(L^1\) gives
\(|E_0\cap R^{-1}E_0|>0\). Apply the previous construction on this
intersection to the original directions and their images under
\(R^{-1}\), with the corresponding weights. Put
\[
                         g(x)=h(x)-h(R^{-1}x).
\]
Both terms have zero directional pushforwards in the original
list, since each weight is radial. Both are supported in \(E_0\),
so all sectorwise radial pushforwards of \(g\) vanish by rotation
invariance. If \(g=0\), then \(h\) is invariant under \(R\) in
\(L^2\). The powers of \(R\) are dense in the rotation group, whose
action is strongly continuous in \(L^2\), so \(h\) would be radial.
No nonzero radial function is supported in \(B\). Hence \(g\ne0\).

Any finite additional scalar constraints can now be imposed by
splitting \(E\) into more disjoint positive-measure pieces than
there are constraints, constructing a \(g\) on each, and taking
a nonzero linear combination in the kernel of the constraint
matrix. The disjoint supports ensure that this combination is
nonzero. Thus every positive-measure fractional zone of two labels
supports a bounded nontrivial perturbation preserving the desired
data when it is added to one label and subtracted from the other.

Consider the feasible set of tuples on \(D\) with the prescribed
data, nonnegative coordinates and coordinate sum one. It is convex
and weak-* compact in \(L^\infty(D)^k\), by Banach–Alaoglu.
Pushforward constraints are weak-* closed by testing with bounded
measurable functions of their projected variable. The remaining
test functions belong to \(L^1(D)\). The starting fractional
tuple makes the set nonempty, and Krein–Milman gives an extreme
point. If the point were fractional on positive measure, some two
labels would both lie in \([\varepsilon,1-\varepsilon]\) on a common
positive-measure set \(E\). A sufficiently small multiple of the
perturbation just constructed, and its negative, would give two
different feasible tuples with the extreme point as their mean.
Consequently every extreme point is an ordinary assignment.

## The annular examples and their score winner

Choose an annulus whose outer radius is sufficiently small for the
sign calculation of the odd perturbation in (CD7). The three radial
cutoffs used for the odd and even perturbations can all be supported
inside that annulus, with the two even cutoffs disjoint. Use the
uniform fractional tuple on the annulus and the exact simplex
elsewhere as the baseline. Its pair contrasts are still odd in the
relevant normal coordinate, since the annulus is radial. Use the
same odd and even perturbations as in (CD7), with label coefficients
\((1,2,-3)\). Their sum is zero and each is nonzero. Thus the same
strict half-mass witnesses apply to all three labels, and the
nonzero even angular mean now gives three distinct angular masses
wherever it is nonzero.

Take the three reference projection axes with
\(\alpha=1/(2(1-r^2))\), and divide the disk containing the annulus
into the six angular sectors cut by the axes. The directional
pushforwards give (CD1). The sectorwise radial data preserve standard
Gaussian mass, all three spectator imbalances, and each angular
mass by radius. On every unmodified circle the reference simplex
uses just one label in each of the six sectors, so preservation of
the sectorwise masses also gives the simplex assignment there up
to null sets. This leaves the exact simplex germ at the origin and
the exact exterior. The annulus can be chosen within any prescribed
positive outer radius.

For a sufficiently small annulus, the zero-price scores of the
ordinary realization still choose the simplex at almost every
point. In fact the compact pair contrast change has smoothed value
and first derivatives bounded uniformly by a constant depending on
\(r\) times the area of the annulus. This follows by bounding the
conditional Gaussian density and its first derivatives in supremum
norm. The smoothed change vanishes on the reference pair line by
(CD1), so its absolute value is bounded by that same constant
times area times \(\min\{\operatorname{dist}(x,L_{ij}),1\}\).
The reference score margin on \(T_i\) is at least a positive
constant depending on \(r\) times this distance factor. For small
area its sign is unchanged. Thus the projection identities and
genuine exterior and central score geometry coexist with a defect
that fails the winner rule on its changed cells.
[The saturated reference score margin](POSITIVE_CORRELATION_CONTINUATION.md#relative-exterior-control-below-every-fixed-correlation-bound)
