# Affine posterior prices and the rank of the model lift

For an ordinary balanced three-cell partition in any finite
Gaussian dimension and \(0<s<1\), the centered simplicial
price map \(\beta(T_s1_A(x))\) is affine in \(x\) exactly when
\(A\) is a regular centered simplex, up to cylinder extension,
labels and null sets. Equivalently its model lift depends on at
most two Gaussian directions exactly in this case. The price
rigidity itself uses no maximization or stationarity assumption.
[Definition of the model lift](GAUSSIAN_WEIGHTED_CENTROID_INEQUALITY.md#a-sharp-posterior-entropy-and-centroid-tradeoff)

## Affine price rigidity

Let \(T_i\) be the centered regular partitions of the model
label plane, and suppose the centered price map is \(Bx+b\)
there. Write \(\sigma=\sqrt{1-s^2}\). The posterior is
\[
 T_s1_{A_i}(x)=\mathbb P(Z-Bx-b\in T_i).
                                                               \tag{PA1}
\]
It is invariant in directions in \(\ker B\), as is \(1_{A_i}\)
by Gaussian convolution injectivity. Indeed translating the
input by \(u\) corresponds to translating the endpoint by
\(s u\) under \(T_s\). One can pass from translation invariance
modulo null sets to a cylindrical representative by smoothing
and letting the smoothing variance tend to zero.
Reduce to the \(\operatorname{rank}B\) active input directions.
Rank zero is impossible since the posterior would be constant
and the input an indicator of mean \(1/3\).

For rank one, take \(u\) unit along \(\operatorname{ran}B\), and
write the label input as \(y u-b_\perp\), where \(y\) ranges over
\(\mathbb R\). Put
\[
 g_i(y)=\mathbb P(yu-b_\perp+Z_\perp\in T_i)
\]
with \(Z_\perp\) standard Gaussian on \(u^\perp\).
After the affine change of the active real variable, (PA1)
takes the form \(P_D f_i=P_1 g_i\), where \(P_D\) is Gaussian
convolution of variance \(D=\sigma^2|B|^2/s^2\) and \(f_i\) are
indicator functions. At least one sector's open interior
contains neither of the two directions on \(u^\perp\).
Its projection on \(u\) then has just one sign, since the
sector is connected and misses the perpendicular line in its
interior. For this \(i\), \(g_i\) vanishes on one open half-line
and is positive on the other. Also \(g_i<1\) at its positive
points: a pointed \(120^\circ\) cone contains no complete
transverse line. If \(D>1\), Gaussian injectivity would imply
\(g_i=P_{D-1}f_i\), which is real analytic and cannot vanish on
a half-line. If \(D\le1\) it would give \(f_i=P_{1-D}g_i\).
The latter is fractional on a positive-measure set, including
when \(D=1\), also impossible.

For rank two the affine change \(y=-Bx-b\) of observation
coordinates yields
\[
 P_D f_i=P_I1_{T_i},\qquad
 D=\frac{\sigma^2}{s^2}BB^{\mathsf T},\qquad
 f_i(y)=1_{A_i}(-sB^{-1}(y+b)).                                 \tag{PA2}
\]
Here \(P_D\) denotes centered Gaussian convolution of covariance
\(D\). We claim \(D\preceq I\). To see this, let \(v_1,v_2\)
generate the two boundary rays of \(T_i\). Off the lines
\(\xi\cdot v_j=0\) its indicator has Fourier transform
\[
      \widehat{1_{T_i}}(\xi)
        =\frac{|\det(v_1,v_2)|}
               {(i\xi\cdot v_1)(i\xi\cdot v_2)} .
\]
For example, multiply the sector indicator in its nonnegative
ray coordinates by \(e^{-\delta(a_1+a_2)}\) and take
\(\delta\downarrow0\) in the distributional transform.
The products of one-dimensional integrals give the displayed
formula away from the two exceptional lines.

Equation (PA2) forces \(\widehat f_i\), away from those lines,
to equal this expression multiplied by
\(e^{\xi^{\mathsf T}(D-I)\xi/2}\). If \(D-I\) has a positive
direction, choose a ray of strict positivity avoiding both
exceptional lines. Pair this smooth expression with translates
of a fixed nonnegative smooth bump at distance \(n\) along
that ray. The sign is constant on these supports, and the
pairings grow at least as \(c n^{-2}e^{c'n^2}\).
This contradicts the polynomial bound for translates of a
fixed test function under a tempered distribution.
Thus \(D\preceq I\).

Gaussian injectivity now gives \(f_i=P_{I-D}1_{T_i}\).
If \(I-D\) is nonzero, a nonzero shift in its support changes
the pointed sector on a set with nonempty interior. There is
then an open set of points where neighborhoods of two shifts
in that support yield both outcomes. The conditional Gaussian
assigns positive probability to both neighborhoods, making
\(f_i\) fractional there. Hence \(D=I\), and (PA2) makes \(A\)
the inverse similarity of a regular model simplex. Its
endpoint masses are the model masses with prices \(\sigma b\),
by homogeneity of the cones and
\(BB^{\mathsf T}=(s^2/\sigma^2)I\).
Equal masses and uniqueness of the centered model prices give
\(b=0\). The converse follows from the Gaussian conditional
rule for a centered simplex.

## Two-dimensionality of the lift

The lifted partition always has two-dimensional model fibers
in the added \(z\)-plane. No nonzero pure-\(z\) shift can
preserve the partition, since each such fiber is a translate
of a pointed fan. If the lift depends on at most two
directions, its inactive subspace in the product space has
dimension at least the dimension of the input. Projection
of the inactive subspace to the input is injective by the
pure-\(z\) observation, hence this subspace is the graph of
a linear map \(V\) from the input into the label plane.

For every fixed input shift \(u\), invariance under
\((u,Vu)\) means that the translated model assignments
at \(x+u\) and \(x\) agree after the \(z\) shift by \(Vu\).
Their centered prices must differ by \(Vu\), by uniqueness
of the translation of a pointed model fan. The equality
holds for almost every \(x\) and extends to all \(x\) by
continuity of the posterior prices. Thus they are affine.
Conversely the lift of a centered simplex is itself
a centered simplex, depending on two directions.

## Variable orientations cannot restore maximality

One can rotate the model fiber by an arbitrary measurable
\(O(x)\in O(H)\), using scores \((O(x)z)_i-\beta_i(T_s1_A(x))\).
This preserves the conditional masses. If the resulting
ordinary partition depends on at most two directions, the
same graph-of-invariances argument applies. For a fixed
input shift \(u\), agreement of the two fibers after a
\(z\)-translation forces equality of their oriented \(i,j\)
normal vectors for every pair. These vectors are
\(O(x+u)^{\mathsf T}d_{ij}/\sqrt2\) and
\(O(x)^{\mathsf T}d_{ij}/\sqrt2\). They agree because two
equal labeled fans have the same open facet between each
pair and the same two sides there. The label differences
span \(H\), so \(O(x+u)=O(x)\) almost everywhere for every
\(u\). By smoothing each matrix entry this gives a constant
orientation almost everywhere. After removing this
constant rotation the affine price argument applies.

In particular, fix \(\rho,s\in(0,1)\) and an actual balanced
three-cell maximizer \(A\) that is not a simplex. The
supremum of noise stability over all such measurable
fiber orientations is strictly below \(S_\rho(A)\).
Indeed dimension reduction gives the same globally optimal
value in the lifted ambient dimension. If the supremum
attained that value along a sequence, weak \(L^2\)
compactness of the simplex-valued assignments and compactness
of \(T_\rho\) would give a balanced relaxed maximizer.
Every relaxed maximizer is ordinary, by strict convexity
on the fixed-mass directions. The \(L^2\) norms of the
indicator tuples are one, so the convergence is then
strong. Take a subsequence converging in Gaussian fiber
\(L^1\) for almost every \(x\). For each such \(x\), the
rotation orbit of the model assignment with fixed masses
\(p(x)\) is compact by dominated convergence off its null
facets. The limit is therefore in the same fiber orbit.
A measurable rotation can be recovered from its
\(z\)-moment matrix, which on the orbit equals
\(L_{p(x)}O(x)/\sqrt2\), by inverting \(L_{p(x)}\) on \(H\).
Thus the limit is itself an oriented model lift and,
as a maximizer, depends on at most two directions.
The preceding rigidity contradicts the choice of \(A\).
[Dimension reduction and ordinary relaxed maximizers](GAUSSIAN_SPECTRAL_DIMENSION_REDUCTION.md#optimization-and-volume-repair)
