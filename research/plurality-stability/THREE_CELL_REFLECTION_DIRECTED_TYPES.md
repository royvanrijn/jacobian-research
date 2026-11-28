# Directed reflection types at equal Gaussian posterior scores

Fix \(0<r<1\) and let \(A\) be an ordinary balanced planar
three-cell Gaussian noise-stability maximizer. Write
\(F=(1_{A_i})_{i=1}^3\) and \(U_i=T_r1_{A_i}\). For a
centered reflection \(R\) and one of its open half-planes \(H\),
a directed mismatch \(i\to j\) is present if the set
\(\{x\in H:x\in A_i,\ Rx\in A_j\}\) has positive Gaussian
measure.

If there is an off-mirror observation \(p\in H\) with
\(U_i(p)=U_i(Rp)\) for all three labels, then at least five
of the six directed mismatches are present. In particular
at least two unordered pairs occur in both orientations
at this prescribed reflection. This strengthens the
three-unordered-type condition obtained by ruling out
two-type sorting.
[Two-type sorting, DG12](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md#exact-comparisons-within-reflection-orbits)
[Targeted reflection at equal-radius triple points](THREE_CELL_POSTERIOR_DEGREE_AND_TRIPLE_BISECTOR.md#the-targeted-reflection-at-equal-vertex-radii)

More generally, let \(p,q\) be any two distinct planar observations
with the same posterior vector, and set \(e=(p-q)/|p-q|\).
For at least five ordered pairs \(i\to j\), a positive-measure
family of lines parallel to \(e\) has a positive two-dimensional
measure of ordered pairs \(t>t'\) with label \(i\) at coordinate
\(t\) and label \(j\) at \(t'\). Thus two missing ordered
slice relations along \(e\) preclude an equal-posterior pair
in that direction. The two opposite orders for an unordered
pair may be witnessed on different lines.

## A cut of the directed mismatch graph

On \(H\), put \(\delta(x)=F(x)-F(Rx)\) and
\(D=\{x\in H:\delta(x)\ne0\}\). The odd Mehler kernel
\(W(x,y)=M_r(x,y)-M_r(x,Ry)\) is strictly positive on
\(H\times H\). Under a reflection-pair swap, the even part of the
stability is unchanged and the odd part is
\[
                  \frac12\iint_{H^2}W(x,y)
                              \delta(x)\cdot\delta(y)
                                    \,d\gamma(x)d\gamma(y).
\]
A swap on a measurable subset \(C\subset D\) changes the sign of
\(\delta\) there. It preserves every cell mass because the labels
within each reflected pair are exchanged and \(R\) preserves
Gaussian measure. The exact stability change is
\[
               -2\iint_{C\times(D\setminus C)}
                           W(x,y)\delta(x)\cdot\delta(y)
                                      \,d\gamma(x)d\gamma(y).
                                                               \tag{DR1}
\]

The six nonzero possible values of \(\delta\) are the roots
\(e_i-e_j\), \(i\ne j\), written in the standard label basis.
They have common norm \(\sqrt2\) and form a regular hexagon
in the centered label plane. Distinct roots have positive dot
product exactly when they are neighbors on the hexagon. All
other dot products of distinct roots are negative.

The directed types present on positive measure must induce
a connected subgraph of this hexagon. Otherwise choose \(C\)
as the orbits with root in one connected component. Both sides
have positive measure, and every cross dot product is strictly
negative. Strict positivity of \(W\) in (DR1) gives a finite
improvement. The same argument works in any Gaussian dimension
and requires no score information.

For an equal-posterior point \(p\in H\),
\[
                   0=(U_i(p)-U_i(Rp))_{i=1}^3
                       =\int_D W(p,y)\delta(y)\,d\gamma(y).
                                                               \tag{DR2}
\]
This is a positive weighted sum of all present roots. In the
plane \(D\) has positive measure, since reflection invariance of
all three labels would place their centroids in the mirror line,
whereas the two active centroid directions span the plane.
The same reasoning applies to any reflection whose normal has
a nonzero projection on the active plane of the cylinder used
below.
[Full active-plane centroid rank](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible)

A nonempty connected proper subset of the hexagon induces
a path of consecutive vertices. With one or two vertices the
roots lie inside some open half-plane through the origin. With
three or four consecutive vertices they lie in a closed half-plane
through the origin and at least one lies in its interior. A
positive weighted sum of the roots cannot vanish in any of
these cases. Equations (DR1) and (DR2) therefore require at
least five directed types. An observation on the mirror is
outside this argument, since the odd kernel there vanishes.

## Ordered slices for arbitrary equal-posterior observations

The cylinder \(\widehat A_i=A_i\times\mathbb R\) in \(\mathbb R^3\)
is a global maximizer by Gaussian dimension reduction, with
scores independent of the added coordinate. Choose \(0<\alpha<1\),
set \(b=\sqrt{1-\alpha^2}\), and reflect through the centered
plane with normal \(n=(\alpha e,b)\). For \(m=(p+q)/2\) and
\(a=|p-q|/2\), put
\[
 \widehat m=(m,-\alpha(m\cdot e)/b),\qquad
 \widehat p=\widehat m+(a/\alpha)n,\qquad
 \widehat q=\widehat m-(a/\alpha)n .
\]
The midpoint lies in the reflection plane, and the observations
\(\widehat p,\widehat q\) project to \(p,q\). Hence they are
exchanged by the reflection and have the same posterior vector.
Apply the directed type bound in the half-space containing
\(\widehat p\). Its nonempty mismatch support follows from the
centroid rank, since all cylinder centroids are planar and
\(n\) projects onto a nonzero multiple of \(e\).
[Heilman–Tarter, Theorem 1.9](https://doi.org/10.1017/fms.2021.56)

Write a general point of the cylinder as \((z+te,h)\), where
\(z\in e^\perp\), and denote the projected reflected coordinate
by \(t'\). Then
\[
                        t'=(1-2\alpha^2)t-2\alpha b h .
\]
The reflection keeps \(z\) fixed, and the chosen half-space is
equivalent to \(t>t'\). Because \((z,t,h)\mapsto(z,t,t')\)
is invertible, the pushed-forward Gaussian measure has a
positive density. For each present directed type, Fubini
gives a positive-measure set of \(z\) with positive
two-dimensional measure of ordered \((t,t')\) realizing the
corresponding labels. Thus the claim strengthens the
unordered overlap witness using the same cylindrical reflection.
[Unordered conditional slice witness](THREE_CELL_REFLECTION_SCHUR_MISMATCH.md#a-conditional-slice-witness-without-equal-radii)


## Critical posterior directions and ordered half-lines

For any fixed unit vector \(e\), the ordered slice relations
in this direction form a connected subgraph of the hexagon.
Indeed the same tilted reflection of the cylinder realizes those
relations whether or not a posterior anchor is given. In either
direction of this equivalence, positive ordered two-dimensional
measure on a positive-measure set of \(z\) is equivalent by Fubini
and the invertible coordinate change to positive reflected
mismatch measure. The set of relations is nonempty, since
an empty set would make the cylinder reflection invariant,
contrary to centroid rank.

In fact all directional derivatives of the posterior are positive
combinations of the present ordered roots. Write
\(\sigma^2=1-r^2\), let \(x_e=x\cdot e\) and
\(x_\perp=x-x_e e\), and denote the scalar Gaussian density with
mean \(m\) and variance \(\sigma^2\) by \(g_{m,\sigma}\).
Then
\[
 \begin{split}
 \partial_e U(x)=\frac{r}{\sigma^2}
 \int_{e^\perp}\int_{t>t'} &(t-t')
      [F(z+te)-F(z+t'e)]\\
    &\qquad g_{rx_e,\sigma}(t)\,g_{rx_e,\sigma}(t')\,
            g_{r x_\perp,\sigma}(z)\,dt\,dt'\,dz .      \tag{DR3}
 \end{split}
\]
To derive the formula, differentiate the
conditional density of the longitudinal coordinate given
\(X=x\). Its score factor is \(r(t-rx_e)/\sigma^2\).
For each transverse \(z\), draw independent \(t,t'\) from
\(g_{rx_e,\sigma}\). The expectation of
\((t-t')[F(z+te)-F(z+t'e)]\) is twice the covariance of
\(t\) and \(F(z+te)\). This integrand is symmetric in
\(t,t'\), so integrating only over \(t>t'\) gives that
covariance. Gaussian integrability justifies both the
differentiation and the double integral.

The coefficient of each present ordered root in (DR3) is
strictly positive, because the densities and \(t-t'\) are
positive on its ordered set. If \(\partial_e U(x)=0\) at
even one observation, the closed-half-plane argument for
(DR2) now requires at least five ordered relations along \(e\).

Conversely, at most four ordered relations occur along \(e\)
exactly when some label occupies an upper or a lower half-line
in almost every \(e\)-parallel fiber, with a common orientation
and with empty and full fibers allowed. For the forward
implication, extend the connected path of present roots to a
path of four consecutive vertices, by enlarging its endpoints
so that some present vertex is interior to the four-vertex
path. For a path of length four in vertices, its endpoints
represent opposite orientations of one unordered pair. The
two interior roots both have the remaining label as their
source, or both have it as their target. Consequently the two
opposite orders involving that remaining label are absent.
On almost every fiber, absence of both such orders means
that this label is an upper or lower half-line, modulo a null
set. For upper-half-line ordering, if the essential
infimum of the label set were below the essential supremum
of its complement, a threshold between them would give
positive label measure below it and positive complementary
measure above it. This gives the forbidden ordered pairs
in that fiber. Otherwise the set is an upper half-line
up to a null set. Lower-half-line ordering reverses the
inequalities. Conversely
a common half-line orientation excludes two oriented
relations, one for each of the other labels. The chosen path
has at least one present interior root, so the half-line
label's posterior has strict longitudinal derivative of the
corresponding sign at every point, by (DR3). In particular
there is both injectivity along individual \(e\)-parallel
lines and absence of any Jacobian kernel in direction \(e\).

Now let \(C\) be a bounded open strict-winner component
and use the two-score map
\(\mathcal F=(U_1-U_3,U_2-U_3)\). The componentwise degree
formula gives opposite Jacobian signs among the preimages
of any regular value in \(\mathcal F(C)\). A path in \(C\)
between such preimages contains a point with singular
Jacobian. If \(\det D\mathcal F\) vanishes identically on
\(C\), every point is already singular. These cases cover
all possibilities. Otherwise a point with nonsingular
Jacobian has an open inverse chart, and the image of the
critical set has zero planar measure. Indeed the critical
set of the nonzero analytic determinant is null on \(C\),
and \(\mathcal F\) is locally Lipschitz. The inverse chart
then supplies a regular value in \(\mathcal F(C)\).

At any such singular point there is a unit null direction
\(e\) with \(\partial_e U_i=0\) for all labels, since the
two contrasts have derivative zero and \(\sum_i U_i=1\).
Thus every bounded strict-winner component forces at least
five ordered relations along one of its interior posterior
null directions. If every direction has a common half-line
label on its slices, the posterior map is globally injective
and has a nonsingular Jacobian everywhere. Such a partition
has no bounded strict-winner component.
[Componentwise degree](THREE_CELL_POSTERIOR_DEGREE_AND_TRIPLE_BISECTOR.md#degree-zero-inside-each-bounded-winner-component)


## A third posterior preimage for compact simplex discrepancies

Suppose in addition that \(A\) agrees almost everywhere
outside a disk with a centered regular simplex partition \(T\).
View \(U=(U_1,U_2,U_3)\) as a map into the open probability
triangle, and orient its target by identifying \(u\) with
\(\sum_i u_i v_i\), where \(v_i\) are the unit centers of the
three sectors of \(T\). Every globally regular posterior
value has signed preimage count \(+1\) in this orientation.
Moreover, for any bounded strict-winner component \(C\)
there are globally regular posterior values in \(U(C)\).
Each such value has opposite-sign preimages within \(C\)
and at least one further preimage in an unbounded
strict-winner component of the same label.

For properness, the three closed reference sectors intersect
only at the origin. Thus
\[
         c=\min_{|\omega|=1}\max_i\operatorname{dist}(\omega,T_i)>0.
\]
For each observation \(x\ne0\) at least one reference sector is
at distance \(rc|x|\) or greater from \(rx\). Its reference
posterior probability tends uniformly to zero as \(|x|\to\infty\).
The difference between any label's posterior for \(A\) and \(T\)
is bounded by the conditional probability of the fixed
exceptional disk and also tends uniformly to zero. It follows
that the linear homotopy of the two posterior maps is uniformly
proper into the open probability triangle. Every compact set
in that triangle is bounded away from a zero coordinate,
so all its preimages under the homotopy stay in a common disk.

For the reference map let
\[
          \Psi(x)=\mathbb E\max_i v_i\cdot (rx+\sigma Z),
          \qquad Z\sim\gamma_2 .
\]
Differentiation of the Lipschitz maximum shows
\(\nabla\Psi=r\sum_i v_i U_i^T\). If \(N_{ij}\) is the
unit normal from \(T_i\) to \(T_j\) and \(\Gamma_{ij}\) their
common ray, differentiation of the piecewise constant gradient
of the maximum in distributions gives, for a vector \(h\),
\[
 D^2\Psi(x)[h,h]
    =r^2\sum_{i<j}|v_i-v_j|(h\cdot N_{ij})^2
                     \int_{\Gamma_{ij}}\phi_\sigma(y-rx)\,ds(y) .
\]
Here \(\phi_\sigma\) is the centered planar Gaussian density
of variance \(\sigma^2 I\). The interface term is the jump
of the gradient in the normal direction, followed by Gaussian
convolution. The three normals span the plane and each
ray integral is strictly positive. Hence the Hessian is
positive definite at every \(x\). The reference posterior in
the stated target coordinates has positive Jacobian and is
injective by integrating its derivative on line segments.
It is a local diffeomorphism, with image both relatively open
and relatively closed in the target, the latter by properness.
The target is connected, so the reference map is a
diffeomorphism onto it.

For completeness the degree comparison follows directly from
Stokes' theorem. In target coordinates take a smooth two-form
of compact support in the open triangle with integral one.
Pull it back by the homotopy on \([0,1]\times\overline B_R\),
choosing \(R\) large enough that all preimages of the support
avoid the lateral boundary. The two-form is closed, so the
integrals of its pullbacks by the two endpoint maps agree.
The reference integral is one. For a globally regular value
of \(U\), choose the support in a sufficiently small ball
around that value. Properness and local invertibility make
its fiber finite and ensure that nearby preimages stay in
inverse charts for that fiber. The integral for \(U\) is then
the sum of the Jacobian signs in those charts, proving the
count \(+1\).

The same integral identity shows that the Jacobian of \(U\)
cannot vanish identically. It is analytic on the plane, so its
critical set has measure zero, as does the image of that set
by local Lipschitz continuity. Each nonempty open component
\(C\) contains a point of nonsingular Jacobian, giving an open
target patch in \(U(C)\). Hence globally regular values
occur in \(U(C)\). For such a value the adjusted winner label
is determined by the common posterior vector, so all preimages
are strict winners of the same label. Each bounded winner
component with a preimage contributes signed count zero by
the componentwise degree formula. Only finitely many components
meet the fiber, since the fiber is finite. The full count
\(+1\) therefore comes from unbounded components, and \(C\)
itself contributes both signs. The resulting preimages give
the directed slice restrictions above along each of their
pairwise displacements.
