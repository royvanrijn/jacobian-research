# A weighted perimeter criterion at a right contact

Fix \(0<r<1\). Let \(A\) be an ordinary balanced three-cell Gaussian
maximizer in its active plane, and let \(T\) be an oriented, labelled
centered regular simplex. An explicit weighted perimeter of \(A\),
constructed from the simplex scores, must fall strictly below
\(S'_r(T)\) if \(A\) is a nonsimplex tie approached from the right
by simplex-optimal correlations. The criterion uses the full
second variation in one extra Gaussian coordinate.

## A direction supremum for the correlation derivative

Orient each interface \(e=(i,j)\) from the first label to the second.
Write \(d_e=e_i-e_j\), \(N_e\) for this normal, and \(\nu_e\) for Gaussian
surface measure on the interface. For a measurable vector field \(u_e\)
on each interface, put
\[
 m_A^u=\sum_e d_e\otimes u_e\,\nu_e,\qquad
 \mathscr E_r(m)=\|T_{\sqrt r}m\|_{L^2(\gamma_2)}^2,
\]
where the noise operator on a measure denotes integration against its
Mehler kernel and the squared norm sums the tensor coordinates. The
surface energy can also be computed by integrating \(M_r(x,y)\) against
the two copies of the tensor measure. For the actual maximizer,
\[
 S'_r(A)=\sup_{|u_e(x)|\le 1}\mathscr E_r(m_A^u) .             \tag{WC1}
\]
It suffices to use fields with two vector coordinates. Apply the first
extra-coordinate variation (G1) to each scalar coordinate of \(u\)
and sum. The interaction is at most \(r^{-1}\int a|u|^2\,d\nu\), where
\(a_e=|\nabla(T_r1_{A_i}-T_r1_{A_j})|\). This is at most
\(r^{-1}\int a\,d\nu=S'_r(A)\) by (G3). The field \(u=N\) attains equality
by the translation equations (G2). These formulas extend to bounded
measurable fields by the product integrability established with (G3).
Equivalently the identity at \(u=N\) follows by Gaussian integration by
parts in the correlation derivative, which also proves
\(\mathscr E_r(m_T^N)=S'_r(T)\) for the simplex.
[The extra-coordinate form and the translation equations](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)
[The simplex derivative by boundary integration](THREE_CELL_DEFICIT_AND_PERIMETER.md#exact-correlation-derivative-for-the-three-cell-simplex)

Equality in (WC1) requires \(|u|=1\) almost everywhere and each scalar
coordinate to lie in the translation space, by (G17). Conversely a
field satisfying both conditions is an equality case, by (G2).
[The first extra-coordinate gap at any actual maximizer](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#negative-index-and-the-higher-hermite-blocks-at-a-maximizer)

## The simplex weight and its equality case

Let \(U_i^T=T_r1_{T_i}\), set \(B_i=-r^{-1}\nabla U_i^T\), and define
\[
 W_r(A;T)=\sum_{e=(i,j)}\int_{\Sigma_e(A)}|B_i-B_j|\,d\nu_e .
\]
Gaussian divergence through the simplex cells gives
\(B_i=T_r\sum_f(d_f)_i N_f^T\nu_f^T\), using the simplex normals and
surface measures on the right. Thus, for any field \(u\),
\[
 \langle T_{\sqrt r}m_A^u,T_{\sqrt r}m_T^N\rangle
      =\sum_{e=(i,j)} \int u_e\cdot(B_i-B_j)\,d\nu_e .       \tag{WC2}
\]
Choose \(u_e\) in direction \(B_i-B_j\). The positive score derivative
in the pair normal into \(T_i\) shows \(B_i-B_j\ne0\) at every point in
the active plane. In fact the conditional score contrast is the
Gaussian average of \(\operatorname{sgn}(v)\Phi(b+|v|/\sqrt3)\) in
the normal variable \(v\), with a scaled tangential coordinate \(b\).
This function is increasing, with a strictly positive jump at zero.
Its Gaussian average therefore has positive derivative in the normal
mean for every finite observation.
[The reflected pair score](THREE_CELL_CORRELATION_CONTACT.md#pointwise-score-growth)

Cauchy–Schwarz in (WC2) and (WC1) prove
\[
                   W_r(A;T)^2\le S'_r(A)S'_r(T).            \tag{WC3}
\]
At \(A=T\) equality follows because \(B_i-B_j\) is a positive multiple
of the outward pair normal along the corresponding simplex interface.
Suppose \(S'_r(A)\le S'_r(T)\) and \(W_r(A;T)\ge S'_r(T)\). Then all
terms in (WC3) agree, and Cauchy–Schwarz gives identical smoothed tensor
measures. Injectivity of Gaussian smoothing on finite signed measures
gives \(m_A^u=m_T^N\). For injectivity, write the Mehler integral as the
Laplace transform of the measure after multiplication by a strictly
positive Gaussian. Analytic continuation and uniqueness of the Fourier
transform give the claim componentwise. The surface measures are finite
by the Gaussian perimeter bound for maximizers.
[Uniform Gaussian perimeter of maximizers](FINITE_REFLECTION_TEST_FOR_MAXIMIZERS.md#uniform-gaussian-perimeter-and-a-hermite-tail)

The pair surfaces of \(A\) are disjoint up to surface-null sets, and
\(u_e\) is unit everywhere. Equality of the measures therefore makes
every interface of \(A\) lie in the simplex interface set up to surface
measure zero. All three label indicators are distributionally constant
on each open simplex sector. Balance makes their labels a permutation
of the simplex labels. The tensor equality on each simplex ray also
identifies the unordered pair of labels on that ray, so the permutation
fixes every pair and is the identity.

At a tie with simplex-optimal correlations approaching \(r\) from above,
\(S'_r(A)\le S'_r(T)\). Thus a nonsimplex candidate must satisfy
\(W_r(A;T)<S'_r(T)\) for every orientation and labelling. A lower bound
\(\sup_T W_r(A;T)\ge S'_r(T)\) for actual maximizers would exclude these
right contacts. The supremum is attained by continuity on the compact
simplex orbit, since the score gradients are uniformly bounded and
\(A\) has finite Gaussian perimeter. A rotated or relabelled simplex
different from the chosen \(T\) has strict inequality in (WC3) when
its derivative equals \(S'_r(T)\), by the equality argument above.

This deficit is uniform for nonsimplex maximizers with
\(S'_r(A)\le S'_r(T)\) on a compact interior correlation band.
Otherwise choose a sequence for which \(W_r(A;T)\) tends to
the simplex derivative. Rotate the sequence so the reference
simplex is fixed and take a convergent correlation subsequence.
Write \(u_e\) for the unit field in the direction of \(B_i-B_j\).
By (WC1)–(WC3), the squared distance of the two smoothed tensor
measures is at most \(2(S'_r(T)-W_r(A;T))\), which tends to zero.

The individual surface measures have uniformly bounded total mass
on the band. Extract local weak limits of the nonnegative measures
for all three pairs. The fields \(u_e\) converge uniformly on compact
sets by the nonvanishing simplex score derivative. The Mehler kernel
for each fixed observation decays uniformly at source infinity on
the correlation band. Hence the tensor measure obtained from the
weak limits has the same smoothing as the simplex tensor, and is
equal to it by injectivity. On an open set disjoint from the simplex
rays, contract this tensor identity with the matrix of the vectors
\(B_i\) at the limit correlation. The resulting sum of the nonnegative
measures \(|B_i-B_j|\nu_{ij}\) is zero. Thus every weak surface limit
is supported on the simplex rays. Strong maximizer compactness and
the distributional derivative identity put the boundaries of the
ordinary limit partition on these same rays. Balance makes it a
simplex. Uniform local growth contradicts the nonsimplex sequence.
[Compactness and local measurable stability](GAUSSIAN_SIMPLEX_DIMENSION_FREE_LOCAL_GROWTH.md)

## The exact deficit for the weighted perimeter

For the simplex regret \(L_r(A;T)=\sum_i\langle 1_{T_i}-1_{A_i},U_i^T\rangle\),
Gaussian divergence gives
\[
 W_r(A;T)-S'_r(T)
   =\sum_{e=(i,j)}\int
        \big[|B_i-B_j|-(B_i-B_j)\cdot N_e\big]\,d\nu_e
    -\partial_r L_r(A;T).                                \tag{WC4}
\]
Indeed \(\operatorname{div}_\gamma B_i=-r^{-1}L U_i^T=\partial_r U_i^T\),
with \(L=\Delta-x\cdot\nabla\). Its integrals over the cells of \(A\) sum
to the boundary pairing in (WC4). The formula is justified by Gaussian
cutoffs, bounded score gradients, and finite perimeter. The first term
in (WC4) is nonnegative and measures misalignment of the interface
normals with the simplex score gradients. The regret derivative is
also nonnegative by (CT3). Ruling out a right contact through (WC3)
requires controlling this derivative by the normal misalignment after
choosing a reference orientation.
[Pointwise nonnegative growth of simplex score gaps](THREE_CELL_CORRELATION_CONTACT.md#pointwise-score-growth)

## The unrestricted weighted perimeter comparison fails near one

Let \(A\) consist of three concentric Gaussian regions, a disk, an
annulus, and the complement of the outer disk. Take the two radii so
that the three masses equal \(1/3\). For a simplex ray with outer
normal \(n_f\) and outward radial tangent \(\tau_f\), direct integration
over the ray gives
\[
 T_r\nu_f^T(x)=\frac1{\sigma}
       \phi\!\left(\frac{r x\cdot n_f}{\sigma}\right)
       \Phi\!\left(\frac{r x\cdot\tau_f}{\sigma}\right),
       \qquad \sigma=\sqrt{1-r^2}.
\]
Bound the last factor by one and use the representation of \(B_i\)
from (WC2). On a circle of fixed radius \(R>0\), each term is controlled
by
\[
 \frac{R}{\sigma}\int_0^{2\pi}
       \phi\!\left(\frac{rR\sin\theta}{\sigma}\right)d\theta
       \le \frac{C}{r}.
\]
This uses the linear lower bound \(\sin\theta\ge 2\theta/\pi\) on the
first quadrant and is independent of the ray orientation. The Gaussian
surface density is constant on the circle. Thus \(\sup_T W_r(A;T)\)
stays bounded as \(r\uparrow1\), whereas (D2) gives
\(S'_r(T)\sim 3/(2\pi\sigma)\). The unrestricted inequality
\(\sup_T W_r(A;T)\ge S'_r(T)\) fails for this balanced partition.
[The exact simplex slope (D2)](THREE_CELL_DEFICIT_AND_PERIMETER.md#exact-correlation-derivative-for-the-three-cell-simplex)

More generally, if the Gaussian surface measure of a fixed
finite-perimeter partition assigns zero mass to every line through
the origin, then \(\sup_T W_r(A;T)=o(\sigma^{-1})\). In fact the
surface mass of a strip of width \(2\varepsilon\) about a centered
line tends to zero uniformly over its directions. A contrary
sequence of directions has a convergent subsequence, and the
upper limit of the strip indicators is supported on the limiting
line. Finite measure gives the uniform claim. After multiplying
the displayed ray formula by \(\sigma\), the terms outside a fixed
strip vanish uniformly as \(r\uparrow1\). The strip-mass bound
controls the remaining terms. This includes strongly nonradial
partitions and isolates the geometric limitation of an
unrestricted weighted-perimeter inequality.
