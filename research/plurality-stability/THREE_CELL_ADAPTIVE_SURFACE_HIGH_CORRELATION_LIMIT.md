# The high-correlation limit of the adaptive Gaussian surface energy

The optimal adaptive surface energy for balanced planar three-cell
partitions has the same leading coefficient as the simplex as
\(r\uparrow1\). Any sequence attaining this coefficient converges
in Gaussian label distance to the simplex orbit. The lower bound
already holds for scalar directions. It uses compactness at the
Mehler scale and the planar decomposition of a finite-perimeter
boundary into Jordan curves.
[The two-direction functional](THREE_CELL_ADAPTIVE_SURFACE_VARIATIONAL_PROBLEM.md)
[The unrestricted direction functional](THREE_CELL_ADAPTIVE_SURFACE_UNRESTRICTED_DIRECTIONS.md)

For a planar partition \(A\), let
\(\mu_{ij}=\mathcal H^1|_{\partial^*A_i\cap\partial^*A_j}\),
\(\nu_{ij}=\phi_2\mu_{ij}\), and \(d_{ij}=e_i-e_j\), where
\(\phi_2\) is the standard planar Gaussian density. Write
\(\mu=\sum_{i<j}\mu_{ij}\) and \(P_\gamma(A)=\sum_{i<j}\nu_{ij}(\mathbb R^2)\).
For \(D\in\mathbb N\cup\{\infty\}\), let \(F_{D,r}(A)\) be the
supremum of \(\|T_{\sqrt r}\sum_{i<j}d_{ij}\otimes u_{ij}\nu_{ij}\|_2\)
over measurable fields \(|u_{ij}|\le1\) in \(\mathbb R^D\), or in
\(\ell^2\) for \(D=\infty\). Assign infinite energy to partitions of
infinite Gaussian perimeter.

Put \(\sigma_n=\sqrt{1-r_n^2}\) for \(r_n\uparrow1\). If
\[
             \sup_n\sigma_n F_{D_n,r_n}(A^{(n)})^2<\infty ,
                                                               \tag{HC1}
\]
then \(A^{(n)}\) is precompact in Gaussian label measure. Every limit
\(A\) has finite Gaussian perimeter and along the subsequence converging
to \(A\),
\[
 \liminf_n \sigma_n F_{D_n,r_n}(A^{(n)})^2
                  \ge \sqrt{\frac2\pi}\,P_\gamma(A).
                                                               \tag{HC2}
\]
The direction ranks \(D_n\) can vary and include \(1\). The balanced simplex
attains the coefficient on the right of (HC2).

## Compactness at the Mehler scale

Let \(\sigma=\sqrt{1-r^2}\). The joint Gaussian Lebesgue density is
\[
 j_r(x,y)=\frac{1}{(2\pi)^2\sigma^2}
   \exp\!\left(-\frac{|x|^2+|y|^2}{2(1+r)}
               -\frac{r|x-y|^2}{2\sigma^2}\right).
                                                               \tag{HC3}
\]
It equals \(M_r(x,y)\phi_2(x)\phi_2(y)\) for the relative
Mehler kernel \(M_r\). The Mehler semigroup therefore gives the
\(L^2(\gamma_2)\) pairing of smoothed Gaussian surface measures by
integrating \(j_r\) against their Euclidean surface measures.

Put \(h_{ij}=T_{\sqrt r}\nu_{ij}\), which is in \(L^2\) at finite
energy by testing the tag alone. Choose independent uniform scalar signs
on the tags. The expected squared norm of the resulting tensor is
\(2\sum_{i<j}\|h_{ij}\|_2^2\), so
\[
 F_{D,r}(A)^2\ \ge\ \frac23\left\|\sum_{i<j}h_{ij}\right\|_2^2
                =\frac23\iint j_r(x,y)\,d\mu(x)\,d\mu(y).
                                                               \tag{HC4}
\]

On a grid of side length \(\sigma_n\), assign each square a label
of largest area in \(A^{(n)}\). For an open grid square \(Q\),
each other label has area at most \(|Q|/2\). The relative
isoperimetric inequality in the square bounds the sum of their
areas by \(C\mu^{(n)}(Q)^2\). On the squares of any fixed bounded
region, \(j_{r_n}\ge c/\sigma_n^2\) on \(Q\times Q\), by (HC3).
Summing (HC4) on these disjoint products bounds the total area
error of the grid approximation by
\(C\sigma_n^2 F_{D_n,r_n}(A^{(n)})^2\), which tends to zero by (HC1).

For two adjacent squares with different largest labels, take the
open rectangle \(R\) that is the interior of the union of their
closures. Each of the two labels has area at least \(\sigma_n^2/3\)
in \(R\), so relative isoperimetry gives \(\mu^{(n)}(R)\ge c\sigma_n\).
The products \(R\times R\) of adjacent rectangles have bounded
overlap. The lower bound \(j_{r_n}\ge c/\sigma_n^2\) also holds
on these products in a fixed bounded region. Hence the Euclidean
perimeter there of the grid partition is bounded by
\[
 C\sigma_n^{-1}\sum_R \mu^{(n)}(R)^2
                    \ \le\ C\sigma_n F_{D_n,r_n}(A^{(n)})^2 .
\]
Here the sum can include the adjacent rectangles in a fixed bounded
enlargement of the region. Local \(BV\) compactness of the grid
partitions and the vanishing area error give a subsequence of
\(A^{(n)}\) with a locally finite-perimeter limit in local label
measure. Exhaustion and the Gaussian tail of planar label measure
give Gaussian label convergence. Prescribed cell masses pass to
the limit.

## A separating path inside a ball

Work in coordinates centered at \(p\), let \(H\) be the upper
half-plane, and suppose the locally finite-perimeter set \(E\) obeys
\[
       |(E\triangle H)\cap B(p,a)|\le \eta^2a^2,
                           \qquad 0<\eta<1/20 .
                                                               \tag{HC5}
\]
The same argument works for the lower half-plane. Select a rectangle
\(Q\) with left and right side distances between \((1-2\eta)a\) and
\((1-\eta)a\) and each height between \(\eta a/2\) and \(\eta a\).
The full range of rectangles lies inside the ball since
\((1-\eta)^2+\eta^2<1\). On almost every coordinate line the two
\(BV\) traces of \(E\) agree with the sliced indicator almost
everywhere, and the line has zero perimeter measure. Integrate
the trace mismatch with \(H\) over the full ranges of possible
side positions, using the longest possible side segments.
Fubini and (HC5) allow each vertical side to be chosen with
mismatch length at most \(2\eta a\), and each horizontal side
with mismatch at most \(4\eta a\). These lengths are twice
the respective upper bounds from averaging. Thus the total
mismatch length on \(\partial Q\) is at most \(m=12\eta a\).

Keep \(E\) inside \(Q\) and set it equal to \(H\cap(-2a,2a)^2\)
outside \(Q\) in the same coordinates. The patched set is bounded
and has finite perimeter. The perimeter on \(\partial Q\) is
the trace-mismatch length by \(BV\) gluing. Its boundary
outside \(\overline Q\) is one open polygonal arc whose
two endpoints are the intersections of \(\partial H\) with the
vertical sides of \(Q\).

In the planar Jordan decomposition of the measure-theoretic
boundary, some rectifiable Jordan curve meets this exterior arc.
Every decomposing curve is contained in the support of the
perimeter measure, since it lies on the measure-theoretic boundary
almost everywhere in arclength and its own arclength measure
has full support on it. At each exterior point the local
perimeter support is the polygonal arc. Local injectivity of
the Jordan curve parametrization makes its intersection with
the open arc relatively open. The curve is closed, so the
intersection is also relatively closed and contains the whole
exterior arc. The complementary closed path on this curve
connects the two endpoints in \(\overline Q\). It has length
at least \(L=2(1-2\eta)a\). All but at most \(m\) of its
arclength lies on \(\partial^*E\cap Q\), because the part
on \(\partial Q\) is charged to the gluing perimeter and
the measure-theoretic and reduced boundaries agree up to a
length-zero set.
[Ambrosio, Caselles, Masnou and Morel, Corollary 1](https://math.univ-lyon1.fr/homes-www/masnou/fichiers/publications/connected.pdf)

Parametrize a length \(L\) portion of the path by arclength \(s\).
Its chords satisfy \(|x(s)-x(t)|\le|s-t|\). By (HC3),
\[
 j_r(x(s),x(t))\ \ge\
 \frac{\exp(-(|p|+a)^2/(1+r))}{(2\pi)^2\sigma^2}
                     \exp\!\left(-\frac{|s-t|^2}{2\sigma^2}\right).
                                                               \tag{HC6}
\]
The integral of the last kernel on \([0,L]^2\), divided by
\(\sigma\), tends to \(\sqrt{2\pi}L\). Discarding the pairs with
either coordinate in the missing path length reduces this
divided integral by at most \(2\sqrt{2\pi}m\), uniformly in
the positions of the missing portions. If \(E_n\) satisfy
(HC5) in a fixed ball and \(\mu_{E_n}=\mathcal H^1|_{\partial^*E_n}\),
we conclude
\[
 \liminf_n \sigma_n\iint_{B(p,a)^2}
       j_{r_n}(x,y)\,d\mu_{E_n}(x)d\mu_{E_n}(y)
 \ \ge\ \frac{\phi_2(p)}{\sqrt{2\pi}}\,
           e^{-(2|p|a+a^2)/2}\,2a(1-14\eta).
                                                               \tag{HC7}
\]
This uses \(L-2m\ge2a(1-14\eta)\). The good rectangles
for different \(n\) can be chosen independently.

## Localizing the label tensor

In a fixed ball \(B\) with target labels \(i,j\), write \(k\)
for the third label and \(h_{ab}=T_{\sqrt r}(\nu_{ab}|_B)\)
for the smoothed tagged measures, regardless of index order.
Choose scalar signs in the ball orienting the tag vectors as
\(d_{ij},d_{ik},d_{kj}\). The squared norm is
\[
 \begin{aligned}
 \|d_{ij}h_{ij}+d_{ik}h_{ik}+d_{kj}h_{kj}\|_2^2
   &=\|h_{ij}+h_{ik}\|_2^2+\|h_{ij}+h_{kj}\|_2^2\\
   &\hspace{15mm}+\|h_{ik}-h_{kj}\|_2^2 .
 \end{aligned}                                                   \tag{HC8}
\]
The first two terms are the scalar boundary energies of cells
\(i,j\) in the ball. For finitely many disjoint balls, choose
these signs locally and multiply them in each ball by an
independent uniform sign. Use zero directions outside the
balls. Averaging the squared norm cancels the cross terms
between balls. Thus \(F_{D,r}^2\) bounds the sum of the first
two terms of (HC8) over all the balls, already at rank one.

Suppose \(A^{(n)}\) converges in local label measure to a
locally finite-perimeter partition \(A\). At
\(P_\gamma(A;\cdot)\)-almost every interface point \(p\), two
cells have complementary half-plane blowups and
\(P_\gamma(A;B(p,a))=(1+o(1))2a\phi_2(p)\) as \(a\downarrow0\).
Choose balls there so that each of these cells has area
mismatch at most \(\eta^2a^2/2\) with its half-plane, the
relative perimeter error is at most \(\eta\), and
\(e^{-(2|p|a+a^2)/2}\ge1-\eta\). For each fixed ball, local
convergence gives (HC5) for the two approximating cells.

The balls can also be chosen with perimeter-null boundary
and \(P_\gamma(A;B(p,5a))\le10P_\gamma(A;B(p,a))\) by the
density formula at the centers. Fix a bounded open region
and restrict the fine family to balls with closure there.
The \(5r\) covering lemma gives a disjoint countable
selection whose fivefold expansions cover the eligible
centers. The doubling bound gives a finite subfamily
containing at least a fixed fraction of the perimeter
measure in the region. Repeat with balls whose closures
are in the open complement of the already selected closed
balls within the same region. The boundaries are null
and the region has finite perimeter measure, so finitely
many rounds leave arbitrarily small uncovered measure.

For each finite family, use (HC8) and then (HC7) on both
cells in every ball. This bounds the liminf in (HC2) by
\[
 \sqrt{\frac2\pi}\,
 \frac{(1-\eta)(1-14\eta)}{1+\eta}\sum_B P_\gamma(A;B)
\]
from below. Approximate the perimeter on bounded regions
by such families, send \(\eta\) to zero, and exhaust the plane.
This proves (HC2) also when its right side is infinite.
Combined with (HC1), it gives finite Gaussian perimeter of
every limit.

## Simplex recovery and minimizing sequences

The balanced simplex \(T\) has three rays meeting at
\(120^\circ\). On a ray \(\mathbb R_+v_e\), change variables
\(u=t+\sigma v\) in (HC3) to obtain
\[
 \begin{aligned}
 \sigma\int_0^\infty\!\int_0^\infty j_r(tv_e,uv_e)\,du\,dt
 &\longrightarrow (2\pi)^{-2}\int_0^\infty e^{-t^2/2}\,dt
                                 \int_{\mathbb R}e^{-v^2/2}\,dv
 =\frac{1}{4\pi}.
 \end{aligned}
\]
The substituted integrand is dominated by a constant times
\(e^{-t^2/4-v^2/4}\) for \(r\ge1/2\). For distinct rays,
\(|tv_e-uv_f|^2=t^2+u^2+tu\ge t^2+u^2\), and their pair
energy is bounded before multiplication by \(\sigma\).
For any directions, the same-tag coefficient is at most
two in absolute value, and the cross-tag coefficient is
at most one. This yields an upper bound matching (HC2),
uniformly over the direction ranks:
\[
 \lim_{r\uparrow1}\sqrt{1-r^2}\,F_{D,r}(T)^2
          =\frac{3}{2\pi}=\sqrt{\frac2\pi}\,P_\gamma(T).
                                                               \tag{HC9}
\]

The Gaussian multi-bubble theorem gives \(P_\gamma(A)\ge
P_\gamma(T)\) for balanced partitions, with equality only
on the simplex orbit. Its dimension hypothesis holds for
three cells in the plane. Using compactness, (HC2), and
(HC9) on almost-minimizing sequences gives
\[
 \lim_{r\uparrow1}\sqrt{1-r^2}\inf_{A\ {\rm balanced}}
                          F_{D,r}(A)^2=\frac{3}{2\pi}.
                                                               \tag{HC10}
\]
Every balanced sequence with scaled energy at most
\(3/(2\pi)+o(1)\) has Gaussian label distance to the
simplex orbit tending to zero. Otherwise a subsequence at
positive distance would have a limit by (HC1), and (HC2)
would make that limit a perimeter minimizer. This applies
to adaptive minimizers and to planar actual noise-stability
maximizers along right-contact sequences approaching \(r=1\).
At such a contact the maximizer's unrestricted adaptive
energy equals its correlation derivative, which is at
most the simplex derivative.
[Milman and Neeman, Theorems 1.1–1.2](https://arxiv.org/html/1805.10961v3)
[Normal saturation at maximizers and at the simplex](THREE_CELL_ADAPTIVE_SURFACE_UNRESTRICTED_DIRECTIONS.md#the-unique-scalar-diagonal-and-quantitative-weight-rigidity)

Exact optimality for the adaptive minimum near \(r=1\)
still requires a strict comparison in shrinking neighborhoods
of the simplex orbit.
