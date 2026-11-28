# The low-correlation slope of the optimal adaptive Gaussian surface energy

For balanced planar three-cell partitions, the squared adaptive
surface minimum and the simplex correlation derivative have the
same value and first derivative at \(r=0\) from the right.
The first-order statement requires control of the tagged
surface moments. Gaussian perimeter convergence alone does
not control these unbounded moments. Cutting off the tail and
repairing the three masses gives a perimeter tail estimate,
which combines with one Mehler test to give first-moment
convergence.
[The adaptive surface variational problem](THREE_CELL_ADAPTIVE_SURFACE_VARIATIONAL_PROBLEM.md)
[The unrestricted direction problem](THREE_CELL_ADAPTIVE_SURFACE_UNRESTRICTED_DIRECTIONS.md)

Let \(F_{D,r}\) denote the adaptive surface functional with
\(D\) direction coordinates, and write \(D=\infty\) for the
unrestricted value. In this note \(2\le D\le\infty\). Let \(T\)
be the balanced planar simplex and set \(P_*=P_\gamma(T)\),
counting each interface once. Then
\[
 \inf_{A\ {\rm balanced}} F_{D,r}(A)^2
                =P_*^2+\frac{9}{8\pi^2}r+o(r)
                       \qquad(r\downarrow0).
                                                               \tag{LC1}
\]
The asymptotic is uniform over these ranks. More generally,
let \(A_r\) be any sequence of balanced ordinary partitions
with \(F_{D,r}(A_r)^2\le F_{\infty,r}(T)^2+o(r)\), allowing
the rank to vary. If \(p_{e,r}\) is the Gaussian perimeter of
tag \(e\), then
\[
 P_\gamma(A_r)-P_*=o(r),\qquad
    \sum_e\bigl(p_{e,r}-P_\gamma(A_r)/3\bigr)^2=o(r).
                                                               \tag{LC2}
\]
Along every subsequence converging in label measure to an
aligned simplex, the three tagged Gaussian surface first
moments converge to the corresponding simplex moments.

## Near-minimal perimeter and tagged convergence

For cyclic tags \(e=(1,2),(2,3),(3,1)\), select constant unit
directions forming an equilateral triple. Write \(\nu_e\)
for the Gaussian surface measures, \(p_e=\nu_e(\mathbb R^2)\),
\(\nu=\sum_e\nu_e\), \(P=\sum_e p_e\), and
\(h_e=T_{\sqrt r}\nu_e\). The label inner product is
\(-1\) across tags and the direction inner product is
\(-1/2\). Thus the squared equilateral energy is
\[
 \frac12\left\|\sum_e h_e\right\|_2^2+\frac32\sum_e\|h_e\|_2^2
                 \ \ge\ \left\|\sum_e h_e\right\|_2^2 .
                                                               \tag{LC3}
\]
Its constant part is \(P^2+\frac32\sum_e(p_e-P/3)^2\).
In particular \(F_{D,r}\ge P\).

The simplex is normally saturated for the unrestricted
functional, so \(F_{\infty,r}(T)^2=S'_r(T)\). By the simplex
ray formula the expansion is
\[
 F_{\infty,r}(T)^2=P_*^2+\frac{9}{8\pi^2}r+O(r^2).
                                                               \tag{LC4}
\]
Indeed the derivative is \(6I_1(r)+3I_{-1/2}(r)\), where the
ray integral has \(I_t(0)=1/(8\pi)\) and \(I_t'(0)=t/(4\pi^2)\).
The bound \(P_\gamma(A_r)\ge P_*\) and its simplex equality case
follow from the Gaussian multi-bubble theorem.
[Normal saturation at the simplex](THREE_CELL_ADAPTIVE_SURFACE_UNRESTRICTED_DIRECTIONS.md#the-unique-scalar-diagonal-and-quantitative-weight-rigidity)
[The simplex ray formula (D1)–(D2)](THREE_CELL_DEFICIT_AND_PERIMETER.md#exact-correlation-derivative-for-the-three-cell-simplex)
[Milman and Neeman, Theorems 1.1–1.2](https://arxiv.org/html/1805.10961v3)

The hypothesis, (LC3), and (LC4) give \(P-P_*=O(r)\).
Gaussian perimeter compactness and the equality statement
give subsequential label convergence, with rotations or
label permutations to align the limit with \(T\). Fix such
a subsequence. Its tagged measures converge locally weakly
to the corresponding tagged measures of \(T\). To verify
this, extract any local weak limits \(\zeta_{ij}\) of the
tags. Lower semicontinuity of each cell's perimeter as a
measure gives \(\sum_{j\ne i}\zeta_{ij}\ge P_\gamma(T_i;\cdot)\).
The sum of these three limit cell measures has total
mass at most \(2P_*\) by compact exhaustion, and the sum
of the lower bounds has mass exactly \(2P_*\). Each measure
inequality is therefore equality. The three cell measures
determine the three pair measures by inverting the incidence
sums, proving the asserted tagged convergence.

## Uniform first-moment tails

Use two smooth compactly supported ambient vector fields
whose volume derivatives at \(T\) span the centered
three-label mass plane. One can choose bump fields normal
to two distinct simplex interfaces away from the vertex.
They give independent derivatives proportional to the
two corresponding label differences. For their composed
flows, Gaussian masses are integrals on the original
cells of the Gaussian Jacobian. Their first derivatives
are continuous under label convergence and their second
derivatives are bounded uniformly for small parameters.
Since the two-dimensional volume derivative remains
invertible for \(A_r\) with uniformly bounded inverse,
the inverse function theorem gives mass corrections of
parameter size at most \(C\) times the error for all
sufficiently small prescribed errors, uniformly in \(r\).
The Gaussian tangential Jacobian differs from one by at
most \(C\) times the parameter size on the support of the
flows. The perimeter cost is bounded by the same order,
since the original local perimeters are bounded.

For \(R\) larger than a fixed \(R_0\), replace the exterior
of \(B_R\) by one label. All original interfaces outside
the disk are deleted. The extra perimeter is at most
the Gaussian length \(R e^{-R^2/2}\) of its circle, using
the \(BV\) gluing formula at almost every \(R\).
The volume error is \(O(e^{-R^2/2})\). The replacement
does not affect the compact support of the chosen flows,
so they correct the mass with perimeter cost
\(O(e^{-R^2/2})\), uniformly. Comparing the repaired
partition with the perimeter minimum gives
\[
 \nu(\{|x|>R\})\le P-P_*+C(1+R)e^{-R^2/2},
                        \qquad R\ge R_0.
                                                               \tag{LC5}
\]
The bound extends from almost every radius by taking
good radii decreasing to the desired one.

For \(h=T_{\sqrt r}\nu\), \(\int h\,d\gamma_2=P\).
By (LC3) and (LC4), \(\|h-P\|_2=O(\sqrt r)\).
For any unit vector \(v\), positivity and the Mehler
kernel give
\[
 \int \cosh(\sqrt r\,x\cdot v)\,d\nu(x)
      = e^{-(1-r)/2}\int h(z)\cosh(z\cdot v)\,d\gamma_2(z).
\]
The right side equals \(P+O(\sqrt r)\) by Cauchy–Schwarz
and \(\int\cosh(z\cdot v)\,d\gamma_2(z)=e^{1/2}\).
Since \(\cosh u-1\ge u^2/2\), testing both coordinate
directions gives \(\int |x|^2\,d\nu=O(r^{-1/2})\).

Set \(L=r^{-3/4}\). For fixed \(R\ge R_0\), integrating
(LC5) up to \(L\) bounds the first moment on
\(R<|x|\le L\) by \(C rL+C(R+1)^2 e^{-R^2/2}\).
Beyond \(L\), the mass bound \(O(r)\) from (LC5) and the
second-moment bound give first moment \(O(r^{1/4})\)
by Cauchy–Schwarz. Consequently
\[
  \lim_{R\to\infty}\limsup_{r\downarrow0}
                          \int_{|x|>R}|x|\,d\nu(x)=0
                                                               \tag{LC6}
\]
along the aligned subsequence. Since each tagged measure
is at most \(\nu\), (LC6) and the local tagged convergence
show that \(m_e=\int x\,d\nu_e(x)\to(2\pi)^{-1}v_e\)
for the unit ray directions \(v_e\) of \(T\).

## The first-degree contribution

The degree-one Hermite coefficient vector of \(h_e\)
is \(\sqrt r\,m_e\). This follows from Fubini and
the conditional Gaussian mean in the Mehler kernel.
The first-moment estimate justifies absolute
integrability for these tests. Keeping the constant and
degree-one parts of the nonnegative expression in
(LC3) now gives
\[
 F_{D,r}(A_r)^2\ \ge\
 P^2+\frac32\sum_e(p_e-P/3)^2+
       r\left(\frac12\left|\sum_e m_e\right|^2+
                              \frac32\sum_e|m_e|^2\right).
                                                               \tag{LC7}
\]
For the simplex moments \(\sum_e v_e=0\) and \(|v_e|=1\).
The last parenthesis therefore tends to \(9/(8\pi^2)\).
Combine (LC7) with (LC4) and \(P\ge P_*\) to get
(LC2) and the lower bound in (LC1). This argument applies
on every subsequence after extracting a further aligned
one. To obtain the infimum statement, take competitors
within \(o(r)\) of the infimum, which is bounded above
by \(F_{\infty,r}(T)^2\) for every rank in the stated
range.

This determines the first-order value of the auxiliary
minimum. A strict simplex comparison at positive
correlations would require control of the remainder
relative to the distance from the simplex orbit.
