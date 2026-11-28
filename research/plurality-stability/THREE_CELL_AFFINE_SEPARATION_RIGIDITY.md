# A separating line and connectedness of three-cell maximizers

Let \(A\) be an actual balanced planar three-cell Gaussian
noise-stability maximizer at \(0<r<1\). Suppose that an affine line
separates two of its labels, up to null sets. In orthonormal
coordinates write the line as \(x=c\), with \(A_1\subseteq\{x>c\}\)
and \(A_2\subseteq\{x<c\}\), and let \(\sigma=\sqrt{1-r^2}\).

After possibly reversing the other coordinate, the third cell is
an epigraph with an extended-real height. Every component of every
strict-winner set is unbounded, and there is at most one genuine
triple junction. For \(r\ge 2/3\), the separating line passes
through the Gaussian origin, every height is finite, and all
three strict-winner sets are connected. Semialgebraic jump
networks in this class therefore give the regular centered simplex.

## Rearrangement of the perpendicular fibers

Let \(u_x=1_{A_3}(x,\cdot)\). On the fiber at \(x\ne c\), the
other label is \(\ell(x)=1\) or \(2\), according to the side of
the line. Conditional on two horizontal observations \(x,s\),
the part of conditional agreement that depends on the shapes
of the vertical fibers is
\[
 (1+1_{\{\ell(x)=\ell(s)\}})
             \langle u_x,T_r^{[1]}u_s\rangle_{\gamma_1}.
                                                        \tag{ASR1}
\]
The remaining term depends only on the two fiber masses. Replace
the fibers of \(A_3\) by upper half-lines of the same Gaussian
masses. The masses of both other labels remain fixed. The
two-set Gaussian comparison maximizes every term in (ASR1)
simultaneously. Its equality case says that two nontrivial
fibers achieving equality are same-direction half-lines.
[Mossel–Neeman, Theorems 1.1 and 1.3](https://ems.press/content/serial-article-files/32017)

The two horizontal observations have positive joint density and
the coefficient in (ASR1) is positive. Thus maximality forces
equality for almost every pair of nontrivial fibers. There is a
positive measure of such fibers, since otherwise the whole
partition depends on \(x\) alone, contrary to full active rank.
By Fubini, the nontrivial fibers are parallel half-lines with a
common orientation up to a null set of fibers. Reversing \(y\)
if needed gives
\[
 A_3=\{y>h(x)\},\qquad
 A_1=\{x>c,\ y<h(x)\},\qquad
 A_2=\{x<c,\ y<h(x)\}                                    \tag{ASR2}
\]
up to null sets. Here \(h\) may take the values \(\pm\infty\).
[Full active rank](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible)

Write \(V_i=T_r1_{A_i}-\beta_i\) for the adjusted stationary
scores. For each \(i=1,2\), the difference \(V_3-V_i\) has
strictly positive vertical derivative throughout the plane.
The source function \(1_{A_3}-1_{A_i}\) is nondecreasing in
the vertical coordinate on every fiber in (ASR2) and is
nonconstant on the positive-measure collection of nontrivial
fibers. Its vertical derivative as a distribution is a
nonnegative nonzero measure. Convolution with the positive
Gaussian gives the strict derivative.

On either open side of the line, take \(h(x)\) to be the
unique zero of \(V_3-V_{\ell(x)}\) in \(y\) when it exists.
Define \(h(x)=+\infty\) when the difference is always negative,
and \(h(x)=-\infty\) when it is always positive. This agrees
with (ASR2), since the other non-3 label is absent from that
side and score ties form a null set. The finite-height set on each
side is open, and its graph is analytic and lies in the closure
of the actual jump interfaces. The derivative in \(y\) is
nonzero, and each of the two labels wins in an open region
arbitrarily near the graph on its appropriate side.

The component of any strict winner of label 3 contains its
entire upward ray, since both pair contrasts increase
strictly. For labels 1 and 2 the open part of the corresponding
side where \(V_i>V_3\) is occupied by \(i\) almost everywhere.
The equalities with the absent label inside that open part
are isolated by the dormant-tie Hessian bound (FC4), applied
in a ball contained in the open part. A downward segment can
therefore be followed within the strict winner by detouring
around finitely many such points in a compact neighborhood
of the segment. Taking endpoints arbitrarily far down proves
unboundedness of every strict-winner component.
[Dormant ties and (FC4)](THREE_CELL_FIRST_CONTACT_VARIATIONS.md#point-transfers-at-a-dormant-interior-score-tie)

## Scores and height distributions along the separating line

The positive \(1/2\) interface supplied by complete adjacency
lies on \(x=c\). Hence \(V_1(c,y)=V_2(c,y)\) for every real
\(y\), by analyticity along the line. The common difference
from \(V_3\) is strictly decreasing in \(y\), so it has at most
one zero. Every genuine triple junction must occur on the line.
Wherever \(V_1(c,y)>V_3(c,y)\), the \(1/2\) interface is
active. In a neighborhood label 3 loses, and labels 1 and 2
occupy their respective sides almost everywhere. Wherever
\(V_1(c,y)<V_3(c,y)\) there is no jump along the line. The
\(1/2\) wall is therefore either the whole line or a downward ray.
[Complete adjacency and analytic interfaces](THREE_CELL_GLOBAL_OBSTRUCTIONS.md)

Let \(p_\sigma\) denote the centered Gaussian density of
standard deviation \(\sigma\). The constant value of
\(T_r(1_{A_1}-1_{A_2})(c,y)\) and Gaussian injectivity in
the vertical variable show that
\[
 \int_{s>c}1_{\{t<h(s)\}}p_\sigma(s-rc)\,ds
  -\int_{s<c}1_{\{t<h(s)\}}p_\sigma(s-rc)\,ds
                    =\beta_1-\beta_2
 \quad\text{for almost every }t.                           \tag{ASR3}
\]
The distributional derivative in \(t\) shows that the
pushforwards of \(p_\sigma(s-rc)\,ds\) restricted to the two
sides and finite values of \(h\) agree on \(\mathbb R\).
In particular both sides have a nonempty open finite-height set.

## Limiting vertical profiles and transient widths

Assume \(r\ge2/3\). We show that every height is finite without a
global assumption on the interface network. Fix
\(\tau\in\{-1,1\}\), use \(y=\tau t\), and let \(t\to+\infty\).
Write \(F_i(t,x)=V_3(x,\tau t)-V_i(x,\tau t)\) for \(i=1,2\).
By (ASR2) and dominated convergence its limit \(F_i^\infty(x)\)
is the one-dimensional Mehler smoothing of the limiting
slice contrast, minus the same pair price. Uniform bounds on the horizontal
score derivatives give locally uniform convergence. The limiting
slice gives label 3 where \(h<+\infty\) in the upper direction,
and where \(h=-\infty\) in the lower direction. The other
points have their side label. The function \(F_i^\infty\)
is entire in \(x\) and is not identically zero. If it were
zero, Gaussian injectivity would make its source contrast
constant. On side \(i\) the source contrast has values in
\(\{-1,1\}\), so that constant would be a price difference
of absolute value one, which is impossible.

Let \(q_i(t,z)\) indicate the points of side \(i\) where
\(h(z)\) is finite and \(t<\tau h(z)\). Thus \(q_i(t,z)\)
records a discrepancy between the slice at \(y=\tau t\)
and its limiting label. For each \(z\), it is nonincreasing
in \(t\) and tends to zero. Define
\(D_i(t,x)=\tau(F_i^\infty(x)-F_i(t,x))\).
The vertical conditional law \(S\sim N(rt,\sigma^2)\)
gives the exact formula
\[
 D_i(t,x)=\sum_{\ell=1}^2 (1+1_{\{i=\ell\}})
      \int_{\mathbb R}p_\sigma(z-rx)\,
            \mathbb E[q_\ell(S,z)]\,dz .                    \tag{ASR4}
\]
In both directions the source difference of the two slices
is \(\tau(e_3-e_\ell)\) on the defect of side \(\ell\),
which gives the positive coefficients. The formula is
strictly positive for every finite \(t,x\), since the set
of finite heights has positive measure. On side \(i\), away
from the null set of zeros of \(F_i^\infty\), the defect is
exactly the set where \(0<\tau F_i^\infty(x)<D_i(t,x)\).
This follows from the strictly increasing vertical score
and the definition of its height as the unique zero.
We call \(\tau F_i^\infty(x)>0\) the relevant sign for the
chosen direction.

Suppose \(F_i^\infty(b)=0\) and the relevant sign occurs
in a one-sided interval of side \(i\) abutting \(b\).
Then \((F_i^\infty)'(b)\ne0\). To see this, let \(I\)
be a sufficiently short such interval and \(m(t)\) its
Lebesgue measure of defect. By strict positivity in (ASR4)
and continuity near \(b\), \(m(t)>0\) for every finite \(t\).
It decreases to zero by dominated convergence. In (ASR4)
keep only the contribution of \(I\) itself. On \(I\),
\(D_i(t,x)\ge\kappa\,\mathbb E m(S)\) with a fixed
\(\kappa>0\) that can be kept positive as \(I\) shrinks, by
the positive diagonal density \(p_\sigma((1-r)b)\).
If the derivative at \(b\) were zero, shrink \(I\) so that
\(|F_i^\infty(x)|\le L|x-b|\) on \(I\) with \(L<\kappa/2\).
The defect then contains an initial portion of \(I\) of
length \(\kappa\mathbb E m(S)/L\) for sufficiently large
\(t\). Indeed this length tends to zero and the defining
inequality for the defect holds there except possibly at
its endpoints. But
\(\mathbb E m(S)\ge\Pr(S\le t)m(t)\) and \(\Pr(S\le t)\to1\).
This contradicts \(m(t)>0\), proving the derivative claim.

## The one-dimensional switch matrix

Index the ordinary switches of the limiting slice inside either
open side of \(c\) by \(j\), with offsets \(b_j\) and incidence
\(d_j\) equal to the left label vector minus the right label
vector. These switches are simple zeros by the preceding
claim. They have at most quadratic counting growth in
\([-R,R]\). In fact \(F_i^\infty\) is a Gaussian convolution
of a bounded source, with complex growth bounded by
\(C\exp(C|\operatorname{Im} x|^2)\). Jensen's formula in a disk
of radius \(2(R+|x_0|+1)\) centered at a real nonzero point
\(x_0\) for \(F_i^\infty\) bounds the number of real zeros
in \([-R,R]\) by \(O(1+R^2)\).
[Gaussian complex growth and Jensen counting](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#polynomial-growth-of-the-network)

At an ordinary switch, choose the one-sided interval with
the relevant sign. If \(F_i^\infty(c)=0\) and a side has
the relevant sign next to \(c\), add a switch with offset
\(b_j=c\) and choose the one-sided interval within that side.
Both limiting contrasts have the same value at \(c\) by
(ASR3) and the adjusted pair equality there. Add the fixed
axis incidence \(d_0=e_2-e_1\) at \(c\) exactly when the
\(1/2\) interface is present for all large finite \(t\).
The following rules ensure that the incidences at \(c\)
sum to the true limiting slice incidence there.

- If the limiting contrast at \(c\) is nonzero, both sides
  agree nearby with their finite-\(t\) winners for large
  \(t\). There is no added moving switch at \(c\); the fixed
  axis is present exactly when the contrast is negative.
- If it is zero in the upper direction, the finite score
  difference at \(c\) is negative because (ASR4) is
  positive. The fixed axis is present. A relevant left
  interval has incidence \(e_3-e_2\), and a relevant right
  interval has incidence \(e_1-e_3\). Adding these to
  \(d_0\) gives the transition from the actual limiting
  left label to the actual limiting right label, including
  the case in which the sum is zero.
- If it is zero in the lower direction, the finite score
  difference at \(c\) is positive and the fixed axis is
  absent. A relevant left interval has incidence \(e_2-e_3\),
  and a relevant right interval has incidence \(e_3-e_1\).
  Their sum again gives the limiting incidence, including
  the possibility of a limiting \(2/1\) switch.

Let \(J\) contain all switches other than the fixed axis.
For each \(j\in J\), let \(i(j)\) be its non-3 label.
There is a sign \(\eta_j\in\{-1,1\}\) with
\(\eta_j d_j=e_{i(j)}-e_3\). The sign of
\(d_j\cdot(V^\infty)'(b_j)\) is negative, where
\(V^\infty\) is the vector of limiting adjusted scores.
At an interior switch this follows from its orientation
and simplicity. At \(c\) the relevant side has strict
sign \(\tau\) for \(F_i^\infty\) and a simple zero at
the endpoint, so the listed incidences give the same
negative derivative in each case.

Write \(w_j=p_1(b_j)\), where \(p_1\) is the standard normal
density, and put
\[
 K_{j\ell}=(d_j\cdot d_\ell)p_\sigma(b_\ell-rb_j),\quad
 D_j^{\rm row}=\sum_{\ell\ \mathrm{all}}K_{j\ell},\quad
 \widehat K_{j\ell}=|K_{j\ell}|\quad(j,\ell\in J).
                                                        \tag{ASR5}
\]
Here \(\mathrm{all}\) includes the fixed axis if present.
The distributional derivative of the limiting slice label
vector has atoms \(-d_\ell\) at \(b_\ell\), with the stated
decomposition at \(c\). Gaussian convolution therefore
gives \(-rD_j^{\rm row}=d_j\cdot(V^\infty)'(b_j)\), so
\(a_j=rD_j^{\rm row}=|(F_{i(j)}^\infty)'(b_j)|>0\).
The kernel on \(J\) is explicitly
\[
 \widehat K_{j\ell}
      =(1+1_{\{i(j)=i(\ell)\}})p_\sigma(b_\ell-rb_j).
\]
All these sums converge absolutely, by polynomial
counting growth and Gaussian decay.

For every \(j\in J\),
\[
                 D_j^{\rm row}
                       \le \sum_{\ell\in J}\widehat K_{j\ell}.
                                                        \tag{ASR6}
\]
If the fixed axis is absent or its term is nonpositive,
this follows by taking absolute values. Otherwise consider
a right-side switch, reflecting \(x\) and interchanging
labels for a left-side one. The positive pairing with
\(d_0=e_2-e_1\) occurs at a switch from 3 to 1. There is a
switch \(\ell\) from 1 to 3 between \(c\) and switch \(j\).
If the limiting third cell reaches \(c\) on that side,
this earlier switch is the added upper-direction switch
at \(c\). In the lower direction the fixed axis can be
present only when the limiting contrast at \(c\) is
negative, so the earlier switch is interior. In either
case \(b_\ell\in[c,b_j]\) and \(d_\ell=-d_j\).
The surplus \(|K_{j\ell}|-K_{j\ell}=4p_\sigma(b_\ell-rb_j)\)
dominates the axis contribution \(p_\sigma(c-rb_j)\).
Indeed the maximum squared-distance excess over
\([c,b_j]\) is at an endpoint, and
\[
 (b_j-rb_j)^2-(c-rb_j)^2
   =(1-2r)b_j^2+2rcb_j-c^2
   \le \frac{c^2(1-r)^2}{2r-1}.
\]
Balance gives \(|c|\le1\), since each side contains a
cell of mass \(1/3\). The ratio of the two densities is
therefore at least
\(\exp[-c^2(1-r)/(2(1+r)(2r-1))]\ge e^{-3/10}>1/4\).
This proves (ASR6).

## The finite matrix contradiction

Take a finite subset \(J_0\) of \(J\) and small pairwise
disjoint one-sided intervals \(I_j\) of relevant sign at
\(b_j\), \(j\in J_0\). Coincident offsets at \(c\) use
opposite intervals. Let \(m_j(t)\) be the Lebesgue measure
of defect on \(I_j\). Each \(m_j(t)>0\) and decreases to
zero. For \(0<\delta<1\), the intervals can be made short
enough that
\[
 |F_{i(j)}^\infty(x)|\le(1+\delta)a_j|x-b_j|
 \quad(x\in I_j)
\]
and, for \(x\in I_j,\ z\in I_\ell\),
\(p_\sigma(z-rx)\ge(1-\delta)p_\sigma(b_\ell-rb_j)\).
Keeping only the selected intervals in (ASR4) then gives
\[
 D_{i(j)}(t,x)\ge
        (1-\delta)\sum_{\ell\in J_0}\widehat K_{j\ell}
                                               \mathbb E m_\ell(S).
\]
The right side tends to zero. For large \(t\), the defect
therefore includes the initial interval of length equal to
the right side divided by \((1+\delta)a_j\), by the slope
bound. Also \(\mathbb E m_\ell(S)\ge(1-\delta)m_\ell(t)\)
for all sufficiently large \(t\). Consequently
\[
 a_j m_j(t)\ge
       \frac{(1-\delta)^2}{1+\delta}
       \sum_{\ell\in J_0}\widehat K_{j\ell}m_\ell(t)
 \qquad(j\in J_0).                                        \tag{ASR7}
\]

Suppose \(J\) is nonempty. Gaussian detailed balance makes
\(w_j\widehat K_{j\ell}=w_\ell\widehat K_{\ell j}\).
The double sum of these positive terms over \(J^2\) is finite.
The factor \(p_1(b_j)p_\sigma(b_\ell-rb_j)\) is the joint
normal density, bounded by a constant times
\(\exp[-(b_j^2+b_\ell^2)/(2(1+r))]\), and the counts grow
polynomially. By (ASR6), \(\sum_J w_jD_j^{\rm row}\) is
positive and at most \(\sum_{J^2}w_j\widehat K_{j\ell}\).
For a finite exhaustion of \(J\), the ratio of these two
sums restricted to \(J_0\) and \(J_0^2\), respectively,
therefore gives some finite \(J_0\) with
\[
 \frac{\sum_{j,\ell\in J_0}w_j\widehat K_{j\ell}}
      {\sum_{j\in J_0}w_j D_j^{\rm row}}>r .                \tag{ASR8}
\]
This is the Rayleigh quotient of the constant vector for
the finite matrix \(\operatorname{diag}(D^{\rm row})^{-1}
\widehat K|_{J_0^2}\) in the positive inner product with
weights \(w_jD_j^{\rm row}\). Its largest eigenvalue exceeds
\(r\). But (ASR7) supplies a strictly positive supersolution
and bounds that eigenvalue from above by
\(r(1+\delta)/(1-\delta)^2\). For example conjugating the
nonnegative matrix by the diagonal matrix of \(m_j(t)\)
bounds its spectral radius by its maximum row sum.
Taking \(\delta\) sufficiently small contradicts (ASR8).
Hence \(J\) is empty in both directions.

## Finite heights and centered separation

Each open side has some finite heights by (ASR3). Its
limiting contrast is therefore positive somewhere in the
upper direction and negative somewhere in the lower one.
Since \(J\) is empty, a change of sign on that side is
impossible. A zero touching the relevant sign from both
sides would also have contradicted the derivative argument
following (ASR4). Thus the limiting contrasts are strictly
positive in the upper direction and strictly negative in
the lower direction on both open sides. They have the
same strict signs at \(c\) by continuity and the absence
of a switch with offset \(c\). The unique score root against
label 3 is therefore finite everywhere, and \(h\) is
continuous across the separating line.

On \(x=c\), the limit \(y\to+\infty\) in the adjusted pair
equality makes both unadjusted scores tend to zero. Thus
\(\beta_1=\beta_2\). Sending \(y\to-\infty\) gives equal
masses of the two half-lines under \(N(rc,\sigma^2)\), so
\(c=0\). For \(x>0\), the horizontal conditional density
relative to that at \(x=0\) is
\(\exp[-r^2x^2/(2\sigma^2)+rxs/\sigma^2]\).
The factor \(\exp(rxs/\sigma^2)\) is strictly greater
than one on \(A_1\) and strictly less than one on \(A_2\),
and the two scores at \(x=0\) agree. Hence \(V_1(x,y)>V_2(x,y)\).
Their order reverses for \(x<0\). The strict winners are
therefore exactly the three open regions in (ASR2), with finite
continuous \(h\) and \(c=0\), and all are connected.
At the junction \((0,h(0))\), the difference \(V_1-V_2\)
has zero vertical derivative and a strictly positive
horizontal derivative. The latter is the conditional
Gaussian integral of a positive constant times
\(s(1_{A_1}-1_{A_2})\) at horizontal observation zero.
Combined with \(\partial_y(V_3-V_1)>0\), this makes the
unique triple junction nondegenerate.

If the jump network is semialgebraic, the
[connected semialgebraic rigidity theorem](THREE_CELL_CIRCULAR_INTERFACE_OBSTRUCTION.md#connected-semialgebraic-maximizers)
identifies the maximizer as the regular centered simplex.
