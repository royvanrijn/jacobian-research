# A two-block comparison at the simplex junction

Let \(k\ge3\), \(d=k-1\), \(E=\mathbf1^\perp\subset\mathbb R^k\),
and let \(T\) be the regular simplex partition, with label the largest
coordinate. Split the labels into two nonempty blocks \(C,D\).
In each block use the reference winner within that block. A competitor
can choose either block arbitrarily at each point, then assigns the
winner within the chosen block. For the tetrahedron split as
\(2+2\), this means only that labels 1 and 2 lie on their respective
sides of \(x_1=x_2\), and labels 3 and 4 lie on their respective
sides of \(x_3=x_4\). No reflection symmetry or volume constraint
is imposed.

Each such competitor that differs from \(T\) on a set of finite
Lebesgue measure strictly increases the ordered cut for Euclidean
unit Gaussian convolution, unless it agrees almost everywhere.
There are constants \(c,C,\delta>0\) depending only on \(k\) such
that the Gaussian deficit satisfies
\[
 S_\rho(T)-S_\rho(U)\ge c\,e^{-C(R+1)^2}L_\rho(U;T),
 \qquad 1\le R\le\delta\sqrt{\log(1/e)},\qquad e=\sqrt{2(1-\rho)},
 \tag{BC1}
\]
for \(\rho\) sufficiently close to one and changes supported in the
physical active ball \(B_{eR}\). Here \(L_\rho\) is the exact oriented
OU transfer cost. On each fixed \(R\) one also has
\(S_\rho(T)-S_\rho(U)\ge c_R |\{U\ne T\}|^2/e^d\), with physical
Lebesgue measure. The local-step argument (BC10) below extends this comparison to
a fixed physical ball with a polynomial relative margin.

The argument uses the sign calibration from
[Pagliari, Theorem 2.3 and Lemma 2.4](https://arxiv.org/pdf/1905.00623).
That theorem applies to scalar nonlocal perimeter with prescribed
exterior data and zero calibration divergence. For the folded kernel here, we verify cancellation of the row
potential directly. The [Lawlor–Morgan tetrahedral calibration,
Section 1](https://msp.org/pjm/1994/166-1/pjm-v166-n1-p04-p.pdf)
concerns local interfacial area.

## Sorting the finite reflection orbit

Use Lebesgue measure and a symmetric pair kernel \(K(p,q)\) on
\(E^2\) with bounded marginals, invariant under simultaneous
orthogonal transformations. Write \(m\) for the indicator of
block \(C\) in \(U\) and \(\eta\) for that in \(T\).
Let \(\chi_C(p,q)\) indicate agreement of the two internal
\(C\)-winners, and similarly for \(\chi_D\). With
\[
 H(p,q)=K(p,q)(\chi_C(p,q)+\chi_D(p,q)),\qquad
 u(p)=\int K(p,q)(\chi_C(p,q)-\chi_D(p,q))dq,
\]
the ordered cut relative to \(T\) is
\[
 \mathcal J_K(m)
 =\tfrac12\iint H(p,q)(|m(p)-m(q)|-|\eta(p)-\eta(q)|)dp\,dq
   -\int u(p)(m(p)-\eta(p))dp .
 \tag{BC2}
\]
Indeed agreement at \((p,q)\) is
\(\chi_C m(p)m(q)+\chi_D(1-m(p))(1-m(q))\).
For binary values, \(2m(p)m(q)=m(p)+m(q)-|m(p)-m(q)|\),
and symmetry yields (BC2).

Let \(W=S_C\times S_D\) act by within-block permutations and
\(N=|W|\). Sort the \(N\) orbit copies of the block indicator into
\[
 m_j(p)=1_{\{\sum_{w\in W}m(wp)\ge j\}},\qquad 1\le j\le N .
\]
These indicators are individually invariant under \(W\).
The sum over \(j\) of \(|m_j(p)-m_j(q)|\) is exactly the difference
in absolute value between the two orbit sums. It is at most
\(\sum_w|m(wp)-m(wq)|\). The linear terms add exactly.
By invariance this gives \(\mathcal J_K(m)\ge
N^{-1}\sum_j\mathcal J_K(m_j)\). Also
\[
 \sum_j|m_j(p)-\eta(p)|=\sum_w |m(wp)-\eta(wp)|.
 \tag{BC3}
\]
Thus the mean defect measure is preserved. The mean transfer cost
is preserved as well: its integrand is \(|m-\eta|\) times the
absolute difference between the two reference internal-winner scores,
a \(W\)-invariant weight. Reference score ordering follows by
reflecting across the hyperplane exchanging the two labels, so this
weighted defect is precisely the oriented score regret. Ball support
is preserved, and all relative integrals converge absolutely by
the bounded marginal and finite defect hypotheses.

## Folding and the translated score tie

Take the fundamental chamber \(Q\) with one fixed coordinate
ordering in each block, and let \(i\in C,j\in D\) be the winners
there. For an invariant profile on \(Q\) define
\[
 A(p,q)=\sum_{w:w i=i}K(p,wq),\qquad
 B(p,q)=\sum_{w:w j=j}K(p,wq),\qquad
 U(p)=\int_Q(A-B)(p,q)\,dq .
\]
Both kernels are symmetric by inversion in their respective
stabilizers and are at least \(K(p,q)\). Folding (BC2) over the
first variable into \(Q\) gives
\[
 \frac{\mathcal J_K(m)}{N}
  =\frac12\iint_{Q^2}(A+B)
     (|m(p)-m(q)|-|\eta(p)-\eta(q)|)\,dp\,dq
      -\int_Q U(p)(m(p)-\eta(p))\,dp .
 \tag{BC4}
\]
Let \(g(p)=p_i-p_j\), so \(\eta=1_{\{g>0\}}\) on \(Q\), and let
\(v\in E\) be constant within each block with \(g(v)=1\).
Translation by \(c v\) preserves \(Q\) and turns the block
indicator of \(T\) there into \(1_{\{g>c\}}\).
Set \(s(p,q)=\operatorname{sgn}(g(p)-g(q))\). At \(p\in Q\)
put \(c=g(p)\). The row identity
\[
 \int_Q(A+B)(p,q)s(p,q)\,dq-U(p)
 =2\left[
   \int_{\{q\in Q:g(q)<c\}} B(p,q)\,dq
    -\int_{\{q\in Q:g(q)>c\}} A(p,q)\,dq
   \right]
 \tag{BC5}
\]
expresses the right side as twice the score of label \(j\) minus
that of \(i\) for \(T+c v\) at \(p\).
For \(K(p,q)=\phi_1(q-p)\) these scores tie. Indeed \(p-c v\)
lies in the pair hyperplane \(\{g=0\}\), and transposition of
\(i,j\) reflects the two reference cells while fixing this point.

Since \(|\eta(p)-\eta(q)|=s(p,q)(\eta(p)-\eta(q))\),
antisymmetry of \(s\) makes (BC4) the sum of the row-error
term, given by the left side of (BC5) multiplied by \(m-\eta\)
and integrated, and the nonnegative integral
\[
 \tfrac12\iint_{Q^2}(A+B)
 \bigl[|m(p)-m(q)|-s(p,q)(m(p)-m(q))\bigr]\,dp\,dq .
 \tag{BC6}
\]
In the Euclidean case the row error vanishes. A wrong point
\(p\) with \(g(p)>0\) has strictly positive slack against all
correct points in \(Q\) with \(0<g(q)<g(p)\).
This region has infinite measure because it contains a strip
along the positive-dimensional internal order chamber, and
only finite defect volume was removed. For \(g(p)<0\) take
\(g(p)<g(q)<0\). Thus a defect of positive measure has positive
slack. Equations (BC2)–(BC3) give the same strict conclusion
for the original unsorted competitor.

## Bounds for scaled Ornstein–Uhlenbeck kernels

Set \(a=\sqrt{(1+\rho)/2}\). For \(X=e p,Y=e q\), their joint
density in the scaled variables equals \(e^d K_e(p,q)\), with
\[
 K_e(p,q)=\phi_a(e(p+q)/2)\phi_1(q-p)
          =\phi_1(e p)\phi_a(q-\rho p),
 \tag{BC7}
\]
where \(\phi_a\) is the centered Gaussian density of covariance
\(a^2 I_d\). The first formula follows by writing
\((X+Y)/2\) and \((Y-X)/e\) as independent Gaussians.
The second is conditional on \(p\).
Thus the actual stability deficit and \(L_\rho\) are \(e^d\)
times the relative cut and score regret computed with \(K_e\).

The two translated scores in (BC5) are now calculated using the
conditional mean \(\rho p-cv\) and variance \(a^2 I_d\), with
row factor \(\phi_1(e p)\). Their difference vanishes whenever
the mean is on \(\{g=0\}\), by reflection. The difference has
a uniform Lipschitz constant as a function of the conditional mean,
from the \(L^1\) bound on the gradient of the Gaussian density
for \(a^2\ge1/2\). Since \(g(\rho p-cv)=(\rho-1)g(p)\), the
row error in (BC5) is bounded on \(B_R\) by
\[
 C e^2(R+1)\min(|g(p)|,1).
 \tag{BC8}
\]
The unshifted reference score difference is bounded by
\(C\min(|g(p)|,1)\) on all \(Q\), by the same reflection and
Lipschitz bound, together with bounded scores. Hence the folded
score regret divided by \(N\) is at most \(C\int_{\{m\ne\eta\}\cap Q}
\min(|g(p)|,1)\,dp\).

To bound (BC6), let \(V=v^\perp\cap E\) and write
\(q=z+(h-g(z))v\) for \(z\in V\) and \(h=g(q)\).
The measure in these coordinates is a constant times \(dz\,dh\).
The internal order chamber in \(V\) has an interior unit vector
\(z_*\), since \(k\ge3\). Take \(z\) in a small fixed-radius
ball about \((R+3)z_*\) inside that chamber and \(|h|<1\).
The resulting points stay in \(Q\) because the \(v\)-direction
preserves each internal ordering. They have \(|q|\ge|z|>R\)
and \(|q|\le C(R+1)\), and therefore have the reference label.
The \(q\)-volume for \(h\) in any interval inside \((-1,1)\)
is at least a positive constant times its length. If \(g(p)>0\)
and \(p\) is wrong, take \(0<h<\min(g(p),1)\). If \(g(p)<0\)
and \(p\) is wrong, take \(\max(g(p),-1)<h<0\).
On each of these pairs the bracket in (BC6) equals 2.
For \(p\in B_R\) and all these \(q\), (BC7) and \(A+B\ge K_e\)
give a lower bound \(c\exp(-C(R+1)^2)\) on the kernel.
Consequently the slack is at least
\[
 c\exp(-C(R+1)^2)\int_{\{m\ne\eta\}\cap Q}\min(|g(p)|,1)\,dp .
\]
For \(1\le R\le\delta\sqrt{\log(1/e)}\), choosing \(\delta\)
small and then \(e\) small permits absorption of (BC8) in half
this lower bound. This proves (BC1) for invariant profiles,
and (BC3) gives it for the original profile.

For \(R\) fixed, the \(g\)-strip of width \(2s\) in \(B_R\)
has volume at most \(C_R s\). For any defect subset of volume
\(M\), taking \(s=\min(1,M/(2C_R))\) shows that
\(\int\min(|g|,1)\) over the defect is at least \(c_R M^2\).
The folded bound therefore controls the square of the unsorted
defect volume by (BC3) and Jensen's inequality. On restoring the
physical scale the bound is \(c_R|\{U\ne T\}|^2/e^d\).


The Gaussian score-relative comparison also holds in ambient
dimension (nge k-1) when the support condition is imposed on the
active projection. Conditional on the two correlated transverse
points, positivity of the active OU operator bounds each mixed
interaction of the active defects by the mean of their two self
interactions. The score cost integrates exactly over active
slices. This is the same active-space lift as for
[the separated exterior comparison](THREE_CELL_EXTERIOR_TUBES.md#changes-away-from-four-label-junctions).

For four cells this verifies the noise-scale core comparison for
the two-pair ordering class. The unrestricted core can exchange
labels within either ordered pair. Those exchanges and the
region between the logarithmically growing noise-scale ball and
the physical localization radius in (CT19) still require control.

## Reflected tetrahedral pairs

For \(k=4\), (BC1) also applies when, for either pair, the
ordering assumption is replaced by exact equivariance under the
reflection exchanging that pair and fixing the other two labels.
The Euclidean finite-volume conclusion has the same extension.

For a reflected pair \(i,j\) with normal \(n\), put
\(H^+=\{n\cdot p>0\}\) with \(i\) its reference-side label, and
\(b=1_{U_i}-1_{U_j}\) on \(H^+\). Then \(b\) extends oddly and
the pair union is even. Replace \(b\) by \(|b|\) on \(H^+\) and
continue it oddly. The pair now has the required ordering.
Expanding the pair interaction into its sum and difference gives
the gain in agreement, in the scaled coordinates,
\[
 \iint_{(H^+)^2} [K(p,q)-K(p,R_nq)]
 \bigl(|b(p)b(q)|-b(p)b(q)\bigr)\,dp\,dq ,
 \tag{BC9}
\]
where \(R_n\) is the reflection. The kernel bracket is positive.
Distinct disjoint-pair reflections commute, and the sign of
\(n\cdot p\) is invariant under the other reflection, so sorting
one pair preserves the hypothesis on the other. The competitor
already agrees with the reference outside its defect support,
which the sorting does not enlarge.

For \(K=K_e\) in (BC7), on a fixed-radius reference anchor
\(q\) in cell \(i\) at distance \(R+O(1)\) outside \(B_R\),
one has \(b(q)=1\) and \(n\cdot q\ge c>0\).
For \(p\in B_R\cap H^+\), the ratio
\(K_e(p,R_nq)/K_e(p,q)=\exp(-2\rho(n\cdot p)(n\cdot q)/a^2)\)
and the kernel lower bound on this anchor give
\[
 K_e(p,q)-K_e(p,R_nq)
 \ge c \exp(-C(R+1)^2)\min(n\cdot p,1).
\]
The decrease in exact score cost due to the sorting is at most
\(C\int_{\{b<0\}\cap H^+}\min(n\cdot p,1)\,dp\) in scaled units.
Indeed it is the sum of the reference-score differences between
the exchanged labels at the two reflected points. Those scores
tie on the pair plane and have uniformly bounded scaled
gradients by (BC7). Thus the gain in (BC9) controls the score
cost removed by the sorting, with the margin in (BC1).
Applying (BC1) to the resulting ordered partition controls
its remaining score cost, giving the comparison for the
original partition. The Euclidean gain is strictly positive
for a wrong-side region of positive measure: \(b(q)=1\) on
a positive-measure portion of the infinite reference cell
on \(H^+\), as the defect has finite measure.

## A local-step comparison up to a fixed physical radius

The two-block domain admits a polynomial comparison. In Gaussian
notation there are \(c,\delta>0\), depending only on \(k\), such
that the ordered two-block comparison holds with
\[
 S_\rho(T)-S_\rho(U)\ \ge\ c R^{-2} L_\rho(U;T),
 \qquad 1\le R\le \delta/e,\qquad \{U\ne T\}\subset B_{eR}.
 \tag{BC10}
\]
No balance condition is needed. In particular a fixed physical
ball of sufficiently small radius is admissible, with a
transfer-relative margin of order \(e^2\) at that radius.
The same conditional comparison holds in extra Gaussian
dimensions under the active support condition. For the reflected
tetrahedral pair class the margin in this range is \(cR^{-3}\).

Use the chamber and linear function \(g\) from (BC4)–(BC6).
Write \(Q\) in coordinates \(p=z+(h-g(z))v\), with
\(z\in Z=Q\cap v^\perp,\ h=g(p)\). Here \(Z\) is a full-dimensional
open convex cone in \(v^\perp\), of dimension \(d-1\ge1\), and the
volume element is a constant multiple of \(dz\,dh\). Set
\(D_\pm=\{m\ne\eta,\ \pm h>0\}\) in \(Q\).
The following argument uses sequences of short random steps to
lower bound (BC6).

Choose a unit vector \(z_*\in Z\) and a small ball \(W_0\Subset Z\)
with \(\inf_{w\in W_0}w\cdot z_*>0\).
Starting with density \(|h|1_{D_\pm}(p)\,dp\), move by independent
steps
\[
 z\longmapsto z+w,\qquad
 h\longmapsto(1-\theta/N)h,\qquad
 w\sim\mathrm{Unif}(W_0),\quad \theta\sim\mathrm{Unif}([1,2]).
\]
Take \(N=\lceil C_0(R+1)\rceil\). The cone property keeps all points
in \(Q\). The final \(z\)-projection along \(z_*\) exceeds \(R\)
if \(C_0\) is large enough, so the last point lies outside the
change ball. All points stay within \(C_k R\) of the origin, since
the initial point lies in \(B_R\) and \(|h|\) decreases.
The steps have length bounded by a fixed small constant, by
choosing \(W_0\) small and \(C_0\) large. The sign of \(h\) is
preserved.

For fixed increments up to any time \(j\le N\), the map in the
coordinates \((z,h)\) is a translation in \(z\) and a contraction
of \(h\) by a factor at least \((1-2/N)^N\ge1/16\) for \(N\ge4\).
Thus the marginal at time \(j\) has density at most \(C|h|\).
The transition density is at most \(CN/|h|\), because the
conditional \(h\)-interval has length \(|h|/N\). Therefore every
pair of successive points has joint density at most \(CN\).
The initial point is wrong and the terminal point is correct, so
there is a transition from wrong to correct. For \(h>0\) such a
transition decreases \(g\) and switches \(m\) from zero to one.
For \(h<0\) it increases \(g\) and switches \(m\) from one to zero.
In either case the bracket in (BC6) equals two. Summing over the
\(N\) possible steps bounds the initial mass by \(CN^2\) times
the measure of these violating pairs. The pairs used here have
both points in \(B_{C_kR}\) and bounded small separation.

For (BC7) the kernel on these pairs is bounded below by
\(c\exp(-C e^2R^2)\): the displacement Gaussian has bounded
argument and only the midpoint Gaussian sees points of size \(R\).
Therefore (BC6) is at least
\[
 c R^{-2}\exp(-C e^2R^2)
          \int_{\{m\ne\eta\}\cap Q}|g(p)|\,dp .
\]
The translated scores in (BC5) have difference in absolute value
at most \(C e^2 |g(p)|\), directly from the scaled conditional mean
\(\rho p-g(p)v\), which lies at distance
\((1-\rho)|g(p)|/\sqrt2\) from the pair plane, and the uniform
Lipschitz score bound there. Choose \(\delta\) small enough that
the row error is at most half the positive term for \(eR\le\delta\).
The result bounds the exact score cost, by the reference
Lipschitz estimate used in (BC8). The scaling factor \(e^d\)
cancels in the comparison. Orbit sorting preserves the average
score cost by (BC3), so the inequality applies to unsorted
competitors. The usual positive-operator argument for transverse
Gaussian slices supplies the ambient lift.

### Weyl-symmetric kernels

The Euclidean statement holds for every integrable, nonnegative,
even kernel \(k(q-p)\) invariant under the site permutation group
and bounded below by a positive constant on a ball about zero.
The translated row identity in (BC5) is exact by translation and
pair reflection. The short-step comparison therefore gives
\(c R^{-2}\int |g|1_{\{m\ne\eta\}}\) after folding, for defects
supported in \(B_R\), \(R\ge1\). For a finite-volume defect with
unbounded support, start the same chain on a bounded subset
of positive mass of one \(D_\pm\). The terminal law leaves \(B_R\),
has uniformly bounded density there and has bounded \(|h|\)
by the initial bound. Its mass of wrong terminal points
consequently tends to zero as \(R\to\infty\). Some finite choice
of \(R\) then forces a positive measure of violating pairs.
This gives strictness without requiring long-range interactions.

For bounded support one can also allow singular kernels satisfying
\(\int \min(1,|z|)\,k(z)\,dz<\infty\), with the same symmetries
and local lower bound. Compare local cuts with pairs restricted
to those having at least one endpoint in \(B_R\). The reference
cut is finite, by the strip bound along the finitely many
reference supporting hyperplanes, and the competitor cut
may take value \(+\infty\). Apply the preceding bound with
\(k_M=\min(k,M)\), then use monotone convergence on both
local cuts. The lower bound for local interacting pairs is
uniform in sufficiently large \(M\). This includes fractional
kernels \(c|z|^{-d-s}\) with \(0<s<1\).

A nontranslation-invariant extension uses
\(K_\epsilon(p,q)=\psi(\epsilon(p+q)/2)k(q-p)\).
Here \(\psi\ge0\) is site-permutation invariant, bounded and
\(C^2\) with bounded second derivative, and bounded below on
a neighborhood of zero. Assume also \(\int |z|k(z)\,dz<\infty\).
At \(c=g(p)\), change variables from the translated comparison
partition \(T+cv\) back to \(T\). Write \(b=p-cv,\ a=p+cv\).
For the pair transposition the normal component of \(b\)
vanishes, so the reflected displacement factor
\(k(z-b)\) agrees exactly. Along the pair normal the midpoint
factor is an even function of that normal coordinate.
The difference between its values at
\(\epsilon(a_N+z_N)/2\) and \(\epsilon(a_N-z_N)/2\)
is at most \(C\epsilon^2|a_Nz_N|\), by the double integral
of the second derivative and its vanishing when either input
is zero. Since \(a_N=\sqrt2 c\) and \(z_N=(z-b)_N\), the first
moment of \(k\) bounds the (BC5) row error by
\(C\epsilon^2|g(p)|\). On the local-step pairs the kernel has a
positive uniform lower bound when \(\epsilon R\le\delta\).
The same absorption therefore gives \(cR^{-2}\int |g|1_{\{m\ne\eta\}}\)
after folding for \(1\le R\le\delta/\epsilon\).

### Reflected tetrahedral pairs at fixed physical radius

The polynomial comparison extends the reflected-pair step (BC9).
Work in the positive halfspace of the pair \(i,j\), where \(i\)
is the correct internal winner. Write \(b=1_{U_i}-1_{U_j}\)
as in (BC9). After both sorts, let \(L_{\rm ord}\) be the
score cost of the ordered partition in scaled units. Sorts
preserve block membership, so the final ordered labels
can be used when charging defects from either sort.
For the current sort, the removed score cost is at most
a constant times the integral of \(\min(n\cdot p,1)\) on
\(b=-1\) in the positive halfspace. This follows from the
scaled Lipschitz score bound and reflection, as in (BC9).

Take a fixed small \(\eta>0\). When \(p_i<p_l-\eta\) for some
label \(l\) outside this pair, the final ordered partition
still assigns the wrong block, and its exact regret at this
point is bounded below by a positive constant in scaled units.
Indeed the winning reference score has gap from \(i\) at least
\(c\min(p_{\rm win}-p_i,1)\) by the saturated transfer estimate.
For completeness that estimate follows by integrating on
\(\rho p+a Z\) with \(Z\) in a fixed compact patch inside the
winning cone, bounded away from the plane of the two sites.
The Gaussian mass of the patch is fixed, the translated patch
stays in the cone, and the difference from the reflected
density is at least a constant times the minimum of the
plane distance of \(p\) and one. The row factor
\(\phi_1(e p)\) is bounded below on the fixed-radius domain.
Charge all such initial points directly to \(L_{\rm ord}\).

For the other initial points \(p\) let the starting measure be
\(\min(n\cdot p,1)1_{\{b=-1\}}dp\), restricted to this
approximate winning region. Choose independent increments
uniform on a fixed small ball in the interior of \(T_i\),
with every coordinate gap \(w_i-w_l\ge2\eta\), shrinking
\(\eta\) if necessary. After the first step the points are
in \(T_i\) with all gaps at least \(\eta\), and with
\(n\cdot q\ge c>0\). The \(n\)-coordinate increases.
After \(N=O(R)\) steps the point lies outside \(B_R\)
in \(T_i\); thus \(b=+1\). Pairs have bounded separation
and both points lie within \(CR\). At each step the
marginal density is bounded by \(C\min(n\cdot p,1)\):
for each sequence of increments the map is a translation
increasing this weight. The conditional step density is bounded.

The first transition away from \(b=-1\) either goes to \(+1\)
or to \(0\). In the former case, charge its joint mass to
the gain (BC9): for this pair the reflected kernel difference
is at least \(c\min(n\cdot p,1)\), since \(n\cdot q\ge c\)
and \(eR\le\delta\). In the latter case the next point lies
in the deep \(i\)-cone but belongs to the other block, so
the final ordered regret there is at least a fixed positive
constant. Charge this mass to \(L_{\rm ord}\) using the
bounded marginal. Summing over \(O(R)\) steps shows that
the score cost removed by either sort is at most
\(CR(\text{agreement gain of that sort}+L_{\rm ord})\).
The preceding ordered comparison has deficit at least
\(cR^{-2}L_{\rm ord}\). Add the nonnegative gains of the
required sorts to obtain the claimed \(cR^{-3}\) bound for
the original reflected-pair class.

### One ordered block by Euclidean fiber rearrangement

Suppose one block is a pair, \(C=\{i,j\}\), and require only
reference ordering within its complement \(D\). The labels
\(i,j\) may occur on either side of their plane. For the
Euclidean statement, the comparison remains strict for this
larger class for every integrable radial nonincreasing kernel
which is positive almost everywhere and has a finite first
moment. In particular the Euclidean tetrahedral comparison
only needs the separation of a single prescribed pair.

Slice in the direction of the normal for \(i,j\), with transverse
coordinate \(s\). All labels within \(D\) have scores whose
relative order is independent of the normal coordinate. Thus
each slice has at most three possible labels: the two pair
labels and the reference winner \(M(s)\) in \(D\). The reference's
two tails are half-lines and its middle interval, of label
\(M(s)\), is centered. Apply (J11)–(J13) of
[the Euclidean three-label comparison](EUCLIDEAN_THREE_CELL_JUNCTION.md#strict-comparison-by-fiber-rearrangement)
to the two tails. Their signed lengths are well-defined using
finite-volume change. For each pair of transverse fibers the
conditional displacement kernel is even in the normal
coordinate, so averaging the two half-line lower bounds gives
the corresponding centered-tail bound. For a pair of fibers
with the same \(D\)-winner, the one-dimensional Riesz
rearrangement of their \(D\) sets increases agreement under the
nonincreasing normal kernel. For fibers whose \(D\)-winners
differ, agreement between their \(D\) labels is zero before and
after the rearrangement, and the sum of these two binary
cuts is just their summed lengths. The reduction therefore
centers the middle interval and fills the two tails without
increasing cut. It gives a competitor of finite-volume change
in the ordered two-block class.

In the equality case the two integrated tail inequalities
must saturate for almost every pair of fibers and almost
every normal displacement, by positivity of the kernel.
Choose a sequence of good transverse fibers escaping along
the open cone where the average site of \(i,j\) strictly
wins against all sites of \(D\), with the pair-tail defects
converging to zero in \(L^1\). Such fibers exist by finite
volume and Fubini. On them the reference middle interval
is empty and both reference tails have threshold zero.
Continuity in normal displacement, followed by the same
translation argument as in (J12), identifies the original
pair labels on almost every fiber as half-lines. Strict
convexity of the even averaged normal cost makes their
midpoints constant across fibers. Finite-volume change
makes that constant zero. Thus the original competitor
already agrees with the centered rearrangement, and the
ordered comparison gives equality only at the reference.

## Binary and cyclic changes without ordering assumptions

The same polynomial comparison applies to arbitrary transfers
between a fixed pair of labels, with all other labels unchanged
from the simplex. It also applies to changes whose directed
label transfers contain at most one arrow with any given source
and at most one arrow with any given destination. There is no
symmetry or orientation assumption on these changes. In the
Euclidean case the assertion holds for the Weyl-symmetric
convolution kernels specified above. For the Gaussian kernel
the support range and \(cR^{-2}L_\rho\) bound are as in (BC10).

For the binary step work on \(\Omega=T_i\cup T_j\). The cut with
an endpoint outside \(\Omega\) is unchanged. Use the signed
unit pair coordinate \(h=(p_i-p_j)/\sqrt2\) and tangential
coordinate \(z\). Outside \(\Omega\) each normal fiber is a
centered interval, up to endpoints, because the larger of
the two incident site scores is their mean plus
\(|h|/\sqrt2\), while all other scores are independent of \(h\).
With \(\eta=1_{\{h>0\}}\) the scalar sign proof in (J1) applies,
keeping also the nonnegative slack as in (BC6).

For \(h>0\), the integral of \(\operatorname{sgn}(h-q_h)\)
against \(K_e(p,q)\) on the excluded fibers is nonnegative.
In the conditional distribution of \(q\) given \(p\), replace
the sign threshold \(h\) by the smaller \(\rho h\).
On every excluded centered interval, the conditional normal
density is an even Gaussian centered at \(\rho h\ge0\).
Reflection about \(\rho h\) shows the portion below the mean
has at least the mass of the portion above it. The full-space
sign integral is
\(\phi_1(ep)[2\Phi((1-\rho)h/a)-1]\le C e^2 h\).
Thus the harmful contribution of the row integral on \(\Omega\)
is bounded by \(C e^2|h|\), and the negative \(h\) case follows
by reflection. For the centered Euclidean displacement kernel,
the full-space integral is zero and conditional evenness
gives the same sign on excluded intervals.

To bound the slack, take the local-step chain above with one
modification: the tangential increments lie in a small ball
where the increase of the average site score of \(i,j\)
exceeds that of every other site, by a fixed amount.
Such increments exist along the average pair direction.
The half-width of the excluded normal interval therefore
decreases by at least a fixed amount at each step until it is
zero. Choose the normal contraction so that the decrease
in \(|h|\) is less than this fixed amount.
A point of \(\Omega\) then stays in \(\Omega\) throughout
the chain. Tangential escape, pair congestion and the lower
bound on the kernel are unchanged. This gives the same
polynomial slack, which absorbs the row error in the
Gaussian physical ball of radius \(\delta\).
In Euclidean unbounded finite-volume changes, start on a
bounded subset of positive wrong mass and use the same
vanishing bound on incorrect terminal points outside \(B_R\).

For a general directed change of the stated class, write
\(h_U=1_U-1_T=\sum_{i\to j}(e_j-e_i)f_{ij}\), with \(f_{ij}\)
the transfer indicators. The exact deficit is \(2L-Q\) and
the linear score cost adds. The coefficient of a cross term
in \(Q\) is the inner product of two label-transfer vectors.
Distinct arrows have nonpositive inner product when they
share neither source nor destination. Since the kernel is
nonnegative, discarding those cross terms only lowers the
deficit. Each remaining term is the deficit of its single
transfer and obeys the binary bound just proved. This extends
the cyclic comparison (J2) to every cell count and to
the stated near-one Gaussian radius.

## A fixed-radius junction bound with arbitrary incoming transfers

Let \(T\) be the regular \(k\)-cell simplex partition in
\(E=\mathbf1^\perp\), \(k\ge3\), and let \(U\) be an ordinary competitor
whose changes are supported in the active rescaled ball \(B_R\),
\(R\ge1\). Write
\[
 F_{ij}=\{p:T(p)=i,\ U(p)=j\},\qquad
 h_{ij}(p)=(p_i-p_j)/\sqrt2,\qquad
 H=\sum_{i\ne j}\int_{F_{ij}}h_{ij}(p)\,dp.
\]
For Ornstein--Uhlenbeck noise put \(e=\sqrt{2(1-\rho)}\) and
\(a=\sqrt{(1+\rho)/2}\), and use the rescaled joint kernel
\[
 K_e(p,q)=\phi_1(ep)\phi_a(q-\rho p)
         =\phi_a(e(p+q)/2)\phi_1(q-p).
\]
Let \(\Delta_K\) denote the *ordered* relative cut
\(\iint K(p,q)(1_{U(p)\ne U(q)}-1_{T(p)\ne T(q)})\,dp\,dq\).
Then \(S_\rho(T)-S_\rho(U)=e^{k-1}\Delta_{K_e}\).
Define the unpaid shared-source interaction
\[
 B=\sum_i\ \sum_{\substack{j<\ell\\j,\ell\ne i}}
 \int_{F_{ij}\times F_{i\ell}} K(p,q)
  1_{\{h_{ij}(p)<h_{ij}(q),\ h_{i\ell}(q)<h_{i\ell}(p)\}}\,dp\,dq .
 \tag{BR1}
\]
There are \(c,\delta>0\), depending only on \(k\), such that for
\(eR\le\delta\)
\[
 \Delta_{K_e}\ge cR^{-2}H-2B.
 \tag{BR2}
\]
The same bound holds for the unit Gaussian Euclidean convolution kernel
on every \(R\ge1\). In particular the hypotheses permit arbitrarily many
transfer arrows into a label whenever each source has at most one arrow
out. They also permit outgoing branching whenever the interaction in
(BR1) vanishes or is sufficiently small. In the Gaussian case,
\(L_\rho(U;T)\le C e^{k-1}H\), so the zero-\(B\) class has the same
fixed-physical-ball relative \(cR^{-2}\) margin as (BC10). All assertions
allow unbalanced competitors.

### Single-arrow row credit

For each arrow \(i\to j\), let \(\Delta_{ij}\) be the ordered deficit of
the competitor that performs only this arrow. Changes in the cut with
one endpoint outside \(\Omega_{ij}=T_i\cup T_j\) vanish. On this domain
write \(\eta=1_{\{h>0\}},\,m=\eta-1_{F_{ij}}\), \(h=h_{ij}\) and
\(s(p,q)=\operatorname{sgn}(h(p)-h(q))\). The binary sign identity is
\[
 \Delta_{ij}
  =\iint_{\Omega_{ij}^2}K(p,q)
       [\,|m(p)-m(q)|-s(p,q)(m(p)-m(q))\,]\,dp\,dq
       -2\int_{F_{ij}}\int_{\Omega_{ij}} K(p,q)s(p,q)\,dq\,dp.
 \tag{BR3}
\]
Indeed the reference bracket vanishes, while antisymmetry converts
the signed change into twice its row integral. If \(p\in F_{ij}\)
and \(q\in T_i\setminus F_{ij}\) has \(h(q)<h(p)\), the bracket is two
at both \((p,q)\) and \((q,p)\). Therefore the first integral is at
least \(4A_{ij}\), where
\[
 A_{ij}=\int_{F_{ij}\times(T_i\setminus F_{ij})}
             K(p,q)1_{\{h(q)<h(p)\}}\,dp\,dq .
\]

For the second part, resolve \(q=z+r n\),
\(n=(e_i-e_j)/\sqrt2\). At fixed \(z\), the part outside \(\Omega_{ij}\)
is a centered interval \((-a_z,a_z)\), possibly empty. If
\(r=h(q)\) is in this interval, and \(b=T(q)\), then
\(r+a_z=\sqrt2(q_b-q_j)\). This follows since \(q_b\) and
\((q_i+q_j)/2\) do not depend on \(r\) and
\(a_z=\sqrt2(q_b-(q_i+q_j)/2)\).
Fix \(p\in F_{ij}\), and write \(u=h(p)>0\).
With \(q\) conditioned on \(p\), its normal density is Gaussian
centered at \(\rho u\). Each \(r>u\) inside \((-a_z,a_z)\) reflects
about \(u\) to \(2u-r\) in the same interval, below \(u\).
The density at \(2u-r\) is at least that at \(r\), because its
distance to the center \(\rho u\) is smaller or equal.
The remaining below-\(u\) strip is exactly
\((-a_z,\min(a_z,2u-a_z))\).
Consequently the excluded part of the sign row is bounded below by
\[
 \int_{E\setminus\Omega_{ij}}
     K(p,q)1_{\{h_{T(q),j}(q)<h_{ij}(p)\}}\,dq .
 \tag{BR4}
\]
The full-space sign row is
\(\phi_1(ep)[2\Phi((1-\rho)u/a)-1]\le C e^2 u\).
The Euclidean identity uses center \(u\), so the pairing is exact and
the full-space row vanishes. Combining with (BR3) gives
\[
 \Delta_{ij}\ \ge\ 4 A_{ij}+2 I_{ij}
                   -Ce^2\int_{F_{ij}}h_{ij},
 \qquad
 I_{ij}=\sum_{b\notin\{i,j\}}
       \int_{F_{ij}\times F_{bj}}K(p,q)
       1_{\{h_{bj}(q)<h_{ij}(p)\}}\,dp\,dq .
 \tag{BR5}
\]
In the Euclidean case the last term is zero.

### Paying shared destinations and retaining a margin

Write \(1_U-1_T=\sum_{i\ne j}(e_j-e_i)1_{F_{ij}}\).
The exact expansion \(\Delta_K=2L-Q\), with the linear transfer cost
\(L\) and quadratic interaction \(Q\), shows that two distinct arrows
with a shared source or shared destination cost
\(2\int_{F_{ij}\times F_{b\ell}}K\) once in the expansion of the
deficit relative to the sum of all single-arrow deficits. All other
distinct-arrow contributions to the deficit are nonnegative. This uses
\[
 (e_j-e_i)\cdot(e_\ell-e_b)
 =\delta_{j\ell}+\delta_{ib}-\delta_{jb}-\delta_{i\ell}.
\]
For a pair with common destination \(j\), the interaction appears exactly
once in \(\sum_{ij}I_{ij}\), according to which of the two transfer gaps
is larger. Gap equalities have zero measure, since the difference
is a nonzero linear function of \((p,q)\). Thus the incoming costs are
paid by the \(2I_{ij}\) terms. For shared source \(i\), the pair appears
in \(A_{ij}\) if \(h_{ij}(p)>h_{ij}(q)\), or in \(A_{i\ell}\)
if \(h_{i\ell}(q)>h_{i\ell}(p)\). If either holds, the full expense can
be charged to \(2\sum_{ij} A_{ij}\). The other pairs cost \(2B\).
In particular
\[
 \Delta_{K_e}\ge 2\sum_{ij} A_{ij}-Ce^2 H-2B.
 \tag{BR6}
\]

For completeness, the remaining binary margin in (BR6) can be proved
by a chain in the single reference cell. Fix \(i,j\), and write
\(p=z+h n\) using the normal coordinate above. Start with density
\(h1_{F_{ij}}\,dp\). In tangential space choose a small open ball \(W\)
of increments \(w\) with \(w_i=w_j\) and
\(w_i-w_b>\epsilon\) for every \(b\notin\{i,j\}\), bounded away
from zero in the direction \(w_*=(e_i+e_j-2\mathbf1/k)/
|e_i+e_j-2\mathbf1/k|\). Such a ball exists about a small positive
multiple of \(w_*\). At each step add an independent uniform \(w\in W\)
to \(z\) and contract \(h\) by an independent \(1-\theta/N\) with
\(\theta\) uniform on \([1,2]\), where \(N=\lceil C_0(R+1)\rceil\).
For \(C_0\) sufficiently large, \(h\) stays positive and decreases,
while \(p_i-p_b\) strictly increases for \(b\notin\{i,j\}\).
The path stays in \(T_i\), ends outside \(B_R\), stays throughout
\(B_{C_kR}\), and every step has uniformly bounded length.

For fixed increments, the maps translate \(z\) and contract \(h\) by
a factor bounded below by \((1-2/N)^N\ge1/16\). The marginal density
at each step is therefore at most \(C h\), and the conditional
one-step density at most \(CN/h\). Thus every successive pair has
joint density at most \(CN\). The path must leave \(F_{ij}\)
and its first such step is among the pairs of \(A_{ij}\).
On all paths \(K_e\ge c\exp(-C e^2R^2)\), using the midpoint form
and bounded displacements. A union bound for the crossing step
gives \(A_{ij}\ge cR^{-2}\int_{F_{ij}} h\) when \(eR\le1\).
The same applies to the Euclidean unit Gaussian. Summing in (BR6)
and decreasing \(\delta\) absorbs the row drift.
Lastly the reference \(i,j\) scores tie when \(h=0\) by reflection.
Their uniform rescaled Lipschitz bound and the bounded Gaussian row
factor yield the stated upper bound on \(L_\rho\). The Gaussian transverse lift for a fixed global outdegree-one class,
or for slices satisfying the same uniform zero or small-\(B\) condition,
uses positivity of the transverse OU operator as in the two-block
comparison. Formula (BR2) itself uses active-space interaction.

### Scope of the remaining expense

For a genuine shared-source event, let \(k=4\) and take \(i=1,j=2,\ell=3\).
Choose small balls about the centered projections of
\((2,1,0,-1)\) and \((2,0,1,-1)\), after common scaling to place
both inside the support ball. Both balls lie in \(T_1\). Assign labels
2 and 3 to the two balls, respectively. Their cross interaction
lies entirely in (BR1), with strict inequalities throughout for
sufficiently small balls. Hence its positive lower bound cannot be
paid by the two corresponding binary chord contributions described
above. This is a statement about those terms; the ordinary
competitor's total deficit is a separate quantity. A general outgoing
bound needs an additional curvature, row, or joint-label estimate.

The terms discarded in (BR6) can also be kept. Let \(G\) be the
same-source interaction counted once over pairs complementary to
(BR1). Let \(N\ge0\) be the contribution to the deficit from pairs
of distinct arrows whose label-transfer vectors have negative
inner product. The unused excluded-row credit is
\[
 D_{\rm row}=\sum_{i\ne j}
    \int_{F_{ij}}\int_{E\setminus(T_i\cup T_j)}
       K(p,q)1_{\{h_{T(q),j}(q)<h_{ij}(p),\ U(q)\ne j\}}\,dq\,dp .
                                                               \tag{BR7}
\]
The binary bracket (BR3) gives exactly \(4A_{ij}\), because it can
be nonzero only when the flipped point and the other endpoint
have reference label \(i\). The excluded-row reflection is exact
for Euclidean convolution by the unit Gaussian. Thus the Euclidean
identity and Gaussian lower bound are respectively
\[
 \begin{split}
 \Delta_K&=4\sum_{ij}A_{ij}-2G-2B+2D_{\rm row}+N,\\
 \Delta_{K_e}&\ge4\sum_{ij}A_{ij}-2G-2B+2D_{\rm row}+N-Ce^2H .
 \end{split}                                                     \tag{BR8}
\]
All quantities in each line use its own kernel. The unused strip
and the negative interaction terms remain available to pay for
\(B\).

Two classes allow arbitrary incoming and outgoing degrees and still give
the full relative bound.

First, for each source \(i\) fix an allowed target set \(J_i\) and
constants \(b_{ij}\). At every changed point of \(T_i\) assign a
maximizer of \(h_{ij}(p)+b_{ij}\) over \(j\in J_i\). The region
where changes occur can be arbitrary. For \(p\in F_{ij}\) and
\(q\in F_{i\ell}\) the two maximization inequalities add to
\[
     h_{ij}(p)-h_{ij}(q)+h_{i\ell}(q)-h_{i\ell}(p)\ge0 .
\]
The opposed strict inequalities in (BR1) are then impossible and
\(B=0\).

Second, the margin also holds if every pair counted in (BR1) has
separation at least \(D_R=C_k\sqrt{\log(2+R)}\), up to a null set.
For such pairs \(K_e(p,q)\le C\exp(-cD_R^2)\phi_{\sqrt2}(p-q)\)
by the midpoint formula. The same bound holds for the Euclidean
kernel. The thick/thin section estimate in
[(ET4)](EUCLIDEAN_SIMPLEX_JUNCTION_TRANSFER.md) bounds the
resulting absolute interaction by \(C_k e^{-cD_R^2}H\).
Increasing the constant in \(D_R\) lets (BR2) absorb it.
In particular only the pairs with opposed gaps need separation.
For an ordinary active-space tie or improvement, (BR2) requires
\(B\ge (c/2)R^{-2}H\).

The \(B\) term can exceed the margin reserved in (BR2) by a growing
factor even in an ordered example. Fix a small \(\tau>0\), and
inside \(B_R\) lower the score of one source \(i\) by \(\tau\)
while keeping every other score unchanged. This is an ordered
two-block change on the split \(\{i\}\) and its complement.
Each outgoing transfer lies in a normal strip of width
\(\tau/\sqrt2\), so \(H\le C_k\tau^2 R^{d-1}\).
For two distinct rivals \(j,\ell\), choose a tangential patch inside
the \((i,j,\ell)\) triple face of volume at least \(c_k R^{d-2}\),
with all spectator gaps at least \(c_k R\). For \(d=2\) the
tangential space is a point and the spectator condition is vacuous.
Such a patch for \(d\ge3\) is obtained by scaling an interior ball
of the triple-face cone inside \(B_{R/2}\). Take two normal pairs
with raw gaps
\[
 (p_i-p_j,p_i-p_\ell)\in[\tau/4,\tau/2]\times[1,2],
 \qquad
 (q_i-q_j,q_i-q_\ell)\in[1,2]\times[\tau/4,\tau/2].
\]
Their nearest rivals are distinct, their transfers are prescribed
by the score decrease, and their gaps have the opposing order.
The two gap coordinates are nonsingular on the normal plane.
Restrict the tangential separation to at most one, taking the first
endpoint inside a smaller concentric patch. Both endpoints stay
in \(B_R\) for \(R\) sufficiently large. The kernel is bounded below
on these pairs for \(eR\le\delta\), using bounded displacement and
midpoints of size \(O_k(R)\). It follows that
\(B\ge c_k\tau^2R^{d-2}\). Consequently \(B/H\ge c_k/R\) on this
family. A general estimate that tries to absorb \(B\) into just the
\(R^{-2}H\) term of (BR2) cannot close the argument. In this example
the ordered theorem supplies the full comparison.

The starting scalar sign identity is the binary calibration in
[Pagliari, Theorem 2.3 and Lemma 2.4](https://arxiv.org/pdf/1905.00623).
The strip credit above uses the exact excluded-fiber geometry of the
simplex and the ordinary transfer expansion. The
[Lawlor--Morgan local tetrahedral calibration](https://msp.org/pjm/1994/166-1/pjm-v166-n1-p04-p.pdf)
addresses local area and does not yield this finite-noise row estimate.

## Balanced fractional obstruction for the small physical core

Fix \(k\ge3\) and work in its active space of dimension \(d=k-1\).
For any fixed physical radius \(r>0\), for all \(\rho\) sufficiently
close to one the simplex indicators fail to minimize the convex
componentwise total-variation cut among balanced simplex-valued
functions that agree with those indicators outside \(B_r\).
Here the cut is
\[
 \mathcal C_\rho(u)=\frac12\sum_{i=1}^k
       \mathbb E\,|u_i(X)-u_i(Y)|,\qquad |u|_1=1,\quad u_i\ge0 ,
\]
for the correlated Gaussian pair in the active space. Balance is
\(\int u_i\,d\gamma_d=1/k\). The corresponding Euclidean unit
Gaussian relative cut has a strict fractional improvement with
bounded support for every \(k\ge3\), as well. This extends the
planar translation example (J3)–(J4) of the
[Euclidean three-label comparison](EUCLIDEAN_THREE_CELL_JUNCTION.md#a-fractional-obstruction-to-convex-total-variation-calibration).
It concerns the convex cut, which agrees with the ordinary cut
for hard labels.

### Volume-preserving translations with cutoff

Use scaled coordinates \(p=X/e\) and the kernel (BC7), so that
Gaussian cut differences are \(e^d\) times scaled cut differences.
Write \(\omega_e(p)=\exp(-e^2|p|^2/2)\), polarize the active space
with orthonormal coordinates, and fix a vector \(v\) of bounded
length in that space. Let \(\chi_R\) equal one on \(B_R\), vanish
outside \(B_{2R}\), and have derivatives of order \(j\le3\) bounded
by \(C_j R^{-j}\). Define
\[
 A_{ab}(p)=\omega_e(p)\chi_R(p)
              \frac{v_a p_b-p_a v_b}{d-1},\qquad
 V_a(p)=\omega_e(p)^{-1}\sum_b\partial_b A_{ab}(p).
\]
Antisymmetry gives \(\operatorname{div}(\omega_e V)=0\), so the
smooth compactly supported flow \(F_t\) of \(V\) preserves
Gaussian volume in scaled coordinates. The same definition
with \(e=0\) preserves Lebesgue measure. Since
\(\nabla\log\omega_e=-e^2p\), direct differentiation gives
\[
 \|V\|_\infty\le C,\quad
 \|DV\|_\infty\le C/R,\quad
 \|D^2V\|_\infty\le C/R^2,\quad
 |V(p)-v|\le C(eR)^2\quad(p\in B_R),
\]
provided \(eR\le\delta_0\). For the first term with no cutoff
derivative, the identity
\(\sum_b\partial_b[(v_a p_b-p_a v_b)/(d-1)]=v_a\) was used.
Terms with \(\nabla\chi_R\) have respectively size \(O(1)\),
\(O(R^{-1})\), \(O(R^{-2})\). Terms from
\(\nabla\log\omega_e\) have bounds \(O(e^2R^2)\),
\(O(e^2R)\), \(O(e^2)\), respectively.

Put \(T(t)=F_tT\) and let \(J(t)\) be its ordinary scaled
cut minus that of \(T\). For fixed sufficiently small \(|t|\)
and \(R\ge1\),
\[
                         |J(t)|\le C t^2 R^{d-3}.             \tag{BC11}
\]
The first derivative is zero. Differentiating the two
cell-indicator integrals gives the sum of normal speeds
against the differences of the two incident reference scores
on the facets. Reflection fixes each facet point and exchanges
those cells, so the differences vanish. The volume derivatives
are justified by the divergence theorem on the polyhedral
cells and the compactly supported smooth flow.

For the second derivative pull back both variables and integrate
over the fixed differing-label pairs of \(T\).
Let \(Z_t=F_t(p)-F_t(q)\) and \(W_t=e(F_t(p)+F_t(q))/2\).
For bounded \(|t|\), \(|Z_t|\) is comparable to \(|p-q|\) by
the Lipschitz bound and invertibility of the flow. Differentiating
in \(t\) gives
\[
 |\dot Z_t|\le C |p-q|/R,\quad
 |\ddot Z_t|\le C |p-q|/R^2,\quad
 |\dot W_t|\le Ce,\quad |\ddot W_t|\le Ce/R .
\]
The second displacement estimate uses the Lipschitz norm of
\(DV\cdot V\), bounded by \(\|D^2V\|\|V\|+\|DV\|^2\).
The Jacobians have bounded size, and their first and second
time derivatives are \(O(R^{-1})\) and \(O(R^{-2})\), from
the log-Jacobian equation with \(\operatorname{div}V\).
In the kernel \(\phi_a(W_t)\phi_1(Z_t)\), the midpoint
derivatives are bounded uniformly as functions of \(W_t\),
with \(a^2=(1+\rho)/2\ge1/2\). The first and second time
derivatives of this factor are therefore \(O(R^{-1})\) and
\(O(R^{-2})\) for \(eR\le\delta_0\). Differentiating the
displacement factor gives a Gaussian times a polynomial in
\(|p-q|\), with the same powers of \(R^{-1}\).
Thus the second derivative of the pulled-back kernel, including
the two Jacobians, is bounded by
\[
 C R^{-2}(1+|p-q|^4)e^{-c|p-q|^2}.
\]
It vanishes if both points are outside \(B_{2R}\).
For fixed displacement \(p-q=z\), the differing-label pairs
meet a strip of width \(2|z|\) about one of the finitely many
pair planes, since the joining segment crosses such a plane.
The moving endpoint condition imposes membership of \(p\) in
\(B_{2R}\) or its translate by \(z\). The corresponding
volume in the strips is at most \(CR^{d-1}|z|\), uniformly
in the translate. Integrating in \(z\) gives
\(|J''(t)|\le C R^{d-3}\) and proves (BC11).
The same derivation for \(e=0\) uses a constant midpoint factor.

### A strict Jensen gap at an entire triple face

Choose distinct labels \(i,j,l\) and a unit \(v\) such that
\((v_i-v_j)(v_i-v_l)<0\). In the active space their triple
tie subspace \(H\) has dimension \(d-2\); the relatively
open face where the three sites beat the others is
nonempty. If \(d>2\), take \(z_*\) relatively interior and
unit in that face and tangential points \(z\) in a
small fixed-fraction ball in \(H\) about \((R/8)z_*\).
Its volume is \(cR^{d-2}\). Choose its fraction so that
the three sites beat all others by \(cR\) there and
\(|z|\le R/4\), and shrink it by a fixed factor to allow
each second tangential point in a unit ball about the
first. If \(d=2\) take the single tangential point,
with zero-dimensional volume one.

In \(H^\perp\) choose two bounded strips of width \(ct\)
and fixed positive length along the two open pair rays
\(i/j\) and \(i/l\) of the normal tripod. The strips
are separated from the third ray and from the junction.
They can be centered on their respective pair lines.
For \(0<t\le t_0\) fixed sufficiently small, on one strip
the indicator of label \(i\) in \(T(t)\) is one and
in \(T(-t)\) is zero, while the reverse holds on the
other strip, after interchanging the strips if necessary.
Indeed in this region
\(F_{\pm t}^{-1}(p)=p\mp tv+O(t\delta_0^2)\) by the
uniform bound on \(V-v\) and by remaining inside \(B_R\).
The opposite normal projections of \(v\) make the flips
robust to the latter error when \(\delta_0\) is small.
The order-one gaps on the strips to the third site,
and the order-\(R\) gaps to the remaining sites,
prevent any other label from winning. This holds for
large \(R\) uniformly in the indicated tangential points.

Let \(u=(1_{T(t)}+1_{T(-t)})/2\).
On a pair of points from the opposite strips with
tangential coordinates at distance at most one,
the label-\(i\) absolute difference is one for each
deformation and zero for \(u\). The scaled kernel
is bounded below there: the pair displacement is
bounded and \(eR\le\delta_0\). Such pairs have product
volume at least \(c t^2R^{d-2}\). Convexity for every
other pair and (BC11) consequently give
\[
 \mathcal C_e(u)-\mathcal C_e(1_T)
 \le C t^2 R^{d-3}-c t^2R^{d-2}<0
\]
for sufficiently large \(R\), writing \(\mathcal C_e\)
for the scaled cut (or its relative version).
Taking \(R=\min(\delta_0,r)/(2e)\) produces the claimed
Gaussian support and preserves balance by construction.
The Euclidean conclusion follows by letting \(R\) be large
with \(e=0\).

### Sharpness of the ordered-ball exponent

The power \(R^{-2}\) in (BC10) is sharp up to constants, even
for balanced competitors. In the flow for (BC11), choose \(v\)
constant on each of the two blocks, with their two constants
distinct, and choose the cutoff radial. The formula for \(V\)
then has the form
\(V(p)=a(p)v+b(p)(v\cdot p)p\) with scalar radial coefficients.
Under the flow every within-block coordinate difference
satisfies a scalar linear differential equation, so its sign
is preserved. Thus \(F_tT\) remains in the ordered class.
It is balanced by preservation of \(\omega_e\)-volume.

For \(R\) large, choose a cross-block reference facet patch
of \((d-1)\)-volume \(cR^{d-1}\) inside \(B_{R/2}\), away
from every other facet by a fixed positive fraction of \(R\).
The direction \(v\) crosses its pair plane. Since
\(V=v+O((eR)^2)\) on \(B_R\), for small fixed \(t>0\)
the flow transfers a strip of thickness \(ct\) on this patch.
Integrating the saturated score margin over a substrip at
normal distances between \(c't\) and \(ct\), and using
\(eR\le\delta\), gives \(L\ge c t^2 R^{d-1}\) in scaled units.
Bound (BC11) gives \(\Delta\le C t^2R^{d-3}\) in the same
units. The support is \(B_{2R}\), which only changes the
constants and the admissible small physical radius by fixed
factors. This yields \(\Delta/L\le C R^{-2}\).
In Euclidean coordinates the same example uses the weighted
defect \(\int |g|1_{\{m\ne\eta\}}\) in place of \(L\).
