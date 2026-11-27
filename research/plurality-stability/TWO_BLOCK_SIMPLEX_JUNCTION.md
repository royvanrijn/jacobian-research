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
Lebesgue measure. The support range in (BC1) grows in noise
lengths. The comparison still requires the two block orientations.

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
