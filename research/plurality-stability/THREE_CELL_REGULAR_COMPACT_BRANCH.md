# Analytic reduction near a regular compact three-cell maximizer

Fix \(0<r_0<1\). Let \(A^0\) be a balanced actual planar Gaussian
maximizer, with a pure exterior label, numbered 3, so labels 1
and 2 are supported in a disk. Use the normalized adjusted scores
\(V_i^0=U_i^0-\beta_i^0\), \(U_i^0=T_{r_0}1_{A_i^0}\),
\(\beta_3^0=0\). Assume that their maximal ties form a finite
compact embedded network of real-analytic pair edges, with the
following regularity.

- At every open edge exactly its incident pair is maximal and
  its pair-score difference has nonzero gradient.
- Every vertex is a triple equality where
  \((V_1^0-V_3^0,V_2^0-V_3^0)\) has invertible planar
  differential.
- There are no other maximal ties.

Closed pair edges without vertices are allowed. The conditions
exclude dormant maximal ties and collinear score gradients at a
triple. The local stationary score-winner equations in this class reduce to a
finite-dimensional real-analytic system, with no nondegeneracy
assumption on the constrained shape Hessian.

In particular suppose that \(A^0\) ties the centered regular
simplex \(T\), and that balanced actual maximizers \(A^j\to A^0\)
in planar Gaussian label distance at \(r_j\to r_0\) satisfy
\(S_{r_j}(A^j)>S_{r_j}(T)\). Suppose the \(r_j\) lie on one
fixed side of \(r_0\). Then there are \(c>0\) and positive
integers \(m,n\) such that throughout a sufficiently short interval
on that side,
\[
            V(r)-S_r(T)\ \ge\ c|r-r_0|^{n/m},              \tag{RB1}
\]
where \(V(r)\) is the balanced planar maximum. In particular
strict failure cannot appear only intermittently or with an
infinitely flat value gap along this regular compact branch.
The constructed competing branch need not maximize. The
statement does not exclude the tie.

If every nonsimplex balanced planar maximizer at the tied
correlation \(r_0\) satisfies these compact regular hypotheses,
then the exact global value gap has a convergent one-sided
Puiseux expansion on each side of \(r_0\). Equation (RB5)
states this stronger finite-chart consequence.

## Score convergence and coordinates

Pure-exterior compactness puts labels 1 and 2 of all sufficiently
late \(A^j\) in one common disk. The adjusted score of label 3
for \(A^0\) strictly dominates far outside that disk: its
unadjusted score tends uniformly to one there, the other two
tend to zero, and every price difference is strictly less than one
in magnitude. Choose a disk \(D\) containing the common support
with a fixed collar of strict label 3 for the reference scores.
[Pure-exterior compactness](THREE_CELL_PURE_EXTERIOR_COMPACTNESS.md)
[Score prices and complete adjacency](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#complete-adjacency)

Choose a bounded, conjugation-invariant complex neighborhood
\(\Omega\subset\mathbb C^2\) of \(\overline D\). Use the
real Banach space \(X\) of pairs of holomorphic functions on
\(\Omega\), continuous on \(\overline\Omega\) and real on the
real disk, with the supremum norm. Its complexification uses
the same holomorphic space without the reality condition.
The pair of unadjusted scores \((U_1,U_2)\) determines
\(U_3=1-U_1-U_2\). The label prices use the two coordinates
\((\beta_1,\beta_2)\) and \(\beta_3=0\).

In this space the actual scores of \(A^j\) converge uniformly
on \(\overline\Omega\) to those of \(A^0\). For real \(t\) in
a compact band in \((0,1)\) and complex observations in a
bounded set, the relative Mehler kernel and its \(t\) derivative
are bounded uniformly in the real source variable. The
negative quadratic source exponent absorbs the linear exponent
in the complex observation and the polynomial from
differentiation. Strong Gaussian label convergence and
\(r_j\to r_0\) then give convergence in the required norm.
The normalized prices of \(A^j\) are bounded because each
label wins on a positive-measure set and the unadjusted
scores lie in \((0,1)\). Any subsequential price limit describes
the winners of \(A^0\), by continuity of the scores and passage
to a subsequence with almost-everywhere indicator convergence.
The normalized prices of \(A^0\) are unique. Indeed any two
choices have the same pair differences by the equalities on
the regular active pair interfaces and complete adjacency.
The prices of \(A^j\) therefore converge to \(\beta^0\).

## Analytic winner integrals

For real data \(p=(U,\beta)\) near \(p_0=(U^0,\beta^0)\) in
\(X\times\mathbb R^2\), let \(B(p)\) be the ordinary winner
partition of \(D\) for \(U_i-\beta_i\), completed by label 3
outside \(D\). Near the boundary of \(D\), label 3 remains the
unique winner by the strict collar margin. Away from a small
neighborhood of the reference network, strictness and
compactness keep the reference winners. Inside the neighborhoods
the reference regularity gives the same labeled graph with
perturbed edge curves, as follows.

On regular pair portions use coordinates \((u,v)\) for which
the transverse derivative in \(v\) of the pair difference is
nonzero. Solving that difference equal to zero gives an analytic
graph in \(u\), depending analytically on \(p\). The spectator
remains strictly below it there. At a vertex use the two
contrast coordinates \((V_1-V_3,V_2-V_3)\). Their local inverse
depends analytically on the data by the implicit function
theorem. In these coordinates the winner interfaces are three
fixed rays meeting at zero, so the inverse gives the junction
and its three outgoing edges. The local inverse also ensures
that no additional winning edges occur in the vertex chart.
The graphical descriptions hold for nearby real data because
the uniform holomorphic norm controls real first derivatives
on smaller neighborhoods by Cauchy estimates.

More explicitly, cut off a short segment of each outgoing
ray in its vertex coordinates at a fixed positive parameter.
Its image under the local inverse supplies an analytic cut
point in the overlapping regular pair chart. Subdivide the
remaining regular arcs into finitely many graphical segments
with fixed transverse sections at intermediate cuts. Cut
points on these sections depend analytically on \(p\) by the
scalar implicit function theorem. For a closed edge insert
such cuts along the edge. In each graph chart the integration
endpoints are the projections of its two cut points onto the
graph coordinate. The endpoints therefore vary analytically.
After affine rescaling of the graph coordinate, each perturbed
segment is parametrized on a fixed compact interval, with
its orientation inherited from the reference network.
The same implicit functions extend to complex data and
complex parameters in neighborhoods of the segment intervals.
The evaluation \((U,z)\mapsto U(z)\) needed here is jointly
holomorphic on a smaller spatial chart by Cauchy estimates,
including when the variable \(U\) belongs to the Banach space.
These are finite chart constructions, so the neighborhoods of
the input data can be chosen common to all segments.

Write \(M_i(p)=\gamma_2(B_i(p))\) and
\(G_{r,i}(p)=T_r1_{B_i(p)}\) for \(i=1,2\). Both maps are
real analytic in the relevant inputs. The following boundary
formula proves this without differentiating an indicator in
a normed space. For an entire density \(F(y)\), integrate its
one-form
\[
                -\left[\int_0^{y_2}F(y_1,z)\,dz\right]dy_1
\]
along the oriented boundary of \(B_i(p)\). The exterior
derivative of this form is \(F(y)\,dy_1\wedge dy_2\).
Thus Green's formula gives \(\int_{B_i(p)}F(y)\,dy\) as the
sum of the finitely many line integrals. The vertex points
have no separate term. This applies to Gaussian mass density
and to the conditional Gaussian density
\[
  F_{r,w}(y)=\frac{1}{2\pi(1-r^2)}
             \exp\!\left[-\frac{\sum_{\ell=1}^2(y_\ell-rw_\ell)^2}
                               {2(1-r^2)}\right].
\]
The squares here do not use complex conjugation. The latter
density is entire in \(w,y\) and holomorphic in \(r\) near
\(r_0\). For complex data define the extension by the very
same line integrals on the chart curves, without a complex
winner rule. The parametrizations and their derivatives
are uniformly holomorphic in the data on the fixed integration
intervals, using the finite charts and Cauchy estimates.
Integration on these intervals yields holomorphic maps
of the complexified input data, with uniform bounds on
any fixed larger observation neighborhood
\(\Omega^+\supset\overline\Omega\). In particular \(G_r\)
is analytic into the holomorphic supremum-norm space for
\(\Omega^+\). Its derivative as a map into \(X\) is compact:
restrict the bounded derivative into the larger space to
\(\Omega\), using the Cauchy derivative bounds and
Arzela–Ascoli on \(\overline\Omega\).

## Finite-dimensional stationary equations with exact masses

The system is
\[
 \mathcal F(r,p)=
    \big(U-G_r(p),\ (M_1(p),M_2(p))-(1/3,1/3)\big)=0,
                                                              \tag{RB2}
\]
with values in \(X\times\mathbb R^2\). A real root produces
an ordinary partition with all masses \(1/3\), as the
third follows from the first two. The scores in the disk
are its actual unadjusted scores. The exterior label can
also be checked directly on the roots: the two compact cells
remain supported in a fixed smaller disk, and enlarging \(D\)
initially makes their actual scores uniformly small outside
\(D\) for \(r\) near \(r_0\). The price margins for label 3
there persist from \(\beta^0\). Only exact balance and
ordinary partition status are needed for the value
comparison below. No translation or rotation quotient is
taken, and the masses are equations in the system rather than
first-order constraints.

At the reference point the partial derivative
\(L=D_p\mathcal F(r_0,p_0)\) is identity minus compact on
\(X\times\mathbb R^2\). The first component has the
compact score derivative just proved, and all additional
price or mass blocks are finite rank. Thus \(L\) is Fredholm
of index zero. Split its finite-dimensional kernel from a
closed complement and its range from a finite-dimensional
cokernel. The range projection of (RB2) is invertible
in the complement variable at the reference point, because
\(L\) restricts there to an isomorphism onto its closed range.
The analytic implicit function theorem solves this
component as an analytic function of \(r\) and of the
kernel variable \(q\). The remaining equations are
\[
                        \Psi(r,q)=0,                         \tag{RB3}
\]
a finite-dimensional real-analytic system. This construction
represents all sufficiently close real roots of (RB2),
with \(p\) an analytic function of \((r,q)\).

On a real root the stability is
\[
       S_r(B)=\frac13+
           \sum_{i=1}^2\int_{B_i}
                     (G_{r,i}-(1-G_{r,1}-G_{r,2}))\,d\gamma_2 .
                                                              \tag{RB4}
\]
Indeed the full integral of the third score is the third
mass \(1/3\), by invariance of Gaussian measure. Subtract
its integrals over \(B_1,B_2\) and add the first two
self-score integrals. The displayed expression is analytic
on the chart before imposing the root equations as well,
using the same boundary integrals since the score outputs
are entire in the observation. Hence the deficit relative
to the simplex, with sign reversed, is analytic in the
finite variables on the stationary chart. Write
\(E(r,q)=S_r(B(r,q))-S_r(T)\) for this analytic gain,
using (RB4) on the roots. The simplex stability is analytic
in \(r\) by its Hermite expansion.

## Selecting an improving branch

The convergent improving maximizers in the statement are
represented by real roots of (RB3), because of the common
support disk, score convergence, and price convergence.
Let \(\pm\) be the fixed sign of \(r_j-r_0\). Their reduced
parameters place \((r_0,0)\) in the closure of the
finite-dimensional semianalytic set
\[
                 \{\Psi=0,\ \pm(r-r_0)>0,\ E>0\}.
\]
The real-analytic curve selection lemma gives a parameter
\(u\) and a real-analytic arc \((r(u),q(u))\) from
\((r_0,0)\) into this set for all sufficiently small \(u>0\).
Since the reference partition ties, both the signed
correlation increment and the gain vanish at \(u=0\).
Their first nonzero coefficients are positive by the
strict inequalities on the selected arc. Thus
\[
 \pm(r(u)-r_0)=c_0 u^m+O(u^{m+1}),\qquad
 E(r(u),q(u))=c_1 u^n+O(u^{n+1}),
\]
for \(c_0,c_1>0\) and positive integers \(m,n\).
The first quantity is strictly increasing for small
positive \(u\) and fills an interval. The balanced
partitions on this arc are legitimate competitors, so
their gains imply (RB1).
[Real-analytic curve selection for an analytic set intersected with strict inequalities, Garcia, Fact 12](https://repositorio.usp.br/directbitstream/139e3878-bcbf-42f9-b201-00df4369ce1f/841054.pdf)

If the convergent sequence instead consists of exact ties,
use the analytic equation \(E=0\) in place of \(E>0\).
It gives a stationary branch attaining the simplex value
on a whole short interval on that side. Its members are
co-maximizers wherever the simplex is globally optimal.
They need not be simplex partitions, since the chart
partitions approach the fixed nonsimplex reference.
The reduction leaves open arbitrary networks with
dormant maximal ties, nontransverse junctions, or more
than one exterior label.

## Exact one-sided germs of the global value

Suppose now that \(T\) is globally optimal at \(r_0\) and that
every nonsimplex balanced planar maximizer at \(r_0\) satisfies
the compact regular network hypotheses, after relabeling if
needed. On each side of \(r_0\) there is a positive integer
\(m_\pm\) and a real-analytic function \(F_\pm\) defined near
zero such that the exact global gap satisfies
\[
  V(r_0\pm h)-S_{r_0\pm h}(T)=F_\pm(h^{1/m_\pm})
                      \qquad(0\le h<h_0).                  \tag{RB5}
\]
Thus a side has either equality throughout a short interval
or a convergent Puiseux series with a positive leading
coefficient on a whole short interval of strict failure.

The family of nonsimplex maximizers at \(r_0\) is compact
in planar Gaussian label distance. Indeed strong maximizer
compactness gives subsequential limits, and the uniform
measurable simplex comparison on a compact correlation
band gives a fixed positive distance of any nonsimplex
maximizer from the simplex orbit. This latter comparison
also excludes convergence to the simplex of nonsimplex
maximizers at correlations tending to \(r_0\): any such
maximizer within the local neighborhood would have strictly
smaller stability than the admissible simplex.
[Strong maximizer compactness](GAUSSIAN_NEAR_MAXIMIZER_COMPACTNESS.md)
[Uniform measurable simplex comparison, B11](POSITIVE_CORRELATION_CONTINUATION.md#an-interface-gap-implies-a-uniform-measurable-local-comparison)
[Three-cell interface gap](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem)

Take finitely many regular compact charts covering the
nonsimplex maximizing family at \(r_0\). More precisely, in
each chart take a closed bounded box of its reduced variables
\((r,q)\), with closure inside its analytic domain and
containing the base point in its interior. There is a smaller
label-distance neighborhood of that chart's base partition
whose nearby actual maximizers have their parameters in
the interior of this box. Otherwise a sequence of such
maximizers converging to the base in labels and correlation
would contradict pure-exterior confinement and the
score-price convergence already proved. Use these smaller
neighborhoods for the finite cover. After restricting to
one common short correlation interval \(I\) about \(r_0\),
every nonsimplex actual maximizer belongs to one of the
boxes. If this failed, compactness applied to a violating
sequence would give a nonsimplex maximizer at \(r_0\) in
one of the smaller neighborhoods. The sequence eventually
lies in that neighborhood, which places its actual-maximizer
parameters in the box by the choice of neighborhood. If there is no nonsimplex maximizer at \(r_0\),
the same compactness argument gives a simplex-only interval.

Within \(I\), every root in a chart box represents a balanced
ordinary competitor, and every nonsimplex global value is
captured by one of these boxes. Let
\(\Delta(r)=V(r)-S_r(T)\). The union of the baseline
rectangle \(\{r\in I,\ -1\le y\le0\}\) and the
\((r,y)\)-projections of the finitely many sets
\[
  \{(r,q,y):(r,q)\text{ in chart box }j,\ r\in I,\
             \Psi_j(r,q)=0,\ -1\le y\le E_j(r,q)\}
\]
is exactly \(\{r\in I,\ -1\le y\le\Delta(r)\}\).
Here all boxes and \(I\) are closed, lie inside the respective
analytic domains for the evaluations, and the sets being
projected are compact semianalytic. The equality of the
projected union with this hypograph band uses both facts:
a chart competitor never exceeds the global value, and a
global value above the simplex comes from an actual maximizer
that is covered by a chart. The global value and simplex
value are continuous by the uniform correlation derivative
bound of the Ornstein--Uhlenbeck operators on a compact
positive band, so the upper boundary of the union near
\((r_0,0)\) is the graph of \(\Delta\).

By Definition 3.1 and Theorem 6.1(2) of Bierstone--Milman,
these compact semianalytic projections form a subanalytic,
hence semianalytic, subset of the plane. Its boundary is
semianalytic by Corollary 2.8. On either side of \(r_0\)
the continuous graph is a one-dimensional semianalytic
set whose punctured germ is locally connected, since the
first coordinate identifies it with an interval. Apply
their Lemma 6.3(1) to parametrize each one-sided graph near
\((r_0,0)\) by a real-analytic arc \(u\mapsto(r(u),y(u))\).
For \(u>0\) small its signed correlation increment has form
\(\pm(r(u)-r_0)=u^{m_\pm}a(u)\), where \(a\) is analytic
and \(a(0)>0\). The analytic coordinate \(u a(u)^{1/m_\pm}\)
has nonzero derivative at zero, so composing \(y(u)\) with
its analytic inverse gives (RB5).
[Bierstone and Milman, Semianalytic and subanalytic sets, Definition 3.1, Corollary 2.8, Theorem 6.1(2), and Lemma 6.3(1)](https://www.numdam.org/item/10.1007/BF02699126.pdf)

The exact expansion also yields a branch of actual maximizers
on any side of strict failure. Choose a sequence \(h_j\downarrow0\)
on that side and put \(z_j=h_j^{1/m_\pm}\). A subsequence of
maximizers is represented in one fixed closed chart box and has
parameters \(q_j\to q_*\) within that box. The finite analytic
system in \((z,q)\)
\[
 \Psi(r_0\pm z^{m_\pm},q)=0,\qquad
 E(r_0\pm z^{m_\pm},q)=F_\pm(z),\qquad z>0
\]
therefore accumulates at \((0,q_*)\). Curve selection gives an
analytic arc out of that point into the set. The chart is
defined on a neighborhood of the closed box, so all roots
on the arc give ordinary balanced partitions. The value
equation makes them actual maximizers. The correlation
increment along the arc again fills a one-sided interval.
At each of its correlations \(r=r_0\pm z^{m_\pm}\), every
actual maximizer \(A\) has the same fixed-partition slope:
\[
 S'_r(A)-S'_r(T)
       =\pm\frac{F_\pm'(z)}{m_\pm z^{m_\pm-1}}.               \tag{RB6}
\]
Indeed \(V(t)-S_t(A)\) has a local minimum zero at \(t=r\),
and both functions are differentiable there. By continuity
of the slopes of strongly convergent partitions, if the
leading gap order exceeds one, the selected limiting tie
has the simplex slope at \(r_0\). The lower bound of one on
the leading order follows from the uniform Lipschitz bound
for the gap.
