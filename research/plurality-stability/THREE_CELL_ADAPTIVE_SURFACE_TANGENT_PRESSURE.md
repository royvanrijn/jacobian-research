# Pressure selection at an adaptive surface tangent

For an adaptive surface minimizer, the pressure mixture
can be chosen among the attaining tests that are sharp
at any one fixed tangent partition. In particular, if a
two-phase tangent line carries a collapsed third-phase
layer, a common pressure law can be chosen whose tests
place the third score between the other two scores at
the base point. Its averaged triangle deficit has zero
value and zero gradient there.

The selection follows by applying an ambient variation
of size \(h t\) after replacing the cells on the tangent
scale \(t\). Sending \(t\) to zero first and then \(h\) to
zero restricts the pressure tests to the equality set of
the tangent comparison.
[Exact first-order tangent patch and equality tests (TC1–TC2)](THREE_CELL_ADAPTIVE_SURFACE_TANGENT_CLASSIFICATION.md#the-comparison-at-a-tangent-scale)
[The pressure mixture for the original minimizer (ASV1–ASV2)](THREE_CELL_ADAPTIVE_SURFACE_STATIONARITY.md#the-envelope-first-variation)

## Combining a tangent patch and a volume-tangent flow

Fix a balanced planar minimizer \(A\) at \(0<r<1\),
and let \(\mathcal M\) be its compact set of normalized
attaining tests. For a fixed base point \(p\), take
a sequence \(t_n\downarrow0\) along which the rescaled
cells and tagged arclength measures converge to \(C\)
and \(\eta_e\), respectively. Choose a continuity radius
\(R\) for the limiting measures. As in (TC1), write
\(w_e(G)=|T_{\sqrt r}G_i(p)-T_{\sqrt r}G_j(p)|\) and
\[
 d_R(G)=\sum_e w_e(G)
                [\eta_e(B_R)-P_e(C;B_R)]\ge0,\qquad
 \mathcal M_R=\{G\in\mathcal M:d_R(G)=0\}.
\]
Nonnegativity comes from (TC2), and \(\mathcal M_R\)
is compact and nonempty. Let \(X\) be a compact smooth
field with zero cell mass first variations \(M_i(X)\).
Fix a real \(h>0\).

Use the replacement of the original cells inside the
shrinking ball by the conical cells \(C\) in the proof
of (TC1), with radius \(R_n\to R\) and seam length
\(o(t_n)\). At this stage leave the \(O(t_n^2)\) volume
errors uncorrected. Call the resulting partition
\(B_n\). Apply the ambient flow \(\Phi_{h t_n}\) of
\(X\) to \(B_n\), then correct the masses exactly
with the two exterior flows used in (TC1). Denote
the final balanced partition by \(D_n(h,X)\).

The final correction times are \(O_{h,X}(t_n^2)\).
Indeed \(\Phi_{h t_n}A\) has a mass error of this
order because \(M(X)=0\). The image of the local
change from \(A\) to \(B_n\) still has Gaussian
volume \(O(t_n^2)\), by the uniformly bounded
Jacobian on the fixed compact neighborhood. The
patch and its ambient image avoid the support of
the exterior correction fields for large \(n\).
On those supports the mass response is that of
\(\Phi_{h t_n}A\), so the same invertible volume
derivative matrix gives the stated estimate by the
implicit function theorem.

For every unit-ball dual test \(G\), the surface
contribution of the final correction is
\(O_{h,X}(t_n^2)\). Before correction, compare the
flowed patched partition with the flowed original
partition. The deleted and inserted measures are
supported on the patch and have total length
\(O(t_n)\). The flow changes the weighted Jacobian
on this neighborhood by \(O_{h,X}(t_n)\), and the
test weights change by \(O_{h,X}(t_n)\), uniformly
over the unit-ball tests. Thus the difference of
their weighted perimeters has the same
first-order coefficient as the unflowed patch.
That coefficient at a fixed base attaining test
is \(-\phi_2(p)d_R(G)\), by weak tagged
convergence and the seam estimate.

The change from \(A\) to \(D_n(h,X)\) has the
following exact first-order formula:
\[
 \lim_{n\to\infty}
   \frac{\mathcal F_r(D_n(h,X))-\mathcal F_r(A)}{t_n}
    =\max_{G\in\mathcal M}
          [\,h L_G(X)-\phi_2(p)d_R(G)\,].                  \tag{TPR1}
\]
For the lower bound evaluate each fixed attaining
test of \(A\) on the final partition. The ordinary
flow contribution is \(h L_G(X)\) by (ASV1), and
the patch and repair contributions were computed
above.

For the upper bound take an attaining test \(G_n\)
at the final partition. The compact-source norm
estimate for changes of size \(O(t_n)\) implies
that any subsequential limit of these tests in
norm belongs to \(\mathcal M\), and that a
convergent subsequence exists. To verify this,
start from a tensor attaining the final norm.
Transport its directions back on the unchanged
interfaces under the two flows and use zero
directions on the removed patch of \(A\).
The surface edits and the compact ambient
transport change the resulting tensor by
\(O_{h,X}(t_n)\) in \(L^2\). Use compactness of
the reference tensors and the minimizing
inequality, just as in (TC1).

Along a subsequence \(G_n\to G\in\mathcal M\),
the ordinary flow difference quotient for those
tests converges to \(h L_G(X)\). Indeed their
smoothed fields converge in \(C^1\) on the compact
transport support. The limiting pair differences
are nonzero at almost every point of the
corresponding tagged measure by (ASV4).
The difference quotients of their Euclidean norms
therefore converge to the derivative of the
limiting norm at those points. The common
Lipschitz bound gives dominated convergence,
including the tangential Gaussian Jacobian term.
The patch coefficient also converges to
\(-\phi_2(p)d_R(G)\). Using
\(W_{G_n}(A)\le\mathcal F_r(A)\) for the
upper bound proves (TPR1).

## A pressure law supported on equality tests

Minimality makes the right side of (TPR1)
nonnegative for every \(h>0\). If \(G_h\) attains
that maximum, then \(d_R(G_h)=O_X(h)\), since
\(L_G(X)\) is bounded on \(\mathcal M\).
Also \(L_{G_h}(X)\ge0\). By compactness and
continuity, sending \(h\downarrow0\) therefore
gives
\[
     \max_{G\in\mathcal M_R}L_G(X)\ge0
                         \qquad\text{whenever }M(X)=0.
                                                               \tag{TPR2}
\]
For any finite list of such fields, separation
in the space of their evaluations now puts zero
in the convex hull of the evaluation vectors of
tests in \(\mathcal M_R\). Compactness of the
probability measures on \(\mathcal M_R\) gives
one law annihilating these variations for every
field with \(M(X)=0\). Surjectivity of the volume
derivatives then gives its centered pressure
vector, exactly as in (ASV2).

The law can be chosen sharp on all continuity
radii of this fixed tangent subsequence. Indeed
the defect in (TC2) is a nonnegative measure for
each test, so \(d_R(G)\) is nondecreasing with
\(R\). The equality sets \(\mathcal M_R\) are
therefore nested. Apply the finite intersection
argument simultaneously to the support conditions
on these sets and to the annihilation of the
balanced fields. A finite list of the support
conditions only requires the largest radius in
that list. This selection concerns one fixed
tangent subsequence.

## The pressure metric of a collapsed line

Suppose the tangent is the line between labels
\(i,j\) and has a nonzero excess measure of the
form (TC4), with \(k\) the absent label.
Every test in the selected law satisfies
\[
 T_{\sqrt r}G_k(p)
   =(1-\lambda_G)T_{\sqrt r}G_i(p)
        +\lambda_G T_{\sqrt r}G_j(p),
               \qquad 0<\lambda_G<1 .
\]
This is the equality case of the triangle
inequality from (TC4); its strict interior
condition follows from positivity on every
tag present in the excess. The selected tests
form a compact set of positive weights at \(p\).
Consequently their Euclidean norms are smooth
on a common neighborhood, with locally uniform
derivative bounds. Thus the averaged metric
satisfies
\[
 \begin{split}
 (\bar w_{ik}+\bar w_{kj}-\bar w_{ij})(p)&=0,\\
 \nabla(\bar w_{ik}+\bar w_{kj}-\bar w_{ij})(p)&=0.
 \end{split}                                                  \tag{TPR3}
\]
The second equality follows from nonnegativity
of the triangle deficit at every observation.

Its second derivative is also explicit. For
one selected test write \(F=T_{\sqrt r}(G_j-G_i)\),
\(E=T_{\sqrt r}(G_k-G_i)\), and let \(P_F^\perp\)
be orthogonal projection perpendicular to
\(F(p)\). In any planar observation direction
\(v\), the second derivative of that test's
triangle deficit at \(p\) is
\[
 \left.D_v^2\big(|E|+|F-E|-|F|\big)\right|_p
  =\frac{\left|P_F^\perp
            (D_vE(p)-\lambda_G D_vF(p))\right|^2}
         {\lambda_G(1-\lambda_G)|F(p)|}.                   \tag{TPR4}
\]
Indeed the terms pairing the second derivatives
of the vectors with their common unit direction
cancel. The remaining Hessians of the Euclidean
norm give the displayed square. The averaged
deficit in (TPR3) has Hessian obtained by
integrating this formula. The selection and
its second-order identity do not by themselves
exclude a collapsed layer.

## A second-order condition for a smooth touching layer

Assume that near \(p\), in tangential and normal
coordinates \((x,y)\), the two sheets of the layer
are \(C^2\) graphs \(y=h_-(x)\) and \(y=h_+(x)\).
The phase \(i\) lies below, \(j\) above and \(k\)
between them. Assume \(h_-<h_+\) for all
sufficiently small nonzero \(x\), with both values
and first derivatives zero at \(x=0\). Write
\(\kappa_-=h_-''(0)\), \(\kappa_+=h_+''(0)\) and
\(\Delta\kappa=\kappa_+-\kappa_-\).
For a sharp test \(G\), let
\(\Lambda_G(x,y)=w^G_{ik}(x,y)+w^G_{kj}(x,y)-w^G_{ij}(x,y)\).
All the weights without arguments below are
evaluated at \(p\). For every fixed \(0\le\gamma\le1\)
there is a sharp attaining test \(G\) such that
\[
 \partial_{xx}\Lambda_G(p)
    \ \le\ (\Delta\kappa)^2
           [\,\gamma^2 w^G_{ik}+(1-\gamma)^2 w^G_{kj}\,].
                                                               \tag{TPR5}
\]
In particular, taking \(\gamma=1/2\) gives
\(\partial_{xx}\Lambda_G(p)/w^G_{ij}\le
(\Delta\kappa)^2/4\). With equal curvatures, a
sharp test must have vanishing tangential
transverse separation in (TPR4).

To prove this, fix \(0<c<1/4\). Choose an even
Lipschitz function \(k_c\) on \([-1,1]\), zero
near the endpoints, with
\[
 k_c(s)=-s^2\ \ (|s|\le c),\qquad
 k_c(s)>-s^2\ \ (c<|s|\le1),\qquad
 \int_{-1}^1 k_c(s)\,ds=0.
\]
It can be chosen so that its derivative energy
outside \([-c,c]\) is \(O(c^4)\). For example,
join \(-c^2\) at \(c\) to zero at \(1/3\) by a
linear segment, reflect evenly, and add an even
nonnegative bump supported in
\(1/2<|s|<3/4\) with coefficient \(O(c^2)\) to
set the mean to zero.
Set \(\chi_c(s)=1+k_c(s)/s^2\) for \(s\ne0\)
and \(\chi_c(0)=0\). This is an even Lipschitz
cutoff for fixed \(c\), nonnegative and equal to
zero on \([-c,c]\), and equal to one near the
outer endpoints.

Write \(\Delta h=h_+-h_-\) and
\(h_0=h_-+\gamma\Delta h\). In \(|x|\le t\)
replace the two graphs by
\[
 \widetilde h_-=h_0-\gamma\chi_c(x/t)\Delta h,\qquad
 \widetilde h_+=h_0+(1-\gamma)\chi_c(x/t)\Delta h .
\]
The sheets merge in \(|x|\le ct\), giving a
single \(i,j\) wall there, and agree with the
original sheets near \(|x|=t\). In rescaled
coordinates \(x=ts\), the leading height
changes divided by \(t^2\) are
\(-\gamma\Delta\kappa\,k_c(s)/2\) and
\((1-\gamma)\Delta\kappa\,k_c(s)/2\),
respectively, uniformly in \(s\).
The corresponding slopes divided by \(t\)
converge to the derivatives of these profiles
almost everywhere, with uniform bounds for
fixed \(c\). This follows from the \(C^2\)
expansions of the original sheets.

Each leading height change has zero integral.
The exact Gaussian volume changes are thus
\(o(t^3)\): integrate the height transfers with
the Gaussian density, whose leading factor
is \(\phi_2(p)\). The exterior flows correct
the volumes with cost \(o(t^3)\) in the weighted
perimeter of every unit-ball test. Denote the
resulting balanced partition by \(A_t\).

Take an attaining test \(G_t\) for \(A_t\).
The perimeter of the changed part is \(O(t)\),
so the compact transport estimate from (TPR1)
gives subsequential norm convergence to some
base attaining test \(G\). At order \(t\) the
weighted perimeter change along this subsequence
is \(-2c t\,\phi_2(p)\Lambda_G(p)+o(t)\).
Indeed only the merged interval changes the
pair types to leading order. Minimality forces
\(\Lambda_G(p)=0\). The two original sheets have
positive limiting measure at \(p\), so (AT1)
gives both their frozen weights positive.
The merged pair weight is then positive too.
The weights along the subsequence converge in
\(C^2\) on a fixed neighborhood of \(p\).

For the order \(t^3\) comparison, first count
both new sheets over the whole interval as
virtual sheets with their old pair tags, even
where they coincide. The changes of the linear
normal-position terms in the two weighted
perimeters are \(o(t^3)\) because the leading
height changes have zero integral. The changes
in squared slopes have zero cross term to
leading order: the original scaled slopes are
\(\kappa_-s,\kappa_+s\) and
\(\int_{-1}^1 s k_c'(s)\,ds=-\int_{-1}^1 k_c(s)\,ds=0\).
Since the arclength expansion is
\(\sqrt{1+y'^2}=1+y'^2/2+O(t^4)\)
for these slopes, the virtual comparison
contributes
\[
 \phi_2(p)t^3\frac{(\Delta\kappa)^2}{8}
        [\,\gamma^2 w^G_{ik}+(1-\gamma)^2 w^G_{kj}\,]
                    \int_{-1}^1 k_c'(s)^2\,ds +o(t^3).
\]
The use of frozen weights in the squared slope
terms has error \(o(t^3)\) by their common
local continuity. The use of the first normal
derivative of the weighted Gaussian perimeter
density in the height terms has the same error
by \(C^1\) convergence and the zero integrals.

On the merged interval one must subtract
\(\Lambda_{G_t}(x,h_0(x))\) times the Gaussian
arclength factor of the common graph. This
factor is bounded below by
\(\phi_2(p)-O(t)>0\), and the deficit is
nonnegative. Taylor expand the deficit at
\(p\) before integration on the symmetric
interval. The constant term is nonnegative,
the first tangential term integrates to zero,
and the normal linear term is \(o(t^3)\)
since \(\partial_y\Lambda_{G_t}(p)\to0\)
and \(h_0=O(t^2)\). The quadratic leading
term is thus bounded below by
\[
 \frac{\phi_2(p)t^3}{2}
       \partial_{xx}\Lambda_G(p)\int_{-c}^c s^2\,ds +o(t^3).
\]
Together with \(W_{G_t}(A)\le\mathcal F_r(A)\)
and the exact volume correction, minimality
implies
\[
 \partial_{xx}\Lambda_G(p)
 \ \le\ \frac{(\Delta\kappa)^2}{4}
        [\,\gamma^2 w^G_{ik}+(1-\gamma)^2 w^G_{kj}\,]
       \frac{\int_{-1}^1 k_c'(s)^2\,ds}{\int_{-c}^c s^2\,ds}.
\]
The sharp tests form a compact set, and the
ratio of the two integrals tends to \(4\) as
\(c\downarrow0\), by the \(O(c^4)\) outer
derivative energy. Taking a subsequence of
these tests proves (TPR5).

If the ratio \(w_{ik}/(w_{ik}+w_{kj})\) is the
same across the sharp tests, \(\gamma\) can be
chosen to minimize the bracket in (TPR5).
Its minimum is \(w_{ik}w_{kj}/(w_{ik}+w_{kj})\).
For a test whose row differences are aligned
with the oriented sheet normals throughout
the punctured graphs, (TPR4) gives equality
in this optimized coefficient. Indeed the
two unit normals have derivatives along
the common tangent differing in magnitude
by \(|\Delta\kappa|\). Orthogonal projection
in (TPR4) cancels the derivatives of the
scalar weights and leaves
\[
 \partial_{xx}\Lambda_G(p)
        =\frac{w^G_{ik}w^G_{kj}}{w^G_{ik}+w^G_{kj}}\,
                      (\Delta\kappa)^2 .
\]
This equality uses alignment along the sheets,
not just collinearity at \(p\). The comparison
does not settle such an aligned layer.

## A shared pressure law along common scales

The selection works simultaneously at finitely many
distinct points \(p_j\) when their tangents and tagged
limits are taken along the same sequence \(t_n\).
Choose one continuity radius for each of these
tangents. The corresponding shrinking patches are
disjoint for large \(n\), and the seam slices can be
chosen separately within their respective annuli.
The two volume-correcting fields can be taken outside
fixed disjoint disks about the finitely many points.
Choose the disks small enough that each label still
has positive mass in their open complement. The
complement is connected, so the same volume
annihilator argument gives surjectivity there.

Perform all the patches, apply one macro flow with
velocity \(X\) and time \(h t_n\), and use a single
exterior mass correction of size \(O(t_n^2)\).
For a finite number of patches, the estimates
in (TPR1) are unchanged except for adding the
first-order patch contributions. The resulting
maximum is
\[
       \max_{G\in\mathcal M}
       \left[hL_G(X)-\sum_j \phi_2(p_j)d_{R_j,j}(G)\right],
                                                               \tag{TPR6}
\]
where \(d_{R_j,j}\) is the defect of the chosen
tangent at \(p_j\). Each summand of the defect is
nonnegative, and each Gaussian density is positive.
Since the maximum is nonnegative for every \(h>0\),
the limit as \(h\downarrow0\) selects tests that are
sharp at all the indicated points and radii. Applying
the separation argument of (TPR2) gives one stationary
pressure law supported on their common equality set.

For countably many distinct prescribed points one
may first choose a common diagonal subsequence of
scales along which all their cell and tagged limits
exist. This uses the local compactness at each
fixed point. Along such a chosen subsequence the
same finite intersection argument on probability
laws gives a single law sharp at all these tangents
and all their continuity radii. A finite set of
conditions only uses finitely many points, the
largest indicated radius at each, and finitely
many volume-tangent fields. The construction
requires common scales for the chosen tangents.
[Compactness of cells and tagged measures on tangent scales](THREE_CELL_ADAPTIVE_SURFACE_TANGENTS.md#homogeneity-and-nonconstancy-of-cell-limits)

## A limitation of local metric conditions

Constant positive tensions with triangle equality
already permit collapsed layers in the Euclidean
weighted perimeter problem. Take two disjoint open
disks tangent at \(p\), assign them labels \(i,j\),
and assign their exterior label \(k\). For any
partition with the two prescribed disk areas and
constant tensions satisfying \(w_{ij}=w_{ik}+w_{kj}\),
the weighted perimeter is exactly
\[
                       w_{ik}P(A_i)+w_{kj}P(A_j),
\]
where \(P\) denotes Euclidean perimeter. The two
tangent disks minimize this sum by the Euclidean
isoperimetric inequality for each disk separately.
Their tangent at \(p\) is a line between \(i,j\),
and both \(i,k\) and \(k,j\) interfaces converge to
that line. The triangle deficit and its
derivatives vanish identically. Thus even
minimization for a single positive degenerate
metric can allow a collapsed line. Excluding
collapse for the adaptive Gaussian problem calls
for an additional restriction on its Mehler
tests or their interaction.
