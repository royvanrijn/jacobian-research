# Finite transfers and motion at a three-cell contact

Let \(A\) be a balanced planar Gaussian maximizer at \(0<r<1\).
Use the signed pair incidence vectors, the reference Gaussian
surface measure, and the quadratic form \(Q_0\) from the
[maximizer analysis](THREE_CELL_GLOBAL_OBSTRUCTIONS.md).
A tied maximizer has the stability of the regular simplex \(T\)
at \(r\). The first two results below constrain how such a
competitor can improve when the correlation or its labels change.
The last concerns the creation of a component of one label
inside a different cell of an actual maximizer.

## Motion from a tangential endpoint tie

For the fixed partition \(A\) put \(D(s)=S_s(T)-S_s(A)\).
Suppose \(D(r)=D'(r)=0\) and \(D''(r)>0\). Let \(A_h\) be any
balanced partition beating \(T\) at \(r+h\), for sufficiently
small nonzero \(h\), in the same plane and label coordinates.
If \(E_h=\{A_h\ne A\}\), then
\[
 \int_{E_h}(1+|x|)\,d\gamma_2\ \ge c|h|,\qquad
 \gamma_2(E_h)\ \ge
         c'\frac{|h|}{\sqrt{\log(e/|h|)}}.                   \tag{FC1}
\]
The constants can be uniform over nonsimplex reflected
tangential endpoint ties in a compact positive correlation band,
by their uniform positive curvature bound.
[Reflected contact curvature](THREE_CELL_FIRST_CONTACT_RIGIDITY_TESTS.md#reflected-ties-have-at-most-second-order-contact)

For any indicator \(f\), differentiating its conditional Gaussian
density gives \(|\partial_s T_s f(x)|\le C(1+|x|)\) on a
compact correlation band. The covariance and mean terms in the
derivative are polynomials of degrees two and one in a standard
conditional Gaussian, with the mean term linear in \(x\).
By self-adjointness and polarization, the difference of the
two fixed-partition stability derivatives is thus bounded by
\(C\int_{E_h}(1+|x|)d\gamma_2\), uniformly in \(s\) between
\(r\) and \(r+h\). At \(r\) the difference
\(S_r(A_h)-S_r(A)\) is nonpositive. At \(r+h\) it exceeds
\(D(r+h)\ge c_0h^2\). Integration proves the first bound.
For a measurable set of mass \(\delta\), splitting at
\(R=\sqrt{2\log(1/\delta)}\) and evaluating the radial
Gaussian tail gives
\(\int_E(1+|x|)d\gamma_2\le C\delta\sqrt{\log(e/\delta)}\).
The function on the right is increasing for \(0<\delta\le1\).
Substitution of a sufficiently small constant times
\(|h|/\sqrt{\log(e/|h|)}\) into it proves the second bound.
Uniform Taylor remainders follow from the bounded third
correlation multipliers on a compact positive band.

There is an exact equality requirement for a smoothly moving
improving competitor. Suppose \(r\) is approached by distinct
simplex-optimal correlations, \(D(r)=D'(r)=0\), and a balanced
improving path has first-order normal speed \(u\) per unit
change of correlation at \(r\). Assume the joint Taylor
expansion below is valid for the path and for balanced compact
test speeds. This includes compact smooth deformations with
exact volume correction. Write \(g\) for the mixed pair-score
derivative in (G10). Then
\[
 Q_0(u,u)=D''(r)/2,\qquad
 Q_0(u,f)=\langle g,f\rangle
     \quad\hbox{for every balanced compact regular }f.       \tag{FC2}
\]
Here in addition the expansion and the nonnegativity test
are assumed applicable to \(u+\lambda f\) for real
\(\lambda\) and all balanced compact regular \(f\).
Assume also that \(u\) lies in the closure of the balanced
compact speeds in a norm in which \(Q_0\) and the mixed
pairing are continuous. For compact smooth \(u\) the
perturbed directions are just other compact smooth
balanced speeds and the closure condition is automatic.

For a signed correlation increment \(h\) and speed \(f\)
multiplied by \(h\), the deficit divided by \(h^2\) tends to
\[
 D''(r)/2-2\langle g,f\rangle+Q_0(f,f).
\]
This is the same Taylor expansion as in (G10), with the
deformation scaled by the signed increment. Along the
simplex-optimal sequence the expression is nonnegative for
all such speeds. Along the improving path its value for
\(f=u\) is nonpositive. Hence \(u\) saturates the quadratic
inequality and minimizes the expression. Variation in
\(\lambda\) gives the second identity of (FC2), and its
specialization to \(u\) by the stated continuity gives
the first.

These are necessary conditions at second order. For example,
the smooth two-variable deficit model
\((t-h)^2-h^3\) has positive fixed-shape curvature at \(t=0\),
nonnegative left-side optimum, and an improving right-side
path \(t=h\) saturating its mixed quadratic bound. With a
degenerate shape Hessian the model \(t^4-ht^2+c h^2\) for
\(0<c<1/4\) instead has an improving path of scale
\(\sqrt{h}\). The models concern the variational argument,
with no assertion that they are Gaussian partitions.

## Balanced cuts using only two labels

At a nonsimplex tie \(A,T\), put \(e=1_A-1_T\) in vector
notation. Let \(E\) be a measurable cut with
\(\int 1_Ee\,d\gamma_2=0\), and suppose the transfers on
\(E\) use a single unordered pair. Let \(B=T+1_Ee\) denote
the patched partition in indicator notation. Uniformly on
a compact positive correlation band,
\[
 \sum_i\langle 1_Ee_i,T_r(1_{E^c}e_i)\rangle
   \ \ge\ c L_r(B;T)\ \ge\ c'\gamma_2(B\ne T)^2             \tag{FC3}
\]
with positive constants independent of the tie, the orientation
of \(T\) and the cut.

For small \(\gamma_2(B\ne T)\) use the compact flux proof
of (B11) in the fixed orientation \(T\). The only possible
limiting signed flux is supported on the facet of the one
pair. This is a closed subspace of balanced speeds. It
meets the rotation kernel trivially, since a nonzero
planar rotation has nonzero speed on each tripod ray.
Its angle from the finite-dimensional kernel is positive.
Thus the interface form has a uniform positive gap on this
subspace, by (S37) and compactness of the band. The proof
of (B13)–(B14) applies there without the rotation gauge:
the same transfer estimates give the flux, its regret
bound and its interaction limit, and the subspace gap
gives the fixed-orientation score-relative comparison.
Consequently \(S_r(T)-S_r(B)\ge c_1 L_r(B;T)\) for small
label displacement. The identity (G12) expresses twice the
cross energy in (FC3) as this deficit plus that of the
complementary patch. The latter is nonnegative at the tie.
The saturated margin implies the last inequality of (FC3).
[Measurable local comparison](POSITIVE_CORRELATION_CONTINUATION.md#an-interface-gap-implies-a-uniform-measurable-local-comparison)

For cuts of displacement bounded below, if the cross energy
tended to zero the two patched partitions would both become
maximizing. Fix the orientation and pair by rotations and
passing to a subsequence. On the fixed active plane and at a
positive limiting correlation, near-maximum compactness gives
strongly convergent subsequences of \(A\) and \(B\), with
maximizing limits. The limit of \(A\) is nonsimplex by uniform
local isolation of the simplex. The limit of \(B\) preserves
the untouched simplex cell. By (G16) it is the simplex with
either the original labels or a transposition of the two
changed labels. The first case contradicts the displacement
bound. In the second case \(B\) differs from \(T\) throughout
the two transposed cells. The patched label comes from either
\(A\) or \(T\), a property preserved by subsequential
almost-everywhere convergence of these discrete labels.
Thus \(A\) has the transposed labels on those two full
cells, and balance forces its third cell as well. This
contradicts the nonsimplex limit. The cross energy consequently
has a positive uniform lower bound for cuts of displacement
bounded below. Since \(L_r(B;T)\le1\), this finishes (FC3).
[Single simplex cell rigidity and balanced patch identity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md)

One can always orient \(T\) so that all transfers decompose
into at most three balanced pair cuts, with at least one third
of the label mass unchanged. For a rotation by \(\theta\)
write \(p_{ij}(\theta)=\gamma_2(T_i^\theta\cap A_j)\).
Equal row and column sums imply
\(p_{12}-p_{21}=p_{23}-p_{32}=p_{31}-p_{13}\).
With \(\omega=e^{2\pi i/3}\), set
\(z(\theta)=\sum_{i,j}p_{ij}(\theta)\omega^{j-i}\).
This is continuous. Its imaginary part is zero exactly when
the common antisymmetric difference vanishes, and its real
part is \((3\sum_i p_{ii}-1)/2\). A rotation through \(2\pi/3\)
multiplies it by a fixed primitive cube root. If it vanishes
at some rotation, use that rotation. Otherwise its winding
number over \(2\pi\) is congruent to \(1\) or \(-1\) modulo
three, so it visits the positive real axis. The corresponding
matrix is symmetric with the stated unchanged label mass.
For a nonsimplex tie, at least two pairs in this orientation
are active by single-cell rigidity. Inequality (FC3) then
requires each active pair to have strictly positive cross
energy with the sum of the others.

More precisely let \(f_e=1_{\{T_i\to A_j\}}-
1_{\{T_j\to A_i\}}\) for each cyclically oriented pair
\(e=(i,j)\), and put \(F=\sum_e f_e\). Write
\(q_r(f)=\langle f,T_rf\rangle\) for scalar interaction.
The cyclic label-incidence products give the identity
\[
 Q_r(e)=3\sum_{\text{pairs }l}q_r(f_l)-q_r(F) ,
\]
where on the left \(e=1_A-1_T\) is the full vector
transfer. Summing (FC3) over the pair cuts gives
\(\sum_l q_r(f_l)-q_r(F)\ge c L_r(A;T)\).
Using \(Q_r(e)=2L_r(A;T)\) at the tie first gives
\(L_r(A;T)\ge\sum_l q_r(f_l)\), and consequently
\[
                  q_r(F)\le (1-c)\sum_l q_r(f_l).          \tag{FC7}
\]
Thus the signed sum of forward and reverse transfers in
this orientation must have a strict smoothing deficit
relative to the separate pair fields.

Positivity of every balanced patch test in this orientation
can occur at strictly suboptimal partitions. To see this,
let \(v_i\) be the three sector-center unit vectors and
choose two disjoint subsets of equal small Gaussian mass
\(p\) in each ball \(B_1(Rv_i)\), one for each output
label \(j\ne i\). Take \(R\) large enough that the balls
lie in their respective reference sectors. Transfer each
chosen set from \(i\) to its assigned \(j\); all masses
are preserved and the overlap matrix is symmetric.
For an arbitrary balanced cut let \(u_i\in[0,2]\) be the
fractional mass it selects from the two source-\(i\) sets,
in units of \(p\), and \(U=\sum_i u_i\). Balance gives
\(u_i\le U/2\), because its incoming flux must come from
the other two labels. Swap the cut with its complement if
needed so that \(0<U\le3\). Then
\[
 \sum_i u_i(2-u_i)\ge U/2,\qquad
 \sum_{i\ne j}u_i(2-u_j)\le4U .
\]
Incidence dot products for transfers from the same source
are at least one, and other dot products have absolute
value at most two. Let \(m\) be the minimum Mehler kernel
on the same-source balls and \(m'\) the maximum between
different balls. Here \(m'/m\to0\) as \(R\to\infty\):
the quadratic cross exponent loses
\(3rR^2/[2(1-r^2)]\) and moving either argument inside
its unit ball changes it by \(O(R)\). Hence the cross
energy is at least
\(p^2(m-16m')\sum_i u_i(2-u_i)>0\).
Taking \(p\) small makes the constructed partition strictly
suboptimal by the local simplex gap and comparison. The
sign tests by themselves therefore permit a symmetric
overlap with all three positive pair couplings.

There is also an exclusion for near-tie competitors.
At any \(r_0\) where \(T\) is globally optimal, there are
\(\varepsilon,\eta,c>0\) such that for \(|r-r_0|<\eta\)
any balanced planar partition \(A\) with one cell within
\(\varepsilon\) in Gaussian symmetric difference of a
central \(120^\circ\) wedge satisfies
\[
 S_r(T)-S_r(A)\ge c\,d_3(A)^2 ,
\]
where \(d_3\) denotes distance to the simplex orbit.
Use the uniform local comparison for small \(d_3\).
Otherwise, a sequence violating the claim as the band
width, cell distance and coercivity constant tend
to zero has near-maximal value at \(r_0\). Rotate and relabel
its distinguished cell to fix the reference wedge. Strong
maximum compactness at \(r_0\) gives a limit maximizing
partition with that cell, hence a simplex by (G16).
This contradicts the uniform lower bound on \(d_3\) outside
the local comparison neighborhood.

Finally, assume a nonsimplex \(A\) ties \(T\) at a correlation
approached by simplex-optimal correlations and agrees with
\(T\) outside a bounded set. If \(D'(r)=0\) then \(D''(r)>0\),
with no reflection or connectedness hypothesis.
If the curvature were zero, (G10) would make the linear
functional \(\langle g,f\rangle\) vanish on the kernel
of the mass flux. Flux is onto the centered label space
using two regular bumps on a spanning tree, so \(g\) on
all regular interfaces is given by constant label-price
differences. On each of the three shared exterior
reference rays, the \(T\) correlation derivative of the
pair contrast is zero by reflection. Its difference
for \(A\) is the smoothed derivative of a bounded-support
function, which decays Gaussianly on that ray by the
conditional Gaussian derivative formula. Hence all three
price differences vanish. The dilation identity (G6)
and positivity of the score gradient then force the
dilation speed to vanish. Centered conical rigidity
gives the simplex. Nonnegative curvature at the original
tangency follows by Taylor expansion on the approaching
optimal correlations.

For the exterior-agreement hypothesis, a geometric
restriction also holds at any maximizer, without derivative
conditions. Each adjusted pair contrast vanishes on the
whole supporting line of the matching simplex facet, by
analytic continuation from the outside ray. The price
differences are zero, since the contrast of the reference
is zero on the line and the score of the bounded perturbation
decays along the outside ray. No strict winner of label \(i\)
can cross either of its two reference lines, where it ties
another label. An unbounded strict-winner component of label
\(i\) meets the exterior \(T_i\) and is consequently
confined to that sector. Thus every changed point, up to
null sets, lies in a bounded strict-winner component.

The origin in this situation is a genuine nondegenerate
triple junction. All scores coincide there. The gradients
of the three cyclic contrasts lie along the three cyclic
pair normals and sum to zero, so they are a common multiple
of those normals. If the multiple were zero, (G20) and
local stratification would preclude actual jump arcs
approaching the origin. A whole ball would then lie in
one label, with the other scores still tied at its center,
which is excluded by (FC5) below. The multiple is therefore
nonzero. The stationary zero lines and the implicit
function theorem then identify the partition locally as
the reference simplex or its spatial reflection through
the origin.

## Point transfers at a dormant interior score tie

Suppose a ball about \(p\) lies in cell \(i\) up to null sets
and the adjusted \(i,j\) scores tie at \(p\), where \(j\ne i\).
Call this a dormant interior tie. Write
\(H=\nabla^2[U_i-U_j](p)\), \(U_i=T_r1_{A_i}\).
Then \(H\) is positive definite and
\[
 \frac{\sqrt{\det H}}{2\pi\phi_2(p)}
  \ \ge\ 2M_r(p,p)+(r^{-1}-1)\|s_v\|_a^2,                 \tag{FC4}
\]
where \(\phi_2\) is standard planar Gaussian density and
\(v\) is the unique velocity whose first-order volume change
compensates a unit transfer from \(i\) to \(j\). In particular
no dormant interior tie lies in
\[
 |p|^2\le\frac{2(1+r)}{1-r}\bigl(1-2\log r\bigr).          \tag{FC5}
\]

The adjusted \(i,j\) contrast \(W\) is nonnegative in the
ball and vanishes at \(p\). Thus \(\nabla W(p)=0\) and
\(H\succeq0\). Insert an ellipse \(E_t=p+P B_{R_t}\) of
Gaussian mass \(t\) into label \(j\), with \(\det P=1\).
Translate all original cells by \(tv+O(t^2)\) before this
insertion so that the two changes preserve the masses exactly.
The full moment rank makes the derivative of the translation
volume map invertible onto the centered label space, giving
this translation by the implicit function theorem. For small
\(t\) the shrinking ellipse still lies in label \(i\) of the
translated cells.

Use the augmented objective \(S_r-2\sum_i\beta_i\gamma(A_i)\)
with the stationary prices. Its translation gain is
\((r^{-1}-1)\|s_v\|_a^2 t^2+o(t^2)\) by the translation
identity (G2) and the constrained second-variation formula
for the augmented objective. The first derivative vanishes
by stationarity, even though the translation itself changes
volumes. The linear cost of the insertion is \(2\int_{E_t}W\)
up to \(o(t^2)\): translating by \(w\) changes the score at
\(x\) to the original score at \(x-w/r\), and \(\nabla W(p)=0\).
Taylor expansion on this ellipse gives
\[
 \lim_{t\downarrow0}\frac{2\int_{E_t}W\,d\gamma_2}{t^2}
   =\frac{\operatorname{tr}(P^{\mathsf T}HP)}{4\pi\phi_2(p)}.
\]
Indeed \(t\sim\pi\phi_2(p)R_t^2\) and the uniform disk
second moments are \(R_t^2/4\) per coordinate. The quadratic
self-interaction of the insertion divided by \(t^2\) tends to
\(2M_r(p,p)\), by continuity of Mehler and the squared incidence
norm two. Maximality makes the total change nonpositive.
Minimizing the trace over \(P\) gives \(2\sqrt{\det H}\),
including an infimum of zero for singular \(H\). This gives
(FC4) and strict positive definiteness.

For the explicit exclusion, differentiation of the conditional
Gaussian gives
\(\operatorname{tr}H\le (r^2/(1-r^2))\mathbb E||Z|^2-2|
=4r^2/[e(1-r^2)]\) for a standard planar \(Z\). The
last expectation follows by integrating the centered
exponential variable \(|Z|^2\) of mean two. Also
\(M_r(p,p)=(1-r^2)^{-1}e^{r|p|^2/(1+r)}\).
The nonzero balancing translation has strictly positive
\(\|s_v\|_a^2\) by the full active rank. Combining these
identities with \(\sqrt{\det H}\le\operatorname{tr}H/2\)
gives (FC5). At a dormant tie outside this disk, the
fixed-shape correlation derivative of \(U_i-U_j\) at \(p\)
is \(-\operatorname{tr}H/r<0\) by the OU generator identity
and the vanishing gradient.

For each fixed ordered pair, the dormant interior points
are uniformly separated on compact sets and compact
correlation bands. Indeed (FC4) gives a uniform positive
determinant bound on the compact set. Gaussian
differentiation bounds the second and third derivatives
of the contrasts uniformly on the band. The Hessian at
each point therefore has a positive lower eigenvalue
bound and remains positive definite in a ball of
uniform radius around it. A function with positive
definite Hessian on a ball can have only one critical
point there, by strict monotonicity of its gradient.

The same argument applies simultaneously at finitely many
distinct dormant points \(p_l\). Let \(d_l=e_{j_l}-e_{i_l}\)
be the change of volume label per unit inserted mass and
\(H_l\) the corresponding Hessian. For arbitrary nonnegative
amplitudes \(u_l\) choose the translation velocity \(v\)
compensating \(\sum_l u_l d_l\). Inserting disjoint ellipses
of masses \(tu_l\), and minimizing their individual trace
costs, gives the inequality
\[
 \sum_l \frac{\sqrt{\det H_l}}{2\pi\phi_2(p_l)}u_l^2
 \ \ge\
 (r^{-1}-1)\|s_v\|_a^2+
 \sum_{l,m}u_l u_m(d_l\cdot d_m)M_r(p_l,p_m).              \tag{FC6}
\]
The cross interaction converges to the displayed Mehler
matrix by continuity on the shrinking ellipses. The linear
translation error is still \(o(t^2)\) at every point because
all the score gradients vanish there. In particular (FC6)
is a restriction on simultaneous point insertions in addition
to the individual interior exclusion (FC5).
