# Moving maximizers and the junction scale near correlation one

Balanced three-cell maximizers converge in Gaussian measure to the simplex
orbit as $\rho\uparrow1$, with a dimension-independent conclusion after
dimension reduction. Their rescaled position-and-direction measures also
converge to the simplex interface measures. The corresponding score
prices are smaller than the noise length. At the noise length, every
three-label minimizing limit is a Euclidean tripod. Actual maximizers
have a unique triple score tie and uniformly aligned interface graphs
on bounded physical sets. The open-ray graph comparison and relative
exterior estimate give additional control of shrinking deformations.

## Compactness and the sharp perimeter lower bound

The Gaussian perimeter comparison uses Milman–Neeman, Theorems 1.1
and 1.2, for positive prescribed masses and $k\le d+1$. Here the active
dimension is two and there are three cells. Corollary 1.11 and Section
12 of Heilman's *Stable Gaussian Minimal Bubbles* provide a
fixed-competitor estimate with an error depending on the competitor.
Related Gaussian Gamma-convergence results concern fractional energies.
The following proof works directly at a single OU correlation on
varying partitions. [Milman–Neeman](https://arxiv.org/html/1805.10961v3),
[Heilman](https://link.springer.com/article/10.1007/s00526-026-03408-z),
[fractional Gaussian perimeter comparison, Theorem 1](https://arxiv.org/pdf/2011.07544).

Let $e=\sqrt{2(1-\rho)}$, $a=\sqrt{(1+\rho)/2}$, and let $w$ be the
midpoint of the two observations and $z$ their rescaled difference.
The joint law of $(w,z)$ is $\phi_a(w)dw\,d\gamma_d(z)$, where
$\phi_a$ is the density of $N(0,a^2I)$. Write $u_i=1_{A_i}$ and

\[
 d\mu_{e,i}(w,z)=\frac{|u_i(w+ez/2)-u_i(w-ez/2)|}{e}
                    \phi_a(w)dw\,d\gamma_d(z).
 \tag{P2}
\]

If $N_\rho(A)=1-S_\rho(A)$ is the probability of disagreement, then
$E_e(A):=N_\rho(A)/e=\frac12\sum_i\mu_{e,i}(\mathbb R^{2d})$.

Suppose $e_j\downarrow0$ and $E_{e_j}(A^{(j)})$ remains bounded in a
fixed dimension. On each fixed spatial ball and for $|z|\le1$ the
two densities in (P2) have a uniform positive lower bound. Translate
the midpoint to an endpoint to obtain, for every fixed $R$,

\[
 \int_{|z|\le1}\int_{B_R}|u_{j,i}(x+e_jz)-u_{j,i}(x)|dx,dz
                                  \le C_R e_j E_{e_j}(A^{(j)}).
\]

Mollify $u_{j,i}$ at scale $e_j$ by an even smooth unit-integral bump
supported in $B_1$. The $L^1(B_R)$ error is $O_R(e_j E_{e_j})$, while
the $L^1(B_R)$ gradient of the mollification is $O_R(E_{e_j})$. In the
second bound use the zero integral of each derivative of the bump to
subtract $u_{j,i}(x)$ before applying the displayed estimate. Local BV
compactness and a diagonal subsequence give local $L^1$ convergence to
indicators $u_i$ of an ordinary partition $A$. The Gaussian tails of
bounded functions upgrade it to $L^1(\gamma_d)$, preserving the masses.

Let $D u_i$ be the distributional gradient of the limiting indicator,
and write $D u_i=n_i\,|D u_i|$, so $n_i$ is the measure-theoretic
normal up to orientation. For fixed $z$, the signed difference
quotient in (P2), tested against any compact smooth function of $w$,
tends distributionally to $z\cdot D u_i$. To see this, shift the
smooth test to the two endpoints, take its difference quotient, and
use local $L^1$ convergence of $u_{j,i}$. The weights $\phi_{a_j}$ and
their derivatives converge uniformly on each compact set to $\phi_1$.
Lower semicontinuity of total variation against any nonnegative compact
continuous test in $w$ gives, for fixed $z$, the liminf with measure
$\phi_1(w)|z\cdot n_i(w)|d|D u_i|(w)$. More explicitly, the weighted
variation is the supremum over compact smooth signed tests bounded in
absolute value by the nonnegative continuous weight, with smooth
approximation from below on compact subsets where the weight is positive.
This argument applies to each $w$ slice of a nonnegative compact
continuous test in $(w,z)$. Fatou in $z$ then yields the joint measure
liminf with

\[
 d\mu_i(w,z)=\phi_1(w)|z\cdot n_i(w)|d|D u_i|(w)d\gamma_d(z).
\]

Taking increasing compact cutoffs gives

\[
 \liminf_j E_{e_j}(A^{(j)})\ge
      \sqrt{2/\pi}\,\mathcal P_\gamma(A),
 \qquad \mathcal P_\gamma(A)=\tfrac12\sum_i P_\gamma(A_i).
 \tag{P3}
\]

Here $\mathbb E|z\cdot n_i|=\sqrt{2/\pi}$ for unit $n_i$.

Apply this in $d=2$ to balanced sequences with
$N_{\rho_j}(A^{(j)})\le N_{\rho_j}(T)+o(e_j)$. The
[simplex expansion](THREE_CELL_DEFICIT_AND_PERIMETER.md#exact-correlation-derivative-for-the-three-cell-simplex)
gives $E_{e_j}(T)\to\sqrt{2/\pi}\,\mathcal P_\gamma(T)$.
Milman–Neeman therefore forces every subsequential limit to be a simplex.
Since the total mass in (P2) converges to the sum of the $\mu_i$ masses,
the individual measure liminfs imply narrow convergence $\mu_{e_j,i}\to
\mu_i$ after alignment. Indeed any vague cluster point dominates
$\mu_i$, and the total upper mass bound forces equality for each index.
This also gives tightness and convergence of the total mass in each
component.

In particular, for every fixed $\delta>0$ there is an $\eta(\delta)>0$
such that each planar balanced $A$ with $d_3(A)\ge\delta$ satisfies
$S_\rho(A)<S_\rho(T)$ for $1-\eta(\delta)<\rho<1$. Here $d_3$ is the
infimum of the total Gaussian symmetric differences from the simplex
orbit. A violating sequence would have an aligned convergent subsequence
as above. Applying this to global maximizers via Heilman–Tarter gives
convergence in any ambient dimension, with the same distance threshold
in all dimensions. The comparison for arbitrary fixed-distance
competitors is stated in the plane. This compactness result gives no
rate for $d_3(A^{(j)})$ in terms of $e_j$.

## Score prices are smaller than the noise length

Consider planar balanced maximizers $A^{(j)}$ at $\rho_j\uparrow1$,
aligned to $T$ by the preceding result. Take centered score prices
$\beta_{j,i}$ so that the winning scores are $T_{\rho_j}1_{A_i^{(j)}}-
\beta_{j,i}$. Such prices exist by the
[linearized assignment argument](MAXIMIZER_PERTURBATIVE_RIGIDITY.md#assignment-and-the-exact-information-supplied-by-global-maximality).
For a compact smooth vector field $V$, pull the joint density back by
its flow in both variables. With $x=w-ez/2$ and $y=w+ez/2$, the
first variation of $N_\rho/e$ is the integral of

\[
 D_e(V)=\operatorname{div}V(x)+\operatorname{div}V(y)
 -\frac{w\cdot(V(x)+V(y))}{2a^2}
 -\frac{z\cdot(V(y)-V(x))}{e}
\]

against $\frac12\sum_i\mu_{e,i}$. On compacts this tends uniformly to
$2\operatorname{div}V(w)-w\cdot V(w)-z^TDV(w)z$.
For the tails, if $D_e(V)$ is nonzero then at least one endpoint
lies in a fixed ball $B_R$. Hence $|w|\le R+e|z|/2$ and
$|D_e(V)|\le C_V(1+|z|^2)$. Compare the joint density at noise length
$e$ to the density at $2e$ in the original $(x,y)$ coordinates. Their
ratio on these pairs is bounded by $C_R\exp(-c|z|^2)$ for small $e$.
Indeed the rescaled difference exponent gives $-3|z|^2/8$, and the
midpoint exponent changes by $O(e^2)|w|^2$. If $\rho'=1-4(1-\rho)>0$,
then $N_{\rho'}(A)\le4N_\rho(A)$ for every partition. In its Hermite
expansion this is the inequality
$1-(\rho')^m\le4(1-\rho^m)$, obtained by writing the quotient as a
ratio of two geometric sums. These estimates bound the integral of
$|D_e(V)|$ on $|z|>L$ by a vanishing Gaussian tail times a uniform
bound for $E_e(A^{(j)})$. The midpoint coordinate is bounded for bounded
$z$ where $D_e(V)$ is nonzero. Thus the narrow convergence just proved
applies to the first variations.

At the limit the integral equals $\sqrt{2/\pi}\,\mathcal P_\gamma'(T)[V]$.
To check it, integrate the identity
$\mathbb E[|z\cdot n|zz^T]=\sqrt{2/\pi}(I+nn^T)$,
leaving $\operatorname{div}_{n^\perp}V-w\cdot V$ on each facet.
For the simplex rays the first variation is zero, by weighted
integration along each ray and the sum-zero junction tangents. Hence
$N'_{\rho_j}(A^{(j)})[V]=o(e_j)$ for every fixed compact smooth $V$.

On the other hand the linearized adjusted score is maximized over all
partitions, so its flow derivative is zero. It follows that
$N'_{\rho_j}[V]=-2\beta_j\cdot m'_j[V]$, where $m'_j[V]$ is the vector
of volume derivatives. Two fixed bumps in distinct open simplex rays
give independent volume columns at $T$. Their columns converge for
$A^{(j)}$ by the bulk formula $m'_{j,i}[V]=\int_{A_i^{(j)}}
(\operatorname{div}V-x\cdot V)d\gamma_2$. Invert these columns on the
centered label space. This proves

\[
                         |\beta_j|=o(e_j).
 \tag{P4}
\]

In particular there is no order-$e_j$ bias from the Gaussian mass
constraints in an interior near-one scaling limit of the stationarity
equations. The conclusion is a statement about actual maximizers.

## A compact Euclidean junction test

Let $Y$ be the $120$-degree partition of Euclidean $\mathbb R^2$.
For any measurable labeled partition $B$ equal to $Y$ outside a fixed
ball, define

\[
 \mathcal J(B)=\iint \phi_1(v-u)
       [\mathbf1_{b(u)\ne b(v)}-\mathbf1_{Y(u)\ne Y(v)}]du\,dv,
 \tag{P5}
\]

where here $\phi_1$ is the density of $N(0,I)$ in two dimensions.
Only pairs with at least one endpoint in the ball contribute, so the
integral is absolutely convergent. A negative value in (P5) implies
that the balanced Gaussian simplex fails to maximize stability for
all $\rho$ sufficiently close to one.

To prove this implication, replace the partition inside the ball
scaled by $e$ by the corresponding scaled copy of $B$. Changing
variables $(x,y)=(eu,ev)$ in the joint kernel yields an exact
factor $e^2\phi_a(e(u+v)/2)\phi_1(v-u)$ in dimension two. Thus the
disagreement change is $e^2\phi_1(0)\mathcal J(B)+o(e^2)$ by dominated
convergence. The mass changes are $O(e^2)$. Correct them by two fixed
smooth bump flows supported at positive distance from the origin on
independent simplex rays. Their volume derivative matrix is invertible
on the centered mass directions, so the correcting parameters are
$O(e^2)$ by the implicit function theorem. The inner change does
not affect those supports, and the correction achieves exact balance.

For the outer correction alone, the signed loss is at most twice its
score regret relative to $Y$, using positivity of $T_\rho$ in the exact
quadratic expansion. Every changed point is within $O(e^2)$ of a
regular facet, its score difference vanishes on that facet, and the
Gaussian convolution gradient has bound $O(e^{-1})$. The correction
therefore costs at most $O(e^3)$. The mixed inner/outer interaction is
exponentially small by fixed positive separation of their supports in
the joint kernel. Adding both changes gives a balanced competitor with
strictly lower disagreement whenever $\mathcal J(B)<0$.

**Small-support Euclidean junction comparison.** Put
$c_*=(1-\exp(-1))\exp(-25/32)/32$. If $B$ differs from $Y$ only in a disk of radius
$0<R<c_*\pi/24$, then $\mathcal J(B)>0$ unless it agrees almost
everywhere. The disk can have any center and the comparison has no mass
constraint.

For $x\in Y_i$, let $t$ be its inward distance to the pair line with
$Y_j$, and let $P_i=\phi_1*1_{Y_i}$. Reflection gives

\[
 P_i(x)-P_j(x)=\int_{Y_i-x}\phi_1(z)
                             (1-e^{-2t(t+z\cdot n_{ij})})dz.
\]

On the radius-$1/4$ ball about the unit bisector of $Y_i$, the quantity
$z\cdot n_{ij}$ is at least $1/2$. This ball lies in $Y_i-x$ by cone
addition and has standard Gaussian mass at least $\exp(-25/32)/32$.
Restricting the nonnegative integral proves
$P_i(x)-P_j(x)\ge c_*\min(t,1)$.
For a transfer $Y_i\cap B_j$ of area $m_{ij}$ in the disk, the strip
$0\le t<u$ has area at most $2Ru$. Choose $u=m_{ij}/(4R)\le1$,
using $m_{ij}\le\pi R^2$ and $R\le1$. The integral of $\min(t,1)$
over the transfer is at least $m_{ij}^2/(8R)$. Hence the score regret
$L$ is at least $c_* m^2/(48R)$, where $m$ is the total transferred
area. With $\delta_i=1_{B_i}-1_{Y_i}$, the interaction $Q$ is bounded
above by $\phi_1(0)\sum_i\|\delta_i\|_1^2\le m^2/\pi$, since each
norm is at most $m$ and their sum is $2m$. The exact finite-difference
identity $\mathcal J(B)=2L-Q$ now gives strict positivity.

The Euclidean interface form is also strictly positive for arbitrary
$L^2$ speeds supported in $[0,R]$ on the three rays, at each finite
$R$. The kernels of (S14) become
$k_t(r,s)=(2\pi)^{-1}\exp(-(r^2+s^2)/2+trs)$ on the measure $dr$.
Write $K_+=2k_1+k_{-1/2}$ and $K_-=2(k_1-k_{-1/2})$. Direct
completion of the square gives $K_+1=a$ and $K_-r=ar$, where $a$ is
the common local coefficient
$\phi_{\mathrm{1D}}(0)[2\Phi(r)+e^{-3r^2/8}\Phi(-r/2)]$.
Thus the sector forms are the squared-difference identities for $h$
in the standard sector and $h/r$ in the sign sector. Restrict one
integration variable to $[0,R]$ and the other to $[R+1,R+2]$, where
the speed is zero. The integrands then have coefficients $K_+(r,s)$
and $K_-(r,s)s/r$ multiplying $h(r)^2$, respectively. Both have
positive lower bounds there, including the limit at $r=0$ for the
second coefficient. This proves a strictly positive gap depending on
$R$, without a mass constraint. A compact improvement would therefore
require a finite deformation with nontrivial spatial extent in noise
units.

## Blow-ups of moving maximizers minimize the Euclidean cut

Let $A^{(j)}$ be the aligned planar maximizers above and let $x_j$ tend
to any finite $x_*$. Along a subsequence the rescaled indicators
$1_{A_i^{(j)}}(x_j+e_j u)$ converge in $L^1_{\rm loc}(du)$ to an
ordinary partition $U$. Every such limit minimizes the Euclidean
Gaussian-kernel cut against compact replacements, with arbitrary
replacement volumes.

Start with a local weak-star limit $u$ in $L^\infty$, which is
simplex-valued. Replace the labels within any fixed rescaled ball by
an ordinary partition there. This changes physical volumes by
$O(e_j^2)$. Two fixed tangent-normal bump fields on different open
simplex rays away from $x_*$ have a uniformly invertible volume
derivative matrix for the original $A^{(j)}$, by the bulk volume
formula and $L^1$ convergence. They correct the mass by flow parameters
$O(e_j^2)$, since the second volume derivatives computed by density
pullback are uniformly bounded. Their supports avoid the changed ball.

For the outer correction, pull back the joint density by the fixed
smooth bump flow family. Its first parameter derivatives vanish when
both endpoints avoid the supports. Otherwise they are bounded by
$C(1+|z|^2)$ times the joint density at the deformed endpoints, by
boundedness of the flow derivatives and Lipschitz dependence on the
endpoints. The flow and its inverse are $C^1$ close to the identity
for small parameters. Thus its short-distance exponent still improves
by a fixed Gaussian factor over the widened noise length $2e_j$.
For the midpoint factor use $|w|\le R+e_j|z|/2$ where a derivative is
nonzero. A bounded displacement of $w$ in the Gaussian exponent gives
at most an additional $C_R+C e_j|z|$, and the change in midpoint
variance gives $O(e_j^2)|w|^2$. The negative difference exponent
absorbs these and the polynomial. The parameter derivatives of
disagreement are therefore bounded by $C N_{2e_j}(A^{(j)})=O(e_j)$,
uniformly for small parameters. The outer correction costs $O(e_j^3)$.
Its mixed interaction with the inner change is exponentially small by
fixed positive separation of their supports.

Divide optimality by $e_j^2$. The rescaled joint density is, after
division, $\phi_{a_j}(x_j+e_j(u+v)/2)\phi_1(v-u)$. In the energy
difference one endpoint lies in the replacement ball. Truncate the
other to a large ball. The kernel converges uniformly there and can
be approximated uniformly by finite sums of products of continuous
functions of the two variables. Weak-star convergence of the two
label factors then gives convergence of the truncated pair integrals.
The remaining Gaussian tail is uniformly small. This proves the
limiting comparison for the fractional profile and ordinary
replacements, interpreting agreement as the dot product of the label
vectors in the finite energy difference.

On the replacement ball, any fractional profile is a weak-star limit
of ordinary profiles. Subdivide fine grid cells with volume proportions
equal to its label averages in each grid cell. The same kernel argument
extends the limiting comparison to these fractional replacements. If
two components of $u$ are positive on a set of positive measure, choose
a bounded nonzero compactly supported variation $h$ transferring between
them, on a subset where both have a uniform positive lower bound.
Both signs of a sufficiently small $h$ stay in the simplex. The sum
of the cut changes for $u+h$ and $u-h$ is
$-2\sum_i\langle h_i,\phi_1*h_i\rangle<0$, by the strictly positive
Gaussian Fourier multiplier. This contradicts minimality. Hence $u$
is ordinary, and local weak-star convergence of indicators to indicators
gives strong local $L^1$ convergence. Its cut comparison is (P5) with
$Y$ replaced by this limiting partition.

## Uniform phase density and two-label blow-ups on open rays

On every bounded physical ball there are constants $c,r_0,R_0>0$
independent of $j$ such that every point $x_0$ where label $i$ is a
score winner satisfies

\[
 |A_i^{(j)}\cap B_s(x_0)|\ge c s^2
                       \quad (R_0e_j\le s\le r_0).
 \tag{P6}
\]

Use the noise-length variables about $x_0$, abbreviate $e=e_j$, and let
$V(r)$ be the area of phase $i$ in the rescaled ball $B_r$. The pair
weight after dividing the joint law by $e^2$ is
$W(u,v)=\phi_a(x_0+e(u+v)/2)\phi_1(v-u)$. It satisfies $W\le C\phi_1(v-u)$
globally and $W\ge c_L>0$ when $|u|\le r\le r_0/e$ and $|v-u|\le1$,
uniformly on the given bounded physical ball.

Recolor $E=A_i^e\cap B_r$ completely into one other phase and correct
the masses by the outer bump construction. The two bumps may be
chosen from finitely many fixed, well-separated patches on two
different simplex rays, so each chosen support stays a fixed positive
distance from $B_{r_0}(x_0)$. Their volume matrices are uniformly
invertible for large $j$. The correcting flow has size $O(e^2V(r))$
and costs $O(e^3V(r))$ by the uniform pullback derivative above. The
mixed interaction, divided by $e^2$, is bounded by
$O(e V(r))$ via fixed physical separation and the Gaussian factor.
For either other label $l$, the disagreement variation before repair,
divided by $2e^2$, is $L(r)-\int_{E\times A_l^e}W$ with
$L(r)=\int_{E\times(A_i^e\setminus B_r)} W$.
Optimality bounds each such cross integral by $L(r)+CeV(r)$. Summing
gives $P_W(E):=\int_{E\times E^c}W\le3L(r)+CeV(r)$.

The lower bound is $P_W(E)\ge c_L'\min(V(r),\sqrt{V(r)})$. To see
it, the lower bound on $W$ controls the averaged unit-scale difference
quotients of $1_E$. A unit-scale smooth mollification $v$ therefore
satisfies $\|v-1_E\|_1+\|\nabla v\|_1\le C P_W(E)$. If the first
term is small compared to $V(r)$, then $v$ is at least $1/2$ on a set
of area at least $V(r)/2$, and $\|v\|_2\ge c\sqrt{V(r)}$. The planar
bound $\|v\|_2^2\le\|\partial_1 v\|_1\|\partial_2 v\|_1$ follows
by bounding $|v|$ by the integral of the directional derivative along
each coordinate line, multiplying the bounds and integrating. Otherwise
$P_W(E)\ge cV(r)$. Choose $r_0$ small, then $e$ small, to absorb the
repair term, since $e\sqrt{V(r)}\le\sqrt\pi r_0$.

For $R\le r_0/e$, integrate the resulting bound between $R/2$ and
$R$. Each pair with $u$ inside and $v$ outside the ball contributes
for a range of radii at most $|v-u|$. The Gaussian first absolute moment
is finite, so $\int_{R/2}^R L(r)dr\le C V(R)$. Thus

\[
 V(R)\ge c_1 R\min(V(R/2),\sqrt{V(R/2)}).
\]

At the center, being a winner and (P4) imply $T_\rho1_{A_i}(x_0)\ge1/4$
for large $j$. A fixed ball of radius $R_0$ in noise units captures all
but $1/8$ of the conditional Gaussian uniformly over bounded centers,
and its density in those units is bounded above. Thus $V(R_0)\ge m_0>0$.
On larger balls, $\min(V,\sqrt V)\ge\min(1,\sqrt{m_0})\sqrt V$.
Apply the recurrence at dyadic radii to get $V(R)\ge c_2R^2$, since
the ratio $V(R)/R^2$ is bounded below at each step by a positive
multiple of the square root of the preceding ratio. Monotonicity of
$V$ covers intermediate radii. Scaling back proves (P6).

The Gaussian measure convergence and (P6) imply that the closed
winner sets converge locally in distance to the corresponding closed
simplex cells. A winner point at a fixed positive distance from its
limiting cell would place a fixed positive amount of that phase inside
a region where its mass tends to zero. The reverse distance follows
from measure convergence on balls centered in the interior of a
limiting cell, then approximation of its boundary by interior points.
In particular a third phase cannot be present in a fixed ball around
a point of an open simplex ray that stays away from the other rays
and junction.

At centers on a two-label interface converging to such a point, both
labels survive in every microscopic limit by (P6). No third phase
survives. Every limit is a two-phase compact-perturbation minimizer
for the Euclidean Gaussian kernel. Cinti–Serra–Valdinoci's Corollary
1.23 implies its boundary is a line. Their hypotheses (1.3)–(1.7)
hold after multiplying the kernel by a constant to get $K\ge1$ on
$B_2$. Symmetry and integrability are direct, and the derivative
envelope is bounded by $C(1+|z|^4)e^{-|z|^2/8}$ because their second
derivative is taken at points within $|z|/2$ of $z$. The plane passes
through the origin of the rescaled coordinates: the two scores at
the interface center are equal up to their price difference, and
strong local convergence with Gaussian domination passes these scores
to the convolutions of the limiting halfplanes, which tie only on
the interface. This identifies the possible shapes along open rays.
The orientation and the three-label junction require further control.
[Cinti–Serra–Valdinoci, Corollary 1.23](https://cvgmt.sns.it/media/doc/paper/2917/Quantita-Nonloca-Perim10.pdf).

## Intermediate scales and localization of triple ties

Consider the aligned planar maximizers $A^{(j)}$, with $x_j\to x_*$
and scales $\ell_j\downarrow0$ satisfying $\delta_j=e_j/\ell_j\to0$.
Their labels in the coordinates $x_j+\ell_j u$ are precompact in
$L^1_{\rm loc}$. Each limit minimizes Euclidean perimeter against
compact changes of arbitrary volumes. Every nonempty phase has
infinite area, so the possibilities are a constant label, two
halfplanes, or three $120^\circ$ sectors with a translated junction.
Along each convergent subsequence the local normalized disagreement
energies converge at balls whose boundary carries no limiting
perimeter.

For the proof, write $G_\delta(v-u)=\delta^{-2}\phi_1((v-u)/\delta)$
in two dimensions. The local rescaled energy $J_j(U;B_R)$ is
$1/\delta_j$ times the integral of the label disagreement on pairs
with at least one endpoint in $B_R$, weighted by

\[
 W_j(u,v)=\phi_{a_j}(x_j+\ell_j(u+v)/2)G_{\delta_j}(v-u).
\]

The corresponding physical normalization is $e_j\ell_j$.
A compact change in $B_R$ gives a physical volume error $O_R(\ell_j^2)$.
Two fixed volume bumps away from $x_*$ correct this by the uniform
pullback argument above, costing $O_R(e_j\ell_j^2)$. The mixed
interaction with the change is exponentially small because the
physical supports are separated. The repair thus costs $o(1)$ after
normalization. Filling $B_R$ by one label leaves only crossing pairs,
whose unnormalized integral is bounded by $C_R\delta_j$. Indeed for
each displacement $h$ the crossing area is at most the perimeter of
$B_R$ times $|h|$, and the Gaussian has first moment of order
$\delta_j$. Optimality bounds $J_j$ on each fixed ball.

On fixed compact sets the midpoint factors converge uniformly to
$\phi_1(x_*)>0$. The local mollification and directional difference
proof of (P2)–(P3), applied on smaller balls with both endpoints in
$B_R$, gives local $L^1$ compactness and

\[
 \liminf_j J_j(A^{(j)}_{x_j,\ell_j};B_R)
 \ge \sqrt{2/\pi}\,\phi_1(x_*)\mathcal P(U;B_R).
 \tag{P7}
\]

Here $A^{(j)}_{x_j,\ell_j}$ denotes the rescaled labels and
$\mathcal P$ is the unweighted interface perimeter counted once.
Take exhaustion by balls to obtain a locally BV limit. For any fixed
locally finite-perimeter competitor $V$ agreeing with $U$ outside
$B_R$, glue $V$ to the rescaled maximizers in $B_t$, with $t\in(a,b)$
and $R<a<b$. The possible increase of its local cut over that of
$V$ in $B_t$ lies in pairs crossing $\partial B_t$. It is bounded
above by

\[
 2\int_{|x|<t<|y|} W_j(x,y)
                 1_{\{A^{(j)}_{x_j,\ell_j}(y)\ne V(y)\}}dx\,dy.
\]

Averaging in $t$ restricts the interval of radii for a pair to one
of length at most $|x-y|$. On $y\in B_{b+1}\setminus B_a$ the
averaged integral is at most $C\delta_j$ times the local $L^1$
mismatch with $U$, which is $o(\delta_j)$. For $y$ outside $B_{b+1}$
it is a Gaussian tail since $|x-y|\ge1$ and $|x|\le b$. Choose
$t_j\in(a,b)$ with error $o(\delta_j)$. Compare with the repaired
glued partition, and use nonnegativity to restrict the left to
$B_a$. This gives

\[
 \limsup_j J_j(A^{(j)}_{x_j,\ell_j};B_a)
                         \le \limsup_j J_j(V;B_b).
\]

For fixed $V$, the right side is at most
$\sqrt{2/\pi}\,\phi_1(x_*)\mathcal P(V;\overline{B_b})$.
In fact each displacement $\delta_j z$ counted in the integral
has one endpoint in $B_b$. The scalar directional BV estimate
bounds the integrated label-$i$ difference by
$\delta_j |z\cdot D1_{V_i}|(B_{b+2\delta_j|z|})$. This estimate follows
from integrating the distributional directional derivative along
segments, first for smooth mollifications and then by BV convergence
on a slightly larger ball. For $|z|\le1/\delta_j$ dominate the
normalized bound by $|z|$ times the variation in $B_{b+2}$, and use
the uniformly convergent midpoint weight on fixed compacts. For
$|z|>1/\delta_j$ bound the pair difference by one and the possible
first-endpoint area by $2|B_b|$, and use a Gaussian tail. This proves
the stated upper bound. Taking $b\downarrow a$ at radii with no
limiting interface mass on the sphere and using (P7) proves perimeter
minimality. Taking $V=U$ also gives the claimed local energy convergence.

Every phase of $U$ of positive area has infinite area. Choose winner
points of that label in bounded rescaled coordinates by $L^1_{\rm loc}$
convergence. For any fixed rescaled radius, (P6) applies for large $j$
since $R_0e_j/\ell_j\to0$. Pass the estimate to the limit. As the
rescaled radius can be arbitrarily large, the phase has infinite area.
The planar classification then follows from Novaga–Paolini–Tortorelli,
Theorems 4.2 and 4.5. Those theorems concern compact volume-preserving
comparisons, which are included here. Their all-infinite-area cases
with up to three phases give precisely the three possibilities above.
[Novaga–Paolini–Tortorelli, arXiv:2312.13709](https://arxiv.org/pdf/2312.13709).

For all sufficiently large $j$, every fixed small disk about the
limiting junction contains at least one triple score tie. To prove
existence, take three points on the boundary circle in the interiors
of the limiting sectors. Along each successive arc between such
points only the two corresponding labels can win by (P6) and local
distance convergence. At the endpoints the winning labels are unique.
Write $s_i$ for the scores minus their prices. If there were no triple
equality in the disk, $\max_i s_i-\min_i s_i$ would have a positive
minimum on it. For smaller $\varepsilon>0$, the map with components
$(s_i-\max s+\varepsilon)_+$, divided by their sum, would be a
continuous map of the disk to the boundary of the probability triangle.
On each of the three boundary arcs the excluded label has a strictly
positive minimum gap. Decrease $\varepsilon$ so the three arcs map to
successive edges of the triangle and their endpoints to the distinct
corners. This boundary map has winding number one, which prevents
the asserted extension over the disk.

In any fixed compact disk around the limiting junction, the set of
triple ties has diameter $O(e_j)$. All such points tend to the
origin by local distance convergence of the winner sets. If two
triple ties had separation $\ell_j\gg e_j$, take the intermediate
scale $\ell_j$ centered at one tie. Pass to a subsequence on which
the other tie tends to a point at distance one. Each of the three
labels is a winner at both centers, so (P6) and strong convergence
give positive area proportions for all three phases in every ball
centered at either limiting point. But among the classified limiting
partitions only the three-sector partition has such a point, and it
has exactly one. This contradiction proves the diameter bound.

## Uniform gradients at winning ties near one

On any fixed compact set of the aligned active plane there is $c>0$
such that, for all sufficiently large $\rho<1$, every winning tie
between labels $i,j$ of a balanced maximizer satisfies

\[
 e\,|\nabla(U_i-U_j)(x)|\ge c.
 \tag{P8}
\]

Otherwise take a sequence of ties $x_j$ in the compact set and pass
to a noise-length blow-up limit as above, with the same pair of
labels. The rescaled conditional means and standard deviations are
$\rho_j u-e_jx_j/2$ and $\sqrt{(1+\rho_j)/2}$ in coordinates
$x=x_j+e_ju$. Strong local label convergence and Gaussian tails
therefore give local $C^1$ convergence of the rescaled scores to
$P_i=G*1_{A_i}$ for the Euclidean minimizing limit $A$.
For derivatives this follows by differentiating the conditional
density, with its uniformly integrable Gaussian tail times a linear
polynomial. By (P4) the price differences tend to zero, so the limiting
scores have a winning tie at $u=0$. A vanishing scaled gradient would
contradict the [Euclidean compact recoloring argument](EUCLIDEAN_THREE_CELL_JUNCTION.md#nonzero-gradients-at-winning-score-ties).

## The unique noise-length junction

Every noise-length blow-up at bounded centers is a constant
partition, a two-phase halfplane partition, or a regular
three-sector partition. The limit is a Euclidean compact-perturbation
minimizer by the earlier comparison. In the two-phase case apply the
two-phase classification above. In the three-phase case apply
the [strict Euclidean comparison and classification](EUCLIDEAN_THREE_CELL_JUNCTION.md#strict-comparison-by-fiber-rearrangement).

Let $p_j$ be a triple score tie of aligned near-one maximizers.
At the noise scale centered at $p_j$, the three limiting scores
tie at the origin by (P4) and convergence of the convolutions.
The limit therefore has three phases and is a tripod centered at
zero. Indeed in a constant or two-phase limit the absent label
has zero score and cannot be a simultaneous winner. The tripod's
only triple tie is its vertex, by its strict interior score
margins and the two-way tie on each open ray. The conclusions hold
along every subsequence after extracting orientations.

Write $F_j(u)$ for the two centered score differences
$(s_1-s_2,s_2-s_3)$ at $p_j+e_j u$. They converge locally in
$C^k$ for every fixed $k$ along a subsequence, by differentiating
the rescaled Gaussian densities as in (P8). The limiting map is
invertible to first order at zero. Each incident pair difference
has nonzero normal derivative there by the linear score margin
of (P5), and the two supporting lines are distinct. Thus, for
some small ball, $DF_j$ throughout that ball is arbitrarily close
to the invertible limiting derivative at zero for large $j$.
Integrating the derivative on line segments gives injectivity on
the ball. Every other triple tie in the fixed physical disk is
at bounded noise-scale distance by the diameter bound. A
subsequence of such ties would converge to a tie of the
limiting tripod and hence to its vertex. Injectivity then
rules it out. This proves uniqueness of the triple tie in any
fixed bounded disk containing the limiting origin in its
interior, for all sufficiently large $j$.

The inverse score coordinates give three analytic arcs through
that point, with incident angles tending to $120^\circ$.
In those coordinates the winner regions are the three polyhedral
cells of the largest of three scores. Away from zero in any
fixed noise-scale ball, the implicit function theorem for the
incident pair and the strict price-adjusted score gap for the
third give convergence of the interfaces to the corresponding
rays, with all derivatives on compact open-ray portions.

There is also a separation of the two local models. If $x_j$
are bounded winning pair ties with
$d_j=|x_j-p_j|\gg e_j$, every microscopic limit at $x_j$ is
a halfplane pair. For $x_j$ tending to a nonzero point, use
local class-distance convergence to the physical tripod. If
$x_j\to0$, apply (P7) centered at $p_j$ at scale $d_j$.
The triple tie and (P6) force a tripod centered at zero in the
limit. The two labels at $x_j$ survive at a point at distance
one, which must lie on an open ray of the rescaled tripod.
There is then no third winner in a ball of radius $c d_j$
around $x_j$ for large $j$, with $c>0$ by the class-distance
consequence of density at that scale. Rescaling further at
$e_j$ yields the two-phase halfplane classification.

## A bounded physical defect cannot maximize

For every fixed $R<\infty$, if $\rho$ is sufficiently close to one,
a balanced planar global maximizer that agrees with a labeled
central simplex off $B_R$ is that simplex almost everywhere.

Suppose otherwise along a sequence $\rho_j\uparrow1$. Use the
exterior simplex as reference. The compactness conclusion from
(P3) aligns the maximizers to this same reference in the limit:
any limit simplex agrees with it off $B_R$, which fixes the
labeled sectors. On each open exterior ray the incident adjusted
scores agree. The unadjusted score difference is real analytic
on the supporting line, so it equals its pair price difference
on the entire line. At points tending to infinity along the
ray, the unadjusted difference tends to zero. The central
simplex scores tie identically there by reflection, and the
conditional probability of hitting the bounded defect tends
to zero for each fixed $\rho_j>0$. Thus all pair prices are equal
and the three unadjusted differences vanish on their supporting
lines. This use of line analyticity is the same mechanism as
in the flat-interface argument of Heilman–Mossel–Neeman.
[Corollary 5.1](https://www.stevenheilman.org/~heilman/papers/PluralityV5.pdf)

There is therefore a triple tie at the origin, which is the
only one in any fixed larger disk by the preceding section.
On each positive simplex ray the incident labels are tied
everywhere and win near the outer endpoint, where the partition
agrees with the exterior simplex. On the open ray their common
score cannot cross the third score without another triple tie.
They therefore win along the entire ray in the disk.
By the bounded graph result below, each pair interface is
a unique normal graph on bounded outgoing ray segments.
Since each reference ray already consists of winning pair ties
from the same vertex, those graphs coincide with the reference
rays. All other points of the disk have the reference labels
by continuity and the exterior assignment. This is a
contradiction. The comparison excludes bounded defect support
for actual global maximizers, including supports much larger
than a fixed number of noise lengths.

## Tangent alignment across bounded physical regions

Along any sequence of aligned planar balanced maximizers as $\rho_j\uparrow1$,
the oriented normals at winning pair ties converge uniformly on fixed
physical compact sets to the corresponding simplex normals. Center
the branch coordinates at the unique nearby triple tie. On each fixed
bounded interval of its outward coordinate the corresponding interface
is a single analytic normal graph $g_j$, with $g_j(0)=0$,
$\|g_j\|_{C^1}\to0$ and $e_j^{k-1}\|g_j^{(k)}\|_\infty\to0$ for
each fixed $k\ge2$.

The proof compares stress fluxes across transverse cuts of increasing
length. It uses the fact that every microscopic limit is flat on its
open interfaces and that the force on an intermediate region tends to
zero with its physical radius.

### Stress and its force density

Use coordinates of length $e$ centered at a bounded physical
point $x_*$, so $u_i(w)=1_{A_i^{(j)}}(x_*+ew)$, and write
$q(x,y)$ for label disagreement in these coordinates. Suppress
the index $j$. Put
\[
 w_j(w)=\phi_{a_j}(x_*+ew),\qquad
 \gamma_j(w)=\phi_1(x_*+ew),\qquad
 b_j(w)=\sum_i\beta_{j,i}u_i(w).
\]
Here $w_j$ denotes the scalar weight, and $w$ its argument.
The rescaled disagreement, after dividing the physical disagreement
by $e^2$, uses the ordered pair density
$G(y-x)w_j((x+y)/2)$.
Define
\[
 \begin{aligned}
 S_j(w)
 &=2I\int G(z)w_j(w+z/2)q(w,w+z)\,dz\\
 &\quad-\int G(z)z\otimes z\int_0^1
    w_j(w+(1/2-t)z)q(w-tz,w+(1-t)z)\,dt\,dz,\\
 f_j^0(w)
 &=\int G(z)\nabla w_j(w+z/2)q(w,w+z)\,dz .
 \end{aligned}
\]
In these coordinates $\nabla w_j$ includes the factor $e$.
Compact pullback by $w\mapsto w+tV(w)$ yields
$\int S_j:DV+\int f_j^0\cdot V$. The midpoint term uses
$[V(x)+V(y)]/2$ and the two endpoint contributions are
equal by pair symmetry. The difference kernel derivative is
$-z\cdot(V(y)-V(x))$. Write the endpoint difference as its
integral along the segment to obtain the stress above.
By (P4), physical $N_\rho+2\sum_i\beta_{j,i}m_i$ is stationary
under these flows, with $m_i$ the Gaussian cell masses.
The volume derivative divided by $e^2$ is
$\int u_i\operatorname{div}(\gamma_j V)$. Hence in distributions
\[
 \operatorname{div}\mathsf T_j=f_j,\qquad
 \mathsf T_j=S_j+2\gamma_j b_j I,\qquad
 f_j=f_j^0+2b_j\nabla\gamma_j .
 \tag{WS1}
\]
The divergence is row-wise. No power of $e$ multiplies the pressure
tensor, since both physical volumes and physical cut have the Jacobian
$e^2$ under the rescaling.

The following bounds are uniform for bounded physical centers.
For $1\le R\le r_0/e$, with $r_0$ small and fixed,
\[
 \int_{B_R}|f_j(w)|\,dw\le C(eR+e|\beta_j|R^2).
 \tag{WS2}
\]
In fact, replacing the labels inside $B_R$ by one label leaves
only crossing pairs and their cut is at most $CR$, by the
Gaussian first moment and the disk perimeter. It changes the
physical volumes by $O(e^2 R^2)$. Two remote fixed flow bumps,
chosen as in (P6), correct the volumes, with uniformly invertible
volume derivative and physical cut cost $O(e^3 R^2)$. The
inner-outer repair interaction is exponentially small because
the supports are physically separated. Thus the cut in the
original $B_R$ is at most $CR+CeR^2\le C'R$ in scaled units.
The midpoint weight is bounded below uniformly in this range for
pairs with $|z|\le1/e$ and at least one endpoint in the ball.
The remaining Gaussian tail is at most $CR^2e^{-c/e^2}$.
It follows that the unweighted order-one kernel disagreement
with one endpoint in $B_R$ is at most $CR$. Since
$\|\nabla w_j\|_\infty+\|\nabla \gamma_j\|_\infty\le Ce$,
(WS2) follows. The same proof works in a ball centered elsewhere
in the common bounded physical region. The radius and remote bumps
can be chosen uniformly by using the finite list from (P6).

If $B_d(w)$ lies in one phase up to null sets, then
\[
 |\mathsf T_j(w)|\le Ce^{-c d^2}+C|\beta_j|
 \tag{WS3}
\]
almost everywhere. All disagreement terms in $S_j$ vanish for
$|z|<d$; absorb the stress polynomial weight in the Gaussian
tail. To take ordinary line fluxes, convolve $\mathsf T_j$ and
$f_j$ with one fixed smooth nonnegative kernel of integral one
supported in the unit ball. This preserves (WS1), gives (WS2)
on $B_R$ after enlarging the center ball by one, and gives
(WS3) with $d$ decreased by one. Use these mollified fields
in the rest of the argument.

For an oriented straight cut $\Gamma$ in the comparison region,
write its flux as
\[
 \mathfrak F_j(\Gamma)=\int_\Gamma \mathsf T_j J^\top\,dw ,
\]
where $J$ is counterclockwise quarter-turn and $dw$ is the
oriented tangent line element. If two cuts are both oriented from
the $k$ side to the $i$ side and their corresponding endpoints
can be joined in deep regions of their respective phases, close
the two cuts with these joining segments. Stokes' theorem applied
row-by-row introduces the integral of $f_j$ over the oriented
polygon. A crossing polygon can be triangulated into two oriented
triangles and has the same bound up to a constant. In a ball of
radius $CR$, with joining segments of total length $O(R)$
at phase depth at least $cR$, (WS2) and (WS3) give
\[
 |\mathfrak F_j(\Gamma)-\mathfrak F_j(\Gamma')|
 \le C\left(eR+e|\beta_j|R^2+|\beta_j|R+R e^{-cR^2}\right).
 \tag{WS4}
\]
All comparisons below can be arranged with this geometry.

If translations of the labels about tie points converge locally
to a binary halfplane with inward normal $n$ towards $i$,
and the physical points tend to $p$, translated stresses
converge locally in $L^1$ to $\phi_1(p)$ times the Euclidean
halfplane stress. To check this first truncate $|z|$. In each
stress summand, local $L^1$ convergence of the label difference
at both endpoints suffices after shifting the $w$ integration
variable. The tail is bounded uniformly by Gaussian decay.
Convolution with the fixed smooth kernel makes this convergence
locally uniform. The whole normal flux for the halfplane is
\[
 \int_\mathbb R \phi_1(p) S_H(sn)J^\top n\,ds
   =\phi_1(p)\sqrt{2/\pi}\,J^\top n .
 \tag{WS5}
\]
This is (J8). Mollification does not affect the full integral,
since the stress is independent of tangential position. A cut
of fixed large half-length $R_0$ makes only
$O(e^{-cR_0^2})$ error in (WS5).

### Uniform flux comparison at shrinking physical scales

First center the coordinates at the triple tie $x_{0,j}$. Pass
to a subsequence with unit-scale limit a tripod $Y_*$, which
exists by unit-scale compactness and has vertex zero by score
convergence. Every centered scale $R$ with $1\le R\le L_j$,
where $L_j\to\infty$ and $e L_j\to0$, is uniformly close to
some labeled tripod with vertex zero. More precisely, on
$B_{4R}$ its closed winner sets, after division of distances
by $R$, have local distance discrepancy tending to zero from
such a tripod. Otherwise a contradicting sequence has either
bounded $R$, where the unit limit and score convergence apply,
or $R\to\infty$, where (P7) applies. All three labels meet the
limiting center by (P6), so the (P7) limit is a tripod with
its vertex at zero. The passage from measures to winner sets at
intermediate scales follows from (P6). At bounded noise radii,
use uniform score convergence and the strict limiting score gaps.

For each pair ray of the approximating tripod at scale $R$,
take a cut centered at radius $R$ on that ray, with endpoints
a distance $\eta R$ from the ray on its opposite phase sides,
for fixed sufficiently small $\eta>0$. It is oriented from
phase $k$ to phase $i$. The approximating normals at $R$
and $2R$ are close by overlap of the corresponding scaled
tripod approximations with the same labels. Thus corresponding
cut endpoints join by deep straight segments: both endpoints
are within the same reference wedge with depth at least
$cR$, which is preserved along the joining segment by
convexity. The winner-set errors, for large $j$, are smaller
than $cR/2$. The same holds for neighboring radii of ratio
between one and two. Sum (WS4) over dyadic radii from a fixed
large $R_0$ to any $R\le L_j$. The pressure terms sum to
$O(|\beta_j| L_j+e|\beta_j| L_j^2)=o(1)$ by (P4). The
inhomogeneous stress term sums to $O(eL_j)=o(1)$.
The depth tails sum to $O(e^{-c'R_0^2})$.
At the starting radius, $Y_*$ can be used as approximating
tripod, and its cut flux tends to the Euclidean halfplane
flux in (WS5) with error $O(e^{-cR_0^2})$. To see this,
at radius $R_0$ the Euclidean tripod agrees with the incident
halfplane in a ball of radius comparable to $R_0$ around the
cut. Contributions to the stress from outside that ball and
the whole-line halfplane flux tails are exponentially small.
The physical weight there tends to $\phi_1(0)$.

Consider any active tie $z_j$ of this pair with
$R=|z_j|\to\infty$ and $R\le L_j$, still in these unit
coordinates. The third label is absent in $B_{c_1R}(z_j)$
by the centered tripod approximation, for a fixed small
$c_1>0$. Every unit translation limit at $z_j$ is therefore
a two-phase compact-perturbation minimizer. Its halfplane passes
through the center by score convergence. Denote its normal by
$n$ along a convergent subsequence. Its intermediate
approximations on every scale $r$ with $R_0\le r\le c_2R$
are also uniformly close on the comparison balls to labeled halfplanes through $z_j$,
where $c_2>0$ is smaller than $c_1/10$. To check the uniform
assertion, use the starting microscopic halfplane limit for bounded $r$.
For $r\to\infty$ with $r/R\to0$, use (P7) at $z_j$;
the third phase is absent from every fixed rescaled ball.
For $r/R$ bounded below, use the centered tripod approximation
at $R$, restricted to the small ball about $z_j$. The
density estimate makes these approximations comparisons of
winner sets, with both phases meeting their center in all
limits.

Start at $z_j$ with the cut along the limiting normal $n$
and half-length $R_0$, and then take cuts along successive
halfplane normals up to $r$ comparable to $c_2R$. The
starting sphere uses the unit limit with $R_0$ fixed.
Neighboring cuts again have joining segments at depth $cr$:
the normals at neighboring dyadic scales are close by overlap,
and both endpoints of each joining segment are on the same
side of the approximating line by at least a fixed multiple
of the smaller radius. Join the last cut to the canonical cut
at radius $R$ of the centered tripod. Both pairs of endpoints
are at phase depth at least a fixed multiple of $R$, and their
segments stay in the correct approximate convex wedges.
The extra comparison costs the right side of (WS4) at $R$.
Thus the micro cut at $z_j$ agrees with the micro tripod
flux up to $O(e^{-c'R_0^2})+o(1)$. Taking $j\to\infty$
then $R_0\to\infty$ in (WS5) shows that its microscopic
normal must agree with the corresponding $Y_*$ normal.
At bounded $z_j$ in these units the same conclusion follows
directly from $C^1$ score convergence to the unit tripod.

To identify $Y_*$ with the aligned macroscopic simplex
$T$, choose $r_j=eL_j\downarrow0$ slowly enough that
$r_j\gg |x_{0,j}|$, $r_j\gg e$, and the $L^1$ error
of physical $A^{(j)}$ versus $T$ on a fixed ball is
$o(r_j^2)$. If a particular bounded physical test sequence
tends to the origin, $r_j$ can also be chosen much larger
than its distance to the origin. Then rescaling about
$x_{0,j}$ at $r_j$ gives the aligned $T$, locally.
If a pair normal $n_*$ of $Y_*$ differed from the
corresponding normal $n_T$, choose a vector $v$ with
$v\cdot n_*>0$ and $v\cdot n_T<0$. A small segment
parallel to $v$ crossing the open pair ray of $T$ in
this last rescaling lies where the third phase is excluded
and has opposite score signs at its endpoints, from
positive to negative. Every intervening zero of the pair
score is an active tie. All such ties have positive
directional derivative along $v$ for large $j$, by
the flux comparison and microscopic $C^1$ score convergence.
This is inconsistent with a positive-to-negative crossing.
The segment can be chosen inside the scaled ball $B_{L_j}$
and away from the junction. This establishes the desired
normal alignment for sequences of active ties with physical
limit zero.

If instead active tie points $x_j$ tend to a point $p\ne0$
on an open simplex ray, no third phase is present in a fixed
physical neighborhood by (P6). Choose $r_j\to0$ slowly enough
that $r_j/e\to\infty$, $|x_j-p|=o(r_j)$, and local physical
$L^1$ error is $o(r_j^2)$. The outer scaled neighborhood
about $x_j$, of radius $L_j=r_j/e$, now approaches a
halfplane with the macroscopic normal. Take a unit translation
subsequence at $x_j$ and let its halfplane normal be $n_*$.
Apply the same binary dyadic flux comparison from the
microscopic cut at $x_j$ to a canonical cut of length
$\eta L_j$ in the outer halfplane. At any other tie point
in a fixed small fraction of the outer ball the identical
argument compares that tie's microscopic cut to the same
canonical cut. For this comparison the top binary cut can be
centered at the new tie, and its deep endpoints join those
of the canonical cut by straight segments on the same
halfplane sides. Compactness gives the required intermediate
binary approximations by the same cases of bounded rescaled
radius, diverging radius $o(L_j)$, or comparable to $L_j$.
In (WS5) the micro weight factor is $\phi_1(p)$ for all
these ties. Therefore the micro normals agree with $n_*$
throughout the smaller outer ball. The same positive-to-negative
crossing test forces $n_*=n_T$.

### From tangent alignment to normal graphs

The sequential result gives uniform tangent convergence at all
active ties in every fixed physical ball. The unit-score
gradient magnitude bound follows by subsequences: at bounded
distance from the triple tie in noise units it is the
strict simplex pair margin derivative on the appropriate
supporting line, and on all other subsequences it is the
binary halfplane derivative. For a fixed bounded arc of a
reference ray, the actual interface can be described as a
unique normal graph starting at the shifted vertex. For
distances much larger than the noise length, the centered
tripod approximations give the appropriate labels on opposite
sides of each narrow tube, with the third phase absent.
Any limit of these approximations has the aligned orientation
by the same crossing argument and uniform tangent estimate.
At bounded distances in noise units the unit tripod score
convergence and the inverse function theorem at the unique
triple tie give the three graphs. Normal-line uniqueness
further out follows because every active pair zero has
positive normal derivative. The physical graph offsets and
their slopes tend uniformly to zero on bounded regions.
For each fixed derivative order $k\ge2$, the normalized
derivative $e^{k-1}g_j^{(k)}$ tends uniformly to zero there,
by applying the implicit function theorem in the local unit
coordinates and using $C^k$ score convergence. In the unit
limits the pair tie line of either a tripod or a halfplane
is exactly straight.

## Exact comparison on bounded noise-scale supports

For each fixed $R>0$ there are $c_R>0$ and $\rho_R<1$ such that
for $\rho_R<\rho<1$, a planar ordinary partition $B$ agreeing
with the central Gaussian simplex $T$ outside $B_{Re}(0)$ satisfies

\[
 S_\rho(T)-S_\rho(B)\ \ge\
            c_R\,\frac{|\{B\ne T\}|^2}{e^2}.
 \tag{P9}
\]

Here $|\cdot|$ is planar Lebesgue area and $e=\sqrt{2(1-\rho)}$.
The estimate needs no cell mass constraint. It is a uniform
comparison on each fixed noise-scale support radius.

In coordinates scaled by $e$, let $b$ be the vector label
function of $B$, $u=1_Y$, and $m=|\{b\ne u\}|$, with equality
to $u$ off $B_R$. The stability difference divided by $e^2$ is
the relative cut $\mathcal C_e(b)$ for the symmetric pair weight
$W_e(x,y)=\phi_a(e(x+y)/2)G(y-x)$. The limiting weight $W_0$ is
$\phi_1(0)G$. For $e_j\to0$ these relative cuts converge along
weak-star convergence of simplex-valued profiles in $B_R$,
with the cut of a soft profile given by one minus the dot
product of its label vectors. To see this, expand the relative
cut as $2L_e-Q_e$ in $b-u$. The regret integrates a score
converging uniformly on $B_R$, by uniform Gaussian tail bounds.
The quadratic term converges by uniform kernel convergence
on $B_R^2$ and approximation of the continuous limiting kernel
there by finite sums of products.

The limiting relative cut is nonnegative for every soft
profile. Indeed approximate the soft profile by ordinary
partitions with the same exterior labels, allocating the label
averages in each small grid cell by measurably subdividing it,
and use the Euclidean comparison. If a zero-cost soft profile
has two components positive on a set of positive area, choose
a nonzero bounded transfer $h$ between them supported where
both components are bounded below. The average of the cuts
of $b+h$ and $b-h$ drops by
$\phi_1(0)\sum_i\langle h_i,G*h_i\rangle>0$ because the Gaussian
Fourier multiplier is positive. This contradicts nonnegativity.
An ordinary zero-cost profile is $u$ by strict Euclidean
comparison. Thus $u$ is the sole zero-cost soft profile.

If (P9) fails, there exist $e_j\to0$ and $b_j$ distinct from
$u$ with $\mathcal C_{e_j}(b_j)/m_j^2\le1/j$. After a weak-star
subsequence, the cuts converge to zero or less since $m_j$ is
bounded by $|B_R|$. Their limit must be $u$. Weak-star
convergence of indicators to those of $u$ gives $m_j\to0$.
It remains to control the normalized cuts for these small
perturbations.

Let $E_{i\ell}^{(j)}$ be the set transferred from $i$ to $\ell$
in the $j$-th partition.
On its supporting line use tangential coordinate $r$ measured
outward along the ray and nonnegative inward coordinate $t$
on the side of $i$. The pointwise regret for this transfer is
the difference of the $W_{e_j}$ scores of $Y_i$ and $Y_\ell$.
It vanishes on the pair line by reflection. The difference
converges in $C^1(B_R)$ to $\phi_1(0)$ times the Euclidean
score difference, by Gaussian domination of derivatives.
Consequently (P5) gives a lower bound $c\,t$ on $Y_i\cap B_R$,
uniformly for large $j$. In fact the error relative to the
limiting score difference is at most $o(1)t$, by integration
of the derivative error along the normal segment to the line.
The upper bound $\mathcal C_{e_j}(b_j)\le m_j^2/j$ and the compact
kernel bound $|Q_{e_j}|\le C m_j^2$ now give
$\sum_{i\ell}\int_{E_{i\ell}^{(j)}} t\,dx=O(m_j^2)$.

For each transfer let $f_{i\ell}^{(j)}(r)$ be its section length in
$t$, divided by $m_j$. The inequality
$\int_E t\,dt\ge |E|^2/2$ for subsets of $[0,\infty)$ bounds these
functions in $L^2([-R,R])$. Pass to weak limits $f_{i\ell}\ge0$.
The measures $1_{E_{i\ell}^{(j)}}dx/m_j$ converge weakly to these
densities on the projected lines because their average normal
distance tends to zero. Their supports in the limit are on
the active rays, since the $Y_i$ side meets the pair line
only along that ray. Also $\sum_{i\ell}\int_0^R f_{i\ell}\,dr=1$.
For a fixed $\eta>0$ the part of a projected density from
$t>\eta$ has integral $O(m_j/\eta)$ and stays bounded in
$L^2$, so discarding it does not change the weak limit. For
$r\ge0$ and $0\le t\le\eta$, the pointwise regret is at least
$(\phi_1(0)a(r)-o_j(1)-\omega(\eta))t$, with
$\omega(\eta)\to0$ by smoothness of the Euclidean score and
$a(r)$ the common interface coefficient of (P5).
Use the section inequality again and weak lower
semicontinuity of the weighted $L^2$ norms, then send
$\eta$ to zero, to get

\[
 \liminf_j\frac{2L_{e_j}}{m_j^2}
       \ \ge\ \phi_1(0)\sum_{ij}\int_0^R a(r) f_{ij}(r)^2\,dr.
 \tag{P10}
\]

Orient the three rays and write $v_{ij}=f_{ij}-f_{ji}$ on each
oriented pair. Weak measure convergence on the compact support
gives $Q_{e_j}/m_j^2\to\phi_1(0)\langle v,\mathsf K v\rangle$
for the Euclidean interface interaction in (P5). This follows
from the vector signed facet measure of the transfer and
uniform convergence of the continuous kernel. Subtracting this
limit from (P10), the deficit bound would force

\[
 \sum_{i<j}\int_0^R a v_{ij}^2\,dr-\langle v,\mathsf K v\rangle
       +2\sum_{i<j}\int_0^R a f_{ij}f_{ji}\,dr \ \le\ 0 .
\]

The first two terms form the strictly coercive Euclidean
quadratic form on speeds supported in $[0,R]$, by (P5).
The last term is nonnegative. It follows that $v=0$ and
$f_{ij}f_{ji}=0$ almost everywhere for every pair, forcing all
the transfer densities to be zero. This contradicts their
unit total mass and proves (P9).

## A uniform comparison for graph deformations away from the vertex

Fix $r_0,M>0$. On each simplex ray replace the line
locally by a normal graph with height $g_i$, supported on
$[r_0,M/e]$ in its tangential coordinate. Reassign only the
signed strips between the reference ray and these graphs to the
incident labels. There are constants $c,\varepsilon>0$ and
$\rho_0<1$, depending only on $r_0,M$, such that if the zero-extended heights
have $\|g_i\|_\infty+\operatorname{Lip}(g_i)\le\varepsilon$, then
for $\rho_0<\rho<1$ this competitor $B$ satisfies

\[
 S_\rho(T)-S_\rho(B)
             \ \ge\ \frac{c}{e}\sum_i
                   \langle g_i,(I-T_\rho^{(1)})g_i\rangle_{\gamma_1}.
 \tag{P11}
\]

Here $T_\rho^{(1)}$ is the scalar Ornstein–Uhlenbeck operator
on the full tangential line. This comparison has no mass
constraint, and its amplitude threshold is independent of $e$.
The supports may cover arbitrarily many noise lengths.

First compare a binary graph $y=g(s)$ with the full flat line.
Write $\sigma=\sqrt{1-\rho^2}$ and let $C_\rho(u,v)$ be the
disagreement probability for two correlated normal coordinates
thresholded at $u,v$, respectively. Write $\phi,\Phi$ for the scalar
Gaussian density and distribution function, $p_\rho(a,b)$ for the
joint density and $D_\rho(u)=C_\rho(u,u)-C_\rho(0,0)$. Then

\[
 C_\rho(u,v)-C_\rho(0,0)=
 \tfrac12(D_\rho(u)+D_\rho(v))
   +\int_{\min(u,v)}^{\max(u,v)}
        \int_{\min(u,v)}^{\max(u,v)}p_\rho(a,b)\,da\,db .
 \tag{P12}
\]

This is the identity obtained by writing the disagreement as
the sum of the two marginal probabilities minus twice the lower
rectangle probability. Its diagonal derivative is
$D_\rho'(u)=2\phi(u)[1-2\Phi((1-\rho)u/\sigma)]$.
Hence $D_\rho(u)\ge-2\phi(0)^2(1-\rho)u^2/\sigma$, by the
Gaussian density bound.

For $(s,t)$ the correlated tangential observations, take
$u=g(s),v=g(t)$. The normal-density exponent relative to its
value at zero is $(a+b)^2/(4(1+\rho))+(a-b)^2/(4(1-\rho))$.
On the square in (P12) this is at most
$C\varepsilon^2(1+|s-t|^2/e^2)$. Using $e^{-x}\ge1-x$, the
edge term in (P12) is bounded below by
$\phi(0)^2 (g(s)-g(t))^2/\sigma$ minus this same term times
$C\varepsilon^2(1+|s-t|^2/e^2)$.
For pairs having one nonzero endpoint of $g$ the polynomial
factor times the tangential joint density is bounded by $C_R$
times the joint density with noise length $2e$, with a constant
depending only on $M$. Indeed one endpoint is at most $M/e$
in absolute value, so the change in midpoint exponent is
bounded by $C e^2(M/e+|s-t|/2)^2$. The widened difference
Gaussian absorbs the polynomial and the part depending on
$|s-t|$, leaving a bounded factor. The widened tangential Dirichlet form is at most
four times the original one by the Hermite multiplier bound
$1-(1-4(1-\rho))^n\le 4(1-\rho^n)$. Thus the binary deficit is
at least

\[
 \frac{2\phi(0)^2}{\sigma}
  \left[(1-C_M\varepsilon^2)
                 \langle g,(I-T_\rho^{(1)})g\rangle_{\gamma_1}
                        -(1-\rho)\|g\|_{L^2(\gamma_1)}^2\right].
 \tag{P13}
\]

For all sufficiently large $\rho$ the tangential Dirichlet form
on functions supported in $[r_0,\infty)$ is bounded below by
$(1+\eta)(1-\rho)\|g\|_2^2$ for some $\eta>0$ depending on $r_0$.
Otherwise normalize $\|g\|_2=1$ along a sequence $\rho\uparrow1$.
The eigenvalues $(1-\rho^n)/(1-\rho)$ increase in $n$ and
converge to $n$ for each fixed $n$. A bound on the normalized
Dirichlet forms gives a uniform Hermite tail bound and an
$L^2$-strong subsequential limit supported in $[r_0,\infty)$.
Its Gaussian derivative energy is bounded above by the liminf.
For this limit use its odd extension across zero. Its first
Hermite projection squared is at most $b$ times its squared
norm, where $b=\int_{|x|\ge r_0}x^2d\gamma_1(x)<1$. The other
possible odd Hermite degrees are at least three. Thus its
derivative energy is at least $(3-2b)$ times its squared norm.
Choosing $0<\eta<2(1-b)$ contradicts the failure of the asserted
Dirichlet bound. Now take $\varepsilon$ small enough that (P13)
bounds below a fixed positive multiple of the edge term divided
by $\sigma$.

To pass to the tripod, the signed strips lie in disjoint
tubes for sufficiently small $\varepsilon$. On the strip at
outward coordinate $s\ge r_0$, the two-phase indicators are
identical to the incident tripod indicators on a ball of
radius $cs$. The difference between the two score differences
vanishes on the whole supporting line by reflection. Its
normal derivative in the strip is bounded by
$C \sigma^{-1}e^{-c s^2/\sigma^2}$, by integrating the Gaussian
derivative outside the agreement region. The conditional mean
$\rho x$ stays at comparable distance from that region of
indicator difference. Integrating from the line along the strip
and then against Gaussian measure bounds the regret error by
$C \sigma^{-1}e^{-c r_0^2/\sigma^2}\|g_i\|_2^2$. The self
quadratic term of each strip is identical in the binary
and tripod comparisons.

For distinct strips with tangential parameters $s,t\ge r_0$,
their points are separated by at least $c(s+t)$. Conjugate
the Mehler operator to Lebesgue measure with the square root
of the Gaussian density. Its kernel is bounded by
$C\sigma^{-2}e^{-c(s^2+t^2)/\sigma^2}$ on these strip pairs.
Integrating the normal coordinates against their bounded
Gaussian density factors leaves $|g_i(s)|\sqrt{\phi(s)}$
and $|g_j(t)|\sqrt{\phi(t)}$. Cauchy–Schwarz in each tangential
coordinate bounds the mixed term by
$C\sigma^{-1}e^{-c r_0^2/\sigma^2}\|g_i\|_2\|g_j\|_2$.
Both types of error are absorbed by the scalar bound because
the tangential Dirichlet form controls $(1-\rho)\|g_i\|_2^2$.
This proves (P11).

## Relative control in the extreme Gaussian tail

Let $D=\{A\ne T\}$ and $D_R=D\cap\{|x|\ge R\}$. For a
measurable partition $A$ write $h_R=1_{D_R}(1_A-1_T)$ for
its vector label change, let $Q_R=\langle h_R,T_\rho h_R\rangle$
be its summed interaction and let $L_R$ be its linear simplex
regret. Put $L=L_0$. Uniformly for $\rho$ sufficiently close
to one and every $R\ge0$,

\[
                    Q_R\le C e^{-c e^2R^2} L_R .
 \tag{P14}
\]

For aligned planar balanced maximizers there are also absolute
$K_0,C,c>0$ such that for $R\ge K_0/e$,

\[
 L_R \le C(e^{-c e^2 R^2}+\tau_R)L,\qquad
 \tau_R=\sigma\,\gamma_1\{|s|\ge\sqrt{(R^2-\sigma^2)_+}\}
                       +\gamma_2\{|x|\ge R\},\qquad
 \sigma=\sqrt{1-\rho^2}.
 \tag{P15}
\]

The simplex score margin for a transfer from $i$ to $j$ is
at least $c\min(t/\sigma,1)$, with $t$ the inward distance to
the pair line. Pair each $y\in T_i$ with its reflection and
write $y=\rho x+\sigma z$ in the conditional integral.
The density difference yields the factor
$1-\exp[-2\rho t(\rho t+\sigma z_N)/\sigma^2]$. The integral
contains a fixed ball in $T_i$ with $z_N\ge1/2$, by cone
addition, and the remaining integrand is nonnegative.

On conjugating Mehler to Lebesgue $L^2$, its symmetric kernel is

\[
 (2\pi\sigma^2)^{-1}
 \exp\left[-\frac{\rho|x-y|^2}{2\sigma^2}
           -\frac{1-\rho}{4(1+\rho)}(|x|^2+|y|^2)\right].
 \tag{P16}
\]

On pairs exterior to $B_R$, the second exponential is at
most $e^{-c e^2 R^2}$. The Gaussian convolution left over
has uniformly bounded $L^2$ norm. For the part of a transfer
with $t\ge\sigma$, this bounds self-interaction by the discount
times a constant times its Gaussian area, which its regret controls.
For $0\le t\le\sigma$, drop the normal factor of the convolution,
leaving a tangential convolution with norm $O(1/\sigma)$.
If $v(s,t)$ indicates the transfer in each section, then

\[
 \left(\int_0^\sigma v(s,t)\sqrt{\phi(t)}\,dt\right)^2
                     \le C\int_0^\sigma t\,v(s,t)\phi(t)\,dt .
\]

Here $\phi$ is the scalar Gaussian density, bounded above and
below on this interval, and a set of length $\ell$ in the
positive half-line has first moment at least $\ell^2/2$. After
the tangential Gaussian square-root factors, the self-interaction
is at most the discount times the regret. Sum all pieces by
Cauchy–Schwarz for the positive quadratic form to prove (P14).

Put $m_R=\gamma_2(D_R)$ and $m=m_0$. The same sections show
$m_R^2\le C\tau_R L_R$. For $t\le\sigma$, any exterior point
has $|s|\ge\sqrt{(R^2-\sigma^2)_+}$. The section bound and
Cauchy–Schwarz in $s$ bound the squared mass of this part by
$C\sigma\,\gamma_1\{|s|\ge\sqrt{(R^2-\sigma^2)_+}\}L_R$.
For $t\ge\sigma$ its mass is at most both a constant times its
regret and $\gamma_2\{|x|\ge R\}$. In particular $m^2\le CL$.

Use the fixed compact volume bumps from (P4). The centered
score prices of an aligned planar maximizer obey $|\beta|\le Cm$
uniformly near one. Indeed, the difference between its
disagreement derivative and the simplex derivative against a
fixed bump is bounded by $Cm$. In the pullback formula from
(P4), the change of the disagreement indicator requires a
changed endpoint, and the derivative weight is bounded by
$C(1+|z|^2)$ with at least one endpoint in the fixed support.
Its conditional absolute integral in the other endpoint is
uniformly bounded. For a bounded conditioning endpoint the
standardized displacement has bounded moments. For a distant
endpoint to contribute, its partner must be in the support,
and the conditional Gaussian tail dominates the polynomial.
The simplex derivative is zero. Invert the volume derivative
matrix for the maximizer to obtain the price estimate.

Replace $A$ by $T$ outside $B_R$ to obtain $A_R$, and restore
its volumes using the fixed bumps inside $B_R$. Their volume
map has the same derivatives as for $A$, and the correction
parameter is $O(m_R)$ by the uniform implicit function theorem.
The disagreement derivatives of $A_R$ differ from those of $A$
by $O(m_R)$ by the same conditional estimate. The fixed-bump
second derivatives are uniformly bounded for small parameters.
Differentiation of the pullback density gives a polynomial
bound in the deformed midpoint and standardized displacement.
Derivatives of the flow differences are bounded by a constant
times the endpoint distance. Thus there are no inverse noise
factors after standardizing the difference, and changing back
to pair variables gives uniformly finite Gaussian moments.
By stationarity of $A$ the volume repair changes disagreement
by at most $C(|\beta|m_R+m_R^2)$ in absolute value.

By the [exact deficit expansion (B1)](POSITIVE_CORRELATION_CONTINUATION.md)
the disagreement change from truncation alone is
$-2L_R+Q_R+2\langle h_0-h_R,T_\rho h_R\rangle$.
Global optimality makes the sum of this change and the repair
cost nonnegative. By (P14) and Cauchy–Schwarz for the positive
quadratic form, the mixed term is bounded by
$C\sqrt{L e^{-c e^2R^2}L_R}$, using the $R=0$ bound on the
interior transfers. The repair term is at most
$C(\sqrt{L\tau_R L_R}+\tau_R L_R)$. For $R\ge K_0/e$ with
$K_0$ large and $\rho$ sufficiently close to one, absorb the
terms linear in $L_R$. Divide by $\sqrt{L_R}$ if it is
positive and square to obtain (P15).

Ordinary Gaussian perimeter convergence suppresses the
$e^2$ corrections that this test resolves. The gap in the
[spectral theorem](SIMPLEX_INTERFACE_SYMMETRY.md) is of order $e$, and
the fixed-correlation measurable comparison has a neighborhood depending
on $\rho$. A uniform comparison still has to couple the moving
junction to the open-ray deformations, including those for which the
vertex displacement is large compared to $e$. Extending the
bounded-physical-region graph control to the radii needed by (P11)
and combining it with (P15) are further obligations.

The arguments are written analytic proofs with internal review.
