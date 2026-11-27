# Moving maximizers and the junction scale near correlation one

Balanced three-cell maximizers converge in Gaussian measure to the simplex
orbit as $\rho\uparrow1$, with a dimension-independent conclusion after
dimension reduction. Their rescaled position-and-direction measures also
converge to the simplex interface measures. The corresponding score
prices are smaller than the noise length. A compact improvement of the
Euclidean three-sector junction at the noise length would nevertheless
produce a balanced Gaussian improvement. This leaves a specific local
rescaling problem at the endpoint.

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
has exactly one. This contradiction proves the diameter bound. It
leaves the number and geometry of ties within the noise scale open.

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

Ordinary Gaussian perimeter convergence suppresses the
$e^2$ corrections that this test resolves. The gap in the
[spectral theorem](SIMPLEX_INTERFACE_SYMMETRY.md) is of order $e$, and
the fixed-correlation measurable comparison has a neighborhood depending
on $\rho$. The remaining near-one problem is to control all shrinking
defects of actual maximizers, including changes on the scale in (P5).

The arguments are written analytic proofs with internal review.
