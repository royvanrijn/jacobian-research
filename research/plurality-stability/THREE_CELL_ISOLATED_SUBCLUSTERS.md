# Restrictions on bounded components surrounded by one cell

Take an actual planar three-cell global maximizer at $0<\rho<1$.
Write $O_j$ for the open strict-winner set of label $j$ under
the adjusted stationary scores. A bounded component $S$ of
$\mathbb R^2\setminus\overline{O_j}$ cannot be isolated from
the remaining non-$j$ cells by an arbitrarily wide region of
label $j$. In a fixed central disk, such components also have
a positive minimum diameter under the quantitative condition
below. These restrictions allow either or both non-$j$ labels
within $S$. A single-label component also needs a strictly positive
amount of exterior interaction. Angular variations impose additional
conditions on the position of either kind of component. Every
centered disk pocket of a balanced maximizer gives a third
negative direction. The circular-interface test also excludes
some translated circles and persists for pockets with boundaries
trapped in sufficiently thin annuli. The number of surrounded
components contained in a fixed disk is uniformly bounded, with
an explicit disk in which there are none.

Use the normals, the interface measure, the weight $a$, and the
quadratic forms from the [maximizer analysis](THREE_CELL_GLOBAL_OBSTRUCTIONS.md).
Write $A=(1_{A_i})_{i=1}^3$ for its indicator vector and $B$ for
the indicator vector obtained by filling $S$ with label $j$.
Set $H=A-B$ and let $D$ be the set of affected interfaces,
including the internal interfaces within $S$ and its exterior
boundary. For a bounded component as above, the boundary of
$S$ is adjacent to $j$ modulo surface-null singular points.
This follows from the local analytic stratification of the
strict-winner sets and their closures. On each affected jump
the jump of $H$ is the original label jump. Equivalently,
the identities below apply to any bounded measurable fill
with $D_vH=-\sum_{e\in D}d_e(v\cdot N_e)\,d\mathcal H^1$ in
distributions for every constant $v$. Here the sum means
restriction of each oriented edge measure to its portion in $D$.

## Covariance bound from the first extra-coordinate block

Let $f_v=1_D(v\cdot N)$. The exact trace identity is
\[
 Q^{(1)}(f_v)
 =-\sum_i\int H_i(x)
       [\partial_{vv}-(x\cdot v)\partial_v]T_\rho B_i(x)\,d\gamma_2(x).
 \tag{IC1}
\]
Indeed the differentiated convolution of $H$ gives
$d_e\cdot\partial_vT_\rho H=-\rho Kf_v$ on every original
edge. The derivative of the original score is
$d_e\cdot\partial_vT_\rho A=-a_e(v\cdot N_e)$ there.
Subtract the two derivatives, multiply by $-f_v$, and
integrate. The resulting surface integral on $D$ equals
(IC1) by Gaussian integration by parts with the derivative
of $H$. Bounded support and local finite perimeter suffice
for this calculation.

Put $\sigma^2=1-\rho^2$, $c_\rho=(1+\rho^2)/(2\rho)$ and
$\alpha_\rho=\sigma^2/(2\rho)$. Let $M$ be the total joint
Gaussian mass weighted by
\[
 1_S(x)\bigl[H(x)\cdot B(y)+1\bigr] .
\]
The weight is nonnegative, equal to zero on $B_j$ in the
second endpoint, to two when the two non-$j$ labels match,
and to one when they are different. Suppose the remaining
non-$j$ set has positive Gaussian mass. Then $M>0$ and we
can normalize the weight to a probability law, with
expectation denoted by $\mathbb E_*$. Adding the constant
one to the integrand does not affect (IC1), as the
directional differential operator annihilates constants.

Conditional on $X=x$, differentiation of the Gaussian
density of $Y$ gives the kernel factor
\[
 \frac{\rho^2}{\sigma^4}
 \left[\sigma^2|v|^2-\bigl((y-\rho x)\cdot v\bigr)^2
              +\frac{\sigma^2}{\rho}(x\cdot v)((y-\rho x)\cdot v)\right]
\]
for the negative differential operator in (IC1). Completing
the square and using $Q^{(1)}\ge0$ gives
\[
 \mathbb E_*[(Y-c_\rho X)(Y-c_\rho X)^{\mathsf T}]
 \preceq \sigma^2 I+\alpha_\rho^2\mathbb E_*[XX^{\mathsf T}].
 \tag{IC2}
\]
The difference of the two quadratic forms in a direction $v$
is $\sigma^4 Q^{(1)}(f_v)/(\rho^2 M)$. In particular (G18)
gives a lower bound on it from the two normal Gram matrices
of the interface cut. Taking a trace in (IC1) also yields
\[
 \sum_{v\in\{e_1,e_2\}}Q^{(1)}(f_v)
 =\rho\,\partial_\rho \sum_i\langle H_i,T_\rho B_i\rangle_{\gamma_2}
 \ge0,
 \tag{IC3}
\]
where the derivative varies the correlation with the cells fixed.

As one consequence, if $|x|\le R$ on $S$, the Euclidean
distance from $S$ to the remaining non-$j$ labels, with
null sets discarded, cannot exceed
\[
 (c_\rho-1)R+\sqrt{2\sigma^2+\alpha_\rho^2 R^2} .
 \tag{IC4}
\]
On the normalized law a larger separation would make
$|Y-c_\rho X|$ strictly larger than the displayed square
root, contradicting the trace of (IC2). The directional
version follows from an analogous bound on $|(Y-c_\rho X)\cdot v|$
and uses $\sigma^2$ in place of $2\sigma^2$ for unit $v$.
On a compact interval of positive correlations, both
thresholds can be evaluated at its lower endpoint since
$c_\rho-1$, $\sigma$, and $\alpha_\rho$ decrease with $\rho$.

## Minimum diameter in a central disk

Let $[\rho_-,\rho_+]\subset(0,1)$ be compact and $R>0$. For
balanced maximizers on this correlation interval, bounded
components $S$ contained in $B_R$ have a uniformly positive
minimum diameter whenever
\[
 \exp\left[-\frac{1-\rho_-}{2(1+\rho_-)}R^2\right]
       >\frac{4\rho_+^2}{3\exp(1)} .
 \tag{IC5}
\]
If $S$ uses only one label, the right side can be replaced
by $\rho_+^2/\exp(1)$.

Put $W(x)=\max_{i\ne j}(U_i(x)-\beta_i-U_j(x)+\beta_j)$ on
$S$. It is nonnegative and vanishes on the entire $\partial S$:
boundary points are limits of strict wins of $j$ and of
points outside $\overline{O_j}$. Its positive maximum is
attained in the interior, and the smooth score branch
attaining that maximum has zero gradient there. Its value
on $\partial S$ is at most zero, so the uniform score Hessian
bound gives $\|W\|_\infty\le C\operatorname{diam}(S)^2$.

Apply Gaussian divergence in each phase of $S$ to
$\nabla(U_i-U_j)$. On the outer boundary the trace is $-a$,
and the traces on the two sides of an internal interface sum
to $-a$. It follows that
\[
 \int_D a\,d\nu
  =-\sum_{i\ne j}\int_{S\cap A_i}\Delta(U_i-U_j)\,d\gamma_2
     -\int_S(2-|x|^2)W(x)\,d\gamma_2(x).
 \tag{IC6}
\]
In the drift term use $\nabla W=\nabla(U_i-U_j)$ almost
everywhere in $S\cap A_i$ and zero outer trace for $W$.
The Euclidean Laplacian of the conditional convolution is
bounded in absolute value by
$(\rho^2/\sigma^2)\mathbb E||Z|^2-2|=4\rho^2/(\exp(1)\sigma^2)$,
where $Z$ is standard planar Gaussian and the test contrast
has absolute value at most one. The moment formula follows
by splitting the mean absolute deviation of the exponential
variable $|Z|^2/2$ at one.

Now take $f=1$ on the outer boundary oriented away from $S$,
and zero on the remaining interfaces. The sign is interpreted
in the fixed edge orientation. If $\mu_i$ denotes the Gaussian
surface area of the outer boundary against phase $i\ne j$,
then
\[
 \langle f,Kf\rangle\ge
 \min_{\overline S\times\overline S}M_\rho\,
 \left(2\sum_{i\ne j}\mu_i^2+2\prod_{i\ne j}\mu_i\right)
 \ge \tfrac32 \min_{\overline S\times\overline S}M_\rho\,
                   (\mu_{i_1}+\mu_{i_2})^2 .
 \tag{IC7}
\]
Indeed the two outward label vectors have squared norm two
and mutual dot product one. If only one is present, the last
coefficient is two. The negative-index estimate from (G17)
gives $Q_0(f)\ge-(1-\rho)\|P_Vf\|_a^2/\rho$.
The projection is bounded in squared norm by
$C(\int_{\partial S}a\,d\nu)^2$ since the translation Gram
matrix has uniformly bounded inverse on the compact family
of balanced maximizers on this band, as follows from (G24)
and strong compactness. By (IC6) this cost is $O(|S|^2)$
uniformly in the displayed central disk.

Write $d=\operatorname{diam}(S)$ and let $p\in S$. The
Gaussian density and the Mehler kernel on $\overline S$
are within $1+O(d)$ of their values at $p$ and $(p,p)$,
respectively, uniformly on the band and disk. The planar
Euclidean isoperimetric inequality gives
$\operatorname{Per}(S)^2\ge4\pi|S|$. Combine this with
(IC6), (IC7), and the bound on $Q_0(f)$. After dividing by
$\phi_2(p)|S|$ and using $|S|=O(d^2)$, the result is
\[
 \exp\left[-\frac{1-\rho}{2(1+\rho)}|p|^2\right]
 \le \frac{4\rho^2}{3\exp(1)}+O_{\rho_-,\rho_+,R}(d) .
 \tag{IC8}
\]
Here $\phi_2(p)M_\rho(p,p)=(2\pi\sigma^2)^{-1}$ times the
exponential on the left. A strict inequality in the opposite
direction in (IC5) rules out a sequence with $d\downarrow0$.
When there is only one inside label, the stronger coefficient
in (IC7) gives the stated $\rho^2/\exp(1)$ threshold.

## A count of bounded surrounded components

Fix an exterior label \(j\) and take any finite collection
\(S_1,\ldots,S_N\) of bounded components of
\(\mathbb R^2\setminus\overline{O_j}\). Write
\(\Gamma_\ell=\partial^* S_\ell\), \(L_\ell=\mathcal H^1(\Gamma_\ell)\),
and \(\Gamma=\bigcup_\ell\Gamma_\ell\). Set
\[
 \Delta_\Gamma=\sup_{x,y\in\Gamma}|x-\rho y|,\qquad
 c=\frac32 .
\]
One can use \(c=2\) when all the outer interfaces meet the
same inside label. The collection satisfies
\[
 \sum_{\ell=1}^N\frac{L_\ell^2}{|S_\ell|}
     <\frac{8\pi\rho^2}{c\,e}
                \exp\!\left(\frac{\Delta_\Gamma^2}{2\sigma^2}\right).
                                                               \tag{PC1}
\]
In particular, the number of such components contained in a
disk \(B_R(p)\) is finite and obeys
\[
 N<\frac{2\rho^2}{c\,e}
       \exp\!\left(\frac{((1-\rho)|p|+(1+\rho)R)^2}{2\sigma^2}\right)
                                                               \tag{PC2}
\]
when \(N>0\). If the bound is at most one, there is no
such component. For one component the stronger shape
exclusion is
\[
 \frac{\Delta_{\partial^*S}^2}{2\sigma^2}>
       \log\!\left(\frac{c\,e}{8\pi\rho^2}
                              \frac{L^2}{|S|}\right).           \tag{PC3}
\]

The local analytic stratification of the strict winners
and their closures gives locally finite perimeter for
these sets. The open set \(\mathbb R^2\setminus\overline{O_j}\)
is locally decomposed into components of this stratification.
On a one-dimensional boundary stratum of \(S_\ell\), the
component occupies one open side. The other side has strict
wins of \(j\), since the stratum belongs to \(\overline{O_j}\).
Refine the strata to the actual adjacent winners so that
the interior side has a fixed label \(i\ne j\). These are
regular reduced pair arcs. Points of the boundary outside
these arcs have zero \(\mathcal H^1\)-measure, by local
analytic stratification and the regular-wall description.
The outer boundaries of distinct components are disjoint
in surface measure because each such boundary arc has
label \(j\) on its exterior side.
[Regular interfaces and local analytic stratification](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

Apply Euclidean divergence to \(\nabla(V_i-V_j)\) on each
inside phase \(S_\ell\cap A_i\), and sum over \(i\ne j\).
The outer normal derivative on \(\Gamma_\ell\) is \(-a\)
and the sum of the two derivatives at each internal
interface is also \(-a\). The absolute Euclidean Laplacian
of each pair contrast is at most
\(4\rho^2/(e\sigma^2)\) by the conditional Gaussian
calculation in (IC6). The bound is strict at every point
because the pair indicator contrast vanishes on the third
cell, which has positive conditional probability and
\(\big||Z|^2-2\big|>0\) almost everywhere. It follows that
\[
 \int_{\Gamma_\ell}a\,ds<\frac{4\rho^2}{e\sigma^2}|S_\ell|,
 \qquad
 \int_{\Gamma_\ell}a^{-1}ds>
                   \frac{e\sigma^2}{4\rho^2}
                             \frac{L_\ell^2}{|S_\ell|}.          \tag{PC4}
\]
The second inequality uses Cauchy–Schwarz.

On \(\Gamma\) choose the speed \(f=(a\phi_2)^{-1}\),
oriented out of the corresponding component, and extend
it by zero to the other interfaces. The reciprocal-gradient
bound (G20) and the positive minimum of \(\phi_2\) on the
compact family of boundaries put \(f\) in \(L^2(a\,d\nu)\).
Compact speeds supported in regular interface interiors
are dense in this space, and the interaction has a
continuous extension there by (G1). The weighted pairing
of \(f\) with any translation \(s_v\) is
\(\sum_\ell\int_{\Gamma_\ell}v\cdot N\,ds=0\) by Euclidean
divergence. Thus \(Q_0(f)\ge0\) by the translation-complement
index bound. This use of the index bound applies to the
weighted limit of regular speeds.
[The translation complement and finite negative index](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#negative-index-and-the-higher-hermite-blocks-at-a-maximizer)
[Reciprocal-gradient control, G20](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

Write \(d\mu=a^{-1}ds\) on \(\Gamma\). The oriented
label vectors are \(e_i-e_j\) for \(i\ne j\). Their
projection on the common direction
\((e_{i_1}+e_{i_2})/2-e_j\) is constant with squared norm
\(3/2\), and the orthogonal part contributes nonnegative
Mehler energy. For a single pair type the full squared
norm is two. Consequently \(Q_0(f)\ge0\) gives
\[
 \int_\Gamma\phi_2^{-1}\,d\mu
      \ge c\iint_{\Gamma^2} M_\rho(x,y)\,d\mu(x)d\mu(y).
\]
The factor \(\phi_2(x)M_\rho(x,y)\) is the conditional
Gaussian density at \(x\) given \(y\), hence is at least
\((2\pi\sigma^2)^{-1}\exp(-\Delta_\Gamma^2/(2\sigma^2))\).
Substitution and cancellation bound \(\mu(\Gamma)\) by
\((2\pi\sigma^2/c)\exp(\Delta_\Gamma^2/(2\sigma^2))\).
Now (PC4) gives (PC1). The planar isoperimetric inequality
and the disk bound on \(|x-\rho y|\) yield (PC2), first
for every finite subcollection and thus for the full
family. Taking one component gives (PC3).

There is a different radius bound when \(\partial S\)
is a single rectifiable Jordan curve of length \(L\).
Put \(m_2=L^{-1}\int_{\partial S}|x|^2ds\),
\(\kappa=(1-\rho)/(2(1+\rho))\), and \(P(z)=2\Phi(z)-1\).
Then
\[
          \kappa m_2>
          \log\!\left(\frac{c\,P(e\sqrt{\pi/(2\rho)})}
                            {2\sqrt\rho}\right).               \tag{PC5}
\]
Indeed use \(f=\phi_2^{-1/2}\) on the boundary in
\(Q^{(1)}\ge0\) with the same orientation and label
projection. This gives \(\int a\,ds\ge c\rho\iint
k_\rho\,ds_xds_y\) for the normalized joint kernel in
(NT4). The arclength chord estimate there applies to
any rectifiable Jordan curve. Use the first bound in
(PC4) and \(\int a\,ds<\sqrt{2/\pi}\rho L/\sigma\),
strict by the missing third label in the normal derivative
bound. The ensuing two bounds for \(e^{-\kappa m_2}\)
are \(2\sqrt\rho/[cP(\ell)]\) and
\(2\rho\sqrt{2/\pi}\,\ell/[ceP(\ell)]\) with
\(\ell=L\sqrt\rho/(2\sigma)\). Their minimum is maximized
at their intersection \(\ell=e\sqrt{\pi/(2\rho)}\),
which gives (PC5).
[The arclength chord bound, NT4](THREE_CELL_CONNECTED_NETWORK_TOPOLOGY.md#a-radius-bound-for-compact-dormant-walls)

## The interior challenger equality in a bounded cell component

Let \(S\) be any bounded connected component of one strict
winner \(O_i\), with the other labels denoted by \(j,k\).
Write \(\Gamma=\partial^*S\), \(L=\mathcal H^1(\Gamma)\) and
\(\Delta_\Gamma=\sup_{x,y\in\Gamma}|x-\rho y|\) as before.
Inside \(S\), let \(D\) consist of the regular stratum arcs
where the sign of \(V_j-V_k\) changes. Neither challenger
wins on these interior arcs. At a regular boundary point
let \(a\) be the outward slope of the actual pair in the
sense of (G1). Let \(\widetilde a\) be the outward slope
against the larger challenger on the inside, with the
same sign convention as \(a\). The slopes agree away
from dormant three-score ties. The following bound allows
either or both other labels on the exterior:
\[
 \int_D |\nabla(V_j-V_k)|\,ds+\int_\Gamma(a-\widetilde a)\,ds
 \ \ge\
   \frac{3L^2}{4\pi\sigma^2}
       \exp\!\left(-\frac{\Delta_\Gamma^2}{2\sigma^2}\right)
          -\frac{4\rho^2}{e\sigma^2}|S| .                       \tag{PC6}
\]
When the dormant ties on the outer wall have zero
arclength, a positive right-hand side requires an
interior challenger equality with positive weighted length.

Use the speed \(f=(a\phi_2)^{-1}\) on \(\Gamma\).
Its translation pairing vanishes by Euclidean divergence,
and its oriented label vectors now have a common inside
label. Their common projection still has squared length
\(3/2\). Thus the argument for (PC1) through the reciprocal
gradient bound gives
\[
 \int_\Gamma a^{-1}ds\le
       \frac{4\pi\sigma^2}{3}
                     \exp(\Delta_\Gamma^2/(2\sigma^2)).
\]
By Cauchy–Schwarz, \(\int_\Gamma a\,ds\) is at least the first
term on the right of (PC6).

For the divergence step, split \(S\) into the regions
\(S_j\), \(S_k\) where the indicated challenger has the
higher score. Integrate the ordinary Laplacian of
\(V_i-V_j\) on \(S_j\) and that of \(V_i-V_k\) on \(S_k\).
The summed outward derivative on \(\Gamma\) is
\(-\widetilde a\). If both challengers are tied on the
boundary, the one with smaller outward derivative of
its score has the larger score on the inside to first
order. The actual exterior winner has the larger
outward derivative to first order. Both slopes are
measured relative to \(V_i\), so
\(\widetilde a\le a\). If the first derivatives coincide,
the choice has no effect on the formula. On \(D\), the outward derivative from the
\(S_j\) side plus the outward derivative from the \(S_k\)
side is \(+|\nabla(V_j-V_k)|\). To see the sign, the outward
normal from \(S_j\) goes from positive to negative values
of \(V_j-V_k\), and the difference of the two integrated
contrasts is \(V_k-V_j\). The terms vanish wherever this
pair gradient is zero. Lower-dimensional junction sets
have no surface contribution by local analytic
stratification. Consequently
\[
 \int_\Gamma \widetilde a\,ds
  =\int_D |\nabla(V_j-V_k)|\,ds
        -\int_{S_j}\Delta(V_i-V_j)\,dx
        -\int_{S_k}\Delta(V_i-V_k)\,dx .
\]
Each Laplacian has absolute bound \(4\rho^2/(e\sigma^2)\)
as in (PC4). Adding the outer slope deficit yields (PC6).
These two contributions are the terms absent from the
surrounded-component area estimate.

## Exterior interaction required by a single-label pocket

Suppose $S$ contains just one label different from $j$. Orient $D=\partial^*S$
outward and set $u=1_D$ in this orientation, so $K=2M_\rho$ on $D^2$. Define
\[
 I=\iint_{D^2} K\,d\nu^2,\qquad
 I_N=\iint_{D^2} K\,N(x)\cdot N(y)\,d\nu^2,\qquad
 X_D=\sum_{v\in\{e_1,e_2\}}Q^{(1)}(1_D v\cdot N).
\]
Tracing the definition of $Q^{(1)}$ gives $X_D=\|u\|_a^2-\rho I_N$.
Let $P_V$ be the weighted projection onto translations. By (G2), the
translation part of $Q_0$ is orthogonal to the remaining part, and
\[
 Q_0(u-P_Vu)=X_D-I+\rho I_N+\frac{1-\rho}{\rho}\|P_Vu\|_a^2.
 \tag{IC9}
\]
A negative right-hand side gives a third negative direction. The trace
identity and (G17) at $m=1$ give
\[
 X_D=\rho(I-I_N)+Q^{(1)}(u)
 \ge 2\rho m_0\left(P_\gamma(S)^2-\left|\int_S x\,d\gamma_2\right|^2\right)
       +(1-\rho)\|u-P_Vu\|_a^2,\qquad
 m_0=\min_{\overline S\times\overline S}M_\rho>0.
 \tag{IC10}
\]
Here $P_\gamma(S)=\nu(D)$ and $\int_D N\,d\nu=-\int_S x\,d\gamma_2$.
The difference in parentheses is
$\iint_{D^2}(1-N(x)\cdot N(y))d\nu^2$. It is bounded below by
$(\min_{\overline S}\phi_2)^2\operatorname{Per}(S)^2>0$, since
$\int_D N\,d\mathcal H^1=0$ by Euclidean divergence and $S$ is bounded
with positive area.

Let $\delta$ be the essential infimum of $|y-c_\rho x|$ on the
interaction-weighted law of (IC2), and suppose $S\subset B_R$. Taking
the trace of the identity behind (IC2) and applying (IC10) gives
\[
 M(2\sigma^2+\alpha_\rho^2R^2-\delta^2)
 \ge \frac{\sigma^4}{\rho^2}X_D
 \ge \frac{2\sigma^4 m_0}{\rho}
     (\min_{\overline S}\phi_2)^2\operatorname{Per}(S)^2>0.
 \tag{IC11}
\]
This gives a positive minimum exterior interaction for each fixed
pocket geometry. The general sign in (IC9) still requires control of
$X_D$ and the translation projection.

## Angular variations of an enclosed component

Let $J$ be counterclockwise rotation by $\pi/2$, put $b=\rho/\sigma^2$,
and define
\[
 W_{\mathrm{in}}(x,y)=H(x)\cdot H(y),\qquad
 W_{\mathrm{out}}(x,y)=H(x)\cdot B(y)+1_S(x).
\]
On $S$, $H$ has the form $e_i-e_j$ for $i\ne j$. The weights are
nonnegative, supported respectively on $S^2$ and on $S\times E$, where
$E$ is the set of non-$j$ labels remaining after filling $S$. On each
support the weight is two for matching non-$j$ labels and one for
distinct labels.

For any smooth function $\eta$ put $\eta_x=\eta(|x|^2)$. The speed
$f_\eta=1_D\eta_x Jx\cdot N$ has zero first variation of all masses.
It satisfies the following two identities, with $d\gamma^2$ denoting
the product of the planar Gaussian measures:
\[
\begin{split}
0={}&\iint\left[\tfrac12 W_{\mathrm{in}}(\eta_x-\eta_y)
                    +W_{\mathrm{out}}\eta_x\right]
            (Jx\cdot y)M_\rho(x,y)\,d\gamma^2,\\
0\le Q_0(f_\eta)={}&b\iint\left[\tfrac12 W_{\mathrm{in}}(\eta_x-\eta_y)^2
                    +W_{\mathrm{out}}\eta_x^2\right]
       [x\cdot y-b(Jx\cdot y)^2]M_\rho(x,y)\,d\gamma^2.
\end{split}
\tag{IC12}
\]

To prove this, write $R_x=Jx\cdot\nabla_x$ and $L_x=\eta_x R_x$.
The vector field of $L$ has zero Gaussian divergence. The distributional
jump formula gives $LH=-\sum_{e\in D}d_e f_\eta\,d\mathcal H^1$.
Integrating against Gaussian density yields zero in every label, as
$H$ has bounded support. This proves balance. Approximate the bounded
speed on $D$ by smooth compact speeds on the regular interfaces in
$L^2(a\,d\nu)$ and in $L^1(d\nu)$ on a fixed compact disk. These
approximations use local finiteness of the interface measure and
restriction to regular patches. Subtract two fixed regular mass
bumps to make the approximations balanced. The correction coefficients
tend to zero by flux convergence. Now (G1) ensures continuity of $Q_0$
in the weighted norm, so ordinary stability applies to $f_\eta$.

The first variation is zero because the stationary pair-score constants
are differences of cell prices. With Gaussian integration by parts this
gives $\langle H,LT_\rho A\rangle=0$. The same integration by parts with
the normal-score equation gives $\|f_\eta\|_a^2=-\langle H,L^2T_\rho A\rangle$.
Integrating the two derivatives of $H$ against Mehler's kernel yields
\[
 \langle f_\eta,Kf_\eta\rangle
  =\iint W_{\mathrm{in}}\eta_x\eta_y R_xR_y M_\rho\,d\gamma^2.
\]
Simultaneous rotation invariance gives $R_y M_\rho=-R_x M_\rho$.
The kernels $R_x M_\rho$ and $R_x^2 M_\rho$ are respectively
antisymmetric and symmetric in $x,y$, and $L_x^2=\eta_x^2R_x^2$.
Substitute $A=B+H$ into the pairings with $T_\rho A$ and symmetrize
the $H$ terms. Adding $1_S(x)$ to each exterior term integrates to zero
since $\int M_\rho(x,y)d\gamma_2(y)=1$. Use
$R_x M_\rho=b(Jx\cdot y)M_\rho$ and
$R_x^2 M_\rho=b[b(Jx\cdot y)^2-x\cdot y]M_\rho$ to obtain (IC12).
The derivatives and integrals are justified by bounded support in the
differentiated pocket variables and the Gaussian conditional moments.

If $E$ has positive mass, $\eta=1$ gives conditions on the same
normalized law as in (IC2):
\[
 \mathbb E_*[JX\cdot Y]=0,\qquad
 \mathbb E_*[X\cdot Y]\ge b\,\mathbb E_*[(JX\cdot Y)^2]>0.
 \tag{IC13}
\]
The signed angular displacement must have both signs. In particular,
it is impossible that $x\cdot y\le b(Jx\cdot y)^2$ almost everywhere on
$S\times E$, which also excludes pair angles all at least $\pi/2$.
The equality set has product measure zero, so the displayed pointwise
condition makes the balanced speed in (IC12) strictly unstable. If
$E$ has zero mass, the inequality in (IC12) still applies to a radial
cut $\eta=1_{\{r_1<|x|<r_2\}}$ by bounded approximation in its integral.
It excludes the same pointwise condition across the cut if both sides
have positive pocket mass.

For a finite rotation $R_\theta$ satisfying
$\gamma_2(R_\theta S\cap E)=0$, the partition $B+H\circ R_{-\theta}$
has the original Gaussian masses. The self-interaction of $H$ is
rotation invariant. Thus the stability change is twice the change
in its interaction with $B$. If $E$ has positive mass, maximality gives
\[
 \mathbb E_*\exp\!\left(b\big[(\cos\theta-1)X\cdot Y+
                             \sin\theta\,JX\cdot Y\big]\right)\le1 .
 \tag{IC14}
\]
The exponential is the ratio of rotated to original Mehler kernels.
The $1_S$ part of the cross weight integrates to $\gamma_2(S)$ at
both angles, so it cancels in the comparison.

## Circular interfaces and pockets in a thin annulus

At a planar maximizer, a whole circle $D=p+\partial B_R$ with
$R>0$ cannot separate the same two labels with consistent
orientation under either condition
\[
 \begin{array}{ll}
 R\le1,&e^{-\kappa_\rho |p|^2}\ge\rho,\\
 R^2\le3,&e^{-\kappa_\rho |p|^2}\ge\rho^{1-R^2/3},
 \end{array}
 \qquad
 \kappa_\rho=\frac{1-\rho}{2(1+\rho)}.
 \tag{IC15}
\]
The first condition includes all $|p|\le2$ when $R\le1$.
Indeed $\log(1/\rho)>2(1-\rho)/(1+\rho)$ for $0<\rho<1$.
The difference vanishes at one and its derivative is
$-(1-\rho)^2/[\rho(1+\rho)^2]$. The second condition includes
$|p|^2/4+R^2/3\le1$.

Orient $D$ outwards from its bounded side and let $F=U_i-U_j$,
where $i$ is the inside label along $D$ and $j$ the outside
label. The score-gradient identity gives $a=-\partial_N F$ almost
everywhere on $D$. Use the speed $f=\phi_2^{-1/2}$ on $D$
and zero elsewhere, with its sign changed if the fixed edge
orientation is the reverse. Its local energy is $\int_D a\,ds$,
where $ds=d\mathcal H^1$. This compactly supported bounded speed
is admissible by the extension of (G1).

Put $\sigma^2=1-\rho^2$. Averaging $F(p+r\theta)$ over
unit directions $\theta$ radializes the input contrast about
$\rho p$, because the conditional Gaussian mean is
$\rho p+\rho r\theta$. If $h(s)$ is that radialized contrast,
then $|h|\le1$ and the average is $\int_0^\infty h(s)q_r(s)\,ds$,
with conditional radius density
\[
 q_r(s)=\frac{s}{\sigma^2}
 e^{-(s^2+\rho^2 r^2)/(2\sigma^2)} I_0(\rho rs/\sigma^2).
 \tag{IC16}
\]
Here $I_0(z)=\int_0^{2\pi}e^{z\cos\theta}\,d\theta/(2\pi)$ and
$I_1=I_0'$. For $b=I_1/I_0$, differentiation gives
$b'=\operatorname{Var}_z(\cos\theta)>0$, using the probability
density proportional to $e^{z\cos\theta}$. Periodic integration
of $(\sin\theta\,e^{z\cos\theta})'$ gives
$b'=1-b/z-b^2$ for $z>0$. At $z_R=\rho R^2/\sigma^2$,
this implies $b(z_R)<\rho$: the expression $1-u/z_R-u^2$
decreases for $u\ge0$, is positive at $u=b(z_R)$, and is
nonpositive at $u=\rho$ when $R\le1$.

At $r=R$, the derivative $\partial_r q_r(s)$ changes sign
from negative to positive at the unique $t>0$ given by
$t b(\rho Rt/\sigma^2)=\rho R$. This follows by differentiating
(IC16), since $s b(\rho Rs/\sigma^2)$ increases strictly from
zero to infinity. Also $t>\rho R$ since $b<1$, and $t>R$
if $R\le1$ by the preceding bound. The same differentiation, with
$I_1'=I_0-I_1/z$, gives
\[
 \partial_r q_r(s)=-\frac{\rho}{\sigma^2}\,
   \partial_s\!\left[s e^{-(s^2+\rho^2r^2)/(2\sigma^2)}
                                      I_1(\rho rs/\sigma^2)\right].
\]
Thus the radial derivative of the average contrast is bounded
in absolute value by
\[
 \int_0^\infty |\partial_r q_r(s)|_{r=R}\,ds
  =\frac{2\rho^2 R}{\sigma^2}
       e^{-(t^2+\rho^2R^2)/(2\sigma^2)} I_0(\rho Rt/\sigma^2)
  =2\rho^2 q_R(R)\frac{F_0(t)}{F_0(R)},\qquad F_0(s)=q_R(s)/s .
 \tag{IC17}
\]
The boundary terms at zero and infinity vanish and the equality
uses the equation for $t$. Differentiation is justified by
Gaussian domination locally in $r$. The function $F_0$ is
strictly decreasing on $s\ge R$, since its logarithmic derivative
is $(\rho R b(\rho Rs/\sigma^2)-s)/\sigma^2<0$. Thus $R\le1$
gives $F_0(t)/F_0(R)<1$.

For arbitrary $R>0$, we also have $F_0(t)/F_0(R)<\rho^{-R^2/3}$.
It suffices to consider $t<R$. Monotonicity of $b$ and the
equation for $t$ give
\[
 \begin{split}
 \log\frac{F_0(t)}{F_0(R)}
 &\le \frac{1}{\sigma^2}\int_t^R (s-\rho^2R^2/t)\,ds\\
 &=\frac{R^2}{\sigma^2}
      \left[\frac{1-\lambda^2}{2}-\frac{\rho^2}{\lambda}+\rho^2\right],
 \qquad \lambda=t/R.
 \end{split}
\]
The bracket is maximal at $v=\rho^{2/3}$ on $(0,\infty)$,
where it equals $(1-v)^2(1+2v)/2$. After division by
$\sigma^2=1-v^3$ this is strictly less than
$\tfrac12\log(1/v)=\tfrac13\log(1/\rho)$, since
$(1+v)(1+2v)<2(1+v+v^2)$ and
$\log(1/v)\ge 2(1-v)/(1+v)$.

For the interaction let $G_\rho$ be the joint Lebesgue density
and $k_\rho(x,y)=G_\rho(x,y)/\sqrt{\phi_2(x)\phi_2(y)}$.
On a centered circle, $\int_{\partial B_R}k_\rho(x,y)\,ds_y
=q_R(R)$ for every $x\in\partial B_R$, since the two marginal
densities are equal there. The same-pair label factor is two.
Write $J_\circ=2(2\pi R)q_R(R)>0$ for the centered interaction.
The joint Gaussian density gives
\[
 k_\rho(p+x,p+y)
  =k_\rho(x,y)e^{-\kappa_\rho(|p|^2+p\cdot(x+y))}.
 \tag{IC18}
\]
Under simultaneous inversion of $x,y$ on the centered circle,
$k_\rho$ is unchanged. Averaging the two integrals replaces the
factor $e^{-\kappa_\rho p\cdot(x+y)}$ by a hyperbolic cosine,
so $\langle f,Kf\rangle\ge e^{-\kappa_\rho|p|^2}J_\circ$.
By (IC17) the local energy is at most
$\rho^2 J_\circ F_0(t)/F_0(R)$.
Either condition in (IC15) gives $Q^{(1)}(f)<0$, contrary to
(G1). More explicitly, let $h=f-P_V f$. By (G2) translations
are in the bilinear kernel of $Q^{(1)}$, so
$Q^{(1)}(h)=Q^{(1)}(f)<0$. Then
$Q_0(h)=[Q^{(1)}(h)-(1-\rho)\|h\|_a^2]/\rho<0$.
The two translations and $h$ span a negative three-dimensional space.

If an enclosed component is a disk, its boundary carries one pair
even when both non-enclosing labels occur deeper inside. The
analytic strict-winner closures stratify its circle into
finitely many arcs and exceptional points. Each arc has the enclosing
label $j$ outside and one of the other labels inside. If both types
occur, choose a point where the inside label changes from $i$ to
$k$. Both adjusted contrasts against $j$ vanish on the entire
circle by analyticity. The tangential derivative of $V_i-V_k$
is zero. Its inward normal derivative is nonnegative on the $i$
arcs and nonpositive on the $k$ arcs by the winning-score rule.
By continuity the two gradients coincide at the change point. The
three winners accumulate there, contrary to the distinct gradients
at genuine triple junctions proved after (G20). Thus one pair
covers the circle with consistent orientation. For a balanced
maximizer a centered disk therefore satisfies (IC15). Its mass
is at most $2/3$ since the surrounding cell has mass $1/3$,
so $R^2\le 2\log 3<3$.

The exclusion persists for pockets whose whole outer boundary is
trapped in a sufficiently thin annulus. Fix a compact set of
parameters $(\rho,p,R)$, with $0<\rho<1$ and $R>0$, satisfying
one of the conditions in (IC15) at each parameter. There is a
uniform $\varepsilon>0$ with the following property. Let $S$ be a
bounded set of finite perimeter and $D=\partial^* S$. Suppose
\[
 B(p,R-\varepsilon)\subset S\subset B(p,R+\varepsilon)
\]
up to volume-null sets, and $D$ lies almost everywhere on one
pair interface with its normal from $i$ to $j$ agreeing with
the outer normal of $S$. This configuration cannot occur at
a maximizer.

Indeed, the same speed on $D$ has local energy
$-\int_S\Delta F\,dx$ by Gauss–Green. Its difference from the
flux $-\int_{B(p,R)}\Delta F\,dx$ is bounded by
$4\pi R\varepsilon\|\Delta F\|_\infty$, which is uniformly
$O(\varepsilon)$ by the convolution derivative bound.
The virtual circular flux is bounded by (IC17), regardless of
whether the comparison circle is itself an interface. For the
interaction, write $r=|x-p|$, $e_r=(x-p)/r$, and let
$g\ge0$ be smooth on the unit circle. Apply Gauss–Green to
$g(e_r)e_r/r$, with a radial cutoff equal to one on the
annulus and zero near $p$ and beyond the exterior ball.
Away from the cutoffs its divergence is zero. The inner
cutoff changes from zero to one inside the included ball and
contributes $\int g\,d\theta$ to the integral of divergence
over $S$. The outer cutoff contributes zero. It follows that
\[
 \int_D g(e_r)\,ds
 \ \ge\ (R-\varepsilon)\int_D g(e_r)(N\cdot e_r)/r\,ds
 \ =\ (R-\varepsilon)\int g\,d\theta.
\]
Thus the radial pushforward of $ds$ dominates
$(R-\varepsilon)d\theta$. On $D\times D$, the positive
kernel $k_\rho$ is at least $1-O(\varepsilon)$ times its
value at the respective projections onto $p+\partial B_R$.
This follows from the uniform bound on the derivatives of
$\log k_\rho$ on the compact spatial and correlation range.
The two projected measures and the positive lower bound on $R$ give an
interaction at least $1-O(\varepsilon)$ times the circular
interaction. The resulting circular upper bound
$\rho J_\circ[\rho F_0(t)/F_0(R)-e^{-\kappa_\rho|p|^2}]$
for $Q^{(1)}$ is continuous and strictly negative at every
parameter in the compact set. Here $t$ is continuous by
strict monotonicity in its defining equation. Taking $\varepsilon$
sufficiently small therefore again gives $Q^{(1)}(f)<0$.

### Mixed inside labels in a thin annulus

For balanced maximizers the thin-annulus exclusion also holds with
a common outside label $j$ and either other label on the inside
of each reduced boundary arc. Suppose it fails for a sequence with
annulus widths tending to zero and parameters in a compact set as
above. Pass to convergent parameters $(p,R,\rho)$ and relabel so
the outside label is fixed. Strong maximizer compactness gives a
limiting actual balanced maximizer. The normalized prices have a
convergent subsequence because price differences are bounded by
the score range. The scores and their spatial gradients converge
uniformly on compact sets by the strong convergence and the
conditional Gaussian kernel bounds on compact correlation
intervals. All limiting pair contrasts are nonconstant by
injectivity of $T_\rho$, so ties have zero Gaussian measure.

Write $F_i=U_i-U_j-\beta_i+\beta_j$ for the two limiting
adjusted contrasts, and denote their prelimit versions by $F_i^n$.
The angular flux estimate for the annulus shows that $ds$ on
the prelimit reduced boundary projects to at least
$(R_n-\varepsilon_n)d\theta$. Almost everywhere on that boundary
the outside label $j$ ties with one of the two remaining labels,
and both contrasts are nonpositive. Consequently the same holds
at every point of the limiting circle by uniform score convergence
and positive projected boundary measure on every angular interval.
Use (G20) and a lower bound for the Gaussian density on the
bounding disk to obtain a uniform bound
$\int_{\partial^*S_n}a^{-1}ds\le C$. Suppose that on an angular
interval $L$ the only possible active pairs for large $n$
are among $(i,j)$ for $i$ in a set of labels $E$. Along
the corresponding prelimit boundary, the local $a$ is bounded
above by $m(\theta)+o(1)$, uniformly, where
\[
 m(\theta)=\max_{i\in E}|\nabla F_i(p+R\theta)|.
\]
This follows from the uniform score-gradient convergence and
the annular trapping. The angular measure bound and Fatou's
lemma yield $\int_L m(\theta)^{-1}\,d\theta\le C/R$. The
function $m$ is Lipschitz, so it cannot vanish at any interior
point of $L$.

At least one of the two restricted contrasts vanishes
identically on the circle, since otherwise their analytic zero
sets on the circle would be finite and could not cover it.
If only $F_i$ vanishes identically, the other contrast is
strictly negative away from finitely many points. On every
closed subarc avoiding those points, take $E=\{i\}$ in the
preceding bound. Then $\nabla F_i\ne0$ there and the limit
circle is an active $i,j$ interface on that subarc. At the
finitely many exceptional points the nonvanishing result for
closure points of active pair interfaces again gives
$\nabla F_i\ne0$. Thus the circle carries the same oriented
pair almost everywhere.

If both contrasts vanish identically on the circle, take $E$
to consist of both labels. At every angle the three adjusted
score derivatives in the inward normal direction have distinct
maximum and minimum values, since $m>0$, and tangential
derivatives of their contrasts vanish. In the signed inward
distance $t$ from the circle, for each pair of labels the
first nonidentically zero Taylor coefficient of its contrast
is an analytic function of the angle. Such a coefficient exists
because the pair contrast cannot vanish in an open set. Exclude
the finitely many angular zeros of these three coefficients.
On each remaining arc, at each angle and for sufficiently small
$t>0$ or $t<0$, the winner is a single label fixed on that
side of the arc. The inside and outside winners differ, since
the maximum and minimum first derivatives differ.

If the winner on one side changes at an exceptional angle,
the two changing labels have the same extremal first derivative
there. The third has a strictly separated derivative, hence
loses in a neighborhood on that side for all sufficiently
small nonzero $t$. In arbitrarily small connected collar
neighborhoods of the exceptional point, both changing labels
have strict-winner sets of positive area. Their reduced pair
interface therefore accumulates there, since otherwise their
indicators would have zero distributional gradient and be
constant in such a neighborhood. But their contrast has zero
tangential and normal derivatives at that point, again contrary
to nonvanishing on the closure of an active pair interface.
No change of winners along either side is possible. The limit
circle is therefore one consistently oriented pair interface
in this case as well.

In either case the limiting actual maximizer has an interface
excluded by (IC15). This proves the uniform thin-annulus
exclusion for a fixed outside label.
