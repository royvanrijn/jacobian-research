# Correlation contact and score derivatives for three cells

Let $A$ be an ordinary balanced three-cell Gaussian partition
and $T$ a central regular three-sector partition, cylindrical
when the ambient dimension exceeds two. Write
\[
 D_\rho=S_\rho(T)-S_\rho(A),\qquad
 L_\rho=\mathbb E\!\left[T_\rho1_{T_{I_T(X)}}(X)
                           -T_\rho1_{T_{I_A(X)}}(X)\right],
 \quad Q_\rho=\sum_i\langle e_i,T_\rho e_i\rangle,
 \quad e_i=1_{A_i}-1_{T_i}.
\]
Here the scores are evaluated at $X$ after the labels
$I_T(X),I_A(X)$ are selected. Thus $D_\rho=2L_\rho-Q_\rho$.
The regret $L_\rho$ depends on the orientation of $T$, whereas
$D_\rho$ is invariant under a change of this orientation.

## A bound on the order and number of contacts

Let $W_m(A)=\sum_i\|(1_{A_i})_m\|_2^2$ be the total Hermite
weight of degree $m$. If $A$ is outside the simplex orbit, its
first-moment loss $\delta_1=W_1(T)-W_1(A)$ is strictly positive
by the [first-moment rigidity theorem](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md).
Balance and Parseval give $\sum_{m\ge1}W_m(A)=\sum_{m\ge1}W_m(T)=2/3$.
Consequently
\[
 F_A(z)=D_z/z=\sum_{m\ge1}(W_m(T)-W_m(A))z^{m-1}
\]
is analytic in the unit disk, satisfies $|F_A|\le4/3$ there,
and has $F_A(0)=\delta_1$. Jensen's formula gives, for
$0<b<1$, the following bound on zeros in $|z|\le b$,
with multiplicities:
\[
    \#\{z:|z|\le b,\ F_A(z)=0\}
       \le \frac{\log((4/3)/\delta_1)}{\log(1/b)}.             \tag{CT1}
\]
Indeed apply Jensen at radii $R>b$ avoiding zeros on their
circles and take $R\uparrow1$. A particular positive tie at
$r$ has order at most the right side with $b=r$.

For actual nonsimplex tied maximizers on a compact positive
correlation band, the bound is uniform. The
[measurable local comparison](POSITIVE_CORRELATION_CONTINUATION.md#an-interface-gap-implies-a-uniform-measurable-local-comparison)
and the [all-correlation three-cell gap](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem)
give a uniform open neighborhood of the simplex orbit in the
active plane containing no nonsimplex tie. To obtain uniformity
cover the band by finitely many correlation neighborhoods
in (B11). By dimension reduction the maximizers can be compared
in their active plane. The first-moment distance bound (G2)
then gives a uniform positive lower bound on $\delta_1$ for
these candidates. Interior ties at correlations in an open
global-optimality interval have even contact order and positive
leading coefficient. For a first transition the maximum may
change with $\rho$, so this finite zero count for each fixed
partition leaves open the existence of a first transition.

## Pointwise score growth

For $x$ in $T_i$ write $a\ge0$ for the signed distance to
the line separating $i$ and $j$, and $b\ge-a/\sqrt3$ for
the coordinate along their incident ray. The pair scores
at $x$ differ by
\[
 G(u)=\int_0^\infty
 [\phi(v-ua)-\phi(v+ua)]\Phi(ub+v/\sqrt3)\,dv,
 \qquad u=\frac{\rho}{\sqrt{1-\rho^2}} .                       \tag{CT2}
\]
Here $\phi,\Phi$ denote the one-dimensional standard Gaussian
density and distribution. This follows by conditioning the
point $ux+Z$ on its normal coordinate. The upper restriction
of each of the two competing sectors is the same after
reflection, namely vertical coordinate at least $-|v|/\sqrt3$.
The gap satisfies
\[
 0\le uG'(u)\le\tfrac54 G(u),\qquad
 0\le \rho(1-\rho^2)L'_\rho\le\tfrac54L_\rho .                 \tag{CT3}
\]

For $b\ge0$, write $G=D(A,B)$ with $A=ua$, $B=ub$.
Differentiating the expectation of
$\operatorname{sgn}(V)\Phi(B+|V|/\sqrt3)$ for $V\sim N(A,1)$
gives
\[
 D_A=2\phi(A)\Phi(B)+\tfrac1{\sqrt3}
                          \mathbb E\phi(B+|V|/\sqrt3).
\]
This is positive and decreases with $A\ge0$, because $|V|$
is stochastically increasing in $A$ and the density factor
decreases with $|V|$. Thus $A D_A\le D$. For $A>0$, the
probability density on $v>0$ proportional to
$\phi(v-A)-\phi(v+A)$ stochastically dominates the density
$v e^{-v^2/2}$, by the increasing likelihood ratio
$\sinh(Av)/(Av)$. The expectation of $e^{-v^2/6}$ under the
latter density is $3/4$. Using $\Phi(B+v/\sqrt3)\ge\Phi(B)$
and $\phi(B+v/\sqrt3)\le\phi(B)e^{-v^2/6}$ shows
$0\le D_B/D\le (3/4)\phi(B)/\Phi(B)$. Since
$B\phi(B)\le\phi(1)<1/4$ and $\Phi(B)\ge1/2+B\phi(B)$,
we have $B D_B\le D/4$. The case $A=0$ is the zero score
gap.

For $-a/\sqrt3\le b\le0$, use the Gaussian boundary
derivative on each sector ray. Put
$c_i=(a+\sqrt3b)/2$, $c_j=(-a+\sqrt3b)/2$,
$t_i=(\sqrt3a-b)/2$ and $t_j=(-\sqrt3a-b)/2$. Then
\[
 G'=2a\phi(ua)\Phi(ub)+c_i\phi(uc_i)\Phi(ut_i)
                                       -c_j\phi(uc_j)\Phi(ut_j).
\]
All three terms are nonnegative. In $G''$ the terms from
differentiating $\Phi$ sum to $3ab\,\phi(ua)\phi(ub)\le0$,
since $c_i^2+t_i^2=c_j^2+t_j^2=a^2+b^2$ and
$c_it_i-c_jt_j=ab$. The terms from differentiating $\phi$
are nonpositive by $c_i\ge0\ge c_j$ and $a\ge0$.
Thus $G$ is nondecreasing and concave, with $G(0)=0$.
The two cases prove the pointwise bound. Integrate it on
the changed-label sets and use $\rho(1-\rho^2)u'=u$.
On compact correlation intervals differentiation under the
integral follows from the same bound and $0\le G\le1$.

Let $W_m(e)=\sum_i\|(e_i)_m\|_2^2$. At a nonsimplex tie
where $D'_\rho\ge0$ we have $Q_\rho=2L_\rho>0$ and
$Q'_\rho\le 2L'_\rho$. As the terms of degree $m\ge2$ in
$\rho Q'_\rho$ count at least twice those in $Q_\rho$, (CT3)
gives the orientation-independent necessary condition
\[
 \frac{\rho W_1(e)}{Q_\rho}
       \ge 2-\frac{5}{4(1-\rho^2)}.                          \tag{CT4}
\]
The inequality holds for every oriented reference $T$,
with weights and interaction computed in that orientation.
It is positive for $\rho<\sqrt{3/8}$. It applies to
tangential first ties and to endpoints at which globally
optimal correlations accumulate from above.

## A holomorphic estimate relative to score regret

Fix a compact interval $I=[a_0,b_0]\subset(0,1)$. There are
$\eta_I,C_I>0$ such that for every $r\in I$,
\[
 |D_z|\le C_I L_r\quad (|z-r|<\eta_I),\qquad
 |D_r^{(n)}|\le C_I n!\eta_I^{-n} L_r .                      \tag{CT5}
\]
The constants here are for the fixed ambient dimension.
For the nearest-moment alignment of $T$ they also give
$|D_z|\le C_I'\delta_1$ by (M3) of the first-moment note.

To prove the estimate, choose
$b_0/[2(1+b_0)]<\alpha<1/4$. Write $u_z=z/\sqrt{1-z^2}$
on the branch positive along $I$. Shrink $\eta_I$ so that
$u_z$ stays in a bounded set with positive real part and
$\operatorname{Re}(u_z^2)\ge c|u_z|^2$ for a fixed $c>0$.
The same bounds hold on the line segment to $u_r$.
The Euclidean Gaussian smoothing $P_i(w)$ of the sector
indicator extends to an entire function of $w\in\mathbb C^2$,
using the quadratic form without conjugation in its density.
For real $w$, differentiation and integration over its two
rays give
\[
 \nabla P_i(w)=-\sum_{\text{rays }\tau}
                      n\,\phi(w\cdot n)\Phi(w\cdot\tau),
\]
where $n$ is the outward normal and $\tau$ the outgoing unit
tangent. The formula extends to complex $w$.
At $w=u_z x$ for real $x$ the scalar entire function $\Phi$
is bounded, by integrating $\phi$ along a ray with
$\operatorname{Re}(u_z^2)\ge c|u_z|^2$. Both
$\phi(w\cdot n)$ and $(w\cdot n)\phi(w\cdot n)$ are also
bounded. Hence both the gradient and Euler derivative
$w\cdot\nabla P_i(w)$ are bounded. Integrating
$\partial_u P_i(ux)=(w\cdot\nabla P_i(w))/u$ from
$u_r$ gives boundedness of the score itself. The pair
contrast vanishes on its reflection line for all $z$, so
the gradient bound yields
$|P_i(u_z x)-P_j(u_z x)|\le C_I\min(a,1)$ with $a$ from
(CT2). At the real anchor $r$, restricting (CT2) to $v\ge u_r a$
gives $G(u_r)\ge(\Phi(2u_r a)-1/2)/2\ge c_I\min(a,1)$.
Thus $|L_z|\le C_I L_r$.

The complex Mehler exponent is
$z|x+y|^2/[4(1+z)]-z|x-y|^2/[4(1-z)]$. After shrinking
the same disks, its real part is at most
$\alpha(|x|^2+|y|^2)$. This gives an absolute bound on
$Q_z$ by a constant times the squared weighted masses of
the transfer sets. Apply the half-space mass-cost estimate
(B5), including its full-space version with finite constant,
to each directed transfer set. The preceding real score
lower bound gives $|Q_z|\le C_I L_r$. The majorants
justify holomorphy of the integrals on the disks, and the
real identity for $D$ extends there. The derivative bound
is Cauchy's estimate. In the closest-moment alignment the
real score gaps are also bounded by a constant times the
linear gap $a$. Integrating and using (M3) gives the stated
bound relative to $\delta_1$.

## The relative deficit can increase near zero

Neither monotonic direction holds universally for $D_\rho/L_\rho$.
Take a bounded positive-measure $E\subset T_i$ in the pair
coordinates of (CT2) with $a>0$ and $0<b<\phi(0)$. Exchange
its labels with its reflection across the pair line.
The resulting partition is exactly balanced. Set
$q_1=\int_E a\,d\gamma$, $q_2=\int_E ab\,d\gamma$, and
$m=q_2/q_1\in(0,\phi(0))$. The transfer first and second
Hermite weights are $8q_1^2$ and $8q_2^2$: the two transfer
labels have opposite signed moments, and the only second
degree change is in the product of the two coordinates.
Expansion of (CT2) gives
\[
 G(u_\rho)=\ell a\rho+c\,ab\rho^2+O_E(\rho^3),\qquad
 \ell=3\phi(0)/2,\quad c=3/(4\pi),\quad c/\ell=\phi(0).
\]
Indeed $D_A(0,0)=\ell$, $D_{AB}(0,0)=c$, and $D$
is odd in its first coordinate. Thus the ratio of the second
and first coefficients in $L_\rho$ is $\phi(0)m$, while in
$Q_\rho$ it is $m^2$. This makes $Q_\rho/L_\rho$ strictly
decrease near zero, so $D_\rho/L_\rho$ strictly increases
there. Choose $E$ small enough that its limiting relative
deficit at zero, $2-4q_1/\ell$, is positive. At the other
endpoint the strong Ornstein–Uhlenbeck limit gives
$D_\rho/L_\rho\to0$, as $L_\rho\to2\gamma(E)>0$.
The example therefore also excludes universal increase.

## Anisotropic dilation and mass-tangent directions

For a planar global maximizer at $r$, keep standard Gaussian
marginals and vary the symmetric cross-covariance as
$rI+\epsilon B$, where $B=B^{\mathsf T}$. On interface
$e=(i,j)$ write $g_B$ for the corresponding derivative
of $U_i-U_j$, and $b_B(x)=(Bx)\cdot N_e(x)$ for the linear
deformation speed. Use the weighted normal and score
conventions of the [maximizer analysis](THREE_CELL_GLOBAL_OBSTRUCTIONS.md).
The exact identities are
\[
 g_B=\tfrac1r\big[(Bx)\cdot\nabla(U_i-U_j)
                    -\operatorname{tr}(B\nabla^2(U_i-U_j))\big],
 \qquad
 K b_B=\tfrac{a b_B}{r^2}+\tfrac{1-r^2}{r}g_B .               \tag{CT6}
\]
For the first identity, the conditional mean varies by
$Bx$ and the conditional covariance by $-2rB$. For the
second identity, moving the cells by $I+tB$ changes the
law used to integrate their unmoved indicators by the
inverse linear map. The derivatives of that law's
conditional mean and covariance are $-rBx,-2(1-r^2)B$,
respectively. Use the stationary normal gradient to obtain
(CT6). This is the matrix version of the scalar dilation
identity (G6). The weighted polynomial integrability
in (G6) justifies these formulas for the affine speeds.

In particular (G1) implies
\[
 \int g_B b_B\,d\nu\le
          -\frac{\|b_B\|_a^2}{r(1+r)} .                      \tag{CT7}
\]
If $b_B$ is mass balanced and nonzero, it therefore detects
a nonzero mixed derivative on the constrained tangent space.
For a nonbalanced speed, let $s_v$ be its unique balancing
translation as in (G21), and $h_B=b_B-s_v$. For any speed
$h$ in the weighted space with finite pairing one has
\[
 -(1-r^2)\langle g_B,h\rangle
   =Q^{(1)}(h_B,h)+(r^{-1}-1)\langle b_B,h\rangle_a .
                                                               \tag{CT8}
\]
Here $Q^{(1)}(\cdot,\cdot)$ is the bilinear version of the
first extra-coordinate form. This follows from (CT6) and
the translation kernel. In particular the pairing extends
continuously to weighted speeds by the right side.

Suppose now the maximizer is nonsimplex with a genuine triple
junction whose three tangent lines are distinct. Every
nonzero symmetric $B$ has $b_B\ne0$. Otherwise the
distributional equation $(Bx)\cdot D1_{A_i}=0$ gives
invariance of the cells under the linear flow $\exp(tB)$.
The isolated genuine junction is fixed by this flow, and
the three tangent lines through it are invariant, since the
flow is continuous and permutes this finite set. Thus $B$
must be scalar. A nonzero scalar would give a central cone,
excluded for a nonsimplex maximizer by (G15). The mass flux
map on symmetric matrices takes values in the two-dimensional
centered label space, so has a nonzero kernel. Hence there
is always a nonzero mass-balanced linear deformation with
strict (CT7) under the stated junction hypothesis. In
particular the argument applies to a nonsimplex maximizer with
the [single central reflection](THREE_CELL_REFLECTION_REDUCTION.md).

These mixed derivative constraints are anisotropic. The
fixed-shape stability derivative for the same variation is
$r^{-1}\int a\,N^{\mathsf T} B N\,d\nu$. To see this,
differentiate the joint Gaussian density in the cross
covariance, integrate by parts once in each variable, and
use the normal translation identity (G2). The coefficient
matrix in this derivative is positive definite by full
active normal span. In particular every nonzero inward
direction $B\preceq0$ has strictly negative derivative.
The scalar first-transition comparison constrains the
mixed derivative for $B=I$ at tangential contact. It gives
no corresponding vanishing claim for general $B$ in
the mass-flux kernel.
