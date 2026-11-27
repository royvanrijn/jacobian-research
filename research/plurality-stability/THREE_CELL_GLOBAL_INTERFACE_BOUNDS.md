# Global interface and volume-flux bounds for three cells

The stationary scores give polynomial growth of the planar
interface length. At actual maximizers the translation index
gives bounds on volume flux and Gaussian perimeter, as well as
a stronger form of the normal-flux cut inequality.

Use the score, orientation, and weighted quadratic form
conventions of the [maximizer analysis](THREE_CELL_GLOBAL_OBSTRUCTIONS.md).
In particular $N_e$ points from cell $i$ to cell $j$, and
$d_e=e_i-e_j$ with $e=(i,j)$ denotes its label vector.

## Polynomial growth of the network

Consider a planar partition of positive cell masses with
locally finite perimeter. Assume every pair interface lies
up to surface-null sets in a stationary pair level
$F_{ij}=T_\rho(1_{A_i}-1_{A_j})-c_{ij}=0$ at $0<\rho<1$.
The union of interfaces in $B_R$ has Euclidean length
$O_{A,\rho}(1+R^4)$. For actual balanced maximizers the
constant is uniform on compact intervals of correlations.

The Gaussian convolution extends $F_{ij}$ to an entire
function of $z\in\mathbb C^2$ by using the quadratic form
without complex conjugation in the conditional density.
Taking the absolute value inside the convolution bounds it
by $C\exp(\rho^2|\operatorname{Im}z|^2/(2(1-\rho^2)))$,
uniformly in the real part. No $F_{ij}$ vanishes identically,
by injectivity of $T_\rho$ and positivity of the disjoint
cell masses. Choose three noncollinear real points $p_1,p_2,p_3$
off the zero sets of all three functions. Along any complex
line $p_l+z e$, with $e$ a real unit direction, Jensen's
formula in the disk of radius $2(R+|p_l|+1)$ bounds the number
of zeros in the disk of half that radius by $C(1+R^2)$.
The constants use the nonzero values at the three base
points and are independent of $e$. One may take a limiting
sequence of disk radii if a boundary contains a zero.

For a unit tangent $\tau$ of a rectifiable interface at $x$,
the angular Jacobian along the interface of the projection
$x\mapsto \arg(x-p_l)$ is
$|J(x-p_l)\cdot\tau|/|x-p_l|^2$, where $J$ rotates by a right
angle. The maximum of the three numerators is bounded below
by a positive constant depending only on the base points:
subtract two numerators before taking absolute values and
use that the two base-point differences span the plane.
Consequently the sum of their Jacobians is at least
$c/(1+R^2)$ on $B_R$. Coarea for each angular projection bounds
its integral by $C(1+R^2)$, since all points in a fiber lie
on the corresponding radial line and are zeros of a pair
contrast. Combining gives the length estimate.

For actual balanced maximizers at correlations in a compact
interval, take any sequence of partitions and correlations.
Strong maximizer compactness gives a subsequential limit
partition and normalized score prices converge along a
further subsequence. The prices are bounded since each
label wins on a set of positive mass and the scores lie
between zero and one. Score convergence at the three base
points chosen for the limit preserves a common positive
lower bound on all their contrast moduli. The same points
and Jensen estimates work along that subsequence, proving
uniformity by contradiction. In particular the Euclidean
surface measure is integrable under every positive quadratic
Gaussian weight, also after multiplication by a polynomial.

## A global matrix bound on volume flux

For an actual maximizer and a compact regular speed $f$, let
$Bf=\sum_e d_e\int f_e\,d\nu_e$ be its volume flux. Gaussian
Mehler factorization on the vector measure
$m_f=\sum_e d_e f_e\nu_e$ gives
\[
 \langle f,K f\rangle
   =\left\|\int M_{\sqrt\rho}(\cdot,y)\,dm_f(y)\right\|_{L^2(\gamma_2)}^2
   \ge |Bf|^2 .
\]
The inequality is the contribution from the constant
projection. Thus (G1) implies $\rho|Bf|^2\le\|f\|_a^2$.
Extend $B$ by density to $L^2(a\,d\nu)$. Testing one edge
by truncated reciprocals shows $w_e=\int a_e^{-1}d\nu_e<\infty$.
Cauchy–Schwarz with these weights also identifies the extension
of $B$ with the absolutely convergent flux integral. Consequently
\[
 \rho BB^*=\rho\sum_e w_e d_e d_e^{\mathsf T}\prec I
             \quad\text{on the label sum-zero space}.           \tag{IB1}
\]
For strictness, suppose that $\rho BB^*\lambda=\lambda$ for
some nonzero centered $\lambda$. At $f=B^*\lambda$, equality
would then hold both in (G1) and in the constant projection
bound. The signed vector surface measure $m_f$ therefore
has constant Mehler smoothing. The identity extends from
compact speeds by the weighted norm bound on $K$ in (G1)
and convergence of $m_f$ in total variation. For each fixed
real $x$, the kernel $M_{\sqrt\rho}(x,y)$ is bounded in $y$,
so the pointwise formula agrees with the Hilbert-space limit.
Injectivity on finite signed measures follows by writing the
kernel integral as a Laplace transform after Gaussian damping
in $y$, then using uniqueness of the Fourier transform.
Thus $m_f=(Bf)\gamma_2$, which contradicts its surface support
since $Bf=\lambda/\rho\ne0$.

Taking traces in (IB1) gives $\rho\sum_e w_e<1$. Summing
Cauchy–Schwarz over the interfaces, and using (G3), yields
\[
 P_\gamma(A)^2<S'_\rho(A) .                                    \tag{IB2}
\]
Here $P_\gamma(A)=\sum_e\nu_e(\Sigma_e)$ is Gaussian cluster
perimeter, with each interface counted once. At the simplex,
its square is the first Hermite weight and the bound becomes
an equality in the limit $\rho\downarrow0$.

There is also a bound on flux that retains its coupling to the
translations. Let $L v=s_v=v\cdot N$, $G=L^*L$ in the weighted
space, $C=BL$, and $P$ the weighted orthogonal projection
onto $V=\operatorname{ran}L$. The map $C$ is invertible from
the active plane onto the label sum-zero space by full
cell-moment rank. Put $\tau_f=L C^{-1}Bf$. Then $f-\tau_f$
is balanced. The ordinary form $Q_0$ is nonnegative on this
balanced speed, by compact approximation with fixed mass
bumps and continuity of the flux. By (G2) the splitting
into $V$ and $V^{\perp_a}$ is $Q_0$-orthogonal, and
\[
 Q_0(f-Pf)\ge \frac{1-\rho}{\rho}\|\tau_f-Pf\|_a^2,
 \qquad Q^{(1)}(f)\ge(1-\rho)
          [\|f-Pf\|_a^2+\|\tau_f-Pf\|_a^2].                  \tag{IB3}
\]
Since $K$ is positive semidefinite, the first inequality
implies $\|\tau_f-Pf\|_a\le\sqrt{\rho/(1-\rho)}\|f-Pf\|_a$.
The triangle inequality and weighted orthogonality then give
$\|\tau_f\|_a^2\le \|f\|_a^2/(1-\rho)$, equivalently
\[
 (1-\rho)BB^*\preceq C G^{-1} C^{\mathsf T} .                  \tag{IB4}
\]
Applied to $f=1_Ds_v$, (IB3) adds a nonnegative term to the
right-hand side of (G18). Its extra term compares the
translation with the same cell fluxes as the cut to the
weighted projection of that cut onto translations. The
first inequality of (IB3) also implies that every vector in
the unrestricted $Q_0$ operator kernel is balanced: it is
orthogonal to the translation space by (G2), hence its
balancing translation has norm zero.

## The same bounds for ordinary stable stationary partitions

Let a planar three-cell partition with positive masses have
locally finite perimeter and stationary pair scores on the
interfaces. Assume the regular interface interiors carry all
surface measure and the usual constrained second variation
$Q_0(f)=\int a f^2d\nu-\langle f,Kf\rangle\ge0$ holds for
smooth compact speeds there with $Bf=0$. Here
$a_e=-\partial_{N_e}(U_i-U_j)$ is initially signed. Assume
the interface normals span the plane in surface measure, so
no nonzero translation has speed zero almost everywhere.
Then $a>0$ almost everywhere, the ordinary form has index
exactly two on unrestricted speeds, and (IB1)–(IB4) hold
without assuming global maximality or stability in an added
coordinate.

First, on any compact part of the interior of one regular
interface, the second variation extends to bounded measurable
speeds of zero surface mean by approximation in $L^2(d\nu)$
using one fixed bump to correct the mean. The kernel and $a$
are bounded there. A positive-measure subset on which $a\le0$
would support a nonzero bounded speed of mean zero. Its kernel
quadratic form is strictly positive by the Mehler factorization
and injectivity for finite signed measures used above. Thus
second-order stability forces the positive orientation.

Polynomial length growth and bounded score gradients give
integrability of all translation terms. For the double integral,
$M_\rho(x,y)\phi_2(x)\phi_2(y)$ has a strictly positive
Gaussian quadratic decay rate in $(x,y)$, so bounded speeds
are integrable against it on the interface product. Gaussian
integration by parts therefore gives (G2) and (G3) as above.
Approximation by bounded compact speeds, with fixed independent
mass bumps, extends the balanced second variation to bounded
balanced speeds. The bumps exist since the interface graph
is connected by the finite-perimeter argument for (G1).

For a nonzero translation $s_v$, (G2) gives
$Q_0(s_v)=-(1-\rho)\|s_v\|_a^2/\rho<0$. Thus $C=BL$ is
injective by constrained stability, and it is an isomorphism
onto the centered two-dimensional label space. The negative
space of translations has dimension two, while the balanced
subspace has codimension two. This proves the index assertion.
For a compact speed $f$ we can now use its balancing
translation $\tau_f$, and the proof of (IB3) gives
$Q^{(1)}(f)\ge0$. The weighted-space extension, the strict
matrix inequality and the subsequent bounds then follow by
the same arguments. Under the planar second-variation
hypothesis, the first extra-coordinate block is therefore
nonnegative with kernel exactly the translations.

## A perimeter bound that retains the score pressures

Let $c=\phi_1(0)$, $\sigma=\sqrt{1-\rho^2}$, and take
$U_i-U_j=c_{ij}$ on the interface $e=(i,j)$.
For the maximizers and stable stationary partitions above,
\[
 \sigma S'_\rho(A)\le 2cP_\gamma(A)
  -\frac{1}{4c}\sum_{e=(i,j)}\int_{\Sigma_e}
             (c_{ij}^2+U_k(x)^2)\,d\nu_e(x),
           \qquad \{i,j,k\}=\{1,2,3\}.                       \tag{IB5}
\]
To see this, differentiate the conditional Gaussian
representation of the scores in direction $N_e$. The
largest possible negative derivative of the score difference
at the conditional masses $U_i,U_j$ is bounded by
$\rho[I(U_i)+I(U_j)]/\sigma$, where
$I(p)=\phi_1(\Phi^{-1}(p))$. For each label this follows by
selecting the corresponding upper or lower tail of one
standard normal coordinate. Since $I''=-1/I\le-1/c$ and
$I'(1/2)=0$, one has
\[
 I(U_i)+I(U_j)
  \le 2c-\frac{(U_i-U_j)^2+(1-U_i-U_j)^2}{4c}.
\]
Integrate and use (G3). Together, (IB2) and (IB5) bound the
pressure and third-score integral strictly by
$4cP_\gamma(A)(2c-\sigma P_\gamma(A))$.
