# The simultaneous first angular moment envelope for arbitrary circular masses

For any prescribed circular mass profile, the norm caps and closure are
the full feasible set of first moment triples. Every feasible field has
a measurable realization with at most thirteen intervals in each circle
cut at angle zero. Consequently an integrable radial kernel admits an
exact moment optimization with one common moment field across all radii.
This is an internal analytic result.

## Feasible moments on a single circle

Use $d\sigma=d\theta/(2\pi)$ and write $a_i(p)=\sin(\pi p_i)/\pi$.
The feasible first moments at prescribed $p\in\Delta_3$ are precisely
\[
 \mathcal K_p=\{(u_1,u_2,u_3)\in(\mathbb R^2)^3:
                  \sum_i u_i=0,\quad |u_i|\le a_i(p)\} .         \tag{SM1}
\]
Necessity is the cap inequality and closure. For sufficiency first
allow fractional labels, measurable $g_i\ge0$ with $\sum_i g_i=1$.
The feasible moment set is then compact and convex by weak compactness
in $L^2(\sigma)^3$. Each extreme point of $\mathcal K_p$ is realized by
ordinary labels as in [GAC2](THREE_CELL_GLOBAL_ANGULAR_COMPATIBILITY.md):
two caps are saturated in the interior mass case, and the two centered
arcs with those moments are disjoint. The complement has the required
third moment. On the two-label boundary use complementary arcs for
the extreme points of the moment disk. On the one-label boundary the
moment body is a singleton. Taking convex hulls shows that every point
of $\mathcal K_p$ is fractionally feasible.

Elimination of fractional labels for finitely many moments is a case of
the atomless setting in [Dvoretzky, Wald and Wolfowitz, Theorem 4](https://msp.org/pjm/1951/1-1/pjm-v1-n1-p06-s.pdf).
Here there is also a direct argument that bounds the number of intervals.
Let $C$ be the set of all fractional masses and first moments together,
allowing the masses to vary. This is a compact convex set in finite
dimensions. At a point of its relative interior, maximize the additional
linear functional $\sum_i i\int g_i(\theta)\cos(2\theta)d\sigma$ over
profiles with the given moments. The value as a function on $C$ is
finite and concave, so it admits a supporting affine majorant at this
point. Subtract its linear part from the auxiliary functional. The
chosen fractional maximizer then maximizes the resulting linear
functional with no moment constraints. Its three scores are linear
combinations of $1,\cos\theta,\sin\theta$ plus $i\cos(2\theta)$.
Every pair contrast therefore has at most four zeros. Off the at most
twelve zeros the winner is unique. Thus the fractional maximizer is
ordinary almost everywhere, piecewise constant between the zeros.

For boundary points, take convex mixtures with a fixed relative
interior point and pass to a limit in the interval parameters. With
a cut at zero, the ordinary realizers have at most thirteen intervals.
Their labels range over a finite set of sequences, and the ordered
endpoints lie in a compact simplex, allowing zero-length intervals.
The map from these parameters to masses and moments is continuous.
Its image is compact, so every boundary point also has an ordinary
realizer of the same complexity. This proves (SM1) with the claimed
interval bound.

The construction permits simultaneous measurable realization. Use the
compact interval parameter space just described. At level one cover
it by finitely many compact sets of diameter at most $1/2$, and at
each next level refine each compact set by a finite compact cover of
diameter at most $2^{-n}$. For each prescribed joint mass and moment
vector, choose the first set in the current cover whose image contains
the vector, starting from the full space. The tests are Borel because
all the images are compact. The unique point in the nested sets is
therefore a Borel selection from the fiber. Composing with a measurable
field of masses and moments gives measurable endpoints and hence an
ordinary planar partition with these exact data.

## Exact shared envelope and finite radial approximations

For $W\in L^1(\nu\otimes\nu)$, even signed, define
\[
 M_W(p)=\max_{u(t)\in\mathcal K_{p(t)}}
       \iint W(t,s)\sum_i u_i(t)\cdot u_i(s)d\nu(t)d\nu(s).       \tag{SM2}
\]
The maximum has exactly the value of the original first moment problem
over ordinary circular partitions with masses $p(t)$. Indeed (SM1)
gives the two feasible directions, and the bounded closed convex set
of feasible fields is weakly compact in $L^2(\nu;(\mathbb R^2)^3)$.
Approximate $W$ in $L^1$ by finite linear combinations of rectangle
indicators. On the bounded fields the objective is a uniform limit
of weakly continuous quadratic functions, since every rectangle
interaction uses only the integrated moments on its two sides.
It follows that the maximum exists, and the selection above gives
an ordinary maximizer. With two equal masses, [EP1](THREE_CELL_EQUAL_PAIR_ANGULAR_MOMENTS.md)
reduces this problem further when the kernel is nonnegative.

The objective is stable both in the kernel and the profile:
\[
 |M_W(p)-M_{\widetilde W}(q)|\le\frac3{\pi^2}\|W-\widetilde W\|_1
   +\frac1\pi\int \bigg(\int \big[|W(t,s)|+|W(s,t)|\big]d\nu(s)\bigg)
                              |p(t)-q(t)|_1d\nu(t).             \tag{SM3}
\]
For the kernel term use the bound $|u_i|\le1/\pi$. For the profile term,
starting from an ordinary realizer of $u\in\mathcal K_p$, reassign surplus
mass along each circle to the deficient labels. This uses angular measure
$|p-q|_1/2$ and changes the sum of moment distances by at most
$|p-q|_1$. The reassignment can be done measurably: order the labels
and use intervals starting at the angle-zero cut to select the required
amounts from the surplus sets. Their cumulative masses are continuous
in the cut and measurable in radius, so the selected quantiles are
measurable as well. Change one argument of the product at a time and
use $|u_i|\le1/\pi$ again. Taking the maximum over profiles in both
directions proves (SM3).

If $W$ and $p$ are constant on the rectangles and cells, respectively,
of one finite radial partition, we can average $u$ on each cell. This
preserves feasibility by convexity and preserves every interaction.
The problem is then exactly a finite quadratic optimization with
three Euclidean caps and a closure constraint per cell. Conditional
averages on refining radial intervals converge to $p$ almost everywhere
and to $W$ in $L^1$ for the product partitions. Thus (SM3) gives convergence
of these finite optimizations to the exact simultaneous comparison for
each specified profile. This replaces the use of the separate pair
envelope as an upper bound, including at large amplitudes. It does not
control the sign of the comparison when the profile is allowed to vary,
nor does it capture higher angular modes by itself.

## The exact angular Fourier head with a uniform tail bound

Let $P_n$ be the orthogonal projection onto angular Fourier orders
$0,1,\ldots,n$ in the plane. This commutes with the Ornstein–Uhlenbeck
operator $T_\rho$ because $T_\rho$ commutes with rotations. The space
orthogonal to $P_n$ is also orthogonal to every polynomial of degree
at most $n$: in polar coordinates such a polynomial has angular
frequencies of magnitude at most its degree. Thus the Hermite norm
of $T_\rho$ on the complementary space is at most $\rho^{n+1}$.
Every ordinary balanced planar partition consequently satisfies
\[
 0\le S_\rho(A)-\sum_i\langle P_n1_{A_i},T_\rho P_n1_{A_i}\rangle
   \le \rho^{n+1}\sum_i\|1_{A_i}-P_n1_{A_i}\|_2^2
   \le \tfrac23\rho^{n+1}.                                      \tag{SM4}
\]
For the last inequality, the projection contains constants, and the
sum of the three squared centered norms is $2/3$.

The projected objective has an exact finite-angular-moment formulation.
For a slice retain the masses $p_i$ and $u_{i,j}=\int_{A_i}
(\cos j\theta,\sin j\theta)d\sigma$ for $1\le j\le n$. The full
body $C_n$ of such data, including the masses as variables, is compact
and convex. Repeat the proof of (SM1) with all the moment constraints
and an auxiliary $\sum_i i\int g_i\cos((n+1)\theta)d\sigma$.
Every pair contrast now has at most $2(n+1)$ roots. This shows that
the fractional convex body equals the ordinary moment body and every
point has an ordinary realizer with at most $6(n+1)+1$ intervals after
cutting at angle zero. The same measurable selection applies. The
support of $C_n$ in any specified direction is the ordinary integral
of the largest of three prescribed trigonometric polynomials of order
at most $n$, by the pointwise linear assignment. This characterizes
$C_n$ without angular endpoint choices in the optimization.

Put $d=1-\rho^2$, $v=2\rho\sqrt{ts}/d$, and
$W_{\rho,j}(t,s)=d^{-1}e^{-\rho^2(t+s)/d}I_j(v)$, with $I_j$ the
modified Bessel function. The projected objective in (SM4) is
\[
 \iint W_{\rho,0}(t,s)\sum_i p_i(t)p_i(s)d\nu^2
  +2\sum_{j=1}^n\iint W_{\rho,j}(t,s)
                 \sum_i u_{i,j}(t)\cdot u_{i,j}(s)d\nu^2 .        \tag{SM5}
\]
Maximizing (SM5) over one measurable field in $C_n$ with the prescribed
profile, or over all fields in $C_n$ with integrated masses $1/3$,
therefore gives a lower bound on the corresponding true maximum that
differs from it by at most $(2/3)\rho^{n+1}$. Compactness and attainment
follow as for (SM2), since all the moment components are bounded and
all the kernels belong to $L^1$.

For the optimization over all balanced profiles, replace the kernels
by their conditional averages on a common finite radial partition.
Averaging a field on each cell preserves $C_n$, the integrated masses,
and every step-kernel interaction. This gives a finite quadratic
optimization over one copy of $C_n$ per cell and linear mass constraints.
The uniform error in (SM5), over all feasible fields, is at most
\[
 \|W_{\rho,0}-W_{\rho,0}^{\rm step}\|_1
       +2\sum_{j=1}^n\|W_{\rho,j}-W_{\rho,j}^{\rm step}\|_1.       \tag{SM6}
\]
Indeed $\sum_i p_i(t)p_i(s)\le1$ and $|u_{i,j}(t)|\le p_i(t)$.
For each fixed $n$, the error tends uniformly to zero as the radial
partition refines, on every compact correlation interval. Relative
to $dt\,ds$ the bound $I_j(v)\le e^v$ gives an integrable envelope
$(1-b^2)^{-1}e^{-(t+s)/(1+b)}$ for $\rho\le b<1$.
The kernels are jointly continuous in $\rho,t,s$ on compact radial
boxes. This proves $L^1$ continuity in $\rho$, and uniform convergence
then follows from compactness in $L^1$ and contraction of conditional
expectation. For a concrete mesh, one can also bound the error by the
kernel oscillations in the box and its integral outside the box.

The errors in (SM4) and (SM6) are additive. The finite optimizations
therefore approximate the all-profile Gaussian optimum uniformly
on compact intervals. An additive approximation does not establish
exact optimality at the simplex. A comparison separating the local
simplex neighborhood from possible distant maximizers is still needed.

## A finite strict comparison outside the local simplex neighborhood

Fix a compact interval $I=[a,b]\subset(0,1)$ and a number $\delta>0$
such that the measurable local simplex theorem applies whenever
$d(A)<\delta$, uniformly on $I$. Here
$d(A)=\min_{T'}\sum_i\gamma_2(A_i\triangle T'_i)$ and $T'$ ranges
over the simplex orbit. Such a number exists by the
[measurable local comparison and maximizer compactness](THREE_CELL_GLOBAL_OBSTRUCTIONS.md).
Given any admissible $\delta$, the finite moment optimizations above
can be restricted by linear constraints to cover every partition with
$d(A)\ge\delta$. This supplies a sufficient finite strict comparison
for global optimality throughout $I$. Conversely, if all planar
maximizers throughout $I$ are simplex partitions, then this finite
comparison succeeds at some finite orders and meshes.

To construct the constraints, take a grid of $q$ equally spaced
rotations of the labeled simplex and include every label permutation.
For each grid representative $T'$ keep the overlap
$O_n(A,T')=\sum_i\langle P_n1_{A_i},P_n1_{T'_i}\rangle$.
These overlaps are linear in the integrated angular moment field since
the simplex is conical. For one simplex cell the squared Fourier tail
past $n\ge1$ is at most $2\sum_{j>n}(1/(j\pi))^2\le2/(\pi^2 n)$.
The sum of the three squared tail norms for $A$ is at most $2/3$,
so Cauchy–Schwarz bounds the overlap error by $2/(\pi\sqrt n)$.
The distance between two labeled simplex partitions rotated through
an angle $\theta$ is at most $3|\theta|/\pi$. Hence
\[
 \left|d(A)-2\big(1-\max_{T'\text{ in grid}}O_n(A,T')\big)\right|
       \le \epsilon_{n,q}:=4/(\pi\sqrt n)+3/q.                   \tag{SM7}
\]
If $d(A)\ge\delta$, every $O_n(A,T')$ is at most
$1-(\delta-\epsilon_{n,q})/2$. Conversely, if these linear constraints
hold then $d(A)\ge\delta-2\epsilon_{n,q}$. Averaging the retained
moment field over radial cells preserves the linear constraints.
All feasible step fields have ordinary realizers by the finite
moment result, so the reverse distance implication also holds
for the finite optimization.

For each correlation, upper bound the stability of the distant
partitions by the constrained finite optimum in (SM5) with the
step kernels, plus the errors (SM4) and (SM6). To make this comparison
on a correlation grid, the stability of every balanced partition
and of the simplex varies at rate at most $\tfrac23(1-b)^{-2}$
on $I$, by the Hermite expansion. The upper bound at a grid point
therefore controls the deficit at the neighboring correlations,
with a further error at most $\tfrac43(1-b)^{-2}$ times the distance
to that grid point. Strict negativity after adding these errors
proves global uniqueness using the local comparison inside the
distance threshold. By planar reduction for maximizers, the same
conclusion would hold in every $d\ge2$.

For the converse, if every global maximizer throughout $I$ is a
simplex, there is a uniform positive stability gap over all balanced
partitions with $d(A)\ge\delta/2$. Otherwise take a sequence with
stability tending to the simplex value and pass to convergent
correlations and weakly convergent indicators in $L^2(\gamma_2)^3$.
The set of balanced fractional assignments is weakly compact and the
positive-correlation Ornstein–Uhlenbeck operator is compact, so the
fractional maximum is attained. Strict convexity of the stability
functional makes every such optimizer ordinary: a nonordinary
assignment can be split into two distinct assignments of the same
balanced masses by perturbing two fractional labels oppositely on
subsets of a set where both have a positive margin. It follows that
the fractional and ordinary maxima agree. Compactness of the operator
and operator-norm continuity in $\rho$ now give a maximizing limit
of the sequence. This limit is ordinary and is a simplex by hypothesis.
The indicators then converge strongly by preservation of their
squared norms at the ordinary limit. This contradicts $d(A)\ge\delta/2$.

Choose $n,q$ with $2\epsilon_{n,q}<\delta/2$ and with the error in
(SM4) smaller than the uniform gap. Next choose sufficiently fine
radial and correlation meshes. The finite step objective at a feasible
point differs from the stability of its distant ordinary realizer
by at most the displayed errors. Thus a strict margin is present.
It can be verified by finite outer constraints on each $C_n$: their
support values are integrals of maxima of finitely many smooth
trigonometric polynomials and can be enclosed to arbitrary accuracy.
Refining the support directions gives nested compact outer polytopes
converging to $C_n$. Continuity of the quadratic objective on the
bounded coefficient boxes, and compactness, allow rigorous upper
bounds by refinement to preserve any strict margin. The radial
kernel integrals can likewise be enclosed by controlled quadrature
on boxes and the displayed exponential tails. This establishes
existence of a successful finite comparison under global uniqueness.
No such comparison is evaluated here on the unresolved interval.

## A linear mesh error for the mass profiles of maximizers

For an ordinary planar partition with finite Gaussian perimeter,
write the mass profile as a function of the radius $r>0$ and use
$d\mu(r)=r e^{-r^2/2}dr$. With the outward normal $N_i$ of each cell,
the radial distributional derivative of its profile is
\[
 Dp_i=-|x|_\#\left(\frac{N_i(x)\cdot (x/|x|)}{2\pi |x|}
              \mathcal H^1\!\restriction_{\partial^* A_i}\right).
                                                               \tag{SM8}
\]
Indeed test against a smooth compactly supported scalar function
$\psi(r)$ away from zero. The radial vector field
$X(x)=\psi(|x|)(x/|x|)/|x|$ has divergence $\psi'(r)/r$.
Integrating over the cell in polar coordinates and applying
Gauss–Green gives (SM8). In particular the profile is locally of
bounded variation on the positive radius axis and
\[
 \sum_i\int_0^\infty r e^{-r^2/2}d|Dp_i|(r)
              \le 2\sum_e\int |N_e\cdot(x/|x|)|d\nu_e(x)
              \le 2P_\gamma(A).                               \tag{SM9}
\]
Here $\nu_e$ is Gaussian surface measure, and the perimeter counts
each interface once.

Divide $[0,R]$ into intervals of length at most $h$ and include a
final tail interval $[R,\infty)$. On each interval replace $p$ by
its $\mu$-average, giving $p^{\rm step}$ with the same global masses.
Then
\[
        \int |p-p^{\rm step}|_1d\mu
                      \le h P_\gamma(A)+2e^{-R^2/2}.            \tag{SM10}
\]
To prove the coefficient of $h$, let $m(r)=r e^{-r^2/2}$ and consider
one finite interval $I=(a,b)$. Put $L(t)=\mu(a,t)$ and $H(t)=\mu(t,b)$.
By averaging differences and applying the one-dimensional variation
bound between the two arguments,
\[
 \int_I|f-f_I|d\mu\le\int_I\frac{2L(t)H(t)}{L(t)+H(t)}d|Df|(t).
\]
The integrals allow an interval starting at zero since the pointwise
variation estimate only uses pairs of positive radii. For every $t>0$
in $I$, the log-concavity of $m$ gives $m(r)\le m(t)e^{k(r-t)}$,
where $k=(\log m)'(t)$. Substitute the resulting upper bounds on
$L,H$ in their increasing harmonic-mean expression. When $k\ne0$,
as $t-a$ and $b-t$ range over nonnegative numbers with their sum
fixed, the expression is bounded by
$m(t)(2/|k|)\tanh(|k|(b-a)/4)\le m(t)(b-a)/2$.
For example for $k>0$ this follows on writing the numerator after
integration as proportional to
$1+e^{k(b-a)}-e^{k(t-a)}-e^{k(b-t)}$, which is maximized when the
two lengths are equal. Reflect the interval for $k<0$. For $k=0$
the same bound is the arithmetic-geometric mean inequality.
Sum the estimate over intervals and labels and use (SM9).
On the tail, the distance of two mass vectors is at most two.

In particular, any global planar maximizer has in (SM10) the
upper bound $P_\gamma(A)<2\phi(0)/\sqrt{1-\rho^2}$ from (IB2)
and (IB5) of the [interface bounds](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md).
The integrated exact circular pair envelope changes by at most
the right side of (SM10) under this replacement. In fact, for
the probability-normalized circular kernel, changing one mass
vector can be implemented by reassigning an angular mass equal
to half its $\ell^1$ change. The joint angular marginal is uniform,
so the pair value changes by at most that angular mass. Combine
both mass changes and use the two exact radial marginals $\mu$.
This gives a uniform explicit radial mesh error for the exact
pairwise envelope of actual maximizers. The simultaneous envelope
can instead use (SM3) or the full finite Fourier comparison.
