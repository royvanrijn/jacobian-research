# A shared first angular moment comparison for two equal circular masses

Suppose two labels have equal masses on almost every circle. The maximum
first angular interaction over all partitions with these masses is attained
by a single family of arcs, with two arcs and their complement at each radius.
One scalar function parametrizes the optimization, and its square has a
concave objective. The result applies to every nonnegative integrable radial
kernel, including the interaction subtracted in the global three-cell
comparison. It is an internal analytic result.

## Exact simultaneous maximization

Write $t=r^2/2$, $d\nu(t)=e^{-t}dt$, and $d\sigma(\theta)=d\theta/(2\pi)$.
Prescribe circular masses $p(t)=(x(t),(1-x(t))/2,(1-x(t))/2)$ for a
measurable $x$ valued in $[0,1]$. Set
\[
 \beta(t)=\frac{\cos(\pi x(t)/2)}{\pi},\qquad
 b(t)=\sin^2(\pi x(t)/2),\qquad
 K(t,s)=2W(t,s)\beta(t)\beta(s),
\]
where $W$ is a nonnegative symmetric kernel in $L^1(\nu\otimes\nu)$.
For slices with these masses let $u_i(t)=\int_{A_i(t)}(\cos\theta,\sin\theta)
d\sigma(\theta)$. Then
\[
 \max_A\iint W(t,s)\sum_{i=1}^3 u_i(t)\cdot u_i(s)d\nu(t)d\nu(s)
 =\max_{0\le y\le b}\iint K(t,s)
    \big[3\sqrt{y(t)y(s)}+\sqrt{(1-y(t))(1-y(s))}\big]d\nu(t)d\nu(s).
                                                               \tag{EP1}
\]
One may fix $y=0$ at $x=1$, where the kernel has zero weight.

To prove this, use the orthonormal label coordinates $(2,-1,-1)/\sqrt6$
and $(0,1,-1)/\sqrt2$ of $H=\mathbf1^\perp$. Let $G(t)$ be the Gram on
$H$ of the closing moment triple at $t$. For any pair, its dot product
is at most the matrix fidelity
$F(G(t),G(s))=\operatorname{tr}[(G(t)^{1/2}G(s)G(t)^{1/2})^{1/2}]$,
as in [AM2](THREE_CELL_EXACT_ANGULAR_MOMENT_ENVELOPE.md). The infimum
formula there makes $F$ jointly concave, including singular Grams.
Conjugating both Grams by the label swap $2\leftrightarrow3$ leaves $F$
unchanged. Averaging each Gram with its conjugate therefore increases
the upper bound for every pair simultaneously. The averaging preserves
the squared norm caps because the two swapped caps are equal.

The averaged Gram is diagonal, with entries $e,f\ge0$. The caps on the
three physical vectors are equivalent to
\[
 e\le\tfrac32(\sin(\pi x)/\pi)^2=6\beta^2 b,
 \qquad e/6+f/2\le\beta^2.
\]
At a pair of diagonal Grams, $F=\sqrt{e(t)e(s)}+\sqrt{f(t)f(s)}$.
We can increase $f$ to $2\beta^2-e/3$ and write $e=6\beta^2y$.
This gives the right side of (EP1) as an upper bound.

For the converse, use the same axes at every radius. At $x<1$, put
$\vartheta=\arccos\sqrt y$ and take the arcs of labels $2,3$, each of
angular width $\pi(1-x)$, centered at $\pi+\vartheta$ and
$\pi-\vartheta$, respectively. They are disjoint since
$\pi(1-x)/2\le\vartheta\le\pi/2$. Give their complement to label $1$.
The first moment of label $1$ is $(2\beta\sqrt y,0)$, and the two other
moments are $(-\beta\sqrt y,\mp\beta\sqrt{1-y})$. Their dot product
between radii is exactly the integrand of (EP1) divided by $W$.
The construction is measurable and also works at $x=0$. At $x=1$, all
moments vanish.

For existence, the objective on the right is concave in $y$ by the
concavity of the geometric mean on nonnegative arguments. It is
continuous on the feasible set in the strong $L^2(\nu)$ topology:
from every convergent sequence extract a subsequence converging pointwise
and use the integrable bound $4K$. It is consequently weakly upper
semicontinuous on that convex set, which is weakly compact by its
$L^2$ bound and its closed pointwise constraints. Thus the maximum is
attained, without having to make choices of different arc orientations
for different partner radii.

## When the pair bounds can be attained simultaneously

For two circular masses parameterized by $x,z\in[0,1)$, put
$a=\min(\pi x/2,\pi z/2)$ and $c=\max(\pi x/2,\pi z/2)$. The exact
pair maximum $\mathfrak F$ from (AM2) is
\[
 2\beta_x\beta_z
 \begin{cases}
  3\sin a\sin c+\cos a\cos c,&\tan c\le3\tan a,\\
  \sqrt{1+8\sin^2 a},&\tan c>3\tan a.
 \end{cases}                                                   \tag{EP2}
\]
If either parameter is $1$, the value is zero. To see (EP2), the
same Gram averaging and arc construction reduce the problem to maximizing
$3\sin u\sin v+\cos u\cos v$ on $0\le u\le a$, $0\le v\le c$.
At a given $u$ the unique maximizing $v$ is
$\min(c,\arctan(3\tan u))$, which is at least $u$. The derivative in
$u$ at this choice is $3\cos u\sin v-\sin u\cos v$, strictly positive
for $u>0$. Thus the unique choice is $u=a$ and the stated $v$,
also when $a=0$. In the first case of (EP2), the exact angular loss
relative to the sum of the three individual norm products is
\[
 4\beta_x\beta_z\sin^2\big(\pi(x-z)/4\big).                     \tag{EP3}
\]

Assume $W>0$ almost everywhere. On the set $D=\{t:x(t)<1\}$ of
positive weight, write $\tau(t)=\tan(\pi x(t)/2)$. If $D$ is null,
all interactions are zero. Otherwise the shared maximum (EP1) equals
the integral of the separate pair maxima (EP2) precisely when
\[
        \operatorname*{ess\,sup}_{D}\tau
           \ \le\ 3\operatorname*{ess\,inf}_{D}\tau .           \tag{EP4}
\]
Under (EP4), $y=b$ attains all the pair maxima. For necessity, take a
maximizer of (EP1). If its value equals the integrated pair maximum,
almost every pair must attain the pointwise optimum. For almost every
$t$ with $0<\tau(t)<\infty$, there is a positive-measure set of
partners $s$ with $\tau(t)/3<\tau(s)<3\tau(t)$. This follows by taking
$\tau(t)$ in the support of the pushforward of $\nu|_D$, and using
Fubini for the pair equality. For one of those pairs, (EP2) forces
$y(t)=b(t)$. The same equality is forced by feasibility when $\tau(t)=0$.
With $y=b$ almost everywhere on $D$, the unique pair optimizer forces
$\tau(t)\le3\tau(s)$ for almost every $(t,s)\in D^2$, yielding (EP4).
Thus whenever (EP4) fails the integrated separate pair bound is
strictly loose. This includes any profile with an origin interval on
which the singled label is absent and a positive set of three-label
circles.

There is also a weaker criterion for when the shared optimum itself
still uses $y=b$. Under the same positivity assumption, it does so if
and only if for almost every $t$ with $0<x(t)<1$,
\[
 \tan(\pi x(t)/2)\le 3\,
 \frac{\int W(t,s)\beta(s)\sin(\pi x(s)/2)d\nu(s)}
      {\int W(t,s)\beta(s)\cos(\pi x(s)/2)d\nu(s)} .            \tag{EP5}
\]
At $y=b$, the derivative coefficient of the concave objective is
$3\int K(t,s)\sqrt{b(s)}d\nu(s)/\sqrt{b(t)}-
\int K(t,s)\sqrt{1-b(s)}d\nu(s)/\sqrt{1-b(t)}$.
Nonnegativity is necessary by decreasing $y$ on measurable subsets
where $b$ is bounded away from $0,1$ and the kernel row integral is
bounded. For sufficiency use the tangent bound for each geometric
mean. The variations at $b=0$ are zero, and the $b=1$ rows have zero
weight. The positive part of the integrated tangent bound is finite:
$y-b\le0$, and $\beta(t)/\sqrt{1-b(t)}=1/\pi$ on $D$. Thus when the
derivative coefficient is nonnegative the integrated tangent bound is
nonpositive. Cancelling the common row factor gives (EP5).

## Consequence for the Gaussian comparison

Apply (EP1) to $W=V_\rho$ in (GAC2) of the
[global angular comparison](THREE_CELL_GLOBAL_ANGULAR_COMPATIBILITY.md).
If $a_i(t)=\sin(\pi p_i(t))/\pi$ and $M_{V_\rho}(x)$ is the right
side of (EP1), every measurable partition with the stated circular
masses satisfies
\[
 S_\rho(A)-S_\rho(T)\le
   U_\rho(p)-\iint V_\rho(t,s)\sum_i a_i(t)a_i(s)d\nu(t)d\nu(s)
                  +M_{V_\rho}(x).                             \tag{EP6}
\]
This comparison improves the bound obtained by integrating (AM2)
whenever (EP4) fails. For example, it applies to the mass profiles of
all genuinely reflected partitions, including the reflected maximizers
in [RF1](THREE_CELL_REFLECTION_REDUCTION.md). No reflection assumption
on the actual slices is needed. The scalar concave optimization computes
the full shared first-mode envelope for these profiles. Proving the
sign of (EP6) for all balanced profiles, adding other angular harmonics,
and treating profiles without a common equal pair remain open.
