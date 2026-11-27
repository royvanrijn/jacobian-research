# The three-cell deficit and its perimeter-angle normalization

The universal nondecreasing law for
$(S_\rho(T)-S_\rho(A))/\arccos\rho$ is false, already for a balanced
three-slab partition in $\mathbb R^2$. The obstruction is an explicit
quadratic angle correction for the simplex junction. This note derives
that coefficient and the counterexample analytically. It does not
produce a partition improving the simplex. The
[all-correlation three-cell interface gap](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem)
is a separate local result; global tied maximizers and exact optimality
near correlation one remain open. These arguments have internal checks
only; external review and formal verification have not occurred.

## Exact correlation derivative for the three-cell simplex

Let $T$ be the central partition into three $120$-degree sectors in
$\mathbb R^2$, and set $\theta=\arccos\rho\in(0,\pi/2)$. Write
$N_A(\theta)=1-S_{\cos\theta}(A)$ and let
$\mathcal P_\gamma$ denote the ordinary Gaussian interface sum, each
interface counted once. For the simplex,
$\mathcal P_\gamma(T)=3\phi_1(0)/2$.

The simplex value has the known closed form

\[
 S_\rho(T)=\frac13+
       \frac3{4\pi^2}\big[\arccos(-\rho)^2-\arccos(\rho/2)^2\big].
 \tag{D0}
\]

It is printed in Heilman,
[*Three Candidate Plurality is Stablest for Correlations at most 1/10*,
arXiv:2306.03312v3](https://arxiv.org/pdf/2306.03312v3), Section 1.3,
the display after Conjecture 1.7 (printed p. 6), and Section 10 (p. 42),
where it is credited to KPW04. Those passages were inspected; the KPW04
original was not. Formula (D4) below is a direct expansion of this known
value, not a claim of a new simplex-value formula. The boundary calculation
also derives its derivative without importing (D0).

Use the ray kernels $M_t$ and measure $\nu$ of (S14) in the linked
note. A Gaussian quadrant integral gives

\[
 I_t(\rho):=\iint M_t(r,s)d\nu(r)d\nu(s)
   =\frac{1}{2\pi\sqrt{1-\rho^2t^2}}
        \left[\frac14+\frac{\arcsin(\rho t)}{2\pi}\right].
 \tag{D1}
\]

This follows by completing the square to standardize a centered
two-dimensional Gaussian of correlation $\rho t$; its quadrant
probability is $1/4+\arcsin(\rho t)/(2\pi)$, also obtained by polar
integration. The two outward facet normals of any one cell have inner
product $1/2$, while the ray directions have inner product $-1/2$.
Consequently

\[
 \frac{d}{d\rho}S_\rho(T)=6I_1(\rho)+3I_{-1/2}(\rho).
 \tag{D2}
\]

To verify the derivative without a stationarity or maximality
assumption, the joint Gaussian density satisfies
$\partial_\rho p_\rho(x,y)=\sum_a\partial_{x_a}\partial_{y_a}p_\rho(x,y)$.
Integrate once in each cell variable. The result is the boundary-pair
integral with coefficient $n_i(x)\cdot n_i(y)$. There are six diagonal
ray terms and six ordered off-diagonal terms of normal inner product
$1/2$, giving (D2). Truncation removes the boundary at infinity by
Gaussian decay at every fixed $\rho<1$. The junction is surface-null.

Since $\arcsin(\cos\theta)=\pi/2-\theta$, (D1)--(D2) yield the exact
angle derivative

\[
\begin{split}
 N_T'(\theta)={}&\frac3{2\pi}-\frac{3\theta}{2\pi^2}\\
 &+\frac{3\sin\theta}{2\pi\sqrt{1-\cos^2\theta/4}}
   \left[\frac14-\frac{\arcsin(\cos\theta/2)}{2\pi}\right].
\end{split}
 \tag{D3}
\]

The bracket and square-root factor in the second line are analytic
even functions near zero. The whole second line is odd, with linear
coefficient $1/(2\pi\sqrt3)$. Thus

\[
 N_T(\theta)=\frac3{2\pi}\theta-\beta_3\theta^2+O(\theta^4),
 \qquad
 \boxed{\displaystyle\beta_3=\frac3{4\pi^2}
                           -\frac1{4\pi\sqrt3}>0.}
 \tag{D4}
\]

The leading coefficient is
$\sqrt{2/\pi}\,\mathcal P_\gamma(T)$. The second coefficient is
nonzero; it cannot be discarded in a comparison normalized by the
perimeter angle. Positivity follows from $\pi<3\sqrt3$.

## A balanced three-slab counterexample

Put $q=\Phi^{-1}(2/3)>0$ and let

\[
 A_1=\{x_1<-q\},\qquad A_2=\{-q\le x_1\le q\},\qquad
 A_3=\{x_1>q\}.
 \tag{D5}
\]

These cells have Gaussian measure $1/3$, with two smooth disjoint
interfaces. They have ordinary interface perimeter
$\mathcal P_\gamma(A)=2\phi_1(q)$.

For a single threshold $q$, let $C_q(\theta)$ be the probability that
the two correlated scalar observations lie on opposite sides. Integrating
$\partial_\rho p_\rho=\partial_x\partial_y p_\rho$ on the lower
quadrant gives the elementary formula

\[
 C_q(\theta)=\frac1\pi\int_0^\theta
          \exp\left(-\frac{q^2}{1+\cos u}\right)du.
 \tag{D6}
\]

The number of thresholds between the two observations is one when
adjacent slab labels differ and two when the extreme labels differ.
Writing $J(\rho)=\Pr(X<-q,Y>q)$ gives the exact identity

\[
 N_A(\theta)=2C_q(\theta)-2J(\cos\theta).
 \tag{D7}
\]

The extreme-label event implies $Y-X>2q$, so a Gaussian Chernoff bound
gives $J(\rho)\le\exp(-q^2/(1-\rho))$. This is smaller than every
power of $\theta$ as $\theta\downarrow0$. The integrand in (D6) is
analytic and even near zero. Hence

\[
 N_A(\theta)=\frac2\pi e^{-q^2/2}\theta+O(\theta^3).
 \tag{D8}
\]

Unlike (D4), there is no quadratic term. Subtracting (D4) from (D8)
therefore gives

\[
 \frac{S_\rho(T)-S_\rho(A)}{\arccos\rho}
   =c+\beta_3\arccos\rho+O((\arccos\rho)^2),\qquad
 c=\sqrt{2/\pi}\,[\mathcal P_\gamma(A)-\mathcal P_\gamma(T)].
 \tag{D9}
\]

For all sufficiently large $\rho<1$ the expression is greater than
its limit $c$ at one. It therefore cannot be nondecreasing in $\rho$.
The same balanced example works in every ambient dimension $d\ge2$
by cylindrical extension.

The limit is positive here: $q<1/2$, since
$\int_0^{1/2}\phi_1(t)dt>1/6$ (use $\phi_1(1/2)>1/3$).
Then $e^{-q^2/2}>7/8$, so $c>1/(4\pi)>0$. As the normalized deficit
tends to zero at $\rho=0$, it is not nonincreasing either. This is a
counterexample to both directions of universal monotonicity for this
particular normalization. It does not exclude a different monotone
quantity or a comparison restricted to actual maximizers.

## Consequence for the near-one programme

The correction in (D4) is of order $\theta^2$, or $1-\rho$, whereas
smooth separated interfaces have no term at that order in this example.
Thus even the fixed simplex prevents a general replacement of its
heat content by just its perimeter with an $o(\theta^2)$ error. A
deformation-relative expansion for moving competitors must account for
the junction contribution. The calculation supplies that obstruction,
not an error estimate uniform over competitors or exact near-one
optimality.

The formulas are analytic. The companion
`check_simplex_interface_symmetry.py` records finite normalization and
coefficient controls; its floating quadrature does not prove the
asymptotic statement or the unrestricted counterexample. No stationary
partition classification is assumed. In particular the existing
[hyperstability source audit](MAXIMIZER_PERTURBATIVE_RIGIDITY.md)
retains the extra mixed-variation condition in that theorem; ordinary
second-variation stability alone is still the open classification target.
