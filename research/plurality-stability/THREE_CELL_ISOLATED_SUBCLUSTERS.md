# Restrictions on bounded components surrounded by one cell

Take an actual planar three-cell global maximizer at $0<\rho<1$.
Write $O_j$ for the open strict-winner set of label $j$ under
the adjusted stationary scores. A bounded component $S$ of
$\mathbb R^2\setminus\overline{O_j}$ cannot be isolated from
the remaining non-$j$ cells by an arbitrarily wide region of
label $j$. In a fixed central disk, such components also have
a positive minimum diameter under the quantitative condition
below. These restrictions allow either or both non-$j$ labels
within $S$.

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
 \left[\sigma^2|v|^2-igl((y-\rho x)\cdot v\bigr)^2
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
