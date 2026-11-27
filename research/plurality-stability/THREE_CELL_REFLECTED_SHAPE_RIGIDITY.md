# Shape restrictions for a maximizer with one central reflection

Use the notation of the [reflection reduction](THREE_CELL_REFLECTION_REDUCTION.md).
In particular a balanced global maximizer with a single
central reflection has the fixed cell above the graph
$y=h(|x|)$, with a unique junction $(0,h_0)$ and
$h_0=h(0)$. Its two other cells occupy the right and
left parts below the graph. Let $m_*=1/\sqrt3$ be the
slope of the balanced simplex graph.

## The graph crosses the simplex ray through its junction

If the reflected maximizing partition is nonsimplex, its
height is strictly above $h_0+m_*x$ somewhere and strictly
below it somewhere on $x>0$.

Write $m=h'=f/w$ and introduce the probability kernel
\[
 P_p(dq)=\frac{(2M+M_R)w(q)d\nu(q)}{D(p)},\qquad
 D(p)=\int(2M+M_R)w\,d\nu=a(p)w(p)/\rho
\]
on the right graph. It is reversible with respect to the
finite measure $a w^2d\nu$, because the kernel times this
measure is symmetric in the two arguments. The two translation
equations (RF3) give
\[
             m-Pm=\frac{J_0-2I_R}{D},\qquad
             I_R=\int M_R f\,d\nu.                          \tag{RS1}
\]
Thus the right side has zero integral against $a w^2d\nu$.

Suppose $h(s)\le h_0+m_*s$ for all $s\ge0$. Balance forces
$h_0\ge0$, and $h_0=0$ can occur in this comparison only at
the simplex: otherwise the fixed cell would strictly
contain the reference simplex cell in Gaussian measure.
Write $k=h^{-1}$ and set, at $p=(x,h(x))$,
$z=\rho x/\sigma$, $W=(\rho h(x)-h_0)/\sigma$, and
$\bar\Phi=1-\Phi$. Parametrizing the reflected source
interface by its height $h_0+\sigma v$ gives
\[
 \frac{I_R}{J_0}=\int_0^\infty
    \frac{\phi(z+k(h_0+\sigma v)/\sigma)}{\phi(z)}
    \frac{\phi(v-W)}{\bar\Phi(W)}\,dv .                     \tag{RS2}
\]
Here the denominator $J_0$ integrates the conditional
Gaussian on the downward axis ray. The assumed comparison
gives $k(h_0+\sigma v)/\sigma\ge v/m_*$ and $W\le m_*z$,
with the latter strict for a nonsimplex. For $v\ge0$ the
second quotient in (RS2) is increasing in $W$, since it
equals the increasing Gaussian hazard
$\phi(W)/\bar\Phi(W)$ times $e^{vW-v^2/2}$. The first
quotient decreases in the nonnegative inverse height.
Thus the ratio is strictly below its value for the simplex,
which equals $1/2$ by (RS1) with constant slope. This
contradicts the invariant mean of (RS1). The same comparison
reverses if $h(s)\ge h_0+m_*s$: balance requires $h_0<0$
unless the graph is the simplex, the inverse height is
smaller and $W>m_*z$. This proves both crossings.

## Global convexity or concavity forces the simplex

Suppose first that $h$ is globally convex. Its bounded
derivative increases to a limit $m_\infty>0$. The
discrete homogeneity and bounded error in (RF5) then give
the linear asymptotic profile $m_\infty x$. In particular
\[
       h(x)=m_\infty x+b+v(x),\qquad v(x)\downarrow0,
\]
since $h(x)-m_\infty x$ is decreasing and bounded.
For a concave $h$ the same argument gives $v(x)\uparrow0$.

Let $\beta=U_3-U_1$ on the right graph. On the tail,
conditioning on the horizontal source coordinate and
using the affine asymptote shows
$\beta=\lambda(m_\infty,b)$, where
\[
 \lambda(m,c)=2\Phi\!\left(\frac{-(1-\rho)c}
                                {\sigma\sqrt{1+m^2}}\right)-1 .
\]
This is the full binary affine score at observation
$(x,mx+c)$. To compute the correction for the folded
affine problem, use $z=\rho x/\sigma$,
$d=(1-\rho)c/\sigma$ and define
\[
 A(z)=\int_0^\infty\phi(z+t)\bar\Phi(mz-d+mt)\,dt,\qquad
 B(z)=\int_0^\infty\phi(z+t)\bar\Phi(mz-d-mt)\,dt .
\]
Only the negative horizontal source contributes a
correction. There the fixed cell uses $m|x|+c$, and the
right-hand label is absent. Subtracting the full binary
score gives $2A(z)-B(z)$.

This correction is strictly negative for large $z$ when
$m>m_*$, and strictly positive for large $z$ when
$0<m<m_*$. The signs are uniform for $c$ in a compact
interval at fixed $m,\rho$. To see this, scale $t=u/z$
and divide by $\phi(z)\bar\Phi(mz-d)/z$. The respective
limiting integrands are $e^{-(1+m^2)u}$ and
$e^{-(1-m^2)u}$, since the Gaussian hazard at a large
positive $s$ is $s+O(1/s)$. The first integral converges
by domination by $e^{-u}$. For $m<1$ the second converges
as well: monotonicity of the hazard bounds its tail
ratio by $\exp[(mu/z)\mathcal H(mz-d)]$, and the product
is dominated by $e^{-(1-m^2)u/2}$ for large $z$.
For $m\ge1$ its limit inferior is infinite by Fatou.
For $m<1$ the sign follows from
$2/(1+m^2)-1/(1-m^2)$.

For convex $h$, compare the true score at $(x,h(x))$
to the folded affine graph of slope $m=m_\infty$ and
intercept $b+v(x)$. For $0\le s\le x$ the true
source threshold lies above this affine graph, so its
contribution to $U_3-U_1$ is no larger. For $s>x$ the
possible upward score error is at most $C v(x)$
times the conditional tail probability of $|Y_x|>x$,
because $v(s)\ge0$ and one-dimensional Gaussian
probabilities are Lipschitz in the vertical threshold.
This gives
\[
 \beta\le\lambda(m,b+v(x))+2A(z)-B(z)
                                 +C v(x)e^{-c_1 x^2}.      \tag{RS3}
\]
If $m>m_*$, the correction is negative, and
$\lambda(m,b+v(x))\le\beta-c_2v(x)$ with $c_2>0$.
This contradicts (RS3), so $m\le m_*$. Convexity
would then give $h(x)\le h_0+m_*x$ everywhere,
excluded for a nonsimplex by (RS1)–(RS2). For
concave $h$ the threshold comparison reverses, with
error bounded by $C|v(x)|e^{-c_1x^2}$. A slope
$m<m_*$ would have positive correction and
$\lambda(m,b+v(x))\ge\beta+c_2|v(x)|$, giving the
opposite contradiction. Thus $m\ge m_*$, and
$h(x)\ge h_0+m_*x$ is likewise excluded. Hence a
nonsimplex reflected maximizing graph has curvature
of both signs.

## Local affine limits at infinity

Every translated tail limit of a reflected maximizing graph
is affine. In particular $h''(x)\to0$ and
$h'(x)-h'(\rho x)\to0$. The affine slopes of local limits
can depend on the sequence used to take the limit.

By (RF4)–(RF5) and the slope bounds,
$|h(\rho x)-\rho h(x)|\le C$ on the tail. In (RF3)
restrict the vertical-translation integral to
$|s-\rho x|\le1$. Both conditional Gaussian factors are
then bounded below, giving $a w\ge c>0$. The corresponding
lower bound on the vertical score derivative, together
with uniform smooth Gaussian score bounds, gives bounded
derivatives of $h$ of every fixed order by the implicit
function theorem.

Take any $x_j\to\infty$. Along one subsequence, for every
integer $n\ge0$ the profiles and constants
\[
 H_n(t)=\lim_j[h(\rho^n x_j+t)-h(\rho^n x_j)],\qquad
 d_n=\lim_j[\rho h(\rho^n x_j)-h(\rho^{n+1}x_j)]
\]
exist in local smooth convergence for the profiles.
Use a diagonal subsequence with the uniform derivative
bounds above. Set $m_n=H_n'$. In the two translation
equations the reflected terms and the axis term vanish
in these local limits. Their conditional Gaussian weights
are exponentially small when the observation has large
positive horizontal coordinate; the bounded slope controls
integration in source coordinates. Dividing the limiting
horizontal identity by the vertical identity gives
\[
 m_n(t)=\int m_{n+1}(s)\,P_{n,t}(ds),\qquad
 P_{n,t}(ds)\ \propto\
 p_\sigma(s-\rho t)
 p_\sigma(H_{n+1}(s)-\rho H_n(t)-d_n)\,ds .                 \tag{RS4}
\]
Passing to the limit in the unreflected factors is justified
by bounded slopes and Gaussian domination in the horizontal
argument. The bound $|h(\rho x)-\rho h(x)|\le C$ implies
$|H_{n+1}(\rho t)-\rho H_n(t)-d_n|\le C$ for all $n,t$.
Together with the uniform slope bound this gives a
normalizing integral in (RS4) bounded below uniformly:
restrict to $|s-\rho t|\le1$. It is bounded above by
$p_\sigma(0)$. Thus $P_{n,t}\le C_1p_\sigma(s-\rho t)\,ds$
and $\int|s|P_{n,t}(ds)\le\rho|t|+C_2$.
For any fixed $R$, on $|t|\le R$ the kernel also
dominates a fixed positive multiple of the uniform
probability on $[-1,1]$, by the same residual and slope
bounds.

Run two chains with these transitions starting at
$t$ and $t'$. At each simultaneous visit to
$|t|+|t'|\le R$ use the common uniform component
to couple them, then keep them equal. The sum of their
absolute values has expected next value at most its
current value times $\rho$ plus $2C_2$, for every coupling.
For $R$ large enough this is a drift of at least one
toward the compact set outside it. Stopping on entry
and taking expectations shows that the expected entry
time from an outside starting state is at most the
starting sum of absolute values. This holds at every
starting time. Thus after each failed coupling attempt
the chains return to the compact set almost surely,
and eventually couple almost surely by the uniform
success probability. Iterate (RS4). Boundedness of
all $m_n$ shows that the difference of the starting
values of $m_0$ is bounded by twice that bound times
the probability of no coupling by step $N$. It tends
to zero. Hence $H_0$ is affine, and the same proof
starting at any $n$ shows every $H_n$ is affine.
Equation (RS4) makes their slopes equal. The stated
tail derivative limits follow by subsequence compactness.

## A constant correlation derivative constrains the tail

Write $g(x)=\partial_\rho(U_3-U_1)(x,h(x))$, taking
the derivative at fixed shape and fixed observation.
If $g$ is bounded, the asymptotic profile in (RF5) is
linear and $h(x)=m x+O(1)$, $h'(x)=m+O(1/x)$ for
some $m>0$. If $g$ is constant on the graph, its value
is forced by the price $\beta$:
\[
                  g=-\,\frac{2 z\phi(z)}{1-\rho^2},
                  \qquad z=\Phi^{-1}((1+\beta)/2).         \tag{RS5}
\]

For the proof set $S=\rho x+\sigma Z$ for $Z\sim N(0,1)$
and $L=(\rho h(x)-h(S))/\sigma$, using an arbitrary
smooth bounded-slope extension for negative $S$.
The negative-source branch and its correlation
derivative are exponentially small for large $x$.
The binary part of the score is $2\mathbb E\Phi(L)-1$.
Write $A_x=\mathbb E\phi(L)\ge c>0$, using (RF4) and the
slope bounds on $|S-\rho x|\le1$. Equation (RF3), with
$f=h'w$, gives
$\mathbb E[\phi(L)h'(S)]=A_x h'(x)+O(e^{-c_1x^2})$.
Indeed the difference is given by the reflected and axis
terms, which have that Gaussian bound. Differentiating
the conditional expression for the full score therefore
gives
\[
 g(x)=\frac{2A_x}{\sigma}(h(x)-x h'(x))
       +\frac{2\rho}{\sigma^2}
                \mathbb E[\phi(L)(Z h'(S)+L)]
       +O(e^{-c_2x^2}).                                    \tag{RS6}
\]
All terms in the second expectation are bounded in
absolute mean, by bounded slope and boundedness of
$|L|\phi(L)$. Thus bounded $g$ gives
$|(h(x)/x)'|\le C/x^2$, proving the assertions about
$m$ and the linear asymptote.

Let $b(x)=h(x)-m x$. It is bounded with $b'(x)=O(1/x)$,
so in the conditional score replacing $b(S)$ by $b(\rho x)$
has error $O(1/x)$. Gaussian integration of the remaining
affine expression gives
$\rho b(x)-b(\rho x)=\delta+O(1/x)$ with
$\delta=\sigma\sqrt{1+m^2}\,z$. Iterating outward using
boundedness of $b$ shows
$b(x)=b_*+O(1/x)$ with $b_*=\delta/(\rho-1)$.
Consequently $L$ tends in the conditional expression to
$\delta/\sigma-mZ$, $A_x\to A_*=\phi(z)/\sqrt{1+m^2}$,
and the second expectation in (RS6) tends to
$\delta A_*/\sigma$. If $g$ is constant, (RS6) now
forces $x b'(x)$ to have a finite limit. Boundedness of
$b$ requires that limit to be zero. Substituting $b_*$,
$\delta$ and $A_*$ in (RS6) gives (RS5).

## A pressure-free higher-order reflected tie is impossible

Suppose a reflected global maximizer ties the simplex and
its fixed-partition deficit has $D'=D''=0$ at a correlation
approached by distinct simplex-optimal correlations.
The [mixed tie inequality (G10)](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-first-crossing-obstruction-at-the-level-of-hermite-weights)
makes $g$ constant along the graph: one can test with any
compact regular speed of zero mean on this one interface.
Thus (RS5) applies. Furthermore, for $q(x)=x h'(x)-h(x)$
the scalar dilation identity (G6) and the first equation
of (RF3) give
\[
        q=\rho Pq-\frac{\rho(1-\rho^2)}{a w}\,g,            \tag{RS7}
\]
with $P$ from (RS1). The axis has zero dilation speed,
so there is no axis term in this even block. Boundedness
of $q$ follows from (RS6); $a w$ is bounded below on
the whole graph. Invert $I-\rho P$ on bounded functions
by its positive series. If $\beta=0$, then $g=0$ by (RS5),
so $q=0$ and balance makes $h(x)=x/\sqrt3$. If
$\beta\ne0$, then $q$ has the strict sign of $\beta$
by (RS5)–(RS7). In that case $h(0)$ has the opposite sign
and $h(x)/x$ is strictly monotone. If $\beta>0$ the
graph is below its limiting unshifted ray $mx$, so
balance requires $m>1/\sqrt3$. If $\beta<0$, it is
above the ray and $m<1/\sqrt3$. The pressure-free
conclusion does not cover these nonzero-price cases.
