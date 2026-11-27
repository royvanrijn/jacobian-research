# The sharp positive coefficient in the aligned-arc bound

The sum of the three positive Taylor remainders in the aligned-arc
comparison has an exact optimal quadratic coefficient. It is attained
at circular mass sums $(0,1,1)$, up to label permutation. This gives a
stronger scalar radial operator than the one in the
[unrestricted comparison through $53/100$](THREE_CELL_GLOBAL_SCALAR_THROUGH_HALF.md).
The directed estimates below extend the unrestricted comparison to
$29/50$.

## Sharp angular inequality

Use $k_a(v)=e^{a\cos(\pi v)}$ for $a\ge0$ and the even second primitive
$G_a$ with $G_a(0)=G_a'(0)=0$. Put $b=1/3$ and

\[
 P_\sharp(a)=3[2G_a(1)-3G_a(2/3)].                         \tag{SP1}
\]

For every $s\in[0,2]^3$ with $\sum_i s_i=2$,

\[
 \sum_i[G_a(s_i)-G_a(2/3)]
      \le \frac{P_\sharp(a)}2\sum_i(s_i-2/3)^2.           \tag{SP2}
\]

The coefficient is best possible for each $a$. At $s=(0,1,1)$,
the squared sum of displacements is $2/3$ and (SP2) is an equality.

To prove (SP2), use the normalized Taylor remainder

\[
 R_a(\delta)=2\int_0^1(1-u)k_a(2/3+u\delta)\,du.
\]

Write $P_0(a)=R_a(-2/3)$, $P_1(a)=R_a(1/3)$ and
$E_0(a)=e^{-a/2}$. Evaluating the Taylor remainders at $(0,1,1)$
gives $P_\sharp=(2P_0+P_1)/3$. Also $P_0\ge E_0$ and

\[
             P_0+2P_1\ge3E_0.                            \tag{SP3}
\]

For (SP3), the left side divided by three is an average of
$e^{ac}$, where $c=\cos(\pi v)$. Its mean $c$ is the derivative
at $a=0$. Direct integration gives

\[
 P_0'(0)=\frac{3\sqrt3}{2\pi}-\frac{27}{4\pi^2},\qquad
 P_\sharp'(0)=-\frac{3}{2\pi^2},\qquad
 P_0'(0)+2P_1'(0)=\frac{45}{4\pi^2}-\frac{9\sqrt3}{2\pi}>-\frac32.
\]

The last inequality follows from
$6\pi^2-18\sqrt3\pi+45=6(\pi-3\sqrt3/2)^2+9/2>0$.
Jensen's inequality proves (SP3), which implies
$P_\sharp\ge(P_0+E_0)/2\ge E_0$.

Suppose first there is one negative displacement $-d$ in
$s-(2/3)\mathbf1$. The other two, totaling $d$, lie in $[0,2/3]$.
On that interval the function

\[
 H(\delta)=G_a(2/3+\delta)-G_a(2/3)-G_a'(2/3)\delta
                         -\tfrac12P_\sharp(a)\delta^2
\]

is concave because $k_a(2/3+\delta)\le E_0\le P_\sharp$.
Equalizing the two nonnegative displacements can only increase
their summed $H$. For the resulting triple $(-d,d/2,d/2)$, the
normalized remainder coefficient is
$(2R_a(-d)+R_a(d/2))/3$. This is nondecreasing for $0\le d\le2/3$.
Indeed, for each $u\in[0,1]$, set $\theta=\pi(2/3-ud)$.
The two angles in the integrand are $\theta$ and $\pi-\theta/2$.
On $0\le\theta\le2\pi/3$ we have
$\sin\theta\ge\sin(\theta/2)$ and
$\cos\theta\ge\cos(\pi-\theta/2)$. The positive derivative from
the negative displacement dominates the negative derivative from
the positive displacement, including their factors $2$ and $1/2$.
The coefficient at $d=2/3$ is $P_\sharp$.

Suppose instead there are two negative displacements and a positive
one $D$. Each negative remainder coefficient is at most $P_0$ by
monotonicity. If $D\le2/3$, the positive coefficient is at most $E_0$
and the negative displacements carry at most half the squared sum.
If $D>2/3$, let $r=(2/3)/D\in[1/2,1]$. The negative squared share
is at most

\[
 \omega=\frac{r^2+(1-r)^2}{1+r^2+(1-r)^2},
\]

since their magnitudes each lie in $[0,2/3]$ and add to $D$.
For the positive remainder, reflect the part beyond $4/3$ about $1$.
Its density on $[v_0,2/3]$, $v_0=4/3-D$, is $2(v-v_0)/D^2$.
Relative to the density $(9/2)v$ of $P_0$, this is at most
$r-r^2$, by evaluating the increasing ratio at $v=2/3$.
Subtracting $E_0$ before bounding this part gives
$R_a(D)\le E_0+(r-r^2)(P_0-E_0)$. Therefore the full coefficient
of $P_0$ is at most $\omega+(1-\omega)(r-r^2)=1/2$.
Both cases are bounded by $(P_0+E_0)/2\le P_\sharp$. This proves
(SP2), since the linear Taylor terms sum to zero.

For $a>0$, equality in (SP2) requires $s=(2/3,2/3,2/3)$ or a
permutation of $(0,1,1)$. In the one-negative case the concavity
used to equalize the other coordinates is strict and the resulting
coefficient strictly increases with $d$. In the two-negative case
the bound $(P_0+E_0)/2<P_\sharp$ is strict by the strict mean
inequality in (SP3).

## Consequences for the radial operator

In (GH2)--(GH6) of the [unrestricted comparison](THREE_CELL_GLOBAL_SCALAR_THROUGH_HALF.md),
one can replace $P_*$ by $P_\sharp$. The negative coefficient
$M(a)=2\int_0^1(1-v)k_a(v)\,dv$ stays the same.
The probability measure defining $P_\sharp$ has density $3v$ on
$[0,2/3]$ and $6(1-v)$ on $[2/3,1]$. Its integral against
$\cos(l\pi v)$ for integer $l>0$ is

\[
 \frac{9\cos(2l\pi/3)-3-6(-1)^l}{(l\pi)^2}.              \tag{SP4}
\]

The integral is one for $l=0$. The binomial cosine expansion used
in (GH6) therefore computes the new moments by the same exact method.

Write $M_n,P_n$ for the moments of $\cos^n(\pi v)$ in the measures
of $M$ and $P_\sharp$. Folding at $1/2$ gives the following signed
density coefficients.

| Moment | On $[0,1/3]$ | On $[1/3,1/2]$ |
| --- | ---: | ---: |
| $M_n-P_n$, $n$ even | $2-9v$ | $-1$ |
| $M_n+P_n$, $n$ odd | $2-7v$ | $2v-1$ |
| $M_n-P_n$, $n$ odd | $2-v$ | $5(1-2v)$ |

The first density has one sign change, from positive to negative,
and integral zero. The second has the same sign pattern and
positive integral. Weighting by the decreasing function
$\cos^n(\pi v)$ preserves a positive moment for $n>0$ in both
cases. The third density is nonnegative. The even sum is positive
from the original measures. Thus the kernel has a positive rank-one
series and the full row weight has a nonnegative series as before.

In the even difference, the positive mass on $[0,1/10]$ is
$31/200$ and has $\cos^2(\pi v)\ge9/10$. All negative mass starts
at $2/9$, where $\cos^2(2\pi/9)<3/5$. The latter follows, for
example, from the alternating Taylor bound for cosine and
$\pi>157/50$. The remaining positive mass has squared cosine
at least as large as on the negative support. Hence

\[
 M_{2l}-P_{2l}\ge\frac{31}{200}[(9/10)^l-(3/5)^l],\qquad l\ge1.
\]

With the notation of (GH3), summing the even series as in (GH9)
gives

\[
 J(t)\ge\frac{31}{800}e^{-\kappa x^2}
                   [e^{9\kappa x^2/10}-e^{3\kappa x^2/5}]
       \ge\frac{31}{3200}e^{-\kappa x^2/10}
       \quad(\rho\ge1/2,\ x\ge2).
\]

The last inequality uses $1-e^{-2/5}>1/4$. It gives a trace
tail ratio $2\rho/(1+9\rho^2/10)$ in the analogue of (GH10).

## Strictness at a nonpositive radial bound

For $a>0$ the negative Taylor ratio strictly exceeds $M(a)$ when
$0<|p_i-q_i|<1$. Therefore equality in (GH2) with the sharp
coefficient requires one of the following pairs of mass profiles.

- Both are uniform.
- Both are the same permutation of $(0,1/2,1/2)$.
- Both are pure labels, with different labels in the two profiles.

The list follows by combining the equality cases for the sum
$p+q$ with the requirement $p_i-q_i\in\{0,1,-1\}$. A balanced
radial profile that attains equality for almost every pair of radii
must be uniform almost everywhere. Indeed, the list first forces
its values into the displayed finite collection. Any pure label
of positive radial measure violates the allowed pairs on the
product of that set with itself. A nonuniform half-mass profile
of positive measure permits only its own type at other radii and
cannot be balanced. Thus a nonpositive sharp scalar radial form
already implies strict comparison for every nonuniform balanced
mass profile, even when the form has null directions. Uniform
mass profiles have the simplex as their only equality case by DG6.

## Directed extension through correlation 11/20

For every ordinary balanced Gaussian three-cell partition $A$ in
$\mathbb R^d$, $d\ge2$, and every $0<\rho\le11/20$,

\[
                           S_\rho(A)\le S_\rho(T),
\]

where $T$ is a regular simplex partition. On $[53/100,11/20]$,
equality requires a simplex partition up to rotation, relabeling
and null sets. The interval through $53/100$ is supplied by the
[earlier unrestricted comparison](THREE_CELL_GLOBAL_SCALAR_THROUGH_HALF.md).

For the finite block (GH6) of the new scalar operator, take
$m=24$ and $N=96$. The square root of the block sum bounds the
quadratic form relative to $\int J_m z_i^2d\nu$ and then to
$\int J z_i^2d\nu$. In the trace bound for the remaining modes,
retain the small-radius argument of (GH10) for $x\le2$ and use
$J(t)\ge(31/3200)e^{-\kappa x^2/10}$ for $x\ge2$. With
$\eta=1+19\kappa/10$ and $r=\gamma/\eta$, this gives the tail

\[
 \frac{16(4\gamma)^{N+1}}
      {b_1(N+1)![1-4\gamma/(N+2)]}
       +\frac{3200r^{N+1}}{31\sigma\eta(1-r)}.              \tag{SP5}
\]

The [directed checker](check_three_cell_sharp_positive_arc_envelope.py)
encloses the sum of the finite norm and (SP5) by $99/50$.
It uses the radial moment enclosure (GH8) with radius $12$,
16384 equal radial intervals and 100-bit Arb precision. It covers
$[53/100,27/50]$ by 40 correlation balls of width $1/4000$ and
$[27/50,11/20]$ by 80 of width $1/8000$. The entry point reuses
the radial moment enclosure in the
[through-half checker](check_three_cell_global_scalar_through_half.py).

The planar stability excess is consequently bounded by
$-(1/50)\sum_i\int J(t)|p_i(t)-1/3|^2d\nu(t)$ on this new
interval. Equality and reduction to the plane for global
maximizers follow as in the earlier unrestricted comparison.
The directed run has internal verification only.

### Discrete plurality

For every fixed $0<\rho\le11/20$ and $\varepsilon>0$, there exist
an integer $D$ and $\tau>0$ such that for every $m\ge1$ and
$f:[3]^m\to\Delta_3$ with $\mathbb E f=\mathbf1/3$ and
$\operatorname{Inf}_i^{\le D}(f_j)\le\tau$ for all $i,j$,

\[
                   S_\rho(f)\le
      \lim_{k\to\infty}S_\rho(\operatorname{PLUR}_{k,3})+\varepsilon.
\]

The inputs are uniform and each coordinate is independently retained
with probability $\rho$ or resampled uniformly. The influence is
the squared Fourier mass in degrees at most $D$ on supports
containing $i$. The discrete transfer in the earlier note applies
with the larger Gaussian interval by the same invariance theorem
and plurality limit. [Isaksson--Mossel, Theorem 7.1 and Lemma 7.3](https://arxiv.org/pdf/0903.3362)

## A row-sum comparison through correlation 29/50

The [second directed checker](check_three_cell_sharp_arc_row_extension.py)
uses the same quadratic comparison on $[11/20,29/50]$, with $m=24$
and $N=112$. For the finite block, take the Gram matrix
\[
 C_{ij}=\sqrt{k_i k_j}\big(T_{i+j}-U_iU_j/W\big),\qquad 0\le i,j\le N.
\]
It is the Gram of the weighted centered features in (GH6), so
its largest eigenvalue equals the norm of the positive finite
block. For any strictly positive vector $v$, each eigenvalue lies
in absolute value below
\[
        \max_i v_i^{-1}\sum_j |C_{ij}|v_j.                       \tag{SP6}
\]
Indeed, at a coordinate maximizing $|x_i|/v_i$ of a corresponding
eigenvector $x$, the eigenvector equation gives the bound. A short
floating point power iteration of entrywise upper bounds selects
the weights. Only the ensuing inequality (SP6), evaluated in Arb
with the resulting exact dyadic weights, is used for certification.

For the radial moments, on a panel $[l,r]$ interpolate
$F(x)=8\sigma e^{-\eta x^2}/Q(x)$ linearly between the endpoints.
Here $\eta$ runs over $\beta,1,1-\kappa$, and $Q(x)$ is the
positive-coefficient polynomial below (GH7). In bounding the
interpolation error use $F(x)\le F(l)$, and set
$D_1=Q'(r)/Q(l)$, $D_2=Q''(r)/Q(l)$. The logarithmic derivative is
$-2\eta x-Q'(x)/Q(x)$. Differentiating once more and using positivity
of the coefficients gives
\[
 \sup_{[l,r]} |F''|\le F(l)[(2\eta r+D_1)^2+2\eta+D_2+D_1^2].    \tag{SP7}
\]
The error in linear interpolation is at most $(r-l)^2/8$ times
this bound. Multiply by $x^j$ and integrate the interpolant and
constant error bound exactly. In particular its two endpoint weights
are $(r\int_l^r x^jdx-\int_l^r x^{j+1}dx)/(r-l)$ and the complementary
weight $(\int_l^r x^{j+1}dx-l\int_l^r x^jdx)/(r-l)$. The radial
enclosure uses 2048 panels on $[0,12]$ and (GH8) outside. Combine
the resulting bound (SP6) with (SP5) for the infinite positive tail.

The correlation cover starts with thirty intervals of width $1/1000$.
Subdivide each into exact dyadic children until the sum of (SP6)
and (SP5) is strictly below $1999/1000$. Each accepted child is
verified by 100-bit Arb enclosures on the entire correlation ball.
The full replay accepts 129 balls and bounds the norm sum from above
by $1.998766$.
The norm bound on this cover gives the Gaussian comparison on
$[11/20,29/50]$, with weighted stability deficit bounded below by
$\frac1{1000}\sum_i\int J(t)|p_i(t)-1/3|^2d\nu(t)$ in the plane.
Equality rigidity and the passage to all $d\ge2$ follow as in the
previous interval. The adjoining intervals give unrestricted
optimality for every $0<\rho\le29/50$. The same discrete plurality
statement holds for every fixed correlation in this interval with
the stated balance and low-degree influence conditions.
