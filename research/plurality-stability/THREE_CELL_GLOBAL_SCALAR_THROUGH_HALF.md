# An unrestricted three-cell comparison through correlation 53/100

Let $T$ be a regular simplex partition into three equal Gaussian
cells. For every ordinary measurable balanced three-cell partition
$A$ in $\mathbb R^d$, $d\ge2$, and every $0<\rho\le53/100$,

\[
                         S_\rho(A)\le S_\rho(T).
\]

For $11/25\le\rho\le53/100$ equality holds only for a simplex partition
up to rotation, relabeling and null sets. The new argument compares
all planar partitions without a circular mass restriction. It uses a
directed scalar estimate with analytic tails. The result has internal
verification only. [Directed checker](check_three_cell_global_scalar_through_half.py)

## An angular inequality for all circular masses

Write $a\ge0$, $k_a(v)=e^{a\cos(\pi v)}$ and let $G_a$ be its even
second primitive with $G_a(0)=G_a'(0)=0$. Centered arcs of masses
$p,q$ have unnormalized angular pairing
$B_a(p,q)=[G_a(p+q)-G_a(p-q)]/2$. Circular rearrangement bounds the
pairing for arbitrary measurable sets of masses $p,q$ by $B_a(p,q)$.
The primitive identity follows by differentiating in $p,q$ and
using the zero boundary values. [Circular rearrangement and its application in DG6](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md)

Put $b=1/3$ and define

\[
 \begin{split}
 P_0(a)&=\frac92\int_0^{2/3}v k_a(v)\,dv,\qquad E_0(a)=e^{-a/2},\\
 P_*(a)&=\frac23 P_0(a)+\frac13E_0(a),\qquad
 M(a)=2\int_0^1(1-v)k_a(v)\,dv.
 \end{split}                                                   \tag{GH1}
\]

For $p,q\in\Delta_3$, write $z=p-b\mathbf1$ and $w=q-b\mathbf1$.
Then

\[
 \sum_i B_a(p_i,q_i)-3B_a(b,b)
 \le\frac14\big[P_*(a)|z+w|^2-M(a)|z-w|^2\big].             \tag{GH2}
\]

To prove this, for $\delta\in[-2/3,4/3]$ let

\[
 \mathcal R_a(\delta)=2\int_0^1(1-u)k_a(2/3+u\delta)\,du.
\]

The Taylor remainder of $G_a$ at $2/3$ is
$\delta^2\mathcal R_a(\delta)/2$. On negative displacements,
$\mathcal R_a(\delta)\le\mathcal R_a(-2/3)=P_0(a)$ by monotonicity
of the kernel on $[0,1]$. On $0\le\delta\le2/3$ the cosine in
the integrand is at most $-1/2$, giving $\mathcal R_a(\delta)\le E_0(a)$.

For $2/3<\delta\le4/3$, put $v_0=4/3-\delta$. In the part of
the integral where $2/3+u\delta>4/3$, reflect the kernel argument
about $1$ to obtain a density $2(v-v_0)/\delta^2$ on $[v_0,2/3]$.
On this interval

\[
 1-v_0/v\le1-3v_0/2\le(1-3v_0/4)^2=(3\delta/4)^2,
\]

so the density is at most $(9/8)v$. Subtracting $E_0(a)$ before
estimating this part gives
$\mathcal R_a(\delta)\le E_0(a)+(P_0(a)-E_0(a))/4$.
The rest of the integral contributes a nonpositive difference
from $E_0(a)$.

Apply these bounds to $\delta=z+w$, a vector of sum zero. If it
has one negative coordinate, every positive coordinate is at most
$2/3$ and the negative coordinate carries at most two thirds of
$|\delta|^2$. If it has two negative coordinates, they carry at
most half of $|\delta|^2$. Including the possible contribution from
the positive coordinate then gives the coefficient $P_0(a)$ a
weight at most $1/2+(1/2)/4=5/8$. In both cases
$\sum_i\delta_i^2\mathcal R_a(\delta_i)\le P_*(a)|\delta|^2$
because $P_0(a)\ge E_0(a)$. The linear terms at $2/3$ cancel.
For the other remainder, $|p_i-q_i|\le1$ and monotonicity gives
$G_a(p_i-q_i)\ge (p_i-q_i)^2 M(a)/2$. These inequalities prove (GH2).

## Positive radial series and mean-zero centering

Write $\sigma=1-\rho^2$, $\kappa=\rho^2/\sigma$,
$\gamma=2\rho/\sigma$, $\beta=1+\kappa$ and $d\nu(t)=e^{-t}dt$.
Use $x=\sqrt t$, $y=\sqrt s$ and $a=\gamma xy$. A circular
mass profile $p(t)$ of a balanced planar partition gives
$z_i(t)=p_i(t)-b$ with $\int z_i\,d\nu=0$. Define the kernels
and row weight

\[
 \begin{split}
 K(x,y)&=\frac{e^{-\kappa(x^2+y^2)}}{2\sigma}[M(a)+P_*(a)],\\
 J(t)&=\int \frac{e^{-\kappa(x^2+y^2)}}{4\sigma}[M(a)-P_*(a)]\,d\nu(s).
 \end{split}                                                    \tag{GH3}
\]

Expanding the squares in (GH2) and integrating in the two radii
bounds the stability excess by

\[
 \sum_i\left[\iint K(\sqrt t,\sqrt s)z_i(t)z_i(s)d\nu(t)d\nu(s)
                        -2\int J(t)z_i(t)^2d\nu(t)\right].   \tag{GH4}
\]

Let $M_n,P_n$ be the moments of $c^n$, with $c=\cos(\pi v)$, in
the two probability measures in $M$ and $P_*$. They are the
coefficients of $a^n/n!$ in (GH1). Their densities on $v$ are
$2(1-v)$ on $[0,1]$ for $M_n$, and $3v$ on $[0,2/3]$ plus an
atom of mass $1/3$ at $2/3$ for $P_n$. Then $M_n+P_n>0$ for
all $n\ge0$, and $M_n-P_n>0$ for $n\ge1$.

To verify the signs, fold at $v=1/2$ using $c(1-v)=-c(v)$.
For even $n$, positivity of the sum is immediate. For the
difference, the folded density is $2-3v$ on $[0,1/3]$ and $-1$
on $[1/3,1/2]$, with an atom of mass $-1/3$ at $1/3$. The
positive mass is $1/2$, equal to the total negative mass, so
weighting by $c(v)^n$ gives a strict positive sign for even $n>0$.
For odd $n$, the sum has density $2-v$ on $[0,1/3]$, density
$2v-1$ on $[1/3,1/2]$ and atom $-1/3$ at $1/3$. Its moment
is at least $(1/2)^n(11/18-1/36-1/3)>0$. The odd difference
has density $2-7v$ on $[0,1/3]$, density $5(1-2v)$ on
$[1/3,1/2]$ and atom $1/3$ at $1/3$. The first density has
integral $5/18$ and one sign change from positive to negative.
Weighting it by the decreasing function $c(v)^n$ preserves its
positive sign, and the remaining contributions are nonnegative.

Consequently $K$ is a sum of positive rank-one kernels
$k_n f_n(x)f_n(y)$, with

\[
 f_n(x)=e^{-\kappa x^2}x^n,\qquad
 k_n=\frac{\gamma^n(M_n+P_n)}{2\sigma n!}.
\]

Truncating the nonnegative series of $J$ after $m$ terms gives

\[
 J(t)\ge J_m(t)=\frac{e^{-\kappa x^2}}{4\sigma}
                    \sum_{n=1}^m b_n x^n,\qquad
 b_n=\frac{\gamma^n(M_n-P_n)\Gamma(1+n/2)}{n!\,\beta^{1+n/2}}.
                                                               \tag{GH5}
\]

Put

\[
 T_j=\int\frac{x^j e^{-2\kappa x^2}}{J_m(t)}d\nu(t),\quad
 U_j=\int\frac{x^j e^{-\kappa x^2}}{J_m(t)}d\nu(t),\quad
 W=\int\frac{1}{J_m(t)}d\nu(t).
\]

For the mean-zero functions in (GH4), replacing $f_n$ by
$f_n-U_n/W$ preserves the quadratic form. The squared
Hilbert--Schmidt norm of the first $N+1$ terms on the space
weighted by $J_m$, is

\[
 \sum_{i,j=0}^{N}k_i k_j(T_{i+j}-U_iU_j/W)^2.              \tag{GH6}
\]

For $\rho\le1/2$ take $m=12$ and $N=48$. The positive tail has
operator norm at most its trace. Its centered variance for mode $n$
is at most $T_{2n}$. Using the highest term
of (GH5), $k_n\le\gamma^n/(\sigma n!)$ and
$\gamma/\beta=2\rho\le1$ gives the trace bound

\[
 \begin{split}
 \sum_{n>N} k_n T_{2n}
 &\le\frac{4\beta^{m/2-1}}{b_m}
                 \sum_{n>N}\frac{\Gamma(n+1-m/2)}{n!}\\
 &=\frac{4\beta^{m/2-1}}{b_m}
             \frac{(N-m/2+1)!}{(m/2-1)N!}.
 \end{split}                                                    \tag{GH7}
\]

The last sum telescopes. These bounds ensure convergence of the
centered operator and justify the series calculation on the stated
weighted space.

## Directed bounds through correlation one half

The angular moments are evaluated using
$\cos^n\theta=2^{-n}\sum_{j=0}^n\binom nj\cos((n-2j)\theta)$.
For a positive integer $l$, the cosine integral in the measure of
$M$ is $2(1-(-1)^l)/(l\pi)^2$ and the cosine integral in the
measure of $P_0$ is

\[
 \frac{3\sin(2l\pi/3)}{l\pi}
 +\frac{9(\cos(2l\pi/3)-1)}{2(l\pi)^2}.
\]

Both integrals equal one at $l=0$. The exact angles in these
formulas reduce to rational numbers and $\sqrt3$.

For the radial moments in (GH6), write
$Q(x)=\sum_{n=1}^m b_n x^{n-1}$. Each integral has the form
$8\sigma\int_0^\infty e^{-\eta x^2}x^j/Q(x)\,dx$ with
$\eta=\beta,1$ or $1-\kappa$. The checker encloses the integral
on $[0,10]$ over 8192 equal intervals. For $L=10$ and
$\alpha=\max(0,(j-1)/2)$, the bound on the omitted part is

\[
 \frac{4\sigma L^{j-1} e^{-\eta L^2}}
      {b_1(\eta-\alpha/L^2)}.                                  \tag{GH8}
\]

Here $Q(x)\ge b_1$. The change of variable $u=x^2$ and the
bound $u^{(j-1)/2}\le L^{j-1}e^{\alpha(u-L^2)/L^2}$ for
$u\ge L^2$ give (GH8).

At 100-bit Arb precision the checker bounds the sum of the square
root of (GH6) and (GH7) by $193/100$ at $\rho=1/2$ and by
$194/100$ on $[11/25,1/2]$. The interval run covers 240 adjacent
correlation balls of width $1/4000$. Thus (GH4) gives

\[
 S_\rho(A)-S_\rho(T)
 \le -\frac3{50}\sum_i\int J_m(t)|p_i(t)-1/3|^2d\nu(t)
 \qquad (11/25\le\rho\le1/2).
\]

If equality holds, the masses are balanced on almost every circle.
In this case the first cosine term in DG6 has strictly positive
radial weight. Equality forces each first moment to saturate its
circular cap and the corresponding first moments to align across
almost every pair of radii. Saturation gives single arcs, and
alignment fixes the labeled arcs across radii. Three equal disjoint
arcs give the simplex. [DG6](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md)

## Exponential control above correlation one half

For $1/2\le\rho\le53/100$, use $m=16$ and $N=80$ in the
finite block (GH6). Its bound is relative to $\int J_m z_i^2d\nu$
and therefore to $\int J z_i^2d\nu$. For the tail, use a lower
bound on the full weight $J$.

In the folded signed measure for $M_{2l}-P_{2l}$, the positive
mass on $[0,1/6]$ is $7/24$ and has $c(v)^2\ge3/4$. The rest
of the positive mass has $c(v)^2\ge1/4$ and the negative mass
has $c(v)^2\le1/4$. Hence

\[
 M_{2l}-P_{2l}\ge\tfrac7{24}[(3/4)^l-(1/4)^l],\qquad l\ge1.
\]

The even terms of the full row series, with coefficients as in (GH5), and the bound
$l!/(2l)!\ge1/(4^l l!)$ from $\binom{2l}{l}\le4^l$ give

\[
 J(t)\ge\tfrac7{96}e^{-\kappa x^2}
             [e^{3\kappa x^2/4}-e^{\kappa x^2/4}].         \tag{GH9}
\]

Here $\sigma\beta=1$ and $\gamma^2/\beta=4\kappa$.
If $x\ge2$ and $\rho\ge1/2$, then $\kappa x^2/2\ge2/3$
and $1-e^{-\kappa x^2/2}\ge1/3$. Therefore
$J(t)\ge (7/288)e^{-\kappa x^2/4}$ on this range.

The trace of the uncentered tail relative to $J$ is bounded by
$\sum_{n>N}k_n\int f_n^2/J\,d\nu$. On $0\le x\le2$ use
$J\ge J_1$ from (GH5) and $k_n\le\gamma^n/(\sigma n!)$.
The contribution of a mode is at most
$16(4\gamma)^n/(b_1 n!)$. On $x\ge2$ use (GH9) and integrate
the resulting Gaussian monomial over $[0,\infty)$. With
$\eta=1+7\kappa/4$ and $r=\gamma/\eta$, the contribution is
at most $288r^n/(7\sigma\eta)$. Summing the two series bounds
the trace by

\[
 \frac{16(4\gamma)^{N+1}}
      {b_1 (N+1)![1-4\gamma/(N+2)]}
       +\frac{288r^{N+1}}{7\sigma\eta(1-r)}.              \tag{GH10}
\]

The `--extension` checker mode verifies the denominators in (GH10)
and bounds the sum of the finite norm from (GH6) and this tail by
$199/100$. It uses 16384 equal radial intervals on $[0,12]$, then
(GH8) with $L=12$, at 100-bit precision. It covers $[1/2,13/25]$
by 80 correlation balls of width $1/4000$ and $[13/25,53/100]$
by 80 of width $1/8000$. Combining the finite norm against $J_m$
and the tail against $J$ bounds the stability excess by
$-\frac1{100}\sum_i\int J(t)|p_i(t)-1/3|^2d\nu(t)$.
The same equality argument applies.

For $d\ge2$, existence and dimension reduction for balanced
maximizers reduce the global comparison on the new interval to
the planar result, as in DG4. The earlier unrestricted interval
through $11/25$ from DG10 completes the stated correlation range.
[Heilman--Tarter, Lemma 2.3 and Theorem 1.9](https://doi.org/10.1017/fms.2021.56)
[DG4 and DG10](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md)

## Discrete plurality

Fix $0<\rho\le53/100$ and $\varepsilon>0$. There are an integer
$D$ and a number $\tau>0$ such that the following holds for every
$m\ge1$. Let $X$ be uniform on $[3]^m$ and obtain $Y$ from $X$
by independently retaining each coordinate with probability $\rho$
and resampling it uniformly with probability $1-\rho$. Let
$f:[3]^m\to\Delta_3$ have $\mathbb E f(X)=\mathbf1/3$ and
$\operatorname{Inf}_i^{\le D}(f_j)\le\tau$ for all $i,j$.
Here the low-degree influence is the sum of squared product Fourier
coefficients of degrees at most $D$ whose supports contain $i$.
Then

\[
 \mathbb E\langle f(X),f(Y)\rangle
   \le \lim_{k\to\infty} S_\rho(\operatorname{PLUR}_{k,3})+\varepsilon.
\]

Isaksson and Mossel's invariance result produces an ordinary balanced
Gaussian partition in $\mathbb R^{2m}$ whose stability differs from
that of $f$ by at most $\varepsilon$. The Gaussian comparison here
bounds its stability by that of the simplex. The Gaussian limit of
discrete plurality has exactly this simplex stability.
[Isaksson--Mossel, Theorem 7.1 and Lemma 7.3](https://arxiv.org/pdf/0903.3362)
