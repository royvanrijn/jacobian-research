# A planar circular mass cube near correlation one half

Suppose an ordinary balanced Gaussian partition of the plane into
three cells has circular masses $p_i(t)$ satisfying

\[
 \left|p_i(t)-\frac13\right|\le\frac{11}{50}
 \qquad (i=1,2,3)
\]

for almost every $t=|x|^2/2$. For every
$11/25\le\rho\le53/100$, its noise stability is at most the simplex
stability, with equality only for a rotated or relabeled simplex up to
null sets. The comparison allows arbitrary angular shapes on each
circle. The circular mass hypothesis restricts the theorem.

The [directed checker](check_three_cell_near_half_mass_cube.py) gives
the scalar estimate used below. It has internal verification only.

## The exact arc interaction on the cube

Write $d=1-\rho^2$, $\kappa=\rho^2/d$, $\gamma=2\rho/d$ and
$d\nu(t)=e^{-t}dt$. We use $x=\sqrt t$ and $y=\sqrt s$ as scalar
radial variables and set $a=\gamma xy$. The unnormalized angular
kernel at these radii is $k_a(v)=e^{a\cos(\pi v)}$. Let $G_a$ be
its even second primitive with $G_a(0)=G_a'(0)=0$. The pairing of
centered arcs of masses $p,q$ for $k_a$ is

\[
 B_a(p,q)=\tfrac12[G_a(p+q)-G_a(p-q)].                    \tag{MC1}
\]

Indeed, its mixed derivative is
$[k_a(p+q)+k_a(p-q)]/2$ and it vanishes if either mass is zero.
Circular rearrangement bounds the angular pairing for each cell by
this expression. It applies to every measurable slice by approximation
in circular measure, as in [DG6](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md).

Put $b=1/3$, $z=p-b\mathbf1$, $w=q-b\mathbf1$ and $D=11/25$.
Both $|(z+w)_i|$ and $|(z-w)_i|$ are at most $D$ on the cube.
For a displacement $\delta$ from $2b$, Taylor's formula gives

\[
 G_a(2b+\delta)-G_a(2b)-G_a'(2b)\delta
 =\delta^2\int_0^1(1-u)e^{a\cos(\pi(2b+u\delta))}\,du.
\]

For $v\in[0,D]$, put $F(v)=\cos(2\pi/3-\pi v)+1/2$. On
$[0,1/6]$, convexity and the values at its endpoints give $F(v)\le3v$.
On $[1/6,D]$, concavity bounds $F$ by its tangent at $v=1/4$.
The tangent slope is $\pi\sin(5\pi/12)<76/25$, and its value
at $1/6$ is less than $(76/25)/6$. Thus $F(v)\le(76/25)v$
throughout. For positive $\delta$, the argument $2b+u\delta$
lies in $[2/3,4/3]$, where the cosine is at most $-1/2$.
Together these facts give

\[
 \cos(\pi(2b+u\delta))\le -\tfrac12+(\tfrac{419}{500}+\tfrac12)u,
 \qquad |\delta|\le D,
\]

since $-1/2+(76/25)D<419/500$. Bound the exponential of this
affine majorant by its chords on $[0,1/2]$ and $[1/2,1]$.
Integration against $2(1-u)du$ gives the upper mixture

\[
 E_+(a)=\sum_{j=0}^2 w_j e^{a c_j},\qquad
 (w_0,w_1,w_2)=(5/12,1/2,1/12),\quad
 (c_0,c_1,c_2)=(-1/2,169/1000,419/500).
\]

For the Taylor remainder at zero with $|\delta|\le D$,
Jensen's inequality for the exponential and monotonicity of the
cosine on $[0,\pi D]$ give

\[
 2\int_0^1(1-u)e^{a\cos(\pi u\delta)}du
 \ge \exp\left[a\,\frac{2(1-\cos(\pi D))}{(\pi D)^2}\right]
 \ge E_-(a):=e^{17a/20}.
\]

The linear terms at $2b$ cancel when summed over the cells. Applying
these bounds to (MC1) yields

\[
 \sum_i B_a(p_i,q_i)-3B_a(b,b)
 \le\tfrac14\big[E_+(a)|z+w|^2-E_-(a)|z-w|^2\big].      \tag{MC2}
\]

With $c_-=17/20$, put

\[
 R_c(x,y)=d^{-1}e^{-\kappa(x^2+y^2)+\gamma cxy},\qquad
 K(x,y)=\tfrac12\big(\sum_j w_j R_{c_j}(x,y)+R_{c_-}(x,y)\big),
\]

and $J(t)=\tfrac14\int(R_{c_-}(\sqrt t,\sqrt s)-
\sum_j w_jR_{c_j}(\sqrt t,\sqrt s))d\nu(s)$. Polar integration of (MC2)
bounds the stability excess by

\[
 \sum_i\left[\iint K(\sqrt t,\sqrt s)z_i(t)z_i(s)d\nu(t)d\nu(s)
                   -2\int J(t)z_i(t)^2d\nu(t)\right].    \tag{MC3}
\]

Each $z_i$ has mean zero for $\nu$ by Gaussian balance.

## A directed one-dimensional estimate

Since $|c_j|<c_-$ and $\sum_j w_j=1$, the power series of $K$ is
a sum of positive rank-one kernels. Center its feature vectors at
$x_0=14/25$. The resulting positive kernel has diagonal

\[
 H(x)=K(x,x)-2K(x,x_0)+K(x_0,x_0).
\]

Centering does not change the quadratic form on mean-zero functions.
The same bounds on $c_j$ show that every power series coefficient
of $R_{c_-}-\sum_jw_jR_{c_j}$ is nonnegative. With
$\beta=1+\kappa$, truncating after degree twelve gives

\[
 J(t)\ge J_{12}(t)=\frac{e^{-\kappa x^2}}{4d}P(x),\qquad
 P(x)=\sum_{n=1}^{12}b_nx^n,\qquad
 b_n=\frac{\gamma^n(c_-^n-\sum_jw_jc_j^n)\Gamma(1+n/2)}{n!\,\beta^{1+n/2}}.
\]

Weighted Cauchy--Schwarz for the centered feature vectors bounds
the first term for each cell in (MC3) by
$C_\rho\int J_{12}(t)z_i(t)^2d\nu(t)$, where

\[
 C_\rho=\int_0^\infty
          \frac{8d\,e^{-(1-\kappa)x^2}H(x)}{P(x)/x}\,dx.    \tag{MC4}
\]

The denominator at zero is its continuous limit $b_1>0$.
For interval evaluation, set
$m_c(x)=\exp[(\gamma c/2-\kappa)x^2]$. With the same weights as
in $2K$ (weight one for $c_-$), the centered diagonal is

\[
 2d H(x)=\sum_{c\in\{c_0,c_1,c_2,c_-\}} w(c)
 \left[(m_c(x)-m_c(x_0))^2+2m_c(x)m_c(x_0)
             (1-e^{-\gamma c(x-x_0)^2/2})\right].
\]

The checker evaluates (MC4) on $[0,6]$ by upper interval Riemann
sums of 4096 equal panels with 100-bit Arb arithmetic. It verifies
$K(x_0,x_0)<2$, $g=\gamma c_--2\kappa>0$ and
$M=1+\kappa-\gamma c_->0$ on each correlation ball. For $x\ge6$,
$H(x)\le(d^{-1}+2)e^{g x^2}$ and $P(x)/x\ge b_{12}x^{11}$.
Since $d\le1$ and $x^{-11}\le x/6^{12}$ on the tail, the remaining
integral is at most

\[
 \frac{12 e^{-36M}}{M b_{12} 6^{12}}.
\]

The full bound is $C_{1/2}<44/25$ and $C_\rho<48/25$ on
$[11/25,53/100]$, covered by 90 adjacent correlation balls of
width $1/1000$. In particular

\[
 S_\rho(A)-S_\rho(T)
 \le -\frac2{25}\sum_i\int J_{12}(t)|p_i(t)-1/3|^2d\nu(t).
\]

Equality requires balanced masses on almost every radius. In that
case the first cosine term separated in DG6 has a strictly positive
coefficient at every pair of positive radii. Equality forces each
first moment to saturate its circular cap and the corresponding
moments to align for almost every pair of radii. Saturation gives a
single arc in each cell, and alignment fixes the labeled arcs across
radii. Three disjoint arcs of mass $1/3$ give the simplex.
