# An unrestricted three-cell comparison through correlation three fifths

For every ordinary balanced Gaussian three-cell partition $A$ in
$\mathbb R^d$, $d\ge2$, and every $0<\rho\le3/5$,

\[
                         S_\rho(A)\le S_\rho(T),
\]

where $T$ is a regular simplex partition. On the new interval
$[29/50,3/5]$, equality requires a simplex up to rotation, relabeling
and null sets. The proof handles arbitrary planar circular masses
and angular shapes. Comparison through $29/50$, the reduction to the
plane for global maximizers, and the exact scalar kernel are supplied
by the preceding results. [Sharp aligned-arc comparison](THREE_CELL_SHARP_POSITIVE_ARC_ENVELOPE.md)

The [directed checker](check_three_cell_quartic_extension.py) covers
the new interval. The result has internal verification only.

## A quartic cost for differing circular masses

Use the notation of GH2 with the sharp positive coefficient
$P_\sharp$ of SP1. Write $d=1-\rho^2$, $\kappa=\rho^2/d$,
$\gamma=2\rho/d$, $x=\sqrt t$ and $d\nu(t)=e^{-t}dt$.
Let $K,J$ denote the scalar kernel and row weight formed from
$M,P_\sharp$ as in GH3, and let $I_0(a)=\int_0^1 k_a(v)dv$.

For $0\le s\le1$, the function
$2G_a(\sqrt s)/s=2\int_0^1(1-u)k_a(u\sqrt s)du$ is convex, with its
continuous value at zero. Indeed, $s\mapsto\cos(\pi u\sqrt s)$
is convex on this interval. The sign of its second derivative
is the sign of $\sin y-y\cos y\ge0$ for $0\le y\le\pi$.
Composition with the increasing convex exponential preserves convexity.
At $s=1$ the function has value $M(a)$ and derivative $I_0(a)-M(a)$.
Its tangent yields

\[
 2G_a(h)\ge h^2 M(a)+h^2(1-h^2)[M(a)-I_0(a)],\qquad |h|\le1.
\]

The even moments of the measures in $M$ and $I_0$ agree by folding
at $v=1/2$. The odd moments of the measure in $M$ are positive and
its first moment is $4/\pi^2$. Consequently
$M(a)-I_0(a)\ge4a/\pi^2$. In the three-cell comparison,
$h_i=p_i(t)-p_i(s)$ has sum zero and
$\sum_i h_i^4=(\sum_i h_i^2)^2/2$. Thus the extra negative term
is bounded below by a rank-one radial pairing. If $\mu_\rho$ has
density $\sqrt t e^{-\kappa t}/C_\rho$ relative to $\nu$, with
$C_\rho=\sqrt\pi\,d^{3/2}/2$, then

\[
 \begin{split}
 S_\rho(A)-S_\rho(T)&\le \sum_i q_\rho(z_i)-\mathcal L_\rho(p),\\
 q_\rho(z)&=\langle z,Kz\rangle-2\int Jz^2d\nu,\qquad
 z_i(t)=p_i(t)-1/3,\\
 \mathcal L_\rho(p)&\ge \frac{\rho d}{2\pi}
    \mathbb E\left[|p-q|^2(1-|p-q|^2/2)\right].
 \end{split}                                                   \tag{QE1}
\]

Here $p,q$ in the expectation are independent mass vectors drawn
from the same profile under $\mu_\rho$. The scalar coefficient
follows by multiplying $(1/4d)e^{-\kappa(t+s)}$ by
$4\gamma\sqrt{ts}/\pi^2$, then integrating the two rank-one factors.
Since $p\cdot q\ge0$ and $\mathbb E[|p-q|^2\mid p]=
|p-\mathbb E q|^2+\operatorname{Var}(q)$, we have

\[
 \mathbb E[|p-q|^2(1-|p-q|^2/2)]
 \ge (1-\mathbb E|p|^2)\operatorname{Var}_{\mu_\rho}(p).    \tag{QE2}
\]

The inequality follows by using the disjoint contributions
$(1-|p|^2)/2$ and $(1-|q|^2)/2$ inside $1-|p-q|^2/2$.

## A fixed finite operator on a correlation slab

Take a rational slab $[l,u]\subset[29/50,3/5]$ and write
$d_u=1-u^2$, $d_l=1-l^2$. Let $m=24$ and $N=136$ be the row and
kernel truncation orders. Use the weight $J_\mathrm{mix}$ formed
at correlation $u$ by multiplying the coefficient $b_n$ in (GH5)
by $(l/u)^n$. Let $K_{N,u}$ be the first $N+1$ terms of the
kernel at $u$.

For any $\rho\in[l,u]$, change variables to $r=\sqrt{t/d}$ and
put $t_0=d_u r^2$. Transport a scalar profile to
$\widehat z(t_0)=\sqrt{d/d_u}\,z(t)$ and set
$v(t_0)=\sqrt{J_\mathrm{mix}(t_0)}\widehat z(t_0)$.
The mixed norm is at most $\int J(t)z(t)^2d\nu(t)$.
To see the two comparisons, put

\[
 Q_l(r)=\sum_{n=1}^m
       \frac{(2l)^n(M_n-P_n)\Gamma(1+n/2)}{n!}r^{n-1}.
\]

In the physical variable $t=dr^2$ the corresponding auxiliary
weight is $\frac14 e^{-\rho^2 r^2}rQ_l(r)$, which is bounded
by $J_m(t)$ from monotonicity of each row coefficient in $\rho$
in the scaled variable. Its weighted measure in $r$ has density
$(d/2)r^2 e^{-r^2}Q_l(r)$. In these same coordinates, the finite
kernel pairing on $v$ is a sum of positive rank-one forms with
fixed features $e^{-r^2/2}r^n/\sqrt{Q_l(r)}$. The coefficient of
the $n$th form is $4(2\rho)^n(M_n+P_n)/n!$, which is at most
the coefficient at $u$. Hence the finite form is bounded by the
fixed operator $B$ on $L^2(\nu(dt_0))$ with kernel
$K_{N,u}/\sqrt{J_\mathrm{mix}(t_0)J_\mathrm{mix}(s_0)}$.

The omitted kernel terms are bounded against the full $J$ by the
trace estimate (SP5). Let $\tau$ be its upper bound over the slab.
Since $2-\tau>0$, the preceding comparisons give
$q_\rho(z)\le\langle v,Bv\rangle-(2-\tau)\|v\|^2$.

The physical mean-zero condition becomes orthogonality of $v$ to
$e^{-(d/d_u-1)t_0}/\sqrt{J_\mathrm{mix}(t_0)}$. Put
$w=1/\sqrt{WJ_\mathrm{mix}}$ with $W=\int J_\mathrm{mix}^{-1}d\nu$
and let $P$ project off $w$. By $1-e^{-x}\le x$ for $x\ge0$,

\[
 |\langle v,w\rangle|\le L\|v\|,
 \qquad L=\frac{u^2-l^2}{d_u}\sqrt{\frac{\int t_0^2/J_\mathrm{mix}\,d\nu}{W}}.
                                                               \tag{QE3}
\]

Writing $b_w=\|PBw\|$ and $k_w=\langle w,Bw\rangle$ gives
$\langle v,Bv\rangle\le\langle Pv,BPv\rangle+
(2b_wL+k_wL^2)\|v\|^2$.

## A certified adverse direction

At the anchor $u$ write $f_n(t_0)=e^{-\kappa_u t_0}t_0^{n/2}$,
$U_n=\int f_n/J_\mathrm{mix}\,d\nu$ and
$g_n=(f_n-U_n/W)/\sqrt{J_\mathrm{mix}}$.
Their Gram matrix $C$ has entries $T_{i+j}-U_iU_j/W$ as in GH6.
The trial uses the fixed rational coefficient vector $\alpha$ in
the checker, supported in degrees at most 24. Put
$D=\alpha^TC\alpha$ and $v_*=\sum\alpha_i g_i/\sqrt D$.
It is a unit vector in the anchored hyperplane. Write
$h=C\alpha$, and denote the kernel coefficients at $u$ by $k_i$.
The Rayleigh quotient and a cross term bound are

\[
  a=\frac{\sum_i k_i h_i^2}{D},\qquad
  R^2=\frac{b^TCb}{D},\qquad
  b_i=k_i h_i-\lambda\alpha_i,
\]

where $\lambda$ is the midpoint of the enclosure of $a$. The square
root of $\sum_{i,j}k_i k_j(C_{ij}-h_i h_j/D)^2$ bounds the
complement norm, denoted by $\Lambda$. To compute $b_w$ and $k_w$,
use the coefficient vector $\beta_i=k_i U_i/\sqrt W$.
Then $b_w^2=\beta^TC\beta$ and $k_w=\sum k_i U_i^2/W$.

For the three cell profiles decompose $v_i=\alpha_i'v_*+g_i'$,
with $g_i'\perp v_*$, and put $A=|\alpha'|$ and
$G^2=\sum_i\|g_i'\|^2$. Write $\ell=2b_wL+k_wL^2$.
The complement bound applies to $Pg_i'$, whose norm is at most
$\|g_i'\|$. Therefore

\[
 \sum_i q_\rho(z_i)\le e A^2+2RAG-\delta G^2,
 \qquad e=a-2+\tau+\ell,\qquad \delta=2-\tau-\ell-\Lambda . \tag{QE4}
\]

## Closing the finite-amplitude comparison

The physical reference is
$\phi_\rho(t)=\sqrt{d_u/d}\,v_*(t_0)/\sqrt{J_\mathrm{mix}(t_0)}$,
so $z=\alpha'\phi_\rho+\epsilon$ with the residual norm $G$
in the auxiliary weighted measure. Its pointwise pole is controlled
by the anchor relation
$\sqrt t\,\phi_\rho(t)=\sqrt{t_0}\,v_*(t_0)/\sqrt{J_\mathrm{mix}(t_0)}$.
For $x\le X=1/4$, the checker bounds this below by $c-Cx$,
using $C=13/10$ and the fact that $\sqrt{t_0}\le x$. The auxiliary
weight times $e^{-x^2}$ is at least $j x$ on this band, where
$j=b_1^\mathrm{mix}e^{-X^2/d_u}/(4\sqrt{d_ud_l})$.

The density of $\mu_\rho$ relative to the auxiliary weighted
measure is bounded by $B_0=16\pi/(11l d_u)$, using $M_1-P_1=11/(2\pi^2)$
from the sharp coefficient and the first term of $Q_l$. Under the
change to $r$ the measure $\mu_\rho$ has density
$(4/\sqrt\pi) r^2 e^{-r^2}dr$. Hence the reference variance has
lower bound $V$ from the anchor variance scaled by $d_u/d_l$.
The checker uses the anchor core $0\le\sqrt{t_0}\le3$.
For a core of mass $m_0$, the variance over all radii is at least
$m_0\int_\mathrm{core}\phi^2d\mu_u-(\int_\mathrm{core}\phi d\mu_u)^2$
by the two-sample variance identity.

If the scalar upper bound can be nonnegative and $A>0$, put
$r_e=G/A$. With $\delta>0$, (QE4) bounds $r_e$ by the positive
root of $e+2Rr-\delta r^2=0$. The probability simplex about
$\mathbf1/3$ has circumradius $R_c=\sqrt{2/3}$. On the origin
band we have

\[
 \frac{G^2}{A^2}\ge2j\int_0^X
                  [c-(C+R_c/A)x]_+^2dx.
\]

For an upper cap $r$ on $r_e$, the resulting bound on $A$ is

\[
 A\le A_r=\frac{3R_c r^2}{2jc^3-3Cr^2},                  \tag{QE5}
\]

provided the denominator is positive and $3r^2/(2jc^2)\le X$.
At $A_r$ the integral of the squared affine function up to its zero
equals $r^2/(2j)$. The integrand increases with $A$, so any larger
$A$ is excluded by the same origin window.

For $r_e\le r$, auxiliary weighted orthogonality gives

\[
 \begin{split}
 1-\mathbb E_{\mu_\rho}|p|^2
     &\ge 2/3-B_0 A_r^2(1+r^2),\\
 \operatorname{Var}_{\mu_\rho}(p)
     &\ge A^2(\sqrt V-\sqrt{B_0}\,r)_+^2.
 \end{split}
\]

The purity bound follows from $\mathbb E|p|^2=1/3+\mathbb E|z|^2$.
For the variance, subtract the $\mu_\rho$ mean and use the triangle
inequality with residual norm at most $\sqrt{B_0}G$.

On each correlation slab the checker bounds $e,R$ above and
$\delta$ below. It tests 64 subintervals for $r_e$ up to the
positive root. On a subinterval $[s,r]$ it bounds the scalar
coefficient by $e+2Rr-\delta s^2$. If this is nonnegative, it
checks that the product of the two lower bounds above, multiplied by
$u d_u/(2\pi)$, is positive and strictly larger than the scalar
coefficient. Here $u d_u\le\rho d$ on the stated interval.
This proves strict comparison for any nonuniform mass profile.
When $A=0$, (QE4) is already strictly negative unless $G=0$.
The DG6 equality argument applies to uniform mass profiles.

The checker uses 140-bit Arb enclosures, the piecewise linear
radial moment enclosure of the row extension checker on $[0,12]$
with 8192 panels, and the analytic moment tails in GH8. The
correlation cover begins with 20 equal slabs and bisects any slab
whose directed inequalities require refinement. The `--point-only`
option checks $\rho=3/5$ independently of the contiguous cover.

## Discrete plurality

For each fixed $0<\rho\le3/5$ and $\varepsilon>0$, there are an
integer $D$ and $\tau>0$ such that for all $m\ge1$, any balanced
$f:[3]^m\to\Delta_3$ with $\operatorname{Inf}_i^{\le D}(f_j)\le\tau$
for all $i,j$ satisfies

\[
 S_\rho(f)\le\lim_{k\to\infty}S_\rho(\operatorname{PLUR}_{k,3})+\varepsilon.
\]

Here the noise independently retains each coordinate with probability
$\rho$ and otherwise resamples it uniformly. Balance is
$\mathbb E f=\mathbf1/3$, and low-degree influence sums the squared
product Fourier coefficients of degrees at most $D$ on supports
containing $i$. The Gaussian realization in $\mathbb R^{2m}$
preserves balance, and the discrete plurality limit is the simplex
stability, as in the preceding transfer. [Isaksson--Mossel, Theorem 7.1 and Lemma 7.3](https://arxiv.org/pdf/0903.3362)
