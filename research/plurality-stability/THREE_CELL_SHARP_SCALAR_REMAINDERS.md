# Quantitative remainders in the sharp circular arc bound

The sharp angular Taylor coefficients admit two explicit finite-amplitude
remainders. One detects that the sum of the two circular mass vectors
stays away from having an absent label. The other is an exact quartic
penalty for the difference in their mass vectors. Both are available
before using any first angular moment compatibility. These results are
internal analytic bounds.

Use the normalized angle variable of the [sharp positive arc estimate](THREE_CELL_SHARP_POSITIVE_ARC_ENVELOPE.md).
For $a\ge0$ put $k_a(v)=e^{a\cos(\pi v)}$ and write $G_a$ for its even
second primitive with $G_a(0)=G_a'(0)=0$. On the two mass vectors
$p,q\in\Delta_3$ write $z=p-\mathbf1/3$, $w=q-\mathbf1/3$,
$s=z+w$, and $d=z-w$. Use $P_0$, $P_1$ and
$P_\sharp=(2P_0+P_1)/3$ from the sharp estimate and put
\[
 M(a)=2\int_0^1(1-v)k_a(v)dv,\qquad
 I_0(a)=\int_0^1 k_a(v)dv,
 \qquad m=\min_i(p_i+q_i).
\]
Here $M-I_0\ge0$: reflect at $1/2$ and use that $k_a$ decreases.
Set $\eta=1/2+15/(4\pi^2)-3\sqrt3/(2\pi)>0$ and
$c=\min(\pi\sqrt3/24,3\eta/4)>0$.
The positivity of $\eta$ follows for example from
$3.14<\pi<22/7$ and $\sqrt3<7/4$.

## A quantitative positive remainder

Let $\mathcal R_a(u)=2\int_0^1(1-v)k_a(2/3+vu)dv$.
Then for the centered sum of any two circular mass vectors,
\[
 \sum_i s_i^2\mathcal R_a(s_i)
   \le \big[P_\sharp(a)-c a e^{-a/2}m\big]|s|^2 .               \tag{SR1}
\]
If $s=0$ this is immediate. Suppose first that there is exactly
one negative coordinate $-h$ with $0<h\le2/3$. For
$P(h)=(2\mathcal R_a(-h)+\mathcal R_a(h/2))/3$, we have
$P(2/3)=P_\sharp(a)$ and $P(0)=e^{-a/2}$. Writing
$\theta=\pi(2/3-hv)$, differentiating under the integral gives
\[
 P'(h)=\frac{a\pi}{3}\int_0^1 v(1-v)
           \big[4\sin\theta\,e^{a\cos\theta}
               -\sin(\theta/2)\,e^{-a\cos(\theta/2)}\big]dv.
\]
On $0\le\theta\le2\pi/3$, one has
$\sin\theta\ge\sin(\theta/2)$ and
$\cos\theta\ge-\cos(\theta/2)$, so the bracket is at least
$3\sin\theta\,e^{a\cos\theta}$. On $0\le v\le1/2$ this is
at least $(3\sqrt3/2)e^{-a/2}$, and the integral of
$v(1-v)$ over that subinterval is $1/12$. It follows that
$P'(h)\ge(\pi\sqrt3/24)a e^{-a/2}$.

To bound the ratio at a possibly unequal split of the two positive
displacements, subtract $P(h)u^2/2$ from the Taylor remainder of
$G_a$ at $2/3$. The resulting function is concave on the allowed
positive displacements $[0,h]$, since $k_a(2/3+u)\le e^{-a/2}$ there
and $P(h)\ge P(0)$. At fixed negative coordinate $-h$ its sum over
the positive displacements is maximized by splitting $h$ equally.
This makes the sum over all three coordinates zero by the definition
of $P(h)$. Thus the ratio is at most $P(h)$. Integrating its derivative
bound up to $2/3$ gives (SR1) since $m=2/3-h$.

If there are two negative coordinates, the sharp arc proof bounds
the ratio by $(P_0+e^{-a/2})/2$. The difference from $P_\sharp$ is
$\tfrac12[(P_0+2P_1)/3-e^{-a/2}]$. The probability measure in
$(P_0+2P_1)/3$ has mean cosine $-1/2+\eta$, by integrating the
two elementary linear densities. Jensen's inequality and
$e^{a\eta}-1\ge a\eta$ bound the difference below by
$(\eta/2)a e^{-a/2}$. Since $m\le2/3$, this implies (SR1) with
the stated constant.

## The quartic difference remainder and the integrated comparison

For $u\in[-1,1]$,
\[
 G_a(u)\ge \frac{u^2}{2}\big[M(a)+(1-u^2)(M(a)-I_0(a))\big].    \tag{SR2}
\]
For $u\ne0$, write $x=u^2$. The Taylor coefficient $2G_a(u)/u^2$
is $2\int_0^1(1-v)e^{a\cos(\pi v\sqrt x)}dv$. It is convex in
$x\in[0,1]$. Indeed the second derivative of $\cos(\pi v\sqrt x)$
has the sign of $\sin y-y\cos y\ge0$, where $0\le y\le\pi$.
Continuity handles zero and the exponential is increasing and
convex. At $x=1$ the coefficient equals $M(a)$ and its derivative is
\[
 -a\pi\int_0^1 v(1-v)\sin(\pi v)e^{a\cos(\pi v)}dv
                         =-\big(M(a)-I_0(a)\big)
\]
by integration by parts. The tangent bound gives (SR2), which holds
also at $u=0$. For a vector $d$ of zero sum with three entries,
$\sum_i d_i^4=|d|^4/2$. This follows by squaring
$d_1d_2+d_2d_3+d_3d_1=-|d|^2/2$ and using $\sum_i d_i=0$.

Consequently arbitrary circular partitions satisfy the aligned-arc
upper estimate (GH2) with $P_\sharp$, and with the following extra
terms subtracted:
\[
 \frac14\left[ c a e^{-a/2}m |z+w|^2
       +(M(a)-I_0(a))|z-w|^2\big(1-|z-w|^2/2\big)\right].      \tag{SR3}
\]
The first arises from (SR1), the second from (SR2) and the three-label
identity. Both terms are nonnegative since $|p-q|^2\le2$ on the mass
simplex. The quartic factor vanishes at identical profiles and at
two distinct pure labels. Away from these extremal differences it
strengthens the coefficient of the radial Dirichlet term. The first
factor gives an additional loss when two profiles together give every
label positive mass. Integrating (SR3) with the positive radial factor
from (GH3) gives a global sharper comparison for all $0<\rho<1$.
Controlling the sign for every feasible amplitude beyond the scalar
spectral range remains open.

## Retaining the first angular moment as well

The sharp positive Taylor estimate applies to any angular kernel
$g(\cos(\pi v))$ with $g$ increasing, convex and differentiable on
$[-1,1]$. In the proof, monotonicity gives the same negative and
positive displacement comparisons. Convexity and the cosine means
give $P_\sharp(g)\ge g(-1/2)$ and
$(P_0(g)+2P_1(g))/3\ge g(-1/2)$ by Jensen's inequality. The
derivative formula for $P(h)$ now has $g'(\cos\theta)$ and
$g'(-\cos(\theta/2))$ in place of $a e^{a\cos\theta}$ and
$a e^{-a\cos(\theta/2)}$. Nonnegative monotone $g'$ makes the same
argument work, with derivative lower bound
$(\pi\sqrt3/24)g'(-1/2)$. The second case of (SR1) has the same
bound with $g'(-1/2)$ by the tangent inequality for Jensen's gap.
Also $x\mapsto g(\cos(\pi v\sqrt{x}))$ is convex with the same
endpoint derivative identity, giving (SR2) with the linear
functionals $M(g)$ and $I_0(g)$ in place of $M(a)$ and $I_0(a)$.

Separate $\lambda_a\cos(\pi v)$ from the exponential kernel,
where $\lambda_a=ae^{-a}$. The remaining kernel
$g(c)=e^{ac}-\lambda_a c$ is increasing and convex, so it admits
both sharp scalar estimates and the remainders. It is also
nonnegative, so circular rearrangement applies to it. The scalar
coefficient for the cosine alone is
$P_\sharp(\mathrm{id})=-3/(2\pi^2)$ and $M(\mathrm{id})=4/\pi^2$,
where $\mathrm{id}(c)=c$. These follow by integrating the densities
for $P_\sharp$ and $M$ against the cosine.

For the two actual partitions write $u_i,v_i\in\mathbb R^2$ for
their first moments and $a_i(p)=\sin(\pi p_i)/\pi$, $a_*=a_i(1/3)$.
Define
\[
 \begin{split}
 \Phi(p,q)&=3a_*^2-\sum_i a_i(p)a_i(q)
              -\frac{3}{8\pi^2}|z+w|^2-\frac{1}{\pi^2}|z-w|^2,\\
 \Delta(u,v)&=\sum_i a_i(p)a_i(q)-\sum_i u_i\cdot v_i .
 \end{split}                                                     \tag{SR4}
\]
Here $\Delta\ge0$ by the sine caps. To see $\Phi\ge0$ apply the
sharp Taylor estimate to $g(c)=c$, where the aligned cosine pairing
is exactly $a_i(p)a_i(q)$. In (GH2) with $P_\sharp$, we can then
subtract the combined correction
\[
 \begin{split}
 \lambda_a(\Phi(p,q)+\Delta(u,v))+\tfrac14\big[&
  c\,a(e^{-a/2}-e^{-a})m|z+w|^2\; +\\[-2pt]
 &\big(M(a)-I_0(a)-4\lambda_a/\pi^2\big)
           |z-w|^2(1-|z-w|^2/2)\big].                           \tag{SR5}
 \end{split}
\]
Indeed the cosine part has its exact value $\sum_i u_i\cdot v_i$.
Apply (SR1) and (SR2) to the increasing convex remainder, noting
$g'(-1/2)=a(e^{-a/2}-e^{-a})$ and the linearity of the Taylor
coefficients. Each of the terms in (SR5) is nonnegative. In the
quartic term this follows from $M(g)-I_0(g)\ge0$ for increasing $g$.

For an explicit rank-one version of the scalar part $\Phi$, let
$d_\rho=1-\rho^2$, $\gamma=2\rho/d_\rho$, $\kappa=\rho^2/d_\rho$,
and $h(t)=\sqrt t\,e^{-(\kappa+\gamma/2)t}$. The weight multiplying
$\Phi$ after radial integration is at least
$(\gamma/d_\rho)h(t)h(s)$ by $2\sqrt{ts}\le t+s$, as in the
[global first moment comparison](THREE_CELL_GLOBAL_ANGULAR_COMPATIBILITY.md).
Let $\mu$ have density $h/\int h\,d\nu$ relative to $\nu$. The exact
average against $\mu\otimes\mu$ is
\[
 3a_*^2-\sum_i(\mathbb E_\mu a_i(p))^2
          -\frac{11}{4\pi^2}\mathbb E_\mu |z|^2
          +\frac{5}{4\pi^2}|\mathbb E_\mu z|^2 .                \tag{SR6}
\]
Thus SR5 permits retaining both the exact global first-mode
comparisons and additional scalar and quartic terms in the same
upper bound. Neither the separate positivity nor the rank-one
formula resolves the sign of this bound at all correlations.

## Rank-one bounds for the quartic penalty

The integrated quartic remainder in (SR3) alone has a direct
probability expression. The even coefficients in the exponential
series for $M(a)$ and $I_0(a)$ coincide by reflection at $1/2$.
The odd coefficients of their difference are positive, since the
folded density of $M$ is $2(1-2v)$ on $[0,1/2]$ for odd orders.
The coefficients of $a$ and $a^3/6$ are $4/\pi^2$ and
$28/(9\pi^2)$, by the cosine and triple-angle integrals. Thus
$M(a)-I_0(a)\ge4a/\pi^2$. After retaining the first mode as in
(SR5), the remaining factor obeys
$M(a)-I_0(a)-4\lambda_a/\pi^2\ge [28/(9\pi^2)]a^3/6$.

Let $\mu_j$ have density proportional to $t^{j/2}e^{-\kappa t}$
relative to $\nu$, for $j=1,3$. Under two independent draws of
the mass profile from $\mu_j$, put $b=\mathbb E p$,
$s_2=\mathbb E |p-b|^2$, $s_4=\mathbb E|p-b|^4$, and
let $\Sigma$ be the covariance matrix of $p$. The exact average
of the quartic factor is
\[
 \mathcal Q_j=\mathbb E\big[|p-q|^2(1-|p-q|^2/2)\big]
    =2s_2-s_4-s_2^2-2\operatorname{tr}(\Sigma^2)
    \ge s_2\,\mathbb E(1-|p|^2).                                 \tag{SR7}
\]
The identity follows by expanding the difference of two independent
centered vectors. For the inequality use $p\cdot q\ge0$ to bound
$1-|p-q|^2/2\ge[(1-|p|^2)+(1-|q|^2)]/2$ and condition on $p$.
The conditional mean of $|p-q|^2$ is $|p-b|^2+s_2$.

In the upper bound with (SR3), the integrated quartic penalty is
at least $\frac{\rho d_\rho}{2\pi}\mathcal Q_1$. In the combined
bound (SR5), it is at least
$\frac{7\rho^3 d_\rho}{12\pi}\mathcal Q_3$. These coefficients
follow from the separable factors $a^j=\gamma^j(ts)^{j/2}$ and
the normalization
$\int t^{j/2}e^{-\kappa t}d\nu=\Gamma(1+j/2)d_\rho^{1+j/2}$.
This takes account of radial variation and label purity together.
In particular the positive lower bound in (SR7) survives any
pure slices as long as they do not exhaust the measure and the
profile is not constant. Its domination of the positive adverse
radial term beyond the certified ranges is still unresolved.

## A one-mode amplitude criterion using the quartic bound

Write $J_m$ for the positive truncated row weight of the sharp
scalar operator and let $\phi$ be a mean-zero radial trial in
$L^2(J_m\,d\nu)$ of norm one. Decompose a feasible centered
profile as $z=\alpha\phi+g$, where each component of $g$ is
weighted-orthogonal to $\phi$. Set $A=|\alpha|$ and
$G^2=\int J_m|g|^2d\nu$. Suppose the scalar quadratic upper
bound before the quartic correction has the estimate
\[
 q(z)\le e A^2+2RAG-\delta G^2,
          \qquad e,R\ge0,\quad\delta>0.                          \tag{SR8}
\]
The tail of the scalar series may be bounded against the full $J$
before replacing the remaining negative weight by $J_m$. For
$q(z)\ge0$ and $A>0$, (SR8) gives $G/A\le H$ with
$H=(R+\sqrt{R^2+\delta e})/\delta$.

Put $B=\sup(d\mu_1/(J_m,d\nu))$. The first moment difference
is $M_1-P_{\sharp,1}=11/(2\pi^2)$. The first term of $J_m$ and
the normalization of $\mu_1$ then give
\[
    B\le\frac{16\pi}{11\rho(1-\rho^2)} .                       \tag{SR9}
\]
Given a cap $A\le A_0$ for the possible nonnegative profiles,
their purity and variance satisfy
\[
 \begin{split}
 \mathbb E_{\mu_1}(1-|p|^2)&\ge 2/3-B A_0^2(1+H^2),\\
 \mathbb E_{\mu_1}|p-\mathbb E_{\mu_1}p|^2
       &\ge A^2(\sqrt{\operatorname{Var}_{\mu_1}\phi}-\sqrt B H)_+^2.
 \end{split}                                                     \tag{SR10}
\]
For the first, use $|p|^2=1/3+|z|^2$ and weighted orthogonality.
For the second, subtract the $\mu_1$-mean, use the reverse
triangle inequality, and bound the remainder norm by $\sqrt B G$.
Consequently (SR7) makes the global comparison strict whenever
\[
 \frac{\rho(1-\rho^2)}{2\pi}
       [2/3-B A_0^2(1+H^2)]_+
       (\sqrt{\operatorname{Var}_{\mu_1}\phi}-\sqrt B H)_+^2
                         > e+2RH .                              \tag{SR11}
\]
This suffices for every profile: $q(z)<0$ is already strict,
and if $A=0$ a nonnegative scalar upper bound forces $G=0$.
At $z=0$, equality reduces to the usual equal-circle rigidity.

One way to obtain $A_0$ is by a pole lower bound. Suppose
$\phi(x^2)\ge L/x-C$ and $J_m(x^2)e^{-x^2}\ge jx$ for
$0<x\le X$, with $L,j>0$ and $C\ge0$. Every point in the mass
simplex is within $r_*=\sqrt{2/3}$ of its center, so
$|g(x^2)|\ge(A\phi(x^2)-r_*)_+$. If
$x_0=A_0 L/(r_*+A_0 C)\le X$, then any $A\ge A_0$ obeys
\[
  G^2/A^2\ge 2j\int_0^{x_0}x^2(L/x-C-r_*/A_0)^2dx
                  =\frac{2j L^3 A_0}{3(r_*+A_0 C)} .            \tag{SR12}
\]
If this exceeds $H^2$, such $A$ are excluded by $q(z)\ge0$.
The criterion (SR8)–(SR12) is a sufficient route to certifying
a nonlinear interval beyond scalar nonpositivity. It does not
by itself verify the required spectral and trial constants.

The criterion can also keep the dependence on $r=G/A$. With
the same pole hypothesis, every possible pair $(A,r)$ satisfies
\[
 r^2\ge \Psi(A):=
 \frac{2j}{3(C+r_*/A)}
       \big[L^3-(L-(C+r_*/A)X)_+^3\big].                         \tag{SR13}
\]
This is the same clipping calculation integrated over the whole
window, including when the pole would predict a wider clipped set.
The function $\Psi$ is increasing, since its defining integrand
$2j(L-(C+r_*/A)x)_+^2$ is increasing in $A$ at each $x>0$.
Let $A_*(r)$ be any verified upper bound on feasible amplitudes
under (SR13) and the bounded mass-simplex constraint. In (SR11)
use $A_*(r)$ instead of $A_0$ and $r$ instead of $H$, then compare
with $e+2Rr-\delta r^2$ for each $0\le r\le H$. This retains
the whole negative complement term and allows a one-dimensional
verification. In the range where the predicted clipped interval
fits in the window and $2jL^3-3Cr^2>0$, the bound from (SR13) gives
$A\le3r_*r^2/(2jL^3-3Cr^2)$, with constraint $L/(C+r_*/A)\le X$.

One can strengthen the purity factor in (SR11). With
$r_*=\sqrt{2/3}$, $M_2=\mathbb E_{\mu_1}|p-p_*|^2$ and
$m_0=|\mathbb E_{\mu_1}p-p_*|$, we have
$s_2=M_2-m_0^2$ and $|p-\mathbb E_{\mu_1}p|\le r_*+m_0$.
The identity in (SR7) and $\operatorname{tr}(\Sigma^2)\le s_2^2$
give
\[
 \mathcal Q_1\ge [2-(r_*+m_0)^2-3s_2]s_2 .                     \tag{SR14}
\]
If only $M_2$ is known, maximizing the expression being subtracted
over $0\le m_0\le\sqrt{M_2}$ gives the lower coefficient
$2-(r_*+\sqrt{M_2})^2$ for $M_2\le1/6$ and $1-3M_2$ otherwise.
This can be used together with the original purity coefficient
in (SR7) by taking their maximum. A uniform upper bound on $M_2$
can be substituted because the two-branch lower coefficient is
continuous and nonincreasing in $M_2$.

The variance factor also improves using the weighted orthogonality
of the remainder. Let $T$ be the operator on the weighted space
of (SR8) with quadratic form $\operatorname{Var}_{\mu_1}$.
Then $0\preceq T\preceq B I$. Put $V=\langle\phi,T\phi\rangle$.
We have
$\|P_{\phi^\perp}T\phi\|^2=\langle\phi,T^2\phi\rangle-V^2
\le V(B-V)$. For $A,V>0$, project $T^{1/2}z$ onto the unit vector
$(\alpha/A)T^{1/2}\phi/\sqrt V$ in the label-valued space.
It follows that
\[
 s_2\ge A^2\big(\sqrt V-\sqrt{B-V}\,r\big)_+^2.                \tag{SR15}
\]
The bound is trivial at $V=0$. Finally, the physical mass balance
$\mathbb E_\nu z=0$ bounds the tilted mean by
$m_0^2\le\chi s_2$, where
\[
 \chi=\int (d\nu/d\mu_1-1)^2d\mu_1
      =\frac{\pi d_\rho^{3/2}}{2\sqrt{1-\kappa}}-1
                     \qquad(\kappa<1).                         \tag{SR16}
\]
Indeed $-\mathbb E_{\mu_1} z$ is the pairing of the centered
profile with $d\nu/d\mu_1-1$. Cauchy–Schwarz gives the bound,
and the integral is elementary in $t$. Combining this estimate
with $m_0^2+s_2=M_2$ allows use of the sharper constraint
$m_0^2\le\chi M_2/(1+\chi)$ when optimizing (SR14).

The cosine scalar penalty in (SR6) can replace the quartic penalty
when the profile has very little mixed mass. For the tilt $\mu$ in
(SR6) set $P=\mathbb E_\mu(1-|p|^2)$. Convexity of
$x\mapsto\cos(\pi\sqrt x/2)$ on $[0,1]$ bounds it by the chord
$1-x$. With $x=(2p_i-1)^2$ this gives
$a_i(p)\le(4/\pi)p_i(1-p_i)$. Also
$p_i(1-p_i)\le(1-|p|^2)/2$, as seen by summing the products
$p_jp_k$ over the unordered pairs of labels. Thus
$\max_i\mathbb E_\mu a_i\le (2/\pi)P$ and
$\sum_i\mathbb E_\mu a_i\le(4/\pi)P$. In (SR6), using
$\mathbb E_\mu|z|^2=2/3-P$ gives
\[
 \mathbb E_{\mu\otimes\mu}\Phi(p,q)
           \ge\frac1{\pi^2}\left(\frac5{12}+\frac{11}4P-8P^2\right).
                                                               \tag{SR17}
\]
In particular this is at least $5/(12\pi^2)$ for $P\le11/32$.
The two tilted measures obey $d\mu/d\mu_1\le(1+\rho)^{3/2}$ by
their exact gamma normalizations. A small value of the purity
factor under $\mu_1$ therefore guarantees a positive cosine loss,
while larger values can be used in the quartic estimate. The
normalization multiplying (SR17) in the radial penalty is
$\rho\pi(1-\rho)/(2(1+\rho)^2)$.

In the one-mode decomposition there is a separate purity lower
bound from the direction transverse to $\alpha$. Let $n$ be a
unit centered label vector perpendicular to $\alpha$ and write
$b_i=\mathbb E_{\mu_1}p_i$. The variance of $n_i$ when a label is
sampled with probabilities $p_i$ is at most $1-|p|^2$: write it as
$\sum_{i<j}p_i p_j(n_i-n_j)^2$ and use $(n_i-n_j)^2\le2$.
Since $n\cdot p=n\cdot g$ and $\max_i n_i^2\le2/3$, this gives
\[
 \mathbb E_{\mu_1}(1-|p|^2)
    \ge \sum_i b_i n_i^2 - BG^2
    \ge 2/3-u_\rho-BG^2,
 \qquad u_\rho=\frac{e^{-1/2}}
          {3\Gamma(3/2)d_\rho^{3/2}\sqrt{2\kappa}}.             \tag{SR18}
\]
For the last inequality, balance under $\nu$ gives $b_i\le u_\rho$
using the maximum of the density $d\mu_1/d\nu$. Therefore
$\sum_i b_i n_i^2\ge u_\rho-(2/3)\sum_i(u_\rho-b_i)=2/3-u_\rho$.
This lower bound uses only the transverse residual to charge
imperfection of the one-mode representation. It can replace the
first inequality of (SR10) whenever it is stronger.

There is also a sharper raw second-moment bound if the trial's
raw tilted moment $M_\phi=\mathbb E_{\mu_1}\phi^2$ is bounded.
Let $W$ be the multiplier $d\mu_1/(J_m\,d\nu)$ on the weighted
space, extended to label-valued functions. Thus $0\preceq W\preceq B I$.
For $r=G/A\le\sqrt{(B-M_\phi)/M_\phi}$ with $M_\phi>0$, one has
\[
          \mathbb E_{\mu_1}|z|^2
               \le A^2(\sqrt{M_\phi}+r\sqrt{B-M_\phi})^2 .      \tag{SR19}
\]
Indeed take the weighted unit vectors $v=(\alpha/A)\phi$ and
$w=g/G$ when $G>0$. If $t=\langle w,Ww\rangle$, positivity of
$W$ and $BI-W$ and orthogonality show
$|\langle v,Ww\rangle|^2\le\min(M_\phi t,(B-M_\phi)(B-t))$.
The upper bound for the quadratic expression
$M_\phi+2r|\langle v,Ww\rangle|+r^2t$ increases up to
$t=B-M_\phi$. On the remaining interval it decreases under the
stated condition on $r$, by differentiating the square-root
bound there. Substitution gives (SR19). The case $G=0$ is an
identity. The right hand side is increasing in $M_\phi$ up to
$B/(1+r^2)$, so a certified upper bound in that range can be used.

For the low-purity alternative (SR17), the exact density ratio
gives a stronger bound than its supremum. Let $F$ be the gamma
distribution function of shape $3/2$ and rate one. If
$P_1=\mathbb E_{\mu_1}(1-|p|^2)$ and $F(x_0)\ge3P_1/2$, then
\[
                     P\le\frac23 F((1+\rho)x_0) .              \tag{SR20}
\]
The function $1-|p|^2$ ranges in $[0,2/3]$. Its expectation under
$\mu$ with a specified expectation under $\mu_1$ is maximized by
placing it at $2/3$ on an initial radial interval, because
$d\mu/d\mu_1$ decreases in the radius. This follows by subtracting
the density ratio at the cutoff and checking the sign on both
sides. In the variable $(1+\kappa)t$ the measure $\mu_1$ has
distribution function $F$, and $\mu$ has effective rate $1+\rho$.
This gives (SR20). Combining exact tilt concentration with the
two possible nonlinear penalties is a finite-amplitude route beyond
the verified interval.
