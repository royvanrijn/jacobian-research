# Simplex feasibility in the finite-amplitude radial comparison

This note gives an exact amplitude obstruction from circular
mass feasibility. It strengthens the existing one-direction
finite-amplitude argument. It also gives a finite-dimensional
sufficient test when several radial directions could be adverse.
Neither statement assumes regularity of the angular slices.

Use the weighted radial notation of
[the finite-amplitude comparison](THREE_CELL_FINITE_AMPLITUDE_COMPARISON.md).
Write $p_*=\mathbf1/3$, $z=p-p_*$, $x=\sqrt t$,
and $d\nu(t)=e^{-t}dt$. In a one-direction decomposition
put $z=\alpha\phi+\epsilon$, $A=|\alpha|$,
$G^2=\int J_+|\epsilon|^2d\nu$, and
$v=\alpha/A\in H=\mathbf1^\perp$ for $A>0$.

## Exact projection onto the circular mass simplex

At every radius the residual has size at least
$\operatorname{dist}(p_*+A\phi v,\Delta_3)$.
At the pole this distance gives a sharp amplitude cost,
with a closed formula depending on the direction $v$.
Sort $v_1\ge v_2\ge v_3$, write $a=-v_3$ and
$b=v_1-v_2$. Then for $u\ge0$,
\[
 \operatorname{dist}\big(u p_*+v,u\Delta_3\big)^2
       =\tfrac16(3a-u)_+^2+\tfrac12(b-u)_+^2,
 \qquad
 \mathcal C(v):=\int_0^\infty\operatorname{dist}
            \big(u p_*+v,u\Delta_3\big)^2du
       =\tfrac32a^3+\tfrac16b^3.                            \tag{AP1}
\]
Indeed for $u\ge3a$ all coordinates are already
nonnegative. For $b\le u<3a$, the projection sets
the last coordinate to zero and subtracts
$(a-u/3)/2$ from each of the first two. Then the
second projected coordinate is $(u-b)/2$. For $u<b$
it must also be set to zero, so the solution is
$(u,0,0)$, and the squared residual gains an
additional $(b-u)^2/2$. These projections can
also be verified directly by the variational
inequality for nearest points in the simplex.

The constraints on the unit centered vector are
$3a^2+b^2=2$ and $0\le b\le3a$. The minimum of
$\mathcal C$ is $1/\sqrt6$, attained exactly at
permutations of $(2,-1,-1)/\sqrt6$. For example,
minimize $x^3/(2\sqrt3)+y^3/6$ for
$x=\sqrt3 a$, $y=b$ on $x^2+y^2=2$,
$0\le y\le\sqrt3 x$. The function decreases
until $y=\sqrt3 x$, by differentiation on
the circle.

Suppose $J_+(x^2)e^{-x^2}=j x+O(x^2)$ and
$\phi(x^2)=c/x+O(1)$ near zero with $j,c>0$.
Then the integral of the squared pointwise distance
on any sufficiently small fixed origin band is
\[
 \int_{\sqrt t\le x_0}J_+(t)
         \operatorname{dist}(p_*+A\phi(t)v,\Delta_3)^2d\nu(t)
     =2jc^3A^3\mathcal C(v)+o(A^3)                          \tag{AP2}
\]
as $A\downarrow0$, uniformly in unit $v\in H$.
Indeed after $x=Ac u$ the clipping region is
uniformly bounded in $u$ since $p_*$ is interior.
Multiply the distance by $u$ before passing
to the limit. Its square becomes the bounded
distance expression in (AP1), and the Gaussian
weight and reference errors converge uniformly
on this rescaled band. This is the sharp
coefficient for pointwise feasibility on the
origin band. The true residual can also have
additional requirements from balance and projection.

For exact adverse eigenmodes of the compressed weighted
operator (RMA2), each radial eigenprofile with nonzero
eigenvalue has the form $\beta/x+O(1)$ near zero. In
fact its equation is $\mathcal K_+z=\lambda J_+z+d_0$.
As functions of $x=\sqrt t$, both $\mathcal K_+z(x^2)$
and $J_+(x^2)$ are analytic near zero, and $J_+$ has
a simple zero with positive first derivative there.
The kernel integrals can be differentiated locally in
complex $x$: their absolute majorant is a polynomial
in $\sqrt s$ times $e^{-\kappa s+C\sqrt s}$ after
canceling $I_0$, and Cauchy–Schwarz uses the lower
weight bound (ST3). Thus the pole coefficient is a
single linear functional on any finite adverse space.
Its kernel consists of profiles bounded near the origin.
In particular origin clipping alone only controls the
pole part of a general combination of adverse modes.

## A stronger analytic version of the one-mode calculation

Under the existing bounds (FA8), the circumradius
$R_c=\sqrt{2/3}$ of $\Delta_3$ about $p_*$ gives a
simple finite-amplitude bound. If $A\ge A_0$
and $X_0=A_0c_* /(R_c+A_0 C_*)\le1/10$, then
\[
 \frac{G^2}{A^2}
      \ge \frac{2j_*c_*^3 A_0}{3(R_c+A_0 C_*)}.             \tag{AP3}
\]
Indeed the norm residual divided by $A$ is at least
$(|\phi|-R_c/A)_+$, which increases with $A$.
Use the lower pole envelope in (FA8) and integrate
$2j_*[c_*-(R_c/A_0+C_*)x]^2$ from zero to $X_0$.
This is six times the lower bound (FA11) from one
coordinate. With the retained constants and $A_0=1/200$,
the right side exceeds $.0028$, giving
$A<1/200$ in the nonnegative-scalar case (FA7).

Use the exact first cosine penalty of
[the angular envelope](THREE_CELL_EXACT_ANGULAR_MOMENT_ENVELOPE.md),
denoted $\Gamma_{\mathrm{ex}}(p,q)$.
On the existing annulus of (FA13) the reference vectors
$p^0=p_*+\alpha\phi$ now lie in the cube of (AM5)
with coordinate margin greater than $.03$. For any
actual pair on that annulus one has
\[
 \Gamma_{\mathrm{ex}}(p(t),p(s))\ge
  \frac{3}{160}A^2(\phi(t)-\phi(s))^2
        -\frac1{20}(|\epsilon(t)|^2+|\epsilon(s)|^2).        \tag{AP4}
\]
When both profiles lie in the cube, apply (AM5)
and $|u+v|^2\ge |u|^2/4-|v|^2/3$.
If one leaves the cube, its residual norm is at
least $.03$; meanwhile $A|\phi(t)-\phi(s)|$
is at most $2(97/20)/200=.0485$. The right side
is then nonpositive, proving the bound by
nonnegativity. Integrating (AP4) and using
(FA13) and (FA15) gives
\[
 \mathcal L_{\mathrm{ex}}(p)
       \ge \frac9{320}A^2-\frac{20}{27}G^2.                 \tag{AP5}
\]
This provides a larger radial deficit margin
on the certified band of (FA5) using its retained
spectral and radial inputs. The angular estimate
itself is analytic, and no mass-cube box search
is needed for (AP4). The estimates here impose
no new correlation cutoff.

## A finite-dimensional test with an optimized angular penalty

Fix a positive correlation and an adverse spectral
space $E_\rho$ of finite dimension as in (RMA4).
Use $\sqrt J z_i=\alpha_i+g_i$ with $\alpha_i\in E_\rho$
and $g_i\perp E_\rho$, and put
$\bar p_i=p_{*,i}+\alpha_i/\sqrt J$. Write
$C_V=\sup_t J(t)^{-1}\int V(t,s)d\nu(s)<\infty$
as in (FA15), and let $\eta>0$ be the negative
scalar margin on the remaining space. For $\ell>0$
define the finite-profile angular cost
\[
 B_\ell(a,b)=\min_{p,q\in\Delta_3}
   [\Gamma_{\mathrm{ex}}(p,q)+\ell(|p-a|^2+|q-b|^2)] .
\]
It is continuous and nonnegative, and the minimum
is on a compact finite-dimensional set.
For $\ell<\eta/(2C_V)$ the exact deficit comparison
implies
\[
\begin{split}
 S_\rho(A)-S_\rho(T)\ \le{}&
  \sum_i\langle\alpha_i,(B_\rho-2)\alpha_i\rangle
  -\iint V(t,s)B_\ell(\bar p(t),\bar p(s))d\nu(t)d\nu(s)\\
 &-(\eta-2\ell C_V)\int J(t)
                  \operatorname{dist}(\bar p(t),\Delta_3)^2d\nu(t).
\end{split}                                                   \tag{AP6}
\]
Indeed, applying the minimum at the actual mass profiles
costs at most $2\ell C_V\sum_i\|g_i\|^2$ after integration.
Use the remaining scalar margin for pointwise feasibility.
The row bound also ensures integrability of all terms.
The directions $\alpha$ lie in an explicitly bounded
subset of $E_\rho^3$, since $|p-p_*|^2\le2/3$ and $J$
is bounded. A sufficient comparison for this correlation
is therefore nonpositivity of the right side of (AP6)
on that finite-dimensional set with $\sum_i\alpha_i=0$.
Additional angular modes or the simultaneous compatibility
loss may still be needed to make this test pass.
