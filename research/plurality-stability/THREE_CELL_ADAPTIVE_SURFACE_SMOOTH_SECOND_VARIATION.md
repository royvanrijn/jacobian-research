# The smooth second variation of the three-cell surface energy

At every \(0<r<1\), the centered Gaussian simplex has
positive second variation of the equilateral pair surface
energy along smooth, exactly balanced deformations with
nonzero normal speed, modulo rotation. This includes
deformations that move its common junction. The same
calculation gives strict quadratic growth along these paths
for the adaptive direction supremum. It concerns smooth
paths preserving the three-arm geometry near the simplex.

The proof uses the weighted radial term of the common
junction comparison. Exact balance gives a common
first-order mean for the three normal speeds. A covariance
estimate bounds the term that couples the junction
displacement to the remaining speed on each arm.
[Weighted radial comparison (OA11)](THREE_CELL_ORIGIN_ARM_SURFACE_COMPARISON.md#quadratic-loss-at-a-common-displaced-junction)

## Smooth deformations and the radial term

Let \(T\) be the centered regular simplex in the plane.
Consider a smooth family \(A_\epsilon\) obtained by
transporting \(T\) through diffeomorphisms that equal the
identity outside a fixed compact set, with \(A_0=T\).
Assume all three Gaussian cell masses stay equal to
\(1/3\). A common rotation can be composed with the
diffeomorphisms without changing the energies. Write
\(v_i\) for the three outward unit tangents of the simplex
rays, and use normal orientations \(n_i=Jv_i\) for the
same quarter turn \(J\). Denote the common vertex by
\(p_\epsilon\), where \(p_0=0\) and
\(p=\partial_\epsilon p_\epsilon|_0\). Each arm has a
smooth arclength parametrization \(x_i^\epsilon(s)\)
from \(p_\epsilon\) for \(s\ge0\). Put
\[
 f_i(s)=\left.\partial_\epsilon x_i^\epsilon(s)\right|_0\cdot n_i,
       \qquad q_i=p\cdot n_i=f_i(0).
\]
The first variation of the cell masses is the sum of the
normal fluxes with the cyclic label-transfer signs.
The three means \(\int f_i\,d\mu\) therefore have a common
value \(b\), where \(\mu\) is the standard half-Gaussian
probability measure on \(s>0\).

Let \(u_i^\epsilon(s)=x_i^\epsilon(s)-p_\epsilon\).
At \(\epsilon=0\) the unit-speed identity
\(|\partial_s u_i^\epsilon|=1\) gives
\[
 \left.\partial_\epsilon u_i^\epsilon(s)\right|_0
                 =(f_i(s)-q_i)n_i,\qquad
 v_i\cdot\left.\partial_\epsilon^2 u_i^\epsilon(s)\right|_0
                 =-\int_0^s f_i'(t)^2\,dt .
\]
For the first equality, the tangential derivative in
\(s\) of the initial shape speed is zero and its
relative value at \(s=0\) vanishes. For the second,
differentiate the unit-speed identity twice and integrate
from \(s=0\). Consequently \(s^2-|u_i^\epsilon(s)|^2\)
has quadratic coefficient
\[
                 D_i(s)=s\int_0^s f_i'(t)^2\,dt
                                 -(f_i(s)-q_i)^2\ge0 .
\]
Using \(\sum_i v_i=0\) and
\(\sum_i q_i^2=3|p|^2/2\), the sum of radial deficits
about the Gaussian origin satisfies
\[
 \sum_i(s^2-|x_i^\epsilon(s)|^2)
   =\epsilon^2\sum_i\bigl[D_i(s)-2q_i f_i(s)\bigr]
          +o(\epsilon^2).                                  \tag{AV1}
\]
In fact the terms from the translated vertex are
\(-3|p|^2-2\sum_i q_i(f_i-q_i)\), and the linear and
acceleration terms paired with \(s\sum_i v_i\) vanish.
Taylor expansion under the Gaussian integrals below is
justified by compact support of the transport. Beyond
that support the only arclength change is a smooth
constant offset on each ray. There is thus polynomial
domination in \(s\) for the necessary derivatives.

Use the reference marginal \(w(s)\) in (OA11), and write
\(w(s)\,ds=\alpha(s)\,d\mu(s)\). Put
\[
 m=\int s\,d\mu(s),\qquad W=\int_0^\infty w(s)\,ds,\qquad
 M=\int_0^\infty s\,w(s)\,ds,\qquad C=M/m-W .
\]
The function \(\alpha\) is strictly increasing. To see
this, integrate the same-arm and cross-arm reference
kernels in the second coordinate. With
\(a=r/\sqrt{1-r^2}\), the result for \(\alpha\) is a
positive constant times
\[
     2\Phi(as)+e^{-3a^2s^2/8}\Phi(-as/2).
\]
For \(x=as>0\), the derivative of the displayed
expression with respect to \(x\) is
\(\tfrac32\phi_1(x)-\tfrac{3x}{4}e^{-3x^2/8}\Phi(-x/2)>0\)
by the Gaussian upper tail bound
\(y\Phi(-y)<\phi_1(y)\) at \(y=x/2\).
In particular
\[
 C=\frac1m\int_{s>t}(s-t)(\alpha(s)-\alpha(t))
                                  \,d\mu(s)d\mu(t)>0.       \tag{AV2}
\]

## A covariance bound for the balanced arm speeds

Subtract the affine function and put
\[
 k_i(s)=f_i(s)-q_i(1-s/m)-(b/m)s .
\]
Then \(k_i(0)=0\) and \(\int k_i\,d\mu=0\).
The arclength deficit is unchanged by this substitution:
\(D_i(s)=s\int_0^s k_i'(t)^2\,dt-k_i(s)^2\).
For a single such function \(k\), write \(g(s)=k(s)/s\)
and \(D(s)=s\int_0^s k'(t)^2\,dt-k(s)^2\).
For \(s>t>0\),
\[
   |k(s)-(s/t)k(t)|^2\le(s/t-1)D(s).
\]
Indeed \(k'-k(s)/s\) has integral zero on \([0,s]\).
Pair it with \(1_{[0,t]}-t/s\), whose squared
\(L^2([0,s])\) norm is \(t(1-t/s)\), and apply
Cauchy–Schwarz. Therefore, by symmetrizing against
the zero mean of \(k\),
\[
\begin{split}
 \int w k&=\frac1m\int_{s>t} st(\alpha(s)-\alpha(t))
                       (g(s)-g(t))\,d\mu(s)d\mu(t),\\
 \left|\int w k\right|^2
 &\le C\left[\frac1m\int_{s>t}
         t(\alpha(s)-\alpha(t))D(s)\,d\mu(s)d\mu(t)\right]\\
 &\le C\int wD .                                            \tag{AV3}
\end{split}
\]
For the first inequality split the factor
\(\sqrt{t(s-t)D(s)}\bigl(\alpha(s)-\alpha(t)\bigr)\)
by Cauchy–Schwarz and use (AV2). For the last, use
\(\alpha(s)-\alpha(t)\le\alpha(s)\) and
\(\int_0^s t\,d\mu(t)\le m\).

As \(\sum_i q_i=0\), the quadratic coefficient of the
weighted radial term in (AV1) is
\[
 \begin{split}
 R
   &=\sum_i\left[\int wD_i+2Cq_i^2-2q_i\int w k_i\right]\\
   &\ge \sum_i\left[
            \bigl(\sqrt{\int wD_i}-\sqrt C\,|q_i|\bigr)^2
                                  +Cq_i^2\right].            \tag{AV4}
 \end{split}
\]
Equality at zero requires \(q_i=0\) and \(D_i=0\) on
every arm. The equality case of Cauchy–Schwarz in the
definition of \(D_i\) then gives \(f_i(s)=a_i s\).
The common mean condition makes all \(a_i\) equal,
so this speed is a rotation. With compact support it
is zero.

## Coercivity modulo rotation

On every compact interval \(I\subset(0,1)\), the lower
bound in (AV4) controls the Gaussian \(H^1\) distance to
the rotation speeds:
\[
 R\ \ge\ c_I\inf_{a\in\mathbb R}
                  \sum_i\|f_i-a s\|_{H^1(\mu)}^2.           \tag{AV5}
\]
Here \(\|h\|_{H^1(\mu)}^2=\int(h^2+h'^2)\,d\mu\).
The formula for \(\alpha\) makes it bounded below by a
positive constant uniform on \(I\). Also \(C>0\) is
continuous in \(r\) by (AV2) and dominated convergence,
so it is bounded below there. Since
\((x-y)^2+y^2\ge(x^2+y^2)/3\), (AV4) bounds \(R\)
below by a positive multiple of
\(\sum_i(q_i^2+\int D_i\,d\mu)\), uniformly on \(I\).

Put \(g_i=f_i-q_i\) and write \(g_i=a_i s+h_i\), where
\(\int s h_i\,d\mu=0\). Both \(g_i\) and \(h_i\) vanish
at zero. If \(d\mu(s)=\rho(s)\,ds\), then
\(\int_t^\infty s\,d\mu(s)=\rho(t)\). Integrating the
definition of \(D_i\) in the opposite order gives
\[
 \int D_i\,d\mu
   =\|g_i'\|_{L^2(\mu)}^2-\|g_i\|_{L^2(\mu)}^2
   =\|h_i'\|_{L^2(\mu)}^2-\|h_i\|_{L^2(\mu)}^2 .
\]
The second equality follows by Gaussian integration by
parts, since \(\int h_i'\,d\mu=\int s h_i\,d\mu=0\)
and \(\int s^2\,d\mu=1\).
The odd extension of \(h_i\) has no degree-one Hermite
component. Thus the next possible degree is three,
giving \(\|h_i'\|_2^2\ge3\|h_i\|_2^2\) and
\(\int D_i\,d\mu\ge\tfrac12\|h_i\|_{H^1(\mu)}^2\).

Write \(\bar a=\tfrac13\sum_i a_i\). The common flux
mean \(b=q_i+m a_i+\int h_i\,d\mu\) and \(\sum_i q_i=0\)
give
\[
 m(a_i-\bar a)
    =-q_i-\int h_i\,d\mu+\tfrac13\sum_j\int h_j\,d\mu .
\]
Consequently \(\sum_i(a_i-\bar a)^2\) is bounded by a
constant times \(\sum_i(q_i^2+\|h_i\|_2^2)\). Apply
this estimate to \(f_i-\bar a s=q_i+h_i+(a_i-\bar a)s\)
to obtain (AV5).

## Consequence for the adaptive bound

Apply (OA11) for each \(\epsilon\) and retain its
radial term. The angular term is nonnegative. It follows
that
\[
 \mathcal J_r(A_\epsilon)-\mathcal J_r(T)
          \ge \frac{R}{1+r}\epsilon^2+o(\epsilon^2).        \tag{AV6}
\]
The equilateral energy is twice differentiable along
the smooth path, by differentiating the Gaussian
surface integrals in their arclength coordinates.
The two signs of (AV6) force its first derivative to
vanish. Its second derivative is at least
\(2R/(1+r)>0\) whenever the normal speed is nonzero
modulo rotation.

The adaptive squared energy \(\mathcal F_r^2\)
dominates \(\mathcal J_r\) and the two agree at \(T\).
Thus (AV6) also gives its strict quadratic lower
growth along the indicated paths, without requiring
twice differentiability of the supremum. This
argument does not give a uniform local comparison
in Gaussian label measure, which would have to
control nonsmooth changes and changes of the
interface topology.
[Adaptive energy of the simplex](THREE_CELL_ADAPTIVE_SURFACE_VARIATIONAL_PROBLEM.md#a-dual-formula-with-metric-surface-weights)
