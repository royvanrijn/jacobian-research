# Normal saturation of the adaptive Gaussian surface energy

At an actual three-cell noise-stability maximizer, using the interface
normals as direction fields attains the adaptive surface energy.
For a general planar partition this saturation has a converse
description. With three direction coordinates, it is equivalent to
normal score gradients and nonnegative first extra-coordinate
variation. It also forces complete adjacency. A separate bound shows
that any fixed finite family of dual tests loses the scale of the
simplex energy as the correlation tends to one.

## Normal saturation and the scalar interface form

Let \(A\) be a planar three-cell partition of finite Gaussian
perimeter, fix \(0<r<1\), and let \(\mathcal F_{r,3}(A)\) be the
adaptive functional with directions in \(\mathbb R^3\). Use the
notation \(e=(i,j)\), \(d_e=e_i-e_j\), \(N_e\) pointing out of
\(A_i\) and into \(A_j\), and \(\nu_e\) for the tagged Gaussian
interface measure. Embed the two-dimensional normal fields in
the first two direction coordinates. The definitions of the
adaptive functional and its pairing use the same formulas as
in (AS1) with the norm summing three direction coordinates.
[The adaptive functional and (AS1)](THREE_CELL_ADAPTIVE_SURFACE_VARIATIONAL_PROBLEM.md#a-dual-formula-with-metric-surface-weights)

Assume \(\mathcal F_{r,3}(A)<\infty\). Taking the three distinct
unit coordinate directions on the three tagged interfaces shows
that each \(T_{\sqrt r}\nu_e\) lies in \(L^2(\gamma_2)\). In
particular the scalar kernel products below are absolutely
integrable for bounded directions and speeds. Define
\[
 \begin{split}
 (Kv)_e(x)&=\sum_f(d_e\cdot d_f)
      \int M_r(x,y)v_f(y)\,d\nu_f(y),\\
 B_e(x)&=(KN)_e(x)
   =-\frac1r\nabla T_r(1_{A_i}-1_{A_j})(x).
 \end{split}                                                     \tag{NS1}
\]
The last identity follows from Gaussian distributional integration
by parts and the derivative commutation formula for the Mehler
operator. In detail, the Gaussian distributional gradient of
\(1_{A_i}\) is \(-\sum_e(d_e)_i N_e\nu_e\). Its smoothing by \(T_t\),
multiplied by \(t\), equals \(\nabla T_t1_{A_i}\). For fixed
observations the kernel and its source derivatives have Gaussian
decay at infinity, so smooth cutoff approximations justify this
identity for finite Gaussian perimeter. The gradient on the right
in (NS1) is uniformly bounded, since the conditional Gaussian
derivative has an integrable first moment independent of \(x\).

Write \(m^v=\sum_e d_e\otimes v_e\,\nu_e\) and use \(\langle\, ,\,\rangle_\nu\)
for the summed interface pairing. The semigroup identity gives
\(\|T_{\sqrt r}m^v\|_2^2=\langle v,Kv\rangle_\nu\). Moreover,
writing \(S_r(A)=\sum_i\langle 1_{A_i},T_r1_{A_i}\rangle\), one has
\[
 S'_r(A)=r^{-1}\sum_i\|\nabla T_{\sqrt r}1_{A_i}\|_2^2
        =\|T_{\sqrt r}m^N\|_2^2.
\]
For the first equality, the degree \(n\) Hermite summand on each
side is \(nr^{n-1}\) times the squared norm of that component.
The second follows from the same derivative commutation formula.

The equality \(\mathcal F_{r,3}(A)^2=S'_r(A)\) holds if and only
if there exist bounded nonnegative \(\lambda_e\) such that
\[
 B_e=\lambda_e N_e\quad \nu_e\text{-almost everywhere},\qquad
 \langle f,Kf\rangle_\nu\le\sum_e\int\lambda_e f_e^2\,d\nu_e
       \quad\text{for all bounded real speeds } f.                  \tag{NS2}
\]
When these conditions hold, \(\lambda_e>0\) almost everywhere
on each nonzero tagged interface.

For the forward implication, \(N+t(v-N)\) is admissible for
\(0\le t\le1\) and every field \(v\) with \(|v|\le1\). Its
one-sided derivative at zero gives
\(\langle B,v-N\rangle_\nu\le0\). Pointwise maximization of the
left pairing over unit balls gives
\(B_e\cdot N_e=|B_e|\) almost everywhere. This proves the first
condition with \(\lambda_e=|B_e|\). For any bounded real speed
take the unit field
\[
 u_e(t)=\sqrt{1-t^2 f_e^2}\,N_e+t f_e\,e_3,
\]
where \(e_3\) denotes the direction orthogonal to the active
plane. The coefficient of \(t^2\) in its energy is
\(\langle f,Kf\rangle_\nu-\sum_e\int\lambda_e f_e^2\,d\nu_e\).
The expansion follows from the uniform Taylor expansion for
bounded \(f\) and absolute product integrability. Maximality
of \(N\) makes this coefficient nonpositive. If \(\lambda_e\)
vanished on a set of positive \(\nu_e\)-measure, a speed one
on that set and zero on all other interfaces would give
strictly positive interaction by \(d_e\cdot d_e=2\) and strict
positivity of the Mehler kernel. This proves positivity
almost everywhere. For the converse, apply the scalar bound
in (NS2) to each coordinate of an admissible field. The
energy is at most \(\sum_e\int\lambda_e\,d\nu_e\), with equality
for \(N\) by \(KN=\lambda N\).

On every connected smooth open arc where the pair interface
has positive surface density, (NS1) and (NS2) imply that
\(T_r(1_{A_i}-1_{A_j})\) is constant along the arc. The
continuous tangential derivative vanishes almost everywhere
and hence everywhere there. This conclusion allows different
constants on different components of the same pair. Actual
maximizers satisfy (NS2) by the translation equations and
the first extra-coordinate form. The converse here concerns
the surface energy.
[The translation equations and extra-coordinate form](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

## Complete adjacency at saturation

If all three cells have positive mass and the normal fields
saturate the surface energy, each pair interface has positive
Gaussian measure. Indeed, if the adjacency graph were
disconnected, some cell would have Gaussian distributional
gradient zero, which would make it constant almost everywhere.
If only two edges are present, orient both out of their common
cell. Reversing an edge reverses both \(d_e\) and \(N_e\) and
does not change the functional or (NS2). For the two edges,
all four products \(d_e\cdot d_f\) are strictly positive.

For each \(v\in\mathbb R^2\), the scalar speed \(h_e=N_e\cdot v\)
satisfies \(Kh=\lambda h\). If \(h\) took both signs on sets of
positive tagged measure, replacing \(h\) by \(|h|\) would keep
the local side of (NS2) and strictly increase the interaction.
This contradicts (NS2). A countable dense set of \(v\) shows
that the oriented normals on the two edges agree with one
unit vector almost everywhere: any two distinct normal
directions in the support of their distribution can be strictly
separated in sign by one of these \(v\).

The distributional gradients of all the indicators are
therefore parallel to that vector. Each cell depends only on
the corresponding scalar coordinate, by the vanishing
tangential derivative. The common cell has a distributional
derivative of one sign in that coordinate, so its indicator
is monotone and it is a half-space up to null sets. Inside the
open complementary half-space the other two indicators have
zero gradients, since no mutual edge exists and their only
interfaces with the common cell are on its boundary. Each is
constant there. Both cannot have positive mass. This proves
the claim.

## A uniform bound for fixed dual tests

Let \(G=(G_1,G_2,G_3)\) have square-integrable vector-valued
rows with summed \(L^2(\gamma_2)\) norm at most one, and write
\(V_i=T_{\sqrt r}G_i\). For an ordinary partition write
\[
 W_G(A)=\sum_{i<j}\int |V_i-V_j|\,d\nu_{ij}.
\]
For every \(r\) and every such test, there exist balanced
ordinary partitions with \(W_G(A)\le 2\). More generally, for
\(m\) such tests \(G^{(1)},\ldots,G^{(m)}\) one has the bound
\[
       \inf_{A\ {\rm balanced}}\max_{1\le j\le m} W_{G^{(j)}}(A)
          \le 2\sqrt m .                                            \tag{NS3}
\]
The estimate is independent of the number of vector
coordinates.

To prove it, let \(z(x)=1-\exp(-|x|^2/2)\), which is uniform on
\([0,1]\) under \(\gamma_2\). For \(0<p<1/3\), label the radial
shells by the three consecutive intervals of length \(1/3\)
on the \(z\)-circle starting at \(p\). This gives a balanced
partition \(A(p)\). Its interfaces are the circles with \(z\)
equal to \(p,p+1/3,p+2/3\), with a different tagged pair at
each cut. On each circle, use
\[
 \max_j |V^{(j)}_i-V^{(j)}_k|
       \le \sqrt2\left(\sum_{j,\ell}|V^{(j)}_\ell|^2\right)^{1/2}.
\]
Averaging the sum of these three surface integrals over
\(p\) with density three puts one integral on each radius.
With \(R=|x|\), the radial Gaussian weight
\(\omega(R)=R\exp(-R^2/2)\) satisfies \(dz=\omega(R)dR\).
The coarea formula for \(R\) then bounds the average of
\(\max_j W_{G^{(j)}}(A(p))\) by
\[
 3\sqrt2\int \omega(|x|)
        \left(\sum_{j,\ell}|V^{(j)}_\ell(x)|^2\right)^{1/2}
       d\gamma_2(x)
 \le 3\sqrt2\left(\int_0^\infty R^3e^{-3R^2/2}\,dR\right)^{1/2}
       \sqrt m
 =2\sqrt m.
\]
The operator \(T_{\sqrt r}\) is an \(L^2\) contraction. For a
general \(L^2\) test its smoothed values are well defined
pointwise by Cauchy–Schwarz against the Mehler kernel.
The coarea estimate in particular supplies finite surface
integrals for almost every \(p\). This proves (NS3).

At the simplex \(\mathcal F_{r,3}(T)=\sqrt{S'_r(T)}\) by the
simplex extra-coordinate gap and the normal identity. Its value
diverges on the scale \((1-r^2)^{-1/4}\). Even a finite family
of tests chosen separately at each correlation, with its
cardinality uniformly bounded, cannot certify a sharp
adaptive minimum near one. The estimate concerns families
fixed before the competing partition is chosen. The adaptive
supremum allows the test to depend on the partition.
[The simplex correlation derivative](THREE_CELL_DEFICIT_AND_PERIMETER.md#exact-correlation-derivative-for-the-three-cell-simplex)
[The simplex extra-coordinate gap](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem)
