# Paired analytic halves give a uniform factorial transfer

Consider
\[
 B_j(z;S)=\sum_{k=0}^j\binom jk e^{-S k(j-k)/(2j)}z^k
\]
on a compact positive parameter interval \(I\) of \(S\). Suppose, for
some \(R>1\), that \(B_j=I_j(z)+z^jI_j(1/z)\) where \(I_j\)
has nonnegative coefficients. Suppose that for all sufficiently large
\(j\) there is a conjugation-symmetric logarithm
\[
          H_j=\log I_j=jF+C+u_j\qquad (|z|\le R).
\]
Here \(F,C\) can depend on \(S\) and are common to all \(j\);
\(\|F\|,\|C\|\) are bounded and \(\|u_j\|\to0\) uniformly in
\(S\in I\). Assume \(\|zF'(z)\|\le s<1/2\) on this disk.

**Theorem.** There are a finite symbol radius \(T>1\), a positive
tolerance \(\omega_0\), and a terminal lower bound \(n_0\),
uniform in \(S\in I\), with the following property. For \(n\ge n_0\),
\(N=n+1\), multiply the Gaussian coefficient at index \(k\) of order \(j\)
by \(e^{E_j(k/j)}\), where the functions \(E_j\) are holomorphic
on \(|v|\le T\), have conjugation symmetry and vanish at \(0,1\).
The two resulting polynomials have disjoint root arguments if
\[
 \|E_j''\|\le j\omega,\qquad
 \|E_N''-E_n''\|\le\omega,\qquad 0\le\omega\le\omega_0 .
\]
Consequently there is a positive order-independent tolerance on
\(S\eta\) for exact ordinary factorial profiles with \(S\in I\)
and \(n\ge n_0\). If the bounds on \(F,C,u_j\) have explicit
uniform constants and a computable rate of convergence for \(u_j\),
then the tolerance and order bound can be chosen effectively.

This abstracts the [curvature-eight](MULTIVARIATE_DIFFUSE_S8_CENTERED_ANGLE_LIFT.md)
and [curvature-six](MULTIVARIATE_DIFFUSE_S6_CENTERED_ANGLE_LIFT.md)
transfers. The absolute size of \(F\) does not enter the
radial confinement or the final phase margin.

## Uniform centered multiplier

Fix \(1<R_1<R_2<R\) and a positive \(h_0\) with \(R_2e^{h_0}<R\).
Write \(D=z\,d/dz\) and \(\sigma=DF\). The bound on \(\sigma\)
gives a uniform bound \(M\) on \(D\sigma\) in \(|z|\le R_2e^{h_0}\)
by Cauchy's estimate. For either \(\Phi_j(v)=E_j(v)\) or
\(E_j(1-v)\), choose the fixed point
\[
                   w_j=\sigma(z e^{\Phi_j'(w_j)/j})
                      \qquad (|z|\le R_1).
\]
Work in \(|w|\le(1+s)/2\). The endpoint zeros give
\(|\Phi_j'(w)|/j\le2\omega\); with \(\omega\) small the
dilation remains in \(R_2\) and the fixed-point map has
contraction constant at most \(M\omega\).
For \(Z=z e^{h_j}\), \(h_j=\Phi_j'(w_j)/j\), set
\({\cal R}_j(v)=\Phi_j(v)-\Phi_j(w_j)-\Phi_j'(w_j)(v-w_j)\).
The actual half is then exactly
\[
 \widetilde I_j(z)=
 e^{\Phi_j(w_j)-w_j\Phi_j'(w_j)}
                  [e^{{\cal R}_j(D/j)} I_j](Z).
\]

For \(w=\sigma(Z)\) define
\(D_{j,m}=[(D/j-w)^m I_j](Z)/I_j(Z)\).
There are uniform constants \(K,L,P\) with
\[
       |D_{j,m}|\le K (L\sqrt{m/j})^m\quad(m\le j),
       \qquad |D_{j,m}|\le K P^m\quad(m>j/4).
\]
The first estimate is Cauchy's bound applied to
\(e^{-jwh}I_j(Ze^h)/I_j(Z)\) at radius \(h_0\sqrt{m/j}\).
The exponent \(j(F(Ze^h)-F(Z)-wh)\) has zero first derivative
and second derivative bounded by \(jM\).
One can take \(L=e^{Mh_0^2/2}/h_0\).
For the second bound, coefficient positivity yields
\(|D_{j,m}|\le (1+s)^m I_j(|Z|)/|I_j(Z)|\);
\(j<4m\) gives \(P=(1+s)e^{8\|F\|}\).
Both bounds permit \(K=\exp(2\|C\|+2\sup_j\|u_j\|)\), using
the respective uniform norms.

Choose \(A\ge2\max\{L,P,1\}\) and \(T=A+2\).
On \(|v-w_j|\le A\), \(|{\cal R}_j(v)|\le B|v-w_j|^2\)
with \(B=j\omega/2\).
For the series \(e^{{\cal R}_j}-1=\sum_{m\ge2}c_m(v-w_j)^m\):
 
- If \(m\le j/4\) and \(m\le A^2j\omega\), the Cauchy radius
  \(\sqrt{m/(j\omega)}\) gives
  \(|c_mD_{j,m}|\le K(L\sqrt{e\omega})^m\).
- If \(m\le j/4\) and \(m>A^2j\omega\), use radius \(A\)
  and \(|e^v-1|\le|v|e^{|v|}\), with \(A^2B<m/2\).
  Then \(|c_mD_{j,m}|\le(K e L^2/2)m\omega
  (e^{1/2}/4)^{m-2}\).
- If \(m>j/4\), the same circle gives
  \(|c_mD_{j,m}|\le (K A^2/2) j\omega(3/4)^m\)
  once \(e^{2A^2\omega}\le3/2\).

The first two sums are \(O(\omega)\). The third uses \(j<4m\)
and is also \(O(\omega)\). Here and below these constants are
uniform over the given disk and parameter interval, with
effectively calculable bounds from this construction.
After taking its logarithm and moving \(C+u_j\) back from
\(Z\) to \(z\), the actual half logarithm on \(R_1\) is
\(H_j+W_j+e_j\), with \(e_j=O(\omega)\) and
\[
 W_j=j(F(z e^{h_j})-F(z))+\Phi_j(w_j)-w_j\Phi_j'(w_j).
\]
It satisfies \(W_j=O(j\omega)\) and \(W_N-W_n=O(\omega)\).
Indeed at fixed \(w\) the shifts \(\Phi_j'(w)/j\) change by
\(O(\omega/N)\). The uniform contraction then bounds both
\(w_N-w_n\) and \(h_N-h_n\) on the same scale. The remaining
term \(\Phi_j(w)-w\Phi_j'(w)\) vanishes at zero and has derivative
\(-w\Phi_j''(w)\), giving its asserted cross-order bound.

## Phase margin from \(s<1/2\)

Choose \(q_0>0\) with \(e^{q_0}<R_1\). Let \(J_j(q,\theta)\)
be the difference of the two actual half logs at
\(e^{q+i\theta}\) and \(e^{-q-i\theta}\). Use \(J_j^0\) for
the same expression with the Gaussian half logs. Cauchy bounds
on the interior disk give
\[
 \|\partial_q J_j\|\le 2sj+O(1)+O(j\omega),\qquad
 \|\partial_q\partial_\theta J_j\|=O(j).
\]
Since \(2s<1\), the first bound is at most \((1-\delta)j\)
for a fixed \(\delta>0\) after increasing \(n_0\) and decreasing
\(\omega_0\). At \(q=0\) the real Gaussian difference is zero
by conjugation, so \(|\Re J_j(0,\theta)|=O(j\omega)\).
Rouché at \(q=q_0\), then on the reciprocal, uses
\(\Re J_j(q_0,\theta)<jq_0\). It confines all roots to
\(|q|<q_0\) for sufficiently small \(\omega\).
The real root equation then gives \(|q_j|=O(\omega)\).

At two hypothetical roots sharing angle \(\theta\), compare the
real equations at \(q_N,q_n\). At \(q_n\), the real baseline order
difference is \(O(|q_n|)\), by Cauchy's estimate on \(H_N-H_n\)
and its conjugation symmetry at \(q=0\). The \(W\) and \(e\)
order differences cost \(O(\omega)\).
Absorbing \((1-\delta)N|q_N-q_n|\) from moving the upper radius
gives \(N|q_N-q_n|=O(\omega)\).
In the imaginary equation, the same radius change costs
\(O(\omega\theta)\) by the mixed derivative bound.
Resetting the baseline radius to zero and comparing the two
half corrections also cost \(O(\omega\theta)\), by the
second- and first-derivative Cauchy bounds, respectively.
All imaginary errors vanish at \(\theta=0\).
The baseline imaginary comparison is
\(2\Im(F+u_N-u_n)(e^{i\theta})\), bounded by
\((2s+o(1))\theta\) uniformly. Shrinking \(\omega_0\) and
increasing \(n_0\) gives \(|\Im(J_N-J_n)|<\theta\).
Thus the difference of the two root phases would force an
integer multiple of \(2\pi\) strictly between zero and \(2\pi\)
for \(0<\theta\le\pi\), proving separation.

For ordinary coefficients use the
[exact Gamma interpolation](MULTIVARIATE_DIFFUSE_S8_CENTERED_ANGLE_LIFT.md#exact-factorial-interpolation).
On \(|v|\le T\), if \((T+1)\eta\le1/2\), then
\(|x_i-L_i|\le(T+1)|d_i|\) and \(\Re x_i\ge L_i/2\).
The same kernel estimates give
\(\|E_j''\|\le2(T+3)jS\eta\) and
\(\|E_N''-E_n''\|\le2(T+3)S\eta\).
Choosing the tolerance on \(S\eta\) small enough uniformly on
the compact positive \(I\) proves the ordinary statement.

## Explicit application down to curvature \(16/3\)

The
[Gaussian half estimate at \(16/3\)](MULTIVARIATE_DIFFUSE_EXPLICIT_S16_OVER_3_INTERLACING.md)
satisfies the theorem's assumptions on
\(16/3\le S\le6\), with \(R=1001/1000\) and \(s=3/20\).
In this interval the exact ordinary factorial moments have disjoint
arguments for
\[
                         n\ge10^7,\qquad S\eta\le5\cdot10^{-24}.
\]
Here are explicit choices and bounds for this application.

The relative half error is less than \(1/50\), with adjacent
difference less than \(10^{-6}\). Thus \(\|u_j\|\le1/49\)
and \(\|u_N-u_n\|<2\cdot10^{-6}\).
The amplitude estimate \(|2aA_*-1|<91/100\) in the saddle proof
also gives \(\|C\|\le\tfrac12\log(100/9)<5/4\).
Use \(R_1=5004/5000\), \(R_2=100081/100000\),
\(h_0=1/20000\), \(A=65536\), and \(\omega\le10^{-18}\).
The log shifts and the Cauchy contours remain inside \(R\).
Cauchy's bound gives \(M\le1600\) on \(R_2e^{h_0}\).
In the centered derivative estimates take \(K=32,\ L=20010,\ P=7\).
The three multiplier sums are bounded respectively by
\(36\cdot10^9\omega,\ 62\cdot10^9\omega,\ 12000\omega\).
For the first bound use \(L\sqrt e<33400\). For the second, the
geometric ratio \(L\sqrt e/(2A)\) is below \(13/50\).
For the third, \(P e^{2A^2\omega}/A<1/7000\).
The logarithmic remainder including the shift of \(C+u_j\)
is bounded by \(\varepsilon=2\cdot10^{11}\omega\), and the leading
terms satisfy \(\|W_j\|\le j\omega,\ \|W_N-W_n\|\le2\omega\).
The latter estimates use \(|w_j|\le3/20\) in the same
fixed-point comparison.

Use the radial test at \(q_0=1/5000\). The argument moduli are
below \(10003/10000\), where Cauchy's estimates give
\(|DH_j|\le j/6,\ |D^2 H_j|\le216j\). The corrections cost
\(2001M_0\) and \((2001+2\cdot2001^2)M_0\), respectively,
with \(M_0=j\omega+\varepsilon\).
Thus \(|\partial_q J_j|\le j/2\),
\(|\partial_q\partial_\theta J_j|\le440j\).
At \(q_0\), the baseline real half difference is at most \(jq_0/3\).
The comparison with the corrected halves confines all roots
to \(|q|<q_0\), and the individual modulus equation gives
\(|q_j|\le4(\omega+\varepsilon/j)\).

For \(H_N-H_n=F+u_N-u_n\) on the same root annulus,
the baseline difference has \(q\)-derivative bounded by one and
mixed derivative bounded by \(500\).
Hence the two modulus equations at roots with common argument give
\(N|q_N-q_n|\le24\omega+9\varepsilon\).
The imaginary comparison to the baseline at \(q=0\) costs at most
\[
 \bigl[440(24\omega+9\varepsilon)+2000(\omega+\varepsilon/n)
                         +8004\omega+8004\varepsilon\bigr]\theta
 <\theta/100 .
\]
The baseline itself is bounded by
\(2(3/20+1000\cdot2\cdot10^{-6})\theta\), leaving the strict
phase margin. Finally \(T=65538\) gives
\(\omega=2(T+3) S\eta\le10^{-18}\) under the stated ordinary
tolerance and verifies the exact factorial interpolation.


## The saddle branch loses its fixed analytic radius at criticality

For the Gaussian saddle with \(a=S/2>2\), the inverse saddle map
is \(z=q\exp(a-2aq/(1+q))\). The derivative vanishes at
\(q_c=a-1-\sqrt{a(a-2)}=e^{-r}\), with \(r=\operatorname{arcosh}(a-1)\).
It is positive on the segment \(0\le q<q_c\), and its first zero
at \(q_c\) is simple. The small branch from \(q=0\) thus has a
square-root singularity at \(z_c\), where
\[
          z_c=\exp(\sinh r-r),\qquad
          z_c-1\sim \frac{\sqrt{2}}{3}(a-2)^{3/2}
                                \quad(a\downarrow2).
\]
The formula uses \(a\tanh(r/2)=\sinh r\).
Accordingly a disk of radius \(R>1\) for this same analytic
branch can only be uniform at positive distance from critical
coupling. Separately, the
[half-log obstruction](MULTIVARIATE_GAUSSIAN_HALF_LOG_OBSTRUCTION.md)
locates actual canonical half zeros around the critical regime.
