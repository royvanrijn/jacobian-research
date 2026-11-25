# Degree-independent factorial angular separation down to curvature eight

Let \(a,b\) be distinct ordinary exponent columns. Discard equal rows
and put \(d_i=b_i-a_i\), \(L_i=\max(a_i,b_i)\),
\(S=\sum_i d_i^2/L_i\), and \(\eta=\max_i |d_i|/L_i\).

**Theorem.** The two exact ordinary factorial moment polynomials of
orders \(n,N=n+1\) have disjoint root arguments when
\[
         n\ge 600,\qquad 8\le S\le13,\qquad
                         S\eta\le 5\cdot10^{-12}.              \tag{1}
\]
The conclusion holds under independent positive rescalings. Rational
Bézout then gives coprimality over every characteristic-zero field.

The estimate uses the [Gaussian saddle at curvature eight](MULTIVARIATE_DIFFUSE_EXPLICIT_S8_INTERLACING.md).
Its small complex disk still admits an exact centered multiplier
estimate. In the imaginary part of the consecutive-order equation,
the radial displacement is bounded using the real part of the
same equation.

## Analytic hypotheses

Take \(8\le S\le13\). Write the symmetric Gaussian reference as
\(I_j(z)+z^j I_j(1/z)\), where \(I_j\) has nonnegative coefficients.
On \(|z|\le6/5\) the saddle theorem gives
\[
 I_j(z)=\exp(H_j(z)),\qquad
 H_j=jF+C+u_j,\quad
 \|F\|\le1/20,\quad\|zF'\|\le1/26,\quad
 \|C\|\le1/3,\quad\|u_j\|\le1/11,\quad
                   \|u_N-u_n\|<7/200.                        \tag{2}
\]
for \(j\ge600\). The stationary value \(F\) and the logarithm \(C\)
of the Gaussian amplitude are common to the two orders. The
last inequality uses the paired saddle error, which differentiates
the central integral after the Gaussian rescaling. All functions
in (2) have conjugation symmetry.

More generally than (1), multiply the Gaussian coefficient at
index \(k\) of order \(j\) by \(e^{E_j(k/j)}\), with \(E_j\)
holomorphic on \(|v|\le128\), real-symmetric, and zero at \(v=0,1\).
Assume
\[
 \|E_j''\|\le j\omega,\qquad
 \|E_N''-E_n''\|\le\omega,\qquad \omega\le10^{-9}.              \tag{3}
\]
We prove separation for these arrays. In each half set
\(\Phi_j(v)=E_j(v)\) or \(E_j(1-v)\). The derivative norms in (3)
are valid for \(\Phi_j\) on \(|v|\le127\). In either case
\(\Phi_j(0)=\Phi_j(1)=0\), hence
\[
 |\Phi_j'(w)|\le j\omega,\qquad
 |\Phi_N'(w)-\Phi_n'(w)|\le\omega\quad (|w|\le1/20).
\]
Indeed the corresponding derivative at zero is bounded by
one half the second-derivative norm by integrating to the other
endpoint. Integration from zero to \(w\) costs at most another
\(1/20\) of the same norm.

## Centered multiplier on the smaller saddle disk

Let \(D=z\,d/dz\) and \(\sigma=DF\). For \(|z|\le23/20\) solve
\[
 w_j=\sigma(z e^{h_j}),\qquad h_j=\Phi_j'(w_j)/j .
\]
On \(|w|\le1/20\) the right side is bounded by \(1/26\),
with \(|z e^{\Phi_j'(w)/j}|\le29/25\). For \(|Z|\le59/50\),
Cauchy's estimate on radius \(6/5\) gives
\(|D\sigma(Z)|\le (59/50)/(26(6/5-59/50))<3\).
Thus the fixed point is holomorphic and unique by iteration of
a contraction with constant at most \(3\omega\).
Set \(Z_j=z e^{h_j}\) and
\[
 {\cal R}_j(v)=\Phi_j(v)-\Phi_j(w_j)-\Phi_j'(w_j)(v-w_j).
\]
The exact coefficient multiplier for the half satisfies
\[
 \widetilde I_j(z)=
 e^{\Phi_j(w_j)-w_j\Phi_j'(w_j)}
             [\, e^{{\cal R}_j(D/j)} I_j\,](Z_j).              \tag{4}
\]
Here \(w_j,h_j\) are held fixed for the diagonal action in \(D\).
On the coefficient of \(z^k\), the exponent in (4), including
the dilation, is exactly \(\Phi_j(k/j)\).

For \(|Z|\le29/25\), \(w=\sigma(Z)\), and \(m\ge1\), put
\[
 D_{j,m}(Z)=\frac{[(D/j-w)^m I_j](Z)}{I_j(Z)} .
\]
For \(m\le j\) apply Cauchy's formula in \(h\) to
\(e^{-jwh}I_j(Ze^h)/I_j(Z)\) on
\(|h|=\sqrt{m/j}/64\).
The analytic points satisfy \(|Ze^h|\le59/50\).
The first derivative of \(F(Ze^h)-F(Z)-wh\) is zero at \(h=0\),
and the second derivative is bounded by three.
The amplitude quotient is at most
\(\exp(2(1/3+1/11))<3\). With \(m!\le m^m\) this gives
\[
               |D_{j,m}(Z)|\le 3(65\sqrt{m/j})^m .
\]
For \(m>j/4\) use the nonnegative coefficients and (2) instead.
The ratio \(I_j(|Z|)/|I_j(Z)|\) is at most \(3e^{j/10}\), giving
\[
               |D_{j,m}(Z)|\le 3(8/5)^m
\]
since \((27/26)e^{2/5}<8/5\).

Write \(B=j\omega/2\). On \(|v-w_j|\le100\) one has
\(|{\cal R}_j(v)|\le B|v-w_j|^2\). Let \(c_m\) be the Taylor
coefficients of \(e^{{\cal R}_j(v)}-1\) in \(v-w_j\), starting
at \(m=2\). The case \(\omega=0\) is exact with no error.
For \(\omega>0\) the relative error in (4) is bounded by summing
\(|c_m D_{j,m}(Z_j)|\) in three ranges.

- If \(2\le m\le j/4\) and \(m\le10000 j\omega\), the Cauchy
  radius \(\sqrt{m/(j\omega)}\) gives
  \(|c_mD_{j,m}|\le3(108\sqrt\omega)^m\). The sum is below
  \(36000\omega\).
- If \(2\le m\le j/4\) and \(m>10000j\omega\), use radius \(100\)
  and \(|e^u-1|\le |u|e^{|u|}\).
  Since \(10000B<m/2\), \(\sqrt{m/j}\le1/2\) and
  \(e^{1/2}<5/3\), the sum is bounded by
  \[
     3\cdot 5000(13/20)^2\cdot3\,\omega
             \sum_{m\ge2}m(13/24)^{m-2}<133000\omega .
  \]
- If \(m>j/4\) and \(m\ge2\), the same radius and the coefficient
  positivity bound give \(15000j\omega(1/50)^m\).
  Here \(j<4m\) and \(e^{5000j\omega}\le(5/4)^m\).
  Summing gives less than \(50\omega\).

Thus the total relative error is less than \(170000\omega\).
Its logarithm is bounded by \(171000\omega\).
Moving \(C+u_j\) from \(Z_j\) back to \(z\) costs less than
\(15\omega\) by Cauchy's bound on its first derivative along the
dilation. On \(|z|\le23/20\) the perturbed half therefore has
the analytic log
\[
 \log\widetilde I_j = H_j+W_j+e_j,\qquad
 W_j=j(F(Z_j)-F(z))+\Phi_j(w_j)-w_j\Phi_j'(w_j),\qquad
 \|e_j\|\le\varepsilon:=180000\omega.                         \tag{5}
\]
The logs have the same conjugation symmetry as the coefficients.

For each half the leading corrections satisfy
\[
                  \|W_j\|\le j\omega,\qquad
                  \|W_N-W_n\|\le2\omega .                    \tag{6}
\]
To see this, the logarithmic dilations at fixed \(w\) differ by
at most \(2\omega/N\). Using \(|D\sigma|<3\) and the contraction
gives \(|w_N-w_n|\le7\omega/N\) and \(|h_N-h_n|\le3\omega/N\).
The contribution from \(F\) in \(W_N-W_n\) is at most
\((\omega+3\omega)/26\). The remaining expression
\(K_j(w)=\Phi_j(w)-w\Phi_j'(w)\) obeys
\(K_j'(w)=-w\Phi_j''(w)\) and \(K_j(0)=0\).
At fixed \(w\) its absolute value is at most \(j\omega/800\)
and its order difference at most \(\omega/800\). Shifting \(w\)
costs at most \(7\omega^2/20\) in that difference.
These inequalities also prove the first bound in (6).

## Comparison at the actual root radii

By (2), (5) and (6) the half logs are bounded by \(j/18\)
strictly on their common disk. On \(|z|=e^{1/9}<9/8\),
\(z^j\widetilde I_j^-(1/z)\) dominates
\(\widetilde I_j^+(z)\) because the difference in the logarithmic
moduli is positive. The first polynomial has all its roots
inside the circle, since the half has no zeros in the disk.
Rouché, then the same comparison of the reciprocal polynomial,
puts every root in \(|\log|z||<1/9\).

For real \(q,\theta\) let \(J_j(q,\theta)\) be the plus half log
at \(e^{q+i\theta}\) minus the minus half log at
\(e^{-q-i\theta}\). Use \(J_j^0\) for the Gaussian halves with
both logs equal to \(H_j\). On \(|q|\le1/9\) the log arguments
have modulus at most \(9/8\). For a root \(e^{q_j+i\theta}\)
with \(0<\theta\le\pi\),
\[
 j(q_j+i\theta)=J_j(q_j,\theta)+(2k_j+1)i\pi .
                                                               \tag{7}
\]
Positive real coefficients exclude \(\theta=0\).
At either log argument, (2) and Cauchy's bounds give
\(|D H_j|\le j/20\) and \(|D^2H_j|\le19j/20\).
For example the corresponding \(F\) bounds are \(1/26\) and
\(15/26\), and the amplitude term uses \(\|C+u_j\|<17/40\).
For the corrections in (5), the Cauchy radius difference is
at least \(1/40\), giving the bounds \(45M\) and \(4095M\)
with \(M=j\omega+\varepsilon\).
Using \(j\ge600\) now gives
\[
        |\partial_q J_j|\le j/9,\qquad
        |\partial_q\partial_\theta J_j|\le 2j .
\]
Since \(\Re J_j^0(0,\theta)=0\), equation (7) implies
\(|q_j|\le3(\omega+\varepsilon/j)\).

If two roots share the argument \(\theta\), the real parts of
(7) control their different radii. At a common \(q\),
\(\Re(J_N^0-J_n^0)\) vanishes at \(q=0\) and its \(q\)-derivative
has modulus below three. This follows from
\(\|H_N-H_n\|<17/200\) on \(|z|\le6/5\).
The difference of the real equations, first at \(q_n\) and
then at \(q_N\), yields
\[
 \begin{split}
 N|q_N-q_n|
 &\le \tfrac{9}{8}\bigl(4|q_n|+4\omega+4\varepsilon\bigr)\\
 &\le20\omega+8\varepsilon .
 \end{split}
\]
In the imaginary part, changing \(q_N\) to \(q_n\) costs at
most \((40\omega+16\varepsilon)\theta\).
For the baseline order difference, changing \(q_n\) to zero
costs at most \(100|q_n|\theta\), by Cauchy's second-derivative
estimate and the bound on \(H_N-H_n\).
The variation in the two leading corrections costs at most
\(180\omega\theta\), and that of the errors at most
\(180\varepsilon\theta\). All these estimates integrate in
\(\theta\) from zero, where the imaginary parts vanish by
conjugation symmetry. The total error is at most
\[
   (520\omega+(196+300/n)\varepsilon)\theta <\theta/25 .
\]
At \(q=0\) the baseline difference is \(2\Im(H_N-H_n)(e^{i\theta})\).
Using Cauchy's radius difference \(1/5\) there bounds its
modulus by \(17\theta/20\).
Hence \(|\Im(J_N(q_N,\theta)-J_n(q_n,\theta))|<89\theta/100\).
Subtracting (7) places a multiple of \(2\pi\) strictly between
zero and \(189\theta/100<2\pi\), yielding separation.

## Exact factorial interpolation

Under (1) all unequal rows have both exponents positive.
Write \(x_i(v)=a_i+d_iv\) and use the analytic normalized Gamma
chord \(A_j(v)\) of
[matched log transport](MULTIVARIATE_HOLOMORPHIC_MATCHED_LOG_TRANSPORT.md).
For the Gaussian reference set
\[
             E_j(v)=j\bigl[A_j(v)+\tfrac S2 v(1-v)\bigr].
\]
Its sampled exponent \(E_j(k/j)\) gives exactly the ordinary
endpoint-normalized factors relative to the Gaussian.
On \(|v|\le128\), \(|x_i-L_i|\le129|d_i|\) and
\(\Re x_i\ge (99/100)L_i\). Each \(1/L_i\le\eta\) by integrality.
With \(h(u)=u/(e^u-1)\), the curvature kernel gives
\[
 A_j''(v)=\sum_i d_i^2\int_0^\infty e^{-x_i(v)u}h(u/j)\,du .
\]
The bound \(0\le1-h(u/j)\le u/(2j)\) yields
\[
 \|A_j''-S\|
 \le \left(\frac{129}{99/100}+\frac{1}{2j(99/100)^2}\right)S\eta
 <131S\eta .
\]
Also \(0\le h(u/N)-h(u/n)\le u/(2nN)\), so
\(N\|A_N''-A_n''\|\le S\eta/(2n(99/100)^2)\).
Since \(E_j''=j(A_j''-S)\), these estimates give (3) with
\(\omega=132S\eta\le10^{-9}\).
Endpoint normalization only rescales \(z\) by a positive factor.
The analytic comparison proves the theorem.

## A signed matched-tanh criterion in the same curvature range

There is also a sufficient condition that retains cancellation
between row orientations. Use the [matched tanh factors](MULTIVARIATE_MATCHED_TANH_MOMENTS.md)
at orders \(n,N\) and put
\[
 r_i=\frac{n}{n-1}T_i,\qquad A=\sum_i r_i,\qquad
 P=\sum_i\epsilon_i r_i^2,\qquad V_p=\sum_i r_i^p .
\]
The index ranges over Gauss factors, with orientation \(\epsilon_i\).
For \(n\ge600\) the two exact factorial moment polynomials have
disjoint root arguments if
\[
         8\le2A\le13,\qquad
        514|P|+176200V_3+4V_2/n\le10^{-9}.                    \tag{8}
\]

Apply the [midpoint kernel computation](MULTIVARIATE_MATCHED_TANH_ORDINARY_ANGLE_LIFT.md#exact-curvature-of-the-assembled-symbol)
with reference parameter \(S=2A\).
In (8), \(176200V_3\le10^{-9}\) forces every lower midpoint
parameter below \(1/2056\), and the exact matched identity
\(r_i'=r_i/(1-\xi_i)\), with \(|\xi_i|\le r_i/(nN)\),
gives \(r_i'<2r_i\). On \(|v|\le128\), \(g=2v-1\) has
modulus at most \(257\), so \(\Re(1+\epsilon_i r g)\ge3/4\)
on the entire path between the two matched midpoint parameters.
Writing \(d=3/4\) for this bound, the reciprocal expansion and
the kernel error estimate in the midpoint computation give
\[
 \|A_n''-2A\|
   \le 514|P|+\left(\frac{2\cdot257^2}{d}+
                           \frac{2}{3n^2d^3}\right)V_3 .
\]
The signed quadratic term is summed before taking its modulus.
The change of midpoint parameters is bounded by
\(\sum_i|r_i'-r_i|\le V_2/(nN-1/1028)\).
Differentiating the exact kernel in the parameter gives the
bound \(2/d^2\) for its derivative. For fixed parameter, the
change of grid from \(n\) to \(N\) costs at most
\(2 V_3(n^{-2}-N^{-2})/(3d^3)\).
Thus, multiplying the difference of the normalized curvatures
by \(N\) gives
\[
 N\|A_N''-A_n''\|
    \le 4 V_2/n +\frac{2V_3}{3d^3}\frac5{12}.
\]
Since \(2\cdot257^2/d+2(1/4+5/12)/(3d^3)<176200\),
\(E_j(v)=j[A_j(v)+A v(1-v)]\) satisfies (3) under (8).

For example, use the reflected rows \((\ell,\ell+d_0)\),
\((\ell+d_0,\ell)\) with \(\ell=8\cdot10^{15}\) and \(d_0=2\cdot10^8\).
Here \(P=0\),
\(d_0^2/(\ell+d_0)\le A\le d_0^2/\ell\),
\(V_2\le1/16000000\), and \(V_3\le1/1280000000000000\).
These verify (8) at all \(n\ge600\). For these rows the
relative-gap quantity \(S\eta=2d_0^3/(\ell+d_0)^2\) exceeds
the tolerance in (1).
