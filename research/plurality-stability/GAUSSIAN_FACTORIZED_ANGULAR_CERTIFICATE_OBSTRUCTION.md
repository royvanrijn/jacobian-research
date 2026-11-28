# Obstruction to factoring the mixed-radius angular certificate

In polar coordinates for the Gaussian plane, write \(d\mu(r)=
r e^{-r^2/2}\,dr\) and let \(t\in[0,\pi]\) be the smaller angular
separation. At \(0<\rho<1\) the relative Mehler kernel is
\[
 M_\rho(r,s,t)=C_\rho e^{-\kappa_\rho(r^2+s^2)}
                   e^{\omega_\rho rs\cos t},
 \qquad
 C_\rho=\frac1{1-\rho^2},\quad
 \kappa_\rho=\frac{\rho^2}{2(1-\rho^2)},\quad
 \omega_\rho=\frac{\rho}{1-\rho^2}.
\]
The exact one-circle triangle certificate accepts any decreasing
angular weight. It also accepts a distribution of radii when the
two endpoint radii are sampled independently from that same
distribution. The Mehler weight cannot be represented on the
whole active angle interval by a positive superposition of such
independent-radius certificates. The obstruction persists under
radial centering and addition of row potentials.
[The exact circle certificate](GAUSSIAN_AGREEMENT_TRIANGLE_RELAXATION.md#the-angular-inequality-and-an-exact-positive-certificate)

## A negative centered radial test on the active interval

Set \(a=2\pi/3\), the angular support endpoint of the averaged
simplex agreement. Consider a factorization on \(0<t<a\) of
the form
\[
 M_\rho(r,s,t)= A(r,t)+A(s,t)
          +\int k_z(t)\,v_z(r)v_z(s)\,d\eta(z),             \tag{FA1}
\]
where \(\eta\) is positive and each \(k_z\) is nonincreasing.
The radial weights for an independent-sampling certificate
are nonnegative, but the obstruction even allows signed
weights. Assume the pairings are defined and the integration
may be interchanged on bounded compactly supported radial
tests, as needed for the proposed certificate. The terms
\(A\) can depend on angle. One may also add any term supported
on \(t\ge a\), since it does not affect this interval.

For each bounded compactly supported \(f\) with
\(\int f\,d\mu=0\), (FA1) would make
\[
          t\longmapsto
          \iint f(r)M_\rho(r,s,t)f(s)\,d\mu(r)d\mu(s)
\]
nonincreasing. The row terms vanish, and each \(k_z\) is
paired with the nonnegative square \((\int f v_z\,d\mu)^2\).
This argument uses only monotonicity, without differentiating
the proposed factorization. Since the displayed Mehler
pairing is smooth, its negative derivative would have to
be nonnegative on every such centered test.

Fix \(t_0\in(\pi/2,a)\), set \(\omega=-\omega_\rho\cos t_0>0\),
and put \(K=-\partial_t M_\rho(\,\cdot\,,\,\cdot\,,t_0)\).
With a positive constant \(c\),
\[
       K(r,s)=c\,rs\exp[-\kappa_\rho(r^2+s^2)-\omega rs].
\]
The diagonal \(K(r,r)\) is positive for \(r>0\), vanishes
in the limits \(r\downarrow0\) and \(r\to\infty\), and has
a single maximum. Choose distinct \(r_1,r_2>0\) with
\(K(r_1,r_1)=K(r_2,r_2)>0\). Then
\[
 \frac{K(r_1,r_2)}{\sqrt{K(r_1,r_1)K(r_2,r_2)}}
       =\exp\!\left(\frac{\omega}{2}(r_1-r_2)^2\right)>1 .
\]
Taking the difference of the two normalized radial indicator
bumps on small intervals about \(r_1,r_2\) gives a mean-zero
test. Its quadratic form tends to
\(K(r_1,r_1)+K(r_2,r_2)-2K(r_1,r_2)<0\) as the intervals
shrink. Thus the Mehler pairing has positive angular derivative
at \(t_0\) for this test, contrary to (FA1).

This also excludes such a factorization of
\(P_{1^\perp}M_\rho(t)P_{1^\perp}\), the radially centered
kernel, even with arbitrary angle-dependent row terms. If the
factorization is stated almost everywhere in \(t\), the tested
smooth pairing would still agree almost everywhere with a
nonincreasing function, which gives the same contradiction.

## An exact cumulative bound for angles below \(\pi/2\)

For \(0<t<\pi/2\) the same derivative has a positive
independent-radius expansion. In fact
\[
 -\partial_t M_\rho(r,s,t)
 =C_\rho\omega_\rho\sin t
      \sum_{n\ge0}\frac{(\omega_\rho\cos t)^n}{n!}
               v_n(r)v_n(s),\qquad
 v_n(r)=r^{n+1}e^{-\kappa_\rho r^2}.                        \tag{FA2}
\]
This can be combined with a triangle comparison without a
conditional row assumption. Let \(F\) be an even rotational
agreement kernel in the triangle relaxation and put
\[
    P(r,s)=\frac1\pi\int_0^\pi F(r,s,u)\,du,\qquad
    L_t(r,s)=\int_0^t\big(F(r,s,u)-F_T(u)\big)\,du .
\]
Then for \(0<t<\pi/2\),
\[
 \iint (-\partial_tM_\rho)L_t\,d\mu^2
       \le \frac{\pi t^2}{a^2}
           \iint(-\partial_tM_\rho)(P-1/3)\,d\mu^2.        \tag{FA3}
\]

To prove this, average the radii independently against any
probability density \(q\) with respect to \(\mu\). Write
\(F_q(t)=\iint q(r)q(s)F(r,s,t)\,d\mu(r)d\mu(s)\) and
\(p_q=\pi^{-1}\int_0^\pi F_q\). The triangle inequality passes
to \(F_q\) by sampling three radii independently. For
\(H_q(t)=\int_0^t(1-F_q)\), the integrated triangle inequality
gives \(tH_q'(t)\le2H_q(t)\) almost everywhere on \((0,\pi)\).
Thus \(H_q(t)/t^2\) is nonincreasing. The box condition
gives
\[
      H_q(a)\ge a-\pi p_q=\frac a2-\pi(p_q-1/3).
\]
For \(t\le a\) it follows that the averaged \(L_t\) is at
most \((\pi t^2/a^2)(p_q-1/3)\). Apply this to the normalized
weights \(v_n\) of (FA2), and sum with their nonnegative
coefficients. The absolute pairings converge by boundedness
of \(L_t,P\) and integrability of the derivative kernel.

The contribution in (FA3) is nonpositive when all mixed
angular means are \(1/3\). In general it contains the radial
mean correction, and above \(\pi/2\) the derivative no longer
has a positive independent-radius expansion. These are
specific constraints on extending the one-circle certificate
to the full pair objective. The calculation does not rule out
a mixed-radius dual whose triangle multipliers involve all
three radii jointly.
