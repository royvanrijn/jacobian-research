# Next target: uniform exclusion across the quartic transition

This is an OPEN research target, not a proved extension of the uniform exact-cusp theorem.

The completed exact-cusp route has a=s^2,d=2s+1. Its next extension should let the quadratic unfolding move simultaneously with the complex coefficient field. A fixed-support theorem and a theorem uniform on the equality curve do not imply this two-parameter statement.

For EVERY integer gap d>=2 put s=(d-1)/2 (a positive integer or half-integer), e=a-s^2 with integer a>=1, C=a+s, L=sqrt(C), kappa=(s+1)/s, w0=C^d*kappa, q=w/w0, and t=(S-C)/L. The actual normalized local phase is

    F(t,q)=a log(1+t/L)-L t
          +log[((1+t/L)^d+kappa*q)/(1+kappa*q)].

At q=1,t=0 the point is stationary for EVERY e, and direct differentiation gives

    F_t=0, F_tt=-e/C.

The higher raw jets are F_ttt=2e/C^(3/2) and
F_tttt=-[2(s(s+1))^2+6e]/C^2; hence at e=0 they are0 and-2.
The exact mixed first derivative is

    F_tq(0,1)=-s(s+1)/(d*L),

which equals -L/d only on the equality curve. It must not be frozen at its cusp value when e varies.

The natural bounded transition variables are

    u=-sqrt(m)*e/C,
    h=m^(3/4)*[-s(s+1)/(d*L)]*(q-1).

The candidate leading local integral is a real-quadratic deformation of the quartic Fourier-Laplace kernel. The retained Newman/Lee-Yang bridge already supplies relevant model-zero results; inspect its EXACT range of real u and simplicity before importing them. The raw moving-center cubic jet also needs a uniform bound. A pointwise formal Taylor expansion is insufficient.

An integrality advantage is available in the bounded thermal window. If e is nonzero then |e|>=1/4, so |u|<=U implies C>=sqrt(m)/(4U) and L>=m^(1/4)/(2sqrt(U)). The own-order radial displacement of the scaled complex field consequently need not be small. This observation is only a prospective argument: its constants, lower-endpoint capture and complete uniform approximation still have to be proved.

For growing field set rho=|beta_e(q-1)|^(1/3), k=m*rho^4. The thermal deformation seen by the blown-up phase is T=-(e/C)/rho^2=u/sqrt(k). For k large compared with U^2 the saddle geometry is a small perturbation of the exact-cusp geometry. The leading amplitude may now depend on its angular field at order T: its q-derivative can include O(T/rho^3), rather than only the exact-cusp O(rho^-2). This adds O(U/k^(3/2)) to the relative transport budget and must be explicitly retained.

The immediate action is to prove complete actual-period expansions with errors, derivative collars, and lower-endpoint capture UNIFORM on bounded real u and bounded/growing h. Under adjacent-order comparison e,C,s stay fixed but u_m varies as sqrt(m); include that thermal drift in the exact log-ratio derivative. Match unbounded u to the regular/simple-turn or two-well packets on the correct sides of the transition. The half-integer s convention includes even gaps; exact integer cusps occur only when e=0 and s is integral.

The relevant classical model is already identified in the retained Curie-Weiss bridge. The new task is zero-sensitive transfer to actual Gamma periods while m/d may be unbounded. An exact finite-Gamma circle law is false and is not a candidate premise.

Success would extend uniform control from a single critical support curve into a neighborhood of that curve. It would still leave the global remaining support region and the finitely many low terminal orders. No global finite-bad-triple result, all-order theorem, numerical cutoff, external review, or formal verification is asserted here.
