# A uniform quartic large-lower-exponent Strong Factorial tail

The uniform analytic hypotheses(U1)--(U13) are proved in UNIFORM_ACTUAL_PHASE_AND_CONTOUR_ESTIMATES.md. The independent global root-strip input below is the retained ITERATED_GRACE_POLYNOMIAL_ROOT_STRIP.md. This proof is separate from the quadratic Bell/homotopy route.



The independently derived movable-centre Grace lemma gives, for all d>=2,N>=1,a>=100d^2, the actual monic moment roots in the normalization(U1):

|Re z|<7*d^2/sqrt(a), |Im z|^2<17*d.                     (U14)

Its source is ITERATED_GRACE_POLYNOMIAL_ROOT_STRIP.md. The second bound is no longer needed for the following argument, because Section4 covers every outside height.

**Theorem UniformQuarticStrongFactorialTail.** For all integers d>=2, a>=10^20*d^4 and N>=2, the actual consecutive ordinary factorial moments of a nonzero X^a(A+B*X^d) cannot both vanish, over every characteristic-zero field.

Proof using the uniform analytic lemmas and the independently proved root strip. The actual monic polynomial is

A_(m;a,d)(Z)=sum_(k=0)^m binom(m,k)*(a*m+1)_(d*k)*Z^(m-k).

After nonzero scalars its moment integral is

I_m^orig(w)=integral_0^infinity [S^a(w+S^d)]^m exp(-m*S)dS,
w=Z/m^d=a^d*(-1+delta*z),
I_m^orig(w)=sqrt(a)*[a^(a+d)*delta*exp(-a)]^m*I_m(z).

The exact moment identity is(am)!*A_m(Z)=m^((a+d)*m+1)*I_m^orig(Z/m^d), with I_m the normalized actual F-period. Thus no fixed-Gamma or common-shape surrogate is substituted.

The exact large-lower criterion already covers a>=d*(2N-3), including equality. In the complementary hard range N>a/(2d)+3/2, n=N-1 and N both exceed4*10^20. The proposed threshold gives

delta=d/sqrt(a)<=10^(-10)/d<=5*10^(-11),
|Re z|<7*10^(-10)

by(U14). If |Im z|>=2.1, the upper actual period is nonzero by Section4, regardless of its height. Conjugation reduces the rest to the upper half-plane.

The exact actual lower-order argument is

z_n=z+alpha*(z-delta^(-1)), alpha=(N/n)^d-1.              (U15)

In the hard range d/n<2delta^2; since log(N/n)<1/n,

alpha<exp(2delta^2)-1<3delta^2.

The same bound holds along the whole coupled radial segment. For |z|<=3, its real shift is at most3delta+9delta^2 and its positive imaginary part is multiplied by at most1+3delta^2. Thus an upper bulk point |Im z|<=1.9 stays inside the intermediate bulk rectangle of Section5, and an upper endpoint point1.9<Im z<2.1 stays inside the full expanded endpoint rectangle. In both cases the real part is far below the required1e-7/2e-7 or1e-5 allowances.

The central case cannot give simultaneous zeros by(U10) and its complete packet bounds. In the endpoint case, |zeta_z+i|<1/4 along(U15) gives

Im[zeta(z_n)-zeta(z)]>=alpha/(2delta)
 >=d/(2delta*N)=sqrt(a)/(2N).

But(U13) at the two own actual parameters bounds its modulus by

5200*(1/N+1/n)<15600/N.

The lower bound is at least(10^10*d^2)/(2N), a contradiction. This excludes all hard-order common roots. The soft exact criterion covers every remaining N, so every consecutive pair of rational monic polynomials is coprime over Q. Its rational Bezout identity transfers to every characteristic-zero field; monomial cases never have a zero positive factorial moment. This proves the theorem.

