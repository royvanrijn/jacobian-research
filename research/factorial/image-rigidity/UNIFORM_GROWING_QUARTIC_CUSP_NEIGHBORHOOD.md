# Uniform consecutive exclusion in one physical quartic-cusp neighborhood

Scratch mathematical proof for independent review. No canonical source or
status is changed. The conclusion is local in the dimensionless coefficient
q, uniform in the INTEGER cusp parameter. It supplies no numerical cutoff
or assertion about common roots outside that neighborhood.

For s>=1 put

    a=s^2, d=2s+1, Sc=s(s+1), L=sqrt(Sc),
    wc=Sc^d*(s+1)/s, beta_s=-L/d,
    A_m(x)=sum_(j=0)^m binom(m,j)*(am+1)_(dj)*x^(m-j).

**Theorem UniformPhysicalCuspNeighborhood.** There are r>0 and N0,
independent of s, such that for all integers s>=1,N>=N0 the two ACTUAL
polynomials A_N,A_(N-1) have no common complex root x satisfying

    |x/(N^d*wc)-1|<r.                                     (V1)

The fixed-cusp shrinking-neighborhood theorem alone does not give this
uniform statement. We re-establish its complete estimates uniformly and
replace its support-dependent radial small-step argument by a two-endpoint
cone argument. A separate positive-side catch is essential: the lower
own-order parameter is not assumed local merely because the upper one is.

## 1. A compact analytic family and uniform complete real tails

Write nu=1/d in (0,1/3], include nu=0 as a limiting parameter, and set

    kappa=(s+1)/s=(1+nu)/(1-nu), t=(S-Sc)/L, p=1+t/L.

Subtract the central phase value. The EXACT local phase is

    F_nu(t,q)=a*Log(p)-L*t
                 +Log[(p^d+kappa*q)/(1+kappa*q)].           (V2)

It equals Phi(Sc+L*t,wc*q)-Phi(Sc,wc*q), where
Phi=a Log S+Log(w+S^d)-S. At nu=0 its limit is

    F_0(t,q)=-t-t^2/2+Log[(exp(2t)+q)/(1+q)].               (V3)

This is a compact holomorphic local family in (t,q), continuous and
analytic through nu=0. To see the apparent cancellation explicitly,
a/L-L=-sqrt((1-nu)/(1+nu)), a/L^2=(1-nu)/(1+nu), and
L^(-1)=2nu/sqrt(1-nu^2). Thus the power series of a Log(1+t/L)-Lt has
bounded analytic coefficients through nu=0. Likewise
d Log(1+t/L) has analytic coefficients, starting with
2t/sqrt(1-nu^2). A common sufficiently small complex t-disc and q-disc
about1 therefore work for the entire compact nu interval. In particular

    F_t=F_tt=F_ttt=0, F_tttt=-2 at (t,q)=(0,1),
    F_tq(0,1)=beta_s=-sqrt(1-nu^2)/2,
    sqrt(2)/3<=|beta_s|<=1/2.                              (V4)

All fixed local derivatives are uniformly bounded on smaller discs.

Complete tails require more than these local jets. Use the normalized
real v-family G_nu of UNIFORM_BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md,
with v=1+sqrt(kappa)*t. Its exact derivative and real-tail inequalities
(U3)--(U6) prove the following at q=1, uniformly in nu: t=0 is the unique
global quartic maximum; on a fixed local interval F_nu(t,1)<=-c*t^4;
outside any smaller fixed interval there is a fixed negative phase gap;
and both whole real tails have uniform exponential integral bounds.
For v<=-1/2 the derivative G_v>=1/2, and for v>=6 it is <=-3/7. These
inequalities include the receding left endpoint and the nu=0 whole-line
limit. The compact middle v-interval has p bounded away from zero, so
its strict gap follows from uniqueness and compactness. Fixed powers of
the phase are integrable with the same gap, as proved there.

For real t in its original domain (-L,infinity) and Re q>=1/2,

    |F_q(t,q)|<=4,
    |F_nu(t,q)-F_nu(t,1)|<=4|q-1|.                         (V5)

Indeed differentiate the two logarithms in(V2); each denominator has
positive real part and each ratio is at most2. Near t=0 the stronger
bound C|q-1||t| follows because F_q(0,q)=0 and F_tq is uniformly bounded.
These facts preserve the complete nonlocal gap BEFORE multiplying by m.
They also preserve it for a fixed small sector |arg m|<epsilon_m, since
the base phase is real and the imaginary perturbation is uniformly
bounded. Normalizing physical periods by their common central factor
and L removes all support-dependent scalar growth.

The branches used for actual integrals are legitimate. Choose a narrow
S-sector with |d arg S|<c<pi/2 and Re q>0. Both S and wc*q+S^d are
zero-free there. Every finite local deformation with |t| small fits in
this sector uniformly, since d/L is bounded. The positive-ray tails are
unchanged. Exact split periods are thus holomorphic in (m,q) for the
required complex-m sector. No logarithm on an unjustified full right
half-plane is asserted.

## 2. Uniform blown-up packets and their Cauchy collars

Set eta=beta_s*(q-1)=rho^3*omega, rho=|eta|^(1/3). For t=rho*u,

    rho^(-4)*F_nu(rho*u,q)=-u^4/12+omega*u+O(rho),          (V6)

with all fixed derivatives on bounded complex u-discs, uniformly in nu
and in a complex neighborhood of compact omega sets. The family extends
analytically at rho=0; the modulus rho is only a center size for estimates.

The fixed-cusp proof's complete horizontal contour construction now has
UNIFORM constants. Here are the details supplying that uniformity.
On the real local interval, the base quartic bound and(V5) give

    Re(F_nu(rho*u,q)/rho^4)<=-c*u^4+C|u|,
    |Im(F_nu(rho*u,q)/rho^4)|<=C|u|.                       (V7)

These bounds control the ENTIRE intermediate range from fixed |u|=R to
the fixed physical cut, not just a bounded u-interval. The far original
tails have the fixed gap from Section1. Relative to a Gaussian packet,
their normalized contribution is at most
C*sqrt(|m|)*rho*exp(-c|m|), hence O(1/|k|) for k=m*rho^4>=K.

Choose a small fixed epsilon>0. On the model outer-angle sets
Re omega>=sin epsilon or Re omega<=-sin epsilon, its selected saddle
horizontal has one strict global maximum, a uniform quadratic curvature
and a uniform gap elsewhere. These assertions follow from the explicit
quartic horizontal polynomial in canonical(C5). Compact parameter
continuation, fixed large u vertical cuts, (V7), and the complete tails
then give one nonzero actual Gaussian packet for all large k, uniformly
in nu. Thus every local actual zero with k>=K lies in

    |arg eta+pi/2|<epsilon                                (V8)

when Im q>0. A wider cone of width2epsilon is reserved for collars.

On that wider central cone, the two selected saddles have normalized
locations approaching 3^(1/3)*(+/-sqrt(3)/2-i/2). The exact midpoint split
defines P_L from0 to S_M and P_R from S_M to infinity; their sum is the
actual I_m. Define the paths by a fixed real split plus or minus one
finite path to S_M. They are exact holomorphic functions of (m,q).
The two actual saddle horizontals and short midpoint connectors have
strict gaps by the model identities, c=3^(1/3)/2, x_i=+/-sqrt(3)c:

    Re[f(u_i)-f(x-ic)]=(x^2-3c^2)^2/12,
    Im[f(x-ic)-f(u_i)]=(c/3)*(x-x_i)^2*(x+2x_i).

On the appropriate halfline the absolute imaginary difference is at
most five times the real loss. A fixed small complex-m sector therefore
preserves all losses. The local cubic Gaussian correction is odd on a
symmetric interval; quartic and squared-cubic remainders have uniform
integrable majorant C(|k|y^4+|k|^2 y^6)exp(-c|k|y^2). Combining this with
(V7) and the whole tails gives exact uniform expansions

    P_i=exp(m Phi_i)*sqrt(2pi/m)*c_i*(1+R_i),
    |R_i|<=C/|k|, c_i=(-Phi_SS)^(-1/2), i=L,R.              (V9)

The individual c_i have size L/rho. The common L/rho cancels in c_R/c_L;
it must not be bounded separately. Formula(V9) holds on product discs
|m'-m|<c*m, |q'-q|<c*rho^3 with the CENTER rho held fixed. Their constants
are uniform by the angular/physical collars just proved. Cauchy gives

    |R_(i,m)|<=C/(m*k), |R_(i,q)|<=C/(k*rho^3).            (V10)

Here the comparison K is enlarged by a fixed factor from the basic
Gaussian existence threshold, so every product disc remains above that
basic threshold. All subsequent uses retain this margin.

The Puiseux ratio c_R/c_L is a nonzero analytic unit in eta^(1/3), with
unit-modulus limit, uniformly in nu. Thus its local logarithm b0 obeys
|Re b0|<=C*rho, |b0_q|<=C/rho^2. For K large put

    D=Phi_R-Phi_L,
    b=b0+log(1+R_R)-log(1+R_L), E=m*D+b.

Uniformly on the wider cone,

    |Re b|<=C(rho+1/k), |b_m|<=C/(m*k),
    |b_q|<=C(rho^(-2)+1/(k*rho^3)),
    Re DeltaS>=c*L*rho, |DeltaS|<=C*L*rho, |D|<=C*rho^4,
    d*q*D_q=-DeltaS, DeltaS=S_R-S_L.                       (V11)

The last equality is exact stationarity. At every actual zero Re E=0.
All estimates come from exact functions, not derivatives of O-symbols.

## 3. A uniform positive-side catch for the lower endpoint

Fix the packet q-disc from Section2, of radius r_p>0. For any fixed
Q>2 there are delta_p>0 and M_p, uniform in s, such that the ACTUAL
period is zero-free for m>=M_p whenever

    1+r_p/2<=Re q<=Q, |Im q|<delta_p.                      (V12)

To prove this, use the canonical normalized real family on the compact
set nu in[0,1/3], real q in[1+r_p/2,Q]. Its unique maximizer v0 lies
in[c0,1) by(U3)--(U5). Formula(U7), or its exact second derivative,
gives uniformly negative curvature: q-1 is bounded below, so(U8) bounds
1-v0 below by a positive constant. Whole real tails obey(U6), and the
remaining compact real cycle has a uniform gap outside a fixed saddle
interval. The analytic implicit-function theorem continues this saddle
to one complex q-strip of uniform width. Move only a bounded local
interval to the saddle horizontal and join its ends back to the real
cycle below the fixed gap. The global q-perturbation bound from the
positive denominators preserves the real tails. A complete Gaussian
estimate, with uniform nonzero leading coefficient and O(1/m) relative
error, proves(V12). This is an ordinary positive-real regular saddle
argument; it does NOT apply the bounded-H theorem at a growing H.

Exclude soft pairs first. The OLD Gamma cor:linear-lower-exponent gives
coprimality if a>=d(2N-3). At a cusp its strict complement implies
d/(N-1)<8; hence lambda_1=(N/(N-1))^d<3^8. Put Q=2*3^8 in(V12).

Choose the claimed upper-disc radius r MUCH smaller than r_p. If q_N is
in |q_N-1|<r, its lower own-order parameter q_n=lambda_1*q_N has positive
real part, bounded above by Q, and |Im q_n|<=3^8*r. Choose this smaller
than both r_p/4 and delta_p. If |q_n-1|>=r_p, then Re q_n>=1+r_p/2:
the alternative left-hand side of that circle is impossible because
Re q_n>=1-r>1-r_p/2. Formula(V12) would make the lower period nonzero.
Consequently a common root forces BOTH endpoints into the packet disc.
This catch is logically prior to the lower endpoint's angular label.

## 4. Endpoint cones replace a support-dependent small-step argument

Choose K for Section2 and consider k_N=N|eta_N|^(4/3)>=8K. The upper zero
has(V8). Put n=N-1, lambda_theta=(N/(N-theta))^d, and

    eta_theta=beta_s*(lambda_theta*q_N-1)
             =lambda_theta*eta_N-|beta_s|*(lambda_theta-1).

If y_N=-Im eta_N>0, then -Im eta_n=lambda_1*y_N. BEFORE assigning a
lower angle this already implies

    k_n>= (n/N)*cos(epsilon)^(4/3)*k_N>=K.                 (V13)

Take epsilon small and reserve a fixed threshold margin. Section3 put
q_n in the packet disc, so its zero also satisfies(V8).

Let z_theta=Re eta_theta/(-Im eta_theta). EXACTLY

    z_theta=z_N-|beta_s|*(1-1/lambda_theta)/y_N.            (V14)

It decreases monotonically between the two endpoint ratios. Since both
lie in[-tan epsilon,tan epsilon], the ENTIRE path stays in that central
cone; no assumption of a small d* rho/k step is used. Moreover

    1-1/lambda_1<=2*tan(epsilon)*y_N/|beta_s|.

Choosing r small gives lambda_1<=4/3 and rho_theta comparable to rho_N,
and in fact lambda_1-1<=C*epsilon*y_N, uniformly in s,N.
In particular k_theta>=K throughout the path and a
fixed physical collar remains available. The common-root hypothesis at
BOTH integer endpoints, not an a priori bound on the lower scaled field,
is what supplies this conclusion.

## 5. Exact log-ratio comparison and the bounded core

Put m=N-theta and q_theta=lambda_theta*q_N. By(V11),

    dE/dtheta=-D-DeltaS-b_m+(d*q_theta/m)*b_q.             (V15)

At the upper zero, |Re D|<=C(rho_N/N+1/(N*k_N)). Along the path the exact
identity dD/dtheta=-DeltaS/m gives

    |Re D(q_theta)|<=C(rho/N+1/(N*k)+L*rho/N),

where comparability has absorbed fixed constants. Since d/L is uniformly
bounded and L>=sqrt(2), dividing ALL errors in(V15) by L*rho gives

    C*(1/N+rho/k+1/k^2).                                 (V16)

For example the b_q terms give rho/k and1/k^2 exactly, using m=k/rho^4;
the upper-zero 1/(N*k) term divided by L*rho is bounded by C/k^2. Choose
K large, the physical radius small, and N large with constants uniform
in s. Then d Re E/dtheta<=-c*L*rho/2<0, whereas the two actual endpoint
zeros require Re E=0 at both ends. This contradiction excludes every
growing-field common root in(V1).

If k_N<=8K, the field h_N=N^(3/4)*eta_N satisfies
|h_N|<=H=(8K)^(3/4). The canonical UniformBoundedQuarticCuspExclusion
theorem supplies a single N_0(H), independent of s. This is a FIXED H
after K has been chosen. It excludes this remaining bounded core.
Positive real q near1 has no actual zero by positivity, and the lower
half-plane follows by conjugation. Taking the maximum of the uniform
thresholds proves(V1).

## Scope and exact antecedents

The analytic imports are the compact normalized family, whole real
tail/concentration estimates, and bounded-H theorem of
UNIFORM_BOUNDED_QUARTIC_CUSP_CONSECUTIVE_EXCLUSION.md; the fixed-cusp
complete construction C1--C15 is revalidated here with its constants
uniform. The old Gamma soft corollary is used only to bound lambda_1
before the lower endpoint catch. No new sharp lower annulus, global
q-root localization, or identification of every complex stationary
branch is required. The proof gives no explicit r,N0, no all-order
claim, and no information about common roots outside(V1).
