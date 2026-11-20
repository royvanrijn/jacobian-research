# Large-thermal matching through the actual lower quartic fold

This is a scratch proof for adversarial review. No canonical file or status is
changed. The statement is sequential near the critical support curve; it gives
no numerical cutoff, no all-order assertion, and no uniform theorem on all
ordinary supports.

For integers a>=1,d>=2 put

    s=(d-1)/2, e=a-s^2, C=a+s, L=sqrt(C), Sc=s(s+1),
    r=e/C, beta=-Sc/(dL), kappa=(s+1)/s, w0=C^d*kappa,
    q_m=x/(m^d*w0), eta_m=beta*(q_m-1), u_m=-sqrt(m)*r.

The original moment polynomials and actual positive-cycle integral are the
ones in UNIFORM_BOUNDED_THERMAL_TRANSITION_EVENTUAL_SFC.md. The raw phase is

    F(t,q)=a Log(1+t/L)-Lt
          +Log[((1+t/L)^d+kappa*q)/(1+kappa*q)].

A nonzero detuning has |e|>=1/4, because e is quarter-integral.

**Sequential large-thermal claim.** There is no sequence of actual consecutive
common roots with N tending to infinity, nonzero r tending to0, and
sqrt(N)*|r| tending to infinity.

All limits below refer to a putative counterexample sequence. They are not
coefficient-relative approximations of the actual integral.

## 1. Scales, actual phase, and global capture

Write

    tau=sqrt(|r|), sigma=-sign(r) in {+1,-1},
    gamma=L*tau=sqrt(|e|)>=1/2,
    M_m=m*tau^4=u_m^2,
    omega_m=eta_m/tau^3.                                (L1)

After a subsequence sigma is fixed. Since |e|>=1/4 and r->0, C,s,L tend to
infinity, d/L->2 and beta->-1/2. In a fixed bounded v-disc, t=tau*v gives

    f(v,omega)=tau^(-4)*F(tau*v,1+tau^3*omega/beta)
       -> f_sigma(v,omega):=-v^4/12+sigma*v^2/2+omega*v.   (L2)

The convergence has every fixed analytic derivative and is compact through
tau=0,nu=1/d=0. To see the parameter regularity, take r=-sigma*tau^2 and
C=Sc/(1+sigma*tau^2). The exact local jets are

    F_tt=-r, F_ttt=2r/L,
    F_tttt=-2(1-r)^2-6r/C,
    F_tq=beta, F_ttq=0.

Their convergent Taylor series show (L2). In particular the cubic correction
at this scale is O(tau/L), not an uncontrolled term involving M*tau.

The complete original real tails satisfy, on a fixed local real interval,

    Re(F(t,q))<=-c*t^4+C*|r|*t^2+C*|eta|*|t|.

Consequently in the entire intermediate range R<=|v|<=t0/tau,

    Re f(v,omega)<=-c*v^4+C*v^2+C*|omega|*|v|.           (L3)

For bounded omega this supplies integrable quartic tails uniformly. Outside
the fixed physical interval |t|<t0 the original real phase retains a fixed
negative gap: the r=0 phase has its unique quartic maximum at0; the compact
middle interval is stable, and the exact derivative identity

    partial_t Log Psi=(-s-Lt+dH)/(L p), 0<=H<=1,

controls BOTH whole original tails including t=-L. Small r and q-1 preserve
that gap before multiplying by m. Finite vertical connectors at a large
fixed v-cut R also have loss by the -v^4 term. Throughout, the leading
exponent uses the ACTUAL critical value.

The global entire-phase compactification and annulus in Section2 of the
bounded-thermal proof use only r->0,s->infinity, not bounded u. Thus they
apply here without change: both actual endpoint roots have q_m in one
fixed compact annulus, and joint s,m->infinity circle capture gives

    q_N->q_*, q_n->q_*, |q_*|=1,
    lambda=(N/n)^d->1, d/N->0, n=N-1.                   (L4)

If q_*!=1, the uniform regular two-packet log-ratio excludes the sequence.
Therefore q_N,q_n->1. By conjugation assume Im q_N>0. Then Im omega_N,
Im omega_n<0. Real q near1 is positive and cannot be a root.

## 2. Unbounded omega is already the small-temperature quartic regime

Suppose |omega_N|->infinity. Put rho=|eta|^(1/3) and k=m*rho^4.
The residual temperature in the field blowup is

    T=-r/rho^2=sigma/|omega|^(2/3),
    k=M_m*|omega|^(4/3).                               (L5)

Thus T->0 and k->infinity. The complete small-T construction of the
bounded-thermal proof Section4 is valid independently of a bound on u:
its estimates are uniform for small |T|,rho,nu, and its transported error
term is |T|/k. The representation T=u/sqrt(k) does not impose an upper
bound on u when T is directly known small.

Every upper zero is consequently in the central eta cone. Since
Im eta_n=lambda*Im eta_N, this forces |omega_n|->infinity as well and
k_n>=c*k_N before assigning a lower cone. The lower zero is therefore
central too. The exact endpoint ratio monotonicity keeps the whole path
central, with comparable rho,k. The comparison errors divided by the
main L*rho are

    O(1/N+rho/k+|T|/k+1/k^2)=o(1).

The exact log ratio has a strict negative real drop, contradicting the two
integer endpoint zeros. Hence a remaining sequence has bounded omega_N.

The lower omega_n is bounded too. Its imaginary part is lambda times the
bounded upper imaginary part. If its modulus were unbounded, (L5) would
put that lower zero in a central eta cone, whose imaginary part is a fixed
fraction of its modulus; this contradicts the bounded imaginary part.
Thus both endpoint omega fields lie in a fixed compact set.

## 3. Direct complete contour geometry of the two thermal models

For either sigma and v=x+iy, the x-dependent real part of f_sigma is

    V(x)=-x^4/12+(y^2+sigma)*x^2/2+Re(omega)*x.          (L6)

If Re omega>0, reflection strictly favors x>0. On that halfline,

    V'(x)=-x^3/3+(y^2+sigma)*x+Re omega

has exactly one positive zero. If y^2+sigma<=0 it is strictly decreasing.
If y^2+sigma>0 it first increases and then strictly decreases to minus
infinity. Its zero has V''<0 and is the unique global maximum.
Reflection gives the analogous left maximum for Re omega<0.

For each right-lower omega quadrant the cubic stationary map

    omega=v^3/3-sigma*v                               (L7)

has an analytic selected root with Re v>0,Im v<0. Begin on the positive
real omega axis at its favored positive real maximum. Critical values
are on the quadrant boundary (real for sigma+, imaginary for sigma-).
Cubic properness prevents escape on compact target sets. The selected
root cannot cross Re v=0 or Im v=0, since (L7) maps each corresponding
axis into an axis. It therefore continues simply throughout the quadrant.
The horizontal through this root has the unique favored maximum of(L6).

Original real tails and finite vertical cuts have strict quartic loss by
(L3). Uniform analytic perturbation on a fixed finite v-rectangle transfers
the unique horizontal maximum to the ACTUAL f. A symmetric local Gaussian
calculation gives a nonzero leading term and relative O(1/M) error. Hence
every bounded-omega actual-zero limit satisfies Re omega=0.

On omega=-iz, z>=0, the active model pair, when present, is

    v_R=x-ik, v_L=-x-ik,
    x^2=3(k^2+sigma), z=(8/3)k^3+2sigma*k.              (L8)

For sigma=+, every z>=0 has a unique k>=0. Both simple maxima are present
on one horizontal, with curvature -2(k^2+1)<0. Their midpoint loss is
x^4/12>0. The real model folds omega=+-2/3 lie on the DISFAVORED side of a
tilted real contour; its favored maximum is simple and strictly higher.
They do not require an Airy contribution to the original period.

For sigma=-, 0<=z<2/3 has the one dominant saddle v=-iy,0<=y<1, determined
by z=y-y^3/3. Its real horizontal curvature is y^2-1<0 and its modulus has
one global maximum. Thus this interval is zero-free away its endpoint.
For z>2/3, (L8) has unique k>1 and two strict simple maxima, curvature
-2(k^2-1)<0. The sole remaining model turn is

    sigma=-1, omega_f=-2i/3, v_f=-i.                    (L9)

All statements are complete horizontal-contour facts, not a claim that a
high saddle is active merely because its action is large.

Both bounded endpoint fields now have limiting real part0. Their exact
difference is

    omega_n-omega_N=beta*q_N*(lambda-1)/tau^3,
    omega_theta'=-gamma*(1-r)*q_theta/M_theta.          (L10)

Because beta is real negative, q_N->1, and lambda->1, equality of both
limiting real parts forces

    gamma/M_N->0, omega_n-omega_N->0.                  (L11)

Here Bernoulli's lambda-1>=d/n gives the converse lower scale. Both
endpoints therefore converge to the SAME bounded imaginary omega point.

## 4. Regular bounded omega and the two-well side

Away from(L9), a limiting imaginary-axis point is either zero-free by the
single packet just proved or has the two simple maxima of(L8). In the
latter case choose a fixed small omega-disc. The actual limiting contours
and full tails transfer uniformly, giving two exact midpoint pieces

    P_i=lambda_i^m*m^(-1/2)*c_i*(1+R_i), |R_i|<=C/M.

The common L/tau scale in c_i cancels from the amplitude ratio; its omega
derivative is bounded on that fixed disc. Exact stationarity gives

    d*q*D_q=-DeltaS, DeltaS=gamma*(v_R-v_L).             (L12)

Only integer endpoint values of R_i are needed here. At a zero,

    Re[mD+b0]=O(1/M), b0=Log(c_R/c_L).

By(L11) the whole actual omega path remains in this fixed disc, with
Re(v_R-v_L)>=c>0. The derivative of the bare expression mD+b0 is

    -D-DeltaS+(d*q/m)*partial_q b0.

Since omega_q=beta/tau^3, the last term divided by gamma is O(1/M).
At the upper zero Re D=O(1/m); its path change is O(gamma/m).
Thus the main real drop is at most -c*gamma/2 for all large members.
Both endpoint values are O(1/M), whereas gamma>=1/2. This is impossible.

This handles ALL bounded omega on the two-well side sigma+, including its
coexistence point0. On the single-well side only the fold(L9) remains.

## 5. The actual lower-fold cycle and the inactive high saddle

At the model fold the exact finite horizontal identity is

    f_-(x-i,-2i/3)=-1/4+i*x^3/3-x^4/12.                (L13)

The third stationary point v=2i has action2, but it is NOT on this contour.
For any fixed epsilon>0, the pieces |x|>=epsilon on the horizontal have
a strict modulus loss epsilon^4/12 relative to the folded value. Large
fixed vertical cuts and the original real tails have loss by(L3).
First choose the fixed cuts and epsilon, then take tau,nu small.
The same complete outside loss persists for the ACTUAL fold and nearby
fields. This proves inactivity of the high third saddle on this original
cycle; it is not assumed from an action ordering or a formal local model.

The analytic equations f_v=f_vv=0 have a unique actual fold
(v_f,omega_f) near(-i,-2i/3). Their Jacobian in(v,omega) is invertible at
the model because f_vomega=1 and f_vvv=2i. Uniform local analytic cubic
normalization consequently gives

    f(v,omega)=A(omega)+z^3/3-zeta(omega)*z,
    v=V(z,omega), V_z near i,
    zeta(omega_f)=0, zeta_omega near -i.                (L14)

The parameters tau,nu range on one small compact analytic neighborhood;
the radii, nonzero derivative bounds and all fixed analytic derivative
bounds are COMMON. This follows by the same explicit local analytic
inverse construction as the retained simple-turn proof A1--A3; its
invertible linear operator and derivative bounds are continuous near the
single model germ above. The new complete-cycle assertion is supplied
by(L13), not imported from a fixed-support cutoff.

At the model fold write v=-i+x. With the branch V_z=i the canonical map is

    z=-i*x+x^2/12+O(x^3).                               (L15)

Indeed substituting into z^3/3 gives i*x^3/3-x^4/12+O(x^5).
Both real local endpoints therefore enter the RIGHT two cubic decay
sectors, of arguments near +-pi/2 but bent toward the right half-plane.
Those sectors deform to the Airy rays +-pi/3. Their orientation and
dv=i dz give the leading scalar 2pi*Ai, not a Scorer function.
The scalar is uniformly nonzero under the actual compact perturbation.
The fixed outer pieces in(L13) lie below the fold action and cannot add
another leading local contour class.

For this step DLMF2.4(v),36.12(ii) give classical cubic/Airy frameworks,
but do not select an actual contour automatically;2.4 explicitly allows
Scorer approximants for other contours. The actual right-exit selection,
the whole original tails and the exclusion of the high saddle have just
been established by the finite horizontal deformation.



## 6. Uniform complete cubic packets near the actual fold

Choose one fixed local cubic radius in(L14), with the outer connectors
from Section5 already below the folded height by a fixed amount. Shrink
the omega-disc and parameter neighbourhood. The following estimates have
COMMON constants and refer to the actual map and actual action values.

Put

    zeta=zeta(omega), rhoA=|zeta|^(1/2), K=M*rhoA^3.

For |arg zeta-pi|>=epsilon and large K the selected Airy contour has one
nonzero Gaussian leading term. This statement includes canonical Stokes
directions away the negative real ray. A precise contour construction is
available without assuming that a bounded real-phase band is compact:
in the scaled cubic P(z)=z^3/3-chi*z, |chi|=1, continue each arm by the
proper inverse equation P(z)=P(z_i)-h, h>=0. A bounded complex P value has
bounded preimages; arm velocities are O(h^-1/2) at a simple saddle and
O(h^-2/3) at infinity. When a high descending arm meets the other lower
saddle, stop at one fixed intermediate sublevel, use a fixed small Morse
disc detour, and continue within that lower component. The detour is
exponentially smaller than the high leading packet. Compact angular
subsets give a finite cover and common constants.

The selected branch and its nonzero coefficient are the ordinary Ai
class just identified in Section5. For constant amplitude its sectorial
coefficient can also be read from DLMF9.7.5, valid on every closed sector
away the negative ray. A uniformly analytic amplitude multiplies the
local leading coefficient by its nonzero value at that saddle; Taylor
expansion and the proper-arm loss give O(1/K) relative error. Thus this
model fact is not being promoted to an actual-cycle assertion: the actual
class, map, outer connectors and complete tails were separately verified.

Every actual zero with K large must therefore satisfy

    |arg zeta-pi|<epsilon.                              (L16)

For a wider negative cone, choose the analytic rhoA=sqrt(-zeta), positive
when zeta is negative real, only for defining contours. The common
canonical midpoint z_M=-rhoA has phase -4rhoA^3/3 at zeta=-rhoA^2,
strictly below BOTH critical phases +-2i rhoA^3/3. This gives a uniform
cubic-size gap, unlike the weak quartic-size midpoint loss on the original
horizontal. To identify the shared connector explicitly, at zeta=-rhoA^2
scale z=rhoA*y. The cubic y^3/3+y has critical values +-2i/3 and no
critical value in the open left value-half-plane. Its branch containing
the negative real y-axis is an unbranched proper degree-one sheet there;
its boundary limits at +-2i/3 are respectively +-i, as is seen from
y=i*v, v-v^3/3 running from0 to2/3. Thus the two left descending arms
belong to the SAME left sheet. On the fixed value-line Re(phi)=-4/3, both arms join to y=-1 by lifted
vertical segments still strictly below both critical heights.
Compact angular perturbation preserves this construction. Lift the
two proper cubic paths through the actual inverse V, connecting them
to their original finite exit paths. Splitting the
original integral exactly at V(z_M,omega) defines two periods P_L,P_R
with I=P_L+P_R. The left and right original tails remain unchanged.

The local branches and their logarithms lie in one bounded zero-free
S-sector about the large positive centre C: d*arg S=O(tau), and
q=1+O(tau^3). Hence these exact periods are holomorphic for a fixed small
complex-m sector and in the wider omega cone. At original real endpoints
the powers have their positive-ray branches; only integer m identities
with the factorial polynomial are needed. The same complete-tail loss
holds in that sector, and no original-cycle identity for a path crossing
a phase zero at arbitrary noninteger m is assumed.

Uniform local Gaussian expansion on the proper cubic paths gives

    P_i=exp(m Phi_i)*sqrt(2pi/m)*c_i*(1+R_i),
    |R_i|<=C/K,
    |partial_m R_i|<=C/(m*K),
    |partial_omega R_i|<=C/(K*rhoA^2).                  (L17)

The last two inequalities are genuine Cauchy estimates on product discs
|m'-m|<c*m and |omega'-omega|<c*|zeta|; support r,s remains FIXED.
They follow from the uniform complete O(1/K) estimate just established,
reserving angular, physical and threshold margins first. Analytic
amplitudes and their fixed derivatives are bounded on the fixed cubic
disc. Their exterior tails are O(exp(-cM)); choose |zeta| small enough
that this remains below both actual selected packet scales before
multiplication by M. At a Stokes direction the fixed-sublevel detour has
the same strict loss. Thus no unestimated nonlocal term occurs in(L17).

Let D=Phi_R-Phi_L. The common singular factor in c_i cancels in their
ratio. By the analytic cubic map, that ratio tends to a constant of
modulus1 at zeta=0, uniformly in the compact actual fold parameters.
For its logarithm b0 and the exact remainder-corrected b one has

    |Re b0|<=C*rhoA, |b0_omega|<=C/rhoA,
    b=b0+Log(1+R_R)-Log(1+R_L),
    |Re b|<=C(rhoA+1/K), |b_m|<=C/(mK),
    |b_omega|<=C[rhoA^(-1)+1/(K*rhoA^2)].               (L18)

With labels oriented on the negative cone,

    Re DeltaS>=c*gamma*rhoA, |DeltaS|<=C*gamma*rhoA,
    d*q*D_q=-DeltaS.                                   (L19)

The real-part inequality follows from V_z near i and
sqrt(zeta) near i*rhoA; its leading separation is real nonzero.
All statements concern actual stationary points, not the model positions
inserted into the actual exponent.

## 7. Lower-fold endpoint capture and the growing-field path

At a sequence approaching the fold, (L11) already puts BOTH endpoints
in the same fixed actual fold disc. The omega path is an exact straight
segment as a set, because it is affine in the positive real variable
lambda_theta. Equation(L10) shows its oriented tangent tends to the
negative real direction. Since zeta_omega->-i, its zeta increment lies
in a fixed narrow cone about the POSITIVE imaginary direction:

    Delta zeta=zeta_n-zeta_N,
    Im Delta zeta>=c|Delta omega|>=0,
    |Re Delta zeta|=o(1)*Im Delta zeta.                 (L20)

The estimates follow by integrating the uniformly continuous derivative
on the shrinking straight segment. They hold before assigning the
lower canonical angular label.

Suppose K_N->infinity. The upper root lies in(L16). A vector in a small
cone about the negative real ray plus a vector in a small cone about
the positive imaginary ray cannot cancel its norm. Thus

    |zeta_n|>=c|zeta_N|, K_n>=c*K_N->infinity.           (L21)

This is prior to the lower cone assignment. The lower root now has(L16).
Both endpoints are negative-cone points, whereas their difference has
the almost-vertical orientation(L20). It follows that their moduli are
comparable and |Delta zeta|<=C*epsilon*|zeta_N|. Invertibility of zeta
gives |Delta omega|<=C*epsilon*|zeta_N|.

The whole zeta path remains in the wider negative cone. Indeed a straight
chord between the endpoints lies in the narrower convex negative cone.
The image of the straight omega segment differs from that chord by at
most C|Delta omega|^2. Relative to |zeta_N| this is
O(epsilon^2*|zeta_N|)=o(1), by the uniform second derivative bound.
All moduli and K values along the actual adjacent path are consequently
comparable, with fixed angular/physical collars.

Set E=mD+b. At each actual integer zero Re E=0. Exact stationarity gives

    E_theta=-D-DeltaS-b_m+omega_theta*b_omega.           (L22)

By(L10),(L17)--(L19), the transported amplitude terms divided by the main
gamma*rhoA obey

    C[1/(M*rhoA^2)+1/(M*K*rhoA^3)]
                 =C[rhoA/K+1/K^2].                    (L23)

The m-derivative term contributes
C*tau^4*rhoA^2/(gamma*K^2), also harmless since gamma>=1/2.
At the upper zero |Re D|<=C(rhoA+1/K)/m. Its path variation is at most
C*gamma*rhoA/m. The total error/main ratio is therefore bounded by

    C[1/N+rhoA/K+1/K^2]=o(1).                           (L24)

Choose the packet threshold first and physical fold radius next, then
pass sufficiently far along the sequence. Equation(L22) has strictly
negative real derivative throughout the path, while both exact endpoint
zeros impose Re E=0. This contradiction excludes growing canonical fold
fields. No adjacent step was assumed small before both endpoint captures.

## 8. Bounded canonical fold fields

It remains that H_N=M_N^(2/3)*zeta_N is bounded. If H_n were unbounded,
then |zeta_n|>>M_N^(-2/3) whereas |zeta_N|=O(M_N^(-2/3)).
Thus zeta_n has the almost-positive-imaginary direction of(L20), in
contradiction with(L16) because K_n=|H_n|^(3/2)->infinity.
Hence BOTH endpoint H fields are bounded.

The complete local Airy class has the uniform two-term expansion

    J_M(H)=A0*[Ai(H)+M^(-1/3)*c1*Ai'(H)+O(M^(-2/3))],   (L25)

after removing its nonzero action/scale scalar. The coefficients A0,c1
depend on the fixed actual fold parameters, are uniformly bounded, and
A0 is uniformly nonzero. To derive this on the COMPLETE cycle, expand
the analytic inverse amplitude V_z at z=M^(-1/3)y,
zeta=M^(-2/3)H through its linear z term. The quadratic and zeta
remainders are bounded by C M^(-2/3)(1+|y|^2+|H|) on fixed Airy decay
arms. Their exponential cubic envelope is integrable uniformly for
bounded H. The fixed physical cuts have exp(-cM) loss by Section5.
This proves(L25), including derivatives on a larger fixed H-disc.
It is not a compact-local replacement for an unevaluated whole period.

The leading Ai has only negative real zeros, all simple: the standard
location is DLMF9.9(i), and simultaneous Ai,Ai' zero contradicts uniqueness
for the Airy ODE9.2.1. Uniform convergence and Rouche exhaust the finitely
many zeros in any bounded field window. Their scaled physical branches
therefore have the COMMON two-term form

    omega_j(M)=omega_f+alpha_j*M^(-2/3)
                             +beta_j*M^(-1)+O(M^(-4/3)), (L26)

uniformly in the actual compact fold parameters. The coefficients
alpha_j,beta_j are fixed for the support while its order changes.
This follows by the simple-root implicit expansion in(L25) and the
analytic inverse of zeta; a derivative of an unspecified remainder is
not used.

Both bounded endpoint zeros approach real Ai zeros. Equations(L10),(L20)
give

    Im(H_n-H_N)>=c*gamma*M_N^(-1/3)+O(1/N).              (L27)

The harmless O(1/N) comes from M_n/M_N=n/N in the bounded field.
As M_N<=N for small tau and gamma>=1/2, it is negligible compared with
gamma M_N^(-1/3). If that main scale failed to tend to0, both endpoints
could not have imaginary part tending to0. Therefore it tends to0.
The full field difference tends to0 as well, since the increment
orientation is almost vertical and all derivatives are bounded.
Uniform separation of the finitely many real model zeros now forces
the SAME root label j at both orders.

But the actual radial step satisfies the exact lower bound

    |omega_n-omega_N|
       =|beta*q_N|*(lambda-1)/tau^3 >=c*gamma/M_N,       (L28)

using lambda-1>=d/n and d|beta|=L(1-r).
On the other hand M_N-M_n=tau^4, so(L26) at the two endpoints gives

    |omega_j(M_n)-omega_j(M_N)|
         <=C*tau^4*M_N^(-5/3)+C*M_N^(-4/3)
         =o(gamma/M_N),                                (L29)

uniformly because gamma>=1/2 and M_N tends to infinity. This contradicts
(L28). The bounded canonical fold fields are excluded.

## 9. Exhaustion and the exact prospective consequence

The preceding cases exhaust a putative sequence: unbounded omega is the
small-T quartic regime; bounded omega away the lower simple fold has
strict complete Gaussian packets; near that fold growing and bounded
canonical fields are handled in Sections7--8. Conjugation covers lower q.
This proves the sequential large-thermal claim, subject to whole audit
of the actual contour construction and its uniform estimates above.

A useful quantified form is equivalent by a diagonal argument: there exist
UNEVALUATED delta0>0,U0,N0 such that for every ordinary integer support and
N>=N0 with

    0<|r|<delta0, sqrt(N)*|r|>U0,

the actual consecutive pair has no common root. If no such triple of
constants existed, choose counterexamples with |r|<1/j, sqrt(N)|r|>j
and N>=j, contradicting the sequential statement.

Combining this large-thermal result with the separately proved bounded-
thermal theorem at the ONE fixed U0 gives one unevaluated terminal cutoff
for ALL ordinary supports in a fixed relative near-critical strip

    |a-(d-1)^2/4|/[a+(d-1)/2]<delta0.                  (L30)

The exact equality curve is included through that bounded-thermal theorem.
This implication does not require proving all real temperatures at once,
and it does not give a numerical width or cutoff. Nor does it settle
supports outside that strip or full all-order binomial Strong Factorial.
An all-field coefficient version follows only after the actual rational
coprimality statement, by the unchanged factorial normalization and
rational Bezout descent.

## 10. Classical inputs and audit boundary

The actual current NIST pages were inspected:
https://dlmf.nist.gov/36.12.ii (36.12.11--12),
https://dlmf.nist.gov/2.4.v (2.4.17--21),
https://dlmf.nist.gov/9.7.ii (9.7.5--6),
https://dlmf.nist.gov/9.9.i, and https://dlmf.nist.gov/9.2.i .

They supply the classical cubic/Airy framework and model facts. They do
NOT certify the actual Gamma contour, exclude its high inactive saddle,
or automatically give the present uniform parameter bounds. Those
obligations are the explicit finite-horizontal and proper-cubic arguments
in Sections5--8. The analytic normal-form and proper-arm mechanisms are
the same ones written in SIMPLE_TURN_SHRINKING_WINDOW_CONSECUTIVE_EXCLUSION.md;
here the common model germ, finite contour gaps, compact parameter
domain, actual tails and lattice lower bound are reverified.

No finite experiment can establish these complete analytic estimates.
Exact symbolic checks can support(L6)--(L15),(L23)--(L24), but the
sequential theorem remains a written mathematical proof requiring whole
review. All work is scratch; no canonical, external-review, independent-
replay or formal-verification status is upgraded.

