# One support-uniform eventual cutoff for individual ordinary-binomial squarefreeness

For integers a>=0,d>=1,m>=1 define the ACTUAL monic polynomial

    A_(m;a,d)(x)=sum_(j=0)^m binom(m,j)*(am+1)_(dj)*x^(m-j).

**Theorem SupportUniformEventualOrdinaryBinomialSquarefreeness.**
There is an unevaluated integer M_star such that for EVERY ordinary
support a>=0,d>=1 and EVERY m>=M_star,

    gcd_Q[x](A_(m;a,d), partial_x A_(m;a,d))=1.           (Q1)

The same polynomial remains squarefree after extension to every
characteristic-zero field. The cutoff is independent of the support,
the coefficients and the field. No numerical M_star, all-order statement
or emptiness of the remaining repeated-root set is asserted.

This is a within-order proof. Consecutive coprimality, whether
support-fixed or support-uniform, is never substituted for (Q1).
The retained fixed-support INDIVIDUAL squarefreeness theorem is used only
for finitely many bounded supports. The proof uses the already proved
actual contour estimates as analytic lemmas, rechecking their one-order
hypotheses and holomorphic parameter collars below.

## 1. Exact normalization and the two derivative mechanisms

For a>=1 set D=a+d and

    I_m(w)=int_0^infinity [S^a(w+S^d)exp(-S)]^m dS
           =(am)!*A_(m;a,d)(m^d*w)/m^(Dm+1).            (Q2)

The scalar and affine variable change are nonzero. All other
normalizations used below multiply this actual period by a locally
holomorphic nonzero scalar and make a locally invertible analytic
change of variable. They preserve root multiplicity. For integer
a,d,m, finite paths may cross zeros of the integrand because the
integrand itself is entire. Estimating paths may vary nonanalytically;
the EXACT partial periods are defined by analytic moving endpoints
and are holomorphic for the same fixed support and fixed integer m.

Here are the two elementary mechanisms used throughout.

**Exact two-packet derivative.** Suppose on a complex parameter chart z
the original integer period is EXACTLY P_L+P_R, where

    P_i=e^(m Phi_i)*c_i*m^(-1/2)*(1+R_i),

the factors are nonzero and holomorphic, and |R_i|<1/2.
Locally put

    E=m(Phi_R-Phi_L)+log(c_R/c_L)
                         +log(1+R_R)-log(1+R_L).

At an actual zero, P_R=-P_L!=0 and

    partial_z(P_L+P_R)=P_R*E_z.                         (Q3)

Thus a nonzero E_z proves within-order simplicity. A support-dependent
common Gaussian width cancels in c_R/c_L; it does not enlarge its
parameter derivative. Local logarithms suffice.

Whenever a uniform bound |R_i|<=C/k is proved on a LARGER complex
z-chart of radius r_z, Cauchy gives

    |(R_i)_z|<=C'/(k*r_z)                              (Q4)

on a reserved smaller chart. This is a derivative of the exact
holomorphic remainder. We never differentiate an unspecified real-axis
O-symbol or change the support while taking this derivative.

**Compact simple-model limit.** Suppose normalized actual periods are
holomorphic in h on a fixed enlarged bounded h-domain and converge
there uniformly to a nonzero model Z(h,t), uniformly for a REAL
parameter t in a compact set. If every zero of each such Z is simple,
a sequence of ACTUAL double roots with bounded h is impossible.
Indeed pass to h_j->h_*,t_j->t_*; Cauchy's formula gives convergence
of the h-derivatives as well. The double-root equations would give
Z(h_*,t_*)=Z_h(h_*,t_*)=0. Equivalently, finite compact covers and
Rouche give exactly one root COUNTED WITH MULTIPLICITY near every
relevant model root. An enlarged collar handles roots on a boundary.
No all-index model separation theorem is needed.

The exact stationary envelope, at two distinct actual stationary
points S_L,S_R of the phase in (Q2), is

    d*w*partial_w(Phi_R-Phi_L)=-(S_R-S_L).               (Q5)

It follows directly by stationarity from
w Phi_w=(D-S)/d. With w=Wq, (Q5) is d*q*D_q=-DeltaS.
This identity, rather than consecutive-order transport, supplies the
leading within-order derivative.

## 2. A uniform high-ratio individual theorem is already available

The new scratch proofs UNIFORM_ACTUAL_BULK_ROOT_SIMPLICITY.md and
UNIFORM_ACTUAL_ENDPOINT_ROOT_SIMPLICITY.md prove actual root simplicity
for delta=d/sqrt(a)<=10^(-10), m>=10^16 in the central box and both
endpoint boxes respectively. Their union covers

    |Re z|<=10^(-7), |Im z|<2.1,
    x=m^d*a^d*(-1+delta*z).

The bulk proof replaces the nonholomorphic stronger-packet assignment
by exact analytic-midpoint periods and proves BOTH center connectors
have a fixed loss. Its Cauchy argument gives Re E_z>m/3.
The endpoint proof uses an OPEN common zeta collar and the complete
piecewise Airy vector remainder; its derivative has modulus greater
than .7*m^(2/3)*|Ai'|. Neither conclusion follows from coprimality.

Lemma UniformActualMomentRootRectangle in
UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md, (R6), captures EVERY
individual root whenever d>=2,m>=10^16,0<delta<=10^(-12).
Its Section2 explicitly imposes NO adjacent-order hard-range condition
on this within-order capture. The homotopy counts multiplicities and is
used solely for LOCATION. The two independent derivative proofs give

    a>=10^24*d^2, d>=2, m>=10^16
       ==> A_(m;a,d) is squarefree.                    (Q6)

This is the high-ratio input. It covers only the displayed large orders;
no soft coprimality theorem is used to claim small-order simplicity.

## 3. The small-ratio sector: conditional individual-root capture

Put c=a/(d-1)^2, J=d^2/(a+d), R=(a+d)/a for a>=1,d>=2.
We prove that some unevaluated epsilon_s>0,M_s works for every support
with c<=epsilon_s and every m>=M_s.

It suffices to exclude a sequence of repeated actual roots with
c->0,m->infinity. The RAW sparse theorem
lower_exponent_initial_and_boundary.tex, thm:sparse, states that
J>=4log(8m) makes the first m actual moment polynomials INDIVIDUALLY
SQUAREFREE as well as pairwise coprime. Thus a repeated root forces

    J<4log(8m).                                        (Q7)

As c->0, J>=1/(c+sqrt(c))->infinity. Also R<=4J for J>=1.
Consequently log R=O(log log m), and every fixed exponential loss
in m absorbs sqrt(Rm), uniformly under (Q7).

Use the exact normalization of
UNRESTRICTED_SMALL_RATIO_GAMMA_TWO_WELL_EXCLUSION.md, (U5):

    f_ell(S)=(S/ell)^ell exp(ell-S), D=a+d,
    W=D^D exp(-d)/a^a, S0=W^(1/d),
    I_m(Wq)=R0^m*J_m(q),
    J_m(q)=int_0^infinity (q f_a+f_D)^m dS,
    R0=D^D exp(-D).

Its complete one-packet constructions (U9)--(U17) are within-order
statements. For |q|<=1/2 they are uniformly zero-free. For |q|>=2,
divide by q and use the left Gamma well. Its only extra nonlocal factor
is sqrt(m)(1+sqrt(R))exp(-delta_0*m), which (Q7) absorbs.
Therefore EVERY repeated root along the proposed sequence has

    1/2<|q|<2,     q=x/(m^d W).                         (Q8)

This is conditional capture for the one repeated-root order, not an
unconditional all-order lower annulus. Applying thm:sparse to an
individual repeated root is essential.

On the LARGER annulus 1/3<=|q|<=3, split the exact integer period at the
SUPPORT-FIXED identric seam S0. Source (U18)--(U20) proves

    P_L=sqrt(a)*lambda_L^m*sqrt(2pi/m)c_L*(1+R_L),
    P_R=sqrt(D)*lambda_R^m*sqrt(2pi/m)c_R*(1+R_R),
    |R_i|<=C/m.                                       (Q9)

These are exact holomorphic q-functions. The source proves the complete
left cross-L1 estimate in its OWN width, so (Q9) is uniform even when
R is unbounded or a stays bounded. Its compact analytic Gamma family
H_b, b=1/sqrt(a) in[0,1], includes such bounded a.
No Gaussian-limit assumption a->infinity is made.

Use q-discs of radius1/12 inside the larger annulus, with finitely
many local log charts. Cauchy gives (R_i)_q=O(1/m) at (Q8).
The actual Morse amplitudes and their q-derivatives are uniformly
bounded, and c_i stay away from0. The support-fixed sqrt(R) in their
ratio has zero q-derivative. Finally source (U15) and (Q5) give

    Re DeltaS/d->1,  D_q=-DeltaS/(d q),
    |D_q|>=1/4                                       (Q10)

eventually, uniformly on (Q8). Thus the exact ratio in (Q3) has

    E_q=mD_q+O(1),       |E_q|>=m/4-C>0.

The actual root is simple, contradicting its choice.
If no epsilon_s,M_s existed, choose repeated roots with c<=1/j,m>=j
and obtain the excluded sequence. This proves the small-ratio
individual theorem. No adjacent endpoint, d/m-smallness, model-power
replacement, or all-support coprimality theorem entered the argument.

## 4. Compact subcritical ratio intervals

Fix 0<c_-<=c_+<1/4. A sequence of repeated roots with m->infinity
and c in this interval has either bounded d, leaving finitely many
integer supports, or d->infinity after a subsequence. The bounded
case is excluded ONLY by the retained408 fixed-support INDIVIDUAL
squarefreeness theorem and a finite maximum of its cutoffs.

For unbounded d use the ACTUAL compactification (S1)--(S8) of
COMPACT_SUBCRITICAL_RATIO_UNIFORM_EVENTUAL_SFC.md. Its nonzero scalars
preserve multiplicity; its original integer integrand is entire on
finite deformations. The whole original real-tail derivative bounds
and the large fixed vertical cuts hold on fixed complex q-discs.
The old Gamma UPPER disk contains all normalized roots. The q=0
Gaussian chart is zero-free. The favored simple whole-horizontal
maximum excludes compact sets off the unit circle. Thus any actual
root sequence has a subsequential limit on |q|=1, without a lower
annulus or a circle-law transfer.

For beta=1/(4c)>1 every point of that circle, including q=1 and-1,
has two strict nondegenerate global maxima on one horizontal and
a strict middle gap. On a small fixed COMPLEX q-chart the exact
analytic-midpoint periods have the complete source estimate (S15),

    P_i=lambda_i^m*m^(-1/2)c_i*(1+R_i), |R_i|<=C/m.

Actual critical values stay in their exponents; no m/d restriction is
used. On a reserved smaller fixed chart, Cauchy yields
(R_i)_q=O(1/m), and partial_q log(c_R/c_L)=O(1).
By (Q5) and source (S14),

    D_q=-L*(t_R-t_L)/(d q), Re(t_R-t_L)>=b0>0,
    L/d ->2c.

Consequently |D_q| is uniformly bounded below by a positive constant
on every member of a finite cover of the compact circle/c-parameter
set. The exact derivative E_q=mD_q+O(1) cannot vanish at large m.
Equation (Q3) excludes the repeated root. The complete one-packet
neighborhoods, two-packet charts and finite bounded-support fallback
prove one individual-squarefreeness cutoff on the fixed interval.


## 5. One-order capture near the critical support curve

Put

    s=(d-1)/2, Sc=s(s+1), e=a-s^2, C=a+s, L=sqrt(C),
    r=e/C, beta=-Sc/(dL), kappa=(s+1)/s, w0=C^d*kappa,
    q=x/(m^d*w0), eta=beta*(q-1).

Here s may be a half-integer. The exact actual normalized factor is

    Psi(t,q)=p^a exp(-Lt)(kappa*q+p^d), p=1+t/L,
    I_m(w0q)=L[C^(a+d)exp(-C)]^m
                         int_(-L)^infinity Psi(t,q)^m dt.       (Q11)

All scalars are nonzero. Along r->0,d->infinity the factor and every
fixed analytic derivative on finite rectangles converge to

    Psi_0(t,q)=exp(-t^2/2)*(q exp(-t)+exp(t)).                    (Q12)

The exact real-tail derivative, for positive q and
H=p^d/(kappa*q+p^d) in[0,1], is

    partial_t log Psi=(-s-Lt+dH)/(Lp).                           (Q13)

As in the actual bounded-thermal proof, t>=4 has derivative at most-2/5
and t<=-3 has derivative at least1, eventually; the actual left endpoint
t=-L is included. The inequality |Psi(t,q)|<=Psi(t,|q|) controls complex
q. First choose a large fixed real cut, then use (Q12) on the finite
rectangle to control both vertical cuts. Thus complete original tails,
including the long left interval, have fixed uniform loss relative to
selected compact positive critical moduli. No Taylor approximation is
used on unbounded original tails.

For clarity a sharp lower annulus is NOT needed for the one-order
argument. The OLD Gamma upper bound gives

    |q|<=((a+d)/C)^d/kappa ->exp(2)                              (Q14)

when r->0,d->infinity. At q=0 the limiting factor is the Gaussian
exp(t-t^2/2), with unique maximum t=1 and curvature-1. Its strict
finite-horizontal gaps, the complete tails (Q13), and actual analytic
perturbation prove a fixed zero-free complex q-disc about0. Off the
unit circle the favored single whole-horizontal maximum of (Q12)
gives the same complete one-packet zero-free charts. Finite covers of
the bounded upper disk minus that q0 chart therefore force every
subsequential actual root limit onto |q|=1.

On every closed circle arc avoiding1 there are two strict simple
maxima on one horizontal with a strict middle loss. Exact analytic-
midpoint pieces, the entire integer integrand, and (Q13) give complete
actual Gaussian packets on a fixed complex q-chart, with R_i=O(1/m).
Their q-Cauchy bounds and (Q5), with L/d->1/2 and the nonzero saddle
separation, give E_q=mD_q+O(1)!=0 at large m, exactly as in Section4.
A repeated-root sequence approaching the critical ratio must consequently
satisfy

    q->1.                                                       (Q15)

This is a SINGLE-order localization and derivative statement.
No second actual endpoint, radial order step, d/m-smallness or
consecutive noncollision is used. A real q near1 is positive and
cannot be a root; conjugation lets us take Im q>0 and Im eta<0.

If d is bounded in a sequence with r->0, a is also bounded:
for fixed s, r=(a-s^2)/(a+s)->1 as a->infinity.
Only finitely many supports remain, handled by408 individual
squarefreeness. Thus all new near-critical analysis may take d->infinity.
It allows e=0; integrality is not needed to remove the equality curve
in the bounded-thermal argument.

## 6. Bounded thermal parameter and a complete compact quartic limit

Fix finite U. Suppose a sequence of repeated roots has m->infinity and

    |u|<=U, u=-sqrt(m)*r.

After the finite bounded-support fallback, r->0,d->infinity, and
Sections5's actual contour capture applies. Hence q->1. Put

    h=m^(3/4)*eta.

When h is bounded, set t=m^(-1/4)v in (Q11), removing the nonvanishing
central factor. The exact raw jets in the bounded-thermal proof give

    mF(m^(-1/4)v,q)
          =-v^4/12+u*v^2/2+h*v+o(1),                           (Q16)

uniformly for bounded v,h and real |u|<=U. The complete-period
domination, not just a local Taylor series, is essential: on a fixed
physical real interval,

    Re F(t,q)<=-c*t^4+C|r|t^2+C|q-1||t|.

After scaling this gives the integrable envelope
exp(-c*v^4+C_U*v^2+C_H|v|). The finite outside interval has a strict
gap, preserved by small r,q-1; (Q13) controls the complete far tails.
The same argument is uniform on a LARGER complex h-disc, because
|h| remains bounded there. Therefore the normalized original period
is holomorphic in h and converges locally uniformly to

    Z(h,u)=int_R exp(-v^4/12+u*v^2/2+h*v)dv.                     (Q17)

The temperature u is REAL and is held fixed when h varies at this
one fixed support/order. Actual higher jets may have complex corrections;
no assertion of an exactly real canonical temperature is made.

Only Sections1--2 of QUARTIC_LEE_YANG_ZERO_AND_ROBIN_BRIDGE.md are used:
Newman's actual quartic-measure theorem places all Z zeros on the
imaginary axis for EVERY real u, and the exact ODE plus heat identity
proves all those zeros simple. The proof of simplicity uses

    F_xxx+3u F_x-3xF=0,  F_u=-F_xx/2, F(x,u)=Z(ix,u).

A zero cannot have multiplicity at least3 by ODE uniqueness. A double
zero would split off the real x-axis under a small negative real heat
step, contradicting Newman's location theorem. No HJN zero-motion,
monotonicity, quotient bound, or all-index Robin estimate is imported.

Thus Section1's compact simple-model lemma excludes a repeated actual
root with bounded h, uniformly for real u in[-U,U]. There are finitely
many model zeros in each enlarged bounded h-domain; alternatively a
hypothetical actual double root would converge, with its derivative,
to Z=Z_h=0. The nonzero affine h_q=m^(3/4)*beta and the nonvanishing
normalization preserve multiplicity.

## 7. A general small-residual-temperature growing-quartic lemma

The remaining bounded-thermal case has |h|->infinity. More generally
consider ANY actual near-critical sequence with q->1 and

    rho=|eta|^(1/3)->0, k=m*rho^4->infinity,
    T=-r/rho^2->0.                                             (Q18)

Use t=rho*v. The complete ACTUAL phase family, compact through
rho=T=nu=1/d=0, has limit

    -v^4/12+T*v^2/2+omega*v, omega=eta/rho^3.

The complete intermediate real envelope is
-c*v^4+C|T|v^2+C|v|. Whole far tails retain the fixed physical loss
from Section5; after choosing finite scaled cuts, compact perturbation
transfers strict model horizontal and midpoint gaps. These are precisely
the revalidated one-order constructions (T16)--(T19), not their
subsequent adjacent-order comparison.

Outside the central negative-imaginary eta cone there is one strict
whole-horizontal maximum, so no zero. Inside a slightly wider central
cone the exact analytic-midpoint pieces have

    P_i=e^(mPhi_i)*sqrt(2pi/m)c_i*(1+R_i),
    |R_i|<=C/k,
    |(R_i)_q|<=C/(k*rho^3),                            (Q19)
    |partial_q log(c_R/c_L)|
             <=C[rho^(-2)+|T|rho^(-3)].

The Cauchy radius is c*rho^3 in the complex q-plane. During this
variation the support r,beta and the CENTER rho are fixed; hence
T=-r/rho^2 is fixed too. The T/rho^3 term is indispensable.
The uniform value estimate holds on a larger such chart, so the
derivative in (Q19) is legitimate. Entire integer paths, strict
midpoint losses and whole real tails define actual holomorphic
periods, not only formal local saddle integrals.

The actual stationary points satisfy |DeltaS|>=c*L*rho.
Since L/d->1/2 and q->1, (Q5) gives

    |mD_q|>=c*m*rho.

At an actual zero the exact derivative of the ratio consequently obeys

    |E_q-mD_q|/|mD_q|
       <=C[rho/k+|T|/k+1/k^2]->0.                              (Q20)

Thus E_q!=0 and (Q3) proves simplicity. The ACTUAL critical actions
remain in their exponentials; no condition m*rho, m/d or
m times a phase-perturbation error is imposed.

For the bounded-thermal sequence, k=|h|^(4/3)->infinity and
T=u/sqrt(k)->0, so (Q18)--(Q20) apply. Together with Section6 and
the compact regular-circle analysis, this excludes EVERY repeated-root
sequence with m->infinity, |sqrt(m)r|<=U. Therefore one unevaluated
within-order cutoff M_b(U) works uniformly throughout each such tube,
including the equality curve r=0.

## 8. Large thermal parameters: actual thermal models and regular fields

We next exclude every repeated-root sequence with

    m->infinity, 0<|r|->0, sqrt(m)|r|->infinity.                (Q21)

Only the ANALYTIC whole-contour modules of
LARGE_THERMAL_QUARTIC_FOLD_MATCHING.md are used. Its adjacent endpoint
relations and gamma/M->0 consequence are not premises here.

Since e=a-s^2 is quarter-integral and nonzero, |e|>=1/4.
Put

    tau=sqrt(|r|), sigma=-sign(r),
    gamma=L*tau=sqrt(|e|)>=1/2,
    M=m*tau^4=(sqrt(m)|r|)^2->infinity,
    omega=eta/tau^3.                                            (Q22)

The fixed actual support determines tau,sigma,beta for the present
omega-derivatives. After t=tau*v, the actual scaled phase is a compact
holomorphic perturbation of

    f_sigma(v,omega)=-v^4/12+sigma*v^2/2+omega*v,
    sigma in{+1,-1}.                                           (Q23)

All leading exponentials are actual. The complete REAL intermediate
range satisfies -c*v^4+C*v^2+C|omega||v|. The outer original tails
retain a fixed physical gap by (Q13). Finite vertical cuts at large
fixed scaled v have the quartic loss. For bounded omega these give
complete relative O(1/M) Gaussian errors uniformly, not a Taylor
series on an unbounded contour. In comparing the far original pieces
to the scaled leading width, the common L cancels and the remaining
factor tau^(-1)*sqrt(M)=sqrt(m)*tau is at most sqrt(m).
Shrinking tau first keeps the bounded actual scaled critical heights
from consuming the fixed physical gap. Thus the far pieces are
exponentially small uniformly; no gamma/M assumption is hidden here.

If |omega|->infinity, the field blowup (Q18) has
k=M|omega|^(4/3)->infinity and
T=sigma/|omega|^(2/3)->0. Section7 excludes a repeated root.
It remains to pass to bounded omega and a subsequential limit.

The complete horizontal constructions of the two thermal models give:
off the imaginary omega-axis one favored unique simple global maximum,
hence zero-free actual neighborhoods. On omega=-iz,z>=0 the active pair is

    v_R=x-iy, v_L=-x-iy,
    x^2=3(y^2+sigma), z=(8/3)y^3+2sigma*y.

For sigma=+ both simple peaks exist for every z>=0, including z=0;
the midpoint loss is x^4/12>0. For sigma=-, 0<=z<2/3 is a zero-free
single-peak regime, z>2/3 has the two strict peaks, and the only
remaining original-cycle fold is omega_f=-2i/3,v_f=-i.
The other high saddle is inactive, as addressed in Section9.

At any regular two-peak limiting point choose a fixed complex
omega-disc and exact analytic-midpoint pieces. The complete actual
contours give R_i=O(1/M) on a larger disc; Cauchy gives
(R_i)_omega=O(1/M), and the common singular width cancels from the
ratio amplitude, whose logarithmic omega-derivative is O(1).
If Delta_v denotes their actual scaled separation, exact (Q5) yields

    D_omega=tau^4*Delta_v/[(1-r)q],
    mD_omega=M*Delta_v/[(1-r)q].                               (Q24)

Here d*beta=-L(1-r) and DeltaS=L*tau*Delta_v.
The scaled separation is bounded away from0 on this compact regular
chart. Therefore E_omega=mD_omega+O(1) is nonzero for large M.
All regular bounded omega repeated roots are excluded.

## 9. The actual lower thermal fold: complete Airy simplicity

Near sigma=-1,omega_f=-2i/3, the model's finite horizontal has the EXACT
identity

    f_-(x-i,-2i/3)=-1/4+i*x^3/3-x^4/12.                        (Q25)

Its outside pieces have strict loss; the high third saddle v=2i,
with action2, is NOT on this original contour. The original real tails
and fixed finite vertical cuts were controlled in Section8.
Choose these cuts and a fixed local interval first; compact analytic
perturbation preserves their losses for the ACTUAL fold.

The uniform analytic cubic map is

    f=A(omega)+z^3/3-zeta(omega)z,
    v=V(z,omega), V_z near i, zeta_omega near-i.

The model horizontal has z=-ix+x^2/12+O(x^3), so its two actual
endpoints lie in the RIGHT two cubic decay sectors. They select one
nonzero Ai cycle, not a Scorer cycle. Whole nonlocal pieces are strictly
below it. These facts and the inactive third saddle are the explicit
full-period statements (L13)--(L17), not consequences of a model
zero-location theorem.

Let rho_A=|zeta|^(1/2), K=M*rho_A^3.
For K->infinity outside a fixed negative-real zeta cone, the proper
cubic inverse-arm and sublevel-detour constructions give one nonzero
complete Gaussian packet. Inside a wider negative cone they split
the ORIGINAL period at V(-sqrt(-zeta),omega), using the shared left
cubic sheet. The two pieces are actual holomorphic functions, with

    |R_i|<=C/K, |(R_i)_omega|<=C/(K*rho_A^2),
    |partial_omega log(c_R/c_L)|<=C/rho_A.                      (Q26)

These are Cauchy bounds on SUPPORT-FIXED complex omega-collars of
radius c*rho_A^2, with center rho_A fixed. The fixed physical exit
loss and proper arms control all nonlocal paths. No bounded-real-phase
compactness assertion about the cubic is used.

The actual scaled separation is at least c*rho_A. By (Q24),

    |mD_omega|>=c*M*rho_A,
    |E_omega-mD_omega|/|mD_omega|
          <=C[rho_A/K+1/K^2]->0.                              (Q27)

Thus growing canonical fields have only simple actual roots.

If K is bounded, H=M^(2/3)zeta is bounded. After removing its analytic
nonzero action/scale, the COMPLETE original period has the uniform
expansion

    A0*[Ai(H)+M^(-1/3)c1*Ai'(H)+O(M^(-2/3))],                  (Q28)

on an ENLARGED bounded H-domain, with A0 uniformly nonzero.
The inverse amplitude and original full contours, including the fixed
outer loss in (Q25), prove this estimate and its holomorphy in H.
The coordinate derivative zeta_omega is uniformly nonzero.
Every Ai zero is simple by its ODE, so Section1's compact simple-model
argument excludes a repeated actual root. No comparison of root labels
at two orders is needed.

This exhausts (Q21). By the diagonal argument, some unevaluated
delta_t>0,U_t,M_t exclude repeated roots whenever

    m>=M_t, 0<|r|<delta_t, sqrt(m)|r|>U_t.                       (Q29)

Otherwise choose counterexamples with |r|<1/j, sqrt(m)|r|>j,m>=j.
Together with Section7's bounded-thermal theorem at the ONE fixed U_t,
this proves one within-order simplicity cutoff on the entire fixed
relative critical strip |r|<delta_t, including r=0.


## 10. Compact supercritical ratio intervals

Fix 1/4<c_-<=c_+<infinity. Bounded d leaves finitely many supports.
Their fallback is408 INDIVIDUAL squarefreeness, not407 or the
coprimality conclusion of the compact-supercritical theorem.
Consider a repeated-root sequence with d,m->infinity and c in this
fixed interval.

Use the actual normalization and whole tails (A1)--(A7) of
COMPACT_SUPERCRITICAL_RATIO_UNIFORM_EVENTUAL_SFC.md.
Here C=a+(d-1)/2, L=2C/d, w0=C^d*(d+1)/(d-1), and

    Psi_(d,c)(t,q)
       ->exp[-t^2/(2beta)](q exp(-t)+exp(t)),
    beta=1/(4c) in a compact subset of(0,1).

The old Gamma upper root disk is bounded uniformly on this interval.
The q0 Gaussian chart and every off-circle favored whole-horizontal
chart are zero-free, using their actual critical values and complete
original tails. Thus a subsequential root limit is on |q|=1.
Selected-saddle existence uses precisely Kabluchko Theorem2.3 with
the retained substitution, not a circle-law or Theorem2.4 transfer.

On the circle gap around1 the stable saddle is a unique simple global
maximum, hence zero-free. On each closed two-peak arc away from its
edges, the exact integer analytic-midpoint pieces have holomorphic
complete errors R_i=O(1/m) on fixed complex q-charts. Their Cauchy
derivatives are O(1/m), and the amplitude-ratio derivative is O(1).
The exact identity

    D_q=-L*(t_R-t_L)/(d q)

has modulus bounded below by a positive constant on each compact
two-peak chart. Consequently E_q=mD_q+O(1) is nonzero. Repeated roots
can only approach an actual edge.

The edge hypotheses must refer to the actual cycle. They are supplied
by (A14)--(A20): on the model fold horizontal there is a unique
quartic real-part maximum with strict fixed outside loss; finite
vertical connectors and the whole original real tails share that loss.
The actual nearby folds admit a uniform analytic cubic inverse.
The positive quadratic bend of the horizontal image gives exactly the
RIGHT pair of cubic exits and one nonzero Ai cycle, not a Scorer cycle.
The explicit pole/zero-free strip contains the finite deformations;
the actual paths and their q-logarithms persist on genuine complex
parameter collars. These are one-period assertions, not consequences
of adjacent noncollision.

Write the actual canonical coordinate zeta=zeta(q),
rho=|zeta|^(1/2), K=m*rho^3. For K->infinity outside the negative cone,
proper cubic arms and their fixed-sublevel detours give a nonzero
single packet. In the larger negative cone the actual shared-left-sheet
midpoint split gives

    |R_i|<=C/K, |(R_i)_q|<=C/(K*rho^2),
    |partial_q log(c_R/c_L)|<=C/rho.

The support-fixed q-collar has radius c*rho^2. The actual separation
is at least c*L*rho and L/d is bounded above and below on the ratio
interval. Therefore

    |mD_q|>=c*m*rho,
    |E_q-mD_q|/|mD_q|<=C[rho/K+1/K^2]->0.                      (Q30)

Equation (Q3) proves growing-field simplicity.

For bounded K, H=m^(2/3)zeta is bounded. The COMPLETE period expansion
(A26), on a larger H-disc, is

    A0*[Ai(H)+m^(-1/3)c1*Ai'(H)+O(m^(-2/3))],

with uniformly nonzero A0. Its physical fixed-cut and original-tail
losses were established before the local amplitude expansion.
Its normalized period is holomorphic in H and zeta_q is uniformly
nonzero. The compact simple-model argument excludes double roots,
using only the simple zeros of Ai. There is no adjacent root label,
radial transversality or second-root-capture premise.

Thus every unbounded-gap repeated-root sequence in the fixed compact
supercritical interval is excluded. The finite bounded-support
fallback supplies one cutoff for the whole interval.

## 11. Boundary supports without heavy finite certificates

For d=1 and every a>=0, put gamma=am+1. The exact polynomial
G_m(x)=A_(m;a,1)(x) satisfies

    xG_m''+(1-gamma-m-x)G_m'+mG_m=0.

Its constant coefficient (gamma)_m is nonzero. If x0 were a root of
multiplicity at least2, x0!=0 and the lowest vanishing-order term
x0*G_m'' could not cancel with either other term. Hence every order
has simple roots. This is an exact ODE criterion, not an inference
from adjacent Padé coprimality.

For a=0, A_m has the EXACT Appell identity

    A_m'=m*A_(m-1).                                            (Q31)

The following RAW proofs give one uniform cutoff m>=10000:
for d=2,3,4,5 use PRIME_POWER_TWO_TERM_SFC.tex, thm:sfc at shape1,
together with (Q31); for 6<=d<=sqrt(m), use
ALL_UNIVARIATE_TWO_TERM_SFC.tex, thm:analytic, whose conclusion is
individual simplicity. For d>sqrt(m), put t=sqrt(m)>=100.
The elementary function 2^t/(4t^4) is increasing there, since
log2-4/t>=1/2-4/100>0, and 2^100>4*100^4. Therefore

    2^d>=4m^2.

The LITERAL conclusion of RAW thm:large in
ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex is that

    F_(m,d)(t)=sum_(j=0)^m binom(m,j)(dj)!*t^j

is squarefree. Exactly A_m(x)=x^m*F_(m,d)(1/x), and both constant and
leading coefficients are nonzero. Reversal preserves multiplicities,
so A_m is squarefree DIRECTLY. No Appell round trip through adjacent
coprimality is needed for this large-gap criterion. Appell (Q31)
is used for the small prime-power thm:sfc transfer above.

The raw exponent-six analytic entropy source may be authenticated for
the credit/crosslink in the analytic proof; its finite prefix is not
a premise. No 92935-pair,9998-order, lower-one, relative-gap or other
heavy boundary certificate is required. The general invalid inference
from consecutive coprimality to individual simplicity is not used.

## 12. A finite cover of ALL positive supports

Sections3,4,7--10 prove the following ONE-ORDER regional statements:

1. Some epsilon_s>0,M_s give simplicity for all c<=epsilon_s,m>=M_s.
2. Every fixed compact interval inside(0,1/4) has a uniform cutoff.
3. Some delta_t>0 and a uniform cutoff give simplicity for |r|<delta_t.
4. Every fixed compact interval inside(1/4,infinity) has a uniform cutoff.
5. Region (Q6) has the explicit individual cutoff m>=10^16.

There is no hidden sequence-exclusion shortcut in the critical region:
Section9's large-thermal contradiction gives a fixed strip and ONE
fixed U_t by the stated diagonal argument; Section7 is then used at
that finite U_t. It is not necessary to take a maximum over all real
temperatures or all bounded-field windows.

For a>=1,d>=2 the exact identity is

    r=(c-1/4)/(c+1/[2(d-1)]).

Choose

    epsilon=min(epsilon_s,1/16),
    h=min(1/8,delta_t/16), K=4*10^24.

If 1/4-h<=c<=1/4+h, the denominator in the r-identity is at least1/8.
Thus |r|<=8h<=delta_t/2<delta_t. If c>=K, then

    a=c(d-1)^2>=10^24*d^2,

because 4(d-1)^2-d^2=(d-2)(3d-2)>=0 for d>=2.
Consequently the FIVE closed regimes

    (0,epsilon],
    [epsilon,1/4-h],
    [1/4-h,1/4+h],
    [1/4+h,K],
    [K,infinity)

are covered by the five individual-simplicity inputs. The middle
two compact intervals stay strictly away from their excluded endpoints.
Take the maximum of their finitely many unevaluated cutoffs, the
small/critical cutoffs,10^16 and10000. Boundary Section11 covers a=0
or d=1; every other support has c>0. This proves (Q1) with one M_star.

For completeness the sequential arguments in each region give the
claimed uniform quantifiers: failure of a uniform cutoff supplies
repeated-root orders m_j>=j. Compact ratios admit parameter subsequences;
bounded d leaves only finitely many supports and is handled by408.
The small and critical strips use their explicit diagonal sequences.
Every possible such sequence is excluded by the within-order derivative
or compact-simple-model arguments, so no estimate is inferred backwards.

## 13. Fields, coefficient phases and exact remaining scope

For every support and m>=M_star, no complex repeated root means that
the rational polynomials A_m,A_m' have a rational Bezout identity.
It remains valid over every characteristic-zero field. No embedding
of an arbitrary such field into C is assumed.

The actual homogeneous moment form is

    Lambda((X^a(A+B X^d))^m)=B^m*(am)!*A_m(A/B)

on the B!=0 coefficient chart. At B=0 its value is (am)!*A^m,
nonzero for A!=0; at A=0,B!=0 it is ((a+d)m)!*B^m.
Thus every projective zero of this moment form is simple in its
coefficient-ratio coordinate, uniformly at all orders m>=M_star.
The statement is about individual moment forms; it is not an
all-order Strong Factorial assertion or a new inference about
consecutive moments.

The proof does NOT claim a numerical M_star, individual simplicity
at every smaller order, or emptiness of the remaining bad set.
Finiteness of the GLOBAL repeated-root triple set additionally
requires a separate fixed-terminal finite-support argument.
That finiteness is not inferred from a soft coprimality bound here.

## 14. Exact source boundary for the one-order adaptations

The logical inputs are:

- UNIFORM_ACTUAL_BULK_ROOT_SIMPLICITY.md and
  UNIFORM_ACTUAL_ENDPOINT_ROOT_SIMPLICITY.md, plus401 (R6) for location
  only, supplying the individual high-ratio theorem.
- UNRESTRICTED_SMALL_RATIO_GAMMA_TWO_WELL_EXCLUSION.md (U5)--(U20):
  exact seam periods, own-width L1 tails, compact Gamma Morse family
  and uniform COMPLEX-q error estimates. RAW thm:sparse supplies its
  INDIVIDUAL simplicity alternative for the width condition.
- COMPACT_SUBCRITICAL_RATIO_UNIFORM_EVENTUAL_SFC.md (S1)--(S16):
  actual compactification, whole tails, q0 chart, circle geometry and
  holomorphic midpoint pieces. Its adjacent-order conclusion is unused.
- UNIFORM_BOUNDED_THERMAL_TRANSITION_EVENTUAL_SFC.md (T3)--(T19):
  actual detuned whole contours, compact raw quartic limit and
  support-fixed singular q-collars; the T/rho^3 amplitude term is kept.
- LARGE_THERMAL_QUARTIC_FOLD_MATCHING.md (L1)--(L9),(L13)--(L19),(L25):
  actual thermal compactification, complete finite horizontals,
  inactive third saddle, proper cubic paths/shared sheet, full
  holomorphic errors and nonzero Ai amplitude. No adjacent-path
  assumption gamma/M->0 is imported.
- COMPACT_SUPERCRITICAL_RATIO_UNIFORM_EVENTUAL_SFC.md (A1)--(A7),
  (A9)--(A21),(A26): its complete actual root localization, regular
  charts, selected right Ai class, zero-free log lifts and full
  cubic collars. Bounded-support fallback is changed to408.
- QUARTIC_LEE_YANG_ZERO_AND_ROBIN_BRIDGE.md Sections1--2 ONLY:
  Newman's all-real-temperature model location and the independent
  ODE/heat simplicity proof. HJN and the later Robin/velocity bounds
  are not premises.
- EVENTUAL_SQUAREFREENESS_OF_ORDINARY_BINOMIAL_MOMENTS.md (408):
  used only for finitely many bounded integer supports.
- GAMMA_PRODUCT_ROOT_MODULUS.tex, thm:product-annulus, for the OLD
  upper root disks. No new lower annulus is required by this proof.
- The RAW constant-family analytic/large-exponent/prime-power proofs
  and the elementary Gamma ODE for the boundary supports.

In every use of a regional account the new argument names the ACTUAL
whole-period estimate and its open parameter collar, then applies a
within-order derivative or a compact simple-root count. None of the
regional coprimality theorem statements,407 or418, is a logical
substitute for that work. All leading exponentials are actual critical
values. No finite test establishes a universal analytic estimate,
and no new canonical status is assigned by this scratch synthesis.


## 15. Separate finite-exception theorem and genuine six-way reverse

This section adds the independently proved theorem
FIXED_TERMINAL_REPEATED_ROOT_SUPPORT_FINITENESS.md (847):
for each FIXED m>=2 there are only finitely many ordinary supports
whose actual A_m has a repeated root. It is NOT an input to the
five-region uniform cutoff proved in Sections1--13.

Define

    R={(a,d,m): a>=0,d>=1,m>=2,
                      gcd_Q(A_(m;a,d),A'_(m;a,d))!=1},
    S={(a,d): some (a,d,m) belongs to R}.

Let H_(m;a,d)(A,B)=Lambda((X^a(A+B X^d))^m) be the degree-m
HOMOGENEOUS MOMENT FORM IN THE COEFFICIENTS A,B. Define P to be
the intrinsic projective phases (a,d,[A:B]) over C such that, for
some m>=2,

    H_m(A,B)=partial_A H_m(A,B)=partial_B H_m(A,B)=0.             (Q32)

The support a,d is part of the phase. This singularity concerns
the COEFFICIENT projective moment form, not roots of the binomial
as a polynomial in X.

Every point in P has AB!=0: the exact values at either coefficient
axis were shown in Section13 to be nonzero. On the B!=0 chart,

    H_m=B^m*(am)!*A_m(x), x=A/B,
    partial_A H_m=B^(m-1)*(am)!*A_m'(x).

Consequently (Q32) is equivalent to x being a root of
gcd_Q(A_m,A_m'). The reverse uses Euler's identity
A*H_A+B*H_B=mH to obtain H_B=0. Such x is nonzero and algebraic
over Q. Each bad triple contributes finitely many such points,
at most floor(m/2) distinct repeated roots.

**Corollary UniformSquarefreenessFiniteExceptionEquivalence.**
Under the written fixed-terminal finiteness theorem847 and the
retained408 fixed-support INDIVIDUAL eventual squarefreeness theorem,
the following six assertions are equivalent:

1. There is one M_star such that every rational A_(m;a,d) is squarefree
   for all ordinary supports and all m>=M_star.
2. The global repeated-root triple set R is finite.
3. There is one cutoff with gcd_K(A_m,A_m')=1 for every ordinary
   support and EVERY characteristic-zero field K.
4. There is one cutoff after which no complex projective moment
   form has a singular coefficient phase (Q32).
5. The intrinsic singular-projective-phase set P is finite.
6. The bad-support projection S is finite.

Here are all directions and the exact reverse premises.

For1=>2, theorem847 makes each fixed-order fiber finite. All bad
triples lie at m<M_star, so a finite union proves R finite.
THIS implication uses847, and does not use408 or a soft pair theorem.
For2=>1 choose2 when R is empty, otherwise one plus its largest order.

For1=>3 use the rational Bezout identity for A_m,A_m'. The assertion
persists over every characteristic-zero field, including all its
extensions. For3=>4 restrict to C and use the coefficient-chart
identities above. For4=>1 a nontrivial rational gcd has a complex
root x, producing the forbidden singular phase [x:1].
No converse for only one prescribed nonclosed field is substituted.

For2=>5 take the finite union of the finitely many repeated roots
at each member of R. For5=>6 project its INTRINSIC support; every
bad support supplies at least one such complex algebraic phase.
For6=>2, apply408 to each of the finitely many supports in S.
Each has only finitely many bad orders, so their union is finite.
WITHIN THESE reverse implications,408 is used only for the
support/projective-projection-to-finite-triples direction.
It is not needed for1=>2. The uses of408 for finitely many bounded
supports in the analytic proof were already specified separately.

The theorem in this note independently proves1;847 now supplies2
and hence all six assertions. In particular the global repeated-root
triple set and its support and intrinsic projective phase projections
are finite. This does NOT prove any of them empty, compute a
containing box, or evaluate M_star. Finiteness is projective:
arbitrary nonzero coefficient scaling can still produce infinitely
many unscaled binomials from one exceptional phase.

The contexts in the reverses are essential. In an abstract pattern,
infinitely many supports bad only at order2 can have a uniform tail;
without fixed-terminal finiteness this need not give finite R.
One support bad at infinitely many orders has a finite support
projection; without408 that need not give finite R.
These are logical countermodels, not actual factorial examples.
No estimate is inferred backwards and no bad triple is claimed to
replicate into infinitely many bad supports or orders.
