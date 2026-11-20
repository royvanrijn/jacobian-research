# A bounded-thermal transition tube for ordinary binomial moments

Internal mathematical proof, pending canonical registration and source/code audit. This extends the equality-curve
argument only to a transition tube whose width depends on the terminal
order. It is not an all-support uniform terminal theorem, an effective
cutoff, or an all-order assertion. No canonical files or status are changed.

For integers a>=1,d>=2 put

    s=(d-1)/2>=1/2, e=a-s^2, C=a+s, L=sqrt(C),
    Sc=s(s+1), kappa=(s+1)/s, w0=C^d*kappa,
    r=e/C, beta=-Sc/(d*L),
    A_m(x)=sum_(j=0)^m binom(m,j)*(am+1)_(dj)*x^(m-j).

The parameter s may be a half-integer. For the two actual own orders use

    q_m=x/(m^d*w0), u_m=-sqrt(m)*r,
    eta_m=beta*(q_m-1), h_m=m^(3/4)*eta_m.                 (T1)

**Theorem BoundedThermalTransitionUniformExclusion.** For every
finite U>=0 there is an integer N0(U), independent of a,d and the phase
coefficients, such that if N>=N0(U) and

    sqrt(N)*|a-(d-1)^2/4|/[a+(d-1)/2]<=U,                (T2)

then A_N,A_(N-1) are coprime over Q. Equivalently, the same support/order
condition excludes simultaneous factorial-moment vanishing over EVERY
characteristic-zero coefficient field. No numerical N0(U) is supplied.

## 1. Exact jets, integrality, and the contradiction setup

Set t=(S-C)/L, p=1+t/L and subtract the central phase. Exactly

    F(t,q)=a Log p-Lt+Log[(p^d+kappa*q)/(1+kappa*q)].     (T3)

At t=0,q=1 the derivatives are

    F_t=0, F_tt=-e/C,
    F_ttt=2e/C^(3/2),
    F_tttt=-(2Sc^2+6e)/C^2,
    F_tq=-Sc/(dL)=beta, F_ttq=0.                        (T4)

Thus freezing beta at -L/d is valid only when e=0. The exact identities
C=Sc+e=Sc/(1-r), a=C-s are useful below.

If e=0, integrality of a implies s is an integer: a square of a strict
half-integer has fractional part1/4. The separately proved uniform
integer-cusp eventual theorem excludes those pairs at one fixed cutoff.
If U=0 this already proves the assertion. Hence suppose U>0 and e!=0.
Because e is quarter-integral, |e|>=1/4. Under(T2),

    C>=sqrt(N)/(4U), L>=N^(1/4)/(2sqrt(U)).              (T5)

To disprove the proposed theorem there would be common actual roots at
a sequence N->infinity satisfying(T2). After discarding the excluded
e=0 cases, every such sequence has

    L->infinity, s->infinity, r->0,
    s/L->1, d/L->2, a/C->1, beta->-1/2.                (T6)

Indeed C=Sc/(1-r), |r|<=U/sqrt(N), and(T5) imply these statements. The
proof proceeds by excluding every such sequence. Along each adjacent
comparison the integers a,d and hence e,C,r,beta STAY FIXED; only m and
q_m vary. The thermal variables satisfy u_n=sqrt(n/N)*u_N, n=N-1.

## 2. Global root capture and actual regular contours survive detuning

The exact integer positive-cycle normalization is

    J_m(q)=integral_(-L)^infinity Psi(t,q)^m dt,
    Psi(t,q)=p^a*exp(-Lt)*(kappa*q+p^d).                (T7)

It differs from the original I_m(w0q) by the nonzero scalar
L*[C^(a+d)*exp(-C)]^m. Since a,d,m are integers, the finite-contour
integrand is entire. Zero-crossing finite middle paths are legitimate;
global division by q+1 is not.

On each fixed complex t-disc, all fixed derivatives of Psi converge,
locally uniformly in q, under(T6) to

    Psi_0(t,q)=exp(-t^2/2)*(q*exp(-t)+exp(t)).            (T8)

For example a Log(1+t/L)-Lt has linear coefficient -s/L, quadratic
coefficient -a/(2C), and higher coefficients O(L^(2-j)); meanwhile
d Log(1+t/L) tends to2t. This convergence transfers only STRICT contour
geometry and gaps. Every leading exponential below uses the ACTUAL
critical value of Psi, so neither m/s nor m*r is assumed small.

For real t and positive q, put H=p^d/(kappa*q+p^d) in[0,1]. The exact
whole-tail derivative is STILL

    partial_t Log Psi=(-s-Lt+dH)/(Lp).                  (T9)

For large members of the sequence, s/L and(s+1)/L are at most3/2 and
L>=1. Thus t>=4 gives derivative <=-t/(2p)<=-2/5; t<=-3 in the
original domain gives derivative >=(-t-3/2)/p>=1. The finite left
endpoint -L is included. Also |Psi(t,q)|<=Psi(t,|q|). Exactly as in the
uniform equality-curve regular proof, first choose a large fixed real
cut R, then use(T8) on its fixed rectangle. Complete original real tails
and both vertical connectors are uniformly below every selected compact
positive critical-modulus threshold. No Taylor series is used on an
unbounded original tail.

The sharp Gamma annulus gives every own-order normalized root the bounds

    (a/C)^d/kappa < |q_m| <= ((a+d)/C)^d/kappa.           (T10)

Along(T6) the lower and upper endpoints tend to exp(-2),exp(2).
Consequently both actual endpoint roots stay in one fixed compact annulus,
for example1/32<|q_m|<16 after passage sufficiently far along the sequence.

The limiting phase(T8) has the exact contour geometry proved in
UNIFORM_CUSP_CURIE_WEISS_REGULAR_COMPLEMENT.md Sections3--6: off the unit
circle there is a complete horizontal with one strict simple maximum;
on a closed unit-circle arc avoiding1 there are two strict simple maxima
on one horizontal, with fixed middle and outer gaps. These are statements
about the limiting entire function and therefore apply unchanged here.
Their transfer mechanism is revalidated by(T8)--(T10): analytic actual
saddles and curvatures persist on fixed q-discs; the real tails and finite
vertical cuts have uniform loss; actual complete Gaussian leading terms
and O(1/m) relative errors follow on these finite horizontal contours.

Hence BOTH actual endpoint roots are asymptotically captured by the unit
circle. Write lambda=(N/n)^d>1, so q_n=lambda*q_N. Their radii force
lambda->1. After taking a subsequence both converge to the same circle
point q_*. If q_*!=1, use the two actual complete packets on its fixed
q-disc. Their exact envelope is

    d*q*D_q=-DeltaS, DeltaS=L*(t_R-t_L),
    Re(t_R-t_L)>=c>0, |b0_q|<=C.                       (T11)

The value estimate at a zero is Re[mD+b0]=O(1/m). Along
m=N-theta, q_theta=(N/m)^d q_N, the SAME calculation as in that regular
proof gives derivative -Re D-Re DeltaS+O(d/m)<=-cL/2. The endpoint values
cannot both be O(1/N). No complex-order identity for a path crossing a zero
is used; only the actual integer endpoints need the complete period split.
It follows that every putative sequence must satisfy

    q_N->1, q_n->1, lambda->1.                         (T12)

The actual polynomials have real positive coefficients. A real q_N near1
is positive and cannot be a root. By replacing a common root by its
complex conjugate if necessary, assume Im q_N>0; then Im q_n>0 and both
eta fields have negative imaginary part. This fixes the cone orientation
used below without excluding any common root.

This is not a blind application of the equality-curve theorem outside its
hypotheses: its limiting horizontal geometry is reused, while(T7)--(T11)
verify all actual detuned whole-cycle hypotheses and exact envelopes.

## 3. Complete bounded thermal/field convergence

For bounded real u_m and bounded complex h_m, put t=m^(-1/4)*v. Equations
(T3)--(T4) give, uniformly on bounded v-discs,

    mF(m^(-1/4)v,q_m)
       =-v^4/12+u_m*v^2/2+h_m*v+o(1).                 (T13)

The moving-center cubic contributes
-u_m*v^3/(3L*m^(1/4)), hence vanishes uniformly. The quartic coefficient
is -1/12+O(r); the other local derivatives are uniformly bounded. Terms
nonlinear in q-1 also vanish at this scale. The exact zero mixed quadratic
jet in(T4) is available, although mere boundedness already suffices for
the leading limit.

The expansion is a COMPLETE original-period limit. On one fixed real
local interval, Taylor's theorem uniformly gives

    Re F(t,q)<=-c*t^4+C|r|t^2+C|q-1||t|.

After rescaling the dominating envelope is
exp[-c*v^4+C_U*v^2+C_H|v|], integrable on R. The actual whole far tails
have a fixed gap by(T9) and compact continuity at r=0,q=1; detuning and
the small coefficient perturbation preserve it before multiplication by m.
Their normalized contribution is exponentially small. Therefore, after
removing a nonvanishing central scalar and L*m^(-1/4), the full period
converges uniformly on compact h/u sets to

    Z(h,u)=integral_R exp(-v^4/12+u*v^2/2+h*v)dv.         (T14)

For EVERY real u, the retained Newman/Lee--Yang theorem says all complex
zeros of Z(h,u) are imaginary. Only this location statement is needed
here, not simplicity, zero labels, velocities or a bound uniform on all
real temperatures. On compact real-u/complex-h sets, (T14) implies:

    an actual bounded-field zero has Re h_m->0.         (T15)

Indeed any contradictory subsequence has a limiting real u and complex h
with Re h!=0, but uniform convergence would make Z(h,u)=0. Actual
canonical temperature may have small complex corrections; (T13)--(T15)
are a leading compact limit and do not assert that it is exactly real.

## 4. Uniform large-field packets with small residual temperature

For eta=beta*(q-1), rho=|eta|^(1/3), k=m*rho^4 define

    T=-r/rho^2=u_m/sqrt(k), t=rho*v.

When k->infinity and |u_m|<=U, T->0. The actual scaled phase extends as a
compact holomorphic family through rho=0,T=0,nu=1/d=0 and has limit

    -v^4/12+T*v^2/2+omega*v, omega=eta/rho^3.            (T16)

For clarity the analytic parameters may be rho,T,nu with r=-T*rho^2;
C=Sc/(1-r). This removes any hidden singularity in the thermal unfolding.
All fixed-parameter complex-m and q collars below hold SUPPORT r fixed,
and CENTER rho fixed. Thus T is fixed on those collars, even though the
notation u_m/sqrt(k) contains m twice. No frozen O(u_m) is differentiated.

The complete small-T packet construction follows the uniform cusp proof
with explicitly changed bounds. On the entire intermediate real range,
the bound replacing its quartic envelope is

    Re(F/rho^4)<=-c*v^4+C|T|v^2+C|v|.

The far original tails retain their fixed loss. On compact model contours
small T preserves the outer-angle unique horizontal maximum, the central
two-maximum horizontals, the midpoint loss and the complex-m sector. The
Gaussian odd cubic term cancels; quartic and squared-cubic majorants give
R_i=O(1/k). These facts prove the SAME angular localization: every local
large-k zero lies in |arg eta+pi/2|<epsilon when Im q>0. All estimates are
uniform for T sufficiently small, rho sufficiently small and nu near0.

In the wider central cone the exact midpoint periods P_L,P_R give
I=P_L+P_R and

    P_i=exp(m Phi_i)*sqrt(2pi/m)*c_i*(1+R_i),
    |R_i|<=C/k, |R_(i,m)|<=C/(mk),
    |R_(i,q)|<=C/(k*rho^3),
    Re DeltaS>=cLrho, |DeltaS|<=CLrho.                  (T17)

The derivatives use genuine complex product collars |m'-m|<cm,
|q'-q|<c*rho^3 and exact functions, as in the equality-curve proof.

One amplitude change MUST be retained. At T=0 the leading c_R/c_L ratio
is a constant unit in omega, but for T!=0 it may depend on omega at
order T. Analytic dependence on(T,omega), with the common L/rho factor
cancelled, therefore gives

    |Re b0|<=C(rho+|T|),
    |b0_q|<=C[rho^(-2)+|T|rho^(-3)].                   (T18)

The T/rho^3 term cannot be discarded. To derive it, keep center rho fixed:
the leading ratio's omega derivative is O(T), its finite-rho correction
is O(rho), and omega_q=beta/rho^3. For
b=b0+log(1+R_R)-log(1+R_L), E=mD+b, consequently

    |Re b|<=C(rho+|T|+1/k), |b_m|<=C/(mk),
    |b_q|<=C[rho^(-2)+|T|rho^(-3)+1/(k*rho^3)].         (T19)

Every actual zero obeys Re E=0. The exact envelope(T11) still holds.

## 5. Bounded upper field: integrality defeats both endpoint zeros

First pass to a subsequence with |h_N| bounded. From(T12),
Im h_n=lambda*(n/N)^(3/4)*Im h_N is bounded. If |h_n| were unbounded,
then k_n=|h_n|^(4/3)->infinity and rho_n->0. Formula(T16) has T_n->0,
so its outer-angle localization puts the lower root in the central cone.
That cone bounds |h_n| by a fixed multiple of |Im h_n|, a contradiction.
Thus BOTH endpoint fields are bounded after passage to a subsequence.

Use(T15) at both actual own orders, with u_n=sqrt(n/N)u_N also bounded.
Then Re h_N and Re h_n tend to0. Their EXACT relation is

    h_n=lambda*(n/N)^(3/4)*h_N
                         +beta*n^(3/4)*(lambda-1).        (T20)

The first multiplier tends to1. Bernoulli's inequality gives
lambda-1>=d/n, so the modulus of the final real term is at least

    (Sc/L)*n^(-1/4)
       =L*(1-r)*n^(-1/4)>=1/(4sqrt(U))                  (T21)

for all sufficiently large members, by(T5) and |r|<1/2. Its sign is
negative. Thus the real parts in(T20) cannot both tend to0. No same-index
argument, derivative of a model zero, or small adjacent field step was
used. This excludes the bounded upper field case.

## 6. Growing upper field: both endpoint cones control the whole path

It remains |h_N|->infinity, equivalently k_N->infinity. The upper root is
in the smaller central cone. Since -Im eta_n=lambda*(-Im eta_N),

    k_n >= (n/N)*cos(epsilon)^(4/3)*k_N->infinity.

The lower root therefore also has the central cone label, because its
thermal T_n->0 and q_n->1. For m=N-theta,q_theta=(N/m)^d*q_N, let
lambda_theta=(N/m)^d and y_N=-Im eta_N>0. Exactly

    Re eta_theta/(-Im eta_theta)
       =Re eta_N/y_N-|beta|*(1-1/lambda_theta)/y_N.      (T22)

It is monotone BETWEEN its two endpoint ratios. Hence the entire path
is central, rho_theta is comparable to rho_N, and k_theta is comparable
to k_N. Here lambda->1 was already proved globally. In particular the
path and its Cauchy collars remain inside the constructed physical and
angular domains. Also |u_m|<=U and |T_theta|<=U/sqrt(k_theta).

Differentiate the exact log ratio, keeping the support fixed:

    dE/dtheta=-D-DeltaS-b_m+(d*q_theta/m)*b_q.           (T23)

At the upper zero, |Re D|<=C(rho+|T|+1/k)/N. The exact path derivative
dD/dtheta=-DeltaS/m then adds at most CLrho/N along the comparable path.
Using m=k/rho^4, d/L bounded and(T17)--(T19), every error in(T23) divided
by Lrho is bounded by

    C_U*[1/N+rho/k+1/k^(3/2)+1/k^2].                    (T24)

In particular the NEW transported amplitude term |T|rho^(-3) contributes
|T|/k<=U/k^(3/2), not rho/k. The upper-zero |T|/N term contributes at most
U*rho^3/(L*k^(3/2)) and is harmless. All terms in(T24) tend to0. The main
real part is therefore <=-cLrho/2<0 throughout the adjacent path, whereas
both integer endpoint zeros impose Re E=0. This final contradiction
excludes the growing upper field case and hence every proposed sequence.

## 7. Quantifiers, field descent, and the remaining transition problem

If no N0(U) existed, one could choose counterexamples at N tending to
infinity. The exact-cusp theorem removes e=0; the preceding argument
excludes all remaining subsequences. This contradiction proves the uniform
existence conclusion. Rational coprimality follows from absence of a
complex common root. Its rational Bezout identity persists over every
characteristic-zero field; the exact coefficient formula is
Lambda(f^m)=B^m*(am)!*A_m(A/B), and coefficient-zero phases are monomials.

The model input is only the already retained Newman location theorem for
Z(h,u) at REAL u. There is no exact finite-Gamma Lee--Yang claim, no
transfer based only on coefficient-relative errors, and no assumption
m/d small. The local/global equality-curve mechanisms are reused only
after their actual detuned hypotheses and complete tails are verified.

The theorem leaves sequences |sqrt(N)e/C|->infinity untreated. Matching
those regimes requires further uniform control near the thermal folds,
including fields of size |u|^(3/2) and the two-well side. No bound uniform
on all real model temperatures, all-index model zero separation, numerical
cutoff, all-order binomial theorem, or finite GLOBAL bad-triple set is
asserted here. All new work remains scratch and unintegrated.


## 8. Regionwise finite-bad-set and escape equivalences

For integers a>=1,d>=2,N>=2 put

    s=(d-1)/2, e=a-s^2, C=a+s,
    theta(a,d,N)=sqrt(N)*|e|/C.

Let B be the set of triples at which the two ACTUAL rational moment
polynomials A_N,A_(N-1) are not coprime, and define

    B_U={(a,d,N) in B: theta(a,d,N)<=U}, U>=0.

The following fixed-terminal finiteness input is already available WITHOUT
any new finite terminal calculation. Set

    D_N=2(N-1)*ceil(log_2(4N(N-1))).

The OLD gap theorem in lower_exponent_initial_and_boundary.tex excludes
d>=D_N for every a>=1. The OLD Gamma linear-soft theorem excludes
a>=d(2N-3). Consequently every bad triple at a specified terminal N has

    2<=d<D_N, 1<=a<d(2N-3).

This is a finite support set. The exact integer ceiling may equally be
written (4N(N-1)-1).bit_length(). Neither this finiteness observation nor
its use below requires the terminal-six finite certificate, terminal-ten
prefix, all-support eventual theorem, or any a=0/d=1 boundary theorem.

**Corollary BoundedThermalBadSetEquivalence.** For EACH fixed finite U
the following are equivalent:

1. B_U is finite.
2. There exists one terminal N0(U) such that every actual pair is coprime
   whenever N>=N0(U) and theta(a,d,N)<=U.
3. The same one cutoff excludes two vanishing factorial moments over
   every characteristic-zero field and every nonzero coefficient pair
   at all those support/order parameters.

For 1=>2 choose2 if B_U is empty, otherwise1+max of its terminal
coordinates. For 2=>1, all members have N<N0(U), and the union of the
finitely many fixed-terminal finite support sets above is finite.
The rational Bezout identity gives2=>3. Conversely a rational common
factor has a complex root x; the phase X^a*(x+X^d) has both moments zero,
contradicting3. The universal complex coefficient statement is equivalent
as well; no converse over a prescribed nonclosed field is substituted.

**Corollary ThermalEscapeForDistinctBadSequences.** The family of these
properties for EVERY finite U is equivalent
to the following escape statement:

    Every infinite sequence of pairwise DISTINCT bad triples satisfies
    theta(a_j,d_j,N_j)->infinity.

If a bounded-U subsequence existed, it would give infinitely many distinct
members of B_U. Conversely, an infinite B_U can be enumerated distinctly
and violates the escape statement. Fixed-terminal finiteness also forces
N_j->infinity along every distinct bad-triple sequence, so an equivalent
form may explicitly quantify sequences with terminal orders tending to
infinity. Repetition of one triple is deliberately excluded.

Thus the bounded-thermal theorem yields finiteness in
EACH bounded thermal tube. It does not imply that B_U is empty or that B
itself is finite: bad sequences may leave every such tube. The abstract
diagonal triples(a,d,N)=(n,n,n) have theta of order n^(3/2), illustrating
this logical escape without claiming any actual moment counterexample.
No cutoff or finite set is numerically evaluated here.

## 9. Source boundary and exact reproduction

The reviewed mathematical body was frozen at SHA256
`b47c3b457d2fe296442936e5e4fa9645386542c865305aa84a826fb0fa0e1958`.
This account changes only its process-status wording and appends Section8's
regional reverse; the analytic argument is unchanged. Its imports are the
new uniform-cusp raw mechanisms and theorem, the sharp Gamma annulus, the
Newman real-temperature location theorem, and the RAW old high-gap/old
linear-soft bounds for fixed-terminal finiteness. The broad all-support
eventual theorem, fiberwise theorem and old a=0/a=1 finite prefixes are not
premises. Unbounded-temperature matching remains open.

Run the standalone standard-library guard from any directory:

    python3 PATH/verify_bounded_thermal_transition.py --repo-root REPO \
      --support-root UNIFORM_CUSP_PACKAGE --output review.json

In a combined installation the support directory defaults to the checker
directory. The mandatory ledger authenticates all proof inputs and the
fresh upstream receipt. Universal algebra and exact finite arithmetic
support the proof's jets, own-order maps, scale powers and scope. They do
not replace the complete contour proofs or constitute formal verification,
external review or an all-order computational certificate.
