# Uniform eventual noncollision on compact subcritical ratio intervals

Internal mathematical proof. Canonical sources remain unchanged.

**Theorem CompactSubcriticalRatioUniformEventualSFC.** Fix real numbers
0<c_-<=c_+<1/4. There is an integer N0(c_-,c_+), not evaluated here,
such that for every pair of integers a>=1,d>=2 satisfying

    c_- <= c:=a/(d-1)^2 <= c_+,

and every N>=N0, the actual rational polynomials

    A_(m;a,d)(x)=sum_(j=0)^m binom(m,j)(am+1)_(dj)x^(m-j)

at m=N and m=N-1 are coprime. Equivalently, at these terminals no
nonzero X^a(A+B X^d) over any characteristic-zero field has both ordinary
factorial moments zero.

The cutoff is uniform on this fixed ratio interval. It is not evaluated,
is not uniform as c_- approaches0 or c_+ approaches1/4, and is not an
all-order assertion. The theorem's bounded-d fallback uses only the
retained strict-subcritical fixed-support theorem, not the407 all-support
theorem. The new unbounded-d argument uses complete actual contours, not
the model circle law or coefficient-relative convergence.

## 1. Exact normalization and a holomorphic compactification

Write

    s=(d-1)/2, C=a+s, L=2C/d, kappa=(d+1)/(d-1),
    p=1+2t/d, w0=C^d*kappa, S=C+Lt=Cp.

The original positive-cycle integral is

    I_m(w)=integral_0^infinity [S^a(w+S^d)]^m exp(-mS)dS
          =(am)! A_m(m^d*w)/m^((a+d)m+1).               (S1)

For integer m>=1,

    I_m(w0*q)=L [C^(a+d)*exp(-C)]^m J_m(q),
    J_m(q)=integral_(-d/2)^infinity Psi_(d,c)(t,q)^m dt,
    Psi_(d,c)=p^a exp(-Lt)(kappa*q+p^d).                (S2)

The displayed scalar is nonzero. For the ACTUAL integer supports and
integer orders the finite t-integrand is entire. Finite contour
deformations therefore remain valid even if a middle path crosses a zero.
No complex-order identity across such a zero is used.

Put sigma=1/d. On each fixed t-disc and small complex sigma-disc, choose
the logarithm of1+2sigma*t which vanishes at t=0. Cancellation of the
apparently singular linear term gives exactly

    a Log p-Lt
      =c(1-sigma)^2
          [Log(1+2sigma*t)-2sigma*t]/sigma^2
          -(1-sigma)t,
    p^d=exp[Log(1+2sigma*t)/sigma],
    kappa=(1+sigma)/(1-sigma).                          (S3)

The removable values at sigma=0 give the jointly holomorphic family

    Psi_(0,c)(t,q)=exp(-2c*t^2)*(q*exp(-t)+exp(t))
                 =:Psi_beta(t,q), beta=1/(4c)>1.        (S4)

Consequently every fixed derivative converges uniformly on fixed
t/q-polydiscs and c in the prescribed compact interval. This is an
amplitude/phase compactification. Actual leading terms below ALWAYS use
Psi_(d,c) and its ACTUAL critical values. No m/d bound is required.

## 2. Whole original tails, including q near zero

Let q range in a fixed disk |q|<=Q. For real t>-d/2 and r>=0, put

    H=p^d/(kappa*r+p^d) in[0,1].

The exact derivative is

    partial_t Log Psi_(d,c)(t,r)
       =2[-s-Lt+dH]/(d*p).                              (S5)

The case r=0 is included: H=1 on the open original interval. Moreover

    |Psi_(d,c)(t,q)| <= Psi_(d,c)(t,|q|).                (S6)

For all large d let g_-=c_->0. The number L/d=2C/d^2 is at least g_-.
Since2s/d<1 and2(s+1)/d<=3/2, choose a fixed R>=max(1,3/g_-).
Then for t>=R,

    partial_t Log Psi(t,r)
       <=[3/2-2g_-*t]/p
       <=-g_-*t/p <=-g_-/2,

because d>=2 implies t/(1+2t/d)>=t/(1+t)>=1/2.
For -d/2<t<=-R,

    partial_t Log Psi(t,r)
       >=[-1+2g_-*|t|]/p >=g_-*|t|/p >=1.              (S7)

Increase R if necessary. These bounds are uniform for r in[0,Q].
They bound the two complete real tail integrals by fixed constants times
Psi(+-R,Q)^m/m. The finite endpoint -d/2 is included by integration of the
left bound and by the actual vanishing there.

For bounded |Im t| and bounded q, the limiting finite-rectangle estimate is

    |Psi_beta(+-R+iy,q)|
        <=K_Q exp[-2c_-*R^2+R+K_y].                    (S8)

This tends to0 as R tends to infinity. First choose R sufficiently large,
then d large by(S3). The complete real tails AND both fixed finite vertical
connectors can thus be placed uniformly below any fixed positive lower
bound on selected critical moduli. No expansion on an unbounded original
tail and no comparison of exp(mF_actual) to exp(mF_model) is made.

**Lemma ZeroDiscGaussian.** There are fixed r0>0,d0,m0 such that
J_m(q) is nonzero for |q|<r0,d>=d0,m>=m0, throughout the ratio interval.

Indeed at q=0 the limiting factor is exp[t-t^2/(2beta)], with its unique
simple global real maximum at t=beta and second logarithmic derivative
-1/beta. On the compact beta interval these maxima, their nonzero factors,
curvatures and fixed outer losses have common bounds. Use(S7)--(S8) and
(S3) to perturb the selected saddle and a finite real/horizontal contour
for a small complex q-disc. The actual complete Gaussian estimate is

    J_m(q)=lambda_actual(q)^m*m^(-1/2)
              c_actual(q)*(1+O(1/m)),                  (S9)

uniformly, with lambda_actual and c_actual nonzero and uniformly bounded
away from0. On a symmetric local saddle interval the cubic correction is
odd; the quartic and squared-cubic terms are bounded by integrable
Gaussians and give the stated relative O(1/m). The other pieces have a
fixed exponential loss. This proves the lemma. In particular a sharp
LOWER Gamma annulus is not an input.

## 3. Elementary model horizontal dominance

Write q=exp(2h+i theta), with h real. On t=x+iy put X=x-h and
b=cos(2y-theta). Apart from a constant independent of X,

    V(X)=-X^2/(2beta)-hX/beta
                      +(1/2)log(cosh(2X)+b).            (S10)

A zero of the bracket is a zero of the factor and not a maximum.
For h>0 reflection gives V(-z)-V(z)=2hz/beta>0 for z>0.
Thus a global maximum is on the favored side X<0. The derivative of
V(-z) is(h-z)/beta+sinh(2z)/(cosh(2z)+b), positive on0<z<=h.
At a stationary point z>h put v=z-h, A=sinh(2z), B=cosh(2z).
Its equation is beta=v(B+b)/A, and

    V''=-1/beta+2(1+bB)/(B+b)^2 <0.                     (S11)

For if1+bB<=0 the strict sign is immediate; otherwise the numerator
whose positivity is equivalent to(S11) satisfies

    A(B+b)-2v(1+bB)
       > A(B+b)-2z(1+bB)
       >=(B+1)(A-2z)>0.

The middle inequality uses A-2zB<0 and b<=1.
For fixed b the stationary equation has exactly one such root:
T(z)=beta*sinh(2z)/(z-h)-cosh(2z) starts at+infinity and tends to-infinity,
and every crossing T(z)=b has negative derivative by(S11).
Two distinct downward crossings would require an intervening upward or
tangent crossing. This proves a unique strict GLOBAL maximum on every
horizontal. For h<0 use reflection. The same argument with h=0 shows
that each POSITIVE stationary maximum is simple, and there can be at
most one on the positive side.

This is the elementary model lemma GENERAL_CURIE_WEISS_HORIZONTAL_DOMINANCE.md.
It supplies no Gamma conclusion on its own.

## 4. Selected saddles and the entire unit circle at beta>1

The nonzero model stationary points have map

    q=Q_beta(t):=exp(2t)*(beta-t)/(beta+t),
    Q_beta'(t)
      =2exp(2t)*[beta(beta-1)-t^2]/(beta+t)^2.            (S12)

For beta>1 the two critical points of this map are real. On the open strip
-pi/2<Im t<0 there is therefore no critical point. Its top boundary maps
to the real line apart from the pole -beta; its bottom boundary maps to
the lower half-plane. For targets in the upper half-plane an inverse
cannot escape horizontally, to the pole, or to either boundary.
For q=iM large there are exactly two inverses: one near
-beta+2beta*exp(-2beta)/q and one near Log(-q)/2. After deleting a small
pole disk, any further inverse is forced to large positive real part,
where the logarithmic equation is a contraction on its forced
O(1/log M) neighborhood. Thus every upper target has exactly two simple
inverse points, jointly locally analytic in beta,q.

They satisfy Re t_L<0<Re t_R. Neither can cross the imaginary axis:
for0<y<pi/2,

    arg Q_beta(iy)=2y-2atan(y/beta) lies strictly in(0,pi),

since beta>1. Hence Q_beta(-iy) is lower nonreal, whereas q is upper.
The large-target ordering fixes the labels.

At |q|=1 the reciprocal involution

    Q_beta(-conjugate t)=1/conjugate(Q_beta(t))

interchanges the two saddles, so they have the SAME horizontal and
opposite nonzero real parts. At h=0, the value of
beta*sinh(2z)/z-cosh(2z) at0 is2beta-1>1; its level b in[-1,1]
therefore has exactly one positive stationary root by the curvature
argument of Section3 and its negative-infinity limit. The two saddles
are exactly the two strict symmetric global maxima of that horizontal,
and the middle is strictly lower (possibly zero).

This includes q=-1 by the real odd-factor horizontal. It also includes
q=1: the two selected saddles are the nonzero real solutions
t=beta*tanh t, both strict global maxima. The middle t=0 is an unstable
stationary point there and is NOT a selected maximum. Thus no edge or
cusp occurs anywhere on the circle when beta stays in a compact subset
of(1,infinity).

By compactness of this beta interval and the unit circle, the selected
saddle separation, their negative horizontal real Gaussian curvatures,
their nonzero critical moduli, and the middle/outer losses have common
positive margins. Local analytic saddle charts give a finite cover of
this whole compact parameter set.

For off-circle targets choose the favored selected saddle in Section3.
Its existence for upper nonreal targets also follows from the above
two-strip-inverse construction; for positive or negative real targets
its favored nondegenerate real maximum extends directly. The lower
half-plane follows by conjugation. Around q=0 use Section2 instead.

## 5. Complete actual localization: only the old UPPER disk is needed

The retained Gamma product bound thm:product-annulus, with gamma=am+1,
gives every root x of A_m the estimate

    |x| <=(am+1+d(m-1))_d <=[m(a+d)]^d.                 (S13)

Hence every q=x/(m^d*w0) lies in the disk

    |q| <=((a+d)/C)^d/kappa.

As d tends to infinity this upper bound tends to exp[1/(2c)],
uniformly in the ratio interval. Choose one fixed upper disk containing
all large-d actual roots. No new lower bound is needed.

Every model point in that disk off the unit circle has either the favored
single simple whole-horizontal maximum of Section3, or the Gaussian
q=0 chart of Section2. Choose a finite central contour through that
saddle and make the complete real tails and vertical connectors small
by(S7)--(S8). All strict inequalities persist for actual large d by(S3).
The symmetric local Gaussian calculation then proves(S9) on a small
complex target neighborhood, with constants uniform on any compact set
off the circle. Leading factors and amplitudes are ACTUAL throughout.

It follows that for any sequence d_j,m_j tending to infinity in the
ratio interval, every subsequential limit of normalized actual roots
has modulus1. This conclusion is about actual whole periods. It does
not follow by transferring the finite CW circle law.

## 6. Two complete actual packets on all circle charts

Fix a circle point q_* and a limiting c_* in the ratio interval. Section4
gives two strict simple global maxima on one horizontal, a strict middle
gap and strict outer gaps. On slightly larger parameter neighborhoods
the two actual saddles t_i(q), actual factors lambda_i(q) and their
Gaussian amplitudes c_i(q) are analytic, nonzero and uniformly controlled
by(S3), including fixed q-derivatives. Choose the labels with

    Re Delta_t:=Re(t_R-t_L)>=b0>0.                      (S14)

Split the EXACT original integral at the analytic midpoint
M=(t_L+t_R)/2 into P_L+P_R. Outside a fixed real cut keep the original real
tails. On the finite middle rectangle the ACTUAL integer integrand is
entire, so deform each half to its own saddle horizontal, with finite
outer and middle connectors. At the limiting chart all middle
connectors are below BOTH peaks; shrinking the chart preserves this.
The whole real tails and outer connectors have the losses of Section2.

The same complete Gaussian proof gives, for every integer m large,

    J_m=P_L+P_R,
    P_i=lambda_i^m*m^(-1/2)*c_i*(1+R_i),
    |R_i|<=K/m,                                      (S15)

uniformly in large d,c and the complex q-chart. The prefactors are
nonzero. The exact split ensures both are ACTUAL active packets; no
inference from homology rank or a limiting zero-density theorem is used.
The negative-real chart may have a zero on the middle path, which is
harmless for integer orders.

Choose local analytic logs

    D=Log lambda_R-Log lambda_L,
    b=Log c_R-Log c_L.

At any ACTUAL integer-order zero,(S15) implies

    Re[mD(q)+b(q)]=O(1/m), |b_q|<=K, |Re b|<=K.        (S16)

No derivative of R_i or of an unspecified O-symbol will be taken.

## 7. The exact adjacent-order transport

Assume contrary to the theorem a sequence of bad pairs with N->infinity.
If d is bounded, the ratio condition bounds a as well; finitely many
strict-subcritical supports remain. The retained strict-subcritical
eventual theorem excludes all of them beyond their maximum cutoff.
Thus after passing to a subsequence d->infinity, c->c_* in the interval.

Let n=N-1 and x be an actual common root. Its own-order arguments are

    q_N=x/(N^d*w0), q_n=x/(n^d*w0)=lambda*q_N,
    lambda=(N/n)^d>1.                                 (S17)

Both arguments are actual roots. Section5 and the upper disk imply
|q_N|->1 and |q_n|->1, so lambda->1. In particular d/N->0 follows
as a CONSEQUENCE, not as an initial restriction on the Gamma expansion.
Choose a circle chart containing both endpoint arguments and the whole
radial segment, and define for0<=theta<=1

    m=N-theta, q_theta=(N/m)^d*q_N,
    E(theta)=mD(q_theta)+b(q_theta).

The support a,d,c,C,L is FIXED during this interpolation. The action and
leading amplitude are analytic functions of q. The actual stationary
identity gives exactly

    d*q*D_q=-L*(t_R-t_L),
    D_theta'=-L*Delta_t/m,
    E_theta'=-D-L*Delta_t+(d*q_theta/m)*b_q.              (S18)

These are identities for the actual phases and leading amplitudes. No
noninteger exact cycle identity is asserted.

By(S16), Re D(q_N)=O(1/N). Along the path,(S18) and bounded saddles give

    |Re D(q_theta)|<=K/N+K*L/N.

Since d/L is bounded on the ratio interval and q_theta stays in its fixed
chart,(S14)--(S18) imply for all large members

    Re E_theta' <=-b0*L+K*(1+L+d)/N <=-b0*L/2.          (S19)

But(S16) at the TWO actual integer endpoints gives
Re E(0)=O(1/N) and Re E(1)=O(1/n). Integration of(S19) is impossible,
since L is bounded below and in fact L grows proportionally to d.
This excludes the hypothetical sequence and proves the uniform cutoff.

## 8. Field and coefficient quantifiers

The above analytic argument excludes every complex common root of the
rational pair A_N,A_n. The pair is therefore coprime over Q and has a
rational Bezout identity. That identity remains valid over EVERY
characteristic-zero field; no embedding of an arbitrary field into C is
assumed. For B!=0 the exact formula is

    Lambda((X^a(A+B X^d))^m)=B^m(am)! A_m(A/B).

A zero coefficient gives a monomial, whose nonzero factorial moments
are immediate. The same cutoff works for all coefficients and all fields.

## 9. Classical antecedents and exact limits of the primary theorem

The elementary all-beta horizontal lemma is classical-model geometry;
the direct formula(S10)--(S11) is reproduced here. Kabluchko,
arXiv2203.05533v3, Theorem2.3 gives the simple favored saddle via its
upper-half-plane inverse of zeta-beta*tan(zeta). His Theorem2.4 is for
FIXED positive sigma^2 and is locally uniform off the unit circle.
Remark2.6 explicitly identifies circle Plancherel--Rotach asymptotics as
a further expectation. It supplies neither the complete two-packet
circle estimates here nor Gamma-relative errors uniform when m/d is
unrestricted. The exact model integral identity is

    integral_R Psi_beta(t,q)^m dt
      =sqrt(2*pi*beta/m)*exp(beta*m/2)
           *(-1)^m H_m(-q;4beta/m).

Its finite circle law is not a substitute for Sections2,5--7.

Actual inspected primary:
research/factorial/image-rigidity/uniform_cusp_sources/kabluchko_2203_05533v3.pdf;
Theorems2.3--2.4 and Remark2.6 on printed6, plus equations1.5,1.8--1.9.
Official version: https://arxiv.org/abs/2203.05533v3 .
Retained hypotheses used independently:
SUBCRITICAL_ORDERED_PAIR_EVENTUAL_SFC.md for the finite bounded-d fallback;
GAMMA_PRODUCT_ROOT_MODULUS.tex thm:product-annulus for the upper root disk.
No407 theorem,401 cutoff or sharp404 lower-annulus premise is needed.

## 10. Complementary compact supercritical ratios: remaining work

For a compact positive c interval strictly ABOVE1/4, beta<1. The same
normalization, real tails, off-circle one-packet estimates and q0 chart
remain available. The stationary map(S12) now has an interior lower
turning point

    t_*=-i*sqrt(beta*(1-beta)),
    q_*=exp(i theta_beta),
    theta_beta=2acos(sqrt(beta))-2sqrt(beta*(1-beta)).

There is a conjugate lower edge. The two-peak arc is the complement of
the circle gap |arg q|<theta_beta. Inside that gap one simple horizontal
peak is dominant. At the edges a uniform ACTUAL single-Ai cycle,
nonzero leading amplitude, full vector error, and radial adjacent-order
transversality must be proved uniformly in beta and the1/d perturbation.
The existing local Airy geometry by itself is insufficient. This note
does not assert the corresponding supercritical uniform theorem or
control ratios approaching0 or1/4.


## 11. Finite bad sets on each fixed ratio interval: genuine reverse

For the fixed interval I=[c_-,c_+] define

    B_I={(a,d,N): a>=1,d>=2,N>=2, a/(d-1)^2 in I,
                   gcd_Q(A_N,A_(N-1)) is nonconstant}.

**Corollary CompactSubcriticalRatioFiniteBadEquivalence.** The following
are equivalent for this fixed interval:

1. B_I is finite.
2. One terminal cutoff works for every support in I.
3. One terminal cutoff excludes both zero moments for every nonzero
   coefficient pair over every characteristic-zero field on those supports.
4. One terminal cutoff excludes both moments for all nonzero complex coefficient
   pairs on those supports.

For1=>2 take2 when B_I is empty and otherwise one plus its largest terminal.
For2=>1 use only the retained OLD fixed-terminal complements. At a fixed
N, the raw high-gap theorem gives d<D_N, where

    D_N=2(N-1)*ceil(log_2(4N(N-1))).

The OLD linear-soft Gamma theorem gives a<d(2N-3). Thus each fixed-terminal
bad support set is finite; restricting it to I leaves it finite. A finite
union over N below the uniform cutoff is finite. This is a reuse of those
universal written theorems, not a new finite terminal computation.
Rational Bezout gives2=>3;3=>4 is restriction. For4=>2, a complex common
root x supplies the forbidden phase X^a(x+X^d) at the same support/order.
The converse is not asserted for only one prescribed nonclosed field.

The theorem above proves these equivalent properties for every fixed
compact interval I inside(0,1/4). They do not imply a single cutoff on the
whole open interval, an evaluated cutoff, emptiness of B_I, or finiteness
of all ordinary-binomial bad triples. The regimes approaching c=0 or
c=1/4 and the compact supercritical sector remain separate obligations.

## 12. Exact support and source boundary

    python3 verify_compact_subcritical_ratio_uniform.py --repo-root REPO \
      --support-root STRICT_SUBCRITICAL_PACKAGE --output review.json

For an installation containing both packets in one directory, support-root
defaults to the checker directory. The mandatory source ledger authenticates
the written strict-subcritical fixed-support fallback, its existing
guard/input/ledger/receipt and all of its six declared proof inputs, the
old Gamma upper/soft theorem, and the raw fixed-terminal high-gap theorem.
The strict fallback guard is freshly replayed and its receipt is compared
byte for byte. Kabluchko's cached primary is pinned as classical context.

Finite controls support exact normalization, the stationary map, envelope
and transport identities, whole-tail inequalities, actual own-order
moments, fixed-terminal counts and scope. They are not a formal verification
of the universal contour proof and do not provide a numerical cutoff.
No401,404,407,411 or412 theorem is a premise, no heavy finite-prefix
computation is replayed, and no canonical file is changed here.

