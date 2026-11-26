# Eventual Strong Factorial in the strict subcritical regime

Fix integers a>=1,d>=3 such that

    4*a < (d-1)^2.                                      (S1)

Put b=(d-1)/2, D=a+d, A=a+b, q=sqrt(b^2-a)>0, and define

    psi_w(S)=S^a*(w+S^d)*exp(-S),
    I_m(w)=integral_0^infinity psi_w(S)^m dS,
    A_m(x)=sum_(k=0)^m binom(m,k)*(a*m+1)_(d*k)*x^(m-k).

The two meanings of A are distinguished typographically here:
the bare A=a+b is a constant; A_m is the actual polynomial.
The letter h below is used for real critical heights. Orders m,N
are positive integers. Thus every period integrand is entire in S.

**Theorem StrictSubcriticalEventualConsecutiveCoprimality.**
For every fixed integer pair satisfying(S1), there is an
unevaluated integer N_0(a,d) such that the actual polynomials
A_N,A_(N-1) are coprime over Q for every N>=N_0. Equivalently,
they have no common zero over C, or over any characteristic-zero
field.

The result is eventual, with no effective cutoff. No assertion is
made for the remaining finitely many terminal orders or for a
uniform bound as a,d vary.

The proof verifies a globally ordered physical pair in the lower
wedge, then uses complete one- or two-Gaussian contours on a compact
root-containing disk. A positive real coexistence point involves
two NONDEGENERATE maxima, so it does not require a quartic local
theorem. The critical equality case is outside(S1) and is treated
by a separate fixed-cusp proof.

## 1. The physical stationary cover and its critical values

Write Phi=a Log S+Log(w+S^d)-S locally. Its stationary equation is

    w=Omega(S):=S^d*(D-S)/(S-a).                          (S2)

No stationary point in the open lower physical wedge

    C={S: -pi/d<arg S<0}

is a zero or pole of psi_w. Differentiation gives exactly

    Omega'(S)
      =-d*S^(d-1)*[(S-A)^2-q^2]/(S-a)^2.                 (S3)

The two extra critical points are S_-=A-q and S_+=A+q.
They lie strictly between a and D. In fact

    S_+=(b+q)*(b+1-q),
    S_-=(b-q)*(b+1+q).

The corresponding critical values are

    w_+=(b+q)^(d-1)*(b+1-q)^(d+1),
    w_-=(b-q)^(d-1)*(b+1+q)^(d+1).                      (S4)

Both are positive. Treating q as a variable with0<q<b, direct
logarithmic differentiation gives

    d/dq log(w_+/w_-)=4*d*q^2/[a*(a+d)]>0,

and the ratio tends1 as q tends0. Thus 0<w_-<w_+.
All branch values of Omega are real:0,infinity,w_-,w_+.
Consequently every solution of(S2) is simple for Im(w)>0.

There are exactly two roots in C for every Im(w)>0. Here are
boundary and counting details, independent of a cusp-only theorem.
The positive boundary maps to the real w-axis. On the lower
angular ray S=r exp(-i*pi/d), S^d is negative real and

    Im[(D-S)/(S-a)]=-d*Im(S)/|S-a|^2>0.

Therefore that angular boundary maps strictly to Im(w)<0.
For w confined to a compact subset of the upper half-plane,
no root can reach0 or infinity, by the leading and constant
terms of S^(d+1)-D*S^d+w*S-a*w. No root can cross either finite
side boundary by the preceding signs.

At a large upper-half-plane w there is one root near the pole,

    S=a+d*a^d/w+O(w^-2),

which lies in C, and exactly one large root in C,

    S=(-w)^(1/d)*(1+o(1)),
    arg S=(arg w-pi)/d+o(1).

The other large roots lie outside this wedge. Thus the constant
wedge root count is two. The roots are simple and the parameter
upper half-plane is simply connected, so these two roots label
globally holomorphic branches once their labels are fixed below.

## 2. A vertical divider and the global real-part order

Set S_c=b*(b+1), so A<S_c. For S=A+i*y with y>0,

    arg Omega(A+i*y)
      =d*atan(y/A)-atan(y/b)-atan(y/(b+1)).                (S5)

When this point lies in the upper reflected physical wedge,
atan(y/A)<pi/d. The expression in(S5) is strictly positive and
strictly less than pi. To prove positivity, use A<=S_c and put
Y=y/S_c. The expression is at least

    chi(Y)=d*atan(Y)-atan(b*Y)-atan((b+1)*Y).

A direct derivative calculation yields

    chi'(Y)=d*S_c*Y^2*[1+(S_c-1)*Y^2]/
      [(1+Y^2)*(1+b^2*Y^2)*(1+(b+1)^2*Y^2)] >0           (S6)

for d>=3,Y>0, since S_c>=2, while chi(0)=0. The strict upper
bound follows by dropping the two positive subtracted terms
from(S5). Hence Omega(A-i*y) is strictly in the LOWER w
half-plane whenever A-i*y is in C.

The two upper-w branches can therefore never cross the line
Re S=A. At large upper w, the pole root has Re S<A and the
large root has Re S>A. Denote them S_L,S_R. Throughout the
whole upper half-plane,

    Re S_L < A < Re S_R.                                (S7)

Their difference never has zero real part. This ordering is
for these two distinguished physical branches; it does not
compare all stationary points and does not assert that every
stationary point contributes to the original cycle.

The stationary envelope is

    d*w*partial_w(Phi_R-Phi_L)=-(S_R-S_L).                (S8)

Indeed stationarity gives S^d/(w+S^d)=(S-a)/d and thus
w Phi_w=(D-S)/d. Stationary root derivatives drop out before
subtraction. This also proves(S8) across any locally chosen
logarithm.

## 3. Positive real coexistence and strict inactive spinodals

The restriction of Omega to(a,D) decreases from infinity to
w_- on(a,S_-), increases from w_- to w_+ on(S_-,S_+), and
decreases from w_+ to0 on(S_+,D), by(S3). For w in(w_-,w_+)
there are exactly three positive real stationary points

    p_1 in(a,S_-), p_2 in(S_-,S_+), p_3 in(S_+,D).

The outer two are nondegenerate maxima of the positive real
psi_w, and p_2 is a minimum. For example, differentiating
Phi_S(S,Omega(S))=0 and using Phi_Sw<0 shows that
sign(Phi_SS)=sign(Omega') on this real interval.

At w=w_- the points p_1,p_2 coalesce at S_-. The critical
equation has a double zero there, so Phi_S has the same sign
on its two sides: the point is a strict stationary inflection,
not a local maximum. The remaining p_3 is a simple maximum.
It is the unique global maximum, and its value is STRICTLY
larger than the value at S_-. Otherwise an interior point
not a local maximum would attain the global maximum, which
is impossible. The symmetric assertion holds at w_+: p_2,p_3
merge at a strict stationary inflection and p_1 is the unique
strictly dominant maximum.

Thus the height difference between the two outer maxima,

    H_+(w)=Phi(p_3,w)-Phi(p_1,w),

has a positive limit at w_- and a negative limit at w_+.
On(w_-,w_+) equation(S8) gives

    d H_+/d log w=-(p_3-p_1)/d<0.

There is exactly one coexistence point

    w_*=w_*(a,d) in(w_-,w_+)                             (S9)

at which their heights agree. Both maxima there are simple
and separated, and the real minimum between them is strictly
below them.

For every real w_0>=0 other than w_*, the positive cycle has
a unique nondegenerate global maximum. This includes both
spinodal parameters, whose degenerate stationary inflections
are strictly inactive, as proved above. At w_0=0 the integrand
psi_0=S^D exp(-S) has its unique maximum at S=D.

Near w_* split the original positive ray at a fixed real
point between p_1 and p_3 whose phase is strictly below the
two maxima, for example the intervening minimum at w_*.
Each exact endpoint period has one nondegenerate Gaussian
maximum and a strict global gap on its complement. These
two complete packets persist uniformly in a small COMPLEX
w-neighborhood by the whole-cycle estimates in Section6.

## 4. The single physical equal-magnitude arc is Stokes-free

The critical factors are nonzero. Define their holomorphic
logarithmic lifts on Im(w)>0 by

    Phi_i=log(d)+D*Log(S_i)-Log(S_i-a)-S_i.               (S10)

Both arguments Arg S_i and Arg(S_i-a) lie in(-pi,0), with
Arg S_i more precisely in(-pi/d,0); these logarithms are
globally holomorphic on the indicated physical branches.
Let E=Phi_R-Phi_L and H=Re E.

At a fixed angle theta in(0,pi), (S7)--(S8) give

    partial_(log|w|) H=-Re(S_R-S_L)/d<0.                 (S11)

As w tends0 along that ray, S_L tends0 while S_R tends D.
The left factor has magnitude of order |w|^(1+a/d) and
the right tends D^D exp(-D)>0, so H tends plus infinity.
As |w| tends infinity along that ray, S_L tends a and its
factor grows like a positive constant times |w|. Meanwhile
S_R=(-w)^(1/d)(1+o(1)) has positive real part of order
|w|^(1/d), so its factor is only a power times exp(-S_R)
and H tends minus infinity. These limits are uniform near
either endpoint angle0 or pi: cos(pi/d)>0 bounds the large
root's real part from below.

Hence H=0 is exactly one smooth radial graph in each
upper-half-plane ray. It ends at w_* on the positive axis.
For completeness, the boundary identification does not
assume that an arbitrary nonreal saddle is inactive.
The physical branch factors extend continuously to every
positive real w, including their simple branch mergers.
Their difference E has a continuous boundary derivative
given by(S8), because S_R-S_L stays finite and neither
factor vanishes. The divider order remains strict on
this boundary: a limit at A-i*y with y>0 is excluded by
(S5); at S=A, Omega'(A)>0, so its local root for upper
w would have Im S>0 and cannot be a physical branch.
Thus the boundary H is strictly decreasing in log w.
In(w_-,w_+) the two branches limit to the outer real maxima,
so their boundary H vanishes precisely at(S9). Compactness
from the uniform endpoint asymptotics then forces the upper
radial graph to end there.

On the negative real axis, there are unique real stationary
maxima p_L in(0,a) and p_R>D on the intervals separated by
zeta=(-w)^(1/d)>0. The map Omega is strictly decreasing on
both intervals. Their real height difference decreases
strictly with log|w| by(S8), from plus infinity to minus
infinity, so there is exactly one negative crossing.
The upper physical branches limit to these real points:
the root count near either real simple maximum and
Omega'<0 places its upper-w continuation below the real
S-axis, hence inside C.

At the positive coexistence endpoint, (S10) has E tending0.
At the negative crossing, Arg(S_L-a) tends-pi whereas all
other real arguments tend0, so Im E tends EXACTLY -pi.
No unspecified winding is permitted by the explicit logs.

Writing S_R-S_L=P+iQ with P>0 and t=log|w|, the equations
on H=0 are

    H_t=-P/d, H_theta=Q/d,
    t'(theta)=Q/P,
    d(Im E)/dtheta=-(P^2+Q^2)/(d*P)<0.                   (S12)

Therefore -pi<Im E<0 along the complete connected
equal-magnitude arc. In particular a possible Stokes point,
where lambda_R/lambda_L is positive real, NEVER has
equal magnitudes. No claim that both saddles are active at
every parameter follows from this geometric fact.


## 5. Revalidate the actual two-handle argument at these parameters

Theorems stated only for integer cusps are not silently applied
to the present strict parameters. We list and verify the
hypotheses used by their contour proofs.

The normalized vector field and its side signs are unchanged:

    S_dot=-1/Phi_S,   d log|psi_w|/dt=-1,
    B(S)=S Phi_S=a-S+d*S^d/(w+S^d).                     (S13)

On the positive angular side,
Im B=-d*r^d*Im(w)/|w+r^d|^2<0. On the lower side,

    Im B=r*sin(pi/d)+d*r^d*Im(w)/|w-r^d|^2>0.

Thus normalized descent is strictly inward on both sides.
There is one interior zero zeta=(-w)^(1/d), and the ends have

    log|psi_w|=a log|S|+O(1) near0,
    log|psi_w|=log|S-zeta|+O(1) nearzeta,
    log|psi_w|=D log|S|-Re S+O(|S|^-d) at infinity.

Since cos(arg S)>=cos(pi/d)>0, all three ends decay.
Choose disjoint radial end caps. A sufficiently negative
sublevel has exactly these three contractible components:
the remaining compact set has a finite lower bound for
log|psi_w|. Every finite closed band is compact away from
the ends.

The compactified closed wedge is a disk X, with E={0,zeta,
infinity}. On each critical-free band, stop(S13) at its
lower boundary; the stopping time is the exact phase
difference, the flow is inward, and compact ODE continuation
gives a strong deformation retraction. Thus there are no
extra boundary handles. Both stationary points are simple
by Section1. The analytic normal form

    Phi=Phi_i-u^2/2

adds one index-one handle at each. Three starting components
and two handles reach the connected contractible X.
Consequently both handles connect components and their
incidence is a tree, without a loop.

At non-Stokes w the descending arms cannot reach the other
critical point: psi has constant argument on a descent,
whereas the two critical factors have different arguments.
They cannot return to their own point after a phase loss.
Compact finite-band continuation and the stable caps force
each arm to one of the three ends. Thus the actual thimbles
are the two tree edges. Their relative classes are an
integral basis of H_1(X,E), and the original0-to-infinity
path has coefficients0,+/-1 given by its unique tree path.

These relative identities are exact PERIOD identities.
For integer a,m the differential is entire. The small end
cap integrals are O(r^(am+1)) at0 and O(r^(m+1)) atzeta.
The outer sector arc is bounded by a power of R times
exp(-m R cos(pi/d)). These vanish. Truncation followed by
Cauchy's theorem therefore identifies periods of the
endpoint-relative contours; no nonzero finite cap integral
is discarded.

On each compact non-Stokes cell, each complete thimble has

    P_i(m,w)=lambda_i(w)^m*sqrt(2*pi/m)*c_i(w)
                             *(1+O(1/m)),
    c_i=1/sqrt(-Phi_SS(S_i,w)) !=0.                      (S14)

The proof is direct: a symmetric real local Morse arc
gives the Gaussian term and O(m^-1) relative error;
critical-argument separation excludes the other saddle
from all nonlocal arms, a fixed local phase loss excludes
return to the same saddle, and the uniformly stable end
caps have integrable velocities. Finite regular bands are
compact, so the complete tails have strict exponential
loss uniformly on a smaller compact collar.

At a possible Stokes point, Section4 guarantees unequal
critical magnitudes. The finite-sublevel detour from
Section2 of CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md
uses only the hypotheses just verified, not the equations
a=s^2,d=2s+1. To spell out the construction: choose a level
c strictly between the two critical heights, with a uniform
margin. The lower critical point has a complete regular
thimble because it cannot descend to the higher one.
For the higher critical point, stop its two arms at c,
above the lower point. The sublevel at c has two components,
one containing two ends and the other the third. The upper
handle joins these components. Join the two cutpoints
within their respective low-level components to selected
ends, below c plus a small fixed margin. The resulting
higher contour has exactly one Gaussian peak and a strict
loss on every connector and end tail. Together with the
lower thimble it is an endpoint tree basis with fixed
coefficients across the wall. If the original coefficient
of the higher contour is nonzero it dominates; if zero,
the original integral equals the lower period up to sign.
Either case gives a uniformly zero-free neighborhood for
all sufficiently large m. No Stokes jump formula is used.

Thus all regular upper-half-plane points have a neighborhood
where either a complete one-packet estimate excludes zeros
or a complete two-packet estimate is available for the
radial consecutive comparison in Section7. Lower-half-plane
points follow by conjugation. No activation assumption for
all stationary points, or even for both physical saddles
throughout the plane, has entered the argument.

## 6. All real boundary points, including coexistence

The following elementary whole-cycle estimate justifies
the necessary uniform complex neighborhoods. On the
positive real S-axis,

    |psi_w(S)-psi_(w_0)(S)|
        <=|w-w_0|*S^a*exp(-S).                          (S15)

Its right-hand factor has finite supremum and integral,
and integral|psi_w| is uniformly bounded for w near w_0.
Therefore a strict global supremum gap away from a chosen
real maximum persists throughout a small complex w-disc.
For any fixed contour piece with sup|psi_w|<=M<|lambda_i|,

    |integral psi_w^m dS|
       <=M^(m-1)*integral |psi_w(S)| |dS|.               (S16)

This bounds the full tails, not only finite pieces.

For w_0>=0,w_0!=w_*, Section3 gives one nondegenerate strict
global maximum, even at the inactive spinodals and at0.
Deform only a small interval about that maximum to its
analytic local Gaussian arc. Its short endpoint connectors
retain a strict loss, as do all remaining real tails by
(S15)--(S16). The nonzero leading factor is analytic.
Hence I_m is zero-free on a complex neighborhood of every
such w_0 for all sufficiently large m.

At w_* use the exact split from Section3. The left and
right partial integrals have their respective simple
global maxima and a uniform loss on the rest. Local Morse
deformation and(S15)--(S16) give two uniform complete
Gaussian packets on a complex neighborhood. The two
stationary branches remain distinct and Re(p_3-p_1)>0.
A Stokes ratio of the two factors at the real center is
irrelevant: the two original subinterval contours, rather
than complete downward thimbles across a wall, supplied
the two packets directly.

At any w_0<0 split the actual period at the analytic
endpoint zeta(w), a positive real number at w_0:

    I_m=integral_0^zeta(w) psi_w^m dS
              + integral_zeta(w)^infinity psi_w^m dS.   (S17)

On the real intervals separated by zeta, |psi_(w_0)| has
exactly two simple maxima, one p_L in(0,a), the other
p_R>D. This follows from the monotonicity of Omega on
those intervals. There is a strict gap elsewhere.
The same local-Gaussian and full-tail estimates give
two complete packets uniformly on a complex neighborhood,
using a short moving connector at zeta. The leading
factors are nonzero. Their factor ratio is negative real
at w_0, and the local logarithms exist because neither
factor vanishes. Moreover Re(p_R-p_L)>0 throughout a
smaller disc. Thus Section7 applies. No false claim that
each individual I_m is zero-free near the negative
equimodular crossing is needed.

## 7. Consecutive comparison, compact covering and arithmetic

Suppose a neighborhood has two complete Gaussian packets.
Their orientations or tree coefficients can contribute
signs, which do not affect the modulus equation. Put

    E=log(lambda_R/lambda_L),
    B=log|c_R/c_L|.

At an actual zero of I_m, the complete expansions imply

    m Re E(w)+B(w)=O(1/m).                              (S18)

The leading factors and their derivatives are bounded
and bounded away from0 on a compact smaller neighborhood.
By(S7) in the upper physical region, or by the two-real-
maxima separation near either real boundary, there is a
positive lower bound Re(S_R-S_L)>=c. The exact own-order
normalization is

    I_m(w)=(am)! A_m(m^d w)/m^((a+d)m+1).               (S19)

For a common root x at N,n=N-1 set w=x/N^d and

    m=N-theta, w_theta=(N/(N-theta))^d*w, 0<=theta<=1.

For large N this path remains in a chosen compact collar.
The upper zero and(S18) give Re E(w)=O(1/N); (S8) keeps
that bound along the path. The exact derivatives are

    dE(w_theta)/dtheta=-(S_R-S_L)/m,
    d[m Re E(w_theta)]/dtheta
          =-Re E(w_theta)-Re(S_R-S_L),
    dB(w_theta)/dtheta=O(1/N).                           (S20)

Consequently the real zero-equation expression decreases
by at least c+O(1/N) between its two endpoint orders.
Both actual zeros would force its endpoint values to be
O(1/N), a contradiction. Only values of uniform Gaussian
remainders are used; no unspecified O-term is differentiated.

For each w in C choose a controlled neighborhood: a regular
or Stokes neighborhood in the upper or lower half-plane,
a real one-packet neighborhood, or a real two-packet
neighborhood at coexistence or on the negative axis.
The preceding argument yields eventual pair exclusion
on a smaller neighborhood. A finite covering proves
eventual pair exclusion on every compact subset of C.

Only the retained OLD global Gamma upper bound is needed
to contain all possible roots. Theorem thm:product-annulus
in GAMMA_PRODUCT_ROOT_MODULUS.tex gives

    |x| <= ((a+d)m-d+1)_d <= ((a+d)m)^d

for each root of A_m, so every upper-own-order root
w=x/N^d lies in the fixed compact disk |w|<=(a+d)^d.
Apply the finite-cover theorem to that disk. This excludes
all actual consecutive common complex roots for N large,
proving StrictSubcriticalEventualConsecutiveCoprimality.

The coefficients of each actual A_m are integers. The
absence of a common complex root is equivalent to gcd1
over Q, so the Euclidean algorithm supplies a Bezout
identity with rational polynomial coefficients. The
identity remains valid over every characteristic-zero
field under its canonical copy of Q. Conversely such
a universal statement includes C. This proves the
stated arithmetic equivalence without a claim about
zeros only in a specified non-algebraically-closed field.

## 8. The exact square-gap formulation and remaining scope

For positive integers a,d, the inequality

    4a<=(d-1)^2

is equivalent to sqrt(a+d)-sqrt(a)>=1. Indeed the latter
is equivalent to d>=1+2sqrt(a); both sides d-1 and
2sqrt(a) are nonnegative in the regime in question,
so squaring is reversible. Equality for integers forces
d=2s+1 and a=s^2 for some positive integer s, hence
(a,a+d)=(s^2,(s+1)^2). Thus the independent eventual
integer-cusp theorem at equality, together with the
strict result proved here, covers every fixed ordinary
integer pair with square-root gap at least1.

This note proves only the STRICT case. The separate
critical-cusp theorem is a necessary premise for that
combined conclusion; this note does not reprove it.
No effective cutoff, uniform terminal bound, all-order
claim, or closure for the remaining supports follows.

Proof-source context:
the one-dimensional Morse, exact endpoint contours and
Gaussian framework is classical. The actual statements
needed here were revalidated in Section5 rather than
importing a theorem whose hypotheses only said cusp.
CUSP_TWO_THIMBLE_POSITIVE_CYCLE.md supplies the detailed
finite-band, handle, cap and Gaussian proofs, and
CUSP_GLOBAL_REGULAR_COMPLEMENT_AND_STOKES_DETOUR.md
supplies the finite Stokes detour construction repeated
above. Their geometric mechanisms use precisely the
explicit hypotheses verified here. The old global
Gamma theorem is the only root-bound premise; no sharp
new lower-annulus result is needed. No canonical
source or status file was changed by this research note.

