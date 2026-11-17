# Real-seed density bounds and finite period traces

15 November 2026. Internal written proofs with selected exact controls;
no novelty, independent specialist replay, external review or formal
verification is claimed. See the [checker](verify_seed_density_and_trace.py),
[receipt](seed_density_trace_verification.json) and
[source audit](seed_density_trace_source_audit.json).

The actual Haar seed law has an additional property beyond positivity:
on every compact interval inside its range it dominates a positive
multiple of Lebesgue measure, including at critical values. This
excludes the **particular** inward-branch interval law constructed in
the [filled-graph note](POSITIVE_GRAPH_HULLS_AND_INTERIOR_BALAYAGE.md)
from any smooth positive connected source, in every dimension. It does
not exclude other laws on that quartic polynomial graph.

There is also an exact finite trace of the real-seed Stieltjes transform
that computes the pure phase resolvent and the complete seed-power
marker list, with an explicit inverse identity. Finally, a meromorphic
differential pole criterion excludes every degree of polynomial phase
for a broad fixed-algebraic-curve class of seed laws. Polynomial and
hyperelliptic interval densities are examples. Arbitrary Haar seed
laws are not assumed to be in that class. Both exceptional MZ decisions
remain **OPEN**, with full non-torus-invariant and image-interior scopes
retained.

## 1. A real seed on a connected smooth source has locally positive mass

Let M be a connected smooth manifold of positive dimension, without
boundary, with a measure having continuous strictly positive density
in smooth charts. M need not be compact; a positive Gaussian source is
also allowed. Let h:M->R be C^1 and nonconstant. For each compact
interval J contained in the interior of h(M), there is c_J>0 such that

    (h_*measure)(E)>=c_J*length(E) for every Borel E subset J. (1)

Only the source measure on a fixed compact tube is needed, so no
global lower density or finite total volume assumption is used here.

Proof. Choose two points with h values below and above J. Join them
by a smooth embedded compact arc. A sufficiently thin tube around
this arc has coordinates(s,u) in[0,1] times a fixed small transverse
ball, positive volume density at least c>0, and

    h(0,u)<min J, h(1,u)>max J,
    |partial_s h(s,u)|<=L<infinity.

For dimension1 the transverse ball is a point. For every fixed u,
the Lipschitz function s->h(s,u) covers J. Its preimage of E therefore
has length at least length(E)/L: a Lipschitz map sends a measurable
set to a set of length at most L times its length. Integrate over u
and use the positive tube density to prove(1). The change of endpoints
ensures L>0. The argument does not assume that levels in J are regular.

In particular, at every interior range value v,

    measure{v-epsilon<h<v}>=c_v*epsilon,
    measure{v<h<v+epsilon}>=c_v*epsilon                      (2)

for sufficiently small epsilon. If the pushforward has a density,
that density cannot tend to0 from either side of an interior value.
It may diverge or have critical-value singularities.

There is a dimension-free consequence for the existing
[regular outer-boundary cap estimate](POSITIVE_GRAPH_HULLS_AND_INTERIOR_BALAYAGE.md).
For a pure graph phase p(h), if an INTERIOR parameter value v maps
to a regular analytic arc of the outer boundary and p has order k
there, that estimate gives a seed band O(epsilon^k). Equation(2)
forces k=1 in **every source dimension**. Thus an interior critical
polynomial parameter cannot map to a regular outer arc. Endpoint
parameters are outside(1) and retain the earlier dimension-dependent
extremum bounds. Nothing here supplies a regular-arc estimate at a
singular outer vertex or an arbitrary internal graph point.

For the compact exceptional cores, normalized real seeds are regular
and hence C^1, and Haar is positive smooth volume. Thus(1)-(2) hold
on both full three-dimensional cores and on their original sphere
quotients. The same lemma applies to real polynomial Gaussian seeds;
this statement is about their distributions and does not assert any
new Gaussian pure phase.

The regularity assumption matters. On the uniform circle let
h(theta)=sgn(sin theta)*sqrt(|sin theta|). Its law has density

    2*|v|/(pi*sqrt(1-v^4)), -1<v<1,

and one-sided mass arcsin(epsilon^2)/pi at0, which is o(epsilon).
This h is not locally Lipschitz at the zeros of sin theta.

## 2. The constructed quartic interval law has no smooth Haar realization

Recall the exact construction for

    p(t)=t^4-(3/2)t^2+1/4+i*(t-t^3),
    I=[-1-delta,1], 0<delta<1/4 sufficiently small.

Its outer loop is p([-1,1]); the extra arc p([-1-delta,-1)) lies
inside the loop. For a normalized inverse Riemann map psi, the measure
assigned to this extra arc, and pulled back to its parameter, is

    dmu(t)=g(t)dt, g(t)=1-|psi(p(t))|, -1-delta<=t<-1.        (3)

The boundary remainder was positive and gave all pure moments zero.
That abstract interval law remains valid. But continuity of the
Jordan boundary map gives

    g(t)->0 as t increases to-1,
    mu((-1-epsilon,-1))=o(epsilon).                         (4)

Here-1 is an **interior** value of the full seed interval I. Equation(4)
contradicts(2). Therefore this specific full-support interval law cannot
equal h_*positive_volume for any C^1 real seed on a connected smooth
manifold as in Section1, in any positive dimension. In particular it
has no regular Haar-seed realization on either exceptional quotient
or sphere model. Real affine reparameterization does not remove the
obstruction.

This excludes a proposed measure, not the polynomial p or every
positive measure on its graph. Another interior density and boundary
remainder may obey the same sweeping equation. Existence or exclusion
for the actual Haar law on that graph remains OPEN.

## 3. A finite trace computes every real-curve basis resolvent

Let mu be a positive probability of full support a nondegenerate
compact real interval I, and define its Stieltjes germ at infinity

    G(t)=integral_I dmu(u)/(t-u).

Let p(t)=sum_(k=0)^d a_k t^k be nonconstant and define, for0<=j<d,

    F_j(s)=integral_I u^j/(s-p(u)) dmu(u).                    (5)

For large s with distinct roots xi_1,...,xi_d of p(xi)=s,

    F_j(s)=sum_i xi_i^j*G(xi_i)/p'(xi_i).                    (6)

All roots are large as s->infinity, so G is evaluated on its same
principal infinity germ. There is no unspecified continuation across
the seed interval. Formula(6) follows from the partial fraction identity

    u^j/(s-p(u))=sum_i xi_i^j/(p'(xi_i)*(xi_i-u)), j<d.

The basis resolvents determine G in the opposite direction too. Define

    b_j(t)=sum_(k=j+1)^d a_k*t^(k-1-j).

The divided-difference identity is

    (p(t)-p(u))/(t-u)=sum_(j=0)^(d-1) b_j(t)*u^j.

Using(6) and setting s=p(t), all terms at xi_i!=t vanish, while the
term at xi_i=t has the limiting factor p'(t). Hence, near infinity,

    G(t)=sum_(j=0)^(d-1) b_j(t)*F_j(p(t)).                   (7)

This is an explicit finite inverse, not just a necessary trace condition.

For a normalized real seed h and phase P=p(h), the original Haar pure
condition is exactly

    F_0(s)=1/s.                                             (8)

For the fixed marker h^j, its marked tail eventually vanishes exactly
when F_j belongs to s^-1*C[s^-1]. Thus(6)-(7) retain both pure and
marked information; merely composing p with a seed resolvent would
not give the correct phase resolvent.

If(8) holds, some fixed F_j with1<=j<d is not a finite polynomial in
1/s. Otherwise(7) would make G rational. A positive measure with
full interval support cannot have rational G: rationality gives a
linear recurrence on its moments, so after shifting there is a
nonzero polynomial r with integral u^n*r(u)dmu=0 for every n>=0.
Taking a polynomial linear combination equal to bar(r) gives
integral|r|^2 dmu=0, contradicting full interval support. This is the
resolvent form of the finite-curve marker theorem.

More precisely, if G is transcendental over C(t), some F_j with
1<=j<d must be transcendental over C(s); otherwise(7) would make G
algebraic. This last conclusion is conditional on the actual seed
resolvent, not an assertion that every seed has transcendental G.
The [original matrix marker bound](FILTERED_MATRIX_DEGREES_AND_CUBIC_THIN_PHASES.md)
still gives degree<=D-8/12 for a bad seed-power marker. No finite
moment prefix or general period producer is computed by(6)-(7).

## 4. A compact-curve pole criterion for all degrees

Let C be a compact connected Riemann surface, f a nonconstant meromorphic
function, omega a meromorphic differential, and gamma a finite piecewise
analytic path chain whose endpoints and support avoid poles of f and omega.
Suppose

    poles(omega) subset f^-1(infinity).                      (9)

Define the infinity germ

    A(s)=integral_gamma omega/(s-f).

Then

    A(s) rational => A(s)=0.                                (10)

Consequently eventual vanishing of integral_gamma f^m omega forces
vanishing also at m=0. This conclusion does not need positivity.

Here is the local argument, including the finite-endpoint case.
After transporting the integration chain while avoiding the moving points
f^-1(s), the only possible finite singular values are critical values
of f and images of chain endpoints. Outside small disks about their
preimages, the integrals are holomorphic in s. Inside such a disk, put
f=c+u^k. By(9), omega is holomorphic there. Endpoint arcs have the
local Cauchy form

    integral a(u)du/(s-c-u^k), a holomorphic,

and changes of path add only residue terms a(u)/(k*u^(k-1)).
The absolute homology class of the chain in the fixed surface is
unchanged; differences created by moving punctures are small loops
around those punctures, so no other type of finite singular term is
introduced. Fixed handle portions stay outside the disks.

Taylor expansion, or division by u^k followed by scaling, bounds these
terms on each fixed local sector by

    O(|s-c|^(1/k-1)*(1+|log(s-c)|))+O(1).                    (11)

In particular they are o(|s-c|^-1). Higher vanishing order of omega
only improves the bound. If A is rational, it is single-valued and
cannot have a pole of any positive integer order at c, by(11). It has
no other finite poles. Since A(s)->0 at infinity, it must be zero.
If the moments eventually vanish, A is a finite polynomial in1/s,
so(10) applies and its coefficient at1/s, namely integral omega,
vanishes.

For closed chains on an arbitrary compact surface, the pole-compatible
case and its rationality/vanishing mechanism are already described in
[Gavrilov--Pakovich, Section2, equations(5)--(8) and the continuation argument](https://arxiv.org/pdf/1107.3029).
Their stated integration curve is closed. The endpoint analysis above
retains the finite open-chain scope used here. For C=P^1 the rational
case is also the pole mechanism in
[Pakovich, Theorem3.4 and Corollary3.5](https://arxiv.org/pdf/0910.2105).
That published theorem assumes function poles(q) subset poles(f) and
f(infinity)=infinity in its chosen coordinate. The coordinate-free
differential condition(9) is the hypothesis used here. No priority
claim is made for this pole argument.

The same conclusion holds for a **finite sum** of such curve-chain
integrals, even on different surfaces. The sum has the same subpolar
local bounds and tends to0 at infinity. Rationality of the sum does
not require or imply rationality of its individual summands; the local
bound proves directly that the rational sum is zero.

## 5. Conditional exclusions for actual seed laws

Suppose a real-seed functional admits a fixed compact-curve model

    integral r(t)dmu_h(t)=integral_gamma r(pi)omega
                                  for every polynomial r,   (12)

where pi is nonconstant meromorphic on C, integral_gamma omega=1, and

    poles(omega) subset pi^-1(infinity).                     (13)

Finite sums of representations(12)-(13) are allowed, with total mass1.
For every nonconstant polynomial p, f=p(pi) has the same pole set
as pi. Section4 shows that p(h) cannot have all positive pure
moments zero, in any degree. Establishing(12)-(13) for a candidate
Haar seed is essential; they are not automatic from positivity,
real-valuedness, algebraicity of the source, or holonomicity of its
periods.

Useful classes already excluded by this criterion include:

* Polynomial interval densities q(t)dt of nonzero total mass.
  Take C=P^1,pi=t,omega=q(t)dt. The published interval result is also
  [Francoise--Pakovich--Yomdin--Zhao, Proposition3.3 and Theorem3.3](https://www.math.bgu.ac.il/~pakovich/Publications/FPYZ.pdf).
  The finite-sum version also excludes piecewise polynomial densities
  on a finite interval partition, and finite mixtures of the classes
  below, whenever the total mass is nonzero.
* Densities r(t)dt/sqrt(F(t)), where r is polynomial and F is
  nonconstant squarefree, on a finite lifted interval chain avoiding
  poles and of nonzero total mass. On the compact hyperelliptic curve
  y^2=F(t), omega=r(t)dt/y has no finite poles: at a simple branch
  point t-a=u^2, the factor dt cancels the zero of y. Only points over
  infinity can be poles. Positive integrable branches are included.
* Polynomial tilts of Gegenbauer weights
  q(t)*(1-t^2)^(lambda-1/2),lambda in(1/2)Z_{≥0},
  with nonzero total mass. Integer exponents are polynomial weights;
  half-integer exponents are the preceding hyperelliptic form with
  F=1-t^2. The untilted version is the published Proposition4.2 in
  the same source; the tilt and compact-curve scope follow from(13).

No claim is made for every algebraic density. Finite poles of the
representing differential away from pi^-1(infinity) violate the
hypothesis, even if the chosen real branch is smooth and integrable.
Nor is a varying-fibre Abelian period automatically a meromorphic
differential on one fixed finite curve.

## 6. Updated all-degree frontier

The particular quartic abstract law from the preceding construction is
now excluded as a smooth seed law. Other distributions on the same
inward-branch graph remain possible targets. Every proposed actual
Haar seed law must satisfy the interior lower bound(1), the existing
outer-domain sweeping balance, and all original regularity constraints.

For a supplied real seed the finite trace identities(6)-(7) isolate the
pure condition and all complete seed-basis marked resolvents. A
pole-compatible fixed algebraic-curve representation closes its whole
polynomial phase class at once, without a support or degree sweep.
The outstanding issue is whether the remaining actual Haar seed
periods satisfy the required trace identity, or admit a broader
obstruction. The full image-interior/all-marker problem, non-T-invariant
thin phases and the separate symmetrized height gate remain OPEN.
Neither exceptional MZ decision nor the full active goal is complete.
