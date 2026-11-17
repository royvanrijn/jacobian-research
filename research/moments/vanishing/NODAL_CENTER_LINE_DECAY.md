# Nodal critical lines: uniform cut bounds and a zero-fiber exclusion

14 November 2026. Internal written analytic proof with exact auxiliary
algebra. The checker does not verify the analysis. Imported compactness,
semialgebraic triviality and rapid-decay comparison are identified below.
No external review, independent specialist replay, formal verification or
priority is claimed.

## Result and boundary

Retain the full affine Gaussian exponent Q from the
[nodal normal form](NODAL_GAUSSIAN_CRITICAL_PROFILES.md), with m=Q_TT!=0.
Consider its center-line boundary

    P=Tuv+a*u^3+d*v^3+b*u^2+c*v^2+k, a*d!=0.       (1)

All Gaussian means and metric cross terms are allowed. The critical values
of P are among

    k, k+A, k+B, A=4*b^3/(27*a^2), B=4*c^3/(27*d^2). (2)

Omit a nonzero-axis point when its b or c is zero; the repeated value k
then adds no condition.

**Half-plane exclusion.** If these values lie in a closed complex half-plane
through0, some positive Gaussian pure moment of P is nonzero. Equivalently,
the origin must not lie in the interior in R^2 of their convex hull.

In particular **k=0 is always excluded**, for arbitrary b,c and the full
affine Q. This closes the zero-action critical-line boundary of the
irreducible nodal problem. It includes affine translates of homogeneous
irreducible nodal cubics, with the Gaussian mean retained. The older
homogeneous theorem alone did not cover a nonzero Gaussian mean.

If k!=0 and the three values in (2) form a triangle strictly containing0,
this half-plane argument is silent. The subsequent
[Morse growth-dominance theorem](NODAL_MORSE_DOMINANCE.md) closes that
triangle and hence the entire critical-line class. This note also does not close the
finite-critical nodal profiles h!=0 or (e,f)!=(0,0), whose phase-critical
points escape. No numerical uniform pure prefix is supplied.

## 1. A gradient estimate uniform in the parameter

For S_t=tP-Q there are constants C,c,t0>0 such that

    ||grad S_t(x)|| >= c*||x||,
    |t|>=t0, ||x||>=C, x in C^3.                  (3)

These are Hermitian norms on complex affine space; Q is holomorphic.
This is stronger than separate coercivity for each fixed t.

If (3) failed, take |t_j|>=t0, ||x_j||->infinity, with
||grad S_tj(x_j)||/||x_j||->0. Dividing by t_j*||x_j||^2 shows that the
projective limit is the unique leading node. Write x_j=(u_j,v_j,T_j),
so |T_j|->infinity and |u_j|+|v_j|=o(|T_j|). The first two equations,
divided by t_j*T_j, have the form

    v+(3a*u^2+2b*u)/T = D/t+O((|u|+|v|+1)/(|t*T|))+o(1/|t|),
    u+(3d*v^2+2c*v)/T = E/t+O((|u|+|v|+1)/(|t*T|))+o(1/|t|).

Here D,E are the fixed Gaussian cross coefficients. Absorb the errors
linear or quadratic in u,v, using (u,v)=o(T). It follows that

    |u_j|+|v_j|=O(1/|t_j|).

But

    (S_t)_T=t*u*v-m*T-l_T-D*u-E*v=-m*T+O(1),

contradicting the assumed gradient ratio. This argument allows t_j to be
bounded or unbounded and proves the stated uniform estimate. Consequently
all phase-critical points remain in one fixed ball. For each t, S_t is
tame; its [affine critical length is7](NODAL_GAUSSIAN_CRITICAL_PROFILES.md).

Bounded critical points by themselves would not prove (3). For example
P=Z^2 with circular Gaussian Q=ZW+T^2/2 is pure-null and S_t has only
the origin critical, but grad S_t(1,2t,0)=(0,-1,0). Its coercivity constant
degenerates. This negative control is excluded by (3), not by point counting.

## 2. Compactified cuts and a jointly semialgebraic basis

Fix any ray t=R*exp(i*theta), R->infinity. The critical points and values
of S_t are algebraic functions of R, after a finite ramification. Raise
R0 so their branches and multiplicities stabilize. Group identical critical
values; coincident points or non-Morse points are not discarded.

Choose eta with |eta|=1 and Re eta>=1/2 so the cuts

    w=c_i(R)-eta*s, s>=0                         (4)

are disjoint for large R, except when their endpoints are identical and
already grouped. Such eta exists: each nonzero difference of two algebraic
critical-value branches has one leading direction at infinity, and there
are only finitely many differences. The argument theta need not be generic;
eta is chosen independently. This matters on the opposite-action boundary.

The ordinary tame Lefschetz-cut description gives all rapid-decay middle
homology as the sum of local vanishing homologies attached to these cuts,
of total rank7. This is the fixed-phase homology construction in
[Pham, PartI and Appendix](https://numdam.org/item/AST_1985__130__11_0.pdf).
It remains to obtain one controlled family of representatives, rather than
choose unrelated representatives for each R.

Along a cut, normalized complex gradient transport has

    dx/ds=-eta*conj(grad S_t)/||grad S_t||^2,
    dS_t/ds=-eta.

Outside the fixed ball in (3), d||x||^2/ds<=2/c. Local vanishing thimbles
start in a common fixed ball, since all critical points do. Consequently
all the generators can be carried inside

    ||x||^2<=A0+B0*s,                             (5)

with fixed A0,B0 independent of R and of the cut. Choose B0>2/c so the
transport points strictly into this region at its boundary. No gradient
trajectory escapes at a finite cut length; the cuts avoid other critical
values. Small local thimbles may depend on R and may be degenerate; their
supports still lie inside the common starting ball.

Here is a precise semialgebraic replacement of this transport. Embed x
in its spherical compactification x/sqrt(1+||x||^2), and put sigma=s/(1+s).
Take the closure of the family defined by (4),(5). Its only boundary at
spatial infinity has sigma=1, because (5) forbids spatial escape for finite
s. Let G be this rapid boundary. The compact relative pairs (Z_R,G_R)
therefore represent the cut-supported cycles in the constrained carrier.
Equivalence with the stabilized truncated relative homology follows by
semialgebraic conic structure at sigma=1; it is the relative-chain model
used in Pham's Appendix. The gradient transport above proves that the map
from these constrained carriers to the full cut homology contains every
local generator.

Apply compatible Hardt triviality to the projection onto R, including the
rapid boundary, the constrained carrier and its inclusion in the full cut
carrier as distinguished semialgebraic subsets. There is one unbounded
interval R>R0 on which this diagram is trivial. Choose finitely many
integral relative chains at one R and transport them by this trivialization.
The resulting chains form a jointly semialgebraic family, have fixed
integer multiplicities, satisfy (5), and map to a flat rank7 basis of full
rapid-decay homology. Continuous transport on every compact R-interval is
proper on finite phase bands, so it is transport in the actual homology
local system. It does not impose an unproved identification of formal
exponential blocks.

The compatible finite-family version of Hardt triviality is used here.
See [Basu--Vorobjov, Theorem1.1 and its setup](https://www.math.purdue.edu/~sbasu/definable.pdf)
and [Hardt--Lambrechts--Turchin--Volic, Section2.1](https://arxiv.org/pdf/0806.0476).
Compact semialgebraic chains and their Stokes property are as in the
latter paper, Section3. These are applications of those theorems to the
explicit compactified cut pairs, not a claim that a gradient flow itself
is semialgebraic.

## 3. A uniform volume estimate and uniform tails

A fixed semialgebraic family of integral k-chains in a ball of radius L,
with bounded multiplicities, has mass<=C*L^k, with C independent of the
parameters. One elementary proof partitions smooth tangent planes according
to which coordinate k-minor is largest. Some minor has absolute value at
least binomial(n,k)^(-1/2). The area formula bounds each piece by the
volume of its coordinate projection times its generic fiber count. The
fiber counts are uniformly bounded by semialgebraic uniform finiteness.
There are finitely many pieces and multiplicities. The same argument
bounds the boundary chains. This strengthens pointwise finite volume to
the uniform estimate needed here.

Lift a cut chain to (x,s). Its band j<=s<=j+1 lies in a ball of radius
O(1+j), by (5). The joint semialgebraic description has fixed complexity,
so its mass is at most C*(1+j)^M, independently of R. For omega=du*dv*dT,

    |exp(-c_i(R))*integral_(s>=L) exp(S_t)*omega|
      <=C*sum_(j>=floor L) (1+j)^M*exp(-j/2).       (6)

The right side tends to zero independently of R. No factor growing with R
has been hidden in C. Polynomial test forms have the same estimate with
a different M. Thus the estimate also controls any fixed polynomial
de Rham representatives, although the exclusion below only needs omega.

More generally, without (5), a jointly proper semialgebraic cut family
whose finite-band radius is polynomial in R+s gives
C*R^A*(1+j)^M*exp(-j/2). A cutoff L=K*log R then makes the tails smaller
than any assigned inverse power. This general observation is useful for
the escaping profiles, but does not itself identify their formal blocks.

## 4. Fixed bands collapse into a complex level set

Let c_i(R)/t->alpha. For (1), alpha is one of (2). On the lifted chain,

    P(x)-alpha=(c_i(R)-alpha*t+Q(x)-eta*s)/t.       (7)

Fix a finite cut length L. The supports of the lifted truncated chains lie
in a fixed compact subset of C^3 times[0,L], by (5), and their masses and
boundary masses are uniformly bounded by the preceding area estimate.
They are integral3-currents. Along any sequence R->infinity, the
Federer--Fleming compactness theorem gives a subsequence converging to an
integral3-current. We use the precise mass-and-boundary-mass hypothesis;
weak distributional compactness alone would not give a rectifiable limit.
See [Simon, Chapter6 Section3, Theorem3.15](https://math.stanford.edu/~lms/ntu-gmt-text.pdf).

Equation(7) shows that every such limit is supported in

    {P=alpha} times[0,L].                         (8)

Choose generic truncation levels so the boundary slices are integral;
semialgebraic stratification provides an arbitrarily large sequence of
such levels. A complex algebraic level set has a finite complex Whitney
stratification. At almost every point of the rectifiable limiting current,
its tangent plane lies in a stratum of (8). The pullback of the holomorphic
top form omega to each such stratum is zero: its x-projection lies in a
proper complex subvariety of C^3. This also handles singular strata;
strata of real dimension below3 carry no mass, and higher-dimensional
complex strata still annihilate a holomorphic3-form.

Test the weak limit with exp(-eta*s)*omega, multiplied by a fixed smooth
compact cutoff equal to1 on the common support. Every subsequential limit
of its integral is zero. Therefore

    exp(-c_i(R))*integral_(0<=s<=L) exp(S_t)*omega ->0. (9)

This is a statement about currents, not about the full real3-volume tending
to zero. A chain may retain real volume while its holomorphic volume form
vanishes in the limiting complex fiber.

First send R to infinity in (9), then L to infinity in (6). The conclusion
for every member of the flat cut basis is

    exp(-c_i(R))*integral_Gamma_i(R) exp(S_t)*omega ->0. (10)

This treats collisions at (m*T+l_T)*(T^2-4bc)=0 without replacing them by
generic Morse points. It gives no rate and uses no formal stationary-phase
splitting.

## 5. The supporting ray and the moment contradiction

For each bounded critical branch let epsilon=1/t and
F(epsilon)=P(x(epsilon))-epsilon*Q(x(epsilon)). Criticality gives
F'=-Q(x(epsilon)), which is bounded. Hence

    c_i(R)=alpha_i*t+O(1), alpha_i in (2).          (11)

Choose theta with Re(exp(i*theta)*alpha_i)<=0 for every alpha_i. Such a
theta exists precisely under the closed-half-plane condition. Equations
(10),(11) then show that **all seven scalar basis periods tend to zero on
this one ray**. Even when two actions are opposite and purely imaginary
on this ray, (10) applies; the cut angle eta remains independently available.

If every positive pure Gaussian moment vanished, Wick integration by parts
would define a horizontal formal functional with seed value1. On the
differential submodule generated by the seed its values are rational
(1,0,0,... on the successive derivatives). Fiber extension and the dual
ODE extend it to an actual horizontal functional at an ordinary t-point.
The retained perfect rapid-decay comparison then represents it by a flat
cycle. This is the [seed-submodule bridge](NODAL_NEWTON_EXCLUSION.md), which
requires neither a cyclic seed nor convergence of the original Wick series.
The fixed affine coordinate Jacobian only changes a nonzero normalization.

That flat cycle is a constant linear combination of the basis above along
the chosen ray. Its scalar period tends to zero, contradicting its nonzero
constant value. This proves the half-plane exclusion.

For k=0 there are at most two nonzero alpha values A,B. They always lie in
a closed half-plane through0. Explicitly, for nonzero A,B put z1=A/|A|,
z2=B/|B|. If z1+z2!=0, choose exp(i*theta)=-conj(z1+z2)/|z1+z2|.
If z1+z2=0, choose exp(i*theta)=i*conj(z1), giving both real parts zero.
The cases with a zero A or B are immediate. This proves the entire k=0
assertion, including nonzero Gaussian means and arbitrary b,c.

## Evidence, correction and remaining work

Run `.venv/bin/python research/moments/vanishing/verify_nodal_center_line_decay.py`.
The [receipt](nodal_center_line_decay_verification.json) checks the scaled
gradient identities, critical values, phase-band identity, regular-level
top-form annihilation, a Cauchy--Binet volume identity, exact tail majorants,
supporting-ray examples including opposite actions, and retained evidence
hashes. It verifies the shifted nodal prefix trap

    P=(T+tau)*uv+u^3+v^3,
    M1=0, M2=tau^2+31, M3=54*tau.

Thus tau^2=-31 gives two initial zeros and a nonzero third moment. The
Gaussian mean generated by the translation is essential to its normal-form
interpretation. The first checker version used structural equality for an
expanded versus factored expression; its rejected input and receipt remain.

The [source audit](nodal_center_line_source_audit.json) specifies the imported
theorems and reading coverage. The gradient argument, constrained cut-family
application, uniform volume estimate, current-limit argument and Gaussian
moment exclusion are internal written proofs. None is formally certified
by the symbolic checker.

The earlier critical-profile note called b=c=k=0 a completed homogeneous
subcase without explicitly requiring a centered Q. That wording was too
broad: an affine Gaussian mean is not covered by the old homogeneous proof.
This note supplies the stronger result and records the correction.

The method does not cover, for example, critical values
(-1-i,2-i,-1+2i), whose average is0 and whose triangle has positive area.
This is a boundary example for the method, **not** a pure-null polynomial
witness. The subsequent [dominance theorem](NODAL_MORSE_DOMINANCE.md)
closes the strict triangle and the finite Morse/distinct-value locus with
nonzero quartic gamma. Degenerate finite points, coincident finite values
and the quartic resonance remain open. The general cuspidal/reducible
leading cases, numerical pure gates, coefficient schemes and lower-power
GVC remain separate obligations.
