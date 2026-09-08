# An equation-only infinite family of certified determinant1092 seeds

## Result and boundary

**New constructive deduction, independently verified.** Let `t(u)` and
`Q(u)` be the already constructed rational parametrization and elliptic
section of the orbit8044 conic. Put

\[
\boxed{\mathcal B=6437157837635645727368063702636659.}
\]

For **every integer** `n`, the fibre at `t(B*n)` is smooth over `Q`, and
the17 inherited sections together with `Q(B*n)` are linearly independent.
Thus

\[
\boxed{Q(\mathcal B n)\notin H_{t(\mathcal B n)}\otimes\mathbf Q,
\qquad \operatorname{rank}E_{t(\mathcal B n)}(\mathbf Q)\ge18.}
\]

There is no unresolved finite specialization-exception set on this
progression. The map `t(u)` has degree two and no pole on it, so these
give infinitely many distinct rational parameters. All equations and
points are defined over `Q`; no number-field extension is needed.

This is an oracle-free seed producer on an explicit infinite subfamily.
It is **not** a reconstruction of the historical302 seed: the same conic
is nonsplit at302 and all eight unchanged controls. Their other seeds and
ranks are not excluded. No amplifier was run, no pilot was changed, and
no later discovered point or V3 artifact was read as a construction input.

## What is new, and what was already proved?

**Verified application already established.** The
[original conic construction](DET1092_RANK18_BASE_CHANGE_AND_INITIAL_UNLOCK_2026-09-08.md#explicit-rank-18-construction-over-q)
provides the quadratic equation, rational parametrization, elliptic lift,
trace relation and positive18-section function-field height Gram. It
proves rank at least18 over `Q(u)`, but did not certify a rational
specialization or identify its exceptional set. Those calculations are
reused, not presented as new.

**New deduction in this application.** One fixed formula evaluation at
`u=0` gives an independently certified rank18 subgroup. An explicit
congruence argument preserves its entire finite-group certificate for
every `u=B*n`. The exported routine also recognizes input parent parameters
in this image. No parameter is selected using catalogue ranks, exceptional
coordinates, local landscape scores or a point search.

## Explicit construction and parameter criterion

**Verified application.** The existing conic is

\[
C:s^2=q(t)=\alpha t^2+\beta t+\gamma,
\]

with its exact coefficients and rational point `(t_c,s_c)` recorded in
the [original equation certificate](../../artifacts/generated-results/elliptic-curves/det1092_orbit8044_rank18_base_change_v2.json).
These choices are unchanged. Its line parametrization is

\[
t(u)=t_c+\frac{q'(t_c)-2s_cu}{u^2-\alpha},\qquad
s(u)=s_c+u(t(u)-t_c).
\]

The already fixed RR residual equation `a(t)x^2+b(t)x+c(t)=0` and line
`f0(t)+f1(t)x+f2(t)y=0` have discriminant `h(t)^2 q(t)`. Their lift is

\[
x_Q(u)=\frac{-b(t(u))+h(t(u))s(u)}{2a(t(u))},\qquad
y_Q(u)=-\frac{f_0(t(u))+f_1(t(u))x_Q(u)}{f_2(t(u))}.
\]

Use the reduced rational maps in the certificate, which remove apparent
chart cancellations. The independent replay verifies the conic identity,
the RR equation, the line equation and the elliptic identities of all18
points as rational functions, before proving uniform specialization.

**Exact sufficient criterion for an input parent parameter `tau`.** Solve
the degree-at-most-two equation

\[
\operatorname{num}(t(u))-\tau\operatorname{den}(t(u))=0.
\]

If a rational root `u` has nonzero denominator and `u/B` is an integer,
the routine returns the rational elliptic point `Q(u)`, its17 generic
companions, and the uniform rank18 certificate reference. Otherwise it
returns `UNKNOWN_OUTSIDE_PROVED_PROGRESSION_IMAGE`, not a rank upper bound
or a seed-nonexistence claim. This criterion is sufficient, not necessary
for an exceptional direction on the fibre.

The first constructed parameter is

\[
\tau_*=
\frac{1636579823608344775969579426883491310836441501800610644096637979339414901424849081843284973507691000}
{6102275947055890661075259641569164647354846459461782772324094998826318153155853546790385417152874979}.
\]

Its [seed bundle](../../artifacts/generated-results/elliptic-curves/det1092_conic_seed_progression_v1/seed-at-n0-v2.json)
contains the exact elliptic equation,17 inherited points, the constructed
extra point and the conic point. This is a subgroup rank lower bound;
neither the full rank nor catalogue/literature novelty was investigated.

## Independent rank certificate at the fixed address

**Verified computation.** The protocol froze `u=0` and all odd primes at
most1009 before finite rank tests. It made no replacement or retry at a
different family parameter. Of168 prime exposures,156 gave good cubic
Kummer reductions, seven were bad reduction, and five had equation
denominators. Every exposure and skipped reason is retained.

The constructed matrix has generic rank17 and rank18 with the conic point.
The independent checker does not import that constructor or its cubic
character calculation. It enumerates complete finite groups and their
doubled subgroups at the selected proof primes, using elementary integer
addition on `Y^2=X^3+aX^2+bX+c`. It independently recomputes the quotient
coordinates of the18 rational points and verifies matrix rank18.

The proof primes, whose product is `B`, are

```text
19,23,43,47,53,59,67,71,79,83,101,107,109,127,131,149,157,163.
```

At23 the finite group has order33. Good reduction therefore excludes
rational2-torsion. Full18-column rank in the product of finite quotients
by doubling then proves that the18 displayed rational points are
independent, rather than merely bounding rank plus2-torsion.

The [independent certificate](../../artifacts/generated-results/elliptic-curves/det1092_conic_seed_progression_v1/replay.json)
retains the complete finite point sets, doubled subgroups, coset
representatives, point reductions and the independently recomputed matrix
blocks. These finite-group enumerations are proofs, not rational-point
searches.

## Why the certificate holds for every integer n

**New deduction, with exact polynomial checks.** Write each pulled-back
Weierstrass coefficient as `A_i(u)/D_i(u)` with integral polynomial
numerator and denominator. Write every point as a primitive integral
projective triple, and the base map as an integral polynomial pair.
The checker proves, at every selected prime `p`, that:

- every coefficient denominator and the base-map denominator is a unit
  at `u=0`;
- every point triple is nonzero modulo `p` at `u=0`;
- the elliptic equation at `u=0` has good reduction.

For `u=B*n` the polynomial evaluations have exactly those same residues.
Hence every denominator remains nonzero, the curve remains smooth, and
all18 points have exactly the certified reductions. The odd-order
reduction at23 remains unchanged as well. The rank18 certificate therefore
holds for every integer `n`, without invoking an unspecified exceptional
set in a specialization theorem.

The base-map denominator cannot vanish on this progression: a zero would
contradict its nonzero reduction at any proof prime. Since `t(u)` is
nonconstant of degree two, each rational parameter has at most two
preimages. The progression thus gives infinitely many distinct fibres.

**Established literature, applied to verified rational points.** In the
monic cubic model, each seed has Kummer representative `X_Q-theta` with
square norm `Y_Q^2`. Its class is non-generic by the certificate and is
rational, so it lies in Selmer with zero Sha image. This uses the usual
[Kummer–Selmer–Sha exact sequence](https://www.mathe2.uni-bayreuth.de/stoll/papers/Explicit-Descent-III.pdf),
not a computation of a Selmer dimension. No linear map from an RR
Jacobian to this elliptic curve is asserted.

## Controls and remaining302 problem

**Verified application.** Exact integer square-root bounds prove that
`q(tau)` is nonsquare at every address below:

```text
scale-0131232, scale-0257585, scale-0487239, scale-0177036,
scale-0043332, scale-0590501, scale-0290097, scale-0748009,
302 (t=0).
```

Thus this producer has no coverage on the frozen comparison panel. All
nine calls return the explicit noncoverage/UNKNOWN status. This is not
evidence that those elliptic fibres have no exceptional points;302 is
already a counterexample to that inference.

**Unresolved.** This construction supplies prospectively available seeds
on other determinant1092 fibres, but it does not remove the calibration
from the historical302 RR-net class or characterize why302 unlocks.
The broader goal therefore remains open. No claim is made that these
newly seeded fibres will amplify beyond18, or that their height/cost is
competitive with the existing pilot. No amplification is part of this work.

**New verified obstruction.** The
[trace/twist Kummer identity](DET1092_TRACE_TWIST_KUMMER_OBSTRUCTION_2026-09-08.md)
shows that the conic's anti-invariant point, rational on the quadratic
twist even at nonsplit fibres, has inherited2-Kummer class under the
natural torsion-module identification. This class is rational/Selmer with
zero Sha image on302 and all eight controls. It cannot itself distinguish
the first seed. This does not contradict independence of the conic branch
on the proved split progression, or exclude other independent points.

The two sheets are not two extra independent directions: the existing
identity `Q+sigma(Q)=P_w` puts both in one quotient direction modulo the
generic group. No subsequent exceptional direction or combined-cover
claim is made.

## Reproduction and input boundary

**Verified computation.** Constructor, independent proof replay and producer
interface tests each finished in under one second with25-second caps.
There are no detached jobs. No class group, number-field calculation,
new parameter sweep, later exceptional point, V3 artifact or pilot
modification was used.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_conic_seed_progression.sage
sage -python research/elliptic-curves/cas/verify_det1092_conic_seed_api.sage
sage -python research/elliptic-curves/cas/det1092_conic_seed_v2.sage --n 0
```

For a supplied rational parent parameter:

```sh
sage -python research/elliptic-curves/cas/det1092_conic_seed_v2.sage --parameter 'NUMERATOR/DENOMINATOR'
```

The [producer](../cas/det1092_conic_seed_v2.sage) evaluates only the generic
parent and conic maps. It does not load a saved exceptional point as an
execution oracle. The modulus is a proof output from the single
formula-constructed specialization, not a point-selection score. Its
current maps are hash-bound to the independently verified proof. The
initial producer and correct outputs are retained; version2 only makes
that transitive input binding explicit.
