# Determinant 1092: an explicit rational rank-18 base change and the initial unlock

## Result and classification

**Verified application.** The existing norm-10 orbit 8044 bisection gives an
explicit degree-two map `C -> P1_t` over **Q**, with `C` a conic having a
displayed rational point. The pulled-back elliptic family has an independently
checked rank-at-least-18 subgroup. Its new anti-invariant section has Shioda
height 12. All equations, rational maps and the rank-18 height Gram are in the
[construction certificate v2](../../artifacts/generated-results/elliptic-curves/det1092_orbit8044_rank18_base_change_v2.json);
the [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_orbit8044_rank18_base_change_replay_v2.json)
checks the coefficient identities and geometric height hypotheses.

**Verified application / exact obstruction.** The bisection of V3's first
successful generic norm-10 centre, orbit 127449, is nonsplit at `t=0`; the
discovered rational point does not lie on it. A pointed degree-two chart on a
single elliptic fibre and a degree-two multisection of the surface are distinct
objects even when their trace section agrees.

The collection has exactly **two** members: previously constructed orbit 8044,
and one additional 19-by-20 Riemann--Roch calculation for the initial V3 centre.
Only the first is parametrized for a rank-18 construction. No point search,
parameter sweep, class-group computation, subgroup census or modification to
the eight-fibre V3 pilot belongs to this work.

**Later fixed-centre supplement, verified application.** The historical first
centre has now received one further RR calculation, also yielding an exact
nonsplitting obstruction at zero. This brings the examined collection to three;
the original two-member evidence above is unchanged. The separate
[two-cover note](DET1092_TWO_COVER_BRANCH_GATE_2026-09-08.md) supplies the later
parametrization of the V3-centre cover, without making its fibre at zero split.

## Audit of completed results

The following are **verified applications already established in this
repository**, reused here rather than recomputed as new theorems:

| Existing authority | Consequence used here |
|---|---|
| `EC-CURVE302-RECOVERED-MW17-PARENT` | Explicit 17 sections and exact height matrix `G`, determinant 1092; `t=0` is literally 302. |
| `EC-CURVE302-PARENT-GEOMETRIC-PICARD19` | Full geometric MW rank is exactly 17; the elliptic K3 has 24 `I1` fibres. A constant-field extension cannot add an independent section. |
| `EC-CURVE302-PARENT-QUADRATIC-DESCENT-GATE` | The original surface cannot descend from a pointed quadratic source of ranks 1 through 16. This does not obstruct a new base change of the original surface. |
| `EC-CURVE302-PARENT-DEGREE2-MULTISECTION-QUOTIENT` | Complete 131072-class quotient; 40917 geometrically rational norm-10 bisection classes. Geometric genus zero alone does not assert a Q-rational point. |
| `EC-CURVE302-PARENT-CHEAPEST-LATTICE-BISECTION` | Orbit 8044 has an exact RR equation; its fibre at zero is nonsplit. Its conic parametrization and lifted rank-18 group were not previously supplied. |

Canonical inputs are the [parent note](CURVE302_RECOVERED_MW17_PARENT_2026-09-07.md),
[geometric/descent note](CURVE302_PARENT_SEARCH_AND_GEOMETRY_2026-09-07.md), and
[multisection note](CURVE302_LOW_DEGREE_MULTISECTIONS_2026-09-07.md).

**Established literature.** Base change, the intersection height pairing, its
positive definiteness modulo torsion, and degree scaling are used from
[Schuett--Shioda, *Elliptic Surfaces*, §§5 and 11](https://arxiv.org/abs/0907.0298).
The rank calculation below is an application to explicit checked equations.
It is not a claim that the general construction principle is new literature.

## The exact initial unlocks

**Verified application.** The [initial-unlock certificate](../../artifacts/generated-results/elliptic-curves/det1092_initial_unlock_construction_audit_v1.json)
reconstructs both the earliest historical gain and the first gain of the completed
autonomous V3 calibration. This avoids conflating their different frozen policies.

| Run | First gaining chart, one-based | Rational search coordinate | Certified transition |
|---|---:|---:|---|
| Historical generic-only 49-chart wave | 8 | `-1714/2373` | `17 ->18`; the whole fixed wave later reaches 19 |
| Completed autonomous V3 | 17 | `848/2073` | `17 ->18`, followed by the frozen rebuild |

Both start with precisely the 17 specialized generic sections. The literal
specializations are checked by the [second RR audit](../../artifacts/generated-results/elliptic-curves/det1092_initial_unlock_bisection_obstruction_v1.json),
using the original parent formulas and the exact short-model change
`X=36x+15`, `Y=108(2y+x+1)`.

In that ordered generic basis the first V3 centre is

```text
w = (0,1,0,2,0,3,0,-1,-1,0,-2,2,2,-4,-1,2,-1),
w^t G w = 10; orbit 127449, canonical lane.
```

The certificate displays its exact rational `(Xc,Yc)`. Its searched quartic is
`r^2 = sum(c_j n^j d^(4-j))`, with ascending coefficients

```text
4215832226801039816080943064392704
-73822486035349620599009516655116928
277634842813457827614924531846985476
12894706083225824080517670544535148
-326526540600688922861386593383871
```

The exact homogeneous witness is

```text
(n,d,r) = (848,2073,627903986691316782763008).
```

Let `[a,b;c,e]` and `k` be the exact raw-slope matrix and ordinate factor in
the certificate. The verified map is

\[
 \ell=\frac{an+bd}{cn+ed},\qquad
 v=\frac{kr}{(cn+ed)^2},\qquad
 X=\frac{\ell^2-X_c-v}{2},\quad Y=\ell(X-X_c)-Y_c.
\]

It returns

```text
X = -145408491449074952644831800869259116405301/3243601
Y = 19103023468249730422176006737246636039677568252645668448005632/5841725401.
```

The original-302 coordinates are obtained by the inverse change above and
are also explicitly recorded. Independently enumerated finite groups at 18
fixed primes through 197 give binary column rank 18 on these 17 generic
points plus the new point. The 2-division cubic has no root modulo 31, so
infinite descent proves integral independence. This rechecks the first gain
without trusting a floating height matrix or running a rational point search.

Before this gain, the centre, its rational chart and all 32 generic anchor
cosets were already available from generic data. After the gain, the new
point and new parity-extension cosets become available. The later adaptive
landscapes use those exceptional discoveries. They were not inputs to the
first chart.

**Verified application: unchanged control.** The completed native11952 V3
positive control starts at 17, scores 32 anchor cosets and completes all 82
scheduled charts, ending at 17 with independent replay. Thus it never reaches
a positive extension dimension. The 302 initial success and 11952 initial
miss measure bootstrap recovery under finite exposure. Different fibres and
parents prevent a controlled attribution to rank incidence or visibility.
In particular the 11952 miss does not imply that rational multisections or
extra specialized points are absent. Neither search was tuned or rerun here.

## Why the winning chart does not itself provide a rational base change through 302

**Verified application / exact obstruction.** For the V3 word `w` above,
construct `P_{-w}` and solve the existing RR system

\[
 f_0(t)+f_1(t)x+f_2(t)y=0,\qquad
 (\deg f_0,\deg f_1,\deg f_2)\le(9,5,3).
\]

The matrix has rank 19. Removing `P_{-w}` leaves its unique rational bisection,
with trace `P_w`. The stored quadratic at `t=0` has nonsquare discriminant;
direct substitution of the first V3 point into the line has a nonzero rational
value. Both exact witnesses are retained in the second RR audit.

Consequently that particular rational bisection cannot explain the first
rational302 gain. The point chart searches varying slopes on the specialized
elliptic curve; the RR bisection imposes one specific slope as a rational
function of the family parameter. The shared norm/parity label does not make
these conditions equivalent. Higher-genus pencils or other multisections
remain possible explanations; no conclusion about them is established here.

## Historical first-centre obstruction: a limited negative result

**Verified application.** The actual historical first chart has generic word

```text
w = (1,-1,-2,1,0,-1,1,-1,1,-1,1,0,1,0,0,-1,0).
```

Its norm is 10. The new
[constructor](../cas/construct_det1092_historical_unlock_obstruction.sage)
uses only this fixed word and the generic parent/17 sections. The centre is
selected retrospectively from the historical run; this is not claimed as a
prospective selection experiment. Neither the exceptional point nor its search
coordinate is an input to constructing the equation.

The unique RR line in the same degree bounds (9,5,3) gives a residual
`a(t)x^2+b(t)x+c(t)` with

\[
b^2-4ac=h(t)^2q(t),\qquad \deg q=2,\quad \gcd(q,q')=1.
\]

All coefficients are explicit in the
[historical obstruction certificate](../../artifacts/generated-results/elliptic-curves/det1092_historical_unlock_obstruction_v1.json).
Off zeros and poles of `a,f2,h` and the elliptic discriminant, the exact
splitting condition is `q(t)` square in Q. The associated curve
`s^2=q(t)` is a degree-two cover with smooth projective geometric genus zero.
Its rational solubility is not asserted by this supplement. The maps are

\[
x=\frac{-b+hs}{2a},\qquad y=-\frac{f_0+f_1x}{f_2}.
\]

At zero all displayed nondegeneracy conditions hold, but the *certified
normalization* of this quadratic satisfies

\[
q(0)\equiv17\pmod{23},\qquad17^{11}\equiv-1\pmod{23}.
\]

Thus this cover has no Q-rational point above zero. The line does not contain
the elliptic point at infinity there (`f2(0) != 0`), so an omitted affine chart
does not evade the obstruction. The historical first exceptional point cannot
be supplied by this particular residual multisection.

The [independent checker](../cas/verify_det1092_first_unlock_obstructions.sage)
does not import the RR constructor. It independently clears rational-function
denominators, proves the interpolation rank 19, checks the line and trace-removal
identities, and verifies nonsplitting for both historical and V3 first centres.
Its [replay certificate](../../artifacts/generated-results/elliptic-curves/det1092_first_unlock_obstructions_replay_v1.json)
passes. Both bounded calculations finished in seconds, with no detached job.

**Scope of the conclusion.** This met alternative 4 of the earlier goal:
an exact obstruction to a specified proposed source of the first 302 direction.
It does **not** produce an affirmative splitting condition satisfied at zero,
identify the actual arithmetic source, exclude every multisection in a parity
orbit, or distinguish 302 from all eight pilot fibres. These remain unknown.
The separate rank-18 rational base change above misses zero and is not offered
as an explanation of that first point. No new point search or V3 change follows.

The user has since reopened the broader goal of understanding302's mechanism;
these limited obstructions do not complete that goal. The subsequent
[RR-net calculation](DET1092_FIRST_UNLOCK_RR_NET_2026-09-08.md) constructs the
next three-dimensional system and verifies its full pointed-chart restriction.

Replay command from the repository root:

```sh
sage -python research/elliptic-curves/cas/verify_det1092_first_unlock_obstructions.sage
```

## Explicit rank-18 construction over Q

**Verified application.** Reuse orbit 8044 with

```text
w = (0,-1,0,0,0,0,0,0,0,1,0,0,0,-1,0,0,0).
```

Its existing residual equation is `a(t)x^2+b(t)x+c(t)=0`, accompanied by
the RR line `f0+f1*x+f2*y=0`. The coefficients are supplied verbatim in the
new certificate. Exact factorization gives

\[
 b(t)^2-4a(t)c(t)=h(t)^2q(t),\quad q(t)=\alpha t^2+\beta t+\gamma,
\]

with the rational constant square factor retained and

```text
alpha = 65787359279227310322979727123145990962151271221560127726834795993231289
beta  = 5407077334757298104342553746936743632103681815165706562836351116238000
gamma = 127598680603895545533594201587042450373210388891137967376445049000000.
```

Take `C: s^2=q(t)`. Its discriminant is nonzero, so its smooth projective
model has genus zero. It has the exact rational point

```text
t0 = -2502544513254536830634997737000/7142333315568558691744661205247
s0 = 567337156389201179805046395842070887694530746377558132477366600000/7142333315568558691744661205247.
```

The line of slope `u` through this point gives the rational parametrization

\[
 t(u)=t_0+\frac{q'(t_0)-2s_0u}{u^2-\alpha},\qquad
 s(u)=s_0+u\bigl(t(u)-t_0\bigr).
\]

The inverse is `u=(s-s0)/(t-t0)` where defined. Exact substitution verifies
`s(u)^2=q(t(u))`, and `t(u)` has degree two. This constructs infinitely many
rational points on `C` and infinitely many rational family parameters, after
removing the finite excluded set. No number-field extension is required.

The pulled-back family is exactly `E(t(u))`, with the 17 old sections composed
with `t(u)`, and the additional section

\[
 x_Q=\frac{-b(t)+h(t)s}{2a(t)},\qquad
 y_Q=-\frac{f_0(t)+f_1(t)x_Q}{f_2(t)}.
\]

The certificate gives these functions explicitly in `u`, together with a
polynomial Weierstrass model and polynomial coordinates for `Q`. Every
coefficient identity is independently checked.

**Exact splitting condition.** Off the recorded degree-29 excluded-parameter
polynomial and the affine infinity chart, the residual fibre splits over Q
if and only if `q(t)` is a rational square. The two signs of `s` give its two
points. The excluded polynomial contains the surface discriminant, branch
polynomial, trace-removal degeneracies and every displayed denominator.
Ramified fibres require separate evaluation; a rational root of `q` gives
one rational conic point. All rational `u` outside the finite exceptional
set supply a split fibre through the displayed parametrization.

At zero, `q(0)=gamma` remains nonsquare. That fibre lifts only after adjoining
its quadratic square root. This local fact coexists with a parametrization
of the full conic over Q; it does not make the rank-18 construction a
constant-number-field construction.

## Independence proof

**Established literature applied to verified data.** The quadratic branch
polynomial is coprime to the parent's degree-24 discriminant. Both branch
fibres are smooth, so the base-changed minimal surface has 48 `I1` fibres
and Euler characteristic `chi=4`.

Write `D(u)=den(t(u))`. The checker verifies that `D^4*x_Q` and `D^6*y_Q`
are polynomials of degrees at most 8 and 12. On the global minimal model
this proves `Q.O=0`, including at infinity. There are no reducible-fibre
corrections, so the height formula gives

\[
 \langle Q,Q\rangle=2\chi+2Q.O=8.
\]

The deck conjugate is obtained by `s -> -s`. The elliptic group law verifies
coefficientwise that

\[
 Q+\sigma Q=P_w.
\]

The old trace has height 10 on the parent, hence height 20 after the
degree-two base change. Thus

\[
 \langle Q,\sigma Q\rangle=2,\qquad
 \langle Q-\sigma Q,Q-\sigma Q\rangle=12>0.
\]

The inherited 17-dimensional space is deck invariant. Applying `1-sigma`
to a putative dependence forces the coefficient of `Q` to vanish, since
`Q-sigma Q` has positive height; inherited independence then finishes the
proof. Equivalently the explicit 18-section height Gram is

\[
 H=\begin{pmatrix}2G&Gw\\w^tG&8\end{pmatrix},\qquad
 8-(Gw)^t(2G)^{-1}(Gw)=3,
 \quad\det H=3\,2^{17}\,1092=429391872>0.
\]

This proves `rank E(Q(C)) = rank E(t(u))(Q(u)) >=18`. It does not assign rank
18 to every rational specialization or prove the full generic rank is 18.

**Subsequent verified application.** The
[conic seed progression](DET1092_CONIC_SEED_PROGRESSION_2026-09-08.md)
now certifies the fixed `u=0` specialization and every `u=B*n`, `n in Z`,
for an explicit squarefree integer `B`. This gives an equation-only
infinite rank18 seed family with no unresolved exceptional set on that
progression. It still misses302 and all eight controls and does not
strengthen the claim to every rational `u`.

**New deduction, elementary.** Replacing this lifted section by a translate
`Q+P_i` does not provide another independent direction: modulo the inherited
subgroup it is still `Q`. Even the two sheets give only one extra direction,
because `Q+sigma Q=P_w`. This explicit trace relation supplies the relevant
independence obstruction before any attempt to combine covers. No composite
cover or several-direction claim is made in this first theorem.

## Checkpoints, correction, and replay

The construction and checks each completed in small bounded runs; no detached
job is left running. The 45-second per-process cap was used for these checks.
The construction certifies its quadratic cover, trace identity and height
Gram in separate printed checkpoints. The additional RR construction is
limited to one orbit and one 19-by-20 system.

The first local certificate v1 is retained as rejected evidence: it used
Sage's `squarefree_part` for the exclusion polynomial, which keeps odd
valuations rather than all distinct factors. The independent checker caught
the lost exceptional factors. Version 2 uses `f/gcd(f,f')`; all exclusion
divisibility checks pass. The cover, section and height calculations were
unchanged. Only v2 and its independent replay support the theorem here.

From the repository root:

```sh
sage -python research/elliptic-curves/cas/verify_det1092_rank18_base_change.sage
sage -python research/elliptic-curves/cas/audit_det1092_initial_unlock.sage
```

The construction producer
`construct_det1092_rank18_base_change_v2.sage` and RR producer
`audit_det1092_unlock_bisection.sage` refuse to overwrite their certificates.
All four output claims are labelled **verified application**; the general
height tools are labelled **established literature**, and the translated-sheet
independence observation is labelled **new deduction**. No conjecture is
promoted to a result.
