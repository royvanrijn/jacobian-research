# A Euclidean formula for norm-ten bisections; the15-orbit panel is closed

## Result and scope

**New constructive deduction/application.** The norm-ten RR construction
needs neither a19-by20 matrix solve nor polynomial factorization. A single
inverse modulo the square of the trace pole polynomial constructs the
unique line, its branch quadratic, and polynomial maps back to the elliptic
surface. A sixth-order cancellation is forced by the elliptic equation.

**Independently verified application.** The three remaining masks61,107,111
of the frozen panel all give nonsplit fibres over302. Each also misses all
eight unchanged controls:27 exact rational nonsquare certificates. Combining
these with the earlier proofs closes **all fifteen identified orbits** of
that panel at302, including every generic section translation. It does not
close the40,917-orbit atlas. No extra orbit or prime was admitted.

The three new curves have explicit equations and maps, but give **no new
rational point on302 or the controls**. All27 tested values are positive;
their obstruction is not real solubility. No class/unit group, parameter
sweep, point search, later exceptional point, production change or detached
job was used. The construction reuses only saved generic trace points;
the independent checker reconstructs their words exactly over `Q(t)`.

## 1. Audit and relation to existing constructions

**Completed inputs.** The [degree-two census](CURVE302_LOW_DEGREE_MULTISECTIONS_2026-09-07.md)
identifies minimum norm-ten classes with the40,917 smooth rational bisection
orbits. The [index theorem](DET1092_RATIONAL_BISECTION_INDEX_AND_ODD_DIVISOR_CONSTRUCTION_2026-09-09.md)
proves their rationality over `Q`. The old norm-eight pencil already uses
principal-part cancellation, and the existing RR interface constructs the
norm-ten curves by a linear system.

What is new here is the explicit zero-free-parameter Euclidean specialization
for norm ten, a proof of all divisions and the sixth-order cancellation,
and exact resolution of the three previously surviving orbits. This is an
application of the existing pointed-quartic machinery, not a claim that the
general chord construction or polynomial Euclidean algorithm is new.

## 2. Explicit formula from a generic trace

Use the short parent equation

\[
 E_t:Y^2=X^3+A(t)X+B(t),\qquad\deg A\le8,\quad\deg B\le12.
\]

Let `w` be a minimum norm-ten word and write its negative trace as

\[
 P_{-w}=(N_x/h^2,N_y/h^3),\quad h\text{ monic},\quad
 \deg h=3,\ \deg N_x\le10,\ \deg N_y\le15,\quad (N_x,h)=1.
                                                               \tag{1}
\]

If an intersection with `O` is at infinity, first change base chart; the
[four-chart result](DET1092_UNIFORM_SMOOTH_BISECTION_REDUCTION_2026-09-09.md)
supplies a deterministic remedy at its certified primes. The same
intersection-degree argument supplies a chart in characteristic zero.
The three saved traces here already have degree-three finite pole polynomials.

**New explicit construction.** Take the unique degree-below-six remainder

\[
 m=-N_yN_x^{-1}\pmod{h^2},\qquad
 g=\frac{mN_x+N_y}{h^2},\quad
 b=\frac{m^2-N_x}{h^2},\quad
 k=\frac{mb-2g}{h^2}.                                      \tag{2}
\]

Then all these expressions are polynomials. Their degree bounds are
`deg(m,g,b,k)<=(5,9,4,3)`. Define

\[
 \boxed{q(t)=\frac{4mk-3b^2-4A}{h^2}
 =\frac{m^4-6N_xm^2-8N_ym-3N_x^2-4Ah^4}{h^6}.}            \tag{3}
\]

This is a polynomial of degree at most two. The desired curve and maps are

\[
 \boxed{C_w:W^2=q(t),\qquad
 X=\frac{b+hW}{2},\qquad Y=-\frac{hk+mW}{2}.}               \tag{4}
\]

Its line in `H0(3O+9F)` is simply

\[
 -g+mX+hY=0.                                               \tag{5}
\]

To return to the original Weierstrass model, use
`x=X-b2/12`, `y=Y-(a1*x+a3)/2`, with the fixed parent's original invariants.
No exceptional coordinate enters (1)--(5).

## 3. Why every division is exact

**New algebraic deduction.** Put `N=N_x`, `V=N_y`. The trace equation is

\[
 V^2=N^3+ANh^4+Bh^6,\qquad V=-mN+h^2g.                   \tag{6}
\]

Since `N` is a unit modulo `h^2`, (6) first gives `m^2=N mod h^2`, proving
that `b` is polynomial and also `(m,h)=1`. Substituting into (6) gives

\[
 N^2b-2mNg=h^2(AN+Bh^2-g^2).                              \tag{7}
\]

Thus `Nb-2mg` is divisible by `h^2`. Reducing this again modulo `h^2`
gives `m(mb-2g)=0`; invertibility of `m` proves that `k` is polynomial.

For a short proof of the final cancellation set

\[
 c=(Nb-2mg)/h^2-A,\qquad Nc=Bh^2-g^2.
\]

Modulo `h^2`,

\[
 N(b^2+4c)\equiv m^2b^2-4g^2
             =(mb-2g)(mb+2g)\equiv0.
\]

Therefore `b^2+4c` is divisible by `h^2`. Direct substitution identifies
the resulting quotient with both expressions in (3). This proves the
sixth-order cancellation without discovering or factoring any polynomial.
The degree bound follows from numerator degree at most20 minus `deg h^6=18`.

### Why this is the complete RR line, not an ansatz

Any line with coefficient bounds `(9,5,3)` vanishing on the trace satisfies

\[
 f_0h^3+f_1Nh+f_2V=0.
\]

Reducing modulo `h` forces `h|f2`; degree at most three then gives
`f2=lambda*h`. If `lambda=0`, invertibility of `N` modulo `h^2` and
`deg f1<6` force `f1=f0=0`. Otherwise normalize `lambda=1`. The remaining
equation forces exactly `f1=m` and `f0=-g` from (2). Thus (5) is the unique
nonzero RR relation. A maximal-pole trace already forces rank19; no
separate nullspace or rank calculation is necessary.

The same argument works after admissible finite-field reduction. Once the
previous Gauss lifting lemma gives integral `h,N,V` and the full reduced
pole degree, `N` is a unit in `Z_p[t]/(h^2)`. Its inverse and the monic
polynomial divisions lift integrally. The Euclidean construction therefore
also replaces the earlier large-minor kernel certificate; it does not
remove the need to check admissible trace reduction.

### An independently checked polynomial normal form

Equations (2),(3),(6) imply

\[
 A=mk-\frac34b^2-\frac14h^2q,
\qquad
 B=\frac{m^2q+b^3-2mbk+h^2(k^2-bq)}4.                     \tag{8}
\]

Conversely these identities, with
`N=m^2-h^2*b`, `V=-m^3+3*h^2*m*b/2-h^4*k/2`, give the trace equation and
the lift (4) by universal polynomial substitution. The companion checker
verifies them in a polynomial ring over `Q`, not at a parameter sample.
The degree tuple `(h,m,b,k,q)<=(3,5,4,3,2)` gives precisely the degree8/12
K3 coefficient bounds. No new determinant1092 family is inferred by
varying these formal coefficients: its fixed `A,B` and17 sections remain
constraints.

## 4. Exact incidence and independence are separate

**Verified application of completed geometry.** The residual curve of (5)
has class `2O+4F+phi(w)` and is the unique smooth rational bisection of
that class. Its generic irreducibility forces `q` to be a nonsquare with
squarefree degree one or two. All three new examples have degree two.
The [all-prime multisection theorem](DET1092_ALL_PRIME_DIVISION_AND_MULTISECTION_THEOREM_2026-09-09.md)
therefore makes (4) independent of the inherited generic span over
`Q(C_w)`. Each separately gives generic rank at least18 after that base
change; no joint independence or rank18 on a specified fibre follows.

For every finite smooth original fibre, polynomial maps (4) extend even
when `h(t0)=0`. They lose no branch: choose polynomial Bezout coefficients
`U*h+V*m=1`. The inverse on the image is

\[
 W=U(2X-b)-V(2Y+hk).                                      \tag{9}
\]

Hence a rational fibre point from this curve exists exactly when
`q(t0)` is a rational square, including the ramified zero case. If
`q(t0)=0`, the two branches coincide and their inherited trace satisfies
`2Q=P_w(t0)`. Such a point cannot add a rational MW direction. For a
nonzero square, the two branches are distinct but may still be inherited;
the [halving-or-cycle classifier](DET1092_SPLIT_SEED_HALVING_DICHOTOMY_2026-09-08.md)
decides that remaining question under its certified generic-footprint
hypotheses.

For a *fixed conic and fixed parameter*, this incidence problem is an
ordinary rational-square test, not a genus-one Sha obstruction. A positive
nonsquare rational number has an odd valuation at some finite prime, even
if no such prime was factored or sampled. Passing149,151 and the real place
does not amount to everywhere local solubility. This observation does not
resolve the different, positive-genus carrier/Jacobian-class problem.

## 5. Exact outcome on the frozen panel

**Verified application.** The new computation takes exactly61,107,111,
the old unresolved masks; the nine target/control addresses are unchanged.
It constructs (2)--(4) from each saved generic trace, then evaluates `q(t)`
only on that roster. Every one of the27 values is positive and nonsquare.

For a reduced positive value `q=n/d`, each record supplies nonnegative
integers `a,b` such that

\[
 a^2\le n<(a+1)^2,\qquad b^2\le d<(b+1)^2,
 \qquad a^2\ne n\ \text{or}\ b^2\ne d.
\]

These integer inequalities certify nonsquareness without prime factors.
The checker verifies them directly and does not call a rational square-test
routine. It also reconstructs all three generic trace words by its own
reversed-order rational-function group law, verifies the entire map in
`Q(t)[W]/(W^2-q)`, and compares the lines with their old151 reductions.

The closure ledger is:

| Orbits | Evidence at302 |
| --- | --- |
|87,103,109,110,115,117,121,122,123|Earlier independently replayed modular exclusions|
|61,107,111|New exact rational nonsquare certificates|
|8044,47755,103186|Earlier exact conic/control certificates|

All fifteen miss302, including every generic translation. The new three
and the old three explicit-conic regressions also miss the eight controls;
the other nine were not newly evaluated on those controls. Historical
UNKNOWN outcomes remain intact in their original protocols and are
superseded only by this explicit closure record.

This remains a small, coordinate-selected panel. It neither excludes all
smooth bisections at302 nor restricts the known seed to this atlas. It
does not contradict the positive conic seeds at other parameters, the
dependent conic split, or the nine inherited MW16-to17 recoveries already
classified by the existing checker. No subsequent rank direction was used.

## 6. Reproduction and next boundary

The constructor took0.046 seconds internally. The independent rational-word
replay and universal symbolic check completed within their25-second caps.
No long-running process was started. The expensive generic traces were
already saved; this timing is not a promise for computing all40,917 traces.

- [Protocol and frozen roster](../../artifacts/generated-results/elliptic-curves/det1092_euclidean_bisection_survivors_v1/protocol.json)
- [Three equations and27 outcomes](../../artifacts/generated-results/elliptic-curves/det1092_euclidean_bisection_survivors_v1/summary.json)
- [Independent replay and15-orbit closure](../../artifacts/generated-results/elliptic-curves/det1092_euclidean_bisection_survivors_v1/independent-replay.json)
- [Universal normal-form identity](../../artifacts/generated-results/elliptic-curves/det1092_euclidean_bisection_survivors_v1/universal-identity.json)

```sh
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_euclidean_bisection_survivors.sage
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_euclidean_bisection_identity.sage
```

Main checker SHA256:
`c1d0bd03fbe7c8cec6dbe71648cf60567c9e86943148f4d4cec5433236025231`.
Main replay SHA256:
`075339b10b39f2d714678182df81881ee2e2563210ab6c88c9ca4e73e789ad2b`.
The proof is not formally verified or externally reviewed.

The prospective302 seed criterion is still open. The direct formula removes
an equation-construction bottleneck; it does not select a successful
remaining orbit. No full-atlas campaign is authorized or launched by this
checkpoint.

## 7. Cold positive/dependent admission preflight

**Verified application, not a new theorem or new rank result.** A separate
three-address replay now starts from the original generic parent and the
table word for orbit8044, without reading saved trace, conic or seed
coordinates. Its positive and dependent addresses are fixed historical
regressions, not prospectively selected parameters. The reduced-parameter
chart is used only to transport the addresses to the original parent.

The reconstructed trace word is `w=-e2+e10-e14`. Equations (2)--(4) construct
the conic. At each address a finite quotient footprint is frozen using
only the seventeen generic sections, before evaluating the square condition.
When the value is a nonzero square, the nonnegative rational root determines
one branch, which is passed to the existing halving-or-cycle classifier.

| Fixed address | Incidence and exact admission result |
| --- | --- |
|302, original `t=0`|Nonsplit; no branch from this conic|
|Reduced `s=1926/2699`|Rational split; independent of M17 by a nonhalving obstruction|
|Reduced `s=-528/3635`|Rational split; exactly `Q=P10-P14-P15` in the original generic basis|

At the positive control the finite footprint stops at prime151 with rank17.
The branch has a compatible footprint parity
`e1+e4+e6+e12+e14+e16`: **that local compatibility is not a dependence
certificate**. Subtracting this unique parity word gives a rational point
whose duplication quartic, after primitive integral normalization and
monic reduction at157, is

\[
 x^4+94x^3+10x^2+81x+27\pmod {157}.
\]

It has no root in `F157`; its leading coefficient is a unit, so a rational
root cannot escape through infinity. Thus the unique parity-compatible
target has no rational half. The generic mod2 injection proves2-saturation
in the rational span, and the halving-or-cycle theorem proves that this
branch is outside `M17 tensor Q`. The independent checker verifies all157
residues; it does not trust the producer's rational quartic factorization.
Prime157 is used here for a fixed rational-fibre duplication polynomial,
**not** under the global smooth-reduction theorem, which does not apply
to the parent at157.

The dependent control instead reaches an exact cycle after two successful
halves, with multiplier1 and the displayed integral relation. The companion
branch is `-P2+P15`, since the branch sum is `P_w`. There is no second new
direction hidden in the two roots. This explicitly distinguishes a split
incidence from a specialization rank gain, with no exceptional point input.

The independent replay reconstructs the negative trace by a manual,
reversed-order `Q(t)` group law, checks the polynomial maps, specializes all
seventeen generic sections, and rebuilds the finite quotients using Sage's
elliptic group implementation rather than the producer's manual group code.
It also verifies that the generic footprint stops at the first sufficient
prefix. All three generic frames, square tests, and admission outcomes are
retained, including the nonsplit302 exposure.

Construction took0.143 seconds internally; the three address actions took
0.043,0.116,0.120 seconds. Independent replay completed in under one second
wall time including startup. Every action had a25-second process cap. No
point search, parameter sweep, full-atlas job, later cascade point or
production change was used. These timings apply only to this one sparse
three-section word and three old addresses.

- [Frozen cold-input protocol](../../artifacts/generated-results/elliptic-curves/det1092_euclidean_seed_admission_v1/protocol.json)
- [Exact conic and polynomial maps](../../artifacts/generated-results/elliptic-curves/det1092_euclidean_seed_admission_v1/construction.json)
- [Independent admission certificates](../../artifacts/generated-results/elliptic-curves/det1092_euclidean_seed_admission_v1/independent-replay-v2.json)

Checker revision2 replaces a lexicographic list comparison by separate
degree-bound checks. All six actual degrees are at their stated bounds;
no arithmetic outcome changed. The revision1 checker, protocol and replay
are retained, rather than replacing their source hashes silently.

```sh
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_euclidean_seed_admission.sage
```

**Remaining obstruction.** This validates the entire generic-word-to-seed
pipeline, but does not choose a successful cover through302. The known
positive-genus302 carrier remains retrospectively calibrated; no later
point or V3 discovery has been used to replace that missing selection rule.
