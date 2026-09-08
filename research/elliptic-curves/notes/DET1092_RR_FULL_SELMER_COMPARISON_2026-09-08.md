# RR residual Selmer comparison: exact local images, unfinished global descent

## Current mathematical endpoint

**Verified application.** The universal RR sextic, inherited rational
basepoint, first-unlock rational Jacobian class, and source-only control
comparison are complete. The [first class](DET1092_FIRST_UNLOCK_JACOBIAN_CLASS_2026-09-08.md)
and [all nine control classes](DET1092_RR_GENERIC_POINT_SPECIFICITY_2026-09-08.md)
are independent of their rank17 Picard images. They generate certified
rank18 subgroups with zero rational2-torsion and18-dimensional rational
Kummer images. Their Sha images are zero.

**Open.** Full `Sel_2(J)` groups, their quotients by the known rational
images, other Sha classes, and a predictive incidence discriminator are
not computed. The elementary Jacobian signature already fails specificity.
The [torsion-module theorem](DET1092_MARKED_KUMMER_TRANSPORT_2026-09-08.md#broader-obstruction-irreducible-2-torsion)
also excludes the ordinary elliptic-quotient and coefficient-induced
Selmer-transfer routes on these members.

**Verified application, local computation.** A Sage-only construction and
independent replay now certify the complete local 2-Kummer image at 553
of the frozen 640 curve–prime pairs. The inherited generic divisors alone
span each of these local images. This is not a full global Selmer group
or a new rank claim. No licensed backend is required for this calculation.

**Verified application, further local arithmetic.** The
[real/dyadic replay](DET1092_RR_REAL_DYADIC_DESCENT_2026-09-08.md)
now certifies all ten complete real images and the exact true/fake
2-adic image dimensions. The inherited `D0=[K_C-2P0]` is locally nonzero
modulo2 on five controls, and locally divisible by2 on the other four
source-only members and the historical302 member. The subsequent
[complete 2-adic replay](DET1092_RR_COMPLETE_DYADIC_IMAGES_2026-09-08.md)
now constructs full true/fake image generators using inherited divisors
and independently verifies them with finite arithmetic modulo8. Global
Selmer completeness and the remaining unresolved places are still open.

## Cohort and information separation

**Verified application.** Keep the existing nine source-only members
`C_t:s^2=c*q(T;r0(t),0)`, comprising the eight frozen V4 fibres and a
generic-section control on302. The historical first-unlock curve is an
additional, explicitly **retrospective calibration arm**, not a replacement
for the source-only302 member. There is no new selector or parameter.

The [frozen preparation protocol](../../artifacts/generated-results/elliptic-curves/det1092_rr_full_selmer_inputs_v1/protocol.json)
binds those ten equations. Each separate worker input contains exactly
`schema`, `P`, and `Q` for an integral model

\[
y^2+Q(x)y=P(x).
\]

The filenames are equation digests. No exceptional point, rank bound,
outcome label or generic-point divisor is sent to a future descent worker.
The preparation ledger retains the arm labels and exact isomorphisms for
subsequent comparison. Equation-only isolation is not a claim that the
historical first-unlock curve was selected blindly.

## Exact normalization and independent replay

**Verified application.** Primitive normalization and extraction of exactly
verified square factors give an integral sextic `F` and positive rational
`k` with

\[
cq(T)=k^2F(T).
\]

One Cremona--Stoll reduction per case gives integers `a,b,c,d` with
`ad-bc=1`, an integral polynomial `H`, and integral `P,Q`, satisfying

\[
2H=Q,\qquad
(cx+d)^6F\!\left(\frac{ax+b}{cx+d}\right)=P(x)+H(x)^2.
\]

The exact map back to the original RR model is

\[
T=\frac{ax+b}{cx+d},\qquad
s=\frac{k(y+H(x))}{(cx+d)^3}.
\]

The [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_rr_full_selmer_inputs_v1/replay.json)
checks the complete polynomial identities, integral coefficients, primitive
sextics, and every equation-only payload. It does not rerun PARI reduction.
The normalization preserves the curve, its Jacobian and its Selmer group;
it is not a search over different covers. The map is written on a finite
patch and extends to the smooth projective curves.

**Verified application; arithmetic limitation.** The reduced coefficient
sizes are3292--4523bits. Trial division through97 leaves polynomial
discriminant cofactors of31884--44206bits. Each exact cofactor and every
trial valuation is retained, and the replay recomputes the polynomial
discriminant by a Sylvester determinant.

These are **not certified field discriminants**, complete bad-prime lists,
or lower bounds for the cost of every possible descent algorithm. No
unfactored cofactor is assumed prime, squarefree, or irrelevant. No full
factorization, maximal-order computation, class group or unit group was run.

## Exact good-prime descent without a licensed backend

**Verified application.** The
[frozen protocol](../../artifacts/generated-results/elliptic-curves/det1092_rr_good_local_images_v1/protocol.json)
uses precisely the previous 64 odd primes and ten RR equations. It uses
the 17 inherited generic section divisors and the inherited basepoint,
not the marked exceptional point or the marked control point. The
historical equation remains a retrospective calibration arm.

Write the inherited divisor representatives as

\[
\beta_j=(-1)^{d_j}g_j(\theta)/(T_0-\theta)^{d_j}
\quad\text{in } A_p^*/A_p^{*2}\mathbf Q_p^*,
\qquad A_p=\mathbf Q_p[T]/q(T).
\]

Normalize the sextic and each divisor polynomial separately by powers of
the current prime. These rational scalars do not alter their fake Kummer
classes. If one divisor is nonunit, omit only that divisor, not the whole
prime. If the basepoint has negative valuation, use the scalar-equivalent
anchor `1-T/T0`. This recovers 152 pairs previously discarded by the
all-generators-must-be-units test. No new prime or point search is involved.

The normalization verifies a squarefree degree-six reduction and an even
valuation of the adjusted scalar twist before claiming good reduction.
Failure of this test is **not** a proof of intrinsically bad reduction.

**Established literature.** At good odd reduction the Kummer image is
unramified; the even-degree fake map has the explicitly described kernel
in [Poonen–Schaefer, §§11–12](https://math.mit.edu/~poonen/papers/descent.pdf).

**New deduction and verified application.** Let the reduced sextic have
`k` irreducible factors of degrees `d1,...,dk`. Good reduction and unique
2-divisibility of the formal group identify

\[
J(\mathbf Q_p)/2J(\mathbf Q_p)
\simeq J(\mathbf F_p)/2J(\mathbf F_p).
\]

The dimension is the Frobenius-fixed dimension of `J[2]`, represented by
even subsets of six branch points modulo complementation:

\[
d_p=\begin{cases}
k-1,&\text{all }d_i\text{ are even},\\
k-2,&\text{some }d_i\text{ is odd}.
\end{cases}
\]

For this genus-two good-reduction situation the fake kernel is zero.
Indeed an odd-degree factor gives injectivity by the cited kernel theorem;
if all factor degrees are even, the unramified quadratic extension splits
each factor into conjugate halves and gives its second injectivity case.
Equivalently, there is a Frobenius-fixed theta characteristic: a fixed
branch point, an invariant `3+3` partition, or an alternating partition
on the even cycles. The replay explicitly verifies a fixed theta and all
16 torsion classes for every encountered factor pattern.

The unit-squareclass target has `k` factor characters. Its norm-square
condition is that the sum of the bits is zero. Rational scalar units
are quotiented by the vector `(d1 mod 2,...,dk mod 2)`. Thus this ambient
space has dimension exactly `d_p`. The recorded generic classes attain
that dimension at every pair marked complete; consequently they span the
**whole true local Kummer image**, not merely a sampled part.

| Frozen case | Complete local images | Reduction-test unresolved | Basepoint unresolved |
| --- | ---: | ---: | ---: |
| scale-0131232 | 53 | 11 | 0 |
| scale-0257585 | 54 | 10 | 0 |
| scale-0487239 | 57 | 6 | 1 |
| scale-0177036 | 55 | 8 | 1 |
| scale-0043332 | 58 | 6 | 0 |
| scale-0590501 | 56 | 8 | 0 |
| scale-0290097 | 53 | 11 | 0 |
| scale-0748009 | 55 | 8 | 1 |
| 302 generic-section control | 54 | 9 | 1 |
| 302 first-unlock calibration | 58 | 6 | 0 |
| Total | 553 | 83 | 4 |

The [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_rr_good_local_images_v1/replay.json)
uses finite-extension Frobenius norms instead of the constructor's
modular polynomial powering, bit elimination instead of its rank routine,
and branch-subset enumeration instead of assuming the dimension formula.
Source hashes bind the previously proved inherited-divisor identities.
All 640 outcomes, including the 87 unresolved pairs, are retained.

**New deduction, with an important limitation.** At each certified place
separately, no local Kummer class lies outside the inherited local image.
This is not simultaneous global spanning: the combination of inherited
generators required may differ between primes. The already certified
global independence of the marked classes therefore remains consistent
with this result. These computations do not identify new Sha classes.

## Backend audit and remaining local implementation

**Verified environment audit.** The available Sage10.9/PARI2.17.3 instance
has no exposed full Selmer method on its genus-two Jacobian class, and no
`magma` executable is found on `PATH`. The repository's completed small
Jacobian Selmer examples use elliptic quotient gluing and are not a general
backend for these simple Jacobians. No computation is currently running.

The earlier public-page check returned HTTP401. The
[access snapshot](../../artifacts/generated-results/elliptic-curves/det1092_rr_full_selmer_inputs_v1/backend-access-audit.json)
is retained as historical evidence, **not the current endpoint status**.
A later check of Sage's documented XML endpoint successfully evaluated
`1+1`. A toy-curve probe failed immediately at `Alarm(3)` with
`Illegal operation`; its Selmer instruction never executed. No RR equation
was submitted. Thus the webpage failure did not establish calculator
unavailability. Per the user's resource constraint, online Magma is not
the substantial-computation route and no license or paid backend is sought.

**Established implementation documentation.** Magma's
[TwoSelmerGroup for hyperelliptic Jacobians](https://magma.maths.usyd.edu.au/magma/handbook/text/1618)
returns the Selmer group and its map to the descent algebra, with options
for raw and even-degree fake data. Its documentation identifies class and
unit arithmetic as a major part of the computation. The curve
`TwoCoverDescent` routine is a different object and must not substitute for
the Jacobian group.

**Implementation plan, not a completed theorem.** Continue with local
Sage/PARI arithmetic and small independently checkable custom components:
resolve the four anchor cases and the remaining genuinely bad-place
images with certified local squareclass arithmetic. The complete real
and 2-adic images, including their true/fake kernel distinction, are now
certified and spanned by inherited subgroups.
A full global
answer also needs a complete supported squareclass space, including any
unramified class-group contribution; a span of known rational classes or
trial-prime support is not a substitute. There is currently no certified
completeness algorithm for that global step in this repository.

No broad class-group calculation is authorized by these scripts. Any
future run must retain its full/conditional status, field and
unit certificates, complete required local support, raw group/map data,
and all failures. Neither heuristic class-group bounds nor a timeout may
produce an unconditional Selmer dimension. Global point-search entrypoints
such as rank-bound routines that also search points are not part of this
descent-only step.

## What a completed Selmer answer would prove

**New deduction from the certified rational subgroup and Kummer sequence.**
For every case, let `V` be the certified18-dimensional rational Kummer
subspace, containing the17-dimensional inherited image.

- If a full unconditional calculation gives `dim Sel_2(J)=18`, then
  `V=Sel_2(J)`, the Jacobian has exact rank18, and `Sha(J)[2]=0`.
- If it gives a larger dimension, the excess is outside the currently
  certified rational image. It must **not** automatically be called Sha:
  further rational points or an exact obstruction are required to decide.
- Cassels--Tate calculations on only the known rational subspace vanish
  and do not answer that excess-class question.

Only after the raw descent results are frozen should they be compared
with the retained inherited images, first witness and V4 outcome snapshot.
A difference of Selmer dimensions alone would be an observed discriminator,
not an equivalence with elliptic rank incidence or a construction of an
extra elliptic point. No such difference is currently asserted.

## Replay and limits

```sh
sage -python research/elliptic-curves/cas/verify_det1092_rr_full_selmer_inputs.sage
sage -python research/elliptic-curves/cas/verify_det1092_rr_good_local_images.sage
```

Ten checkpointed preparations took about half a second each. Whole-panel
independent replay took about half a second. Each preparation and the replay
had a25-second cap. The bounds were ten fixed cases, one model reduction per
case, and trial primes at most97. No point search, full Selmer run, class
group, full integer factorization, or V3/V4 mutation occurred.

The local-image constructor and independent whole-panel replay each
completed within a 25-second cap, with per-prime construction checkpoints
and per-case replay checkpoints. They use no remote calculations or
class groups. No detached job was started or left running.
