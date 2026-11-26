# Uniform eventual Strong Factorial for two monomials

Let \(a,b\) be distinct nonnegative integer exponent columns and
remove rows on which they agree. Throughout, the factorial moments
are the ordinary ones,
\[
 Q_j(t)=\sum_{k=0}^j \binom jk
            \prod_i ((j-k)a_i+kb_i)!\,t^k.
\]
Write \(S=\sum_i(b_i-a_i)^2/\max(a_i,b_i)\).
A row is private when one entry is zero. A profile is fully shared
when there are no private rows.

There is an absolute terminal cutoff above which the consecutive
factorial moments are coprime for every ordinary two-monomial
profile. The proof of this uniform eventual result appears below.
The [effective absolute theorem](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md)
specifies a terminating computation of an absolute cutoff.
Its numerical value is unevaluated. Coprimality at the earlier
terminals requires a separate argument.

**Theorem.** For each fixed \(s_0>0\) there is a terminal
\(N_0(s_0)\) such that every fully shared profile with \(S\ge s_0\)
has \(\gcd(Q_{N-1},Q_N)=1\) for \(N\ge N_0(s_0)\).
This is uniform in the dimension, the exponents, and the difference
of the total degrees. The threshold is unevaluated.

For fully shared profiles the midpoint curvature \(W\) satisfies
\(S\le W\le2S\). Suppose an unbounded sequence of terminals had
common roots. The [unbounded-curvature theorem](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md)
allows us to pass to a subsequence with \(W\) in a compact positive
band. Put \(\delta=\max\rho-1\), where
\(\rho(q)=q(1-q)\sum_i(b_i-a_i)^2/(a_i+(b_i-a_i)q)\).
The [strictly subcritical](MULTIVARIATE_UNIFORM_SUBCRITICAL_BANDS.md)
and [strictly supercritical](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md)
band theorems force \(\delta\to0\).
The [critical regular comparison](MULTIVARIATE_CRITICAL_REGULAR_COMPLEMENT.md)
then forces both centered fields into the quartic neighborhood.
Bounded \(\sqrt n|\delta|\) is covered by the
[moving quartic theorem](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md).
Otherwise pass to a subsequence with \(\sqrt n|\delta|\to\infty\)
and fixed sign of \(\delta\). The
[supercritical thermal proof](MULTIVARIATE_SHARED_SUPERCRITICAL_THERMAL_APPROACH.md)
excludes the positive sign. The
[subcritical thermal comparison](MULTIVARIATE_SHARED_SUBCRITICAL_THERMAL_REGULAR_FIELDS.md)
places both fields at the same limiting scaled fold, which the
[Airy branch comparison](MULTIVARIATE_SHARED_SUBCRITICAL_THERMAL_FOLD.md)
excludes. This exhausts the sequences and proves the theorem.

**Corollary.** There is one absolute eventual terminal cutoff for
fully shared ordinary profiles. The proof also gives one absolute
eventual cutoff on the union of fully shared profiles and profiles
with private variables on both sides. If degree-unbounded
counterexamples exist, they have a subsequence with bounded \(S\)
and with private variables on exactly one side.

Indeed the [vanishing-curvature theorem](MULTIVARIATE_VANISHING_CURVATURE_EVENTUAL.md)
gives one cutoff throughout a fixed neighborhood of \(W=0\).
Since \(S\le W\le2S\) for shared profiles, this and the theorem
above yield one cutoff on the full shared family. The
[both-private theorem](MULTIVARIATE_BOTH_PRIVATE_UNIFORM_EVENTUAL.md)
handles the second family. Together with the unbounded-curvature
theorem this leaves exactly one private side and bounded \(S\).
There the sum of the private exponents is a bounded positive
integer, since every private row contributes its exponent to \(S\).
The [bounded-private-degree annulus](MULTIVARIATE_BOUNDED_PRIVATE_ROOT_ANNULUS.md)
gives a uniform compact range for endpoint-normalized roots
throughout this family.
For one private linear variable, the
[moving cubic model](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md)
gives the exact linked interpolation near the critical endpoint.

## The remaining large-order regime

**Theorem.** There is an absolute \(N_*\) such that
\(\gcd(Q_{N-1},Q_N)=1\) for all \(N\ge N_*\) and all ordinary
two-monomial profiles. The conclusion holds over every
characteristic-zero field.

Suppose an unbounded sequence of terminals has common roots.
First use the preceding support-type corollary. The
[uniform one-private comparison](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md)
excludes private total degree at least two and all strictly
supercritical limits for private degree one. It localizes critical
limits to the positive endpoint.
The [finite subcritical activation](MULTIVARIATE_UNIFORM_SUBCRITICAL_PRIVATE_ACTIVATION.md)
does the same at strictly subcritical limits, including zero
shared curvature. The
[uniform endpoint comparison](MULTIVARIATE_PRIVATE_ENDPOINT_UNIFORM_THERMAL_TAILS.md)
then excludes all strictly subcritical limits and every escaping
temperature or scaled field at critical curvature. Along a
subsequence the only private row is linear, the shared endpoint
curvature \(\sigma\) tends to one, and the cubic temperatures and
both scaled fields stay bounded.

Use the [moving cubic interpolation](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md).
Write \(n=N-1\), \(\epsilon_j=j^{-1/3}\), and
\[
 T_j=\ell^2(1-\sigma)/\epsilon_j,\qquad
 H_j=\ell\xi_j/\epsilon_j^2.
\]
Both \(T_j\) converge along a subsequence to the same real
\(T_*\). At a common root the exact drift is
\(\xi_N-\xi_n=\Delta\log(1+1/n)\in\mathbb R\), so the two scaled
limits have equal imaginary height. Uniform interpolation makes
both limits zeros at \(T_*\). The
[global Riccati comparison](MULTIVARIATE_HALF_LINE_CUBIC_RICCATI_POSITIVITY.md)
forces the limits to be the same simple zero \(H_*\).

For each actual profile the normalized interpolation near this
simple zero has a unique branch \(h(T,\epsilon)\). On a fixed
small field circle around \(H_*\), uniform convergence to the
model gives that one root by Rouché. On a smaller disk the
field derivatives converge uniformly to \(F_H\), which stays
close to its nonzero value at the center. Thus the exponentially
small interpolation residual at either integer order moves its
field by only an exponentially small amount. The same convergence
for the temperature derivative gives \(h_T\to H_T\) along the
adjacent intervals. The bounded first \(\epsilon\) derivative
of the interpolation gives bounded \(h_\epsilon\) there.
These are the local estimates in the
[simple-branch comparison](MULTIVARIATE_PRIVATE_CUBIC_DOUBLE_LIMITS.md).
In physical field coordinates the interpolation branch is
\(\xi(\epsilon)=\epsilon^2 h(T(\epsilon),\epsilon)/\ell\), with
\(T(\epsilon)=\ell^2(1-\sigma)/\epsilon\). Its derivative is
\[
 \xi'(\epsilon)=\frac{\epsilon}{\ell}
       (2h-T h_T+\epsilon h_\epsilon).
\]
Since \(\ell\) is bounded above and away from zero, integration
from \(\epsilon_n\) to \(\epsilon_N\) gives
\[
 \xi_N-\xi_n=-\frac{2H_*-T_*H_T(T_*)+o(1)}{3\ell n^{5/3}} .
\]
For nonzero \(\Delta\) this contradicts the exact drift of
size at least \(\log(1+1/n)\). For \(\Delta=0\) it forces
\(2H_*-T_*H_T(T_*)=0\).

The Riccati comparison excludes zero speed for \(T\le0\), and
the [extended interval certificate](MULTIVARIATE_HALF_LINE_CUBIC_EXTENDED_TEMPERATURES.md)
excludes it for \(0\le T\le2\).
For \(T\ge2\) the [compact field theorem](MULTIVARIATE_HALF_LINE_CUBIC_SPEED_COMPACT_FIELD.md)
and the [positive speed certificate](MULTIVARIATE_HALF_LINE_CUBIC_POSITIVE_SPEED_CERTIFICATE.md)
exclude it throughout the upper half-plane. Conjugation gives
the result at lower fields, and there are no real zeros.
This contradicts the existence of the initial sequence and
proves the uniform eventual assertion. No estimate of \(N_*\)
follows from the sequential comparison used here.

## The fixed-terminal obligation

The absolute eventual theorem leaves only finitely many possible
terminal orders, starting at fourteen by the
[complete terminal-thirteen theorem](MULTIVARIATE_TERMINAL_THIRTEEN_ALL_PROFILES.md).
The [all-degree both-private theorem](MULTIVARIATE_BOTH_PRIVATE_ALL_DEGREES.md)
closes every profile with private variables on both sides. For one
private side, write \(P\) for its total private exponent and \(S_{\rm sh}\)
for the shared absolute curvature. The
[private entropy theorem](MULTIVARIATE_ONE_PRIVATE_ENTROPY_THIRTEEN.md)
closes \(2P+S_{\rm sh}\ge13\). The
[private-degree-six theorem](MULTIVARIATE_ONE_PRIVATE_DEGREE_SIX_ALL_ANGLES.md)
handles the remaining \(P=6\) profiles, so every \(P\ge6\) is
covered with arbitrary shared rows. The [integer-twist theorem for one private unit](MULTIVARIATE_ONE_PRIVATE_UNIT_ALL_INTEGER_TWISTS.md)
closes \(P=1\) at a computable positive weak shared curvature cap,
with arbitrary signed shared degree difference. The
[Gamma resolvent comparison](MULTIVARIATE_ONE_PRIVATE_SHARED_UNIT_ALL_SHAPES.md)
gives every degree for the exact template \((a,a+1),(0,1)\) at every
positive integer \(a\). Other one-private profiles and the fully
shared profiles still require comparisons
in arbitrary dimension. The [fixed-terminal compactification](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md)
describes the full closure of the normalized coefficient pairs and
the relative center. A certified separation of that family gives
a finite algebraic certificate for the terminal. The closure includes
a signed center parameter as well as a diffuse variance. A separation
test that varies the two moment arrays independently admits spurious
collisions already at six, as shown by the
[coefficient-link calculation](MULTIVARIATE_TERMINAL_SIX_COEFFICIENT_LINKS.md).
A proof at each remaining terminal must retain the links between the
two orders or use the exact arithmetic of the factorial coefficients.
The [nonsemialgebraic boundary theorem](MULTIVARIATE_COMPACTIFICATION_NONSEMIALGEBRAIC_BOUNDARY.md)
finds an infinite discrete algebraic section of the exact compactified
outer-loss projection, already in a fixed positive-curvature band.
The [finite-prime interpolation theorem](MULTIVARIATE_REAL_LOCAL_COEFFICIENT_APPROXIMATION.md)
shows that finite local coefficient data are independent of the real
closure at fixed signed degree difference. The approximating profiles
can be chosen with all-degree angular separation. These statements leave
finite outer certificates and exact profile-dependent arithmetic tests
as applicable approaches. The
[balanced row interpolation](MULTIVARIATE_UNIT_GAP_LOCAL_COEFFICIENT_INTERPOLATION.md)
gives arbitrary adjacent coefficient prescriptions at every finite
set of primes with a bounded number of rows of one fixed gap.
It preserves any prescribed diffuse limit while forcing matched
exterior branches whose permitted valuation spread grows at any
specified rate. Even a single rational quadratic can match all
finite local branches and a prescribed real conjugate pair near
a diffuse contact. The constructed ordinary profiles can be
regular in every degree, so these finite tests alone cannot
give a uniform bound on common-factor denominators.
The [bounded-gap decision theorem](MULTIVARIATE_BOUNDED_GAP_FINITE_DECISION.md)
gives two-sided exact decisions for angular separation and several
physical rescaling questions at fixed total gap. Its affordable twists
measure exactly the signed gap that can escape while consuming the
remaining gap budget. The [one-unbounded-type reduction](MULTIVARIATE_ONE_UNBOUNDED_TYPE_FINITE_DECISION.md)
instead decides the original physical property for profiles with one
row type above a fixed height, even at arbitrary multiplicities.
A hard-row cutoff gives a finite original decision for every profile
with at most two unequal rows.
The [joint curvature measure](MULTIVARIATE_JOINT_FACTORIAL_CURVATURE.md)
uses an exact factorial density for each primitive row curve. It
couples all orders and determines the relative center. Relaxing only
these densities to arbitrary positive measures on the same curves
already gives a collision at orders two and three.
The [uniform curvature closure](MULTIVARIATE_JOINT_CURVATURE_PROFILE_CLOSURE.md)
identifies the actual diffuse part as a diagonal density \(r\,dr\).
Its error is uniform in the two curvature indices, and the relative
centers have one additional signed diffuse parameter.
The [exact row replacement theorem](MULTIVARIATE_EXACT_FACTORIAL_ROW_REPLACEMENTS.md)
characterizes raw moment equality along common translations by
integer projection counts. Matching through the absolute cutoff
transfers all-degree separation to replacement profiles. This
includes replacing a mirrored shared block by one with no reflected
pairs and the same unordered row counts. The
[single-order rigidity theorem](MULTIVARIATE_SINGLE_ORDER_ROW_RIGIDITY.md)
recovers the ordered rows from one normalized array once its degree
exceeds the stated gap-count bound.
The [endpoint comparison](MULTIVARIATE_ORDINARY_ADJACENT_COEFFICIENT_ATTENUATION.md)
quantifies the decline of each coefficient multiplier at fixed
distance from an endpoint and retains four additional cross-order
links at terminal six.
At terminal six, [one-private-side separation](MULTIVARIATE_TERMINAL_SIX_ONE_PRIVATE_SIDE.md)
and the [two-private linear certificate](MULTIVARIATE_TERMINAL_SIX_TWO_PRIVATE_LINEAR.md),
together with the [other private-degree criteria](MULTIVARIATE_TERMINAL_SIX_PRIVATE_ANGLES.md),
cover every profile with any private variable. Thus the remaining
terminal-six family is fully shared. The
[combined curvature certificate](MULTIVARIATE_TERMINAL_SIX_HYBRID_CURVATURE.md)
covers every fully shared profile with \(S\ge4\).
The [effective relative-gap bands](MULTIVARIATE_TERMINAL_SIX_DIFFUSE_BANDS.md)
and [nine linked outer-mean bands](MULTIVARIATE_TERMINAL_SIX_JOINT_RAY_BANDS.md)
cover additional profiles in the intermediate range. The
[extended near-flat certificate](MULTIVARIATE_TERMINAL_SIX_NEAR_FLAT_EXTENSION.md)
and the [intermediate Cayley certificate](MULTIVARIATE_TERMINAL_SIX_INTERMEDIATE_CAYLEY.md)
together cover the order-five outer mean \(u\ge18/25\).
This includes \(S\le4/9\). Indeed, if
\(\eta=\max_i |b_i-a_i|/\max(a_i,b_i)\) then \(\eta\le S\) by
integrality. The coefficient bound in the diffuse-band note gives
\[
 u\ge\exp\!\left(-\frac{2S}{5(1-\eta)}\right)\ge e^{-8/25}
 >\frac{18}{25};
\]
the last inequality follows from
\(e^x\le1+x+x^2/[2(1-x/3)]\) at \(x=8/25\), by the exponential
series. When \(S<4\), the
[retained-row certificate](MULTIVARIATE_TERMINAL_SIX_RETAINED_ROWS.md)
covers every profile with a Gauss reciprocal above \(1/8\) except
when the retained list consists of a single row of unordered
exponents \((1,2)\). The
[marked-row tail certificate](MULTIVARIATE_TERMINAL_SIX_MARKED_ROW_TAIL.md)
covers that list for the whole range \(S<4\), completing the
case of at least one reciprocal above \(1/8\). For reciprocal cap
\(1/8\), the first joint ray band covers \(0<u\le7/25\), the
[paired Bernstein certificate](MULTIVARIATE_TERMINAL_SIX_PAIRED_LIGHT_ROWS.md)
covers \(7/25\le u\le3/5\), and the
[diffuse Cayley certificate](MULTIVARIATE_TERMINAL_SIX_DIFFUSE_CAYLEY.md)
covers \(3/5\le u\le18/25\). Together with the intermediate and
near-flat certificates, this completes terminal six.
At terminal seven, the [private-degree proof](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_ROWS.md)
and the [one-private-linear proof](MULTIVARIATE_TERMINAL_SEVEN_PRIVATE_LINEAR.md)
cover every profile with a private row. For fully shared profiles,
the [joint high-curvature proof](MULTIVARIATE_TERMINAL_SEVEN_SHARED_HIGH_CURVATURE.md)
covers \(S\ge17/4\). At lower curvature the
[retained-row replay](MULTIVARIATE_TERMINAL_SEVEN_RETAINED_ROWS.md)
and the two marked extensions for
[(1,2) rows](MULTIVARIATE_TERMINAL_SEVEN_MARKED_LIGHT_TAIL.md)
and [(2,4) rows](MULTIVARIATE_TERMINAL_SEVEN_MARKED_24_TAIL.md)
cover every profile with some sixth-order reciprocal above \(1/8\).
At or below that cap, the [small-outer-mean proof](MULTIVARIATE_TERMINAL_SEVEN_SMALL_OUTER.md)
covers \(u\le1/10\), and the
[extended cap-one-eighth replay](MULTIVARIATE_TERMINAL_SEVEN_LIGHT_CAP_EXTENSION.md)
covers \(1/10\le u\le12/25\). The
[paired Bernstein bands](MULTIVARIATE_TERMINAL_SEVEN_PAIRED_LIGHT.md)
cover the stated smaller caps from \(12/25\) to \(22/25\).
The [upper Cayley bridge](MULTIVARIATE_TERMINAL_SEVEN_UPPER_BRIDGE.md)
covers \(22/25\le u\le19/20\), and the
[near-flat Cayley proof](MULTIVARIATE_TERMINAL_SEVEN_NEAR_FLAT.md)
covers \(u\ge19/20\), both with unrestricted reciprocals.
Between \(12/25\) and \(22/25\), the
[cross-order skew transfer](MULTIVARIATE_TERMINAL_SEVEN_CROSS_SKEW.md),
[two retained bands](MULTIVARIATE_TERMINAL_SEVEN_RETAINED_BANDS.md)
and [late Cayley bridge](MULTIVARIATE_TERMINAL_SEVEN_LATE_BRIDGE.md)
cover all remaining shared profiles, completing terminal seven.
At terminal eight, the [private-degree proof](MULTIVARIATE_TERMINAL_EIGHT_PRIVATE_ROWS.md)
and the [private linear proof](MULTIVARIATE_TERMINAL_EIGHT_PRIVATE_LINEAR.md)
cover every profile with a private row. The
[shared high-curvature proof](MULTIVARIATE_TERMINAL_EIGHT_SHARED_HIGH_CURVATURE.md)
covers \(S\ge19/4\), including the
[fixed-ray range](MULTIVARIATE_TERMINAL_EIGHT_HIGH_CURVATURE.md)
\(S\ge15/2\). At lower curvature, the
[retained rows](MULTIVARIATE_TERMINAL_EIGHT_RETAINED_ROWS.md)
and the [three marked singletons](MULTIVARIATE_TERMINAL_EIGHT_MARKED_SINGLETONS.md)
cover all fully shared profiles with a seventh-order reciprocal
above \(1/8\). The [small outer-mean](MULTIVARIATE_TERMINAL_EIGHT_SMALL_OUTER.md)
and [paired light-row](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_LIGHT.md)
proofs cover reciprocal cap \(1/8\) through \(u=12/25\), and the
[paired bridge](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_BRIDGE.md)
and its [continuation](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_BRIDGE_72.md)
extend that light-row result to \(u=18/25\) under their stated
reciprocal caps. The
[middle retained bands](MULTIVARIATE_TERMINAL_EIGHT_RETAINED_BANDS.md)
cover all profiles above those caps and at or below \(1/8\)
throughout the same bands. Thus \(u\le18/25\) is covered for
every shared profile. For \(u\ge22/25\), the
[upper bridge](MULTIVARIATE_TERMINAL_EIGHT_UPPER_BRIDGE.md)
and [near-flat theorem](MULTIVARIATE_TERMINAL_EIGHT_NEAR_FLAT.md)
cover unrestricted reciprocals. The
[late shared bridge](MULTIVARIATE_TERMINAL_EIGHT_LATE_BRIDGE.md)
covers \(18/25\le u\le22/25\), also for unrestricted reciprocals.
Together these results complete terminal eight.
The [terminal-nine theorem](MULTIVARIATE_TERMINAL_NINE_ALL_PROFILES.md)
covers both private cases and every fully shared profile. The proof
combines linked matched rays, exact retained rows, and the near-flat
Cayley comparison. Thus all terminals through nine are complete.
The [linked Hermite proof](MULTIVARIATE_FIXED_TERMINAL_HERMITE_SEPARATION.md)
independently gives a computable small-curvature separation radius
at every fixed terminal.
The [terminal-ten theorem](MULTIVARIATE_TERMINAL_TEN_ALL_PROFILES.md)
covers all private profiles and all fully shared profiles. Its upper
middle band uses seventeen exact retained lists and the same cubic
matched model as the lighter bands. The consecutive pairs through
terminal ten are therefore complete. The
[terminal-eleven theorem](MULTIVARIATE_TERMINAL_ELEVEN_ALL_PROFILES.md)
combines the exact private and retained cases with the shared
moving-ray and Cayley bands. All terminals through eleven are
complete. The [terminal-twelve theorem](MULTIVARIATE_TERMINAL_TWELVE_ALL_PROFILES.md)
combines the private and retained cases with the exact light bands,
the weighted Cayley models and the conditional signed-moment bounds.
The [terminal-thirteen theorem](MULTIVARIATE_TERMINAL_THIRTEEN_ALL_PROFILES.md)
combines the private and retained cases with light rays, grouped root
disks and upper Cayley bands. The covered outer-mean intervals and
reciprocal-cap complements exhaust every ordinary profile. Thus the
consecutive pairs through terminal thirteen are complete. The
[terminal-fourteen fixed rays](MULTIVARIATE_TERMINAL_FOURTEEN_HIGH_CURVATURE.md)
cover \(S\ge8\), every total private degree at least four, and the
listed lower private-degree cases above their shared-curvature
thresholds. The [terminal-fourteen retained profiles](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md)
cover every private profile and all fully shared profiles with a
reciprocal above \(1/8\). The [small-outer-mean rays](MULTIVARIATE_TERMINAL_FOURTEEN_SMALL_OUTER.md)
cover all fully shared profiles with \(0<u\le.02\).
The [left flank](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_MID_LEFT.md),
[middle](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_MID.md), and
[right flank](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_MID_FLANKS.md)
give full coverage on \([.25,.50]\). The
[lower-middle disks](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_LOW_MID.md)
and [retained low-middle complements](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_LOW_MID.md)
cover \([.50,.64]\). The [middle disks](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_DISKS.md)
and [middle complements](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_MIDDLE.md),
followed by the [upper disks](MULTIVARIATE_TERMINAL_FOURTEEN_HIGH_U.md),
[retained bridge](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_BRIDGE.md),
and [upper complements](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_UPPER.md),
close the remaining bands of \([.64,1)\).
Thus fully shared profiles are covered on \((0,.02]\) and \([.25,1)\).
Terminal fourteen is the first remaining fixed terminal.
The [linked-skew supplement](MULTIVARIATE_TERMINAL_THIRTEEN_LINKED_CAYLEY.md)
demonstrates a common separating ray that moves with the matched
signed second moment on a narrow upper band.


For a degree-uniform extension of these comparisons, the
[matched tanh transport](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md)
bounds the cross-order moment error by the same excess outer loss
\(W-RU\). The [Cayley operator test](MULTIVARIATE_CAYLEY_OPERATOR_RAY_TEST.md)
carries every grouped moment through the exact change of variable.
Its ray inequality requires a bound at the reference angular gap scale.
The [Gaussian reference comparison](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md)
gives full interlacing at equal heat parameters and the inner inequalities
at the linked parameters. The outer reference inequalities and the
transformed factorial remainder require further estimates across the
remaining range. The [full-kernel transport](MULTIVARIATE_MATCHED_TANH_KERNEL_TRANSPORT.md)
extends the excess-loss estimate to smooth test functions and gives
explicit bounds for the whole even and odd logarithmic corrections at
common grid locations. The [higher signed-moment centering](MULTIVARIATE_CENTERED_HIGHER_SIGNED_TANH.md)
reduces the independent skew remainder by taking out a larger exact
multiple of the matched second moment. The [sharp capped moment bounds](MULTIVARIATE_SHARP_SKEW_CONDITIONED_TANH_MOMENTS.md)
give an attained simultaneous envelope for the full kernels at fixed skew.
The [midpoint-grid comparison](MULTIVARIATE_MATCHED_FRACTIONAL_LOG_TRANSPORT.md)
carries the full-kernel estimate through to interpolants of the normalized
log coefficients. The [holomorphic matched log transport](MULTIVARIATE_HOLOMORPHIC_MATCHED_LOG_TRANSPORT.md)
bounds all complex derivatives of the exact normalized Gamma chord
difference by the excess outer loss and gives the two same-degree
multipliers for the upper moment. At circle-rooted references, an
integrated root-product estimate controls the first angular variation
of the matched chord by the reference root spacing. The [alternating binomial estimate](MULTIVARIATE_ANALYTIC_FRACTIONAL_BINOMIAL_TRANSPORT.md)
carries this analytic transport through exponentiation and the
alternating sum against the upper-grid lift. The [angular sensitivity bound](MULTIVARIATE_ANGULAR_PHASE_CONDITIONING.md)
shows exponential phase sensitivity to general bounded log-coefficient
directions at angles separated from the negative ray. A [Gaussian comparison path](MULTIVARIATE_FACTORIAL_GAUSSIAN_ANNULUS_HOMOTOPY.md)
keeps a narrow normalized root annulus. For the Gaussian references,
the explicit all-degree comparisons apply at all physical time ratios
for [\(s\ge13\)](MULTIVARIATE_DIFFUSE_EXPLICIT_S13_INTERLACING.md) and
[\(0<s\le2\cdot10^{-4}N^{-3/2}\)](MULTIVARIATE_GAUSSIAN_SMALL_COUPLING_INTERLACING.md). A [fixed-order error obstruction](MULTIVARIATE_FIXED_ORDER_COEFFICIENT_ERROR_OBSTRUCTION.md)
shows that independent coefficient remainders can swamp the gap at high
orders even with sharp endpoint attenuation and exact moment parameters.
The [polynomial diffuse lift](MULTIVARIATE_DIFFUSE_POLYNOMIAL_ANGLE_LIFT.md)
preserves cancellation by applying the analytic coefficient multiplier
to each nonvanishing Gaussian half separately. For the ordinary
moments, this gives angular separation at every terminal \(N\) with
\(S\ge13\) and \(S\eta N\le10^{-6}\).
The [degree-independent diffuse lift](MULTIVARIATE_DIFFUSE_ABSOLUTE_ANGLE_LIFT.md)
centers the analytic multiplier at the shifted saddle, applies its
linear part exactly and couples the radial phase comparisons. This gives
the additional ordinary range \(S\ge13,\ S\eta\le10^{-9}\) at every
terminal. The [signed midpoint curvature test](MULTIVARIATE_MIDPOINT_CURVATURE_SYMBOL.md)
improves the small-gap condition when the signed cubic row sum cancels.
The [matched-tanh phase lift](MULTIVARIATE_MATCHED_TANH_ORDINARY_ANGLE_LIFT.md)
gives an actual all-order angular criterion in the matched moments,
with common Gaussian strength \(2\sum_i r_i\ge13\), signed second
moment \(\sum_i\epsilon_i r_i^2\), and unsigned second and third moments.
The [matched tanh test at the exact midpoint clocks](MULTIVARIATE_MIDPOINT_INFLATION_POLYNOMIAL_ANGLE_TRANSFER.md#direct-matched-tanh-criteria)
covers every positive curvature subject to a signed skew and unsigned
third-moment bound. Inflation of the Gauss midpoint rates with a
compensating Gaussian term preserves the midpoint curvature and
controls the full complex motion of each simple root. The explicit
tolerance for a pair is proportional to \(1/[n(n+1)]\). Replacing
\(n(n+1)\) by the square of the effective absolute cutoff gives a
single all-degree test using the tanh moments from any one order.
The [quadratic skew refinement](MULTIVARIATE_QUADRATIC_SKEW_TANH_ANGLE_TRANSFER.md)
projects the moving root chart onto the unit circle. The cubic
signed curvature contributes first to radial motion, so its
angular effect is quadratic in the matched skew.
The homogeneous bounds leave the common part of the adjacent angular
responses uncontrolled at the comparison scale for general rows.
The [curvature-eight centered lift](MULTIVARIATE_DIFFUSE_S8_CENTERED_ANGLE_LIFT.md)
covers \(8\le S\le13\) for \(n\ge600\) with \(S\eta\le5\cdot10^{-12}\).
The [smaller-disk centered lift](MULTIVARIATE_DIFFUSE_S6_CENTERED_ANGLE_LIFT.md)
covers \(6\le S\le8\) for \(n\ge10000\) with \(S\eta\le4\cdot10^{-15}\).
The [abstract paired half lift](MULTIVARIATE_DIFFUSE_ABSTRACT_PAIRED_HALF_LIFT.md)
uses the [Gaussian comparison at \(16/3\)](MULTIVARIATE_DIFFUSE_EXPLICIT_S16_OVER_3_INTERLACING.md)
to reach \(16/3\le S\le6\) for \(n\ge10^7\) and
\(S\eta\le5\cdot10^{-24}\).
The curvature-eight and curvature-six signed matched-moment criteria allow larger individual gaps
when the row orientations cancel. The underlying Gaussian
comparisons use an adjacent saddle estimate before transferring
the exact factorial coefficients.
The [matched curvature hierarchy](MULTIVARIATE_MATCHED_TANH_ANGULAR_HIERARCHY.md)
bounds the higher derivatives of the normalized symbol and of its adjacent
increment relative to their exact midpoint curvatures. It gives a signed
matched-moment criterion on \(16/3\le2A\le6\) and a cubic-order constraint on
the weak reflected factorial samples. A common positive-Laplace symbol can
produce an exact angular collision at orders three and four. Its samples
violate the weak factorial constraint.
The [weighted matched-row collision](MULTIVARIATE_WEIGHTED_TANH_ANGLE_OBSTRUCTION.md)
retains the exact adjacent factorial kernels, their homogeneous transport
bounds, and the integer-profile size inequality in midpoint form. A
certified resultant sign change gives one common unit-circle pair at
orders eight and nine for two positively weighted mirrored rows. A
discrete packed tanh moment inequality excludes those weights.
The [row-cancelled angular lift](MULTIVARIATE_MATCHED_TANH_ROW_CANCELLED_ANGLE_LIFT.md)
pairs opposite shifts within each Gauss row. It replaces the unsigned
second-moment drift in two sufficient angular criteria by the smaller
of that moment and the sum of squared row rates. Combining the lower
strip with quadratic skew control gives an explicit all-degree criterion
near the curvature interval \([16/3,6]\) from a single matched pair.
The [row concentration theorem](MULTIVARIATE_MATCHED_TANH_ROW_CONCENTRATION.md)
gives a convex-order comparison after first-moment normalization when
the rows share one rate. Global normalization in the mixed-rate case can
move the second moment in either direction even in ordinary mirrored
profiles. The [quadratic packing collision](MULTIVARIATE_MATCHED_TANH_QUADRATIC_PACKING_OBSTRUCTION.md)
satisfies the homogeneous transport and sharp quadratic packing at
both orders. The weighted example violates cubic packing and the
whole-row multiplicity constraints.
The [integer multiplicity criterion](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md)
recovers whole Gauss blocks from the joint tanh factors by divisor
inversion. It also recovers exact row counts from the two marginal
inverse moments and detects the fractional weights in the collision.
The [inverse-moment obstruction](MULTIVARIATE_MATCHED_TANH_FINITE_INVERSE_MOMENT_OBSTRUCTION.md)
preserves the weighted collision while matching an ordinary profile
on every aggregate inverse-tanh polynomial of any prescribed finite
degree. Lagrange extrapolation makes the added shape corrections
arbitrarily small, while the local row count still detects the
fractional weight.
The [first-sector comparison](MULTIVARIATE_COMPARABLE_ROWS_FIRST_SECTOR_COMMON_ROOTS.md)
excludes a common physical root with \(|\arg z|\le2\pi/(n+1)\)
at any orders \(n,n+1\) for at most two comparable unequal rows,
and for arbitrarily many homothetic rows. It uses a conditional
chord comparison and a variance monotonicity inequality for a
Gamma product. The [arbitrary-dimensional comparison](MULTIVARIATE_COMPARABLE_FIRST_SECTOR_HCM_RADIAL_COMPARISON.md)
uses a beta–Gamma coupling and hyperbolic complete monotonicity to
cover any comparable profile. In the increasing orientation it also
rules out all first-sector contacts with upper radius at least the
lower radius. An [effective stability theorem](MULTIVARIATE_COMPARABLE_EFFECTIVE_FIRST_SECTOR_MIXED_TAIL.md)
preserves this sector conclusion under sufficiently small shared
tails with aligned total gap, including some truly mixed profiles.
The [finite positive-moment reconstruction](MULTIVARIATE_MATCHED_TANH_FINITE_POSITIVE_CORE.md)
recovers all ordinary rows above any fixed relative gap under a curvature
cap, by polynomial filtering and integer rounding. A
[finite angular certificate](MULTIVARIATE_MATCHED_TANH_FINITE_ANGULAR_CERTIFICATE.md)
uses positive mixed moments from the single reference pair \((2,3)\)
to approximate the exact log coefficients in every degree below the
absolute angular cutoff. Its error sums with the first moment, so
the certificate also controls the contribution of arbitrarily small
factors. It succeeds at sufficient precision whenever the profile
has disjoint consecutive root arguments.
The [opposed-slot matched-tanh ray theorem](MULTIVARIATE_MATCHED_TANH_OPPOSED_AFFINE_ANGLES.md)
uses two lower Beta slots and a positively rescaled upper binomial
endpoint. Its quotient eliminant keeps only opposite half-planes,
which excludes the upper binomial quotients at the endpoint.
In particular, any fixed ordinary core has all-degree angular
separation after adding two opposite moving rows of a common
smaller base, for all sufficiently large integer bases, even
with unequal gaps. The first-slot cohort may have unequal
orientation counts. The [opposed affine trace bound](MULTIVARIATE_OPPOSED_AFFINE_TRACE_MENU.md)
gives a polynomial bound on the number of exceptional bases even
when unrelated traces coincide identically. With two mirrored
unit-gap rows, any list of distinct bases longer than this bound
contains an all-degree witness for any core.
The [mirrored endpoint phase comparison](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_FOURTEEN.md)
gives strict angular interlacing in every degree when the mirrored
curvature plus any common Gaussian strength is at least fourteen.
The [curvature-thirteen refinement](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_THIRTEEN.md)
lowers the threshold by retaining the dependence between curvature
and the endpoint mass in the order derivative and using the squared
cost of the largest mirrored relative gap. The
[curvature-twelve comparison](MULTIVARIATE_MATCHED_TANH_MIRRORED_CURVATURE_TWELVE.md)
uses a squared-curvature endpoint bound, the paired binomial
entropy derivative and a sharper differentiated connected
expansion to lower the threshold again.
The [ratio-balanced extension](MULTIVARIATE_RATIO_BALANCED_ANGULAR_SEPARATION.md)
retains \(S\ge12\) for fully shared profiles with equal total
scale of both orientations at each reduced endpoint ratio.
It accounts for the radial displacement of the actual roots
when the rows themselves are not mirrored.
The [private-unit geometric test](MULTIVARIATE_PRIVATE_UNIT_MATCHED_TANH_GEOMETRIC_AND_INTEGER_TWISTS.md)
uses the unsigned shared first moment to give all-degree separation
for two private unit rows and a weak shared profile of any orientation.
With one comparable shared unit-gap row, its integer-twist result also
gives coprimality at every positive rescaling \(((n+1)/n)^h\) with
integer \(h\), allowing small shared tails with any signed total gap.
The [critical core packing theorem](MULTIVARIATE_MATCHED_TANH_CRITICAL_CORE_PACKING.md)
uses these linked data to obtain every-degree angular separation
near \(W=4\) with small signed skew and at least one relative gap
\(d/(2\ell+d)>1/3\). It also allows an asymmetric residual profile.
The [Gaussian gap obstruction](MULTIVARIATE_GAUSSIAN_MIDPOINT_GAP_EXTREMALITY_OBSTRUCTION.md)
shows that the exact ordinary corrections can reduce the Gaussian
midpoint gap even with arbitrarily small individual relative gaps.
The [two-midpoint-clock test](MULTIVARIATE_MATCHED_MIDPOINT_CLOCK_ANGLE_TEST.md)
uses the exact physical Gaussian clocks at the two orders. Its signed
shape error is the exact midpoint skew plus a row-curvature-weighted
sum of squared relative gaps. The all-coupling physical gap
makes this a single all-degree signed tolerance at unrestricted
curvature, including asymmetric profiles with limiting skew
cancellation.
The [Gaussian physical-clock comparison](MULTIVARIATE_GAUSSIAN_PHYSICAL_CLOCK_COMPARISON.md)
proves a late-clock outer-gap theorem and reduces the intermediate outer sign
to a first-contact test in matched squared-height moments. The
[weak sharp-clock corollary](MULTIVARIATE_UNIFORM_WEAK_COMMON_SYMBOL_INTERLACING.md#sharp-clock-gaussian-corollary)
proves Gaussian interlacing at the smallest admissible clock ratio for
every degree on a fixed computable weak interval. The
[matched squared-height coefficients](MULTIVARIATE_GAUSSIAN_MATCHED_HEIGHT_COEFFICIENTS.md)
give an exact positive integral for the coefficient ratios at the sharp
clock and prove strict Gaussian interlacing at all times for clock
ratios sufficiently close to one. At the sharp ratio, the
[bulk first-contact exclusion](MULTIVARIATE_GAUSSIAN_SHARP_CLOCK_BULK_CONTACT.md)
shows that a large-degree first contact at bounded positive normalized time
must approach the limiting edge, with rank sublinear in the degree.
The same proof gives an exact zero-free sector up to the limiting
subcritical edge at every degree.
The [exact Gaussian saddle phase](MULTIVARIATE_GAUSSIAN_EXACT_SADDLE_PHASE.md)
labels every finite-degree root in a phase interval of length
\(\pi/2\), gives microscopic edge bounds uniform through the critical
crossover, and gives strict sharp-clock interlacing on the part with
phase at least \(\pi/4\) when the normalized time is at most six.
The [radial minimum argument](MULTIVARIATE_GAUSSIAN_HIGH_PHASE_ALL_COUPLINGS.md)
extends the high-phase conclusion to every positive normalized time.
The [subcritical sharp-clock theorem](MULTIVARIATE_GAUSSIAN_SUBCRITICAL_SHARP_CLOCK_INTERLACING.md)
proves the strict adjacent phase sign for all ranks and all degrees at
normalized coupling at most two. The
[quantitative physical gaps](MULTIVARIATE_GAUSSIAN_QUANTITATIVE_PHYSICAL_GAPS.md)
give explicit angular margins at both endpoint clocks and throughout
the interval between them. The
[supercritical collar](MULTIVARIATE_GAUSSIAN_SUPERCRITICAL_COLLAR.md)
gives an explicit sharp-clock margin beyond criticality at each degree
and a uniform collar in terms of the eventual Gaussian threshold.
The [paired supercritical half integral](MULTIVARIATE_GAUSSIAN_PAIRED_HALF_SUPERCRITICAL.md)
gives all-order Gaussian separation at normalized time \(t\ge26/5\).
Its exact rational replay controls the difference of two successive
contour errors and reuses the earlier low-order Gaussian exclusions.
The [all-coupling phase reserve](MULTIVARIATE_GAUSSIAN_SUPERCRITICAL_PHASE_RESERVE.md)
proves strict consecutive Gaussian interlacing at every degree and
every positive physical clock. Its continuous three-parameter
certificate is needed only for \(2<t<4\). The dyadic width estimate
then gives a uniform positive adjacent phase sign in every degree.
The [polynomial physical gap](MULTIVARIATE_GAUSSIAN_POLYNOMIAL_PHYSICAL_GAPS.md)
gives an explicit \(t/[36n(n+1)(1+t)]\) angular margin across the
full physical clock interval.
The [absolute ordinary diffuse transfer](MULTIVARIATE_DIFFUSE_ALL_ORDER_RESULTANT_TRANSFER.md)
uses that phase margin to give one computable positive order-independent
diffuseness tolerance on the entire curvature range. It also gives
one uniform tolerance for independent blends of the two exact symbols.
No exceptional Gaussian curvature set remains.
The [arbitrary-core transfer](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md)
keeps any fixed ordinary core exact. A resultant on root ratios
covers all positive algebraic Gaussian strengths, even with asymmetric
off-circle roots. The ordinary extension theorem works in all
degrees for every sufficiently diffuse tail at each prescribed
positive rational added curvature, and in a full one-sided range
at zero added curvature with a quantitative tolerance. The
[uniform tight-core theorem](MULTIVARIATE_UNIFORM_TIGHT_CORE_DIFFUSE_STABILITY.md)
makes this tolerance uniform over cores with bounded row count
and computably tight curvature tails. A bound on the total absolute
core gap is sufficient, with all baselines and orientations allowed.
The [finite diffuse strength menu](MULTIVARIATE_UNIVERSAL_DIFFUSE_STRENGTH_MENU.md)
removes all core restrictions by allowing the strength to be selected
from a finite universal list in any positive interval. Its rational
strengths and any signed integer gap can be realized exactly by
ordinary rows with arbitrarily small relative gaps.
The [physical diffuse theorem](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md)
retains signed center drift as a second variable. For any fixed core,
positive diffuse limits with a possible physical common root form a
finite explicit set of points in the curvature and signed-drift plane.
The two coefficients in a contact are logarithms of positive algebraic
numbers. In particular every algebraic affine line other than one of
rational slope through the origin avoids the exception set. A quantitative
perturbation bound localizes any ordinary failure to a neighborhood
of one of these points, uniformly in the tail row count.

For \(0<\rho\le1\), retain as the core all rows with \(d/M\ge\rho\),
where \(M\) is the larger exponent and \(d\) the absolute gap.
Each retained row contributes at least \(\rho^2 M\) and at least
\(\rho\) to \(S\). Thus \(S<12\) gives both \(M<12/\rho^2\) and
fewer than \(12/\rho\) retained rows. At each fixed threshold
there are finitely many cores, while the diffuse part can have
unboundedly many rows and unbounded exponents.

The [signed curvature triangle](MULTIVARIATE_STATIONARY_BIAS_CURVATURE_TRIANGLE.md)
is the sharp limiting domain of the next quartic and squared
skew responses at a fixed strength and bias. It allows all
ordinary tails and every prescribed integer signed gap. The
[stationary dichotomy](MULTIVARIATE_STATIONARY_TAIL_DICHOTOMY.md)
uses its three vertices at simple even stationary interior
contacts with the stated skew vanishing. A positive minimum excludes
all sufficiently diffuse ordinary tails near that branch pair.
A negative minimum gives nearby simple transverse Gaussian contacts
with finite ordinary cores. If the constant second-order phase
response is nonzero, one strict test always applies. The exact
integer exclusion at transverse contacts remains open. The
[zero-clock refinement](MULTIVARIATE_STATIONARY_ZERO_CLOCK_DICHOTOMY.md)
gives a compact exact moment test when the constant clock response
vanishes, using three positive semidefinite matrices of size two.
The [lower-clock two-scale theorem](MULTIVARIATE_STATIONARY_LOWER_CLOCK_TWO_SCALE.md)
handles a same-sign lower-order clock response by testing the
two possible rates of strength displacement. The critical rate
charges the additional inequality \(Z\ge |V|^{2/3}\). Strict
positive margins give uniform ordinary exclusion. At interior
bias, negative margins give finite transverse Gaussian references.
The [transverse response polytope](MULTIVARIATE_TRANSVERSE_TAIL_RESPONSE_POLYTOPE.md)
gives a sharp four-vertex interval for an exact physical bias
residual at a simple positive-strength contact. For distinct
simultaneous simple root pairs, projecting their joint response
equations can exclude all unrestricted ordinary tails approaching
that configuration. If all reference roots are simple and all
two-pair tests pass, every surviving common factor is a single
simple irreducible quadratic, uniformly on a compact positive
strength interval for that fixed core. Independent clock velocities
also give uniform exclusion for the indicated moment-balanced
families.

The [signed core replacement theorem](MULTIVARIATE_SIGNED_CORE_REPLACEMENT_CONTACTS.md)
extends the response equations to varying cores by accounting for
both added rows and omitted rows from a countable limit. Under
the stated rank and strict separation tests, the two second-order
sizes must be comparable and can be forced to agree asymptotically.
At any infinite retained core exact integer replacements can make
these responses cancel to second order. There are retained cores
for which equal reciprocal moments give coefficient approximation
through every finite order.

The [common-height tail theorem](MULTIVARIATE_COMMON_HEIGHT_BOUNDED_GAP_TAIL_CLOSURE.md)
excludes exact ordinary common roots in every degree for any fixed
core and a tail of bounded integer gaps at a sufficiently large
shared upper exponent. The height bound holds uniformly over all
multiplicities and orientations, including at every positive
Gaussian exception. The [private-prime density theorem](MULTIVARIATE_SPREAD_OUT_PRIVATE_PRIME_DENSITY.md)
treats independently variable ordinary row heights in consecutive
windows. Its exceptional proportion is exponentially small
under the stated window and retained-prime conditions. This
includes windows approaching any prescribed positive diffuse
strength uniformly across the tuples. The
[dense cluster theorem](MULTIVARIATE_DENSE_CLUSTERED_UNIT_GAP_CLOSURE.md)
gives deterministic closure for sufficiently narrow, densely occupied
bands of unit-gap rows with arbitrary orientation multiplicities.
Each adjacent pair has an isolated large prime among its marked
factors. A fixed retained core or the stated allowance for varying
retained prime supports is permitted, and balanced clusters reach
any prescribed positive diffuse strength. The
[affine-height theorem](MULTIVARIATE_AFFINE_HEIGHT_BOUNDED_GAP_TAIL_CLOSURE.md)
also allows distinct, diverging heights with independent
multiplicities. The
[count-simplex test](MULTIVARIATE_AFFINE_TAIL_COUNT_POLYTOPE.md)
gives a finite rational criterion for any fixed affine template.
It compares the two linked Newton slope multisets over every
possible vector of type fractions. It includes overlapping slope
valuations and supplies a uniform exact height bound whenever the
test passes. For large counts the test can be restricted to the
signed-gap hyperplane, with a polynomial resultant check for
finitely many smaller count vectors. The
[retained-scale theorem](MULTIVARIATE_ARBITRARY_RETAINED_SCALE_DECOUPLING.md)
separates an occupied band from arbitrary independently placed
rows when the retained logarithmic weight is too small to
absorb the band marks. The band's gaps and the number of
different row types may grow. A band carrying positive
curvature excludes retained heights of order \(o(M/D)\) at
height \(M\) and gap scale \(D\). For a certified bounded-gap
affine band the retained height range improves to \(o(M)\).
The [common-slope budget](MULTIVARIATE_RETAINED_COMMON_SLOPE_BUDGET.md)
subtracts the exact linear valuation drift shared by a retained
row in both orders. It charges the maximum depth at each prime
rather than all factors in a long ramp. The resulting matching
test applies to arbitrary selected subsets, and it permits
retained rows with independent large heights and gaps even
when their curvature stays positive.
The [cyclic discrepancy theorem](MULTIVARIATE_RETAINED_CYCLIC_DISCREPANCY.md)
gives the exact retained valuation differences by residue
counts at each prime power. Complete residue packets can have
zero effect on a prime's shared slopes while their individual
heights vary independently. The [envelope mark theorem](MULTIVARIATE_UNIVERSAL_ENVELOPE_MARKS.md)
selects a literal factorial factor in an adjacent pair at every
even terminal for arbitrary ordinary rows. It is absent on at
most two adjacent indices of the selected parent. The theorem
classifies the exact odd-terminal exception to this window rule
and charges an explicit prime-overlap budget to any common factor.
The [aligned unit-gap test](MULTIVARIATE_ALIGNED_UNIT_GAP_CENTRAL_NEWTON.md)
reduces the high-count Newton polygons to fixed central faces
when all smaller exponents satisfy the stated congruence.
The [prepared core theorem](MULTIVARIATE_PREPARED_ALIGNED_TAIL_ALL_DEGREE.md)
constructs arbitrarily weak balanced extensions of any core
whose central slopes separate simultaneously at every low pair.
Every sufficiently numerous aligned unit-gap tail then gives
all-degree coprimality, with independently varying heights
and no restriction on its signed gap.
The [endpoint-aligned preparation theorem](MULTIVARIATE_PREPARED_ENDPOINT_ALIGNED_TAILS.md)
uses a fixed balanced unit-gap extension and permits any tail
count, including the empty tail. At a witness prime it treats
endpoint or neutral ramps of bounded reduced gap with independently
varying row heights and \(p\)-power inflations. A residual oscillation
budget also permits arbitrary additional rows. With opposite
slope orderings at two primes, one arbitrary bounded-reduced-gap
row is allowed at unbounded valuation depth. Arbitrarily many
synchronous marks can be allowed when they aggregate into the
same single coefficient step at both orders. The proof covers
ordinary tails through positive exceptional Gaussian limits for
the prepared core. When synchronous rows span multiple cutoff
pairs, arbitrary prescribed finite-depth screens at finitely
many large primes can retain the same quadratic on both central
faces. A separate witness prime gives all-degree coprimality
for those ordinary profiles. The residue certificate requires
preparation and does not resolve arbitrary unmodified profiles.
The [oriented two-slope budget](MULTIVARIATE_ORIENTED_NATIVE_TWO_SLOPE_BUDGET.md)
gives an exact necessary divisibility on the original rows at
every adjacent pair. The increasing first blocks in either order
can be paid for only by opposite rows in that order and
same-direction rows in the other order. Even positivity permits
a second-largest slope bound. The resulting secant gcd and
adaptive prime budgets certify all degrees when the indicated
factors are insufficient at each low pair, allowing arbitrary
retained heights, orientations and counts. The
[gcd-free Newton residual test](MULTIVARIATE_GCD_FREE_NEWTON_RESIDUAL_SIEVE.md)
recovers the first common residual width at every coefficient
prime without factoring the coefficients. Its integer power
tests use the reduced slope denominators, and its block products
bound projective denominators. On odd blocks the Jacobi test
also constrains possible degree allocations across primes.
The [gcd-free Smith refinement](MULTIVARIATE_GCD_FREE_TRUNCATED_SMITH_SIEVE.md)
computes exact arbitrary-precision truncated Smith caps on
those same coprime blocks, even at wildly ramified slopes.
A separate [monotone upper-cluster theorem](MULTIVARIATE_MONOTONE_GAP_UPPER_CLUSTER_CLOSURE.md)
uses discrete convexity of the literal prime-factor counts.
It proves all-degree coprimality for large-gap upper clusters
with monotone multiplicities at each height and orientation,
unbounded gap ratios, and retained rows within the stated
logarithmic budget.
The [short-window theorem](MULTIVARIATE_INDEPENDENT_HEIGHT_WINDOW_CLOSURE.md)
also allows every exponent in a window of width
\(Ad+d\log(2+d)/(512J^3)\) to vary independently, where \(d\ge2\)
is the least tail gap. It certifies every fixed common-slope
template with gaps at least two and exhibits a signed-balanced
unit-gap template that defeats the single-mark count test.
The
[isolated two-clock theorem](MULTIVARIATE_ISOLATED_TWO_CLOCK_DIFFUSE_EXCLUSION.md)
excludes exact collisions near an isolated positive contact for
fixed reflected affine templates repeated arbitrarily many times.
It includes multiple roots and contacts of arbitrary finite order.
The [balanced affine jet theorem](MULTIVARIATE_BALANCED_AFFINE_JET_DIFFUSE_EXCLUSION.md)
extends this elimination to repeated signed-balanced affine
templates with unequal slopes, unmatched orientations and fixed
signed-gap offsets. Its two algebraic directions are calculated
from the full exact factorial jet.
The [affine divisor-rank theorem](MULTIVARIATE_AFFINE_DIVISOR_RANK_DIFFUSE_EXCLUSION.md)
uses the exact prime valuation directions of the template.
A proper projection of the degree-two common-factor locus
gives exclusion even when the first two limiting directions
have a nonisolated contact. Full divisor rank holds for
the four-type balanced unit-gap obstruction to the single-mark
test, giving all-degree closure of its repeated templates.
The [exact slope-cone theorem](MULTIVARIATE_EXACT_SLOPE_CONE_MULTIPLICITIES.md)
uses good primes to constrain independently varying counts
to a finite union of rational cones. Low-dimensional balanced
recession cones reduce to finitely many fixed counts and
exact balanced rays. Two opposed, separated affine windows
and their stated polynomial analogues close in all degrees,
with each count independent. For unit gaps the affine
constraints are exact flow conservation on a directed graph.
Acyclic graphs exclude common roots uniformly in the counts.
The [multiplicity antichain theorem](MULTIVARIATE_MULTIPLICITY_ANTICHAIN_BOUND.md)
bounds the number of exact collision profiles on any nested chain
whose whole row coefficient vectors have bounded multiplicative
rank. For finitely many row types it confines the exceptional
counts to finitely many lattice cosets in hyperplanes with
positive normals, with a bound independent of the row heights.
The [polynomial-scale theorem](MULTIVARIATE_POLYNOMIAL_SCALE_DIFFUSE_EXCLUSION.md)
extends divisor rank to polynomial heights and polynomial weights
on multiplicities. It applies when several distinct height scales
each contribute positive limiting curvature.
The [polynomial native-mark theorem](MULTIVARIATE_POLYNOMIAL_HEIGHT_NATIVE_MARKS.md)
uses effective bounds on smooth parts of polynomial values.
A separating nonlinear divisor excludes collisions for arbitrary
counts and any fixed ordinary core. The divisor exists generically
when one specified nonlinear type is present.
The [first-block depth theorem](MULTIVARIATE_FIRST_BLOCK_DEPTH_AND_MOVING_RESIDUES.md)
compares the deep part of a literal first block with the maximum
opposite-order depths and the depths of all other rows. It
includes a joint test using every row and a gcd certificate
computed from small remainders. Its all-degree residue classes
allow the number of distinct types and all their heights to
grow, with independent multiplicities at the stated widths.
The [signed jet tests](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md)
exclude further positive exceptions under effective signs on the
signed-bias boundary or at even-order stationary contacts. The
[clock transversality tests](MULTIVARIATE_EXCEPTIONAL_GAUSSIAN_CLOCK_TRANSVERSALITY.md)
give additional exclusion regions for reflected tails and for
some common-height tails with nonzero signed gap.
The [opposed independent-height theorem](MULTIVARIATE_OPPOSED_INDEPENDENT_HEIGHT_RATIONAL_BIAS.md)
rules out every simple transverse positive rational-bias
contact in the strict two-type curvature cone. Fixed signed
gap and integrality force a candidate sequence onto an exact
affine ray. A first-factor depth argument excludes all
sufficiently high two-type affine families in every degree,
uniformly in the two counts.
The [paired prime-mark descent](MULTIVARIATE_PAIRED_PRIME_MARK_AFFINE_DESCENT.md)
forces exact affine alignment for independently varying
opposed heights when both first factors have sufficiently
small prime-power cofactors. A same-order depth estimate
also covers powers of small base primes. The resulting
all-degree closure allows growing gaps. The same interval
argument permits a growing number of clustered types within
the stated logarithmic bound. At a positive two-type contact with both shares positive, one sufficiently
dominant first mark suffices, including tangential and
multiple contacts.
The [zero-bias two-type descent](MULTIVARIATE_OPPOSED_LARGE_GAP_ZERO_BIAS_DESCENT.md)
allows independent counts and growing gaps at transverse contacts.
Equal counts, vanishing projected skew response, sufficiently close
gaps and sufficiently deep single-block marks each give an exclusion.
The same note includes tangential tests under fixed count ratios.

An [interior persistence theorem](MULTIVARIATE_DIFFUSE_INTERIOR_CONTACT_PERSISTENCE.md)
gives a conditional obstruction to removing all the exceptional
neighborhoods by a positive separation margin. At any simple
interior reference contact whose angle difference changes sign,
arbitrary ordinary tails can approach contact with fixed signed
gap, retained weak rows, and prescribed finite prime data.
Real baselines with integer multiplicities can give exact
contact in that setting. General ordinary tails near such a
reference contact still require an exact integer exclusion or
a proof that the reference contact itself is absent.
The [finite truncation theorem](MULTIVARIATE_COMPACTIFIED_INTERIOR_CONTACT_FINITE_REDUCTION.md)
reduces sign-changing simple interior contacts of a countable
retained core to finite ordinary cores. It also unfolds some
stationary contacts using exact cubic, two-clock or quartic
row responses, with the total signed gap fixed.
The [stationary interior unfolding](MULTIVARIATE_INTERIOR_ANGULAR_TANGENCY_UNFOLDING.md)
also covers vanishing cubic response when its variation in
strength has sufficiently low order. A second-order
discriminant test and a reflected-row test use a controlled
transfer of Gaussian strength to ordinary rows. The
[third-order boundary theorem](MULTIVARIATE_THIRD_ORDER_BOUNDARY_AND_INTERIOR_REDUCTION.md)
excludes further exact ordinary tails at saturated-bias
contacts and reduces boundary contacts to finite interior
ones under the specified ordinary-row or angular and radial
response tests.
The [height quadrature theorem](MULTIVARIATE_SATURATED_BOUNDARY_HEIGHT_QUADRATURE.md)
reduces every transverse saturated-bias contact at \(n,N=n+1\)
with \(n\ge3\) to a response on at most \(N\) aligned height
nodes. The positive quadrature preserves moments and bounds
the loss from dispersion. For a finite ordinary core the
one-height sign, when needed, is determined by at most
\(n+N-1\) Taylor coefficients.
The [multiple-root unfolding](MULTIVARIATE_MULTIPLE_ROOT_INTERIOR_HEAT_UNFOLDING.md)
uses exact reflected ordinary rows to compare local Hermite
and quartic heat polynomials. It reduces the stated multiple
interior contacts of countable cores to finite ordinary
cores with simple transverse interior contacts.
The [reverse heat-shift construction](MULTIVARIATE_REVERSE_CLOCK_ROW_REPLACEMENT.md)
replaces one weak retained row with two matched ordinary rows.
It handles the lower-double and upper-triple contact, completing
all four double and triple multiplicity pairs for infinite
retained cores.
The [stationary midpoint theorem](MULTIVARIATE_STATIONARY_MIDPOINT_SHIFT_REDUCTION.md)
uses two third-order angular responses at a simple contact
where both velocities and the leading cubic response vanish.
Its fourth-order test shifts an ordinary midpoint by one unit
when a fixed signed gap occurs infinitely often.
The [angular rank reversal theorem](MULTIVARIATE_MULTIPLE_CONTACT_ANGULAR_RANK_REVERSAL.md)
compares the two sides of the logarithmic heat split at a
multiple root. It reduces arbitrary multiplicities, even
with one simple root, whenever at least one regularized
angular velocity is nonzero. A third-order rank test covers
further cases. Arbitrarily weak ordinary rows make the
resulting finite-core contacts simple and transverse.
The [integer row split](MULTIVARIATE_INTEGER_ROW_SPLIT_UNFOLDING.md)
gives a local differential operator from an ordinary
integer replacement of a retained row. It resolves the
stated multiple contacts even when the angular velocities
vanish. The [paired multiple-root theorem](MULTIVARIATE_PAIRED_MULTIPLE_UNBOUNDED_GAPS.md)
completes the reduction whenever the retained gaps are
unbounded. At each fixed pair the proof only needs infinitely
many gaps above a computable bound, by taking the maximum
of the thresholds over the possible multiplicities.
The [midpoint and cloud theorem](MULTIVARIATE_PAIRED_MULTIPLE_MIDPOINT_CLOUDS.md)
removes the retained-gap threshold altogether. Every
paired multiple positive interior contact of an infinite
ordinary retained core reduces to simple transverse
contacts of finite ordinary cores. A single weak retained
row can have unit gap.
The [first complex jet theorem](MULTIVARIATE_SIMPLE_CONTACT_COMPLEX_JET_UNFOLDING.md)
uses the same midpoint and cloud operators at simple
stationary contacts of infinite cores. Hermite sign control
gives finite-core reduction when the first complex order
is at least four, and when the order-three weighted response
is nonzero. The remaining simple even case has first
complex order three with weighted cubic cancellation.
The [mixed multiplicity jet theorem](MULTIVARIATE_MIXED_MULTIPLICITY_COMPLEX_JET_UNFOLDING.md)
also reduces several mixed simple and multiple interior contacts
to finite transverse references. Its remaining infinite-core
mixed cases have even multiple-root multiplicity and an earlier
even jet or a tied second jet, subject to the explicit double-root
exclusion in that note.

The [prime-height theorem](MULTIVARIATE_DIFFUSE_CLOCK_PRIME_HEIGHT.md)
uses the angular and physical resultants in the same rational
factorial coordinates. If their monomial ratios concentrate
enough logarithmic height on any fixed finite prime set, both
limiting clocks have a common rational power from an explicit
finite list. At other contacts, an actual collision sequence
must place a positive fraction of its height outside every
fixed finite prime set on a resultant monomial ratio.
The arithmetic test allows arbitrary ordinary tail counts,
gaps and heights.
The [projective factor sieve](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md)
tracks a common factor's root ratios independently of the
physical normalization. Its exact denominator is determined
by the spread of the shared Newton slopes and is bounded by
a gcd of scale-invariant adjacent coefficient ratios. Algebraic
contact confinement then requires a growing spread at irrational
limiting shapes. The same shape also constrains the product
of primes with binomial stepped residues, even as those
primes vary with the tail.
The [depth-weighted quotient test](MULTIVARIATE_DEEP_BINOMIAL_PROJECTIVE_QUOTIENTS.md)
charges the full prime-power depth of binomial Newton faces
at several shared valuations simultaneously. Its rational
invariant distinguishes actual torsion ratios from those
that only approach a torsion value at a diffuse contact.
For irreducible common factors the
[short resultant bound](MULTIVARIATE_PROJECTIVE_BINOMIAL_PRIME_RESULTANTS.md)
uses the distinct projective conjugates and their discriminant
to improve the moving-prime product bound. The deep uniform
case admits a single resultant with a prime-power bound.
The [oriented block gcd bounds](MULTIVARIATE_ORIENTED_BLOCK_DENOMINATOR_GCD.md)
give sharper exact denominator budgets from adjacent factorial
blocks, with increasing and decreasing rows accounted for
separately. For two opposed types of independent heights, the
same-orientation contributions at large primes come only from
the short integer strips at their smaller exponents.
The [truncated Newton face theorem](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md)
compares the two faces in their common physical coordinate
above residue-field precision. Smith valuations bound the
common width at each shared root valuation. Binomial faces
give an exact power-ratio congruence at the full available
depth. In the synchronized terminal strips this gives a
further denominator gcd. Explicit congruence classes supply
all-degree examples with unbounded distinct rows approaching
any prescribed positive interior Gaussian parameter.
The [simultaneous contact height theorem](MULTIVARIATE_SIMULTANEOUS_CONTACT_HEIGHT.md)
uses the \(p\)-adic Subspace Theorem on independent limit
coordinates of a common factor. Its necessary denominator
growth estimate allows moving primes and separately charges
approximation to torsion configurations at fixed deep
binomial faces. Its exclusion neighborhood is qualitative.
The [moving-prime trace resultant](MULTIVARIATE_MOVING_PRIME_TRACE_RESULTANT.md)
extends the scale-invariant denominator budget to algebraic
trace targets certified by nonbinomial truncated faces.
One integer resultant charges the moving prime depths and
the real error at a matching diffuse reference. Its bound
continues to apply after removing exact target traces.
The [fixed-field trace orbit theorem](MULTIVARIATE_FIXED_FIELD_TRACE_ORBIT_HEIGHT.md)
relates the denominator of each unordered Galois trace orbit
to the height of an ordered root quotient. A finite field
menu strengthens the approximation exponent, and local
decomposition groups restrict the allowed moving-prime
valuation spread separately for each orbit.
The [Hermite trace height bound](MULTIVARIATE_TRACE_HERMITE_JET_HEIGHT.md)
uses the small successive derivatives when several projective
traces cluster. The lower bound allows the trace field and the
denominator primes to vary, and it can improve the squarefree
discriminant bound when there is more than one cluster. The
[moving-prime projection](MULTIVARIATE_MOVING_PRIME_TRACE_PROJECTION.md)
can charge local depths and real approximation together when
the local target coefficients lie in a fixed rational subspace
and the projected real reference is nonzero. The
[denominator-prime discriminant](MULTIVARIATE_TRACE_DENOMINATOR_PRIME_DISCRIMINANT.md)
uses primes where several traces already have negative valuations.
Their valuation spreads and normalized unit coincidences decrease
the available height budget. It allows a Galois-invariant selection
of trace pairs and needs no algebraicity of their real limits.
The [directional quadratic height theorem](MULTIVARIATE_TRANSVERSE_QUADRATIC_DIRECTIONAL_HEIGHT.md)
uses the two exact contact equations to constrain the trace
and scaled rational norms to one first-order direction.
The resulting smaller cross-direction error sharpens the
height requirement. For a fixed quadratic field, signed angle
coordinates replace the spread denominator by numerator budgets.
The [principal divisor height theorem](MULTIVARIATE_QUADRATIC_PRINCIPAL_DIVISOR_LOG_HEIGHT.md)
uses the split prime support and one ideal-class factor to bound
approximation to an algebraic direction and to a scaled root.
For bounded fields and fixed split prime support, a nonattained
contact needs a trace denominator exceeding every polynomial
bound in reciprocal diffuseness. There is a moving-prime bound
with explicit support dependence.
The [single-order norm gap](MULTIVARIATE_HIGH_ORDER_NORM_DENOMINATOR_GROWTH.md)
uses a signed high-order factorial difference to construct two
scaled rational norms with a strictly positive small difference.
The [paired norm gap](MULTIVARIATE_PAIRED_FACTORIAL_DIFFERENCE_NORM_HEIGHT.md)
compares the two factorial grids and cancels their common leading
term. It yields the effective error scale
\(\eta^{r-2}Z+W\eta^r\) for every \(2\le r\le n\), with
no reference algebraicity or transversality assumption. These
tests charge the sum of the selected local root valuations rather
than their spread. The [evaluated norm bound](MULTIVARIATE_EVALUATED_FACTOR_NORM_HEIGHT.md)
instead evaluates both parents at the algebraic limiting root.
For a common quadratic approaching it, the prime-to-fixed-set
part of the gcd of these two evaluated resultants must grow at a
specified rate. The bound applies in reference degree at least
five, and in degree four under the stated rank or field condition.
Its sharper budget counts only matching prime ideals of the two
evaluations and at most two total residue degrees at each good
prime. The [ordinary norm amplification](MULTIVARIATE_TWO_NORM_BUDGET_AMPLIFICATION.md)
constructs arbitrarily weak balanced extensions that preserve the
starting tail's maximum relative gap. At a fixed prime the same
extensions make every paired denominator exponent and the common
evaluated ideal valuation arbitrarily large, with two common unit
slots surviving any prescribed finite precision. The actual
extensions have disjoint root arguments in every degree. A uniform
upper bound on these finite-screen norm budgets therefore needs
information beyond the fixed real controls and tested precision.

The [uniform bias theorem](MULTIVARIATE_UNIFORM_DIFFUSE_BIAS_SEGMENTS.md)
allows one chosen signed-drift intercept to work throughout a
prescribed compact strength interval, with no row-count restriction.
The [countable core height theorem](MULTIVARIATE_COUNTABLE_CORE_CURVATURE_HEIGHT_SEPARATION.md)
compares omitted curvature to the factorial height of every
truncation at positive algebraic strength. It gives uniform
all-degree transfer for families with the stated decay and
without a retained row-count bound. For tails whose curvature
decays faster than every inverse power of the truncation height,
the exceptional strengths across all degrees are nonalgebraic
and form a set of Hausdorff dimension zero.
The [generic diffuse closure theorem](MULTIVARIATE_DIFFUSE_RATIONAL_SLOPE_GENERIC_CLOSURE.md)
excludes physical contacts on prescribed algebraic curves of the two
clock parameters for generic ordinary cores with sufficiently many
independent base coordinates. In particular it gives nonadaptive
stability along rational signed-drift slopes, including zero.
The [nested Newton clock sieve](MULTIVARIATE_RATIONAL_BIAS_NESTED_NEWTON_CLOCKS.md)
applies at a prescribed rational bias slope for an arbitrary
fixed ordinary core. The strict vertices disappear only once
as the clock valuation varies, so finitely many active slope
comparisons can exclude the entire positive reference family.
The [generic transverse quadratic reduction](MULTIVARIATE_GENERIC_TRANSVERSE_QUADRATIC_DIFFUSE_CONTACTS.md)
shows that sufficiently diffuse collisions on generic wide cores
would have a quadratic gcd, with a transverse reference and a
single simple conjugate pair. Arbitrarily weak balanced integer
rows give such reference contacts near any specified transverse
interior reference. A prescribed proper algebraic locus of
quadratic coefficients can be excluded generically.
The [third-order matched angular response](MULTIVARIATE_MATCHED_TANH_CUBIC_CORE_ANGLE_RESPONSE.md)
resolves some weak-tail contacts of simple core roots where the
quadratic response reaches zero. The
[two-clock quotient resultant](MULTIVARIATE_GAUSSIAN_TWO_CLOCK_RESULTANT_AT_INFINITY.md)
determines its exact leading terms and separates independently
deformed pairs when just one clock is algebraic, or when the
two transcendental clocks satisfy the stated rational power relation.
The [prime-adic clock bounds](MULTIVARIATE_GAUSSIAN_CLOCK_PADIC_EXCEPTIONS.md)
constrain every algebraic exceptional Gaussian parameter at every
prime and reduce rational exceptions to an explicit finite list.
The [Gaussian endpoint transfer](MULTIVARIATE_MIRRORED_CORE_GAUSSIAN_ENDPOINT_TRANSFER.md)
gives separation at every strength for mirrored cores satisfying
the exact endpoint tests, including any fully shared mirrored core
with a pair \(d\ge7\ell\). A uniform angular gap permits small relative
gaps in an arbitrary asymmetric ordinary tail. The
[large mirrored pair theorem with arbitrary rows](MULTIVARIATE_HARD_MIRRORED_PAIR_ARBITRARY_TAIL.md)
allows any additional ordinary rows, including private and unmatched
rows, whenever one shared mirrored pair has \(d\ge7\ell\).
It also allows any common positive Gaussian strength.
A [retained mirrored subset test](MULTIVARIATE_RETAINED_MIRRORED_ENDPOINT_CRITERION.md)
extends this proof to products of exact mirrored endpoint factors.
The [oriented endpoint product theorem](MULTIVARIATE_ORIENTED_ENDPOINT_PRODUCT_CRITERION.md)
uses the Gamma chord of each shared row in either orientation and
permits asymmetric retained profiles. One qualifying shared row can
certify the whole profile, and several smaller rows can meet the
endpoint tests collectively. The tests persist under arbitrary
added ordinary rows. The [opposed-gap theorem](MULTIVARIATE_OPPOSED_SHARED_GAPS_ENTROPY_INTERLACING.md)
requires two shared rows of opposite sign with individual relative
gap at least five. The more general
[binomial entropy criterion](MULTIVARIATE_ENTROPY_ENDPOINT_PHASE_CRITERION.md)
allows any collection of retained shared rows satisfying the endpoint
and derivative bounds, including examples with every relative gap
below five. Both allow arbitrary additional ordinary rows and common
Gaussian strength. The [absolute curvature theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_16.md)
handles the Gaussian-deformed pair for any ordinary core when its
curvature plus Gaussian strength is at least sixteen. The
[weighted endpoint refinement](MULTIVARIATE_DOUBLED_SHARED_CURVATURE_ENTROPY.md)
also covers absolute curvature from fourteen when enough comes from
shared rows with the larger exponent at least twice the smaller.
The [radial comparison for mirrored cores](MULTIVARIATE_MIRRORED_CORE_SIGNED_DIFFUSE_RADIAL_SEPARATION.md)
gives an absolute diffuse-tail tolerance when the exact center
shift has a signed lower bound. For a total-degree-balanced profile
it suffices that the tail log ratio has magnitude at least a fixed
positive multiple of its curvature. This covers arbitrary mirrored
cores in every degree, including cores with private rows.
In the weak-curvature window the [fourth-order angular model](MULTIVARIATE_MATCHED_SECOND_ORDER_ANGLES.md)
preserves the reference interlacing. The remaining exact phase
comparison would follow from a paired angular remainder smaller than
\(y/(32n^2)\) at the neighboring Hermite anchors in its stated window.
The [Gaussian Laplace response](MULTIVARIATE_GAUSSIAN_LAPLACE_ANGLE_RESPONSE.md)
gives a degree-uniform bound on the first angular variation of
the exact matched curvature difference at a Gaussian reference.
The [nonlinear fixed-grid transfer](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md)
treats the positive jump multiplier for that exact difference
and bounds angular motion by \(O(n^{-2})\) at a sufficiently
spaced reference for fixed shared rows under its rate condition.
[Compensated off-circle transport](MULTIVARIATE_OFF_CIRCLE_COMPENSATED_LAPLACE_TRANSFER.md)
extends the same-grid argument to asymmetric factorial increments at
noncircular reference roots, under a thin-annulus or a paired-root
cancellation condition.
The [sharpened mesh estimate](MULTIVARIATE_MATCHED_LAPLACE_MESH_SHARPENING.md)
gives this nonlinear transfer at sufficiently high order for every
fixed fully shared mirrored profile.
For weak mirrored profiles the [symmetric same-grid transfer](MULTIVARIATE_UNIFORM_WEAK_MIRRORED_LAPLACE_TRANSFER.md)
gives a uniform bound of order \(W^{3/2}(1+\log n)/n^2\) for the exact
matched angular displacement at every degree.
The [weak common-symbol theorem](MULTIVARIATE_UNIFORM_WEAK_COMMON_SYMBOL_INTERLACING.md)
compares the two degrees at the same exact upper symbol for fully shared
mirrored profiles below an absolute computable weak cap. It uses a
quantitative gap for the actual weak roots and the symmetric same-grid
bound.
The [asymmetric weak comparison](MULTIVARIATE_WEAK_ASYMMETRIC_COMMON_SYMBOL_ANGLE_SEPARATION.md)
uses the exact midpoint law and a shared analytic order interpolation.
It proves a quantitative all-degree gap for fully shared profiles with
independent blends of the two symbols on either grid. The
[large-curvature blend theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_27_COMMON_SYMBOL.md)
gives the same blend freedom for fully shared profiles at \(S\ge27\)
at every order, with no symmetry requirement.
On compact positive shared bands separated from stationary
criticality, the
[asymmetric fixed-band theorem](MULTIVARIATE_FIXED_SHARED_BAND_COMMON_SYMBOL_ANGLE_SEPARATION.md)
gives a uniform eventual angular gap for all independent blends.
It controls both root radii and the imaginary packet errors.
The [asymmetric cusp theorem](MULTIVARIATE_NEAR_CRITICAL_SHARED_COMMON_SYMBOL_ANGLE_SEPARATION.md)
covers the remaining near-critical band. With the weak and high
curvature theorems it gives one absolute eventual angular threshold
for unrestricted fully shared ordinary profiles, for both actual
symbols and for independent blends on the two grids. In every [fixed subcritical mirrored band](MULTIVARIATE_SUBCRITICAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md)
\(W\le W_+<4\), there is also a uniform large-order threshold for strict
interlacing of the actual and common-upper-symbol pairs.
The [supercritical mirrored comparison](MULTIVARIATE_SUPERCRITICAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md)
gives uniform eventual strict interlacing on each compact band above
\(W=4\), including the entire interpolation between the actual and
upper symbols on the lower grid.
The [moving critical mirror theorem](MULTIVARIATE_CRITICAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md)
covers \(|W-4|\le K/\sqrt n\) for any fixed \(K\) at uniformly
large orders, with the same strict ordering along that interpolation.
The [uniform mirrored theorem](MULTIVARIATE_UNIFORM_EVENTUAL_MIRRORED_COMMON_SYMBOL_INTERLACING.md)
covers the remaining intermediate thermal scale at effective order
\(n|W-4|^2\). It gives one eventual threshold for all ordinary
fully shared mirrored profiles, at unrestricted curvature and for
every lower-symbol interpolation parameter.
The compression analysis quantifies the nonlinear defect introduced by the
usual interlacer,
and the [positive-jump obstruction](MULTIVARIATE_POSITIVE_LAPLACE_ANGLE_OBSTRUCTION.md)
gives an exact inward angular motion and a separate finite loss
of circle location under an unlinked positive direction.
The [linear Cayley compression obstruction](MULTIVARIATE_LINEAR_CAYLEY_COMPRESSION_OBSTRUCTION.md)
excludes universally interlacing degree-lowering linear operators
even for the actual adjacent moments of one overlapping mirrored
unit gap. For a shared mirrored gap-two pair, the
[exact section interlacer](MULTIVARIATE_MIRRORED_GAP_SECTION_INTERLACER.md)
gives the comparison signs at every lower root from a padded
Gegenbauer section. For the ordinary pair \((1,3),(3,1)\) a
[direct differential comparison](MULTIVARIATE_MIRRORED_CRITICAL_GAP_TWO_INTERLACING.md)
gives the actual upper signs in every degree. Its
[angular reserve](MULTIVARIATE_MIRRORED_GAP_TWO_FORCE_RESERVE.md)
bounds the minimum adjacent gap below by \(10^{-10}(n+1)^{-2}\).
The
[adjacent Gaussian smoothing theorem](MULTIVARIATE_GAUSSIAN_ADJACENT_SMOOTHING_OBSTRUCTION.md)
proves universal interlacing preservation for equal angular heat times.
At unequal rates, including those from a common normalized quadratic
symbol on the two grids, abstract positive reciprocal log-convex
inputs can lose their adjacent interlacing.

For a mirrored profile the [polar-half criterion](MULTIVARIATE_MIRRORED_POLAR_HALF_ANGULAR_REDUCTION.md)
gives nonvanishing halves and increasing phases at every curvature.
It reduces the unresolved angular comparison to a relative
signed estimate for one exact matched multiplier at lower roots.
The [Euler factorization](MULTIVARIATE_MIRRORED_POLAR_EULER_FACTORS.md)
separates the simultaneous advance of the Gamma shapes from
the one-sided residual factors, keeping each residual attached
to its own upper Gegenbauer factor. Its exact counterexample
shows why exterior stability of the residual and coefficient
log convexity alone give insufficient control of the sign.
The [two-slot force comparison](MULTIVARIATE_TWO_SLOT_MIRRORED_FORCE_COMPARISON.md)
proves one direction of the adjacent interlacing inequalities
for two shared mirrored unit-gap pairs or one shared gap-two
pair. The other direction requires a bound on the actual
displacement under the Gamma shape advance.

Across all terminals, the [absolute curvature twelve theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
and the [uniform weak-curvature theorem](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md)
confine a possible counterexample to \(\eta_*/2<S<12\), where
\(\eta_*>0\) is a computable absolute rational. The
curvature-twelve proof uses the one-sided phase contribution of
the first coefficient, Stirling prefactors in the tails, and
Gamma argument bounds that retain the sign of each row. The
[curvature thirteen theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_13.md)
supplies the separate cubic curvature budgets, allowing the two
orientations to use different ranges of the linked phase proof. The
[absolute curvature fourteen theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_14.md)
introduces the directed quadratic budget for short halves. The curvature-fifteen
bound combines an improved largest-gap budget, a linked Gamma derivative
including private ending rows, and weighted radial and angular estimates. The
[curvature-twenty argument](MULTIVARIATE_ABSOLUTE_CURVATURE_20.md)
uses the exact endpoint cost to absorb excess Gamma curvature and
sharpens the differentiated connected logarithm bound.
The curvature-sixteen refinement also bounds the product of
the largest relative gap and the endpoint coefficient scale.
Both proofs allow a common Gaussian strength and arbitrary private
and asymmetric rows.

The [square-factorial phase criterion](MULTIVARIATE_MIXED_PRIVATE_SQUARE_FACTORIAL_ANGLES.md)
covers profiles with a private variable on both sides and total
private exponent at least three on either side, with disjoint
root arguments at every order. Two extra private factorial
factors give a squared factorial coefficient envelope, and the
signed sectors use only scalar bounds at the small orders.
A [shorter phase bound](MULTIVARIATE_MIXED_PRIVATE_DEGREE_THREE.md)
treats all lower orders at least twelve. The [prime-power private theorem](MULTIVARIATE_PRIME_POWER_PRIVATE_TANH_ANGULAR_STABILITY.md)
uses the exact modular root quotients for single private rows of
same-prime power exponents, including exponent one. An explicit
cap on the shared matched tanh rates gives all-degree stability
under asymmetric shared tails.

The [private totals one and two theorem](MULTIVARIATE_PRIVATE_ONE_TWO_ALL_DEGREES.md)
gives disjoint consecutive root arguments in every degree for
both realizations of the total two, with unrestricted shared
rows. The [endpoint Euclidean descent](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md)
gives physical coprimality in every degree for any both-private profile
with \(1-P_a\le\Delta\le P_b-1\), where \(P_a,P_b\) are the private
totals. In particular, every balanced profile with private variables
on both sides is covered with no cutoff or curvature hypothesis.
The [first-sector coefficient-ratio theorem](MULTIVARIATE_BOTH_PRIVATE_ALL_DEGREES.md)
excludes the physical first box for every both-private profile.
It uses one side of the exact adjacent Beta quotient and the
single sign change of the imaginary coefficient terms. Together
with [higher-label exclusion](MULTIVARIATE_PRIVATE_UNIT_ALL_HIGHER_LABELS.md)
for two private units and the three larger private-total cases,
this proves ordinary SFC(2) for all profiles with private variables
on both sides. The same first-sector argument applies when only
one column has private rows and it has strictly higher total degree.
The [polar derivative gate](MULTIVARIATE_POLAR_DERIVATIVE_FIRST_SECTOR.md)
also gives first-sector exclusion for all comparable exponent columns
and for mixed shared profiles satisfying its exact coefficient
comparison or its curvature and imbalance bound. The
[direct coefficient ratio](MULTIVARIATE_DIRECT_RATIO_FIRST_SECTOR.md)
also separates the entire first sector for fully shared profiles
with \(V\le16/5\), including balanced columns and arbitrary orientations.
Its [normalized curvature theorem](MULTIVARIATE_DIRECT_RATIO_CURVATURE_AND_SECTOR_COLLAR.md)
gives a single possible change of sign in the exact adjacent log
ratio's second derivative and discrete second differences in one
orientation of every survivor. For \(N\ge14\), \(S+\tau\le12\),
and effective signed gap at most \(-2\) in that orientation,
a two-cut sine comparison excludes an explicit interval on both
sides of \(2\pi/N\). This includes ordinary fully shared profiles
with absolute total degree difference at least two after reversal.
The [joint positive moment obstruction](MULTIVARIATE_JOINT_POSITIVE_MOMENT_SHAPE_OBSTRUCTION.md)
realizes exact adjacent common roots together with this monotone
normalized curvature, strict direct log convexity, and any fixed
number of positive even matched derivatives at sufficiently large
orders. An exact affine Euler factor and increasing companion
ratio also permit simple common roots under a positive joint law.
The [smooth log-concave obstruction](MULTIVARIATE_LOG_CONCAVE_ADJACENT_MOMENT_OBSTRUCTION.md)
uses a single concentrated positive ratio law with rational
moments and the exact derivative relation. For every upper order
at least five it produces one common irreducible quadratic.
These relaxations isolate the need for further factorial
constraints in positive moment and direct shape arguments.
The [exact factorial rectangle](MULTIVARIATE_FACTORIAL_RECTANGLE_FIRST_SECTOR.md)
uses the stationary curvature to obtain explicit first-sector
angular separation on each strictly subcritical fully shared band.
It also treats any relative radii in the first sector for one
private unit with \(V\le4/5\).
The [critical coefficient concentration theorem](MULTIVARIATE_SUBCRITICAL_FIRST_SECTOR_CONCENTRATION.md)
gives a uniform zero-free sector for individual fully shared moments,
including stationary criticality. The
[two-sector coefficient chord](MULTIVARIATE_TWO_SECTOR_COEFFICIENT_CHORD.md)
uses consecutive arithmetic and harmonic ratio bounds to account
for an additional sine sign change. Its all-order ordinary application
covers balanced shared profiles under signed midpoint bounds.
The [joint matched phase comparison](MULTIVARIATE_PRIVATE_UNIT_MATCHED_PHASE_HIGH_LABELS.md)
proves separation under independent positive rescalings in high angle
boxes. It also treats every angle when \(V\le1/n\), or when
\(V\le13/(5n)\) at lower order at least thirteen.

For two private units the
[curvature fifteen-eighths theorem](MULTIVARIATE_BOTH_PRIVATE_UNIT_CURVATURE_15_OVER_8.md)
gives all-degree angular separation at shared absolute
curvature at least \(15/8\), using a triangular Gamma chord
bound and exact first-box ordering of the linked Beta
coefficient quotient. In the weak range the
[attenuated scalar comparison](MULTIVARIATE_PRIVATE_UNIT_WEAK_CURVATURE_SEVENTEEN.md)
gives physical coprimality for all signed gaps when the
shared smaller-endpoint curvature is at most \(17/J\),
where \(J\) is the absolute angular cutoff enlarged to
at least \(10^8\).
The [absolute weak-curvature bound](MULTIVARIATE_PRIVATE_UNIT_ABSOLUTE_WEAK_RADIAL_DRIFT.md)
also proves physical coprimality at \(V\le10^{-11}\), with
\(V\) the shared smaller-endpoint curvature, without a terminal
cutoff. The [signed-skew bound](MULTIVARIATE_PRIVATE_UNIT_NONWEAK_SKEW_RADIAL.md)
permits unrestricted shared curvature when the unmatched signed
skew and its unsigned remainders meet its \(A/V\le10^{-15}\)
criterion. For arbitrary mixed skew, the
[absolute imbalance theorem](MULTIVARIATE_PRIVATE_UNIT_IMBALANCE_TWENTY_FIVE.md)
bounds any common physical root in the two-private-unit class
to a signed degree difference of absolute value at most 24.
Its radial walls require no eventual terminal cutoff.

The [uniform private-unit angular reserve](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_UNIFORM_ANGULAR_RESERVE.md)
gives a cross-degree angular gap of order \(N^{-2}\)
and a derivative lower bound independent of curvature
for every mirrored private-unit core. These bounds
give all-degree stability under polynomial-size weak
asymmetric tails and midpoint-matched Gauss slots.
The [quadratic near-mirror test](MULTIVARIATE_PRIVATE_UNIT_NEAR_MIRROR_QUADRATIC_ANGLES.md)
allows opposite slots matched across different integer row gaps.
A reciprocal deformation cancels first-order angular motion and
a signed Laplace-kernel bound retains further slot cancellation.
[Balanced Gauss blocks](MULTIVARIATE_PRIVATE_UNIT_BALANCED_GAUSS_BLOCKS.md)
give a separate second-order cancellation criterion for
ordinary nonmirrored blocks, including different slopes.
The [signed skew radius comparison](MULTIVARIATE_PRIVATE_UNIT_SIGNED_SKEW_RADIAL_SEPARATION.md)
uses same-order symmetrization and the actual integer center drift
to separate the physical radii in every degree under its matched
skew and drift conditions. It needs no eventual cutoff.

The [paired Newton-binomial theorem](MULTIVARIATE_PAIRED_NEWTON_BINOMIAL_ANGLES.md)
separates arguments using one prime at both orders, even
when the endpoint slopes are horizontal. It also gives
all-degree families through simultaneous first-block hits
on movable shared rows. The [long Newton edge theorem](MULTIVARIATE_LONG_NEWTON_CYCLE_ANGULAR_GATE.md)
uses transitivity and a long tame cycle to give the
unordered-pair angular criterion.
For exactly circular arrays, the
[reciprocal paired-slope test](MULTIVARIATE_RECIPROCAL_PAIRED_SLOPE_CYCLOTOMIC_GATE.md)
accounts for opposite slopes jointly in the common-factor degree.
When the paired allocations exclude all nonzero slopes, it reduces
the remaining obstruction to finitely many cyclotomic polynomials.
The [paired multi-edge gate](MULTIVARIATE_PAIRED_MULTI_EDGE_ANGULAR_GATE.md)
also uses disjoint differences of Newton slopes and the
residue cosets from both orders. Its all-degree example
uses two edges at the same prime.

For a constant and a monomial, the
[endpoint squarefreeness theorem](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_SEVEN.md)
with its [degree-six refinement](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_SIX.md)
and [degree-five refinement](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_FIVE.md)
covers total degrees at least five. The
[degree-four theorem](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_FOUR.md)
completes all-degree SFC(2) for a constant and any nonconstant
ordinary monomial. It uses squarefreeness of the Appell sequence and
a logarithm estimate at the nonconstant endpoint. The
[grouped middle comparison](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_FOUR_MIDDLE.md)
uses the centered exponential jet and exact scalar Gamma majorants.
The
[matched Gauss remainder](MULTIVARIATE_MATCHED_GAUSS_UNIFORM_REMAINDER.md)
also supplies the explicit sufficient angular range
\(S\le10^{-18}N^{-13/6}\). The [earlier Hermite-edge estimate](MULTIVARIATE_EFFECTIVE_WEAK_CURVATURE_HERMITE_EDGE.md)
gives a second explicit weak-curvature range for smaller orders.
For all ordinary profiles the
[absolute eventual angular theorem](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
gives a computable cutoff under independent positive rescaling.
It compares modulus equations to control the difference of real
fields before subtracting phases. At the one-private cubic scale,
the signed imaginary speed gives a uniform strict margin on each
bounded model window. The numerical cutoff is unevaluated.
The [effective stationary-fiber bounds](MULTIVARIATE_EFFECTIVE_STATIONARY_FIBER_LOCALIZATION.md)
give a finite net for bounded-curvature profiles on interior field
rectangles. The [inverse-map estimates](MULTIVARIATE_EFFECTIVE_STATIONARY_INVERTIBILITY.md)
then quantify root geometry away from critical values and localize
the fold at very small curvature. Real-lattice expansions now give computable order cutoffs inside
bounded [shared quartic](MULTIVARIATE_EFFECTIVE_BOUNDED_SHARED_QUARTIC.md)
and [one-private cubic](MULTIVARIATE_EFFECTIVE_BOUNDED_PRIVATE_CUBIC.md)
windows. The one-private proof also has effective
[quadratic](MULTIVARIATE_EFFECTIVE_PRIVATE_QUADRATIC_CORRIDOR.md)
and [escaping cubic-field](MULTIVARIATE_EFFECTIVE_PRIVATE_CUBIC_FIELD.md)
corridors. The [effective thermal endpoint](MULTIVARIATE_EFFECTIVE_PRIVATE_THERMAL_ENDPOINT.md)
and [positive-boundary entry](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_POSITIVE_ENTRY.md)
cover the small positive angle sector for a single private linear row
on a computable band \(0\le\sigma\le1+\delta\). The
[negative-boundary proof](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_NEGATIVE_BOUNDARY.md)
covers the sector near the negative real axis for \(0\le\sigma\le5/4\).
The [effective interior comparison](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_INTERIOR.md)
covers the intervening strip on the latter band, away from its real
boundaries by any prescribed positive angle. These give one
computable all-field cutoff on \(0\le\sigma\le1+\delta\).
For bounded positive shared curvature the
[effective transversality bounds](MULTIVARIATE_EFFECTIVE_SHARED_TRANSVERSALITY.md)
quantify the real-field slopes for separated regular pairs and fold
directions away from the classified singleton geometry. The
[regular contour cover](MULTIVARIATE_EFFECTIVE_SHARED_REGULAR_COVER.md)
uses actual saddle activation to exclude roots in fixed interior
supercritical bands and in the regular part of subcritical bands,
including the inactive singleton rays. The
[effective fold comparison](MULTIVARIATE_EFFECTIVE_SHARED_SUBCRITICAL_FOLD.md)
includes every Airy scale in those fixed subcritical interior bands.
The [shared real-boundary proof](MULTIVARIATE_EFFECTIVE_SHARED_REAL_BOUNDARIES.md)
completes the all-field comparisons for fixed positive shared bands
separated from criticality. The
[vanishing-curvature theorem](MULTIVARIATE_EFFECTIVE_VANISHING_CURVATURE.md)
gives a computable absolute small-curvature cutoff, including the
full physical Hermite fold. The
[both-private](MULTIVARIATE_EFFECTIVE_BOTH_PRIVATE_LINEAR.md)
and [one-private](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_REGULAR_COMPLEMENT.md)
comparisons make the large-order thresholds computable on every
private support type. The
[effective global reduction](MULTIVARIATE_EFFECTIVE_GLOBAL_REDUCTION.md)
therefore confines the remaining global cutoff problem to
fully shared profiles approaching the quartic transition,
with curvature in one fixed positive band.
The [effective critical entry](MULTIVARIATE_EFFECTIVE_SHARED_CRITICAL_ENTRY.md)
then localizes both physical fields to any prescribed neighborhood
of the moving quartic center. The
[escaping-field bound](MULTIVARIATE_EFFECTIVE_SHARED_ESCAPING_QUARTIC_FIELDS.md),
[supercritical thermal comparison](MULTIVARIATE_EFFECTIVE_SHARED_BOUNDED_THERMAL_CUSP.md),
and [subcritical thermal regular](MULTIVARIATE_EFFECTIVE_SHARED_SUBCRITICAL_THERMAL_REGULAR.md)
and [fold](MULTIVARIATE_EFFECTIVE_SHARED_SUBCRITICAL_THERMAL_FOLD.md)
comparisons cover the remaining scales. The
[absolute hierarchy](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md)
combines these estimates with the bounded quartic window.

At prime-power terminals the [balanced unit-gap Eulerian reduction](MULTIVARIATE_BALANCED_UNIT_GAP_EULERIAN_REDUCTION.md)
gives an exact bounded-degree arithmetic test for arbitrarily many
rows in each orientation. The [linked beta Frobenius reduction](LINKED_BETA_FROBENIUS_REDUCTION.md)
gives a fixed-degree modular test at upper orders \(N=hp^q\) with
\(p>h(a+b+1)\) for the unequal linked two-row beta family.
The [fixed-cofactor reduction](MULTIVARIATE_MIRRORED_FIXED_COFACTOR_FROBENIUS.md)
extends the mirrored unit-gap residue test to terminals \(mp^q\).
Its test degree depends on \(m\) and the number of row pairs,
and is independent of \(q\). The
[opposed-unit Jacobi strip](MULTIVARIATE_OPPOSED_UNIT_JACOBI_STRIP.md)
combines an upper-degree numerical range with a linked
fixed-degree remainder. The strip excludes every possible common
root for the original rows \((3,4),(2,1)\) in every degree.
For a single shared mirrored pair of any fixed gap and smaller
exponent, the [bounded Frobenius test](MULTIVARIATE_SINGLE_MIRRORED_GAP_FROBENIUS.md)
gives an explicit infinite progression of prime-power exponents
with coprimality, using an eliminant of degree independent of the
prime power. The [prime-offset descent](MULTIVARIATE_UNIT_GAP_PRIME_OFFSET_DESCENT.md)
reduces the arithmetic question at \(p+t\) and \(2p+t\) to
much smaller terminal \(t\) for same-orientation unit-gap rows
under stated size conditions. A finite residue test also handles
the shared slopes that remain outside the simple descent.
The [unbalanced prime-power criterion](MULTIVARIATE_UNBALANCED_PRIME_POWER_NEWTON.md)
bounds the common-factor degree for arbitrary mixed unit-gap bases
and gives a residue test independent of the prime-power exponent.
For base-one gap-two rows it gives explicit Frobenius exclusions
with arbitrary unequal orientation counts.
The [two-prime unit-gap theorem](MULTIVARIATE_TWO_PRIME_UNIT_GAP_DENSITY.md)
gives direct coprimality certificates of positive lower density for
actual same-orientation profiles with arbitrary bases.
The [arbitrary-gap extension](MULTIVARIATE_TWO_PRIME_ARBITRARY_GAP_DENSITY.md)
also permits one gap of any fixed size under explicit base and order
congruences.

The [growing-profile native-prime theorem](MULTIVARIATE_NATIVE_PRIME_GROWING_PROFILES.md)
proves all-degree angular separation for a proportion tending
uniformly to one in boxes of independent bases with mixed
orientations and growing row counts. The natural squared-gap
bound at the dimension-growing scale suffices, and additional
retained rows are allowed within a stated logarithmic weight.

A prescribed finite set of primes and precisions alone cannot give
a universal exclusion. The
[simultaneous prime-power construction](MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md#a-simultaneous-obstruction-at-finitely-many-primes)
produces fully shared profiles with both normalized moments congruent
to their binomial polynomials at all the prescribed moduli, within
any open positive curvature interval. In particular, the ordinary
Newton and residual screens at those moduli allow both candidate
common-factor degrees at terminal six. The
[bounded-row packing theorem](SFC_PADIC_BOUNDED_ROW_PACKING.md)
preserves all normalized coefficient valuations and units at a
fixed prime and precision through any prescribed finite range of
moments, after a single positive rational diagonal scaling. The
row bound depends on the support size, prime and precision alone.
For binary supports the same records occur at arbitrarily small
curvature, with all-degree angular separation, and the digit graph
makes fixed-prime truncated screen survival decidable in Presburger
arithmetic. The bound concerns one prime and prescribed precision.
The [simultaneous multi-order theorem](SFC_SIMULTANEOUS_MULTIORDER_FACTORIAL_RECORDS.md)
gives an exact compatibility test for constant-normalized records
at several primes across any finite set of moment orders, for any
fixed support size. At each prime the only valuation links are
the prime-power scaling chains, and the linked unit quotients
are determined by a common finite residue-count table. Every
compatible record is attained by a bounded number of ordinary
rows in the original variable. For binary profiles, an arbitrary
core can be retained and its real coefficients approximated
when the signed gap is fixed. The exact same local record can
be realized with all-degree angular separation. A shifted
complementary construction also preserves an existing binary
record at all tested orders and primes with just two weak
opposed rows, using a separate witness prime.
The [minimal fixed-prime obstruction](MULTIVARIATE_MINIMAL_SINGLE_PRIME_SCREEN_OBSTRUCTION.md)
produces arbitrarily deep binomial congruences at one prescribed
prime with a single shared row after one common variable scaling,
or in the original variable with a single mirrored pair. Both
constructions allow any prescribed positive curvature interval.
The [marked-row Eisenstein criterion](MULTIVARIATE_MARKED_ROW_VARIABLE_PRIME.md)
gives examples in every degree using a prime selected from the
exponent data. It can also extend any fixed profile by two balanced
shared rows with arbitrarily small added curvature, preserving
private degrees and total-degree difference. The
[opposite marked-row test](MULTIVARIATE_OPPOSITE_MARKED_ROW_PRIME.md)
uses two rows with overlapping below-prime windows. It
covers balanced, arbitrarily small-curvature marks even
when their common maximum prevents isolating one row at
the chosen prime. The [deep-block degree bound](MULTIVARIATE_CENTRAL_DEEP_PRIME_DIVISOR_GATE.md)
confines shared slopes to a narrow interval of coefficient indices,
even when the other array has nonzero slopes or the prime divides
the marked order. The
[common-upper-exponent theorem](MULTIVARIATE_COMMON_UPPER_EXPONENT_DEEP_MARKS.md)
gives simultaneous native divisor tests for any existing group
with a common upper exponent and both orientations, allowing
arbitrary unequal gaps and an arbitrary retained profile. The
product of central marks gives a power bound for a putative
exception in the retained exponent height. The
[large-gap common-upper theorem](MULTIVARIATE_COMMON_UPPER_LARGE_GAPS_ALL_DEGREES.md)
uses those divisors and an elementary prime-factor count to
certify original profiles in every degree when a group of
comparable large gaps dominates the remaining gap sum. The
larger exponents of the group can lie within a controlled
distance of the global maximum. The
[large-prime product theorem](MULTIVARIATE_LARGE_FACTORIAL_MARK_PRIME_BUDGET.md)
charges only the retained factors above the prime cutoff and permits
arbitrary retained multiplicity below it or in a passive top strip.
The [rational-band theorem](MULTIVARIATE_RATIONAL_UPPER_BANDS_ALL_DEGREES.md)
uses bounded integer multipliers to align separated bands of original
upper exponents. A [bounded-row cluster](MULTIVARIATE_BOUNDED_UPPER_CLUSTER_GAP_SPREAD.md)
can have arbitrary ratios between its gaps, and an
[empty interval between gap scales](MULTIVARIATE_GAP_LACUNA_UNBOUNDED_MULTIPLICITIES.md)
allows unbounded multiplicity within each remaining scale.
The [dense-gap theorem](MULTIVARIATE_DENSE_GAP_MARKS_ALL_DEGREES.md)
instead uses strict convexity of the number of copies of a marked factor
at each coefficient. It applies when one factor can occur everywhere.
The [tilted-band version](MULTIVARIATE_TILTED_DENSE_GAP_MARKS_ALL_DEGREES.md)
allows mixtures of bands with upper offsets proportional to the gaps
and no comparison between their orientation densities.
The [separated lower-endpoint theorem](MULTIVARIATE_SEPARATED_LOWER_ENDPOINT_ALL_DEGREES.md)
uses an alternative first-block mark with an orientation-dependent
retained budget, even for a complete spectrum of gap sizes.
The [simultaneous first-block theorem](MULTIVARIATE_SIMULTANEOUS_FIRST_BLOCK_ANGULAR_SEPARATION.md)
uses a common prime factor at all low orders. It yields binomial
reductions and angular separation for irregular gap distributions in
aligned bands, including rationally scaled bands and offsets bounded
by a multiple of each row's gap. The
[simultaneous midpoint theorem](MULTIVARIATE_SIMULTANEOUS_MIDPOINT_MARKS_ALL_DEGREES.md)
uses an alternative shared mark for the actual even and odd Newton
polygons. Its general comparison controls possible common-factor degree
from the nominal mark-count slopes.
The [universal two-row extension](MULTIVARIATE_UNIVERSAL_OPPOSED_PADIC_ANGLE_MARKING.md)
uses two opposed gaps to separate root quotients modulo a single
prime across all degrees below the absolute cutoff. Every
ordinary core then admits arbitrarily weak two-row extensions
with disjoint consecutive arguments in every degree, at any
prescribed signed tail gap. Equal gaps preserve every endpoint
center exactly. The
[geometric modular rescaling criterion](MULTIVARIATE_GEOMETRIC_MODULAR_RATIONAL_TWIST_STABILITY.md)
gives fully shared families with unrestricted dimension and
curvature that are coprime in every degree under all positive
rational rescalings. For terminal \(p+1\) with \(p\) an odd prime, the
[mixed unit-gap Newton-width criterion](MULTIVARIATE_TERMINAL_SIX_UNIT_GAP_NEWTON_WIDTH.md)
certifies coprimality on additional profiles and gives row changes
and additions that preserve its valuation test. The
[native-prime endpoint test](MULTIVARIATE_NATIVE_PRIME_CHORD_ANGULAR_CERTIFICATE.md)
checks four upper factorial-block valuations with every original
row fixed. Its all-degree witnesses have density one in boxes
with fixed signed gaps, also within Cartesian progressions.
The [primitive logarithmic criterion](MULTIVARIATE_PRIMITIVE_PADIC_LOG_ANGULAR_CERTIFICATE.md)
computes the exact fractional coordinates of the upper root
quotient logarithms, including terms beyond the leading interior
height. A cyclic subfield constraint gives a separate criterion
using only upper and lower tame ramification degrees.
The [global Galois orbit criterion](MULTIVARIATE_GALOIS_ORBIT_ANGULAR_SEPARATION.md)
excludes every adjacent root-angle coincidence when the upper
Galois group is transitive on unordered root pairs, or when the
upper degree is an odd prime and the upper moment is irreducible.
A delayed factorial mark gives an explicit point-stabilizer cycle,
and a native edge of width one less than the upper degree suffices
for physical coprimality. This can certify original profiles whose
endpoint chord is not primitive. The
[prime-power cycle theorem](MULTIVARIATE_PRIME_POWER_CYCLE_ANGULAR_GATE.md)
uses a full Galois cycle in prime-power upper degree. It also
works with a tame primitive full edge at a small prime,
without an interior valuation condition.
The [quadratic Newton trace sieve](MULTIVARIATE_QUADRATIC_NEWTON_TRACE_SIEVE.md)
constrains a possible degree-two common factor using the shared
root valuations and the exact angular boxes. It applies to
multi-edge polygons, including a two-index unit cluster.
The [primitive quadratic norm gate](MULTIVARIATE_PRIMITIVE_QUADRATIC_NORM_GATE.md)
uses the same local valuations with any supplied right-half-plane
angle cap. Its norm equation, exact discriminant index, and
ideal-class test apply to any rational parent pair.
The [joint quadratic budget](MULTIVARIATE_QUADRATIC_JOINT_PARITY_DISCRIMINANT_BUDGET.md)
makes the two parents share their parity and critical-point
conditions at each prime. It bounds the angle numerator and
discriminant complement alongside the spread denominator and
uses one slope allocation for all three. Finite local compatibility
at all primes of a candidate's arithmetic support can still be
prescribed. The [cofactor Smith theorem](MULTIVARIATE_ISOLATED_COFACTOR_SYLVESTER_BUDGET.md)
charges cofactor valuation overlap through the nonzero maximal
Sylvester minors. It also shows that the full Smith form of
rational parent arrays alone cannot bound an isolated quadratic.
The [field discriminant test](MULTIVARIATE_COMMON_FACTOR_FIELD_DISCRIMINANT.md)
filters the possible common factors of general profiles by their
ramification and signed field discriminants. It allows arbitrary
root moduli and intersects the discriminant possibilities with
the actual stepped residual branches. The
[tame field refinement](MULTIVARIATE_TAME_FIELD_DISCRIMINANT_SQUARECLASSES.md)
also determines the admissible tame local completions from a
shared stepped branch and filters field discriminant unit
squareclasses, including at unramified places.
The [CM-field screen](MULTIVARIATE_CM_COMMON_FACTOR_GATE.md)
uses the common root angles and uniform valuations to restrict a
CM common factor with two private sides to four cyclotomic
fields. It allows distinct absolute values for the conjugates.
For two private units the
[degree-growth obstruction](MULTIVARIATE_PRIVATE_UNIT_COMMON_FACTOR_DEGREE_GROWTH.md)
uses the absolute radial wall and a bound on algebraic units to
force the cyclotomic alternative at sufficiently large terminal
relative to the common-factor degree, without assuming the
field is CM.

The [paired Gamma annulus](MULTIVARIATE_PAIRED_GAMMA_ROOT_ANNULUS.md)
also gives explicit coprimality tests in the physical coefficient
ratio. Exactly mirrored row pairs contribute no radial width,
so the comparison can certify asymmetric profiles with an
arbitrary mirrored part. The
[linked Gauss matching criteria](MULTIVARIATE_LINKED_GAUSS_MATCHING_SEPARATION.md)
use both slopes and offsets, for differing shared ratios and
arbitrary private matchings. They also transfer angular separation
from a mirrored reference to profiles with unbalanced ratios.
The [signed degree-imbalance theorem](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
uses the linked radial centers and an integrality bound on shared
curvature to cover every ordinary profile with
\(|\Delta|\ge146J-73\), with arbitrary orientations and private
rows. Here \(J\) is the absolute angular cutoff.

At each terminal, coprimality over \(\mathbb Q\) gives the assertion
over any characteristic-zero field by the rational Bézout identity.
Both pure monomial coefficients are nonzero at every order, so
the univariate reduction in \(t\) covers all projective coefficient
choices. An all-degree proof still requires coverage of the terminals
from fourteen through \(N_*-1\), or an argument that applies directly
at every terminal.

## The discrete row-merger route

For proportional rows the [joint density comparison](MULTIVARIATE_JOINT_FACTORIAL_CURVATURE.md#inflation-and-identification-along-a-primitive-slope)
orders the effect of integer inflation and identification on all
coefficient curvatures. For nonproportional shared rows the
[exact merger sign law](MULTIVARIATE_ROW_MERGER_CURVATURE_SIGN.md)
gives an explicit flattening threshold and a decrease of the leading
stationary curvature.
The [ordered opposite-pair theorem](MULTIVARIATE_ORDERED_OPPOSITE_PAIR_RADIAL_CLOSURE.md)
compares the exact signed center to only the odd coefficient
perturbation around a mirrored core. The resulting all-degree
criterion includes arbitrarily many same-gap opposite pairs
with individually ordered midpoints and very small signed bias.
The [ordered gap-mass theorem](MULTIVARIATE_ORDERED_GAP_MASS_RADIAL_CLOSURE.md)
uses integer gap units to compare different scales and normalized
midpoints. It covers a primitive row inflated on one side
against the corresponding number of uninflated reversed copies
on the other, and also higher-order moment matching of crossing
scale distributions. General row identification and the case of
canceling center contributions remain open.

There is a discrete transfer for arithmetic certificates. The
[initial-form merger](SFC_PADIC_ROW_MERGER_INITIAL_FORMS.md)
combines two ordinary exponent rows by separating their base-\(p\)
digits. One diagonal change of the monomial coefficients preserves
all weighted moment initial forms through any prescribed order,
including for arbitrary support size. Iterating the merger reduces
these finite-window certificates to a single exponent row. The
transfer proves coprimality when the initial forms exclude common
residue zeros and leaves the lifting of a residue collision open.
The [higher-precision inflation theorem](SFC_PADIC_ROW_INFLATION_HIGH_PRECISION.md)
preserves linked factorial units modulo \(p^3\) under arbitrary
rowwise \(p\)-power inflation for primes \(p\ge5\), with increasing
precision at successive layers. The [modular row operations theorem](MULTIVARIATE_MODULAR_ROW_OPERATIONS_ALL_DEGREES.md)
fixes one prime for any all-degree regular core. It permits independent
baseline translations, arbitrarily many prescribed neutral residue
blocks, and subsequent row inflation and digit packing at that prime.
The [complete shift-block deletion theorem](MULTIVARIATE_COMPLETE_SHIFT_BLOCK_DELETION.md)
uses literal boundary products to transfer core resultants
modulo a composite modulus at higher precision. It includes
cases where a first residue test has a common factor, and
allows independent lifts and repetitions giving arbitrary
positive diffuse strength under the complete-block hypotheses.
The [complementary-row comparison](MULTIVARIATE_COMPLEMENTARY_ROW_EULER_REDUCTION.md)
reduces a complementary pair to a weighted Euler operator at every
required degree when the common cap has sufficiently large
\(p\)-adic valuation. It gives an all-degree ordinary family for
linear disjoint core monomials and any two fixed positive weights.
The [complementary Euler extension](MULTIVARIATE_ARBITRARY_CORE_COMPLEMENTARY_EULER_EXTENSION.md)
uses several distinct positive base scales to separate possible
multiple roots of an arbitrary exact binary core. It gives an
ordinary all-degree extension along every sufficiently divisible
total at any chosen prime, with arbitrarily small added curvature
and no change in the total-degree difference.
For exactly mirrored cores, the [mirrored complementary congruence](MULTIVARIATE_MIRRORED_COMPLEMENTARY_ENDPOINT_CONGRUENCE.md)
gives a direct condition on existing complementary groups.
Their Euler factors can annihilate every mixed coefficient
below the absolute cutoff modulo one odd prime, leaving
coprime endpoint binomials and retaining exact mirror symmetry.

The [signed diffuse projective theorem](SFC_SIGNED_DIFFUSE_PROJECTIVE_REGULARIZATION.md)
uses both the diagonal variance and the signed center drift at each
order. The latter orders the pure coefficient costs by reciprocal
moment degree. This gives simultaneous projective regularity for all
order sets of the support size in the algebraic-strength limit, and
ordinary fully shared, total-degree-balanced realizations on any
prescribed finite collection of windows. The arbitrary-core binary
result also supplies a common all-order diffuse tail for any finite
family of chosen ordinary row mergers. The
[integer-twist stability theorem](MULTIVARIATE_ARBITRARY_CORE_TWISTED_TAIL_STABILITY.md)
identifies the exact limiting rescaling \(((n+1)/n)^\Delta\)
for a fixed ordinary core and a small-curvature shared tail of signed
gap \(\Delta\). A rational resultant at each of finitely many possible
twists gives an effective stability threshold outside the exceptional
set. For comparable both-private cores any exception is a negative
twist, with a further prime-power exclusion for unit-gap rows.
The [projective first-block criterion](SFC_SIMULTANEOUS_FIRST_BLOCK_PROJECTIVE_WINDOWS.md)
uses separating high exponent patterns on existing rows and equal
pure factorial coefficients. One prime gives pure-power reductions
and a Vandermonde comparison for every support-size window up to
a prescribed maximum order.
The [triangular row construction](SFC_TRIANGULAR_INTEGER_ROW_TAGGING.md)
forces projective regularity at any prescribed finite family of
windows for arbitrary support size by adding \(s-1\) integer tag
rows. Explicit congruences in the smaller exponents give ordinary
realizations at arbitrarily small added curvature, and one fixed
compensating row preserves total-degree differences.

The [sharp private-degree theorem](SFC_PHYSICAL_ASYMMETRIC_PRIVATE_DEGREES.md)
gives all-order physical binary transfer when both untouched private degree
sums are positive and their total is at least four. The proof applies to
any positive base with finite mixed moments. An exact common root at orders
four and five with private degrees one and two shows that the full
positive-base statement cannot be weakened to that pair of degrees.
For every ordinary profile with private variables on both sides,
the [asymmetric angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
prove simplicity in every degree and exclude common roots in the
closed left half-plane. If there is a common root, it must have the
same label in the two right-half boxes. The
[endpoint Euclidean descent](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md)
excludes every label at relative scales in the exact first-coefficient
radius interval. It replaces the angular restriction of the
[linked first-box comparison](MULTIVARIATE_BOTH_PRIVATE_FIRST_BOX_LINKED_PHASE.md)
by an interlaced-ratio condition preserved under endpoint cancellation.
The [binomial-weighting obstruction](MULTIVARIATE_BINOMIAL_RATIO_DESCENT_OBSTRUCTION.md)
shows that stripped coefficient ratios and positive Mellin quotients
alone do not extend this descent to profiles with missing private
factorials. The [ordinary remainder example](MULTIVARIATE_ORDINARY_BINOMIAL_REMAINDER_OBSTRUCTION.md)
exhibits failure of the proposed inherited ratio inequality on two
balanced shared rows. The
[single-factorial weighting example](MULTIVARIATE_SINGLE_FACTORIAL_WEIGHT_OBSTRUCTION.md)
shows that the bare ratio invariant fails even when one
reciprocal factorial is inserted. The
[residual curvature theorem](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_CURVATURE_INTERLACING.md)
and its [residual-curvature-four refinement](MULTIVARIATE_BOTH_PRIVATE_RESIDUAL_CURVATURE_FOUR.md)
prove strict all-degree angular interlacing at total absolute
curvature \(S\ge6\), regardless of the remaining orientations.
The [seven-halves residual theorem](MULTIVARIATE_BOTH_PRIVATE_RESIDUAL_CURVATURE_SEVEN_HALVES.md)
lowers the sufficient threshold to \(S\ge11/2\), with strict
interlacing at terminals from fourteen onward.
The [conditional monotonicity comparison](MULTIVARIATE_BOTH_PRIVATE_COMPARABLE_SHARED_INTERLACING.md)
proves every-degree coprimality at arbitrary curvature if the
exponents on the shared rows are coordinatewise comparable.
The [balanced tail stability theorem](MULTIVARIATE_COMPARABLE_CORE_BALANCED_TAIL_STABILITY.md)
extends this family to sufficiently small shared tails of either
orientation with zero total gap. The
[effective aligned-tail theorem](MULTIVARIATE_COMPARABLE_CORE_EFFECTIVE_TAIL_STABILITY.md)
gives a tolerance depending only on the absolute cutoff and allows
any nonnegative total tail gap in the core orientation.
The [shifted midpoint reduction](MULTIVARIATE_PRIVATE_UNIT_SHIFTED_TANH_RESOLVENT_REDUCTION.md)
gives exact even and odd equations with arbitrary shared orientations.
The odd equation uniquely specifies a radius at each angle in a right
half-plane box. A folded order comparison would suffice at coincident
candidate radii. The centered convex-order comparison holds unconditionally,
and the stronger translated comparison can fail at arbitrary radii.
The [three-atom construction](SFC_PRIVATE_UNIT_ANGLE_BOX_SHARPNESS.md)
shows that every overlap of two angle boxes can support a common
root for a fixed positive exchangeable base. When the private exponent
total is at least two on each side, the
[reciprocal binomial coefficient cap](MULTIVARIATE_BOTH_PRIVATE_EXPONENT_SUM_TWO.md)
instead gives all-degree angular separation with arbitrary
remaining rows. Its phase proof applies to independent
asymmetric coefficient arrays satisfying the same cap.
The [physical mirrored circle theorem](SFC_PHYSICAL_MIRRORED_ROOT_CIRCLE.md)
gives simple circular roots at every order under arbitrary independent
correlations of mirrored row pairs, with no lower bound on their gaps
or reserved degrees. It preserves disjoint residual radius intervals.
For one private mirrored unit pair it also gives consecutive interlacing
throughout the noncollapsed physical path.
The same note's hypergeometric Riccati comparison gives exact strict
curvature flattening for any correlated pair with opposite equal gaps,
even with unmatched smaller exponents. At common fractional arguments
it orders the normalized correlation derivatives across integer degrees
and identifies their limiting symbol.
The [matched Euler comparison](SFC_MATCHED_EULER_REAL_INTERLACING.md)
gives all-order strict interlacing for a fixed Gamma tilt and a class
of positive shape mixtures. The [merger edge comparison](SFC_PHYSICAL_MIRROR_MERGER_EDGE.md)
gives a uniform angular confinement bound for arbitrary mirrored
blocks and a strictly interlacing scaled endpoint limit for a
single ordinary block of any gap.
The [simultaneous merger limit](SFC_PHYSICAL_MULTIBLOCK_MERGER_MATCHINGS.md)
identifies the multiblock limit as an expected weighted matching
polynomial of a colored permutation graph, or equivalently a finite
additive convolution. It has simple real roots. The adjacent-order
comparison of those matching polynomials remains open.
The [one-private angle boxes](MULTIVARIATE_ONE_PRIVATE_MIRRORED_ANGLE_BOXES.md)
give one of the two interlacing inequalities with any exchangeable
positive shared base. The [root monotonicity theorem](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_ROOT_MONOTONICITY.md)
orders fixed-degree root angles under changes of mirrored rows.
The [linked phase theorem](MULTIVARIATE_ONE_PRIVATE_MIRRORED_LINKED_PHASE_INTERLACING.md)
uses the adjacent folded likelihood-ratio comparison to prove strict
interlacing with a distinguished private unit pair and any additional
mirrored rows. It covers every degree without a curvature restriction,
including Gamma shape changes and common Gaussian augmentation.
The [exact-core stability theorem](MULTIVARIATE_PRIVATE_UNIT_EXACT_CORE_ANGLE_STABILITY.md)
gives a uniform polynomial angular gap for the same reference with
physically matched Gaussian clocks. It permits a tail of ordinary
rows with arbitrary orientations under a universal relative-gap
condition or a signed shape bound. The reference keeps the entire
mirrored core exact, and the tolerance depends polynomially on
the absolute cutoff.
The [large-gap pair theorem](MULTIVARIATE_MIRRORED_HARD_PAIR_INTERLACING.md)
gives strict all-degree interlacing for fully shared mirrored profiles
if one pair has gap at least seven times its smaller exponent.
The [circular angular transport theorem](MULTIVARIATE_MIRRORED_UNIT_CORE_ANGULAR_TRANSPORT.md)
bounds fixed-degree root displacements under circular coefficient
convolution. An explicit Gegenbauer gap then retains all-degree
angular separation near an exact shared unit-gap core.

A proof for arbitrary support size cannot rely on zero-freeness
throughout the physical Gaussian interpolation. An
[ordinary three-monomial example](SFC_PHYSICAL_ROW_MERGE_INTERIOR_COLLISION.md)
has regular endpoints for the window of orders two through four
and a common projective zero at a physical squared correlation
strictly between \(21/25\) and \(17/20\). At the exact coefficient
direction the [arithmetic obstruction](SFC_ORDINARY_SEARCH_NEAR_PHYSICAL_COLLISION.md)
excludes every nonzero rational homogeneous form through degree
eighteen, including an ordinary quadratic moment. The discrete
endpoint transfer and the general binary interpolation remain open.
The [private-unit removal obstruction](SFC_PRIVATE_UNIT_REMOVAL_FOLDED_ORDER_OBSTRUCTION.md)
shows that both folded likelihood-ratio inequalities still allow
an interior adjacent collision under physical correlation or
private Gamma shape inflation. This holds for a fixed exchangeable
core and for one mirrored Gamma pair with real noninteger exponents.
The construction uses no ordinary integer-gap counterexample.


## Exact arithmetic along common translations

The [single-row fiber theorem](MULTIVARIATE_SINGLE_ROW_ANGULAR_EXCEPTIONAL_FIBERS.md)
holds with every other row fixed arbitrarily. The negative-baseline
specialization gives a nonzero root-ratio resultant in the chosen
baseline, leading to a uniform all-degree angular exception bound
depending only on the chosen gap and the absolute cutoff. Varying
both rows of a reflected pair also gives cofinite all-degree
coprimality for an arbitrary core and preserves exact mirroring
when the core is mirrored.

For any fixed terminal \(N\), add an integer \(T\) to both entries
of every unequal row. If \(D\) is the sum of the unsigned row gaps,
the [exact translated resultant degree](MULTIVARIATE_TRANSLATED_RESULTANT_DEGREE.md)
is \(N(N-1)D\) when the total signed gap is nonzero, and
\(N(N-1)(D-1/2)\) in the balanced case. This bounds the number
of exceptional translations at the given terminal for every
profile. The [all-order exceptional polynomial](MULTIVARIATE_ALL_ORDER_EXCEPTIONAL_LOCUS.md)
combines the finite resultant product with its exact highest
homogeneous part in the row bases. It bounds the exceptions
on every positive integer ray and shows that, at fixed
dimension \(r\), the exceptional fraction among all exponent
pairs in \(\{0,\ldots,H\}^{2r}\) is \(O_r(H^{-(r+1)/2})\),
with an \(O(H^{-1})\) bound uniform in \(r\).
Its partial-translation leading face is a scaled resultant
on the smaller support. For at least \(2(J-1)\) unequal rows, the
[codimension-two theorem](MULTIVARIATE_EXCEPTIONAL_CODIMENSION_TWO.md)
places the all-order exceptional locus inside a proper algebraic
set of codimension at least two. From any initial base vector,
some positive integral rowwise direction then avoids every
exception at positive integer translation. The [next homogeneous layer](MULTIVARIATE_RESULTANT_SECOND_HOMOGENEOUS_LAYER.md)
also has positive coefficients on its support in the unbalanced case. The [divisor-prime criterion](MULTIVARIATE_TRANSLATED_DIVISOR_PRIME.md)
gives a coprimality certificate using suitable prime divisors of a
linear function of the translation, including primes below
the factorial arguments. The
[quantitative prime-screen construction](MULTIVARIATE_QUANTITATIVE_PRIME_SCREEN_OBSTRUCTION.md)
shows why even a growing collection of small-prime tests can
leave degree two admissible. These statements leave the
integer points of the all-order exceptional locus to be excluded.
For a linear monomial with a shared variable, the
[intermediate Gamma polynomial](MULTIVARIATE_LINEAR_MONOMIAL_EULER_BRIDGE.md)
reduces consecutive coprimality to a single squarefreeness question.
Its differential equation bounds the maximum length of a vanishing
window. An exact synthetic law for the private factor gives a repeated
quadratic while preserving the linked Gamma operators.
For [prime-power terminals](MULTIVARIATE_LINEAR_MONOMIAL_PRIME_POWER_WINDOWS.md),
valuations and the residue polynomial give exact physical separation
in additional ordinary cases with private prime-power exponents.
The [Gamma numerical-range proof](MULTIVARIATE_SHARED_UNIT_GAP_GAMMA_NUMERICAL_RANGE.md)
gives physical coprimality for adjacent Gamma moment polynomials with
positive lower shape and any nonnegative integer shape increment.
It is a separate all-degree proof for one unequal ordinary unit-gap row.
The [native Newton response test](MULTIVARIATE_NATIVE_NEWTON_ANGULAR_RESPONSE.md)
uses tame edge data to certify first-order splitting of any core ray
contact under a sufficiently weak ordinary shared tail.
The [multiple-root response theorem](MULTIVARIATE_WEAK_TAIL_MULTIPLE_ROOT_RESPONSE.md)
applies the same type of weak perturbation at every root
multiplicity. Scaled Hermite zeros separate all noncentral
branches. At an odd central pair, an isolated close conjugate
pair at a prime certifies angular splitting, and the linked
signed center and cubic skew give a further physical test at
integer twists. The estimates apply with arbitrary weak-tail
row counts and gaps, including controlled two-scale approximation
of the retained body. They give additional finite-height
certificates for branches of the original bounded-gap property.
The [stationary central jet theorem](MULTIVARIATE_STATIONARY_CENTRAL_WEAK_TAIL_JETS.md)
gives the uniform second angular response for pairs of odd
central branches with identical first angular velocity. A
compact description of all weak-tail jet statistics, with
explicit finite outer errors and realizability at fixed signed
gap, gives exact tests for strict margins, even at a stationary
physical twist.
The [single-row progression theorem](MULTIVARIATE_SINGLE_ROW_ALL_ORDER_ARITHMETIC_PROGRESSION.md)
bounds the number of exceptional choices on every ordinary fibre obtained
by varying one smaller exponent. The bound depends on its gap and the
absolute eventual cutoff alone. It constructs an all-order regular
progression inside any prescribed residue class. The
[one-movable-row angular progression](MULTIVARIATE_SINGLE_MOVABLE_ROW_PADIC_PROGRESSION.md)
obtains the stronger separation under independent positive
rescalings by comparing the first nonintegral term in the upper
root quotients with the lower tame ramification lattices.
The [affine-ray theorem](MULTIVARIATE_AFFINE_BASE_RAY_PADIC_ANGLES.md)
applies this comparison to several existing movable rows when
the rightmost upper factor has the stated orientation imbalance.
It also produces Cartesian congruence classes with disjoint
arguments and bounds exceptional fibers in the original exponent grid.
The [all-direction affine-ray theorem](MULTIVARIATE_AFFINE_BASE_RAY_ALL_DIRECTIONS.md)
gives all-degree physical coprimality away from finitely many
parameters on every affine ray. A resultant at the first upper
factor gives an explicit integer divisor test at each bounded
terminal, including balanced and nonprimitive cohorts.
