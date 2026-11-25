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
terminal orders, starting at twelve by the
[complete terminal-eleven theorem](MULTIVARIATE_TERMINAL_ELEVEN_ALL_PROFILES.md).
Each such order still includes arbitrary dimension and exponents. The
[fixed-terminal compactification](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md)
describes the full closure of the normalized coefficient pairs and
the relative center. A certified separation of that family gives
a finite algebraic certificate for the terminal. The closure includes
a signed center parameter as well as a diffuse variance. A separation
test that varies the two moment arrays independently admits spurious
collisions already at six, as shown by the
[coefficient-link calculation](MULTIVARIATE_TERMINAL_SIX_COEFFICIENT_LINKS.md).
A proof at each remaining terminal must retain the links between the
two orders or use the exact arithmetic of the factorial coefficients.
The [joint curvature measure](MULTIVARIATE_JOINT_FACTORIAL_CURVATURE.md)
uses an exact factorial density for each primitive row curve. It
couples all orders and determines the relative center. Relaxing only
these densities to arbitrary positive measures on the same curves
already gives a collision at orders two and three.
The [uniform curvature closure](MULTIVARIATE_JOINT_CURVATURE_PROFILE_CLOSURE.md)
identifies the actual diffuse part as a diagonal density \(r\,dr\).
Its error is uniform in the two curvature indices, and the relative
centers have one additional signed diffuse parameter.
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
complete. At terminal twelve the
[fixed rays](MULTIVARIATE_TERMINAL_TWELVE_HIGH_CURVATURE.md)
cover \(S\ge8\), and the
[private-degree theorem](MULTIVARIATE_TERMINAL_TWELVE_PRIVATE_DEGREES.md)
covers every profile with total private degree at least two.

Across all terminals, the [explicit high-curvature theorem](MULTIVARIATE_EFFECTIVE_ABSOLUTE_CURVATURE_BOUND.md)
and the [uniform weak-curvature theorem](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md)
confine a possible counterexample to \(\eta_*/2<S<64\), where
\(\eta_*>0\) is a computable absolute rational. The
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

A prescribed finite set of primes and precisions alone cannot give
a universal exclusion. The
[simultaneous prime-power construction](MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md#a-simultaneous-obstruction-at-finitely-many-primes)
produces fully shared profiles with both normalized moments congruent
to their binomial polynomials at all the prescribed moduli, within
any open positive curvature interval. In particular, the ordinary
Newton and residual screens at those moduli allow both candidate
common-factor degrees at terminal six.
The [marked-row Eisenstein criterion](MULTIVARIATE_MARKED_ROW_VARIABLE_PRIME.md)
gives examples in every degree using a prime selected from the
exponent data. It can also extend any fixed profile by two balanced
shared rows with arbitrarily small added curvature, preserving
private degrees and total-degree difference. The
[opposite marked-row test](MULTIVARIATE_OPPOSITE_MARKED_ROW_PRIME.md)
uses two rows with overlapping below-prime windows. It
covers balanced, arbitrarily small-curvature marks even
when their common maximum prevents isolating one row at
the chosen prime. For terminal \(p+1\) with \(p\) an odd prime, the
[mixed unit-gap Newton-width criterion](MULTIVARIATE_TERMINAL_SIX_UNIT_GAP_NEWTON_WIDTH.md)
certifies coprimality on additional profiles and gives row changes
and additions that preserve its valuation test.

The [paired Gamma annulus](MULTIVARIATE_PAIRED_GAMMA_ROOT_ANNULUS.md)
also gives explicit coprimality tests in the physical coefficient
ratio. Exactly mirrored row pairs contribute no radial width,
so the comparison can certify asymmetric profiles with an
arbitrary mirrored part.

At each terminal, coprimality over \(\mathbb Q\) gives the assertion
over any characteristic-zero field by the rational Bézout identity.
Both pure monomial coefficients are nonzero at every order, so
the univariate reduction in \(t\) covers all projective coefficient
choices. An all-degree proof still requires coverage of the terminals
from eleven through \(N_*-1\), or an argument that applies directly
at every terminal.

## The discrete row-merger route

For proportional rows the [joint density comparison](MULTIVARIATE_JOINT_FACTORIAL_CURVATURE.md#inflation-and-identification-along-a-primitive-slope)
orders the effect of integer inflation and identification on all
coefficient curvatures. For nonproportional shared rows the
[exact merger sign law](MULTIVARIATE_ROW_MERGER_CURVATURE_SIGN.md)
gives an explicit flattening threshold and a decrease of the leading
stationary curvature. These comparisons leave the corresponding
root transfer open.

There is a discrete transfer for arithmetic certificates. The
[initial-form merger](SFC_PADIC_ROW_MERGER_INITIAL_FORMS.md)
combines two ordinary exponent rows by separating their base-\(p\)
digits. One diagonal change of the monomial coefficients preserves
all weighted moment initial forms through any prescribed order,
including for arbitrary support size. Iterating the merger reduces
these finite-window certificates to a single exponent row. The
transfer proves coprimality when the initial forms exclude common
residue zeros and leaves the lifting of a residue collision open.

There is an [all-order physical binary transfer](SFC_PHYSICAL_BINARY_RESIDUAL_PRIVATE.md)
when untouched independent rows leave at least two units of private
exponent on each side. The proof applies to every positive correlated
base with finite mixed moments, by Hölder's inequality and separation
of the root arguments.

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
