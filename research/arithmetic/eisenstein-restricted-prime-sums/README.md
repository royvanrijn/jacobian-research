# Eisenstein restricted-prime sums: research and verification

This is the canonical home for the branch's supporting research notes,
verification programs, receipts and supplied inputs. The
[exact primitive-norm companion](../../../results/eisenstein-restricted-prime-sums/README.md)
contains the four-discriminant article, PDF, review records and source
packages. It is a published preprint on [Hexagon v1](https://hexagonmath.org/2611.00673v1),
with exact public files in its [deposited edition](../../../results/eisenstein-restricted-prime-sums/published/).
The [pointwise companion](../../../results/pointwise-binary-splitting-semigroups/README.md)
is also published on [Hexagon v1](https://hexagonmath.org/2611.00672v1), with its own
[deposited edition](../../../results/pointwise-binary-splitting-semigroups/published/).
The companion citations do not make the exact classifications depend on the pointwise argument.

## Research focus

The primary target is the exact ordinary balanced classification: every positive
`n = 2 mod 6`, except `8` and `44`, should be a sum of two members of `A`
with both summands in `[n/3, 2n/3]`. The retained finite certificate covers
`n < 10^8`; the eventual balance theorem has not yet been connected to that
range. The new
[fixed-interval reduction](../EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md#exact-balanced-exceptions-prime-families-and-reduced-targets)
covers every `8p` and `44p` with `p = 1 mod 3` prime, and all-inert targets
above `2^60`. The residual task is now to cover the all-inert admissible cores
in `[10^8, 2^60)`, by sharper estimates or an effective finite bridge. This is
equivalent to the full conjecture; the residual interval remains unresolved.

For fixed abelian subgroup restrictions on both summands, the
[general splitting theorem](../SPLITTING_SEMIGROUP_LOCAL_GLOBAL.md) now gives
pointwise leading asymptotics under the explicitly retained zero-free and
weighted-dilation inputs. This includes `C12+C12` on every sufficiently large
target `N = 2 mod 12`. The older
[asymmetric subgroup route](../ASYMMETRIC_SUBGROUP_SPLITTING.md) and its sieve
reductions retain their independent scopes. Arbitrary prime prescriptions and
dependent coupled systems are not covered merely by the subgroup theorem.

## Current mathematics

- [Exact balanced exceptions: prime-family closure and the remaining inert cores](../EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md#exact-balanced-exceptions-prime-families-and-reduced-targets).
- [Pointwise semigroup-scale counts for Eisenstein sums](../EISENSTEIN_POINTWISE_SEMIGROUP_SCALE.md).
- [Pointwise Eisenstein counts in moving summand windows](../EISENSTEIN_MOVING_WINDOWS.md).
- [Power-shrinking Eisenstein windows, with an explicit onset](../EISENSTEIN_POWER_SHRINKING_WINDOWS.md).
- [Pointwise logarithmic lower bounds for distinct Eisenstein value pairs](../EISENSTEIN_POINTWISE_VALUE_COUNTS.md).
- [Every-target Eisenstein sums with bounded numbers of prime factors](../EISENSTEIN_BOUNDED_PRIME_FACTORS.md).
- [Sifted angular distribution and the direct finite-polynomial balance proof](../EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md).
- [Binary sums with prescribed prime factors](../BINARY_SUMS_WITH_PRESCRIBED_PRIME_FACTORS.md).
- [Obstructions to very narrow central windows](../PRESCRIBED_PRIME_CENTRAL_GAPS.md).
- [Semigroup-scale density one](../PRESCRIBED_PRIME_SEMIGROUP_DENSITY_ONE.md).
- [Binary variance and the remaining thin-window energy problem](../PRESCRIBED_PRIME_BINARY_VARIANCE.md).
- [Quantitative full-count variance in power windows](../PRESCRIBED_PRIME_POWER_WINDOW_VARIANCE.md).
- [Eisenstein Hecke transfer audit](../EISENSTEIN_HECKE_TRANSFER_AUDIT.md).
- [Power-saving exceptions for the modulus-12 semigroup lower count](../PRESCRIBED_PRIME_POWER_SAVING_EXCEPTIONS.md).
- [Simultaneous Gaussian--Eisenstein case study](../SIMULTANEOUS_GAUSSIAN_EISENSTEIN_SUMS.md).
- [Every-target weighted asymmetric Gaussian--Eisenstein sums](../ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md).
- [Uniform small coefficients in the asymmetric convolution](../ASYMMETRIC_UNIFORM_COEFFICIENTS.md).
- [Every-target unweighted asymmetric counts at the semigroup scale](../ASYMMETRIC_SEMIGROUP_SCALE.md).
- [Asymmetric theory for every fixed modulus and quadratic character](../ASYMMETRIC_FIXED_MODULUS_THEORY.md).
- [Abelian subgroup splitting with nonnested allowed sets](../ASYMMETRIC_SUBGROUP_SPLITTING.md).
- [Exact sign transfer and its periodic-replacement obstruction](../ASYMMETRIC_SIGN_TRANSFER.md).
- [Cofactor drift and two-large-prime bounds in the unweighted mixed sieve](../COUPLED_GAUSSIAN_EISENSTEIN_SIEVE.md).
- [Well-factorability obstructions and the exact large-prime Rosser switch](../ASYMMETRIC_WELL_FACTORABLE_SWITCH.md).
- [The positive single-index operator information barrier](../EISENSTEIN_SINGLE_INDEX_OPERATOR_BARRIER.md).
- [Averaging auxiliary coefficients and a positive distribution level](../ASYMMETRIC_AUXILIARY_AVERAGE.md).
- [Half-dimensional criterion for two modulus-12 summands](../ASYMMETRIC_SEMILINEAR_SIEVE.md).
- [Critical-level switching and the two-large-prime bound](../ASYMMETRIC_CRITICAL_SWITCHING.md).
- [Fixed periodic factors in divisor boxes](../FIXED_PERIODIC_DIVISOR_BOXES.md).
- [Fixed-prime marginals and raw dilation laws](../SIMULTANEOUS_NORM_DILATION_LAWS.md).
- [Simultaneous-norm pointwise upper bounds and the fixed-prime sieve criterion](../SIMULTANEOUS_NORM_POINTWISE_SIEVE.md).
- [Obstructions to absolute character-divisor truncation](../SIMULTANEOUS_NORM_DIVISOR_OBSTRUCTIONS.md).
- [Weighted simultaneous-norm means in arithmetic progressions](../SIMULTANEOUS_NORM_WEIGHTED_AP_MEANS.md).
- [Raw and fully sifted simultaneous-norm correlations with power exceptions](../SIMULTANEOUS_NORM_POWER_CORRELATIONS.md).

The manuscript proves both exact identities for the positive integers `A`
whose prime divisors are `1 mod 3`, including `1`: `A+A` is precisely the
positive `2 mod 6` class, and `A+(3·A)` is precisely the positive `4 mod 6`
class. The companion also proves the exact primitive norm-value sumsets at
discriminants `-7`, `-8` and `-11`, independently of the pointwise companion.
Its higher-sumset consequences retain their original scope. The current angular
note proves eventual balance. The separate
power-window proof gives every fixed difference exponent `2/3 < theta < 1`, in
both ordinary and mixed cases, with an explicit onset formula and at least
`n^(theta-1) R_o(n)/20` coordinate witnesses per orientation. For example,
difference at most `n^(3/4)` holds for every admissible `n >= 2^(2^40)`.
Distinct value pairs now have matching upper and lower counts of order
`n^theta J(n)/log n`, where `J(n)` is the inert target-prime product,
for every sufficiently large admissible target in the same exponent range.
This also holds for fixed central widths, including balance and unrestricted
pairs. The cofactor proof restores a logarithmic factor to the rough-pair
sieve and treats every target valuation. The earlier higher-moment lower
bound remains available as a separate method. The moving-window theorem
now proves matching counts for arbitrary fixed macroscopic intervals and,
more generally, for `H <= M = min(x,n-x)` with
`H^3 >= n^(1+epsilon) M`. A designated small summand may have size
`n^theta` for every `theta > 1/2`; the nearly equal threshold stays
`theta > 2/3`. These new counts have existential onsets. The joint coset
theorem now compresses the finite-level remainder at every target valuation,
including prescribed mixed orientations and exact prime-square tests. Its
aggregate bound couples the cofactor and sieve budgets without the earlier
residue-space size powers.
A further lower sieve gives `n^theta/(log n)^2` pairs
whose two permitted parts each have a bounded number of prime factors, all
exceeding a fixed positive power of the target. The exact triangle expansion and central-annulus tail yield
the power-window range;
the positive Jackson--Parseval construction also verifies the weaker `theta > 6/7`
route. The fixed-interval argument separately closes the prime-multiplier
families and proves balance on all-inert targets above `2^60`. It reduces the
ordinary exact-exception conjecture to the all-inert interval `[10^8, 2^60)`,
which remains unresolved. The [separate window paper](../../../results/eisenstein-pointwise-counts/)
contains the moving-window value-count results; the exact companion retains
its classification and certificate scope.

The asymmetric theorem now gives a weighted pointwise asymptotic on EVERY
sufficiently large target `N = 2 mod 6`, with the first summand in `C12`
and the second in `A3`. Its positive target-dependent singular series gives
one-sided existence. A separately proved small-coefficient extension and
cofactor construction now give matching unweighted counts of order
`N J_E(N)/(log N)^(5/4)`, uniformly in target valuations. This has no
balance restriction or explicit cutoff, and does not give two `C12`
summands. A half-dimensional sieve criterion and a further critical
switching bound isolate a square-root-level remainder sufficient for the
two-sided problem, up to logarithmic losses. The switching upper bound
holds unconditionally on every target, including all target valuations.
The local model, parity reduction and exact kernel tail are proved;
the critical global remainder is open. The full auxiliary-modulus average
now controls coefficients through `N^xi` before choosing its smaller error
saving, and gives an L1 sieve level at every fixed exponent below `xi`.
The general theory extends asymmetric counts to every fixed modulus and
chosen quadratic splitting condition. A separate clause gives every-target
counts at the full two-sided logarithmic scale while permitting boundedly
many large excluded prime factors in the second summand; eliminating them
remains open.

The subgroup extension permits nonnested first-side and quadratic splitting
sets. Its finite residue density determines admissibility, and its target
product uses primes outside both allowed sets. The resulting C4+C3 count
does not force two C12 summands. The sign-transfer note retains the signed
main model and identifies the full-order cost of a fixed-periodic
replacement of the remaining coefficient. The critical remainder stays open.

For the simultaneous modulus-12 restriction, the new pointwise note proves
upper bounds and a natural-scale large-prime tail. A separate elementary
four-character hyperbola argument now supplies weighted arithmetic-progression
means, including polynomially growing moduli. The new full weighted
correlation theorem proves the natural lower bound and actual bad-prime
marginals outside one power-sized exceptional set, including arbitrary target
valuations and full moving balance intervals. Its fully sifted weighted
proportion has an exact target-dependent Euler product. The all-target versions of these particular weighted power-window estimates
remain open. Unrestricted `C12+C12` coverage and its leading binary count are
supplied separately by the flagship theorem under its retained recent inputs.
The full semigroup variance is proved for every fixed modulus and
`5/6 < theta < 1`, including the original `.99` energy target. Every-target
power-window counts and power-saving full-count variance remain open. The separate
modulus-12 lower-count theorem now gives `H/(log X)^(3/2)` pairs outside
`O(X^(1-10^(-6)))` compatible targets for `H=X^(99/100)`. The general note records the
separate pointwise, angular, variance and proof-audit campaigns.

## Files and provenance

- The root Python/Java programs are the retained binary verifiers, exact
  constant checks and spectral/ramified diagnostics.
- [Research notes](research_notes/) retain the supplied spectral argument,
  bounded literature audit and earlier unfinished weighted route. The linked
  current angular note gives the completed qualitative balance arguments.
- [Receipts](receipts/) retain imported evidence;
  [local replay receipts](receipts/local-replay/) record the previous internal
  executions, with their original scopes and historical commands.
- [Supplied extended-package receipts](receipts/supplied-extension/) accompany
  the byte-preserved [original extended input archive](supplied_extension.zip).
- `supplied/compute_mixed_certificate.py` is the received implementation;
  `compute_mixed_certificate.py` is the separate canonical witness checker.

The paper's [verification and review record](../../../results/eisenstein-restricted-prime-sums/verification.json)
retains the existing calculation provenance and manuscript review. Its current
artifact paths point here. Historical receipts and source hashes were not
rewritten. A paper source archive may contain distribution copies of the code
and evidence needed for reproduction; those copies are generated from here.

## Replay from this directory

The printed standalone verifier needs only Python 3. Run the commands below
from `research/arithmetic/eisenstein-restricted-prime-sums/`:

```sh
mkdir -p ../../../.scratch/eisenstein-replay
python3 verify_binary_certificates.py > ../../../.scratch/eisenstein-replay/binary_certificates_receipt.json
```

It checks all `16,666,667` ordinary and `16,666,666` mixed targets below
`10^8`, with the recorded zero exceptions. The least small summands in those
finite certificates are at most `211` and `199`. These are retained results;
no numerical replay was performed as part of the relocation.

The separate NumPy generators recheck every witness:

```sh
python3 compute_finite_certificate.py --limit 99999999 --output ../../../.scratch/eisenstein-replay/ordinary
python3 compute_mixed_certificate.py --limit 99999999 --output ../../../.scratch/eisenstein-replay/mixed
python3 verify_mixed_constants.py > ../../../.scratch/eisenstein-replay/mixed_constants_receipt.json
python3 check_ramified_counts.py --limit 6000 --output ../../../.scratch/eisenstein-replay/ramified_counts_receipt.json
```

The recorded environment used Python 3.12.3 and NumPy 1.26.4. The ramified
program needs the explicit `6000` limit to reproduce the retained receipt;
its default is `3000`. It checks 444 reduced even targets, exact quotient
calculations and a wrong-sign control. These checks support the written
all-index argument and are not a substitute for it.

Some diagnostics write outputs next to their source. Copy them into ignored
scratch storage before executing, so the imported receipts remain untouched:

```sh
cp compute_elementary_constants.py compute_mixed_odd_checks.py compute_mixed_wheel_checks.py mixed_spectral_check.py ../../../.scratch/eisenstein-replay/
python3 ../../../.scratch/eisenstein-replay/compute_elementary_constants.py
python3 ../../../.scratch/eisenstein-replay/compute_mixed_odd_checks.py
python3 ../../../.scratch/eisenstein-replay/compute_mixed_wheel_checks.py
python3 ../../../.scratch/eisenstein-replay/mixed_spectral_check.py
python3 supplied/compute_mixed_certificate.py --limit 99999999 --output ../../../.scratch/eisenstein-replay/supplied-mixed
python3 compute_mixed_ramified_checks.py --limit 3000 --output ../../../.scratch/eisenstein-replay/supplied-ramified.json
```

The degree-20 weight checker needs SymPy. Its recorded internal replay used
SymPy 1.14.0 in the repository's existing virtual environment, with assertions
enabled. It also writes its receipt beside the copied source:

```sh
cp research_notes/literature_balanced_polynomial.py ../../../.scratch/eisenstein-replay/
../../../.venv/bin/python ../../../.scratch/eisenstein-replay/literature_balanced_polynomial.py
```

The canonical printed verifier uses explicit exceptions and has recorded
invalid-input checks under `python3 -O`. Some supplied diagnostics use
`assert`; run those normally. The polynomial certificate verifies finite
rational inequalities only. The harmonic-theta argument linking that weight
to eventual balance is in the current angular note.

The unweighted mixed sieve has half-dimensional local density, but its bare
Haar main fails at every positive power level through cofactor density
drift. Its unconditional two-large-prime bound leaves a sifted lower bound
open. This is distinct from the correct tau-weighted R2 model. Pure
well-factorable lower weights cannot reproduce the critical Rosser sieve;
the exact large-prime switch leaves a global remainder to estimate.

For Eisenstein windows, a continuous single-index operator can preserve
nonnegative functions and all finite prime-to-index coset means while
avoiding a central band on the two-thirds scale. It is not positive
semidefinite, an arithmetic shell or a multiplicative Hecke family. No
arithmetic exponent or every-target modulus-12 conclusion follows from
these obstruction results.

## Exact finite controls for the new analytic arguments

From this directory, the portable diagnostics print their results by
default and write JSON only when `--output` is supplied:

```sh
python3 verify_joint_coset_distribution.py
python3 verify_balanced_exception_reduction.py
python3 verify_asymmetric_sums.py
python3 verify_asymmetric_uniform.py
python3 verify_asymmetric_cofactors.py
python3 verify_asymmetric_switching.py
python3 verify_asymmetric_auxiliary_average.py
python3 verify_fixed_modulus_asymmetric.py
python3 verify_subgroup_asymmetric.py
python3 verify_asymmetric_sign_transfer.py
python3 verify_coupled_gaussian_eisenstein_sieve.py
python3 verify_asymmetric_well_factorable_switch.py
```

The [balanced-exception controls](balanced_exception_reduction_checks.json)
certify 39 finite coefficient, normalization, interval-envelope and onset
checks. They enumerate 616 inert primes through `10000`, without scanning
any new target range or verifying the remaining inert-core interval. The
written proof supplies the uniform moment input, infinite coefficient tail
and descent. The checker does not certify those analytic statements.

The other retained receipts are [joint coset controls](joint_coset_distribution_checks.json)
and [asymmetric arithmetic controls](asymmetric_sums_checks.json), together
with [uniform-family controls](asymmetric_uniform_checks.json) and
[cofactor controls](asymmetric_cofactor_checks.json), and
[switching controls](asymmetric_switching_checks.json),
[auxiliary-average controls](asymmetric_auxiliary_average_checks.json), and
[fixed-modulus controls](fixed_modulus_asymmetric_checks.json),
[subgroup controls](subgroup_asymmetric_checks.json), and
[sign-transfer controls](asymmetric_sign_transfer_checks.json),
[coupled-sieve controls](coupled_gaussian_eisenstein_sieve_checks.json), and
[well-factorability controls](asymmetric_well_factorable_switch_checks.json). All commands reject
Python `-O`, because their exact finite assertions must run.
The first checks the asymmetric singular recurrences, 240 scalar identities,
and every invertible matrix modulo 49, with wrong-coefficient controls. The
second checks local integrals, gcd and hyperbola identities, including
controls that reject omitted corrections and strata. These support the
written identities. The uniform-family checker covers 89 exact local
cases, 6,120 coefficient reductions, and 1,500 shifted-kernel identities,
including the failure of a mod-12 support assertion for shifted second
kernels. The cofactor checker covers 376 exact local compensations,
1,331 lower-polynomial sign cases and 5,104 injective fibres, including
higher target valuations and an invalid lower-polynomial control.
The switching checker verifies 507 exact local integrals, 507 corrected
Henriot inequalities and 153 Euler identities, with missing-content and
incorrect-expectation controls. These finite checks do not certify an
unrestricted asymptotic, prove the critical distribution hypothesis, or
replace the imported analytic estimates. The auxiliary-average checker adds
60 exact family cases, 40 local prime-power cases, 2,397 CRT checks,
222 kernel/model identities and 120 cutoff cases, including 30 with a cutoff
larger than the divisor level. The fixed-modulus checker includes even and
odd characters, prime 3, fixed-place densities, 21,500 shifted second-kernel
cases, 8,028 strict divisor-pairing identities, 22,612 restricted-cofactor
fibres and nine negative controls. Their receipts identify all finite ranges
and current input bytes.
The subgroup checker adds 13 nested and nonnested configurations, including
odd index, 144 local integrals at 3, 5 and 7, and recoverable directives on
the recorded finite ranges. The sign-transfer checker covers finite
convolutions, exact signs, Rosser support and least-squares projection
identities. These diagnostics do not establish asymptotic means or the
critical remainder. The coupled-sieve receipt contains 11,534 finite exact
controls for coefficient identities, residues, parity, content and local
factors. The well-factorability receipt contains 2,145 finite controls,
including 1,426 exact large-prime switch cases. The operator barrier is
analytic and has no finite checker. These receipts support only their
stated finite ranges, not the unrestricted analytic conclusions.

## Java and assurance limits

The archived mixed Java source and receipt match the recorded Python result,
but no Java runtime was available for a fresh Java execution. With an existing
Java 21 runtime, its documented invocation is:

```sh
java -Xmx512m VerifyMixedEisenstein.java 100000000 ../../../.scratch/eisenstein-replay/mixed-java.json
```

Run Java from this research directory: both verifiers locate their own source
for hashing relative to the working directory. Alternatively copy the source
to scratch and run with that scratch directory as the working directory.
The Python bound is inclusive and this Java bound is exclusive.

Internal execution, imported receipts, source inspection and internal AI
review remain distinct from external specialist review, external replay and
formal verification. No such assurance changed when these files moved.
