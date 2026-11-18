# Eisenstein restricted-prime sums: research and verification

This is the canonical home for the branch's supporting research notes,
verification programs, receipts and supplied inputs. The
[unpublished manuscript](../../../unpublished/eisenstein-restricted-prime-sums/README.md)
contains only paper sources, its PDF, paper records and its generated source
archive. The relocation preserves every moved file's bytes and mathematical
scope; it is not a new calculation or an assurance upgrade.

## Current mathematics

- [Power-shrinking Eisenstein windows, with an explicit onset](../EISENSTEIN_POWER_SHRINKING_WINDOWS.md).
- [Sifted angular distribution and the direct finite-polynomial balance proof](../EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md).
- [Binary sums with prescribed prime factors](../BINARY_SUMS_WITH_PRESCRIBED_PRIME_FACTORS.md).
- [Obstructions to very narrow central windows](../PRESCRIBED_PRIME_CENTRAL_GAPS.md).
- [Semigroup-scale density one](../PRESCRIBED_PRIME_SEMIGROUP_DENSITY_ONE.md).
- [Binary variance and the remaining thin-window energy problem](../PRESCRIBED_PRIME_BINARY_VARIANCE.md).
- [Simultaneous Gaussian--Eisenstein case study](../SIMULTANEOUS_GAUSSIAN_EISENSTEIN_SUMS.md).
- [Simultaneous-norm pointwise upper bounds and the fixed-prime sieve criterion](../SIMULTANEOUS_NORM_POINTWISE_SIEVE.md).

The manuscript proves both exact identities for the positive integers `A`
whose prime divisors are `1 mod 3`, including `1`: `A+A` is precisely the
positive `2 mod 6` class, and `A+(3·A)` is precisely the positive `4 mod 6`
class. Its primitive-norm and higher-sumset consequences retain their original
scope. The current angular note proves eventual balance. The separate
power-window proof gives every fixed difference exponent `theta > 2/3`, in
both ordinary and mixed cases, with an explicit onset formula and at least
`n^(theta-1) R_o(n)/20` coordinate witnesses per orientation. For example,
difference at most `n^(3/4)` holds for every admissible `n >= 2^(2^40)`.
The exact triangle expansion and central-annulus tail yield this range;
the positive Jackson--Parseval construction also verifies the weaker `theta > 6/7`
route. The explicit onsets leave an enormous gap above the finite certificate,
so the complete balanced exception lists are still unresolved. These are research notes; the manuscript
has not been expanded to include the new power-window theorem.

For the simultaneous modulus-12 restriction, the new pointwise note proves
upper bounds and a natural-scale large-prime tail. Its lower-bound and
fixed-prime distribution inputs remain open, as does every-target binary
coverage. The general note records the separate pointwise, angular, variance
and proof-audit campaigns.

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

The paper's [verification and review record](../../../unpublished/eisenstein-restricted-prime-sums/verification.json)
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
