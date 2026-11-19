# Exact Sumsets of Primitive Quadratic Norms

This companion manuscript treats the four discriminants **-3, -7, -8 and
-11**. The current [PDF](paper.pdf), edition hashes and scoped review records
are maintained with the manuscript artifacts below.

The presentation places the elementary identity P_-3 = A disjoint union 3A
beside the opening theorem, groups the Eisenstein argument in one main section,
and summarizes the finite scans, interval certificates, infinite estimates and
descent in an opening coverage table. The theorem scopes, certificate programs
and retained receipts are unchanged by this reorganization.

[LaTeX](main.tex) · [Metadata](metadata.json) ·
[Verification record](verification.json) · [Source archive](source.zip)

## Mathematical scope

For positive primitive numerical values P_d of the displayed form, the paper
states the following exact binary sumsets. Every row has N >= 2. Coordinates
may vanish or have either sign; gcd(x,y) = 1 is required, and different
coordinate representations do not create different numerical summands.

| Discriminant d | Form Q_d(x,y) | Exact condition for N in P_d + P_d |
|---|---|---|
| -3 | x² - xy + y² | N even and N not congruent to 0 or 12 modulo 18 |
| -7 | x² + xy + 2y² | 49 does not divide N, excluding 7, 77 and 203 |
| -8 | x² + 2y² | 16 does not divide N |
| -11 | x² + xy + 3y² | N even and 121 does not divide N |

The stronger split-prime core identities, the existing Eisenstein higher-sumset
corollaries and the averaging criterion retain their stated scopes.
No further discriminant or quadratic order is added. The balanced Eisenstein
exception classification remains an open question beyond its declared finite
scan; the observed failures 8 and 44 do not establish its unrestricted version.
No numerical leading correlation asymptotic is claimed here. The separate
pointwise-correlation manuscript is not a proof dependency of these exact
classifications: this paper uses its own classical-input arguments and finite
certificates.

## Build

From this directory, or the root of the extracted source archive:

```sh
SOURCE_DATE_EPOCH=1795046400 FORCE_SOURCE_DATE=1 latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
```

All included TeX and verification sources are present locally. The standalone
source archive needs no repository checkout. Its certificate files are
snapshots of the canonical arithmetic checkers; canonical development stays
under research/arithmetic in the repository. `build_bundle.py` creates the
archive and checks included paths. To check a fresh extraction against the
frozen PDF, run `python3 build_bundle.py --verify-rebuild`.

## Replay the exact certificates

Use Python3 with assertions enabled, without -O, -OO or PYTHONOPTIMIZE. The
main programs use only the standard library. From the source root:

```sh
python3 certificates/verify_binary_certificates.py > certificates/binary_certificates_receipt.json
python3 certificates/verify_primitive_norm_minus7.py --output certificates/primitive_norm_minus7_checks.json
python3 certificates/verify_primitive_norm_minus8.py --output certificates/primitive_norm_minus8_checks.json
python3 certificates/verify_primitive_norm_minus11.py --output certificates/primitive_norm_minus11_checks.json
python3 certificates/verify_primitive_norm_coordinate_membership.py --output certificates/primitive_norm_coordinate_membership_checks.json
python3 certificates/verify_primitive_norm_quaternary.py --output certificates/primitive_norm_quaternary_checks.json
python3 certificates/compute_elementary_constants.py
python3 certificates/verify_mixed_constants.py
python3 certificates/check_ramified_counts.py --output certificates/ramified_counts_receipt.json
```

The coordinate authenticator uses a separate coordinate/gcd algorithm, with a
standalone path adapter. It loads the -8/-11 support sieves solely for whole-array
comparison and authenticates every finite membership flag. The -7 checker
independently authenticates its full membership arrays too. These are separate
internal implementations; no external human replay is claimed. The larger
arrays require a few hundred megabytes of memory. Retained NumPy and Java
Eisenstein implementations provide additional replay options; their actual
execution scopes remain in the preserved verification record.

## Checks already recorded

`verification.json` records prior internal replays on 19 November 2026. Their
retained receipts have the following scopes; this documentation pass has not
rerun the mathematical computations.

| Existing check | Recorded finite scope |
|---|---|
| Eisenstein exact scan | All 33,333,333 admissible ordinary/mixed targets below 10^8 |
| -7 exact checker | Euler and local quotient checks, 257 fine and 17 higher interval certificates, complete coverage through 2^24 and complete listed exceptions |
| -8 exact checker | Euler checks, 107 interval certificates, complete C+C and C+2C coverage through 2^27 and corruption controls |
| -11 exact checker | Sturm coefficients, Euler/lattice bounds, 35 interval inequalities, complete C+C coverage through 2^25 and corruption controls |
| Coordinate authentication | Whole-array coordinate/gcd comparison through 2^27 at -8 and 2^25 at -11; -7 authentication is in its exact checker |
| Auxiliary arithmetic | Dyadic theta/Sturm and Euler certificates, exact Eisenstein constants, finite ramified/multiplicativity diagnostics |
| Runtime controls | Four assertion-based exact verifiers reject -O, -OO and PYTHONOPTIMIZE=1 before writing output: 12 controls |

These checks certify their declared finite domains. The written spectral
estimates, rational interval coverage and multiplicative descent supply the
extension to all targets. Passing a checker does not check those written
arguments or the imported theorems.

The existing build record reports resolved references and citations, no
overfull boxes, and a byte-identical PDF from a fresh source extraction. Its
layout record specifies which page renders were inspected. An isolated
AI reader reviewed the preceding 36-page PDF; its supplied context, findings and
rechecks are recorded separately. All of this is internal review. External
specialist review, external human certificate replay, formal verification and
historical priority are not established.

The earlier scope-freeze pass checked receipt hashes and preservation of the
mathematical files, and validated its standalone source archive. Its outcome
is recorded under `review_preparation` in `verification.json`. The current
editorial edition requires its own build and layout record.

## Questions for outside review

Please assess the written argument and the declared supplement together:

1. Check the spectral import in Section 2 against the cited Voight,
   Jacquet--Langlands and Deligne statements: all-integral-ideal normalization,
   global norm-character removal, maximal level at the unramified index primes,
   and self-adjointness for the correct invariant mass.
2. Check the arithmetic local adaptations in the Eisenstein proof (Section 3) and the discriminant -7, -8 and -11
   proofs (Sections 4--6): the odd
   and ramified filters, the nonmaximal dyadic level at -8, target-dividing
   branches, and all ideal classes with inverse-unit mass at -11. In particular,
   verify that each local quotient and shell count matches the spectral setup.
3. Check that every rational interval certificate meets its finite scan and
   infinite bound, and that scaling/descent covers every claimed target with
   exactly the stated congruence obstructions and exceptions. Appendix A
   explains the finite reductions and byte conventions.
4. Replay the exact certificates independently, including full membership-array
   authentication and negative controls. Compare outputs and declared domains,
   rather than relying only on the retained hashes.
5. Audit attribution and historical priority against the closest restricted-prime
   and primitive-norm sumset literature. The adaptation from Howe--Huffman--
   Kohlmeyer--Kruer is attributed; a complete priority audit remains open.

The bundle is prepared for these reviews. No outreach, submission or public
deposit is part of this preparation.
