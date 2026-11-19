# Exact Sumsets of Primitive Quadratic Norms

This is the current unpublished companion manuscript. It evolves the earlier
Eisenstein manuscript in this directory and keeps its complete binary proof.
The paper determines the exact numerical primitive norm-value self-sumsets at
discriminants -3, -7, -8 and -11; the additional -7 exceptions are 7,77,203.
It also proves the stronger split-prime core identities used in the descents.

[PDF](paper.pdf) · [LaTeX](main.tex) · [Metadata](metadata.json) ·
[Verification](verification.json) · [Source archive](source.zip)

## Build

From this directory, or the root of the extracted source archive:

```sh
latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
```

All included TeX and verification sources are present locally. The standalone
source archive needs no repository checkout. Its certificate files are
snapshots of the canonical arithmetic checkers; canonical development stays
under research/arithmetic in the repository. `build_bundle.py` creates the
archive and checks included paths.

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

The independent coordinate authenticator has a standalone path adapter; its
membership and gcd algorithms are unchanged. It loads the -8/-11 support
sieves solely for whole-array comparison. It authenticates every finite
membership flag, rather than merely comparing target masks or hashes.
The -7 checker independently authenticates its full membership arrays too.
The larger arrays require a few hundred megabytes of memory. Retained NumPy
and Java Eisenstein implementations provide additional replay options; their
older execution scopes remain in the preserved review record.

All recipes above were replayed in this preparation. The exact rational
certificates, scans and negative controls passed. They certify their declared
finite domains; the written spectral estimates, interval coverage and descent
are what extend the result to all targets.

## Review boundaries

The current 36-page PDF compiles without unresolved references, citations or
overfull boxes. All pages were rendered and visually inspected. The existing
verification record preserves earlier Eisenstein reviews and records current
source and finite replay checks separately. This is internal review; there is
no claimed external specialist review, formal verification or established
historical priority. The PDF-only blind-reader status is recorded explicitly. A fresh extraction of
source.zip rebuilds a byte-identical PDF.

Before submission, obtain external specialist review of the Hecke and local
quotient arguments, independent human certificate replay and a targeted
priority audit. No numerical leading asymptotic is claimed here. The separate
flagship manuscript treats pointwise correlations; these exact identities use
their own classical-input proofs and complete finite certificates.
