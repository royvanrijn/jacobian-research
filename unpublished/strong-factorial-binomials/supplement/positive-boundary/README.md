# Positive-exponent boundary

This is the electronic proof appendix cited in the section “The large
relative gaps” of *The Strong Factorial Conjecture for Univariate
Binomials*. It contains written proofs, exact arithmetic programs, and
the input and result records needed for those arguments. The files are
preserved byte for byte in [the proof snapshot](proof-snapshot.zip).
The [archive checksum](ARCHIVE_SHA256SUMS) verifies the ZIP, and
the [member checksums](SHA256SUMS) verify the extracted files.
The archive contains the sources and replay inputs listed below. Cross-links to earlier work outside these proof routes refer to historical context in the corresponding canonical research sources.
Paths in the following list are relative to the extracted
`research/factorial/image-rigidity/` directory.

## Lower exponent one

- **Prime-power gaps.** The theorem labelled
  `thm:prime-power-lower-one` in `ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex`
  gives squarefreeness for a positive \(p\)-integral shape and gap
  \(p^s\). Here the shape is the integer \(N\). The all-gap algebraic
  presentation used in that proof is in `gap_three_hurwitz_structure.tex`
  and `HIGHER_GAP_ADJOINT_RANK.tex`.
- **Gap six.** The contour estimates in
  `GAP_SIX_RESULTANT_STRICT_POSITIVITY.md` prove squarefreeness at
  the physical shape \(N\) for \(N\ge 10{,}000\). The article gives
  the derivative identities converting this conclusion to
  coprimality at the actual consecutive moment orders.
- **Gaps at least ten.** The effective contour estimates are in
  `ORDINARY_LOWER_ONE_STRONG_FACTORIAL.md`, Sections 1–6, together
  with the deformation estimates in
  `FACTORIAL_UNIFORM_LARGE_GAP_SADDLE_SIMPLICITY.md` and
  `FACTORIAL_GAP_TEN_SADDLE_EXTENSION.md`.

## Lower exponent at least two and gap at least three times it

- **Contour estimates.** `ORDINARY_LOWER_TWO_STRONG_FACTORIAL.md`,
  Sections 2–6, contains the base contour proof and its Bernstein
  certificates for \(a=2\). The scaled bounds and the uniform
  estimate for the order-\(N-1\) amplitude are in
  `RELATIVE_GAP_THREE_LOWER_WEIGHTED_SADDLE_TAIL.md`, Sections 2–6.
  These proofs give the effective onset
  \(N\ge\max(4000,128(a+d))\).
- **Finite domains.** The article uses the ratio bound to restrict
  the finite work to three smaller domains. The [finite boundary verifier](../verify_positive_boundary.cpp)
  excludes all 326,146 supports remaining after the ratio test by
  prime-block and Newton-polygon checks at the two actual orders. Earlier,
  larger finite runs described in the archived written proofs are
  unnecessary for this article.

## Exact replay

Run the following commands from the directory containing the ZIP.
The analytical comparison programs use standard-library Python 3.
Checking the gap-six contour constants through the archived program
also requires python-flint 0.9.0 at import time.

```sh
sha256sum --check ARCHIVE_SHA256SUMS
snapshot_dir=$(mktemp -d)
unzip -q proof-snapshot.zip -d "$snapshot_dir"
cd "$snapshot_dir"
sha256sum --check SHA256SUMS

src=research/factorial/image-rigidity
python3 "$src/verify_factorial_practical_uniform_saddle_tail.py" --repo-root "$PWD" --output /tmp/sfb-positive-one-tail.json
python3 "$src/verify_lower_two_weighted_saddle_tail.py" --repo-root "$PWD" --output /tmp/sfb-positive-two-tail.json
python3 "$src/verify_relative_gap_three_weighted_tail.py" --repo-root "$PWD" --output /tmp/sfb-positive-relative-tail.json

PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$src" python3 -c 'from verify_gap_six_resultant_strict import check_contour_constants; check_contour_constants()'
```

The analytical verifiers check rational inequalities and, where
applicable, the exact Bernstein coefficients. The written contour
proofs give the stated unbounded ranges. The article derives the
remaining finite domains independently and identifies their checker.
