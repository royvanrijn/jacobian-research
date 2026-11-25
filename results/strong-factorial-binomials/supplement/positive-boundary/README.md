# Positive-exponent boundary

The [contour proof](contour-proof.md) gives the three unbounded
estimates used in the article's appendix “The large relative gaps.”
It derives the stationary branches, deformed contour bounds and
terminal-order cutoffs in one place. The prime-power proof and the
fixed-shape to actual-order derivative identity are in the article.

| Range | Proof in this appendix | Exact replay |
| --- | --- | --- |
| $a=1,\ d=6,\ N\ge10000$ | [Gap six](contour-proof.md#gap-six) | [Local rational checker](verify_contour_constants.py) |
| $a=1,\ d\ge10,\ N\ge\max(5000,128(d+1))$ | [Lower exponent one and gaps at least ten](contour-proof.md#lower-exponent-one-and-gaps-at-least-ten) | Gap-ten contour and practical-onset comparisons in the snapshot |
| $a\ge2,\ d\ge3a,\ N\ge\max(4000,128(a+d))$ | [Lower exponent at least two and gap at least three times it](contour-proof.md#lower-exponent-at-least-two-and-gap-at-least-three-times-it) | Corner Bernstein coefficients, support-monotonicity guards, and weighted Gaussian comparisons in the snapshot |
| Orders below those onsets | The article's ratio bound leaves three finite domains | [Finite boundary verifier](../verify_positive_boundary.cpp) |

The finite verifier excludes all 326,146 supports remaining after
the uniform filter $d^2<21(a+d)$, by prime blocks and Newton polygons
at the actual two orders. The local checker independently counts
those supports and encloses the logarithms used to obtain the
cutoffs $d\le19$ and $d\le24$.
See the [supplement overview](../README.md) for the C++ replay.

## Exact replay

Run these commands from this directory with standard-library Python 3,
unzip, and sha256sum.

~~~sh
set -eu
python3 verify_contour_constants.py
sha256sum --check ARCHIVE_SHA256SUMS
snapshot_dir=$(mktemp -d)
unzip -q proof-snapshot.zip -d "$snapshot_dir"
(cd "$snapshot_dir" && sha256sum --check SHA256SUMS)

src="$snapshot_dir/research/factorial/image-rigidity"
python3 "$src/verify_factorial_gap_ten_saddle_extension.py" --repo-root "$snapshot_dir" --output "$snapshot_dir/gap-ten.json"
python3 "$src/verify_factorial_practical_uniform_saddle_tail.py" --repo-root "$snapshot_dir" --output "$snapshot_dir/lower-one-tail.json"
python3 "$src/verify_lower_two_weighted_saddle_tail.py" --repo-root "$snapshot_dir" --output "$snapshot_dir/lower-two-tail.json"
python3 "$src/verify_relative_gap_three_weighted_tail.py" --repo-root "$snapshot_dir" --output "$snapshot_dir/relative-gap-tail.json"
~~~

The written proofs establish the estimates throughout the unbounded
ranges. The scripts certify the rational inequalities and reconstruct
the Bernstein polynomials and their coefficients. The archived scripts
also report domain counts for earlier reductions. The C++ verifier
above performs the finite coprimality checks used in the article.

## Source record

The [checksummed snapshot](proof-snapshot.zip) retains the original
proof sources and exact replay inputs. Its [archive checksum](ARCHIVE_SHA256SUMS)
identifies the ZIP, and the [member checksums](SHA256SUMS) identify
its extracted files. The three contour arguments consolidate the
analytic portions of the following sources:

- The gap-six argument, from “Uniform simplicity for terminal order
  $N\ge10{,}000$” through “Cancellation on the overlap has a simple
  zero” in `GAP_SIX_RESULTANT_STRICT_POSITIVITY.md`.
- The practical lower-one tail, Sections 1–6 of
  `ORDINARY_LOWER_ONE_STRONG_FACTORIAL.md`, with the raw contour
  bounds from Sections 1–7 of
  `FACTORIAL_UNIFORM_LARGE_GAP_SADDLE_SIMPLICITY.md` and Sections 2–5
  of `FACTORIAL_GAP_TEN_SADDLE_EXTENSION.md`.
- The weighted tail, Sections 2–6 of both
  `RELATIVE_GAP_THREE_LOWER_WEIGHTED_SADDLE_TAIL.md` and
  `ORDINARY_LOWER_TWO_STRONG_FACTORIAL.md`.

Those files reside in the snapshot's
`research/factorial/image-rigidity/` directory. Older finite
prefixes, shape-resultant identities and external links in the
snapshot give provenance for the retained estimates.
