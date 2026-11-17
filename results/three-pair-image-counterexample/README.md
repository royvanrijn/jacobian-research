# A Four-Term Counterexample to the Special Image Conjecture in Three Pairs

An explicit four-term polynomial disproves the Special Image Conjecture in three contraction pairs.

**Published preprint · 2026-07-28 · [10.5281/zenodo.21634058](https://doi.org/10.5281/zenodo.21634058)**

[PDF](three-pair-image-counterexample.pdf) · [LaTeX](main.tex) · [Zenodo](https://zenodo.org/records/21634058) · [Submission metadata](zenodo-metadata.json) · [Publication metadata](metadata.json)

The PDF and deposited bundle files are the exact Zenodo edition. [Checksums](deposited-files.json) record their provenance.

## Verification

Written all-power identities and an accompanying finite exact regression checker. The deposited frontier 2 <= r_SIC <= 3 is historical; later local work establishes failure with two pairs. The finite checker alone is not an all-order theorem proof.

- [Evidence 1](verify_three_pair_image_mathieu_counterexample.py)
- [Evidence 2](../../research/foundations/extended-geometry/TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md)

The catalogue corrects the record title's grammatical typo; the downloaded files retain their original bytes.

Build the deposited source from this directory with `latexmk -pdf main.tex`. Formal packages and broader research notes remain under `research/`.
