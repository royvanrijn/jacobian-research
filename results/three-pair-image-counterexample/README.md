# A Four-Term Counterexample to the Special Image Conjecture in Three Pairs

[Current PDF](three-pair-image-counterexample.pdf) · [LaTeX](main.tex) · [Source bundle](source.zip) · [Verification](verification.json) · [Metadata](metadata.json)

The current repository edition gives the four-term three-pair construction and includes a complete proof of the later bidegree-(4,4) two-pair example. The one-pair theorem then identifies two as the minimum failing pair dimension. The proof explains the factorial lift and the endpoint-zero cancellation before drawing the dimension conclusion.

This is a local revision. The exact [published edition](published/) retains the deposited PDF, source, checker, submission metadata and [checksums](published/deposited-files.json) for [Zenodo 21634058](https://zenodo.org/records/21634058). Its statement that the failing dimension is two or three is historical.

Build from this directory with `latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=three-pair-image-counterexample main.tex`.

Run the dependency-free three-pair regression with `python3 verify_three_pair_image_mathieu_counterexample.py`. It writes the current local receipt and never touches the archived deposit. The two-pair checker and canonical proof remain in [research](../../research/foundations/extended-geometry/TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md).

The manuscript proves both identities for every positive exponent. The finite checks only test recorded ranges. No new formalization or independent external review is claimed. Independent mathematical review of the revised article remains a next step.

An isolated internal AI reader read the PDF and rechecked the clarified mixed-selection argument and final pagination. Its answers and limits are recorded in the verification record.
