# Three-pair Image counterexample

This directory contains the short standalone preprint
*A Four-Term Counterexample to the Special Image Conjecture in Three
Pairs*.

Compile it with:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The final repository PDF is copied to:

```text
output/pdf/three-pair-image-counterexample.pdf
```

Verify the mathematical certificate from the repository root with:

```bash
python3 scripts/verify_three_pair_image_mathieu_counterexample.py
```

From the standalone Zenodo bundle, run:

```bash
python3 verify_three_pair_image_mathieu_counterexample.py
```

The proof in the paper establishes the identities for every positive
power. The checker is a finite exact replay and regression test.
