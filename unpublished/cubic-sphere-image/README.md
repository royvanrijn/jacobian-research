# Cubic sphere moments and the two-pair Image problem

[Read the paper](cubic-sphere-image.pdf) · [LaTeX source](main.tex) · [Source bundle](source.zip) · [Verification record](verification.json)

Roy van Rijn. Local manuscript; no public deposit or submission is recorded. External mathematical review and formal verification have not occurred.

The paper classifies all pure phases in the full sixteen-dimensional filtered space of degree at most three on `t²+2xy=1`. It proves the sharp mixed-moment cutoff `m > 3 deg(g)`, transfers it to bihomogeneous bidegree `(3,3)` in two contraction pairs over every characteristic-zero field, and reproduces Long's quartic witness to establish that balanced degree four is the first failure.

The proof has been narrowed to a proper-map volume-period obstruction, explicit boundary resolutions, a four-moment interval lemma, and an eight-moment calculation. The cusp resolution removes the large eleven-parameter algebraic certificate from the essential proof. Every essential new argument and finite algebraic input is printed in the article. The short accompanying program is a reproducibility aid; it does not establish the geometric or all-order conclusions by testing finitely many exponents.

## Build and replay

From this directory, with a standard TeX Live installation:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=cubic-sphere-image main.tex
```

The source files `main.tex`, `geometry.tex` and `finite.tex` compile without the research workspace. The bibliography is included in `main.tex`.

For the exact algebra replay, use Python 3.12 and SymPy 1.14.0:

```sh
python -m pip install -r requirements.txt
python verify.py
```

The verifier regenerates the moment polynomials, checks every displayed elimination, and verifies the small resultant obstruction modulo 19. It uses rational and finite-field arithmetic, with no external polynomial data. Optimized Python (`python -O`) is rejected so that assertions cannot be silently disabled.

## Review and scope

The edition hashes, compilation, arithmetic replay, page inspection, primary-source locators and actual PDF-only review are recorded in `verification.json`. An isolated AI blind-reader check is internal review and does not represent external human review.

Before submission: obtain specialist review of the proper direct-image/Fourier period argument and the exhaustive boundary resolution, and independent replay of the finite calculation. No unrestricted Image conjecture, general mixture of bidegrees, or result with additional contraction pairs is claimed.