# Prime Knots with Augmentation-Rank Gaps at Unbounded Bridge Number

**Published preprint · [Hexagon v1](https://hexagonmath.org/2611.00595v1) · 19 November 2026 · [Deposited files](published/)**

Explicit braid families have complex augmentation rank at most \(b-1\) and
meridional rank and bridge number equal to \(b\) for every \(b\ge5\).
Each family contains infinitely many hyperbolic prime knot types for
\(b=5,6\) and for every sufficiently large \(b\). The thresholds for
hyperbolicity are existential.

[PDF](prime_knot_augmentation_gaps.pdf) · [LaTeX](main.tex) ·
[Computational supplement](prime_knot_certificates.zip) ·
[Metadata](metadata.json) · [Verification](verification.json)

Build with pdfTeX from TeX Live 2023 and latexmk:

~~~sh
SOURCE_DATE_EPOCH=1794960000 latexmk -pdf -jobname=prime_knot_augmentation_gaps -interaction=nonstopmode -halt-on-error main.tex
~~~

The [certificate replay guide](certificate-replay.md) describes the exact
inputs and checks. Supporting notes give the
[uniform algebra](../../research/meridional-rank/representations/UNIFORM_PROTECTED_AUGMENTATION_GAPS.md)
and the [diagram transfer](../../research/meridional-rank/representations/UNBOUNDED_PRIME_AUGMENTATION_GAP_FAMILIES.md).
