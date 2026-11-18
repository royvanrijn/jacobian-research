# Prime Knots with Augmentation-Rank Gaps at Unbounded Bridge Number

[PDF](prime_knot_augmentation_gaps.pdf) · [LaTeX](main.tex) ·
[Computational supplement](prime_knot_certificates.zip) ·
[Metadata](metadata.json) · [Verification](verification.json)

The paper constructs explicit families with complex augmentation rank at
most $b-1$ and meridional rank and bridge number both equal to $b$ for
every $b\ge5$. At bridge numbers five and six and every sufficiently
large bridge number, each family contains infinitely many hyperbolic prime
knot types. The paper gives the uniform rank-one obstruction, the alternating
quotients, and the oriented diagram induction that reduces the geometry to
one fixed three-cusped exterior and two direct small-bridge certificates.

The bridge-four family has a separate geometry and exact rank calculation.
This paper starts at bridge five, where one braid formula and one port
induction apply. The paper's thresholds are non-effective, and its result is
a strict inequality for augmentation rank.

Build in this directory with pdfTeX from TeX Live 2023 and latexmk:

```sh
SOURCE_DATE_EPOCH=1794960000 latexmk -pdf -jobname=prime_knot_augmentation_gaps -interaction=nonstopmode -halt-on-error main.tex
```

All LaTeX sources needed for the build are in this directory. The
[replay instructions](certificate-replay.md) also travel in the computational
supplement, with copies of the exact inputs and checkers. Regenerate the
archive with `python3 bundle-certificates.py`. The canonical research records
give the [uniform algebra](../../research/meridional-rank/representations/UNIFORM_PROTECTED_AUGMENTATION_GAPS.md)
and the [diagram transfer](../../research/meridional-rank/representations/UNBOUNDED_PRIME_AUGMENTATION_GAP_FAMILIES.md).

Independent specialist review should examine the signed crossing reductions,
the peripheral marking in the two fillings, and the uniform Reidemeister
induction with its marked insertion band.
