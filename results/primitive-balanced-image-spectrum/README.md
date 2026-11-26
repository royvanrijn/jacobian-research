# The primitive balanced-degree spectrum in the two-pair Image problem

**Published preprint · 26 November 2026 · [Hexagon 2611.01080v1](https://hexagonmath.org/2611.01080v1)**

The manuscript gives the exact primitive, absolutely nonpower failure
spectrum \(\{4,6,7,8,\ldots\}\) over characteristic-zero fields, and a
density-one ordinarity theorem for abelian varieties with rank-two CM
multiplication and at most two balanced embedding pairs.

[Deposited edition](published/) ·
[Locally reviewed PDF](primitive-balanced-image-spectrum.pdf) ·
[Original source ZIP](primitive-balanced-image-spectrum-hexagon-source.zip) ·
[Verification ZIP](primitive-balanced-image-spectrum-verification.zip) ·
[Edition metadata](metadata.json) ·
[Local checksums](SHA256SUMS.txt)

## Build and replay

Extract the original TeX ZIP and run this from its root with TeX Live
2023 and latexmk 4.83.

    SOURCE_DATE_EPOCH=1795564800 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=primitive-balanced-image-spectrum main.tex

All five TeX files match the deposited sources. The locally reviewed PDF
and the PDF compiled by Hexagon are kept as separate artifacts. The
[verification supplement](supplement/) supplies exact replay commands,
inputs, and the mathematical dependency inventory. The
[verification record](verification.json) identifies the reviewed PDF and
distinguishes finite replay, PDF-only internal review, and external
mathematical review. The [literature comparison](ordinary-reduction-literature.md)
records the precise arithmetic theorem overlap and scope.

[AI provenance](AI-PROVENANCE.md) · [Announcement](ANNOUNCEMENT.md) ·
[Research theorem and proof map](../../research/foundations/extended-geometry/PRIMITIVE_BALANCED_IMAGE_DEGREE_SPECTRUM.md)
