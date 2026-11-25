# The primitive balanced-degree spectrum in the two-pair Image problem

The manuscript gives the exact primitive, absolutely nonpower failure
spectrum \(\{4,6,7,8,\ldots\}\) over characteristic-zero fields, and a
density-one ordinarity theorem for abelian varieties with rank-two CM
multiplication and at most two balanced embedding pairs.

[Manuscript PDF](primitive-balanced-image-spectrum.pdf) ·
[Hexagon source ZIP](primitive-balanced-image-spectrum-hexagon-source.zip) ·
[Verification ZIP](primitive-balanced-image-spectrum-verification.zip) ·
[Submission metadata](metadata.json) ·
[Checksums](SHA256SUMS.txt)

## Build and replay

Extract the TeX ZIP and run this from its root with TeX Live 2023 and
latexmk 4.83.

    SOURCE_DATE_EPOCH=1795564800 FORCE_SOURCE_DATE=1 latexmk -pdf -interaction=nonstopmode -halt-on-error -jobname=primitive-balanced-image-spectrum main.tex

The ZIP contains the complete LaTeX sources. The
[verification supplement](supplement/) supplies exact replay commands,
inputs, and the mathematical dependency inventory. The
[verification record](verification.json) identifies the reviewed PDF
and distinguishes finite replay, PDF-only internal review, and external
mathematical review. The [literature comparison](ordinary-reduction-literature.md)
records the precise arithmetic theorem overlap and scope.

## Hexagon preparation

The TeX archive is the upload candidate. The verification archive is
separate from the TeX upload. The contributor will choose the distribution
license and confirm the deposit declarations at submission.

[AI provenance](AI-PROVENANCE.md) · [Announcement](ANNOUNCEMENT.md) ·
[Mathematical dependency inventory](supplement/DEPENDENCIES.md)
