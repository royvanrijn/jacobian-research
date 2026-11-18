# Counterexamples to the Mathieu Property

*From Explicit Quotient Constructions to a General Classification*

**Published preprint**

- **Zenodo · 17 November 2026.** [PDF](paper.pdf) · [LaTeX](paper.tex) · [Zenodo record](https://zenodo.org/records/24310836) · [DOI 10.5281/zenodo.24310836](https://doi.org/10.5281/zenodo.24310836) · [Metadata](zenodo-metadata.json) · [Checksums](deposited-files.json)
- **Hexagon v1 · 18 November 2026.** [PDF](hexagon-v1/paper.pdf) · [Source bundle](hexagon-v1/source.zip) · [Metadata](hexagon-v1/hexagon-metadata.json) · [Checksums](hexagon-v1/deposited-files.json) · [hexagon:2611.00556v1](https://hexagonmath.org/2611.00556v1)

[Publication metadata](metadata.json) records both deposits.

The [current repository edition](../) matches the archived Zenodo PDF and LaTeX byte for byte. The Hexagon source bundle contains the same LaTeX and the AI-provenance files. Hexagon supplies its own PDF. The archived checksums identify each deposit.

The paper classifies the reductive homogeneous spaces whose normalized Reynolds integration kernel is a Mathieu–Zhao subspace and constructs fixed-multiplier counterexamples when the property fails.

The [verification record](../verification.json), [source audit](../source_audit.json) and [exact checker](../verify.py) describe the internal review and finite controls supporting the repository edition. The finite checks do not certify the all-order proofs or classification. Independent replay, formal verification and external mathematical review are not claimed.

Build the deposited source from this directory with `latexmk -pdf -interaction=nonstopmode -halt-on-error paper.tex`. Use a separate output directory to preserve the archived PDF bytes.
