# The Dimension-Degree Boundary for Gaussian Moment Vanishing

**Unpublished · ready for review.**

Determines the complex Gaussian dimension-degree boundary and proves the sharp ternary-cubic mixed-moment cutoff.

[PDF](manuscript.pdf) · [LaTeX](manuscript.tex) · [Metadata](metadata.json) · [Verification](verification.json) · [Build instructions](BUILD.txt)

The boundary is universal purity-to-safety for n <= 2, for degree <= 2 in any dimension, and for ternary cubics. In the last case, the sharp mixed-moment cutoff is m > 3 deg Q. The manuscript states the imported binary and quadratic results separately from its ternary-cubic classification and reproduces the two minimal counterexamples.

For a first reading, follow Sections 1–7 and 12. Sections 8–11 contain the full formal and analytic reductions, located by the proof map and leading-cubic case table. The all-degree constant-fibre theory and its scalar cancellation example are developed in the separate fibre paper.

From this directory, run `../../.venv/bin/python verify.py --build` to replay the selected exact algebra, compile every local TeX input, and compare every rendered PDF page. See [BUILD.txt](BUILD.txt) for dependencies and the explicit refresh command after source edits. These checks do not certify the general analytic proof. Independent mathematical review and publication remain pending.

The [two-variable revision](../../results/gaussian-moments-two-variables/) supplies a separate elementary proof with a sharp binary cutoff. The [marker-reduction paper](../gaussian-constant-fibre-base-marker-certificates/) combines finite base-marker reduction with its exact specialization obstruction. The [one-base constant-fibre paper](../gaussian-quadratic-fibre-jordan-cancellation/) treats arbitrary Jordan blocks. These all-degree extensions are separate from the dimension-degree proof.

The 17 November review corrections make the glued chain complex and degree-one comparison explicit, isolate the exclusion of finite-target contributions from spatial infinity, display the total/irregular/regular rank ledger, and specify the relative pairs used for the flat cut basis. The cubic table refers to nine geometric types.

Next: Specialist review of these comparisons and the complete qualitative cubic classification. These exposition repairs and the finite algebra checks do not independently certify the global analytic proof.
