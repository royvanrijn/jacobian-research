# AI Math Research

Mathematical research by Roy van Rijn, developed with AI research assistants, especially OpenAI Codex. This repository brings together written proofs, exact computations and Lean formalizations around the Jacobian Conjecture, polynomial maps and moment-vanishing problems.

Start with the released papers below. Each paper directory opens with the best current source and PDF; its `published/` archive preserves the exact Zenodo edition, submission metadata and checksums. Current revisions can differ from Zenodo or Hexagon and do not inherit their public-deposit status. Independent review and the scope of formal verification are stated separately.

## Released papers

| Paper | Main result | Read |
|---|---|---|
| **The Generalized Vanishing Conjecture: The Two-Variable Theorem and the First Failing Dimension** | Proves GVC in two variables and constructs a homogeneous failure in dimension three. | [Paper and verification](results/generalized-vanishing-two-variables/) · [Zenodo](https://zenodo.org/records/21803325) · [Hexagon](https://hexagonmath.org/2610.00580v2) |
| **A Filtered Frobenius Proof of the Two-Variable Gaussian Moments Conjecture** | An alternative proof of the known binary classification and sharp cutoff; Lean coverage has two explicit mathematical inputs. | [Paper and verification](results/gaussian-moments-two-variables/) · [Zenodo](https://zenodo.org/records/21534699) |
| **A Four-Term Counterexample to the Special Image Conjecture in Three Pairs** | An explicit four-term three-pair construction, with the later two-pair boundary proved in the current edition. | [Paper and verification](results/three-pair-image-counterexample/) · [Zenodo](https://zenodo.org/records/21634058) |
| **Sparse Minimality of Gaussian-Moments Counterexamples in Dimension Three** | Excludes degree-at-most-four counterexamples with at most four monomials in Gaussian Witt coordinates. | [Paper and verification](results/sparse-minimality-gaussian-moments-dimension-three/) · [Zenodo](https://zenodo.org/records/21534917) |

## In preparation

[Unpublished papers](unpublished/) include the Gaussian dimension-degree boundary, **Counterexamples to the Mathieu Property** and smaller drafts. Revisions of released papers stay with their current editions in `results/`. Each directory records what it proves, its review status and the remaining work.

[Research topics](research/) contain the current results, open questions, proofs, code and supporting evidence. The [scope audit](research/SCOPES.md) and [dependency graph](research/DEPENDENCIES.md) show exact claim limits, proof routes and unresolved obligations.

## Working on research and papers

The [research and paper working guide](working-guide/README.md) describes repository organization, mathematical checks, readable paper writing and review. Future agents should start there and follow the [repository instructions](AGENTS.md).
