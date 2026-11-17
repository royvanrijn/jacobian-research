# Research

This is the working area for AI-assisted mathematics. Published papers are listed in the [portfolio](../README.md); manuscripts awaiting publication live in [unpublished](../unpublished/README.md).

| Topic | What we study | Status |
|---|---|---|
| [Moments and vanishing](moments/README.md) | GVC, Gaussian moments and Image/Mathieu problems. | GVC(2), GMC(2) and counterexamples are published; cubic coverage and new manuscripts remain open. |
| [Jacobian and Keller maps](jacobian/README.md) | Polynomial maps, boundary geometry and universal fibres. | Restricted constructions and exclusions are established; the Jacobian conjecture remains open. |
| [Hessian and Schur reductions](hessian/README.md) | HC4, polynomial graphs and termination. | Several low-degree branches are closed; global coverage and termination remain open. |
| [Factorial and rigidity](factorial/README.md) | Factorial moments, Image spaces and polynomial rigidity. | Binary proof audits and restricted pencil theorems are complete internally; general FC and SFC remain open. |
| [Dixmier](dixmier/README.md) | Weyl-algebra maps, formal lifts and residue transport. | Partial reductions and obstructions; the full conjecture remains open. |
| [Meridional rank](meridional-rank/README.md) | Knot-group representations and branched covers. | Exact results in known families; the full conjecture and novelty questions remain open. |
| [Arithmetic fibres](arithmetic/README.md) | Prescribed finite etale Keller fibres and arithmetic failures. | Constructive results have draft papers; review and publication remain outstanding. |
| [Foundational proof library](foundations/README.md) | Earlier proofs, formal projects and exact checkers. | Mixed proved, partial and superseded work; consult the ledger for each statement. |

The [mathematical ledger](MATH_STATUS.json) is authoritative for individual claims and their review, replay and formal-verification flags. [STATUS.md](STATUS.md) is its generated overview. Publication is recorded separately in each paper's metadata.

For a specific question, start with [research memory](KNOWLEDGE_BASE.md), the [claim catalogue](index/README.md), or the [work ledger](knowledge/WORK_LEDGER.md). [Historical branches](archive/README.md) preserve closed and superseded routes.

```sh
python3 research/scripts/research.py search "your topic"
python3 research/scripts/research.py show CLAIM-ID
make check-navigation
```
