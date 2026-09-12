# Archived projects

Archived on 12 September 2026 at the user's request. Active research now covers
[elliptic curves](../../elliptic-curves/README.md) and their
[supporting K3 constructions](../../elkies-k3/README.md).

| Archived programme | Preserved entry point | Claim catalogue |
|---|---|---|
| Keller, cancellation and arithmetic | [Verified core](verified/README.md) | [Core claims](../../index/core.md) |
| Gaussian moments, GVC, SIC and geometry | [Geometry map](extended-geometry/README.md) | [Geometry claims](../../index/geometry.md) |
| Hessian, HC4 and Schur reductions | [Canonical sources](../../index/hessian.md) | [Hessian claims](../../index/hessian.md) |
| Plane Jacobian conjecture | [Plane map](plane-jc/README.md) | [Plane claims](../../index/plane-jc.md) |
| Formal verification and papers | [Formal projects](formal/README.md) · [Papers](papers/README.md) | [Formal claims](../../index/formal.md) · [Papers](../../index/papers.md) |

The live [shared authority](../../MATH_STATUS.json) retains all 901 archived
claims with `programme_status: archived`. Their mathematical states, dependency
edges, proof assurances and checker hashes were preserved. Open questions and
unfinished cleanup remain open and unfinished. The local
[MATH_STATUS.json](MATH_STATUS.json) is a frozen pre-archive snapshot, for
historical replay only; it is not a second current authority.

The [manifest](MANIFEST.json) records 2,319 relocated files, original paths,
SHA-256 hashes, frozen navigation copies and explicit links to shared inputs.
Files were moved without rewriting their proof or replay contents. Some older
certificates remain in the shared artifact store because their paths are pinned.

## Reading and replay

From the repository root, use `research.py search QUERY --history`,
`research.py work --history`, or `research.py show ID` through
`python3 research/scripts/research.py`. The archived programme's notes and old
Makefile describe historical work, not current instructions to restart it.

For an explicitly scoped archival replay, first inspect the exact command and
inputs. Old relative project paths now resolve from
`research/archive/non-elliptic/`. Shared EC code, artifact storage and prior
archives are linked, without duplicating those inputs. Frozen commands that
use absolute checkout paths or a historical environment still require their
original environment or a documented path adjustment. No replay was performed
as part of this archive move.

Historical [method records](knowledge/lessons.json), [work proposals](knowledge/work_items.json),
[checklist reconciliation](knowledge/legacy_work_review.json) and
[source reviews](knowledge/partial_reviews.json) preserve the pre-archive scope.
Updated navigation reads the live records in [research/knowledge](../../knowledge/).
