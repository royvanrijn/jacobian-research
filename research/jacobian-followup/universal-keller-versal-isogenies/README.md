# General versal isogeny self-covers

A free field-versal torus base exists for every finite group. Fixed self-isogeny covers then have an exact twisted-H1 criterion. Norm tori give the sharp coprimality condition for every transitive group, while cyclotomic principal ideals give additional primitive covers, including degree8 with AGL1(F8) monodromy.

- [Proofs](PROOF.md): free bases, exact cohomology criterion, general norm boundary, all odd prime powers, cyclotomic atomicity and stable ideal obstruction, and the published all-affine limitation.
- [Literature and remaining gaps](LITERATURE_AND_GAPS.md).
- [Exact controls](verify.py), [receipt](receipt.json), [source pins](sources.json), [integration](integration.json), [focused validation](validation.json).

Replay `python verify.py` with [requirements](requirements.txt). The new self-covers realize specified group torsors and remain outside affine space. Field versality supplies no natural choice of a realizing point and no all-affine scheme universality. Independent proof review remains needed; OP-KDESC remains open.

The [semilinear extension](../universal-keller-semilinear-spectrum/PROOF.md) now fills every remaining odd binary degree and proves the exact prime-power strong-atomic degree spectrum. Its changed ideal lattice and holonomy automorphism bypass the fixed-O principal-prime restriction while preserving that restriction in its original scope.

The [canonical trace-filter compiler](../universal-keller-canonical-filter/PROOF.md) realizes every finite etale algebra as a complete marked fibre of a functorial whole-affine-fibre bundle self-map. This is neither a naturally chosen torsor point nor one fixed scheme receiver: the bundle varies with A and the map is non-atomic.
