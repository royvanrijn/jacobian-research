# All-degree holonomy families

The universal quartic construction extends to degree-n-squared norm self-covers on every Z_m. The covers preserve volume; their specified G_n torsors are field-versal exactly away from multiples of three. Prime-square maps are strong absolute stable atoms except at prime3, where same-variety atomicity survives but a different smooth intermediate gives two degree-three factors.

- [Proofs](PROOF.md): general projection blocks, commutant/outer-automorphism criterion, complete Z_m endomorphisms and factor lengths, norm-torsor family and the all-m canonical quartic exclusion.
- [Literature and limitations](LITERATURE_AND_GAPS.md).
- [Exact controls](verify.py), [receipt](receipt.json), [source pins](sources.json), [integration](integration.json), [focused validation](validation.json).

Replay `python verify.py` with [requirements](requirements.txt). Constant-group torsor universality assumes mu_n in the base; the all-endomorphism classification uses characteristic zero. Independent mathematical review remains needed. OP-KDESC remains open.

General extension: [arbitrary-group versal isogenies](../universal-keller-versal-isogenies/PROOF.md) replace special free-base searches with induced characters, generalize the norm criterion to all transitive groups and separate field versality from the published all-affine scheme obstruction.
