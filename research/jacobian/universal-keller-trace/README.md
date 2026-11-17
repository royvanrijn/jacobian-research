# Trace compression and canonical realization spaces

Further extension: [universal sextics and sparse stable classification](../universal-keller-sextic/README.md) includes rank six in the N-2 bound and closes the characteristic-zero stable-orbit gap on the trace boundary.

For N>=7, an explicit universal determinant-one Keller map now uses N-2 coordinates. A trace-zero primitive generator removes the penultimate coefficient in addition to the previous scaling coordinate. Complete fibre universality, S_N monodromy and absolute stable atomicity remain.

The [proof](PROOF.md) also extends universality to every infinite field of characteristic other than two, and constructs a canonical correspondence whose fibres are open projective spaces of generator lines. Every fixed algebra has dense realizing targets, with compatible weak approximation over number fields.

- UATRACE1: the N-2-dimensional universal map.
- UAODD1: all odd characteristics and a sufficient finite-field bound.
- UALIFT1: rational lift spaces, canonical family correspondence and dense targets.
- UABOUND1: complete pi-zero fibres and a real obstruction to double depression.

[Exact controls](verify.py) · [receipt](receipt.json) · [sources](sources.json) · [mathematical authority](../../MATH_STATUS.json) · [preceding N-1 construction](../universal-keller-compression/README.md).

    .venv/bin/python research/jacobian/universal-keller-trace/verify.py

Written universal proofs and exact representative controls; no formalization, independent whole-proof replay, external review or priority claim.

