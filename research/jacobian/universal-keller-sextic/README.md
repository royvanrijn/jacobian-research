# Universal sextics and stable trace moduli

Further extension: [trace quadrics and the rank-seven minimum](../universal-keller-quadrics/README.md) gives N-3-dimensional maps with an exact arithmetic criterion, classifies field-versality when N is invertible, and proves ktdim_k(S7)=4 over characteristic-zero bases containing sqrt(-14).

Rank six now fits the N-2 theorem: an explicit determinant-one map on A4 realizes every rank-six finite etale algebra over every infinite field of characteristic other than two. It has degree6, monodromy S6 and absolute stable atomicity.

The [proof](PROOF.md) uses the surviving coprime weights5 and7, setting a4=h^2,a5=0,a6=h^3. It also extends the characteristic-zero stable-orbit theorem across all intermediate coefficient-zero boundaries. Every fixed rank-N algebra occurs densely in an N-5-dimensional family of distinct stable classes, N>=6.

- UASEXT1: universal sextics in four coordinates, including a rational canonical lift space.
- UASPARSE1: stable classification for arbitrary intermediate zeros.
- UATRACEMOD1: stable trace families and the surviving N-5-dimensional presentation defect.

[Exact controls](verify.py) · [receipt](receipt.json) · [sources](sources.json) · [mathematical authority](../../MATH_STATUS.json).

    .venv/bin/python research/jacobian/universal-keller-sextic/verify.py

Universality holds in every infinite characteristic except two. Stable inequivalence is proved in characteristic zero. The proofs are not formally or externally verified.

