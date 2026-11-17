# Both nonimmersed r7 classes excluded without source finality

**Universal extension:** [PF2R7QALL1/PF2R7QEX1](../jc2-r7-universal-quartic/README.md)
now exclude the entire quartic-contact family, including every immersed
specialization. This directory preserves the earlier two-class certificate
and the universal polynomial-base bridge. Only the special reciprocal
cubic family remains in the declared single-ramified-prime setting.

14 November 2026. [PF2R7EX2](PROOF.md) excludes both exact cusp classes
under the degree-six Keller and entire single ramified affine-prime
(e,f)=(2,1) hypothesis, retaining every unramified omitted prime.

The cubic exclusions leave one quartic contact escape. Its parameter
must invert to a quartic polynomial. Exact conic-contact minors and four
univariate Bezout certificates exclude every such polynomial for both
classes. This replaces the older finality-dependent geometric bridge;
the other named boundary and positivity inputs remain.

```sh
timeout 60s .venv/bin/python research/jacobian/jc2-r7-quartic-contact/verify.py
```

[Proof](PROOF.md) · [certificates](certificates.json) · [replay](receipt.json) ·
[validation](validation.json) · [cubic reduction](../jc2-r7-cubic-cover/README.md).
[Author-side obligation audit](REVIEW.md) distinguishes written geometry
from the exact algebraic replay.

Special multiple-root reciprocal cases, additional ramification patterns,
normal r5, nonnormal slices and JC2 remain open.
This is a written geometric proof with exact algebra; no formalization
or independent external review is claimed.
