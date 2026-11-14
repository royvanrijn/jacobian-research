# Missing radial blocks: intersection equations and a chain exclusion

The [pointed-pencil budgets](../hc4-pointed-pencil/README.md) now have a
local intersection refinement. On a complete resolution over each omitted
normalization point, the radial multiplicities e, direction degrees f and
exceptional intersection matrix Q must satisfy

\[
 Qe=-f,\qquad \ell_q=e^Tf=f^T(-Q)^{-1}f.
\]

Here ell_q is the omitted local fibre length, and their sum is d-1.
**HC4MYRB1** also gives the canonical-order test and an all-length chain
congruence. In particular a resolved chain with every weight at most -2
cannot have just one degree-one radial horizontal component.

**HC4MYRB2** checks the exact scope using the retained uv=w^n node
models. Their minimal resolution resolves the radial pencil for even n;
odd n needs an additional blowup. That extra -1 curve has radial degree
one, so it must not be excluded by applying the minimal-chain test to an
unresolved pencil. These are local normal-cover controls, not Keller maps.

[Proof and attribution](PROOF.md) · [checker](verify.py) · [replay](receipt.json) ·
[primary-source pins](sources/manifest.json).

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-radial-intersection/verify.py
```

The singleton-zero Keller gate, OP-HC4-SCHUR and JC2 remain open.
Written geometric arguments are supported, not formally certified, by
the finite replay. No external review or novelty priority is claimed.
[MATH_STATUS.json](../../MATH_STATUS.json) is authoritative.

The later [global attachment theorem](../jc2-boundary-attachment/README.md)
excludes the entire collided-node projection family from a plane Keller
finite normalization. The local radial calculations here remain valid;
the global exclusion uses the separate Borisov boundary input.
