# Global boundary attachment: tame-node exclusion and the r7 missing-image gate

Applying Borisov's plane Stein-boundary theorem gives two results beyond
the local radial tests:

- **PF2BATT1:** omitted boundary curves are disjoint and unibranch in the
  affine finite normalization. The collided finite map uv=w^n -> (u,v)
  cannot occur there for any n>=2. A generically unramified omitted curve
  maps with degree one to an immersed polynomial curve.
- **PF2R7BND1:** at geometric degree six, under the retained single
  ramified-row hypothesis (e,f)=(2,1), either classified nonimmersed r7
  curve forces two or three additional unramified boundary primes. Their
  immersed images differ from that cusp curve. In particular a complete
  nonproper-value locus consisting only of that curve is excluded, even
  if arbitrary unramified rows over it were allowed.

The second result uses stratified Euler characteristic, not an assumed
complete one-row Orevkov ledger. It does not rule out additional images,
the immersed r7 family, normal r5, nonnormal slices, or JC2.

The later [adjoint-fibration proof](../r7-adjoint-exclusion/README.md)
excludes both nonimmersed coefficient classes under PF2R7BND1's explicit
single-ramified-row hypothesis, retaining the additional unramified
images until the contradiction. Its extra Borisov source-model input is
identified separately; the earlier attachment calculation remains evidence.

[Proof and exact imported hypotheses](PROOF.md) · [checker](verify.py) ·
[replay](receipt.json) · [primary-source pins](sources/manifest.json).

```sh
timeout 60s .venv/bin/python research/jacobian/jc2-boundary-attachment/verify.py
```

The global boundary theorem is an explicit external input, not certified
by the finite checker. No independent full reproof of Borisov, formal
verification, external review of these deductions or novelty priority is
claimed. [MATH_STATUS.json](../../MATH_STATUS.json) is authoritative.
