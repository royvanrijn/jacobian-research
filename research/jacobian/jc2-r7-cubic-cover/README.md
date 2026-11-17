# Cubic conic-cover gates: one exclusion and one multiple-root restriction

14 November 2026. [PF2R7V1/PF2R7G1](PROOF.md) use the actual field degree
to bound primitive conic coefficients. Without the disputed source-finality
input, they exclude the varying-direction parabola escape (indices2,4)
and force a finite multiple root of p' in the reciprocal-conic escape.

Only infinity indices (1,5) remain. Squarefree p' forces the **quartic**
contact escape; [PF2R7EX2](../jc2-r7-quartic-contact/README.md) subsequently
excludes that escape for both retained cusp classes without source finality.
The [universal quartic theorem](../jc2-r7-universal-quartic/README.md)
subsequently excludes every immersed quartic case too. The remaining
multiple-root reciprocal cubic case is now excluded by the
[exact cubic-contact certificates](../jc2-r7-reciprocal/README.md), an
additional route to PF2OMIT5's declared single-prime exclusion.
Global additional-image data remain open; no whole r7 or JC2 closure is claimed.

```sh
timeout 60s .venv/bin/python research/jacobian/jc2-r7-cubic-cover/verify.py
```

[Proof and two-row frontier](PROOF.md) · [checker](verify.py) ·
[replay](receipt.json) · [validation](validation.json) ·
[preceding three-escape theorem](../jc2-r7-universal-contact/README.md).

The geometric inputs and entire single-ramified-prime hypothesis are
retained explicitly. No formalization or independent external review.

The capped [quartic probe](quartic_probe.py) and its
[one-variable branch data](quartic_branch_discovery.json) are preserved
discovery history. The later quartic checker replays selected exact
certificates directly; it does not rerun that search.
