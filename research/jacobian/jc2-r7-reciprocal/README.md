# The remaining reciprocal cubic case is excluded

14 November 2026. [PF2R7C3ALL1/PF2R7ROW1](PROOF.md) exclude the final
reciprocal cubic configuration. Combined with the earlier quartic and
varying-direction results, this closes the entire degree-six,
**single ramified affine-prime (e,f)=(2,1), gap-21 stratum** without
source finality, retaining all unramified omitted primes.

This supplies an additional algebraic proof of the scoped exclusion
already obtained by [PF2OMIT5](../jc2-omitted-divisors/R7_INFINITY_COMPLETION.md).
The rational cubic contact matrix reduces the marked-root problem to a
centered coefficient scheme and two linear-consistency branches. Exact
scheme memberships, univariate units and a two-variable b^13 certificate
exclude all of them.

PF2R7RAM2 gives the next global gate: any degree-six Keller map containing
such a prime must have exactly one more ramified prime, with a different
target image and one of three explicit Euler ledgers. That second image
and the original terminal passport's complete ramification data remain open.
The current [residue](../jc2-omitted-divisors/RESIDUE_LEDGER.md) and
[local-boundary](../jc2-omitted-divisors/SIMPLE_RAMIFICATION_AND_RADIAL.md)
results already fix all residual residue degrees, local lengths and radial
packets in the immersed case. Their global realization remains the frontier.

```sh
timeout 60s .venv/bin/python research/jacobian/jc2-r7-reciprocal/verify.py
```

[Proof](PROOF.md) · [certificates](certificates.json) · [replay](receipt.json) ·
[validation](validation.json) · [preceding quartic result](../jc2-r7-universal-quartic/README.md).

This is a written proof with exact algebra, without formalization or
independent external review. It does not close whole normal-r7 or JC2.
