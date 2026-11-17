# The entire r7 quartic-contact family is excluded

14 November 2026. [PF2R7QALL1](PROOF.md) excludes conic contacts over every
quartic polynomial parameter for the whole birational degree-(6,10),
gap-21 family. It includes immersed and multiple-derivative-root strata.
Universal minor coefficients force the quartic parameters, and two
literal four-variable unit certificates exclude the remaining branches.

PF2R7QEX1 consequently excludes the entire quartic infinity escape under
the retained degree-six Keller and entire single ramified e2,f1 prime
hypothesis. Its remaining **reciprocal cubic** case is now excluded by
the [exact cubic-contact certificates](../jc2-r7-reciprocal/README.md).
This gives an additional route to PF2OMIT5's entire single-prime exclusion.
No source finality is used.

```sh
timeout 60s .venv/bin/python research/jacobian/jc2-r7-universal-quartic/verify.py
```

[Proof](PROOF.md) · [certificates](certificates.json) · [replay](receipt.json) ·
[validation](validation.json) · [preceding two-class result](../jc2-r7-quartic-contact/README.md).
[Author-side audit](REVIEW.md) separates exact algebra from the inherited
geometric proof obligations.

This is a written proof with exact algebra; no formalization or independent
external review. The complete ramification of the original terminal
passports, global second-image data, other r7 patterns and JC2 remain open.
