# Seven high-contact polynomial quintics and exclusion of contact41

14 November 2026. PF2R7B5RAT1 proves that no member of the closest-quintic
pencil has polynomial normalization degrees(3,5). Such an additional
image must meet C in finite intersection multiplicity at least10.

PF2R7B5REC1 explicitly reconstructs all degree(3,5) polynomial images
with infinity contact at least38: a finite free correspondence of degree7
over every point of each complete first-image family. Seven counts
multiplicity; a nonempty open has seven distinct choices. Their contacts
are38 or40. Actual B is not forced into this chart or contact range.
Compatibility with the Keller covering data remains open.

```sh
timeout 60s .venv/bin/python research/jacobian/jc2-r7-quintic-normalizations/verify.py
```

[Proof](PROOF.md) · [equations and normalization formulas](normalizations.json) ·
[replay](receipt.json) · [author-side audit](REVIEW.md) ·
[validation](validation.json).

Written geometric interpretation with exact monic Groebner and literal
Fraction replays; no formalization, independent external review or JC2
closure. Discovery artifacts are indexed in [discovery.json](discovery.json).

**Subsequent incidence results:** [PF2R7B5GEN1/PF2R7B5CUSP1/PF2R7B5N1](../jc2-r7-quintic-incidence/README.md)
certify generic local compatibility, restrict a nonimmersive second image
to160 first-image coefficient points, and exclude singleton fibres in
the oriented nonimmersed rowN under this degree/contact hypothesis.
The global covering and other degree/contact charts remain open.
