# Nine term Keller collision exclusion and a 49 term example

A normalized three-variable Keller collision with raw coordinate degrees at
most seven requires **at least ten nonlinear monomial occurrences**,
counted with coordinate multiplicity.

All 90 production shards and an independent global traversal finish with
no candidate support through size nine after the proved coordinate gates.
This proof uses the degree-controlled plane reduction and Moh's degree-100
theorem; no residual nine-term coefficient ideal reaches an algebra solver.

There is also an explicit **49-term normalized example** of degree seven.
It is optimal among the specified frames based on the boundary of the known
foundational map. The [extended transverse-frame theorem](../collision-frame-search/README.md)
now supplies a 48-term example, so the current normalized interval is [10,48].
Attainment at ten remains open.

- [Exclusion proof and scope](PROOF.md)
- [Explicit example and sharp minimum in the boundary frame class](NORMALIZED_WITNESS.md)
- [Complete search manifest](results/manifest.json)
- [Independent coverage](coverage-validation.json)
- [Positive controls and checkpoint regression](controls-validation.json)
- [49-term polynomial coefficients](normalized-49-witness.json)
- [Complete frame coefficient factorization](frame-family-certificate.json)
- [Verification](validation.json)
- [Scoped metadata checks and navigation limitation](integration.json)

From the repository root:

```sh
.venv/bin/python research/jacobian/sparse-nine/verify.py --replay-controls
```

That audits the terminal cover and actual replay receipts. A fresh independent
traversal can use a new checkpoint path; the search and resume commands are
in the proof. These results do not determine unrestricted support or degree
minima, nor the census of all dense degree profiles.

