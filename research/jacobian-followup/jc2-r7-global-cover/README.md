# Global Nori contact bound and corrected two-image genus budgets

14 November 2026. PF2R7NORI1 applies Nori Proposition3.27 to the required
S6 finite cover. For a nodal degree(3,5) second image B with mutual
finite singular-locus avoidance, let T be the sum of finite C/B tangency
multiplicities. Then

    I_infinity(C,B)+2*T >=34.

Thus finite-transverse pairs of contacts20 through32 are excluded.
Contact34 is the equality case; the previous38/40 candidates survive
this test. RowA with a distinct H needs the stated avoidance hypothesis
checked; rowC and shared-image rowA force it by their branch counts.

PF2R7LINE1 gives general-target-line genera3--5 in rowA and5--7 in
rowC for a quintic B. The old single-prime genus2 adjoint pencil does
not apply. Neither the actual B degree chart nor its finite incidence
pattern is inferred from these budgets.

```sh
timeout 60s .venv/bin/python research/jacobian-followup/jc2-r7-global-cover/verify.py
```

[Proof and literature comparison](PROOF.md) · [exact budgets](budgets.json) ·
[primary sources](sources/manifest.json) · [replay](receipt.json) ·
[validation](validation.json) · [author-side audit](REVIEW.md).

Authored geometric deduction with imported Nori topology and exact
algebraic checks; no formalization, independent external review or JC2
closure. Global cover/source realization remains open.
