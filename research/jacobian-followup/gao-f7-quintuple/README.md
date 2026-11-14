# F7: all 118 fivefold-root targets and their seven-point fibres

**GAOF7FIVE1:** exactly 118 normalized F7 targets Y have a root of
multiplicity at least five. Each has one fivefold root and seven other
simple roots, so its F7 fibre has exactly seven source points.

[Proof and scope](PROOF.md), [exact parameterization](parameterization.json),
[5490-check receipt](receipt.json), [separate verifier](verify.py),
[integer certificate](certificate.json), [finite witness](finite-witness.json),
[coefficient input](quintuple-input.json), [dependencies](requirements.txt).

The complete set is the known origin plus117 values indexed by a degree 117
integer polynomial with rational recovery formulas. The proof covers every
exceptional denominator locus and proves that every retained root occurs.
On c!=0 these give118 disjoint C* curves in the original target chart.

This excludes three-point types(5,4) and(5,2,2). Combined with the earlier
cube and maximum-multiplicity proofs, only(4,3,2) and(3,2,2,2) remain open.
The verifier takes about 10 seconds and uses finite FREE algebra unit tests
with proved characteristic-zero transport. It runs no Groebner search.

Independent end-to-end review, formal verification and novelty assessment
remain outstanding. Global three-point existence, full F7 strata, F5 and
JC2 are open. The retained probes and generation snapshots are evidence
for reproduction, not a second status ledger.

The subsequent [quadruple-root incidence proof](../gao-f7-quadruple/PROOF.md)
proves generic degree eight and finitely many nonordinary normalized
quadruple-root targets, including the unresolved(4,3,2) case. It also
classifies all24 targets on the slope-coordinate pole chart.
