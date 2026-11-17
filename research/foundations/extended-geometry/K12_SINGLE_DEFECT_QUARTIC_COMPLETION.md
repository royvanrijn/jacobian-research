# Quartic completion on the single-defect \(K_{12}\) graph families

## Status

This is an **exact bounded obstruction** over \(\mathbf Q\), not a
dimension-eleven lower bound. It extends the
[cubic parameterized frontier](K12_Z8_CUBIC_COMPLETION.md) by one target
degree on the five quadratic graph families having a single bad retained
component.

## Certificate

For each pivot

\[
j=7,9,10,11,12,
\]

the general quadratic graph coordinate is retained with all four or five
parameters. Target monomials of degrees one through four in the other ten
raw retained outputs give

\[
10+55+220+715=1000
\]

formal columns. The ten linear columns have no high-degree part, leaving 990
nonzero columns.

The sparse minor-first compiler selects rows at the literal parameter point
modulo \(1{,}000{,}003\), then reconstructs the selected minors over the full
rational parameter ring. For every family:

- the selected \(990\times990\) determinant is a nonzero constant;
- the \(991\times991\) augmented determinant is a nonzero constant multiple
  of it.

Therefore the high-degree defect is outside the full target-completion span
for every parameter value.

> **Single-defect quartic-completion obstruction.** None of the quadratic
> graph families with pivots \(z_7,z_9,z_{10},z_{11},z_{12}\) admits a
> one-stage target completion of degree at most four that restores degree at
> most three.

The remaining case at target degree four is the multi-defect \(z_8\) family.
Its sparse size profile is already substantially larger:

\[
992\text{ nonzero columns},\qquad
447{,}521\text{ source rows},\qquad
8{,}028{,}313\text{ parameter terms}.
\]

That case should use streamed modular row selection and reconstruct only the
selected entries; retaining its whole sparse parameter matrix is no longer
the right memory model.

## Reproduction

Run

```bash
make verify-k12-single-defect-quartic-completion
```

The exact generated record is
[`k12_single_defect_quartic_completion_frontier.json`](../artifacts/generated-results/k12_single_defect_quartic_completion_frontier.json).
