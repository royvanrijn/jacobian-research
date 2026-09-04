# NS0024 direct-QQ Inose source obstruction

Date: 2026-09-04.

<!-- status-consumer: EC-K3-NS0024-DIRECT-QQ-INOSE-OBSTRUCTION e87afc1b3529a07f -->
<!-- status-consumer: EC-K3-NS0024-QQ-MARKING-OBSTRUCTION b7f0cf002c0411fe -->

## Result

The symbolic `2E8/MW1` Inose source is not a direct arithmetic source over
`QQ` for `NS0024`.

For nonisomorphic elliptic curves `E1,E2`, the Mordell--Weil lattice of the
Inose fibration is

```text
Hom(E1,E2)<2>.
```

Thus its required height-`950` generator corresponds to an isogeny of degree
`950/2=475`.  Utsumi's construction produces the rational section over the
field of definition of an isogeny `phi:E1 -> E2`.

The Mazur--Kenku classification says that the complete list of degrees of
cyclic isogenies of elliptic curves over `QQ` is

```text
1,...,19, 21, 25, 27, 37, 43, 67, 163.
```

Since `475=19*25` is not on this list, `X0(475)(QQ)` has no noncuspidal
rational point.  In particular there is no pair `E1,E2/QQ` joined by a
`QQ`-rational cyclic `475`-isogeny.  The previously advertised instruction
to find such a point over `QQ` was therefore not merely unfinished; it was
impossible.

The exact application is replayed by
[`scripts/certify_ns0024_direct_qq_inose_obstruction.py`](scripts/certify_ns0024_direct_qq_inose_obstruction.py),
with compact output in
[`../artifacts/generated-results/elkies-k3-ns0024-direct-qq-inose-obstruction-v1.json`](../artifacts/generated-results/elkies-k3-ns0024-direct-qq-inose-obstruction-v1.json).
The classification theorem itself is an external input, not reproved by the
checker.

## Descent boundary

This direct argument removes only the rational level-structure route. The
subsequent Fricke-quotient argument in
[`NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](NS0024_QQ_MARKING_OBSTRUCTION_2026-09-04.md)
also excludes the rational marking needed by the foundry milestone: Momose's
theorem leaves only cusps and CM points on `X0+(475)(QQ)`, while the rank-one
Inose frame requires a noncuspidal, non-CM point.

A quadratic `Q`-curve or larger-quotient model may still exist with a proper
Galois-invariant sublattice, but it cannot make all nineteen NS classes
rational. A Galois action that exchanges the two `E8` fibres or sends the
generator to its negative is therefore geometric or number-field data, not
an arithmetic source for the stated `QQ(t)` objective.

## Source decision

All `NS0024` source routes are now parked for the arithmetic MW17 milestone.
The semistable `A3+A4+A6/MW4` reconstruction and the `D5+E8/MW4`
completed-core route remain valid geometric investigations, but neither can
produce a full rational NS0024 marking over `QQ`. The live different-NS
objective and its replacement `NS0031` source gate are in
[`DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md`](DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md).

## References

- B. Mazur and M. A. Kenku's complete classification, summarized as Theorem
  1.1 and Table 1.1 in B. Banwait, F. Najman, and O. Padurariu,
  [*Cyclic isogenies of elliptic curves over fixed quadratic fields*](https://arxiv.org/abs/2206.08891).
- K. Utsumi,
  [*The Mordell--Weil lattice of an Inose surface arising from isogenous elliptic curves*](https://arxiv.org/abs/2209.02463),
  especially Proposition 3.1 and Theorems 5.1--5.2.

## Replay

```bash
python3 elkies-k3/scripts/certify_ns0024_direct_qq_inose_obstruction.py
python3 elkies-k3/scripts/certify_ns0024_direct_qq_inose_obstruction.py --check
```
