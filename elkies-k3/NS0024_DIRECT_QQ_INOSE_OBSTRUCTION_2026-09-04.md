# NS0024 direct-QQ Inose source obstruction

Date: 2026-09-04.

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

This removes only the direct level-structure route.  It does not prove that
`X0(475)/<w475>` or a larger Atkin--Lehner quotient has no relevant rational
point, and it does not rule out a quadratic `Q`-curve construction in which
the two elliptic curves or the isogeny are conjugate.

Such a descent would not pass the foundry source gate merely because the
Inose equation descends.  It must also show that the two `E8` configurations,
the fibre and zero, and the height-`950` Mordell--Weil generator descend as
nineteen individual `QQ`-rational divisor classes.  A Galois action that
exchanges the two `E8` fibres or sends the generator to its negative gives a
smaller invariant Neron--Severi rank and is insufficient.

## Source decision

The active source programme is therefore the semistable
`A3+A4+A6/MW4` reconstruction.  It already has the exact resolved component
profiles, a determinant-`95/14` four-section height lattice, a thirteen-edge
marked degree-two corridor, and a prepared first compiler adapter.  Its next
gate is a common characteristic-zero producer for the surface and all four
sections, followed by the rational rank-19 and Picard-rank checks.

The `D5+E8/MW4` completed-core route remains a geometric control.  It becomes
an arithmetic source contender only if an independent explicit equation and
rational marking are found; the abstract `17,13,7` Kneser sequence is not an
equation route.

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
