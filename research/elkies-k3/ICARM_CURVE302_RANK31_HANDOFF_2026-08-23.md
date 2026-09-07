# ICARM curve 302 rank-31 handoff to the Elkies K3 programme

Date: 2026-08-23

Update, 2026-08-31: the published compact rootless `R17` equation and its
seventeen sections are now available.  The exact degree-24 equation
`j_R17(t)=j_302` has an irreducible degree-24 reduction modulo `397`, so curve
302 is not a direct rational specialization of that chart.  The historical
handoff below explains the earlier programme decision; the current exact
boundary is in
[`../elliptic-curves/notes/ELKIES_BISECTION_VISIBILITY_AND_RECORD_CURVES.md`](../elliptic-curves/notes/ELKIES_BISECTION_VISIBILITY_AND_RECORD_CURVES.md).

## New benchmark

ICARM curve 302 is now independently certified in this repository to satisfy

```text
rank E(Q) >= 31
```

unconditionally. The certificate, exact point data, and arithmetic audit are:

```text
elliptic-curves/cas/icarm_curve302.py
elliptic-curves/cas/verify_icarm_curve302_rank31.py
elliptic-curves/cas/check_icarm_curve302_rank31_pinned.py
elliptic-curves/notes/ICARM_CURVE302_RANK31.md
artifacts/generated-results/elliptic-curves/icarm_curve302_rank31_v1.json.gz
```

The 31 public points have full rank in an exact product of finite quotients
`E(F_p)/2E(F_p)`, independently cross-checked by a full-rank
Brumer--Cremona quadratic-character matrix. Rational torsion is trivial. The conditional BSD+GRH
rank-31 closure reported publicly is not used by the repository proof.

## What this changes

The external raw-rank benchmark is now 31, so the next record target is 32.
For the reconstructed generic rootless `MW17` programme, a rank-31
specialization would require fourteen additional directions beyond the generic
sections, one more than the rank-30 benchmark required.

This did **not** change the selected H3 lattice/equation route.  The later
published rootless equation has now supplied the exact test and excludes a
direct rational specialization of curve 302.

## What is not known

No inspected public source gives a family equation or specialization parameter
for curve 302. There is no exact evidence yet that it comes from

```text
H3 E7+E8/MW2 -> ... -> rootless/MW17.
```

Do not merge curve 302 into the H3 construction ledger as a specialization
endpoint. Keep it as an external benchmark and an unresolved provenance target.
The full boundary is documented in
[`../elliptic-curves/notes/ICARM_CURVE302_CONSTRUCTION_INVESTIGATION.md`](../elliptic-curves/notes/ICARM_CURVE302_CONSTRUCTION_INVESTIGATION.md).

## Exact matching data

The curve's reduced `j`-invariant has SHA-256

```text
5939208330113d89ae063d62053f0c8383e18b3a564919b86f86a02a4d13a550.
```

Its bad-fibre fingerprint is

```text
I15@2 + I4@3 + IV@5 + I6@7 + I4@11 + I5@13
+ I2@19 + I2@23 + I3@29 + I2@37 + I2@41 + I2@73
+ I2@131 + I2@167 + 6 I1.
```

These remain cheap rejection filters for any different candidate family.  A
positive `j` match would still require the exact Q-isomorphism and all
seventeen generic section transports.

## Programme priority

The current strategic order is:

1. keep curve 302 as an external rank-31 and finite-quotient regression;
2. do not search for a nonexistent rational parameter in the published chart;
3. test a different construction only when it supplies an explicit invariant
   map or isogeny certificate;
4. search the published `R17` parameter geometry for rank 32 rather than
   spending the main effort on a conjecturally absent 32nd point on curve 302.
