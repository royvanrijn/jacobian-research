# ICARM curve 302 rank-31 handoff to the Elkies K3 programme

Date: 2026-08-23

## New benchmark

ICARM curve 302 is now independently certified in this repository to satisfy

```text
rank E(Q) >= 31
```

unconditionally. The certificate, exact point data, and arithmetic audit are:

```text
elliptic-curves/cas/icarm_curve302.py
elliptic-curves/cas/verify_icarm_curve302_rank31.py
elliptic-curves/notes/ICARM_CURVE302_RANK31.md
artifacts/generated-results/elliptic-curves/icarm_curve302_rank31_v1.json
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

This does **not** change the selected H3 lattice/equation route. It strengthens
the reason to complete it: without an explicit generic rootless equation and
section map, curve 302 cannot be tested as a specialization and its search
neighbourhood cannot be reproduced.

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

These are cheap rejection filters once a generic candidate family is explicit.
A positive `j` match still requires the exact Q-isomorphism and all seventeen
generic section transports.

## Programme priority

The strategic order remains:

1. turn the now-exact q24 horizontal section into the resolved-RR `D12/MW5`
   child equation, binary-quartic Jacobian, and fibre certificate, then continue
   the exact neighbour compilation to the generic rootless `MW17` equation;
2. derive its exact `j`-map and solve against curve 302;
3. if a specialization exists, recover the `17+14` section decomposition;
4. search the nearby parameter geometry for rank 32 rather than spending the
   main effort on a conjecturally absent 32nd point on curve 302 itself.
