# Exact `074d9` cross-fibre carrier-transfer audit

Date: 2026-09-04  
Status: complete rigid census; exact canonical and cheapest-16 norm-eight screen; bounded exact late-label holdout audit

<!-- status-consumer: EC-K3-R17-074D9-RIGID-CROSS-FIBRE-TRANSFER abbedd192865f172 -->
<!-- status-consumer: EC-K3-R17-074D9-NORM8-CROSS-FIBRE-TRANSFER-16 262e405b0adbbb73 -->
<!-- status-consumer: EC-K3-R17-074D9-LATE-POINT-HOLDOUT 284e0f92def23419 -->

## Result

The proposed record-to-record carrier test is negative at both tested layers.
No member of the complete 39,120-character rigid bisection atlas splits at
both rank-at-least-29 fibres 356 and 385.  For one target-blind canonical
norm-eight trace, and again for the first sixteen finite-pole traces in the
same target-blind equation-cost ordering, every frozen off-diagonal
specialization is nonsplit.  Thus these searches produce neither a
record-anchored rank-at-least-18 cover nor a two-character rank-at-least-19
compositum.

The five repeated labels `P18,...,P22` also miss every predeclared low-degree
interpolation template when curves 351, 356, 376, and 377 are used for
reconstruction and curve 385 is held out.  This is a bounded negative result,
not a proof that no higher-degree multisection explains the labels.

## Complete rigid atlas

The five targets do not have rational `j`-preimages on the original published
equation.  Consequently, substituting their `074d9` parameters into the old
stored branch polynomials would test the wrong coordinate.  The certificate
instead transports all 39,120 norm-ten trace characters through the exact
integral MW-lattice isometry and reconstructs every residual quadratic on
the native `norm12-orbit-074d9` equation.

| curve | displayed quotient | splitting covers | exact span in displayed quotient |
|---:|---:|---:|---:|
| 351 | `Z^8` | 7 | 6 |
| 356 | `Z^12` | 2 | 2 |
| 376 | `Z^5` | 28 | 5 |
| 377 | `Z^6` | 7 | 5 |
| 385 | `Z^12` | 2 | 2 |

Every split point is constructed exactly and written in the stated free
quotient.  In the ordered bases `P18,...,P29`, the record-fibre classes are:

| curve | zero-based atlas index | cover | quotient coordinates |
|---:|---:|---|---|
| 356 | 20297 | `074d9-orbit-04b07` | `(0,0,0,0,0,0,0,1,0,0,0,0)` |
| 356 | 21395 | `074d9-orbit-11a44` | `(0,-1,-1,0,-1,0,0,1,0,1,0,0)` |
| 385 | 6986 | `074d9-orbit-11279` | `(0,0,0,1,0,0,0,0,0,0,0,0)` |
| 385 | 33745 | `074d9-orbit-080fa` | `(0,0,0,1,0,0,1,-1,0,0,0,0)` |

The two index sets are disjoint.  Reversing source and target tests the same
fixed-cover intersection, so this closes both requested directions for the
complete rigid inventory.

## Norm-eight transfer matrix

The target-blind canonical trace is priority rank one in the complete 63,925
norm-eight class ordering.  Its equation-basis word is

```text
(-1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0).
```

For each source, one member is fitted through each of the twelve displayed
quotient-basis directions and then frozen.  The exact transfer ranks are:

| source | `T` to 351 | `T` to 356 | `T` to 376 | `T` to 377 | `T` to 385 |
|---:|---:|---:|---:|---:|---:|
| 356 | 0 | 12 | 0 | 0 | 0 |
| 385 | 0 | 0 | 0 | 0 | 12 |

The diagonal entries are fitted controls.  Each off-diagonal entry has zero
splitting covers, not merely zero quotient span.  Extending the same protocol
to the first sixteen finite-pole traces gives 192 fitted covers per source;
the diagonal still spans 12 and all eight off-diagonal source/target entries
still have split count and transfer rank zero.

Curve 12 is a fibre of `norm12-orbit-11952`, not of `074d9`.  It therefore has
no `074d9` base parameter and is recorded as
`NOT_APPLICABLE_DIFFERENT_FIBRATION`, not as a zero.  The same issue prevents
a literal all-69 substitution matrix: cross-chart tests require an explicit
common-K3 birational/intersection transport first.

Because no tested frozen carrier meets both records, the proposed
base-Jacobian calculation has no input at this checkpoint.  No inference
about untested norm-eight traces is made.

## Late-label holdout

For each of `P18,...,P22`, the exact audit reconstructs from curves
351, 356, 376, and 377 and reserves curve 385 as an untouched holdout.  It
tests:

- the unique cubic polynomial `x(u)`;
- every rational `x=N/D` with `deg(N)+deg(D)=3`;
- quadratic and cubic inverse relations `u=u(x)`;
- the all-five quartic interpolant, followed by the exact Weierstrass square
  test for an integral section.

All five labels miss all declared holdout templates, and none of the five
quartics gives an integral section.  A general cubic relation with coefficients
affine-linear in `u` has six free coefficients and is not identifiable from
four reconstruction fibres, so it is deliberately excluded rather than
post-hoc fitted.

## Artifacts and replay

- [Complete rigid transfer certificate](../artifacts/generated-results/elkies-k3-r17-074d9-cross-fibre-bisection-transfer-v1.json)
- [Norm-eight transfer certificate](../artifacts/generated-results/elkies-k3-r17-074d9-norm8-cross-fibre-transfer-v1.json)
- [Late-point interpolation certificate](../artifacts/generated-results/elkies-k3-r17-074d9-late-point-interpolation-v1.json)

```bash
sage -python elkies-k3/scripts/certify_r17_074d9_cross_fibre_bisection_transfer.sage --check
sage -python elkies-k3/scripts/certify_r17_074d9_norm8_cross_fibre_transfer.sage --check
sage -python elkies-k3/scripts/audit_r17_074d9_late_point_interpolation.sage --check
```
