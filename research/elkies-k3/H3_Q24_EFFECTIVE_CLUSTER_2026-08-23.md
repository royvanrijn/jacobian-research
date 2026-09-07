# H3 q24 effective-cluster frontier — 2026-08-23

The q24 `D13/MW4 -> D12/MW5` equation route has now passed three increasingly geometric gates over `GF(100003)`:

1. The global degree-two preflight has ambient dimension 58, smooth-collision rank 48, and post-collision dimension 10.
2. The actual `I9*` germ resolves through 12 point blow-up centres but 13 geometric exceptional components because centre `C10` has rank-two tangent cone and splits. The resulting root graph is D13 with arm lengths `1,1,10`.
3. The deterministic lattice D13 graph matches the geometric graph up to the spinor-arm swap. The earlier coordinate-jet scan checks all 45 coordinate two-planes in the ten-dimensional ordinary-jet space; none produces genus one. Therefore the missing codimension-eight condition is genuinely an infinitely-near/resolved-component condition, not ordinary `u^8` vanishing or any coordinate two-jet complement.

## Effective sign

The independent D13 chamber certificate fixes the sign that graph matching alone cannot see: the effective D13 simple roots are the **negatives** of the first thirteen root-adapted frame vectors.

Thus the equation-frame vertical coefficients

```text
vr = 9,7,1,2,3,4,5,6,7,5,8,4,8
```

mean that, with the global `-7F` twist represented at infinity, the finite `I9*` local divisor contains `-vr[i] E_i`. A section therefore needs resolved divisorial valuation at least `vr[i]` along the corresponding effective component.

Converting these thresholds through the actual blow-up chronology gives the common nonzero infinitely-near point cluster

```text
C01 : +2
C02 : +2
C04 : +2
C06 : +3
```

for both spinor orientations. All other component lower bounds follow from the valuations already forced at earlier centres; the `C10a/C10b` asymmetry is dominated before the split and therefore does not create a branch condition at this stage.

The exact derivation is implemented in

```text
scripts/derive_h92_q24_i9star_effective_cluster_modp.sage
```

and the next equation probe is

```text
scripts/probe_h92_q24_d12_resolved_cluster_rr_modp.sage
```

which applies the four orders successively in the actual strict-transform surface local rings and then compiles any resulting two-dimensional kernel by the binary-quartic route.

## Claim boundary

Neither the dimension `10` nor the desired final dimension `2` is by itself a proof. Promotion to an exact modular q24 RR certificate requires the cluster probe to produce a codimension-eight kernel and a degree-three/four squarefree chord radicand with D12 fibre data. A later independent replay must then verify all redundant component valuations before characteristic-zero lifting.
