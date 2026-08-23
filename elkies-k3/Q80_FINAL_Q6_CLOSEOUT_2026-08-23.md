# Q80 final q6 characteristic-zero closeout — 2026-08-23

## Status

The Q80 low-q terminal neighbour is now closed at both levels that must be kept distinct:

```text
generic lattice:
A1/MW16 --q6(2,3)--> rootless/MW17, det(NS)=948

CM24 characteristic-zero specialization:
A1+A2+A3+A4+A5/MW3 --q6--> 4A2+A3+A5/MW2
```

The specialized exact child has

```text
finite fibres = 4 I3 + I4 + I6 + 2 I1
infinity      = smooth
root rank     = 16
root count    = 66
root det      = 1944
Euler number  = 24
MW rank       = 2
```

The generic endpoint remains the separately certified rootless `MW17` frame. The CM24 child is not the generic rootless equation; it is the exact characteristic-zero specialization shadow of the selected generic divisor.

## What solved the final marking problem

The final generic q6 divisor was already fixed by the lattice search. Replaying the retained neighbour chain with its unimodular basis transports recovered the exact candidate1 horizontal class

```text
P_candidate1 =
(5,1,-44718,-282065,63356,564493,-98198,249323,239104,-1054,
 -22328,-389456,-231271,-641746,-570362,-123785,227276,-186445,89497)
```

with

```text
P^2 = -2
P.F = 1
P.O = 4
D_reduced = O + P - F.
```

Transporting this curve backwards produced very large multisection coordinates in earlier fibrations. That confirmed that the useful small labels `P1,P2,P3` arise only after CM24 specialization and re-chambering; they are not small generic MW coordinates along the whole suffix.

The direct final-section lift at `p=73` was singular. Instead of solving the difficult three-node `P2-P3` section directly, the successful construction recovered easier exact `P.O=0` sections meeting four or five reducible nodes and identified the correct pair by reduction modulo `73`. Their exact elliptic-curve difference is the desired final horizontal.

Thus the final characteristic-zero horizontal is obtained as an exact MW-basis difference, while the historical modular point is used only for identification/regression:

```text
H = P2 - P3
hits = I3, I4, I6
P.O = 0
height = 1.
```

## Exact RR closure

With that exact horizontal, the final resolved Riemann--Roch problem is completely exact over `QQ(sqrt(-3))`:

```text
ambient dimension = 4
whole-A4 rank      = 1
connected-A5 rank  = 1
condition rank     = 2
kernel dimension   = 2
h0(D)              = 2.
```

The whole-A4 row at the exact `I5` fibre reduces to the transported historical pinned row. The connected-A5 quotient at the exact `I6` fibre is recovered by an exact leading-jet toric calculation and reduces to the transported `+/-4` quotient line.

The A5 implementation deliberately uses only the leading local residue needed for the quotient line. It does not divide full Laurent series or expand every rational-function coefficient row; this is mathematically equivalent to the earlier toric calculation and avoids the large symbolic slowdown caused by the exact candidate1 coefficients.

## Final exact child

Compiling the certified two-dimensional pencil gives an exact degree-four binary quartic. Its Jacobian finite-minimizes and classifies as

```text
4 I3 + I4 + I6 + 2 I1,
infinity smooth,
root lattice 4A2+A3+A5,
root rank 16,
root determinant 1944,
MW rank 2.
```

This is certified by

```text
elkies-k3/data/fibrations/q80-final-q6-char0/Q80_CHAR0_FINAL_Q6_CERTIFICATE.md
elkies-k3/data/fibrations/q80-final-q6-char0/q80_char0_final_q6_child.sage
```

and reproduced by

```bash
sage elkies-k3/scripts/trace_q80_candidate1_marked_transport.sage
sage elkies-k3/scripts/recover_q80_final_q6_via_basis_sections.sage
sage elkies-k3/scripts/certify_q80_final_q6_char0_rr_from_basis.sage
sage elkies-k3/scripts/compile_q80_final_q6_char0_child.sage
```

The scripts write the generated intermediate certificates under

```text
artifacts/generated-results/q80-candidate1-marked-transport.json
artifacts/generated-results/q80-final-q6-char0-horizontal-via-basis-sections.json
artifacts/generated-results/q80-final-q6-char0-rr.json
artifacts/generated-results/q80-final-q6-char0-child.json
```

These generated artifacts need not be committed for the proof path to be reproducible; the scripts and pinned equation/model certificates are the repository sources.

## Superseded final-q6 approaches

The following remain useful diagnostics but are no longer the construction path:

- direct two-parameter characteristic-zero resultant/Groebner elimination for the three-node final section;
- digit-by-digit `73`-adic lifting in the singular residue disk;
- the local-`73` singularity probes used to explain why that lift was slow;
- treating the CM24 polynomial section `P2-P3` as the generic definition of the divisor.

The bad `GF(73)` point was genuinely non-transverse in the redundant two-parameter section coordinates. That explains the failure mode; it does not weaken the exact characteristic-zero construction above.

## Durable lesson

For this suffix:

> **Search provenance first; preserve or reconstruct easy markings second; use the CM24 section only as a regression certificate.**

The final q6 reconstruction should now be treated as closed. Future Q80 work should build on the exact child/certificates above rather than reopening the singular final-section lift.