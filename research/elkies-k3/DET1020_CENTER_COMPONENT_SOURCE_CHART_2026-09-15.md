# Determinant1020: a rational chart using the central I10 component

**Subsequent construction:** the [explicit source certificate](DET1020_EXPLICIT_RATIONAL_SOURCE_2026-09-15.md)
now supplies the full split source equation and Picard19. The unresolved
source-equation statements below describe the earlier scope of this note;
physical MW17 transport and section compilation remain open.

The admitted NS has an alternative source frame A2+A2+A3+A9, with a
saturated MW generator of height17/6 and P.O=1. Its component indices
are0,1,0,5. This leads to a rational six-parameter K3 chart in which I10
and one I3 are built in. The other I3 (or IV) and I4 remain unconstructed;
no member has yet been certified as the determinant1020 source.

## Exact frame and arithmetic scope

Adjoin a norm6 vector to the indicated Cartan root lattice, pairing1 with
the first root in the second A2 and the fifth root in A9, and0 elsewhere.
Its projected height is6-2/3-5/2=17/6. The frame determinant is1020.
Since17/6>2, all roots lie in the primitive rank16 root lattice; the exact
check finds114 roots. The cyclic discriminant generator has norm8837/1020,
matched to the admitted T generator-403/1020 by the unit59 modulo2Z.
Indefinite uniqueness and rational nef/divisor descent, as in the
[MW17 proof](DET1020_ARITHMETIC_MW17_EXISTENCE_2026-09-15.md), identify
this as an actual arithmetic source fibration on the admitted K3.
The primitive quotient and height formula give torsion zero and P.O=1.
A9 and A3 mean I10 and I4; each A2 can mean I3 or IV.

The middle I10 component changes the equation cost. In coordinates
Y²=X³+A X²+B X+C² with P=(0,C), conditions t^5|B and t^5|C give

    D=t^10*[A²(b²-4Ac²)+t^5*(18Abc²-4b³)-27t^10*c^4],
    b=B/t^5, c=C/t^5.

For nonzero leading bracket this is I10 and the section meets component5.
There is no long high-order cancellation as in the earlier
[I9 component1 chart](DET1020_LOW_POLE_SOURCE_FRAME_2026-09-15.md).
The tradeoff is the single pole of P, which must be handled globally.

## Rational construction

Place I10 at0, the nonidentity I3 at1, and the single pole of P at infinity.
The pole need not lie on a singular fiber. Use the minimal short equation

    y²=x³+f(t)x+g(t),       P=(A(t)/3,C(t)),
    deg f<=8, deg g<=12,   deg A=6, deg C=9,
    leading(A)=3, leading(C)=1.

These leading terms give P.O=1 at infinity. Start with six parameters
(d,e,a2,a3,a4,a5) and

    A=d²+(e²-d²-a2-a3-a4-a5-3)t+a2*t²+a3*t³+a4*t⁴+a5*t⁵+3*t⁶.

Thus A(0)=d² and A(1)=e² enforce the two split tangents. Initially put
f=-truncation_through_degree4(A²/3)+f5*t⁵+f6*t⁶+f8*t⁸ and C=t⁹.
Successively cancel coefficients17,16,15,14 in C²-A³/27-Af/3 by solving
for C8,C7,C6,C5; each step is linear with coefficient2. Set f7 equal to
the remaining coefficient13, so adding f7*t⁷ to f makes

    g=C²-A³/27-Af/3

have degree at most12. Now C(1)=0 determines f8 linearly, B(1)=0 for
B=A²/3+f determines f5 linearly, and

    B'(1)+2e*C'(1)=0

determines f6 linearly after the previous substitutions. The respective
coefficients are1/2,1,1. Thus this reduction introduces no unresolved
quadratic cover. All coefficients are explicit rational polynomials in
the six parameters, retained as sparse exact data.

The e=0 branch is retained: it has orders(2,2,4) for(f,g,D) at1,
with leading g coefficient C'(1)², hence generically split IV. The
checker verifies this branch as well as the section identity, degrees, split tangent values,
generic discriminant orders10 and3, and degree24. A nondegenerate test
(d,e,a2,a3,a4,a5)=(1,2,0,0,0,0) gives

    f=-2t^8-33t^6+30t^5-1/3,
    g=32t^12-30t^11+t^10+(2/3)t^8+11t^6-10t^5+2/27,
    P=(t^6+1/3, t^9-t^5).

Its other eleven fibers are I1, verified by the squarefree residual
discriminant. This only proves the partial-family open set is nonempty;
this example is not claimed to have determinant1020 or MW17.

## Remaining gate and replay

The residual discriminant after dividing by t^10(t-1)^3 has degree11.
Require an additional identity-component I3 or IV and an identity-component
I4 at distinct rational locations. Their exact type, tangent splitting,
section intersections, saturation and Picard19 all require certificates.
The pole can lie on a singular fiber, so imposing smooth infinity is an
open-chart choice, not a global restriction. The rational source parameter
and physical marking remain unknown.

```sh
sage -python research/elkies-k3/scripts/certify_det1020_center_source_chart.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/center-source-certificate.json)
checks both the frame and the partial chart. Independent implementation,
formal verification, external review and literature novelty are unclaimed.
The [packet](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/README.md)
retains earlier I9 failures and the bounded center-chart checks. Complete
F7 and F11 smooth-infinity windows found no admissible point, and a
separate F7 multiplicative pole-boundary window and the e=0 split-IV
window also missed. None is a
rational-point obstruction; no larger finite-field sweep is prescribed.
