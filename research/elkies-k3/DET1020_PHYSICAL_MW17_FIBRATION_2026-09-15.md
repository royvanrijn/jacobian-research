# Determinant1020: an explicit nef MW17 divisor marking

The [explicit rational K3 source](DET1020_EXPLICIT_RATIONAL_SOURCE_2026-09-15.md)
now has an exact primitive rootless U in its physical Néron–Severi basis,
with a full nef certificate. Its fiber, zero and seventeen saturated
rational section **divisor classes** are supplied in the
[certificate](../artifacts/generated-results/elkies-k3-det1020-equation-preflight-v1/physical-mw17-certificate.json).
The pencil and these sections exist over Q. Their rational functions and
Weierstrass coordinates remain to compute; this is not yet a usable
MW17 equation for specialization experiments.

## Physical coordinates and integral transport

Use the source basis

    f, O+f,
    Theta(I3 at lambda)_1,2,
    Theta(I3 at1)_1,2,
    Theta(I4)_1,2,3,
    Theta(I10)_1,...,9,
    v=O+3f-P.

Orient the I3 at1 so P meets component1; P meets I10 component5.
The other split chains can be oriented arbitrarily. The Gram is
U+F(-1), with F exactly the source certificate's central-frame Gram.
The source section coordinates are (2,1,0,...,0,-1), and O=(-1,1,0,...,0).
These verify P²=O²=-2 and P.O=1; the frame glue is included in v.

The retained unimodular row transport M satisfies

    M NS M^t = U + G(-1),    det(G)=1020,    min(G)>=4.

Exact enumeration finds no roots in G. After reflection into the source
nef chamber, its first row is

    D=(23,20,-7,-10,-6,-11,-5,-9,-10,
       -9,-16,-22,-25,-27,-26,-22,-17,-9,-1).

Thus D.f=20. The second row is its isotropic mate, and the remaining
seventeen rows form the entire negative frame, without a finite-index
loss. The transport was obtained from a bounded root-masking neighbor of a
degree12 A1 parent, minimizing source incidence within each parity coset.
Abstract indefinite uniqueness is no longer needed to locate this U.

## Complete nef check

A finite list of source sections alone would not prove nefness. Instead
use the explicit ample class

    H=20O+42f-sum c_i Theta_i,
    c=(2,2; 2,2; 3,4,3; 9,16,21,24,25,24,21,16,9).

In the displayed basis H has coordinates(22,20,-c,0), and H²=514.
It meets O in2, each nonidentity fiber component in2, and each identity
component in20-2(n-1)>0. For any other irreducible horizontal curve C
of old-fiber degree d, O.C>=0 and the component intersections in each
I_n fiber sum to d. The largest c_i in the four chains sum to33, hence

    H.C >= (42-33)d =9d >0.

Nakai–Moishezon therefore makes H ample. If an irreducible curve C has
D.C<0, then C²=-2. The root reflection preserves the positive cone, so
0<H.(D+(D.C)C), and consequently H.C<H.D. After checking O and every
vertical component, it suffices to test horizontal degrees

    1 <= d <= floor((H.D-1)/9).

For D=(a,b,w) and a degree-d root C=(A,d,v), integrality requires
A=(v.F.v-2)/(2d). Since D²=0,

    D.C = b/(2d) * (v-d*w/b).F.(v-d*w/b) - b/d.

Thus a negative root must lie in the exact ellipsoid
(v-d*w/b).F.(v-d*w/b)<2. All such integer vectors are enumerated by the
positive definite augmented form

    Q(v,k)=b²(v-k*d*w/b).F.(v-k*d*w/b)+2b²k²,
    Q <=4b²-1,    k=+/-1.

The checker applies an integral LLL basis change, enumerates this finite
ellipsoid exactly and filters for integral A. Vectors with k=0 are
irrelevant, and |k|>=2 cannot meet the bound. This checks all possible
negative roots, including horizontal multisections.

The current degree20 fiber requires no reflection: H.D=411, and the
complete final sweep checks all45 required degrees without finding a
negative root. Together with O and all fiber components, this proves
full nefness. An earlier degree183 certificate is retained separately:
it required degree4 and22 horizontal reflections and a507-degree final
sweep, demonstrating why checking source sections alone is insufficient.

## Actual arithmetic pencil and saturated sections

D is primitive, nef and square zero, so its complete linear system is
an elliptic pencil. The source certificate represents every NS class by
an actual Q-divisor; the pencil is therefore over Q with base P1_Q.
Its frame G has no roots, hence every geometric fiber is irreducible.
Put Z=M_1-D, where row indices start at0. Then Z²=-2 and Z.D=1.
Riemann–Roch and nefness make Z effective. Its degree-one horizontal
part is a section; any remaining vertical part is a multiple of D.
The self-intersection-2 forces that multiple to vanish, so Z itself is
the zero section.

For each frame row M_(i+2), put

    S_i=(G_ii/2-1)D+M_1+M_(i+2),    i=0,...,16.

Each class has square-2 and fiber degree1, so the same argument gives
an actual rational section. Their Shioda projections are exactly the
seventeen frame basis rows. Thus their height Gram is G, and they are
a saturated basis of MW=Z^17 over Q(t), with zero torsion. No section
coordinates or rational pencil function are inferred from these classes.

## Replay and retained limits

    sage -python research/elkies-k3/scripts/certify_det1020_physical_mw17.py --check

The checker verifies the full integral transport, rootlessness, the
complete final nef sweep, and all section-class and height identities.
K3 Riemann–Roch, the nef-pencil theorem and the ample bound above are
written geometric inputs. Independent implementation, formal verification,
external review and novelty are unclaimed.

The packet retains the initial192-candidate isotropic miss, the256-candidate
regular-root miss, and the successful253rd masked2-neighbor candidate.
All have explicit seeds, caps and terminal checkpoints. Complete fixed
MW-coordinate windows at source degrees10 and11 found no rootless member.
The degree12 window stopped at its120-second cap after680110 raw choices
and22251 transports; it is incomplete and is not restarted. Its retained
A1 parent has source degree12. A512-candidate masked2-neighbor cost pilot
from this parent found the degree20 rootless fiber above. Coordinate
descent within each parity coset minimizes an incidence quadratic only
heuristically; no minimal-incidence or window-completeness theorem is
claimed. The earlier degree183 certificate remains replayable with
`certify_det1020_physical_mw17_degree183.py --check`.
The overall explicit-equation goal remains OPEN.
