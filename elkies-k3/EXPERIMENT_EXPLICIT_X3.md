# Explicit discriminant-3 K3 anchor

We use two published Jacobian fibrations of the unique singular K3 surface X3
with transcendental lattice [[2,1],[1,2]].

Preferred anchor (Utsumi No.1):

    y^2 = x^3 + t^5 s^5 (t-s)^2

Fibers:

    II* at t=0
    II* at s=0
    IV  at t=s

so ADE = E8 + E8 + A2 and MW=0. The trivial lattice already has discriminant 3.

Independent cross-check (No.2):

    y^2 = x^3
          -3 t^2(t^6-16t^3s^3+16s^6)x
          +2 t^3(t^3-2s^3)(t^6+32t^3s^3-32s^6)

Fibers I12*, I3, 3 I1; ADE=D16+A2; MW torsion Z/2. Shioda-Tate gives
|disc(NS)|=(4*3)/2^2=3.

Run:

    sage elkies-k3/scripts/verify_disc3_k3_fibrations.sage

Then inspect the first-order deformation count:

    sage elkies-k3/scripts/probe_disc3_deformation_tangent.sage

The immediate research direction is to deform No.1 while preserving the rank-19
lattice polarization corresponding to X(6,79), rather than preserving all E8^2A2
fibers (which would keep us at/near the extremal locus).
