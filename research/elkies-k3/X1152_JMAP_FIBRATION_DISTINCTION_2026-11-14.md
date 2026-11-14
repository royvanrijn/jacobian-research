# The Hill and ICARM4239 fibrations have different geometric j-maps

The exact generic equations imported in
[the November transfer](RECENT_PAPERS_X1152_TRANSFER_2026-11-13.md) are not
isomorphic as elliptic fibrations, even after a base automorphism over
Qbar. This distinguishes the fibrations only. Whether their underlying K3
surfaces are isomorphic remains UNKNOWN.

## Intrinsic finite critical-value invariant

Write the short model as y²=x³+A(t)x+B(t), and put

    Delta = 4 A³ + 27 B²,
    j = 6912 A³ / Delta,
    C = 3 A' B - 2 A B'.

Then j' = 6912·27 A² B C / Delta². For both imported models, reduction at197
certifies deg(A,B,Delta,C)=(8,12,24,18), squarefreeness of each polynomial,
gcd(A,B)=1 and gcd(C,A B Delta)=1. These conditions lift to characteristic
zero. The cancellation making deg(C)=18 also proves that infinity is
unramified: j' has order2 there in the t-coordinate.

Consequently the18 roots of C are exactly the critical points lying away
from the special critical values0,1728,infinity. Each is simple. Their
j-value multiset is intrinsic under any precomposition by PGL2(Qbar).
Its characteristic polynomial, and in particular its trace, is rational.

The good-reduction conditions make Delta invertible in F197[t]/C and
make the characteristic-zero trace integral at197. Multiplication by
6912 A³/Delta in this18-dimensional quotient algebra gives:

| Generic model | Critical-value trace mod197 | Smooth K3 F197 points |
|---|---:|---:|
| Hill imported parent |56|42381|
| ICARM4239 generic-only parent |188|42381|

The unequal trace reductions prove unequal rational traces in characteristic
zero. Thus no PGL2(Qbar) base change identifies the two j-maps. This
obstruction also applies before choosing a zero-section marking.

Equal surface point counts at this one prime give neither a surface
identification nor a surface distinction. No Picard bound, full generic
record-parent basis, lattice marking or transport between fibrations follows.

## Certificates and assurance

The bounded [Sage producer](scripts/certify_1152_fibration_distinction.py)
uses a fixed prime,10 CPU seconds and4GiB address space. Its
[frozen certificate](../artifacts/generated-results/x1152-fibration-distinction-20261114-v1.json)
binds both generic sources by SHA-256. No exceptional record parameter is used.

The [standard-library replay](scripts/verify_1152_fibration_distinction.py)
independently reconstructs the finite polynomial algebras, squarefreeness and
coprimality gates, traces and all characteristic-polynomial coefficients.
Newton identities are valid because197>18. It also recomputes the surface
point counts by direct finite-field enumeration. The
[v2 independent receipt](../artifacts/generated-results/x1152-fibration-distinction-independent-20261114-v2.json)
supersedes the arithmetic coverage of the retained v1 receipt, which checked
only the characteristic polynomial's trace coefficient.

    python3 research/elkies-k3/scripts/verify_1152_fibration_distinction.py

The two equations may still describe different fibrations on one K3 surface.
Any future surface comparison must preserve that distinction.
