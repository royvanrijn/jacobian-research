# All-degree HC4 graph obstructions under dual automorphisms

For the scaled Meng--Yang potential, allow **every** polynomial dual
automorphism (p,q) over K[x,y] and every polynomial r=R(x,y,P,Q).
A nonzero constant Hessian determinant of the actual pullback forces
its gradient to be a polynomial automorphism. Degrees, coefficients,
y-dependent traces and alternating word lengths are unrestricted.

The [automorphic-trace proof](AUTOMORPHIC_TRACE.md) (`HC4MYGV0/1/2`)
actually needs only an automorphic x=0 dual trace over K[y]. It exhausts
all affine fibre planes, including those contained in x=0. Transverse
planes yield either a triangular endpoint or a locally nilpotent
Hamiltonian with a constant-Jacobian mate. Global invertibility uses
Tang's reproduced structural alternative and the classical Rentschler
criterion, with a [retained primary source](sources/automorphic-trace-import.json).

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-triangular-graphs/verify_automorphic_trace.py
```

[Actual replay receipt](receipt-automorphic-trace.json). The finite
identities and controls do not replace the written all-plane proof or
independently certify the imported theorems.

Retained results have additional or different strength:

- [First slice](FIRST_SLICE_INVERTIBILITY.md), `HC4MYGF1/2`: arbitrary x-divisible corrections to a base-independent **Keller** trace pair, without assuming a trace inverse. It also supplies an explicit unit-Hessian embedded control with an affine plane.
- [Plane-Keller encoding](FIRST_SLICE_INVERTIBILITY.md#6-sharp-trace-boundary-the-larger-substitution-class-contains-jc2), `HC4MYGF3`: the larger arbitrary-substitution class with merely Keller y-dependent traces contains every plane Keller cotangent lift. The encoding has nonconstant full dual Jacobian and is not asserted to be embedded.
- [Matching two divisors](TWO_DIVISOR_EXTENSION.md), `HC4MYGD1/2`: a stronger no-plane conclusion on its stated subclass.
- [Earlier triangular and base-independent theorem](PROOF.md), `HC4MYGT1/2`: retained no-plane arguments and exact marked-covector formulas.

Nonlinear base changes, other eliminated-coordinate charts, embeddings
with uncovered first traces, and genuinely mixed/coisotropic mechanisms
remain open. No preservation of a Hessian determinant under nonlinear
coordinate change is assumed. No general HC4 or JC2 solution, formal
verification, independent external review or novelty priority is claimed.
Mathematical status remains authoritative in [MATH_STATUS.json](../../MATH_STATUS.json).
