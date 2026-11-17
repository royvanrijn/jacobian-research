# Quadratic trace criterion and an exact rank-seven dimension

Further structural extension: the [full receiver and unified theory map](../universal-keller-receiver/THEORY_MAP.md) proves the proposed root-line quotient and its inertia, effective polynomial bundle gluing, fixed-space restrictions and exact affine-source and Laurent-profile landing rules.

The [proof](PROOF.md) constructs determinant-one maps D_N on A^(N-3), N>=7. A rank-N finite etale algebra over an infinite field of characteristic other than two is a complete fibre exactly when its trace quadratic form on ker Tr has a zero outside its radical. All D_N have degree N, monodromy S_N and absolute stable atomicity.

A binary-block trace-form test characterizes field-versality when N is invertible and remains sufficient in modular ranks. In characteristic zero, a base containing sqrt(-14) passes the rank-seven test. Duncan's external lower bound then gives the exact value ktdim_k(S7)=4. Fixed-field arithmetic criteria alone do not prove field-versality.

- UADOUBLE1: exact full-fibre trace-quadric criterion, all N>=7.
- UABINARY1: the binary-block field-versality criterion and modular sufficient test.
- UAK7MIN1: exact target dimension four in rank seven over the specified bases.
- UADOUBLEMOD1: canonical quadric lift spaces and N-6 stable-class parameters.

[Literature and gap map](LITERATURE_AND_GAPS.md) · [checker](verify.py) · [receipt](receipt.json) · [sources](sources.json) · [authority](../../MATH_STATUS.json).

    .venv/bin/python research/jacobian-followup/universal-keller-quadrics/verify.py

Written proofs with exact controls, not formalization or independent whole-proof review. Global arbitrary-field minima, coordinate-degree minima and presentation-free descent remain open.
