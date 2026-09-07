# Q80 orbit 1222 — exact characteristic-zero Weierstrass model

Status: **PASS_EXACT_MODEL**

The exact critical-value computation determines `mu`, while certified monic polynomials `C` (degree 8) and `S` (degree 12) determine the j-map.  No square-root extension is needed: choose

    c4 = mu*C
    c6 = mu^2*S
    Delta = mu^3*(C^3-mu*S^2)/1728.

Then `c4^3-c6^2=1728*Delta` identically, so the short model is

    y^2 = x^3 + A(V)*x + B(V)
    A = -c4/48
    B = -c6/864.

The model is defined over `QQ(sqrt(-3))`.

Exact fiber-support polynomials were extracted from the gcd tower, without full discriminant factorization:

    Delta ~ P7^7 * P2^2 * P1

with degrees `(deg P7, deg P2, deg P1)=(2, 3, 4)`, hence `2 I7 + 3 I2 + 4 I1` over the algebraic closure.

Both independent p=73 Galois j-maps validate by exact cross-multiplication.
