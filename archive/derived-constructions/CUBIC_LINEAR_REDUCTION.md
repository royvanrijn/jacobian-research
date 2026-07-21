# Explicit cubic-linear counterexample

The 95-dimensional cubic-homogeneous counterexample `f=x+H` has been paired,
in the sense of Gorni--Zampieri, with an explicit cubic-linear map

\[
G(X)=X-(AX)^{*3}:\mathbb C^{451}\longrightarrow\mathbb C^{451}.
\]

Here `*3` denotes coordinatewise cubing.  The matrix `A` is rational, has rank
95, `det DG=1`, and three stored distinct rational points have exactly the same
image.  Consequently `G` is an explicit Druzkowski-form Keller counterexample.

The sparse matrix, pairing certificate, collision points, and common image are
stored in `artifacts/generated-results/cubic_linear_counterexample.json`.

## Pairing construction

The polarization identities

\[
ab^2=\frac{(a+b)^3+(a-b)^3-2a^3}{6},
\]

\[
abc=\frac{(a+b+c)^3+(a-b-c)^3-(a+b-c)^3-(a-b+c)^3}{24}
\]

write every cubic monomial as a linear combination of cubes of linear forms.
After identical linear forms are combined, the 148 terms of `H` require 415
distinct cubes.  This gives matrices `B_0,D_0` with

\[
H(x)=-B_0(D_0x)^{*3}.
\]

The matrix `D_0` already has full column rank 95.  The first implementation put

\[
B=(B_0\mid I_{95}),\qquad D=\binom{D_0}{0},\qquad A=DB,
\]

which gives dimension 510.  Exact row reduction, however, gives
`rank(B_0)=59`.  The audited construction appends only a 36-column coordinate
complement and computes a rational right inverse `C`.  This gives dimension
`415+36=451`, and still

\[
BC=I_{95},\qquad AC=D,
\]

and, because `D` is injective,

\[
\ker A=\ker B.
\]

Moreover

\[
f(x)=BG(Cx).
\]

The Sylvester determinant identity gives

\[
\det DG(X)=\det Df(BX)=1.
\]

Thus the determinant certificate is an exact matrix identity, not a collection
of numerical specializations.

## Collision transport

The points `C p_i` do not necessarily have identical images under `G`; their
image differences lie in `ker B=ker A`.  Since

\[
G(X+k)=G(X)+k\qquad(k\in\ker A),
\]

each point can be translated by the corresponding kernel difference to give a
literal common image.  The artifact stores the three resulting rational points,
and the verifier substitutes them into all 451 coordinates exactly.

## Reproduction

```bash
.venv/bin/python scripts/cubic_linear_reduction.py
.venv/bin/python scripts/verify_cubic_linear_counterexample.py
```

The construction follows Proposition 2.1 of Gorni and Zampieri, *On
cubic-linear polynomial mappings* (2013), which makes Druzkowski's reduction
explicit through the pairing matrices `B` and `C`.

Dimension 451 is minimal for this fixed polarization matrix `B_0`; no
minimality among all cube decompositions is claimed.  See
[C15_INDEPENDENT_AUDIT.md](../../experimental/geometry/C15_INDEPENDENT_AUDIT.md).
