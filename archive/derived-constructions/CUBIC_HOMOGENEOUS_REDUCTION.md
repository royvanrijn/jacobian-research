# Explicit cubic-homogeneous counterexample

The announced degree-`(7,6,4)` counterexample has now been carried through an
explicit Bass--Connell--Wright/Yagzhev degree reduction using polynomial stable
equivalences and one Segre equivalence.  The result is a map

\[
K=I+H:\mathbb C^{95}\longrightarrow\mathbb C^{95},
\]

where every nonzero component term of `H` is homogeneous of degree three,
`det DK=1`, and three stored distinct rational points have exactly the same
image.  Thus `K` is an explicit cubic-homogeneous Keller counterexample.

The full map is stored sparsely in
`artifacts/generated-results/cubic_homogeneous_counterexample.json`.  It has 148 nonzero cubic
terms, so the sparse artifact is a substantially clearer exact specification
than 95 expanded coordinate formulas.

## Construction

First reorder and scale the original target coordinates to

\[
G=(F_3/2,F_2,F_1).
\]

Then `G(0)=0`, `DG(0)=I`, and `det DG=1`.

For a term `ab` of degree greater than three in one component, replace the
current map `(f_1,...,f_n)` by

\[
(f_1-(Y+a)(Z+b),f_2,\ldots,f_n,Y+a,Z+b).
\]

The lower-right Jacobian block in `(Y,Z)` is the identity.  Its Schur
complement is the old Jacobian matrix, so this operation preserves the
Jacobian determinant exactly.  On the graph slice `Y=-a`, `Z=-b`, it also
preserves every fiber and therefore transports the certified collision.

The deterministic implementation uses 22 such steps.  It produces a normalized
degree-three map in 47 variables,

\[
D(X)=X+Q(X)+C(X),
\]

where `Q` and `C` are quadratic and cubic homogeneous respectively.

Introduce one Segre variable `t` and 47 stable variables `Y`.  The final map is

\[
K(X,Y,t)=
\bigl(X-t^2Y+tQ(X),\;Y+C(X),\;t\bigr).
\]

Its nonlinear part is visibly cubic homogeneous.  Moreover

\[
\det D_X(X+tQ+t^2C)=\det DD(tX)=1,
\]

and the last displayed transformation is a composition with polynomial
automorphisms of determinant one.  Hence `det DK=1`.

Since `K=I+H`, `H` is cubic homogeneous, and `det(I+JH)=1`, scaling the source
shows `det(I+sJH)=1` identically in `s`.  Therefore the characteristic polynomial
of `JH` is a pure power and `JH` is nilpotent, as required in Yagzhev form.

If `P` is any transported collision point for `D`, then

\[
(P,-C(P),1)
\]

is a collision point for `K`.  The artifact stores the three resulting rational
points and their common 95-dimensional image.

## Reproduction

```bash
.venv/bin/python scripts/cubic_homogeneous_reduction.py
.venv/bin/python scripts/verify_cubic_homogeneous_counterexample.py
```

The first script reconstructs the map from the original three-dimensional
formula.  The second independently parses the sparse artifact, checks cubic
homogeneity, substitutes all three points into all 95 coordinates, validates the
22-step stable-equivalence certificate, and rechecks the initial determinant.

The degree-reduction formula follows the four-step implementation in L. Andrew
Campbell, *Reduction Theorems for the Strong Real Jacobian Conjecture* (2014),
which in turn implements the Bass--Connell--Wright/Yagzhev reduction.
The independent regeneration and primary-source hypothesis audit are recorded
in [C15_INDEPENDENT_AUDIT.md](../../extended-geometry/C15_INDEPENDENT_AUDIT.md).
