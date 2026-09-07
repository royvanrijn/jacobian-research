# Degree-six sparse equivariant search

> **Archived bounded search.** The exact exclusions remain scoped to the box
> below and are not part of the active theorem chain.

## Result

No new Keller map is claimed here.

For the foundational source weights `(1,-1,-2)` and target weights
`(-2,-1,1)`, the complete total-degree-at-most-six support has component
sizes

\[
8,\qquad 7,\qquad 5.
\]

Write `v=xy`, `s=x^2z` and

\[
F=(x^{-2}A(v,s),x^{-1}B(v,s),xC(v,s)).
\]

After normalizing the linear part to `(z,y,x)`, the Keller equation is

\[
\det\begin{pmatrix}
-2A&A_v&A_s\\
-B&B_v&B_s\\
C&C_v&C_s
\end{pmatrix}=-1.
\]

The points `(1,1,1)` and `(-1,-1,1)` have the same image exactly when

\[
B(1,1)=C(1,1)=0.
\]

The complete degree-six subbox which is linear in `z` has 14 nonlinear
coefficients and 21 determinant/collision equations.  Its reduced
Gröbner basis over `QQ`, computed by `msolve`, is

\[
\boxed{[1]}.
\]

Thus no map in that exact subbox supplies the desired collision.

The three extra coefficients in the full degree-six box are the `s^2` and
`vs^2` coefficients of `A` and the `s^2` coefficient of `B`.  Their three
principal opens cover the complement of the `z`-linear closed box.  Adding
an inverse variable for each coefficient in turn gives, over `QQ`,

\[
\boxed{[1],\qquad[1],\qquad[1].}
\]

Together with the closed-box calculation, this excludes the complete
degree-six equivariant support for a collision on the normalized nonzero
order-two orbit.  In particular it also excludes every ten-monomial
subsupport inside this box.

This does not cover a symmetry-breaking map, nor an equivariant collision
whose invariant point lies on `vs=0`.

A separate seeded modular search uses the fact that, once `B,C` are fixed,
every determinant equation is linear in the seven free coefficients of
`A`.  The first 100,000 full-support samples over `F_65521` produced no
consistent linear system.  This search is only a regression alongside the
exact principal-open proof.

An additional exact pass through every single elementary source or target
shear of degree at most three found no parameter over the algebraic closure
which lowers the foundational map to degree six or below, and no rational
specialization with fewer than its sixteen expanded nonconstant terms.

## Reproduction

Run

```bash
.venv/bin/python archive/tooling/search_degree6_sparse_equivariant.py
.venv/bin/python archive/tooling/search_degree6_sparse_equivariant.py --full-exact
```

All four `[1]` calculations are characteristic zero.  The first command
runs the fast closed chart, the shear pass, and the seeded regression.  The
second also runs the three slower principal opens; on the reference machine
their times were approximately 199, 0.3, and 18 seconds.  The seeded pass is
clearly labelled `SEARCH` in the output.
