# The two-plus-one-plus-one quartic-denominator gate

## Status

This note proves `HC4NHM8`: the clean generic-corank-one partition

\[
P=L_0^2L_1L_2
\]

is empty for every residual-line incidence and every concurrency pattern of
the three denominator lines. The proof first compares the forced residual
root partitions on the normalization of the double line, then uses the same
immutable tangent coefficients as `HC4NHM7` on the only coincident-root
boundaries.

Replay the exact local certificates with

~~~bash
.venv/bin/python scripts/verify_hc4_two_one_one_quartic_denominator_gate.py
~~~

## 1. Local residual possibilities

Normalize the double denominator line to \(x=0\). The global cleared vector
is cubic, so `HC4NHM2--3` again force a constant kernel on this essential
line. A transverse kernel has zero determinant under \(x^4\)-divisibility;
normalize the tangent kernel to \(\partial_Z\), where \((Y,Z)\) are binary
coordinates on \(x=0\).

The tangent coefficient ladder has only the following exact possibilities:

| exact multiplicity of (x) | boundary branch | residual root partition on (x=0) |
|---:|---|---|
| (5) | \(\alpha=0,h_1=0,j_1\ne0\) | (4) |
| (4) | \(\alpha\ne0,j_1\ne0\) | (5) |
| (4) | \(\alpha=0,h_1\ne0\) | (4+1) |

For a residual line coincident with the double component, the target is

\[
\det C=L_0^5L_1^2L_2^2.
\tag{1.1}
\]

If \(L_1,L_2\) meet \(L_0\) in distinct points, the local residual type is
\(2+2\), incompatible with the table. If all three denominator lines are
concurrent, the type becomes \(4\), the sole possible exact-five row.

For residual line \(L_1\) or \(L_2\), exact multiplicity four has local type
\(3+2\) in general and type \(5\) on the concurrent boundary. For a fourth
distinct residual line, the generic local type is \(2+2+1\); its only
degenerations matching the table are the fully concurrent type \(5\) and
the type \(4+1\), where the two simple denominator lines meet \(L_0\) at one
point and the residual line at another.

Thus only the power and \(4+1\) concurrency boundaries require calculation.

## 2. Power boundary

On either the exact-five fourth-power row or the exact-four fifth-power row,
let \(Z\) be a tangent coordinate transverse to the common line pencil. The
complete determinant has the immutable coefficient

\[
\boxed{[x^7YZ]\det C=\frac1{18}j_1^3.}
\tag{2.1}
\]

Here \(j_1\ne0\) is exactly the condition that the displayed multiplicity is
exact. A determinant supported on the concurrent line pencil is independent
of \(Z\), so (2.1) is impossible. This removes every power degeneration.

## 3. The \(4+1\) boundary

The remaining possibility has exact multiplicity four on the
\(\alpha=0,h_1\ne0\) branch. Normalize the fourfold local root to \(Y=0\)
and the distinct residual root to \(Z=0\). As in `HC4NHM7`, the complete
next face contains

\[
\boxed{[x^5Y^2Z^2]\det C=-\frac{9h_1^4}{2\beta}.}
\tag{3.1}
\]

The target determinant has only one factor of the distinct residual line,
so it is at most linear in \(Z\). Equation (3.1) is nonzero and impossible.

> **Theorem `HC4NHM8` -- Two-plus-one-plus-one exclusion.** No clean
> generic-corank-one Hessian--Schur packet with quartic denominator partition
> \(2+1+1\) exists, including all line-concurrency strata.

The only clean quartic-denominator partition not yet treated is the
squarefree partition \(1+1+1+1\).
