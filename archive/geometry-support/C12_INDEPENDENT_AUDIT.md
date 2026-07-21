# Clean-room strengthening of the C12 local singularity

This audit isolates the algebraic heart of
[DEGREE12_LOCAL_SINGULARITY.md](../../extended-geometry/DEGREE12_LOCAL_SINGULARITY.md) and checks it
without SymPy or project algebra helpers.

For one exchanged sixfold root, depress the cubic and write

\[
Q=z^3+uz+v,\qquad T=z^2+az+b.
\]

Coefficient comparison in `Q^2-T^3` gives

\[
\begin{aligned}
c_5&=-3a,\\
c_4&=2u-3b-3a^2,\\
c_3&=2v-a^3-6ab,\\
c_2&=u^2-3a^2b-3b^2,\\
c_1&=2uv-3ab^2,\\
c_0&=v^2-b^3.
\end{aligned}
\]

The first three equations put `a`, `2u-3b`, and `v` in the coefficient
ideal.  Moreover,

\[
4c_2-(2u-3b)(2u+3b)=-12a^2b-3b^2,
\]

so `b^2` also belongs to it.  Conversely, setting
`a=v=0`, `u=3b/2`, and `b^2=0` kills all six coefficients.  Therefore

\[
(c_5,c_4,c_3,c_2,c_1,c_0)=(a,v,2u-3b,b^2).        \tag{1}
\]

This replaces the original black-box Groebner-basis step by an elementary
ideal proof.  One block is a dual number.  The two separated roots give two
independent blocks, hence transverse algebra

\[
k[\epsilon,\eta]/(\epsilon^2,\eta^2)
\]

with basis `1,epsilon,eta,epsilon*eta` and length four.

The clean-room checker also verifies that the two quadratic relations from
the common-tangent quotient restrict to

\[
(5x_1-6x_2)^2,qquad 9(144x_0+125x_1)^2
\]

on successive independent linear quotients.  This reproduces the upper
length bound independently from the lower dual-block construction.

Run

```bash
python3 scripts/audit_c12_independent.py
```

The script uses only `fractions`.  The full tangent-rank matrices are still
most conveniently reproduced by the original exact checker, so an unrelated
CAS reproduction remains desirable before labeling the entire theorem
independently verified.
