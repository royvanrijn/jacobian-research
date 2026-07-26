# Audit of the weighted tangent formula

## Status

**Mathematics: verified.**  **Rendered-formula readability: defective.**

An adversarial reader reconstructed two different maps from text extracted
from the rendered display of the former Proposition 2.1.  Neither
reconstruction had constant Jacobian.  The problem was not the proposition:
the extraction lost the outer divisions by `x^2` and `x`.

The canonical definition should therefore be the pair of intrinsic identities

\[
 BC=H'(W)+c\gamma,
 \qquad
 cAC^2=W\bigl(H'(W)+c\gamma\bigr)-H(W),
\]

with

\[
 u=1+xy,
 \qquad
 \gamma=1+a\,xy+b\,x^2z,
 \qquad
 W=u\gamma,
 \qquad
 C=x\gamma.
\]

Solving for the target coordinates gives unambiguously

\[
 B=\frac{H'(W)+c\gamma}{x\gamma},
 \qquad
 A=\frac{W(H'(W)+c\gamma)-H(W)}{c\,x^2\gamma^2}.
\]

Define

\[
 L_H(W)=\int_0^W\frac{V H''(V)}c\,dV.
\]

Since `H(0)=0`, integration by parts gives the coefficientwise identity

\[
 cL_H(W)=WH'(W)-H(W).
\]

Using `W=u\gamma`, the solved coordinates are exactly

\[
 \boxed{
 A=\frac{u+L_H(W)/\gamma^2}{x^2},
 \qquad
 B=\frac{c+H'(W)/\gamma}{x},
 \qquad
 C=x\gamma.}
\]

The outer quotients are essential.  Expressions such as
`u + (L_H(W)/gamma^2) x^2` or
`(u + L_H(W)/gamma^2) x^2` define different maps and are false readings of
the construction.

## Polynomiality in every degree

The admissibility conditions are

\[
 H(0)=H'(0)=H(1)=0,
 \qquad
 H'(1)=-c,
 \qquad
 \kappa=H''(1)/c\ne-2,
 \qquad
 a=-\frac{1+\kappa}{2+\kappa}.
\]

They clear the apparent denominators coefficientwise:

1. `H'(W)` is divisible by `W`, so `H'(W)/gamma` is polynomial after
   `W=u gamma`.  At `(xy,x^2z)=(0,0)`, its numerator satisfies
   `c+H'(1)=0`; hence it belongs to the ideal `(xy,x^2z)` and is divisible by
   `x` after substitution.
2. `L_H(W)` has a double zero at `W=0`, so `L_H(W)/gamma^2` is polynomial.
   The constant term of `u+L_H(W)/gamma^2` is
   `1+L_H(1)=0`.  Its coefficient linear in `xy` is
   \[
   1+\kappa+a(2+\kappa)=0.
   \]
   Therefore the numerator belongs to `((xy)^2,x^2z)` and is divisible by
   `x^2`.

This is an all-degree argument, not extrapolation from sample degrees.

## Constant Jacobian in every degree

Use the intermediate coordinates

\[
 s=H'(W)+c\gamma,
 \qquad
 t=W s-H(W).
\]

The four exact determinants are

\[
 \det\frac{\partial(x,xy,x^2z)}{\partial(x,y,z)}=x^3,
\]

\[
 \det\frac{\partial(W,\gamma,C)}
          {\partial(x,xy,x^2z)}=b\gamma^2,
\]

\[
 \det\frac{\partial(s,t,C)}{\partial(W,\gamma,C)}=-c^2\gamma,
\]

and

\[
 \det\frac{\partial(s,t,C)}{\partial(A,B,C)}=-cC^3.
\]

As `C=x\gamma`, the chain rule gives

\[
 \boxed{\det D(A,B,C)=bc.}
\]

The former paper used `b=1`, so the determinant is exactly `c`.

## Exact quartic regression

For

\[
 H_4(W)=W^4+W^3-2W^2,
 \qquad c=-3,
\]

one obtains

\[
 \kappa=-14/3,
 \qquad a=-11/8.
\]

Exact expansion verifies:

- both apparent quotients have denominator one;
- the coordinate degrees are `(12,11,4)`;
- the full three-variable Jacobian is identically `-3`;
- both intrinsic identities hold identically;
- `H_4(W)-BCW+cAC^2` vanishes identically on the source.

The executable certificate is

```bash
python scripts/verify_weighted_tangent_formula.py
```

It proves the all-degree factorized calculation without expanding `H` and
then independently expands the complete quartic instance.

## Effect on the former stable-inequivalence theorem

The corrected map is not a new map.  It is exactly the map implemented by
`jcsearch/weighted.py` and used by the weighted boundary and marked-root
certificates.  Consequently the unit-rank, Fitting, and boundary-contact
calculations are not changed by this correction.

There was nevertheless a verification-presentation gap: the former
`verify_common_arithmetic_fibers.py` checked the inverse pencil and quoted the
stable invariant values, but did not itself expand the displayed weighted map.
Other repository certificates did verify the correct map, including the
clean-room coefficientwise construction and Jacobian calculation.  The new
certificate makes the formula-to-map bridge explicit and directly guards the
outer parentheses that the PDF extraction lost.

Any future manuscript using this construction should present the two intrinsic
identities first and either omit the nested-fraction display or label it only
as their solved form.
