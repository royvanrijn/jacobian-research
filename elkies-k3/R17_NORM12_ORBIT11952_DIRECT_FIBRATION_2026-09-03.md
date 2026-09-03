# Direct `norm12-orbit-11952` fibration compilation (2026-09-03)

<!-- status-consumer: EC-K3-R17-NORM12-11952-DIRECT-Q80-EQUATION 077c6409d76cbe63 -->

## Result

The alternate-Q80 rootless fibration is now compiled directly on Elkies's
published R17 equation.  The input divisor is

```text
D = (3,2,w),
w = (-1,-1,3,0,0,1,-1,-1,3,1,1,0,2,-1,1,-2,-2),
```

and the exact identity in `NS=U+R17(-1)` is

```text
D = O_old + P_w - F_old.
```

Thus `D.F_old=2`, `D.O_old=1`, and the published zero is also the zero of the
new pencil.  The direct equation, maps, fibre certificate, alternate-Q80
marking, and seventeen saturated sections are stored in
[`../artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json).
The current deterministic output SHA256 is
`76c54483c93c7090def42a8dad256838eb9510cd8479d07c5e3123eefa5cfe66`.

This replaces the historical degree-11511 Q80 transport as the primary
equation route.  That transport remains provenance for the alternate frame,
but its size is a property of a poor marked representative rather than an
intrinsic fibration distance.

## The two-dimensional linear system

Write the norm-twelve section on the published short Weierstrass model as

```text
P_w = (Nx/h^2, Ny/h^3),       deg(h)=4.
```

The degree-two-neighbour construction of Brandhorst--Elkies specializes to
sections with coefficient bounds

```text
deg(a) <= 7,     deg(b) <= 1,
```

and regularity condition

```text
a*Nx - b*Ny == 0 mod h^2.
```

The eight coefficient conditions have rank eight in the ten unknown
coefficients.  Their exact kernel therefore has dimension two, as required by
`h^0(X,O(D))=2`.  If `(a0,b0),(a1,b1)` are the stored kernel rows, the new base
coordinate is the ratio `u=L1/L0`, where

```text
Li = ai*(x*h^2-Nx) + bi*(y*h^3+Ny).
```

No interpolation or finite-field reconstruction enters this step.

## Quartic and Jacobian

On a member of the pencil, the residual chord through `-P_w` has slope

```text
m = (a1-u*a0)/((u*b0-b1)*h).
```

For the two residual intersections, the square of their `x`-difference is

```text
r = m^4 - 6*x(P_w)*m^2 - 8*y(P_w)*m
    - 3*x(P_w)^2 - 4*A_old.
```

Exact factor removal gives

```text
r = q(t,u) * square(t,u)^2,       deg_t(q)=4,
```

so the genus-one model is `W^2=q(t,u)`.  The pole of `m` gives the rational
point cut out by the shared old zero; its two coordinates are recorded and
their equation is checked exactly.

For

```text
q = a*t^4+b*t^3+c*t^2+d*t+e,
I = 12*a*e-3*b*d+c^2,
J = 72*a*c*e+9*b*c*d-27*a*d^2-27*b^2*e-2*c^3,
```

the raw short Jacobian is `Y^2=X^3-27*I*X-27*J`.  Its sole denominator is
removed by the stored degree-two gauge.  The resulting polynomial equation

```text
Y^2 = X^3 + A_new(u)*X + B_new(u)
```

has degrees `(8,12)`.  Its discriminant has degree 24, is irreducible and
squarefree over `QQ`, and is coprime to `A_new`.  The fibre at infinity is
smooth.  Hence the complete geometric fibre configuration is

```text
24 I1.
```

The artifact stores every coefficient of `q`, `A_new`, `B_new`, and the
discriminant, together with the pointed quartic-to-Weierstrass normalization.

## Frame and sections

Splitting off

```text
<D,D+O_old> = U
```

has determinant `-1`.  Its positive rank-17 complement has determinant 948
and no norm-two vectors.  Exact PARI `qfisom` supplies an integral isometry to
the alternate Q80 frame and rejects an isometry to published R17.

Sixteen old R17 sections meet `D` once and transport by a Möbius base change
to rational points on the new equation.  The degree-one old sections generate
an index-two sublattice.  The rational old bisection `orbit-0adf9` also meets
`D` once and supplies the missing glue class.  The seventeen exported points
satisfy the new equation exactly; their coordinate matrix in the compiled
frame has determinant `-1`, and their height Gram is rootless of determinant
948.  They are therefore a saturated Mordell--Weil basis, not merely seventeen
independent points.

<!-- status-consumer: EC-K3-R17-NORM12-11952-CONTROL-J-PREIMAGES 1ef38474a0d7f629 -->

## Published controls are not rational fibres

The four published-R17 rank-25--28 controls cannot be reused as literal
rational-fibre controls in the alternate chart.  For each published parameter
`t0`, exact cross-multiplication of the two `j`-maps gives

```text
A_alt(u)^3 B_old(t0)^2 - A_old(t0)^3 B_alt(u)^2 = 0.
```

Each resulting polynomial has full projective degree 24, has no linear factor
over `QQ`, and has no root at infinity.  Thus none of the four curves occurs,
even up to geometric isomorphism or quadratic twist, at a rational alternate
parameter.  This is an exact accessibility obstruction, not evidence against
high-rank rational fibres in the alternate family.  Alternate-native controls
must therefore be found before quotient minima and low-degree visibility can
be compared fairly.

The exact records are stored in
[`../artifacts/generated-results/elkies-k3-r17-norm12-11952-control-j-preimages-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-control-j-preimages-v1.json).

```bash
sage -python \
  elkies-k3/scripts/certify_r17_norm12_11952_control_j_preimages.sage

sage -python \
  elkies-k3/scripts/certify_r17_norm12_11952_control_j_preimages.sage --check
```

## Replay

```bash
sage -python \
  elkies-k3/scripts/compile_r17_norm12_orbit11952_qq.sage

sage -python \
  elkies-k3/scripts/compile_r17_norm12_orbit11952_qq.sage --check
```

The replay uses exact rational arithmetic throughout.  It checks the source
published equation and sections, the Riemann--Roch kernel, quartic radical
identity, pointed Jacobian map, discriminant, primitive `U` splitting,
rootlessness, both integral-isometry decisions, all seventeen new section
equations, and saturation.  The result is an equation-level construction over
`QQ`; it does not classify the corresponding fibration modulo surface
automorphisms (`J1`) or prove a specialization rank above 17.
