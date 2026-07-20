# Gradient dynamics at infinity

For

\[
\mathcal L=16\lVert F-(-1/4,0,0)\rVert^2,
\qquad \dot X=-\nabla\mathcal L(X),
\]

the explicit Palais--Smale curve suggests the weighted chart

\[
r=1/x,\qquad Y=xy,\qquad Z=x^2z.
\]

It is cleaner to use the construction variables

\[
U=1+Y,\qquad
\gamma=1-\frac32Y-\frac12Z.
\]

The Palais--Smale curve is exactly `U=gamma=0`, with `r -> 0`.

## First desingularization

Transform the Euclidean negative-gradient vector field to `(r,U,gamma)`.
Its components have poles at worst `r^-6`. Multiplication by `r^6` clears
the poles and only reparametrizes time in the interior `r != 0`. On the
boundary `r=0`, the resulting polynomial field is exactly

\[
(\dot r,\dot U,\dot\gamma)=(0,0,-32\gamma).
\]

Thus the first compactification does **not** give an isolated hyperbolic
critical point at the end of the Palais--Smale curve. It gives the entire
equilibrium line

\[
r=0,\qquad \gamma=0,
\]

parametrized by `U`. The linearization at every point of this line is

\[
\operatorname{diag}(0,0,-32).
\]

There is one contracting normal direction and two center directions. A
second blow-up or a center-manifold calculation is therefore necessary before
making claims about basin separatrices.

## Formal center dynamics

Writing the desingularized field in the same time variable, its first formal
center-manifold jet is

\[
\gamma=h(r,U)
=3U^2r^4+\frac{49}{8}U^3r^4
 +\frac{99}{8}r^6+\frac{99}{4}Ur^6+O_8(r,U),
\]

where `O_8` denotes terms of total degree at least 8. Substitution gives the
leading reduced field

\[
\dot r=-392r^{11}+O_{12}(r,U),
\qquad
\dot U=-264r^6-528Ur^6+O_8(r,U).
\]

In particular, the tangential drift does not vanish at `U=0` for `r != 0`.
Formally,

\[
\frac{dU}{dr}\sim\frac{33}{49}r^{-5}.
\]

An interior orbit cannot remain in a fixed small neighborhood of a finite
point on the boundary line while `r -> 0`: `U` is driven out of that
neighborhood much faster than `r` reaches zero. Thus the endpoint of the
explicit Palais--Smale curve is not behaving like a conventional saddle at
infinity whose local stable manifold supplies a basin separatrix. This is a
local conclusion; trajectories may still escape in a different blow-up chart
where `U` grows as `r -> 0`.

## The Palais--Smale curve is not an orbit

On `U=gamma=0`, the unrescaled chart field is

\[
\dot r=-392r^5,
\qquad
\dot U=-8(49r^4+33),
\qquad
\dot\gamma=-4(343r^4-99).
\]

In particular, near `r=0`, both `U` and `gamma` immediately leave zero. The
Palais--Smale curve certifies loss of compactness but is not itself a gradient
trajectory or a basin boundary.

These identities are checked exactly by
`scripts/analyze_gradient_infinity.py`.

## Next calculation

The reduced equation suggests the scaling `U` of order `r^-4` along a possible
escape orbit. That leaves the present finite-`U` chart. The next step is to
introduce the corresponding secondary chart (or return to affine variables
with that scaling), then test its equilibria and stable sets. Only after that
step can one determine whether stable sets at infinity form the boundaries of
the three finite attraction basins.
