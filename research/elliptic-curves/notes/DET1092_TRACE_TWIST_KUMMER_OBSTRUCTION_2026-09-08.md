# Quadratic trace and twist: an explicit 2-Kummer transfer obstruction

## Result and scope

**New algebraic deduction in this application, independently verified.**
Suppose a quadratic multisection supplies `Q` with inherited rational trace
`Z=Q+sigma(Q)`. Its anti-trace `R=Q-sigma(Q)` is a rational point on a
quadratic twist of the elliptic fibre. Under the natural identification of
the two elliptic 2-torsion modules, its Kummer class is exactly inherited:

\[
\boxed{\delta_{E^{(d)}}(R)=\delta_E(Z).}
\]

The explicit square identity below proves this on the stated affine patch
over every field of characteristic different from two. It does not require
a class group, irreducibility of the cubic, a Selmer dimension, or enumeration
of multisections. Thus the obstruction applies to the RR quadratic points
with inherited trace, not merely to one particular conic.

**Verified application.** The unchanged generic orbit8044 conic gives this
inherited class on302 and all eight frozen controls. Every one of these
classes is already represented by a generic rational point on the original
elliptic fibre. Consequently it lies in its Selmer group with zero Sha
image on all nine cases. Global descent survival of this transferred class
cannot distinguish the302 unlock.

**Important boundary.** This is an obstruction to this particular
2-Kummer class transfer, not to an independent rational point, other
classes or higher descent. An independent point can have inherited
2-Kummer image: if `Q` is rational and independent, `2Q-Z` is independent
but has class `delta(Z)`. The earlier conic seed progression remains valid.
No full-rank upper bound, nonzero Sha class, or absence of control seeds is
claimed. Literature novelty of the general identity was not investigated.

## Exact equations and the square witness

**New deduction, verified as polynomial identities.** Write

\[
E:Y^2=f(X)=X^3+aX^2+bX+c_0,\qquad s^2=d\ne0,
\]

and let

\[
Q=(u+vs,w+zs),\qquad \sigma(Q)=(u-vs,w-zs),\quad v\ne0.
\]

Their common chord is `L(X)=mX+n`, where

\[
m=z/v,\qquad n=w-mu.
\]

Put `g(X)=(X-u)^2-v^2d`. The inherited trace is

\[
Z=(c,e),\qquad c=m^2-a-2u,\qquad e=-L(c),
\]

and the first chord identity is

\[
f(X)-L(X)^2=(X-c)g(X).
\]

The chord from `Q` to `-sigma(Q)` has slope `lambda/s`, where
`lambda=w/v`. Its third intersection gives

\[
R=(r,s\eta),\qquad
r=\frac{\lambda^2}{d}-a-2u,\qquad
\eta=\frac{\lambda(u-r)}d-z.
\]

Thus `d*eta^2=f(r)`. Define

\[
N(X)=\lambda(X-u)+dz.
\]

The second chord identity and its exact polynomial consequence are

\[
d f(X)-N(X)^2=d(X-r)g(X),
\]
\[
\boxed{d(r-X)L(X)^2-(c-X)N(X)^2=d(r-c)f(X).}
\]

In the cubic etale algebra `A=K[theta]/f(theta)`, it follows that

\[
\boxed{d(r-\theta)=(c-\theta)
       \left(\frac{N(\theta)}{L(\theta)}\right)^2.}
\]

This is an explicit square witness, not a comparison of norm-square
conditions or numerical scores.

**Patch.** The displayed formulas require `d*v!=0`, smooth `f`, and units
`c-theta`, `L(theta)` and `N(theta)`. The specialization replay checks all
unit norms and every rational-function denominator exactly. It exports
the finite affine exclusion polynomials. Omitted parameters require
another chart; no undefined value is treated as a class or a seed.

## Why this is the correct twisted Kummer class

**Established literature, applied to the verified equations.** For a monic
cubic elliptic model the Kummer representative is `X-theta` in the
norm-square cubic algebra; see
[Cremona--Fisher--O'Neil--Simon--Stoll, section3](https://www.mathe2.uni-bayreuth.de/stoll/papers/Explicit-Descent-III.pdf).
The twist written as `d*eta^2=f(r)` has the monic model

\[
E^{(d)}:Y_d^2=X_d^3+daX_d^2+d^2bX_d+d^3c_0.
\]

Its rational point is `(X_d,Y_d)=(dr,d^2 eta)`, and its cubic root is
`theta_d=d theta`. Therefore its class, after identifying the 2-torsion
modules, is **`d(r-theta)`**, not `r-theta`. Its norm is

\[
\operatorname{Norm}(d(r-\theta))=d^4\eta^2.
\]

The factor `d` must not be discarded when the conic is nonsplit. The square
identity proves that the class equals `c-theta`, whose norm is `e^2` and
whose original-curve rational representative is the inherited point `Z`.

**New consequence.** This class lies simultaneously in the rational Kummer
images for `E` and `E^(d)`. In both Selmer groups its Sha image is zero, by
the [Kummer--Selmer--Sha exact sequence](https://www.mathe2.uni-bayreuth.de/stoll/papers/Explicit-Descent-III.pdf).
This supplies an exact globally rational class on the controls as well as
302, not merely a candidate passing a finite set of local tests.

## Generic-only determinant1092 application

**Verified application.** The inputs are only the recovered parent equation,
its17 generic sections, the existing orbit8044 RR conic, and the unchanged
nine-address roster. The trace word in the inherited basis is

```text
(0,-1,0,0,0,0,0,0,0,1,0,0,0,-1,0,0,0).
```

The replay recomputes this sum on the original function-field elliptic
curve and checks that it is exactly the chord trace. No exceptional point,
historical unlock coordinate, V3 artifact or rank label enters construction.

The original RR quadratic/line equations and `s^2=q(t)` give `x=x0+x1*s`,
`y=y0+y1*s`. We retain their exact rational functions and convert to the
monic cubic by `X=4x`, `Y=8y+4a1*x+4a3`. The formulas above then produce
`R` on the `q(t)`-twist and a square witness over the cubic algebra.

All nine unchanged original-parameter addresses pass:

```text
scale-0131232, scale-0257585, scale-0487239, scale-0177036,
scale-0043332, scale-0590501, scale-0290097, scale-0748009,
302 (t=0).
```

At each address, `q(t)` is a rational nonsquare, with exact integer
square-root inequalities retained. The quadratic point is not a rational
point of the original fibre. The rational twist point exists, but its
transferred class is the inherited trace class. Hence neither twist-point
existence nor Selmer/Sha status of this class separates this panel.

The [equations and square witnesses](../../artifacts/generated-results/elliptic-curves/det1092_trace_twist_kummer_v1/construction.json)
and [nine-case certificate](../../artifacts/generated-results/elliptic-curves/det1092_trace_twist_kummer_v1/panel.json)
are independently replayable. All nine relative class outcomes are zero;
their other seed existence remains `UNKNOWN`, not excluded.

## Relation to the actual first-unlock problem

**Verified prior result.** The first302 elliptic seed has nonzero class
modulo the inherited image, certified by the
[finite-place first-seed gate](DET1092_SEED_KUMMER_COVER_2026-09-08.md#equation-only-finite-place-gate).
Its marked RR Jacobian class and that elliptic class remain different
objects; no linear transfer between those groups is asserted here.

**New obstruction.** Retaining the anti-invariant point instead of the
unordered RR pair does not repair the proposed transfer at the level of
2-Kummer classes. The symmetric pair gave the centre class by addition;
the anti-invariant twist point gives that same class by the explicit square
identity. Neither can recover the non-generic first elliptic Kummer class
without more information.

**Unresolved next step.** A prospective construction still needs a rational
split branch with a certified independent image, or a different
non-generic rational elliptic class. The universal RR curve and the
calibrated first-seed covering are already available, but generic-input
selection that reconstructs302 has not been found. This theorem rules out
one class-transfer mechanism uniformly; it does not justify enlarging a
search or claiming the broader goal complete.

## Independent verification and limits

**Verified computation.** The constructor proves the universal cleared
polynomial identity and the generic RR identities before evaluating any
of the nine old addresses. The checker imports no constructor. It verifies
an exact polynomial factorization and reconstructs the two chord maps.
At each rational fibre it uses elementary cubic multiplication tables and
matrix determinants, rather than a number-field square test or the
constructor's polynomial inverse, to verify the square witness and unit
norms. The generic trace word is independently checked as well.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_trace_twist_kummer.sage
```

Constructor and independent replay each finish in under one second with
separate25-second caps. The
[frozen protocol](../../artifacts/generated-results/elliptic-curves/det1092_trace_twist_kummer_v1/protocol.json)
and [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_trace_twist_kummer_v1/replay.json)
retain exact inputs, exclusions, outcomes and hashes. No point search,
new parameter, class group, Selmer-dimension calculation, paid backend,
detached job or pilot modification was performed.
