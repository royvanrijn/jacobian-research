# Retaining the intrinsic free variables at the ambient critical value

15 November 2026. A closed quadratic-ambient embedding retaining the two
intrinsic free variables cannot encode a nonautomorphic plane Keller
pair if its zero-free-variable trace is the ambient critical value.
The other three coordinates may have arbitrary polynomial dependence
on all four source variables. Only their zero-section and first free
jet are used. No fibre-affine hypothesis, graph-degree bound, ambient
action lifting or nonlinear source Hessian invariance is assumed.

For a constant trace different from the ambient critical value, the
only remaining helper quadratic has rank one and the base image is
one noncritical plane. A natural universal polynomial substitution in
that remaining chart has unit Hessian for every plane Keller pair,
but its global immersion obligation already forces the pair's inverse.
Thus it does not supply a closed-embedding escape from the criterion.
The full noncritical rank-one case remains open.

Claims HC4MYIU0/1/2 are recorded in
[MATH_STATUS.json](../../MATH_STATUS.json). Both plane endpoints and
unrestricted HC4/JC2 remain open.

## 1. Constant traces restrict all three helper ranks (HC4MYIU0)

Let `Q` be a nondegenerate quadratic polynomial on `A5`, with unique
critical point `o` and critical value `lambda`. Let
`b=(b1,b2)` be independent affine ambient linear functions with
`b(o)=0`. Suppose a genuine closed embedding is written

```
j(z,w)=(w1,w2,X(z,w),Y(z,w),Z(z,w))                (1)
```

in an ambient affine frame starting with `b`, and suppose
`Q composed with j(z,0)=mu` is constant. Base change of the closed
embedding to `b=0` gives an actual closed embedded `A2_z` in the
helper `A3`. This is stronger than a generic immersion or a parameter
substitution.

Center at `o`, subtract `lambda` from `Q`, and put `m=mu-lambda`.
The restriction to the helper three-space `ker b` has rank `r>=1`:
its radical is contained in the two-dimensional orthogonal complement
of that three-space. The image of the zero-section is an irreducible
closed surface in the ternary quadric level `q_r=m`, hence equals
one irreducible component. The possibilities are exact:

- If `r=3,m=0`, that surface is an irreducible singular quadric cone,
  so it cannot be the smooth affine plane.
- If `r=3,m!=0`, write its equation as `XY+Z^2=m`. Localizing
  `X` gives the UFD `C[X,1/X,Z]`, with units `C* X^Z`.
  The removed set has two reduced prime lines `X=0,Z=+/-sqrt(m)`.
  Their only relation is `div X=D_++D_-`, so the class group is
  `Z^2/Z(1,1)=Z`, unlike the affine plane. This elementary divisor
  calculation uses normality of the smooth quadric and the displayed
  localization, not a classification of affine surfaces.
- If `r=2,m!=0`, the level `XY=m` has a nonconstant unit `X`,
  so again is not `A2`.
- If `r=2,m=0`, the zero-section image is one of the planes
  `X=0` or `Y=0`. The other two helper coordinates form an actual
  polynomial coordinate pair of the source `z` plane.
- If `r=1`, write `q_r=X^2`. The image is one plane `X=c`,
  `c^2=m`, and its other two helper coordinates again form a
  source coordinate pair. This includes the reduced image plane
  when the quadric equation is the double equation `X^2=0`.

Equality with a component follows from the prime height-one image
ideal in the helper polynomial ring and the dimension two, not from
surjectivity of an arbitrary dominant map onto a surface. Every
coordinate-pair assertion here follows from the actual closed-image
ring equality with that plane.

## 2. The critical-value trace forces the inverse (HC4MYIU1)

Assume in addition that the actual source potential is

```
Q composed with j=lambda+w1 K1(z)+w2 K2(z),
J(K1,K2) in C*.                                    (2)
```

Then `K` is a polynomial automorphism. In particular the actual
four-variable Hessian is `J(K)^2`, and its gradient is invertible:
recover `z` from its two free derivatives `K`, and then recover `w`
from the remaining derivatives using `DK^(-T)`.

The rank-three case was excluded above. In rank two, constant ambient
changes preserving `b=0`, and an invertible constant row change of
`w`, put the whole nondegenerate quadratic into

```
Q-lambda=XY+w1 Z+gamma w2^2, gamma!=0.              (3)
```

Indeed remove cross terms with the nondegenerate `X,Y` block by
linear helper shifts in `w`. The one remaining helper radical must
pair with one base direction; after choosing it as `w1`, the other
base direction must have a nonzero square. Remaining base quadratic
terms containing `w1` are absorbed by a linear shift of `Z`.

Interchange `X,Y` if needed. On `w=0`, write
`X0=0,Y0=a(z),Z0=d(z)`, where `(a,d)` is a source coordinate
pair by Section 1. Differentiate (2)--(3) at the zero-section:

```
K1=d+a X_(w1)|0,
K2=a X_(w2)|0.                                    (4)
```

In actual polynomial coordinates `(a,d)`, the Keller pair therefore
restricts on the entire straight line `a=0` to `(d,0)`. The retained
Gwozdziewicz one-line theorem HC4MYGL0 gives its polynomial inverse.
Only a two-variable source coordinate change is used for that theorem;
the four-variable Hessian is evaluated from the literal form (2).

In rank one, the analogous constant normalization is

```
Q-lambda=X^2+w1 Y+w2 Z.                            (5)
```

Cross terms with `X` are removed by a helper shift, the two helper
radicals pair nondegenerately with the base two-space, and base
quadratics are absorbed by helper shifts. Section 1 gives
`X0=0,Y0=a(z),Z0=d(z)` with `(a,d)` a coordinate pair.
The first free derivatives of (5) are simply

```
K1=a, K2=d.                                        (6)
```

Thus `K` is already that coordinate pair. Higher free derivatives
of all three embedding helpers are unrestricted and do not enter
either rank argument.

### Intrinsic recognition and polynomial base offsets

The criterion applies to the intrinsic cotangent form reconstructed
by IT2. It suffices that two ambient affine linear functions centered
at `o` pull back to

```
b composed with j=A w+v(z), A in GL2(C),
D(z)-v(z)^T A^(-T)K(z)=lambda.                     (7)
```

Use the actual source coordinates `(z,b)`; the potential is literally
`lambda+b^T A^(-T)K(z)`. Its coefficient pair is Keller, so (2)
applies. This verifies the Hessian in the new coordinates from the
cotangent identity, rather than assuming invariance under the nonlinear
offset. All degrees of `v` are allowed. A base-dependent matrix
`A(z)` cannot be substituted in (7) without separately verifying that
`A(z)^(-T)K(z)` is Keller.

For a constant trace `mu!=lambda`, Section 1 leaves only rank one.
In the normalization (5), its exact zero-section and first jet are

```
X0=c!=0, Y0=a(z), Z0=d(z), c^2=mu-lambda,
K1=a+2c X_(w1)|0, K2=d+2c X_(w2)|0.                (8)
```

Thus the critical-value hypothesis has a precise geometric role:
it removes those normal first-jet corrections. The full noncritical
rank-one chart is not asserted solved by this argument.

## 3. A universal compensation does not evade the immersion gate (HC4MYIU2)

Let `K=(K1,K2)` be any plane Keller pair normalized to `J(K)=1`.
Fix `c!=0`, set `z=(s,t)`, and define

```
A(z)=(K(z)-z)/(2c), H=w^T A(z),
J_(c,K)(z,w)=(w1,w2,c+H,s-H A1,t-H A2).             (9)
```

This polynomial substitution satisfies the exact function identity

```
(X^2+w1Y+w2Z) composed with J_(c,K)=c^2+w^T K(z).    (10)
```

It has actual unit Hessian and its zero-section is a closed affine
plane at the noncritical value `c^2`. Nevertheless, **if (9) is
an immersion at every source point, then K is a polynomial
automorphism**. Closed embedding is stronger, so (9) cannot repair
a hypothetical nonautomorphic pair by a genuine embedding.

Here is the all-degree immersion obstruction. At a fixed `z`, write
`B=DA(z)`. A source tangent `(v,0)` is killed by the derivative
of (9) if

```
w^T Bv=0, (I-H B)v=0.                              (11)
```

If `Bv=theta v` with `theta!=0` and `A(z)` is not parallel
to `v`, choose `w` with `w^T v=0,w^T A=1/theta`.
This is possible in dimension two and satisfies (11). Thus global
immersion forbids such an eigenvector.

Since `det(I+2c B)=J(K)=1`, a nonzero eigenvalue cannot occur
alongside a zero one. If `B` has two distinct eigenvalues, both
are nonzero, and at any point with `A(z)!=0` one eigenvector
is not parallel to `A(z)`. Immersion would fail. The nonvanishing
locus of a nonconstant nonzero `A` is dense. Hence the discriminant
of the characteristic polynomial of `B` is identically zero.
If `A` is constant, the inverse is affine and needs no argument.

Writing `r=tr B`, the determinant-one identity and zero discriminant
give

```
r+2c det B=0, r^2=4 det B,
r(1+c r/2)=0.                                     (12)
```

The polynomial source ring is a domain, so either `r=0,det B=0`,
or `r=-2/c,det B=1/c^2`. In the first case `DA` is nilpotent
everywhere. In the second, `DK` has trace `-2` and determinant
one, so `D(K+z)=DK+I` is nilpotent everywhere. In both cases
the elementary plane lemma below makes `K` invertible.

### The plane nilpotent-Jacobian lemma used here

If a polynomial vector `P:C2->C2` has nilpotent Jacobian, then
`epsilon z+P(z)` is a polynomial automorphism for every nonzero
constant `epsilon`. To see this, `det DP=0` makes its components
algebraically dependent. The retained polynomial generative-function
theorem used in IL0 supplies `P1=p(h),P2=q(h)` for one polynomial
`h`. Constant-component cases are included. The trace-zero identity is

```
p'(h) h_s+q'(h) h_t=0.                             (13)
```

For every scalar `a` outside the finite common-zero set of `p',q'`,
each irreducible component of `h=a` is tangent to the nonzero constant
vector `(p'(a),q'(a))`. The directional derivative of its irreducible
equation is divisible by that equation; its smaller degree makes
the derivative zero. In coordinates adapted to that vector, the
equation is an irreducible one-variable polynomial, hence a line.
Components from different levels are disjoint, so those lines must
all be parallel. Therefore the vector direction is constant on a
dense set of levels, and a fixed constant directional derivative
annihilates `h`. Thus `h` is a polynomial in one linear form `ell`.
Both components of `P` depend only on `ell`; trace zero then gives
`ell(P)=constant`. The map `epsilon z+P(z)` is triangular in
constant coordinates starting with `ell`, with a polynomial inverse.
If both components are constant the assertion is immediate.

This proof uses the retained polynomial generative theorem, not a
nilpotent-Jacobian classification in dimension three or more. In the
second case of (12) apply it to `P=K+z,epsilon=-1`; in the first
apply it to `P=2c A,epsilon=1`.

The implication in IU2 is not a converse. Even a genuine nonlinear
plane automorphism can make (9) fail immersion. For

```
K=(s+t^2, t+(s+t^2)^2/3), c=1,
z=(0,1), w=(3,-3), v=(1,1),
```

one has `A=(1/2,1/6)`, `DA v=v`, `w^T v=0`, and `H=1`.
The tangent `(v,0)` is killed by (9), although (10) still has
Hessian determinant one. A substitution with unit Hessian and an
embedded zero-section is not necessarily an immersion or a closed
embedding of the full source.

## Sources and remaining scope

The one-line theorem is the retained explicit import
[Gwozdziewicz, Theorem 1.1](https://arxiv.org/abs/alg-geom/9305008),
with the unchanged primary PDF pinned in the
[earlier source record](sources/critical-jet-import.json). The
polynomial generative theorem is the retained IL0 input in
[Arzhantsev--Petravchuk](https://arxiv.org/html/math/0608157v2).
All quadric, unit, class-lattice, first-jet and spectral arguments
needed here are written above. No general surface rectification,
ambient action classification, or inverse assertion for higher-dimensional
nilpotent-Jacobian maps is imported.

The [source record](sources/free-projection-trace-context.json),
[checker](verify_free_projection_trace.py),
[receipt](receipt-free-projection-trace.json) and
[validation](validation-free-projection-trace.json) separate exact controls
from the written all-degree proof. Ambient translations and distance
critical-point counts were also investigated, but no assertion that
all ambient translations preserve unit Hessian or realize arbitrary
intrinsic target shifts is used. Source parameter recovery remains
part of every embedding hypothesis. The noncritical rank-one chart
outside (9), arbitrary projection data, actual endpoint representatives,
both plane endpoints and unrestricted HC4/JC2 remain open. No formal
verification, independent replay, external review or novelty priority
is claimed.
