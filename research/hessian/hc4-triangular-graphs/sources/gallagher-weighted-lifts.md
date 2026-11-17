# Weighted lifts from the Jacobian counterexample

Date: 2026-07-20

This note records deductions made while building the visual guide. They are
exactly checked in `verify.py`, but no claim of literature priority is made.

## 1. The announced map as a weighted lift

Put

```text
v = xy,  t = x^2 z.
```

The announced map has the form

```text
F = (alpha(v,t)/x^2, beta(v,t)/x, x gamma(v,t)),
```

where

```text
alpha = (1+v)^3 t + v^2(1+v)(4+3v),
beta  = v + 3(1+v)^2 t + 3v^2(4+3v),
gamma = 2-3v-t.
```

The coordinate change `(x,y,z) -> (x,v,t)` has Jacobian `x^3`.
Consequently the three-dimensional Jacobian is the two-variable determinant

```text
det [ -2 alpha   alpha_v   alpha_t ]
    [ -beta      beta_v    beta_t  ]
    [  gamma     gamma_v   gamma_t ] = -2.
```

If

```text
P = beta gamma,  Q = alpha gamma^2,
```

then, for every weighted lift of this form,

```text
Jac(P,Q) = -det(JF) gamma^2.
```

For the announced map, `u=1+v` satisfies

```text
(P^3-P^2-18PQ+27Q^2+16Q)u^3 + (P^2-12Q)u - 4Q = 0.
```

This gives generic fiber degree three. The exact target `(-8/27,0,1)` has
three preimages, listed and checked in `verify.py`.

## 2. Tangent geometry and the branch at infinity

After harmless constant rescalings of the output invariants, the announced
map is encoded by

```text
p(w) = 2w - 3w^2,
q(w) = w^2 - 2w^3.
```

These satisfy

```text
q'(w) = w p'(w).
```

Consequently the plane curve

```text
Gamma(w) = (p(w),q(w))
```

has tangent slope `dq/dp=w`. Its tangent at parameter `w` passes through an
output point `(P,Q)` precisely when

```text
Q - q(w) = w(P - p(w)).
```

Since

```text
R(w) = w p(w) - q(w) = integral_0^w p(s) ds,
```

this is exactly the hidden inverse equation

```text
R(w) = wP - Q.
```

Thus inverse images correspond to tangents from `(P,Q)` to `Gamma`.

A repeated inverse root must satisfy both

```text
R(w)-wP+Q = 0,
p(w)-P = 0.
```

It follows that `(P,Q)=(p(w),q(w))`. In the weighted lift this is exactly
the locus `gamma=0`, while reconstruction uses

```text
x = C/gamma.
```

For `C != 0`, the would-be collision of inverse roots therefore occurs only
when the corresponding input escapes to infinity. This is the geometric
mechanism by which the three-dimensional map remains etale at every finite
point.

For fixed nonzero `C`, this gives an explicit family of asymptotic values

```text
(A,B,C) = (q(w)/C^2, p(w)/C, C)
```

up to the constant output normalizations used above. Determining the full
closure, including all `C=0` boundary components, is a next problem; equality
with the complete nonproper-value set has not yet been claimed.

## 3. A general one-variable construction

Fix a nonzero constant `c`. Choose a polynomial `p(w)` satisfying

```text
p(0) = 0,
p(1) = -c,
integral_0^1 p(w) dw = 0.
```

Define `q(0)=0` and

```text
q'(w) = w p'(w) / c.
```

Let `kappa=p'(1)/c`, assume `kappa != -2`, and set

```text
a = -(1+kappa)/(2+kappa).
```

For any nonzero `b`, put

```text
v = xy,
t = x^2 z,
u = 1+v,
gamma = 1 + a v + b t,
w = u gamma,
beta = c + p(w)/gamma,
alpha = u + q(w)/gamma^2.
```

Because `p(0)=0`, the expression `p(w)/gamma` is polynomial in `(u,gamma)`.
Because `q` has a double zero at zero, `q(w)/gamma^2` is also polynomial.
The endpoint and integral conditions imply

```text
beta(0,0) = 0,
alpha(0,0) = 0,
partial alpha / partial v (0,0) = 0.
```

It follows that

```text
F_p(x,y,z) = (alpha(v,t)/x^2, beta(v,t)/x, x gamma(v,t))
```

is a polynomial map.

To calculate its Jacobian, use quotient coordinates

```text
P = beta gamma = c gamma + p(w),
Q = alpha gamma^2 = w gamma + q(w).
```

In coordinates `(w,gamma)`,

```text
Jac(P,Q) = -c gamma.
```

The changes `(v,t)->(u,gamma)->(w,gamma)` have determinants `b` and
`gamma`. Hence

```text
Jac(P,Q) with respect to (v,t) = -bc gamma^2,
det(J F_p) = bc.
```

Finally, eliminating `gamma` gives

```text
integral_0^w p(s) ds = wP - cQ.
```

If `p` has degree `d`, the generic fiber degree is `d+1`.

## 4. An explicit generic-degree-four map

For a particularly small collision, take the raw one-variable data

```text
c = 1,
p(w) = w - 2w^3,
q(w) = (1/2)w^2 - (3/2)w^4,
a = -4/3,
b = -1.
```

The raw lift has `u=1+xy`, `gamma=1-(4/3)xy-x^2z`, determinant `-1`,
and rational coefficients. Precompose by the linear change `y -> 3y` and
scale the first output by `2`. This gives the integer-coefficient map

```text
u = 1+3xy,
gamma = 1-4xy-x^2z,
A = (2u+u^2-3u^4 gamma^2)/x^2,
B = (1+u-2u^3 gamma^2)/x,
C = x gamma.
```

The apparent quotients cancel. The resulting polynomial map has:

```text
component degrees: (12,11,4),
Jacobian determinant: -6,
generic fiber degree: 4.
```

It sends both

```text
(1,0,0), (-1,0,2)
```

to `(0,0,1)`. Generic fiber degree distinguishes this map from the announced
degree-three example under left and right polynomial coordinate changes.

For a complete generic-degree certificate, put

```text
P = B C,
Q = A C^2,
w = u gamma.
```

Then

```text
w^4 - w^2 + 2Pw - Q = 0.
```

Conversely, on the generic open set a root reconstructs the unique input by

```text
gamma = P-w+2w^3,
u = w/gamma,
x = C/gamma,
y = (u-1)/(3x),
z = (1-(4/3)(u-1)-gamma)/x^2.
```

At the target `(A,B,C)=(1,0,1)`, the quartic is

```text
w^4-w^2-1 = 0.
```

Its discriminant is `-400`; its four distinct roots reconstruct four distinct
inputs, all checked exactly in `verify.py`.

## 5. Every generic fiber degree at least three

For every integer `n >= 3`, put `d=n-1` and take the single explicit seed

```text
p_d(w) = 2w - 3w^2 + w(1-w)(w^(d-2) - 6/(d(d+1))).
```

For `d=2`, the last term vanishes identically and this is the quadratic seed
of section 2.  For `d>=3`, its leading coefficient is `-1`, so it has degree
`d`.  In every case,

```text
p_d(0)=0,
p_d(1)=-1,
integral_0^1 p_d(w) dw=0,
p_d'(1)=-5+6/(d(d+1)) != -2.
```

Apply the construction of section 3 with `b=c=1`.  It gives an explicit
polynomial map

```text
F_n : C^3 -> C^3,
det JF_n = 1.
```

Its inverse equation has degree `d+1=n`:

```text
integral_0^w p_d(s) ds = wP-Q.
```

Away from its discriminant, this equation has `n` distinct roots.  A repeated
root is exactly where `gamma=P-p_d(w)=0`; thus on the further open set `C!=0`,
each simple root has `gamma!=0` and reconstructs one distinct input.  Hence

```text
[C(x,y,z) : C(F_{n,1},F_{n,2},F_{n,3})] = n,
```

or equivalently a Zariski-generic target has exactly `n` complex preimages.
This proves that every generic fiber degree `n>=3` occurs in the fixed
dimension three.  The proof is uniform in `n`; the finite checks in the
scripts are regression tests, not the basis for the all-`n` conclusion.

## 6. Rational certificates for degrees 3 through 100

The all-`n` seed above is concise, but its generic inverse roots need not be
rational.  For a deliberately redundant finite certificate, one can instead
choose the inverse polynomial first and make all of its roots integers.

For `3 <= n <= 100`, define a set of `n-2` roots by

```text
S_n = {2} union {j,1-j : 3 <= j <= (n+1)/2}              (n odd),
S_n = {3,4} union {j,1-j : 5 <= j <= (n+4)/2}            (n even),
```

where an upper limit below the lower limit means that the second set is empty.
Put

```text
G_n(w) = w(w+1) product_(r in S_n) (w-r),
lambda_n = G_n'(0)-G_n'(1),
H_n(w) = G_n(w)/lambda_n.
```

Every paired pair `{j,1-j}` contributes one to the ratio
`S_n(0)/S_n(1)`.  The unpaired root `2`, or the pair `3,4`, contributes two.
Thus `S_n(0)=2S_n(1)`, and consequently

```text
G_n(1)=G_n'(0).
```

The exact atlas checker verifies, for each `3 <= n <= 100`, that
`lambda_n != 0` and `H_n''(1) != -2`.  Define

```text
p_n(w) = H_n'(w)-H_n'(0),
q_n(w) = w H_n'(w)-H_n(w),
a_n = -(1+H_n''(1))/(2+H_n''(1)).
```

Then

```text
p_n(0)=0,
p_n(1)=-1,
integral_0^1 p_n(w)dw=0.
```

Applying section 3 with `b=c=1` gives a polynomial map `F_n` with
`det JF_n=1` and generic fiber degree `n`.  These maps have an especially
simple exact fiber.  Set

```text
T_n = (0,-H_n'(0),1).
```

For each of the `n` distinct integer roots `r` of `G_n`, put

```text
gamma_r = -H_n'(r),
X_(n,r) = (
    1/gamma_r,
    r-gamma_r,
    (gamma_r-1-a_n(r/gamma_r-1)) gamma_r^2
).
```

All coordinates are rational.  Since the roots are simple, `gamma_r != 0`.
Exact substitution gives

```text
F_n(X_(n,r)) = T_n.
```

Conversely, over `T_n` the inverse equation is exactly `H_n(w)=0`, and the
reconstruction above is unique for every root.  Therefore the displayed fiber
has exactly `n` points, not merely at least `n`.  Because the inverse equation
has degree `n` and this fiber consists of `n` simple points, the generic fiber
degree is also `n`.

At `n=3`, the roots are `0,-1,2` and the construction gives
`p_3(w)=2w-3w^2`, the minimal seed underlying the announced map.  The atlas
uses the determinant-one normalization `b=1`; Alpöge's displayed map uses a
different nonzero `b` and a linear output scaling.

The standard-library script `verify_counterexamples.py` checks all 98 maps and
all 5,047 rational preimages using exact `Fraction` arithmetic.  This finite
atlas is an audit aid; section 5 remains the proof for every `n>=3`, including
values beyond 100.
