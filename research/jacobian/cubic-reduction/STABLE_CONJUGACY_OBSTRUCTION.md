# Twelve active coordinates cannot be removed by stable polynomial conjugacy

KCUB25, 15 November 2026. Fix the integral twelve-variable cubic K of
[KCUB3](PROOF12.md). For every m>=0, if K times identity_m is polynomially
conjugate to G times identity_s, then s<=m and dim G>=12. The same holds
for invertible formal conjugacies at the corresponding fixed points.
No degree bound is placed on the conjugacy or G.

Thus12 is the minimum number of active coordinates in this stable
conjugacy class. This is not a global lower bound for cubic Keller maps.
Independent source and target changes, graph restrictions and general
stable left-right equivalence are outside the theorem. In particular it
does not obstruct the asymmetric stable reductions already used to relate
the known maps. It also does not prevent a conjugate with fewer terms
while retaining twelve active coordinates. Global ordinary-cubic dimension
remains[4,12], and all positive simplicity bounds are unchanged.

## A finite obstruction to any invariant coordinate

Write K(X)=X+H2(X)+H3(X), with homogeneous degrees2,3. If a formal series
h=h0+ell+h2+h3+... satisfies h(K(X))=h(X), its degree2 and3 equations are

```
ell(H2)=0,
ell(H3)+D2(h2)=0,       D2=sum_i H2_i*d/dx_i.
```

The unknown cubic part h3 cancels because the linear part of K is I;
higher parts cannot affect these degrees. Consequently this argument
applies to every formal series, not just bounded-degree polynomials.

The twelve quadratic output parts have rank11. Row5 is zero, and the
11-by11 coefficient minor on rows

```
1,2,3,4,6,7,8,9,10,11,12
```

and columns

```
x7*x12, x7*x11, x6*x7, x2^2, x1*x10, x1*x9,
x1*x8, x1*x5, x1*x3, x1*x2, x1^2
```

has determinant -288. Thus ell must be lambda*x5.

For a homogeneous cubic p define

```
Lambda(p)=[x2^3]p+(2/3)*[x3*x7*x12]p
          +(1/12)*[x7*x8*x12]p-(1/4)*[x7*x9*x12]p.
```

This annihilates D2(q) for every homogeneous quadratic q. One direct
coefficient explanation is to let c3,c8,c9 be the coefficients of
x2*x3,x2*x8,x2*x9 in q. The four sampled coefficients of D2(q) are

```
8*c3+c8-3*c9,   -12*c3,   -12*c8,   -12*c9,
```

respectively. The displayed linear combination is identically zero.
The independent checker also verifies this on all78 quadratic monomials.
But H3_5=x2^3, so Lambda(H3_5)=1. Applying Lambda to the degree3 equation
forces lambda=0. Every invariant formal function therefore has dh(0)=0.
This does not assert that all invariant formal functions are constant;
functions with zero linear jet remain outside this conclusion.

A polynomial coordinate has nonzero gradient at every point, since it is
a row of an invertible Jacobian. Hence no invariant polynomial coordinate
exists, and no polynomial conjugacy can make K fix a coordinate, even
if intermediate degrees are arbitrarily large or the conjugacy does not
fix the original origin. Pullbacks of target coordinate functions would
be invariant coordinates at the original fixed point0.

## Stabilizations do not evade the obstruction

For K(X) times identity_m(Z), write a formal invariant h(X,Z). Its degree2
equation still says that its X-linear part is lambda*x5; its Z-linear part
is arbitrary. Set Z=0 in the degree3 equation. The same Lambda argument
forces lambda=0, because cross terms of h2 containing Z vanish on Z=0.
Thus the entire space of possible first derivatives of formal invariants
is exactly the m-dimensional span of dZ1,...,dZm. These derivatives are
attained by the obvious invariant functions Zj.

If a conjugacy split s identity coordinates, their pullbacks would give
s invariant functions with linearly independent differentials at(0,0).
Therefore s<=m. Ambient dimensions give dimG=12+m-s>=12. The identity
conjugacy attains12, so this is a matching bound in this precise class.
The argument is purely characteristic-zero formal algebra over the base
field; no analytic convergence or embedding in C is needed.

## Independent finite certificate and scope

Run

```
.venv/bin/python research/jacobian/cubic-reduction/verify_stable_conjugacy12.py
```

The [certificate](stable-conjugacy12/certificate.json) specifies the minor
and four functional weights. The checker reconstructs H2,H3 from the pinned
K coefficient list, computes the rational determinant, verifies the kernel
claim and all78 Lie derivatives, and checks the degree3 pullback formula
on all degree<=3 monomials. Identity stabilization controls check the
projection mechanism. These arithmetic statements are independently
replayed; the formal-series and arbitrary-m stabilization arguments above
are written, not formalized in a proof assistant.

- [Validation receipt](stable-conjugacy12/validation.json)
- [Independent checker](verify_stable_conjugacy12.py)
- [Prior input and full determinant/collision verification](integer-validation12.json)

The result adds no smaller witness and no global cubic dimension lower
bound. A dimension improvement from this construction must use unequal
left/right transformations or another mechanism; arbitrary conjugacies
alone cannot create a removable identity factor.

## Coupled affine translation experiment retained separately

[joint_translation12.py](joint_translation12.py) tested97 distinct coupled
carrier-zero graph centres, excluding the old coordinate-axis centres,
with exact inverse-derivative target normalization. The smallest displayed
count was145, so none improved36. These are finite rational samples, not
a theorem about arbitrary affine centres or subsequent linear frames.
The four best stored maps all have count145. The new checker reconstructs
the selected maps and the complete finite count table separately from
the universal conjugacy proof.

These tangent-normalized translated maps use unequal source and target
affine transformations. They are outside the conjugacy obstruction above;
its proof must not be used to exclude further coupled affine searches.
The unchanged97-centre grid need not be rerun. The full12-variable affine
orbit and coordinated nonlinear left-right transformations remain open.
