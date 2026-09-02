# Orbit-103 resolved pencil and Weierstrass equation — 2026-09-02

<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT103-EQUATION 8cfa9387612ac443 -->

<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT103-ARITHMETIC-RANK2 2130bc147519ac6b -->

## Result

The sparse quadratic neighbour

```text
D = P0 + P1 + A3_2
```

on the new determinant-36 `E6+A1` K3 has now been compiled through its full
resolved Riemann--Roch pencil, binary quartic, and Jacobian.  Put

```text
D0     = 3*k^2-4,
c      = 2*k/D0,
lambda = -(k^2-4)*D0/4,
H      = 1-t^2.
```

The source K3 is

```text
Y^2 = X^3
    + H^3*(lambda-3H)*X
    + H^4*(c^2*lambda^2+lambda*H-2H^2),
```

with sections

```text
P0 = (-H^2, c*lambda*H^2),

P1 = (4*lambda^2/D0^2 + k^2*lambda*H/D0 - H^2,
      lambda*t*(8*lambda^2/D0^3
                +2*(k^2+2)*lambda*H/D0^2)).
```

After the new base normalization

```text
r = k*(k^2-4)*z/2-k,
```

the child equation is

```text
y^2 = x^3 - 27*(r^2-4)^2*A4(r)*x + 54*(r^2-4)^3*B6(r),
```

where

```text
A4(r) = k^4*r^4
      + 16*k*(k^2+12)*r^3
      + 8*k^2*(k^2+152)*r^2
      + 64*k*(31*k^2-12)*r
      + 16*k^2*(61*k^2-48),

B6(r) = k^6*r^6
      + 24*k^3*(k^2+12)*r^5
      + 12*(k^6+160*k^4+192*k^2+1152)*r^4
      + k*(3072*k^4+35072*k^2+27648)*r^3
      + (1488*k^6+91776*k^4+4608*k^2-55296)*r^2
      + k*(85632*k^4-13824*k^2-110592)*r
      + k^2*(26560*k^4-4608*k^2-55296).
```

The exact replay is
[`scripts/compile_e6a1_rho19_orbit103_rr_weierstrass.sage`](scripts/compile_e6a1_rho19_orbit103_rr_weierstrass.sage),
with generated certificate
[`../artifacts/generated-results/elkies-k3-e6a1-rho19-orbit103-rr-weierstrass-v1.json`](../artifacts/generated-results/elkies-k3-e6a1-rho19-orbit103-rr-weierstrass-v1.json).
The arithmetic descent and secondary-pencil gate replay in
[`scripts/certify_e6a1_rho19_orbit103_arithmetic_and_orbit96_audit.sage`](scripts/certify_e6a1_rho19_orbit103_arithmetic_and_orbit96_audit.sage),
with certificate
[`../artifacts/generated-results/elkies-k3-e6a1-rho19-orbit103-arithmetic-orbit96-audit-v1.json`](../artifacts/generated-results/elkies-k3-e6a1-rho19-orbit103-arithmetic-orbit96-audit-v1.json).

## The resolved Riemann--Roch basis

On the old generic elliptic fibre, let

```text
m = (y1-y0)/(x1-x0)
```

be the secant slope through `P0` and `P1`.  Its apparent coefficient pole is
at `t=2/k`.  It disappears in

```text
L = t-2/k,

L*m = (2/k)*(t^3+(-k^3/4+k/2)*t^2
                  +(k^2/2-3)*t+k/2).
```

The second pencil function is

```text
z = L*(y+y0+m*(x-x0))/((x-x0)*(x-x1)).
```

It has simple horizontal poles at `P0` and `P1`.  At either old `IV*` fibre,
the exact local orders are

```text
ord(A,B)       = (3,4),
ord(x0,y0)     = (2,2),
ord(x1,y1)     = (0,0),
ord(L,L*m)     = (0,0).
```

For an exceptional divisorial valuation of

```text
y^2=x^3+h^3*a*x+h^4*b,
```

the tropical equality implies

```text
min(ord(y),2ord(h)) >= min(ord(x),2ord(h)).
```

The secant quotient is therefore regular on every `E6` exceptional
component.

At the old `I4` fibre, put `s=1/t` and use minimal coordinates

```text
x_s=s^4*x, y_s=s^6*y.
```

The node is `(-1,0)`.  On either side chart

```text
x_s=-1+s*xi, y_s=s*eta,
```

the leading equation and pencil function are

```text
eta^2+3*xi^2=0,
ord(z)=0.
```

On the middle chart

```text
x_s=-1+s^2*xi, y_s=s^2*eta,
```

the exact calculation gives `ord(z)=-1`.  Thus the only vertical pole is the
middle component `A3_2`, with multiplicity one.  Consequently

```text
1, z in H0(X,O(P0+P1+A3_2)).
```

The lattice certificate already proves that this divisor is primitive, nef,
and isotropic.  Hence `h0=2`, so

```text
H0(X,O(D)) = <1,z>.
```

This is a complete resolved pencil, not merely a generic-fibre chord ansatz.

## Quartic and rational origin

Substitute the pencil relation into the old Weierstrass cubic and remove the
two known pole factors `(x-x0)(x-x1)`.  The result is quadratic in the old
`x`; its discriminant is

```text
(t-2/k)^2 * F4(t;z,k),
```

where `F4` is a squarefree binary quartic.  After replacing `z` by

```text
z = 2*(r+k)/(k*(k^2-4)),
```

the quartic has the two rational points

```text
t= 1, w=(r+2)^2/(k*(k-2)),
t=-1, w=(r-2)^2/(k*(k+2)).
```

Thus the genus-one pencil has a rational origin over `QQ(k,r)` and is
isomorphic to its Jacobian.  Classical binary-quartic invariants followed by
the constant scaling unit `k*(k^2-4)` give the displayed clean equation.

## Fibres and Mordell--Weil lattice

The reduced cubic discriminant factors as

```text
4*A^3+27*B^2
 = 2^11*3^13*k^3*(k^2-4)^2*(k^2+4/3)
   *(r+k)^3*(r-2)^7*(r+2)^7*Q4(r),
```

where `Q4` is the squarefree residual quartic serialized in the certificate.
The exact Kodaira profile is

| place | fibre | root |
|---|---|---|
| `r=2` | `I1*` | `D5` |
| `r=-2` | `I1*` | `D5` |
| `r=-k` | `I3` | `A2` |
| `r=infinity` | `I3` | `A2` |
| four roots of `Q4` | `I1` | none |

Hence

```text
fibres     = 2I1*+2I3+4I1,
root type  = 2D5+2A2,
root data  = (14,92,144),
Euler sum  = 24.
```

Since this is the certified orbit-103 fibration on the same generic
Picard-rank-19 K3, Shioda--Tate gives geometric Mordell--Weil rank three.  The
root saturation certificate gives trivial torsion and free height lattice

```text
[ 1/6   1/12  1/12 ]
[ 1/12 17/12  2/3  ]
[ 1/12  2/3  17/12 ],
```

of regulator `1/4`.  Indeed

```text
144*(1/4)=36=|disc NS|.
```

## Two explicit rational points

The quartic covariants give the following polynomial points on the clean
child equation.  First,

```text
X+ = 3*(r+2)*(
       k^2*r^3+(-2*k^2-16*k+96)*r^2
       +(-20*k^2+320*k+192)*r+232*k^2+192*k),

Y+ = 216*(r+k)*(r+2)^2*(
       (k-2)^2*r^2+(-12*k^2+32*k+48)*r
       +52*k^2+80*k+16).
```

Second,

```text
X- = 3*(r-2)*(
       k^2*r^3+(2*k^2-16*k-96)*r^2
       +(-20*k^2-320*k+192)*r-232*k^2+192*k),

Y- = -216*(r+k)*(r-2)^2*(
       (k+2)^2*r^2+(12*k^2+32*k-48)*r
       +52*k^2-80*k+16).
```

Both identities are checked coefficientwise.  The arithmetic descent is now
also exact.  Each section has `P.O=0`, meets one `I1*` spinor component and
the nonidentity component of both `I3` fibres, and therefore has height

```text
4-5/4-2/3-2/3=17/12.
```

They are neither equal nor negatives of one another.  Since geometric torsion
is trivial, equal positive height would force equality up to sign if they were
dependent.  Thus `Q+` and `Q-` are independent over `QQ(k)(r)`.

## The third geometric direction and arithmetic rank

The leading coefficient of the old-base binary quartic is

```text
-48*(r+k)^2/(k^2*(k^2-4)^2).
```

Consequently its two points at old-base infinity are defined over
`QQ(sqrt(-3))`.  Put `j^2=-3`.  Their covariant image on the clean child is

```text
X_delta = (-6*k^2-12)*r^4 - 96*k*r^3
          +(-48*k^2-672)*r^2 - 1152*k*r
          -480*k^2-192,

Y_delta = j*((-18*k^2-24)*r^6
          +(-36*k^3+144*k)*r^5
          +(648*k^2+3168)*r^4
          +(672*k^3+12672*k)*r^3
          +(14112*k^2+12672)*r^2
          +(4032*k^3+20736*k)*r
          +8064*k^2-1536).
```

Exact substitution proves that `Q_delta=(X_delta,Y_delta)` lies on the child.
Conjugation `j -> -j` fixes its abscissa and negates its ordinate, so

```text
sigma(Q_delta)=-Q_delta.
```

It is a nonzero anti-invariant direction.  Since the geometric MW rank is
three, the invariant subspace has rank at most two; the independent invariant
points `Q+` and `Q-` make its rank exactly two.  Therefore

```text
rank E(QQ(k)(r))                 = 2,
rank E(QQbar(k)(r))              = 3,
third direction field           = QQ(sqrt(-3)),
Galois decomposition over QQ    = 2 invariant + 1 anti-invariant.
```

In particular there is **no** third independent generator over `QQ(k)(r)`.
The former arithmetic-rank-three search target was obstructed by descent, not
by failure to find a sufficiently small formula.

## The orbit-96 `A7+D7` comparison gate

The abstract norm-eight census selects orbit `96` as the secondary target.  It
has root type `A7+D7`, height lattice `diag(3/8,1,3)`, two reducible fibres,
and old-basis complexity `(max_abs,L1)=(2,12)`, versus `(1,3)` for orbit `103`.

However, orbit `96` uses old `E6` simple-root coordinates that have not been
identified with a complete resolved physical `E6` marking on the displayed
equation.  The tempting equation-side tangent trace can be tested exactly.
With

```text
m_tan = ((1-3*k^2/4)/k)*(t-1)*(t+1),

z_96,trace = (t+1)*(y-y0-m_tan*(x-x0))/(x-x0)^2,
```

the opposite tangent orientation gives the same Jacobian.  After exact
elimination and minimalization it is

```text
y^2 = x^3 - 3*(z^2-1)^3*C2(z)*x - 2*(z^2-1)^4*C4(z),

C2(z) = z^2-k^4/4+4*k^2/3-7/3,

C4(z) = z^4+(-3*k^4/8+2*k^2-4)*z^2
              -k^6/8+11*k^4/8-4*k^2+3.
```

Its exact fibre profile is

```text
2IV*+I4+4I1,   roots 2E6+A3,   root data (15,156,36),
```

not `A7+D7`.  This low-coefficient equation is therefore a negative control,
not the secondary competitor.  Comparing it with orbit `103` would be a false
optimality claim.

The valid comparison currently is:

| metric | orbit 103 | abstract orbit 96 |
|---|---:|---:|
| old-basis coefficient maximum | `1` | `2` |
| old-basis `L1` | `3` | `12` |
| reducible fibres | `4` | `2` |
| minimal K3 degree bounds `(A,B)` | `(8,12)` | `(8,12)` |
| exact physical Weierstrass equation | yes | not yet |

Thus orbit `103` is proved source-divisor-optimal between the two declared
targets, while orbit `96` remains fibre-count-optimal.  Weierstrass coefficient
optimality is not yet comparable: it requires first resolving the physical
`E6` marking and then compiling the genuine `A7+D7` pencil.

## Reproduction and boundary

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_e6a1_rho19_orbit103_rr_weierstrass.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_e6a1_rho19_orbit103_rr_weierstrass.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_orbit103_arithmetic_and_orbit96_audit.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_orbit103_arithmetic_and_orbit96_audit.sage --check
```

Pinned SHA-256 values are

```text
equation checker   258da7f13d2c0186bd5b4d32eb2debaaa9c685b136d0350c45dac2e85ca72498
equation artifact  c7f08349b9ab874131285cde40023f40b70e35f0d8f90d735d7af446b7a5419a
descent checker    23988e2a5e608301ebf393fb157c0c03a02ae353760aadf1fbe5ed5a689304c1
descent artifact   1daae6cc371f04b541ac1127748db75924098bc959c30ffbe067bbd3027a41d9
```

The resolved `H0` basis, binary quartic, rational origin, Weierstrass equation,
two rational points, exact arithmetic rank two, anti-invariant third geometric
direction, fibre profile, geometric MW rank, torsion, and height lattice are
exact.  A physically marked `A7+D7` equation and a genuine Weierstrass
coefficient-optimality comparison remain open.
