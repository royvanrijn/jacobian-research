# Section-first normal forms for MW1 and MW2 searches — 2026-09-02

## Outcome

The equation search can now start with the arithmetic marking built in.
The reusable implementation is
[`scripts/compile_section_first_normal_forms.sage`](scripts/compile_section_first_normal_forms.sage),
and its exact control artifact is
[`../artifacts/generated-results/elkies-k3-section-first-normal-form-controls-v1.json`](../artifacts/generated-results/elkies-k3-section-first-normal-form-controls-v1.json),
SHA-256 `6edfa5f3487f020a05772bd2b1a6b5d74586c126d852b87418285ef46f843c34`.

This is a normal-form and regression result.  It does not prove that a
section-first coefficient space contains a desired characteristic-zero
foundry surface, that its displayed sections generate the full Mordell--Weil
group, or that a tuned fibre model has the intended saturated Neron--Severi
lattice.

## MW1: one marked section is part of the equation

Over a function field of characteristic different from two, translate a
non-2-torsion marked point to the origin and shear by its tangent.  The result
has the Tate-style form

```text
y^2 + a1*x*y + a3*y = x^3 + a2*x^2,
P = (0,0).
```

Thus a one-section search has no separate section equation.  For an elliptic
K3 over `P1`, start with degree bounds

```text
deg(a1)<=2,  deg(a2)<=4,  deg(a3)<=6
```

and impose the required fibre jets on the discriminant.  Exact fibre order,
the `c4` unit gate, minimality, splitness, and resolved component incidence
remain explicit checks; discriminant divisibility alone is not a Kodaira or
marking certificate.

For a known short model

```text
y^2=x^3+A*x+B,   P=(xp,yp),   yp != 0,
```

the exact change of variables is

```text
x = X+xp,
y = Y+m*X+yp,
m  = (3*xp^2+A)/(2*yp),

a1 = 2*m,
a2 = 3*xp-m^2,
a3 = 2*yp.
```

The new marked point is `(X,Y)=(0,0)`, its tangent is `Y=0`, and `c4`,
`c6`, and `Delta` are unchanged because this is a unit Weierstrass
transformation.  Translation can introduce rational-function coefficients;
that is why a new search should be parameterized in the Tate chart from the
outset rather than translated only after solving a short-Weierstrass system.

## MW2: compile both points and their intersection

Let `h` be the desired affine intersection divisor of two marked sections.
Choose polynomials `a1,h,r,s,k` with

```text
gcd(r,h)=gcd(r,s)=1
```

and compute a polynomial Bezout relation

```text
alpha*s + beta*r^2 = 1.
```

Put

```text
R  = h*(r^3-h*s^2-a1*r*s),
a3 = alpha*R+k*r^2,
a2 = -beta*R+k*s.
```

Then the Tate-style equation above contains the two points

```text
P=(0,0),
Q=(h*r,h^2*s)
```

identically.  Indeed, after substituting `Q` and dividing by `h^2`, its
equation is exactly

```text
a3*s-a2*r^2 = R.
```

The Bezout formula makes the right-hand relation an identity.  Moreover

```text
gcd(h*r,h^2*s)=h,
```

so the affine coincidence divisor is already `h`.  If `h` is coprime to the
discriminant, `r` is a unit along `h`, and there are no omitted intersections
at poles or infinity, the smooth section intersection number on this chart is
`deg(h)`.  The extra factor `h` in the second coordinate is forced by the
tangent normalization: two sections meeting at `P` agree to first order with
the tangent line, so the tangent-normal coordinate vanishes to second order.

The search variables are now `a1,h,r,s,k`, subject to degree and coprimality
open conditions.  The coefficients `a2,a3`, both section equations, and the
pair-incidence equation are compiled away before any fibre solve.  The
remaining closed equations are discriminant jets; the geometric open gates
are checked afterward.

## Positive controls

The compiler replays two independent marked models.

### Golay `G720-S0128`

The exact rational `3I6+6I1` specialization translates to the two-section
chart.  Its marked intersection polynomial is

```text
h=t^2-(6/5)*t+1/25,
```

which has degree two and is coprime to the discriminant.  The translated
second point factors exactly as `(h*r,h^2*s)`, the two-section parameterization
recovers all Tate coefficients, and the discriminant has exact `I6` orders at
`0,1,infinity` with `c4` a unit at each support.

This is a positive control for the normal form only.  The same rational model
has hidden rational 3-torsion and a rational half-section, so its saturated NS
determinant is 20 rather than 720; the existing saturation rejection remains
in force.

### NS0031 `S001`, model 157

The exact marked `GF(7)` pair likewise translates to the chart.  Its smooth
intersection polynomial is

```text
h=t^2+3*t+2,
```

and the fibre orders are exactly `I2` at zero and `I8` at one and infinity,
again with the required `c4` unit gates.  Both section identities, the
`(h*r,h^2*s)` factorization, and recovery of the Tate coefficients are exact
over `GF(7)(t)`.

This preserves the existing arithmetic boundary: model 157 is a finite-field
marked point with a finite 7-adic lift, not a rational equation or an
algebraized characteristic-zero NS0031 family.

## Literature placement

The one-point chart is the untorsioned version of the Tate/Kubert idea: put a
marked point at `(0,0)` and use the tangent to remove the linear `x` term.
Kubert's classical two-parameter specialization is designed to impose torsion
relations; here the coefficient functions remain free and fibre conditions
replace torsion equations.  See D. Kubert,
[*Universal bounds on the torsion of elliptic curves*](https://doi.org/10.1112/plms/s3-33.2.193),
*Proc. London Math. Soc.* **33** (1976), 193--237.

For two independent points, the global literature naturally uses a `dP2`
cubic and then maps it birationally to Tate/Weierstrass form; see M. Cvetic,
D. Klevers, and H. Piragua,
[*F-Theory Compactifications with Multiple U(1)-Factors: Constructing Elliptic Fibrations with Rational Sections*](https://arxiv.org/abs/1303.6970),
especially the two-point model and its Weierstrass map.  The present
`Q=(h*r,h^2*s)` chart is a narrower affine compiler adapted to polynomial K3
fibre searches and to a prescribed section-intersection divisor.  It is not a
claim that every global two-section fibration admits this one polynomial chart
without changing base charts or line-bundle trivializations.

## Reproduction

Generate and then byte-check the control artifact with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_section_first_normal_forms.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/compile_section_first_normal_forms.sage --check
```

The next application should instantiate degree-bounded MW1 and MW2 coefficient
charts for a chosen source profile, eliminate `a2,a3` by the formulas above,
and solve the discriminant jets.  Every surviving rational reconstruction
must still pass torsion, divisibility, component, Picard, and NS-saturation
gates before it is identified with a foundry lattice.
