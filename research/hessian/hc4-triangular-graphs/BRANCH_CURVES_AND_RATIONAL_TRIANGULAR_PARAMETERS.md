# Intermediate branch curves and rational triangular parameters in the plane

This note proves a plane-only geometric strengthening of the retained
[pencil tests](REDUCED_PENCILS_AND_RATIONAL_MONODROMY.md). Every rational
parameter affine in the second target coordinate remains relatively
algebraically closed after plane Keller pullback, with no exceptional
fibre or multiplicity restriction. It closes rational affine fixed-polynomial
outputs with rational output translation, gives an all-factor inertia
criterion, and completes the cubic `(T,xy-P(T,z))` output class.

Claims `HC4MYII0/1/2/3/4` are recorded in
[MATH_STATUS.json](../../MATH_STATUS.json). This is not a solution of
either plane endpoint or unrestricted HC4/JC2. The source is the
original affine plane throughout. The higher-dimensional multiplicity
criteria remain distinct; the curve argument below is plane-specific.

Let `K=(f,g):A2->A2` be complex Keller, `A=C[T,V]`, `B=C[s,t]`,
`F=C(f,g)` and `L=Frac(B)`.

## 1. Intermediate ramification lies in nonproperness (HC4MYII0)

If `F subset E subset L` is finite, and the finite normalization of the
target in `E` has a ramified divisorial place over an irreducible target
curve `C`, then `C` is a component of the nonproperness set `S_K`.
Consequently `C` admits a nonconstant polynomial map from `A1`.
The intermediate field need not be Galois.

Here is the complete missing-sheet argument. By the retained
[HC4MYOR2](../hc4-relative-row-orbit/PROOF.md#6-properness-and-the-exact-missing-length)
and [HC4MYIB1](NORMAL_QUOTIENT_AND_BOUNDARY.md), the source `A2` is a
dense open subset of the finite normalization `N_L` of the target in
`L`. The normalization `N_E` is finite as well, and `N_L->N_E` is
finite. Choose a prime divisor of `N_L` above a ramified divisor of
`N_E`. Multiplicativity of divisorial ramification indices gives total
index greater than one over `C`. It cannot meet the source open:
every visible source prime has index one by the square-free pullback
proof [HC4MYHW0](ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md). If the
divisor met that open at even one point, openness would make its generic
point visible. Thus the whole divisor is omitted.

Its finite image is `C`. At a general closed point on this omitted
divisor, approach through the dense source open. The images under the
finite normalization map converge to the corresponding point of `C`.
The source coordinates cannot remain bounded: a bounded subsequence
would converge inside `A2`, and continuity in `N_L` would identify a
source point with the omitted point. Hence those target points are
nonproper for `K`.

The nonproperness set is closed and has dimension at most one: clear
denominators in integral equations for the source coordinates to make
`K` finite over a dense principal target open. Therefore `C` is a
component. The explicitly retained Jelonek--Lason Theorem 3.2 supplies
nonconstant polynomial curves through its general points. Such a curve
must lie in `C` after avoiding intersections with the other components.
This proves the assertion. Only the qualitative theorem is imported.

Other sheets above `C` may still be visible and unramified. This
argument establishes nonproperness, not omission of every target point.

### A useful forbidden branch curve

If `u(T),v(T)` are coprime and `u` is nonconstant, the curve

```
u(T)V+v(T)=0                                (1)
```

admits no nonconstant polynomial `A1` map. Its coordinate ring is
`C[T,1/u(T)]`: Bezout for `u,v` makes `u` invertible on the curve,
and the equation recovers `V`. Along a polynomial parametrization,
`u(T(x))` would be an invertible polynomial, hence constant. Then `T`
and `V` would both be constant. Thus (1) cannot be a branch curve of
an intermediate field inside `L`.

## 2. Every affine-in-V rational parameter stays regular (HC4MYII1)

Let

```
q(T,V)=(Q(T)V+C(T))/D(T),
Q,D!=0, gcd(Q,C,D)=1.                       (2)
```

There is no square-freeness, exceptional-value-count or multiplicity
restriction. Then `C(q(f,g))` is relatively algebraically closed in `L`.
Every finite Galois cover of this one-variable field retains its full
group after source extension, as in HC4MYIG1.

If `Q` is constant, every finite numerator member is nonconstant with
unit `V` derivative, and HC4MYIG0 already gives the result. Now suppose
`Q` nonconstant and, contrary to closedness, write

```
q(f,g)=Phi(h), h in L, deg Phi=d>1.          (3)
```

Existence of this rational composition is exactly the explicit retained
Petravchuk--Iena Lemma 2(1). Since `V=(Dq-C)/Q`, we have
`F=C(f,q)`. In particular `f` is transcendental over `C(q)`, hence also
over its finite extension `C(h)`. Thus `E=C(f,h)` is an intermediate
field of degree `d` over `F`, obtained by adjoining the independent
transcendental `f` to the rational curve cover `Phi:P1->P1`.

For finite `lambda`, put

```
G_lambda=gcd(Q,C-lambda D),
L_lambda=(Q V+C-lambda D)/G_lambda.          (4)
```

This is the irreducible primitive horizontal factor, linear in `V`.
The target valuation of `q-lambda` at `L_lambda=0` is one. Its generic
residue field is `C(T)`, since `T` is independent and the `V` coefficient
is nonzero there. A branch point of `Phi` over `lambda` therefore gives
a ramified place of `E/F` over this curve, with the same local index:
the completion is the constant extension by `C(T)` of the rational
curve completion in the uniformizer `q-lambda`. Section 1 applies.

The forbidden-curve observation forces `Q/G_lambda` to be constant at
every finite branch value. Equivalently

```
Q divides C-lambda D.                       (5)
```

At most one finite `lambda` satisfies (5). Two distinct values would
make the nonconstant `Q` divide both `C` and `D`, contradicting
`gcd(Q,C,D)=1`.

Thus `Phi` has at most one finite branch value, with infinity also
allowed. Riemann--Hurwitz excludes a nontrivial cover with at most one
branch value, and with exactly two it forces a single totally ramified
point over each: their contributions must add to `2d-2`. A projective
change on the `h` line therefore gives

```
q-lambda_0=c w^d, c in C*, w in C(h).        (6)
```

But `q(K)-lambda_0` cannot be a `d`th power. Its primitive horizontal
target factor in (4) has multiplicity one and is coprime to the target
denominator and vertical content. HC4MYHW0 gives source prime factors
with valuation exactly one. This contradicts (6). Alternatively
`E/F` would be a nontrivial cyclic Galois intermediate field, forbidden
by HC4MYIB0. This proves closedness.

The finite-cover and filtered class-action consequences are precisely
the HC4MYIG1 field argument, now with unrestricted parameter (2).
Every global filtration and affine height-one fibre hypothesis stays
in force. This statement is restricted to plane Keller maps: a punctured
curve times affine space can contain polynomial curves in the free
directions, so the higher-dimensional multiplicity tests are not retired.

## 3. Rational affine fixed-polynomial outputs (HC4MYII2)

Let `R in C[Y]` be nonconstant, and let

```
A,c in C(T)*, beta,B0 in C(T),
P(T,z)=A(T)R(c(T)z+beta(T))+B0(T) in C[T,z]. (7)
```

Now the output translation `B0` may also be rational. Final
polynomiality is required. For every plane Keller pair,
`P(f,Z)+g` stays irreducible over `L` and preserves its full target
splitting group, in every degree.

Indeed `q=(V+B0(T))/A(T)` has a reduced form (2), after clearing
denominators and their common gcd. Section 2 makes its source field
regular. The source root equation divided by `A(f)` is `R(Y)+q(K)=0`
under the invertible field change `Y=c(f)Z+beta(f)`. The generic
equation `R(Y)+q` is irreducible and its full finite splitting field
retains its group. Since `F=C(T,q)` is purely transcendental over
`C(q)`, this is also its full target group.

Every carried closed stable embedding through `(T,xy-P)`, `r>=1`,
forces the plane Keller inverse by the retained HC4MYHZ2 global
filtered proof. This closes all rational affine-isotrivial polynomial
presentations (7), with arbitrary source, coefficient and embedding
degrees. No source Hessian invariance is inferred from the field change.

## 4. Nonparametrizable branch inertia controls all factors (HC4MYII3)

Let `M/F` be finite Galois, `G=Gal(M/F)`. Let `N_bad` be the normal
subgroup generated by inertia groups over target curves that admit no
nonconstant polynomial `A1` parametrization. Let `H` be the image of
`Gal(ML/L)` in `G`, so `M^H=M intersect L`.

Then

```
N_bad is contained in H.                    (8)
```

The intermediate field `M^H` lies in `L`. By Section 1 it is unramified
over every such target curve, at every place. Inertia must therefore
act trivially on every coset of its subgroup `H`; equivalently each
inertia group is contained in the normal core of `H`. This standard
inertia description is one of the retained Stacks `0BSD` inputs.
Taking normal closures gives (8).

Consequences, with supplied finite-cover or class-action certificates:

- If an irreducible separable target polynomial has splitting field
  `M`, and `N_bad` acts transitively on its roots, it stays irreducible
  over `L`. A source factor would select an `H`-invariant, hence
  `N_bad`-invariant, proper subset.
- In particular a primitive faithful root group and one genuine
  ramified curve of this kind suffice: a nontrivial normal subgroup
  of a primitive permutation group is transitive, since its orbits
  form a block system.
- If `G/N_bad` is nilpotent, the retained HC4MYIA0 preserves that
  quotient cover's full group. Together with (8) this gives `H=G`,
  so the entire splitting group is preserved. Other already certified
  preserved quotients may be used in the same way.
- In the retained filtered output criterion, if the geometric class
  group's `N_bad`-fixed subgroup is zero, the source component class
  is zero. The global filtered proof then forces the carried plane
  inverse. All the original ring and affine fibre hypotheses remain.

Discriminant factors alone are not asserted genuine inertia: use the
normal field cover and actual divisorial ramification. Nor is full
splitting-group preservation inferred merely from transitivity.

### An unbounded nonadjacent-binomial family

For every `r>=1`, `1<=q<m`, and polynomial `d(T)`, the outputs

```
P=(T-c_0)^r z^m+z^q+d(T)                    (9)
```

close for every carried plane embedding, and their full generic root
group is preserved. No primitivity or group classification is needed.
Translate `T` and shear `V` to take `c_0=d=0`. Away from the axes
`T V=0`, every nonzero critical root obeys

```
m T^r z^(m-q)=-q,
V=-(m-q)z^q/m,
T^(rq) V^(m-q)=(-1)^m q^q(m-q)^(m-q)/m^m.   (10)
```

Each irreducible component of the last curve lies in the torus and
admits no polynomial `A1` map: both `T` and `V` would be polynomial
units and hence constant. Every remaining finite branch curve lies
in the coordinate axes; outside the leading coefficient and the root
discriminant the integral distinct roots have trivial inertia.

The quotient by `N_bad` is therefore unramified away from the axes.
Purity makes its finite normalization etale over `(C*)^2`. The retained
finite-etale/topological-cover comparison makes its group a finite
quotient of `pi_1((C*)^2)=Z^2`, hence abelian. The nilpotent-quotient
case above proves full preservation. This includes all `m,q,r`, even
when the generic root group is imprimitive. It is not a claim about
arbitrary coefficient pairs `a(T)z^m+b(T)z^q`.

## 5. Every cubic output closes (HC4MYII4)

Let

```
P=a(T)z^3+b(T)z^2+c(T)z+d(T), a!=0,
Delta=b^2-3ac.
```

For every plane Keller pair, `P(f,Z)+g` is irreducible. Thus every
carried closed stable embedding through `(T,xy-P)`, `r>=1`, forces
the Keller inverse, with all coefficient, leading, source and embedding
degrees unrestricted. Here `Delta` is the critical-point discriminant,
not the root discriminant in `z`.

The retained HC4MYIC2 already excludes rational roots if the nonzero
square-free part of `Delta` has degree at least two. If `Delta=0`,
the polynomial is a translated power over `C(T)` and (7) applies. If
`Delta=H^2!=0`, use

```
Y=(3az+b)/H, A=H^3/(27a^2),
B0=d+(2b^3-9abc)/(27a^2),
P=A(Y^3-3Y)+B0.                             (11)
```

Again Section 3 applies, with no pole restriction on `B0`. The
remaining case is

```
Delta=H(T)^2 ell(T), ell=T-c_0.              (12)
```

Nonzero complex constants may be absorbed in `H`. Suppose a rational
source root exists; a reducible cubic would have one. The same change
as (11) gives

```
P=A(Y^3-3ell Y)+B0.                         (13)
```

### 5.1 The critical-value curve forces integrality

The horizontal branch curve of the cubic root field has normalization
parameter `tau`, with

```
T=tau^2+c_0,
V=-B0(T)+2A(T)ell(T)tau.                    (14)
```

It is genuinely ramified: at the generic critical point the derivative
has a simple zero and the two critical values are distinct, giving
index two. Its function field is `C(T,tau)` because `A ell!=0`.
An embedded root field is intermediate in `L`, so Section 1 gives a
nonconstant polynomial parametrization of this branch curve.

Along it `tau` is a rational function satisfying the monic equation
`tau^2=ell(T)` and is therefore polynomial in the affine parameter.
It is nonconstant. If the rational function in (14) had a finite
pole in `tau`, this nonconstant polynomial would attain that value,
giving a pole in the polynomial target coordinate `V`. Hence (14)
is polynomial in `tau`.

Taking even and odd parts in
`C[tau]=C[T] direct_sum tau C[T]` yields

```
B0 in C[T],       A ell in C[T].            (15)
```

Thus `A` has no possible finite pole except a simple pole at `ell=0`.

### 5.2 Source valuations force a polynomial cube scaling

At a zero `c` of `A` of order `n>0`, take a prime of the square-free
source fibre `f=c`. The numerator `g+B0(f)` has valuation zero, by
dominance of this divisor onto the target line. The rational root
`Y` satisfies

```
Y^3-3ell(f)Y+(g+B0(f))/A(f)=0.              (16)
```

The last term has valuation `-n`. An integral valuation for `Y` must
be negative, and its cubic term is strictly below the linear term.
Cancellation requires `3v(Y)=-n`. Hence every zero order of `A` is
divisible by three.

A simple pole of `A` at `ell=0` is impossible. There the constant
term of (16) has valuation one and `ell` has valuation one. If
`v(Y)<=0`, the cubic term is the unique minimum; if `v(Y)>=1`, the
constant term is the unique minimum. No integral valuation can solve
the equation. In view of (15), `A` has no poles at all. Consequently

```
A=C(T)^3 for a nonzero polynomial C,
```

after absorbing a complex cube root of its scalar coefficient.
Now `W=C(f)Y` is a rational root of the monic equation

```
W^3-3ell(f)C(f)^2 W+g+B0(f)=0.              (17)
```

Normality of `B` makes `W` polynomial. Differentiating gives
`J(f,W)[3W^2-3ell(f)C(f)^2]=-J(f,g)`. The polynomials `f,W` are
algebraically independent, since (17) expresses `g` in them. The
bracket is therefore nonconstant and cannot divide a nonzero constant.
This contradiction completes case (12), and hence all cubics.

## 6. Scope and verification

The plane branch-curve argument does not preserve arbitrary finite
target covers, nor show that arbitrary intermediate fields are Galois.
The full-group assertions above concern the certified parameter covers
or the explicitly preserved quotient criterion. The transitive-inertia
criterion alone gives irreducibility. No parameter-cover source is
replaced by `A2`, and no source Hessian invariance is used.

General higher-degree coefficient variation outside these field and
inertia certificates, branch curves that do admit polynomial
parametrizations, uncovered output geometries, both plane endpoints and
unrestricted HC4/JC2 remain open. The earlier all-dimensional pencil
theorems retain their separate strengths and exact equality boundaries.

Run [the exact replay](verify_branch_parameters.py). Its
[receipt](receipt-branch-parameters.json) checks finite controls for the
unrestricted proofs above; the
[validation](validation-branch-parameters.json) separates scoped checks
from unrelated whole-repository failures. The
[source record](sources/branch-parameter-import.json) lists all retained
proof inputs. Formal verification, independent replay, external review
and novelty priority are not claimed.
