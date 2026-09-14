# Common-quartic singularities and the fixed-parent boundary

The common-quartic equations admit explicit **genus-zero and positive-rank
genus-one controls**, each with two independent new sections. These controls
construct new `24I1` K3 parents. Their inherited ranks are **UNKNOWN**; they
are not a solution on a retained MW17 parent and supply no rank19 claim.

On each of the four retained MW17 equations, there is also an exact negative
result: **every algebraic constant abscissa gives normalization genus at
least4**. The proof covers the entire constant parameter line, not a bounded
search. Nonconstant quartic abscissas remain outside this exclusion.

Thus the [high-rank correlated-cover objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains **OPEN**. No lower bound above one for its whole fixed-parent quartic
chart has been proved. The calculations below concern the different system
from [the branch-first note, section4](BRANCH_FIRST_CORRELATED_QUADRATIC_COVERS_2026-09-14.md#4-a-different-system-with-the-common-quartic-built-in);
no Mestre auxiliary function is varied.

## 1. The full local discriminant and singularity equations

Fix a short K3 equation with `deg A<=8`, `deg B<=12`, squarefree degree24
`Delta=4A^3+27B^2`, and smooth infinity. Write

```
x(t)=u0+u1*t+u2*t^2+u3*t^3+u4*t^4,
F(t)=x(t)^3+A(t)*x(t)+B(t).
```

The image of a nonsplit cover is the curve `y^2=F(t)` in the class
`2O+4F_fibre`, of arithmetic genus5. Its finite singularity incidence is

```
F = 0,
F' = (3*x^2+A)*x' + A'*x+B' = 0.                       (1)
```

An ordinary geometric node has in addition

```
F'' = (3*x^2+A)*x'' + 6*x*(x')^2
      + 2*A'*x' + A''*x+B'' != 0.                     (2)
```

Higher contact is imposed by further vanishing derivatives. These equations
are a compact exact incidence presentation of the discriminant locus;
expanding its five-parameter resultant is not needed to define it.
They do not by themselves impose four nodes or a shared quadratic field.

At a root of multiplicity `m`, over an algebraic closure a unit can be
absorbed into `y`, giving `y^2=z^m`. Its delta invariant is `floor(m/2)`;
it branches on normalization exactly when `m` is odd. This includes nodes
(`m=2`), cusps (`m=3`) and higher contacts. Restore infinity using

```
F_h(T,Z)=Z^12*F(T/Z),          m_infinity=12-deg(F).
```

For a geometrically connected double cover, let `b` count the geometric
roots of `F_h` of odd multiplicity. Then

```
delta_total = (12-b)/2,
g_normalization = (b-2)/2 = 5-delta_total.              (3)
```

Consequently genus1 requires total delta4 and genus0 requires delta5.
A vanishing discriminant alone supplies neither conclusion. When `b=0`,
the geometric cover splits; a nonsquare rational scalar can instead give
a constant-field extension, which is also outside the objective. Equation
(3) uses the [Riemann--Hurwitz formula](https://stacks.math.columbia.edu/tag/0C1B)
and the local normalization calculation, not a count of rational singular
points. The two branches above a rational node need not be rational.

## 2. Repeated quartic factors give a valid genus-zero boundary

Keep the polynomial identities but allow the common binary quartic `D`
to become singular:

```
x_i^3+A*x_i+B = D*r_i^2,
deg x_i<=4, deg r_i<=4, r_i!=0, x_1!=x_2,
gcd(D,Delta)=1.                                        (4)
```

Work with binary forms to retain infinity. Remove the actual square factor,
including its constant squareclass:

```
D=e^2*d,    d squarefree of binary degree b in {2,4}.
```

Now a *single* normalization is `C: w^2=d` and

```
P_i=(x_i,e*r_i*w).                                    (5)
```

Degree4 gives genus1. A quartic multiplicity pattern `(2,1,1)` or `(3,1)`
gives degree2 and genus0. Patterns `(2,2)` and `(4)` give no geometric
quadratic field. In particular, a repeated quartic must not be rejected
merely because the original smooth-quartic sufficient chart excluded it.
Its normalization and its rational-point obligation must be checked instead.

At a place where `ord(r_i)=m` and `ord(D)=n`, the image has local exponent
`2m+n`, delta `m+floor(n/2)` and branch parity `n mod2`. Four simple zeros
of `r_i` away from a squarefree `D` give four nodes. One further double
root of `D`, away from `r_i`, gives a fifth node and genus0.

The height argument extends to this boundary. Branching is only at smooth
parent fibres. The pullback has `chi=4` and only irreducible fibres.
The degree bounds make both points integral in the infinity chart as well
as at finite places, so `P_i.O=0` and each height is8. The
[Shioda height formula](https://arxiv.org/html/0907.0298v3) also makes the
torsion group trivial: any nonzero section has height `8+2(P.O)>0`.
If two height8 points were rationally dependent, they would be equal up to
sign modulo torsion, forcing their abscissas to agree. Thus distinct
abscissas in (4) suffice for independence. Deck conjugation negates both
points and fixes the **full** inherited Mordell--Weil group. Its invariant
and anti-invariant rational subspaces intersect trivially. Independence
therefore holds modulo the full inherited group, without knowing its rank.

For symbolic elimination it is useful to solve for the parent coefficients
first. With `h=x_1-x_2!=0`, equations (4) are equivalent to

```
A = D*(r_1^2-r_2^2)/h - (x_1^2+x_1*x_2+x_2^2),
B = D*(x_1*r_2^2-x_2*r_1^2)/h + x_1*x_2*(x_1+x_2).    (6)
```

Polynomial divisibility, coefficient degrees and the fixed values of `A,B`
are still equations, not optional constraints. Formula (6) parametrizes
incidences while permitting the parent to move; using it to construct a
different parent does not solve a prescribed MW17 equation.

## 3. Two explicit controls with all base and independence gates closed

Set `s=t^4+1` and, for either quartic below, define

```
A = D*(2*s+1)-1,          B = D*s^2.                    (7)
```

This is the exact subfamily `x_1=0`, `x_2=1`, `r_1=s`, `r_2=s+1` of (6).
If `D=e^2*d`, then on the one quadratic cover `w^2=d` the two sections and
their sum are

```
P = (0, e*s*w),
Q = (1, e*(s+1)*w),
P+Q = (D-1, e*(1-D-s)*w).                              (8)
```

The three equations and the chord addition are literal polynomial
identities. All three have height8, so

```
Gram(P,Q) = [[8,-4],[-4,8]],       determinant = 48.    (9)
```

The inherited rank is not needed for this gain of two, but it is needed
before advertising a total lower bound19.

### Four nodes: an elliptic base of positive rank

Take `D=d=t^4+t+1`, `e=1`. The parent is

```
A = 2*t^8+2*t^5+5*t^4+3*t+2,
B = t^12+t^9+3*t^8+2*t^5+3*t^4+t+1.
```

The quartic and both `s,s+1` are squarefree and mutually coprime where
required. Each image has four ordinary geometric nodes and normalization
genus1. The following rational maps identify the base with an elliptic curve:

```
C: w^2=t^4+t+1,
E_C: Y^2=X^3-4*X+1,
X=2*(w+t^2),        Y=2*t*X+1,
t=(Y-1)/(2*X),     w=X/2-t^2.                          (10)
```

The inverse is written on `X!=0`; it extends on the smooth projective curves.
The point `R=(0,1)` on `E_C` satisfies

```
2R=(4,7),          3R=(-7/4,13/8).
```

The nonintegral coordinate of `3R` proves that `R` is nontorsion by
[Lutz--Nagell](https://eudml.org/doc/150016) on the integral short equation. Thus this *same* quadratic
base has infinitely many rational points and infinitely many rational
`t` values. No BSD or numerical rank calculation is used.

### Five nodes: a rational base from a double quartic factor

Take `D=t^2*(t^2+1)`, `e=t`, `d=t^2+1`. The parent is

```
A = 2*t^8+2*t^6+3*t^4+3*t^2-1,
B = t^12+t^10+2*t^8+2*t^6+t^4+t^2.
```

Each image has four nodes at the zeros of its `r_i`, and a fifth at `t=0`.
The two points above that last node lie on the smooth conic normalization.
An explicit parametrization is

```
t=2*z/(1-z^2),       w=(1+z^2)/(1-z^2).                (11)
```

This is degree two over the **original t-line**, not a degree-four
compositum or a renamed intermediate quotient. Equations (8) and (9)
apply unchanged with `e=t`.

### K3, smoothness and specialization checks

For both controls `deg A=8`, `deg B=12`. Their leading discriminant
coefficient is `4*2^3+27*1^2=59`, so infinity is smooth. Exact reduction
modulo1009 proves that the degree24 discriminant is squarefree and coprime
to `D`. Hence the parents are K3 surfaces with precisely `24I1` fibres.
The node factors and residual branch factors pass the same squarefreeness
and coprimality checks. Evaluating `A^3/B^2` at0 and1 proves nonconstant j;
no point or parameter search enters these checks.

The usual specialization theorem then supplies infinitely many rational
base values where the two new directions remain independent together with
any fixed independent inherited subgroup. The total inherited rank of
either new parent has not been computed. These are constructive low-genus
controls and a counterexample to a genus-above-one bound for the system
with **unrestricted moving parent coefficients**. They do not contradict
or settle a bound restricted to the four fixed MW17 parents.

## 4. A complete symbolic exclusion on the four retained constant-x lines

Now keep `A,B` fixed at any of published R17, direct11952 alternate Q80,
the recovered Curve302 parent, or X1092 class1. For any algebraic constant
`c`, set `F_c(t)=c^3+A(t)*c+B(t)`. Define

```
H(t)=B*(A')^3-A*B'*(A')^2-(B')^3,
N(c)=Res_t(H(t), c*A'(t)+B'(t)).                       (12)
```

These are exact characteristic-zero expressions; the large integer
coefficients of the resultants are left unexpanded. Their bounds are
`deg H<=33`, `deg_c N<=33`. The certificate modulo1009 proves on **each**
parent:

```
deg H=33,          gcd(H,A')=gcd(H,H')=1,
deg N=33,          gcd(N,N')=1.                        (13)
```

The leading coefficients remain units. Therefore these degrees and
nonvanishing resultants hold in characteristic zero as well. In particular,
there is no unexamined algebraic constant with a denominator at1009: (13)
certifies whole rational polynomials, not reduction of a guessed parameter.

Here is why this proves a genus bound. If `F_c(t)=F_c'(t)=0`, then `A'`
cannot vanish: simultaneous `A'=B'=0` would contradict `gcd(H,A')=1`.
Hence

```
c=a(t)=-B'(t)/A'(t),            H(t)=0.
```

Conversely each root of `H` gives such a singular fibre of this abscissa
pencil. The norm `N` is a nonzero scalar times
`product_(H(alpha)=0) (c-a(alpha))`. Its squarefreeness proves that all33
critical abscissas are distinct. At a root of `H`, differentiating
`H=(A')^3 F_(a(t))(t)` gives

```
H' = -(A')^2*(3*a(t)^2+A)*F_c''(t).                  (14)
```

By (13), the left side is nonzero. Thus `F_c''!=0`: every multiple root
is an ordinary double root, and no `c` has two such roots. Every `F_c`
has degree12 with leading coefficient `B_12`, independent of `c`, so
there is no extra singularity at infinity. Formula (3) now proves

```
g_normalization(F_c) = 4 at the33 critical values,
g_normalization(F_c) = 5 at every other algebraic c.   (15)
```

This excludes any common-quartic solution on these parents that uses even
one constant `x_i`. It closes the simple control subfamily (7) there,
without excluding any nonconstant polynomial `x_i(t)`.

## 5. Evidence, replay and the remaining fixed-parent problem

The [frozen input](../artifacts/generated-results/elkies-k3-common-quartic-singularities-v1/input.json)
binds the previously retained generic coefficient packet, four parent names,
two control inputs, prime1009, and limits of40 CPU seconds and4GiB.
The [result](../artifacts/generated-results/elkies-k3-common-quartic-singularities-v1/result.json)
has a checkpoint for every parent and control. The
[independent replay](../artifacts/generated-results/elkies-k3-common-quartic-singularities-v1/independent-replay.json)
uses only Python integers and fractions. Its modular resultants use the
Euclidean recurrence, independently of Sage's resultant implementation.
It reconstructs `N mod1009` from34 values with a **proved degree bound33**;
these are polynomial interpolation data, not a finite parameter exclusion.
It also checks literal squareclasses, all three section equations, addition,
infinity bounds, fibre/node guards, both rational base maps, and the
nontorsion point multiples.

The eight regression tests cover corrupt norms and coefficients, incomplete
parent coverage, a changed constant squareclass, duplicate sections,
infinity and higher singularities, split covers, false rank promotions and
resultant signs/leading coefficients. Written local-geometry, height and
specialization arguments are not formal verification.

```
python3 research/elkies-k3/scripts/verify_common_quartic_singularity_locus.py
python3 -m unittest discover -s research/tests -p 'test_common_quartic_singularity_locus.py' -q
```

For discovery replay, run the producer with `freeze --output DIR`, then
`run --output DIR`, both under `sage -python` with
`OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1`. Outputs refuse overwrites.
The [execution receipt](../artifacts/generated-results/elkies-k3-common-quartic-singularities-v1/execution.json)
retains the initial bad-prime probe, successful witness selection, preview
and final runs, and a corrected hand-written test expectation. No full
quartic elimination, bisection census, descent or rational-point campaign ran.

The remaining question is precise: solve (1)--(4) on a **fixed** MW17
equation with both abscissas nonconstant, or prove that the full
coefficient incidence has no admissible rational point. A complete negative
proof must include nonordinary contacts, repeated quartic factors, infinity,
and all coefficient-chart boundaries; counting equations or testing a
bounded set of coefficients cannot supply it. For any positive solution,
the normalization must additionally have a rational parametrization or a
certified rational nontorsion point. None of these fixed-parent obligations
has been declared complete here.
