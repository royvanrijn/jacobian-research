# Integral-point quadratic recognition and a K3 lifting obstruction

Authority: `EC-CURVE302-INTEGRAL-SHELL-K3-LIFT-GATE`.

A calibrated exact quadratic-form test recovers the withheld generic
rank-17 height matrix on a known family. Applied to the existing candidate
spaces of 302, it does not recover a parent. Four selected point sets
cannot lie on any nonzero constant-norm quadratic shell. Two others admit
positive forms, but **all 144 positive even integral forms in the declared
finite boxes fail a K3 lifting condition**.

The obstruction uses the previously proved saturation of 302's displayed
group at 2. It separates compatibility with point coordinates from the
existence of a generic height lattice. Other height forms, point selections,
subspaces and elliptic surfaces with Euler characteristic above two remain
open.

## Exact recovery on the native control

Specialize the certified [11952 family](../../elkies-k3/R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md)
at `t=10^40+7`, and put the resulting curve in Sage's integral model.
The selector receives specialized point coordinates and numerical heights.
It LLL-reduces the seventeen generators, then selects up to 4,000 rays
within `7/5` times the numerical minimum, with a 50,000-ray storage cap.
The actual selected pool has 1,313 rays. Exact group law and denominator
checks retain 572 integral points.

For coefficient vectors v, solve the linear system

`v^t G v = c`,

with symmetric G and an unknown constant c. The resulting kernel has
dimension one and `c != 0`. Normalize to `c=4`. Only then compare the
recovered form with the known generic height matrix, transported to the
selected basis. They agree exactly; its determinant is 948.

This is one native-family calibration. It is not a demonstrated robust
recognizer for every family or for contaminated point sets. Numerical
height boundaries are not interval-certified. Every retained point,
quadratic relation and final matrix comparison is exact.

## Target selections and finite coefficient boxes

The same height multiple and caps are used on the existing primitive
candidate spaces in the displayed group D of 302. The point model is
the public integral Weierstrass model. The selected integral vectors
span their respective candidate spaces.

| Space | Rank | Selected rays | Integral rays | Quadratic-kernel dimension | Nonzero constant possible? |
|---|---:|---:|---:|---:|---|
| Existing core | 17 | 1,001 | 401 | 0 | No |
| 32 | 14 | 225 | 103 | 15 | Yes |
| 58 | 14 | 299 | 113 | 6 | No |
| 59 | 15 | 410 | 161 | 6 | No |
| 60 | 16 | 563 | 216 | 6 | No |
| 86 | 16 | 447 | 157 | 10 | Yes |

In rows 58–60 every kernel element has `c=0`, so the surviving homogeneous
quadratic relations cannot be positive height forms.

For rows 32 and 86, real positive-definite forms exist. A numerical
positivity calculation supplies a centre in each affine kernel. Freeze
the floor/ceiling integer choices for each free coefficient: 16,384 forms
for row 32 and 512 for row 86. These boxes were chosen **adaptively after
the target kernel tests**. Their limits are recorded in the
[portable input](../../artifacts/generated-results/elliptic-curves/curve302_integral_shell_inputs_v1.json).
No completeness outside those coefficient boxes is asserted.

Exhaustive exact positive-definiteness checks leave 142 and two positive
even integral forms respectively. Each gives norm four to every selected
integral vector, and each contains an explicitly exhibited norm-two word.
The following theorem excludes all 144 as generic K3 height forms of the
specified primitive subgroups specializing to 302.
The height-two argument treats K3 parents with reducible fibres as well.

## Height-two obstruction

We use Shioda's height formula, Shioda–Tate and primitive discriminant-group
duality for K3 lattices; see
[Schütt–Shioda, sections 6, 11 and 12](https://arxiv.org/pdf/0907.0298).
The saturation argument below is specific to this lifting problem.

Let an elliptic K3 over Q(t) have r independent rational sections, with
`r >= 14`, and suppose a nonzero section P has height two. Geometric
Picard rank is at most 20, so the total fibre-root rank is at most
`18-r <= 4`. The possible root factors of rank at most four have maximal
diagonal corrections

| Root factor | Rank | Maximum correction |
|---|---:|---:|
| A1 | 1 | 1/2 |
| A2 | 2 | 2/3 |
| A3 | 3 | 1 |
| A4 | 4 | 6/5 |
| D4 | 4 | 1 |

Since `h(P)=4+2(P.O)-sum(corrections)=2`, the only possibility is four
A1 factors, with P disjoint from O and meeting every nonidentity
component. Thus `r=14`, geometric MW rank is exactly 14 and geometric
Picard rank is 20. In particular, height two is impossible at rank 15
or higher. The four A1 fibres may be of type I2 or III; rationality of
their individual base values is not assumed.

The same height bound excludes geometric torsion. Since fourteen rational
sections span the geometric MW group over Q and there is no torsion,
every geometric section is rational: a suitable multiple is rational,
and uniqueness of division forces Galois invariance.

Now let M be the subgroup generated by the proposed fourteen sections.
Assume its specialization image is primitive in D. The established
`D intersect 2E302(Q)=2D` implies that M is saturated at 2 in the generic
MW group. Indeed, a putative half specializes into D, then into the
primitive image of M. Specialization is injective here: it is injective
on the rank-fourteen subgroup, the full generic rank is fourteen, and
generic torsion is zero. The half must already belong to M.

Suppose the proposed height matrix is even integral. Each basis section
has even height. Its number k of intersections with the four nonidentity
A1 components satisfies

`even = 4 + 2(S.O) - k/2`,

so `k=0` or `k=4`. In the lattice L generated by F, O, the four geometric
root components and M, the four root rows of the intersection matrix are
therefore identical modulo two. Their three independent differences lie
in its radical. The matrix is alternating modulo two and has even size
20, so its radical dimension is even and at least four.

The quotient `NS/L` is the full geometric MW group modulo M. Saturation
at two makes its order odd. Hence the full NS discriminant group also
needs at least four generators at two. But a primitive rank-20 NS lattice
in the unimodular rank-22 K3 lattice has the same discriminant group as
its rank-two transcendental complement, which needs at most two
generators. This contradiction excludes the proposed lift.

Thus a primitive image in the 2-saturated displayed group of 302 cannot
lift to an even integral rank-at-least-14 K3 section lattice containing
a height-two word. This statement does not exclude non-K3 surfaces.

## Replay and limitations

The [certificate](../../artifacts/generated-results/elliptic-curves/curve302_integral_shell_lift_v1.json)
contains the recovered control matrix, target kernel outcomes, all 144
positive forms and height-two witnesses, the correction-profile table,
and the universal mod-two radical vectors. It also checks each candidate
embedding's primitive Smith factors and the existing rank-31 binary
certificate underlying saturation of D.

```sh
sage -python elliptic-curves/cas/verify_curve302_integral_shell_lift.sage
```

The replay recomputes the bounded selections, exact integrality, quadratic
kernels and every coefficient box. SageMath 10.9, one worker, 120 seconds.
The initial interface failures and numerical positivity pilots remain
under `artifacts/local/`; they contribute no exclusions.

Specialized integrality alone does not imply generic disjointness from
the zero section. These point pools may mix different generic heights.
Neither failure of a constant-norm test nor the finite integer-box result
excludes all parents of a candidate space. The explicit surface, full
generic basis and specialization-to-302 objective remains open.
