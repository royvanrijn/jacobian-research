# Dissection of the `E6+A1` Picard-rank-19 K3 — 2026-09-02

<!-- status-consumer: EC-K3-RES-QBC-E6A1-RHO19 7103fa2a1a4e7ba2 -->

## Executive result

The new quadratic-base-change family has now been pinned as an integral K3
lattice, not only as a rank and height calculation.  For very general
parameter its Neron--Severi and transcendental lattices are

```text
rho = 19,
disc group NS = Z/3 + Z/12,
|disc NS| = 36,
T(X) = U(3) + <4>.
```

The displayed sections form the full Mordell--Weil lattice.  The only
possible proper even overlattices have index three, and an exact good-fibre
test eliminates both of them.  The resulting determinant-36 positive frame
is included in the generated certificate as a direct input for elliptic-
neighbor searches.

Four parameter values give forced singular K3 surfaces with exact
transcendental lattices of determinants `12`, `15`, `60`, and `72`.  They
come from three different mechanisms: additive-fibre enhancement with
two-torsion, multiplicative collision with three-torsion, and branching over
an `I1` fibre.

The exact replay is
[`scripts/certify_e6a1_rho19_k3_dissection.sage`](scripts/certify_e6a1_rho19_k3_dissection.sage),
and its machine-readable result is
[`../artifacts/generated-results/elkies-k3-e6a1-rho19-k3-dissection-v1.json`](../artifacts/generated-results/elkies-k3-e6a1-rho19-k3-dissection-v1.json).

## Equation and the true moduli parameter

Recall

```text
D      = 3*k^2-4,
c      = 2*k/D,
lambda = -(k^2-4)*D/4,
H      = 1-t^2.
```

The K3 equation is

```text
Y^2 = X^3
    + H^3*(lambda-3*H)*X
    + H^4*(c^2*lambda^2+lambda*H-2*H^2).
```

The unmarked equation depends only on

```text
s = k^2,
c^2 = 4*s/(3*s-4)^2,
lambda = -(s-4)*(3*s-4)/4.
```

Thus the coarse one-dimensional parameter is `s`, while `k` is the double
cover that rationalizes the invariant section.  Under `k -> -k`,

```text
P0 -> -P0,
P1 ->  P1.
```

Consequently the generic Mordell--Weil rank is one over `QQ(s)(t)` and two
over `QQ(k)(t)`.  This field-of-definition distinction is essential for
arithmetic specialization: a rational value of `s` need not make `P0`
rational.

If `z1,z2` are the squared `t`-coordinates of the two pairs of residual
`I1` fibres, they are the roots of

```text
z^2 - S(s)*z + P(s),
```

where

```text
S(s) = (27*s^4-54*s^3-126*s^2+656*s+544)
       /(9*(3*s+4)^2),

P(s) = (2*s+1)*(9*s-20)^2/(9*(3*s+4)^2),

S(s)^2-4*P(s)
     = (s-4)^2*(s^2+4*s/3+16/9)^3
       /(9*(s+4/3)^4).
```

The earlier invariant `P(s)` is sufficient to prove non-isotriviality, but
it is not a complete moduli coordinate: the two distinct singular K3 points
`s=20/9` and `s=-1/2` both have `P(s)=0`.

## Generic fibres and sections

The rational elliptic surface has discriminant

```text
-16*u^2*Q(u),
Q(u) = 27*c^4*u^2 + (54*c^2+4)*u - (108*c^2+9).
```

The quadratic map is branched at the `I2` fibre `u=0` and the generically
smooth fibre `u=lambda`.  Hence the K3 has

```text
t =  1,-1 : IV*,
t = infinity: I4,
four other points: I1,
root lattice: 2E6+A3.
```

The section component data are:

| section | first `IV*` | second `IV*` | `I4` |
|---|---:|---:|---:|
| `P0` | nonzero, correction `4/3` | nonzero, correction `4/3` | component `2`, correction `1` |
| `P1` | identity | identity | component `2`, correction `1` |

Thus

```text
<P0,P0> = 1/3,
<P1,P1> = 3,
<P0,P1> = 0.
```

Deck-character orthogonality gives the cross height.  Since the local cross
correction at `I4` is one, the corresponding integral section curves satisfy
`P0.P1=1`.

## Integral Neron--Severi lattice

Use the divisor basis

```text
O,F,E6a_1,...,E6a_6,E6b_1,...,E6b_6,A3_1,A3_2,A3_3,P0,P1.
```

The full 19-by-19 Gram matrix is serialized in the certificate.  Its exact
invariants are

```text
signature             (1,18),
determinant            36,
Smith invariants       1^17,3,12,
discriminant group     Z/3 + Z/12.
```

In compatible Smith generators its finite quadratic form is

```text
[ 4/3   1/3 ]
[ 1/3  -1/4 ]
```

with diagonal entries read modulo `2Z` and off-diagonal entries modulo `Z`.
The nonzero isotropic elements are

```text
(0,4), (0,8), (1,4), (2,8).
```

They form two order-three subgroups.  Therefore index three is the only
possible saturation defect.  At `(k,t)=(1,3)` the good elliptic fibre is

```text
y^2 = x^3 - 11904*x - 490496,
```

and none of `P0`, `P1`, `P0+P1`, or `P0-P1` is divisible by three in its
rational Mordell--Weil group.  A generic third would extend over this good
fibre, so both possible index-three overlattices are excluded.  This proves
that the two displayed sections generate the full generic MW lattice and
that `|disc NS|=36` is exact.

The torsion subgroup is trivial: torsion injects into the component groups,
whose orders at `IV*` and `I4` have greatest common divisor one.

## Transcendental lattice

The orthogonal complement in the K3 lattice has rank three, signature
`(2,1)`, determinant `-36`, and discriminant form opposite to `NS`.  The
certificate gives an explicit finite-quadratic-module isometry with

```text
T(X) = U(3) + <4>,

Gram(T) =
[ 0 3 0 ]
[ 3 0 0 ]
[ 0 0 4 ].
```

Sage's exact genus enumeration returns one class for this genus, so the
discriminant-form identification determines the isometry type, not merely a
candidate in the same genus.

This also shows that the very general surface is not a Kummer surface of the
most immediate `T(A)(2)` type: the transcendental form is not globally
two-divisible.  No stronger exclusion of Shioda--Inose structures is claimed.

## The deck involution

The involution

```text
iota: (t,X,Y) -> (-t,X,Y)
```

is non-symplectic because it negates `dt wedge dX/Y`.  It fixes `P0`, negates
`P1`, and its quotient is the original rational elliptic surface.

Its fixed locus consists of the smooth genus-one branch fibre at `t=0` and
the two strict transforms of the rational components of the branched `I2`
fibre.  Locally, each node has equation `xy=r^2`; the exceptional curve is
not fixed pointwise.  Hence the fixed locus is one genus-one curve plus two
rational curves.  The standard two-elementary fixed-locus formulas give

```text
r = 12,
a = 8.
```

The parity invariant `delta` has not been independently pinned here.

## Four forced singular K3 surfaces

The following table records geometric lattices.  Rationality of individual
sections is stated separately because the equation descends to `s` but `P0`
usually requires `k=sqrt(s)`.

| `s` | mechanism | fibres | torsion / free heights | `disc NS` | `T(X)` | CM discriminant |
|---:|---|---|---|---:|---|---:|
| `0` | `IV* -> III*`; `P0` becomes 2-torsion | `2III*+I4+2I1` | `Z/2`; `[3]` | `-12` | `[[4,2],[2,4]]` | `-12` |
| `-4/3` | `I2+I1 -> I3`; pullback `I6`; `P0` becomes 3-torsion | `2IV*+I6+2I1` | `Z/3`; `[5/2]` | `-15` | `[[4,1],[1,4]]` | `-15` |
| `20/9` | smooth branch hits an `I1`; `P1` meets the new component | `2IV*+I4+I2+2I1` | trivial; `diag(1/3,5/2)` | `-60` | `[[4,2],[2,16]]` | `-60` |
| `-1/2` | smooth branch hits an `I1`; `P1` meets the identity component | `2IV*+I4+I2+2I1` | trivial; `diag(1/3,3)` | `-72` | `diag(6,12)` | `-72` |

For the last two rows, root rank `16` plus two independent sections forces
`rho=20`.  For the first two, root rank `17` plus the displayed free section
does the same.  The section lattices are saturated by exact divisibility tests
at `t=3`; at `s=-4/3`, saturation is automatic because no nontrivial square
divides `15`.  The binary transcendental forms are selected by the exact
opposite discriminant form, not only by their determinant.

The equations at all four values are defined over `QQ`, since they depend on
`s`.  However:

- at `s=20/9`, `P0` is naturally defined over `QQ(sqrt(5))`;
- at `s=-1/2`, it is defined over `QQ(sqrt(-2))`;
- at `s=-4/3`, the displayed 3-torsion section is defined over
  `QQ(sqrt(-3))`;
- at `s=0`, both the torsion section and the free section are rational.

These points are arithmetic controls for reduction, CM, and neighbor
compilers; they are not generic-rank specializations.

## Other degeneration loci

The equation

```text
9*s^2 + 12*s + 16 = 0
```

makes the two residual `I1` fibres of the rational surface collide at the
zero of both `c4` and the discriminant.  The result is type `II`, not `I2`.
The K3 fibre profile is therefore

```text
2IV* + I4 + 2II.
```

There is no root-rank increase, so Shioda--Tate alone does not force
`rho=20`.  Their Picard status remains unproved here.

At `s=4`, `lambda=0` and the displayed Weierstrass cubic is singular rather
than a K3.  The chart has a pole at `s=4/3`.  The compactified limits at
`s=4/3` and `s=infinity` have not yet been resolved; they must not be silently
identified with members of the generic fibre profile.

## Exact frame and the first chamber obstruction

Splitting the old hyperbolic plane gives a positive-definite rank-17 frame
`M` with

```text
det(M) = 36,
minimum = 2,
156 roots,
root span = 2E6+A3 of rank 15 and determinant 36.
```

The certificate contains its full 17-by-17 Gram matrix.  Every future
Jacobian fibration on this K3 must split another copy of `U`; its frame still
has rank 17, determinant 36, and the same discriminant form.  A rootless
fibration would therefore require an even positive-definite rank-17 MW
lattice of determinant 36.  Existence of such a frame is not asserted.

The complete nominal smallest degree-two layer has also been checked.  In split
coordinates `NS=U+(-M)`, consider

```text
D = e + 2*f + w,
w^2 = 4.
```

There are `17688` norm-four vectors, collapsing to `14` orbits under the
current `2E6+A3` Weyl group.  But with `O=f-e`, every one satisfies

```text
D.O = -1.
```

The zero section is fixed, and removing it changes `D` to a primitive
degree-one isotropic class.  The repeated child root data

```text
(root rank, root count, root determinant) = (15,156,36).
```

therefore belong to section presentations, not genuine quadratic neighbours.
The nominal degree-three shape `e+3*f+w`, `w^2=6`, has the same obstruction:
`D.O=-2`, and removing `2O` again leaves degree one.  The first zero-neutral
shapes are instead

```text
q=2: D=2e+2f-w, w^2=8,
q=3: D=3e+3f-w, w^2=18.
```

The complete genuine `q=2` calculation is now recorded separately in
[`E6A1_RHO19_GENUINE_Q2_NEIGHBORS_2026-09-02.md`](E6A1_RHO19_GENUINE_Q2_NEIGHBORS_2026-09-02.md).

<!-- status-consumer: EC-K3-E6A1-RHO19-GENUINE-Q2-MW3 cd4314040bb028f7 -->

## Reproduction and proof boundary

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_rational_surface_quadratic_rank_search.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_k3_dissection.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_e6a1_rho19_k3_dissection.sage --check
```

Pinned SHA-256 values after this dissection are

```text
base checker       d9cf7346dc8273e4df47db8a30b4575783944cda2fecec3f89466cc0428b8673
base artifact      888c7bba059af4cff18a8baa359b60f227575a4835d28b7764f62c42690f8d8c
dissection checker e18a9015700df415c24e86ddc4ff87790119f6284234c55c2d27accab525d9e2
dissection artifact c5da73a84f77392d3be2a657fbef1a963666da1724ede180d057a8d3dbc60ecd
```

The generic `NS`, its saturation, `T(X)`, the four singular-K3 lattices, and
the nominal-layer chamber obstruction are exact.  The following remain open:

1. the compactified limits at `s=4/3` and `s=infinity`;
2. the Picard status at the two type-`II` collision parameters;
3. a complete elliptic-fibration classification for this `NS` genus;
4. an explicit Weierstrass compilation of a certified nef low-root neighbor;
5. arithmetic rank-jump searches on non-CM rational `k` specializations.

The rational-surface chart and the general quadratic-base-change framework
come from [Kimura](https://arxiv.org/abs/1802.05195).  This note's integral
gluing, transcendental identification, singular-K3 boundary table, and
bounded neighbor result are calculations specific to the new subfamily.
