# The determinant-500 and determinant-750 rootless rows are arithmetically excluded

Date: 2026-09-04.

Status: **PROVED, conditional only on the named standard theorem inputs**.

## Result

The two rootless rank-17 catalogue rows

```text
K3-04b86146cc6b284b   determinant 500
K3-10a14a46c14b3150   determinant 750
```

are both `ARITHMETICALLY_EXCLUDED`. Neither admits a characteristic-zero K3
over `QQ` with the selected saturated Neron--Severi lattice and all nineteen
divisor classes rational. Consequently neither exact rootless frame can give
a saturated arithmetic MW17 group over `QQ(t)`, and no equation work is
authorized for either row.

These are two decisions in the existing rootless-MW17 arithmetic queue, not
only low-determinant rows from the larger `T`-first catalogue.

## Literal lattices and discriminant forms

Both rows belong to the family

```text
T_N = [ 0    0    5 ]
      [ 0  10N    0 ] = U(5) + <10N> = 5 (U + <2N>),
      [ 5    0    0 ]

N=2: det(T_N)=-500,    A_T=Z/5 + Z/5 + Z/20;
N=3: det(T_N)=-750,    A_T=Z/5 + Z/5 + Z/30.
```

In the dual basis `e0/5,e2/5,e1/(10N)`, the discriminant bilinear Gram is

```text
[  0   1/5      0   ]
[ 1/5   0       0   ]    in QQ/ZZ,
[  0    0    1/(10N)]
```

and its diagonal gives the quadratic values modulo `2ZZ`. The factor-five
similarity to `U+<2N>` preserves the rational orthogonal group but not the
integral marking problem.

The catalogue already contains one exact rootless frame for each selected
surface. Their Gram hashes are respectively

```text
60c96c1674d08624622723f7d7775cc74811e7e957723291ef89c2e7efaa3ba1
faf65e8ac34588dd5e4537069bae9b1294c106fee4887166e9029ca1231ee201.
```

Both have rank 17, minimum squared norm four, and zero roots. This geometric
fact is not used to infer arithmetic realizability; it only explains why
these two `T` rows occur in the MW17 candidate classifier.

## Primitive and literal Clifford orders

For the primitive similarity form `H_N=U+<2N>`, an exact split embedding of
the even Clifford order in `M2(QQ)`, in the basis
`1,e0e1,e0e2,e1e2`, is

```text
[1 0]  [0 N]  [1 0]  [0 0]
[0 1], [0 0], [0 0], [1 0].
```

Its reduced trace pairing is

```text
[2  0  1  0]
[0  0  0  N]
[1  0  1  0]
[0  N  0  0],
```

with reduced discriminant `N`. Its integral matrices are exactly those with
`N | B`, and its determinant-one group is `Gamma^0(N)`, conjugate to
`Gamma_0(N)`. Thus the primitive-similarity curves are the genus-zero curves
`X_0(2)` and `X_0(3)`. They are only coarse quotients.

For the literal lattice `T_N`, the embedded basis is instead

```text
[1 0]  [0 5N]  [5 0]  [0 0]
[0 1], [0  0], [0 0], [5 0].
```

The literal reduced discriminant is `125N`: `250` for determinant 500 and
`375` for determinant 750. Membership is equivalent to

```text
5N | B,       5 | C,       A = D (mod 5).
```

For determinant-one matrices these are the primitive-order conditions plus
projective identity modulo five.

## Full discriminant action and stable kernel

Conjugation on the trace-zero Clifford subspace gives the spin action on the
literal `T_N`. For both `N=2,3`, joint exact reduction on Sage generators of
`Gamma_0(N)` proves:

```text
image in SL_2(F_5)             120
spin image in O(A_T)            60
kernel in SL_2(F_5)           {+I,-I}
spin image                    PSL_2(F_5) = A5.
```

The spin-image element-order histogram is

```text
1:1, 2:15, 3:20, 5:24.
```

The norm-one spin subgroup is not all of `O^+(T_N)`. The standard normalizer
description for a split prime-level Eichler order adds the Fricke coset; its
discriminant image enlarges `A5` to `S5`, with histogram

```text
1:1, 2:25, 3:20, 4:30, 5:24, 6:20.
```

Adding the negative reflection in `(1,0,-1)` and central inversion gives the
full image `S5 x C2`, of order 240, with histogram

```text
1:1, 2:51, 3:20, 4:60, 5:24, 6:60, 10:24.
```

The checker tests every action in each of the three non-spin cosets; none is
stable. Hence on the marked period component

```text
O^+(T_N)^* = Gamma_0(N) intersection +/-Gamma(5).
```

This has index 60 over the primitive-similarity curve.

## Exact marked curves

Conjugate by `D=diag(5,1)` and put `M=25N`. The stable group becomes

```text
Gamma_H(M),
H = {a in (Z/MZ)^* : a = +/-1 (mod 5)}.
```

This distinction matters. Since `(Z/5Z)^*` has four elements, the diagonal
condition is not automatic after conjugation: the exact marked curve is a
degree-two cover of `X_0(M)`, not merely `X_0(M)`.

The exact signatures are

| determinant | exact marked curve | PSL2 index | genus | elliptic points | geometric cusps | rational cusps |
|---:|---|---:|---:|---:|---:|---|
| 500 | `X_H(50)` | 180 | 4 | none | 24 | `1/25, 3/50, 2/25, Infinity` |
| 750 | `X_H(75)` | 240 | 9 | none | 24 | `2/75, 1/25, 2/25, Infinity` |

The cusp fields are decided by the exact cyclotomic Galois action on every
`Gamma_H` cusp. In each case exactly the four displayed cusps are fixed over
`QQ`.

The useful forgetful maps are

```text
X_H(50) -> X_0(50),    degree 2;
X_H(75) -> X_0(75),    degree 2.
```

The target curves have genera two and five. A rational noncuspidal point on
either target would give a rational cyclic isogeny of degree 50 or 75. Both
degrees are absent from the complete Mazur--Kenku list

```text
1,...,19, 21, 25, 27, 37, 43, 67, 163.
```

Therefore each target has only its four rational cusps, and any rational
point upstairs is one of the four rational cusps already listed. For each
exact marked curve the point accounting is consequently

```text
cuspidal rational points                   4
noncuspidal rational points                0
rational CM points                         0
rational non-CM points                     0
full-curve points realizing an overlattice 0.
```

The coarse curves `X_0(2)` and `X_0(3)` are genus zero with rational cusps,
so each has infinitely many rational quotient-level points. Those points do
not carry the full discriminant marking, and no noncuspidal one lifts to a
rational point on the exact marked cover. Thus they are kept separate rather
than miscounted as full marked or overlattice points. A cusp is not a K3
period point, so neither row supplies even a candidate for the saturated
rational marking.

## Arithmetic certificates

The required decisions are exactly

```text
determinant 500 / K3-04b86146cc6b284b:  ARITHMETICALLY_EXCLUDED
determinant 750 / K3-10a14a46c14b3150:  ARITHMETICALLY_EXCLUDED
```

The checker
[`scripts/certify_det500_det750_qq_marking_obstructions.sage`](scripts/certify_det500_det750_qq_marking_obstructions.sage)
reconstructs both discriminant forms, both primitive and literal Clifford
orders, every spin-generator action, the `A5`, `S5`, and `S5 x C2` images,
the stable congruence groups, signatures, cusp Galois orbits, and quotient
maps. Its generated output is
[`../artifacts/generated-results/elkies-k3-det500-det750-qq-marking-obstructions-v1.json`](../artifacts/generated-results/elkies-k3-det500-det750-qq-marking-obstructions-v1.json).

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_det500_det750_qq_marking_obstructions.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_det500_det750_qq_marking_obstructions.sage --check
```

The checker does not reprove the Mazur--Kenku classification, the standard
split prime-level Eichler normalizer theorem, or the rank-three
ternary-spin/fully-marked-K3 period correspondence. It makes no statement
about other ternary genera of determinants 500 or 750, or about markings over
larger number fields.

The derivation consumes only the lattice-foundry catalogue and its
transcendental-arithmetic ledger. It does not read or update the curve-356,
curve-385, or frozen prospective-experiment artifacts.

The rational cyclic-isogeny list is summarized in Theorem 1.1 of Banwait,
Najman, and Padurariu,
[*Cyclic isogenies of elliptic curves over fixed quadratic fields*](https://arxiv.org/abs/2206.08891).

## Foundry consequence

The rootless-MW17 classifier now has one realized positive control, five exact
exclusions, and sixty `UNKNOWN` rows. The global arithmetic-first queue has
six exact exclusions, one realized positive control, and 820 unresolved
rows. The equation-agent handoff remains empty.

No equation, carrier-receptivity, or record-search work is triggered: both
selected candidates fail at the full rational-marking/period-curve gate.
