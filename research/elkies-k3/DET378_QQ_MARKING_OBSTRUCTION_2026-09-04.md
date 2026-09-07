# The split determinant-378 row has no full rational marking

Date: 2026-09-04.

Status: **PROVED, conditional only on the named standard theorem inputs**.

## Theorem

Let

```text
T = [ 0   0   3 ]
    [ 0  42   0 ] = U(3) + <42>.
    [ 3   0   0 ]
```

There is no characteristic-zero K3 surface over `QQ` with this
transcendental lattice and a full saturated rational rank-19
Néron--Severi marking. Equivalently, the split arithmetic-first row
`K3-7b71a1cc00b0c6e2` is excluded before any `NS=T^perp`, rootless-frame,
or equation work.

This is not a claim about the other determinant-378 row in the catalogue. Its
ternary lattice is anisotropic and has quaternion discriminant `21`; it does
not have coarse curve `X_0(7)`.

## Literal lattice and discriminant form

The lattice has signature `(2,1)`, determinant `-378`, and

```text
A_T = Z/3 + Z/3 + Z/42.
```

In the dual basis `e0/3, e2/3, e1/42`, its discriminant bilinear form is

```text
[  0   1/3    0  ]
[ 1/3   0     0  ]   in QQ/ZZ,
[  0    0    1/42]
```

and the quadratic values are `0,0,1/42` in `QQ/2ZZ`. Thus the factor-three
similarity

```text
T = 3 H,       H = U + <14>
```

does not preserve the marking problem even though it preserves the rational
orthogonal group.

## Primitive and literal Clifford orders

For `H`, an exact embedding of the even Clifford order in `M2(QQ)`, in the
basis `1,e0e1,e0e2,e1e2`, is

```text
[1 0]  [0 7]  [1 0]  [0 0]
[0 1], [0 0], [0 0], [1 0].
```

Its matrices are exactly the integral matrices with `7 | B`, so its
determinant-one group is `Gamma^0(7)`, conjugate to `Gamma_0(7)`. This recovers
the coarse genus-zero curve `X_0(7)`.

For the literal lattice `T`, the embedded basis is instead

```text
[1 0]  [0 21]  [3 0]  [0 0]
[0 1], [0  0], [0 0], [3 0].
```

The literal reduced trace pairing is

```text
[2   0   3   0]
[0   0   0  63]
[3   0   9   0]
[0  63   0   0],
```

with reduced discriminant `189`, not `7`. Membership is equivalent to

```text
21 | B,       3 | C,       A = D (mod 3).
```

For determinant one these are precisely the primitive-order conditions plus
projective identity modulo three.

## Spin action and stable kernel

Conjugation on the trace-zero Clifford subspace gives the exact spin action on
`T`. On Sage generators of `Gamma_0(7)`, its image in `O(A_T)` has order `12`
and element-order histogram

```text
1:1, 2:3, 3:8,
```

so the norm-one spin image is `A4`. Exact joint reduction shows that this action is
the projective reduction modulo `3`; its kernel is `+/-I` modulo `3`.

The norm-one spin subgroup is not all of the integral orthogonal group. By the
standard normalizer description for the split level-`7` Eichler order, the
remaining cosets are represented by the Fricke isometry, a negative
reflection, and central inversion. Here the Fricke isometry is minus the
reflection in the negative norm-`6` vector `(1,0,-1)`; adding its coset
enlarges the proper image to `A4 x C2`. Adding the reflection and inversion
gives the full `O^+(T)` image
`A4 x C2 x C2`, of order `48`, with histogram

```text
1:1, 2:15, 3:8, 6:24.
```

The Fricke, reflection, and inversion actions all lie outside the norm-one
`A4` image in the required cosets. Therefore none of the non-spin cosets is
stable, and on the marked period component

```text
O^+(T)^* = Gamma_0(7) intersection +/-Gamma(3).
```

This subgroup has index `12` over the coarse curve and index `96` in
`PSL2(ZZ)`.

## Exact marked curve and rational points

Put `D=diag(3,1)`. Direct congruence comparison gives

```text
D^-1 (Gamma_0(7) intersection +/-Gamma(3)) D = Gamma_0(63).
```

Hence the exact marked modular curve is

```text
X_0(63),       genus 5.
```

It has eight geometric cusps. Four are rational (`0`, `1/9`, `1/7`, and
`Infinity`); the pairs above denominators `3` and `21` have cusp rationality
order three and are not individually rational over `QQ`.

A noncuspidal rational point would give a rational cyclic `63`-isogeny.
Degree `63` is absent from the complete Mazur--Kenku list

```text
1,...,19, 21, 25, 27, 37, 43, 67, 163.
```

Thus `X_0(63)(QQ)` consists exactly of its four rational cusps. There are no
noncuspidal rational points to screen for CM, lift to a saturated rational
marking, or pass to `NS=T^perp`. The rootless-MW17 and equation gates remain
closed.

The cyclic-isogeny classification is summarized as Theorem 1.1 of Banwait,
Najman, and Padurariu,
[*Cyclic isogenies of elliptic curves over fixed quadratic fields*](https://arxiv.org/abs/2206.08891),
which records the earlier work of Mazur and Kenku.

## Exact replay and boundary

The checker
[`scripts/certify_det378_qq_marking_obstruction.sage`](scripts/certify_det378_qq_marking_obstruction.sage)
reconstructs the discriminant form, both Clifford multiplication tables and
embeddings, both reduced trace pairings, every generator action, the `A4`,
`A4 x C2`, and full `A4 x C2 x C2` images, the non-spin-coset exclusion, the congruence
conjugation, and the `X_0(63)` signature and cusps. Its output is
[`../artifacts/generated-results/elkies-k3-det378-qq-marking-obstruction-v1.json`](../artifacts/generated-results/elkies-k3-det378-qq-marking-obstruction-v1.json).

The checker does not reprove the Mazur--Kenku classification, the standard
split-Eichler normalizer theorem, or the rank-three
ternary-spin/marked-K3 period correspondence. It makes no exclusion over
larger number fields and no statement about the anisotropic determinant-378
catalogue row.

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_det378_qq_marking_obstruction.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_det378_qq_marking_obstruction.sage --check
```

## Queue consequence

The global arithmetic-first ledger now has four exact exclusions, one realized
positive control, and `822` unresolved rows. The rootless-MW17 subcatalogue
still has one positive, three exclusions, and `62` unknowns because this split
determinant-378 row was not in that subcatalogue. The next exact-coarse
experiments are determinant `256` with coarse `X_0(2)` and determinant `512`
with coarse `X_0(4)`.
