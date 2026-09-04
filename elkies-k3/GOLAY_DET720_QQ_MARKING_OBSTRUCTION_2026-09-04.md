# Determinant 720 cannot carry the required rational rank-19 marking

Date: 2026-09-04.

Status: **PROVED, conditional only on the named standard theorem inputs**.

## Theorem

There is no characteristic-zero K3 surface `X/QQ` with geometric
Néron--Severi lattice equal to the determinant-720 Golay lattice for which all
nineteen Néron--Severi divisor classes are defined over `QQ`. Consequently its
rootless fibration cannot have a saturated rank-17 Mordell--Weil basis over
`QQ(t)`.

The known rational `3A5/MW2` equation is not an exception. Its displayed
determinant-720 sublattice has index six in the full Néron--Severi lattice:
rational `3`-torsion and a rational half-section saturate it to determinant
`20`.

## Split Clifford orders

The relevant ternary transcendental lattice is

```text
T = [ 8  -2   2 ]
    [-2  -4  10 ]
    [ 2  10  -4 ],       det(T)=-720,
```

with

```text
A_T = Z/2 + Z/6 + Z/60.
```

The primitive similarity form is `H=T/2`. In rational isotropic coordinates
it has Gram `U + <90>`. An exact embedding of its even Clifford order in
`M2(QQ)`, in the basis `1,e0e1,e0e2,e1e2`, is

```text
[1 0]  [-2 0]  [ 4  2]  [2 1]
[0 1], [ 0 1], [-5 -3], [5 3].
```

For an integral matrix `[[A,B],[C,D]]`, membership in this order is equivalent
to

```text
5 | C,
B = C/5 (mod 3),
A-D = 2B (mod 3).
```

Its determinant-one group is

```text
K Gamma_0(15) K^-1,       K=[1 0; 5 1].
```

Thus the coarse norm-one curve is `X_0(15)`: index `24`, genus `1`, no
elliptic orbits, and cusp widths `1,3,5,15`. This coarse curve is not yet the
full marking curve.

For the literal lattice `T`, the even Clifford basis maps to

```text
[1 0]  [-4 0]  [  8  4]  [ 4 2]
[0 1], [ 0 2], [-10 -6], [10 6].
```

Its reduced discriminant is `360`, rather than the similarity order's `45`.
The literal order conditions are

```text
2 | B,                         10 | C,
B/2 = C/10 (mod 3),           A-D = 2B-3C/5 (mod 6).
```

For determinant one, these are exactly the primitive-order conditions plus
the identity condition modulo `2`.

## Stable discriminant kernel and exact curve

The spin action of the primitive norm-one group on
`A_T=Z/2+Z/6+Z/60` has order `6`; its element-order histogram is
`1:1, 2:3, 3:2`, so the image is `S3`. Exact comparison on generators shows
that this action is its reduction modulo `2`. Therefore

```text
O^+(T)^* = K (Gamma_0(15) intersection Gamma(2)) K^-1
```

on the marked period component. A negative norm-`4` reflection represents the
determinant-minus-one component of `O^+(T)`; its discriminant action lies
outside the displayed proper `S3` image, so it contributes no further stable
element.

Put `D=diag(2,1)`. Then

```text
D^-1 (Gamma_0(15) intersection Gamma(2)) D = Gamma_0(60).
```

Hence the exact full-marking curve is `X_0(60)` over `QQ`. It has index `144`,
genus `7`, no elliptic orbits, and twelve rational cusps.

The complete Mazur--Kenku classification of rational cyclic isogenies excludes
degree `60`. Thus

```text
X_0(60)(QQ) = {its twelve cusps};
```

in particular it has no rational noncuspidal point. Via the standard
rank-three period/spin correspondence, a full rational determinant-720
Néron--Severi marking would give just such a point. This proves the theorem.

## Useful quotient and saturation check

Forgetting the stable mod-`2` marking gives a degree-six map to `X_0(15)`.
That quotient has four rational noncuspidal points, with `j`-invariants

```text
-5^2/2,
-5^2*241^3/2^3,
-5*29^3/2^5,
 5*211^3/2^15.
```

They are non-CM (each rational `j` is nonintegral), and exact Sage checks give
one rational `3`-isogeny and one rational `5`-isogeny at every point. None
lifts to the stable curve: a lift would supply a rational cyclic
`60`-isogeny.

Independently, the existing equation certificate shows that the rational
`s6=10` `3A5` model has Picard rank `19`, rational `3`-torsion, a rational
half of the displayed height-four section, maximal even-overlattice index
`6`, and full determinant `20`. It is a nonprimitive boundary/control for the
determinant-720 lattice calculation, not a saturated determinant-720 marked
point. No point is sent to an equation agent.

## Exact replay and theorem boundary

The checker
[`scripts/certify_golay_det720_qq_marking_obstruction.sage`](scripts/certify_golay_det720_qq_marking_obstruction.sage)
reconstructs both Clifford embeddings and multiplication tables, derives both
integral matrix descriptions, identifies the congruence groups, computes the
full discriminant action and stable kernel, verifies the `X_0(60)` signature
and rational cusps, checks the four quotient `j`-values, and imports the exact
determinant-20 saturation certificate. Its output is
[`../artifacts/generated-results/elkies-k3-golay-det720-qq-marking-obstruction-v1.json`](../artifacts/generated-results/elkies-k3-golay-det720-qq-marking-obstruction-v1.json).

The checker does not reprove the Mazur--Kenku classification, the complete
determination of `X_0(15)(QQ)`, or the standard marked-K3 ternary-spin period
correspondence. It makes no exclusion over larger number fields and does not
deny complex-geometric determinant-720 K3 surfaces or their exact lattice
corridors.

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_det720_qq_marking_obstruction.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_det720_qq_marking_obstruction.sage --check
```

## References

- B. Mazur and M. A. Kenku's complete classification, summarized as Theorem
  1.1 and Table 1.1 in B. Banwait, F. Najman, and O. Padurariu,
  [*Cyclic isogenies of elliptic curves over fixed quadratic fields*](https://arxiv.org/abs/2206.08891).
- The determinant-20 equation-side saturation proof is in
  [`GOLAY_OCTAD_LATTICE_DESIGN_2026-09-01.md`](GOLAY_OCTAD_LATTICE_DESIGN_2026-09-01.md)
  and its linked exact certificates.
