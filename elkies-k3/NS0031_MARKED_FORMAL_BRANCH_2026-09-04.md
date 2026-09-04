# NS0031 model-157 marked formal branch

Date: 2026-09-04.

Status: **PROVED LOCAL THEOREM**.

<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->
<!-- status-consumer: EC-K3-NS0031-MARKED-RATIONAL-PARAMETER-SCAN ca678e520745dd3c -->

## Statement

Let the normalized short-Weierstrass model, the two marked source sections,
and their component jets be the 52-variable, 59-equation system attached to
the artifact-qualified source

```text
(artifacts/generated-results/
 elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json,
 NS0031-S001)
```

at the square-twist `GF(7)` model 157 point. In the open chart where the
Weierstrass coefficient `A` is a unit at all three marked fibre supports and
the pole denominator `C` is a unit at `t=1`, the complete marked germ is
formally smooth of relative dimension one over `ZZ_7`.

Equivalently, the model-157 point lies on a compatible one-parameter formal
`ZZ_7` family satisfying all 59 displayed equations. This is stronger than
the previously certified lift through `7^8`.

## Exact dependence of the residual equations

For

```text
D = 4 A^3 + 27 B^2,
H = 2 A X + 3 B C^2,
F = X^3 + A X C^4 + B C^6,
```

direct expansion over `ZZ` gives

```text
8 A^3 F = D C^4 (H - B C^2) - 9 B H^2 C^2 + H^3.
```

At the pole-zero section, the fibre and component equations give order at
least two at `t=0`; the corresponding orders at infinity give residual order
at least four there. Since its global residual has degree at most 12, it is
of the form

```text
t^2 (u0 + u1 t + ... + u6 t^6).
```

Thus the retained residual coefficients 2 through 8 kill all seven quotient
coefficients, and coefficients `0,1,9,10,11,12` follow identically.

At the pole-one section, the fibre and component equations give residual
order at least two at `t=1`. After retained coefficients 0 through 16 vanish,
the residual is

```text
v17 t^17 + v18 t^18.
```

Its value and derivative at `t=1` yield the integral system

```text
v17 + v18 = 0,
17 v17 + 18 v18 = 0,
```

whose determinant is one. Hence coefficients 17 and 18 vanish. All eight
equations omitted from the maximal Jacobian minor therefore follow on the
localized marked scheme, rather than merely at the residue point.

## Formal smoothness

The retained 51 equations have a `51 x 51` Jacobian minor equal to `1 mod 7`
at model 157. There are 52 variables, with `m9` as the complementary formal
coordinate. The formal implicit-function theorem therefore gives a
one-parameter formally smooth `ZZ_7` germ for those 51 equations. The exact
dependence above identifies it with the full 59-equation marked germ.

The certificate also rechecks that `A` is a unit at each marked fibre support
and that `C(1)` is a unit, so the localization used in the proof contains the
certified residue point.

## Proof boundary

This theorem does **not** algebraize the formal branch, rationally
parameterize it over `QQ`, or produce a `QQ`-rational point. It does not prove
that any characteristic-zero member has geometric Picard rank 19, and it does
not yet provide nineteen individually `QQ`-rational Neron--Severi classes.
Those are the remaining source gate for the different-NS arithmetic MW17
milestone.

## Bounded rational-coordinate scan

The free formal coordinate is `m9`, with residue `1 mod 7`. An exact bounded
scan tested every reduced `m9=n/d` satisfying

```text
|n| <= 40,  1 <= d <= 40,  gcd(n,d)=1,  7 does not divide d,
n/d = 1 mod 7.
```

There are 247 such values. For each one, the checker fixed `m9` throughout a
lift to `7^40`, attempted simultaneous rational reconstruction of all 52
coordinates, and was prepared to substitute every reconstruction into all 59
equations and the exact `I2+2I8+6I1` open gates. All 247 rows stopped at
`NO_FULL_RR`; there were no lift errors and no claimed rational points.

This is a bounded negative calculation, not an obstruction. In particular,
rational reconstruction failure does not prove that an individual `7`-adic
coordinate is irrational, and the box does not classify rational points on
the algebraic marked curve.

## Reproduction

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_lattice_foundry_ns0031_marked_formal_smoothness.sage \
  --check
```

The checker is
[`scripts/certify_lattice_foundry_ns0031_marked_formal_smoothness.sage`](scripts/certify_lattice_foundry_ns0031_marked_formal_smoothness.sage).
Its SHA-256 is
`24500201e3bfd04cbf6cb6506c165140112ba84e7ad400efed588ed3e50f5d59`.
The generated certificate is
[`../artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-marked-formal-smoothness-v1.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-marked-formal-smoothness-v1.json),
with SHA-256
`8f6ab911eee02c65427dc8202d99c2300da1ec9eca9cdf35902fde52fd9c943b`.

The bounded rational-coordinate scan is replayed by

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_lattice_foundry_ns0031_rational_parameters.sage \
  --numerator-bound 40 --denominator-bound 40 \
  --lift-precision 40 --workers 8 --check
```

Its checker and generated certificate have SHA-256 respectively
`cbafe412e528e736c2dd0c87e196dd2cca0683593846c7485232e1c4937153e3`
and
`5bbef1cfcd14008dd985c790ed2f9f423f08962a05df8f95473949649d2abdf6`.
