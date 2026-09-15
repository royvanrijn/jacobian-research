# The star510 point has no rational Fricke lift

The [retained exceptional star510 point](STAR510_ISOTROPIC_PARENT_DESCENT_GATE_2026-09-15.md)
cannot lift rationally to X0(510)/<w510>. Its proposed isotropic Inose-parent
route is closed. The point and corrected quartic remain valid; the failure
is arithmetic descent, not its coordinates or non-CM status.

More generally, **if51 divides N, X0(N)/<wN> has no non-CM rational point**.
This uses an existing complete lower-level classification and needs no
construction of the degree-eight fibre at510.

## Sourced low-degree obstruction

[Ozman–Siksek, Main Theorem and Table8.6](https://arxiv.org/pdf/1806.08192)
classify the quadratic points on X0(51). Their three representatives, together
with conjugates, all have CM: two have j=8000 and discriminant-8, while the
third has j=-2770550784-671956992*sqrt(17) and discriminant-51. There are no
non-CM quadratic points in that complete table.

The rational case is also empty:51 is absent from the
[retained Mazur–Kenku cyclic-isogeny list](NS0024_DIRECT_QQ_INOSE_OBSTRUCTION_2026-09-04.md).
Thus X0(51) has no non-CM point of degree at most2 over Q.
Classification completeness is a published theorem input; no Chabauty or
Jacobian-group computation is replayed here.

## Pull back once, then forget level

Suppose a rational non-CM point existed on X0(N)/<wN>, with51 dividing N.
The degree-two quotient map from X0(N) has a point above it over a number
field of degree at most2. This point is noncuspidal and non-CM.

A cyclic order-N subgroup has a unique order-51 subgroup, preserved by
Galois. Forgetting the rest of the subgroup gives the Q-defined map

```
X0(N) -> X0(51),       (E,C_N) |-> (E,C_51).
```

It preserves the elliptic curve, hence preserves non-CM status, and cannot
increase the residue degree. The resulting degree-at-most-two non-CM point
on X0(51) contradicts the preceding classification.

For N=510 this proves that the degree-eight Fricke-to-star fibre over the
retained exceptional point has no rational point. It need not be explicitly
factored to determine rational solubility. Its full field decomposition is
not computed. This is compatible with the rational point on the full star
quotient: quotienting by all16 Atkin–Lehner elements permits a larger Galois
orbit than quotienting only by w510.

## Exact checks and remaining boundary

The [checker](scripts/verify_star510_fricke_exclusion.py) verifies51|510,
absence of51 from the rational cyclic-isogeny list, and the reported CM
j-values using the Hilbert class polynomials

```
H_-8(J)=J-8000,
H_-51(J)=J²+5541101568*J+6262062317568.
```

It reconstructs the latter polynomial from the conjugate quadratic j-values
and compares both with PARI's class-polynomial computation. The
[result](../artifacts/generated-results/elkies-k3-star510-fricke-exclusion-v1/result.json)
explicitly records the external classification dependency. These finite
checks do not independently establish exhaustion of the modular points.

The prior rootless-frame mismatch at determinant1020 remains true, but
computing a new frame for this specific proposed source is no longer the
next step: the rational Fricke input itself cannot exist. This does not
exclude every isotropic K3 source, other uses of the higher-degree Q-curve
orbit, or the [full correlated-gain objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md).
It does not affect the separately admitted anisotropic determinant1020 K3.

For another isotropic Inose candidate, check proper divisors of its isogeny
level against complete degree-at-most-two modular-point classifications
before constructing its Fricke fibre or searching for frames.

```
.venv/bin/python research/elkies-k3/scripts/verify_star510_fricke_exclusion.py
```
