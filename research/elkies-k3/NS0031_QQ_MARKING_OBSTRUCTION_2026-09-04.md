# NS0031 cannot carry the required rational rank-19 marking

Date: 2026-09-04.

<!-- status-consumer: EC-K3-NS0031-QQ-MARKING-OBSTRUCTION 8e2dc35cdf9b6bc3 -->

## Theorem

There is no characteristic-zero K3 surface `X/QQ` with

```text
NS(X_Qbar) = NS0031
```

for which all nineteen Néron--Severi divisor classes are defined over `QQ`.
Consequently no rootless fibration on `NS0031` can have a saturated rank-17
Mordell--Weil basis over `QQ(t)`.

This closes the determinant-1184 candidate negatively.  It does not conflict
with the formally smooth `ZZ_7` branch through model 157: that branch has
`QQ_7` points but no rational point satisfying the full marking can exist.

## The exact modular curve

The rank-three transcendental lattice attached to `NS0031` is

```text
T = [ 0   0   4 ]
    [ 0  74   1 ]
    [ 4   1  -2 ],       det(T) = -1184.
```

It is rationally isotropic, with primitive isotropic vector `(1,0,0)` of
divisibility `4`.  Its even Clifford algebra is split.  In the basis

```text
1, e0*e1, e0*e2, e1*e2
```

an exact embedding into `M2(QQ)`, followed by conjugation by `diag(4,1)`,
gives the integral order with basis

```text
[1 0]  [0 4]  [4 0]  [ 1 1]
[0 1], [0 0], [0 0], [37 0].
```

Thus a matrix `[A B; C D]` belongs to this order exactly when

```text
37 | C,
B     = C/37 (mod 4),
A - D = C/37 (mod 4).
```

For determinant one, the mod-4 reduction is the norm-one subgroup of the
unramified non-split Cartan and the mod-37 reduction is upper triangular.
The projective norm-one modular curve is therefore

```text
X_ns(4) x_{X(1)} X_0(37),
```

with congruence group `Gamma_ns(4) intersection Gamma_0(37)`.  Direct coset
enumeration gives

```text
index       304
elliptic-2  0
elliptic-3  4
cusps       4, with widths 4,4,148,148
genus       23.
```

The full rational `NS0031` marking fixes the discriminant gluing.  By the
standard rank-three period/spin description, its moduli curve is a cover of
this norm-one curve.  Hence a rational marked K3 would give a noncuspidal
`QQ`-point on the displayed fibre product and, after forgetting the level-4
structure, a noncuspidal rational point of `X_0(37)`.

## The rational-point obstruction

Vélu determined `X_0(37)(QQ)`.  Apart from its cusps, its two rational points
have elliptic-curve `j`-invariants

```text
j1 = -7*11^3,
j2 = -7*137^3*2083^3.
```

Neither point lifts to `X_ns(4)`.  For each `j`, an exact rational minimal
model has a rational 37-isogeny, good reduction at `19`, and

```text
a_19 = -6,       (a_19 mod 4, 19 mod 4) = (2,3).
```

For the full unramified non-split Cartan in `GL(2,Z/4Z)`, the possible
`(trace,determinant)` pairs are

```text
(0,3), (1,1), (1,3), (2,1), (3,1), (3,3).
```

The pair `(2,3)` is absent.  Frobenius at `19` therefore excludes a mod-4
Galois image inside that Cartan.  This is unchanged by quadratic twisting:
the trace changes sign, while `6` and `-6` are both `2` modulo `4`.  Hence
neither noncuspidal point of `X_0(37)(QQ)` lifts, and the fibre product has no
noncuspidal rational point.

This contradicts the point forced by a full rational `NS0031` marking and
proves the theorem.

## Foundry consequence

The model-157 formal branch, the rational-coordinate scan, and the exact
five-edge route to `NS0031-F017` remain valid local and geometric evidence.
They cannot be upgraded to the requested arithmetic source over `QQ`, so the
different-NS foundry objective must move again.  A replacement candidate must
pass the rational-marking arithmetic gate before another equation-facing
coefficient campaign begins.

## Exact replay and theorem boundary

The checker
[`scripts/certify_ns0031_qq_marking_obstruction.sage`](scripts/certify_ns0031_qq_marking_obstruction.sage)
reconstructs the Clifford embedding and integral order, identifies the two
local congruence conditions, computes the full signature, verifies the two
37-isogenies, and performs both Frobenius exclusions.  Its output is
[`../artifacts/generated-results/elkies-k3-ns0031-qq-marking-obstruction-v1.json`](../artifacts/generated-results/elkies-k3-ns0031-qq-marking-obstruction-v1.json).

The checker does not reprove Vélu's global determination of
`X_0(37)(QQ)` or the general marked-K3 period/Clifford correspondence.  Those
are theorem inputs.  No claim is made about geometric `NS0031` surfaces,
models over larger number fields, or rational models with a proper
Galois-invariant sublattice.

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_ns0031_qq_marking_obstruction.sage
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_ns0031_qq_marking_obstruction.sage --check
```

## Reference

- J. Vélu,
  [*Les points rationnels de X0(37)*](https://doi.org/10.24033/msmf.145),
  Bull. Soc. Math. France, Mémoire 37 (1974), 169--179.
