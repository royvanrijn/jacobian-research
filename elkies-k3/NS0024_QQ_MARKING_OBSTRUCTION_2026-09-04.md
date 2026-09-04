# NS0024 cannot carry the required rational rank-19 marking

Date: 2026-09-04.

<!-- status-consumer: EC-K3-NS0024-QQ-MARKING-OBSTRUCTION b7f0cf002c0411fe -->
<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->

## Theorem

There is no characteristic-zero K3 surface `X/QQ` with

```text
NS(X_Qbar) = NS0024
```

for which all nineteen Neron--Severi divisor classes are defined over `QQ`.
Consequently no rootless fibration on `NS0024` can have a saturated rank-17
Mordell--Weil basis over `QQ(t)`.

This rules out determinant `950` as the arithmetic different-NS milestone.
It does not rule out geometric `NS0024` surfaces or models over larger number
fields.

## Proof

The exact `NS0024` lattice contains the primitive source frame

```text
2E8/MW1,  torsion 1,  height Gram [950].
```

Suppose first that all of `NS(X_Qbar)` is represented by rational divisors.
The fibre, zero, two `E8` configurations, and Mordell--Weil generator of this
primitive `U` embedding are then rational divisor classes.  Riemann--Roch and
the root reflections needed to make the fibre nef take place inside that
rational divisor lattice, so the marked `2II*` Jacobian fibration descends.

Every Jacobian K3 fibration with two `II*` fibres is an Inose fibration
attached over `Qbar` to an unordered pair of elliptic curves `(E1,E2)`.  For
nonisomorphic curves its geometric Mordell--Weil lattice is

```text
Hom(E1,E2)<2>.
```

The primitive height-`950` generator therefore corresponds to a primitive
isogeny of degree `475`.  A primitive isogeny of elliptic curves has cyclic
kernel: a noncyclic kernel would contain a full `E[n]` and make the homomorphism
divisible by `[n]`.  Forgetting the order of `(E1,E2)` and identifying the
isogeny with its dual gives a point of the Fricke quotient

```text
X0+(475) = X0(475)/<w475>.
```

Because the marked fibration is defined over `QQ`, this moduli point is
`QQ`-rational.

Momose's composite-level theorem applies to `N=475`: the prime `19` divides
`475`, satisfies `19>=17`, is not `37`, and is among the primes below `300`
for which the required `J0^-(19)(QQ)` is finite.  Hence every rational point
of `X0+(475)` is a cusp or a CM point.

Our point is neither.  It is noncuspidal because it represents two elliptic
curves and an isogeny.  It is non-CM because an isogenous CM pair has
`Hom(E1,E2)` of rank two over `ZZ`, whereas the displayed Inose frame has
geometric Mordell--Weil rank one.  This contradiction proves the first
claim.

Finally, a rootless arithmetic MW17 fibration over `QQ(t)` with a saturated
height lattice of determinant `950` would itself display the fibre, zero, and
seventeen rational sections generating the full rank-19 Neron--Severi
lattice.  It would therefore give the rational marking just excluded.

## Exact replay and theorem boundary

The compact checker
[`scripts/certify_ns0024_qq_marking_obstruction.py`](scripts/certify_ns0024_qq_marking_obstruction.py)
pins the exact source frame, derives `475=5^2*19`, verifies every numerical
hypothesis of Momose's criterion, and records the rank-one/non-CM gate.  Its
output is
[`../artifacts/generated-results/elkies-k3-ns0024-qq-marking-obstruction-v1.json`](../artifacts/generated-results/elkies-k3-ns0024-qq-marking-obstruction-v1.json).

The checker does not reprove Momose's theorem.  The proof also uses the
standard Inose correspondence and the descent of a fibration from its full
rational divisor marking.  Those are theorem inputs, not bounded-search
inferences.

```bash
python3 elkies-k3/scripts/certify_ns0024_direct_qq_inose_obstruction.py --check
python3 elkies-k3/scripts/certify_ns0024_qq_marking_obstruction.py
python3 elkies-k3/scripts/certify_ns0024_qq_marking_obstruction.py --check
```

## Foundry consequence

The direct source problem on `NS0024` should stop.  Modular reconstruction of
its semistable `A3+A4+A6/MW4` frame remains useful as geometric research, but
it cannot close the arithmetic MW17 milestone over `QQ(t)`.

The next working candidate was `NS0031`, determinant `1184`. Its source key is
the artifact-qualified pair
`(prescribed-root-sources-all-ns-3e8-all-a-v1.json, NS0031-S001)`, not the
shard-local identifier alone. The `A1+2A7/MW2` source has complete-basis pole
profile `[0,1]`, two complete marked pairs in the normalized square-twist
`GF(7)` chart, and a one-parameter formally smooth `ZZ_7` marked branch
through model 157. The same lattice has an exact five-edge physical
degree-two corridor to rootless `NS0031-F017`. The later exact
split-Clifford/`X_0(37)` argument now excludes a full rational `NS0031`
marking as well; see
[`NS0031_QQ_MARKING_OBSTRUCTION_2026-09-04.md`](NS0031_QQ_MARKING_OBSTRUCTION_2026-09-04.md).
The reranked live gate is recorded in
[`DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md`](DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md).

## References

- F. Momose,
  [*Rational points on the modular curves X0+(N)*](https://doi.org/10.2969/jmsj/03920269),
  J. Math. Soc. Japan 39 (1987), Theorem 0.1.
- K. Utsumi,
  [*The Mordell--Weil lattice of an Inose surface arising from isogenous elliptic curves*](https://arxiv.org/abs/2209.02463),
  Proposition 3.1 and Theorems 5.1--5.2.
