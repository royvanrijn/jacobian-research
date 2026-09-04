# Frobenius bounds for the two rank-55 singleton twists

## Result

For the two singleton characters in the shortlist-rank-55 `V4` base, exact
toric Frobenius calculations give

```text
1 <= rank E^(q_0fda0)(QQbar(u)) <= 2,
1 <= rank E^(q_1037d)(QQbar(u)) <= 2.               (1)
```

The lower bounds are the existing exact lifted sections.  The upper bounds
are unconditional Picard-number bounds and do not use the Tate conjecture.
They are the sharpest bounds obtainable here merely by counting Tate roots at
one good prime.  In particular, (1) does **not** prove the desired exact
character decomposition `17+1+1+0`; the two singleton rank-one assertions
remain `UNKNOWN`.

The compact certificates are

- [`elkies-k3-r17-singleton-alternate-orbit-0fda0-p131-toric-frobenius-v1.json`](../artifacts/generated-results/elkies-k3-r17-singleton-alternate-orbit-0fda0-p131-toric-frobenius-v1.json);
- [`elkies-k3-r17-singleton-alternate-orbit-1037d-p157-toric-frobenius-v1.json`](../artifacts/generated-results/elkies-k3-r17-singleton-alternate-orbit-1037d-p157-toric-frobenius-v1.json).

## Cohomological calculation

For either quadratic `d=q_i`, the regular model

```text
d(u)y^2=x^3+A(u)x+B(u)
```

has two geometric `I0*` fibres, the original twenty-four `I1` fibres, and
arithmetic genus three.  Its trivial lattice is `U+2D4`, of rank ten, and the
nontrivial elliptic cohomology has degree twenty-four.

The Newton polytope has vertices

```text
(0,0,0), (12,0,0), (0,3,0), (0,0,2), (2,0,2)
```

and primitive Hodge vector `[2,24,2]`.  If `D={d=0}` and
`Z={d=0,x^3+A*x+B=0}`, the independent verifier reconstructs

```text
P_boundary = P_Z/P_D,                         degree 4,
P_triv = (T-p)^2 P_D P_Z,                     degree 10,
P_ambient = (T-p)^2 P_D^2,                    degree 6,
P_H2 = P_ambient P_toric = P_triv P_E,        degree 34,
P_E = P_toric/P_boundary.                     degree 24.       (2)
```

It checks integrality, the exact reciprocal functional equation, the Weil
circle by certified real-root isolation, independently computed fibrewise
power sums for `n=1,2`, and every cyclotomic polynomial allowed in degree 24.

For `q_0fda0` at `p=131`, the normalized `P_E` has exactly the Tate factor
`(Z-1)^2`.  For `q_1037d`, the reductions at `131`, `137`, and `151` have
Tate degrees ten, six, and four; the next good prime `157` has exactly
`(Z-1)^2`.  Thus in the two useful reductions

```text
rho(Fpbar) <= 10+2 = 12,
rank MW(QQbar(u)) <= 12-2-2*rank(D4) = 2.       (3)
```

Together with the known section, (3) proves (1).

## Exact remaining gate

The multiplicity in `(Z-1)^2` matters: counting only the distinct factor
`Z-1` would incorrectly report upper bound one.  A one-prime Tate-root count
therefore cannot close the odd geometric rank.  An exact rank-one theorem
needs additional arithmetic input, for example:

1. a complete two-descent giving Mordell--Weil rank at most one; or
2. two full-rank reduction lattices with incompatible discriminant
   squareclasses, established without assuming that every Tate root is
   algebraic.

Until one of those gates is crossed, the certified character information for
the connected rank-55 `V4` base is

```text
17 + [1,2] + [1,2] + 0,
```

so its total geometric generic rank lies in `[19,21]`.

## Replay

```bash
sage -python \
  elkies-k3/scripts/audit_r17_singleton_twist_finite_field_bounds.sage
elkies-k3/scripts/run_r17_singleton_toric_frobenius.sh \
  alternate-orbit-0fda0 131
elkies-k3/scripts/run_r17_singleton_toric_frobenius.sh \
  alternate-orbit-1037d 157
```

The runner pins the same open-source controlled-reduction commit as the
product certificate and records the build invocation, executable, raw input
and output, parser, verifier, and runner hashes.

