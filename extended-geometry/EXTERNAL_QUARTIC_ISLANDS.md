# Zhuang's quartic islands: weighted-seed identification and boundary audit

This note integrates three explicit maps published by Juntang Zhuang as
`F4a`, `F4b`, and `F4c`.  Zhuang's repository and accompanying PDF call them
Islands A, B, and C, give their expanded coordinate polynomials, check their
constant Jacobians, and record rational collisions.  The external source is
pinned at commit
[`1ff68e870f66afec8c6611f910fcc8f5522fdbce`](https://github.com/jzkay12/jacobian_conjecture/commit/1ff68e870f66afec8c6611f910fcc8f5522fdbce).

The labels, expanded formulas, and collision certificates are credited to
Zhuang's public compilation.  The PDF says that the note was inspired by
Levent Alpöge's public post.  Neither that statement nor the present audit
settles discovery priority for the individual quartic maps.  The resolvent,
normalization-boundary, nilpotency, monodromy, and cancellation-comparison
calculations below are repository analyses of the pinned examples.  No
license file was present at the pinned commit; accordingly this repository
links to the upstream materials but does not copy their checker or prose.

## 1. Compact reconstruction of the three maps

Write the external target coordinates as `(p,q,r)`.  On the common source put

\[
 A_0=1+xy,
 \qquad D=1-x(\lambda y+xz),
 \qquad T={A_0\over x},
 \qquad W=rT=A_0D.
\]

The three values of `lambda` are

\[
 \lambda_a={5\over4},\qquad
 \lambda_b=2,\qquad
 \lambda_c={13\over10}.
\]

Rewriting the expanded maps in `(A_0,T,r)` and eliminating `A_0` gives the
following exact inverse equations.

| map | quartic inverse equation | primitive `H(W)` |
|---|---|---|
| `F4a` | `W^4-W^3+qrW-pr^2=0` | `W^3(W-1)` |
| `F4b` | `W^4-4W^3+3W^2-qrW+2pr^2=0` | `W^2(W-1)(W-3)` |
| `F4c` | `-2W^4+W^3+W^2-qrW+3pr^2=0` | `-W^2(W-1)(2W+1)` |

For `F4a`, the standard weighted target coordinates are
`(A,B,C)=(p,-q,r)`.  For `F4b` and `F4c`, they are `(p,q,r)`.  Thus every
row is an instance of

\[
 E_H(W)=H(W)-BCW+cAC^2,
 \qquad c=-H'(1).
\]

By the [tangent-map core theorem](../verified/TANGENT_MAP_CORE.md), this same
pencil simultaneously supplies the Jacobian factor, discriminant
normalization, reconstruction pole, and Hessian Fitting divisor.  In
particular it proves generic degree four.  It also identifies the apparent islands:

- `F4a` is the canonical triple-zero seed `W^3(W-1)`;
- `F4b` is a split double-zero seed with extra root `rho=3`;
- `F4c` is a split double-zero seed with extra root `rho=-1/2`.

They are weighted-seed strata, not cancellation-parameter roots.

## 2. Normalized source and target boundaries

Let `Q_H` be the reduced saturated discriminant of `E_H`.  The two target
boundary images are

\[
 Z_\Delta=V(Q_H),\qquad Z_C=V(C).
\]

The normalization of the discriminant incidence is the closure of

\[
 BC=H'(u),
 \qquad
 cAC^2=uH'(u)-H(u).
\]

Over `Z_Delta` there is one normalized source-boundary prime `D_Delta` with
`(e,f)=(2,1)` and sheet loss two in all three cases.  The second normalized
source boundary depends on the root profile:

| map | source prime over `Z_C` | `(e,f)` | sheet loss |
|---|---|---|---:|
| `F4a` | zero-cluster prime `D_0` | `(2,1)` | 2 |
| `F4b` | extra-root prime `D_3` | `(1,1)` | 1 |
| `F4c` | extra-root prime `D_(-1/2)` | `(1,1)` | 1 |

The discriminant ramification index is therefore `e_Delta=2` for every map.
The raw pulled-back discriminants have exact `C`-adic orders three, two, and
two, respectively; those coordinate orders are not the discriminant
ramification index.

## 3. Thick boundary intersections

After exact saturation, restriction of `Q_H` to `C=0` gives

| map | scheme-theoretic target intersection | nilpotency index `mu` |
|---|---|---:|
| `F4a` | `k[A,B]/(B^3)` | 3 |
| `F4b` | `k[A,B]/(B^2-24A) ~= k[B]` | 1 |
| `F4c` | `k[A,B]/(B^2-12A) ~= k[B]` | 1 |

Here `mu=1` means that the intersection is reduced.  Hence the canonical
boundary object separates Island A from Islands B and C.  This numerical
boundary signature does not separate B from C.  The stronger
[decorated-normalization invariant](DECORATED_NORMALIZATION_INVARIANT.md)
does: their unmarked cusp and node-pair configurations are affinely
equivalent, but adding the canonical zero-cluster boundary mark gives the
distinct exact invariants `2` and `3/19`.  Thus Islands B and C are not
polynomially left--right equivalent, even after stabilization.

## 4. Monodromy and cancellation exclusion

The universal pencil theorem applies to every polynomial `H` in the table.
Consequently all three geometric and arithmetic generic monodromy groups are
the natural `S_4`.

A degree-four cancellation normal form must satisfy

\[
 r(m+1)+1=4,
\]

so its type is uniquely `(m,r)=(2,1)`.  That type has

\[
 e_\Delta=2,
 \qquad e_{P=0}=1,
 \qquad \mu=mr(m+1)=6.
\]

Island A is excluded because its second boundary has index two and its thick
intersection has index three.  Islands B and C have the same numerical
ramification profile as cancellation type `(2,1)`, but their intersections
are reduced rather than six-fold nilpotent.  Therefore none of the three maps
is polynomially left--right equivalent, even after stabilization, to any
`(m,r,h)` cancellation normal form.

The resulting conclusion is more precise than calling the examples three
unrelated islands: all three lie outside the cancellation skeleton, but all
three lie inside the already classified weighted-seed construction.  They do
not by themselves introduce a new construction direction.

## 5. Exact reproduction

Run

```bash
.venv/bin/python scripts/verify_external_quartic_islands.py
```

The script independently reconstructs the maps from compact formulas, checks
their ordered rational-coefficient fingerprints against the pinned expanded
formulas, checks polynomiality, Jacobian determinants and Zhuang's reported collisions, derives
the three resolvents by elimination, computes the saturated discriminant
traces and nilpotency indices, and checks the cancellation obstruction.  It
uses the written universal-pencil theorem for the `S_4` conclusion rather
than presenting a finite calculation as a monodromy proof.
