# ICARM curve 273: unconditional rank at least 30

Status: **exact unconditional rank lower bound and independent replay**.
No unconditional exact-rank statement is claimed.

## Public record

The current 2026 public record page attributes the curve to Levent Alpöge and
Ava Howell.  It is the same curve that first appeared on the ICARM leaderboard
under the submitter name `ranksunbounded`:

```text
y^2 + x*y = x^3 + A*x + B

A = -201769035260418549083594900060734240952308696994802735114305555
B = 1151107939141058565733479426024323225135665982951300586808823640527729578307228357301072889377.
```

The source is
[Dujella's maintained rank-30 record page](https://web.math.pmf.unizg.hr/~duje/tors/rk30.html).
The thirty displayed points there agree exactly with
[`icarm_curve273.py`](../cas/icarm_curve273.py).

This attribution was rechecked on 2026-08-23.  At that time the rank-history
page and rank-30 page had SHA-256 hashes
`239ed273283ac47b6cc3dcbd08120fdb2c01bcb8494d14d2447862244e841fcc`
and `0828a3f55f3469119797b923a90d15bf6a47367769b6111ff79c3051ef1b28e6`,
respectively.  These source hashes certify the inspected snapshots, not the
independence claim; the exact local certificate below supplies that proof.

## Exact lower-bound certificate

The primary checker is
[`verify_icarm_curve273_rank30.py`](../cas/verify_icarm_curve273_rank30.py).
It performs the following exact replay.

1. All 30 displayed rational points satisfy the integral Weierstrass equation.
2. The rational change of variables

   ```text
   X = 36*x + 3,
   Y = 108*(2*y+x)
   ```

   transports the curve and points to an integral short Weierstrass model.
3. The short 2-division cubic has no root modulo 23. Hence it is irreducible
   over `Q`, so the curve has no rational 2-torsion.
4. Exact exhaustive finite-group calculations at

   ```text
   11,17,19,29,43,61,73,79,89,101,103,113,127,
   131,137,151,173,193,197,211,223,227,229,241,257
   ```

   give 31 rows in the product of the quotients
   `E(F_p)/2E(F_p)`. Their binary matrix has 30 columns and rank 30.

If an integral relation among the points existed, its reduction in every
displayed quotient would force every coefficient to be even. Dividing the
relation by two introduces a rational 2-torsion point; since rational
2-torsion is trivial, this process repeats indefinitely. Therefore every
coefficient is zero and the points are independent over `Z`.

Thus

```text
rank E(Q) >= 30
```

unconditionally.

The compact pinned manifest is
[`icarm_curve273_rank30_v1.json`](../../artifacts/generated-results/elliptic-curves/icarm_curve273_rank30_v1.json),
with SHA-256

```text
e2a7a322fbd4703af4239f497749a69a68f9d5149aa8a1f696b39ab3941a3284
```

Two consecutive generations produced the same whole-file hash.

For cleanup-only provenance checking, including the exact artifact and every
directly imported proof helper, run:

```text
python3 elliptic-curves/cas/audit_icarm_rank_lower_bound_artifacts.py
```

That command does not enumerate a finite group, recalculate a matrix rank, or
invoke PARI or Sage. The commands below remain the mathematical replays.

## Model diagnostics

PARI/GP 2.17.4 independently reports that the displayed integral model is
already global minimal, with trivial rational torsion, root number `+1`, and
exact conductor

```text
2381958488309327728488641214562148525681586925734398576894288277390640894675828305511266759053188773997283310012808983216600083938367948090232306090.
```

The 120-digit numerical height determinant is nonzero, but it is stored only
as a diagnostic and is not used in the rank proof.

## Independent implementation

[`verify_icarm_curve273_rank30_sage.py`](../scripts/verify_icarm_curve273_rank30_sage.py)
replays the result through Sage's independent elliptic-curve interfaces. For
each certificate prime it computes the finite group invariant factors and
discrete logarithms, extracts the even invariant-factor coordinates, and
again obtains a `31 x 30` binary matrix of rank 30. It separately checks the
global minimal model, conductor, root number, torsion, and the modulo-23
2-division polynomial.

## Claim boundary and continuation

What is proved:

```text
rank E(Q) >= 30.
```

What is not proved unconditionally:

```text
rank E(Q) = 30,
rank E(Q) >= 31.
```

The completed direct 31st-point search covered 1,933 of 1,935 declared charts;
two charts timed out and no nonbasis image was found. A later short-height
lattice pass completed 3,000 charts with no new image. These are bounded
negative computations, not a rank upper bound. The unfinished alternate-cover
and custom 2-Selmer work remains under `artifacts/local/elliptic-curves/` and
must not be promoted until its scripts, parameters, outputs, and completeness
claims are pinned.

The principal research continuation is a certified residual 2-Selmer
calculation in the cubic 2-division field. Exact smooth principal-ideal
relations may either produce explicit residual 2-covers to search for a 31st
point or give the missing upper bound needed for an exact-rank result.

The separate source and family investigation is maintained in
[`ICARM_CURVE273_CONSTRUCTION_INVESTIGATION.md`](ICARM_CURVE273_CONSTRUCTION_INVESTIGATION.md).
