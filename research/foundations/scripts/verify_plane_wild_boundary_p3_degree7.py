#!/usr/bin/env python3
"""Exact normalization scan of the first balanced odd-characteristic row.

In characteristic three the retained-sheet support theorem forces a monic
degree-four polynomial to have the form ``A=T^4+bT+a0``.  Over ``F_3`` there
are six choices with ``a0 != 0``.  Singular verifies their normalizations,
conductors, smoothness, relative differents, fierce boundaries, and exact
point counts.  The two rows surviving over ``F_3`` are then counted over
``F_9`` and ``F_27``; both fail the affine-plane count over ``F_27``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
from pathlib import Path
from typing import Any

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "plane_wild_boundary_p3_degree7_scan.json"
)

PRIME = 3
EXPECTED_COVER_POINTS = {
    (1, 0): 9,
    (1, 1): 12,
    (1, 2): 12,
    (2, 0): 15,
    (2, 1): 9,
    (2, 2): 9,
}
EXPECTED_BOUNDARY_POINTS = 3
EXTENSION_FIELD_CONTROLS = {
    2: {
        "minpoly": "a2+a+2",
        "order": 9,
        "rows": {(1, 1): (90, 9), (1, 2): (90, 9)},
    },
    3: {
        "minpoly": "a3+2a+1",
        "order": 27,
        "rows": {(1, 1): (837, 27), (1, 2): (837, 27)},
    },
}


def singular_program(a0: int, linear: int, cover_points: int) -> str:
    """Return an exact Singular certificate for one retained polynomial."""

    return rf'''
LIB "normal.lib";

proc assertReductionZero(poly f, ideal G, string label)
{{
  if (reduce(f,std(G)) != 0)
  {{
    "FAIL: "+label;
    exit(1);
  }}
}}

proc assertIdealEqual(ideal A, ideal B, string label)
{{
  int i;
  for (i=1; i<=size(A); i++)
  {{
    assertReductionZero(A[i],B,label);
  }}
  for (i=1; i<=size(B); i++)
  {{
    assertReductionZero(B[i],A,label);
  }}
}}

ring r=3,(P,Q,T),dp;
poly A=T4+{linear}*T+{a0};
ideal primitiveOrder=A*(T3-P4)+P2*Q*T;
list N=normal(primitiveOrder,"isPrim","withGens");
intvec normalCheck=norTest(primitiveOrder,N);
if (normalCheck[1] != 1 || normalCheck[2] != 1 || normalCheck[3] != 1)
{{
  "FAIL: normalization";
  exit(1);
}}
assertIdealEqual(
  normalConductor(primitiveOrder),
  ideal(P2,P*T,T2),
  "primitive conductor"
);

def R=N[1][1]; setring R;
poly U=var(1);
poly V=var(2);
ideal C=norid;
ideal fierce=C+ideal(Q,T3-P4);
ideal explicitFierce=C+ideal(Q,P-U*V,T-U^2,U-V^2);
assertIdealEqual(radical(fierce),explicitFierce,"fierce affine line");

int nv=nvars(basering);
int nr=size(C);
int nc=nv-2;
matrix relativeJacobian[nr][nc];
int rowIndex;
int columnIndex;
for (rowIndex=1;rowIndex<=nr;rowIndex++)
{{
  for (columnIndex=1;columnIndex<nc;columnIndex++)
  {{
    relativeJacobian[rowIndex,columnIndex]=
      diff(C[rowIndex],var(columnIndex));
  }}
  relativeJacobian[rowIndex,nc]=diff(C[rowIndex],var(nv));
}}
ideal relativeDifferent=C+minor(relativeJacobian,nc);
ideal residualDifferent=sat(relativeDifferent,fierce);
assertIdealEqual(residualDifferent,C+ideal(1),"different off fierce boundary");
assertIdealEqual(
  radical(relativeDifferent),explicitFierce,"different radical"
);

ideal absoluteSingularLocus=C+minor(jacob(C),nv-2);
assertIdealEqual(absoluteSingularLocus,C+ideal(1),"absolute smoothness");

ideal fieldEquations;
int variableIndex;
for (variableIndex=1;variableIndex<=nv;variableIndex++)
{{
  fieldEquations[variableIndex]=var(variableIndex)^3-var(variableIndex);
}}
int coverPoints=vdim(std(C+fieldEquations));
int boundaryPoints=vdim(std(fierce+fieldEquations));
if (coverPoints != {cover_points})
{{
  "FAIL: cover point count";
  exit(1);
}}
if (boundaryPoints != {EXPECTED_BOUNDARY_POINTS})
{{
  "FAIL: boundary point count";
  exit(1);
}}

"PASS_P3_DEGREE7_{a0}_{linear}";
'''


def extension_point_count_program(
    a0: int,
    linear: int,
    extension_degree: int,
    field_order: int,
    minpoly: str,
    cover_points: int,
    boundary_points: int,
) -> str:
    """Return an exact normalization and point-count control over ``F_(3^n)``."""

    return rf'''
LIB "normal.lib";
ring r=(3,a),(P,Q,T),dp;
minpoly={minpoly};
poly A=T4+{linear}*T+{a0};
ideal primitiveOrder=A*(T3-P4)+P2*Q*T;
list N=normal(primitiveOrder,"isPrim","withGens");
intvec normalCheck=norTest(primitiveOrder,N);
if (normalCheck[1] != 1 || normalCheck[2] != 1 || normalCheck[3] != 1)
{{
  "FAIL: extension normalization";
  exit(1);
}}

def R=N[1][1]; setring R;
ideal C=norid;
ideal fierce=C+ideal(Q,T3-P4);
int nv=nvars(basering);
ideal fieldEquations;
int variableIndex;
for (variableIndex=1;variableIndex<=nv;variableIndex++)
{{
  fieldEquations[variableIndex]=
    var(variableIndex)^{field_order}-var(variableIndex);
}}
int coverPoints=vdim(std(C+fieldEquations));
int boundaryPoints=vdim(std(fierce+fieldEquations));
if (coverPoints != {cover_points})
{{
  "FAIL: extension cover point count";
  exit(1);
}}
if (boundaryPoints != {boundary_points})
{{
  "FAIL: extension boundary point count";
  exit(1);
}}

"PASS_P3_DEGREE7_EXT{extension_degree}_{a0}_{linear}";
'''


def verify_symbolic_row(a0: int, linear: int) -> None:
    """Check squarefreeness and the balanced support identity in SymPy."""

    p, q, t = sp.symbols("P Q T")
    retained = t**4 + linear * t + a0
    assert sp.gcd(
        sp.Poly(retained, t, modulus=PRIME),
        sp.Poly(sp.diff(retained, t), t, modulus=PRIME),
    ).degree() == 0
    support_factor = sp.Poly(
        retained - t * sp.diff(retained, t),
        t,
        modulus=PRIME,
    )
    assert support_factor.degree() == 0
    assert int(support_factor.LC()) % PRIME == a0

    boundary_factor = t**3 - p**4
    order = retained * boundary_factor + p**2 * q * t
    identity = (
        retained * sp.diff(order, t)
        - sp.diff(retained, t) * order
        - p**2 * q * (retained - t * sp.diff(retained, t))
    )
    assert sp.Poly(
        sp.expand(identity),
        p,
        q,
        t,
        modulus=PRIME,
    ).is_zero


def retained_root_count(a0: int, linear: int, field_order: int) -> int:
    """Count roots of the retained polynomial in ``F_q`` exactly."""

    t = sp.symbols("T")
    retained = sp.Poly(t**4 + linear * t + a0, t, modulus=PRIME)
    field_polynomial = sp.Poly(t**field_order - t, t, modulus=PRIME)
    return sp.gcd(retained, field_polynomial).degree()


def expected_counts_from_roots(field_order: int, root_count: int) -> dict[str, int]:
    """Apply the proved normalization stratification count."""

    return {
        "cover": field_order**2 + root_count * field_order,
        "boundary": field_order,
        "open": field_order**2 + (root_count - 1) * field_order,
    }


def retained_factor_degrees(a0: int, linear: int) -> list[int]:
    """Recover the squarefree quartic's factor degrees without factor ordering."""

    linear_degree = retained_root_count(a0, linear, PRIME)
    degree_at_most_two = retained_root_count(a0, linear, PRIME**2)
    quadratic_count = (degree_at_most_two - linear_degree) // 2
    degrees = [1] * linear_degree + [2] * quadratic_count
    remaining = 4 - sum(degrees)
    if remaining:
        degrees.append(remaining)
    assert sum(degrees) == 4
    return degrees


def compile_scan() -> dict[str, Any]:
    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required for the degree-seven scan")

    rows = []
    for (a0, linear), cover_points in sorted(EXPECTED_COVER_POINTS.items()):
        verify_symbolic_row(a0, linear)
        marker = f"PASS_P3_DEGREE7_{a0}_{linear}"
        completed = subprocess.run(
            [singular, "-q"],
            input=singular_program(a0, linear, cover_points),
            text=True,
            capture_output=True,
            check=False,
            timeout=55,
        )
        if completed.returncode != 0 or marker not in completed.stdout:
            raise AssertionError(
                f"Singular failed for a0={a0}, b={linear}\n"
                f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
            )
        open_points = cover_points - EXPECTED_BOUNDARY_POINTS
        prime_field_roots = retained_root_count(a0, linear, PRIME)
        predicted = expected_counts_from_roots(PRIME, prime_field_roots)
        assert predicted == {
            "cover": cover_points,
            "boundary": EXPECTED_BOUNDARY_POINTS,
            "open": open_points,
        }
        prime_field_survivor = open_points == PRIME**2
        rows.append(
            {
                "A": f"T^4+{linear}*T+{a0}",
                "a0": a0,
                "b": linear,
                "conductor": "(P,T)^2",
                "cover_points_F3": cover_points,
                "decision": (
                    "requires_extension_count"
                    if prime_field_survivor
                    else "geometrically_obstructed"
                ),
                "decisive_field": None if prime_field_survivor else "F_3",
                "extension_point_counts": [],
                "factor_degrees_over_F3": retained_factor_degrees(a0, linear),
                "fierce_boundary": "A1",
                "fierce_boundary_points_F3": EXPECTED_BOUNDARY_POINTS,
                "normalization": "smooth",
                "open_points_F3": open_points,
                "prime_field_decision": (
                    "point_count_survivor"
                    if prime_field_survivor
                    else "obstructed_over_F3"
                ),
                "retained_roots_F3": prime_field_roots,
                "relative_different_support": "fierce_boundary_only",
            }
        )

    survivors = [
        row for row in rows if row["decision"] == "requires_extension_count"
    ]
    assert [(row["a0"], row["b"]) for row in survivors] == [(1, 1), (1, 2)]

    for extension_degree, control in EXTENSION_FIELD_CONTROLS.items():
        field_order = control["order"]
        minpoly = control["minpoly"]
        for row in survivors:
            row_id = (row["a0"], row["b"])
            cover_points, boundary_points = control["rows"][row_id]
            marker = (
                f"PASS_P3_DEGREE7_EXT{extension_degree}_"
                f"{row['a0']}_{row['b']}"
            )
            completed = subprocess.run(
                [singular, "-q"],
                input=extension_point_count_program(
                    row["a0"],
                    row["b"],
                    extension_degree,
                    field_order,
                    minpoly,
                    cover_points,
                    boundary_points,
                ),
                text=True,
                capture_output=True,
                check=False,
                timeout=55,
            )
            if completed.returncode != 0 or marker not in completed.stdout:
                raise AssertionError(
                    "Singular extension count failed for "
                    f"a0={row['a0']}, b={row['b']}, q={field_order}\n"
                    f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
                )
            root_count = retained_root_count(
                row["a0"], row["b"], field_order
            )
            predicted = expected_counts_from_roots(field_order, root_count)
            open_points = cover_points - boundary_points
            assert predicted == {
                "cover": cover_points,
                "boundary": boundary_points,
                "open": open_points,
            }
            row["extension_point_counts"].append(
                {
                    "boundary": boundary_points,
                    "cover": cover_points,
                    "field": f"F_{field_order}",
                    "open": open_points,
                    "retained_roots": root_count,
                }
            )

    for row in survivors:
        final_control = row["extension_point_counts"][-1]
        assert final_control["field"] == "F_27"
        assert final_control["open"] != 27**2
        row["decision"] = "geometrically_obstructed"
        row["decisive_field"] = "F_27"

    return {
        "characteristic": PRIME,
        "cover_degree": 7,
        "format": "plane-wild-boundary-p3-degree7-scan-v2",
        "rows": rows,
        "scope": (
            "Complete F_3 coefficient scan inside the monic balanced "
            "support-admissible degree-four retained family, followed by "
            "exact F_9 and F_27 counts for the two prime-field survivors. "
            "Both survivors fail over F_27. No affine-plane source or "
            "constant-Jacobian polynomial map is produced."
        ),
        "summary": {
            "coefficient_rows": len(rows),
            "geometrically_obstructed": len(rows),
            "geometric_survivors": 0,
            "normalization_smooth": len(rows),
            "point_count_obstructed_over_F3": len(rows) - len(survivors),
            "point_count_survivor_ids_over_F3": [
                [row["a0"], row["b"]] for row in survivors
            ],
            "point_count_survivors_over_F3": len(survivors),
            "point_count_survivors_over_F27": 0,
            "relative_different_gate_passed": len(rows),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--emit-json", action="store_true")
    args = parser.parse_args()

    artifact = compile_scan()
    serialized = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    if args.emit_json:
        print(serialized, end="")
        return
    if args.write:
        OUTPUT.write_text(serialized)
    else:
        assert OUTPUT.exists(), (
            f"missing {OUTPUT.relative_to(ROOT)}; regenerate with --write"
        )
        assert OUTPUT.read_text() == serialized, (
            f"{OUTPUT.relative_to(ROOT)} is stale; regenerate with --write"
        )
    digest = hashlib.sha256(serialized.encode()).hexdigest()
    summary = artifact["summary"]
    print(
        "PASS p=3 degree-seven retained-polynomial scan: "
        f"rows={summary['coefficient_rows']}, "
        "prime_field_survivors="
        f"{summary['point_count_survivors_over_F3']}, "
        f"geometric_survivors={summary['geometric_survivors']}"
    )
    print(f"SHA256: {digest}")


if __name__ == "__main__":
    main()
