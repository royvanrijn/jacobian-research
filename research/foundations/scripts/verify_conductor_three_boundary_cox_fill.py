#!/usr/bin/env python3
"""Verify the minimal node/cusp three-boundary Cox-fill obstruction.

The construction replaces the separated conductor localization by

    x*y*z = c(t),

where c=t(t-1) for a node and c=t^2 for a cusp.  The three factors remove
the old unit-group objection and make every descended incidence coordinate
polynomial.  The exact obstruction moves to the conductor/canonical ledger
and affine-space recognition.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
SINGULAR_SCRIPT = (
    ROOT / "scripts" / "verify_conductor_three_boundary_cox_fill.sing"
)
ARTIFACT = (
    ROOT
    / "artifacts"
    / "generated-results"
    / "conductor_three_boundary_cox_fill.json"
)

x, y, z, t, p, q, T, L = sp.symbols("x y z t p q T L")


def run_singular() -> dict[str, object]:
    """Replay elimination, smoothness, and singular-locus ideals."""

    singular = shutil.which("Singular")
    if singular is None:
        raise SystemExit("Singular is required on PATH")
    completed = subprocess.run(
        [singular, "-q", str(SINGULAR_SCRIPT)],
        check=True,
        capture_output=True,
        text=True,
    )
    passes = tuple(
        line
        for line in completed.stdout.splitlines()
        if line.startswith("SINGULAR_PASS")
    )
    assert len(passes) == 3, completed.stdout
    version = subprocess.run(
        [singular, "--version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()[0]
    return {
        "version": version,
        "script": str(SINGULAR_SCRIPT.relative_to(ROOT)),
        "script_sha256": hashlib.sha256(
            SINGULAR_SCRIPT.read_bytes()
        ).hexdigest(),
        "certificates": passes,
    }


def conductor_data() -> dict[str, object]:
    """Check the two conductor algebras and their marked normalizations."""

    node_p = t * (t - 1)
    node_q = t * node_p
    cusp_p = t**2
    cusp_q = t**3

    node_relation = q**2 - p * q - p**3
    cusp_relation = q**2 - p**3
    assert sp.expand(node_relation.subs({p: node_p, q: node_q})) == 0
    assert sp.expand(cusp_relation.subs({p: cusp_p, q: cusp_q})) == 0

    node_root = T**2 - T - p
    cusp_root = T**2 - p
    assert sp.expand(node_root.subs({T: t, p: node_p})) == 0
    assert sp.expand(cusp_root.subs({T: t, p: cusp_p})) == 0
    assert sp.expand((p * T - q).subs({T: t, p: node_p, q: node_q})) == 0
    assert sp.expand((p * T - q).subs({T: t, p: cusp_p, q: cusp_q})) == 0

    return {
        "node": {
            "conductor_algebra": "Q[p,q]/(q^2-p*q-p^3)",
            "normalization": "p=t(t-1), q=t^2(t-1)",
            "conductor": "(p,q), pulling back to p*Q[t]",
            "marked_root_equations": ("T^2-T-p=0", "p*T-q=0"),
            "reconstruction_off_conductor": "T=q/p",
            "conductor_fiber": "T(T-1)=0: the two branches are glued",
            "normalization_quotient": (
                "(B/A) is cyclic on T-bar with annihilator (p,q)"
            ),
        },
        "cusp": {
            "conductor_algebra": "Q[p,q]/(q^2-p^3)",
            "normalization": "p=t^2, q=t^3",
            "conductor": "(p,q), pulling back to p*Q[t]",
            "marked_root_equations": ("T^2-p=0", "p*T-q=0"),
            "reconstruction_off_conductor": "T=q/p",
            "conductor_fiber": "T^2=0: the missing first jet is restored",
            "normalization_quotient": (
                "(B/A) is cyclic on T-bar with annihilator (p,q)"
            ),
        },
    }


def polynomiality_and_volume_ledger() -> dict[str, object]:
    """Verify descent and the conductor pole in the dualizing ledger."""

    product = x * y * z
    cases = {
        "node": {
            "c": t * (t - 1),
            "q": t**2 * (t - 1),
            "F": q**2 - product * q - product**3,
            "F_q": 2 * q - product,
            "boundary_orders": {
                "D_x,0; D_y,0; D_z,0": 1,
                "D_x,1; D_y,1; D_z,1": 1,
            },
        },
        "cusp": {
            "c": t**2,
            "q": t**3,
            "F": q**2 - product**3,
            "F_q": 2 * q,
            "boundary_orders": {"D_x; D_y; D_z": 2},
        },
    }
    output: dict[str, object] = {}
    for name, case in cases.items():
        c = case["c"]
        descended_q = case["q"]
        descended_equation = case["F"]
        derivative_q = case["F_q"]

        assert sp.expand(
            descended_equation.subs(
                {z: c / (x * y), q: descended_q}
            )
        ) == 0
        dz_dt = sp.diff(c, t) / (x * y)
        pulled_residue_coefficient = sp.cancel(
            dz_dt
            / derivative_q.subs(
                {z: c / (x * y), q: descended_q}
            )
        )
        normalized_residue_coefficient = 1 / (x * y)
        assert sp.cancel(
            pulled_residue_coefficient
            / normalized_residue_coefficient
            - 1 / c
        ) == 0
        assert sp.expand(product.subs(z, c / (x * y)) - c) == 0

        output[name] = {
            "normalized_fill": f"x*y*z={sp.sstr(c)}",
            "polynomial_descended_coordinates": {
                "p": "x*y*z",
                "q": sp.sstr(descended_q),
            },
            "descended_hypersurface": sp.sstr(descended_equation),
            "reconstruction": "t=q/(x*y*z) off the conductor",
            "dualizing_pullback": "pi^*Omega_desc=Omega_norm/(x*y*z)",
            "conductor_boundary_orders": case["boundary_orders"],
        }
    return output


def affine_space_recognition() -> dict[str, object]:
    """Compute motivic/Hodge ledgers and the cusp singularity obstruction."""

    zero_fiber = sp.expand(L**3 - (L - 1) ** 3)
    node_class = sp.expand((L - 2) * (L - 1) ** 2 + 2 * zero_fiber)
    cusp_class = sp.expand((L - 1) ** 3 + zero_fiber)
    assert node_class == L**3 + 2 * L**2 - L
    assert cusp_class == L**3

    # Singular independently verifies these gradient conclusions.
    node_c = t * (t - 1)
    node_gradient = (
        y * z,
        x * z,
        x * y,
        -sp.diff(node_c, t),
    )
    # At a hypothetical singular point 2t-1=0, the defining equation makes
    # xyz=-1/4, contradicting the first three gradient equations.
    assert sp.expand(node_c.subs(t, sp.Rational(1, 2))) == -sp.Rational(1, 4)
    assert node_gradient[-1] == 1 - 2 * t

    cusp_gradient = (y * z, x * z, x * y, -2 * t)
    assert all(
        entry.subs({t: 0, y: 0, z: 0}) == 0
        for entry in cusp_gradient
    )

    return {
        "node": {
            "smooth": True,
            "grothendieck_class": "L^3+2*L^2-L",
            "hodge_deligne_polynomial": "(uv)^3+2(uv)^2-uv",
            "affine_space": False,
            "stable_affine_space_after_polynomial_stabilization": False,
            "obstruction": "nonzero Hodge-Deligne defect 2(uv)^2-uv",
        },
        "cusp": {
            "normal": True,
            "smooth": False,
            "singular_locus": (
                "V(t,y,z) union V(t,x,z) union V(t,x,y)"
            ),
            "grothendieck_class": "L^3",
            "affine_space": False,
            "stable_affine_space_after_polynomial_stabilization": False,
            "obstruction": (
                "a codimension-two three-axis singular locus, persisting "
                "after polynomial stabilization"
            ),
        },
    }


def main() -> None:
    singular = run_singular()
    output = {
        "schema": "conductor-three-boundary-cox-fill.v1",
        "status": "exact obstruction for the minimal symmetric three-boundary fill",
        "singular_replay": singular,
        "conductor_and_marked_root": conductor_data(),
        "polynomiality_and_determinant_ledger": (
            polynomiality_and_volume_ledger()
        ),
        "affine_space_recognition": affine_space_recognition(),
        "theorem_boundary": (
            "This excludes the symmetric Cox fill x*y*z=p for the nodal "
            "and cuspidal rational conductor algebras, including ordinary "
            "polynomial stabilization.  It does not exclude asymmetric "
            "affine modifications, a non-Cox three-boundary relation, or "
            "a distributed source/target ledger."
        ),
        "reproducing_command": (
            ".venv/bin/python "
            "scripts/verify_conductor_three_boundary_cox_fill.py"
        ),
    }
    ARTIFACT.write_text(json.dumps(output, indent=2) + "\n")
    print("PASS: marked-root normalization and conductor descent")
    print("PASS: all three-boundary incidence coordinates are polynomial")
    print("PASS: the dualizing ledger retains exactly one conductor pole")
    print("PASS: neither normalized fill is affine three-space")
    print(f"PASS: wrote {ARTIFACT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
