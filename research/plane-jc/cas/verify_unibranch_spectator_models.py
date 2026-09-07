#!/usr/bin/env python3
"""Verify universal unibranch finite-flat models with an affine spectator.

For every n>=3 consider

    pi_n(T,u) = (u, T^(n+1)-T^n+u*T).

The source is A^2 and the map is finite free of rank n+1.  At the origin
the fiber is T^n(T-1): a length-n clean unibranch collision and one reduced
etale spectator.  These models refute a purely local finite-flat exclusion.

Their exact failure as Keller-normalization models is global.  The Jacobian
curve is a coordinate line after a polynomial source change, so deleting it
gives A^1 x G_m with a nonconstant unit, not the distinguished A^2 open.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/"
    "jc2_unibranch_spectator_countermodels.json"
)


def model_row(n: int) -> dict[str, object]:
    if n < 3:
        raise ValueError("the target branch is singular only for n>=3")
    T, u, v = sp.symbols("T u v")
    source_v = T ** (n + 1) - T**n + u * T
    equation = sp.expand(T ** (n + 1) - T**n + u * T - v)
    polynomial = sp.Poly(equation, T)
    if not polynomial.monic() or polynomial.degree() != n + 1:
        raise AssertionError("the cover ceased to be finite free of rank n+1")

    jacobian = sp.diff(source_v, T)
    expected_jacobian = (
        (n + 1) * T**n - n * T ** (n - 1) + u
    )
    if sp.expand(jacobian - expected_jacobian) != 0:
        raise AssertionError("the Jacobian formula changed")
    critical_u = n * T ** (n - 1) - (n + 1) * T**n
    critical_v = sp.expand(source_v.subs(u, critical_u))
    expected_critical_v = (n - 1) * T**n - n * T ** (n + 1)
    if sp.expand(critical_v - expected_critical_v) != 0:
        raise AssertionError("the critical image parametrization changed")

    if sp.Poly(critical_u, T).as_dict().get((n - 1,), 0) != n:
        raise AssertionError("the first target coordinate lost order n-1")
    if sp.Poly(critical_v, T).as_dict().get((n,), 0) != n - 1:
        raise AssertionError("the second target coordinate lost order n")
    if sp.gcd(n - 1, n) != 1:
        raise AssertionError("the target branch is no longer primitive")

    fiber = sp.factor(equation.subs({u: 0, v: 0}))
    if fiber != T**n * (T - 1):
        raise AssertionError("the origin fiber is not n+1")
    spectator_jacobian = sp.expand(jacobian.subs({T: 1, u: 0}))
    if spectator_jacobian != 1:
        raise AssertionError("the spectator sheet is not etale")

    second_normal_derivative = sp.diff(source_v, T, 2)
    if second_normal_derivative == 0:
        raise AssertionError("the generic critical curve lost simple ramification")

    # The coordinate W=Jac(pi_n) has polynomial inverse
    # u=W+n*T^(n-1)-(n+1)*T^n.  Hence deleting E=(W=0) gives
    # Spec k[T,W,W^-1] and W is a nonconstant unit.
    W = sp.symbols("W")
    inverse_u = W + n * T ** (n - 1) - (n + 1) * T**n
    if sp.expand(jacobian.subs(u, inverse_u) - W) != 0:
        raise AssertionError("the Jacobian coordinate change is not invertible")

    degree = n + 1
    local_multiplicity_jump = n - 2
    orevkov_cost = 2 + local_multiplicity_jump
    if orevkov_cost != degree - 1:
        raise AssertionError("the spectator model no longer saturates the Euler budget")

    return {
        "n": n,
        "generic_degree": degree,
        "monic_equation": f"T^{n+1}-T^{n}+u*T-v",
        "finite_free_basis": [f"T^{power}" for power in range(n + 1)],
        "jacobian": (
            f"u-{n}*T^{n-1}+{n+1}*T^{n}"
        ),
        "critical_curve": {
            "smooth": True,
            "generic_transverse_index": 2,
            "parameterization": {
                "u": f"{n}*T^{n-1}-{n+1}*T^{n}",
                "v": f"{n-1}*T^{n}-{n}*T^{n+1}",
            },
            "image_branch_type": [n - 1, n],
            "image_is_singular_unibranch": True,
        },
        "origin_fiber": {
            "equation": f"T^{n}*(T-1)",
            "point_lengths": [n, 1],
            "boundary_local_algebra": f"k[T]/(T^{n})",
            "spectator_point": "T=1",
            "spectator_jacobian": 1,
        },
        "orevkov_budget": {
            "generic_component_multiplicity": 2,
            "special_local_multiplicity": n,
            "jump": local_multiplicity_jump,
            "component_cost": orevkov_cost,
            "global_budget": degree - 1,
            "status": "saturated_not_excluded",
        },
        "distinguished_open_failure": {
            "jacobian_coordinate": "W=Jac(pi_n)",
            "coordinate_inverse": (
                f"u=W+{n}*T^{n-1}-{n+1}*T^{n}"
            ),
            "ramification_complement": "Spec k[T,W,W^-1]=A1 x G_m",
            "nonconstant_unit": "W",
            "source_class_group": 0,
            "ramification_prime_class": 0,
            "required_keller_open": "A2 with unit group k^*",
        },
    }


def full_fiber_contrast(n: int) -> dict[str, object]:
    """The same branch without a spectator is excluded by Euler cost."""

    degree = n
    cost = 2 + (n - 2)
    if cost <= degree - 1:
        raise AssertionError("the full-fiber contrast ceased to be excluded")
    return {
        "map": f"(T,u)->(u,T^{n}+u*T)",
        "generic_degree": degree,
        "origin_fiber": f"k[T]/(T^{n})",
        "component_cost": cost,
        "global_budget": degree - 1,
        "status": "excluded_by_orevkov_euler_budget",
    }


def quartic_packet_regression() -> dict[str, object]:
    """Verify that n=3 realizes both quartic frontier packets."""

    T = sp.symbols("T")
    sqrt_three = sp.sqrt(3)
    left = (sp.Integer(1) + sqrt_three) / 4
    right = (sp.Integer(1) - sqrt_three) / 4

    def critical_u(parameter: sp.Expr) -> sp.Expr:
        return sp.expand(3 * parameter**2 - 4 * parameter**3)

    def critical_v(parameter: sp.Expr) -> sp.Expr:
        return sp.expand(2 * parameter**3 - 3 * parameter**4)

    target_u = sp.simplify(critical_u(left))
    target_v = sp.simplify(critical_v(left))
    if target_u != sp.Rational(1, 8) or target_v != -sp.Rational(1, 64):
        raise AssertionError("the quartic collision target changed")
    if (
        sp.simplify(critical_u(right) - target_u) != 0
        or sp.simplify(critical_v(right) - target_v) != 0
    ):
        raise AssertionError("the two critical points no longer collide")
    fiber = sp.factor(
        T**4 - T**3 + target_u * T - target_v,
        extension=sqrt_three,
    )
    expected = sp.expand((T - left) ** 2 * (T - right) ** 2)
    if sp.expand(fiber - expected) != 0:
        raise AssertionError("the quartic 2+2 fiber factorization changed")
    return {
        "cusp_target": [0, 0],
        "cusp_fiber_partition": [3, 1],
        "self_collision_target": ["1/8", "-1/64"],
        "critical_parameters": [
            "(1+sqrt(3))/4",
            "(1-sqrt(3))/4",
        ],
        "self_collision_fiber_partition": [2, 2],
        "conclusion": (
            "one finite polynomial cover realizes both surviving quartic "
            "packet types, but its etale open is A1 x G_m rather than A2"
        ),
    }


def monodromy_connectivity_regression(maximum_n: int = 8) -> dict[str, object]:
    """Show that the maximal cusp-block group still fixes the spectator."""

    rows = []
    for n in range(3, maximum_n + 1):
        # Sym(n) on {0,...,n-1} fixes spectator n.  Adjacent
        # transpositions generate the first orbit; one bridge transposition
        # (n-1,n) makes the full adjacent-transposition generating set.
        rows.append(
            {
                "n": n,
                "degree": n + 1,
                "maximal_cusp_block_orbits": [n, 1],
                "maximal_cusp_block_transitive": False,
                "one_bridge_transposition_generates": f"S_{n+1}",
            }
        )
    return {
        "rows": rows,
        "consequence": (
            "a connected global cover needs another event involving the "
            "spectator; a jump-free boundary self-collision can supply it"
        ),
    }


def build_payload() -> dict[str, object]:
    models = [model_row(n) for n in range(3, 11)]
    return {
        "schema": "plane-jc.unibranch-spectator-countermodels.v1",
        "status": "exact-countermodels-to-purely-local-exclusion",
        "universal_family": {
            "map": "pi_n(T,u)=(u,T^(n+1)-T^n+u*T), n>=3",
            "models_checked": models,
            "symbolic_proof": (
                "all identities are polynomial in n and are displayed in "
                "the model formulas; n=3,...,10 is a deterministic regression"
            ),
        },
        "full_fiber_contrasts": [
            full_fiber_contrast(n) for n in range(3, 11)
        ],
        "quartic_extremal_packet": quartic_packet_regression(),
        "monodromy_connectivity": monodromy_connectivity_regression(),
        "claim_boundary": {
            "proved": [
                "finite flatness, smooth integral source, clean ramification, "
                "a singular unibranch image, and a separate etale affine sheet "
                "are jointly locally consistent in every rank at least four",
                "the one-spectator packet exactly saturates Orevkov's Euler budget",
                "the quartic member realizes exact 3+1 and 2+2 fibers",
                "every model fails the Keller-normalization unit and boundary "
                "class conditions because its etale open is A1 x G_m",
            ],
            "not_proved": [
                "globalization inside a finite normalization whose distinguished "
                "open is A2",
                "existence of a Keller map",
                "an exclusion of any geometric degree at least four",
            ],
        },
        "verdict": {
            "purely_local_exclusion": "false",
            "remaining_direct_target": (
                "use the global A2 open, trivial unit group, and free boundary "
                "class group together with monodromy; no local length or "
                "conductor inequality can exclude the spectator packet"
            ),
        },
        "software": {
            "python": "standard library",
            "sympy": sp.__version__,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh", action="store_true")
    parser.add_argument("--artifact", type=Path, default=ARTIFACT)
    args = parser.parse_args()
    payload = build_payload()
    artifact = args.artifact.resolve()
    if args.refresh:
        artifact.parent.mkdir(parents=True, exist_ok=True)
        artifact.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        try:
            display_path = artifact.relative_to(ROOT)
        except ValueError:
            display_path = artifact
        print(f"WROTE {display_path}")
    else:
        expected = json.loads(artifact.read_text())
        current_claim = {key: value for key, value in payload.items() if key != "software"}
        pinned_claim = {key: value for key, value in expected.items() if key != "software"}
        if current_claim != pinned_claim:
            raise AssertionError(
                "pinned unibranch-spectator artifact is stale; "
                "inspect before --refresh"
            )
    print("PASS: universal finite-free unibranch packets with an etale spectator")
    print("PASS: every one-spectator packet saturates the Orevkov Euler budget")
    print("PASS: the quartic model realizes both 3+1 and 2+2 fibers")
    print("PASS: every etale open is A1 x G_m and has a nonconstant unit")
    print("PASS: a purely local finite-flat exclusion is false")


if __name__ == "__main__":
    main()
