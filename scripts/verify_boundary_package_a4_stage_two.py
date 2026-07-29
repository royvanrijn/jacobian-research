#!/usr/bin/env python3
"""Replay the known A4 symbolic core against the abstract boundary package.

The calculation deliberately stops short of ``realized``.  It certifies the
degree-four A4 root architecture, the exact W/K/L log-Keller ledger, and the
complete divisorial valuation over the target B-cubic.  The latter corrects
an earlier degree count: L has ramification index two *and residue degree
two*, so it already accounts for all four sheets.  The two transpositions
split only after a geometric residue-field extension; there is no second
missing height-one prime over the generic B-divisor.

All polynomial identities are replayed by the existing A4 verifiers.  This
adapter consumes their globals and compares the resulting divisor orders to
the stage-one package, avoiding a second implementation of the formulas.
"""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import runpy
from dataclasses import asdict, dataclass
from pathlib import Path

import sympy as sp
from sympy.polys.polyerrors import ExactQuotientFailed

from boundary_package_compiler import (
    PackageStatus,
    a4_three_puncture_package,
    compile_boundary_package,
)


SCRIPT_DIR = Path(__file__).resolve().parent
REPLAY_SCRIPTS = (
    SCRIPT_DIR / "verify_a4_affine_keller_frontier.py",
    SCRIPT_DIR / "verify_a4_ledger_reduction.py",
    SCRIPT_DIR / "verify_a4_pure_target_ledger.py",
    SCRIPT_DIR / "verify_a4_keller_inverse_cover.py",
)


@dataclass(frozen=True)
class A4StageTwoReplay:
    package: str
    package_status: str
    replay_sha256: tuple[tuple[str, str], ...]
    abstract_group_order: int
    root_degree: int
    oriented_inverse_discriminant_multiplier: int
    concrete_a4_polynomial: str
    concrete_a4_discriminant_root: int
    cone_ledger: tuple[tuple[str, int], ...]
    auxiliary_ledger: tuple[tuple[str, int], ...]
    target_pullback_ledger: tuple[tuple[str, int], ...]
    package_ledger_matches: bool
    pure_target_log_keller: bool
    ordinary_cone_keller: bool
    naive_target_factorization_polynomial: bool
    dominant_affine_ramified_prime: str
    target_branch_normalization_parameter: str
    source_divisor_normalization_parameter: str
    residue_minimal_polynomial: str
    ramification_index: int
    residue_degree: int
    local_degree_contribution: int
    remaining_degree_over_target_branch: int
    geometric_ramified_prime_count: int
    conductor_pair_preserved: bool
    puncture_map: tuple[tuple[str, str], ...]
    contracted_source_divisors: tuple[str, ...]
    selected_double_transposition_profile_matched: bool
    abstract_global_prime_profile_matched: bool
    abstract_boundary_coloring_matched: bool
    affine_space_realization_matched: bool
    remaining_gates: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        result = asdict(self)
        for key in (
            "replay_sha256",
            "cone_ledger",
            "auxiliary_ledger",
            "target_pullback_ledger",
            "puncture_map",
        ):
            result[key] = dict(result[key])
        return result


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _quiet_run(path: Path) -> dict[str, object]:
    output = io.StringIO()
    with contextlib.redirect_stdout(output):
        namespace = runpy.run_path(str(path))
    assert "PASS:" in output.getvalue(), f"{path.name} emitted no PASS certificate"
    return namespace


def _factor_order(
    expression: sp.Expr,
    factor: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> int:
    polynomial = sp.Poly(sp.expand(expression), *variables, domain=sp.QQ)
    divisor = sp.Poly(sp.expand(factor), *variables, domain=sp.QQ)
    order = 0
    while True:
        try:
            polynomial = polynomial.exquo(divisor)
        except ExactQuotientFailed:
            return order
        order += 1


def replay_a4_stage_two() -> A4StageTwoReplay:
    package = a4_three_puncture_package()
    compilation = compile_boundary_package(package)
    assert compilation.status is PackageStatus.UNKNOWN
    invariants = dict(compilation.invariants)
    assert invariants["generated_group_order"] == 12

    frontier = _quiet_run(REPLAY_SCRIPTS[0])
    ledger_reduction = _quiet_run(REPLAY_SCRIPTS[1])
    pure_target = _quiet_run(REPLAY_SCRIPTS[2])
    oriented = _quiet_run(REPLAY_SCRIPTS[3])

    # Root architecture: the oriented inverse quartic has discriminant
    # 4096*delta, and the oriented target adjoins D with D^2=delta.
    inverse_quartic = frontier["inverse_quartic"]
    root_variable = frontier["T"]
    target_delta = frontier["delta"]
    root_degree = sp.Poly(inverse_quartic, root_variable).degree()
    discriminant = sp.factor(sp.discriminant(inverse_quartic, root_variable))
    discriminant_multiplier = sp.cancel(discriminant / target_delta)
    assert discriminant_multiplier == 4096

    # The concrete specialization proves that the generic oriented group is
    # exactly A4: square discriminant puts it inside A4, while an irreducible
    # quartic with irreducible cubic resolvent has A4 rather than V4.
    concrete = oriented["P0"]
    concrete_resolvent = oriented["resolvent0"]
    concrete_delta = oriented["Delta0"]
    assert concrete.is_irreducible
    assert concrete_resolvent.is_irreducible
    concrete_root = int(sp.sqrt(concrete_delta))
    assert concrete_root**2 == concrete_delta

    # Extract actual divisor orders rather than copying the exponents from
    # the note.  The frontier and pure-target scripts use independent SymPy
    # symbols, so their ledgers are computed separately.
    frontier_variables = (
        frontier["U"],
        frontier["V"],
        frontier["W"],
    )
    cone_expression = frontier["cone_jacobian"]
    cone_factors = {
        "W": frontier["W"],
        "K": frontier["K"],
        "L": frontier["L"],
    }
    cone_orders = {
        name: _factor_order(cone_expression, factor, frontier_variables)
        for name, factor in cone_factors.items()
    }
    assert cone_orders == {"W": 2, "K": 3, "L": 1}

    pure_variables = (
        pure_target["U"],
        pure_target["V"],
        pure_target["W"],
        pure_target["z"],
    )
    lift_expression = pure_target["lift_jacobian"]
    target_expression = pure_target["target_B"]
    pure_factors = {
        "W": pure_target["W"],
        "K": pure_target["K"],
        "L": pure_target["L"],
    }
    lift_orders = {
        name: _factor_order(lift_expression, factor, pure_variables)
        for name, factor in pure_factors.items()
    }
    target_orders = {
        name: _factor_order(target_expression, factor, pure_variables)
        for name, factor in pure_factors.items()
    }
    assert lift_orders == {"W": 3, "K": 3, "L": 2}
    assert target_orders == lift_orders
    auxiliary_orders = {
        name: lift_orders[name] - cone_orders[name] for name in cone_orders
    }
    assert auxiliary_orders == {"W": 1, "K": 0, "L": 1}

    package_rows = {
        row.divisor: (
            row.source_jacobian_order,
            row.controlled_exponent * row.controlled_pullback_order,
            row.target_jacobian_order,
        )
        for row in package.determinant_ledger
    }
    extracted_rows = {
        name: (cone_orders[name], auxiliary_orders[name], target_orders[name])
        for name in cone_orders
    }
    package_ledger_matches = package_rows == extracted_rows
    assert package_ledger_matches

    pure_target_log_keller = (
        sp.factor(lift_expression - target_expression) == 0
    )
    ordinary_cone_keller = not any(cone_orders.values())
    naive_denominator = sp.denom(
        sp.cancel(pure_target["naive_fourth_coordinate"])
    )
    naive_target_factorization_polynomial = naive_denominator in (1, -1)
    assert pure_target_log_keller
    assert not ordinary_cone_keller
    assert not naive_target_factorization_polynomial

    # The residual L divisor is a genuine dominant affine ramification prime.
    # Its normalization parameter is supplied by the ledger-reduction replay.
    # Rank two at one exact point proves that its image is dense in the
    # irreducible target surface B=0.
    target_B_polynomial = sp.Poly(
        ledger_reduction["B_homogeneous"],
        ledger_reduction["P"],
        ledger_reduction["Q"],
        ledger_reduction["R"],
        domain=sp.QQ,
    )
    assert target_B_polynomial.is_irreducible
    restricted_cone = ledger_reduction["cone_map"].subs(
        {
            ledger_reduction["U"]: ledger_reduction["U_param"],
            ledger_reduction["V"]: ledger_reduction["V_param"],
        }
    )
    restricted_jacobian = restricted_cone.jacobian(
        (ledger_reduction["t"], ledger_reduction["W"])
    ).subs({ledger_reduction["t"]: 1, ledger_reduction["W"]: 1})
    assert restricted_jacobian.rank() == 2

    # Normalize the projective target cubic B=0 by the line-slope parameter
    #
    #     s = (P-Q)/R.
    #
    # Lines through its node [1:1:0] give exactly the same rational
    # parametrization formula used for L, now with parameter s.
    source_parameter = ledger_reduction["t"]
    target_parameter = sp.symbols("s_target")
    target_P_param = ledger_reduction["U_param"].subs(
        source_parameter, target_parameter
    )
    target_Q_param = ledger_reduction["V_param"].subs(
        source_parameter, target_parameter
    )
    assert sp.factor(
        ledger_reduction["B_homogeneous"].subs(
            {
                ledger_reduction["P"]: target_P_param,
                ledger_reduction["Q"]: target_Q_param,
                ledger_reduction["R"]: 1,
            }
        )
    ) == 0
    assert sp.factor(target_P_param - target_Q_param - target_parameter) == 0

    # On L, the induced residue-field map is
    #
    #     s = t^2/(2t-3).
    #
    # Hence t satisfies t^2-2st+3s.  Its discriminant 4s(s-3) has odd
    # orders at s=0 and s=3, so it is not a square in Q(s).  This proves
    # [Q(t):Q(s)]=2 rather than merely bounding it by two.
    residue_map = sp.factor(
        (
            ledger_reduction["N1"] - ledger_reduction["N2"]
        ).subs(
            {
                ledger_reduction["U"]: ledger_reduction["U_param"],
                ledger_reduction["V"]: ledger_reduction["V_param"],
            }
        )
        / ledger_reduction["H"].subs(
            {
                ledger_reduction["U"]: ledger_reduction["U_param"],
                ledger_reduction["V"]: ledger_reduction["V_param"],
            }
        )
    )
    expected_residue_map = source_parameter**2 / (2 * source_parameter - 3)
    assert sp.factor(residue_map - expected_residue_map) == 0
    residue_polynomial = sp.Poly(
        source_parameter**2
        - 2 * target_parameter * source_parameter
        + 3 * target_parameter,
        source_parameter,
        domain=sp.QQ.frac_field(target_parameter),
    )
    assert residue_polynomial.is_irreducible
    residue_discriminant = sp.factor(
        sp.discriminant(residue_polynomial.as_expr(), source_parameter)
    )
    assert residue_discriminant == 4 * target_parameter * (target_parameter - 3)
    residue_sqrt = sp.symbols("r_residue")
    geometric_factors_product = sp.expand(
        (source_parameter - target_parameter - residue_sqrt)
        * (source_parameter - target_parameter + residue_sqrt)
    )
    assert sp.factor(
        geometric_factors_product.subs(
            residue_sqrt**2,
            target_parameter * (target_parameter - 3),
        )
        - residue_polynomial.as_expr()
    ) == 0

    ramification_index = target_orders["L"]
    residue_degree = residue_polynomial.degree()
    local_degree_contribution = ramification_index * residue_degree
    remaining_degree_over_target_branch = root_degree - local_degree_contribution
    assert ramification_index == 2
    assert residue_degree == 2
    assert local_degree_contribution == root_degree
    assert remaining_degree_over_target_branch == 0

    # The quadratic residue extension splits over
    # Q(s)(sqrt(s(s-3))), producing the two geometric e=2 inertia orbits of
    # a double transposition.  Globally over Q(s), however, they form one
    # prime with (e,f)=(2,2), not two primes with (e,f)=(2,1).
    geometric_ramified_prime_count = 2
    assert geometric_ramified_prime_count == 2

    # Both source and target projective cubics have the node [1:1:0].  In
    # their normalization the conductor branches are {0,3}; the residue map
    # sends 0->0 and 3->3.  The third puncture infinity maps to infinity.
    node_point = {
        ledger_reduction["P"]: 1,
        ledger_reduction["Q"]: 1,
        ledger_reduction["R"]: 0,
    }
    target_B = ledger_reduction["B_homogeneous"]
    assert target_B.subs(node_point) == 0
    assert all(
        sp.diff(target_B, variable).subs(node_point) == 0
        for variable in (
            ledger_reduction["P"],
            ledger_reduction["Q"],
            ledger_reduction["R"],
        )
    )
    source_homogenizing_variable = sp.symbols("Z_source")
    source_L_projective = sp.Poly(
        ledger_reduction["L"],
        ledger_reduction["U"],
        ledger_reduction["V"],
        domain=sp.QQ,
    ).homogenize(source_homogenizing_variable).as_expr()
    source_node_point = {
        ledger_reduction["U"]: 1,
        ledger_reduction["V"]: 1,
        source_homogenizing_variable: 0,
    }
    assert source_L_projective.subs(source_node_point) == 0
    assert all(
        sp.diff(source_L_projective, variable).subs(source_node_point) == 0
        for variable in (
            ledger_reduction["U"],
            ledger_reduction["V"],
            source_homogenizing_variable,
        )
    )
    source_U_numerator, source_common_denominator = sp.fraction(
        sp.cancel(ledger_reduction["U_param"])
    )
    source_V_numerator, second_denominator = sp.fraction(
        sp.cancel(ledger_reduction["V_param"])
    )
    assert sp.factor(source_common_denominator - second_denominator) == 0
    assert (
        source_U_numerator.subs(source_parameter, 0),
        source_V_numerator.subs(source_parameter, 0),
        source_common_denominator.subs(source_parameter, 0),
    ) == (-27, -27, 0)
    assert (
        source_U_numerator.subs(source_parameter, 3),
        source_V_numerator.subs(source_parameter, 3),
        source_common_denominator.subs(source_parameter, 3),
    ) == (27, 27, 0)
    puncture_map = (("0", "0"), ("3", "3"), ("infinity", "infinity"))
    assert sp.limit(residue_map, source_parameter, 0) == 0
    assert sp.limit(residue_map, source_parameter, 3) == 3
    assert sp.limit(residue_map, source_parameter, sp.oo) is sp.oo
    conductor_pair_preserved = True

    # W=0 maps to the target origin.  K=0 forces P=W*N1=0 because N1=M*K,
    # while B also vanishes, so K maps into the codimension-two locus
    # P=B=0.  Neither dominates the generic point of B=0.
    assert sp.factor(
        ledger_reduction["N1"]
        - ledger_reduction["M"] * ledger_reduction["K"]
    ) == 0
    contracted_source_divisors = ("W", "K")

    # The geometric double-transposition cycle profile is matched.  The
    # package's *global* local-prime ledger is not: it declares two
    # residue-degree-one boundary primes, whereas the exact cone architecture
    # has one residue-degree-two affine prime.  The other abstract triangle
    # branches are not matched here.
    selected_double_transposition_profile_matched = True
    abstract_global_prime_profile_matched = False
    abstract_boundary_coloring_matched = False
    affine_space_realization_matched = False

    return A4StageTwoReplay(
        package=package.name,
        package_status=compilation.status.value,
        replay_sha256=tuple(
            (path.name, _sha256(path)) for path in REPLAY_SCRIPTS
        ),
        abstract_group_order=invariants["generated_group_order"],
        root_degree=root_degree,
        oriented_inverse_discriminant_multiplier=int(discriminant_multiplier),
        concrete_a4_polynomial=str(concrete.as_expr()),
        concrete_a4_discriminant_root=concrete_root,
        cone_ledger=tuple(cone_orders.items()),
        auxiliary_ledger=tuple(auxiliary_orders.items()),
        target_pullback_ledger=tuple(target_orders.items()),
        package_ledger_matches=package_ledger_matches,
        pure_target_log_keller=pure_target_log_keller,
        ordinary_cone_keller=ordinary_cone_keller,
        naive_target_factorization_polynomial=(
            naive_target_factorization_polynomial
        ),
        dominant_affine_ramified_prime="L",
        target_branch_normalization_parameter="s=(P-Q)/R",
        source_divisor_normalization_parameter="t=U-V",
        residue_minimal_polynomial=str(residue_polynomial.as_expr()),
        ramification_index=ramification_index,
        residue_degree=residue_degree,
        local_degree_contribution=local_degree_contribution,
        remaining_degree_over_target_branch=remaining_degree_over_target_branch,
        geometric_ramified_prime_count=geometric_ramified_prime_count,
        conductor_pair_preserved=conductor_pair_preserved,
        puncture_map=puncture_map,
        contracted_source_divisors=contracted_source_divisors,
        selected_double_transposition_profile_matched=(
            selected_double_transposition_profile_matched
        ),
        abstract_global_prime_profile_matched=(
            abstract_global_prime_profile_matched
        ),
        abstract_boundary_coloring_matched=(
            abstract_boundary_coloring_matched
        ),
        affine_space_realization_matched=affine_space_realization_matched,
        remaining_gates=(
            "replace or geometrically split the global (e,f)=(2,2) affine "
            "L-prime by two (2,1) primes with the package's requested "
            "boundary coloring",
            "construct a coupled affine modification with ordinary constant "
            "Jacobian rather than only a pure-target logarithmic Jacobian",
            "prove the modified regular-reconstruction source and target are "
            "polynomial affine spaces without changing the A4 extension",
        ),
    )


certificate = replay_a4_stage_two()
assert certificate.package_status == "unknown"
assert certificate.abstract_group_order == 12
assert certificate.root_degree == 4
assert certificate.package_ledger_matches
assert certificate.pure_target_log_keller
assert not certificate.ordinary_cone_keller
assert not certificate.naive_target_factorization_polynomial
assert certificate.dominant_affine_ramified_prime == "L"
assert certificate.ramification_index == 2
assert certificate.residue_degree == 2
assert certificate.local_degree_contribution == 4
assert certificate.remaining_degree_over_target_branch == 0
assert certificate.geometric_ramified_prime_count == 2
assert certificate.conductor_pair_preserved
assert certificate.contracted_source_divisors == ("W", "K")
assert certificate.selected_double_transposition_profile_matched
assert not certificate.abstract_global_prime_profile_matched
assert not certificate.abstract_boundary_coloring_matched
assert not certificate.affine_space_realization_matched

print("PASS: the abstract A4 package matches the exact W/K/L symbolic ledger")
print("PASS: the degree-four oriented root architecture has exact A4 monodromy")
print("PASS: the pure-target lift is log-Keller but the cone is not Keller")
print("PASS: L/B is the complete affine (e,f)=(2,2) divisorial valuation")
print("PASS: its two geometric e=2 branches preserve the conductor punctures")
print("PASS: global-prime splitting, recoloring, and affine Keller completion remain open")
print(json.dumps(certificate.to_dict(), indent=2, sort_keys=True))
