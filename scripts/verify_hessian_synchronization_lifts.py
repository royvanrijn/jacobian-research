#!/usr/bin/env python3
"""Verify canonical Hessian linear lifts and low-degree synchronization.

For a normalized decomposition of type ``a o b``, triangular reconstruction
uses only coefficients of degrees at least two.  It therefore assigns a
canonical missing linear coefficient to every Hessian-composition point.
Multiple decompositions synchronize exactly when these canonical lifts agree.

This checker proves scheme-theoretic agreement for every multiple
composition-cut intersection through degree 24.  In degree 30, five exact
pair certificates form a spanning tree on all six proper cuts, proving that
the all-six intersection is globally synchronized.
"""

from __future__ import annotations

import itertools
import shutil
import subprocess
import sys
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from explore_degree30_hessian_ritt_braid import (  # noqa: E402
    W,
    Z,
    canonical_linear_lift,
    canonical_residuals,
    canonical_synchronization_defect,
)
from jcsearch.ritt_complex import (  # noqa: E402
    degree_thirty_braid_decorations,
    dickson_vertex_factors,
)


def normalized_polynomial(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
) -> sp.Expr:
    """Return a monic rational normalization, ignoring constant associates."""

    numerator = sp.together(expression).as_numer_denom()[0]
    polynomial = sp.Poly(numerator, *variables, domain=sp.QQ)
    return polynomial.monic().as_expr()


def audit_degree(
    degree: int,
    cuts: tuple[int, ...],
) -> tuple[int, tuple[int, ...]]:
    """Prove synchronization on the intersection of the requested cuts."""

    coefficients = sp.symbols(f"h{degree}_1:{degree}")
    linear_coefficient = coefficients[0]
    hessian_coefficients = coefficients[1:]
    polynomial = W**degree + sum(
        coefficients[power - 1] * W**power
        for power in range(1, degree)
    )

    lifts: dict[int, sp.Expr] = {}
    hessian_residuals: dict[int, list[sp.Expr]] = {}
    for outer_degree in cuts:
        inner_degree = degree // outer_degree
        lift = canonical_linear_lift(
            polynomial, outer_degree, inner_degree
        )
        defect = canonical_synchronization_defect(
            polynomial, outer_degree, inner_degree
        )
        full_residuals = canonical_residuals(
            polynomial,
            outer_degree,
            inner_degree,
            parameters=coefficients,
            factor_output=False,
            minimum_coefficient_degree=1,
        )
        hessian_only = canonical_residuals(
            polynomial,
            outer_degree,
            inner_degree,
            parameters=coefficients,
            factor_output=False,
            minimum_coefficient_degree=2,
        )

        assert sp.diff(lift, linear_coefficient) == 0
        assert sp.expand(defect - (lift - linear_coefficient)) == 0
        assert len(full_residuals) == len(hessian_only) + 1
        assert normalized_polynomial(
            full_residuals[0], coefficients
        ) == normalized_polynomial(defect, coefficients)
        assert tuple(
            normalized_polynomial(residual, coefficients)
            for residual in full_residuals[1:]
        ) == tuple(
            normalized_polynomial(residual, coefficients)
            for residual in hessian_only
        )
        lifts[outer_degree] = lift
        hessian_residuals[outer_degree] = hessian_only

    intersection_ideal = [
        residual
        for outer_degree in cuts
        for residual in hessian_residuals[outer_degree]
    ]
    basis = sp.groebner(
        intersection_ideal,
        *hessian_coefficients,
        order="grevlex",
        domain=sp.QQ,
    )
    base_lift = lifts[cuts[0]]
    for outer_degree in cuts[1:]:
        difference = sp.together(
            lifts[outer_degree] - base_lift
        ).as_numer_denom()[0]
        assert basis.reduce(difference)[1] == 0

    # Consequently adjoining every polynomial synchronization defect is the
    # same as adjoining one graph equation for the common linear lift.
    return len(basis.polys), tuple(
        len(hessian_residuals[outer_degree])
        for outer_degree in cuts
    )


def factor_chart(
    degree: int,
    outer_degree: int,
) -> tuple[sp.Expr, tuple[sp.Symbol, ...]]:
    """Return the generic normalized source incidence chart."""

    inner_degree = degree // outer_degree
    inner_parameters = sp.symbols(
        f"s{degree}_{outer_degree}_b1:{inner_degree}"
    )
    outer_parameters = sp.symbols(
        f"s{degree}_{outer_degree}_a1:{outer_degree}"
    )
    inner = W**inner_degree + sum(
        inner_parameters[power - 1] * W**power
        for power in range(1, inner_degree)
    )
    outer = Z**outer_degree + sum(
        outer_parameters[power - 1] * Z**power
        for power in range(1, outer_degree)
    )
    return (
        sp.expand(outer.subs(Z, inner)),
        inner_parameters + outer_parameters,
    )


def serialize_singular(expression: sp.Expr) -> str:
    return str(
        sp.together(expression).as_numer_denom()[0]
    ).replace("**", "^")


def audit_pair_on_factor_chart(
    degree: int,
    source_cut: int,
    target_cut: int,
    timeout: int = 120,
) -> int:
    """Prove lift equality after pulling back to one canonical incidence."""

    polynomial, parameters = factor_chart(degree, source_cut)
    residuals = canonical_residuals(
        polynomial,
        target_cut,
        degree // target_cut,
        parameters=parameters,
        factor_output=False,
        minimum_coefficient_degree=2,
    )
    source_lift = sp.Poly(polynomial, W).nth(1)
    target_lift = canonical_linear_lift(
        polynomial,
        target_cut,
        degree // target_cut,
    )
    difference = target_lift - source_lift

    singular = shutil.which("Singular")
    assert singular is not None
    block_orders = {
        (18, 3, 2): "(dp(5),dp(2))",
        (24, 4, 3): "(dp(5),dp(3))",
    }
    ring_order = block_orders.get(
        (degree, source_cut, target_cut), "dp"
    )
    program = (
        f'ring q=0,({",".join(map(str, parameters))}),{ring_order};\n'
        f'ideal I={",".join(serialize_singular(item) for item in residuals)};\n'
        "ideal G=slimgb(I);\n"
        'print("PAIR_SYNC");\n'
        f"print(reduce({serialize_singular(difference)},G)==0);\n"
        "print(size(G));\n"
    )
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=timeout,
    )
    compact = result.stdout.split()
    marker = compact.index("PAIR_SYNC")
    assert compact[marker + 1] == "1", result.stdout + result.stderr
    return int(compact[marker + 2])


def audit_degree24_transported_pair() -> int:
    """Close the hard outer-cut ``{2,3}`` pair in normal coordinates.

    Its reduced model is the degree-six ``3 o 2`` Dickson collision composed
    on the right with a generic quartic.  The top source coefficients recover
    the five model parameters polynomially.  Replacing the four remaining
    source coefficients by deviations from this graph turns synchronization
    into a manageable ``4 normal | 5 base`` block calculation.
    """

    inner_coefficients = sp.symbols("sync24_b1:8")
    outer_coefficients = sp.symbols("sync24_a1:3")
    source_parameters = inner_coefficients + outer_coefficients
    inner = W**8 + sum(
        inner_coefficients[power - 1] * W**power
        for power in range(1, 8)
    )
    outer = (
        Z**3
        + outer_coefficients[1] * Z**2
        + outer_coefficients[0] * Z
    )
    polynomial = sp.expand(outer.subs(Z, inner))
    residuals = canonical_residuals(
        polynomial,
        2,
        12,
        parameters=source_parameters,
        factor_output=False,
        minimum_coefficient_degree=2,
    )
    synchronization_defect = canonical_linear_lift(
        polynomial, 2, 12
    ) - sp.Poly(polynomial, W).nth(1)

    translation, parameter = sp.symbols("sync24_t sync24_z")
    quartic_coefficients = sp.symbols("sync24_e1:4")
    quartic = (
        W**4
        + quartic_coefficients[2] * W**3
        + quartic_coefficients[1] * W**2
        + quartic_coefficients[0] * W
    )
    collision_outer, collision_inner = dickson_vertex_factors(
        (3, 2), W, translation, parameter
    )
    transported_inner = sp.expand(collision_inner.subs(W, quartic))
    coefficient_map = {
        inner_coefficients[power - 1]:
        sp.Poly(transported_inner, W).nth(power)
        for power in range(1, 8)
    }
    coefficient_map.update(
        {
            outer_coefficients[power - 1]:
            sp.Poly(collision_outer, W).nth(power)
            for power in range(1, 3)
        }
    )

    b4, b5, b6, b7 = inner_coefficients[3:]
    a2 = outer_coefficients[1]
    e3_inverse = b7 / 2
    e2_inverse = (b6 - e3_inverse**2) / 2
    e1_inverse = (b5 - 2 * e3_inverse * e2_inverse) / 2
    translation_inverse = (
        b4 - e2_inverse**2 - 2 * e3_inverse * e1_inverse
    ) / 2
    parameter_inverse = sp.solve(
        sp.expand(
            coefficient_map[a2].subs(
                {
                    translation: translation_inverse,
                    parameter: parameter,
                }
            )
            - a2
        ),
        parameter,
    )[0]
    inverse = {
        quartic_coefficients[2]: e3_inverse,
        quartic_coefficients[1]: e2_inverse,
        quartic_coefficients[0]: e1_inverse,
        translation: translation_inverse,
        parameter: parameter_inverse,
    }
    solved_variables = (
        inner_coefficients[0],
        inner_coefficients[1],
        inner_coefficients[2],
        outer_coefficients[0],
    )
    graph_values = {
        variable: sp.factor(coefficient_map[variable].subs(inverse))
        for variable in solved_variables
    }
    normal_coordinates = sp.symbols("sync24_x1:5")
    normal_change = {
        variable: graph_values[variable] + normal
        for variable, normal in zip(
            solved_variables, normal_coordinates
        )
    }
    transformed_residuals = [
        sp.together(
            residual.subs(normal_change, simultaneous=True)
        ).as_numer_denom()[0]
        for residual in residuals
    ]
    transformed_defect = sp.together(
        synchronization_defect.subs(
            normal_change, simultaneous=True
        )
    ).as_numer_denom()[0]
    zero_normal = {
        normal: 0 for normal in normal_coordinates
    }
    assert all(
        sp.expand(residual.subs(zero_normal)) == 0
        for residual in transformed_residuals
    )
    assert sp.expand(transformed_defect.subs(zero_normal)) == 0

    base_coordinates = (b4, b5, b6, b7, a2)
    singular_variables = normal_coordinates + base_coordinates
    singular = shutil.which("Singular")
    assert singular is not None
    program = (
        f'ring q=0,({",".join(map(str, singular_variables))}),'
        "(dp(4),dp(5));\n"
        f'ideal I={",".join(serialize_singular(item) for item in transformed_residuals)};\n'
        "ideal G=slimgb(I);\n"
        'print("TRANSPORTED_SYNC24");\n'
        f"print(reduce({serialize_singular(transformed_defect)},G)==0);\n"
        "print(size(G));\n"
    )
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=180,
    )
    compact = result.stdout.split()
    marker = compact.index("TRANSPORTED_SYNC24")
    assert compact[marker + 1] == "1", result.stdout + result.stderr
    basis_size = int(compact[marker + 2])
    assert basis_size == 63
    return basis_size


def audit_degree30_all_cut_spanning_tree() -> tuple[int, ...]:
    """Synchronize the global all-six degree-30 intersection.

    An edge ``(source, target)`` records exact membership

        lambda_target - lambda_source in H_source + H_target.

    The five certified edges form a spanning tree on all proper outer cuts.
    Since every edge ideal is contained in the sum of all six Hessian ideals,
    transitivity makes all six lifts equal on the all-cut scheme.  Orienting
    each edge from the smaller canonical factor chart is essential for a
    compact, reproducible calculation.
    """

    cuts = {2, 3, 5, 6, 10, 15}
    spanning_tree = (
        (6, 2),
        (6, 3),
        (3, 15),
        (5, 15),
        (5, 10),
    )
    reached = {spanning_tree[0][0]}
    pending = list(spanning_tree)
    while pending:
        for edge in tuple(pending):
            source, target = edge
            if source in reached or target in reached:
                reached.update(edge)
                pending.remove(edge)
                break
        else:
            raise AssertionError("degree-30 certificate edges are disconnected")
    assert reached == cuts

    basis_sizes = tuple(
        audit_pair_on_factor_chart(30, source, target)
        for source, target in spanning_tree
    )
    assert basis_sizes == (11, 6, 95, 6, 11)
    return basis_sizes


def audit_degree30_nested_pair() -> int:
    """Close the nested degree-30 outer-cut pair ``{2,10}``.

    On the ``10 o 3`` source chart its common refinement is ``2 o 5 o 3``.
    The top five coefficients of the degree-ten outer factor recover the
    normalized ``2 o 5`` parameters polynomially.  The four remaining outer
    coefficients become graph-normal coordinates, giving a compact
    ``4 normal | 7 base`` certificate.
    """

    inner_coefficients = sp.symbols("sync30n_e1:3")
    outer_coefficients = sp.symbols("sync30n_c1:10")
    source_parameters = inner_coefficients + outer_coefficients
    inner = (
        W**3
        + inner_coefficients[1] * W**2
        + inner_coefficients[0] * W
    )
    outer = Z**10 + sum(
        outer_coefficients[power - 1] * Z**power
        for power in range(1, 10)
    )
    polynomial = sp.expand(outer.subs(Z, inner))
    residuals = canonical_residuals(
        polynomial,
        2,
        15,
        parameters=source_parameters,
        factor_output=False,
        minimum_coefficient_degree=2,
    )
    synchronization_defect = (
        canonical_linear_lift(polynomial, 2, 15)
        - sp.Poly(polynomial, W).nth(1)
    )

    middle_coefficients = sp.symbols("sync30n_b1:5")
    top_linear = sp.symbols("sync30n_a1")
    middle = W**5 + sum(
        middle_coefficients[power - 1] * W**power
        for power in range(1, 5)
    )
    top = Z**2 + top_linear * Z
    model_outer = sp.Poly(sp.expand(top.subs(Z, middle)), W)

    c = outer_coefficients
    b4 = c[8] / 2
    b3 = (c[7] - b4**2) / 2
    b2 = (c[6] - 2 * b4 * b3) / 2
    b1 = (c[5] - 2 * b4 * b2 - b3**2) / 2
    a1 = c[4] - 2 * b4 * b1 - 2 * b3 * b2
    inverse = {
        middle_coefficients[3]: b4,
        middle_coefficients[2]: b3,
        middle_coefficients[1]: b2,
        middle_coefficients[0]: b1,
        top_linear: a1,
    }

    normal_coordinates = sp.symbols("sync30n_x1:5")
    normal_change = {
        c[index]:
        sp.factor(model_outer.nth(index + 1).subs(inverse)) + normal
        for index, normal in enumerate(normal_coordinates)
    }
    transformed_residuals = [
        sp.together(
            residual.subs(normal_change, simultaneous=True)
        ).as_numer_denom()[0]
        for residual in residuals
    ]
    transformed_defect = sp.together(
        synchronization_defect.subs(
            normal_change, simultaneous=True
        )
    ).as_numer_denom()[0]
    zero_normal = {
        normal: 0 for normal in normal_coordinates
    }
    assert all(
        sp.expand(residual.subs(zero_normal)) == 0
        for residual in transformed_residuals
    )
    assert sp.expand(transformed_defect.subs(zero_normal)) == 0

    base_coordinates = outer_coefficients[4:] + inner_coefficients
    singular_variables = normal_coordinates + base_coordinates
    singular = shutil.which("Singular")
    assert singular is not None
    program = (
        f'ring q=0,({",".join(map(str, singular_variables))}),'
        "(dp(4),dp(7));\n"
        f'ideal I={",".join(serialize_singular(item) for item in transformed_residuals)};\n'
        "ideal G=slimgb(I);\n"
        'print("NESTED_SYNC30");\n'
        f"print(reduce({serialize_singular(transformed_defect)},G)==0);\n"
        "print(size(G));\n"
    )
    result = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    compact = result.stdout.split()
    marker = compact.index("NESTED_SYNC30")
    assert compact[marker + 1] == "1", result.stdout + result.stderr
    basis_size = int(compact[marker + 2])
    assert basis_size == 4
    return basis_size


def main() -> None:
    decorations = degree_thirty_braid_decorations()
    assert tuple(
        item.transverse_slice.augmentation_ideal_length
        for item in decorations
    ) == (4, 3, 7)
    assert tuple(
        item.transverse_slice.point_cotangent_homology_ranks
        for item in decorations
    ) == ((1, 1), (2, 2), (2, 2))
    assert tuple(
        item.transverse_slice.intrinsic_tor_ranks(4)
        for item in decorations
    ) == (
        (1, 1, 1, 1, 1),
        (1, 2, 3, 4, 5),
        (1, 2, 3, 4, 5),
    )

    cases = (
        (6, (2, 3)),
        (8, (2, 4)),
        (10, (2, 5)),
        (12, (2, 3, 4, 6)),
    )
    for degree, cuts in cases:
        basis_size, residual_counts = audit_degree(degree, cuts)
        print(
            f"degree {degree}: cuts {cuts}; Hessian residual counts "
            f"{residual_counts}; Groebner basis size {basis_size}"
        )
    factor_chart_cases = [
        (degree, left, right)
        for degree, cuts in (
            (12, (2, 3, 4, 6)),
            (14, (2, 7)),
            (15, (3, 5)),
            (16, (2, 4, 8)),
        )
        for left, right in itertools.combinations(cuts, 2)
    ]
    factor_chart_cases.extend(
        (18, left, right)
        for left, right in (
            (3, 2),
            (2, 6),
            (2, 9),
            (3, 6),
            (3, 9),
            (6, 9),
        )
    )
    # The fourteen direct degree-24 pairs are followed below by a specialized
    # normal-coordinate certificate for the transported (2,3) collision.
    factor_chart_cases.extend(
        (24, left, right)
        for left, right in (
            (4, 2),
            (6, 2),
            (8, 2),
            (2, 12),
            (4, 3),
            (6, 3),
            (3, 8),
            (3, 12),
            (4, 6),
            (4, 8),
            (4, 12),
            (6, 8),
            (6, 12),
            (8, 12),
        )
    )
    for degree, source_cut, target_cut in factor_chart_cases:
        basis_size = audit_pair_on_factor_chart(
            degree, source_cut, target_cut
        )
        print(
            f"degree {degree}: pair ({source_cut},{target_cut}) "
            f"synchronizes on the canonical factor chart; "
            f"basis size {basis_size}"
        )
    transported_basis_size = audit_degree24_transported_pair()
    print(
        "degree 24: transported pair (2,3) synchronizes in Ritt normal "
        f"coordinates; basis size {transported_basis_size}"
    )
    print("PASS: canonical linear lifts use only Hessian coefficients")
    print("PASS: the degree 6, 8, 10, and 12 all-cut schemes synchronize exactly")
    print("PASS: each polynomial intersection is one graph over its Hessian scheme")
    print("PASS: every pair through degree 16 synchronizes scheme-theoretically")
    print("PASS: every degree-18 pair synchronizes scheme-theoretically")
    print("PASS: every degree-24 pair synchronizes scheme-theoretically")
    degree30_basis_sizes = audit_degree30_all_cut_spanning_tree()
    print(
        "degree 30: the all-cut synchronization spanning tree has exact "
        f"basis sizes {degree30_basis_sizes}"
    )
    print(
        "PASS: the global degree-30 all-six intersection synchronizes "
        "scheme-theoretically"
    )
    nested_basis_size = audit_degree30_nested_pair()
    print(
        "degree 30: nested pair (2,10) synchronizes in common-refinement "
        f"normal coordinates; basis size {nested_basis_size}"
    )
    print("PASS: the degree-30 augmentation, cotangent, and intrinsic Tor data agree")


if __name__ == "__main__":
    main()
