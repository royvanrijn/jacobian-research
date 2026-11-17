#!/usr/bin/env python3
"""Exact bounded LNED searches for two principal conductor ideals.

For A=Q[u,v,x,y] and D=u*d/dx+v*d/dy, put w=u*y-v*x.  A target h in
im(D), with primitive P, belongs to D((x)) exactly when P modulo x lies in

    im(Q[u,v,w] -> A/(x)) = Q[u,v,u*y].

Thus every residue monomial u^a*v^b*y^c must satisfy a>=c.  This is an
exact valuation-face membership test with no degree truncation.  Pure
powers through six and mixed powers four through six remain a bounded
counterexample search.

For the rational-root carrier

    p = u*x+w = (u-v)*x+u*y,

normalize A/(p) inside Q[u,v,t] by x=u*t and y=-(u-v)*t.  The invariant
image is Q[u,v,u^2*t].  Thus a primitive enters (p) exactly when every
normalized residue monomial u^a*v^b*t^c has a>=2*c.  This is again exact
and untruncated.  The script also checks the exact local-slice primitive
identity for the two-branch carrier (x*y).
"""

from __future__ import annotations

import sympy as sp

import search_lnd_nonprincipal_plinth as npl


u, v, x, y = npl.u, npl.v, npl.x, npl.y
invariant = npl.invariant
w = sp.symbols("w")
t = sp.symbols("t")
p = sp.symbols("p")
rational_root_carrier = sp.expand(u * x + invariant)


def primitive_enters_x_ideal(primitive: sp.Expr | None) -> bool:
    if primitive is None:
        return False
    residue = sp.Poly(
        sp.expand(primitive.subs(x, 0)), u, v, y, domain=sp.QQ
    )
    return all(
        u_power >= y_power
        for (u_power, _, y_power), coefficient in residue.terms()
        if coefficient
    )


def rational_root_residue(
    primitive: sp.Expr | None,
) -> sp.Poly | None:
    if primitive is None:
        return None
    return sp.Poly(
        sp.expand(
            primitive.subs(
                {
                    x: u * t,
                    y: -(u - v) * t,
                }
            )
        ),
        u,
        v,
        t,
        domain=sp.QQ,
    )


def monomial_root_residue(
    primitive: sp.Expr | None,
    u_exponent: int,
    v_exponent: int,
) -> sp.Poly | None:
    if primitive is None:
        return None
    return sp.Poly(
        sp.expand(
            primitive.subs(
                {
                    x: u * t,
                    y: -(
                        u**u_exponent * v**v_exponent - v
                    )
                    * t,
                }
            )
        ),
        u,
        v,
        t,
        domain=sp.QQ,
    )


def primitive_enters_monomial_root_ideal(
    primitive: sp.Expr | None,
    u_exponent: int,
    v_exponent: int,
) -> bool:
    residue = monomial_root_residue(
        primitive, u_exponent, v_exponent
    )
    if residue is None:
        return False
    return all(
        u_power >= (u_exponent + 1) * t_power
        and v_power >= v_exponent * t_power
        for (u_power, v_power, t_power), coefficient in residue.terms()
        if coefficient
    )


def monomial_root_defects(
    primitive: sp.Expr | None,
    u_exponent: int,
    v_exponent: int,
) -> tuple[int, int] | None:
    residue = monomial_root_residue(
        primitive, u_exponent, v_exponent
    )
    if residue is None:
        return None
    defects = [
        (
            (u_exponent + 1) * t_power - u_power,
            v_exponent * t_power - v_power,
        )
        for (u_power, v_power, t_power), coefficient in residue.terms()
        if coefficient
    ]
    return (
        max((defect[0] for defect in defects), default=0),
        max((defect[1] for defect in defects), default=0),
    )


def primitive_enters_rational_root_ideal(
    primitive: sp.Expr | None,
) -> bool:
    residue = rational_root_residue(primitive)
    if residue is None:
        return False
    return all(
        u_power >= 2 * t_power
        for (u_power, _, t_power), coefficient in residue.terms()
        if coefficient
    )


def rational_root_defect(primitive: sp.Expr | None) -> int | None:
    residue = rational_root_residue(primitive)
    if residue is None:
        return None
    defects = [
        2 * t_power - u_power
        for (u_power, _, t_power), coefficient in residue.terms()
        if coefficient
    ]
    return max(defects, default=0)


def rational_root_margin(primitive: sp.Expr | None) -> int | None:
    residue = rational_root_residue(primitive)
    if residue is None:
        return None
    margins = [
        u_power - 2 * t_power
        for (u_power, _, t_power), coefficient in residue.terms()
        if coefficient
    ]
    return min(margins) if margins else None


def rational_root_local_form_after_u(f: sp.Expr) -> sp.Expr | None:
    if not all(
        u_power >= 1
        for (u_power, _, _, _), coefficient in sp.Poly(
            f, *npl.variables
        ).terms()
        if coefficient
    ):
        return None
    transformed = sp.cancel(
        (f / u).subs(
            {
                x: (p - w) / u,
                y: (v * p + (u - v) * w) / u**2,
            }
        )
    )
    try:
        return sp.Poly(
            transformed, u, v, w, p, domain=sp.QQ
        ).as_expr()
    except sp.PolynomialError:
        return None


def tied_root_local_form_after_uv(f: sp.Expr) -> sp.Expr | None:
    if not all(
        u_power >= 1 and v_power >= 1
        for (u_power, v_power, _, _), coefficient in sp.Poly(
            f, *npl.variables
        ).terms()
        if coefficient
    ):
        return None
    transformed = sp.cancel(
        (f / (u * v)).subs(
            {
                x: (p - w) / (u * v),
                y: (p - w) / u**2 + w / u,
            }
        )
    )
    try:
        return sp.Poly(
            transformed, u, v, w, p, domain=sp.QQ
        ).as_expr()
    except sp.PolynomialError:
        return None


def boundary_defect(primitive: sp.Expr | None) -> int | None:
    if primitive is None:
        return None
    residue = sp.Poly(
        sp.expand(primitive.subs(x, 0)), u, v, y, domain=sp.QQ
    )
    defects = [
        y_power - u_power
        for (u_power, _, y_power), coefficient in residue.terms()
        if coefficient
    ]
    return max(defects, default=0)


def boundary_margin(primitive: sp.Expr | None) -> int | None:
    if primitive is None:
        return None
    residue = sp.Poly(
        sp.expand(primitive.subs(x, 0)), u, v, y, domain=sp.QQ
    )
    margins = [
        u_power - y_power
        for (u_power, _, y_power), coefficient in residue.terms()
        if coefficient
    ]
    return min(margins) if margins else None


def divisible_by(f: sp.Expr, variable_index: int) -> bool:
    return all(
        powers[variable_index] >= 1
        for powers, coefficient in sp.Poly(
            f, *npl.variables
        ).terms()
        if coefficient
    )


def in_u_x_support(f: sp.Expr) -> bool:
    return all(
        u_power + x_power >= 1
        for (u_power, _, x_power, _), coefficient in sp.Poly(
            f, *npl.variables
        ).terms()
        if coefficient
    )


def kernel_x_form_after_u(f: sp.Expr) -> sp.Expr | None:
    if not divisible_by(f, 0):
        return None
    transformed = sp.cancel(
        (f / u).subs(y, (w + v * x) / u)
    )
    try:
        return sp.Poly(
            transformed, u, v, w, x, domain=sp.QQ
        ).as_expr()
    except sp.PolynomialError:
        return None


def verify_square_gate_failure() -> sp.Expr:
    """Return a binary face with zero first two moments, but nonzero third."""
    t, zeta = sp.symbols("t zeta")
    alpha = sp.sqrt(-15) / 3
    legendre_face = sp.expand(
        2 * t - 1 + alpha * (6 * t**2 - 6 * t + 1)
    )
    moments = tuple(
        sp.simplify(sp.integrate(legendre_face**power, (t, 0, 1)))
        for power in range(1, 4)
    )
    assert moments[0] == 0
    assert moments[1] == 0
    assert sp.simplify(
        moments[2] - 32 * sp.sqrt(-15) / 315
    ) == 0

    homogenized = sp.expand(
        w * (2 * zeta - w)
        + alpha * (6 * zeta**2 - 6 * w * zeta + w**2)
    )
    averaged = []
    integration_variable = sp.symbols("integration_variable")
    for power in range(1, 4):
        integrand = homogenized.subs(
            zeta, integration_variable
        ) ** power
        numerator = sp.integrate(
            integrand, (integration_variable, w, zeta)
        )
        average = sp.cancel(numerator / (zeta - w))
        averaged.append(average)
    assert sp.simplify(averaged[0].subs(zeta, 0)) == 0
    assert sp.simplify(averaged[1].subs(zeta, 0)) == 0
    assert sp.simplify(averaged[2].subs(zeta, 0)) != 0
    return legendre_face


def verify_two_branch_identity() -> None:
    """Check the zero-constant primitive identity for the carrier (x*y)."""
    integration_variable = sp.symbols("integration_variable")
    zeta = w + v * x
    multiplier = 1 + u + v + x + zeta / u + x * zeta / u
    primitive = sp.cancel(x * zeta * multiplier / u)
    target = sp.expand(u * sp.diff(primitive, x))
    averaged = sp.cancel(
        sp.integrate(
            target.subs(x, integration_variable),
            (integration_variable, 0, x),
        )
        / x
    )
    assert sp.cancel(averaged - zeta * multiplier) == 0


def verify_invariant_content_identities() -> None:
    """Check aligned and crossed invariant-content carrier identities."""
    multiplier = 1 + u + v + x + y + x * y
    assert sp.expand(
        npl.derivation(u * x * multiplier)
        - u * npl.derivation(x * multiplier)
    ) == 0
    assert sp.expand(
        npl.derivation(u * y * multiplier)
        - u * npl.derivation(y * multiplier)
    ) == 0


def verify_invariant_affine_coordinates() -> None:
    """Check the two-generator inverse chart for ker(D^2) carriers."""
    b0 = 1 + invariant
    b1 = u + invariant
    b2 = v + invariant**2
    carrier = sp.expand(b0 + b1 * x + b2 * y)
    speed = sp.expand(b1 * u + b2 * v)
    assert npl.derivation(carrier) == speed
    assert npl.derivation(speed) == 0
    assert sp.expand(
        u * (carrier - b0) - b2 * invariant - speed * x
    ) == 0
    assert sp.expand(
        v * (carrier - b0) + b1 * invariant - speed * y
    ) == 0


def verify_rational_root_normalization() -> None:
    """Check the rational-root ladder in its quotient normalization."""
    for exponent in range(1, 5):
        carrier = sp.expand(u**exponent * x + invariant)
        substitution = {
            x: u * t,
            y: -(u**exponent - v) * t,
        }
        assert sp.expand(carrier.subs(substitution)) == 0
        assert sp.expand(invariant.subs(substitution)) == (
            -u ** (exponent + 1) * t
        )
        assert npl.derivation(carrier) == u ** (exponent + 1)

    delta = sp.symbols("delta")
    y_in_local_coordinates = v * delta / u**2 + w / u
    second_face = sp.expand(u**2 * y_in_local_coordinates**3)
    assert sp.expand(second_face).coeff(delta, 0) == w**3 / u

    for u_exponent in range(1, 4):
        for v_exponent in range(1, 4):
            tied_carrier = sp.expand(
                u**u_exponent * v**v_exponent * x + invariant
            )
            tied_substitution = {
                x: u * t,
                y: -(
                    u**u_exponent * v**v_exponent - v
                )
                * t,
            }
            assert (
                sp.expand(tied_carrier.subs(tied_substitution)) == 0
            )
            assert sp.expand(
                invariant.subs(tied_substitution)
            ) == -u ** (u_exponent + 1) * v**v_exponent * t
            assert npl.derivation(tied_carrier) == (
                u ** (u_exponent + 1) * v**v_exponent
            )


def main() -> None:
    pure_bound = 6
    mixed_start = 4
    seeds = npl.seed_family()
    multipliers = (sp.Integer(1), u, v, x, y, invariant)

    pure_queries = tuple(
        sp.expand(seed**exponent)
        for seed in seeds
        for exponent in range(1, pure_bound + 1)
    )
    pure_primitives = npl.particular_primitives(pure_queries)
    square_gate_face = verify_square_gate_failure()
    verify_two_branch_identity()
    verify_invariant_content_identities()
    verify_invariant_affine_coordinates()
    verify_rational_root_normalization()
    membership_rows = [
        tuple(
            primitive_enters_x_ideal(
                pure_primitives[
                    seed_index * pure_bound + exponent - 1
                ]
            )
            for exponent in range(1, pure_bound + 1)
        )
        for seed_index in range(len(seeds))
    ]
    survivors = [
        seed
        for seed_index, seed in enumerate(seeds)
        if all(membership_rows[seed_index])
    ]
    survivor_indices = [
        seed_index
        for seed_index, row in enumerate(membership_rows)
        if all(row)
    ]
    survivor_margins = {
        sp.srepr(seeds[seed_index]): tuple(
            boundary_margin(
                pure_primitives[
                    seed_index * pure_bound + exponent - 1
                ]
            )
            for exponent in range(1, pure_bound + 1)
        )
        for seed_index in survivor_indices
    }

    prefix_lengths = [
        (
            next(
                (
                    exponent
                    for exponent, enters in enumerate(row)
                    if not enters
                ),
                pure_bound,
            ),
            seed,
        )
        for seed, row in zip(seeds, membership_rows, strict=True)
        if seed not in survivors
    ]
    best_failed_prefix = max(prefix for prefix, _ in prefix_lengths)
    nearest_failures = [
        seed
        for prefix, seed in prefix_lengths
        if prefix == best_failed_prefix
    ]

    mixed_queries = tuple(
        sp.expand(multiplier * seed**exponent)
        for seed in survivors
        for multiplier in multipliers
        for exponent in range(mixed_start, pure_bound + 1)
    )
    mixed_primitives = (
        npl.particular_primitives(mixed_queries)
        if mixed_queries
        else ()
    )
    tail_obstructions: list[
        tuple[sp.Expr, sp.Expr, tuple[int | None, ...]]
    ] = []
    query_index = 0
    for seed in survivors:
        for multiplier in multipliers:
            tail = mixed_primitives[query_index : query_index + 3]
            query_index += 3
            if tail and all(
                not primitive_enters_x_ideal(primitive)
                for primitive in tail
            ):
                tail_obstructions.append(
                    (
                        seed,
                        multiplier,
                        tuple(boundary_defect(p) for p in tail),
                    )
                )

    rational_membership_rows = [
        tuple(
            primitive_enters_rational_root_ideal(
                pure_primitives[
                    seed_index * pure_bound + exponent - 1
                ]
            )
            for exponent in range(1, pure_bound + 1)
        )
        for seed_index in range(len(seeds))
    ]
    rational_full_prefix_indices = [
        seed_index
        for seed_index, row in enumerate(rational_membership_rows)
        if all(row)
    ]
    rational_survivor_indices = [
        seed_index
        for seed_index, row in enumerate(rational_membership_rows)
        if all(row[mixed_start - 1 : pure_bound])
    ]
    rational_survivors = [
        seeds[seed_index] for seed_index in rational_survivor_indices
    ]
    rational_margins = {
        sp.srepr(seeds[seed_index]): tuple(
            rational_root_margin(
                pure_primitives[
                    seed_index * pure_bound + exponent - 1
                ]
            )
            for exponent in range(1, pure_bound + 1)
        )
        for seed_index in rational_survivor_indices
    }
    rational_prefix_lengths = [
        (
            next(
                (
                    exponent
                    for exponent, enters in enumerate(row)
                    if not enters
                ),
                pure_bound,
            ),
            seed,
        )
        for seed, row in zip(
            seeds, rational_membership_rows, strict=True
        )
        if not all(row)
    ]
    rational_best_failed_prefix = max(
        prefix for prefix, _ in rational_prefix_lengths
    )

    rational_mixed_queries = tuple(
        sp.expand(multiplier * seed**exponent)
        for seed in rational_survivors
        for multiplier in multipliers
        for exponent in range(mixed_start, pure_bound + 1)
    )
    rational_mixed_primitives = (
        npl.particular_primitives(rational_mixed_queries)
        if rational_mixed_queries
        else ()
    )
    rational_tail_obstructions: list[
        tuple[sp.Expr, sp.Expr, tuple[int | None, ...]]
    ] = []
    query_index = 0
    for seed in rational_survivors:
        for multiplier in multipliers:
            tail = rational_mixed_primitives[
                query_index : query_index + 3
            ]
            query_index += 3
            if tail and all(
                not primitive_enters_rational_root_ideal(primitive)
                for primitive in tail
            ):
                rational_tail_obstructions.append(
                    (
                        seed,
                        multiplier,
                        tuple(rational_root_defect(p) for p in tail),
                    )
                )

    tied_membership_rows = [
        tuple(
            primitive_enters_monomial_root_ideal(
                pure_primitives[
                    seed_index * pure_bound + exponent - 1
                ],
                1,
                1,
            )
            for exponent in range(1, pure_bound + 1)
        )
        for seed_index in range(len(seeds))
    ]
    tied_full_prefix_indices = [
        seed_index
        for seed_index, row in enumerate(tied_membership_rows)
        if all(row)
    ]
    tied_survivor_indices = [
        seed_index
        for seed_index, row in enumerate(tied_membership_rows)
        if all(row[mixed_start - 1 : pure_bound])
    ]
    tied_survivors = [
        seeds[seed_index] for seed_index in tied_survivor_indices
    ]
    tied_mixed_queries = tuple(
        sp.expand(multiplier * seed**exponent)
        for seed in tied_survivors
        for multiplier in multipliers
        for exponent in range(mixed_start, pure_bound + 1)
    )
    tied_mixed_primitives = (
        npl.particular_primitives(tied_mixed_queries)
        if tied_mixed_queries
        else ()
    )
    tied_tail_obstructions: list[
        tuple[
            sp.Expr,
            sp.Expr,
            tuple[tuple[int, int] | None, ...],
        ]
    ] = []
    query_index = 0
    for seed in tied_survivors:
        for multiplier in multipliers:
            tail = tied_mixed_primitives[
                query_index : query_index + 3
            ]
            query_index += 3
            if tail and all(
                not primitive_enters_monomial_root_ideal(
                    primitive, 1, 1
                )
                for primitive in tail
            ):
                tied_tail_obstructions.append(
                    (
                        seed,
                        multiplier,
                        tuple(
                            monomial_root_defects(primitive, 1, 1)
                            for primitive in tail
                        ),
                    )
                )

    print(
        "SEARCH:",
        "D=u*d_x+v*d_y; I=(x);",
        f"{len(seeds)} valuation-face seeds",
    )
    print(
        "SQUARE-GATE:",
        "first two moments vanish but the third does not for",
        f"p(t)={square_gate_face}",
    )
    print(
        "TWO-BRANCH:",
        "verified T(D(x*y*a))=(u*y)*a in the local slice",
    )
    print(
        "INVARIANT-CONTENT:",
        "verified D(u*x*a)=u*D(x*a) and D(u*y*a)=u*D(y*a)",
    )
    print(
        "INVARIANT-AFFINE:",
        "verified the inverse chart for h=b0+b1*x+b2*y",
    )
    print(
        "RATIONAL-ROOT:",
        "p=u*x+w; normalization x=u*t, y=-(u-v)*t;",
        "invariant image=Q[u,v,u^2*t]",
    )
    print(
        "SUMMARY:",
        f"pure-prefix survivors={len(survivors)}",
        f"bounded-tail obstructions={len(tail_obstructions)}",
        f"best failed prefix={best_failed_prefix}",
    )
    print(
        "SUPPORT:",
        f"divisible-by-u={sum(divisible_by(f, 0) for f in survivors)}",
        f"divisible-by-x={sum(divisible_by(f, 2) for f in survivors)}",
        f"in-(u,x)={sum(in_u_x_support(f) for f in survivors)}",
    )
    print(
        "BOUNDARY:",
        f"y-free={sum(not f.has(y) for f in survivors)}",
        f"invariant={sum(npl.derivation(f) == 0 for f in survivors)}",
        "slope-zero="
        f"{sum(0 in survivor_margins[sp.srepr(f)] for f in survivors)}",
    )
    print(
        "LOCAL-SLICE:",
        "in-u*ker(D)[x]="
        f"{sum(kernel_x_form_after_u(f) is not None for f in survivors)}",
    )
    for seed in survivors[:20]:
        print(
            "  SURVIVOR:",
            f"f={seed}",
            f"margins={survivor_margins[sp.srepr(seed)]}",
        )
    for seed in (f for f in survivors if f.has(y)):
        print(
            "  Y-DEPENDENT:",
            f"f={seed}",
            f"D(f)={npl.derivation(seed)}",
            f"f/u={kernel_x_form_after_u(seed)}",
            f"margins={survivor_margins[sp.srepr(seed)]}",
        )
    for seed in (
        f
        for f in survivors
        if 0 in survivor_margins[sp.srepr(f)]
    ):
        print(
            "  SLOPE-ZERO:",
            f"f={seed}",
            f"D(f)={npl.derivation(seed)}",
            f"margins={survivor_margins[sp.srepr(seed)]}",
        )
    for seed, multiplier, defects in tail_obstructions[:20]:
        print(
            "  CANDIDATE-ONLY:",
            f"f={seed}",
            f"g={multiplier}",
            f"valuation defects={defects}",
        )
    for seed in nearest_failures[:5]:
        print("  NEAREST-FAILURE:", f"f={seed}")
    print(
        "RATIONAL-SUMMARY:",
        f"full-prefix survivors={len(rational_full_prefix_indices)}",
        f"powers-{mixed_start}-through-{pure_bound} survivors="
        f"{len(rational_survivors)}",
        f"bounded-tail obstructions={len(rational_tail_obstructions)}",
        f"best failed prefix={rational_best_failed_prefix}",
    )
    print(
        "RATIONAL-LOCAL-CONE:",
        "in-u*ker(D)[p]="
        f"{sum(rational_root_local_form_after_u(f) is not None for f in rational_survivors)}",
    )
    for seed in rational_survivors[:30]:
        print(
            "  RATIONAL-SURVIVOR:",
            f"f={seed}",
            f"u^-1 form={rational_root_local_form_after_u(seed)}",
            f"margins={rational_margins[sp.srepr(seed)]}",
        )
    for seed, multiplier, defects in rational_tail_obstructions[:30]:
        print(
            "  RATIONAL-CANDIDATE-ONLY:",
            f"f={seed}",
            f"g={multiplier}",
            f"valuation defects={defects}",
        )
    print(
        "TIED-SUMMARY:",
        "q=u*v*x+w;",
        f"full-prefix survivors={len(tied_full_prefix_indices)}",
        f"powers-{mixed_start}-through-{pure_bound} survivors="
        f"{len(tied_survivors)}",
        f"bounded-tail obstructions={len(tied_tail_obstructions)}",
    )
    print(
        "TIED-LOCAL-CONE:",
        "in-u*v*ker(D)[q]="
        f"{sum(tied_root_local_form_after_uv(f) is not None for f in tied_survivors)}",
    )
    for seed in tied_survivors[:30]:
        print(
            "  TIED-SURVIVOR:",
            f"f={seed}",
            f"(u*v)^-1 form={tied_root_local_form_after_uv(seed)}",
        )
    for seed, multiplier, defects in tied_tail_obstructions[:30]:
        print(
            "  TIED-CANDIDATE-ONLY:",
            f"f={seed}",
            f"g={multiplier}",
            f"valuation defects={defects}",
        )
    assert len(survivors) == 48
    assert len(tail_obstructions) == 0
    assert best_failed_prefix == 1
    assert all(in_u_x_support(seed) for seed in survivors)
    assert sum(not seed.has(y) for seed in survivors) == 43
    assert sum(npl.derivation(seed) == 0 for seed in survivors) == 13
    assert sum(
        0 in survivor_margins[sp.srepr(seed)]
        for seed in survivors
    ) == 4
    assert all(
        kernel_x_form_after_u(seed) is not None
        for seed in survivors
    )
    assert len(rational_full_prefix_indices) == 9
    assert len(rational_survivors) == 17
    assert len(rational_tail_obstructions) == 0
    assert all(
        rational_root_local_form_after_u(seed) is not None
        for seed in rational_survivors
    )
    assert len(tied_full_prefix_indices) == 5
    assert len(tied_survivors) == 8
    assert len(tied_tail_obstructions) == 0
    assert all(
        tied_root_local_form_after_uv(seed) is not None
        for seed in tied_survivors
    )
    print("NOTE: membership is exact; the exponent range is bounded")


if __name__ == "__main__":
    main()
