#!/usr/bin/env python3
"""Classify the 35 missing F2 (75,125) Laurent layers at the B0 boundary.

The output is an exact, compressed polynomial-system certificate.  It uses
source-band variables rather than expanding millions of quadratic terms.
For an original monomial on band ``ell=5*i-j``, the Puiseux translation
``y -> y+lambda/X`` (with lambda scaled to one) gives

    t^ell (1+t)^j z^ell,  j=5*i-ell.

The terminal edge imposes an exact initial-jet vanishing condition.  If
``p=t^ell(1+t)^r`` and ``q=t^m(1+t)^s``, the corrected coordinate bracket
``[t,z]_(X,y)=-z`` gives the compressed quadratic generator

    ell*p*q' - m*p'*q
      = (ell*s-m*r) t^(ell+m) (1+t)^(r+s-1).

Thus every scalar coefficient equation and every Kummer sector can be
reconstructed from the emitted generator records without a heuristic support
choice.  This is the finite B0 envelope plus the universal Keller recurrence;
it is not an exhaustive Newton polygon or gamma-branch theorem.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
from math import comb
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[2]
ARTIFACT = (
    ROOT
    / "artifacts/generated-results/jc2_f2_75_125_character_layers.json"
)


@dataclass(frozen=True)
class Band:
    side: str
    degree: int
    terminal_height: int
    layer: int
    source_i_min: int
    source_i_max: int
    source_j: tuple[int, ...]
    jet_vanishing_order: int
    dimension: int
    t_exponent_min: int
    t_exponent_max: int
    active_t_exponents: tuple[int, ...]

    @property
    def source_count(self) -> int:
        return self.source_i_max - self.source_i_min + 1

    def variable(self, source_i: int) -> str:
        return f"{self.side}_{self.layer}_{source_i}"


def ceil_div(numerator: int, denominator: int) -> int:
    return -((-numerator) // denominator)


def rational_rank(rows: list[list[int]]) -> int:
    if not rows:
        return 0
    return int(sp.Matrix(rows).rank())


def make_band(side: str, degree: int, height: int, layer: int) -> Band:
    source_i_min = max(0, ceil_div(layer, 5))
    source_i_max = (degree + layer) // 6
    if source_i_max < source_i_min:
        raise AssertionError(f"empty {side} band {layer}")
    source_j = tuple(
        5 * source_i - layer
        for source_i in range(source_i_min, source_i_max + 1)
    )
    t_exponent_min = max(layer, ceil_div(17 * layer - height, 12))
    jet_order = t_exponent_min - layer
    if jet_order > len(source_j):
        raise AssertionError(f"terminal edge kills {side} band {layer}")
    jet_rows = [
        [comb(source_power, order) for source_power in source_j]
        for order in range(jet_order)
    ]
    jet_rank = rational_rank(jet_rows)
    if jet_rank != jet_order:
        raise AssertionError(f"jet matrix lost rank on {side} band {layer}")
    active_exponents: list[int] = []
    for order in range(jet_order, max(source_j) + 1):
        row = [comb(source_power, order) for source_power in source_j]
        if rational_rank(jet_rows + [row]) > jet_rank:
            active_exponents.append(layer + order)
    expected_dimension = len(source_j) - jet_order
    if expected_dimension <= 0:
        raise AssertionError(f"unexpected zero-dimensional {side} band {layer}")
    if active_exponents != list(
        range(t_exponent_min, layer + max(source_j) + 1)
    ):
        raise AssertionError(f"noninterval coefficient support on {side} band {layer}")
    return Band(
        side=side,
        degree=degree,
        terminal_height=height,
        layer=layer,
        source_i_min=source_i_min,
        source_i_max=source_i_max,
        source_j=source_j,
        jet_vanishing_order=jet_order,
        dimension=expected_dimension,
        t_exponent_min=t_exponent_min,
        t_exponent_max=layer + max(source_j),
        active_t_exponents=tuple(active_exponents),
    )


def coordinate_bracket_regression() -> dict[str, object]:
    X, y, t, z = sp.symbols("X y t z")
    t_xy = X * y
    z_xy = y**-1
    coordinate_bracket = sp.simplify(
        sp.diff(t_xy, X) * sp.diff(z_xy, y)
        - sp.diff(t_xy, y) * sp.diff(z_xy, X)
    )
    if coordinate_bracket != -z_xy:
        raise AssertionError("the coordinate bracket is not -z")

    p_bands = {2: t**3 + 2 * t, -1: t**2}
    q_bands = {3: t**2 - t, 0: 3 * t}
    P_tz = sum(polynomial * z**layer for layer, polynomial in p_bands.items())
    Q_tz = sum(polynomial * z**layer for layer, polynomial in q_bands.items())
    predicted = sp.expand(
        sum(
            z ** (p_layer + q_layer)
            * (
                p_layer * p_polynomial * sp.diff(q_polynomial, t)
                - q_layer * sp.diff(p_polynomial, t) * q_polynomial
            )
            for p_layer, p_polynomial in p_bands.items()
            for q_layer, q_polynomial in q_bands.items()
        )
    )
    direct = sp.simplify(
        (
            sp.diff(P_tz.subs({t: t_xy, z: z_xy}), X)
            * sp.diff(Q_tz.subs({t: t_xy, z: z_xy}), y)
            - sp.diff(P_tz.subs({t: t_xy, z: z_xy}), y)
            * sp.diff(Q_tz.subs({t: t_xy, z: z_xy}), X)
        ).subs({X: t * z, y: z**-1})
    )
    if sp.simplify(direct - predicted) != 0:
        raise AssertionError("the corrected Laurent recurrence failed")
    return {
        "chart": {"t": "X*y", "z": "y^-1"},
        "coordinate_bracket": "[t,z]_(X,y)=-z",
        "band_recurrence": (
            "J_L(t)=sum_(i+j=L)(i*p_i(t)*q_j'(t)-"
            "j*p_i'(t)*q_j(t))"
        ),
        "bracket_term": (
            "[p_i(t)z^i,q_j(t)z^j]="
            "(i*p_i*q_j'-j*p_i'*q_j)z^(i+j)"
        ),
        "formal_top_layer": 40,
        "zero_layers": [39, 5],
        "zero_layer_count": 35,
        "target_layer": 4,
        "target": "X^4=t^4*z^4",
    }


def band_record(band: Band) -> dict[str, object]:
    return {
        "layer": band.layer,
        "source_i_range": [band.source_i_min, band.source_i_max],
        "source_j": list(band.source_j),
        "source_count": band.source_count,
        "translation_basis": "t^layer*(1+t)^source_j",
        "jet_vanishing_order": band.jet_vanishing_order,
        "dimension": band.dimension,
        "active_t_exponent_interval": [
            band.t_exponent_min,
            band.t_exponent_max,
        ],
        "active_kummer_characters": sorted(
            {exponent % 5 for exponent in band.active_t_exponents}
        ),
    }


def digest_linear_constraints(
    bands: dict[int, Band],
    normalizations: dict[tuple[str, int, int], str],
) -> tuple[str, int, int]:
    digest = hashlib.sha256()
    jet_equation_count = 0
    linear_term_count = 0
    for layer, band in sorted(bands.items()):
        for order in range(band.jet_vanishing_order):
            terms = []
            for source_i, source_j in zip(
                range(band.source_i_min, band.source_i_max + 1),
                band.source_j,
            ):
                coefficient = comb(source_j, order)
                if coefficient:
                    terms.append((band.variable(source_i), coefficient))
            digest.update(
                (
                    f"jet:{band.side}:{layer}:{order}:"
                    + ",".join(f"{variable}:{coefficient}" for variable, coefficient in terms)
                    + "=0\n"
                ).encode()
            )
            jet_equation_count += 1
            linear_term_count += len(terms)
    for (side, layer, exponent), value in sorted(normalizations.items()):
        band = bands[layer]
        order = exponent - layer
        terms = []
        for source_i, source_j in zip(
            range(band.source_i_min, band.source_i_max + 1),
            band.source_j,
        ):
            coefficient = comb(source_j, order)
            if coefficient:
                terms.append((band.variable(source_i), coefficient))
        digest.update(
            (
                f"normalization:{side}:{layer}:{exponent}:"
                + ",".join(f"{variable}:{coefficient}" for variable, coefficient in terms)
                + f"={value}\n"
            ).encode()
        )
        linear_term_count += len(terms)
    return digest.hexdigest(), jet_equation_count, linear_term_count


def count_residue_terms(base: int, degree: int) -> dict[str, int]:
    result = {str(character): 0 for character in range(5)}
    for offset in range(degree + 1):
        result[str((base + offset) % 5)] += 1
    return result


def classify_layer(
    layer: int,
    p_bands: dict[int, Band],
    q_bands: dict[int, Band],
    digest: "hashlib._Hash",
) -> dict[str, object]:
    pair_rows: list[dict[str, object]] = []
    equation_exponents: set[int] = set()
    post_jet_exponents: set[int] = set()
    character_term_counts = {str(character): 0 for character in range(5)}
    raw_variable_pairs = 0
    active_generators = 0
    for p_layer in range(layer - 25, 16):
        q_layer = layer - p_layer
        p_band = p_bands[p_layer]
        q_band = q_bands[q_layer]
        for p_exponent in p_band.active_t_exponents:
            for q_exponent in q_band.active_t_exponents:
                if p_layer * q_exponent - q_layer * p_exponent:
                    post_jet_exponents.add(p_exponent + q_exponent - 1)
        pair_raw = p_band.source_count * q_band.source_count
        pair_active = 0
        pair_degree_min: int | None = None
        pair_degree_max: int | None = None
        raw_variable_pairs += pair_raw
        for p_i, p_power in zip(
            range(p_band.source_i_min, p_band.source_i_max + 1),
            p_band.source_j,
        ):
            for q_i, q_power in zip(
                range(q_band.source_i_min, q_band.source_i_max + 1),
                q_band.source_j,
            ):
                factor = p_layer * q_power - q_layer * p_power
                if factor == 0:
                    continue
                generator_degree = p_power + q_power - 1
                if generator_degree < 0:
                    raise AssertionError("active generator acquired negative degree")
                digest.update(
                    (
                        f"{layer}:{p_layer}:{q_layer}:{p_i}:{q_i}:"
                        f"{factor}:{generator_degree}\n"
                    ).encode()
                )
                pair_active += 1
                active_generators += 1
                pair_degree_min = (
                    generator_degree
                    if pair_degree_min is None
                    else min(pair_degree_min, generator_degree)
                )
                pair_degree_max = (
                    generator_degree
                    if pair_degree_max is None
                    else max(pair_degree_max, generator_degree)
                )
                equation_exponents.update(
                    range(layer, layer + generator_degree + 1)
                )
                residue_counts = count_residue_terms(layer, generator_degree)
                for character, count in residue_counts.items():
                    character_term_counts[character] += count
        pair_rows.append(
            {
                "P_layer": p_layer,
                "Q_layer": q_layer,
                "raw_source_variable_pairs": pair_raw,
                "active_binomial_generators": pair_active,
                "generator_degree_interval": (
                    None
                    if pair_degree_min is None
                    else [pair_degree_min, pair_degree_max]
                ),
            }
        )
    equations_by_character = {
        str(character): sum(
            exponent % 5 == character
            for exponent in equation_exponents
        )
        for character in range(5)
    }
    post_jet_by_character = {
        str(character): sum(
            exponent % 5 == character
            for exponent in post_jet_exponents
        )
        for character in range(5)
    }
    return {
        "layer": layer,
        "rhs": "t^4" if layer == 4 else "0",
        "band_pair_count": len(pair_rows),
        "band_pairs": pair_rows,
        "raw_source_variable_pairs": raw_variable_pairs,
        "active_binomial_generators": active_generators,
        "raw_scalar_coefficient_row_count": len(equation_exponents),
        "raw_scalar_t_exponent_interval": (
            None
            if not equation_exponents
            else [min(equation_exponents), max(equation_exponents)]
        ),
        "raw_scalar_rows_by_kummer_character": equations_by_character,
        "post_jet_structurally_active_coefficient_rows": len(post_jet_exponents),
        "post_jet_t_exponent_interval": (
            None
            if not post_jet_exponents
            else [min(post_jet_exponents), max(post_jet_exponents)]
        ),
        "post_jet_rows_by_kummer_character": post_jet_by_character,
        "expanded_quadratic_terms_by_kummer_character": character_term_counts,
    }


def common_power_regression() -> dict[str, object]:
    t, u, v = sp.symbols("t u v")
    H = sp.Function("H")(t)
    p_top = t**21 * H**3
    q_top = -sp.Rational(9, 5) * t**35 * H**5
    top = sp.simplify(
        15 * p_top * sp.diff(q_top, t)
        - 25 * sp.diff(p_top, t) * q_top
    )
    if top != 0:
        raise AssertionError("the normalized common-power bracket is nonzero")

    p_terminal = t**4
    q_terminal = -t
    target = sp.expand(
        3 * p_terminal * sp.diff(q_terminal, t)
        - sp.diff(p_terminal, t) * q_terminal
    )
    if target != t**4:
        raise AssertionError("the terminal Laurent bracket is not t^4")

    # The source-band exponents are multiples of five after u=1+t.  Put
    # A(u)=(u-1)^2 H(u-1).  The P and Q top-band conditions say that A^3
    # and A^5 are invariant under u -> zeta*u.  Bezout (2*3-5=1) makes A
    # invariant too, hence A=G(u^5).  Divisibility by (u-1)^2 makes
    # G(v) divisible by (v-1)^2.
    a, b = sp.symbols("a b")
    c = sp.Rational(1, 25) - a - b
    R = a * v**2 + b * v + c
    geometric_sum = sum(u**power for power in range(5))
    H_parameterized_u = sp.expand(geometric_sum**2 * R.subs(v, u**5))
    H_parameterized_t = sp.expand(H_parameterized_u.subs(u, 1 + t))
    if sp.expand(H_parameterized_t.subs(t, 0)) != 1:
        raise AssertionError("the common-power H normalization failed")
    if sp.Poly(H_parameterized_t, t).degree() != 18:
        raise AssertionError("the parameterized H lost degree eighteen")
    A_parameterized = sp.expand(
        (u - 1) ** 2 * H_parameterized_t.subs(t, u - 1)
    )
    expected_A = sp.expand((u**5 - 1) ** 2 * R.subs(v, u**5))
    if sp.expand(A_parameterized - expected_A) != 0:
        raise AssertionError("the invariant common root parameterization failed")
    for power in (3, 5):
        exponents = [
            exponent[0]
            for exponent, coefficient in sp.Poly(
                sp.expand(A_parameterized**power), u
            ).terms()
            if coefficient != 0
        ]
        if any(exponent % 5 for exponent in exponents):
            raise AssertionError("a common-power source band left k[u^5]")
    discriminant = sp.expand(b**2 - 4 * a * c)

    # The negative-tail seed has an exact polynomial source lift.  Before the
    # Puiseux translation it is x*(x*y^5-1)^2*R(x*y^5).  This observation is
    # why deleting its negative powers in the later monomial chart is not a
    # legitimate source operation.
    X, y_new, x_source, y_source = sp.symbols(
        "X y_new x_source y_source"
    )
    source_w = x_source * y_source**5
    source_root = sp.expand(
        x_source
        * (source_w - 1) ** 2
        * R.subs(v, source_w)
    )
    translated_source_root = sp.expand(
        source_root.subs(
            {x_source: X**5, y_source: y_new + X**-1}, simultaneous=True
        )
    )
    translated_t = X * y_new
    expected_translated_root = sp.expand(
        X**5
        * ((1 + translated_t) ** 5 - 1) ** 2
        * R.subs(v, (1 + translated_t) ** 5)
    )
    if sp.expand(translated_source_root - expected_translated_root) != 0:
        raise AssertionError("the polynomial source lift of C_top failed")
    return {
        "root_band": "C=t^7*H(t)*z^5",
        "H_degree": 18,
        "H_normalization": "H(0)=1; leading coefficient nonzero",
        "P_top_band": "t^21*H(t)^3*z^15",
        "Q_top_band": "-9/5*t^35*H(t)^5*z^25",
        "top_layer": 40,
        "top_bracket": "0",
        "terminal_pair": {
            "P": "t^4*z^3",
            "Q": "-t*z",
            "layer": 4,
            "bracket": "t^4*z^4",
        },
        "joint_top_band_classification": {
            "u_coordinate": "u=1+t",
            "invariant_polynomial": "A(u)=(u-1)^2*H(u-1)",
            "source_conditions": ["A(u)^3 in k[u^5]", "A(u)^5 in k[u^5]"],
            "bezout_consequence": (
                "A(u)^6/A(u)^5=A(u), so invariance of A^3 and A^5 "
                "forces A in k[u^5]"
            ),
            "parameterization": (
                "H(t)=(1+u+u^2+u^3+u^4)^2*R(u^5), u=1+t"
            ),
            "R": "a*v^2+b*v+(1/25-a-b)",
            "conditions": ["a != 0", "R(1)=1/25"],
            "dimension": 2,
            "natural_root_strata": {
                "distinct_roots": f"{discriminant} != 0",
                "double_root": f"{discriminant} = 0",
                "zero_root": "1/25-a-b=0",
                "claim_boundary": (
                    "these are algebraic root strata, not a proved gamma list"
                ),
            },
        },
        "exact_polynomial_source_lift": {
            "before_translation": "C_R=x*(x*y^5-1)^2*R(x*y^5)",
            "source_total_degree": 25,
            "source_terminal_height": 1,
            "source_band": 5,
            "after_translation": (
                "C_R=t^5*z^5*((1+t)^5-1)^2*R((1+t)^5)"
            ),
            "consequence": (
                "the full negative Laurent tail is the monomial-chart image "
                "of an actual polynomial source band, not a disposable error"
            ),
        },
    }


def band_factor_data(band: Band) -> tuple[int, int, int]:
    """Return ``(u_power, vanishing_power, free_w_degree)`` for a band.

    With ``u=1+t`` and ``w=u^5``, the exact jet-reduced band is

        t^ell*u^u_power*(w-1)^vanishing_power*k[w]_(<=free_w_degree).

    This is just the binomial basis in factored form: all source exponents
    are congruent modulo five, and ``w-1`` is a uniformizer at ``t=0``.
    """

    u_power = band.source_j[0]
    w_degree = (band.source_j[-1] - u_power) // 5
    free_w_degree = w_degree - band.jet_vanishing_order
    if free_w_degree + 1 != band.dimension:
        raise AssertionError("the factored source-band dimension changed")
    return u_power, band.jet_vanishing_order, free_w_degree


def upper_descent_classification(
    p_bands: dict[int, Band], q_bands: dict[int, Band]
) -> dict[str, object]:
    """Classify the *linearized* source-band operator below the top.

    An earlier version inserted ``p=C0^2*U`` and ``q=C0^4*V`` before proving
    those divisibilities.  Exact source bands do not satisfy them.  The
    unrestricted tangent kernel is larger: every P-band deformation is
    followed by a Q-band deformation obtained by multiplication by ``C0^2``.
    """

    A, A_prime, p, p_prime = sp.symbols("A A_prime p p_prime")
    descent = sp.symbols("descent")
    q_scale = -sp.Rational(9, 5)
    p_top = A**3
    p_top_prime = 3 * A**2 * A_prime
    q_top = q_scale * A**5
    q_top_prime = 5 * q_scale * A**4 * A_prime

    # For arbitrary p, not merely p divisible by A^2, this is the tangent to
    # Q=(-9/5)*P^(5/3).  Its band is ten below the P band.
    q_follow = -3 * A**2 * p
    q_follow_prime = -3 * (2 * A * A_prime * p + A**2 * p_prime)
    tangent_bracket = sp.expand(
        (15 - descent) * p * q_top_prime
        - 25 * p_prime * q_top
        + 15 * p_top * q_follow_prime
        - (25 - descent) * p_top_prime * q_follow
    )
    if tangent_bracket != 0:
        raise AssertionError("the unrestricted top tangent identity failed")

    # Retain the old ODE only as an explicitly restricted-slice regression.
    U, U_prime, V, V_prime = sp.symbols("U U_prime V V_prime")
    restricted_p = A**2 * U
    restricted_p_prime = 2 * A * A_prime * U + A**2 * U_prime
    restricted_q = q_scale * A**4 * V
    restricted_q_prime = q_scale * (
        4 * A**3 * A_prime * V + A**4 * V_prime
    )
    restricted_bracket = sp.expand(
        (15 - descent) * restricted_p * q_top_prime
        - 25 * restricted_p_prime * q_top
        + 15 * p_top * restricted_q_prime
        - (25 - descent) * p_top_prime * restricted_q
    )
    W = 5 * U - 3 * V
    W_prime = 5 * U_prime - 3 * V_prime
    restricted_reduction = sp.expand(
        q_scale
        * A**6
        * ((5 - descent) * A_prime * W - 5 * A * W_prime)
    )
    if sp.expand(restricted_bracket - restricted_reduction) != 0:
        raise AssertionError("the restricted upper-descent ODE changed")

    rows: list[dict[str, object]] = []
    source_centralizer_powers: list[int] = []
    formal_negative_powers: list[int] = []
    for descent_value in range(1, 36):
        p_band = p_bands[15 - descent_value]
        q_band = q_bands[25 - descent_value]
        p_u, p_vanishing, p_degree = band_factor_data(p_band)
        q_u, q_vanishing, q_degree = band_factor_data(q_band)

        # C0=t^5*(w-1)^2*R(w), deg R=2.  Therefore C0^2*p lies
        # in the exact Q band.  The following quotient data prove the
        # inclusion directly in the factored source bases.
        if (p_u - q_u) % 5:
            raise AssertionError("P/Q tangent bands lost their Kummer match")
        quotient_w_power = (p_u - q_u) // 5
        quotient_vanishing = p_vanishing + 4 - q_vanishing
        image_degree = (
            quotient_w_power + quotient_vanishing + 4 + p_degree
        )
        if min(quotient_w_power, quotient_vanishing) < 0:
            raise AssertionError("C0^2 no longer maps the P band into Q")
        if image_degree > q_degree:
            raise AssertionError("C0^2 exceeded the Q source-degree bound")

        # The q-only equation is
        #   5*C0*q'-(25-delta)*C0'*q=0.
        # Hence q^5 is a scalar multiple of C0^(25-delta).  The fixed
        # (w-1)^2 factor of C0 proves that a source Laurent polynomial exists
        # exactly for delta=5*j with 1<=j<=5, giving C0^(5-j).
        centralizer_power: int | None = None
        if descent_value % 5 == 0:
            candidate_power = 5 - descent_value // 5
            if candidate_power >= 0:
                centralizer_power = candidate_power
                source_centralizer_powers.append(candidate_power)
            else:
                formal_negative_powers.append(candidate_power)
        centralizer_dimension = int(centralizer_power is not None)
        kernel_dimension = p_band.dimension + centralizer_dimension
        rows.append(
            {
                "descent": descent_value,
                "layer": 40 - descent_value,
                "P_band": 15 - descent_value,
                "Q_band": 25 - descent_value,
                "P_band_dimension": p_band.dimension,
                "Q_band_dimension": q_band.dimension,
                "operator_rank": q_band.dimension - centralizer_dimension,
                "kernel_dimension": kernel_dimension,
                "arbitrary_P_follow_mode_dimension": p_band.dimension,
                "Q_follow": "q=(-3)*C0^2*p",
                "extra_source_centralizer": (
                    None
                    if centralizer_power is None
                    else f"C0^{centralizer_power}"
                ),
                "formal_negative_centralizer": (
                    f"C0^{5-descent_value//5}"
                    if descent_value in (30, 35)
                    else None
                ),
            }
        )

    if source_centralizer_powers != [4, 3, 2, 1, 0]:
        raise AssertionError("the nonnegative centralizer list changed")
    if formal_negative_powers != [-1, -2]:
        raise AssertionError("the formal negative resonance list changed")
    if [rows[index - 1]["kernel_dimension"] for index in range(1, 6)] != [
        6,
        6,
        7,
        7,
        10,
    ]:
        raise AssertionError("the corrected first-five kernel profile changed")

    return {
        "correction": (
            "the former common-root continuation used the unproved "
            "divisibilities p%C0^2=0 and q%C0^4=0; it is only a restricted "
            "slice, not the exact source-band kernel"
        ),
        "exact_source_band_factorization": (
            "V_ell=t^ell*u^j0*(u^5-1)^nu*QQ[u^5]_(<=N-nu)"
        ),
        "unrestricted_tangent_identity": {
            "P_top": "C0^3",
            "Q_top": "(-9/5)*C0^5",
            "arbitrary_P_band": "p",
            "forced_Q_follow": "q=(-3)*C0^2*p",
            "multiplication_closure": (
                "C0^2 has source degree 50, height 2, and band 10; "
                "it maps every exact P band into the paired Q band"
            ),
        },
        "q_only_kernel_theorem": {
            "equation": "5*C0*q'-(25-delta)*C0'*q=0",
            "fifth_power_identity": "q^5=constant*C0^(25-delta)",
            "source_kernel_condition": (
                "delta=5*j with 1<=j<=5; generator C0^(5-j)"
            ),
            "fixed_factor_certificate": (
                "C0 contains (u^5-1)^2, so fifth-power valuations force "
                "5|delta; negative powers cannot be Laurent polynomials"
            ),
        },
        "source_centralizer_modes": ["C0^4", "C0^3", "C0^2", "C0", "1"],
        "formal_but_not_source_kernel_modes": [
            {"descent": 30, "layer": 10, "mode": "C0^-1", "role": "lambda"},
            {"descent": 35, "layer": 5, "mode": "C0^-2"},
        ],
        "r3_rank_profile": rows,
        "first_five_exact_kernel_dimensions": [6, 6, 7, 7, 10],
        "restricted_divisible_slice": {
            "assumption": "p=C0^2*U and q=(-9/5)*C0^4*V",
            "W": "5*U-3*V",
            "equation": "(5-delta)*C0'*W-5*C0*W'=0",
            "warning": "the exact source-band equations do not force this slice",
        },
        "fitting_handoff": (
            "at nonlinear descent delta, let T_delta(q)="
            "5*C0*q'-(25-delta)*C0'*q on the exact Q band; the new P "
            "columns lie in im(T_delta), so solvability of the known forcing "
            "is rank([T_delta|forcing])=rank(T_delta), equivalently the "
            "corresponding maximal-minor/Fitting condition"
        ),
        "outcome": (
            "layers 39 through 36 do not force source-root continuation; "
            "the first five tangent kernels have dimensions 6,6,7,7,10. "
            "The layer-35 extra mode is the commuting C0^4 term, while the "
            "formal lambda*C0^-1 resonance is at layer 10 and is not an "
            "independent source-band kernel"
        ),
        "remaining_unclassified_zero_layers": [39, 5],
        "remaining_unclassified_zero_layer_count": 35,
    }


def family_tangent_theorem() -> dict[str, object]:
    """Record the corrected tangent and resonance formulas for every F2 r."""

    r, delta = sp.symbols("r delta", integer=True, positive=True)
    m = 2 * r - 1
    p_layer = 5 * r - delta
    q_layer = 5 * m - delta
    if sp.expand(q_layer - p_layer - 5 * (m - r)) != 0:
        raise AssertionError("the family follow-band identity changed")
    if sp.expand(r + (m - r) - m) != 0:
        raise AssertionError("the family source-degree identity changed")
    if sp.expand(5 * (r + m) - 10 * r - 5 * (r - 1)) != 0:
        raise AssertionError("the lambda-layer formula changed")

    samples: list[dict[str, object]] = []
    for r_value in range(2, 9):
        m_value = 2 * r_value - 1
        top_layer = 5 * (3 * r_value - 1)
        maximum_zero_descent = 5 * (3 * r_value - 2)
        source_powers = list(range(m_value - 1, -1, -1))
        negative_powers = list(range(-1, -r_value, -1))
        if len(source_powers) != m_value:
            raise AssertionError("the source centralizer count changed")
        if len(negative_powers) != r_value - 1:
            raise AssertionError("the negative resonance count changed")
        if top_layer - 10 * r_value != 5 * (r_value - 1):
            raise AssertionError("a sampled lambda layer changed")
        if top_layer - maximum_zero_descent != 5:
            raise AssertionError("the sampled zero window changed")
        samples.append(
            {
                "r": r_value,
                "m": m_value,
                "top_layer": top_layer,
                "zero_descent_interval": [1, maximum_zero_descent],
                "source_centralizer_powers": source_powers,
                "formal_negative_powers": negative_powers,
                "lambda_descent": 10 * r_value,
                "lambda_layer": 5 * (r_value - 1),
                "target_descent": top_layer - 4,
            }
        )

    return {
        "parameters": "r>=2, m=2*r-1",
        "top_layers": {
            "P": "5*r",
            "Q": "5*m",
            "bracket": "5*(r+m)=5*(3*r-1)",
        },
        "arbitrary_follow_mode": {
            "P_band": "5*r-delta",
            "Q_band": "5*m-delta",
            "formula": "q=-B_r*(m/r)*C0^(m-r)*p",
            "source_closure": (
                "m-r=r-1, so multiplication adds source degree 25*(r-1), "
                "height r-1, and band 5*(r-1)"
            ),
        },
        "q_only_equation": (
            "5*C0*q'-(5*m-delta)*C0'*q=0; equivalently "
            "q^5=constant*C0^(5*m-delta)"
        ),
        "source_centralizer_resonances": (
            "delta=5*j, 1<=j<=m, with generator C0^(m-j)"
        ),
        "formal_negative_resonances": (
            "delta=5*(m+k), 1<=k<=r-1, with formal mode C0^(-k); "
            "none is a homogeneous Laurent-polynomial source kernel"
        ),
        "lambda_location": {
            "mode": "C0^-1",
            "descent": "5*(m+1)=10*r",
            "layer": "5*(r-1)",
            "interpretation": (
                "lambda can only arise through nonlinear forced cancellation "
                "with the F tail, not as an independent exact source-band mode"
            ),
        },
        "target": {
            "layer": 4,
            "descent": "5*(3*r-1)-4=15*r-9",
        },
        "samples": samples,
    }


def nonlinear_first_defect_audit() -> dict[str, object]:
    """Use the next four rows to classify the first non-root r=3 defect."""

    h, C = sp.symbols("h C", nonzero=True)
    p1, p2, p3, p4, p5 = sp.symbols("p1 p2 p3 p4 p5")
    normalized_tail = sum(
        coefficient * h**index / C**3
        for index, coefficient in enumerate((p1, p2, p3, p4, p5), 1)
    )
    formal_q = sp.series(
        -sp.Rational(9, 5)
        * C**5
        * (1 + normalized_tail) ** sp.Rational(5, 3),
        h,
        0,
        6,
    ).removeO()
    q_coefficients = [
        sp.factor(sp.expand(formal_q).coeff(h, index))
        for index in range(1, 6)
    ]
    expected_first_four = [
        -3 * C**2 * p1,
        -3 * C**2 * p2 - p1**2 / C,
        -3 * C**2 * p3 - 2 * p1 * p2 / C + p1**3 / (9 * C**4),
        -3 * C**2 * p4
        - (2 * p1 * p3 + p2**2) / C
        + p1**2 * p2 / (3 * C**4)
        - p1**4 / (27 * C**7),
    ]
    if any(
        sp.cancel(actual - expected) != 0
        for actual, expected in zip(q_coefficients, expected_first_four)
    ):
        raise AssertionError("the first nonlinear fractional-power rows changed")

    # At a multiplicity-two prime pi of C, the only valuation below C^2 that
    # survives rows two and three is v_pi(p1)=3.  Scaling h=pi^3*s leaves the
    # quadratic residue 1+a*s+b*s^2.  Rows four and five have the following
    # primitive numerators.  They have no common solution with a nonzero.
    s, a, b = sp.symbols("s a b")
    local_series = sp.series(
        (1 + a * s + b * s**2) ** sp.Rational(5, 3), s, 0, 6
    ).removeO()
    coefficient_four = sp.factor(sp.expand(local_series).coeff(s, 4))
    coefficient_five = sp.factor(sp.expand(local_series).coeff(s, 5))
    numerator_four = a**4 - 9 * a**2 * b + 27 * b**2
    numerator_five = 7 * a**4 - 60 * a**2 * b + 135 * b**2
    if sp.expand(coefficient_four - sp.Rational(5, 243) * numerator_four) != 0:
        raise AssertionError("the double-prime fourth residue changed")
    if sp.expand(
        coefficient_five + sp.Rational(1, 729) * a * numerator_five
    ) != 0:
        raise AssertionError("the double-prime fifth residue changed")
    resultant = sp.factor(sp.resultant(numerator_four, numerator_five, b))
    if resultant != 1701 * a**8:
        raise AssertionError("the double-prime residue resultant changed")

    # The target descent is 36.  If the first non-root defect has spacing k,
    # all rows 2k,...,5k occur before the target for k<=7, so the same local
    # argument forces it back into C^2.  At k=8, row 4k=32 is still zero but
    # the target intervenes before row 5k=40; only numerator_four is forced.
    target_descent = 36
    killed_first_defects = [
        k for k in range(1, 8) if 5 * k < target_descent
    ]
    if killed_first_defects != list(range(1, 8)):
        raise AssertionError("the pre-target first-defect range changed")
    # A nonnegative centralizer C0^(5-j) starts at absolute descent 5*j.
    # On the exceptional double-prime scaling its leading valuation exceeds
    # the P^(5/3) residue by j*(15/k-2).  This is positive for k<=7, so none
    # of the commuting modes can cancel the two primitive residues above.
    if [15 - 2 * k > 0 for k in range(1, 8)] != [True] * 7:
        raise AssertionError("a pre-target centralizer reached the leading residue")
    if 15 - 2 * 8 >= 0:
        raise AssertionError("the descent-eight threshold changed")
    candidate_ratio = sp.Symbol("candidate_ratio")
    candidate_equation = sp.expand(
        numerator_four.subs({a: 1, b: candidate_ratio})
    )
    if candidate_equation != 27 * candidate_ratio**2 - 9 * candidate_ratio + 1:
        raise AssertionError("the first residual ratio equation changed")
    if sp.discriminant(candidate_equation, candidate_ratio) != -27:
        raise AssertionError("the first residual discriminant changed")

    return {
        "formal_particular_solution": "Q=(-9/5)*P^(5/3)",
        "first_rows": {
            "q1": str(q_coefficients[0]),
            "q2": str(q_coefficients[1]),
            "q3": str(q_coefficients[2]),
            "q4": str(q_coefficients[3]),
        },
        "divisibility_rows": [
            "C0 divides p1^2",
            "C0^4 divides p1*(p1^2-18*C0^3*p2)",
        ],
        "simple_prime_consequence": (
            "at every multiplicity-one prime of C0, rows 2 and 3 force "
            "v(p1)>=2"
        ),
        "double_prime_consequence": {
            "rows_2_and_3": (
                "either v(p1)>=4 or the exceptional valuation is "
                "v(C0)=2, v(p1)=3, v(p2)=0"
            ),
            "row_4_numerator": str(numerator_four),
            "row_5_numerator_after_removing_nonzero_a": str(numerator_five),
            "resultant_in_b": str(resultant),
            "consequence": (
                "when rows four and five are both zero, the exceptional "
                "double-prime valuation is impossible"
            ),
        },
        "nonlinear_root_recovery": (
            "if the first non-root P defect occurs at descent k<=7, all five "
            "required multiples precede the target descent 36, so it is "
            "forced back into the C0^2 root slice"
        ),
        "centralizer_separation": (
            "a C0^(5-j) mode has exceptional-prime valuation gap "
            "j*(15/k-2) over the primitive P^(5/3) residue; this is positive "
            "for k<=7 and changes sign first at k=8"
        ),
        "forced_source_root_through_descent": 7,
        "first_exact_residual_candidate": {
            "descent": 8,
            "P_band": 7,
            "reason": (
                "rows 16,24,32 precede the target, but target descent 36 "
                "intervenes before the fifth-multiple row 40"
            ),
            "support": "only multiplicity-two primes of C0",
            "local_valuations": "v(C0)=2, v(p8)=3, v(p16)=0",
            "normalized_ratio": "y=b/a^2",
            "residual_equation": "27*y^2-9*y+1=0",
            "discriminant": -27,
            "field": "QQ(sqrt(-3))",
            "claim_boundary": (
                "this is an exact local residual seed, not a reconstructed "
                "global source pair or plane counterexample"
            ),
        },
    }


def build_payload() -> dict[str, object]:
    coordinate = coordinate_bracket_regression()
    # Layers 4..40 require P bands -21..15 and Q bands -11..25.
    p_bands = {
        layer: make_band("P", 75, 3, layer)
        for layer in range(-21, 16)
    }
    q_bands = {
        layer: make_band("Q", 125, 5, layer)
        for layer in range(-11, 26)
    }
    p_normalizations = {
        ("P", 3, 4): "1",
        ("P", 15, 21): "1",
    }
    q_normalizations = {
        ("Q", 1, 1): "-1",
        ("Q", 13, 18): "-3",
        ("Q", 25, 35): "-9/5",
    }
    p_linear = digest_linear_constraints(p_bands, p_normalizations)
    q_linear = digest_linear_constraints(q_bands, q_normalizations)

    quadratic_digest = hashlib.sha256()
    layers = [
        classify_layer(layer, p_bands, q_bands, quadratic_digest)
        for layer in range(40, 3, -1)
    ]
    missing_layers = [
        record
        for record in layers
        if 5 <= int(record["layer"]) <= 39
    ]
    if len(missing_layers) != 35:
        raise AssertionError("the missing layer count changed")
    missing_p_layers = range(-20, 16)
    missing_q_layers = range(-10, 26)
    missing_parameter_count = sum(
        p_bands[layer].dimension for layer in missing_p_layers
    ) + sum(q_bands[layer].dimension for layer in missing_q_layers)
    if missing_parameter_count != 978:
        raise AssertionError("the B0 parameter count changed")

    return {
        "schema": "plane-jc.f2-75-125-character-layers.v2",
        "status": (
            "exact-B0-envelope-recurrence-and-corrected-top-tangent;"
            "not-exhaustive-normal-form"
        ),
        "coordinate_and_recurrence": coordinate,
        "support_envelope_theorem": {
            "source_degree_bounds": {"P": 75, "Q": 125},
            "source_band_invariant": "layer=5*i-j",
            "translated_descendants": (
                "(X_power,y_power)=(5*i-k,j-k), 0<=k<=j"
            ),
            "terminal_halfspaces": {
                "P": "5*X_power-17*y_power<=3",
                "Q": "5*X_power-17*y_power<=5",
            },
            "band_translation": (
                "x^i*y^j -> t^(5*i-j)*(t+lambda)^j*z^(5*i-j); "
                "lambda is scaled to one"
            ),
            "jet_rank_proof": (
                "the first r coefficient rows binomial(j,k) form a "
                "falling-factorial Vandermonde matrix of rank r"
            ),
        },
        "bands_for_layers_4_through_40": {
            "P": [band_record(p_bands[layer]) for layer in sorted(p_bands)],
            "Q": [band_record(q_bands[layer]) for layer in sorted(q_bands)],
            "source_variable_count": {
                "P": sum(band.source_count for band in p_bands.values()),
                "Q": sum(band.source_count for band in q_bands.values()),
            },
            "linear_dimension_after_terminal_jets": {
                "P": sum(band.dimension for band in p_bands.values()),
                "Q": sum(band.dimension for band in q_bands.values()),
            },
            "missing_35_layer_window_dimension": missing_parameter_count,
            "missing_window_dimension_after_five_normalizations": (
                missing_parameter_count - 5
            ),
            "missing_window_dimension_after_joint_top_band_and_other_three_normalizations": (
                missing_parameter_count - 18 + 2 - 3
            ),
        },
        "linear_constraints": {
            "P": {
                "digest_sha256": p_linear[0],
                "jet_equation_count": p_linear[1],
                "term_count_including_normalizations": p_linear[2],
            },
            "Q": {
                "digest_sha256": q_linear[0],
                "jet_equation_count": q_linear[1],
                "term_count_including_normalizations": q_linear[2],
            },
            "normalizations": {
                "P_3_t4": "1",
                "P_15_t21": "1",
                "Q_1_t1": "-1",
                "Q_13_t18": "-3",
                "Q_25_t35": "-9/5",
            },
        },
        "common_power_and_terminal_regression": common_power_regression(),
        "family_tangent_theorem": family_tangent_theorem(),
        "upper_descent_classification": upper_descent_classification(
            p_bands, q_bands
        ),
        "nonlinear_first_defect_audit": nonlinear_first_defect_audit(),
        "compressed_quadratic_system": {
            "generator_formula": (
                "(ell*s-m*r)*alpha_(ell,i)*beta_(m,k)*"
                "t^(ell+m)*(1+t)^(r+s-1)"
            ),
            "generator_digest_sha256": quadratic_digest.hexdigest(),
            "layers_40_through_4": layers,
            "missing_zero_layers": [39, 5],
            "missing_zero_layer_count": 35,
            "missing_window_totals": {
                "band_pairs": sum(
                    int(record["band_pair_count"]) for record in missing_layers
                ),
                "active_binomial_generators": sum(
                    int(record["active_binomial_generators"])
                    for record in missing_layers
                ),
                "raw_scalar_coefficient_rows": sum(
                    int(record["raw_scalar_coefficient_row_count"])
                    for record in missing_layers
                ),
                "post_jet_structurally_active_coefficient_rows": sum(
                    int(record["post_jet_structurally_active_coefficient_rows"])
                    for record in missing_layers
                ),
            },
        },
        "classification_boundary": {
            "proved": [
                "the corrected monomial-Jacobian Laurent recurrence",
                "the exact finite support envelope forced by total degree and "
                "the two known edge halfspaces",
                "all 35 zero-layer incidence lists and their five Kummer "
                "character splittings",
                "a deterministic compressed representation of every resulting "
                "quadratic coefficient equation",
                "the joint common-power top-band parameterization by one "
                "normalized quadratic R",
                "the corrected exact source-band tangent profile through all "
                "35 zero layers, including the unrestricted P-follow modes",
                "the proof that the formal C0^-1 and C0^-2 resonances are not "
                "independent homogeneous source-band kernels",
                "the nonlinear first-defect resultant, which forces source-root "
                "continuation through descent 7 and isolates the first exact "
                "local residual at descent 8",
            ],
            "not_proved": [
                "that the B0 envelope is the exhaustive lower Newton support",
                "the gamma-branch list for multiplicities (3,5)",
                "the nonlinear forced Fitting equations below the common top",
                "inconsistency of the resulting system or exclusion of F2",
            ],
            "next_elimination": (
                "retain every exact P-band variable and compile the known "
                "nonlinear forcing against coker(T_delta), beginning with the "
                "descent-8 double-prime ratio 27*y^2-9*y+1=0; the first "
                "actual lambda location is descent 30 (layer 10)"
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
                "pinned F2 layer artifact is stale; inspect before --refresh"
            )
    totals = payload["compressed_quadratic_system"]["missing_window_totals"]
    print("PASS: [t,z]_(X,y)=-z and the formal top Laurent layer is 40")
    print("PASS: exactly 35 zero layers, 39 through 5, are character-classified")
    print(
        "PASS: exact B0 window:",
        "978 jet-reduced parameters (973 normalized),",
        f"{totals['post_jet_structurally_active_coefficient_rows']} "
        "structurally active Keller rows",
    )
    print("PASS: the joint common-power top band is a two-parameter quadratic-R family")
    print(
        "PASS: corrected first-five tangent kernels have dimensions "
        "6,6,7,7,10 (the former 2,2,3,3,6 slice assumed divisibility)"
    )
    print("PASS: lambda*C0^-1 is at layer 10 and is not a source-band kernel")
    print(
        "PASS: nonlinear rows force root continuation through descent 7; "
        "the first local residual is 27*y^2-9*y+1 at descent 8"
    )
    print("PASS: the artifact remains explicitly short of an F2 exclusion")


if __name__ == "__main__":
    main()
