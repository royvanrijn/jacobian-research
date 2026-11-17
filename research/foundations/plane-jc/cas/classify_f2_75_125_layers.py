#!/usr/bin/env python3
"""Classify the complete F2 (75,125) Laurent B0 system.

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
choice.  The record includes the familiar 35 zero layers above the target and
all 204 zero layers below it, down to layer -200.  This is the finite B0
envelope plus the universal Keller recurrence; it is not an exhaustive Newton
polygon or gamma-branch theorem.  On the earliest movable-double-root branch,
the record also eliminates the complete fixed-endpoint Hermite block: one
coordinate is the terminal normalization and the other ten have determinant
75000, leaving thirteen global Hermite coordinates after substitution.
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
    p_layer_min = max(min(p_bands), layer - max(q_bands))
    p_layer_max = min(max(p_bands), layer - min(q_bands))
    for p_layer in range(p_layer_min, p_layer_max + 1):
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
    # Independent coefficient supports in the two jet-reduced bands give an
    # inexpensive structural upper bound, not the exact support of the
    # restricted bilinear map: common-power pairs can cancel identically.
    # Intersecting with the exact pre-jet generator rows removes rows which
    # are impossible even before restriction.
    post_jet_exponents.intersection_update(equation_exponents)
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
        "post_jet_support_row_upper_bound": len(post_jet_exponents),
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


def family_full_support_theorem() -> dict[str, object]:
    """Derive the complete B0 band and bracket ranges for every F2 member."""

    samples: list[dict[str, object]] = []
    for r_value in range(2, 9):
        m_value = 2 * r_value - 1
        slope_numerator = 7 * r_value - 4
        slope_denominator = 5 * r_value - 3
        side_rows: dict[str, object] = {}
        for side, exponent in (("P", r_value), ("Q", m_value)):
            degree = 25 * exponent
            height = exponent
            source_count = 0
            dimension = 0
            jet_count = 0
            for layer in range(-degree, 5 * exponent + 1):
                source_i_min = max(0, ceil_div(layer, 5))
                source_i_max = (degree + layer) // 6
                count = source_i_max - source_i_min + 1
                terminal_t_min = max(
                    layer,
                    ceil_div(slope_numerator * layer - height, slope_denominator),
                )
                jet_order = terminal_t_min - layer
                if not (0 <= jet_order <= count):
                    raise AssertionError("a family B0 band acquired invalid dimension")
                source_count += count
                jet_count += jet_order
                dimension += count - jet_order
            side_rows[side] = {
                "degree": degree,
                "band_interval": [-degree, 5 * exponent],
                "source_coefficients_before_terminal_jets": source_count,
                "terminal_jet_equations": jet_count,
                "dimension_after_terminal_jets": dimension,
            }
        bracket_top = 5 * (r_value + m_value)
        bracket_bottom = -25 * (r_value + m_value)
        if bracket_top != 5 * (3 * r_value - 1):
            raise AssertionError("the family full-support top changed")
        if bracket_bottom != -25 * (3 * r_value - 1):
            raise AssertionError("the family full-support bottom changed")
        post_target_zero_count = 3 - bracket_bottom + 1
        total_zero_count = bracket_top - bracket_bottom
        samples.append(
            {
                "r": r_value,
                "m": m_value,
                "P": side_rows["P"],
                "Q": side_rows["Q"],
                "bracket_layer_interval": [bracket_bottom, bracket_top],
                "target_layer": 4,
                "zero_layer_count": total_zero_count,
                "post_target_zero_layer_count": post_target_zero_count,
            }
        )

    r3 = samples[1]
    if r3["P"]["dimension_after_terminal_jets"] != 653:
        raise AssertionError("the family formula lost the r=3 P dimension")
    if r3["Q"]["dimension_after_terminal_jets"] != 1765:
        raise AssertionError("the family formula lost the r=3 Q dimension")
    if r3["post_target_zero_layer_count"] != 204:
        raise AssertionError("the r=3 post-target tail length changed")

    return {
        "parameters": "r>=2, m=2*r-1",
        "source_degrees": {"P": "25*r", "Q": "25*m"},
        "complete_band_intervals": {
            "P": "-25*r <= ell <= 5*r",
            "Q": "-25*m <= ell <= 5*m",
        },
        "band_source_index_range": (
            "max(0,ceil(ell/5)) <= i <= floor((25*n+ell)/6), "
            "j=5*i-ell, with n=r on P and n=m on Q"
        ),
        "terminal_jet_order": (
            "nu=max(ell,ceil(((7*r-4)*ell-h)/(5*r-3)))-ell, "
            "with h=r on P and h=m on Q"
        ),
        "exact_band_factorization": (
            "t^ell*u^j0*(u^5-1)^nu*K[u^5]_(<=N-nu)"
        ),
        "complete_bracket_layer_interval": [
            "-25*(3*r-1)",
            "5*(3*r-1)",
        ],
        "target_layer": 4,
        "total_zero_layer_count": "30*(3*r-1)",
        "post_target_zero_layer_count": "25*(3*r-1)+4",
        "warning": (
            "the 35 layers above the target are only the upper window; the "
            "corner-derived source envelope contains a genuine lower tail"
        ),
        "samples": samples,
    }


def nonlinear_first_defect_audit(
    p_bands: dict[int, Band],
    q_bands: dict[int, Band],
) -> dict[str, object]:
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
    multiplication_by_ratio = sp.Matrix(
        [[0, -sp.Rational(1, 27)], [1, sp.Rational(1, 3)]]
    )
    characteristic_ratio = sp.expand(
        multiplication_by_ratio.charpoly(candidate_ratio).as_expr()
    )
    if characteristic_ratio != candidate_equation / 27:
        raise AssertionError("the residual Artinian multiplication matrix changed")
    if multiplication_by_ratio.det() != sp.Rational(1, 27):
        raise AssertionError("the residual Artinian norm changed")

    # At descent 3*8=24 the candidate requires Q band 1 to vanish at its
    # supporting double prime.  That band is t*u^4*S(u^5), deg S<=20, and
    # the terminal normalization [t]q_1=-1 is exactly S(1)=-1.  Hence none
    # of the four nontrivial roots of u^5=1 can support the candidate.  A
    # multiplicity-two prime must instead come from a double root of R.
    q_one_u, q_one_vanishing, q_one_degree = band_factor_data(q_bands[1])
    if (q_one_u, q_one_vanishing, q_one_degree) != (4, 0, 20):
        raise AssertionError("the exact Q band-one factorization changed")
    local_t, S_at_one, S_slope = sp.symbols("local_t S_at_one S_slope")
    local_w = (1 + local_t) ** 5
    q_one_series = sp.expand(
        local_t
        * (1 + local_t) ** 4
        * (S_at_one + S_slope * (local_w - 1))
    )
    if q_one_series.coeff(local_t, 1) != S_at_one:
        raise AssertionError("the Q band-one normalization regression failed")

    # The remaining movable-double-root branch is not killed by the target
    # row's terminal P_3/Q_1 summand.  The following degree-one invariant
    # polynomials belong to the exact bands
    #
    #   q_1=t*u^4*S(w),              deg S<=20,
    #   p_3=t^3*u^2*(w-1)*A(w),      deg A<=11.
    #
    # They satisfy both endpoint normalizations.  At a nonzero double root
    # w0!=1, q_1 vanishes simply and the local layer-four bracket is exactly
    # t0^4.  Thus the first target jet is compatible; any exclusion has to
    # use the other layer-four summands and/or lower Laurent bands.
    w, w0, u0, t0 = sp.symbols("w w0 u0 t0", nonzero=True)
    S_w = sp.cancel(-(w - w0) / (1 - w0))
    A_w = sp.cancel(
        sp.Rational(1, 5)
        + (w - 1)
        / (w0 - 1)
        * (sp.Rational(1, 15) / w0**2 - sp.Rational(1, 5))
    )
    if sp.cancel(S_w.subs(w, 1) + 1) != 0:
        raise AssertionError("the Q band-one witness lost its endpoint value")
    if sp.cancel(S_w.subs(w, w0)) != 0:
        raise AssertionError("the Q band-one witness lost its movable zero")
    if sp.cancel(A_w.subs(w, 1) - sp.Rational(1, 5)) != 0:
        raise AssertionError("the P band-three witness lost its endpoint value")
    if sp.cancel(A_w.subs(w, w0) - sp.Rational(1, 15) / w0**2) != 0:
        raise AssertionError("the P band-three witness lost its local value")
    q_one_derivative_at_root = sp.cancel(
        5 * t0 * u0**8 * sp.diff(S_w, w).subs(w, w0)
    )
    p_three_at_root = sp.cancel(
        t0**3 * u0**2 * (w0 - 1) * A_w.subs(w, w0)
    )
    target_jet = sp.cancel(
        (3 * p_three_at_root * q_one_derivative_at_root - t0**4).subs(
            w0, u0**5
        )
    )
    if target_jet != 0:
        raise AssertionError("the movable-double-root target jet changed")

    # Row five for spacing eight is absolute descent 40, i.e. bracket layer
    # zero.  It is not the old primitive E5=0 equation.  The complete corner
    # envelope has new source variables p_-25 and q_-15 at that row, and the
    # target-tail pair p_11/q_-11 also contributes.  In particular q_-15 is
    # a 19-dimensional exact source band, so setting the fifth fractional-
    # power numerator to zero would silently delete a genuine lower tail.
    p_minus_25 = band_factor_data(p_bands[-25])
    q_minus_15 = band_factor_data(q_bands[-15])
    p_eleven = band_factor_data(p_bands[11])
    q_minus_11 = band_factor_data(q_bands[-11])
    if p_minus_25 != (25, 0, 8) or p_bands[-25].dimension != 9:
        raise AssertionError("the descent-forty P tail band changed")
    if q_minus_15 != (15, 0, 18) or q_bands[-15].dimension != 19:
        raise AssertionError("the descent-forty Q tail band changed")
    if p_eleven != (4, 5, 6) or p_bands[11].dimension != 7:
        raise AssertionError("the target-to-row-forty P band changed")
    if q_minus_11 != (11, 0, 19) or q_bands[-11].dimension != 20:
        raise AssertionError("the target Q tail band changed")

    # Exhaust the possible *positions* of a first non-root P defect relative
    # to the target.  The P envelope ends at relative descent 90.  Only
    # multiples strictly below 36 are primitive zero rows; at and below the
    # target, noncentral tail terms enter.  This makes clear that descent 8 is
    # the earliest survivor, not the only remaining branch position.
    first_defect_ledger: list[dict[str, object]] = []
    for spacing in range(1, 91):
        pure_rows = [
            multiple
            for multiple in range(2, 6)
            if multiple * spacing < target_descent
        ]
        target_rows = [
            multiple
            for multiple in range(2, 6)
            if multiple * spacing == target_descent
        ]
        first_defect_ledger.append(
            {
                "spacing": spacing,
                "P_band": 15 - spacing,
                "primitive_zero_multiples_before_target": pure_rows,
                "multiple_on_target": target_rows,
            }
        )
    if [
        row["primitive_zero_multiples_before_target"]
        for row in first_defect_ledger[7:18]
    ] != [
        [2, 3, 4],
        [2, 3],
        [2, 3],
        [2, 3],
        [2],
        [2],
        [2],
        [2],
        [2],
        [2],
        [],
    ]:
        raise AssertionError("the pre-target first-defect ledger changed")
    later_third_bands = []
    for spacing in range(9, 12):
        layer = 25 - 3 * spacing
        u_power, vanishing, free_degree = band_factor_data(q_bands[layer])
        if vanishing != 0:
            raise AssertionError("a later third-multiple Q band gained a fixed zero")
        later_third_bands.append(
            {
                "spacing": spacing,
                "Q_band": layer,
                "factorization": (
                    f"t^{layer}*u^{u_power}*S(u^5), "
                    f"deg(S)<={free_degree}"
                ),
                "dimension": q_bands[layer].dimension,
                "endpoint_normalized": False,
            }
        )

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
        "first_defect_position_ledger": {
            "maximum_spacing_from_full_P_envelope": 90,
            "exact_rows": first_defect_ledger,
            "compressed_regimes": [
                {
                    "spacing": "1..7",
                    "status": "excluded by the E4/E5 resultant",
                },
                {
                    "spacing": "8",
                    "pure_zero_multiples": [2, 3, 4],
                    "status": (
                        "earliest survivor; quadratic ratio and movable "
                        "double-R support only"
                    ),
                },
                {
                    "spacing": "9..11",
                    "pure_zero_multiples": [2, 3],
                    "local_consequence": (
                        "simple-prime defects return to the root slice; a "
                        "double prime retains v(C0)=2,v(p_k)=3,v(p_2k)=0"
                    ),
                    "status": (
                        "later double-prime positions remain because row four "
                        "is on or after the target"
                    ),
                },
                {
                    "spacing": "12..17",
                    "pure_zero_multiples": [2],
                    "local_consequence": "only C0 divides p_k^2 is forced",
                    "status": "later first-defect positions remain",
                },
                {
                    "spacing": "18..90",
                    "pure_zero_multiples": [],
                    "status": (
                        "no nonlinear multiple precedes the target; only the "
                        "full target-and-tail Fitting system can test them"
                    ),
                },
            ],
            "claim_boundary": (
                "descent 8 is the first surviving branch, not an exhaustive "
                "classification of all later first-defect positions"
            ),
            "unnormalized_third_multiple_bands_for_spacing_9_to_11": (
                later_third_bands
            ),
        },
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
            "relative_Artinian_core": {
                "base": "B=QQ[w0,w0^-1,(w0-1)^-1]",
                "double_R": "R(w)=(w-w0)^2/(25*(1-w0)^2)",
                "algebra": "B[y]/(27*y^2-9*y+1)",
                "rank_over_base": 2,
                "multiplication_by_y_matrix_in_basis_1_y": [
                    ["0", "-1/27"],
                    ["1", "1/3"],
                ],
                "norm_of_y": "1/27",
                "etale_discriminant": "-1/27",
                "geometric_components_after_QQ_sqrt_minus_3": 2,
            },
            "kummer_orbit_filter": {
                "Q_band_1": "q1=t*u^4*S(u^5), deg(S)<=20",
                "terminal_normalization": "S(1)=[t^1]q1=-1",
                "fixed_u5_equals_1_double_primes": (
                    "excluded because the local branch requires q1=0"
                ),
                "only_remaining_support": (
                    "R has a double root w0 and S(w0)=0"
                ),
            },
            "remaining_double_R_target_jet_witness": {
                "Q_band_1_invariant": str(S_w),
                "P_band_3_invariant": str(A_w),
                "endpoint_values": "S(1)=-1 and A(1)=1/5",
                "root_values": "S(w0)=0 and A(w0)=1/(15*w0^2)",
                "local_target_identity": (
                    "3*p3(t0)*q1'(t0)=t0^4 when u0=1+t0 and "
                    "w0=u0^5"
                ),
                "consequence": (
                    "the terminal P3/Q1 summand does not kill the remaining "
                    "movable-double-root branch"
                ),
            },
            "post_target_fitting_handoff": {
                "fifth_multiple_descent": 40,
                "bracket_layer": 0,
                "new_P_band": {
                    "layer": -25,
                    "factorization": "t^-25*u^25*A(u^5), deg(A)<=8",
                    "dimension": 9,
                },
                "new_Q_band": {
                    "layer": -15,
                    "factorization": "t^-15*u^15*S(u^5), deg(S)<=18",
                    "dimension": 19,
                },
                "target_cross_pair": (
                    "P band 11 (dimension 7) with Q band -11 "
                    "(dimension 20) also lies on layer 0"
                ),
                "linear_operator": (
                    "T_40(q)=5*C0*q'+15*C0'*q="
                    "5*C0^-2*(C0^3*q)'"
                ),
                "consequence": (
                    "E5=0 is not a valid post-target equation; the exact "
                    "condition is the layer-zero Fitting condition with the "
                    "new lower source bands retained"
                ),
            },
            "claim_boundary": (
                "this is an exact local residual seed, not a reconstructed "
                "global source pair or plane counterexample"
            ),
        },
    }


def target_layer_fitting_audit(
    p_bands: dict[int, Band],
    q_bands: dict[int, Band],
) -> dict[str, object]:
    """Give the complete new-band cokernel on target Laurent layer four."""

    w, w0 = sp.symbols("w w0", nonzero=True)
    E = sp.expand((w - 1) * (w - w0))
    E_prime = sp.diff(E, w)

    if band_factor_data(p_bands[-21]) != (21, 0, 9):
        raise AssertionError("the target P_-21 band changed")
    if band_factor_data(q_bands[-11]) != (11, 0, 19):
        raise AssertionError("the target Q_-11 band changed")
    if p_bands[-21].dimension != 10 or q_bands[-11].dimension != 20:
        raise AssertionError("the target new-band dimensions changed")

    # After q_-11 is replaced by its residual modulo the P follower, write
    # q=t^-11*u^11*S(w) and C=t^5*D(w), with D a base-unit multiple of E^2.
    # Direct cancellation of the t-logarithmic derivatives gives
    #
    #   (15*C^3*q' + 33*C^2*C'*q)/t^4
    #     = 15*w^2*E^5*N(S)
    #
    # up to the cube of that base unit, where N is the operator below.
    S = sp.Function("S")
    normalized_operator = sp.expand(
        5 * w * E * sp.diff(S(w), w)
        + (11 * E + 22 * w * E_prime) * S(w)
    )

    # N maps degree <=19 to degree <=21.  In the monomial bases its j-th
    # column has only three entries:
    #
    #   (5*j+11)*w0 at row j,
    #   -(5*j+33)*(1+w0) at row j+1,
    #   5*(j+11) at row j+2.
    #
    # The first twenty rows form a lower-triangular full-rank minor.  Solve
    # them successively; rows twenty and twenty-one are the two remaining
    # exact Fitting residues.
    matrix = sp.zeros(22, 20)
    for column in range(20):
        image = sp.Poly(
            sp.expand(
                5 * w * E * sp.diff(w**column, w)
                + (11 * E + 22 * w * E_prime) * w**column
            ),
            w,
        )
        for (degree,), coefficient in image.terms():
            matrix[degree, column] = coefficient
        if matrix[column, column] != (5 * column + 11) * w0:
            raise AssertionError("the target triangular diagonal changed")
        if matrix[column + 1, column] != -(5 * column + 33) * (1 + w0):
            raise AssertionError("the target triangular middle row changed")
        if matrix[column + 2, column] != 5 * (column + 11):
            raise AssertionError("the target triangular upper row changed")
    triangular_determinant = sp.factor(matrix[:20, :].det())
    expected_determinant = sp.factor(
        w0**20 * sp.prod(5 * index + 11 for index in range(20))
    )
    if triangular_determinant != expected_determinant:
        raise AssertionError("the target full-rank minor changed")

    forcing = sp.symbols("f0:22")
    solved: list[sp.Expr] = []
    for degree in range(20):
        previous_one = solved[degree - 1] if degree >= 1 else sp.Integer(0)
        previous_two = solved[degree - 2] if degree >= 2 else sp.Integer(0)
        solved.append(
            sp.cancel(
                (
                    forcing[degree]
                    + (5 * degree + 28) * (1 + w0) * previous_one
                    - 5 * (degree + 9) * previous_two
                )
                / ((5 * degree + 11) * w0)
            )
        )
    residue_twenty = sp.cancel(
        forcing[20] + 128 * (1 + w0) * solved[19] - 145 * solved[18]
    )
    residue_twenty_one = sp.cancel(forcing[21] - 150 * solved[19])
    reconstructed = matrix * sp.Matrix(solved)
    if any(
        sp.cancel(reconstructed[row] - forcing[row]) != 0
        for row in range(20)
    ):
        raise AssertionError("the target triangular reconstruction changed")
    if sp.cancel(reconstructed[20] - forcing[20] + residue_twenty) != 0:
        raise AssertionError("the first target residue changed")
    if sp.cancel(reconstructed[21] - forcing[21] + residue_twenty_one) != 0:
        raise AssertionError("the second target residue changed")

    # Before the earlier triangular rows are imposed, old source-band
    # generators already span the whole fourteen-dimensional target
    # cokernel.  Two explicit 34-by-34 minors certify this over the localized
    # double-root base.  Their gcd is a power of w0, hence a unit in B.
    def coefficient_vector(expression: sp.Expr) -> list[sp.Expr]:
        polynomial = sp.Poly(sp.expand(expression), w)
        return [polynomial.nth(degree) for degree in range(34)]

    full_new_columns = [
        coefficient_vector(
            w**2
            * E**5
            * (
                5 * w * E * sp.diff(w**degree, w)
                + (11 * E + 22 * w * E_prime) * w**degree
            )
        )
        for degree in range(20)
    ]

    def prior_generator(metadata: tuple[int, int, int]) -> list[sp.Expr]:
        p_layer, p_free_power, q_free_power = metadata
        q_layer = 4 - p_layer
        p_u, p_vanishing, p_degree = band_factor_data(p_bands[p_layer])
        q_u, q_vanishing, q_degree = band_factor_data(q_bands[q_layer])
        if not (0 <= p_free_power <= p_degree):
            raise AssertionError("a target witness left its P band")
        if not (0 <= q_free_power <= q_degree):
            raise AssertionError("a target witness left its Q band")
        if (p_u + q_u - 1) % 5:
            raise AssertionError("a target witness lost its Kummer character")
        w_shift = (p_u + q_u - 1) // 5
        p_factor = (w - 1) ** p_vanishing * w**p_free_power
        q_factor = (w - 1) ** q_vanishing * w**q_free_power
        normalized_bracket = sp.expand(
            w**w_shift
            * (
                p_layer
                * p_factor
                * (q_u * q_factor + 5 * w * sp.diff(q_factor, w))
                - q_layer
                * (p_u * p_factor + 5 * w * sp.diff(p_factor, w))
                * q_factor
            )
        )
        return coefficient_vector(normalized_bracket)

    witness_one = [
        (-20, 0, 0),
        (-20, 0, 1),
        (-20, 0, 2),
        (-20, 0, 3),
        (-20, 0, 4),
        (-20, 0, 5),
        (-20, 0, 6),
        (-9, 0, 0),
        (-6, 0, 0),
        (-6, 0, 1),
        (-4, 0, 0),
        (-1, 0, 0),
        (-1, 0, 1),
        (3, 0, 0),
    ]
    witness_two = [
        (14, 5, 19),
        (14, 5, 18),
        (14, 5, 17),
        (14, 5, 16),
        (14, 5, 15),
        (14, 5, 14),
        (14, 5, 13),
        (12, 6, 19),
        (10, 8, 19),
        (10, 0, 0),
        (7, 8, 20),
        (5, 10, 20),
        (5, 0, 0),
        (3, 11, 20),
    ]

    def witness_minor(witness: list[tuple[int, int, int]]) -> sp.Poly:
        columns = full_new_columns + [prior_generator(item) for item in witness]
        witness_matrix = sp.Matrix(
            34,
            34,
            lambda row, column: columns[column][row],
        )
        return sp.Poly(
            sp.factor(witness_matrix.det(method="domain-ge")),
            w0,
            domain=sp.QQ,
        )

    minor_one = witness_minor(witness_one)
    minor_two = witness_minor(witness_two)
    minor_gcd = sp.monic(sp.gcd(minor_one, minor_two))
    if minor_gcd.as_expr() != w0**12:
        raise AssertionError("the two target saturation minors lost their unit gcd")
    factors_one = [
        [str(sp.factor(factor.as_expr())), exponent]
        for factor, exponent in sp.factor_list(minor_one)[1]
    ]
    factors_two = [
        [str(sp.factor(factor.as_expr())), exponent]
        for factor, exponent in sp.factor_list(minor_two)[1]
    ]

    # The former P_3/Q_1 witness is necessarily only local.  Its normalized
    # target contribution has a factor w, while the new-band image has a
    # factor w^2.  Evaluation at w=0 therefore gives zero rather than the
    # target value one.  Across the complete older envelope exactly two band
    # pairs can supply that missing constant coefficient.
    constant_pair_coefficients: dict[int, int] = {}
    for p_layer in range(-21, 16):
        q_layer = 4 - p_layer
        p_u, p_vanishing, _ = band_factor_data(p_bands[p_layer])
        q_u, q_vanishing, _ = band_factor_data(q_bands[q_layer])
        if (p_u + q_u - 1) % 5:
            raise AssertionError("a target endpoint pair lost its character")
        w_shift = (p_u + q_u - 1) // 5
        if w_shift:
            continue
        p_factor = (w - 1) ** p_vanishing
        q_factor = (w - 1) ** q_vanishing
        generator = sp.expand(
            p_layer
            * p_factor
            * (q_u * q_factor + 5 * w * sp.diff(q_factor, w))
            - q_layer
            * (p_u * p_factor + 5 * w * sp.diff(p_factor, w))
            * q_factor
        )
        endpoint_value = sp.expand(generator.subs(w, 0))
        if endpoint_value:
            constant_pair_coefficients[p_layer] = int(endpoint_value)
    if constant_pair_coefficients != {-1: -5, 5: 5}:
        raise AssertionError("the target endpoint determinant changed")

    return {
        "target_layer": 4,
        "target_descent": 36,
        "new_bands": {
            "P_-21": {
                "factorization": "t^-21*u^21*A(w), deg(A)<=9",
                "dimension": 10,
            },
            "Q_-11": {
                "factorization": "t^-11*u^11*S(w), deg(S)<=19",
                "dimension": 20,
            },
            "combined_residual": (
                "S_-11+3*w^2*(w-1)^4*R(w)^2*A_-21"
            ),
        },
        "double_R_operator": {
            "E": "(w-1)*(w-w0)",
            "formula_up_to_a_base_unit": (
                "15*w^2*E^5*(5*w*E*S'+(11*E+22*w*E')*S)"
            ),
            "forced_local_factor": "w^2*(w-1)^5*(w-w0)^5",
            "forced_local_jet_conditions": 12,
        },
        "reduced_differential_map": {
            "domain": "B[w]_(<=19)",
            "codomain": "B[w]_(<=21)",
            "matrix_shape": [22, 20],
            "rank": 20,
            "cokernel_rank": 2,
            "column_j": {
                "row_j": "(5*j+11)*w0",
                "row_j_plus_1": "-(5*j+33)*(1+w0)",
                "row_j_plus_2": "5*(j+11)",
            },
            "full_rank_minor": str(expected_determinant),
            "triangular_solution": (
                "s_n=(f_n+(5*n+28)*(1+w0)*s_(n-1)-"
                "5*(n+9)*s_(n-2))/((5*n+11)*w0), 0<=n<=19"
            ),
            "two_final_residues": [
                "rho20=f20+128*(1+w0)*s19-145*s18",
                "rho21=f21-150*s19",
            ],
        },
        "complete_target_fitting_statement": {
            "condition_count": 14,
            "conditions": (
                "after moving all older band pairs and the target 1 to the "
                "forcing side, require divisibility by "
                "w^2*(w-1)^5*(w-w0)^5 and then rho20=rho21=0"
            ),
            "relation_to_previous_witness": (
                "the P_3/Q_1 identity checks one leading movable-root jet; "
                "it is one coordinate of this fourteen-condition cokernel"
            ),
        },
        "raw_B0_target_saturation": {
            "new_operator_columns": 20,
            "old_generator_columns_per_witness": 14,
            "ambient_rank": 34,
            "first_witness": [list(item) for item in witness_one],
            "second_witness": [list(item) for item in witness_two],
            "first_minor_nonconstant_factors": factors_one,
            "second_minor_nonconstant_factors": factors_two,
            "minor_gcd": "w0^12",
            "localized_consequence": (
                "after inverting w0, the old source-basis generators plus "
                "the new-band operator span the complete target row space"
            ),
            "claim_boundary": (
                "bilinear source-basis span is not simultaneous nonlinear "
                "solvability after the earlier triangular rows"
            ),
        },
        "globalization_of_the_local_target_witness": {
            "P3_Q1_contribution_at_w0": 0,
            "new_target_band_image_at_w0": 0,
            "target_value_at_w0": 1,
            "consequence": (
                "P_3/Q_1 plus P_-21/Q_-11 cannot satisfy the full target "
                "polynomial identity"
            ),
            "only_old_pairs_with_nonzero_constant_term": [
                {"P_layer": -1, "Q_layer": 5, "coefficient": -5},
                {"P_layer": 5, "Q_layer": -1, "coefficient": 5},
            ],
            "forced_endpoint_equation": (
                "5*(A_P5(0)*B_Q-1(0)-A_P-1(0)*B_Q5(0))=1"
            ),
            "interpretation": (
                "any globalization of the movable-root target jet must "
                "activate off-grid descents (P,Q)=(10,26) and/or (16,20)"
            ),
        },
    }


def u_zero_edge_audit() -> dict[str, object]:
    """Extract the lowest-u edge and test its sparse formal completion."""

    x, e, v, s = sp.symbols("x e v s", nonzero=True)
    # At u=0 the u-order-zero bands have z-exponents divisible by five,
    # while the u-order-one bands have exponents congruent to -1.  Thus
    #
    #   P=A(x)+u*z^-1*B(x)+...,  Q=C(x)+u*z^-1*D(x)+..., x=z^5.
    #
    # The u^0 bracket is z^-1*5*x*(A'*D-B*C').  Since the target is
    # z^4*(u-1)^4 and its u^0 coefficient is z^-1*x, the exact edge equation
    # is A'*D-B*C'=1/5.
    A = x + e**3 * x**3
    C = -sp.Rational(9, 5) * e**5 * x**5
    B = e / 5
    D = sp.Rational(1, 5) - sp.Rational(3, 5) * e**3 * x**2
    bezout = sp.factor(sp.diff(A, x) * D - B * sp.diff(C, x))
    if bezout != sp.Rational(1, 5):
        raise AssertionError("the u=0 edge witness changed")
    if [sp.degree(poly, x) for poly in (A, B, C, D)] != [3, 0, 5, 2]:
        raise AssertionError("the u=0 edge witness left its degree bounds")

    # If the complete u-order-zero edge retained a nonconstant common root
    # A=E^3 and C=(-9/5)E^5, then both terms in the Bezout left side would be
    # divisible by E^2.  The constant target therefore forces a break from
    # the common-power edge somewhere below its top coefficient.
    root, root_prime, arbitrary_B, arbitrary_D = sp.symbols(
        "root root_prime arbitrary_B arbitrary_D"
    )
    common_power_left = sp.factor(
        3 * root**2 * root_prime * arbitrary_D
        - arbitrary_B * (-9 * root**4 * root_prime)
    )
    if sp.factor(common_power_left / root**2) != sp.factor(
        3 * root_prime * arbitrary_D
        + 9 * arbitrary_B * root**2 * root_prime
    ):
        raise AssertionError("the edge common-power factor changed")

    # The displayed r=3 witness is the specialization of a uniform sparse
    # edge section.  For m=2*r-1 and
    #
    #   beta_r=2*r^2/((r-1)*(2*r-1)),
    #
    # take
    #
    #   A=x+e^r*x^r,             B=(r-1)*e/10,
    #   C=-beta_r*e^m*x^m,       D=(1-r*e^r*x^(r-1))/5.
    #
    # Exact arithmetic samples verify A'*D-B*C'=1/5.  The literal pair
    # A+B*v, C+D*v is not a Keller pair: its next coefficient is already
    # nonzero.  If P is kept linear in v, the forced Q coefficient at v^2
    # is rational with denominator A', so the sparse polynomial completion
    # fails at transverse order two for every r>=2.
    family_samples: list[dict[str, object]] = []
    for r_value in range(2, 9):
        m_value = 2 * r_value - 1
        beta_value = sp.Rational(
            2 * r_value**2,
            (r_value - 1) * m_value,
        )
        family_A = x + e**r_value * x**r_value
        family_B = sp.Rational(r_value - 1, 10) * e
        family_C = -beta_value * e**m_value * x**m_value
        family_D = (
            1 - r_value * e**r_value * x ** (r_value - 1)
        ) / 5
        family_edge = sp.factor(
            sp.diff(family_A, x) * family_D
            - family_B * sp.diff(family_C, x)
        )
        if family_edge != sp.Rational(1, 5):
            raise AssertionError("the all-r sparse edge section changed")
        family_linear_jacobian = sp.factor(
            sp.diff(family_A + family_B * v, x)
            * sp.diff(family_C + family_D * v, v)
            - sp.diff(family_A + family_B * v, v)
            * sp.diff(family_C + family_D * v, x)
        )
        family_defect = sp.factor(
            family_linear_jacobian - sp.Rational(1, 5)
        )
        expected_defect = sp.Rational(
            r_value * (r_value - 1) ** 2,
            50,
        ) * e ** (r_value + 1) * x ** (r_value - 2) * v
        if family_defect != expected_defect:
            raise AssertionError("the sparse edge's next coefficient changed")
        forced_q_two = sp.cancel(
            family_B
            * sp.diff(family_D, x)
            / (2 * sp.diff(family_A, x))
        )
        q_two_numerator, q_two_denominator = sp.fraction(forced_q_two)
        if sp.Poly(q_two_denominator, x).degree() != r_value - 1:
            raise AssertionError("the forced second edge denominator changed")
        if sp.Poly(q_two_numerator, x).degree() != r_value - 2:
            raise AssertionError("the forced second edge numerator changed")
        family_samples.append(
            {
                "r": r_value,
                "m": m_value,
                "literal_next_defect": str(family_defect),
                "forced_Q_v2_denominator_degree": r_value - 1,
                "forced_Q_v2_numerator_degree": r_value - 2,
                "polynomial_if_P_has_no_v2_term": False,
            }
        )

    # There is nevertheless an exact formal completion.  Put
    # kappa=1/(5*B), let H be the unique formal series with
    #
    #   H(A(x))=C(x)+kappa*x,
    #
    # and set P=A+B*v, Q=-kappa*x+H(P).  Its ordinary (x,v) Jacobian is
    # B*kappa=1/5, and differentiating the defining identity recovers D as
    # Q_v(x,0).  This explains why the edge equation survives every local
    # first-jet test.  But H cannot be a polynomial: deg(H(A)) would be a
    # multiple of r, whereas deg(C+kappa*x)=2*r-1.
    kappa = sp.cancel(1 / (5 * B))
    target_composition = sp.expand(C + kappa * x)
    H_trunc = sp.Integer(0)
    h_coefficients: list[sp.Expr] = []
    for degree in range(1, 10):
        residual = sp.expand(H_trunc.subs(s, A) - target_composition)
        coefficient = sp.factor(-residual.coeff(x, degree))
        h_coefficients.append(coefficient)
        H_trunc = sp.expand(H_trunc + coefficient * s**degree)
    expected_h_coefficients = [
        1 / e,
        0,
        -e**2,
        0,
        sp.Rational(6, 5) * e**5,
        0,
        -3 * e**8,
        0,
        10 * e**11,
    ]
    if any(
        sp.factor(actual - expected) != 0
        for actual, expected in zip(
            h_coefficients,
            expected_h_coefficients,
            strict=True,
        )
    ):
        raise AssertionError("the r=3 formal edge shear changed")
    composition_error = sp.expand(
        H_trunc.subs(s, A) - target_composition
    )
    if any(composition_error.coeff(x, degree) != 0 for degree in range(1, 10)):
        raise AssertionError("the truncated edge shear lost triangularity")
    formal_P = A + B * v
    truncated_Q = -kappa * x + H_trunc.subs(s, formal_P)
    truncated_jacobian = sp.factor(
        sp.diff(formal_P, x) * sp.diff(truncated_Q, v)
        - sp.diff(formal_P, v) * sp.diff(truncated_Q, x)
    )
    if truncated_jacobian != sp.Rational(1, 5):
        raise AssertionError("the formal edge shear Jacobian changed")
    forced_q_two_r3 = sp.factor(
        B * sp.diff(D, x) / (2 * sp.diff(A, x))
    )
    expected_q_two_r3 = sp.factor(
        -3 * e**4 * x / (25 * (1 + 3 * e**3 * x**2))
    )
    if forced_q_two_r3 != expected_q_two_r3:
        raise AssertionError("the r=3 second transverse coefficient changed")

    return {
        "coordinate": "x=z^5 at u=1+t=0",
        "edge_expansions": {
            "P": "A(x)+u*z^-1*B(x)+O(u^2)",
            "Q": "C(x)+u*z^-1*D(x)+O(u^2)",
            "degree_bounds": {"A": 3, "B": 3, "C": 5, "D": 5},
        },
        "exact_target_equation": "A'(x)*D(x)-B(x)*C'(x)=1/5",
        "constant_coefficient": (
            "5*(A_P5(0)*B_Q-1(0)-A_P-1(0)*B_Q5(0))=1"
        ),
        "common_power_consequence": (
            "A=E^3 and C=(-9/5)*E^5 with nonconstant E is impossible, "
            "because the left side is divisible by E^2"
        ),
        "exact_edge_witness": {
            "top_parameter": "e=-R(0), a unit on the movable double-R base",
            "A": str(A),
            "B": str(B),
            "C": str(C),
            "D": str(D),
            "verified_left_side": "1/5",
            "claim_boundary": (
                "this solves the complete lowest-u edge equation, not the "
                "higher-u coefficients or the earlier Laurent layers"
            ),
        },
        "all_r_sparse_edge_section": {
            "parameters": "r>=2, m=2*r-1",
            "beta_r": "2*r^2/((r-1)*(2*r-1))",
            "A": "x+e^r*x^r",
            "B": "(r-1)*e/10",
            "C": "-beta_r*e^(2*r-1)*x^(2*r-1)",
            "D": "(1-r*e^r*x^(r-1))/5",
            "verified_edge_identity": "A'*D-B*C'=1/5",
            "literal_linear_v_defect": (
                "r*(r-1)^2*e^(r+1)*x^(r-2)*v/50"
            ),
            "forced_Q_v2_if_P_is_linear_in_v": (
                "B*D'/(2*A')=-r*(r-1)^2*e^(r+1)*x^(r-2)"
                "/(100*(1+r*e^r*x^(r-1)))"
            ),
            "consequence": (
                "the literal four-term sparse witness fails at transverse "
                "order two for every r>=2"
            ),
            "samples": family_samples,
        },
        "formal_shear_completion": {
            "auxiliary_coordinate": "v=u*z^-1",
            "definition": (
                "kappa=1/(5*B), H(A(x))=C(x)+kappa*x, "
                "P=A(x)+B*v, Q=-kappa*x+H(P)"
            ),
            "exact_auxiliary_jacobian": "d(P,Q)/d(x,v)=1/5",
            "first_jet": "Q(x,0)=C(x) and Q_v(x,0)=D(x)",
            "r3_H_through_degree_9": str(H_trunc),
            "r3_forced_Q_v2": str(forced_q_two_r3),
            "nonpolynomiality": (
                "if deg(H)=q then r*q=2*r-1, impossible for r>=2; "
                "therefore H has an infinite tail"
            ),
            "exact_failure_mode": (
                "the edge escape is formally integrable but not a finite "
                "polynomial completion with P kept linear in v"
            ),
            "claim_boundary": (
                "a general completion may add a P coefficient at v^2 and "
                "further off-grid bands; excluding those is exactly the "
                "remaining simultaneous Fitting problem"
            ),
        },
    }


def transverse_edge_completion_audit() -> dict[str, object]:
    """Classify the first polynomial repairs of the sparse edge witness.

    The auxiliary edge variables are ``x=z^5`` and ``v=u*z^-1``.  At
    transverse orders below five the Kummer-return monomial ``x*v^5`` has
    not yet appeared, so the coefficient bounds are simply ``deg(P_i)<=r``
    and ``deg(Q_i)<=2*r-1``.  This makes the first repairs a small exact
    Bezout recursion.
    """

    x, e, v = sp.symbols("x e v", nonzero=True)
    tau_zero, tau_one = sp.symbols("tau_zero tau_one")
    family_samples: list[dict[str, object]] = []

    # At order v, adding P_2*v^2 and Q_2*v^2 gives
    #
    #   A'*Q_2-C'*P_2=B*D'/2.
    #
    # A particular solution is uniform in r.  Since gcd(A',C')=1, every
    # other solution is obtained by adding (A'*T,C'*T).  The degree bounds
    # force deg(T)<=1.  The resulting quadratic truncation is nevertheless
    # never a Keller pair; the displayed low coefficients prove this without
    # a search.
    for r_value in range(2, 9):
        m_value = 2 * r_value - 1
        beta_value = sp.Rational(
            2 * r_value**2,
            (r_value - 1) * m_value,
        )
        A = x + e**r_value * x**r_value
        B = sp.Rational(r_value - 1, 10) * e
        C = -beta_value * e**m_value * x**m_value
        D = (1 - r_value * e**r_value * x ** (r_value - 1)) / 5
        A_prime = sp.diff(A, x)
        C_prime = sp.diff(C, x)
        common_scale = (
            sp.Rational(r_value * (r_value - 1) ** 2, 100)
            * e ** (r_value + 1)
        )
        P_two_particular = (
            -sp.Rational(r_value * (r_value - 1) ** 3, 200)
            * e ** (r_value + 2)
            * x ** (r_value - 2)
        )
        Q_two_particular = sp.expand(
            -common_scale * x ** (r_value - 2)
            + common_scale
            * r_value
            * e**r_value
            * x ** (2 * r_value - 3)
        )
        T = tau_zero + tau_one * x
        P_two = sp.expand(P_two_particular + A_prime * T)
        Q_two = sp.expand(Q_two_particular + C_prime * T)
        if sp.factor(
            A_prime * Q_two
            - C_prime * P_two
            - B * sp.diff(D, x) / 2
        ) != 0:
            raise AssertionError("the all-r second edge repair changed")
        if sp.Poly(P_two, x).degree() > r_value:
            raise AssertionError("the second P repair left its edge degree")
        if sp.Poly(Q_two, x).degree() > m_value:
            raise AssertionError("the second Q repair left its edge degree")

        quadratic_P = A + B * v + P_two * v**2
        quadratic_Q = C + D * v + Q_two * v**2
        quadratic_jacobian = sp.expand(
            sp.diff(quadratic_P, x) * sp.diff(quadratic_Q, v)
            - sp.diff(quadratic_P, v) * sp.diff(quadratic_Q, x)
        )
        if quadratic_jacobian.coeff(v, 0) != sp.Rational(1, 5):
            raise AssertionError("the repaired edge lost its target constant")
        if quadratic_jacobian.coeff(v, 1) != 0:
            raise AssertionError("the repaired edge lost its first cancellation")
        quadratic_tail = sp.Poly(quadratic_jacobian.coeff(v, 2), x)
        if r_value == 2:
            if sp.factor(quadratic_tail.nth(2) - sp.Rational(12, 5) * e**4 * tau_one) != 0:
                raise AssertionError("the r=2 quadratic top coefficient changed")
            if sp.factor(
                quadratic_tail.nth(1).subs(tau_one, 0)
                - sp.Rational(12, 5) * e**4 * tau_zero
            ) != 0:
                raise AssertionError("the r=2 quadratic middle coefficient changed")
            terminal_conflict = sp.factor(
                quadratic_tail.nth(0).subs({tau_zero: 0, tau_one: 0})
            )
            if terminal_conflict != -sp.Rational(3, 250) * e**6:
                raise AssertionError("the r=2 quadratic conflict changed")
            coefficient_cascade = [2, 1, 0]
        else:
            if sp.factor(quadratic_tail.nth(0) - tau_one / 5) != 0:
                raise AssertionError("the family quadratic constant changed")
            if sp.factor(
                quadratic_tail.nth(r_value - 2).subs(tau_one, 0)
                - sp.Rational(3 * r_value * (r_value - 1), 5)
                * e**r_value
                * tau_zero
            ) != 0:
                raise AssertionError("the family quadratic pivot changed")
            terminal_conflict = sp.factor(
                quadratic_tail.nth(2 * r_value - 4).subs(
                    {tau_zero: 0, tau_one: 0}
                )
            )
            expected_conflict = (
                -sp.Rational(3 * r_value**2 * (r_value - 1) ** 4, 1000)
                * e ** (2 * r_value + 2)
            )
            if terminal_conflict != expected_conflict:
                raise AssertionError("the family quadratic conflict changed")
            coefficient_cascade = [0, r_value - 2, 2 * r_value - 4]
        if terminal_conflict == 0:
            raise AssertionError("a quadratic edge truncation unexpectedly survived")
        family_samples.append(
            {
                "r": r_value,
                "m": m_value,
                "P2_degree_bound": r_value,
                "Q2_degree_bound": m_value,
                "repair_kernel_dimension": 2,
                "quadratic_truncation_pivot_exponents": coefficient_cascade,
                "quadratic_terminal_conflict": str(terminal_conflict),
            }
        )

    # Specialize to r=3 and continue the exact recursion through order four.
    # At these orders every homogeneous repair is still
    #
    #   (P_n,Q_n) -> (P_n+A'*T_n,Q_n+C'*T_n), deg(T_n)<=1.
    #
    # The cubic truncation is excluded by five displayed coefficients.  The
    # quartic truncation is excluded by a nine-generator Groebner reduction
    # of the v^4 row followed by three tiny v^5 remainders.
    r_value = 3
    A = x + e**3 * x**3
    B = e / 5
    C = -sp.Rational(9, 5) * e**5 * x**5
    D = (1 - 3 * e**3 * x**2) / 5
    A_prime = sp.diff(A, x)
    C_prime = sp.diff(C, x)
    inverse_A_mod_x4 = 1 - 3 * e**3 * x**2
    P_coefficients: list[sp.Expr] = [A, B]
    Q_coefficients: list[sp.Expr] = [C, D]
    repair_parameters: list[sp.Symbol] = []

    def edge_jacobian_coefficient(power: int) -> sp.Expr:
        result = sp.Integer(0)
        for p_index, p_coefficient in enumerate(P_coefficients):
            q_index = power + 1 - p_index
            if not (0 <= q_index < len(Q_coefficients)):
                continue
            q_coefficient = Q_coefficients[q_index]
            result += (
                q_index * sp.diff(p_coefficient, x) * q_coefficient
                - p_index * p_coefficient * sp.diff(q_coefficient, x)
            )
        return sp.expand(result)

    def adjoin_repair(
        order: int,
        constant_parameter: sp.Symbol,
        linear_parameter: sp.Symbol,
    ) -> None:
        forcing = sp.expand(-edge_jacobian_coefficient(order - 1) / order)
        q_particular = sp.rem(
            sp.Poly(sp.expand(forcing * inverse_A_mod_x4), x),
            sp.Poly(x**4, x),
        ).as_expr()
        p_particular = sp.cancel(
            (A_prime * q_particular - forcing) / C_prime
        )
        if sp.Poly(sp.denom(p_particular), x).degree() != 0:
            raise AssertionError("an r=3 edge repair became nonpolynomial in x")
        kernel = constant_parameter + linear_parameter * x
        P_coefficients.append(sp.expand(p_particular + A_prime * kernel))
        Q_coefficients.append(sp.expand(q_particular + C_prime * kernel))
        repair_parameters.extend([constant_parameter, linear_parameter])
        if edge_jacobian_coefficient(order - 1) != 0:
            raise AssertionError("an r=3 edge repair failed its pivot row")

    a_two, b_two, a_three, b_three, a_four, b_four = sp.symbols(
        "a_two b_two a_three b_three a_four b_four"
    )
    adjoin_repair(2, a_two, b_two)
    quadratic_tail = sp.Poly(edge_jacobian_coefficient(2), x)
    if sp.factor(quadratic_tail.nth(0) - b_two / 5) != 0:
        raise AssertionError("the r=3 quadratic first pivot changed")
    if sp.factor(quadratic_tail.nth(1).subs(b_two, 0) - 18 * e**3 * a_two / 5) != 0:
        raise AssertionError("the r=3 quadratic second pivot changed")
    if quadratic_tail.nth(2).subs({a_two: 0, b_two: 0}) != -sp.Rational(54, 125) * e**8:
        raise AssertionError("the r=3 quadratic final conflict changed")

    adjoin_repair(3, a_three, b_three)
    cubic_tail = sp.Poly(edge_jacobian_coefficient(3), x)
    if cubic_tail.nth(7) != 108 * e**8 * b_two**2:
        raise AssertionError("the cubic b_two square pivot changed")
    if sp.factor(cubic_tail.nth(5).subs(b_two, 0) - 108 * e**8 * a_two**2) != 0:
        raise AssertionError("the cubic a_two square pivot changed")
    if cubic_tail.nth(0).subs({a_two: 0, b_two: 0}) != b_three / 5:
        raise AssertionError("the cubic b_three pivot changed")
    cubic_a_three = sp.Rational(6, 125) * e**6
    if sp.factor(
        cubic_tail.nth(1).subs(
            {a_two: 0, b_two: 0, b_three: 0, a_three: cubic_a_three}
        )
    ) != 0:
        raise AssertionError("the cubic a_three solve changed")
    cubic_conflict = sp.factor(
        cubic_tail.nth(3).subs(
            {a_two: 0, b_two: 0, b_three: 0, a_three: cubic_a_three}
        )
    )
    if cubic_conflict != sp.Rational(108, 125) * e**12:
        raise AssertionError("the cubic terminal conflict changed")

    adjoin_repair(4, a_four, b_four)
    quartic_row = sp.Poly(edge_jacobian_coefficient(4), x)
    quartic_equations = [
        coefficient
        for coefficient in quartic_row.all_coeffs()
        if coefficient != 0
    ]
    quartic_basis = sp.groebner(
        quartic_equations,
        *repair_parameters,
        order="grevlex",
        method="f5b",
    )
    if len(quartic_basis.polys) != 9:
        raise AssertionError("the quartic edge-row Groebner length changed")
    quartic_relation = (
        69375 * a_two**2 * e**3
        - 4375 * b_four
        - 2763 * e**10
    )
    if quartic_basis.reduce(quartic_relation)[1] != 0:
        raise AssertionError("the quartic b_four relation changed")
    quintic_row = sp.Poly(edge_jacobian_coefficient(5), x)
    quintic_remainders = {
        exponent: sp.factor(quartic_basis.reduce(quintic_row.nth(exponent))[1])
        for exponent in (0, 1, 2)
    }
    expected_quintic_remainders = {
        0: -sp.Rational(36, 3125)
        * e**4
        * (9 * a_two * e**5 + 25 * a_four),
        1: -sp.Rational(216, 578125)
        * e**4
        * (1875 * b_four + 608 * e**10),
        2: -sp.Rational(54, 3125)
        * e**7
        * (249 * a_two * e**5 + 275 * a_four),
    }
    if any(
        sp.factor(quintic_remainders[exponent] - expected) != 0
        for exponent, expected in expected_quintic_remainders.items()
    ):
        raise AssertionError("the quartic-to-quintic edge remainders changed")
    elimination_determinant = 9 * 275 - 25 * 249
    if elimination_determinant != -3750:
        raise AssertionError("the quartic a_two/a_four determinant changed")
    quartic_conflict = sp.factor(
        1875 * (-sp.Rational(2763, 4375) * e**10)
        + 608 * e**10
    )
    if quartic_conflict != -sp.Rational(4033, 7) * e**10:
        raise AssertionError("the quartic terminal conflict changed")

    return {
        "all_r_second_transverse_repair": {
            "equation": "A'*Q2-C'*P2=B*D'/2",
            "particular_solution": {
                "P2": (
                    "-r*(r-1)^3*e^(r+2)*x^(r-2)/200"
                ),
                "Q2": (
                    "-r*(r-1)^2*e^(r+1)*x^(r-2)/100+"
                    "r^2*(r-1)^2*e^(2*r+1)*x^(2*r-3)/100"
                ),
            },
            "general_solution": (
                "P2=P2_part+A'*(tau0+tau1*x), "
                "Q2=Q2_part+C'*(tau0+tau1*x)"
            ),
            "kernel_dimension": 2,
            "degree_bounds": "deg(P2)<=r, deg(Q2)<=2*r-1",
            "quadratic_truncation": (
                "impossible for every r>=2 by the recorded three-coefficient "
                "cascade"
            ),
            "samples": family_samples,
        },
        "r3_transverse_termination": {
            "base": "QQ[e,e^-1]",
            "quadratic_truncation": {
                "status": "unit coefficient ideal",
                "pivot_exponents": [0, 1, 2],
                "final_conflict": "-54*e^8/125",
            },
            "cubic_truncation": {
                "status": "unit coefficient ideal",
                "pivot_exponents": [7, 5, 0, 1, 3],
                "forced_values": (
                    "b2=a2=b3=0, a3=6*e^6/125"
                ),
                "final_conflict": str(cubic_conflict),
            },
            "quartic_truncation": {
                "status": "unit coefficient ideal",
                "v4_groebner_generator_count": len(quartic_basis.polys),
                "v5_remainder_exponents": [0, 1, 2],
                "a2_a4_elimination_determinant": elimination_determinant,
                "final_conflict": str(quartic_conflict),
            },
            "minimum_unexcluded_transverse_order": 5,
            "first_kummer_return": "w=x*v^5",
            "consequence": (
                "any polynomial completion of the movable-double-root edge "
                "must activate the v^5 band, where the terminal binomial-jet "
                "relations first couple previously separated edge coefficients"
            ),
            "claim_boundary": (
                "this excludes transverse truncations through degree four; "
                "it neither constructs nor excludes a completion beginning "
                "at the first Kummer return v^5"
            ),
        },
    }


def first_kummer_return_audit(
    p_bands: dict[int, Band],
    q_bands: dict[int, Band],
) -> dict[str, object]:
    """Compile the first *source-compatible* transverse return.

    The ordinary edge recursion is closed only while the total transverse
    order is less than five.  At total order five, terms coming from five
    powers of ``t=u-1`` return to the same Kummer character.  This function
    reconstructs that wrap directly from the exact source bands and proves
    that its first-appearance matrix is surjective on the movable-double-R
    branch.
    """

    x, u, e, w, w0 = sp.symbols("x u e w w0", nonzero=True)

    # Edge coefficient i has degree 3+floor(i/5) on P and
    # 5+floor(i/5) on Q.  Through the first return, every coefficient is the
    # first-appearance coordinate of exactly one source band, except that
    # the i=5 coefficients in positive residue-zero bands reuse their K(w)
    # coefficient of order one.
    p_edge_symbols = [
        sp.symbols(f"p_{index}_0:{4 + index // 5}")
        for index in range(6)
    ]
    q_edge_symbols = [
        sp.symbols(f"q_{index}_0:{6 + index // 5}")
        for index in range(6)
    ]
    p_edge = [
        sum(
            coefficient * x**degree
            for degree, coefficient in enumerate(coefficients)
        )
        for coefficients in p_edge_symbols
    ]
    q_edge = [
        sum(
            coefficient * x**degree
            for degree, coefficient in enumerate(coefficients)
        )
        for coefficients in q_edge_symbols
    ]

    def band_u_jet(
        edge: list[sp.Expr],
        bands: dict[int, Band],
        layer: int,
    ) -> sp.Expr:
        u_power, _, _ = band_factor_data(bands[layer])
        if u_power > 5:
            return sp.Integer(0)
        edge_degree = (layer + u_power) // 5
        leading = sp.Poly(edge[u_power], x).nth(edge_degree)
        result = sp.Integer(0)
        for increment in range(6 - u_power):
            if u_power == 0 and increment == 5:
                # This is the first coefficient of K(u^5) not already seen
                # at u=0.  It is precisely the returned edge coefficient.
                coefficient = sp.Poly(edge[5], x).nth(edge_degree + 1)
            else:
                coefficient = (
                    (-1) ** increment
                    * sp.binomial(layer, increment)
                    * leading
                )
            result += coefficient * u ** (u_power + increment)
        return sp.expand(result)

    p_jets = {
        layer: band_u_jet(p_edge, p_bands, layer)
        for layer in p_bands
        if band_factor_data(p_bands[layer])[0] <= 5
    }
    q_jets = {
        layer: band_u_jet(q_edge, q_bands, layer)
        for layer in q_bands
        if band_factor_data(q_bands[layer])[0] <= 5
    }
    if sorted(p_jets) != list(range(-5, 16)):
        raise AssertionError("the first-return P source-band window changed")
    if sorted(q_jets) != list(range(-5, 26)):
        raise AssertionError("the first-return Q source-band window changed")

    # The u^4 coefficient has total transverse order five.  Its character
    # zero packet consists of Laurent layers -5,0,5,...,40.  Layer -5
    # vanishes identically, and the remaining packet is encoded as a
    # polynomial whose x^d row is Laurent layer 5*d.
    source_row = sp.Integer(0)
    layer_rows: dict[int, sp.Expr] = {}
    for total_layer in range(-5, 41, 5):
        coefficient = sp.Integer(0)
        for p_layer, p_jet in p_jets.items():
            q_layer = total_layer - p_layer
            if q_layer not in q_jets:
                continue
            q_jet = q_jets[q_layer]
            coefficient += (
                p_layer * p_jet * sp.diff(q_jet, u)
                - q_layer * sp.diff(p_jet, u) * q_jet
            )
        coefficient = sp.expand(coefficient).coeff(u, 4)
        layer_rows[total_layer] = coefficient
        if total_layer >= 0:
            source_row += coefficient * x ** (total_layer // 5)
    source_row = sp.expand(source_row)
    if layer_rows[-5] != 0:
        raise AssertionError("the first return acquired a layer-minus-five row")

    ordinary_row = sp.expand(
        sum(
            q_index * sp.diff(p_edge[p_index], x) * q_edge[q_index]
            - p_index
            * p_edge[p_index]
            * sp.diff(q_edge[q_index], x)
            for p_index in range(6)
            for q_index in range(6)
            if p_index + q_index == 5
        )
    )

    # If A=sum a_d*x^d and C=sum c_j*x^j, the missing wrap is the exact
    # fifth-binomial correction below.  It is characterized by the fact that
    # a pure t^(5*d)/t^(5*j) pair has zero Laurent bracket.
    correction = sp.Integer(0)
    for p_degree in range(4):
        p_coefficient = sp.Poly(p_edge[0], x).nth(p_degree)
        for q_degree in range(6):
            q_coefficient = sp.Poly(q_edge[0], x).nth(q_degree)
            correction += (
                25
                * p_coefficient
                * q_coefficient
                * (
                    p_degree * sp.binomial(5 * q_degree, 5)
                    - q_degree * sp.binomial(5 * p_degree, 5)
                )
                * x ** (p_degree + q_degree)
            )
    correction = sp.expand(correction)
    if sp.expand(source_row - 5 * ordinary_row - correction) != 0:
        raise AssertionError("the exact fifth-binomial wrap formula changed")

    # On the double-root branch use e=-R(0).  This sign is forced by
    # t=u-1=-1 at u=0.  The exact top bands determine, rather than free, the
    # x^4 coefficient of P_5 and the x^6 coefficient of Q_5.
    R_edge = -e * (w - w0) ** 2 / w0**2
    p_top_factor = sp.expand((w - 1) ** 6 * R_edge**3)
    q_top_factor = sp.expand((w - 1) ** 10 * R_edge**5)
    p_top_zero = sp.factor((-1) ** 15 * p_top_factor.subs(w, 0))
    q_top_zero = sp.factor(
        -sp.Rational(9, 5) * (-1) ** 25 * q_top_factor.subs(w, 0)
    )
    p_top_return = sp.factor(
        sp.binomial(15, 5) * p_top_factor.subs(w, 0)
        - sp.diff(p_top_factor, w).subs(w, 0)
    )
    q_top_return = sp.factor(
        -sp.Rational(9, 5)
        * (
            sp.binomial(25, 5) * q_top_factor.subs(w, 0)
            - sp.diff(q_top_factor, w).subs(w, 0)
        )
    )
    expected_p_top_return = -e**3 * (3009 + 6 / w0)
    expected_q_top_return = (
        sp.Rational(9, 5) * e**5 * (53140 + 10 / w0)
    )
    if p_top_zero != e**3:
        raise AssertionError("the signed P edge-top parameter changed")
    if q_top_zero != -sp.Rational(9, 5) * e**5:
        raise AssertionError("the signed Q edge-top parameter changed")
    if sp.factor(p_top_return - expected_p_top_return) != 0:
        raise AssertionError("the P first-return top coefficient changed")
    if sp.factor(q_top_return - expected_q_top_return) != 0:
        raise AssertionError("the Q first-return top coefficient changed")

    A = x + e**3 * x**3
    C = -sp.Rational(9, 5) * e**5 * x**5
    specialization: dict[sp.Symbol, sp.Expr] = {}
    for degree, coefficient in enumerate(p_edge_symbols[0]):
        specialization[coefficient] = sp.Poly(A, x).nth(degree)
    for degree, coefficient in enumerate(q_edge_symbols[0]):
        specialization[coefficient] = sp.Poly(C, x).nth(degree)
    specialization[p_edge_symbols[5][4]] = p_top_return
    specialization[q_edge_symbols[5][6]] = q_top_return
    specialized_row = sp.Poly(sp.expand(source_row.subs(specialization)), x)
    if sp.factor(specialized_row.nth(8)) != 0:
        raise AssertionError("the exact top bracket failed to cancel at v^5")

    free_return_variables = [
        *p_edge_symbols[5][:4],
        *q_edge_symbols[5][:6],
    ]
    return_matrix, _ = sp.linear_eq_to_matrix(
        [specialized_row.nth(degree) for degree in range(8)],
        free_return_variables,
    )
    pivot_columns = (1, 2, 4, 5, 6, 7, 8, 9)
    return_minor = sp.factor(return_matrix[:, pivot_columns].det())
    expected_return_minor = 3**5 * 5**16 * e**13
    if return_minor != expected_return_minor:
        raise AssertionError("the first-return Fitting minor changed")
    if return_matrix.rank() != 8:
        raise AssertionError("the first-return edge matrix lost surjectivity")

    # Check that the source-band coefficient map really supplies the claimed
    # conditional freedom through the *second* return.  A coefficient of
    # K(w) of index k enters triangularly with nonzero diagonal as long as
    # k is within the recorded free degree.  Only the common-power top bands
    # are fixed rather than free.
    def conditional_source_freedom(
        bands: dict[int, Band],
        transverse_order: int,
        edge_degree: int,
        fixed_top_layer: int,
    ) -> int:
        result = 0
        for degree in range(edge_degree + 1):
            layer = 5 * degree - transverse_order
            u_power, _, free_degree = band_factor_data(bands[layer])
            if (transverse_order - u_power) % 5:
                raise AssertionError("a returned coefficient changed character")
            return_index = (transverse_order - u_power) // 5
            if layer == fixed_top_layer:
                continue
            if not 0 <= return_index <= free_degree:
                raise AssertionError("a returned source coefficient lost freedom")
            result += 1
        return result

    return_profile: list[dict[str, int]] = []
    for transverse_order in range(5, 11):
        p_degree = 3 + transverse_order // 5
        q_degree = 5 + transverse_order // 5
        p_freedom = conditional_source_freedom(
            p_bands, transverse_order, p_degree, 15
        )
        q_freedom = conditional_source_freedom(
            q_bands, transverse_order, q_degree, 25
        )
        return_profile.append(
            {
                "transverse_order": transverse_order,
                "P_degree_bound": p_degree,
                "Q_degree_bound": q_degree,
                "conditional_P_freedom": p_freedom,
                "conditional_Q_freedom": q_freedom,
            }
        )
    expected_profile = [
        (4, 6),
        (5, 7),
        (5, 7),
        (5, 7),
        (5, 7),
        (5, 7),
    ]
    if [
        (row["conditional_P_freedom"], row["conditional_Q_freedom"])
        for row in return_profile
    ] != expected_profile:
        raise AssertionError("the first two source-return blocks changed")

    # For orders six through nine, and for the non-top part of order ten,
    # the new row is the same Bezout map A'*Q-C'*P with degree bounds 4 and
    # 6.  A unit 9-by-9 minor proves that every forcing of degree at most
    # eight can be absorbed.  At order ten the degree-nine top coefficient
    # cancels identically because the exact common-power top bracket is zero.
    p_free = sp.symbols("return_p_0:5")
    q_free = sp.symbols("return_q_0:7")
    bezout_return = sp.Poly(
        sp.expand(
            sp.diff(A, x)
            * sum(coefficient * x**degree for degree, coefficient in enumerate(q_free))
            - sp.diff(C, x)
            * sum(coefficient * x**degree for degree, coefficient in enumerate(p_free))
        ),
        x,
    )
    bezout_matrix, _ = sp.linear_eq_to_matrix(
        [bezout_return.nth(degree) for degree in range(9)],
        [*p_free, *q_free],
    )
    bezout_pivots = (3, 4, 5, 6, 7, 8, 9, 10, 11)
    bezout_minor = sp.factor(bezout_matrix[:, bezout_pivots].det())
    if bezout_minor != 81 * e**10:
        raise AssertionError("the post-return Bezout minor changed")
    if bezout_matrix.rank() != 9:
        raise AssertionError("the post-return Bezout map lost surjectivity")

    sparse_correction = sp.factor(
        correction.subs(specialization)
    )
    expected_sparse_correction = (
        -sp.Integer(2390625) * e**5 * x**6
        - sp.Integer(6496875) * e**8 * x**8
    )
    if sp.expand(sparse_correction - expected_sparse_correction) != 0:
        raise AssertionError("the sparse first-wrap correction changed")

    return {
        "coordinate_return": "w=x*v^5=u^5",
        "signed_top_parameter": {
            "definition": "e=-R(0)",
            "normalized_double_R_value": "e=-w0^2/(25*(1-w0)^2)",
            "unit_on_base": True,
        },
        "exact_first_wrap": {
            "source_row": "[u^4] of Laurent layers 0,5,...,40",
            "ordinary_edge_row": "5*[v^4] det(d(P,Q)/d(x,v))",
            "correction": (
                "Omega_5(A,C)=25*sum_(d,j) a_d*c_j*"
                "(d*binomial(5*j,5)-j*binomial(5*d,5))*x^(d+j)"
            ),
            "identity": "source_row=ordinary_edge_row+Omega_5(A,C)",
            "sparse_correction": str(sparse_correction),
            "layer_minus_five_row": "identically zero",
            "layer_40_row": "identically zero after exact top returns",
        },
        "fixed_top_returns": {
            "P5_x4": str(p_top_return),
            "Q5_x6": str(q_top_return),
        },
        "first_return_fitting_statement": {
            "base": "B=QQ[w0,w0^-1,(w0-1)^-1]",
            "rank_two_extension": "B[y]/(27*y^2-9*y+1)",
            "free_source_return_variables": 10,
            "equation_layers": [5 * degree for degree in range(8)],
            "matrix_shape": [8, 10],
            "rank": 8,
            "unit_maximal_minor": "3^5*5^16*e^13",
            "cokernel": 0,
            "solution_relative_dimension": 2,
            "consequence": (
                "base change to the rank-two descent-eight algebra preserves "
                "surjectivity; both conjugate local branches survive v^5"
            ),
        },
        "propagation_through_second_return": {
            "conditional_source_freedom": return_profile,
            "common_Bezout_matrix_shape": [9, 12],
            "common_Bezout_rank": 9,
            "unit_minor": "81*e^10",
            "consequence": (
                "the edge staircase also absorbs every first-appearance "
                "forcing through v^10; standalone transverse truncation "
                "cannot close the movable-double-R branch"
            ),
        },
        "claim_boundary": (
            "this is the complete first-appearance Kummer packet at u=0. "
            "It does not solve the remaining target and layer-zero Hermite "
            "conditions at w=1 and w=w0, nor the full lower Laurent tail"
        ),
    }


def source_return_lift_and_global_quotient_audit(
    p_bands: dict[int, Band],
    q_bands: dict[int, Band],
) -> dict[str, object]:
    """Decide whether the return pivots are support artifacts.

    Every conditional return coefficient has a triangular source lift

        t^ell*u^j0*(u^5-1)^nu*(u^5)^k.

    The unit minors used through ``v^10`` select particular such lifts.  We
    trace them back to original ``x^i*y^j`` monomials and verify that they
    are strictly below both certified supporting edges.  The same calculation
    gives closed unit-minor and source-lift formulas for every ``r>=2``.  We
    then remove the two controllable jets at ``w=0`` from each of the target
    and layer-zero obstruction modules by an exact Chinese-remainder
    calculation in the ``r=3`` candidate.
    """

    w, u, w0, x, e = sp.symbols("w u w0 x e", nonzero=True)

    # These are exactly the columns of the unit minors in
    # first_kummer_return_audit: the exceptional first return uses an 8-row
    # minor, while orders 6 through 10 use the common 9-row Bezout minor.
    pivot_profile: dict[int, dict[str, list[int]]] = {
        5: {"P": [1, 2], "Q": list(range(6))},
        **{
            order: {"P": [3, 4], "Q": list(range(7))}
            for order in range(6, 11)
        },
    }
    side_data = {
        "P": {"bands": p_bands, "degree": 75, "height": 3},
        "Q": {"bands": q_bands, "degree": 125, "height": 5},
    }
    endpoint_normalized_bands = {"P": {3, 15}, "Q": {1, 13, 25}}
    lift_records: list[dict[str, object]] = []
    lift_keys: set[tuple[str, int, int]] = set()
    for transverse_order, side_profile in pivot_profile.items():
        for side, edge_degrees in side_profile.items():
            bands = side_data[side]["bands"]
            source_degree = int(side_data[side]["degree"])
            terminal_height = int(side_data[side]["height"])
            for edge_degree in edge_degrees:
                layer = 5 * edge_degree - transverse_order
                band = bands[layer]
                u_power, vanishing, free_degree = band_factor_data(band)
                if (transverse_order - u_power) % 5:
                    raise AssertionError("a pivot lift changed Kummer character")
                return_index = (transverse_order - u_power) // 5
                if not 0 <= return_index <= free_degree:
                    raise AssertionError("a pivot left its exact source band")
                key = (side, layer, return_index)
                if key in lift_keys:
                    raise AssertionError("two return pivots reused one source lift")
                lift_keys.add(key)

                # If the band has a fixed terminal endpoint, replace w^k by
                # (w-1)*w^k.  This preserves both all previously exposed
                # K-coefficients and K(1), while its w^k coefficient is the
                # unit -1; hence it is a normalization-preserving lift of
                # the same return pivot.  The two affected pivots are
                # Q_7[x^4] on layer 13 and Q_9[x^2] on layer 1.
                normalization_preserving = (
                    layer in endpoint_normalized_bands[side]
                )
                effective_vanishing = vanishing + int(normalization_preserving)
                if normalization_preserving and return_index + 1 > free_degree:
                    raise AssertionError("a normalized pivot lacks a follower")

                source_terms: list[dict[str, object]] = []
                source_u_combination = sp.Integer(0)
                maximum_total_degree = 0
                for offset in range(effective_vanishing + 1):
                    source_i = band.source_i_min + return_index + offset
                    source_j = 5 * source_i - layer
                    coefficient = (
                        (-1) ** (effective_vanishing - offset)
                        * comb(effective_vanishing, offset)
                    )
                    if not band.source_i_min <= source_i <= band.source_i_max:
                        raise AssertionError("a pivot lift left its source band")
                    if min(source_i, source_j) < 0:
                        raise AssertionError("a pivot source lift became Laurent")
                    source_u_combination += coefficient * u**source_j
                    total_degree = source_i + source_j
                    maximum_total_degree = max(maximum_total_degree, total_degree)
                    if total_degree >= source_degree:
                        raise AssertionError("a pivot reached the total-degree edge")
                    source_terms.append(
                        {
                            "coefficient": coefficient,
                            "source_monomial": [source_i, source_j],
                            "total_degree": total_degree,
                        }
                    )
                expected_source_u_combination = sp.expand(
                    u**u_power
                    * (u**5 - 1) ** effective_vanishing
                    * (u**5) ** return_index
                )
                if sp.expand(
                    source_u_combination - expected_source_u_combination
                ) != 0:
                    raise AssertionError("a pivot failed exact source recomposition")

                # The largest terminal weight occurs at the smallest t power.
                # Strict positive slack proves that every translated descendant
                # lies below, rather than on, the terminal supporting edge.
                terminal_t_min = layer + effective_vanishing
                terminal_slack = (
                    terminal_height
                    - (17 * layer - 12 * terminal_t_min)
                )
                if terminal_slack <= 0:
                    raise AssertionError("a pivot touched the terminal edge")
                lift_records.append(
                    {
                        "side": side,
                        "transverse_order": transverse_order,
                        "edge_coefficient": f"{side}_{transverse_order}[x^{edge_degree}]",
                        "source_band": layer,
                        "source_basis": (
                            f"t^{layer}*u^{u_power}*(u^5-1)^"
                            f"{effective_vanishing}*(u^5)^{return_index}"
                        ),
                        "endpoint_normalization_preserved": normalization_preserving,
                        "source_terms": source_terms,
                        "maximum_total_degree": maximum_total_degree,
                        "total_degree_edge": source_degree,
                        "terminal_t_min": terminal_t_min,
                        "terminal_halfspace_slack": terminal_slack,
                        "strictly_below_both_certified_edges": True,
                    }
                )

    if len(lift_records) != 53:
        raise AssertionError("the return-pivot lift count changed")

    # The same minors have closed forms for the all-r sparse edge section
    #
    #   A=x+e^r*x^r,
    #   C=-2*r^2*e^(2*r-1)*x^(2*r-1)/((r-1)*(2*r-1)).
    #
    # At the first return, after the inhomogeneous fifth-binomial correction
    # is moved to the forcing side, the linear map is
    # 25*(A'*Q_5-C'*P_5).  Choose P degrees 1..r-1 and Q degrees 0..2r-1.
    # At orders 6 through 10 choose P degrees 3..r+1 and Q degrees 0..2r.
    # Direct block elimination gives the two determinant formulas recorded
    # below.  Exact matrices for r=2..8 guard the formulas.
    family_minor_samples: list[dict[str, object]] = []
    for r_value in range(2, 9):
        m_value = 2 * r_value - 1
        beta_value = sp.Rational(
            2 * r_value**2,
            (r_value - 1) * m_value,
        )
        family_A = x + e**r_value * x**r_value
        family_C = -beta_value * e**m_value * x**m_value
        c_value = sp.Rational(2 * r_value**2, r_value - 1)

        first_columns = [
            sp.expand(-25 * sp.diff(family_C, x) * x**degree)
            for degree in range(1, r_value)
        ] + [
            sp.expand(25 * sp.diff(family_A, x) * x**degree)
            for degree in range(2 * r_value)
        ]
        first_dimension = 3 * r_value - 1
        first_matrix = sp.Matrix(
            first_dimension,
            first_dimension,
            lambda row, column: sp.Poly(
                first_columns[column], x
            ).nth(row),
        )
        first_determinant = sp.factor(first_matrix.det())
        expected_first_determinant = sp.factor(
            (-1) ** (r_value - 1)
            * 25 ** (3 * r_value - 1)
            * r_value
            * c_value ** (r_value - 1)
            * e ** (2 * r_value * (r_value - 1) + 1)
        )
        if first_determinant != expected_first_determinant:
            raise AssertionError("the all-r first-return minor changed")

        post_columns = [
            sp.expand(-sp.diff(family_C, x) * x**degree)
            for degree in range(3, r_value + 2)
        ] + [
            sp.expand(sp.diff(family_A, x) * x**degree)
            for degree in range(2 * r_value + 1)
        ]
        post_dimension = 3 * r_value
        post_matrix = sp.Matrix(
            post_dimension,
            post_dimension,
            lambda row, column: sp.Poly(
                post_columns[column], x
            ).nth(row),
        )
        post_determinant = sp.factor(post_matrix.det())
        expected_post_determinant = sp.factor(
            (-1) ** (r_value - 1)
            * c_value ** (r_value - 1)
            * e ** ((r_value - 1) * (2 * r_value - 1))
        )
        if post_determinant != expected_post_determinant:
            raise AssertionError("the all-r post-return minor changed")
        family_minor_samples.append(
            {
                "r": r_value,
                "first_return_shape": [first_dimension, first_dimension],
                "first_return_minor": str(first_determinant),
                "post_return_shape": [post_dimension, post_dimension],
                "post_return_minor": str(post_determinant),
            }
        )

    # The selected all-r columns also have exact corner-compatible source
    # lifts.  The family terminal halfspace is
    #
    #   (7*r-4)*ell-(5*r-3)*a <= h,
    #
    # where a is the t exponent and h is r on P or 2*r-1 on Q.  The only
    # selected normalized bands are Q_(5*r-2) at order 7, degree r+1, and
    # Q_1 at order 9, degree 2; multiplying those directions by w-1
    # preserves their endpoint values.
    family_source_lift_samples: list[dict[str, object]] = []
    for r_value in range(2, 13):
        family_pivots = {
            5: {
                "P": range(1, r_value),
                "Q": range(2 * r_value),
            },
            **{
                order: {
                    "P": range(3, r_value + 2),
                    "Q": range(2 * r_value + 1),
                }
                for order in range(6, 11)
            },
        }
        family_lift_count = 0
        maximum_total_degrees = {"P": 0, "Q": 0}
        minimum_terminal_slack: int | None = None
        normalized_pivots: list[dict[str, int | str]] = []
        for order, side_pivots in family_pivots.items():
            for side, degrees in side_pivots.items():
                exponent = r_value if side == "P" else 2 * r_value - 1
                source_degree = 25 * exponent
                height = exponent
                normalized_layers = (
                    {3, 5 * r_value}
                    if side == "P"
                    else {
                        1,
                        5 * r_value - 2,
                        5 * (2 * r_value - 1),
                    }
                )
                for edge_degree in degrees:
                    layer = 5 * edge_degree - order
                    source_i_min = max(0, ceil_div(layer, 5))
                    source_i_max = (source_degree + layer) // 6
                    u_power = 5 * source_i_min - layer
                    terminal_t_min = max(
                        layer,
                        ceil_div(
                            (7 * r_value - 4) * layer - height,
                            5 * r_value - 3,
                        ),
                    )
                    vanishing = terminal_t_min - layer
                    free_degree = (
                        source_i_max - source_i_min - vanishing
                    )
                    if (order - u_power) % 5:
                        raise AssertionError("an all-r pivot changed character")
                    return_index = (order - u_power) // 5
                    normalized = layer in normalized_layers
                    if not 0 <= return_index <= free_degree:
                        raise AssertionError("an all-r pivot left its band")
                    if normalized and return_index + 1 > free_degree:
                        raise AssertionError(
                            "an all-r normalized pivot lost its follower"
                        )
                    effective_vanishing = vanishing + int(normalized)
                    largest_i = (
                        source_i_min + return_index + effective_vanishing
                    )
                    largest_j = 5 * largest_i - layer
                    maximum_total_degree = largest_i + largest_j
                    if maximum_total_degree >= source_degree:
                        raise AssertionError(
                            "an all-r pivot reached the degree edge"
                        )
                    maximum_total_degrees[side] = max(
                        maximum_total_degrees[side],
                        maximum_total_degree,
                    )
                    terminal_slack = height - (
                        (7 * r_value - 4) * layer
                        - (5 * r_value - 3)
                        * (layer + effective_vanishing)
                    )
                    if terminal_slack <= 0:
                        raise AssertionError(
                            "an all-r pivot reached the terminal edge"
                        )
                    minimum_terminal_slack = (
                        terminal_slack
                        if minimum_terminal_slack is None
                        else min(minimum_terminal_slack, terminal_slack)
                    )
                    family_lift_count += 1
                    if normalized:
                        normalized_pivots.append(
                            {
                                "side": side,
                                "transverse_order": order,
                                "edge_degree": edge_degree,
                                "source_layer": layer,
                            }
                        )
        if family_lift_count != 18 * r_value - 1:
            raise AssertionError("the all-r source-lift count changed")
        if maximum_total_degrees != {
            "P": 13 * r_value + 8,
            "Q": 26 * r_value - 5,
        }:
            raise AssertionError("the all-r source-degree maxima changed")
        if minimum_terminal_slack != 1:
            raise AssertionError("the all-r minimum terminal slack changed")
        expected_normalized_pivots = [
            {
                "side": "Q",
                "transverse_order": 7,
                "edge_degree": r_value + 1,
                "source_layer": 5 * r_value - 2,
            },
            {
                "side": "Q",
                "transverse_order": 9,
                "edge_degree": 2,
                "source_layer": 1,
            },
        ]
        if normalized_pivots != expected_normalized_pivots:
            raise AssertionError("the all-r normalized pivots changed")
        family_source_lift_samples.append(
            {
                "r": r_value,
                "pivot_source_lift_count": family_lift_count,
                "maximum_total_degrees": maximum_total_degrees,
                "total_degree_bounds": {
                    "P": 25 * r_value,
                    "Q": 25 * (2 * r_value - 1),
                },
                "minimum_terminal_halfspace_slack": minimum_terminal_slack,
                "normalized_pivots": normalized_pivots,
            }
        )

    # Confluent evaluation is the exact CRT map for the three pairwise
    # disjoint points 0,1,w0.  Its determinant is a unit in
    # B=QQ[w0,w0^-1,(w0-1)^-1].
    def confluent_matrix(
        multiplicities: list[tuple[sp.Expr, int]],
    ) -> sp.Matrix:
        dimension = sum(multiplicity for _, multiplicity in multiplicities)
        matrix = sp.zeros(dimension, dimension)
        row = 0
        for point, multiplicity in multiplicities:
            for derivative in range(multiplicity):
                for degree in range(dimension):
                    matrix[row, degree] = sp.diff(
                        w**degree, w, derivative
                    ).subs(w, point)
                row += 1
        return matrix

    target_crt = confluent_matrix([(0, 2), (1, 5), (w0, 5)])
    target_crt_determinant = sp.factor(target_crt.det())
    expected_target_crt = 2**10 * 3**4 * w0**10 * (w0 - 1) ** 25
    if target_crt_determinant != expected_target_crt:
        raise AssertionError("the target CRT determinant changed")

    layer_zero_crt = confluent_matrix([(0, 3), (1, 6), (w0, 6)])
    layer_zero_crt_determinant = sp.factor(layer_zero_crt.det())
    expected_layer_zero_crt = (
        2**17 * 3**6 * 5**2 * w0**18 * (w0 - 1) ** 36
    )
    if layer_zero_crt_determinant != expected_layer_zero_crt:
        raise AssertionError("the layer-zero CRT determinant changed")

    # On target layer four, after the common factor t^4 is removed, the
    # local character-zero remainder is G(w).  Its constant and linear
    # coefficients occur at u-orders 0 and 5.  In the edge staircase these
    # are source transverse orders 1 and 6: the exact edge Bezout equation
    # controls the former, and the order-six unit Bezout block controls the
    # latter.  (The order-five packet itself belongs to Laurent layers
    # divisible by five.)
    g0, g1 = sp.symbols("g0 g1")
    target_two_jet = sp.expand((g0 + g1 * w).subs(w, u**5))
    if target_two_jet.coeff(u, 0) != g0:
        raise AssertionError("the target edge row lost G(0)")
    if target_two_jet.coeff(u, 5) != g1:
        raise AssertionError("the order-six target return lost G'(0)")

    # Layer zero has J_0=dH/dt=dH/du.  Its two character-zero return rows are
    # exactly the two nonconstant jets at w=0, reached at source transverse
    # orders five and ten.
    h0, h1, h2 = sp.symbols("h0 h1 h2")
    H_two_jet = h0 + h1 * w + h2 * w**2
    layer_zero_u_jet = sp.expand(
        sp.diff(H_two_jet.subs(w, u**5), u)
    )
    if layer_zero_u_jet.coeff(u, 4) != 5 * h1:
        raise AssertionError("the v^5 row lost H'(0)")
    if layer_zero_u_jet.coeff(u, 9) != 10 * h2:
        raise AssertionError("the v^10 row lost H''(0)")

    return {
        "support_decision": {
            "pivot_source_lift_count": len(lift_records),
            "all_lifts_are_original_polynomials": True,
            "all_lifts_strictly_below_total_degree_edge": True,
            "all_lifts_strictly_below_terminal_edge": True,
            "consequence": (
                "the columns making the v^5-through-v^10 return maps "
                "surjective are genuine source-polynomial directions, not "
                "independent-coefficient artifacts; no refinement using only "
                "the two certified supporting edges can delete them"
            ),
            "claim_boundary": (
                "an additional lower Newton edge not supplied by the corner "
                "chain could still set some interior coefficients to zero"
            ),
        },
        "pivot_source_lifts": lift_records,
        "all_r_return_no_go": {
            "parameters": "r>=2, m=2*r-1, e invertible",
            "edge_section": {
                "A": "x+e^r*x^r",
                "C": (
                    "-2*r^2*e^(2*r-1)*x^(2*r-1)/"
                    "((r-1)*(2*r-1))"
                ),
            },
            "first_return": {
                "linear_map": "25*(A'*Q5-C'*P5)",
                "pivot_degrees": {
                    "P": "1..r-1",
                    "Q": "0..2*r-1",
                },
                "shape": "(3*r-1)-by-(3*r-1)",
                "minor": (
                    "(-1)^(r-1)*25^(3*r-1)*r*"
                    "(2*r^2/(r-1))^(r-1)*e^(2*r*(r-1)+1)"
                ),
            },
            "orders_6_through_10": {
                "linear_map": "A'*Q-C'*P",
                "pivot_degrees": {"P": "3..r+1", "Q": "0..2*r"},
                "shape": "(3*r)-by-(3*r)",
                "minor": (
                    "(-1)^(r-1)*(2*r^2/(r-1))^(r-1)*"
                    "e^((r-1)*(2*r-1))"
                ),
            },
            "source_lift_count_through_order_10": "18*r-1",
            "source_degree_maxima": {"P": "13*r+8", "Q": "26*r-5"},
            "source_degree_gaps": {"P": "12*r-8", "Q": "24*r-20"},
            "minimum_terminal_halfspace_slack": 1,
            "consequence": (
                "for every r>=2 the standalone sparse-edge staircase is "
                "surjective through the second Kummer return by genuine "
                "strict-interior polynomial source directions"
            ),
            "minor_regression_samples": family_minor_samples,
            "source_lift_regression_samples": family_source_lift_samples,
        },
        "three_point_CRT": {
            "base": "B=QQ[w0,w0^-1,(w0-1)^-1]",
            "target_multiplicities": {"0": 2, "1": 5, "w0": 5},
            "target_determinant": str(target_crt_determinant),
            "layer_zero_multiplicities": {"0": 3, "1": 6, "w0": 6},
            "layer_zero_determinant": str(layer_zero_crt_determinant),
        },
        "certified_w0_controllable_block": {
            "target": {
                "controlled_conditions": 2,
                "rows": [
                    "G(0): the u^0 edge equation A'*D-B*C'=1/5 "
                    "at source transverse order 1",
                    "G'(0): the [u^5]J_4=0 target-return equation "
                    "at source transverse order 6",
                ],
                "triangular_control": (
                    "the exact edge witness controls order 1 and the "
                    "9-by-12 unit Bezout block controls order 6"
                ),
                "original_obstruction_rank": 14,
                "remaining_rank": 12,
            },
            "layer_zero": {
                "controlled_conditions": 2,
                "rows": ["H'(0)=0 at v^5", "H''(0)=0 at v^10"],
                "original_obstruction_rank": 14,
                "remaining_rank": 12,
            },
        },
        "smallest_global_residual_module": {
            "coordinate_module_over_B_rank": 24,
            "rank_over_rank_two_candidate_algebra": 24,
            "underlying_rank_over_B_after_candidate_base_change": 48,
            "coordinate_module_target_summands": [
                "B[epsilon_1]/(epsilon_1^5)",
                "B[epsilon_w0]/(epsilon_w0^5)",
                "B^2 (two triangular residues)",
            ],
            "coordinate_module_layer_zero_summands": [
                "B[eta_1]/(eta_1^6)",
                "B[eta_w0]/(eta_w0^6)",
            ],
            "candidate_base_change": (
                "tensor the displayed rank-24 B-module with "
                "A8=B[y]/(27*y^2-9*y+1)"
            ),
            "next_fitting_problem": (
                "substitute the triangular Laurent solutions into these 24 "
                "global coordinates; a unit Fitting ideal excludes k=8, "
                "while a proper ideal is only the exact pre-lower-tail "
                "residual and must still satisfy the remaining Laurent rows"
            ),
        },
    }


def fixed_endpoint_hermite_elimination_audit(
    p_bands: dict[int, Band],
    q_bands: dict[int, Band],
) -> dict[str, object]:
    """Eliminate the complete fixed-endpoint Hermite block at ``w=1``.

    The rank-24 global quotient splits off five target and six layer-zero
    coordinates at the fixed point ``w=1``.  One target coordinate is the
    already normalized terminal identity.  The remaining ten coordinates
    are affine-linear in ten endpoint Taylor coefficients from the exact
    ``P_3`` and ``P_-1`` source bands, with constant determinant 75000.

    The chosen Taylor directions have degree-seven followers divisible by
    ``w^2``.  They therefore preserve the already eliminated value/first-jet
    block at ``w=0`` while leaving the required endpoint coefficients
    unchanged.  This is an exact ideal elimination, not a raw dimension
    count.
    """

    delta, w = sp.symbols("delta w")
    local_w = 1 + delta

    expected_band_data = {
        ("P", 3): (2, 1, 11),
        ("P", -1): (1, 0, 12),
        ("Q", 1): (4, 0, 20),
        ("Q", 5): (0, 2, 18),
        ("Q", -3): (3, 0, 20),
    }
    for (side, layer), expected in expected_band_data.items():
        bands = p_bands if side == "P" else q_bands
        if band_factor_data(bands[layer]) != expected:
            raise AssertionError("a fixed-endpoint pivot band changed")

    def target_pair_expression(
        p_layer: int,
        p_k: sp.Expr,
        q_layer: int,
        q_k: sp.Expr,
    ) -> sp.Expr:
        if p_layer + q_layer != 4:
            raise AssertionError("a target endpoint pair left layer four")
        p_u, p_vanishing, _ = band_factor_data(p_bands[p_layer])
        q_u, q_vanishing, _ = band_factor_data(q_bands[q_layer])
        shift = (p_u + q_u - 1) // 5
        p_factor = delta**p_vanishing * p_k
        q_factor = delta**q_vanishing * q_k
        return sp.expand(
            local_w**shift
            * (
                p_layer
                * p_factor
                * (
                    q_u * q_factor
                    + 5 * local_w * sp.diff(q_factor, delta)
                )
                - q_layer
                * (
                    p_u * p_factor
                    + 5 * local_w * sp.diff(p_factor, delta)
                )
                * q_factor
            )
        )

    def layer_zero_product_expression(
        p_layer: int,
        p_k: sp.Expr,
        q_layer: int,
        q_k: sp.Expr,
    ) -> sp.Expr:
        if p_layer + q_layer != 0:
            raise AssertionError("a layer-zero endpoint pair left layer zero")
        p_u, p_vanishing, _ = band_factor_data(p_bands[p_layer])
        q_u, q_vanishing, _ = band_factor_data(q_bands[q_layer])
        shift = (p_u + q_u) // 5
        return sp.expand(
            p_layer
            * local_w**shift
            * delta ** (p_vanishing + q_vanishing)
            * p_k
            * q_k
        )

    def delta_order(expression: sp.Expr) -> int:
        polynomial = sp.Poly(sp.expand(expression), delta)
        nonzero_degrees = [
            degree[0]
            for degree, coefficient in polynomial.terms()
            if coefficient
        ]
        return min(nonzero_degrees) if nonzero_degrees else 10**9

    # Determine exactly which band pairs and local Taylor coefficients can
    # reach target orders 0..4 or layer-zero orders 0..5.  The zero-layer
    # pair ell=0 has weight zero in H=sum ell*P_ell*Q_-ell and is omitted.
    target_pair_records: list[dict[str, int]] = []
    layer_zero_pair_records: list[dict[str, int]] = []
    jet_requirements: dict[str, dict[int, int]] = {"P": {}, "Q": {}}

    def require(side: str, layer: int, maximum_jet: int) -> None:
        jet_requirements[side][layer] = max(
            jet_requirements[side].get(layer, -1), maximum_jet
        )

    for p_layer in range(-21, 16):
        q_layer = 4 - p_layer
        _, p_vanishing, _ = band_factor_data(p_bands[p_layer])
        _, q_vanishing, _ = band_factor_data(q_bands[q_layer])
        if p_vanishing + q_vanishing > 5:
            continue
        test_orders = [
            delta_order(
                target_pair_expression(
                    p_layer,
                    local_w**p_power,
                    q_layer,
                    local_w**q_power,
                )
            )
            for p_power in range(2)
            for q_power in range(2)
        ]
        minimum_order = min(test_orders)
        if minimum_order >= 5:
            raise AssertionError("a target dependency pair became inactive")
        maximum_jet = 5 - p_vanishing - q_vanishing
        require("P", p_layer, maximum_jet)
        require("Q", q_layer, maximum_jet)
        target_pair_records.append(
            {
                "P_layer": p_layer,
                "Q_layer": q_layer,
                "P_vanishing": p_vanishing,
                "Q_vanishing": q_vanishing,
                "minimum_target_order": minimum_order,
                "maximum_required_K_taylor_jet": maximum_jet,
            }
        )

    for p_layer in range(-25, 16):
        if p_layer == 0:
            continue
        q_layer = -p_layer
        _, p_vanishing, _ = band_factor_data(p_bands[p_layer])
        _, q_vanishing, _ = band_factor_data(q_bands[q_layer])
        total_vanishing = p_vanishing + q_vanishing
        if total_vanishing > 5:
            continue
        maximum_jet = 5 - total_vanishing
        require("P", p_layer, maximum_jet)
        require("Q", q_layer, maximum_jet)
        layer_zero_pair_records.append(
            {
                "P_layer": p_layer,
                "Q_layer": q_layer,
                "product_vanishing": total_vanishing,
                "maximum_required_K_taylor_jet": maximum_jet,
            }
        )

    if len(target_pair_records) != 22:
        raise AssertionError("the fixed target dependency-pair count changed")
    if len(layer_zero_pair_records) != 25:
        raise AssertionError("the fixed layer-zero dependency-pair count changed")
    local_coordinate_counts = {
        side: sum(maximum_jet + 1 for maximum_jet in requirements.values())
        for side, requirements in jet_requirements.items()
    }
    if local_coordinate_counts != {"P": 81, "Q": 81}:
        raise AssertionError("the fixed-endpoint Taylor-coordinate count changed")
    fixed_local_normalizations = {
        "P_3_K(1)": "1/5",
        "Q_1_K(1)": "-1",
        "Q_13_K(1)": "-3/5^5",
    }
    free_local_coordinate_count = (
        sum(local_coordinate_counts.values()) - len(fixed_local_normalizations)
    )
    if free_local_coordinate_count != 159:
        raise AssertionError("the fixed-endpoint free-coordinate count changed")

    # Derive the triangular Laurent closure from the endpoint dependency
    # bands.  A direct P band first occurs in the top-down recurrence when it
    # meets Q_25; a direct Q band first occurs when it meets P_15.  Descend to
    # the smallest such entry layer, then close every pair on all intervening
    # equations.  This derives layer 3 and the enlarged band intervals rather
    # than inserting them as a support mask.
    direct_endpoint_layers = {
        "P": set(jet_requirements["P"]),
        "Q": set(jet_requirements["Q"]),
    }
    if direct_endpoint_layers != {
        "P": set(range(-13, 13)),
        "Q": set(range(-12, 14)),
    }:
        raise AssertionError("the direct endpoint band intervals changed")
    p_top_layer = max(p_bands)
    q_top_layer = max(q_bands)
    if (p_top_layer, q_top_layer) != (15, 25):
        raise AssertionError("the corner-derived top band layers changed")
    entry_layers = {
        "P": {
            layer: layer + q_top_layer
            for layer in sorted(direct_endpoint_layers["P"])
        },
        "Q": {
            layer: p_top_layer + layer
            for layer in sorted(direct_endpoint_layers["Q"])
        },
    }
    closure_floor = min(
        entry_layer
        for side_entries in entry_layers.values()
        for entry_layer in side_entries.values()
    )
    closure_ceiling = p_top_layer + q_top_layer
    if (closure_ceiling, closure_floor) != (40, 3):
        raise AssertionError("the fixed-endpoint Laurent closure changed")
    closure_pairs = {
        (p_layer, total_layer - p_layer)
        for total_layer in range(closure_floor, closure_ceiling + 1)
        for p_layer in p_bands
        if total_layer - p_layer in q_bands
    }
    closure_p_layers = {p_layer for p_layer, _ in closure_pairs}
    closure_q_layers = {q_layer for _, q_layer in closure_pairs}
    if closure_p_layers != set(range(-22, 16)):
        raise AssertionError("the endpoint closure P interval changed")
    if closure_q_layers != set(range(-12, 26)):
        raise AssertionError("the endpoint closure Q interval changed")

    # Use generic local Taylor polynomials for the only occurrences of the
    # two selected P bands in the target and layer-zero endpoint blocks.
    p_three = sp.symbols("ep_p3_0:8")
    p_minus_one = sp.symbols("ep_pm1_0:8")
    q_one = sp.symbols("ep_q1_0:8")
    q_five = sp.symbols("ep_q5_0:8")
    q_minus_three = sp.symbols("ep_qm3_0:8")

    def taylor(coefficients: tuple[sp.Symbol, ...]) -> sp.Expr:
        return sum(
            coefficient * delta**order
            for order, coefficient in enumerate(coefficients)
        )

    selected_target = sp.expand(
        target_pair_expression(3, taylor(p_three), 1, taylor(q_one))
        + target_pair_expression(
            -1, taylor(p_minus_one), 5, taylor(q_five)
        )
    )
    selected_layer_zero = sp.expand(
        layer_zero_product_expression(
            -1, taylor(p_minus_one), 1, taylor(q_one)
        )
        + layer_zero_product_expression(
            3, taylor(p_three), -3, taylor(q_minus_three)
        )
    )
    endpoint_normalization = {
        p_three[0]: sp.Rational(1, 5),
        q_one[0]: sp.Integer(-1),
    }
    target_leading_value = sp.expand(selected_target).coeff(delta, 0).subs(
        endpoint_normalization
    )
    if target_leading_value != 1:
        raise AssertionError("the normalized target endpoint identity changed")
    if sp.expand(selected_layer_zero.subs(delta, -1)) != 0:
        raise AssertionError("the endpoint pivots unexpectedly changed H(0)")

    endpoint_rows = [
        sp.expand(selected_target).coeff(delta, order)
        for order in range(1, 5)
    ] + [
        sp.expand(selected_layer_zero).coeff(delta, order)
        for order in range(6)
    ]
    pivot_variables = (
        *p_three[1:5],
        *p_minus_one[:6],
    )
    for row in endpoint_rows:
        if sp.Poly(row, *pivot_variables).total_degree() > 1:
            raise AssertionError("the endpoint pivot block lost affinity")
    endpoint_matrix = sp.Matrix(
        [
            [
                sp.diff(row, variable).subs(endpoint_normalization)
                for variable in pivot_variables
            ]
            for row in endpoint_rows
        ]
    )
    endpoint_determinant = sp.factor(endpoint_matrix.det())
    if endpoint_matrix.shape != (10, 10):
        raise AssertionError("the fixed-endpoint matrix shape changed")
    if endpoint_determinant != 75000:
        raise AssertionError("the fixed-endpoint unit determinant changed")
    if endpoint_determinant != 2**3 * 3 * 5**5:
        raise AssertionError("the endpoint determinant factorization changed")

    # A Taylor pivot delta^j can be lifted without changing its coefficients
    # through delta^5 and while vanishing to order two at w=0 (delta=-1):
    #
    #   L_j=delta^j+(j-7)(-1)^j delta^6+(j-6)(-1)^j delta^7.
    #
    # The exact P_3 and P_-1 K-polynomial degree bounds are 11 and 12, so
    # every required follower lies in the original polynomial source band.
    follower_lifts: list[dict[str, object]] = []
    for side, layer, orders in (
        ("P", 3, range(1, 5)),
        ("P", -1, range(6)),
    ):
        _, _, free_degree = band_factor_data(p_bands[layer])
        if free_degree < 7:
            raise AssertionError("an endpoint follower left its source band")
        for order in orders:
            lift = sp.expand(
                delta**order
                + (order - 7) * (-1) ** order * delta**6
                + (order - 6) * (-1) ** order * delta**7
            )
            if sp.expand(lift.subs(delta, -1)) != 0:
                raise AssertionError("an endpoint follower changed a w=0 value")
            if sp.expand(sp.diff(lift, delta).subs(delta, -1)) != 0:
                raise AssertionError("an endpoint follower changed a w=0 jet")
            for lower_order in range(6):
                expected = int(lower_order == order)
                if sp.Poly(lift, delta).nth(lower_order) != expected:
                    raise AssertionError("an endpoint follower changed a local pivot")
            quotient, remainder = sp.div(lift, local_w**2, delta)
            if remainder != 0 or sp.degree(quotient, delta) != 5:
                raise AssertionError("an endpoint follower lost w^2 divisibility")
            follower_lifts.append(
                {
                    "side": side,
                    "source_layer": layer,
                    "endpoint_taylor_order": order,
                    "K_variation": str(lift),
                    "vanishing_at_w0": 2,
                    "source_band_K_degree_bound": free_degree,
                }
            )
    if len(follower_lifts) != 10:
        raise AssertionError("the endpoint follower count changed")

    # The toroidal terminal coordinate explains why the fixed point w=1 and
    # the residue-cover five-cycle should not be conflated.  With
    # s=t^17*z^12 and q=t^10*z^7, one has t=s^-7*q^12 and z=s^10*q^-17.
    residue_s, normal_q = sp.symbols("residue_s normal_q", nonzero=True)
    terminal_t = residue_s**-7 * normal_q**12
    terminal_z = residue_s**10 * normal_q**-17
    terminal_p = sp.factor(
        terminal_t**4 * terminal_z**3
        + terminal_t**21 * terminal_z**15
    )
    terminal_minus_q = sp.factor(
        terminal_t * terminal_z
        + 3 * terminal_t**18 * terminal_z**13
        + sp.Rational(9, 5) * terminal_t**35 * terminal_z**25
    )
    expected_terminal_p = (
        normal_q**-3 * residue_s**2 * (1 + residue_s)
    )
    expected_terminal_minus_q = (
        normal_q**-5
        * residue_s**3
        * (1 + 3 * residue_s + sp.Rational(9, 5) * residue_s**2)
    )
    if sp.factor(terminal_p - expected_terminal_p) != 0:
        raise AssertionError("the terminal P toroidal conversion changed")
    if sp.factor(terminal_minus_q - expected_terminal_minus_q) != 0:
        raise AssertionError("the terminal Q toroidal conversion changed")
    kummer_uniformizer = sp.expand((1 + terminal_t) ** 5 - 1)
    leading_normal_coefficient = sp.limit(
        kummer_uniformizer / normal_q**12,
        normal_q,
        0,
    )
    if sp.factor(leading_normal_coefficient - 5 * residue_s**-7) != 0:
        raise AssertionError("the Kummer endpoint normal order changed")

    return {
        "fixed_endpoint": "w=1, delta=w-1",
        "terminal_toroidal_interpretation": {
            "residue_coordinate": "s=t^17*z^12",
            "normal_uniformizer": "q=t^10*z^7",
            "inverse_monomials": {
                "t": "s^-7*q^12",
                "z": "s^10*q^-17",
            },
            "terminal_P": "q^-3*s^2*(1+s)",
            "terminal_minus_Q": "q^-5*s^3*(1+3*s+(9/5)*s^2)",
            "w_minus_one_leading_term": "5*s^-7*q^12",
            "consequence": (
                "w=1 is the terminal divisor's transverse endpoint; its "
                "Kummer return has normal order 12, whereas the residue "
                "five-cycle acts in the independent s-cover"
            ),
        },
        "dependency_cone": {
            "target_pair_count": len(target_pair_records),
            "layer_zero_pair_count": len(layer_zero_pair_records),
            "target_pairs": target_pair_records,
            "layer_zero_pairs": layer_zero_pair_records,
            "required_K_taylor_jets_by_band": {
                side: {
                    str(layer): maximum_jet
                    for layer, maximum_jet in sorted(requirements.items())
                }
                for side, requirements in jet_requirements.items()
            },
            "local_coordinate_counts": local_coordinate_counts,
            "fixed_local_normalizations": fixed_local_normalizations,
            "free_local_coordinate_count": free_local_coordinate_count,
            "triangular_laurent_closure": {
                "direct_endpoint_band_intervals": {
                    "P": [
                        min(direct_endpoint_layers["P"]),
                        max(direct_endpoint_layers["P"]),
                    ],
                    "Q": [
                        min(direct_endpoint_layers["Q"]),
                        max(direct_endpoint_layers["Q"]),
                    ],
                },
                "corner_derived_top_bands": {
                    "P": p_top_layer,
                    "Q": q_top_layer,
                },
                "entry_layers_by_direct_band": {
                    side: {
                        str(layer): entry_layer
                        for layer, entry_layer in side_entries.items()
                    }
                    for side, side_entries in entry_layers.items()
                },
                "layers": [closure_ceiling, closure_floor],
                "target_layer_retained": 4,
                "P_band_interval": [
                    min(closure_p_layers),
                    max(closure_p_layers),
                ],
                "Q_band_interval": [
                    min(closure_q_layers),
                    max(closure_q_layers),
                ],
                "reason": (
                    "Q_-12 first appears with P_15 on layer 3; the same "
                    "descent introduces P_-22 against Q_25"
                ),
            },
        },
        "exact_endpoint_elimination": {
            "original_coordinates": {
                "target_at_w1": 5,
                "layer_zero_at_w1": 6,
            },
            "automatic_coordinate": (
                "the target delta^0 coefficient is exactly 1 from the "
                "normalized P_3/Q_1 terminal pair"
            ),
            "equations_eliminated": 10,
            "pivot_variables": [str(variable) for variable in pivot_variables],
            "matrix_shape": list(endpoint_matrix.shape),
            "unit_determinant": str(endpoint_determinant),
            "unit_determinant_factorization": "2^3*3*5^5",
            "partner_jet_independence": (
                "the determinant is independent of every Q_1, Q_5, and "
                "Q_-3 Taylor jet and of all other forcing terms"
            ),
            "w0_preserving_source_lifts": follower_lifts,
            "remaining_global_hermite_coordinate_count": 13,
            "remaining_coordinates": {
                "target_at_w0": 5,
                "target_triangular_residues": 2,
                "layer_zero_at_w0": 6,
            },
            "claim_boundary": (
                "the ten endpoint equations solve ten exact source Taylor "
                "variables; their substitutions must still be carried into "
                "the Laurent equations on layers 40 through 3 and into the "
                "remaining thirteen global coordinates"
            ),
        },
    }


def layer_zero_fitting_audit(
    p_bands: dict[int, Band],
    q_bands: dict[int, Band],
) -> dict[str, object]:
    """Integrate the first genuine post-target row on the double-R branch.

    On Laurent layer zero every pair is ``(P_ell,Q_-ell)``.  Its bracket is

        ell*(P_ell*Q_-ell' + P_ell'*Q_-ell)
          = ell*(P_ell*Q_-ell)'.

    The exact band factorizations make every product a polynomial in
    ``w=(1+t)^5``.  The layer-zero equation is therefore one finite
    polynomial membership problem.  On the movable double-root stratum its
    cokernel is the rank-fourteen module displayed below.
    """

    t, u, w, w0 = sp.symbols("t u w w0", nonzero=True)
    generic_p = sp.Function("generic_p")(t)
    generic_q = sp.Function("generic_q")(t)
    generic_layer = sp.Symbol("generic_layer")
    product_rule = sp.expand(
        generic_layer * generic_p * sp.diff(generic_q, t)
        + generic_layer * sp.diff(generic_p, t) * generic_q
        - generic_layer * sp.diff(generic_p * generic_q, t)
    )
    if product_rule != 0:
        raise AssertionError("the layer-zero first integral changed")

    pair_records: list[dict[str, object]] = []
    maximum_degree = 0
    for p_layer in range(-25, 16):
        q_layer = -p_layer
        p_u, p_vanishing, p_degree = band_factor_data(p_bands[p_layer])
        q_u, q_vanishing, q_degree = band_factor_data(q_bands[q_layer])
        if (p_u + q_u) % 5:
            raise AssertionError("a layer-zero pair lost its Kummer character")
        w_power = (p_u + q_u) // 5
        vanishing_power = p_vanishing + q_vanishing
        product_degree = w_power + vanishing_power + p_degree + q_degree
        maximum_degree = max(maximum_degree, product_degree)
        pair_records.append(
            {
                "P_layer": p_layer,
                "Q_layer": q_layer,
                "weight": p_layer,
                "product_factor": (
                    f"w^{w_power}*(w-1)^{vanishing_power}"
                ),
                "free_product_degree_bound": p_degree + q_degree,
                "total_degree_bound": product_degree,
            }
        )
    if maximum_degree != 33:
        raise AssertionError("the layer-zero first-integral degree changed")

    # Put R on its remaining nonzero double-root stratum.  Scalar units do
    # not affect the membership problem, but retaining the normalization
    # checks the exact top-band constants and the P_-25/Q_-15 follower.
    R = (w - w0) ** 2 / (25 * (1 - w0) ** 2)
    K = sp.factor(w**3 * (w - 1) ** 6 * R**3)
    K_monic = sp.expand(w**3 * (w - 1) ** 6 * (w - w0) ** 6)
    if sp.Poly(K, w).degree() != 15 or sp.Poly(K_monic, w).degree() != 15:
        raise AssertionError("the double-R layer-zero modulus changed")

    A_tail, S_tail = sp.symbols("A_tail S_tail")
    p_top = t**15 * (w - 1) ** 6 * R**3
    q_top = -sp.Rational(9, 5) * t**25 * (w - 1) ** 10 * R**5
    p_tail = t**-25 * w**5 * A_tail
    q_tail = t**-15 * w**3 * S_tail
    new_tail_integral = sp.factor(
        15 * p_top * q_tail - 25 * p_tail * q_top
    )
    residual_tail = sp.expand(
        S_tail + 3 * w**2 * (w - 1) ** 4 * R**2 * A_tail
    )
    if sp.factor(new_tail_integral - 15 * K * residual_tail) != 0:
        raise AssertionError("the layer-zero P/Q follower cancellation changed")

    # The equation dH/dt=0 permits one integration constant.  Since the new
    # residual tail is K*S with deg(S)<=18, a prior forcing H_<40 of degree
    # at most 33 is solvable precisely when
    #
    #     H_<40 in B + K*B[w]_(<=18).
    #
    # Equivalently its class in B[w]/(K) is constant.  The quotient by the
    # constants has B-rank 15-1=14.  The Hermite conditions below are an
    # explicit basis-free form of the same membership statement.
    ambient_dimension = maximum_degree + 1
    source_columns = 1 + 19
    cokernel_rank = ambient_dimension - source_columns
    if (ambient_dimension, source_columns, cokernel_rank) != (34, 20, 14):
        raise AssertionError("the layer-zero Fitting dimensions changed")
    hermite_conditions = [
        "H'(0)=H''(0)=0",
        "H(1)=H(0) and H^(j)(1)=0 for 1<=j<=5",
        "H(w0)=H(0) and H^(j)(w0)=0 for 1<=j<=5",
    ]
    hermite_condition_count = 2 + 6 + 6
    if hermite_condition_count != cokernel_rank:
        raise AssertionError("the layer-zero Hermite count changed")

    # The pure 8-step fractional-power staircase recovers E5 as only the
    # first movable-root Hermite residue.  It is not the full layer-zero
    # equation because off-grid pairs also enter H_<40.
    h, C, pi, a, b = sp.symbols("h C pi a b", nonzero=True)
    p1, p2, p3, p4, p5 = sp.symbols("p1 p2 p3 p4 p5")
    P_series = C**3 + sum(
        coefficient * h**index
        for index, coefficient in enumerate((p1, p2, p3, p4, p5), 1)
    )
    Q_series = sp.series(
        -sp.Rational(9, 5)
        * C**5
        * (P_series / C**3) ** sp.Rational(5, 3),
        h,
        0,
        6,
    ).removeO()
    P_coefficients = [C**3, p1, p2, p3, p4, p5]
    Q_coefficients = [
        sp.factor(sp.expand(Q_series).coeff(h, index))
        for index in range(6)
    ]
    pure_integral = sp.factor(
        sum(
            (15 - 8 * index)
            * P_coefficients[index]
            * Q_coefficients[5 - index]
            for index in range(6)
        )
    )
    if pure_integral != 0:
        raise AssertionError("the pure eight-step first integral changed")
    pure_new_tail = sp.factor(
        15 * C**3 * (Q_coefficients[5] + 3 * C**2 * p5)
    )
    pure_prior = sp.factor(-pure_new_tail)
    E5 = 7 * a**4 - 60 * a**2 * b + 135 * b**2
    pure_local_first_jet = sp.factor(
        pure_prior.subs(
            {C: pi**2, p1: a * pi**3, p2: b, p3: 0, p4: 0}
        )
    )
    if sp.factor(pure_local_first_jet + a * pi * E5 / 27) != 0:
        raise AssertionError("the pure staircase E5 Hermite residue changed")

    # At the unrestricted B0 level the residue module is already saturated
    # by old (not newly introduced) layer-zero pairs.  With the exact Q_1
    # witness S=(w-w0)/(w0-1), the ell=-1 pair spans
    # w*(w-w0)*B[w]_(<=12).  The ell=1 pair spans
    # w*(w-1)*B[w]_(<=30).  Their difference gives w^j for 1<=j<=13, and
    # w^14=w^13+w*(w-1)*w^12.  Hence together they span every nonconstant
    # residue class.  Earlier triangular equations may restrict this span;
    # the layer-zero row by itself cannot kill the branch.
    saturation_identity = sp.expand(
        w * (w - w0) - w * (w - 1) - (1 - w0) * w
    )
    if saturation_identity != 0:
        raise AssertionError("the raw layer-zero saturation identity changed")

    return {
        "first_integral": {
            "identity": "J_0=d/dt(sum_ell ell*P_ell*Q_-ell)",
            "w_coordinate": "w=(1+t)^5",
            "integral_degree_bound": maximum_degree,
            "band_pair_records": pair_records,
        },
        "movable_double_R": {
            "base": "B=QQ[w0,w0^-1,(w0-1)^-1]",
            "R": "(w-w0)^2/(25*(1-w0)^2)",
            "modulus_up_to_a_unit": "K=w^3*(w-1)^6*(w-w0)^6",
            "degree": 15,
            "local_lengths": {"w=0": 3, "w=1": 6, "w=w0": 6},
        },
        "new_tail_columns": {
            "P_band": -25,
            "Q_band": -15,
            "combined_residual": (
                "S_-15+3*w^2*(w-1)^4*R(w)^2*A_-25"
            ),
            "free_degree_bound": 18,
            "integral_image": "15*K*B[w]_(<=18)",
            "P_follower": "q_-15=-3*C0^2*p_-25",
        },
        "artinian_fitting_statement": {
            "membership": "H_<40 belongs to B + K*B[w]_(<=18)",
            "artinian_algebra": "A_40=B[w]/(K)",
            "artinian_length": 15,
            "obstruction_module": "A_40/(B*1)",
            "obstruction_rank_over_B": 14,
            "presentation_matrix_shape": [34, 20],
            "augmented_matrix_shape_for_one_forcing": [34, 21],
            "fitting_condition": (
                "rank([1,K,w*K,...,w^18*K,H_<40])="
                "rank([1,K,w*K,...,w^18*K])=20"
            ),
            "hermite_conditions": hermite_conditions,
            "hermite_condition_count": hermite_condition_count,
        },
        "pure_spacing_eight_projection": {
            "first_movable_root_residue": "-a*pi*E5/27",
            "E5": "7*a^4-60*a^2*b+135*b^2",
            "interpretation": (
                "E5 is the first Hermite coordinate of the pure 8-step "
                "staircase, not the full post-target equation"
            ),
        },
        "raw_B0_saturation": {
            "ell_minus_1_with_fixed_Q1": (
                "w*(w-w0)*B[w]_(<=12) up to a base unit"
            ),
            "ell_plus_1": "w*(w-1)*B[w]_(<=30)",
            "span_modulo_constants": 14,
            "constructive_identity": (
                "w*(w-w0)-w*(w-1)=(1-w0)*w; multiply through "
                "w^12, then use w^14=w^13+w*(w-1)*w^12"
            ),
            "claim_boundary": (
                "this proves saturation only before imposing the earlier "
                "triangular Keller rows; those rows must be carried into "
                "the specialized residue"
            ),
        },
    }


def build_payload() -> dict[str, object]:
    coordinate = coordinate_bracket_regression()
    # The degree and terminal-halfspace bounds give the complete B0 band
    # intervals.  The formerly emitted bands -21..15 and -11..25 are only
    # the subintervals which can contribute at or above the target layer 4.
    p_bands = {
        layer: make_band("P", 75, 3, layer)
        for layer in range(-75, 16)
    }
    q_bands = {
        layer: make_band("Q", 125, 5, layer)
        for layer in range(-125, 26)
    }
    if sum(band.source_count for band in p_bands.values()) != 706:
        raise AssertionError("the full P source-band count changed")
    if sum(band.dimension for band in p_bands.values()) != 653:
        raise AssertionError("the full P jet-reduced dimension changed")
    if sum(band.source_count for band in q_bands.values()) != 1901:
        raise AssertionError("the full Q source-band count changed")
    if sum(band.dimension for band in q_bands.values()) != 1765:
        raise AssertionError("the full Q jet-reduced dimension changed")
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
        for layer in range(40, -201, -1)
    ]
    if len(layers) != 241:
        raise AssertionError("the full Laurent layer count changed")
    upper_layers = [record for record in layers if int(record["layer"]) >= 4]
    missing_layers = [
        record
        for record in upper_layers
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
        "schema": "plane-jc.f2-75-125-character-layers.v9",
        "status": (
            "exact-full-B0-envelope-recurrence-and-corrected-top-tangent;"
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
            "full_band_intervals": {"P": [-75, 15], "Q": [-125, 25]},
            "full_bracket_layer_interval": [-200, 40],
            "target_layer": 4,
            "zero_layer_intervals": [[40, 40], [39, 5], [3, -200]],
            "zero_layer_count": 240,
        },
        "bands_for_layers_4_through_40": {
            "P": [band_record(p_bands[layer]) for layer in range(-21, 16)],
            "Q": [band_record(q_bands[layer]) for layer in range(-11, 26)],
            "source_variable_count": {
                "P": sum(p_bands[layer].source_count for layer in range(-21, 16)),
                "Q": sum(q_bands[layer].source_count for layer in range(-11, 26)),
            },
            "linear_dimension_after_terminal_jets": {
                "P": sum(p_bands[layer].dimension for layer in range(-21, 16)),
                "Q": sum(q_bands[layer].dimension for layer in range(-11, 26)),
            },
            "missing_35_layer_window_dimension": missing_parameter_count,
            "missing_window_dimension_after_five_normalizations": (
                missing_parameter_count - 5
            ),
            "missing_window_dimension_after_joint_top_band_and_other_three_normalizations": (
                missing_parameter_count - 18 + 2 - 3
            ),
        },
        "complete_B0_band_envelope": {
            "P": [band_record(p_bands[layer]) for layer in sorted(p_bands)],
            "Q": [band_record(q_bands[layer]) for layer in sorted(q_bands)],
            "source_variable_count_before_terminal_jets": {
                "P": sum(band.source_count for band in p_bands.values()),
                "Q": sum(band.source_count for band in q_bands.values()),
            },
            "dimension_after_terminal_jets": {
                "P": sum(band.dimension for band in p_bands.values()),
                "Q": sum(band.dimension for band in q_bands.values()),
                "total": sum(band.dimension for band in p_bands.values())
                + sum(band.dimension for band in q_bands.values()),
            },
            "dimension_after_five_endpoint_normalizations": 2413,
            "first_post_target_tail_bands": {
                "P_minus_22_through_minus_25": [
                    band_record(p_bands[layer]) for layer in range(-22, -26, -1)
                ],
                "Q_minus_12_through_minus_15": [
                    band_record(q_bands[layer]) for layer in range(-12, -16, -1)
                ],
            },
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
        "family_full_support_theorem": family_full_support_theorem(),
        "upper_descent_classification": upper_descent_classification(
            p_bands, q_bands
        ),
        "nonlinear_first_defect_audit": nonlinear_first_defect_audit(
            p_bands, q_bands
        ),
        "target_layer_fitting_audit": target_layer_fitting_audit(
            p_bands, q_bands
        ),
        "u_zero_edge_audit": u_zero_edge_audit(),
        "transverse_edge_completion_audit": transverse_edge_completion_audit(),
        "first_kummer_return_audit": first_kummer_return_audit(
            p_bands, q_bands
        ),
        "source_return_lift_and_global_quotient_audit": (
            source_return_lift_and_global_quotient_audit(p_bands, q_bands)
        ),
        "fixed_endpoint_hermite_elimination_audit": (
            fixed_endpoint_hermite_elimination_audit(p_bands, q_bands)
        ),
        "layer_zero_fitting_audit": layer_zero_fitting_audit(
            p_bands, q_bands
        ),
        "compressed_quadratic_system": {
            "generator_formula": (
                "(ell*s-m*r)*alpha_(ell,i)*beta_(m,k)*"
                "t^(ell+m)*(1+t)^(r+s-1)"
            ),
            "generator_digest_sha256": quadratic_digest.hexdigest(),
            "layers_40_through_minus_200": layers,
            "target_layer": 4,
            "zero_layer_intervals": [[40, 40], [39, 5], [3, -200]],
            "zero_layer_count": 240,
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
                "post_jet_support_row_upper_bound": sum(
                    int(record["post_jet_support_row_upper_bound"])
                    for record in missing_layers
                ),
            },
            "full_B0_totals": {
                "band_pairs": sum(
                    int(record["band_pair_count"]) for record in layers
                ),
                "active_binomial_generators": sum(
                    int(record["active_binomial_generators"])
                    for record in layers
                ),
                "raw_scalar_coefficient_rows": sum(
                    int(record["raw_scalar_coefficient_row_count"])
                    for record in layers
                ),
                "post_jet_support_row_upper_bound": sum(
                    int(record["post_jet_support_row_upper_bound"])
                    for record in layers
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
                "the complete lower B0 tail through layer -200, including "
                "all 204 post-target zero layers and their source bands",
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
                "the Q-band-one filter excluding the fixed u^5=1 double "
                "roots, leaving only the movable double root of R",
                "an exact normalized P3/Q1 target-jet witness on that final "
                "local branch",
                "the complete rank-fourteen target-layer cokernel, split "
                "into twelve local jets and two triangular residues",
                "the exact u=0 edge Bezout equation and a source-degree "
                "compatible witness showing that edge alone survives",
                "the all-r sparse edge section and its exact formal shear "
                "completion, whose unavoidable infinite tail rules out the "
                "literal four-term polynomial escape at transverse order two",
                "the all-r classification of polynomial second-order edge "
                "repairs and the exclusion of every quadratic transverse "
                "truncation, together with exact r=3 unit ideals excluding "
                "cubic and quartic truncations",
                "the exact first Kummer-return packet, including its fifth-"
                "binomial correction and a unit maximal minor proving that "
                "both rank-two descent-eight branches survive v^5",
                "the source-return and Bezout ranks proving that the isolated "
                "edge staircase continues to absorb every forcing through v^10",
                "exact original-source lifts for all 53 unit-minor pivots, "
                "strictly below both certified supporting edges, proving that "
                "their repair freedom is genuine within the corner envelope",
                "the all-r unit return minors and 18*r-1 strict-interior "
                "source lifts proving that the sparse edge staircase is "
                "surjective through order ten for every r>=2",
                "the three-point CRT quotient removing the certified w=0 "
                "control block and leaving a rank-24 global Hermite module at "
                "w=1 and w=w0 over the rank-two candidate algebra",
                "the exact fixed-endpoint dependency cone, involving 22 "
                "target pairs and 25 nonzero-weight layer-zero pairs, with "
                "159 free local Taylor coordinates and Laurent closure only "
                "through layers 40 down to 3",
                "the normalized leading endpoint identity and the remaining "
                "10-by-10 affine-linear endpoint block of determinant 75000, "
                "which eliminates all w=1 Hermite coordinates and leaves 13 "
                "global coordinates after exact w=0-preserving source lifts",
                "the proof that the descent-40 fifth residue is a lower-tail "
                "Fitting row, not the primitive equation E5=0",
                "the exact layer-zero first integral and its rank-fourteen "
                "Artinian/Hermite cokernel on the movable double-R branch",
            ],
            "not_proved": [
                "that the B0 envelope is the exhaustive lower Newton support",
                "the gamma-branch list for multiplicities (3,5)",
                "the nonlinear forced Fitting equations below the common top",
                "inconsistency of the resulting system or exclusion of F2",
            ],
            "next_elimination": (
                "the subsequent endpoint-reduction replay carries these ten "
                "solutions into layers 40 through 3 and the thirteen "
                "residual coordinates, eliminates the endpoint-disjoint "
                "power block through layer 29, and locates the remaining "
                "coupled Schur/Fitting boundary at descent 12 (layer 28)"
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
        f"{totals['post_jet_support_row_upper_bound']} "
        "post-jet support-row upper bound",
    )
    full = payload["compressed_quadratic_system"]["full_B0_totals"]
    print(
        "PASS: complete B0 tail has 240 zero layers through -200, "
        "2418 jet-reduced parameters, and",
        f"{full['active_binomial_generators']} exact compressed generators",
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
    print(
        "PASS: fixed Kummer double roots are excluded; the remaining "
        "earliest double-R branch has a rank-two Artinian core and passes "
        "the first exact target jet"
    )
    print(
        "PASS: the complete double-R target operator has 12 forced jets "
        "and 2 triangular Fitting residues"
    )
    print(
        "PASS: the forced u=0 target determinant extends to an exact "
        "degree-(3,5) edge Bezout witness"
    )
    print(
        "PASS: the all-r sparse edge escape completes formally but its "
        "literal polynomial truncation fails at transverse order two"
    )
    print(
        "PASS: all-r quadratic edge completions and r=3 cubic/quartic "
        "terminations are excluded"
    )
    print(
        "PASS: the exact v^5 Kummer packet has an 8-by-10 surjective "
        "return matrix over the rank-two branch; the edge recursion "
        "continues through v^10"
    )
    print(
        "PASS: the all-r return minors have 18*r-1 strict-interior source "
        "lifts; at r=3 quotienting w=0 leaves a rank-24 Hermite module"
    )
    print(
        "PASS: the w=1 endpoint block is one normalized identity plus a "
        "10-by-10 unit elimination; 13 global Hermite coordinates remain"
    )
    print(
        "PASS: later first-defect spacings 9..90 remain explicitly "
        "classified by their pre-target row count"
    )
    print(
        "PASS: layer zero has an exact length-15 Artinian first integral "
        "and a rank-14 Fitting residue modulo constants"
    )
    print("PASS: the artifact remains explicitly short of an F2 exclusion")


if __name__ == "__main__":
    main()
