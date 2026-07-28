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
    }


def upper_descent_classification() -> dict[str, object]:
    """Classify the first five zero layers after the common-power top."""

    A, A_prime, U, U_prime, V, V_prime = sp.symbols(
        "A A_prime U U_prime V V_prime"
    )
    q_scale = -sp.Rational(9, 5)
    descent = sp.symbols("descent")
    p = A**2 * U
    p_prime = 2 * A * A_prime * U + A**2 * U_prime
    q = q_scale * A**4 * V
    q_prime = q_scale * (4 * A**3 * A_prime * V + A**4 * V_prime)
    p_top = A**3
    p_top_prime = 3 * A**2 * A_prime
    q_top = q_scale * A**5
    q_top_prime = 5 * q_scale * A**4 * A_prime
    bracket = sp.expand(
        (15 - descent) * p * q_top_prime
        - 25 * p_prime * q_top
        + 15 * p_top * q_prime
        - (25 - descent) * p_top_prime * q
    )
    W = 5 * U - 3 * V
    W_prime = 5 * U_prime - 3 * V_prime
    reduced = sp.expand(
        q_scale
        * A**6
        * ((5 - descent) * A_prime * W - 5 * A * W_prime)
    )
    if sp.expand(bracket - reduced) != 0:
        raise AssertionError("the upper-descent operator reduction failed")

    rows: list[dict[str, object]] = []
    for descent_value in range(1, 5):
        root_parameters = 2 if descent_value <= 2 else 3
        vanishing_power = 2 if descent_value <= 2 else 1
        required_order = (6, 5, 3, 2)[descent_value - 1]
        constructed_order = 5 - descent_value + vanishing_power
        if constructed_order != required_order:
            raise AssertionError("the upper root-band terminal order changed")
        if root_parameters + vanishing_power - 1 != 3:
            raise AssertionError("the upper root-band degree bound changed")
        layer = 40 - descent_value
        rows.append(
            {
                "layer": layer,
                "descent": descent_value,
                "homogeneous_ode": (
                    f"{5-descent_value}*C0'*W-5*C0*W'=0"
                ),
                "valuation_at_t_zero": (
                    f"ord(W)=7*{5-descent_value}/5"
                ),
                "consequence": "W=0 because the required valuation is nonintegral",
                "new_root_band": (
                    f"D_{descent_value}=t^{5-descent_value}*u^{descent_value}*"
                    f"(u^5-1)^{vanishing_power}*S_{descent_value}(u^5)"
                ),
                "new_root_band_dimension": root_parameters,
            }
        )

    # At descent five the homogeneous equation is W'=0.  Its constant mode
    # is represented by C0^2 in P band 10 and commutes with Q_top=C0^5.
    extra_mode = sp.expand(
        10 * A**2 * q_top_prime
        - 25 * (2 * A * A_prime) * q_top
    )
    if extra_mode != 0:
        raise AssertionError("the layer-35 C0^2 mode is not closed")
    rows.append(
        {
            "layer": 35,
            "descent": 5,
            "homogeneous_ode": "-5*C0*W'=0",
            "valuation_at_t_zero": "ord(W)=0",
            "consequence": (
                "W is constant; besides a five-parameter root correction "
                "D_5 in k[u^5] of degree at most 20, one C0^2 mode survives"
            ),
            "new_root_band": "D_5=S_5(u^5), deg(S_5)<=4",
            "new_root_band_dimension": 5,
            "non_common_power_mode": "lambda*C0^2 in P band 10",
            "non_common_power_mode_dimension": 1,
        }
    )
    return {
        "operator_reduction": {
            "substitution": "p=C0^2*U, q=(-9/5)*C0^4*V",
            "W": "5*U-3*V",
            "descent_r_equation": (
                "(5-r)*C0'*W-5*C0*W'=known lower-layer forcing"
            ),
            "top_root_order": "ord_t(C0)=7",
        },
        "classified_layers": rows,
        "outcome": (
            "layers 39 through 36 force common-root continuation; layer 35 "
            "is the first genuine mode and contains lambda*C0^2"
        ),
        "remaining_unclassified_zero_layers": [34, 5],
        "remaining_unclassified_zero_layer_count": 30,
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
        "schema": "plane-jc.f2-75-125-character-layers.v1",
        "status": "exact-B0-envelope-and-recurrence-not-exhaustive-normal-form",
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
        "upper_descent_classification": upper_descent_classification(),
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
                "the exact upper descent through layer 35, including the "
                "first surviving C0^2 mode",
            ],
            "not_proved": [
                "that the B0 envelope is the exhaustive lower Newton support",
                "the gamma-branch list for multiplicities (3,5)",
                "the lower-band coefficient relations below layer 35",
                "inconsistency of the resulting system or exclusion of F2",
            ],
            "next_elimination": (
                "propagate the layer-35 C0^2 mode into layer 34 on the "
                "quadratic-R root strata; pivot if its forced de Rham class "
                "does not yield an early inconsistency"
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
    print("PASS: layers 39..36 continue the common root; layer 35 has a C0^2 mode")
    print("PASS: the artifact remains explicitly short of an F2 exclusion")


if __name__ == "__main__":
    main()
